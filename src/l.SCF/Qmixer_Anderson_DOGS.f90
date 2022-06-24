! copyright info:
!
!                             @Copyright 2022
!                           Fireball Committee
! Hong Kong Quantum AI Laboratory, Ltd. - James P. Lewis, Chair
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek
! Arizona State University - Otto F. Sankey

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! California Institute of Technology - Brandon Keith
! Czech Institute of Physics - Prokop Hapala
! Czech Institute of Physics - Vladimír Zobač
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Synfuels China Technology Co., Ltd. - Pengju Ren
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman
!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! Qmixer.f90
! Function Description
! ===========================================================================
!>      The input charges are passed into this routine. New charges are
!> calculated and mixed in with the old input charges.  The resulting charges
!> become output charges.
!
! The mixing is done using the Anderson mixing algorithm. Details can be found
! in the following:
!
!       Donald G. Anderson, J. Assoc. Computing Machinery, 12, 547 (1965)
!       Generalized by:  V. Eyert, J. Comp. Phys. 124, 271 (1996)
!       Using Eq. 7.7, 5.31 of Eyert, we set all betas to be equal for
! simplicity. Our routine computes input vector for the next iteration
! and the order is equal to m + 1 in Eyert. This is the same as the best
! Broyden method (see Eyert).
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Qmixer (t, iscf_iteration, sigma)

! /SYSTEM
        use M_configuraciones

! /SCF
        use M_charges

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iscf_iteration

        ! This is the pointer to the current structure
        type(T_structure), target :: t

! Output
        real, intent (inout) :: sigma

! Local Parameters and Data Declaration
! ============================================================================
        integer, parameter :: max_order = 6 ! order of iterated polynomial

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                     !< counter over the atoms
        integer iloop, jloop              !< counter for amatrix
        integer imix
        integer in1                       !< species number for iatom
        integer issh                      !< counter over shells
        integer istep                     !< short-cut notation
        integer logfile                   !< writing to which unit
        integer mix_order                 !< in case iscf .lt. max_order

        real dqrms, dqmax                 !< rms and max of charge differences
        real Qin, Qout, Qneutral          !< input and output charges
        real renorm
        real zcheck, ztotal_out

! mixing charges arrays
        real, allocatable :: Qinmixer (:)
        real, allocatable :: Qoutmixer (:)

! Needed for evaluating optimal Qout
        double precision, allocatable :: F_dot_delF (:)     ! <delF|F> in Eq. 5.31

        ! checks to see if structure has changed
        type(T_structure), pointer, save :: current

        real, allocatable, save :: delF (:, :)     ! difference between Fv's
        real, allocatable, save :: delX (:, :)     ! difference between Xv's

        ! difference between old and in charges
        real, allocatable, save :: Fv (:, :)

        real, allocatable, save :: sigma_saved (:) ! old sigma values
        real, allocatable, save :: Xv (:, :)! input charges

        double precision, allocatable, save :: amatrix (:, :)

! Needed for the dsyev LAPACK call
        integer info                   ! error information
        integer lwork                  ! size of the working array

        integer, allocatable :: ipiv (:)

        double precision, allocatable :: work (:)   ! working vector

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Check to see if the structure has changed
	    if (.not. associated (current)) then
          current => t
        else if (.not. associated(current, target=t) .and. allocated (Fv)) then
          deallocate (Fv)
          deallocate (Xv)
          deallocate (delF)
          deallocate (delX)
          deallocate (sigma_saved)
          current => t
        end if

! Initialize logfile
        logfile = t%logfile

! Check to see if input charges and output charges are within tolerance.
! If they are within tolerance, then perform a motional time step. If
! they are not within tolerence, then perform another iteration to make
! charges self-consistent.  This also gets the size of imix.
        dqrms = 0.0d0
        dqmax = -99.0d0
        imix = 0
        do iatom = 1, t%natoms
          in1 = t%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            imix = imix + 1
            Qin = t%atom(iatom)%shell(issh)%Qin
            Qout = t%atom(iatom)%shell(issh)%Qout
            dqmax = max(abs(Qin - Qout), dqmax)
            dqrms = dqrms + (Qin - Qout)**2
          end do
        end do
        dqrms = sqrt(dqrms)/(2*t%natoms)
        write (logfile,*)
        write (logfile,100) dqrms
        write (logfile,101) dqmax

! ===========================================================================
!                              Anderson mixing
! ===========================================================================
! What order interpolation do we use this time
        mix_order = min(iscf_iteration, max_order)
        allocate (Qinmixer (imix))
        allocate (Qoutmixer (imix))

        if (.not. allocated (Fv)) then
          allocate (Fv (imix, max_scf_iterations_set)); Fv = 0.0d0
          allocate (Xv (imix, max_scf_iterations_set)); Xv = 0.0d0
          allocate (delF (imix, max_scf_iterations_set)); delF = 0.0d0
          allocate (delX (imix, max_scf_iterations_set)); delX = 0.0d0
          allocate (sigma_saved (max_scf_iterations_set)); sigma_saved = 0.0d0
        end if

! Store all the charges into one dimensional arrays for easier manipulation.
        imix = 0
        do iatom = 1, t%natoms
          in1 = t%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            imix = imix + 1
            Qinmixer(imix) = t%atom(iatom)%shell(issh)%Qin
            Qoutmixer(imix) = t%atom(iatom)%shell(issh)%Qout
          end do
        end do

! Mix Qinmixer and Qoutmixer to get a NEW Qinmixer, which is the NEW input for
! the next run.
        Xv(:, iscf_iteration) = Qinmixer(:)
        Fv(:, iscf_iteration) = Qoutmixer(:) - Qinmixer(:)

! Calculate the new sigma
        sigma = dot_product(Fv(:,iscf_iteration),Fv(:,iscf_iteration))
        sigma = sigma/imix

! Initially - only perform simple extrapolation
        if (mix_order .eq. 1 .or. sigma .lt. scf_tolerance_set) then
          write (logfile,*) ' Doing simple mixing! '
          Qinmixer(:) = Qinmixer(:) + beta_set*Fv(:,iscf_iteration)
          sigma_saved(iscf_iteration) = sigma
        else

! Evaluate new terms
          sigma_saved(iscf_iteration) = sigma
          delF (:,iscf_iteration - 1) =                                       &
     &                           Fv(:,iscf_iteration) - Fv(:,iscf_iteration - 1)
          delX (:,iscf_iteration - 1) =                                       &
     &                           Xv(:,iscf_iteration) - Xv(:,iscf_iteration - 1)

! Make sure step with lowest sigma value is always used (i.e. not lost)
          istep = iscf_iteration - max_order
          if (iscf_iteration .gt. max_order .and. max_order .ge. 6) then
            if (sigma_saved(istep) .lt.                                       &
      &         minval(sigma_saved(istep + 1:iscf_iteration))) then

!             throw away second oldest step instead
              sigma_saved(istep + 1) = sigma_saved(istep)
              Fv(:,istep + 1) = Fv(:,istep)
              Xv(:,istep + 1) = Xv(:,istep)
              delX(:,istep + 1) = Xv(:,istep + 2) - Xv(:,istep + 1)
              delF(:,istep + 1) = Fv(:,istep + 2) - Fv(:,istep + 1)
            end if
          end if

! Build a_matrix Eq. 5.17
          info = 99
          do while (info .ne. 0)
            istep = iscf_iteration - mix_order
            allocate (amatrix(istep + 1:iscf_iteration - 1, istep + 1:iscf_iteration - 1))
            do iloop = istep + 1, iscf_iteration - 1
              do jloop = istep + 1, iscf_iteration - 1
                amatrix(iloop,jloop) = dot_product(delF(:,iloop),delF(:,jloop))
              end do
            end do

! Build delF_F Eq. 5.31
            if (allocated (F_dot_delF)) deallocate (F_dot_delF)
            allocate (F_dot_delF(istep + 1:iscf_iteration - 1))
            do iloop = istep + 1, iscf_iteration - 1
              F_dot_delF(iloop) = dot_product(delF(:,iloop), Fv(:,iscf_iteration))
            end do

! Solve for gammas Eq. 5.31, 7.4 (We move a-inverse to other side: a * gamma = <|>)
            ! first call dsysv with lwork = -1
            lwork = 1
            allocate (work(lwork))
            allocate (ipiv(mix_order - 1)); ipiv = 0
            call dsysv('U', mix_order - 1, 1, amatrix, mix_order - 1,          &
     &                 ipiv, F_dot_delF, mix_order - 1, work, -1, info)

            ! allocate according to result of first call
            lwork = work(1)
            deallocate (work)
            allocate (work(lwork))
            call dsysv('U', mix_order - 1, 1, amatrix, mix_order - 1,          &
     &                 ipiv, F_dot_delF, mix_order - 1, work, lwork, info)
            do iloop = istep + 1, iscf_iteration - 1
            end do

! If there is an error, then just use the simple mixing
            if (info .ne. 0) then
              write (logfile,*) ' Error in Anderson, info = ', info
              if (mix_order .le. 2) stop
              mix_order = mix_order - 1
            end if
            deallocate (work, ipiv, amatrix)
          end do

! Generate new guess at charges Eq. 7.7  (F_dot_delF is now gamma)
          Qinmixer(:) = Qinmixer(:) + beta_set*Fv(:,iscf_iteration)  ! First-order term
          istep = iscf_iteration - mix_order
          do iloop = istep + 1, iscf_iteration - 1
            Qinmixer(:) =                                                        &
     &        Qinmixer(:) - F_dot_delF(iloop)*(delX(:,iloop) + beta_set*delF(:,iloop))
          end do
          deallocate (F_dot_delF)
        end if

! ===========================================================================
!                                 end mixing
! ===========================================================================
! Check the total charge
        ztotal_out = 0
        imix = 0
        do iatom = 1, t%natoms
          in1 = t%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            imix = imix + 1
            t%atom(iatom)%shell(issh)%Qout = t%atom(iatom)%shell(issh)%Qin
            t%atom(iatom)%shell(issh)%Qin = Qinmixer(imix)
            ztotal_out = ztotal_out + t%atom(iatom)%shell(issh)%Qin
          end do
        end do
        renorm = (ztotal_out - t%ztot)/imix
        write (logfile,*)
        write (logfile,201) renorm

! Reset new charges, Qin
        zcheck = 0.0d0
        do iatom = 1, t%natoms
          in1 = t%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            Qneutral = species(in1)%shell(issh)%Qneutral

            t%atom(iatom)%shell(issh)%Qin =                                  &
     &        t%atom(iatom)%shell(issh)%Qin - renorm
            Qin = t%atom(iatom)%shell(issh)%Qin
            t%atom(iatom)%shell(issh)%dQ = Qin - Qneutral

            zcheck = zcheck + t%atom(iatom)%shell(issh)%Qin
          end do
        end do

! Write out results
        write (logfile,*) ' (Before renormalization) zouttot = ', ztotal_out
        write (logfile,*) ' (After  renormalization)  zcheck = ', zcheck
        write (logfile,*) ' (What it must be)           ztot = ', t%ztot
        write (logfile,*)
        write (logfile,202) sigma, scf_tolerance_set, iscf_iteration

! Deallocate Arrays
! ===========================================================================
        deallocate (Qinmixer, Qoutmixer)

! Format Statements
! ===========================================================================
100     format (2x, ' Deviation (rms) of input/output charges = ', f9.4)
101     format (2x, ' Deviation (max) of input/output charges = ', f9.4)
201     format (2x, ' Renormalization of Qin:   renorm = ', 4x, f14.8)
202     format (' =====> sigma = ', e14.7,                                   &
     &          ' Must be less than', e14.7, ' SCF step = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine Qmixer


