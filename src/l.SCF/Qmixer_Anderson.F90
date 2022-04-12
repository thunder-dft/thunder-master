#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

! copyright info:
!
!                             @Copyright 2013
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
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
#ifdef DOGS
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
#else
!>      This is a dummy routine. There is no charge mixing for Harris.
#endif
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Qmixer (t, iscf_iteration, sigma)
#ifdef DOGS
        use M_charges
#ifdef GRID
        use M_density_matrix
        use M_grid
#endif
#endif
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iscf_iteration

        ! This is the pointer to the current structure
        type(T_structure), target :: t

! Output
        real, intent (inout) :: sigma

#ifdef DOGS
! Local Parameters and Data Declaration
! ===========================================================================
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
#ifdef GRID
        integer jatom, ineigh
        integer in2
        integer norb_mu, norb_nu         !< size of the block for the pair
        integer num_neigh                !< number of neighbors
        integer inu,imu
#else
        real dqrms, dqmax                 !< rms and max of charge differences
        real Qin, Qout, Qneutral          !< input and output charges
        real renorm
#endif
        real zcheck, ztotal_out

! mixing charges arrays
        real, allocatable, save :: Qinmixer (:)
        real, allocatable, save :: Qoutmixer (:)

! Needed for evaluating optimal Qout
        double precision, allocatable :: F_dot_delF (:)     ! <delF|F> in Eq. 5.31

        ! checks to see if structure has changed
        type(T_structure), pointer, save :: current

        real, allocatable, save :: delF (:, :)    ! difference between Fv's
        real, allocatable, save :: delX (:, :)    ! difference between Xv's

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

#ifdef GRID
        character (len = 25) :: slogfile

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors_old
        type(T_assemble_neighbors), pointer :: pdenmat_old
        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
#endif
#endif

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================

#ifdef DOGS
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

#ifdef GRID
        inpfile = t%inpfile

! Get the size of imix.
        imix = 0
        do iatom = 1, t%natoms
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = t%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              do inu = 1, norb_nu
                imix = imix + 1
              end do ! inu
            end do ! imu
          end do ! End loop over neighbors
        end do ! End loop over atoms
#else

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
            imix = imix +1
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

#endif

! ===========================================================================
!                              Anderson mixing
! ===========================================================================
! What order interpolation do we use this time
        mix_order = min(iscf_iteration, max_order)
        allocate (Qinmixer (imix))
        allocate (Qoutmixer (imix))

        if (.not. allocated (Fv))then
          allocate (Fv (imix, max_scf_iterations))
          allocate (Xv (imix, max_scf_iterations))
          allocate (delF (imix, max_scf_iterations))
          allocate (delX (imix, max_scf_iterations))
          allocate (sigma_saved (max_scf_iterations))
        end if

! Store all the charges into one dimensional arrays for easier manipulation.
        imix = 0
        do iatom = 1, t%natoms
#ifdef GRID
          ! cut some lengthy notation
          pdenmat=>t%denmat(iatom)
          pdenmat_old=>t%denmat_old(iatom)
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = t%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over neighbors
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            pRho_neighbors_old=>pdenmat_old%neighbors(ineigh)
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              do inu = 1, norb_nu
                imix = imix + 1
                Qinmixer(imix)  = pRho_neighbors_old%block(imu,inu)
                Qoutmixer(imix) = pRho_neighbors%block(imu,inu)
              end do
            end do ! end loop over matrix elements
          end do ! end loop over neighbors
#else
          in1 = t%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            imix = imix + 1
            Qinmixer(imix) = t%atom(iatom)%shell(issh)%Qin
            Qoutmixer(imix) = t%atom(iatom)%shell(issh)%Qout
          end do
#endif
        end do ! end loop over atoms

! Mix Qinmixer and Qoutmixer to get a NEW Qinmixer, which is the NEW input for
! the next run.
        Xv(:, iscf_iteration) = Qinmixer(:)
        Fv(:, iscf_iteration) = Qoutmixer(:) - Qinmixer(:)

! Calculate the new sigma
        sigma = dot_product(Fv(:,iscf_iteration),Fv(:,iscf_iteration))
#ifdef GRID
        sigma = sqrt(sigma)
#else
        sigma = sigma/imix
#endif

! Initially - only perform simple extrapolation
        if (mix_order .eq. 1 .or. sigma .lt. scf_tolerance_set) then
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
            lwork = (mix_order - 1)**2
            allocate (work(lwork))
            allocate (ipiv(mix_order - 1))
            info = 0
            call dsysv('U', mix_order - 1, 1, amatrix, mix_order - 1,         &
     &                 ipiv, F_dot_delF, mix_order - 1, work, lwork, info)

! If there is an error, then just use the simple mixing
            if (info .ne. 0) then
              write (logfile,*) ' Error in Anderson, info = ', info
              if (mix_order .le. 2) stop
              mix_order = mix_order - 1
            end if
            deallocate (work, ipiv, amatrix)
          end do
          istep = iscf_iteration - mix_order

! Generate new guess at charges Eq. 7.7  (F_dot_delF is now gamma)
          Qinmixer(:) = Qinmixer(:) + beta_set*Fv(:,iscf_iteration)  ! First-order term
          do iloop = istep + 1, iscf_iteration - 1
            Qinmixer(:) =                                                        &
     &        Qinmixer(:) - F_dot_delF(iloop)*(delX(:,iloop) + beta_set*delF(:,iloop))
          end do
          deallocate (F_dot_delF)
        end if

! ===========================================================================
!                                 end mixing
! ===========================================================================

#ifdef GRID

! Check the total charge - before setting the density
        ztotal_out = 0.0d0
        do iatom = 1, t%natoms
          ! cut some lengthy notation
          pdenmat=>t%denmat(iatom)
          poverlap=>t%overlap(iatom)
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = t%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            pS_neighbors=>poverlap%neighbors(ineigh)
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              do inu = 1, norb_nu
                ztotal_out = ztotal_out                                      &
     &                      + pRho_neighbors%block(imu,inu)*pS_neighbors%block(imu,inu)
              end do
            end do ! end loop over matrix elements
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Set the new density matrix
        imix = 0
        do iatom = 1, t%natoms
          ! cut some lengthy notation
          pdenmat=>t%denmat(iatom)
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = t%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over neighbors
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              do inu = 1, norb_nu
                imix = imix + 1
                pRho_neighbors%block(imu,inu) = Qinmixer(imix)
              end do
            end do ! end loop over matrix elements
          end do ! end loop over neighbors
        end do ! end loop over atoms

        t%ztot = 0.0d0
        do iatom = 1, s%natoms
          in1 = t%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            t%ztot = t%ztot + species(in1)%shell(issh)%Qneutral
          end do
        end do

! Check the Ztot and we will rescale the density matrix elements so that we
! can ensure the charge is conserved.
        zcheck = 0.0d0
        do iatom = 1, t%natoms
          ! cut some lengthy notation
          pdenmat=>t%denmat(iatom)
          poverlap=>t%overlap(iatom)
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = t%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            pS_neighbors=>poverlap%neighbors(ineigh)
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              do inu = 1, norb_nu
                zcheck = zcheck                                              &
     &                  + pRho_neighbors%block(imu,inu)*pS_neighbors%block(imu,inu)
              end do ! inu
            end do ! imu
          end do ! neighbors.
        end do ! atoms

! Renormalize and set denmat_old - write the density matrix to a file
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.rho'
        open (unit = inpfile, file = slogfile, status = 'unknown', form = 'unformatted')
        do iatom = 1, t%natoms
          ! cut some lengthy notation
          pdenmat=>t%denmat(iatom)
          pdenmat_old=>t%denmat_old(iatom)
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = t%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            pRho_neighbors_old=>pdenmat_old%neighbors(ineigh)
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              do inu = 1, norb_nu
                pRho_neighbors%block(imu,inu) = pRho_neighbors%block(imu,inu)*t%ztot/zcheck
                pRho_neighbors_old%block(imu,inu) = pRho_neighbors%block(imu,inu)
              end do
            end do ! end loop over matrix elements
            write (inpfile) pRho_neighbors_old%block
          end do ! end loop over neighbors
        end do ! end loop over atoms
        close (inpfile)

! Check normalization
        zcheck = 0.0d0
        do iatom = 1, t%natoms
          ! cut some lengthy notation
          pdenmat=>t%denmat(iatom)
          poverlap=>t%overlap(iatom)
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = t%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            pS_neighbors=>poverlap%neighbors(ineigh)
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              do inu = 1, norb_nu
                zcheck = zcheck                                              &
     &                  + pRho_neighbors%block(imu,inu)*pS_neighbors%block(imu,inu)
              end do
            end do ! End loop over matrix elements
          end do ! End loop over neighbors
        end do ! End loop over atoms

#else

! Check the total charge
        ztotal_out = 0
        imix = 0
        do iatom = 1, t%natoms
          in1 = t%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            imix = imix + 1
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

#endif

! Write out results
        write (logfile,*) ' (Before renormalization) zouttot = ', ztotal_out
        write (logfile,*) ' (After  renormalization)  zcheck = ', zcheck
        write (logfile,*) ' (What it must be)           ztot = ', t%ztot
        write (logfile,*)
        write (logfile,202) sigma, scf_tolerance_set, iscf_iteration
#else
! Set sigma to something crazy in order to kick out of the scf loop.
        sigma = -999.0d0
        if(.FALSE.) t%volume=0
        if(.FALSE.) sigma=iscf_iteration
#endif

! Deallocate Arrays
! ===========================================================================
#ifdef DOGS
        deallocate (Qinmixer, Qoutmixer)

! Format Statements
! ===========================================================================
100     format (2x, ' Deviation (rms) of input/output charges = ', f9.4)
101     format (2x, ' Deviation (max) of input/output charges = ', f9.4)
201     format (2x, ' Renormalization of Qin:   renorm = ', 4x, f14.8)
202     format (' =====> sigma = ', e14.7,                                   &
     &          ' Must be less than', e14.7, ' SCF step = ', i3)
#endif

! End Subroutine
! ===========================================================================
        return
        end subroutine Qmixer


