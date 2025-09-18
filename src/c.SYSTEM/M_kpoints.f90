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

! M_kpoints
! Module Description
! ===========================================================================
!>       This is a module which will either read in the kpoints from a file
!! provided by the user - necessary for doing periodic boundary conditions.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Module Declaration
! ===========================================================================
        module M_kpoints

! Type Declaration
! ===========================================================================
! If we are interested in transitions, such as for the evolving the
! coefficients in order to calculate non-adiabatic couplings or for calculating
! absorption. We consider only a range of states near the valence and conduction
! band edges. The number of states that we consider is defined by ntransitions.
        type T_transition
          integer imap                 ! to which band does this transition map

          ! these are the time-dependent Schrodinger coefficients - these will
          ! be evolved in time according to Verlet or Runge-Kutta as desired
          ! these correspond to the non-adiabatic coupled states
          complex cna_old
          complex cna

          ! these are the coefficients of the adiabatic eigenstates - these will
          ! be evolved in time similarly to cna and cna_old
          complex, allocatable :: c_wf (:)

          ! this is the non-adiabatic coupling belonging to the transition state
          complex, allocatable :: djk (:)
          complex, allocatable :: djk_old (:)
          complex, allocatable :: ddjk (:)
        end type

        type T_kpoint
          real weight                          ! weight of kpoint

          integer, pointer :: ioccupy (:)     ! integer occupation number

          real, pointer :: eigen (:)          ! eigenvalues for k
          real, pointer :: eigen_old (:)      ! previous eigenvlues for k
          real, pointer :: deigen (:)         ! interpolted eigen values
          real, pointer :: foccupy (:)        ! occupation real value for k

          real, dimension (3) :: k            ! kpoint vector

          ! k-dependent S12 matrix
          complex, pointer :: S12matrix (:, :)

          ! the real and imaginary components of the eigenvectors
          complex, pointer :: c (:, :)

          ! the real and imaginary components of the Lowdin transformed eigenvectors
          complex, pointer :: c_Lowdin (:, :)

! If we are interested in transitions, such as for the evolving the
! coefficients in order to calculate non-adiabatic couplings or for calculating
! absorption. We consider only a range of states near the valence and conduction
! band edges. The number of states that we consider is defined by ntransitions.
          integer ntransitions                ! number of transitions

          type (T_transition), pointer :: transition (:)
          type (T_transition), pointer :: atransition (:, :)
        end type

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! read_fermie.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>    Allows for occupations beyond just the ground state - read in from file.
!
! ===========================================================================
        subroutine read_fermie (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iband
        integer ikpoint
        integer imu
        integer inpfile                     !< reading from which unit
        integer logfile                     !< writing to which unit
        integer noccupy

        real fband

        character (len = 25) :: slogfile

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! Open the file and read information.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.OCCUPATION'
        open (unit = inpfile, file = slogfile, status = 'old')

        write (logfile,*)
        write (logfile,*) ' Reading from the OCCUPATION file! '
        do ikpoint = 1, s%nkpoints
          s%kpoints(ikpoint)%ioccupy = 0
          s%kpoints(ikpoint)%foccupy = 0.0d0
          read (inpfile,*) noccupy
          if (noccupy .gt. s%norbitals) then
            write (*,*) ' noccupy > norbitals: from OCCUPATION file. '
            stop
          end if
          do imu = 1, noccupy
            read (inpfile,*) iband, fband
            s%kpoints(ikpoint)%ioccupy(imu) = iband
            s%kpoints(ikpoint)%foccupy(imu) = fband
          end do
        end do
        close (unit = inpfile)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_fermie


! ===========================================================================
! fermie.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the fermie energy and fermi occupations.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine fermie (s, qstate, efermi)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        real, intent (in) :: qstate

! Output
        real, intent (out) :: efermi

! Local Parameters
! ===========================================================================
        integer, parameter :: imax = 1000 !< maximum sc iterations
        integer, parameter :: nmax = 5000 !< cutoff for degeneracy check

        real, parameter :: tol = 1.0d-12

! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer imu
        integer inu
        integer iteration
        integer jkpoint

        real delta
        real emin
        real emax
        real qcharge
        real qztot
        real temperature

        real, pointer :: peigen_mu_k
        real, pointer :: peigen_nu_k

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Add in the qstate to the total charge
        qztot = s%ztot + qstate

! The subroutine fermie needs a temperature to calculate the occupations of
! the states so set temperature to some low value (1 eV = 11604 K).
        temperature = efermi_T/P_kconvert

! Find emin and emax. Also make sure degenerate eigenvalues are truly
! degenerate. However, if norbitals*nkpts*norbitals*nkpts is larger than nmax,
! then skip the degeneracy checking. Otherwise, the checking can take a while.
        if (s%norbitals**2*s%nkpoints**2 .lt. nmax) then
          emin = s%kpoints(1)%eigen(1)
          emax = s%kpoints(1)%eigen(s%norbitals)
          do ikpoint = 1, s%nkpoints
            do imu = 1, s%norbitals
              peigen_mu_k=>s%kpoints(ikpoint)%eigen(imu)
              if (peigen_mu_k .lt. emin) emin = peigen_mu_k
              if (peigen_mu_k .gt. emax) emax = peigen_mu_k
              do jkpoint = ikpoint, s%nkpoints
                do inu = imu, s%norbitals
                  peigen_nu_k=>s%kpoints(ikpoint)%eigen(inu)
                  if (abs(peigen_mu_k - peigen_nu_k) .lt. tol) then
                    s%kpoints(ikpoint)%eigen(inu) =                            &
     &                (peigen_mu_k + peigen_nu_k)/2.0d0
                    s%kpoints(ikpoint)%eigen(imu) = peigen_nu_k
                  end if
                end do
              end do
            end do
          end do
        else
          open (11, file = 'WARNINGS', status = 'unknown', position = 'append')
          write (11,*) '  '
          write (11,*) ' ************ WARNING ******* WARNING *********** '
          write (11,*) '          skipping the degeneracy checking  '
          write (11,*) '               in subroutine fermie'
          write (11,*) ' ************************************************ '
          close (11)
          emin = s%kpoints(1)%eigen(1)
          emax = s%kpoints(1)%eigen(s%norbitals)
          do ikpoint = 1, s%nkpoints
            do imu = 1, s%norbitals
              peigen_mu_k=>s%kpoints(ikpoint)%eigen(imu)
              if (peigen_mu_k .lt. emin) emin = peigen_mu_k
              if (peigen_mu_k .gt. emax) emax = peigen_mu_k
            end do
          end do
        end if

! The value of efermi must be between emin and emax
        iteration = 0
        qcharge = 0.0d0
        do while (abs(qcharge - qztot) .gt. tol .and. iteration .le. imax)
          iteration = iteration + 1

! Make a guess at efermi
          efermi = (emax + emin)/2.0d0
          qcharge = 0.0d0
          do ikpoint = 1, s%nkpoints
            do imu = 1, s%norbitals
              delta = (s%kpoints(ikpoint)%eigen(imu) - efermi)/temperature

! Skip exponential for big -/+ delta
              if (delta .gt. 10.0d0) then
                s%kpoints(ikpoint)%foccupy(imu) = 0.0d0
                s%kpoints(ikpoint)%ioccupy(imu) = 0
              else if (delta .lt. -10.0d0) then
                s%kpoints(ikpoint)%foccupy(imu) = 1.0d0
                s%kpoints(ikpoint)%ioccupy(imu) = 1
              else
                s%kpoints(ikpoint)%foccupy(imu) = 1.0d0/(1.0d0 + exp(delta))
                if (s%kpoints(ikpoint)%foccupy(imu) .gt. 1.0d-5) then
                  s%kpoints(ikpoint)%ioccupy(imu) = 1
                else
                  s%kpoints(ikpoint)%ioccupy(imu) = 0
                end if
              end if
              qcharge = qcharge                                              &
     &          + P_spin*s%kpoints(ikpoint)%foccupy(imu)                     &
     &                  *s%kpoints(ikpoint)%weight
            end do
          end do
! Narrow the range that efermi can fall into
          if (qcharge .gt. qztot) then
            emax = efermi
          else
            emin = efermi
          end if
        end do

! Print warning for going over maximum iterations.
        if (iteration .gt. imax) then
          open (11, file = 'WARNINGS', status = 'unknown', position = 'append')
          write (11,*)
          write (11,*) ' ************ WARNING ******* WARNING *********** '
          write (11,*) '        not under tolerance (toll) after ',imax
          write (11,*) '          iterations in subroutine fermie'
          write (11,*)
          write (11,*) '          qcharge = ', qcharge
          write (11,*) '          qztot = ', qztot
          write (11,*) '          emax = ', emax
          write (11,*) '          emin = ', emin
          write (11,*) ' ************************************************ '
          write (11,*)
          close (11)
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine fermie

! End Module
! ===========================================================================
        end module M_kpoints
