! copyright info:
!
!                             @Copyright 2025
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

! M_nonadiabatic_mdet.f90
! Program Description
! ===========================================================================
!>       This routine calculates the density matrices pieces that are related
!> building the nonadiabatic coupling vectors for caclulating the pieces of
!> the time dependent Schrodinger equation for nonadiabatic molecular dynamics.
!
! ===========================================================================
! Code written by:
! James P. Lewis (with Zhaofa Li at Synfuels China Technology)
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ============================================================================
! Module declaration
! ============================================================================
        module M_nonadiabatic_mdet

! /SYSTEM
        use M_configuraciones

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! initialize_mdet.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine gives the initial state for molecular dynamics with
! electronic transitions (MDET) (nonadiabatic calculation)
! ===========================================================================
! Code written by:
!> @author James P. Lewis (with Zhaofa Li at Synfuels China Technology)
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
        subroutine initialize_mdet (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters
! ===========================================================================
        real, parameter :: tol = 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                      !< counter over atoms
        integer iband, ikpoint             !< counter of band and kpoint
        integer iband_in                   !< counter over transitions
        integer in1                        !< species number
        integer inpfile                    !< reading from which unit
        integer issh
        integer logfile                    !< writing to which unit
        integer nfermi                     !< Fermi level state
        integer ipop

        real qztot                         !< total number of electrons
        real qcharge                         !< total number of electrons
        real foccupy

        character (len = 25) :: slogfile

        type(T_transition), pointer :: piband
        type(T_kpoint), pointer :: pkpoint

! Procedure
! ===========================================================================
! Initialize logfile
        write(s%logfile, *) "initialize_mdet.f"

        logfile = s%logfile
        inpfile = s%inpfile

! Loop over the atoms.
! Total charge - ztot
        s%ztot = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            s%ztot = s%ztot + species(in1)%shell(issh)%Qneutral
          end do
        end do

! Initialize the fermi occupations
        qztot = s%ztot
        nfermi = int(qztot)/2
        do ikpoint = 1, s%nkpoints      
          do iband = 1, nfermi
            s%kpoints(ikpoint)%foccupy(iband) = 1.0d0
            s%kpoints(ikpoint)%ioccupy(iband) = 1
          end do
          do iband = nfermi + 1, s%norbitals
            s%kpoints(ikpoint)%foccupy(iband) = 0.0d0
            s%kpoints(ikpoint)%ioccupy(iband) = 0         
          end do           
        end do

! Read information from .mdet.input file
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.mdet.inp'
        open (unit = inpfile, file = slogfile, status = 'old')
        write (logfile,*)
        write (logfile,*) ' Reading from mdet.input file! '
        do ikpoint = 1, s%nkpoints       

          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
        
! Allocate transition type and initialize imap
          read (inpfile,*) pkpoint%nbands
          allocate (pkpoint%transition(pkpoint%nbands))
          do iband = 1, pkpoint%nbands
          
            ! cut some more lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)
            allocate (pkpoint%eigen_old (s%norbitals)); pkpoint%eigen_old = 0.0d0
            read (inpfile,*) iband_in, foccupy, ipop
! FIXME! We might need to read in foccupy and set ioccupy to 1 when foccupy
! is not 0 (since including all states no fixed relationship between ioccupy
! and foccupy)

            ! initialize imap
            piband%imap = iband_in

            ! initialize mdet occupations
            pkpoint%foccupy(piband%imap) = foccupy

            ! NAC initialize the dij, dij_old and c_mdet
            allocate (piband%dij_old(3, s%natoms, pkpoint%nbands))
            piband%dij_old = 0.0d0
            allocate (piband%dij(3, s%natoms, pkpoint%nbands))
            piband%dij = 0.0d0
            allocate (piband%c_mdet(s%norbitals)); piband%c_mdet = cmplx(0.0d0, 0.0d0)

            ! MDET initialize the c_na
            allocate (piband%c_na(pkpoint%nbands))
            piband%c_na = cmplx(0.0d0, 0.0d0)
            piband%c_na(iband) = cmplx(1.0d0, 0.0d0)  
  
            if (foccupy .ge. 0.5d0) then
              pkpoint%ioccupy(piband%imap) = 1
            else
              write (logfile,*) ' The MDET input has some problems: foccupy=', foccupy
              stop
            end if
            write (logfile,*) ' Check Occupations After read mdet input: : '
            write (logfile,*) ' Band ', piband%imap, ' Occupation = ',       &
                                        pkpoint%foccupy(piband%imap)
          end do   ! end loop over bands

! Check the occupation
          qcharge = 0.0d0          
          do iband = 1, s%norbitals
            if (pkpoint%ioccupy(iband) .ne. 0) then
              qcharge = qcharge                                              &
     &          + P_spin*pkpoint%foccupy(iband)*pkpoint%weight
            end if
          end do
          if (abs(qcharge - qztot) .gt. tol) then
            write (logfile,*) ' Check Occupations After read mdet input: : '         
            write (logfile,*) '          qcharge = ', qcharge
            write (logfile,*) '          qztot = ', qztot
            write (logfile,*) ' Must stop in subroutine initialize_mdet'
            stop
          end if

        end do   ! end loop over kpoints
        close (unit = inpfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_mdet


! ===========================================================================
! initialize_nac.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>      This routine gives the initial state for 
!>      nonadiabatic coupliong vectors dij
!
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis (with Zhaofa Li at Synfuels China Technology)
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
        subroutine initialize_nac (s)
        implicit none


! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iband, ikpoint             !< counter of band and kpoint

        type(T_transition), pointer :: piband
        type(T_kpoint), pointer :: pkpoint

! Procedure
! ===========================================================================
! Read information from .mdet.input file
        do ikpoint = 1, s%nkpoints       

          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
        
          do iband = 1, pkpoint%nbands
          
            ! cut some more lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            piband%dij = 0.0d0
            piband%c_mdet = cmplx(0.0d0, 0.0d0)
          end do   ! end loop over bands

        end do   ! end loop over kpoints

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_nac

! ===========================================================================
! density_matrix_nac
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine store the density matix rho
!>       for nonadiabatic 
!
! ===========================================================================
        subroutine density_matrix_nac (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iband, ikpoint           !< counter of band and kpoint
        integer logfile                  !< writing to which unit

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! ****************************************************************************
!
!       C O M P U T E    N O N A D I A B A T I C   C O E F F I C I E N T S
! ****************************************************************************
        write (logfile,*)
        write (logfile,*) ' Calculating density matrix elements for '
        write (logfile,*) ' nonadiabatic coupling vectors based on transitions '
        write (logfile,*) ' from iband to jband. '

! Loop over the special k points.
        do ikpoint = 1, s%nkpoints

          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)

! Loop over all bands
          do iband = 1, pkpoint%nbands
            ! cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            ! set the coefficient for the iband transition
            piband%c_mdet = pkpoint%c(:, piband%imap)

          end do   ! end loop over bands
        end do   ! end loop over kpoints
        
! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i4, 2(2x, f10.6))
200     format (4(2x, f12.4))

! End Subroutine
! ===========================================================================
        return
        end subroutine density_matrix_nac

! ===========================================================================
! writeout_density_nac.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine is a utility to write out the density matrix.
!
! ===========================================================================
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
        subroutine writeout_density_nac (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint                !< counter of band and kpoint
        integer iband                  !< counter of band
        integer inu                    !< counters for mu, nu

        character (len = 25) :: slogfile

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Formatted file needed for Multimwfn
! Only for gamma
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.cdcoeffs-mdet-mwfn'
        open (unit = 22, file = slogfile, status = 'replace')
        do ikpoint = 1, s%nkpoints
          write (22,"('Kpoint=',i10)") ikpoint

          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands
            ! cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            write (22,*)
            write (22,"('Index=',i10)") piband%imap
            write (22,"('Type=',i2)") 0
            write (22,"('Energy=',1PE16.8)") pkpoint%eigen(piband%imap)
            write (22,"('Occ=',f12.8)") pkpoint%foccupy(piband%imap) * P_spin
            write (22,"('Sym= ?')")
            write (22,"('$Coeff')")

! write out the coefficient
            write (22,"(5(1PE16.8))") (real(piband%c_mdet(inu)), inu = 1, s%norbitals)
          end do   ! end loop over bands
          write (22,*)
        end do   ! end loop over kpoints
        close (unit = 22)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (75('*'))
102     format (1x, ' Matrices connected to neighbor ', i3, ',',               &
     &          ' jatom = ', i4, ', ', ' mbeta = ', i4, ', ', ' d = ', f6.3)
103     format (75('='))
104     format (9f8.3)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_density_nac


! ===========================================================================
! build_dij_nac
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine build nonadiabatic coupliongs dij.
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
! Subroutine Declaration
! ===========================================================================
        subroutine build_dij_nac (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
!NAC MODIFIED
        real, parameter :: tolnac = 0.0001d0

! Variable Declaration and Description
! ===========================================================================
        integer ikpoint                    !< counter of band and kpoint
        integer iband, jband               !< counter of transitions

        real eigen_i, eigen_j              !< eigen values in band i and j
        real diff, sgn

        real, dimension (:, :), allocatable :: temp

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
        ! Initialize temp used in antisymmetry of dij
        allocate (temp (3, s%natoms))

        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands - 1

            ! Cut some lengthy notation
            nullify (piband); piband=>pkpoint%transition(iband)

            ! set the iband eigenvalue
            eigen_i = pkpoint%eigen(piband%imap)

            do jband = iband + 1, pkpoint%nbands

              ! Cut some lengthy notation
              nullify (pjband); pjband=>pkpoint%transition(jband)

              ! set the jband eigenvalue
              eigen_j = pkpoint%eigen(pjband%imap)

              diff = abs(eigen_i - eigen_j )
              if (diff .lt. tolnac) then
                write (s%logfile,*) ' TWO EIGENVALUES VERY CLOSE'
                write (s%logfile,*) ' iband', piband%imap, eigen_i
                write (s%logfile,*) ' jband', pjband%imap, eigen_j
                write (s%logfile,*) ' The nonadiabatic coupling is'
                write (s%logfile,*) ' NOT CALCULATED'
                sgn = sign(1.0_dp, eigen_i - eigen_j)
                piband%dij(:, :, jband) = piband%dij(:,:,jband)/(tolnac*sgn)

                ! NAC force anti-symmetry for NAC
                temp = - piband%dij(:, :, jband)
                pjband%dij(:, :, iband) = temp
              else 
                piband%dij(:, :, jband) = piband%dij(:,:,jband)/(eigen_i - eigen_j)

                ! NAC force anti-symmetry for NAC
                temp = - piband%dij(:, :, jband)
                pjband%dij(:, :, iband) = temp
              end if ! end if check for degeneracy
            end do ! end loop over jband
          end do ! end loop over iband
        end do  ! end loop over kpoints
        deallocate (temp)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine build_dij_nac


! ===========================================================================
! writeout_dij_nac
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine writes out nonadiabatic couplings dij.
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
! Subroutine Declaration
! ===========================================================================
        subroutine writeout_dij_nac (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom              !< counter over atoms and neighbors
        integer ikpoint                    !< counter of band and kpoint
        integer iband, jband                     !< counter of transitions
        integer logfile                     !< writing to which unit

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
        logfile = s%logfile
        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands - 1

            ! Cut some lengthy notation
            nullify (piband); piband=>pkpoint%transition(iband)

            do jband = iband + 1, pkpoint%nbands

              ! Cut some lengthy notation
              nullify (pjband); pjband=>pkpoint%transition(jband)

              write (logfile,*)
              write (logfile,103) ' The nonadibaitc couplings: '
              write (logfile,*) piband%imap, pjband%imap
              write (logfile,100)
              write (logfile,101)
              write (logfile,100)
              do iatom = 1, s%natoms
                write (logfile,102) 'dij', iatom,                            &
     &                                     s%atom(iatom)%species%symbol,     &
        &                                  piband%dij(:, iatom, jband)
              end do 
              write (logfile,100)
              write (logfile,*)
              write (logfile,*)
              write (logfile,103) ' The nonadiabatic couplings: '
              write (logfile,*) pjband%imap, piband%imap
              write (logfile,100)
              write (logfile,101)
              write (logfile,100)
              do iatom = 1, s%natoms
                write (logfile,102) 'dij', iatom,                            &
        &                                  s%atom(iatom)%species%symbol,     &
        &                                  pjband%dij(:, iatom, iband)
              end do 
              write (logfile,100)
              write (logfile,*)
            end do ! end loop over jband
          end do ! end loop over iband
        end do  ! end loop over kpoints

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (4x, 70('='))
101     format (4x, 'dij ', 'Atom # ', 2x, ' Type ', 5x,   &
     &              ' x ', 9x, ' y ', 9x, ' z ')
102     format (4x, A,  i5, 7x, a2, 3(2x,ES10.3))
103     format (4x, A)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_dij_nac

! ===========================================================================
! nonadiabatic_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine implements state switching in nonadiabatic dynamics.
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
! Subroutine Declaration
! ===========================================================================
        subroutine nonadiabatic_mdet (s, itime_step)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.
        integer itime_step                       !< counter of nuclear step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ikpoint                    !< counter of kpoint
        integer iband                      !< counter of bands

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Procedure
! ===========================================================================
        write (s%logfile, *) "Before integration"
        call writeout_c_na_mdet (s)

! Update the dij_old, eigen_old, and vatom_old
        if (itime_step .eq. 1) call save_stuff_mdet (s)
          call evolve_ks_states_mdet(s)
          call save_stuff_mdet (s)
!         call dij_phase (s)

          write (s%logfile, *) ' After integration '
          call writeout_c_na_mdet (s)
          call fewest_switches_mdet(s)

          write (s%logfile, *) ' Time Step ', itime_step
          do ikpoint = 1, s%nkpoints

            ! cut some lengthy notation
            nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
            do iband = 1, pkpoint%nbands

              ! cut some lengthy notation
              nullify (piband); piband => pkpoint%transition(iband)
              write (s%logfile,*) ' Band ', piband%imap, ' Occupation = ',   &
                                            pkpoint%foccupy(piband%imap)
            end do
          end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine nonadiabatic_mdet

! ===========================================================================
! writeout_c_na_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine writes out nonadiabatic wavefunction coefficient c_na
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
! Subroutine Declaration
! ===========================================================================
        subroutine writeout_c_na_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ikpoint                    !< counter of band and kpoint
        integer iband, jband                     !< counter of transitions

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Procedure
! ===========================================================================
        do ikpoint = 1, s%nkpoints

          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
          write (s%logfile,*) ' Coeffients c_na of nonadiabatic wavefunctions '

          do iband = 1, pkpoint%nbands

            ! cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            write (s%logfile,"(6(1PE16.8))") (piband%c_na(jband), jband = 1, pkpoint%nbands)
          end do
          write (s%logfile,*)
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (4x, 70('='))
101     format (4x, 'dij ', 'Atom # ', 2x, ' Type ', 5x,   &
     &              ' x ', 9x, ' y ', 9x, ' z ')
102     format (4x, A,  i5, 7x, a2, 3(2x,ES10.3))
103     format (4x, A)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_c_na_mdet

! ===========================================================================
! evolve_ks_states_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This routine integrates the TD equations for the coefficients c_na of the
!> TD-wfs: phi(i) = \Sum_ij c_na(i,j)*psi(j)
!> use Runge-Kutta 4th order
!>
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
! Subroutine Declaration
! ===========================================================================
        subroutine evolve_ks_states_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: nddt = 100          !< number of electron steps

! Variable Declaration and Description
! ===========================================================================
        integer ikpoint            !< counter of band and kpoint
        integer iband              !< counter of transitions
        integer iteration          !< counter of electron step

        real ddt                   !< length of electron step    
        real step                  !< length of rk4 step    

        complex, dimension(:), allocatable :: dc_na, dc_tmp, c_tmp

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          allocate (dc_na (pkpoint%nbands))
          allocate (dc_tmp (pkpoint%nbands))
          allocate (c_tmp (pkpoint%nbands))

          ddt = dt/nddt
          do iteration = 1, nddt
            do iband = 1, pkpoint%nbands

              ! Cut some lengthy notation
              nullify (piband); piband => pkpoint%transition(iband)

! step 1      
! Interpolation       
              c_tmp = piband%c_na
              step = real((iteration - 1.0d0)/nddt)
              call couplings  (s, ikpoint, pkpoint%nbands, step, c_tmp, dc_tmp)
              dc_na = dc_tmp/6.0d0
! step 2      
! Interpolation 
              c_tmp = piband%c_na + dc_tmp*ddt*0.5d0
              step = real((iteration - 0.5d0)/nddt)
              call couplings (s, ikpoint, pkpoint%nbands, step, c_tmp, dc_tmp)
              dc_na = dc_na + dc_tmp/3.0d0     
! step 3      
! Interpolation 
              c_tmp = piband%c_na + dc_tmp*ddt*0.5d0
              step = real((iteration - 0.5d0)/nddt) 
              call couplings (s, ikpoint, pkpoint%nbands, step, c_tmp, dc_tmp)
              dc_na = dc_na + dc_tmp/3.0d0       
! step 4      
! Interpolation 
              c_tmp = piband%c_na + dc_tmp*ddt
              step = real(iteration/nddt)
              call couplings (s, ikpoint, pkpoint%nbands, step, c_tmp, dc_tmp)
              dc_na = dc_na + dc_tmp/6.0d0 

! Integrate coefficients c_na
              piband%c_na = piband%c_na + dc_na*ddt    

            end do ! end loop over iband
          end do ! end time loop
          deallocate (c_tmp, dc_tmp, dc_na)
        end do ! end loop over kpoints

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine evolve_ks_states_mdet

! ===========================================================================
! save_stuff_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine saves some variables needed for the time
!> integration of the TD-wfs
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
! Subroutine Declaration
! ===========================================================================
        subroutine save_stuff_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                           !< counter over atoms
        integer ikpoint                         !< counter of kpoint
        integer iband                           !< counter of transitions

! NAC stuff
        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Procedure
! ===========================================================================
! Save all the velocities
        do iatom = 1, s%natoms
          s%atom(iatom)%vatom_old = s%atom(iatom)%vatom
        end do

        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          pkpoint%eigen_old = pkpoint%eigen

          do iband = 1, pkpoint%nbands

            ! Cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            piband%dij_old = piband%dij
          end do ! end loop over iband
        end do ! end loop over kpoints

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine save_stuff_mdet

! ===========================================================================
! dij_phase
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine modify phase of dij.
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
! Subroutine Declaration
! ===========================================================================
        subroutine dij_phase(s)

        implicit none

        include '../include/constants.h'
! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer :: iatom
        integer :: ikpoint
        integer :: iband, jband

        real :: dij_old_norm, dij_norm
        real :: dij_doc_product
        real :: cos_angle

        real, dimension (:, :), allocatable :: temp

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Allocate Arrays
! ===========================================================================
        allocate (temp (3, s%natoms))

! Procedure
! ===========================================================================
        do ikpoint = 1, s%nkpoints ! begin of kpoints loop

          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands - 1 ! begin of iband loop

            ! cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            do jband = iband + 1, pkpoint%nbands ! begin of jband loop
              dij_old_norm = 0.d0; dij_norm = 0.d0; dij_doc_product = 0.d0
              do iatom = 1, s%natoms ! begin of atoms loop
! Calculate the norm of dij_old
                dij_old_norm = dij_old_norm                                  &
              &   + dot_product(piband%dij_old(:,iatom,jband), piband%dij_old(:,iatom,jband))
! Calculate the norm of dij
                dij_norm = dij_norm                                          &
              &   + dot_product(piband%dij(:,iatom,jband), piband%dij(:,iatom,jband))
! Calculate the dot product of two dij

                dij_doc_product = dij_doc_product                            &
              &   + dot_product(piband%dij_old(:,iatom,jband), piband%dij(:,iatom,jband))
              end do ! end of atoms loop
              dij_old_norm = sqrt(dij_old_norm); dij_norm = sqrt(dij_norm)

              if ((dij_old_norm .eq. 0) .or. (dij_norm .eq. 0)) then
                cos_angle = 1.0
              else
                cos_angle = dij_doc_product/(dij_old_norm*dij_norm)
              end if

              if (cos_angle .ge. 0) then
                piband%dij(:,:,jband) = piband%dij(:,:,jband)
              else
                piband%dij(:,:,jband) = -piband%dij(:,:,jband)
              end if

              ! cut some lengthy notation
              nullify (pjband); pjband => pkpoint%transition(jband)

              temp = piband%dij(:,:,jband)
              pjband%dij(:,:,iband) = temp

            end do ! end loop of jband
          end do ! end loop of iband
        end do ! end loop of kpoints
        deallocate (temp)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        end subroutine dij_phase


! ===========================================================================
! couplings
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine gives the derivative wrt time for the coefficients
!>       c_wf of the TD-wfs
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
! Subroutine Declaration
! ===========================================================================
        subroutine couplings (s, ikpoint, nbands, step, c_tmp, dc_tmp)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        integer, intent(in) :: ikpoint, nbands   !< number of bands

        real, intent(in) :: step                 !< step of rk4

        !< coefficients of wavefunctions
        complex, intent(in), dimension(nbands) :: c_tmp

! Output
        complex, intent(out), dimension (nbands) :: dc_tmp   !< derivative of coefficients

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom              !< counter over atoms
        integer jband, kband       !< counter of transitions
        integer ix                 !< counter over spatial dimension

        real djkov                 !< dot product of dij and veolocity
        real deigen, eigen_tmp     !< step and temporary eigen value in rk4
        real dv, v_tmp             !< step and temporary velocity in rk4
        real ddij, dij_tmp         !< step and temporary NAC in rk4

! NAC stuff
        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
        nullify (pkpoint); pkpoint => s%kpoints(ikpoint)

        ! Initialize dc_tmp
        dc_tmp = cmplx(0.0d0, 0.0d0)

        do jband = 1, pkpoint%nbands

          ! Cut off length
          nullify (pjband); pjband => pkpoint%transition(jband)

          deigen = pkpoint%eigen(pjband%imap) -  pkpoint%eigen_old(pjband%imap)
          eigen_tmp = pkpoint%eigen_old(pjband%imap) + deigen*step

          do kband = 1, pkpoint%nbands

            djkov = 0.0d0
            do iatom = 1, s%natoms
              do ix = 1, 3
                dv = s%atom(iatom)%vatom(ix) - s%atom(iatom)%vatom_old(ix)
                v_tmp = s%atom(iatom)%vatom_old(ix) + dv*step

                ddij = pjband%dij(ix,iatom,kband) - pjband%dij_old(ix,iatom,kband)
                dij_tmp = pjband%dij_old(ix,iatom,kband) + ddij*step

                ! dij_tmp index (jband, kband)
                djkov = djkov + v_tmp*dij_tmp
              end do ! end loop over ix
            end do ! end loop over iatom

            ! c_tmp index (iband, kband) and  dc_tmp index (iband, jband)
            dc_tmp(jband) = dc_tmp(jband) - djkov*c_tmp(kband)
          end do ! end loop over kband

          ! eigen_tmp index (pjband%imap), c_tmp index (iband, jband)
          ! and dc_tmp index (iband, jband)
          dc_tmp(jband) =                                                    &
     &      dc_tmp(jband) - (cmplx(0.0d0, 1.0d0)/hbar)*eigen_tmp*c_tmp(jband)
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine couplings


! ===========================================================================
! fewest_switches_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This routine determines the hoppings between Kohn-Sham state.
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
! Subroutine Declaration
! ===========================================================================
        subroutine fewest_switches_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                           !< counter over atoms and neighbors
        integer ikpoint                         !< counter of band and kpoint
        integer iband, jband                    !< counter of transitions
        integer ix                              !< counter of spatial dimension
        integer iswitch                         !< switch band

        real djiov                              !< inner product of NAC and velocity

        !< probability of switches
        real, dimension(:), allocatable :: probability

        complex aii, aji, bji                   !< coefficients of FSSH algorithm

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
! Calculate by fewest switches algorithm gives hopping probabilities from
! the current state we follow the possible transitions associated with states
        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)
          allocate (probability (pkpoint%nbands))

          do iband = 1, pkpoint%nbands - 1
 
            ! Cut some lengthy notation
            nullify (piband); piband=>pkpoint%transition(iband)

            aii = real(conjg(piband%c_na(iband))*piband%c_na(iband))

            do jband = iband + 1, pkpoint%nbands
              ! Cut some lengthy notation
              nullify (pjband); pjband=>pkpoint%transition(jband)
   
              aji = piband%c_na(jband)*conjg(piband%c_na(iband))
              djiov = 0.0d0
              do iatom = 1, s%natoms
                do ix = 1, 3
                  djiov = djiov + s%atom(iatom)%vatom(ix)*pjband%dij(ix,iatom,iband)
                end do
              end do

              bji = -2.0d0*real(conjg(aji)*djiov)
! Notice eq(29) in  J. Chem. Phys. 15 September 1994; 101 (6): 4657–4667
! Probability of the iband ---> jband transition
              probability(jband) = bji*dt/aii
              write (s%logfile,*) ' probability = ', piband%imap, pjband%imap, probability(jband)
              if (probability(jband) .lt. 0.0d0) probability(jband) = 0.0d0
            end do ! do jband = 1, nbands

! This only works for gamma kpoint
            if (s%nkpoints .eq. 1) then      
              ! Adjudge possible transition between current band <---> switch band 
              call mc_switch (s, iband, pkpoint%nbands, probability, iswitch)

              if (iswitch .ne. 0) then      
                ! Change configuration    
                ! Check energy conservation
                ! Rescaling velocity
                ! Only one transition is permitted 
                call transition (s, iband, iswitch)
              else
                write(s%logfile, *) ' Impossible Transition '
              end if
            else
              write(s%logfile, *) ' Cannot execute FSSH for non-gamma '  
            end if ! end kpoint if

          end do ! end loop over ibands
          deallocate (probability)
        end do ! end loop over kpoints          

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format ( 1i4, 4f8.4)
102     format ( 3i4, 5f8.4)

! End Subroutine
! ===========================================================================
        return
        end subroutine fewest_switches_mdet

! ===========================================================================
! mc_switch
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This routine determines the hoppings switches between Kohn-Sham states
!> based on a Monte-Carlo approach (see J.C. Tully, JCP 93, 1061 (1990).
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
! Subroutine Declaration
! ===========================================================================
        subroutine mc_switch (s, iband, nbands, probability, iswitch)
        implicit none     

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        integer, intent(in) :: iband
        integer, intent(in) :: nbands

        real, dimension (nbands), intent(in) :: probability

! Output
        integer, intent(out) :: iswitch

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer :: jband
        integer :: ikpoint

        real :: aux
        real ::  xrand        

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
! Check that sum of probabilities is smaller than 1
        ikpoint = 1
        ! Cut some lengthy notion
        nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
        nullify (piband); piband => pkpoint%transition(iband)

        write(s%logfile, *) ' Welcome to mc_switch.f '

        ! Cut some lengthy notion
        nullify (pkpoint); pkpoint => s%kpoints(1)
        nullify (piband); piband => pkpoint%transition(iband)

        aux = 0.0d0
        ! NAC MODIFIED IT SO THAT FOR TRANSITIONS NOT FOR BANDS        
        do jband = iband + 1, nbands

          ! cut some lenghty notation
          nullify (pjband); pjband => pkpoint%transition(jband)

          ! Consider only allowed transitions
          if (pkpoint%ioccupy(piband%imap) .ne. 0 .or.                       &
     &        pkpoint%ioccupy(pjband%imap) .ne. 0) then
            aux = aux + probability(jband)
          end if

        end do ! end loop over jband

        if (aux .gt. 1.0d0) then
          write (s%logfile, *) ' Sum of probabilities greater than 1 in mc_switches.f90 '
          write (s%logfile, *) ' total probabilty', aux, ' for state', piband%imap
        end if

        iswitch = 0
        aux = 0.0d0
        ! NAC MODIFIED IT SO THAT FOR TRANSITIONS NOT FOR BANDS        
        do jband = iband + 1, nbands

          ! cut some lengthy notation
          nullify (pjband); pjband => pkpoint%transition(jband)

          ! Random numbers for Monte-Carlo
          call random_number(xrand)
          write (s%logfile, *) 'random', xrand

          ! Consider only allowed transitions
          if (pkpoint%foccupy(piband%imap) .gt. 0.25 .or.                    &
     &        pkpoint%foccupy(pjband%imap) .lt. 0.75) then
              aux = aux + probability(jband)
              if (aux .gt. xrand) then
                iswitch = jband
                write(s%logfile, *) ' Possible Transition ', piband%imap, '<--->', pjband%imap
                exit
              end if
          end if          
        end do
!----------------------------------------------------------
! If iswitch = 0, no switch (hopping) between states
! Otherwise, switch from current state (piband%imap in
! fewest_switches subroutine) to state "iswitch"
!----------------------------------------------------------

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine mc_switch


! ===========================================================================
! transition
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>     This routine performs the electronic transition
!>     current band --> switch band
!>     and rescales the velocities after the transition 
!>     to conserve total energy
!>     The velocities are re-scaled along the direction of the nonadiabatic
!>     coupling (see JCP 101, 4657 (1994) )
!
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
! Subroutine Declaration
! ===========================================================================
        subroutine transition (s, iband, iswitch)
        implicit none     

        include '../include/constants.h'
        
! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        integer, intent(in) :: iband             !< current band
        integer, intent(in) :: iswitch           !< switch band

! Parameters and Data Declaration
! ===========================================================================
        real, parameter :: tolaa = 1.0d-08

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ikpoint                   !< counter of atomn and kpoint
        integer ix                               !< counter of spatial dimension
        integer in1                              !< counter of species
        integer kband                            !< counter of bands                  
        integer band_current                     

        real aa, bb, cc, alfa                    !< coefficients of quadratic equation
        real ejump                               !< energy difference in eV
        real energy                    !< potential and kinetics energy
        real xrand                               !< random number
        real dij_norm                            !< norm of nonadiabatic coupling        
        real etot_before, etot_after             !< For other theories

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband
        type(T_transition), pointer :: pkband

! Procedure
! ===========================================================================
! Switch from current band ---> switch band
        ikpoint = 1
        ! cut some lengthy notation
        nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
        nullify (piband); piband => pkpoint%transition(iband)
        nullify (pjband); pjband => pkpoint%transition(iswitch)

        ! Calculate energy jump
        ! Use Energy gap between before and after hoppings
        call recalculation (s, etot_before)

        ! Update configureation after surface hopping
        if (pkpoint%foccupy(piband%imap) .lt. 0.75 .and.                     &
     &      pkpoint%foccupy(piband%imap) .gt. 0.25 .and.                     &
     &      pkpoint%foccupy(pjband%imap) .lt. 0.75 .and.                     &
     &      pkpoint%foccupy(pjband%imap) .gt. 0.25 ) then
          write(s%logfile, *) ' Single Occupied '
        ! Ensure that both single occupied bands are equivalently excited   
          call random_number(xrand)

          ! Upwards Transition 
          if (xrand .lt. 0.5d0) then 
            write(s%logfile, *) ' Upwards Transition '
            write (s%logfile,*) ' Possible Transition ', piband%imap, '--->', pjband%imap
            band_current = iswitch
            pkpoint%foccupy(piband%imap) = pkpoint%foccupy(piband%imap) - 0.50d0
            pkpoint%foccupy(pjband%imap) = pkpoint%foccupy(pjband%imap) + 0.50d0
        !     ejump = pkpoint%eigen(pjband%imap) - pkpoint%eigen(piband%imap)
          ! Downwards Transition  
          else 
            write(s%logfile, *) ' Downwards Transition '
            write (s%logfile,*) ' Possible Transition ', pjband%imap, '--->', piband%imap
            band_current = iband
            pkpoint%foccupy(piband%imap) = pkpoint%foccupy(piband%imap) + 0.50d0
            pkpoint%foccupy(pjband%imap) = pkpoint%foccupy(pjband%imap) - 0.50d0
        !     ejump = pkpoint%eigen(piband%imap) - pkpoint%eigen(pjband%imap)
          end if   

        ! Update configureation after surface hopping
        ! Upwards Transition 
        else if (pkpoint%foccupy(piband%imap) .gt. 0.25 .and. &     
        &        pkpoint%foccupy(pjband%imap) .lt. 0.75) then        
          write(s%logfile, *) "Upwards Transition"            
          write (s%logfile,*) ' Possible Transition ', piband%imap, '--->', pjband%imap
          band_current = iswitch          
          pkpoint%foccupy(piband%imap) = pkpoint%foccupy(piband%imap) - 0.50d0
          pkpoint%foccupy(pjband%imap) = pkpoint%foccupy(pjband%imap) + 0.50d0
        !   ejump = pkpoint%eigen(pjband%imap) - pkpoint%eigen(piband%imap)

        ! Update configureation after surface hopping
        ! Downwards Transition  
        else if (pkpoint%foccupy(piband%imap) .lt. 0.75 .and. &     
        &        pkpoint%foccupy(pjband%imap) .gt. 0.25) then 
          write(s%logfile, *) "Downwards Transition"           
          write (s%logfile,*) ' Possible Transition ', pjband%imap, '--->', piband%imap
          band_current = iband          
          pkpoint%foccupy(piband%imap) = pkpoint%foccupy(piband%imap) + 0.50d0
          pkpoint%foccupy(pjband%imap) = pkpoint%foccupy(pjband%imap) - 0.50d0
        !   ejump = pkpoint%eigen(piband%imap) - pkpoint%eigen(pjband%imap)

        ! Update configureation after surface hopping
        end if  
        
        write (s%logfile,*) ' Occupation of Possible Transition '
        do kband = 1, pkpoint%nbands

          ! cut some lengthy notation
          nullify (pkband); pkband => pkpoint%transition(kband)

          write (s%logfile,*) ' Band ', pkband%imap, ' Occupation = ',       &
                                        pkpoint%foccupy(pkband%imap)
        end do
        
        ! Update Integer Occupation
        do kband = 1, pkpoint%nbands

          ! cut some lengthy notation
          nullify(pkband); pkband => pkpoint%transition(kband)
          if (pkpoint%foccupy(pkband%imap) > -0.25d0 .and.                   &
     &        pkpoint%foccupy(pkband%imap) <  0.25d0 ) then
            pkpoint%ioccupy(pkband%imap) = 0
          else if (pkpoint%foccupy(pkband%imap) > 0.25d0 .and.               &
     &             pkpoint%foccupy(pkband%imap) < 1.25d0 ) then
            pkpoint%ioccupy(pkband%imap) = 1
          else
            pkpoint%ioccupy(pkband%imap) = 0
          end if
        end do

        ! Update forces and energies after the configuration change
        call recalculation (s, etot_after)
        ejump = etot_after - etot_before

! Transform energy shift from eV to atomic units ( amu*(angs/fs)**2 )
        write (s%logfile,*) 'ejump (eV) = ', ejump
        energy = ejump*P_fovermp

        write (s%logfile,*) 'energy (dynamical units) = ', energy

! ===========================================================================
! Find out if transition current band --> switch band is accesible 
! (i.e. if there is enough kinetic energy) 
! Calculation of norm of dij
        dij_norm = 0.0d0
        do iatom = 1, s%natoms
          dij_norm = dij_norm &
  &     + dot_product (piband%dij(:,iatom,iswitch), piband%dij(:,iatom,iswitch))
        end do
        dij_norm = sqrt(dij_norm)

        aa = 0.0d0
        bb = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          do ix = 1, 3
            aa = aa + 0.50d0*(piband%dij(ix,iatom,iswitch) / dij_norm        &
        &                   * piband%dij(ix,iatom,iswitch) / dij_norm) / species(in1)%xmass
            bb = bb + s%atom(iatom)%vatom(ix)*piband%dij(ix,iatom,iswitch) / dij_norm
          end do
        end do
        write (s%logfile,*) ' aa, 4*energy*aa =', aa, 4.0d0*aa*energy
        write (s%logfile,*) ' bb, bb**2 = ', bb, bb**2
        cc = bb**2 - 4.0d0*aa*energy
        write (s%logfile,*) ' cc = ', cc

        if (aa .gt. tolaa) then
          if (cc .ge. 0.0d0) then 
            if (band_current .eq. iswitch) then
              write (s%logfile,*) ' Successful Transition ', piband%imap, '--->', pjband%imap   
            else if (band_current .eq. iband) then
              write (s%logfile,*) ' Successful Transition ', pjband%imap, '--->', piband%imap               
            end if  

            if (bb .ge. 0.0d0) then
              alfa = (bb - sqrt(cc))/(2.0d0*aa)
            else
              alfa = (bb + sqrt(cc))/(2.0d0*aa)
            end if
            write (s%logfile, *) ' alfa = ', alfa   

          else !(cc .gt. 0)
            ! Revert Occupation
            if (band_current .eq. iswitch) then
              write (s%logfile,*) ' Frustrated Transition ', piband%imap, '--->', pjband%imap
              pkpoint%foccupy(piband%imap) = pkpoint%foccupy(piband%imap) + 0.50d0
              pkpoint%foccupy(pjband%imap) = pkpoint%foccupy(pjband%imap) - 0.50d0          
            else if (band_current .eq. iband) then
              write (s%logfile,*) ' Frustrated Transition ', pjband%imap, '--->', piband%imap              
              pkpoint%foccupy(piband%imap) = pkpoint%foccupy(piband%imap) - 0.50d0
              pkpoint%foccupy(pjband%imap) = pkpoint%foccupy(pjband%imap) + 0.50d0     
            end if                
! a) Do nothing
            
! b) Reflection  pag. 4664 in J. Chem. Phys. 15 September 1994; 101 (6): 4657–4667
            alfa = bb/aa 
            write (s%logfile, *) ' alfa = ', alfa

          end if  !(cc .gt. 0)
        else  !(aa .gt. tolaa)
          alfa = 0.0d0
          write (s%logfile, *) ' alfa = ', alfa    
        end if  !(aa .gt. tolaa)

        ! Write out Occupation after Energy Conservation
        do kband = 1, pkpoint%nbands
          nullify (pkband); pkband => pkpoint%transition(kband)
          write (s%logfile,*) ' Band ', pkband%imap, ' Occupation = ',       &
          pkpoint%foccupy(pkband%imap)
        end do

        ! RESCALING VELOCITIES
        write (s%logfile,*) ' RESCALING VELOCITIES '
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass             
          do ix = 1, 3
            s%atom(iatom)%vatom(ix) = s%atom(iatom)%vatom(ix)                &
     &       - alfa*(piband%dij(ix,iatom,iswitch)/dij_norm)/species(in1)%xmass
          end do
        end do         

        ! Update Integer Occupation
        do kband = 1, pkpoint%nbands

          ! cut some lengthy notation
          nullify(pkband); pkband => pkpoint%transition(kband)

          if (pkpoint%foccupy(pkband%imap) > -0.25d0 .and.                   &
     &        pkpoint%foccupy(pkband%imap) <  0.25d0) then
            pkpoint%ioccupy(pkband%imap) = 0
          else if (pkpoint%foccupy(pkband%imap) > 0.25d0 .and.               &
     &             pkpoint%foccupy(pkband%imap) < 1.25d0 ) then
            pkpoint%ioccupy(pkband%imap) = 1
          else
            pkpoint%ioccupy(pkband%imap) = 0
          end if
        end do

        ! Update energy and force
        call recalculation (s, etot_after)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine transition

! recalculation
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine recalculates the energy and force in MDET.
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
! Subroutine Declaration
! ===========================================================================
        subroutine recalculation (s, etot)
! /SYSTEM
        use M_vdW

 ! /ASSEMBLERS
        use M_assemble_2c
        use M_assemble_3c
        use M_assemble_ewald
        use M_assemble_usr
        use M_assemble_vxc
        use M_assemble_PP_2c
        use M_assemble_PP_3c

! /DASSEMBLERS
        use M_Dassemble_2c
        use M_Dassemble_PP_2c
        use M_Dassemble_3c
        use M_Dassemble_PP_3c
        use M_Dassemble_vxc
        use M_Dassemble_usr
        use M_Dassemble_ewald
        use M_build_forces
! /SOLVESH
        use M_kspace
        use M_density_matrix        
! /SCF
        use M_charges
! /MPI
        use M_mpi
     
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        real, intent(out) ::  etot               !< total energy        

! Parameters and Data Declaration
! ===========================================================================
        include '../include/constants.h'

! Variable Declaration and Description
! ===========================================================================
        integer iatom                     !< counter over atoms and neighbors
        integer in1
        integer iscf_iteration

        real sigma                          !< difference for SCF
        real rms                            !< RMS of the forces

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_scf_begin, time_scf_end
        real time_forces_begin, time_forces_end
        real time_initial, time_final

! Energies
        real ebs                                 ! band-structure energy
        real efermi                              ! Fermi energy
        real uii_uee, uxcdcc                     ! short-range energies
        real vdW                                 ! van der Waal's energy

        real getot                               ! grand total energy
        
        interface
           subroutine Qmixer (t, iscf_iteration, sigma)
             use M_configuraciones
             use M_charges
             implicit none
             integer, intent (in) :: iscf_iteration
             type(T_structure), target :: t
             real, intent (inout) :: sigma
           end subroutine Qmixer

           subroutine writeout_energies (t, ebs, uii_uee, uxcdcc)
             use M_assemble_blocks
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t         ! the structure to be used
             real, intent (in) :: ebs               ! band-structure energy
             real, intent (in) :: uii_uee, uxcdcc   ! short-range energies
           end subroutine writeout_energies

           subroutine writeout_xyz (t, ebs, uii_uee, uxcdcc)
             use M_species
             use M_configuraciones
             implicit none
             type(T_structure), target :: t         ! the structure to be used
             real, intent (in) :: ebs               ! band-structure energy
             real, intent (in) :: uii_uee, uxcdcc   ! short-range energies
           end subroutine writeout_xyz
        end interface

! Allocate Arrays
! ===========================================================================
! None

! ===========================================================================
! Procedure
        write (s%logfile, '(A)') ' ---------------------------------------------------- '
        write (s%logfile, *) ' Begin energy and force recalculation for MDET '
        write (s%logfile, '(A)') ' ---------------------------------------------------- '

! ===========================================================================
! ---------------------------------------------------------------------------
!              S C F   L O O P
! ---------------------------------------------------------------------------
! Note that self-consistency is preformed regardless of method used.
! But, in Harris, we just do a single scf cycle.
! ===========================================================================
            call destroy_denmat (s)
            call destroy_assemble_ewald (s)
            call destroy_assemble_vxc (s)
            call destroy_assemble_vna (s)

            sigma = 999.0d0
            iscf_iteration = 1
            do while (sigma .gt. scf_tolerance_set .and.                     &
      &               iscf_iteration .le. max_scf_iterations_set - 1)
!             call cpu_time (time_scf_begin)
!             write (s%logfile, *)
!             write (s%logfile, '(A, I5, A7, I5, A1)') 'Self-Consistent Field step: ', &
!                  & iscf_iteration, ' (max: ', max_scf_iterations_set, ')'
!             write (s%logfile, '(A)') '----------------------------------------------------'
!             write (s%logfile, *)

!             write (s%logfile, *) ' Two-center charge dependent assemblers. '
              call assemble_vna_2c (s)
              call assemble_ewaldsr (s)
              call assemble_ewaldlr (s)

!             call cpu_time (time_initial)
!             write (s%logfile, *)
!              write (s%logfile, *) ' Three-center charge dependent assemblers. '
              call assemble_vna_3c (s)

!             write (s%logfile, *)
!             write (s%logfile, *) ' Exchange-correlation assemblers. '
              call assemble_vxc (s)

!             call cpu_time (time_final)
!             write (s%logfile, *) ' vna_3c, vxc time: ', time_final - time_initial

! ===========================================================================
! ---------------------------------------------------------------------------
!                         D I A G O N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Calculating the overlap matrix in K space
!             write (s%logfile, *)
!             write (s%logfile, *) ' Kspace '
!             call cpu_time (time_initial)
              call driver_kspace (s, iscf_iteration)
              call density_matrix (s, efermi)            
!             if (iwriteout_density .eq. 1) call writeout_density (s)

!             call cpu_time (time_final)
!             write (s%logfile, *) ' kspace time: ', time_final - time_initial
              if (ifix_CHARGES .ne. 1) then
                call calculate_charges (s)
                call Qmixer (s, iscf_iteration, sigma)
              end if
              if (iwriteout_charges .eq. 1) call writeout_charges (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!                       T O T A L   E N E R G I E S
! ---------------------------------------------------------------------------
! ===========================================================================
! short-range interactions (double-counting interactions)
              call calculate_ebs (s, ebs)
              uii_uee = 0.0d0; uxcdcc = 0.0d0
              call assemble_uee (s, uii_uee)
              call assemble_uxc (s, uxcdcc)
              ! Evaluate total energy
              etot = ebs + uii_uee + uxcdcc

              if (sigma .gt. scf_tolerance_set .and.                      &
      &           iscf_iteration .le. max_scf_iterations_set - 1 .and.    &
      &           ifix_CHARGES .ne. 1) then
!               write (s%logfile, *) ' Destroy some SCF arrays... '
                call destroy_denmat (s)
                call destroy_assemble_ewald (s)
                call destroy_assemble_vxc (s)
                call destroy_assemble_vna (s)
              end if
!             call cpu_time (time_scf_end)
!             write (s%logfile, *)
!             write (s%logfile, *) ' SCF ENERGY time: ', time_scf_end - time_scf_begin

! After building the density matrix, then we can free up ewald and denmat arrays
! - we reallocate these during the next SCF cycle anyways.
! We also free up the vna and vxc arrays if this is not converged.
              if (sigma .gt. 0.0d0) then
                iscf_iteration = iscf_iteration + 1
              else
                exit
              end if
              if (ifix_CHARGES .eq. 1) exit
            end do

!           call writeout_energies (s, ebs, uii_uee, uxcdcc)

! ===========================================================================
! ---------------------------------------------------------------------------
!                               F O R C E S
! ---------------------------------------------------------------------------
! ===========================================================================
!           call cpu_time (time_forces_begin)
            call initialize_forces (s)

! NAC Initialize NAC dij and NAC density            
            call initialize_nac (s)
            call density_matrix_nac (s)

            call writeout_density_nac (s)            
            call densityPP_matrix (s)
            call cape_matrix (s)

! After building the density matrix, then we can free up the kspace memory
            call destroy_kspace (s)

!           write (s%logfile, *)
!           write (s%logfile,'(A)') 'Forces '
!           write (s%logfile,'(A)') '------ '

! Assemble the derivative blocks needed for forces
!           write (s%logfile, *) ' Two-center non-charge dependent Dassemblers.'
            call Dassemble_S (s)
            call Dassemble_T (s)
            call Dassemble_dipole_z (s)
            call Dassemble_svnl (s)
            call Dassemble_vnl_2c (s)

!           write (s%logfile, *) ' Two-center charge dependent Dassemblers.'
            call Dassemble_vna_2c (s)
            call Dassemble_ewaldsr (s)
            call Dassemble_ewaldlr (s)

!           write (s%logfile, *)
!           write (s%logfile,*) ' Three-center non-charge dependent assemblers.'
            call Dassemble_vnl_3c (s)

!           call cpu_time (time_initial)
!           write (s%logfile,*) ' Three-center charge dependent Dassemblers.'
            call Dassemble_vna_3c (s)

!           write (s%logfile, *)
!           write (s%logfile, *) ' Exchange-correlation Dassemblers. '
            call Dassemble_vxc (s)
!           write (s%logfile, *) ' Three-center charge dependent Exchange-correlation Dassemblers. '

! NAC For three center part of vxc for dij
            call Dassemble_vxc_3c (s)
!           call cpu_time (time_final)
!           write (s%logfile, *) ' vxc, vna_3c forces time: ', time_final - time_initial

! short-range interactions (double-counting interactions)
            call Dassemble_uee (s)
            call Dassemble_uxc (s)

            call build_forces (s, rms)

            if (iwriteout_forces .eq. 1) call writeout_forces (s)
            
            write (s%logfile,*)
            write (s%logfile,*) ' Total Forces:'
            do iatom = 1, s%natoms
              write (s%logfile, 512)  iatom, s%forces(iatom)%ftot
            end do

!           call cpu_time (time_forces_end)
!           write (s%logfile, *)
!           write (s%logfile, *) ' FORCES time: ', time_forces_end - time_forces_begin

            write (s%logfile, *)
            write (s%logfile, '(A)') ' Grand Total Energy '
            write (s%logfile, '(A)') ' ------------------ '
 
            s%md%tkinetic = 0.0d0
            do iatom = 1, s%natoms
              in1 = s%atom(iatom)%imass            
              s%md%tkinetic = s%md%tkinetic                                                  &
         &               + (0.5d0/P_fovermp)*species(in1)%xmass                   &
         &                *(s%atom(iatom)%vatom(1)**2 +  s%atom(iatom)%vatom(2)**2 &
         &                                            +  s%atom(iatom)%vatom(3)**2)
            end do
            s%md%T_instantaneous = (2.0d0/3.0d0)*(s%md%tkinetic/s%natoms)*P_kconvert

            write (s%logfile,600) etot
            write (s%logfile,601) s%md%tkinetic
            write (s%logfile,602) vdW
            getot = etot + s%md%tkinetic + vdW
            write (s%logfile,603) getot
            write (s%logfile,604) getot/s%natoms
            write (s%logfile, *)

            write (s%logfile, '(A)') ' ---------------------------------------------------- '
            write (s%logfile, *) ' End energy and force recalculation for MDET"
            write (s%logfile, '(A)') ' ---------------------------------------------------- '

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
512     format (2x, 'f_total =',i6 ,3(2x,f15.6))

600     format (2x, '                                 Potential Energy = ', f18.8)
601     format (2x, '                           Nuclear Kinetic Energy = ', f18.8)
602     format (2x, '                           van der Waal''s Energy = ', f18.8)
603     format (2x, ' Grand Total Energy (Nuclear Kinetic + Potential) = ', f18.8)
604     format (2x, '                      Grand Total Energy per Atom = ', f18.8)
605     format (2x, '                               deltaE/atom  (meV) = ', f18.8)

! End Subroutine
! ===========================================================================
        return
        end subroutine recalculation

! destroy_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing denmat_mdet
!>       information.
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
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ikpoint                    !< counter of band and kpoint
        integer iband                     !< counter of transitions

! Procedure
! ===========================================================================
        ! destroy the density matrix pieces - forces are already evaluated
        do ikpoint = 1, s%nkpoints
          do iband = 1, s%kpoints(ikpoint)%nbands
            deallocate (s%kpoints(ikpoint)%transition(iband)%dij)
            deallocate (s%kpoints(ikpoint)%transition(iband)%c_mdet)          
            deallocate (s%kpoints(ikpoint)%transition(iband)%c_na)
            deallocate (s%kpoints(ikpoint)%transition(iband)%dij_old)
          end do
          deallocate (s%kpoints(ikpoint)%transition)
          deallocate (s%kpoints(ikpoint)%eigen_old)
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_mdet

! End Module
! ===========================================================================
        end module M_nonadiabatic_mdet
