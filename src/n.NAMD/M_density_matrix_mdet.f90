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

! M_density_matrix_mdet.f90
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
        module M_density_matrix_mdet

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
! None

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
        real foccupy

        character (len = 25) :: slogfile

        type(T_transition), pointer :: piband
        type(T_kpoint), pointer :: pkpoint

! Procedure
! ===========================================================================
! Initialize logfile
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
<<<<<<< HEAD
            nullify (piband); piband => pkpoint%transition(iband)
=======
            nullify (piband)
            piband => pkpoint%transition(iband)
>>>>>>> e6480782416e4cda261d48e447c394509396c483

            read (inpfile,*) iband_in, foccupy, ipop
! FIXME! We might need to read in foccupy and set ioccupy to 1 when foccupy
! is not 0 (since including all states no fixed relationship between ioccupy
! and foccupy)

            ! initialize imap
            piband%imap = iband_in
            pkpoint%foccupy(iband) = foccupy

            ! NAC Zhaofa Li initialize the dij and c_mdet
<<<<<<< HEAD
            allocate (piband%c_mdet(s%norbitals)); piband%c_mdet = 0.0d0
            allocate (piband%dij(3, s%natoms, pkpoint%nbands))
=======
            allocate (piband%c_mdet(s%norbitals))
            piband%c_mdet = 0.0d0
            allocate (piband%dij(3, pkpoint%nbands))
>>>>>>> e6480782416e4cda261d48e447c394509396c483
            piband%dij = 0.0d0
  
            if (foccupy .ge. 0.5d0) pkpoint%ioccupy(iband) = 1
            write (logfile,*) ' testing imaps reach '
            write (logfile,*) piband%imap
<<<<<<< HEAD
=======
            nullify (piband)
>>>>>>> e6480782416e4cda261d48e447c394509396c483
          end do   ! end loop over bands
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
! density_matrix_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine calculates the density matix rho and stores it in the
! structure given in M_assemble_block.f90 for molecular dynamics with electronic
! transition (MDET)
!
! ===========================================================================
        subroutine density_matrix_mdet (s)
        implicit none

        include '../include/constants.h'

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

        character (len = 25) :: slogfile

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
<<<<<<< HEAD

            ! cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            ! set the coefficient for the iband transition
            piband%c_mdet = pkpoint%c(:, piband%imap)

          end do   ! end loop over bands
        end do   ! end loop over kpoints
=======
            nullify (piband)
            piband => pkpoint%transition(iband)

            allocate (piband%c_mdet(s%norbitals))
            piband%c_mdet = pkpoint%c(:, piband%imap)
            nullify (piband)
! Finish loop over bands.
          end do
          nullify (pkpoint)
! Finish loop over k-points.
        end do
>>>>>>> e6480782416e4cda261d48e447c394509396c483
        
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
        end subroutine density_matrix_mdet


! ===========================================================================
! writeout_density_mdet.f90
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
        subroutine writeout_density_mdet (s)
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
          nullify (pkpoint)
          pkpoint => s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands

<<<<<<< HEAD
            ! cut some lengthy notation
            nullify (piband)
            piband => pkpoint%transition(iband)
=======
            nullify (piband)
            piband => pkpoint%transition(iiband)
>>>>>>> e6480782416e4cda261d48e447c394509396c483

            write (22,*)
            write (22,"('Index=',i10)") piband%imap
            write (22,"('Type=',i2)") 0
            write (22,"('Energy=',1PE16.8)") 0
            write (22,"('Occ=',f12.8)") 0
            write (22,"('Sym= ?')")
            write (22,"('$Coeff')")

! write out the coefficient
            write (22,"(5(1PE16.8))") (real(piband%c_mdet(inu)), inu = 1, s%norbitals)
<<<<<<< HEAD
          end do   ! end loop over bands
=======
            nullify (piband)
          end do
>>>>>>> e6480782416e4cda261d48e447c394509396c483
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
        end subroutine writeout_density_mdet

! ===========================================================================
! build_dij_mdet
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
        subroutine build_dij_mdet (s)
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

        real eigen_i, eigen_j             !< eigen values in band i and j
        real diff, tolnac

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
        tolnac = 0.0001d0
        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands

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
                write (s%logfile,*) ' NOT CALCULATED '
              else 
                piband%dij(:,:,jband) = piband%dij(:,:,jband)/(eigen_i - eigen_j)

                ! TESTING
                ! do iatom = 1, s%natoms
                !   write(s%logfile,'(5(A,I6))') "ikpoint=", ikpoint, ", iband=", iband, ", jband=", jband, ", iatom=", iatom
                !   write(s%logfile,*) piband%dij(:, iatom, jband)
                ! end do

                pjband%dij(:,:,iband) = pjband%dij(:,:,iband)/(eigen_i - eigen_j)

                ! TESTING
                ! do iatom = 1, s%natoms
                !   write(s%logfile,'(5(A,I6))') "ikpoint=", ikpoint, ", jband=", jband, ", iband=", iband, ", iatom=", iatom
                !   write(s%logfile,*) pjband%dij(:, iatom, iband)  
                ! end do 

              end if ! end if check for degeneracy
            end do ! end loop over jband
          end do ! end loop over iband
        end do  ! end loop over kpoints

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine build_dij_mdet

! ===========================================================================
! writeout_dij_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine writes out nonadiabatic coupliongs dij.
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
        subroutine writeout_dij_mdet (s)
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

          do iband = 1, pkpoint%nbands

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
        end subroutine writeout_dij_mdet

! ===========================================================================
! destroy_denmat_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing denmat_mdet
!! information.
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
        subroutine destroy_denmat_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh              !< counter over atoms and neighbors
        integer ikpoint                    !< counter of band and kpoint
<<<<<<< HEAD
        integer iband                     !< counter of transitions
=======
        integer iiband                     !< counter of transitions
        integer nbands                     !< number of transitions
>>>>>>> e6480782416e4cda261d48e447c394509396c483

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Procedure
! ===========================================================================
        ! destroy the density matrix pieces - forces are already evaluated
        do ikpoint = 1, s%nkpoints
          do iband = 1, s%kpoints(ikpoint)%nbands
            deallocate (s%kpoints(ikpoint)%transition(iband)%c_mdet)
            deallocate (s%kpoints(ikpoint)%transition(iband)%dij)
          end do
          deallocate (s%kpoints(ikpoint)%transition)
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
        end subroutine destroy_denmat_mdet

! End Module
! ===========================================================================
        end module M_density_matrix_mdet
