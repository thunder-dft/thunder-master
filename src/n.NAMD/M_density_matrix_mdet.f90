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
        use M_assemble_blocks
        use M_configuraciones
        use M_kspace

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
        integer ioccupy                    !< input occupation number
        integer in1                        !< species number
        integer inpfile                    !< reading from which unit
        integer issh
        integer logfile                    !< writing to which unit
        integer nfermi                     !< Fermi level state
        integer ipop

        real qztot                         !< total number of electrons
        real foccupy

        character (len = 25) :: slogfile

        type(T_transition), pointer :: ptransition
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
          nullify (pkpoint)
          pkpoint => s%kpoints(ikpoint)
        
! Allocate transition type and initialize imap
          read (inpfile,*) pkpoint%nbands
          allocate (pkpoint%transition(pkpoint%nbands))
          do iband = 1, pkpoint%nbands

            ! cut some more lengthy notation
            nullify (ptransition)
            ptransition => pkpoint%transition(iband)

            read (inpfile,*) iband_in, foccupy, ipop
! FIXME! We might need to read in foccupy and set ioccupy to 1 when foccupy
! is not 0 (since including all states no fixed relationship between ioccupy
! and foccupy)

            ! initialize imap
            ptransition%imap = iband_in
            pkpoint%foccupy(iband) = foccupy

            ! initialize population
!           ptransition%cna = ipop
            ! (Note: might be errors with type diff. between ipop and cna)

            ! NAC Zhaofa Li initialize the dij
            allocate (ptransition%dij(3, pkpoint%nbands))
            ptransition%dij = 0.0d0
  
            if (foccupy .ge. 0.5d0) pkpoint%ioccupy(iband) = 1
            write (logfile,*) ' testing imaps reach '
            write (logfile,*) ptransition%imap
            nullify (ptransition)
          end do   ! end loop over bands
          nullify (pkpoint)
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
        integer iatom, ineigh              !< counter over atoms and neighbors
        integer iband, jband, ikpoint      !< counter of band and kpoint
        integer ihomo                      !< highest occupied level
        integer imu, inu, jnu
        integer in1, in2                   !< species numbers
        integer issh
        integer jatom                      !< neighbor of iatom
        integer logfile                    !< writing to which unit
        integer mmu, nnu

        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor
        integer norb_mu, norb_nu         !< size of the block for the pair
        integer nbands             !< number of bands

        real dot                         !< dot product between K and r
        real gutr                        !< real part of density matrix

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        logical read_occupy

        character (len = 25) :: slogfile

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: ptransition

! Allocate Arrays
! ===========================================================================
        allocate (s%denmat_mdet (s%natoms))

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
          nullify (pkpoint)
          pkpoint => s%kpoints(ikpoint)

! Loop over all bands
          nbands = pkpoint%nbands
          do iband = 1, nbands
            nullify (ptransition)
            ptransition => pkpoint%transition(iband)

            allocate (ptransition%c_mdet(s%norbitals))
            ptransition%c_mdet = pkpoint%c(:, ptransition%imap)
            nullify (ptransition)
! Finish loop over bands.
          end do
          nullify (pkpoint)
! Finish loop over k-points.
        end do
        
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
        integer iatom, ineigh          !< counters for atom, neighbor loops
        integer in1, in2               !< species numbers
        integer imu, inu               !< counters for mu, nu
        integer jatom                  !<
        integer logfile                !< writing to which unit
        integer mbeta                  !< the cell containing neighbor of iatom
        integer num_neigh              !< number of neighbors
        integer norb_mu, norb_nu       !< block size for the H and S block

        real z                         !< distance between two atoms

        real, dimension (3) :: r1, r2  !< positions for iatom, jatom
        real, dimension (3) :: sigma   !< direction along sigma bond

        character (len = 25) :: slogfile

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: ptransition

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================

! Formatted file needed for Multimwfn
! Only for gamma
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.cdcoeffs-mwfn_mdet'
        open (unit = 22, file = slogfile, status = 'replace')
        do ikpoint = 1, s%nkpoints
          write (22,"('Kpoint=',i10)") ikpoint

          nullify (pkpoint)
          pkpoint => s%kpoints(ikpoint)

          do iiband = 1, pkpoint%nbands

            nullify (ptransition)
            ptransition => pkpoint%transition(iiband)

            write (22,*)
            write (22,"('Index=',i10)") iiband
            write (22,"('Type=',i2)") 0
            write (22,"('Energy=',1PE16.8)") 0
            write (22,"('Occ=',f12.8)") 0
            write (22,"('Sym= ?')")
            write (22,"('$Coeff')")
! write out the coefficient
            write (22,"(5(1PE16.8))") (real(ptransition%c_mdet(inu)), inu = 1, s%norbitals)
            nullify (ptransition)
          end do
          write (22,*)
          nullify (pkpoint)
        end do
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
        integer iiband                !< counter of transitions
        integer nbands               !< number of transitions

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: ptransition

! Procedure
! ===========================================================================
        ! destroy the density matrix pieces - forces are already evaluated
        do ikpoint = 1, s%nkpoints
          do iiband = 1, s%kpoints(ikpoint)%nbands
            deallocate (s%kpoints(ikpoint)%transition(iiband)%c_mdet)
            deallocate (s%kpoints(ikpoint)%transition(iiband)%dij)
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
