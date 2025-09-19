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
        integer itransition                !< counter over transitions
        integer ioccupy                    !< input occupation number
!       integer in1                        !< species number
        integer inpfile                    !< reading from which unit
        integer issh
        integer logfile                    !< writing to which unit
!       integer nfermi                     !< Fermi level state
        integer ntransitions               !< number of transitions
        integer ipop

!       real qztot                         !< total number of electrons
        real foccupy

        character (len = 25) :: slogfile

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! Initialize occupations
        do ikpoint = 1, s%nkpoints
          do ioccupy = 1, s%norbitals
            s%kpoints(ikpoint)%ioccupy(ioccupy) = 0
            s%kpoints(ikpoint)%foccupy(ioccupy) = 0.0d0
          end do
        end do

! Loop over the atoms.
! Total charge - ztot
!       s%ztot = 0.0d0
!       do iatom = 1, s%natoms
!         in1 = s%atom(iatom)%imass
!         do issh = 1, species(in1)%nssh
!           s%ztot = s%ztot + species(in1)%shell(issh)%Qneutral
!         end do
!       end do

!       qztot = s%ztot
!       nfermi = int(qztot)/2

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
        
! Allocate transition type and initialize imap
          read (inpfile,*) ntransitions
          s%kpoints(ikpoint)%ntransitions = ntransitions
          allocate (s%kpoints(ikpoint)%transition(ntransitions))

          do itransition = 1, ntransitions
            read (inpfile,*) iband, foccupy, ipop
! FIXME! We might need to read in foccupy and set ioccupy to 1 when foccupy
! is not 0 (since including all states no fixed relationship between ioccupy
! and foccupy)

            ! initialize imap
            s%kpoints(ikpoint)%transition(itransition)%imap = iband
            s%kpoints(ikpoint)%foccupy(iband) = foccupy

            ! initialize population
            s%kpoints(ikpoint)%transition(itransition)%cna = ipop
            ! (Note: might be errors with type diff. between ipop and cna)

            if (foccupy .ge. 0.5) s%kpoints(ikpoint)%ioccupy(iband) = 1
          end do   ! end loop over transitons
          write (logfile,*) ' testing imaps reach '
          write (logfile,*) s%kpoints(ikpoint)%transition(1)%imap
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
        integer itransition, jtransition   !< counter of transitions
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
        integer ntransitions             !< number of transitions

        real dot                         !< dot product between K and r
        real gutr                        !< real part of density matrix

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        logical read_occupy

        character (len = 25) :: slogfile

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat_mdet

! Allocate Arrays
! ===========================================================================
        allocate (s%denmat_mdet (s%natoms))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Loop over the atoms.
! Total charge - ztot
!       s%ztot = 0.0d0
!       do iatom = 1, s%natoms
!         in1 = s%atom(iatom)%imass
!         do issh = 1, species(in1)%nssh
!           s%ztot = s%ztot + species(in1)%shell(issh)%Qneutral
!         end do
!       end do

! Modify the total charge by the charge state
!       s%ztot = s%ztot + qstate

! ****************************************************************************
!
!                      C O M P U T E    D E N S I T I E S
! ****************************************************************************
        write (logfile,*)
        write (logfile,*) ' Calculating density matrix elements for '
        write (logfile,*) ' nonadiabatic coupling vectors based on transitions '
        write (logfile,*) ' from iband to jband. '

! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          nullify (pdenmat_mdet)
          pdenmat_mdet=>s%denmat_mdet(iatom)

! Allocate arrays
          allocate (s%denmat_mdet(iatom)%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            nullify (pRho_neighbors)
            pRho_neighbors=>pdenmat_mdet%neighbors(ineigh)

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (s%denmat_mdet(iatom)%neighbors(ineigh)%block(norb_mu, norb_nu))
            pRho_neighbors%block = 0.0d0

! Loop over the special k points.
            do ikpoint = 1, s%nkpoints

              ! Find the phase which is based on k*r
              vec = r2 - r1
              sks = s%kpoints(ikpoint)%k
              dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
              phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
              ntransitions = s%kpoints(ikpoint)%ntransitions

! Loop over all bands
              do itransition = 1, ntransitions
                iband = s%kpoints(ikpoint)%transition(itransition)%imap
                do jtransition = itransition + 1, ntransitions
                  jband = s%kpoints(ikpoint)%transition(jtransition)%imap
                  phase = phasex
                  do imu = 1, norb_mu
                    mmu = imu + s%iblock_slot(iatom)
                    step1 = phase*conjg(s%kpoints(ikpoint)%c(mmu,iband))
                    do jnu = 1, norb_nu
                      nnu = jnu + s%iblock_slot(jatom)
                      step2 = step1*s%kpoints(ikpoint)%c(nnu,jband)
                      gutr = real(step2)

! Finally the density matrix:
                      pRho_neighbors%block(imu,jnu) =                        &
     &                   pRho_neighbors%block(imu,jnu) + gutr
                    end do
                  end do
                end do
! Finish loop over bands.
              end do
! Finish loop over k-points.
            end do

! Finish loop over atoms and neighbors.
            nullify (pRho_neighbors)
          end do
          nullify (pdenmat_mdet)
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

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
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

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat_mdet

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*) ' '
        write (logfile,*) ' In writeout_density_mdet.f90 '
        write (logfile,*) ' Writing out pieces of the density matrix. '

! Open the file and write information.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.denmat_mdet'
        open (unit = 22, file = slogfile, status = 'unknown')
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          write (22,*)
          write (22,101)
          num_neigh = s%neighbors(iatom)%neighn
          write (22,*) ' There are ', num_neigh, ' neighbors to atom ', iatom
          write (22,101)

          ! cut some lengthy notation
          pdenmat_mdet=>s%denmat_mdet(iatom)
          do ineigh = 1, num_neigh
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            in2 = s%atom(jatom)%imass
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a

            sigma = r2 - r1
            z = distance (r1, r2)
            write (22,*) '  '
            write (22,102) ineigh, jatom, mbeta, z

! Density matrix elements
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat_mdet%neighbors(ineigh)
            write (22,103)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (22,104) (pRho_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do
          end do
        end do
        write (22,101)

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

! Procedure
! ===========================================================================
        ! destroy the density matrix pieces - forces are already evaluated
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%denmat_mdet(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%denmat_mdet(iatom)%neighbors)
        end do
        deallocate (s%denmat_mdet)

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
