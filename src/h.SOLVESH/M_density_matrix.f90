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

! M_density_matrix.f90
! Program Description
! ===========================================================================
!       This routine calculates the density matrices and the band-structure
! energy, as well as similar density matrices.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ============================================================================
! Module declaration
! ============================================================================
        module M_density_matrix

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones
        use M_kpoints
        use M_kspace
        use M_rotations

! /FDATA
        use M_Fdata_2c

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! density_matrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine calculates the density matix rho and stores it in the
! structure given in M_assemble_block.f90
!
! ===========================================================================
        subroutine density_matrix (s, efermi)
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
        integer iband, ikpoint             !< counter of band and kpoint
        integer ihomo                      !< highest occupied level
        integer imu, inu
        integer in1, in2                   !< species numbers
        integer issh
        integer jatom                      !< neighbor of iatom
        integer logfile                    !< writing to which unit
        integer mmu, nnu

        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor
        integer norb_mu, norb_nu         !< size of the block for the pair

        real dot                         !< dot product between K and r
        real efermi
        real gutr                        !< real part of density matrix

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        logical read_occupy

        character (len = 25) :: slogfile

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat

! Allocate Arrays
! ===========================================================================
        allocate (s%denmat (s%natoms))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Loop over the atoms.
! Total charge - ztot
        s%ztot = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            s%ztot = s%ztot + species(in1)%shell(issh)%Qneutral
          end do
        end do

! Modify the total charge by the charge state
        s%ztot = s%ztot + qstate

        write (logfile,*)
        write (logfile,*) ' Calculating density matrix elements here. '
        write (logfile,*) ' The current charge state is, qstate = ', qstate
        write (logfile,*) ' Total number of electrons in the system is, ztot = ', s%ztot

! Calculate the Fermi energy.
! Inquire here regarding the occupations file
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.OCCUPATION'
        inquire (file = slogfile, exist = read_occupy)
        if (read_occupy) then
          call read_fermie (s)
        else
          call fermie (s, qstate, efermi)
        end if
        write (logfile,*) ' Fermi Level = ', efermi

! write out eigenvalues to eigen.dat
        if (modulo(int(s%ztot), 2) .eq. 0) then
          ihomo = int(s%ztot)/2
        else
          ihomo = int(s%ztot)/2 + 1
        end if
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.eigen'
        open (unit = 22, file = slogfile, status = 'replace')
        do ikpoint = 1, s%nkpoints
          write (22,100) ikpoint, efermi,                                    &
     &      s%kpoints(ikpoint)%eigen(ihomo+1) - s%kpoints(ikpoint)%eigen(ihomo)

! write out the eigenvalues
          write (22,200) (s%kpoints(ikpoint)%eigen(imu), imu = 1, s%norbitals)
        end do
        close (unit = 22)

! ****************************************************************************
!
!                      C O M P U T E    D E N S I T I E S
! ****************************************************************************
! Write out the coefficients to a file - *.cdcoeffs
        if (iwriteout_cdcoeffs .eq. 1) then

! Unformatted file
          slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
          slogfile = trim(slogfile)//'.cdcoeffs'
          open (unit = 22, file = slogfile, status = 'replace', form = 'unformatted')

          do ikpoint = 1, s%nkpoints
            do iband = 1, s%norbitals_new
              write (22) iband

! write out the coefficient
              write (22) (s%kpoints(ikpoint)%c(inu,iband), inu = 1, s%norbitals_new)
            end do
          end do
          close (unit = 22)

! Formatted file needed for Multimwfn
! Only for gamma
          slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
          slogfile = trim(slogfile)//'.cdcoeffs-mwfn'
          open (unit = 22, file = slogfile, status = 'replace')

          do ikpoint = 1, s%nkpoints
            write (22,"('Kpoint=',i10)") ikpoint
            do iband = 1, s%norbitals_new
              write (22,*)
              write (22,"('Index=',i10)") iband
              write (22,"('Type=',i2)") 0
              write (22,"('Energy=',1PE16.8)") s%kpoints(1)%eigen(iband)/P_Hartree
              write (22,"('Occ=',f12.8)") s%kpoints(1)%foccupy(iband) * P_spin
              write (22,"('Sym= ?')")
              write (22,"('$Coeff')")
! write out the coefficient
              write (22,"(5(1PE16.8))") (real(s%kpoints(ikpoint)%c(inu,iband)), inu = 1, s%norbitals_new)
            end do
            write (22,*)
          end do
         close (unit = 22)
        end if

! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          nullify (pdenmat)
          pdenmat=>s%denmat(iatom)

! Allocate arrays
          allocate (s%denmat(iatom)%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            nullify (pRho_neighbors)
            pRho_neighbors=>pdenmat%neighbors(ineigh)

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (s%denmat(iatom)%neighbors(ineigh)%block(norb_mu, norb_nu))
            pRho_neighbors%block = 0.0d0

! Loop over the special k points.
            do ikpoint = 1, s%nkpoints

              ! Find the phase which is based on k*r
              vec = r2 - r1
              sks = s%kpoints(ikpoint)%k
              dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
              phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight*P_spin

! Loop over all bands
              do iband = 1, s%norbitals_new
                if (s%kpoints(ikpoint)%ioccupy(iband) .ne. 0) then
                  phase = phasex*s%kpoints(ikpoint)%foccupy(iband)
                  do imu = 1, norb_mu
                    mmu = imu + s%iblock_slot(iatom)
                    step1 = phase*conjg(s%kpoints(ikpoint)%c(mmu,iband))
                    do inu = 1, norb_nu
                      nnu = inu + s%iblock_slot(jatom)
                      step2 = step1*s%kpoints(ikpoint)%c(nnu,iband)
                      gutr = real(step2)

! Finally the density matrix:
                      pRho_neighbors%block(imu,inu) =                        &
     &                  pRho_neighbors%block(imu,inu) + gutr
                    end do
                  end do
                end if

! Finish loop over bands.
              end do

! Finish loop over k-points.
            end do

! Finish loop over atoms and neighbors.
            nullify (pRho_neighbors)
          end do
          nullify (pdenmat)
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
        end subroutine density_matrix


! ===========================================================================
! cape_matrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine calculates the cape matix rho and stores it in the
! structure given in M_assemble_block.f90
!
!       The cape matrix is a*(iband,imu)*a(iband,inu)*eigen(iband) (see Eq. 47
! in Sankey_Niklewski paper). The coeficient is called a(iband,inu) but in this
! routine it is called c(nnu,iband)
! ===========================================================================
        subroutine cape_matrix (s)
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
        integer iband, ikpoint             !< counter of band and kpoint
        integer imu, inu
        integer in1, in2                   !< species numbers
        integer jatom                      !< neighbor of iatom
        integer logfile                    !< writing to which unit
        integer mmu, nnu

        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor
        integer norb_mu, norb_nu         !< size of the block for the pair

        real dot                         !< dot product between K and r
        real gutr                        !< real part of density matrix

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pCape_neighbors
        type(T_assemble_neighbors), pointer :: pcapemat

! Allocate Arrays
! ===========================================================================
        allocate (s%capemat (s%natoms))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! ****************************************************************************
!
!                 C O M P U T E   C A P E   D E N S I T I E S
! ****************************************************************************
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          nullify (pcapemat)
          pcapemat=>s%capemat(iatom)

! Allocate arrays
          allocate (s%capemat(iatom)%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            nullify (pCape_neighbors)
            pCape_neighbors=>pcapemat%neighbors(ineigh)

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (s%capemat(iatom)%neighbors(ineigh)%block(norb_mu, norb_nu))
            pCape_neighbors%block = 0.0d0

! Loop over the special k points.
            do ikpoint = 1, s%nkpoints
              ! Find the phase which is based on k*r
              vec = r2 - r1
              sks = s%kpoints(ikpoint)%k
              dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
              phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight*P_spin

! Loop over all bands
              do iband = 1, s%norbitals_new
                if (s%kpoints(ikpoint)%ioccupy(iband) .ne. 0) then
                  phase = phasex*s%kpoints(ikpoint)%foccupy(iband)
                  do imu = 1, norb_mu
                    mmu = imu + s%iblock_slot(iatom)
                    step1 = phase*conjg(s%kpoints(ikpoint)%c(mmu,iband))!

                    do inu = 1, norb_nu
                      nnu = inu + s%iblock_slot(jatom)
                      step2 = step1*s%kpoints(ikpoint)%c(nnu,iband)
                      gutr = real(step2)

! Finally the density matrix:
                      pCape_neighbors%block(imu,inu) =                         &
     &                  pCape_neighbors%block(imu,inu) + s%kpoints(ikpoint)%eigen(iband)*gutr
                    end do
                  end do
                end if

! Finish loop over bands.
              end do

! Finish loop over k-points.
            end do

! Finish loop over atoms and neighbors.
            nullify (pCape_neighbors)
          end do
          nullify (pcapemat)
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
        end subroutine cape_matrix


! ===========================================================================
! calculate_ebs.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine calculates the band-structure energy.
!
! ===========================================================================
        subroutine calculate_ebs (s, ebs)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Output
        real, intent (out) :: ebs      !< band-structure energy

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iband, ikpoint           !< counter of band and kpoint

        real ztest

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ****************************************************************************
!
!  C O M P U T E    B A N D - S T R U C T U R E    E N E R G Y
! ****************************************************************************
! Compute ebs, the band structure energy.
        ebs = 0.0d0
        ztest = 0.0d0
        do ikpoint = 1, s%nkpoints
          do iband = 1, s%norbitals_new
            if (s%kpoints(ikpoint)%ioccupy(iband) .eq. 1) then
              ebs = ebs + s%kpoints(ikpoint)%weight*P_spin*s%kpoints(ikpoint)%eigen(iband)  &
     &                                             *s%kpoints(ikpoint)%foccupy(iband)
              ztest = ztest + s%kpoints(ikpoint)%weight*P_spin*s%kpoints(ikpoint)%foccupy(iband)
            end if
          end do
        end do

! Test to make sure we get the proper number of states.
        if (abs(ztest - s%ztot) .gt. 1.0d-02) then
          write (*,*) ' *************** error *************** '
          write (*,*) ' ztest = ', ztest, ' ztot = ', s%ztot
          write (*,*) ' In denmat.f - ztest .ne. ztot! '
          stop
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
        end subroutine calculate_ebs


! ===========================================================================
! writeout_density.f90
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
        subroutine writeout_density (s)
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
        type(T_assemble_neighbors), pointer :: pdenmat

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*) ' '
        write (logfile,*) ' In writeout_density.f90 '
        write (logfile,*) ' Writing out pieces of the density matrix. '

! Open the file and write information.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.denmat'
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
          pdenmat=>s%denmat(iatom)
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
            pRho_neighbors=>pdenmat%neighbors(ineigh)
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
        end subroutine writeout_density


! ===========================================================================
! destroy_denmat
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing denmat
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
        subroutine destroy_denmat (s)
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
        ! destroy the coefficients now at the end of their use
        do ikpoint = 1, s%nkpoints
          deallocate (s%kpoints(ikpoint)%eigen)
          deallocate (s%kpoints(ikpoint)%c)
          deallocate (s%kpoints(ikpoint)%c_Lowdin)
        end do

        ! destory the Hamiltonian - we rebuild it
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%Hamiltonian(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%Hamiltonian(iatom)%neighbors)
        end do
        deallocate (s%Hamiltonian)

        ! destroy the density matrix pieces - forces are already evaluated
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%denmat(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%denmat(iatom)%neighbors)
        end do
        deallocate (s%denmat)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_denmat


! ===========================================================================
! destroy_denmat_PP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing denmat
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
        subroutine destroy_denmat_PP (s)
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

! Procedure
! ===========================================================================
        ! destroy the density matrix pieces - forces are already evaluated
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%capemat(iatom)%neighbors(ineigh)%block)
          end do
          do ineigh = 1, s%neighbors_PPp(iatom)%neighn
            deallocate (s%denmat_PP(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%denmat_PP(iatom)%neighbors)
          deallocate (s%capemat(iatom)%neighbors)
        end do
        deallocate (s%denmat_PP)
        deallocate (s%capemat)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_denmat_PP

! End Module
! ===========================================================================
        end module M_density_matrix
