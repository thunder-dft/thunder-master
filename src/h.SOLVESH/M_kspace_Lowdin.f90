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

! M_kspace.f90
! Module Description
! ===========================================================================
!      This is a version of kspace.f90 that uses the blas library.
!      It contains the following subroutines within the module:
!
!      driver_kspace - driver for the kspace routines.
!      writeoutME_SandH - writeout the overlap and Hamiltonina matrix elements
!      build_H - builds the Hamiltonian matrix according to neighbor
!                mapping representation.
!      kspace_Smatrix - calculates and diagonalizes the overlap matrix in
!                           k-space
!      kspace_Lowdin - calculates and diagonalizes S^1/2*H*S^1/2 (Lowdin)
!                          transformed matrix in k space.
! ============================================================================
        module M_kspace

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_species
        use M_configuraciones

! /SOLVESH
        use M_diagonalization

! Type declaration
! =========================================================================
! None

! module procedures
        contains


! ===========================================================================
! driver_kspace
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This subroutine is the driver for all of the kspace routines.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine driver_kspace (s, iscf_iteration)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

        integer, intent (in) :: iscf_iteration   !< which scf iteration?

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer logfile                     !< writing to which unit

        interface
          function phase (dot)
            real phase
            real dot
          end function
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        do ikpoint = 1, s%nkpoints
          call build_Hmatrix (s)

          if (iwriteout_ME_SandH .eq. 1) then
            write (logfile,*)
            write (logfile,*) ' Writing out overlap and Hamiltonian matrices: '
            call writeout_ME_SandH (s)
          end if

          call diagonalization_initialize (s, iscf_iteration, ikpoint)
          if (iscf_iteration .eq. 1) call kspace_Smatrix (s, ikpoint)
          call kspace_Lowdin (s, iscf_iteration, ikpoint)
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine driver_kspace


! ===========================================================================
! build_Hmatrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This subroutine builds the Hamiltonian by adding matrices of kinetic, Vna,
! ... etc. and stores to a T_neighbors type variable called Hmatrix(:).
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine build_Hmatrix (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer in1, in2                !< species number
        integer jatom, num_neigh        !< counters over neighbors
        integer norb_mu, norb_nu        !< size of the (mu, nu) block for pair

        type(T_assemble_block), pointer :: pH_neighbors
        type(T_assemble_neighbors), pointer :: pHamiltonian
        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic
        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna
        type(T_assemble_block), pointer :: pSR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldsr
        type(T_assemble_block), pointer :: pLR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldlr
        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc

! Allocate Arrays
! ===========================================================================
        allocate (s%Hamiltonian(s%natoms))

! Procedure
! ===========================================================================
! loop over atoms in central cell
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          nullify (pHamiltonian)
          pHamiltonian=>s%Hamiltonian(iatom)

          nullify (pkinetic, pvna, pewaldsr, pewaldlr, pvxc)
          pkinetic=>s%kinetic(iatom)
          pvna=>s%vna(iatom)
          pewaldsr=>s%ewaldsr(iatom)
          pewaldlr=>s%ewaldlr(iatom)
          pvxc=>s%vxc(iatom)

          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate(s%Hamiltonian(iatom)%neighbors(num_neigh))

! Now loop over all neighbors ineigh of iatom.
          do ineigh = 1, num_neigh
            ! cut some lengthy notation
            nullify (pH_neighbors)
            pH_neighbors=>pHamiltonian%neighbors(ineigh)

            nullify (pK_neighbors, pvna_neighbors)
            nullify (pSR_neighbors, pLR_neighbors, pvxc_neighbors)
            pK_neighbors=>pkinetic%neighbors(ineigh)
            pvna_neighbors=>pvna%neighbors(ineigh)
            pSR_neighbors=>pewaldsr%neighbors(ineigh)
            pLR_neighbors=>pewaldlr%neighbors(ineigh)
            pvxc_neighbors=>pvxc%neighbors(ineigh)

            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            allocate(s%Hamiltonian(iatom)%neighbors(ineigh)%block(norb_mu, norb_nu))
            pH_neighbors%block = pK_neighbors%block + pvna_neighbors%block     &
        &                       + pvna_neighbors%blocko + pvxc_neighbors%block &
        &                       - pSR_neighbors%block + pLR_neighbors%block
            nullify (pH_neighbors)
            nullify (pK_neighbors, pvna_neighbors)
            nullify (pSR_neighbors, pLR_neighbors, pvxc_neighbors)
          end do
          nullify (pHamiltonian)
          nullify (pkinetic, pvna, pewaldsr, pewaldlr, pvxc)
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine build_Hmatrix


! ===========================================================================
! writeoutME_SandH.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine is a utility to write out the components of the
! Hamiltonian for all the Harris interactions only.
!
! ===========================================================================
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine writeout_ME_SandH (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh          !< counters for atom, neighbor loops
        integer in1, in2               !< species numbers
        integer imu, inu               !< counters for mu, nu
        integer jatom                  !<
        integer jneigh                 !< neighbor counter for vnl
        integer katom, kbeta           !< PP neighbor atom and cell number
        integer logfile                !< writing to which unit
        integer mbeta                  !< the cell containing neighbor of iatom
        integer num_neigh              !< number of neighbors
        integer num_neighPPp           !< number of neighbors for vnl elements
        integer norb_mu, norb_nu       !< block size for the H and S block

        real z                         !< distance between two atoms

        real, dimension (3) :: r1, r2  !< positions for iatom, jatom
        real, dimension (3) :: sigma   !< direction along sigma bond

        ! for writing out the vna block combining the atom and ontop terms
        real, dimension (:, :), allocatable :: vna_block

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic
        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna

        ! exchange-correlation interactions
!       type(T_assemble_block), pointer :: prho_in_neighbors
!       type(T_assemble_neighbors), pointer :: prho_in
!       type(T_assemble_block), pointer :: prho_bond_neighbors
!       type(T_assemble_neighbors), pointer :: prho_bond
!       type(T_assemble_block), pointer :: pWrho_in_neighbors
!       type(T_assemble_neighbors), pointer :: pWrho_in
!       type(T_assemble_block), pointer :: pWrho_bond_neighbors
!       type(T_assemble_neighbors), pointer :: pWrho_bond
        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc

        type(T_assemble_block), pointer :: pvnl_neighbors
        type(T_assemble_neighbors), pointer :: pvnl

        ! long-range interactions
        type(T_assemble_block), pointer :: pSR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldsr
        type(T_assemble_block), pointer :: pLR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldlr

        type(T_assemble_block), pointer :: pH_neighbors
        type(T_assemble_neighbors), pointer :: pHamiltonian

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*)
        write (logfile,*) ' In writeoutME_SandH.f90 '
        write (logfile,*) ' Writing out pieces of the Hamiltonian. '
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          write (logfile,*) 
          num_neigh = s%neighbors(iatom)%neighn
          write (logfile,'(A, I5, A, I5, A1)') '### Atom:',  iatom, '  (Neighbors:', num_neigh,')'

          ! cut some lengthy notation
          poverlap=>s%overlap(iatom)
          pkinetic=>s%kinetic(iatom)
          pvna=>s%vna(iatom)

          ! vxc interactions
!         prho_in=>s%rho_in(iatom)
!         prho_bond=>s%rho_bond(iatom)
!         pWrho_in=>s%rho_in_weighted(iatom)
!         pWrho_bond=>s%rho_bond_weighted(iatom)
          pvxc=>s%vxc(iatom)

          ! long-range interactions
          pewaldsr=>s%ewaldsr(iatom)
          pewaldlr=>s%ewaldlr(iatom)

          pvnl=>s%vnl(iatom)
          pHamiltonian=>s%Hamiltonian(iatom)

          do ineigh = 1, num_neigh
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            in2 = s%atom(jatom)%imass
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a

            sigma = r2 - r1
            z = distance (r1, r2)
            write (logfile,*)
            write (logfile,102) ineigh, jatom, mbeta, z

! Overlap matrix elements
            write (logfile,*) 
            write (logfile,'(4x,A)') 'Interactions: s, t, vna, rho, and vxc '
            write (logfile,'(4x,A)') '===================================== '
            write (logfile,*) 
            write (logfile,'(4x,A)') 's: overlap '
            write (logfile,'(4x,A)') '---------- '
            ! cut some more lengthy notation
            pS_neighbors=>poverlap%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (logfile,104)                                            &
     &         (pS_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Kinetic energy matrix elements
            write (logfile,*) 
            write (logfile,'(4x,A)') 't: kinetic '
            write (logfile,'(4x,A)') '---------- '
            ! cut some more lengthy notation
            pK_neighbors=>pkinetic%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (logfile,104)                                            &
     &         (pK_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Neutral atom matrix elements
            write (logfile,*)
            write (logfile,'(4x,A)') 'vna: Hartree interactions '
            write (logfile,'(4x,A)') '-------------------------------- '
            ! cut some more lengthy notation
            pvna_neighbors=>pvna%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            allocate (vna_block (norb_mu, norb_nu))
            vna_block = pvna_neighbors%block + pvna_neighbors%blocko
            do imu = 1, norb_mu
              write (logfile,104) (vna_block(imu,inu), inu = 1, norb_nu)
            end do
            deallocate (vna_block)

! Exchange-correlation matrix elements - density
!           write (logfile,*)
!           write (logfile,'(4x,A)') 'rho: input density matrices (for vxc) '
!           write (logfile,'(4x,A)') '------------------------------------- '
            ! cut some more lengthy notation
!           prho_in_neighbors=>prho_in%neighbors(ineigh)
!           norb_mu = species(in1)%norb_max
!           norb_nu = species(in2)%norb_max
!           do imu = 1, norb_mu
!             write (logfile,104)                                            &
!    &         (prho_in_neighbors%block(imu,inu), inu = 1, norb_nu)
!           end do

! Exchange-correlation matrix elements - bond density
!           write (logfile,*)
!           write (logfile,'(4x,A)') 'rho: input (bond) density (for vxc)   '
!           write (logfile,'(4x,A)') '------------------------------------- '
            ! cut some more lengthy notation
!           prho_bond_neighbors=>prho_bond%neighbors(ineigh)
!           norb_mu = species(in1)%norb_max
!           norb_nu = species(in2)%norb_max
!           do imu = 1, norb_mu
!             write (logfile,104)                                              &
!    &         (prho_bond_neighbors%block(imu,inu), inu = 1, norb_nu)
!           end do

! Exchange-correlation matrix elements - input weighted density
!            write (logfile,*)
!            write (logfile,'(4x,A)') 'rhoS: input density matrices (for vxc) '
!            write (logfile,'(4x,A)') '------------------------------------- '
            ! cut some more lengthy notation
!            pWrho_in_neighbors=>prho_in%neighbors(ineigh)
!            norb_mu = species(in1)%norb_max
!            norb_nu = species(in2)%norb_max
!            do imu = 1, norb_mu
!              write (logfile,104)                                              &
!     &         (pWrho_in_neighbors%block(imu,inu), inu = 1, norb_nu)
!            end do

! Exchange-correlation matrix elements - bond density
!            write (logfile,*)
!            write (logfile,'(4x,A)') 'rhoS: input (bond) density (for vxc)   '
!            write (logfile,'(4x,A)') '------------------------------------- '
!            ! cut some more lengthy notation
!            pWrho_bond_neighbors=>prho_bond%neighbors(ineigh)
!            norb_mu = species(in1)%norb_max
!            norb_nu = species(in2)%norb_max
!            do imu = 1, norb_mu
!              write (logfile,104)                                            &
!     &         (pWrho_bond_neighbors%block(imu,inu), inu = 1, norb_nu)
!            end do

! Exchange-correlation matrix elements
            write (logfile,*)
            write (logfile,'(4x,A)') 'vxc: exchange-correlation interactions '
            write (logfile,'(4x,A)') '-------------------------------------- '
            ! cut some more lengthy notation
            pvxc_neighbors=>pvxc%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (logfile,104)                                            &
     &         (pvxc_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Short-range ewald matrix elements
            write (logfile,*)
            write (logfile,'(4x,A)') 'ewaldsr: short-range ewald interactions '
            write (logfile,'(4x,A)') '--------------------------------------- '
            ! cut some more lengthy notation
            pSR_neighbors=>pewaldsr%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (logfile,104)                                            &
     &         (pSR_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Long-range ewald matrix elements
            write (logfile,*)
            write (logfile,'(4x,A)') 'ewaldlr: long-range ewald interactions '
            write (logfile,'(4x,A)') '-------------------------------------- '
            ! cut some more lengthy notation
            pLR_neighbors=>pewaldlr%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (logfile,104)                                            &
     &         (pLR_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Complete Hamiltonian (without vnl) matrix elements
            write (logfile,*)
            write (logfile,'(4x,A)') 'Total Hamiltonian - without vnl piece '
            write (logfile,'(4x,A)') '------------------------------------- '
            ! cut some more lengthy notation
            pH_neighbors=>pHamiltonian%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (logfile,104)                                            &
     &         (pH_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! The vnl terms have a different neighbor mapping.  So, for the
! pseudo-potential part - we need to do a comparison - find the piece
! of the vnl matrix element that belongs to iatom, ineigh.
            num_neighPPp = s%neighbors_PPp(iatom)%neighn
            do jneigh = 1, num_neighPPp
              katom = s%neighbors_PPp(iatom)%neigh_j(jneigh)
              kbeta = s%neighbors_PPp(iatom)%neigh_b(jneigh)
              if (katom .eq. jatom .and. kbeta .eq. mbeta) then
                write (logfile,*)
                write (logfile,'(4x,A)') 'vnl: non-local pseudo-potential '
                write (logfile,'(4x,A)') '------------------------------- '
                ! cut some more lengthy notation
                pvnl_neighbors=>pvnl%neighbors(jneigh)
                norb_mu = species(in1)%norb_max
                norb_nu = species(in2)%norb_max
                do imu = 1, norb_mu
                  write (logfile,104)                                        &
     &             (pvnl_neighbors%block(imu,inu), inu = 1, norb_nu)
                end do
              end if
            end do
          end do
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (75('*'))
102     format (4x, 'Matrices connected to neighbor ', i3, ',',             &
     &          ' jatom = ', i4, ', ', ' mbeta = ', i4, ', ', ' d = ', f6.3)
103     format (75('='))
104     format (4x, 9f8.3)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_ME_SandH


! ===========================================================================
! kspace_Smatrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This subroutine calculates overlap matix S in k space:


!! s(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) s((0,iatom), (l,jatom))
!! and diagonalizes. The eigenvalues are stores in 1D array structure:
!! k_slam%values(:)
!
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
        subroutine kspace_Smatrix (s, ikpoint)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer in1, in2                 !< species numbers
        integer iatom, ineigh            !< counter over atoms and neighbors
        integer logfile                  !< writing to which unit
        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor
        integer imu, inu                 !< size of the block
        integer jatom                    !< neighbor of iatom
        integer norb_mu, norb_nu         !< size of the block for the pair
        integer jmu, jnu

        integer, intent(in) :: ikpoint
        real dot

        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! =========================================================================
! Initialize logfile
        logfile = s%logfile

! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms ! natoms passed
          poverlap=>s%overlap(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

! Now loop over all neighbors ineigh of iatom.
          do ineigh = 1, num_neigh ! <===== loop over i's neighbots
            pS_neighbors=>poverlap%neighbors(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

! Find the phase which is based on k*r
            vec = s%xl(mbeta)%a + s%atom(jatom)%ratom - s%atom(iatom)%ratom
            sks = s%kpoints(ikpoint)%k
            dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)

! So this matrix element goes in the i, j slot
            do inu = 1, norb_nu
              jnu = inu + s%iblock_slot(jatom)
              do imu = 1, norb_mu
                jmu = imu + s%iblock_slot(iatom)
                Smatrix(jmu,jnu) =                                          &
     &            Smatrix(jmu,jnu) + phase(dot)*pS_neighbors%block(imu,inu)
              end do ! do imu
            end do ! do inu
          end do ! do ineigh
        end do ! do iatom

! TESTING - Set Smatrix to identity
!       Smatrix = 0.0d0
!       do jmu = 1, s%norbitals
!         do jnu = 1, s%norbitals
!           if (jmu .eq. jnu) Smatrix(jmu,jnu) = 1.0d0
!         end do
!       end do

! DIAGONALIZE THE OVERLAP MATRIX
        call diagonalize_S (s)

        write (logfile,*)
        write (logfile,'(4x, A)') 'S: eigenvalues '
        write (logfile,'(4x, A)') '-------------- '
        write (logfile,200) (eigen(imu), imu = 1, s%norbitals_new)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (75('='))
200     format (4(2x, f12.4))

! End Subroutine
! ===========================================================================
        return
        end subroutine kspace_Smatrix


! ===========================================================================
! kspace_Lowdin
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This subroutine performs Lowdin transformation on the Hamiltonian
!! matix H in k space:
!! (S^-1/2)*H*(S^-1/2)
!! and diagonalizes it, and stores the eigenvalues in the structure
!! eigen_k(ikpoint)%values(:) or eigen%values(:)
!! @author Kh. Odbadrakh
!! @date March 18, 2008
! ===========================================================================
        subroutine kspace_Lowdin (s, iscf_iteration, ikpoint)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

        integer, intent (in) :: iscf_iteration   !< which scf iteration?
        integer, intent (in) :: ikpoint

! Parameters and Data Declaration
! ===========================================================================
! None

! Local variables declarations
! ===========================================================================
        integer in1, in2                !< species numbers
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer num_neigh               !< number of neighbors
        integer imu, inu                !< size of the block
        integer jmu, jnu                !< address of mu:nu block in big matrix
        integer jatom                   !< neighbor of iatom
        integer logfile                 !< writing to which unit
        integer mbeta
        integer norb_mu, norb_nu       !< size of the block for the pair

        real, dimension (3) :: sks  !< k point value
        real, dimension (3) :: vec
        real dot

        ! pointer to k-point dependent H matrix
        type(T_assemble_block), pointer :: pH_neighbors
        type(T_assemble_neighbors), pointer :: pHamiltonian

        type(T_assemble_block), pointer :: pvnl_neighbors
        type(T_assemble_neighbors), pointer :: pvnl

! Allocate Arrays
! ===========================================================================
! Arrays for the eigenvalues and eigenvectors
        allocate (s%kpoints(ikpoint)%eigen(s%norbitals))

        allocate (s%kpoints(ikpoint)%c(s%norbitals, s%norbitals))
        allocate (s%kpoints(ikpoint)%c_Lowdin(s%norbitals, s%norbitals))

! Procedure
! =========================================================================
! Initialize logfile
        logfile = s%logfile

! First, we build the n x n Hamiltonian matrix to be diagonalized based
! on all the neighbors interactions.  Second, we add in the vnl part.
! We do not add in the vnl part with everything else because the vnl
! interactions have a different neighbor map.

! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms ! natoms passed
          pHamiltonian=>s%Hamiltonian(iatom)
          pvnl=>s%vnl(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

! Now loop over all neighbors ineigh of iatom.
! Add in all the interactions into the n x n HUGE matrix for
! diagonalization.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh ! <===== loop over i's neighbots
            pH_neighbors=>pHamiltonian%neighbors(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! Find the phase which is based on k*r
            vec = s%xl(mbeta)%a + s%atom(jatom)%ratom - s%atom(iatom)%ratom
            sks = s%kpoints(ikpoint)%k
            dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)

            do inu = 1, norb_nu
              jnu = inu + s%iblock_slot(jatom)
              do imu = 1, norb_mu
                jmu = imu + s%iblock_slot(iatom)
                Hmatrix(jmu,jnu) =                                           &
     &           Hmatrix(jmu,jnu) + phase(dot)*pH_neighbors%block(imu,inu)
              end do ! do imu
            end do ! do inu
          end do ! do ineigh

! VNL PART:
          num_neigh = s%neighbors_PPp(iatom)%neighn
          do ineigh = 1, num_neigh
            pvnl_neighbors=>pvnl%neighbors(ineigh)
            jatom = s%neighbors_PPp(iatom)%neigh_j(ineigh)
            mbeta = s%neighbors_PPp(iatom)%neigh_b(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! Find the phase which is based on k*r
            vec = s%xl(mbeta)%a + s%atom(jatom)%ratom - s%atom(iatom)%ratom
            sks = s%kpoints(ikpoint)%k
            dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)

! So this matrix element goes in the i,j slot.
            do inu = 1, norb_nu
              jnu = inu + s%iblock_slot(jatom)
              do imu = 1, norb_mu
                jmu = imu + s%iblock_slot(iatom)
                Hmatrix(jmu,jnu) =                                           &
      &           Hmatrix(jmu,jnu) + phase(dot)*pvnl_neighbors%block(imu,inu)
              end do ! do imu
            end do ! do inu
          end do ! do ineigh
        end do ! do iatom

        call diagonalize_H_Lowdin (s, iscf_iteration, ikpoint)
        s%kpoints(ikpoint)%eigen = eigen

        write (logfile,*)
        write (logfile,'(4x, A)') 'H(k): eigenvalues '
        write (logfile,'(4x, A)') '----------------- '

        write (logfile,200) (s%kpoints(ikpoint)%eigen(imu), imu = 1, s%norbitals)

! We did a symmetric orthogonalization followed by a diagnalization of the
! Hamiltonian in this "MO" basis set. This yields a net canonical
! diagnolization with k-dependent coefficients s%kpoints(ikpoint)%c.

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (70('='))
200     format (4(2x, f12.4))

! End Subroutine
! ===========================================================================
        return
        end subroutine kspace_Lowdin


! ===========================================================================
! destroy_kspace
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the kspace matrix
!! information.
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
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_kspace (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh             !< counter over atoms and neighbors
        integer ikpoint                   !< counter over kpoints

! Procedure
! ===========================================================================
        deallocate (eigen)
        deallocate (Smatrix)

        do ikpoint = 1, s%nkpoints
          deallocate (s%kpoints(ikpoint)%S12matrix)
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
        end subroutine destroy_kspace

! End Module
! ===========================================================================
        end module M_kspace
