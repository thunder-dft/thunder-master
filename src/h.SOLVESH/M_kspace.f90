! copyright info:
!                             @Copyright 2013
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel JelinekmNZxbnmb

!
! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Washington University - Pete Fedders
! West Virginia University - Khorgolkhuu Odbadrakh, Ning Ma and Hao Wang
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
!       This is a version of kspace.f90 that uses the blas library.
!       It contains the following subroutines within the module:
!
!       driver_kspace - driver for the kspace routines.
!       writeoutME_SandH - writeout the overlap and Hamiltonina matrix elements
!       build_H - builds the Hamiltonian matrix according to neighbor
!                 mapping representation.
!       kspace_Smatrix - calculates and diagonalizes the overlap matrix in
!                            k-space
!       kspace_Hmatrix - calculates and diagonalizes the non-orthogonal
!                            Hamiltonian matrix in k space.
! ============================================================================
        module M_kspace
        use M_assemble_blocks
        use M_species
        use M_configuraciones
        use M_diagonalization

!       use M_assemble_vxc

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
!       This subroutine is the driver for all of the kspace routines.
!
! ===========================================================================
! Code written by:
! Prokop Hapala
! Department of Thin Films
! Institute of Physics
! Czech Academy of Sciences
! Prague, Czech Republic

! with modifications by:
! James P. Lewis
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
        integer logfile                     ! writing to which unit

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

          call diagonalization_initialize (s, iscf_iteration)
          call kspace_Smatrix (s, ikpoint)
          call kspace_Hmatrix (s, ikpoint)
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
          pHamiltonian=>s%Hamiltonian(iatom)
          pkinetic=>s%kinetic(iatom)
          pvna=>s%vna(iatom)
          pvxc=>s%vxc(iatom)

          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate(pHamiltonian%neighbors(num_neigh))

! Now loop over all neighbors ineigh of iatom.
          do ineigh = 1, num_neigh
            ! cut some lengthy notation
            pH_neighbors=>pHamiltonian%neighbors(ineigh)
            pK_neighbors=>pkinetic%neighbors(ineigh)
            pvna_neighbors=>pvna%neighbors(ineigh)

            pvxc_neighbors=>pvxc%neighbors(ineigh)

            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            allocate(pH_neighbors%block(norb_mu, norb_nu))

            pH_neighbors%block = pK_neighbors%block      & 
        &                      + pvna_neighbors%block    &
        &                      + pvxc_neighbors%block
          end do
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
! Code written by:
! Prokop Hapala
! Department of Thin Films
! Institute of Physics
! Czech Academy of Sciences
! Prague, Czech Republic

! with modifications by:
! James P. Lewis
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

        character (len = 25) :: slogfile

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic
        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna
        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc
        type(T_assemble_block), pointer :: pvnl_neighbors
        type(T_assemble_neighbors), pointer :: pvnl

        type(T_assemble_block), pointer :: pH_neighbors
        type(T_assemble_neighbors), pointer :: pHamiltonian

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*) ' '
        write (logfile,*) ' In writeoutME_SandH.f90 '
        write (logfile,*) ' Writing out pieces of the Hamiltonian. '

! Open the file and write information.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.SandH'
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
          poverlap=>s%overlap(iatom)
          pkinetic=>s%kinetic(iatom)
          pvna=>s%vna(iatom)
          pvxc=>s%vxc(iatom)
          pvnl=>s%vnl(iatom)
          pHamiltonian=>s%Hamiltonian(iatom)

          do ineigh = 1, num_neigh
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            in2 = s%atom(jatom)%imass
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a

            sigma = r2 - r1
            z = distance (r1, r2)
            write (22,*)
            write (22,102) ineigh, jatom, mbeta, z

! Overlap matrix elements
            write (22,'(4x,A)') 'Interactions: s, t, vna, and vxc  '
            write (22,'(4x,A)') '================================  '
            write (22,*) 
            write (22,'(4x,A)') 's: overlap  '
            write (22,'(4x,A)') '----------  '
            ! cut some more lengthy notation
            pS_neighbors=>poverlap%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (22,104) (pS_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Kinetic energy matrix elements
            write (22,*)
            write (22,'(4x,A)') 't: kinetic  '
            write (22,'(4x,A)') '----------  '
            ! cut some more lengthy notation
            pK_neighbors=>pkinetic%neighbors(ineigh)
            write (22,103)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (22,104) (pK_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Neutral atom matrix elements
            write (22,*)
            write (22,'(4x,A)') 'vna: Hartree interactions  '
            write (22,'(4x,A)') '-------------------------  '
            ! cut some more lengthy notation
            pvna_neighbors=>pvna%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (22,104) (pvna_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Exchange-correlation matrix elements
            write (22,*)
            write (22,'(4x,A)') 'vxc: exchange-correlation interactions  '
            write (22,'(4x,A)') '--------------------------------------  '
            ! cut some more lengthy notation
            pvxc_neighbors=>pvxc%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (22,104) (pvxc_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! Complete Hamiltonian (without vnl) matrix elements
            write (22,*)
            write (22,'(4x,A)') 'Total Hamiltonian  '
            write (22,'(4x,A)') '-----------------  '
            ! cut some more lengthy notation
            pH_neighbors=>pHamiltonian%neighbors(ineigh)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (22,104) (pH_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do

! The vnl terms have a different neighbor mapping.  So, for the
! pseudo-potential part - we need to do a comparison - find the piece
! of the vnl matrix element that belongs to iatom, ineigh.
            num_neighPPp = s%neighbors_PPp(iatom)%neighn
            do jneigh = 1, num_neighPPp
              katom = s%neighbors_PPp(iatom)%neigh_j(jneigh)
              kbeta = s%neighbors_PPp(iatom)%neigh_b(jneigh)
              if (katom .eq. jatom .and. kbeta .eq. mbeta) then
                write (22,*)
                write (22,'(4x,A)') 'vnl: non-local pseudo-potential  '
                write (22,'(4x,A)') '-------------------------------  '
                ! cut some more lengthy notation
                pvnl_neighbors=>pvnl%neighbors(jneigh)
                norb_mu = species(in1)%norb_max
                norb_nu = species(in2)%norb_max
                do imu = 1, norb_mu
                  write (22,104) (pvnl_neighbors%block(imu,inu), inu = 1, norb_nu)
                end do
              end if
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
102     format (1x, ' Matrices connected to neighbor ', i3, ','              &
     &          ' jatom = ', i4, ', ', ' mbeta = ', i4, ', ', ' d = ', f6.3)
103     format (75('='))
104     format (9f8.3)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_ME_SandH


! ===========================================================================
! kspace_Smatrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine calculates overlap matix S in k space:
! s(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) s((0,iatom), (l,jatom)).
!
! ===========================================================================
! Code written by:
! Prokop Hapala
! Department of Thin Films
! Institute of Physics
! Czech Academy of Sciences
! Prague, Czech Republic

! with modifications by:
! James P. Lewis
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
! kspace_Hmatrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine calculates Hamiltonian matrix H in k space:
! h(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) h((0,iatom), (l,jatom)).
! ===========================================================================
        subroutine kspace_Hmatrix (s, ikpoint)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Local variables declarations
! ===========================================================================
        integer in1, in2               !< species numbers
        integer iatom, ineigh          !< counter over atoms and neighbors
        integer logfile                !< writing to which unit
        integer num_neigh              !< number of neighbors
        integer imu, inu               !< size of the block
        integer jmu, jnu               !< address of mu:nu block in big matrix
        integer jatom                  !< neighbor of iatom
        integer mbeta
        integer norb_mu, norb_nu       !< size of the block for the pair

        integer, intent(in) :: ikpoint

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
! arrays for eigenvalues and occupation numbers
        allocate (s%kpoints(ikpoint)%eigen(s%norbitals))
        allocate (s%kpoints(ikpoint)%foccupy(s%norbitals))
        allocate (s%kpoints(ikpoint)%ioccupy(s%norbitals))

! array for the eigenvectors
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
          num_neigh = s%neighbors(iatom)%neighn

! Now loop over all neighbors ineigh of iatom.
! Add in all the interactions into the n x n HUGE matrix for
! diagonalization.
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
                Hmatrix(jmu,jnu) = Hmatrix(jmu,jnu) + phase(dot)*pH_neighbors%block(imu,inu)
              end do ! do imu
            end do ! do inu
          end do ! do ineigh

! VNL PART:
          do ineigh = 1, s%neighbors_PPp(iatom)%neighn
            pvnl_neighbors=>pvnl%neighbors(ineigh)
            jatom = s%neighbors_PPp(iatom)%neigh_j(ineigh)
            mbeta = s%neighbors_PPp(iatom)%neigh_b(ineigh)
            in2 = s%atom(jatom)%imass

            ! Find the phase which is based on k*r
            vec = s%xl(mbeta)%a + s%atom(jatom)%ratom - s%atom(iatom)%ratom
            sks = s%kpoints(ikpoint)%k
            dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)

! So this matrix element goes in the i,j slot.
            do inu = 1, species(in2)%norb_max
              jnu = inu + s%iblock_slot(jatom)
              do imu = 1, species(in1)%norb_max
                jmu = imu + s%iblock_slot(iatom)
                Hmatrix(jmu,jnu) = Hmatrix(jmu,jnu) + phase(dot)*pvnl_neighbors%block(imu,inu)
              end do ! do imu
            end do ! do inu
          end do ! do ineigh
        end do ! do iatom

        call diagonalize_H (s, ikpoint)
        s%kpoints(ikpoint)%eigen = eigen
        write (logfile,*)
        write (logfile,*) ' H(k): eigenvalues '
        write (logfile,100)
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
        end subroutine kspace_Hmatrix


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
        integer iatom                             !< counter over atoms
        integer ineigh                            !< counter over neighbors
        integer ikpoint                           !< counter over kpoints

! Procedure
! ===========================================================================
        deallocate (eigen)
        deallocate (Smatrix)
        deallocate (Hmatrix)

        do iatom=1, s%natoms
          do ineigh=1, s%neighbors(iatom)%neighn
            deallocate (s%Hamiltonian(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%Hamiltonian(iatom)%neighbors)
        end do
        deallocate (s%Hamiltonian)

        do ikpoint=1, s%nkpoints
          deallocate (s%kpoints(ikpoint)%eigen)
          deallocate (s%kpoints(ikpoint)%foccupy)
          deallocate (s%kpoints(ikpoint)%ioccupy)
          deallocate (s%kpoints(ikpoint)%c)
          deallocate (s%kpoints(ikpoint)%c_Lowdin)
        end do
        deallocate (s%kpoints)

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
