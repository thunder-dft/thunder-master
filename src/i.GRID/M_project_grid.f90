! copyright info:
!
!                             @Copyright 2014
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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

! M_project_grid
! Module Description
! ===========================================================================
!       This is a module containing all information related to the developing
! grids - either for the grid version of the code to do Kohn-Sham calculations
! or for putting the charge densities on the grid for viewing. The following
! subroutines are called here:
!
!       assemble_Gmatrix.f90 - assemble the Kohn-Sham potentials vna and vxc
!       project_density_grid.f90 - project the density onto the grid
!       project_vna_grid.f90 - project the neutral atom potential onto the grid
!
! ===========================================================================
! Code written by:
! Prokop Hapala
! Pavel Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
!
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
! Module Declaration
! ===========================================================================
        module M_project_grid
        use M_species
        use M_configuraciones
        use M_grid
        use M_neighbors
        use M_atom_functions

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains


! ===========================================================================
! assemble_Gmatrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calculates the Kohn-Sham potentials for Hartree and
! exchange-correlation interactions.

!                 + X0 (iatom)
!                /|\     u1X = g1 - X0
! uX0 = X0- g0  / | \
!              /  |  + g1 (nearest grid point to iatom)
!             /   | /
!            /    |/
!           /     + Y0
!          +
!          g0 (origin of grid coords)
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    uX0 = X0 - g0
!    u1X = g1 - X0
!    r21 = Y0 - X0
!    u1Y = g1 - Y0 = g1 - Y0 - X0 + X0 = u1X - r21
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ============================================================================
! Module declaration
! ============================================================================
        subroutine assemble_Gmatrix (t)
        use M_species
        use M_configuraciones
        use M_atom_functions
        use M_grid
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh               !< counter over atoms and neighbors
        integer imesh                       !< counter over the mesh
        integer i0, j0, k0                  !< grid index points
        integer index0, index1, index2      !< different indexing counters
        integer in1, in2                    !< species numbers
        integer jatom                       !< neighbor of iatom
        integer logfile                     !< writing to which unit
        integer num_neigh                   !< number of neighbors
        integer mbeta                       !< the cell containing neighbor of iatom

        integer n1, l1, m1                  !< quantum numbers
        integer n2, l2, m2                  !< quantum numbers

        integer imu, inu                    !< loop over shells
        integer norb_mu, norb_nu            !< size of the block for the pair

!       integer, dimension (3) :: ipiv     !< index points on grid
        integer, dimension (3) :: nr

        real density
        real r1pin, r2pin                    ! distance from centers to integration point
        real r1max, r2max                    !< maximum extent of cutoff radius

        ! parameters for lda_ceperley_alder
        real exc, muxc, dexc, d2exc, dmuxc, d2muxc

        real, dimension (3) :: r1, r2, r21 !< positions of iatom and jatom
        real, dimension (3) :: r1p, r2p    !< shifted positions of iatom and jatom
        real, dimension (3) :: g1
        real, dimension (3) :: u, u1x
        real, dimension (3) :: x0
        real, dimension (3, 3) :: amatrix
        real, dimension (3, 3) :: ainverse

        real, allocatable, dimension (:) :: psi1
        real, allocatable, dimension (:) :: psi2

        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna
        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc

        interface
          function Ylm (r, l, m)
            integer, intent (in) :: l, m
            real, intent (in), dimension (3) :: r
            real Ylm
          end function Ylm
        end interface

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile
        write (logfile,*) ' Assemble Kohn-Sham potentials - Vxc(G) and Vna(G) '

! set nr(:)
        nr(1) = irm1
        nr(2) = irm2
        nr(3) = irm3

! We need to solve this linear eq.
!
!  | a1x  a2x  a3x |   |n1|   |x|
!  | a1y  a2y  a3y | x |n2| = |y|
!  ! a1z  a2z  a3z |   |n3|   |z|
!
        ! copy and invert original elementary lattice vectors
        ! to get form written above
        ainverse(1,:) = grid(1)%a
        ainverse(2,:) = grid(2)%a
        ainverse(3,:) = grid(3)%a
        amatrix = transpose(ainverse)

        ! inverse A: solving A*n=x -> n=A-1*x
        call invert3x3 (amatrix, ainverse)

! Loop over atoms
        do iatom = 1, t%natoms
          ! cut some lengthy notation
          pvxc=>t%vxc(iatom)
          pvna=>t%vna(iatom)

          r1 = t%atom(iatom)%ratom
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

! vector between the iatom (not centered in the unit cell yet) and
! the initial grid point
          u = t%atom(iatom)%ratom - g0
! get n-vector
          call mult3x1 (ainverse, u)

! round coefficients to get the position of the nearest grid point g1 to the iatom X1
! i,j,k can be positive or negative exceeding rmX (it means not centered in the unit cell)
          i0 = nint(u(1))
          j0 = nint(u(2))
          k0 = nint(u(3))

! find the vector u1 between the iatom X1 and the nearest point g1
          u1x(1) = u(1) - real(i0)
          u1X(2) = u(2) - real(j0)
          u1X(3) = u(3) - real(k0)

! check if the nearest grid point is located within the unit cell of the grid coords
! if not, let's map it within
!i0
          if (u(1) .lt. 0.0d0) then
            i0 = i0 + irm1*(int(abs(i0/irm1)) + 1)
          else
            i0 = i0 - irm1*int(i0/irm1)
          end if
!j0
          if (u(2) .lt. 0.0d0) then
            j0 = j0 + irm2*(int(abs(j0/irm2)) + 1)
          else
            j0 = j0 - irm2*int(j0/irm2)
          end if
!k0
          if (u(3) .lt. 0.0d0) then
            k0 = k0 + irm3*(int(abs(k0/irm3)) + 1)
          else
            k0 = k0 - irm3*int(k0/irm3)
          end if

! find the coordinates of the nearest point g1 witihin the grid coords
          g1(1) = i0*grid(1)%a(1) + j0*grid(2)%a(1) + k0*grid(3)%a(1)
          g1(2) = i0*grid(1)%a(2) + j0*grid(2)%a(2) + k0*grid(3)%a(2)
          g1(3) = i0*grid(1)%a(3) + j0*grid(2)%a(3) + k0*grid(3)%a(3)

! evaluate coordinates of the iatom in the grid coords
          x0(1) = g1(1) + u1x(1)*grid(1)%a(1) + u1x(2)*grid(2)%a(1) + u1x(3)*grid(3)%a(1)
          x0(2) = g1(2) + u1x(1)*grid(1)%a(2) + u1x(2)*grid(2)%a(2) + u1x(3)*grid(3)%a(2)
          x0(3) = g1(3) + u1x(1)*grid(1)%a(3) + u1x(2)*grid(2)%a(3) + u1x(3)*grid(3)%a(3)

! vector pointing from g1 to x0
          u1X = g1 - x0

! save iatom coordinate within the grid unit cell
          ratom2g(:,iatom) = x0(:)

! find index of the gX point within the extended mesh
          index0 = 1 + (i0 + iemx1) + iem1*(j0 + iemx2) + iem1*iem2*(k0 + iemx3)

! Loop over the neighbors
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pvxc_neighbors=>pvxc%neighbors(ineigh); pvxc_neighbors%block = 0.0d0
            pvna_neighbors=>pvna%neighbors(ineigh)

            mbeta = t%neighbors(iatom)%neigh_b(ineigh)
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            r2 = t%atom(jatom)%ratom + t%xl(mbeta)%a
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            r21 = r2 - r1

! Loop over points in the atomic mesh
            do imesh = 1, nam
              index1 = index0 + iam2rc(imesh)
              index2 = ie2r(index1) - 1
              r1p = ram2rc(:,imesh) + u1x
              r1pin = sqrt(r1p(1)**2 + r1p(2)**2 + r1p(3)**2)

! evaluate the vector between the jatom and the mesh point
              r2p = r1p - r21

! distance between the mesh point and jatom
              r2pin = sqrt(r2p(1)**2 + r2p(2)**2 + r2p(3)**2)

! check if jatom overlap with the gP mesh point
              if (r2pin .gt. rc_max) cycle

! Build wavefunction corresponding to iatom
! Loop over the orbital states of the atom - build wavefunction
              allocate (psi1 (species(in1)%norb_max)); psi1 = 0.0d0
              do imu = 1, norb_mu
                n1 = species(in1)%orbital(imu)%issh
                r1max = species(in1)%shell(n1)%rcutoffA
                l1 = species(in1)%orbital(imu)%l
                m1 = species(in1)%orbital(imu)%m

                psi1(imu) = psi1(imu) + psiofr(r1pin, r1max, in1, n1)*Ylm(r1p, l1, m1)
              end do

! Build wavefunction of ineigh, the neighbor of iatom
! Loop over the orbital states of the neighbor - build wavefunction
              allocate (psi2 (species(in2)%norb_max)); psi2 = 0.0d0
              do inu = 1, norb_nu
                n2 = species(in2)%orbital(inu)%issh
                r2max = species(in2)%shell(n2)%rcutoffA
                l2 = species(in2)%orbital(inu)%l
                m2 = species(in2)%orbital(inu)%m

                psi2(inu) = psi2(inu) + psiofr(r2pin, r2max, in2, n2)*Ylm(r2p, l2, m2)
              end do

              density = drhoG(index2) + rhoG0(index2)
              call lda_ceperley_alder (density, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
              vxcG(index2) = muxc

              do inu = 1, species(in2)%norb_max
                do imu = 1, species(in1)%norb_max
                  pvna_neighbors%block(imu,inu) =                            &
     &              pvna_neighbors%block(imu,inu)                            &
     &              + psi1(imu)*vcaG(index2)*psi2(inu)*dvolume
                  pvxc_neighbors%block(imu,inu) =                            &
     &              pvxc_neighbors%block(imu,inu)                            &
     &              + psi1(imu)*vxcG(index2)*psi2(inu)*dvolume
                end do ! do inu
              end do ! do imu
              deallocate (psi1, psi2)
            end do ! end loop over mesh
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Format Statements
! ===========================================================================
! None

        return
        end subroutine assemble_Gmatrix


! ===========================================================================
! project_density_grid
! ===========================================================================
! Subroutine Description
! ===========================================================================
! Project the density onto the grid.
!
!                 + x0 (iatom)
!                /|\     u1X = g1 - x0
! ux0 = x0- g0  / | \
!              /  |  + g1 (nearest grid point to iatom)
!             /   | /
!            /    |/
!           /     + Y0
!          +
!          g0 (origin of grid coords)
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    ux0 = x0 - g0
!    u1X = g1 - x0
!    r21 = Y0 - x0
!    u1Y = g1 - Y0 = g1 - Y0 - x0 + x0 = u1X - r21
! ===========================================================================
! Code written by:
! Prokop Hapala
! Pavel Jelenik
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
!
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
! Subroutine Declaration
! ===========================================================================
        subroutine project_density_grid (t)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh               !< counter over atoms and neighbors
        integer imesh                       !< counter over the mesh
        integer i0, j0, k0                  !< grid index points
        integer i, j, k
        integer index0, index1, index2      !< different indexing counters
        integer in1, in2                    !< species numbers
        integer jatom                       !< neighbor of iatom
        integer logfile                     !< writing to which unit
        integer num_neigh                   !< number of neighbors
        integer mbeta                       !< the cell containing neighbor of iatom

        integer n1, l1, m1                  !< quantum numbers
        integer n2, l2, m2                  !< quantum numbers

        integer imu, inu                    !< loop over shells
        integer norb_mu, norb_nu            !< size of the block for the pair

!        integer, dimension (3) :: ipiv     !< index points on grid
        integer, dimension (3) :: nr

        real adensity, density
        real r1pin, r2pin                   ! distance from centers to integration point
        real r1max, r2max                   !< maximum extent of cutoff radius
        real renorm
        real residual

        real, dimension (3) :: r1, r2, r21 !< positions of iatom and jatom
        real, dimension (3) :: r1p, r2p    !< shifted positions of iatom and jatom
        real, dimension (3) :: r_dipole    !< dipole value
        real, dimension (3) :: g1
        real, dimension (3) :: u, u1x
        real, dimension (3) :: x0
        real, dimension (3, 3) :: amatrix
        real, dimension (3, 3) :: ainverse

        real, target, allocatable, dimension (:) :: drhoG_save
        real, target, allocatable, dimension (:) :: rho

        real, allocatable, dimension (:) :: psi1
        real, allocatable, dimension (:) :: psi2

        character (len=25) xsfname
        character (len=25) message

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat

        real, dimension (:), pointer   :: pmat
        interface
          subroutine writeout_xsf (t, xsf, xsfname, message)
            use M_configuraciones
            use M_species
            use M_grid
            implicit none

            type(T_structure), target :: t           !< the structure to be used
            real, dimension (:), pointer :: xsf
            character (len=25), intent (in) :: xsfname
            character (len=25), intent (in) :: message

          end subroutine writeout_xsf
        end interface

        interface
          function Ylm (r, l, m)
            integer, intent (in) :: l, m
            real, intent (in), dimension (3) :: r
            real Ylm
          end function Ylm
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (drhoG_save (0:nrm-1))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile

! save previous step; reset drhoG to zero
        do index1 = 0,nrm-1
          drhoG_save(index1) = drhoG(index1)
        end do
        drhoG = 0.0d0

! set nr(:)
        nr(1) = irm1
        nr(2) = irm2
        nr(3) = irm3

! We need to solve this linear eq.
!
!  | a1x  a2x  a3x |   |n1|   |x|
!  | a1y  a2y  a3y | x |n2| = |y|
!  ! a1z  a2z  a3z |   |n3|   |z|
!
        ! copy and invert original elementary lattice vectors
        ! to get form written above
        ainverse(1,:) = grid(1)%a
        ainverse(2,:) = grid(2)%a
        ainverse(3,:) = grid(3)%a
        amatrix = transpose(ainverse)

        ! inverse A: solving A*n=x -> n=A-1*x
        call invert3x3 (amatrix, ainverse)

! Loop over atoms
        do iatom = 1, t%natoms
          ! cut some lengthy notation
          pdenmat=>t%denmat(iatom)

          r1 = t%atom(iatom)%ratom
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

! vector between the iatom (not centered in the unit cell yet) and
! the initial grid point
          u = t%atom(iatom)%ratom - g0
! get n-vector
          call mult3x1 (ainverse, u)

! round coefficients to get the position of the nearest grid point g1 to the iatom X1
! i,j,k can be positive or negative exceeding rmX (it means not centered in the unit cell)
          i0 = nint(u(1))
          j0 = nint(u(2))
          k0 = nint(u(3))

! find the vector u1 between the iatom X1 and the nearest point g1
          u1x(1) = u(1) - real(i0)
          u1X(2) = u(2) - real(j0)
          u1X(3) = u(3) - real(k0)

! check if the nearest grid point is located within the unit cell of the grid coords
! if not, let's map it within
!i0
          if (u(1) .lt. 0.0d0) then
            i0 = i0 + irm1*(int(abs(i0/irm1)) + 1)
          else
            i0 = i0 - irm1*int(i0/irm1)
          end if
!j0
          if (u(2) .lt. 0.0d0) then
            j0 = j0 + irm2*(int(abs(j0/irm2)) + 1)
          else
            j0 = j0 - irm2*int(j0/irm2)
          end if
!k0
          if (u(3) .lt. 0.0d0) then
            k0 = k0 + irm3*(int(abs(k0/irm3)) + 1)
          else
            k0 = k0 - irm3*int(k0/irm3)
          end if

! find the coordinates of the nearest point g1 witihin the grid coords
          g1(1) = i0*grid(1)%a(1) + j0*grid(2)%a(1) + k0*grid(3)%a(1)
          g1(2) = i0*grid(1)%a(2) + j0*grid(2)%a(2) + k0*grid(3)%a(2)
          g1(3) = i0*grid(1)%a(3) + j0*grid(2)%a(3) + k0*grid(3)%a(3)

! evaluate coordinates of the iatom in the grid coords
          x0(1) = g1(1) + u1x(1)*grid(1)%a(1) + u1x(2)*grid(2)%a(1) + u1x(3)*grid(3)%a(1)
          x0(2) = g1(2) + u1x(1)*grid(1)%a(2) + u1x(2)*grid(2)%a(2) + u1x(3)*grid(3)%a(2)
          x0(3) = g1(3) + u1x(1)*grid(1)%a(3) + u1x(2)*grid(2)%a(3) + u1x(3)*grid(3)%a(3)

! vector pointing from g1 to x0
          u1X = g1 - x0

! save iatom coordinate within the grid unit cell
          ratom2g(:,iatom) = x0(:)

! find index of the gX point within the extended mesh
          index0 = 1 + (i0 + iemx1) + iem1*(j0 + iemx2) + iem1*iem2*(k0 + iemx3)

! Loop over the neighbors
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)

            mbeta = t%neighbors(iatom)%neigh_b(ineigh)
            jatom = t%neighbors(iatom)%neigh_j(ineigh)
            r2 = t%atom(jatom)%ratom + t%xl(mbeta)%a
            in2 = t%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            r21 = r2 - r1

! Loop over points in the atomic mesh
            do imesh = 1, nam
              index1 = index0 + iam2rc(imesh)
              index2 = ie2r(index1) - 1
              r1p = ram2rc(:,imesh) + u1x
              r1pin = sqrt(r1p(1)**2 + r1p(2)**2 + r1p(3)**2)

! evaluate the vector between the jatom and the mesh point
              r2p = r1p - r21

! distance between the mesh point and jatom
              r2pin = sqrt(r2p(1)**2 + r2p(2)**2 + r2p(3)**2)

! check if jatom overlap with the gP mesh point
              if (r2pin .gt. rc_max) cycle

! Build wavefunction corresponding to iatom
! Loop over the orbital states of the atom - build wavefunction
              allocate (psi1 (species(in1)%norb_max)); psi1 = 0.0d0
              do imu = 1, norb_mu
                n1 = species(in1)%orbital(imu)%issh
                r1max = species(in1)%shell(n1)%rcutoffA
                l1 = species(in1)%orbital(imu)%l
                m1 = species(in1)%orbital(imu)%m

                psi1(imu) = psi1(imu) + psiofr(r1pin, r1max, in1, n1)*Ylm(r1p, l1, m1)
              end do

! Build wavefunction of ineigh, the neighbor of iatom
! Loop over the orbital states of the neighbor - build wavefunction
              allocate (psi2 (species(in2)%norb_max)); psi2 = 0.0d0
              do inu = 1, norb_nu
                n2 = species(in2)%orbital(inu)%issh
                r2max = species(in2)%shell(n2)%rcutoffA
                l2 = species(in2)%orbital(inu)%l
                m2 = species(in2)%orbital(inu)%m

                psi2(inu) = psi2(inu) + psiofr(r2pin, r2max, in2, n2)*Ylm(r2p, l2, m2)
              end do

! Assemble density
              density = 0.0d0
              do inu = 1, species(in2)%norb_max
                do imu = 1, species(in1)%norb_max
                  density = density                                          &
     &                     + pRho_neighbors%block(imu,inu)*psi1(imu)*psi2(inu)
                end do ! do inu
              end do ! do imu
              deallocate (psi1, psi2)

! store variation of density at given point
              drhoG(index2) = drhoG(index2) + density
            end do ! end loop over mesh
          end do ! end loop over neighbors
        end do ! end loop over atoms

! calculate the total denstity
        density = 0.0d0
        do index2 = 0, nrm - 1
          density = density + drhoG(index2)*dvolume
        end do
        write (logfile,*)
        write (logfile,*) ' -- Total atomic density before renormalization = ', density
        renorm = t%ztot/density

! check total density after renormalization
        density = 0.0d0
        do index2 = 0, nrm - 1
          drhoG(index2) = drhoG(index2)*renorm
          density = density + drhoG(index2)*dvolume
        end do
        write (logfile,*) ' -- Total atomic density after renormalization = ', density

! get drho (rest atomic charge)
        drhoG = drhoG - rhoG0

! evaluate residual of drhoG
        residual = 0.0d0
        do index2 = 0, nrm - 1
          residual = max(residual, abs(drhoG(index2) - drhoG_save(index2)))
        end do
        write (logfile, *) ' -- residual drhoG = ', residual, residual*dvolume

! also check if delta density goes to zero
        density = 0.0d0
        adensity = 0.0d0
        do index2 = 0, nrm - 1
          adensity = adensity + abs(drhoG(index2)*dvolume)
          density  = density + drhoG(index2)*dvolume
        end do
        write (logfile, *) ' -- Check sum drho should go to zero =', density
        write (logfile, *) ' -- |drho| =', adensity

! Calculate the dipole with the unit cell
        index2 = 0
        r_dipole = 0.0d0
        do k = 0, irm3 - 1
          do j = 0, irm2 - 1
            do i = 0, irm1 - 1
              r_dipole = r_dipole                                            &
     &                  + drhoG(index2)*(i*grid(1)%a + j*grid(2)%a + k*grid(3)%a)
              index2 = index2 + 1
            end do ! i
          end do ! j
        end do ! k
        r_dipole = r_dipole*dvolume
        write (logfile, *)
        write (logfile, 100) r_dipole/P_Debye

! Visualization in XCrysDen
! ===========================================================================
! Write out the xsf files (densities on 3D-plot grid)
        allocate (rho (0:nrm - 1))
        rho = drhoG + rhoG0

! write out total rho into xsf file
        pmat => rho
        xsfname = 'rho.xsf'
        message = '3D-plot of density'
        call writeout_xsf (t, pmat, xsfname, message)
        deallocate (rho)

! write out delta rho into xsf file
        pmat => drhoG
        xsfname = 'drho.xsf'
        message = '3D-plot of delta density'
        call writeout_xsf (t, pmat, xsfname, message)

! Deallocate Arrays
! ===========================================================================
        deallocate (drhoG_save)

! Format Statements
! ===========================================================================
100     format (2x, ' Dipole of system [D] = ', 3f14.6)

! End Subroutine
! ===========================================================================
        return
        end subroutine project_density_grid


! ===========================================================================
! project_vna_grid
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       Project the neutral atomic potential onto the mesh.
!
!                 + x0 (iatom)
!                /|\     u1X = g1 - x0
! ux0 = x0- g0  / | \
!              /  |  + g1 (nearest grid point to iatom)
!             /   | /
!            /    |/
!           /     + Y0
!          +
!          g0 (origin of grid coords)
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    ux0 = x0 - g0
!    u1X = g1 - x0
!    r21 = Y0 - x0
!    u1Y = g1 - Y0 = g1 - Y0 - x0 + x0 = u1X - r21
! ===========================================================================
! Code written by:
! Prokop Hapala
! Pavel Jelenik
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
!
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
! Subroutine Declaration
! ===========================================================================
        subroutine project_vna_grid (t)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over atoms and neighbors
        integer imesh                       !< counter over the mesh
        integer i0, j0, k0                  !< grid index points
!        integer i, j, k
        integer index0, index1, index2      !< different indexing counters
        integer in1                         !< species numbers
        integer issh                        !< counter over shells
        integer logfile                     !< writing to which unit

!        integer, dimension (3) :: ipiv     !< index points on grid
        integer, dimension (3) :: nr

        real r1pin                           ! distance from centers to integration point
!        real r1max                           !< maximum extent of cutoff radius
        real vna0

        real, dimension (3) :: r1          !< positions of iatom and jatom
        real, dimension (3) :: r1p         !< shifted positions of iatom and jatom
        real, dimension (3) :: g1
        real, dimension (3) :: u, u1x
        real, dimension (3) :: x0
        real, dimension (3, 3) :: amatrix
        real, dimension (3, 3) :: ainverse

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile
        write (logfile,*) ' Project vna onto the grid '

! re-initialize vnaG to zero
        vnaG = 0.0d0

! set nr(:)
        nr(1) = irm1
        nr(2) = irm2
        nr(3) = irm3

! We need to solve this linear eq.
!
!  | a1x  a2x  a3x |   |n1|   |x|
!  | a1y  a2y  a3y | x |n2| = |y|
!  ! a1z  a2z  a3z |   |n3|   |z|
!
        ! copy and invert original elementary lattice vectors
        ! to get form written above
        ainverse(1,:) = grid(1)%a
        ainverse(2,:) = grid(2)%a
        ainverse(3,:) = grid(3)%a
        amatrix = transpose(ainverse)

! Loop over atoms
        do iatom = 1, t%natoms
          r1 = t%atom(iatom)%ratom
          in1 = t%atom(iatom)%imass

! vector between the iatom (not centered in the unit cell yet) and
! the initial grid point
          u = t%atom(iatom)%ratom - g0
! get n-vector
          call mult3x1 (ainverse, u)

! round coefficients to get the position of the nearest grid point g1 to the iatom X1
! i,j,k can be positive or negative exceeding rmX (it means not centered in the unit cell)
          i0 = nint(u(1))
          j0 = nint(u(2))
          k0 = nint(u(3))

! find the vector u1 between the iatom X1 and the nearest point g1
          u1x(1) = u(1) - real(i0)
          u1X(2) = u(2) - real(j0)
          u1X(3) = u(3) - real(k0)

! check if the nearest grid point is located within the unit cell of the grid coords
! if not, let's map it within
!i0
          if (u(1) .lt. 0.0d0) then
            i0 = i0 + irm1*(int(abs(i0/irm1)) + 1)
          else
            i0 = i0 - irm1*int(i0/irm1)
          end if
!j0
          if (u(2) .lt. 0.0d0) then
            j0 = j0 + irm2*(int(abs(j0/irm2)) + 1)
          else
            j0 = j0 - irm2*int(j0/irm2)
          end if
!k0
          if (u(3) .lt. 0.0d0) then
            k0 = k0 + irm3*(int(abs(k0/irm3)) + 1)
          else
            k0 = k0 - irm3*int(k0/irm3)
          end if

! find the coordinates of the nearest point g1 witihin the grid coords
          g1(1) = i0*grid(1)%a(1) + j0*grid(2)%a(1) + k0*grid(3)%a(1)
          g1(2) = i0*grid(1)%a(2) + j0*grid(2)%a(2) + k0*grid(3)%a(2)
          g1(3) = i0*grid(1)%a(3) + j0*grid(2)%a(3) + k0*grid(3)%a(3)

! evaluate coordinates of the iatom in the grid coords
          x0(1) = g1(1) + u1x(1)*grid(1)%a(1) + u1x(2)*grid(2)%a(1) + u1x(3)*grid(3)%a(1)
          x0(2) = g1(2) + u1x(1)*grid(1)%a(2) + u1x(2)*grid(2)%a(2) + u1x(3)*grid(3)%a(2)
          x0(3) = g1(3) + u1x(1)*grid(1)%a(3) + u1x(2)*grid(2)%a(3) + u1x(3)*grid(3)%a(3)

! vector pointing from g1 to x0
          u1X = g1 - x0

! save iatom coordinate within the grid unit cell
          ratom2g(:,iatom) = x0(:)

! find index of the gX point within the extended mesh
          index0 = 1 + (i0 + iemx1) + iem1*(j0 + iemx2) + iem1*iem2*(k0 + iemx3)

! Loop over points in the atomic mesh
          do imesh = 1, nam
            index1 = index0 + iam2rc(imesh)
            index2 = ie2r(index1) - 1
            r1p = ram2rc(:,imesh) + u1x
            r1pin = sqrt(r1p(1)**2 + r1p(2)**2 + r1p(3)**2)

! get the vna potential
            do issh = 1, species(in1)%nssh
              vnaG(index2) = vnaG(index2) + vnaofr (r1pin, in1, issh)
            end do

            vnaG(index2) = vnaG(index2) + vna0
          end do ! do imesh
        end do ! end loop over atoms

! Convert to proper units
        vnaG = vnaG*P_eq2

! Format Statements
! ===========================================================================
! None

        return
        end subroutine project_vna_grid

! End Module
! ===========================================================================
        end module M_project_grid
