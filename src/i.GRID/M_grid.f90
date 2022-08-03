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

! M_grid
! Module Description
! ===========================================================================
!       This is a module containing all information related to the developing
! grids - either for the grid version of the code to do Kohn-Sham calculations
! or for putting the charge densities on the grid for viewing. The following
! subroutines are called here:
!
!       initialize_grid.f90 - build the mesh for the grid
!       intialize_project_grid.f90 - initialize the density projection
!                                   onto the grid
!       destroy_grid.f90 - destroy the arrays associated with the grid
!
! ===========================================================================
! Code written by:
! Prokop Hapala
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
        module M_grid
        use M_species
        use M_configuraciones
        use M_atom_functions

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        ! the regular mesh (spread out over the unit cell)
        integer :: irm1
        integer :: irm2
        integer :: irm3
        integer :: nrm                    !< total number mesh points

        ! the extended mesh (defined by the overlap of the
        ! atoms in neighbors periodic cells)
        integer :: iem1
        integer :: iem2
        integer :: iem3
        integer :: nem

        ! offset of extended mesh:  emi = cmi+ 2*mexi
        integer :: iemx1
        integer :: iemx2
        integer :: iemx3

        ! Map the extended mesh to the regular mesh
        integer, dimension (:), allocatable :: ie2r

        ! initial point of the grid
        real, dimension (3) :: g0

        ! unit (elementary) lattice vector of the submesh
        type (T_vector), dimension (3) :: grid
        type (T_vector), dimension (3) :: xsfgrid

        real dvolume                      !< elementary volume
        real rc_max                       !< maximum cutoff of all shells

        ! atomic mesh information
        integer nam
        integer, dimension (:), allocatable :: iam2rc
        real, dimension (:,:), allocatable :: ram2rc
        real, dimension (:,:), allocatable :: ratom2g

        ! Grid Interactions
        ! neutral atomic potential
        real, target, dimension (:), allocatable :: vnaG
        ! Hartree potential (related to density0)
        real, target, dimension (:), allocatable :: vcaG
        ! xc potential
        real, target, dimension (:), allocatable :: vxcG
        ! variation of density
        real, target, dimension (:), allocatable :: drhoG
        ! atomic density
        real, target, dimension (:), allocatable :: rhoG0

! module procedures
        contains


! ===========================================================================
! initialize_grid
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine initializes the grid mesh.
!
! ===========================================================================
! Code written by:
! Prokop Hapala
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
        subroutine initialize_grid (t)
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
        integer iatom                       !< counter over atoms
        integer in1                         !< species numbers
        integer ispecies                    !< counter over the species
        integer ix                          !< counter over x, y, and z
        integer logfile                     !< writing to which unit

        ! indexing point on mesh
        integer i, j, k
        integer ii, jj, kk
        integer index0, index1

        integer, dimension (3) :: ipiv     !< index points on grid
        integer, dimension (3) :: np

        real distance                       !< distance between mesh points
        real drmax

        real xmin, xmax
        real ymin, ymax
        real zmin, zmax

        real, dimension (3, 3) :: avec      !< temporary holding vector
        real, dimension (3) :: cvec         !< temporary holding vector
        real, dimension (3) :: r1           !< positions of iatom and jatom

        interface
          function a_cross_b (a, b)
            real, dimension (3) :: a_cross_b
            real, intent(in), dimension (3) :: a, b
          end function a_cross_b
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile

! find rc_max
        rc_max = 0.0d0
        do ispecies = 1, nspecies
          if (wf(ispecies)%rcutoffA_max .gt. Rc_max) Rc_max = wf(ispecies)%rcutoffA_max
        end do

! initialize the grid lattice vectors
        grid(1)%a = t%lattice(1)%a
        grid(2)%a = t%lattice(2)%a
        grid(3)%a = t%lattice(3)%a

! Calculate the geometric center of the system
        t%gcm = 0.0d0
        do iatom = 1, t%natoms
          t%gcm = t%gcm + t%atom(iatom)%ratom
        end do
        t%gcm = t%gcm/t%natoms
        write (logfile, *) ' System geometric center, shift according to: '
        write (logfile, 101) t%gcm

! Determine the maximum and minimum coordinate values.
! If we are doing a cluster, then we simply take the maximum extent
! of the atoms wavefunctions and integrate over this maximum extent.
! If we are doing a periodic cell, then we choose the minimum
! points as the coordinate of the first atom in the cell and then
! integrate along the lattice vectors.
        if (t%icluster .eq. 1) then
          ! Initialize to some very large numbers
          xmin = 999.0d0; xmax = -999.0d0
          ymin = 999.0d0; ymax = -999.0d0
          zmin = 999.0d0; zmax = -999.0d0

          do iatom = 1, t%natoms
            r1 = t%atom(iatom)%ratom
            in1 = t%atom(iatom)%imass
            xmin = min(xmin, r1(1) - wf(in1)%rcutoffA_max)
            xmax = max(xmax, r1(1) + wf(in1)%rcutoffA_max)
            ymin = min(ymin, r1(2) - wf(in1)%rcutoffA_max)
            ymax = max(ymax, r1(2) + wf(in1)%rcutoffA_max)
            zmin = min(zmin, r1(3) - wf(in1)%rcutoffA_max)
            zmax = max(zmax, r1(3) + wf(in1)%rcutoffA_max)
          end do  ! end loop over atoms

          ! Add small buffer
          xmin = xmin - 0.1d0; xmax = xmax + 0.1d0
          ymin = ymin - 0.1d0; ymax = ymax + 0.1d0
          zmin = zmin - 0.1d0; zmax = zmax + 0.1d0

          write (logfile, *)
          write (logfile, *) ' grid box size = '
          write (logfile, *) ' xmin, xmax = ', xmin, xmax
          write (logfile, *) ' ymin, ymax = ', ymin, ymax
          write (logfile, *) ' zmin, zmax = ', zmin, zmax

          ! setup lattice vector
          grid(1)%a = 0.0d0; grid(2)%a = 0.0d0; grid(3)%a = 0.0d0
          grid(1)%a(1) = xmax - xmin
          grid(2)%a(2) = ymax - ymin
          grid(3)%a(3) = zmax - zmin

          ! find the reciprocal lattice vectors related to grid and new cell volume
          ! First get a2 X a3.
          cvec = a_cross_b (grid(2)%a, grid(3)%a)

          ! Next find the volume of the cell.
          ! NOTE: t%volume actually has a sign below. At this point the sign is
          ! important since we form g vectors by dividing by a1 dot (a2 X a3).
          ! Oh, you say. what difference does it make if we change the sign of g.
          ! it makes no difference in principle.
          t%volume = grid(1)%a(1)*cvec(1) + grid(1)%a(2)*cvec(2)         &
     &              + grid(1)%a(3)*cvec(3)
          t%g(1)%a = 2.0d0*pi*cvec(:)/t%volume

! Next we get a3 X a1, and g2.
          cvec = a_cross_b (grid(3)%a, grid(1)%a)
          t%g(2)%a = 2.0d0*pi*cvec(:)/t%volume

! Finally we get a1 X a2, and g3.
          cvec = a_cross_b (grid(1)%a, grid(2)%a)
          t%g(3)%a = 2.0d0*pi*cvec(:)/t%volume
          t%volume = abs(t%volume)
          write (logfile,*)
          write (logfile,*) ' Changed cell volume [cubic Angstroms] = ', t%volume

          ! set the initial grid point
          g0(1) = xmin
          g0(2) = ymin
          g0(3) = zmin
        else
          ! set the initial grid point
          g0(1) = t%gcm(1) - 0.5d0*(grid(1)%a(1) + grid(2)%a(1) + grid(3)%a(1))
          g0(2) = t%gcm(2) - 0.5d0*(grid(1)%a(2) + grid(2)%a(2) + grid(3)%a(2))
          g0(3) = t%gcm(3) - 0.5d0*(grid(1)%a(3) + grid(2)%a(3) + grid(3)%a(3))
        end if

! Gauss diagonalization with pivoting
        write (logfile,*)
        write (logfile,*) ' Lattice vectors for grid integration: '
        do ix = 1, 3
          write (logfile,*) grid(ix)%a
          xsfgrid(ix)%a = grid(ix)%a
        end do

        write (logfile,*)
        write (logfile,*)  ' Reciprocal lattice vector [inverse Angstrom]: '
        do ix = 1, 3
          write (logfile,*) t%g(ix)%a
        end do

! Calculate the number of mesh points along an axis
        drmax = 0.0d0
        do ix = 1, 3
          cvec(ix) = sqrt(t%g(ix)%a(1)**2 + t%g(ix)%a(2)**2 + t%g(ix)%a(3)**2)
          ! estimate the number of points along an axis
          np(ix) = int (2.0d0*sqrt(Ecut)/(cvec(ix)*P_abohr) + 1)
          ! iterate until the number of points is an even multiplier
          do
            ! number of points must be multiply of 2!!
            if (mod(np(ix),2) .eq. 0) exit
            np(ix) = np(ix) + 1
          end do

          ! in reciprocal space one has to multiply by np
          grid(ix)%a = (t%g(ix)%a*np(ix))/(2.0d0*pi)

          ! the elementary distance at given axis
          cvec(ix) = cvec(ix)/np(ix)
          if (drmax .gt. cvec(ix)) drmax = cvec(ix)
        end do ! end loop over coordinates
        write (logfile,*)
        write (logfile,*)  ' Elementary lattice vector [Angstrom]:'
        do ix = 1, 3
          write (logfile,*) grid(ix)%a
        end do

! Find direct elementary vectors
        cvec = a_cross_b (grid(2)%a, grid(3)%a)
        dvolume = grid(1)%a(1)*cvec(1) + grid(1)%a(2)*cvec(2) + grid(1)%a(3)*cvec(3)
        avec(1,:) = cvec/dvolume

        cvec = a_cross_b (grid(3)%a, grid(1)%a)
        avec(2,:) = cvec/dvolume

        cvec = a_cross_b (grid(1)%a, grid(2)%a)
        avec(3,:) = cvec/dvolume

        write (logfile,*)
        write (logfile,*) ' Elementary grid lattice vector :'
        grid(1)%a = avec(1,:)
        grid(2)%a = avec(2,:)
        grid(3)%a = avec(3,:)
        write (logfile,101) grid(1)%a
        write (logfile,101) grid(2)%a
        write (logfile,101) grid(3)%a

! store the division of the regular mesh along the axis
        irm1 = np(1)
        irm2 = np(2)
        irm3 = np(3)
        ! total number of points on the regular mesh
        nrm = irm1*irm2*irm3
        write (logfile,*)
        write (logfile,*) ' Regular mesh information: ', irm1, irm2, irm3, nrm

! allocate arrays
        write (logfile,*) ' Allocating arrays for grid interactions. '
        allocate (vnaG (0: nrm - 1)); vnaG = 0.0d0
        allocate (vxcG (0: nrm - 1)); vxcG = 0.0d0
        allocate (vcaG (0: nrm - 1)); vcaG = 0.0d0
        allocate (rhoG0 (0:nrm - 1)); rhoG0 = 0.0d0
        allocate (drhoG (0:nrm - 1)); drhoG = 0.0d0

! elementary volume (inverse of previous calculated dvolume)
        dvolume = 1/dvolume

! Now evaluate the extended mesh at each direction
        do ix = 1, 3
          ! get the minimal distance
          distance = sqrt(grid(ix)%a(1)**2 + grid(ix)%a(2)**2 + grid(ix)%a(3)**2)
          ! calculate the threshold
          ipiv(ix) = rc_max/distance + 1
          ! total points in extended mesh
          np(ix) = np(ix) + 2*ipiv(ix)
        end do

! save offset of the extended mesh
        iemx1 = ipiv(1)
        iemx2 = ipiv(2)
        iemx3 = ipiv(3)

! save division of the extended mesh
        iem1 = np(1)
        iem2 = np(2)
        iem3 = np(3)
        ! total number of points on the extended mesh
        nem = iem1*iem2*iem3
        write (logfile,*) ' Extended mesh information: ', iem1, iem2, iem3, nem

! ***************************************************************************
!                 Map the Extended mesh to the Regular grid
! ***************************************************************************
        write (logfile,*) ' Map the extended mesh to the regular grid  '
        allocate (ie2r(nem))

! Loops over each axis
        do k = 0, iem3 - 1
          do j = 0, iem2 - 1
            do i = 0, iem1 - 1

              ! rest the offset
              ii = i - iemx1
              jj = j - iemx2
              kk = k - iemx3

              ! find the point in the regular mesh
              ii = mod(ii + 100*irm1, irm1)
              jj = mod(jj + 100*irm2, irm2)
              kk = mod(kk + 100*irm3, irm3)

              ! calculate index of the extended mesh point
              index1 = 1 + i + iem1*j + iem1*iem2*k
              ! calculate index of the regular mesh point
              index0 = 1 + ii + irm1*jj + irm1*irm2*kk

              ! save the value
              ie2r(index1) = index0
            end do ! do i
          end do ! do j
        end do ! do k

! ***************************************************************************
!                            Set up the atomic mesh
! ***************************************************************************
        write (logfile,*) ' Calculate the atomic mesh '

! first determine the array sizes
        nam = 0
        do k = -iemx3, iemx3
          do j = -iemx2, iemx2
            do i = -iemx1, iemx1
              ! distance to the mesh cell
              cvec = grid(1)%a*i + grid(2)%a*j + grid(3)%a*k
              distance = sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2)
              ! distance from point to mesh cell
              if ((rc_max + drmax) .gt. distance) nam = nam + 1
            end do
          end do
        end do

! allocate arrays related to atomic mesh
        allocate (iam2rc (nam))
        allocate (ram2rc (3, nam))
        allocate (ratom2g (3, t%natoms))

        index0 = 0
        do k = -iemx3, iemx3
          do j = -iemx2, iemx2
            do i = -iemx1, iemx1
              ! distance to the mesh cell
              cvec = grid(1)%a*i + grid(2)%a*j + grid(3)%a*k
              distance = sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2)
              if ((rc_max + drmax) .gt. distance) then
                index0 = index0 + 1
                iam2rc(index0) = i + iem1*j + iem1*iem2*k
                ram2rc(:,index0) = cvec
              end if
            end do
          end do
        end do
        write (logfile,*) ' Atomic mesh: ', 2*iemx1 + 1, 2*iemx2 + 1, 2*iemx3 + 1, nam

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (2x, 3f10.4)

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_grid


! ===========================================================================
! initialize_project_grid
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine initializes the density to the grid.
!
!                 + x0 (iatom)
!                / \     u1X = g1 - x0
! ux0 = x0- g0  /   \
!              /     + g1 (nearest grid point to iatom)
!             /
!            /
!           +
!          g0 (origin of grid coords)
!
! g1 ... point where the origin of atomic mesh (sphere) is placed
!  vectors we need:
!    x0g0
!    u1X = g1 - x0
! ===========================================================================
! Code written by:
! Prokop Hapala
! Pavel Jelenik
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_project_grid (t)
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
        integer iatom                       !< counter over the atoms
        integer imesh                       !< counter over the mesh
        integer i0, j0, k0
        integer in1                         !< species numbers
        integer index0, index1, index2      !< different indexing points
        integer logfile                     !< writing to which unit

        integer n1, l1, m1                  ! quantum numbers

        integer imu                         ! loop over shells
        integer norb_mu                     ! number atomic orbitals

        integer, dimension (3) :: nr

        real density
        real Qneutral              ! neutral charge on atom
        real qtot
        real r1pin                 ! distance from centers to integration point
        real r1max                          !< maximum extent of cutoff radius
        real renorm

        real, dimension (3) :: r1           !< position of iatom
        real, dimension (3) :: r1p          !< shifted positions of iatom - mesh

        real, dimension (3) :: g1
        real, dimension (3) :: u, u1x
        real, dimension (3) :: x0
        real, dimension (3, 3) :: xlmatrix
        real, dimension (3, 3) :: xlinverse

        interface
          function Ylm (r, l, m)
            integer, intent (in) :: l, m
            real, intent (in), dimension (3) :: r
            real Ylm
          end function Ylm
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile

! re-initialize variables to zero
        rhoG0 = 0.0d0
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
        xlinverse(1,:) = grid(1)%a
        xlinverse(2,:) = grid(2)%a
        xlinverse(3,:) = grid(3)%a
        xlmatrix = transpose(xlinverse)

        ! inverse A: solving A*n=x -> n=A-1*x
        call invert3x3 (xlmatrix, xlinverse)

! integration checking
        renorm = 0.0d0

! Loop over atoms
        do iatom = 1, t%natoms
          r1 = t%atom(iatom)%ratom
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

! vector between the iatom (not centered in the unit cell yet) and
! the initial grid point
          u = t%atom(iatom)%ratom - g0
! get n-vector
          call mult3x1 (xlinverse, u)

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

! save iatom coord within the grid unit cell
          ratom2g(:,iatom) = x0(:)

! find index of the gX point within the extende mesh
          index0 = 1 + (i0 + iemx1) + iem1*(j0 + iemx2) + iem1*iem2*(k0 + iemx3)

! skip certain atoms to avoid the overcounting; it means all jatoms having
! the identification number less then iatom because those pairs
! have been counted previously

! Loop over points in the atomic mesh
          do imesh = 1, nam
            density = 0.0d0
            index1 = index0 + iam2rc(imesh)
            index2 = ie2r(index1) - 1
            r1p = ram2rc(:,imesh) + u1x
            r1pin = sqrt(r1p(1)**2 + r1p(2)**2 + r1p(3)**2)

            ! loop over shells of iatom
            do imu = 1, norb_mu
              n1 = species(in1)%orbital(imu)%issh
              r1max = species(in1)%shell(n1)%rcutoffA
              l1 = species(in1)%orbital(imu)%l
              m1 = species(in1)%orbital(imu)%m
              Qneutral = t%atom(iatom)%shell(n1)%Qneutral/(2.0d0*l1 + 1.0d0)
              density = density + Qneutral*(psiofr(r1pin, r1max, in1, n1)*Ylm (r1p, l1, m1))**2
              vnaG(index2) = vnaG(index2) + vnaofr(r1pin, in1, n1)
            end do ! end loop over shells
            rhoG0(index2) = density + rhoG0(index2)
          end do ! end loop over mesh
        end do ! end loop over atom
        vnaG = vnaG*P_eq2

! calculate the total denstity
        density = 0.0d0
        do index2 = 0, nrm - 1
          density = density + rhoG0(index2)*dvolume
        end do
        write (logfile,*)
        write (logfile,*) ' -- Total atomic density before renormalization = ', density

! the renormalization factor
        qtot = 0.0d0
        do iatom = 1, t%natoms
          qtot = qtot + t%atom(iatom)%Q0
        end do
        renorm = qtot/density

! check total density after renormalization
        density = 0.0d0
        do index2 = 0, nrm - 1
          rhoG0(index2) = rhoG0(index2)*renorm
          density = density + rhoG0(index2)*dvolume
        end do
        write (logfile,*) ' -- Total atomic density after renormalization = ', density

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_project_grid


! ===========================================================================
! destroy_grid
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the grid
!! information - these arrays were read in by read_species.
!
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
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_grid ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
! None

! Deallocate Arrays
! ===========================================================================
       deallocate (ie2r)

       deallocate (iam2rc)
       deallocate (ram2rc)
       deallocate (ratom2g)

       deallocate (vnaG, vcaG, vxcG)
       deallocate (rhoG0, drhoG)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_grid


! End Module
! ===========================================================================
        end module M_grid
