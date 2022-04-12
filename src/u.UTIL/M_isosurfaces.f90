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

! M_isosurfaces
! Module Description
! ===========================================================================
!       This is a module containing all information related to projecting
! the various different isosurfaces onto the grids - either for putting the
! charge densities on the grid for viewing, etc. The following subroutines
! are called here:
!
!       project_orbitals_grid.f90 - driver program which chooses how to
!                                   project the molecular orbitals onto the grid
!       band_plots.f90 - projects chosen band levels to the grid
!       Fukui.f90 - projects the nucleophilic and electrophilic Fukui
!                   functions to the grid
!       Fukui_radical.f90 - projects the radical attacking Fukui function
!                           to the grid
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
! Module Declaration
! ===========================================================================
        module M_isosurfaces
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
! project_orbitals_grid
! ===========================================================================
! Subroutine Description
! ===========================================================================
! Project the molecular orbitals onto the grid.
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
        subroutine project_orbitals_grid (t, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

        integer itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
        if (iwriteout_ewf .eq. -2) then
          call Fukui (t, itime_step)
          call Fukui_radical (t, itime_step)
        end if
        if (iwriteout_ewf .eq. -1) call Fukui_radical (t, itime_step)
        if (iwriteout_ewf .eq. 1) call band_plots (t, itime_step)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine project_orbitals_grid


! ===========================================================================
! band_plots
! ===========================================================================
! Subroutine Description
! ===========================================================================
!
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
        subroutine band_plots (t, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

        integer itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over atoms
        integer imesh                       !< counter over the mesh
        integer i0, j0, k0                  !< grid index points
        integer index0, index1, index2      !< different indexing counters
        integer in1                         !< species numbers
        integer inpfile                     !< reading from which unit
        integer logfile                     !< writing to which unit

        integer n1, l1, m1                  !< quantum numbers

        ! parameters for plotting charge densities
        integer iband, ipband
        integer npbands, npbands_min        !< read in from .ewf.inp file

        integer imu                         !< loop over shells
        integer mmu                         !< location within the coefficient vector
        integer norb_mu                     !< size of the block for the pair

        integer, dimension (3) :: nr

        real density
        real r1pin                          !< distance from centers to integration point
        real r1max                          !< maximum extent of cutoff radius

        real, dimension (3) :: r1           !< positions of iatom and jatom
        real, dimension (3) :: r1p          !< shifted positions of iatom and jatom
        real, dimension (3) :: g1
        real, dimension (3) :: u, u1x
        real, dimension (3) :: x0
        real, dimension (3, 3) :: amatrix
        real, dimension (3, 3) :: ainverse

        real, target, allocatable, dimension (:) :: ewfrho
        real, allocatable, dimension (:,:) :: positive
        real, allocatable, dimension (:,:) :: negative

        real, allocatable, dimension (:) :: psi1

        character (len=25) xsfname
        character (len=25) message
        character (len = 25) :: slogfile

        logical read_ewf

        real, dimension (:), pointer   :: pmat
        interface
          subroutine writeout_xsf (t, xsf, xsfname, message)
            use M_configuraciones
            type(T_structure), target :: t
            character (len = 25), intent (in) :: xsfname
            character (len = 25), intent (in) :: message
            real, pointer :: xsf (:)
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
        allocate (ewfrho (0:nrm - 1))
        allocate (positive (irm1, irm2))
        allocate (negative (irm1, irm2))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile
        inpfile = t%inpfile

! Calculate the absorption probabilities.
        write (logfile,*)
        write (logfile,*) ' Plotting the densities for certain bands. '

! Determine if the abs.inp files exists - if not, then exit.
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.ewf.inp'
        inquire (file = slogfile, exist = read_ewf)
        if (read_ewf) then
! Read from input file - gives dos options
          write (logfile,*) ' Reading from ewf.inp file! '
          open (unit = inpfile, file = slogfile, status = 'old')
          read (inpfile,*) npbands, npbands_min      ! energy grid
          close (inpfile)
        else
          write (logfile,*) ' There is no input file for plotting the charge densities. '
          return
        end if
        write (logfile,*)

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

! Loop over the bands that we are plotting
        do ipband = 0, npbands - 1
          iband = npbands_min + ipband

          ! reset ewfrho
          ewfrho = 0.0d0

! Loop over atoms
          do iatom = 1, t%natoms
            r1 = t%atom(iatom)%ratom
            in1 = t%atom(iatom)%imass
            norb_mu = species(in1)%norb_max

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

! Assemble density
              density = 0.0d0
              do imu = 1, species(in1)%norb_max
                mmu = imu + t%iblock_slot(iatom)
                density = density + t%kpoints(1)%c(mmu,iband)*psi1(imu)
              end do ! do imu
              deallocate (psi1)

! store variation of density at given point
              ewfrho(index2) = ewfrho(index2) + density
            end do ! end loop over mesh
          end do ! end loop over atoms

! Visualization in XCrysDen
! write out charge plot for this band into xsf file
          pmat => ewfrho
          slogfile = t%basisfile(:len(trim(t%basisfile))-4)
          write (xsfname, '(".",i4.4,"-",i4.4,".xsf")') iband, itime_step
          xsfname = trim(slogfile)//xsfname
          message = '3D-plot of band density '
          call writeout_xsf (t, pmat, xsfname, message)

! Visualization in ppm

        end do ! end loop over bands

! Deallocate Arrays
! ===========================================================================
        deallocate (ewfrho)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine band_plots


! ===========================================================================
! Fukui_radical
! ===========================================================================
! Subroutine Description
! ===========================================================================
!
!       This routine calculates the radical attack from the Fukui functions.
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
        subroutine Fukui_radical (t, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

        integer itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over atoms
        integer imesh                       !< counter over the mesh
        integer i0, j0, k0                  !< grid index points
        integer index0, index1, index2      !< different indexing counters
        integer in1                         !< species numbers
        integer inpfile                     !< reading from which unit
        integer logfile                     !< writing to which unit

        integer n1, l1, m1                  !< quantum numbers

        ! parameters for plotting charge densities
        integer ihomo
        integer ipband

        integer imu                         !< loop over shells
        integer mmu                         !< location within the coefficient vector
        integer norb_mu                     !< size of the block for the pair

        integer, dimension (3) :: nr

        real density
        real fsign                          !< which sign (ihomo should be negative)
        real r1pin                          !< distance from centers to integration point
        real r1max                          !< maximum extent of cutoff radius

        real, dimension (3) :: r1           !< positions of iatom and jatom
        real, dimension (3) :: r1p          !< shifted positions of iatom and jatom
        real, dimension (3) :: g1
        real, dimension (3) :: u, u1x
        real, dimension (3) :: x0
        real, dimension (3, 3) :: amatrix
        real, dimension (3, 3) :: ainverse

        real, target, allocatable, dimension (:) :: ewfrho
        real, allocatable, dimension (:,:) :: positive
        real, allocatable, dimension (:,:) :: negative

        real, allocatable, dimension (:) :: psi1

        character (len=25) xsfname
        character (len=25) message
        character (len = 25) :: slogfile

        real, dimension (:), pointer   :: pmat
        interface
          subroutine writeout_xsf (t, xsf, xsfname, message)
            use M_configuraciones
            type(T_structure), target :: t
            character (len = 25), intent (in) :: xsfname
            character (len = 25), intent (in) :: message
            real, pointer :: xsf (:)
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
        allocate (ewfrho (0:nrm-1))
        allocate (positive (irm1, irm2))
        allocate (negative (irm1, irm2))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile
        inpfile = t%inpfile

! Calculate the absorption probabilities.
        write (logfile,*)
        write (logfile,*) ' Plotting the radical attack Fukui function. '

! Loop over the bands that we are plotting
        if (modulo(int(t%ztot), 2) .eq. 0) then
          ihomo = int(t%ztot)/2
        else
          ihomo = int(t%ztot)/2 + 1
        end if

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

! ===========================================================================
! radical attack
! ===========================================================================
        ! initialize the density here because we are summing the HOMO
        ! and LUMO densities
        ewfrho = 0.0d0
        do ipband = ihomo, ihomo + 1
          if (ipband .eq. ihomo) fsign = -1.0d0

! Loop over atoms
          do iatom = 1, t%natoms
            r1 = t%atom(iatom)%ratom
            in1 = t%atom(iatom)%imass
            norb_mu = species(in1)%norb_max

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

! Assemble density
              density = 0.0d0
              do imu = 1, species(in1)%norb_max
                mmu = imu + t%iblock_slot(iatom)
                density = density + t%kpoints(1)%c(mmu,ipband)*psi1(imu)
              end do ! do imu
              deallocate (psi1)

! store variation of density at given point
              ewfrho(index2) = ewfrho(index2) + fsign*density
            end do ! end loop over mesh
          end do ! end loop over atoms
        end do ! end loop over bands

! Visualization in XCrysDen
! write out charge plot for this band into xsf file
        pmat => ewfrho
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        write (xsfname, '(".",i4.4,"-",i4.4,".xsf")') 0, itime_step
        xsfname = trim(slogfile)//xsfname
        message = '3D-plot of Fukui function'
        call writeout_xsf (t, pmat, xsfname, message)

! Deallocate Arrays
! ===========================================================================
        deallocate (ewfrho)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Fukui_radical


! ===========================================================================
! Fukui
! ===========================================================================
! Subroutine Description
! ===========================================================================
!
!       This routine calculates the nucleophilic and electrophilic attack
! from the Fukui functions.
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
        subroutine Fukui (t, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

        integer itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over atoms
        integer imesh                       !< counter over the mesh
        integer i0, j0, k0                  !< grid index points
        integer index0, index1, index2      !< different indexing counters
        integer in1                         !< species numbers
        integer inpfile                     !< reading from which unit
        integer logfile                     !< writing to which unit

        integer n1, l1, m1                  !< quantum numbers

        ! parameters for plotting charge densities
        integer ihomo
        integer ipband

        integer imu                         !< loop over shells
        integer mmu                         !< location within the coefficient vector
        integer norb_mu                     !< size of the block for the pair

        integer, dimension (3) :: nr

        real density
        real r1pin                          !< distance from centers to integration point
        real r1max                          !< maximum extent of cutoff radius

        real, dimension (3) :: r1           !< positions of iatom and jatom
        real, dimension (3) :: r1p          !< shifted positions of iatom and jatom
        real, dimension (3) :: g1
        real, dimension (3) :: u, u1x
        real, dimension (3) :: x0
        real, dimension (3, 3) :: amatrix
        real, dimension (3, 3) :: ainverse

        real, target, allocatable, dimension (:) :: ewfrho
        real, target, allocatable, dimension (:) :: ewfdrho
        real, allocatable, dimension (:,:) :: positive
        real, allocatable, dimension (:,:) :: negative

        real, allocatable, dimension (:) :: psi1

        character (len=25) xsfname
        character (len=25) message
        character (len = 25) :: slogfile

        real, dimension (:), pointer   :: pmat
        interface
          subroutine writeout_xsf (t, xsf, xsfname, message)
            use M_configuraciones
            type(T_structure), target :: t
            character (len = 25), intent (in) :: xsfname
            character (len = 25), intent (in) :: message
            real, pointer :: xsf (:)
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
        allocate (ewfrho (0:nrm-1))
        allocate (ewfdrho (0:nrm - 1))
        allocate (positive (irm1, irm2))
        allocate (negative (irm1, irm2))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile
        inpfile = t%inpfile

! Calculate the absorption probabilities.
        write (logfile,*)
        write (logfile,*) ' Plotting the nucleophilic and electrophilic Fukui functions. '

! Loop over the bands that we are plotting
        if (modulo(int(t%ztot), 2) .eq. 0) then
          ihomo = int(t%ztot)/2
        else
          ihomo = int(t%ztot)/2 + 1
        end if

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

! ===========================================================================
! calculate full density
! ===========================================================================
        ! initialize the density
        ewfrho = 0.0d0
        do ipband = 1, ihomo

! Loop over atoms
          do iatom = 1, t%natoms
            r1 = t%atom(iatom)%ratom
            in1 = t%atom(iatom)%imass
            norb_mu = species(in1)%norb_max

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

! Assemble density
              density = 0.0d0
              do imu = 1, species(in1)%norb_max
                mmu = imu + t%iblock_slot(iatom)
                density = density + t%kpoints(1)%c(mmu,ipband)*psi1(imu)
              end do ! do imu
              deallocate (psi1)

! store variation of density at given point
              ewfrho(index2) = ewfrho(index2) + density
            end do ! end loop over mesh
          end do ! end loop over atoms
        end do ! end loop over bands

! ===========================================================================
! nucleophilic attack
! ===========================================================================
! initialize the density
        ewfdrho = 0.0d0

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
          r1 = t%atom(iatom)%ratom
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

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

! Assemble density for the lumo only!
            density = 0.0d0
            do imu = 1, species(in1)%norb_max
              mmu = imu + t%iblock_slot(iatom)
              density = density + t%kpoints(1)%c(mmu,ihomo + 1)*psi1(imu)
            end do ! do imu
            deallocate (psi1)

! store variation of density at given point
            ewfdrho(index2) = ewfdrho(index2) + density
          end do ! end loop over mesh
        end do ! end loop over atoms

! Visualization in XCrysDen
! write out charge plot for this band into xsf file
        ewfdrho = ewfdrho - ewfrho
        pmat => ewfdrho
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        write (xsfname, '(".",i4.4,"-",i4.4,".xsf")') ihomo + 1, itime_step
        xsfname = trim(slogfile)//xsfname
        message = '3D-plot of Fukui function'
        call writeout_xsf (t, pmat, xsfname, message)

! ===========================================================================
! electrophilic attack
! ===========================================================================
! initialize the density
        ewfdrho = 0.0d0

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
          r1 = t%atom(iatom)%ratom
          in1 = t%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

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

! Assemble density
            density = 0.0d0
            do imu = 1, species(in1)%norb_max
              mmu = imu + t%iblock_slot(iatom)
              density = density + t%kpoints(1)%c(mmu,ihomo)*psi1(imu)
            end do ! do imu
            deallocate (psi1)

! store variation of density at given point
            ewfdrho(index2) = ewfdrho(index2) + density
          end do ! end loop over mesh
        end do ! end loop over atoms

! Visualization in XCrysDen
! write out charge plot for this band into xsf file
        ewfdrho = ewfrho - ewfdrho
        pmat => ewfdrho
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        write (xsfname, '(".",i4.4,"-",i4.4,".xsf")') ihomo, itime_step
        xsfname = trim(slogfile)//xsfname
        message = '3D-plot of Fukui function'
        call writeout_xsf (t, pmat, xsfname, message)

! Deallocate Arrays
! ===========================================================================
        deallocate (ewfrho, ewfdrho)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Fukui


! End Module
! ===========================================================================
        end module M_isosurfaces
