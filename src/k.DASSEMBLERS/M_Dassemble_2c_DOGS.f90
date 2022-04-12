! copyright info:
!
!                             @Copyright 2016
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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

! M_Dassemble_2c
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the matrix elements for the two-center interactions for
!! the Harris interactions.
!! It contains the following subroutines within the module:
!!
!!       Dassemble_S.f90 - assemble the overlap matrix derivatives
!!       Dassemble_T.f90 - assemble the kinetic matrix derivatives
!!       Dassemble_vna.f90 - assemble neutral atom potential matrix derivatives
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_Dassemble_2c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_2c
        use M_neighbors
        use M_rotations
        use M_Drotations
        use M_assemble_2c
        use M_density_matrix

! Type Declaration
! ===========================================================================
! two-center interactions arrays
       type(T_assemble_neighbors), pointer :: dipole_z (:)

! module procedures
        contains

! ===========================================================================
! Dassemble_S.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the overlap matrix interactions.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
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
        subroutine Dassemble_S (s)
        implicit none

        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh             !< counter over atoms and neighbors
        integer in1, in2, in3             !< species numbers
        integer imu, inu				  !< counter over orbitals
        integer jatom                     !< neighbor of iatom
        integer interaction, isorp        !< which interaction and subtype
        integer num_neigh                 !< number of neighbors
        integer mbeta                     !< the cell containing neighbor of iatom

        integer norb_mu, norb_nu          !< size of the block for the pair

        real z                            !< distance between r1 and r2

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2     !< positions of iatom and jatom
        real, dimension (3) :: sighat     !< unit vector along r2 - r1

! sm = overlap matrix in molecular coordinates
! sx = overlap matrix in crystal coordinates
! dsm = derivative of overlap matrix in molecular coordinates
! vdsm = vectorized derivative of overlap matrix in molecular coordinates
! vdsx = vectorized derivative of overlap matrix in crystal coordinates
        real, dimension (:, :), allocatable :: sm
        real, dimension (:, :), allocatable :: sx
        real, dimension (:, :), allocatable :: dsm
        real, dimension (:, :, :), allocatable :: vdsm
        real, dimension (:, :, :), allocatable :: vdsx

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          poverlap=>s%overlap(iatom)

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            pS_neighbors=>poverlap%neighbors(ineigh)

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (pS_neighbors%Dblock(3, norb_mu, norb_nu))
            pS_neighbors%Dblock = 0.0d0

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)
	        call Depsilon_2c (r1, r2, eps, deps)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
! For these interactions, there are no subtypes and isorp = 0
            isorp = 0
            interaction = P_overlap
            in3 = in2

! sm = overlap matrix in molecular coordinates
! sx = overlap matrix in crystal coordinates
! dsm = derivative of overlap matrix in molecular coordinates
! vdsm = vectorized derivative of overlap matrix in molecular coordinates
! vdsx = vectorized derivative of overlap matrix in crystal coordinates
            allocate (sm (norb_mu, norb_nu)); sm = 0.0d0
            allocate (sx (norb_mu, norb_nu)); sx = 0.0d0
            allocate (dsm (norb_mu, norb_nu)); dsm = 0.0d0
            allocate (vdsm (3, norb_mu, norb_nu)); vdsm = 0.0d0
            allocate (vdsx (3, norb_mu, norb_nu)); vdsx = 0.0d0
            call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,          &
     &                             norb_mu, norb_nu, sm, dsm)

! Apply epsilon, the direction of the bondcharge.
! ****************************************************************************
!
! FORCES
! ****************************************************************************
! dsm is the "scalar" derivative of the matrix; vdsm is the "vector" of the matrix
! When we are done, we get: vdsx as the vector derivative of the matrix.

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
            do inu = 1, norb_nu
              do imu = 1, norb_mu
                if (z .gt. 1.0d-3) vdsm(:,imu,inu) = - eta(:)*dsm(imu,inu)
              end do
            end do

! Drotate then puts the vectors in coordinates alone the bond-charge.
            call Drotate (in1, in2, eps, deps, norb_mu, norb_nu, sm, vdsm, vdsx)

! Store the derivitive, rotate vector matrix.
			pS_neighbors%Dblock = vdsx
            deallocate (sm, sx, dsm, vdsm, vdsx)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_S


! ===========================================================================
! assemble_T.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the rate of change of the kinetic matrix
! interactions and then the kinetic forces from that. The dF/dx of the
! interactions is stored in Kinetic_neighbors%Dblock, alongside the %block
! stored by M_assemble_2c, the forces are stored in a Type Forces%T
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
! Program Declaration
! ===========================================================================
        subroutine Dassemble_T (s)
        implicit none

        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh            !< counter over atoms and neighbors
        integer in1, in2, in3            !< species numbers
        integer imu, inu				  !< counter over orbitals
        integer jatom                    !< neighbor of iatom
        integer interaction, isorp       !< which interaction and subtype
        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor

        integer norb_mu, norb_nu         !< size of the block for the pair

        real z                           !< distance between r1 and r2

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2     !< positions of iatom and jatom
        real, dimension (3) :: sighat     !< unit vector along r2 - r1

! tm = kinetic matrix in molecular coordinates
! tx = kinetic matrix in crystal coordinates
! dtm = derivative of kinetic matrix in molecular coordinates
! vdtm = vectorized derivative of kinetic matrix in molecular coordinates
! vdtx = vectorized derivative of kinetic matrix in crystal coordinates
        real, dimension (:, :), allocatable :: tm
        real, dimension (:, :), allocatable :: tx
        real, dimension (:, :), allocatable :: dtm
        real, dimension (:, :, :), allocatable :: vdtm
        real, dimension (:, :, :), allocatable :: vdtx

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          pkinetic=>s%kinetic(iatom)

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            pK_neighbors=>pkinetic%neighbors(ineigh)

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pK_neighbors%Dblock(3, norb_mu, norb_nu))
            pK_neighbors%Dblock = 0.0d0

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)
	        call Depsilon_2c (r1, r2, eps, deps)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in tm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in tx, where x means crytal
! coordinates.
! For these interactions, there are no subtypes and isorp = 0
            isorp = 0
            interaction = P_kinetic
            in3 = in2

! sm = overlap matrix in molecular coordinates
! sx = overlap matrix in crystal coordinates
! dsm = derivative of overlap matrix in molecular coordinates
! vdsm = vectorized derivative of overlap matrix in molecular coordinates
! vdsx = vectorized derivative of overlap matrix in crystal coordinates
            allocate (tm (norb_mu, norb_nu)); tm = 0.0d0
            allocate (tx (norb_mu, norb_nu)); tx = 0.0d0
            allocate (dtm (norb_mu, norb_nu)); dtm = 0.0d0
            allocate (vdtm (3, norb_mu, norb_nu)); vdtm = 0.0d0
            allocate (vdtx (3, norb_mu, norb_nu)); vdtx = 0.0d0
            call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,         &
     &                             norb_mu, norb_nu, tm, dtm)

! Apply epsilon, the direction of the bondcharge.
! ****************************************************************************
!
! FORCES
! ****************************************************************************
! dtm is the "scalar" derivative of the matrix; vdtm is the "vector" of the matrix
! When we are done, we get: vdtx as the vector derivative of the matrix.

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
            do inu = 1, norb_nu
              do imu = 1, norb_mu
                if (z .gt. 1.0d-3) vdtm(:,imu,inu) = - eta(:)*dtm(imu,inu)
              end do
            end do

! Drotate then puts the vectors in coordinates alone the bond-charge.
 			call Drotate (in1, in2, eps, deps, norb_mu, norb_nu, tm, vdtm, vdtx)

! Store the derivitive, rotate vector matrix.
			pK_neighbors%Dblock = vdtx
            deallocate (tm, tx, dtm, vdtm, vdtx)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_T


! ===========================================================================
! Dassemble_dipole_z.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the derivative of the dipole_z matrix
! interactions.
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
! Program Declaration
! ===========================================================================
        subroutine Dassemble_dipole_z (s)
        implicit none

        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer in1, in2, in3           !< species numbers
        integer imu, inu				  !< counter over orbitals
        integer jatom                   !< neighbor of iatom
        integer interaction, isorp      !< which interaction and subtype
        integer num_neigh               !< number of neighbors
        integer mbeta                   !< the cell containing neighbor of iatom

        integer norb_mu, norb_nu        !< size of the block for the pair

        real z                          !< distance between r1 and r2

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

! dipm = kinetic matrix in molecular coordinates
! ddipm = derivative of kinetic matrix in molecular coordinates
! vddipm = vectorized derivative of kinetic matrix in molecular coordinates
! vddipx = vectorized derivative of kinetic matrix in crystal coordinates
        real, dimension (:, :), allocatable :: dipm
        real, dimension (:, :), allocatable :: ddipm
        real, dimension (:, :, :), allocatable :: vddipm
        real, dimension (:, :, :), allocatable :: vddipx

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pdip_neighbors
        type(T_assemble_neighbors), pointer :: pdipole_z

! Allocate Arrays
! ===========================================================================
        allocate (dipole_z(s%natoms))

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          pdipole_z=>s%dipole_z(iatom)

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            pdip_neighbors=>pdipole_z%neighbors(ineigh)

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (pdip_neighbors%Dblock(3, norb_mu, norb_nu))
            pdip_neighbors%Dblock = 0.0d0

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)
            call Depsilon_2c (r1, r2, eps, deps)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
! For these interactions, there are no subtypes and isorp = 0
            isorp = 0
            interaction = P_dipole_z
            in3 = in2

! dipm = overlap matrix in molecular coordinates
! ddipm = derivative of overlap matrix in molecular coordinates
! vddipm = vectorized derivative of overlap matrix in molecular coordinates
! vddipx = vectorized derivative of overlap matrix in crystal coordinates
            allocate (dipm (norb_mu, norb_nu)); dipm = 0.0d0
            allocate (ddipm (norb_mu, norb_nu)); ddipm = 0.0d0
            allocate (vddipm (3, norb_mu, norb_nu)); vddipm = 0.0d0
            allocate (vddipx (3, norb_mu, norb_nu)); vddipx = 0.0d0
            call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,         &
     &                             norb_mu, norb_nu, dipm, ddipm)

! Apply epsilon, the direction of the bondcharge.
! ****************************************************************************
!
! FORCES
! ****************************************************************************
! dtm is the "scalar" derivative of the matrix; vdtm is the "vector" of the matrix
! When we are done, we get: vdtx as the vector derivative of the matrix.

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
            do inu = 1, norb_nu
              do imu = 1, norb_mu
                if (z .gt. 1.0d-3) vddipm(:,imu,inu) = - eta(:)*ddipm(imu,inu)
              end do
            end do

! Drotate then puts the vectors in coordinates alone the bond-charge.
 			call Drotate (in1, in2, eps, deps, norb_mu, norb_nu, dipm,       &
     &                    vddipm, vddipx)

! Store the derivitive, rotate vector matrix.
			pdip_neighbors%Dblock = vddipx
            deallocate (dipm, ddipm, vddipm, vddipx)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_dipole_z


! ===========================================================================
! Dassemble_vna.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates derivatives of the matrix elements for the
! Hartree interactions. This is the self-consistent version of the code.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
!> @author Barry Haycock
! Box 6315, 209 Hodges Hall
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
        subroutine Dassemble_vna_2c (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer in1, in2, in3           !< species numbers
        integer imu, inu				 !< counter over orbitals
        integer jatom                   !< neighbor of iatom
        integer interaction, isorp      !< which interaction and subtype
        integer num_neigh               !< number of neighbors
        integer matom                   !< matom is the self-interaction atom
        integer mbeta                   !< the cell containing neighbor of iatom
        integer issh                    !< counter over shells

        integer norb_mu, norb_nu        !< size of the block for the pair

        real dQ                         !< net charge on atom
        real rcutoff1_min, rcutoff2_min, rend  !< for smoothing
        real smooth                     !< smoothing value
        real Dsmooth                    !< derivative smoothing value
        real xsmooth                    !< for smoothing function
        real z                          !< distance between r1 and r2

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

! bcnam = Hartree matrix in molecular coordinates
! bcnax = Hartree matrix in crystal coordinates
! dbcnam = derivative of Hartree matrix in molecular coordinates
! vdbcnam = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcnax = vectorized derivative of Hartree matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcnam
        real, dimension (:, :), allocatable :: bcnax
        real, dimension (:, :), allocatable :: dbcnam
        real, dimension (:, :, :), allocatable :: vdbcnam
        real, dimension (:, :, :), allocatable :: vdbcnax

! emnpl = monopole term
! vdemnpl = vectorized derivative of monopole term
        real, dimension (:, :), allocatable :: emnpl
        real, dimension (:, :, :), allocatable :: vdemnpl

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        interface
           real function smoother (r, rend, x)
             real, intent (in) :: r
             real, intent (in) :: rend
             real, intent (in) :: x
           end function smoother
        end interface

        interface
           real function Dsmoother (r, rend, x)
             real, intent (in) :: r
             real, intent (in) :: rend
             real, intent (in) :: x
           end function Dsmoother
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors

        type(T_forces), pointer :: pfi

! Allocate Arrays
! ===========================================================================
        do iatom = 1, s%natoms
          pfi=>s%forces(iatom)
          num_neigh = s%neighbors(iatom)%neighn
!         allocate (pfi%vna_atom (3, num_neigh)); pfi%vna_atom = 0.0d0
          allocate (pfi%vna_ontop (3, num_neigh)); pfi%vna_ontop = 0.0d0
        end do

! Procedure
! ===========================================================================
! Here we assemble only the atom cases first, then assemble the ontop cases.
! This is so that we can get the correct allocation size for the different
! blocks.  We calculate the atom cases in a separate loop.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          pvna=>s%vna(iatom)
          pdenmat=>s%denmat(iatom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            pvna_neighbors=>pvna%neighbors(ineigh)
            pRho_neighbors=>pdenmat%neighbors(ineigh)

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pvna_neighbors%Dblock (3, norb_mu, norb_mu))
            pvna_neighbors%Dblock = 0.0d0

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)
            call Depsilon_2c (r1, r2, eps, deps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! CALL GetDMES AND GET VNA FOR ONTOP CASE
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.

            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcnam). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in vbcnax, the vectorized matrix
! elements; x means crytal coordinates.

! FORCES - ONTOP LEFT CASE
! ****************************************************************************
! For the vna_ontopL case, the potential is in the first atom - left (iatom):
! dbcnam is the "scalar" derivative of the matrix; vdbcnam is the "vector"
! derivative of the matrix in molecular coordinates.  When we are done, we get:
! vdtx as the vector derivative of the matrix in crystal coordinates.
              interaction = P_vna_ontopL
              in3 = in2

! bcnam = Hartree matrix in molecular coordinates
! bcnax = Hartree matrix in crystal coordinates
! dbcnam = derivative of Hartree matrix in molecular coordinates
! vdbcnam = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcnax = vectorized derivative of Hartree matrix in crystal coordinates
              allocate (bcnam (norb_mu, norb_nu)); bcnam = 0.0d0
              allocate (dbcnam (norb_mu, norb_nu)); dbcnam = 0.0d0
              allocate (vdbcnam (3, norb_mu, norb_nu)); vdbcnam = 0.0d0
              allocate (vdbcnax (3, norb_mu, norb_nu)); vdbcnax = 0.0d0

! Neutral atom piece
              isorp = 0
              call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,         &
     &                               norb_mu, norb_nu, bcnam, dbcnam)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  if (z .gt. 1.0d-3) vdbcnam(:,imu,inu) = - eta(:)*dbcnam(imu,inu)
                end do
              end do

 			  call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcnam,    &
     &	                    vdbcnam, vdbcnax)

! Notice the explicit negative sign, this makes it force like.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfi%vna_ontop(:,ineigh) = pfi%vna_ontop(:,ineigh)          &
     &             - pRho_neighbors%block(imu,inu)*vdbcnax(:,imu,inu)*P_eq2
                end do
              end do

! Charged atom case
              do isorp = 1, species(in1)%nssh
                dQ = s%atom(iatom)%shell(isorp)%dQ

! Reinitialize
                bcnam = 0.0d0; dbcnam = 0.0d0
                vdbcnam = 0.0d0; vdbcnax = 0.0d0
                call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,      &
     &                                 norb_mu, norb_nu, bcnam, dbcnam)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdbcnam(:,imu,inu) = - eta(:)*dbcnam(imu,inu)
                  end do
                end do

 			    call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcnam,  &
     &	                      vdbcnam, vdbcnax)

! Notice the explicit negative sign, this makes it force like.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfi%vna_ontop(:,ineigh) = pfi%vna_ontop(:,ineigh)        &
     &               - pRho_neighbors%block(imu,inu)*dQ*vdbcnax(:,imu,inu)*P_eq2
                  end do
                end do
              end do ! end loop over isorp

! FORCES - ONTOP RIGHT CASE
! ****************************************************************************
! For the vna_ontopR case, the potential is in the first atom - left (iatom):
! dbcnam is the "scalar" derivative of the matrix; vdbcnam is the "vector"
! derivative of the matrix in molecular coordinates.  When we are done, we get:
! vdtx as the vector derivative of the matrix in crystal coordinates.
              interaction = P_vna_ontopR
              in3 = in2

! Reinitialize
              bcnam = 0.0d0; dbcnam = 0.0d0
              vdbcnam = 0.0d0; vdbcnax = 0.0d0

! Neutral atom piece
              isorp = 0
              call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,         &
     &                               norb_mu, norb_nu, bcnam, dbcnam)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  if (z .gt. 1.0d-3) vdbcnam(:,imu,inu) = - eta(:)*dbcnam(imu,inu)
                end do
              end do

              call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcnam,    &
     &                      vdbcnam, vdbcnax)

! Notice the explicit negative sign, this makes it force like.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfi%vna_ontop(:,ineigh) = pfi%vna_ontop(:,ineigh)          &
     &             - pRho_neighbors%block(imu,inu)*vdbcnax(:,imu,inu)*P_eq2
                end do
              end do

! Charged atom case
              do isorp = 1, species(in3)%nssh
                dQ = s%atom(jatom)%shell(isorp)%dQ

! Reinitialize
                bcnam = 0.0d0; dbcnam = 0.0d0
                vdbcnam = 0.0d0; vdbcnax = 0.0d0
                call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,      &
     &                                 norb_mu, norb_nu, bcnam, dbcnam)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdbcnam(:,imu,inu) = - eta(:)*dbcnam(imu,inu)
                  end do
                end do

                call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcnam,  &
     &                        vdbcnam, vdbcnax)

! Notice the explicit negative sign, this makes it force like.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfi%vna_ontop(:,ineigh) = pfi%vna_ontop(:,ineigh)        &
     &               - pRho_neighbors%block(imu,inu)*dQ*vdbcnax(:,imu,inu)*P_eq2
                  end do
                end do
              end do ! end loop over isorp
              deallocate (bcnam, dbcnam, vdbcnam, vdbcnax)
            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms

! FORCES - ATM CASE
! ****************************************************************************
! For the vna_atom case, the potential is in the first atom - left (iatom):
! dbcnam is the "scalar" derivative of the matrix; vdbcnam is the "vector"
! derivative of the matrix in molecular coordinates.  When we are done, we get:
! vdtx as the vector derivative of the matrix in crystal coordinates.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          matom = s%neigh_self(iatom)

          ! cut some lengthy notation
          pvna=>s%vna(iatom)
          poverlap=>s%overlap(iatom); pS_neighbors=>poverlap%neighbors(matom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            pvna_neighbors=>pvna%neighbors(ineigh)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)
	        call Depsilon_2c (r1, r2, eps, deps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! Find the smoothing quantity - here we calculate the long-range effective
! monopole.  This term is included so that we obtain no discontinuities when
! atoms leave or enter the rcutoff_1 + rcutoff_2 range criteria.
! Therefore, "close" two-center interactions are exact, while more distant
! two-center integrals go to effective monopoles.  The monopoles are effective
! in the sense that the two atoms in the matrix element, each has a different
! charge.  Since they are separated, this gives a monopole contribution at long
! range.

! The smoothing function is found by calling smoother(r,rbegin,rend).
! We define our final matrix element answer as
! smoother(r)*exact_piece + (1 - smoother(r))*longrange.  The distance r is the
! distance of the third center from the "effective" center of the bondcharge.
! The effective center of the bondcharge is (d + rc1 - rc2)/2 from r1 except in
! weird cases (see below). The distance rbegin is the distance at which we
! include only exact answers and do not smooth. The distance rend is the
! distance past which smooth(r) is zero, so that the result is long-range only.
! We skipped self-interaction terms.
            xsmooth = 0.8d0  ! parameter for smoothing

            rcutoff1_min = 99.0d0
            do issh = 1, species(in1)%nssh
              rcutoff1_min = min(rcutoff1_min, species(in1)%shell(issh)%rcutoffA)
            end do

            rcutoff2_min = 99.0d0
            do issh = 1, species(in2)%nssh
              rcutoff2_min = min(rcutoff2_min, species(in2)%shell(issh)%rcutoffA)
            end do

            rend = rcutoff1_min + rcutoff2_min
            smooth = smoother (z, rend, xsmooth)
            Dsmooth = Dsmoother (z, rend, xsmooth)

! For these interactions, there are no subtypes and isorp = 0
            interaction = P_vna_atom
            in3 = in1

! Allocate block size
            norb_nu = species(in3)%norb_max

! bcnam = Hartree matrix in molecular coordinates
! dbcnam = derivative of Hartree matrix in molecular coordinates
! vdbcnam = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcnax = vectorized derivative of Hartree matrix in crystal coordinates
            allocate (bcnam (norb_mu, norb_nu)); bcnam = 0.0d0
            allocate (bcnax (norb_mu, norb_nu)); bcnax = 0.0d0
            allocate (dbcnam (norb_mu, norb_nu)); dbcnam = 0.0d0
            allocate (vdbcnam (3, norb_mu, norb_nu)); vdbcnam = 0.0d0
            allocate (vdbcnax (3, norb_mu, norb_nu)); vdbcnax = 0.0d0

! Neutral atom piece
            isorp = 0
            call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,          &
     &                             norb_mu, norb_nu, bcnam, dbcnam)

! dbcnam is the "scalar" derivative of the matrix; vdbcnam is the "vector" of
! the matrix; vdbcnax as the vector derivative of the matrix in crystal.

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
            do inu = 1, norb_nu
              do imu = 1, norb_mu
                if (z .gt. 1.0d-3) vdbcnam(:,imu,inu) = - eta(:)*dbcnam(imu,inu)
              end do
            end do

            call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcnam,    &
     &                    vdbcnam, vdbcnax)
			pvna_neighbors%Dblock = pvna_neighbors%Dblock + vdbcnax*P_eq2

! Charged atom case
            do isorp = 1, species(in2)%nssh
              dQ = s%atom(jatom)%shell(isorp)%dQ

! Reinitialize
              bcnam = 0.0d0; bcnax = 0.0d0; dbcnam = 0.0d0
              vdbcnam = 0.0d0; vdbcnax = 0.0d0
              call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,        &
     &                               norb_mu, norb_nu, bcnam, dbcnam)

! emnpl = monopole term
! vdemnpl = vectorized derivative of monopole term
              allocate (emnpl (norb_mu, norb_nu)); emnpl = 0.0d0
              allocate (vdemnpl (3, norb_mu, norb_nu)); vdemnpl = 0.0d0

! Note the minus sign. d/dr1 = - eta * d/dd.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  if (z .gt. 1.0d-4) then
                    vdbcnam(:,imu,inu) = - eta(:)*dbcnam(imu,inu)
                    emnpl(imu,inu) = pS_neighbors%block(imu,inu)/z
                    vdemnpl(:,imu,inu) = pS_neighbors%Dblock(:,imu,inu)/z   &
     &                                  + eta(:)*pS_neighbors%block(imu,inu)/z**2
                  end if
                end do
              end do

              call rotate (in1, in3, eps, norb_mu, norb_nu, bcnam, bcnax)
              call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcnam,    &
     &                      vdbcnam, vdbcnax)

              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pvna_neighbors%Dblock(:,imu,inu) = pvna_neighbors%Dblock(:,imu,inu) &
     &             + P_eq2*(- eta(:)*Dsmooth*bcnax(imu,inu) + smooth*vdbcnax(:,imu,inu) &
     &             + eta(:)*Dsmooth*emnpl(imu,inu) + (1 - smooth)*vdemnpl(:,imu,inu))*dQ
                end do
              end do
              deallocate (emnpl, vdemnpl)
            end do ! end loop over isporp
            deallocate (bcnam, bcnax, dbcnam, vdbcnam, vdbcnax)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_vna_2c


! ===========================================================================
! destroy_Dassemble_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the Dassemble_2c
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
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_Dassemble_2c (s)
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

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          do ineigh=1, s%neighbors(iatom)%neighn
            deallocate (s%overlap(iatom)%neighbors(ineigh)%Dblock)
            deallocate (s%kinetic(iatom)%neighbors(ineigh)%Dblock)
            deallocate (s%vna(iatom)%neighbors(ineigh)%Dblock)
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
        end subroutine destroy_Dassemble_2c


! End Module
! ===========================================================================
        end module M_Dassemble_2c
