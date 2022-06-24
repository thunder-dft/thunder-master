! copyright info:
!
!                             @Copyright 2009
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

! M_assemble_2c
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the matrix elements for the two-center interactions for
!! the Harris interactions.
!! It contains the following subroutines within the module:
!!
!!       assemble_S.f90 - assemble the overlap matrix
!!       assemble_T.f90 - assemble the kinetic matrix
!!       assemble_vna_DOGS.f90 - assemble charged atom potential matrix
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_assemble_2c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_2c
        use M_neighbors
        use M_rotations

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains


! ===========================================================================
! assemble_S.f90
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
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_S (s)
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
        integer jatom                   !< neighbor of iatom
        integer interaction, isorp      !< which interaction and subtype
        integer num_neigh               !< number of neighbors
        integer mbeta                   !< the cell containing neighbor of iatom

        integer norb_mu, norb_nu        !< size of the block for the pair

        real z                          !< distance between r1 and r2

        real, dimension (3, 3) :: eps   !< the epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

        real, dimension (:, :), allocatable :: sm
        real, dimension (:, :), allocatable :: sx

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap

! Allocate Arrays
! ===========================================================================
        allocate (s%overlap(s%natoms))

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          poverlap=>s%overlap(iatom)

          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (poverlap%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pS_neighbors=>poverlap%neighbors(ineigh)

            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (pS_neighbors%block(norb_mu, norb_nu))
            pS_neighbors%block = 0.0d0

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

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
! For these interactions, there are no subtypes and isorp = 0
            isorp = 0
            interaction = P_overlap
            in3 = in2

            allocate (sm (norb_mu, norb_nu))
            allocate (sx (norb_mu, norb_nu))
            call getMEs_Fdata_2c (in1, in2, interaction, isorp, z,           &
     &                            norb_mu, norb_nu, sm)
            call rotate (in1, in3, eps, norb_mu, norb_nu, sm, sx)
            pS_neighbors%block = sx
            deallocate (sm)
            deallocate (sx)
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
        end subroutine assemble_S


! ===========================================================================
! assemble_T.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the kinetic matrix interactions.
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
        subroutine assemble_T (s)
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
        integer jatom                    !< neighbor of iatom
        integer interaction, isorp       !< which interaction and subtype
        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor

        integer norb_mu, norb_nu         !< size of the block for the pair

        real z                           !< distance between r1 and r2

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sighat    !< unit vector along r2 - r1

        real, dimension (:, :), allocatable :: tm
        real, dimension (:, :), allocatable :: tx

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic

! Allocate Arrays
! ===========================================================================
        allocate (s%kinetic(s%natoms))

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pkinetic=>s%kinetic(iatom)

          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pkinetic%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pK_neighbors=>pkinetic%neighbors(ineigh)

            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pK_neighbors%block(norb_mu, norb_nu))
            pK_neighbors%block = 0.0d0

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

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in tm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in tx, where x means crytal
! coordinates.
! For these interactions, there are no subtypes and isorp = 0
            isorp = 0
            interaction = P_kinetic
            in3 = in2
            allocate (tm (norb_mu, norb_nu))
            allocate (tx (norb_mu, norb_nu))
            call getMEs_Fdata_2c (in1, in2, interaction, isorp, z,           &
     &                            norb_mu, norb_nu, tm)
            call rotate (in1, in3, eps, norb_mu, norb_nu, tm, tx)
            pK_neighbors%block = tx
            deallocate (tm)
            deallocate (tx)
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
        end subroutine assemble_T


! ===========================================================================
! assemble_vna_2c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates matrix elements for the Hartree interactions.
! This is the self-consistent version of the code.
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
        subroutine assemble_vna_2c (s)
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
        integer jatom                   !< neighbor of iatom
        integer interaction, isorp      !< which interaction and subtype
        integer issh                    !< counting over shells
        integer num_neigh               !< number of neighbors
        integer matom                   !< matom is the self-interaction atom
        integer mbeta                   !< the cell containing neighbor of iatom

        integer norb_mu, norb_nu        !< size of the block for the pair

        real rcutoff1_min, rcutoff2_min, rend  !< for smoothing
        real stinky                     !< smoothing value
        real xsmooth                    !< for smoothing function
        real z                          !< distance between r1 and r2

        real, dimension (3, 3) :: eps   !< the epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

        real, dimension (:, :), allocatable :: bcnam, bcnax, emnpl

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        interface
          function smoother (z, rend, x)
          	real smoother
            real, intent(in) :: z, rend, x
          end function smoother
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap

        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna

! Allocate Arrays
! ===========================================================================
        allocate (s%vna (s%natoms))

! Procedure
! ===========================================================================
! Here we assemble only the ontop cases first, then assemble the atom cases.
! This is so that we can get the correct allocation size for the different
! blocks.  We calculate the atom cases in a separate loop.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pvna=>s%vna(iatom)

          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pvna%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pvna_neighbors=>pvna%neighbors(ineigh)

            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pvna_neighbors%block(norb_mu, norb_nu))
            pvna_neighbors%block = 0.0d0

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

! CALL DOSCENTROS AND GET VNA FOR ONTOP CASE
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.
            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcnam). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in bccax, where x means crytal
! coordinates.
! For the vna_ontopL case, the potential is in the first atom - left (iatom):
              interaction = P_vna_ontopL
              in3 = in2

! Allocate block size
              allocate (bcnam (norb_mu, norb_nu))
              allocate (bcnax (norb_mu, norb_nu))

! Neutral atom case
              isorp = 0
              call getMEs_Fdata_2c (in1, in3, interaction, isorp, z, norb_mu,&
     &                              norb_nu, bcnam)
              call rotate (in1, in3, eps, norb_mu, norb_nu, bcnam, bcnax)
              pvna_neighbors%block = pvna_neighbors%block + bcnax*P_eq2

! For the vna_ontopR case, the potential is in the second atom - right (iatom):
              interaction = P_vna_ontopR
              in3 = in2

! Neutral atom case
              isorp = 0
              call getMEs_Fdata_2c (in1, in2, interaction, isorp, z, norb_mu,&
     &                              norb_nu, bcnam)
              call rotate (in1, in3, eps, norb_mu, norb_nu, bcnam, bcnax)
              pvna_neighbors%block = pvna_neighbors%block + bcnax*P_eq2

              deallocate (bcnam)
              deallocate (bcnax)
            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms


! CALL DOSCENTROS AND GET VNA FOR ATOM CASE
! ****************************************************************************
! The vna two-center terms are: ontop (L), ontop (R), and atom.
! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pvna=>s%vna(iatom)
          poverlap=>s%overlap(iatom)

          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pvna_neighbors=>pvna%neighbors(matom)
            pS_neighbors=>poverlap%neighbors(matom)

            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

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
            stinky = smoother (z, rend, xsmooth)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcnam). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
            interaction = P_vna_atom
            in3 = in1

! Allocate block size
            norb_nu = species(in3)%norb_max
            allocate (bcnam (norb_mu, norb_nu))
            allocate (bcnax (norb_mu, norb_nu))
            allocate (emnpl (norb_mu, norb_nu))

! Set value for emnpl
            emnpl = 0.0d0
            if (z .gt. 1.0d-4) emnpl = pS_neighbors%block/z

! Neutral atom case
            isorp = 0
            call getMEs_Fdata_2c (in1, in2, interaction, isorp, z, norb_mu,  &
     &                            norb_nu, bcnam)
            call rotate (in1, in3, eps, norb_mu, norb_nu, bcnam, bcnax)
            pvna_neighbors%block = pvna_neighbors%block + bcnax*P_eq2

            deallocate (bcnam)
            deallocate (bcnax)
            deallocate (emnpl)
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
        end subroutine assemble_vna_2c


! ===========================================================================
! destroy_assemble_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the assemble_2c_DOGS
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
        subroutine destroy_assemble_2c (s)
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
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%overlap(iatom)%neighbors(ineigh)%block)
            deallocate (s%kinetic(iatom)%neighbors(ineigh)%block)
            deallocate (s%vna(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%overlap(iatom)%neighbors)
          deallocate (s%kinetic(iatom)%neighbors)
          deallocate (s%vna(iatom)%neighbors)
        end do
        deallocate (s%overlap)
        deallocate (s%kinetic)
        deallocate (s%vna)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_assemble_2c

! End Module
! ===========================================================================
        end module M_assemble_2c
