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

! Program Description
! ===========================================================================
!       This is a module containing all of the assembler programs required
! to assemble all of the matrix elements for the two-center interactions for
! the Harris interactions.
! It contains the following subroutines within the module:
!
!       M_assemble_3c.f90 - assemble a general 3c matrix
!
! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
! located in the Fdata directory.  This list will change depending on
! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_assemble_3c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_3c
        use M_rotations

! Type Declaration
! ===========================================================================

! module procedures
        contains

! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calculates the three center matrix interactions
! for (interaction, isorp), and add them to corresponding ME2c.
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
! Program Declaration
! ===========================================================================
        subroutine assemble_vna_3c (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_3c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.
        
! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha, iatom, jatom     ! the three parties involved
        integer ibeta, jbeta             ! cells for three atoms
        integer ineigh, mneigh           ! counter over neighbors
        integer in1, in2, indna          ! species numbers
        integer interaction, isorp       ! which interaction and subtype

        integer norb_mu, norb_nu         ! size of the block for the pair

        real z                           ! distance between r1 and r2
        real x, cost                     ! dnabc and angle

        real, dimension (3, 3) :: eps    ! the epsilon matrix
        real, dimension (3) :: r1, r2, r3, r12  ! positions
        real, dimension (3) :: sighat    ! unit vector along r2 - r1
        real, dimension (3) :: rhat      ! unit vector along bc - r3

        real, dimension (:, :), allocatable :: bcnam
        real, dimension (:, :), allocatable :: bcnax

        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          indna = s%atom(ialpha)%imass
          r3 = s%atom(ialpha)%ratom

          ! loop over the common neigbor pairs of ialpha
          do ineigh = 1, s%neighbors(ialpha)%ncommon
            mneigh = s%neighbors(ialpha)%neigh_common(ineigh)
            if (mneigh .ne. 0) then
              iatom = s%neighbors(ialpha)%iatom_common_j(ineigh)
              ibeta = s%neighbors(ialpha)%iatom_common_b(ineigh)
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              in1 = s%atom(iatom)%imass
              norb_mu = species(in1)%norb_max

              jatom = s%neighbors(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%norb_max

              ! cut some lengthy notation
              pvna=>s%vna(iatom); pvna_neighbors=>pvna%neighbors(mneigh)

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

! ****************************************************************************
! Find rnabc = vector pointing from center of bondcharge to r3
! This gives us the distance dnabc (or x value in the 2D grid).
              r12 = 0.5d0*(r1 + r2)
              x = distance (r12, r3)

              ! unit vector in rnabc direction.
              if (x .lt. 1.0d-05) then
                rhat(1) = 0.0d0
                rhat(2) = 0.0d0
                rhat(3) = 0.0d0
              else
                rhat = (r3 - 0.5d0*(r1 + r2))/x
              end if
              cost = dot_product(sighat, rhat)
              call epsilon_function (rhat, sighat, eps)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_bcna
              isorp = 0

! Allocate block arrays
              allocate (bcnam(norb_mu, norb_nu))
              allocate (bcnax(norb_mu, norb_nu))

              call getMEs_Fdata_3c (in1, in2, indna, interaction, isorp, x,   &
     &                              z, norb_mu, norb_nu, cost, bcnam)

              ! Rotate into crystal coordinates
              call rotate (in1, in2, eps, norb_mu, norb_nu, bcnam, bcnax)

              ! Add this piece into the total
              pvna_neighbors%block = pvna_neighbors%block + bcnax*P_eq2
              deallocate (bcnam, bcnax)
            end if
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
        end subroutine assemble_vna_3c

! End Module
! ===========================================================================
        end module M_assemble_3c
