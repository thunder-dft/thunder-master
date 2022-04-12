! copyright info:
!
!                             @Copyright 2008
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
!>       This is a module containing all of the assembler programs required
!! to assemble all of the matrix elements for the two-center interactions for
!! the Harris interactions.
!! It contains the following subroutines within the module:
!!
!!       M_assemble_vna_3c.f90 - assemble Hartree 3c matrix elements
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
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
!>       This routine calculates the three center matrix interactions
!! for (interaction, isubtype), and add them to corresponding ME2c.
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
        integer ialpha, iatom, jatom     !< the three parties involved
        integer ibeta, jbeta             !< cells for three atoms
        integer ineigh, mneigh           !< counter over neighbors
        integer in1, in2, in3            !< species numbers
        integer interaction, isorp       !< which interaction and subtype
        integer issh                     !< counting over shells
        integer norb_mu, norb_nu         !< size of the block for the pair

        real distance_13, distance_23  !< distance from 3rd atom
        real dQ                          !< net charge on atom
        real rcutoff1_min, rcutoff2_min, rcutoff3_min, rend  !< for smoothing
        real stinky                      !< smoothing value
        real xsmooth                     !< for smoothing function
        real z                           !< distance between r1 and r2
        real x, cost                     !< dnabc and angle

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2, r3, r12!< positions
        real, dimension (3) :: sighat    !< unit vector along r2 - r1
        real, dimension (3) :: rhat      !< unit vector along bc - r3

        real, dimension (:, :), allocatable :: bcnam
        real, dimension (:, :), allocatable :: bcnax
        real, dimension (:, :), allocatable :: dterm
        real, dimension (:, :), allocatable :: sterm
        real, dimension (:, :), allocatable :: emnpl

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
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
        type(T_assemble_block), pointer :: pdip_neighbors
        type(T_assemble_neighbors), pointer :: pdipole_z

        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          in3 = s%atom(ialpha)%imass
          r3 = s%atom(ialpha)%ratom

          ! loop over the common neigbor pairs of ialp
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
              poverlap=>s%overlap(iatom); pS_neighbors=>poverlap%neighbors(mneigh)
              pdipole_z=>s%dipole_z(iatom); pdip_neighbors=>pdipole_z%neighbors(mneigh)

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

! Find other distances -
              distance_13 = distance (r3, r1)
              distance_23 = distance (r3, r2)

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

              rcutoff3_min = 99.0d0
              do issh = 1, species(in3)%nssh
                rcutoff3_min = min(rcutoff3_min, species(in3)%shell(issh)%rcutoffA)
              end do

              rend = rcutoff1_min + rcutoff3_min
              stinky = smoother (distance_13, rend, xsmooth)

              rend = rcutoff2_min + rcutoff3_min
              stinky = stinky*smoother (distance_23, rend, xsmooth)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_bcna

! Allocate block arrays
              allocate (bcnam (norb_mu, norb_nu))
              allocate (bcnax (norb_mu, norb_nu))
              allocate (sterm (norb_mu, norb_nu))
              allocate (dterm (norb_mu, norb_nu))
              allocate (emnpl (norb_mu, norb_nu))

! Set value for emnpl
              sterm = pS_neighbors%block/2.0d0
              dterm = pdip_neighbors%block/z
              emnpl = (sterm - dterm)/distance_13 + (sterm + dterm)/distance_23

! Neutral atom case
              isorp = 0
              call getMEs_Fdata_3c (in1, in2, in3, interaction, isorp, x,    &
     &                              z, norb_mu, norb_nu, cost, bcnam)

              ! Rotate into crystal coordinates
              call rotate (in1, in2, eps, norb_mu, norb_nu, bcnam, bcnax)

              ! Add this piece into the total
              s%vna(iatom)%neighbors(mneigh)%block =                         &
     &          s%vna(iatom)%neighbors(mneigh)%block + bcnax*P_eq2

! Charged atom cases
              do isorp = 1, species(in3)%nssh
                call getMEs_Fdata_3c (in1, in2, in3, interaction, isorp, x,  &
     &                                z, norb_mu, norb_nu, cost, bcnam)

                ! Rotate into crystal coordinates
                call rotate (in1, in2, eps, norb_mu, norb_nu, bcnam, bcnax)

                ! Add this piece into the total
                dQ = s%atom(ialpha)%shell(isorp)%dQ
                pvna_neighbors%block = pvna_neighbors%block                  &
     &            + dQ*(stinky*bcnax + (1.0d0 - stinky)*emnpl)*P_eq2
              end do

              deallocate (bcnam)
              deallocate (bcnax)

              deallocate (sterm)
              deallocate (dterm)
              deallocate (emnpl)
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
