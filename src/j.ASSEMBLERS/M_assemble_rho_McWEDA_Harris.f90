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

! M_assemble_rho_McWEDA
! Module Description
! ===========================================================================
!>       This is a module containing  the assembler programs required
!! to calculate the matrix elements for the densities rho_in and rho_at
!! used in the McWEDA-Harris interactions (e.g. Sankey-Niklewski).
!! It contains the following subroutines within the module:
!!
!!       assemble_rho_2c.f90 - assembles  two center part for rho_in and rho_at
!!       assemble_rho_3c.f90 - three center part for rho_in
!!       assemble_rho_average.f90 - calculates the final result for average
!!                                  densities
!!       assemble_rho_weighted_2c.f90 - assembles two center part
!!                                      for Wrho, Wrho_bond
!!       assemble_rho_weighted_3c.f90 - three center part for Wrho
!!       assemble_S_weighted.f90 - assembles overlap_weighted (and averaged)
!!                                 PRB 71, 235101 (2005):
!!                                 denominators in Eqs. (19), (22) and (25)
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
         module M_assemble_rho_McWEDA
         use M_assemble_blocks
         use M_configuraciones
         use M_Fdata_2c
         use M_Fdata_3c
         use M_rotations

! Type Declaration
! ===========================================================================
! None

! module procedures
         contains

! ===========================================================================
! assemble_rho_2c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates neutral atom  2-center input density matrix
!! interactions (for rho_in and rho_bond)
!
! ===========================================================================
! Code written by:
!> @author Daniel G. Trabada
!! @author Jose Ortega (JOM)
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_rho_2c (s)
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
        integer interaction, isubtype   !< which interaction and subtype
        integer num_neigh               !< number of neighbors
        integer matom                   !< matom is the self-interaction atom
        integer mbeta                   !< the cell containing neighbor of iatom

        integer norb_mu, norb_nu        !< size of the block for the pair

        real Qneutral                   ! charge
        real z                          !< distance between r1 and r2

        real, dimension (3, 3) :: eps   !< the epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: prho_in_neighbors
        type(T_assemble_neighbors), pointer :: prho_in
        type(T_assemble_block), pointer :: prho_bond_neighbors
        type(T_assemble_neighbors), pointer :: prho_bond

! Allocate Arrays
! ===========================================================================
        allocate (s%rho_in(s%natoms))
        allocate (s%rho_bond(s%natoms))

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          prho_in=>s%rho_in(iatom)
          prho_bond=>s%rho_bond(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          allocate (prho_in%neighbors(num_neigh))
          allocate (prho_bond%neighbors(num_neigh))
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            prho_in_neighbors=>prho_in%neighbors(ineigh)
            prho_bond_neighbors=>prho_bond%neighbors(ineigh)

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (prho_in_neighbors%block(norb_mu, norb_nu))
            allocate (prho_bond_neighbors%block(norb_mu, norb_nu))
            prho_in_neighbors%block = 0.0d0
            prho_bond_neighbors%block = 0.0d0

! SET-UP STUFF
! ***************************************************************************
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

! CALL DOSCENTROS AND GET rho_in FOR ONTOP CASE (OFF-SITE matrix elements)
! ***************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.

            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in bcxcx, where x means crytal
! coordinates.
! For the rho_in_ontopL case, the potential is in the first atom - left (iatom):
! <mu|(rho_mu + rho_nu)|nu> -> (left) <mu|(rho_mu)|nu>
              interaction = P_rho_ontopL
              in3 = in2

              allocate (bcxcm (norb_mu, norb_nu))
              allocate (bcxcx (norb_mu, norb_nu))
              do isubtype = 1, species(in1)%nssh
                Qneutral = species(in1)%shell(isubtype)%Qneutral
                call getMEs_Fdata_2c (in1, in3, interaction, isubtype, z,    &
     &                              norb_mu, norb_nu, bcxcm)
                call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)
                prho_in_neighbors%block =                                    &
     &            prho_in_neighbors%block + bcxcx*Qneutral
                prho_bond_neighbors%block =                                  &
     &            prho_bond_neighbors%block + bcxcx*Qneutral
              end do

! For the rho_in_ontopR case, the potential is in the second atom
! - right (iatom): <mu|(rho_mu + rho_nu)|nu> -> (right) <mu|(rho_nu)|nu>
              interaction = P_rho_ontopR
              in3 = in2
              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getMEs_Fdata_2c (in1, in2, interaction, isubtype, z,    &
     &                              norb_mu, norb_nu, bcxcm)
                call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)
                prho_in_neighbors%block =                                    &
     &            prho_in_neighbors%block + bcxcx*Qneutral
                prho_bond_neighbors%block =                                  &
     &            prho_bond_neighbors%block + bcxcx*Qneutral
              end do
              deallocate (bcxcm)
              deallocate (bcxcx)
            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms

! CALL DOSCENTROS AND GET rho_in FOR ATOM CASE
! ***************************************************************************
! The rho_in two-center terms are: ontop (L), ontop (R), and atom.
! First, do rho_in_atom case. Here we compute <i | v(j) | i> matrix elements.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          prho_in=>s%rho_in(iatom); prho_in_neighbors=>prho_in%neighbors(matom)
          prho_bond=>s%rho_bond(iatom); prho_bond_neighbors=>prho_bond%neighbors(matom)

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! SET-UP STUFF
! *************************************************************************
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

            if (iatom .eq. jatom .and. mbeta .eq. 0) then
! one center case : calculate both rho_in and rho_bond
              interaction = P_rho_atom
              in3 = in1

! Allocate block size
              norb_nu = species(in3)%norb_max
              allocate (bcxcm (norb_mu, norb_nu))
              allocate (bcxcx (norb_mu, norb_nu))
              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getMEs_Fdata_2c (in1, in2, interaction, isubtype, z,    &
     &                            norb_mu, norb_nu, bcxcm)
                call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)
                prho_in_neighbors%block =                                    &
     &            prho_in_neighbors%block + bcxcx*Qneutral
                prho_bond_neighbors%block =                                  &
     &            prho_bond_neighbors%block + bcxcx*Qneutral
              end do
              deallocate (bcxcm)
              deallocate (bcxcx)
            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_rho_atom
              in3 = in1

! Allocate block size
              norb_nu = species(in3)%norb_max
              allocate (bcxcm (norb_mu, norb_nu))
              allocate (bcxcx (norb_mu, norb_nu))
              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getMEs_Fdata_2c (in1, in2, interaction, isubtype, z,    &
     &                                norb_mu, norb_nu, bcxcm)
                call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)
                prho_in_neighbors%block =                                    &
     &            prho_in_neighbors%block + bcxcx*Qneutral
              end do
              deallocate (bcxcm)
              deallocate (bcxcx)
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
        end subroutine assemble_rho_2c


! ===========================================================================
! assemble_rho_3c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates neutral atom 3-center matrix interactions
!! for rho_in
!
! ===========================================================================
! Code written by:
!> @author Daniel G. Trabada
!! @author Jose Ortega (JOM)
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_rho_3c (s)
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
        integer interaction, isubtype    ! which interaction and subtype

        integer norb_mu, norb_nu         !< size of the block for the pair

        real Qneutral                    !< charge
        real z                           !< distances between r1 and r2
        real x, cost                     !< dnabc and angle

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2, r3, r12   !< positions
        real, dimension (3) :: sighat    !< unit vector along r2 - r1
        real, dimension (3) :: rhat      !< unit vector along bc - r3

        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

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

! SET-UP STUFF
! ***************************************************************************
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

! ***************************************************************************
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
              interaction = P_rho_3c

! Allocate block arrays
              allocate (bcxcm(norb_mu, norb_nu))
              allocate (bcxcx(norb_mu, norb_nu))
              do isubtype = 1, species(in3)%nssh
                Qneutral = species(in3)%shell(isubtype)%Qneutral
                call getMEs_Fdata_3c (in1, in2, in3, interaction, isubtype, x,&
     &                                z, norb_mu, norb_nu, cost, bcxcm)

                ! Rotate into crystal coordinates
                call rotate (in1, in2, eps, norb_mu, norb_nu, bcxcm, bcxcx)

                ! Add this piece into the total
                s%rho_in(iatom)%neighbors(mneigh)%block =                    &
     &            s%rho_in(iatom)%neighbors(mneigh)%block + bcxcx*Qneutral
              end do
              deallocate (bcxcm)
              deallocate (bcxcx)
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
        end subroutine assemble_rho_3c


! ===========================================================================
! assemble_rho_weighted_2c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine assembles  two center part for rho_in_weighted and
!! rho_bond_weighted.
!!       rho_in_weighted: numerator in Eq. (19): PRB 71, 235101 (2005)
!!       rho_bond_weighted : numerator in Eqs. (22), (25):PRB 71, 235101 (2005)
!
! ===========================================================================
! Code written by:
!> @author Daniel G. Trabada
!! @author Jose Ortega (JOM)
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_rho_weighted_2c (s)
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
        integer iatom, ineigh            !< counter over atoms and neighbors
        integer in1, in2, in3            !< species numbers
        integer jatom                    !< neighbor of iatom
        integer interaction, isubtype    !< which interaction and subtype
        integer num_neigh                !< number of neighbors
        integer matom                    !< matom is the self-interaction atom
        integer mbeta                    !< cell containing neighbor of iatom

        integer nssh_i, nssh_j         !< size of the block for the pair

        real z                           !< distance between r1 and r2
        real Qneutral

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom

        real, dimension (:, :), allocatable :: bcxcm

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pWrho_in_neighbors
        type(T_assemble_neighbors), pointer :: prho_in_weighted
        type(T_assemble_block), pointer :: pWrho_bond_neighbors
        type(T_assemble_neighbors), pointer :: prho_bond_weighted

! Allocate Arrays
! ===========================================================================
        allocate (s%rho_in_weighted(s%natoms))
        allocate (s%rho_bond_weighted(s%natoms))

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          nssh_i = species(in1)%nssh

          ! cut some lengthy notation
          prho_in_weighted=>s%rho_in_weighted(iatom)
          prho_bond_weighted=>s%rho_bond_weighted(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          allocate (prho_in_weighted%neighbors(num_neigh))
          allocate (prho_bond_weighted%neighbors(num_neigh))
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            pWrho_in_neighbors=>prho_in_weighted%neighbors(ineigh)
            pWrho_bond_neighbors=>prho_bond_weighted%neighbors(ineigh)

! Allocate block size
            nssh_j = species(in2)%nssh
            allocate (pWrho_in_neighbors%block(nssh_i, nssh_j))
            allocate (pWrho_bond_neighbors%block(nssh_i, nssh_j))
            pWrho_in_neighbors%block = 0.0d0
            pWrho_bond_neighbors%block = 0.0d0

! Calculate the distance between the two centers.
            z = distance (r1, r2)

! CALL DOSCENTROS AND GET rho_in_weighted FOR ONTOP CASE
! (i.e. OFF-SITE matrix elements)
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.

            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm).
! No rotation requiered for this case
! (weights are spherical, see definition Eq. (18) McWEDA paper)
!
! For the rho_in_weighted_ontopL case, the potential is in the first atom -
! left (iatom): <mu|(rho_mu + rho_nu)|nu> -> (left) <mu|(rho_mu)|nu>
              interaction = P_rhoS_ontopL
              in3 = in2
              allocate (bcxcm (nssh_i, nssh_j))

              do isubtype = 1, species(in1)%nssh
                Qneutral = species(in1)%shell(isubtype)%Qneutral
                call getMEs_Fdata_2c (in1, in3, interaction, isubtype, z,    &
     &                                nssh_i, nssh_j, bcxcm)

                pWrho_in_neighbors%block =                                   &
     &            pWrho_in_neighbors%block + bcxcm*Qneutral
                pWrho_bond_neighbors%block =                                 &
     &            pWrho_bond_neighbors%block + bcxcm*Qneutral
              end do

! For the rho_in_weighted_ontopR case, the potential is in the second atom -
! right (iatom): <mu|(rho_mu + rho_nu)|nu> -> (right) <mu|(rho_nu)|nu>
              interaction = P_rhoS_ontopR
              in3 = in2
              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getMEs_Fdata_2c (in1, in2, interaction, isubtype, z,    &
     &                                nssh_i, nssh_j, bcxcm)

                pWrho_in_neighbors%block =                                   &
     &            pWrho_in_neighbors%block + bcxcm*Qneutral
                pWrho_bond_neighbors%block =                                 &
     &            pWrho_bond_neighbors%block + bcxcm*Qneutral
              end do
              deallocate (bcxcm)

            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms

! CALL DOSCENTROS AND GET rho_in_weighted FOR ATOM CASE
! ****************************************************************************
! The rho_in_weighted two-center terms are: ontop (L), ontop (R), and atom.
! First, do rho_in_weighted_atom case.
! Here we compute <i|v(j)|i> matrix elements.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          nssh_i = species(in1)%nssh
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          prho_in_weighted=>s%rho_in_weighted(iatom)
          pWrho_in_neighbors=>prho_in_weighted%neighbors(matom)
          prho_bond_weighted=>s%rho_bond_weighted(iatom)

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation

            pWrho_bond_neighbors=>prho_bond_weighted%neighbors(matom)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Calculate the distance between the two centers.
            z = distance (r1, r2)

            if (iatom .eq. jatom .and. mbeta .eq. 0) then
! one center case
              interaction = P_rhoS_atom
              in3 = in1

! Allocate block size
              allocate (bcxcm (nssh_i, nssh_i))
              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getMEs_Fdata_2c (in1, in2, interaction, isubtype, z,    &
     &                                nssh_i, nssh_i, bcxcm)

                pWrho_in_neighbors%block =                                   &
     &            pWrho_in_neighbors%block + bcxcm*Qneutral
                pWrho_bond_neighbors%block =                                 &
     &            pWrho_bond_neighbors%block + bcxcm*Qneutral
              end do
              deallocate (bcxcm)
            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). No rotations here
              interaction = P_rhoS_atom
              in3 = in1

! Allocate block size
              allocate (bcxcm (nssh_i, nssh_i))
              do isubtype = 1, species(in2)%nssh
                Qneutral = species(in2)%shell(isubtype)%Qneutral
                call getMEs_Fdata_2c (in1, in2, interaction, isubtype, z,    &
     &                                nssh_i, nssh_i, bcxcm)
                pWrho_in_neighbors%block =                                   &
     &            pWrho_in_neighbors%block + bcxcm*Qneutral
              end do
              deallocate (bcxcm)
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None - 

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_rho_weighted_2c


! ===========================================================================
! assemble_rho_weighted_3c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine assembles  three center part for rho_weighted
!!       rho_weighted : numerator in Eq. (19): PRB 71, 235101 (2005)
! ===========================================================================
! Code written by:
!> @author Daniel G. Trabada
!! @author Jose Ortega (JOM)
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_rho_weighted_3c (s)
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
        integer in1, in2, in3
        integer interaction, isubtype    ! which interaction and subtype

        integer norb_mu, norb_nu         !< size of the block for the pair

        real z                           !< distances between r1 and r2
        real x, cost                     !< dnabc and angle
        real Qneutral

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2, r3, r12     !< positions
        real, dimension (3) :: sighat    !< unit vector along r2 - r1
        real, dimension (3) :: rhat      !< unit vector along bc - r3

        real, dimension (:, :), allocatable :: bcxcm

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
              norb_mu = species(in1)%nssh

              jatom = s%neighbors(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%nssh

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
! coordinates (stored in sm). No rotations requiered here
              interaction = P_rhoS_3c

! Allocate block arrays
              allocate (bcxcm(norb_mu, norb_nu))

              do isubtype = 1, species(in3)%nssh
                Qneutral = species(in3)%shell(isubtype)%Qneutral
                call getMEs_Fdata_3c (in1, in2, in3, interaction, isubtype, x,&
     &                                z, norb_mu, norb_nu, cost, bcxcm)

                ! Add this piece into the total
                s%rho_in_weighted(iatom)%neighbors(mneigh)%block =           &
     &             s%rho_in_weighted(iatom)%neighbors(mneigh)%block          &
     &              + bcxcm*Qneutral
              end do
              deallocate (bcxcm)
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
        end subroutine assemble_rho_weighted_3c


! ===========================================================================
! destroy_assemble_rho
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the rho_McWEDA
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
        subroutine destroy_assemble_rho (s)
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
            deallocate (s%rho_in_weighted(iatom)%neighbors(ineigh)%block)
            deallocate (s%rho_bond_weighted(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%rho_in_weighted(iatom)%neighbors)
          deallocate (s%rho_bond_weighted(iatom)%neighbors)
        end do
        deallocate (s%rho_in_weighted)
        deallocate (s%rho_bond_weighted)

        do iatom = 1, s%natoms
          do ineigh=1, s%neighbors(iatom)%neighn
            deallocate (s%rho_in(iatom)%neighbors(ineigh)%block)
            deallocate (s%rho_bond(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%rho_in(iatom)%neighbors)
          deallocate (s%rho_bond(iatom)%neighbors)
        end do
        deallocate (s%rho_in)
        deallocate (s%rho_bond)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_assemble_rho

! End Module
! ===========================================================================
        end module M_assemble_rho_McWEDA
