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

! M_assemble_vxc
! Module Description
! ===========================================================================
!>       This is a module containing all of the  programs required
!! to assemble all of the matrix elements for the Horsfield
!! exchange-correlation. See Horsfield AP, 1997,
!! Efficient ab initio tight binding, Phys. Rev. B, Vol: 56, P
!! Pages: 6594-6602, ISSN: 1098-0121
!!
!! It contains the following subroutines within the module:
!!
!!       assemble_vxc_2c : for calculating two-center interactions
!!       assemble_vxc_3c : for calculating three-center interactions
!!       destroy_assemble_vxc : for deallocating all interactions
!!
! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
! located in the Fdata directory.  This list will change depending on
! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_assemble_vxc

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones
        use M_rotations

! /FDATA
        use M_Fdata_1c
        use M_Fdata_2c
        use M_Fdata_3c

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! assemble_vxc.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>  This is the main module for assembling the Vxc matrix element
!! interactions (McWEDA).  Subroutines from M_assemble_rho_McWEDA_rho
!! are used. The results are stored in vxc (potential).
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
        subroutine assemble_vxc (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer logfile                    !< writing to which unit

! Allocate Arrays
! ===========================================================================
        allocate (s%vxc(s%natoms))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        write (logfile,*) ' Calculating Horsfield interactions here. '

! Calculate rho_in (density) matrix elements
        write (logfile,*) ' Calling two-center assemblers. '
        call assemble_vxc_2c (s)

        write (logfile,*) ' Calling three-center assemblers. '
        call assemble_vxc_3c (s)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_vxc


! ===========================================================================
! assemble_vxc_2c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates two-center vxc matrix interactions.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_vxc_2c (s)
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
        integer iatom, ineigh, matom    !< counter over atoms and neighbors
        integer in1, in2, in3           !< species numbers
        integer jatom                   !< neighbor of iatom
        integer interaction             !< which interaction and subtype
        integer num_neigh               !< number of neighbors
        integer mbeta                   !< the cell containing neighbor of iatom

        integer imu, inu
        integer issh, jssh
        integer norb_mu, norb_nu         !< size of the block for the pair

        ! different derivative cases
        integer ideriv

        real z                          !< distance between r1 and r2

        real, dimension (3, 3) :: eps   !< the epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

! bcxcm = Hartree matrix in molecular coordinates
! bcxcx = Hartree matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc

! Allocate Arrays
! ===========================================================================
        allocate (s%vxc (s%natoms))

! Procedure
! ===========================================================================
! Here we assemble only the ontop cases first, then assemble the atom cases.
! This is so that we can get the correct allocation size for the different
! blocks.  We calculate the atom cases in a separate loop.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          pvxc=>s%vxc(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pvxc%neighbors(num_neigh))
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            pvxc_neighbors=>pvxc%neighbors(ineigh)

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pvxc_neighbors%block(norb_mu, norb_nu))
            pvxc_neighbors%block = 0.0d0

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

! CALL DOSCENTROS AND GET VXC FOR ONTOP CASE
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.
            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in bcxcx, where x means crytal
! coordinates.

! For the vxc_ontop case, the potential is in the first atom - left (iatom):
              ideriv = 0
              interaction = P_vxc_ontop
              in3 = in2

              allocate (bcxcm (norb_mu, norb_nu)); bcxcm = 0.0d0
              allocate (bcxcx (norb_mu, norb_nu)); bcxcx = 0.0d0
              call getMEs_Fdata_2c (in1, in2, interaction, ideriv, z,         &
     &                              norb_mu, norb_nu, bcxcm)
              call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)

              pvxc_neighbors%block = pvxc_neighbors%block + bcxcx
              deallocate (bcxcm, bcxcx)
            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms


! CALL DOSCENTROS AND GET VXC FOR ATOM CASE
! ****************************************************************************
! The vxc two-center terms are: ontop (L), ontop (R), and atom.
! First, do vxc_atom case. Here we compute <i | v(j) | i> matrix elements.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some more lengthy notation
          pvxc=>s%vxc(iatom); pvxc_neighbors=>pvxc%neighbors(matom)

! Allocate block size
          allocate (pvxc_neighbors%block(norb_mu, norb_mu))
          pvxc_neighbors%block = 0.0d0

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
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

! CALL DOSCENTROS AND GET VXC FOR ONTOP CASE
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! We calculate here the "on-site" one-center term: mu and nu are in the same
! atom "iatom": vxc (mu,nu) = < mu | V_xc (rho_i) | nu >

! Transform from l-block to mu-block and get the matrix element from the data
! files - from Fdata_1c (onecenter:uncentro)
              do inu = 1, norb_mu
                jssh = species(in1)%orbital(inu)%issh
                do imu = 1, norb_mu
                  issh = species(in1)%orbital(imu)%issh
                  if (species(in1)%orbital(imu)%l                             &
     &                .eq. species(in1)%orbital(inu)%l .and.                  &
     &                species(in1)%orbital(imu)%m                             &
     &                .eq. species(in1)%orbital(inu)%m) then
                    pvxc_neighbors%block(imu,inu) =                           &
     &                pvxc_neighbors%block(imu,inu) + vxc_1c(in1)%V(issh,jssh)
                  end if
                end do
              end do
            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in bcxcx, where x means crytal
! coordinates.

! For these interactions, there are no subtypes and ideriv = 0
              interaction = P_vxc_atom
              in3 = in1
              ideriv = 0

! Allocate block size
              norb_nu = species(in3)%norb_max
              allocate (bcxcm (norb_mu, norb_nu)); bcxcm = 0.0d0
              allocate (bcxcx (norb_mu, norb_nu)); bcxcx = 0.0d0
              call getMEs_Fdata_2c (in1, in2, interaction, ideriv, z,         &
     &                              norb_mu, norb_nu, bcxcm)
              call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)

              pvxc_neighbors%block = pvxc_neighbors%block + bcxcx
              deallocate (bcxcm, bcxcx)
            end if ! end if for r1 .eq. r2 case
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
        end subroutine assemble_vxc_2c


! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calculates the three center matrix interactions
! for (interaction, ideriv), and add them to corresponding ME2c.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_vxc_3c (s)
        implicit none

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
        integer in1, in2, indna          !< species numbers
        integer interaction              !< which interaction and subtype

        integer norb_mu, norb_nu         !< size of the block for the pair

        ! different derivative cases
        integer ideriv

        real z                           !< distance between r1 and r2
        real x, cost                     !< dnabc and angle

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2, rna  !< positions - iatom, jatom, ialpha
        real, dimension (3) :: r21, rnabc   !< vectors
        real, dimension (3) :: sighat    !< unit vector along r2 - r1
        real, dimension (3) :: rhat      !< unit vector along bc - r3

        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx

        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc

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
          indna = s%atom(ialpha)%imass
          rna = s%atom(ialpha)%ratom

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
              pvxc=>s%vxc(iatom); pvxc_neighbors=>pvxc%neighbors(mneigh)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
              r21 = r2 - r1
              z = distance (r1, r2)

              ! unit vector in sigma direction.
              if (z .lt. 1.0d-05) then
                sighat(1) = 0.0d0
                sighat(2) = 0.0d0
                sighat(3) = 1.0d0
              else
                sighat = r21/z
              end if

! ****************************************************************************
! Find rnabc = vector pointing from center of bondcharge to rna
! This gives us the distance dnabc (or x value in the 2D grid).
              rnabc = rna - (r1 + r21/2.0d0)
              x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)

              ! unit vector in rnabc direction.
              if (x .lt. 1.0d-05) then
                rhat(1) = 0.0d0
                rhat(2) = 0.0d0
                rhat(3) = 0.0d0
              else
                rhat = rnabc/x
              end if
              cost = dot_product(sighat, rhat)
              call epsilon_function (rhat, sighat, eps)

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_xc3c

! Allocate block arrays
              allocate (bcxcm(norb_mu, norb_nu))
              allocate (bcxcx(norb_mu, norb_nu))

! Neutral atom case
              ideriv = 0
              call getMEs_Fdata_3c (in1, in2, indna, interaction, ideriv, x,  &
     &                              z, norb_mu, norb_nu, cost, bcxcm)

              ! Rotate into crystal coordinates
              call rotate (in1, in2, eps, norb_mu, norb_nu, bcxcm, bcxcx)

              ! Add this piece into the total
              pvxc_neighbors%block  = pvxc_neighbors%block + bcxcx
              deallocate (bcxcm, bcxcx)
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
        end subroutine assemble_vxc_3c


! ===========================================================================
! destroy_assemble_vxc_McWEDA
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the vxc
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
        subroutine destroy_assemble_vxc (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                            !< counter over atoms
        integer ineigh                           !< counter over neighbors

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%vxc(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%vxc(iatom)%neighbors)
        end do
        deallocate (s%vxc)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_assemble_vxc

! End Module
! ===========================================================================
        end module M_assemble_vxc
