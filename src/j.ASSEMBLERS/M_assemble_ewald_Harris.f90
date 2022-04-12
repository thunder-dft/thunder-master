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

! M_assemble_ewald
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the interactions for the short range and long range
!! ewald interactions.
!!
!! It contains the following subroutines within the module:
!!
!!       assemble_ewaldsr.f90 - assemble the overlap matrix
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_assemble_ewald
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_2c

! Type Declaration
! ===========================================================================
! two-center interactions arrays
! Put all the neighbor group belonging to the atom
! so in the end we have something like overlap(mu, nu, ineigh, iatom)
        type(T_assemble_neighbors), pointer :: ewaldsr (:)
        type(T_assemble_neighbors), pointer :: ewaldlr (:)

! module procedures
        contains


! ===========================================================================
! assemble_ewaldsr.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is a dummy routine.
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
        subroutine assemble_ewaldsr (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh, matom    !< counter over atoms and neighbors
        integer in1, in2                !< species numbers
        integer jatom                   !< neighbor of iatom
        integer num_neigh               !< number of neighbors

        integer norb_mu, norb_nu        !< size of the block for the pair

        type(T_assemble_block), pointer :: pSR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldsr

! Allocate Arrays
! ===========================================================================
        allocate (s%ewaldsr (s%natoms))

! Procedure
! ===========================================================================
! First loop over the atoms in the central cell and allocate.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pewaldsr=>s%ewaldsr(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pewaldsr%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pSR_neighbors=>pewaldsr%neighbors(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pSR_neighbors%blocko(norb_mu, norb_nu))
            pSR_neighbors%blocko = 0.0d0
          end do ! end loop over neighbors
          matom = s%neigh_self(iatom)

          ! cut some lengthy notation
          pSR_neighbors=>pewaldsr%neighbors(matom)
          allocate (pSR_neighbors%block(norb_mu, norb_mu))
          pSR_neighbors%block = 0.0d0
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
        end subroutine assemble_ewaldsr


! ===========================================================================
! assemble_ewaldlr.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is a dummy routine.
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
        subroutine assemble_ewaldlr (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh, matom    !< counter over atoms and neighbors
        integer in1, in2                !< species numbers
        integer jatom                   !< neighbor of iatom
        integer num_neigh               !< number of neighbors

        integer norb_mu, norb_nu        !< size of the block for the pair

        type(T_assemble_block), pointer :: pLR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldlr

! Allocate Arrays
! ===========================================================================
        allocate (s%ewaldlr (s%natoms))

! Procedure
! ===========================================================================
! First loop over the atoms in the central cell and allocate.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pewaldlr=>s%ewaldlr(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pewaldlr%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pLR_neighbors=>pewaldlr%neighbors(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pLR_neighbors%blocko(norb_mu, norb_nu))
            pLR_neighbors%blocko = 0.0d0
          end do ! end loop over neighbors
          matom = s%neigh_self(iatom)

          ! cut some lengthy notation
          pLR_neighbors=>pewaldlr%neighbors(matom)
          allocate (pLR_neighbors%block(norb_mu, norb_mu))
          pLR_neighbors%block = 0.0d0
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
        end subroutine assemble_ewaldlr


! ===========================================================================
! destroy_assemble_ewald
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is a dummy routine.
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
        subroutine destroy_assemble_ewald (s)
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
            deallocate (s%ewaldsr(iatom)%neighbors(ineigh)%blocko)
            deallocate (s%ewaldlr(iatom)%neighbors(ineigh)%blocko)
          end do
          matom = s%neigh_self(iatom)
          deallocate (s%ewaldsr(iatom)%neighbors(matom)%block)
          deallocate (s%ewaldlr(iatom)%neighbors(matom)%block)
          deallocate (s%ewaldsr(iatom)%neighbors)
          deallocate (s%ewaldlr(iatom)%neighbors)
        end do
        deallocate (s%ewaldsr)
        deallocate (s%ewaldlr)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_assemble_ewald

! End Module
! ===========================================================================
        end module M_assemble_ewald
