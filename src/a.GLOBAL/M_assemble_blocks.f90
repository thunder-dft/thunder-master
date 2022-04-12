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

! M_assemble_blocks
! Module Description
! ===========================================================================
!>       This module simply declares the types which will store the matrix
!! elements - all matrices in these assemblers follow the same format.
! ===========================================================================
        module M_assemble_blocks

! Type Declaration
! ===========================================================================
! interactions arrays
! array for a mu by nu block
        type T_assemble_block
          real, pointer :: block (:, :)      ! regular matrix element block
          real, pointer :: blocko (:, :)     ! ontop matrix element block

          real, pointer :: Dblock (:, :, :)  ! atom piece/three-center piece
          real, pointer :: Dblocko (:, :, :) ! ontop piece

          real, pointer :: Dblocka (:, :, :)
          real, pointer :: Dblockb (:, :, :)
          real, pointer :: Dblockc (:, :, :)

          ! We use this in case we are needing actual energies rather than
          ! matrix elements such as in the Coulomb double-counting correction
          ! energies.
          real E
        end type T_assemble_block

! put each mu by nu block into a neighbor group
        type T_assemble_neighbors
          type(T_assemble_block), pointer :: neighbors (:)
        end type T_assemble_neighbors

! End Module
! ===========================================================================
        end module M_assemble_blocks
