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

! M_assemble_ewald
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the interactions for the short range and long range
!! ewald interactions.
!!
!! It contains the following subroutines within the module:
!!
!!       Dassemble_ewaldsr.f90 - assemble the short-range ewald matrix
!!       Dassemble_ewaldlr.f90 - assemble the long-range ewald matrix
!!
!!       These routines really do nothing because this is Harris which
!! has no charge dependency.  These routines just set the interactions
!! to zero.
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_Dassemble_ewald

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones
        use M_rotations
        use M_Drotations

! /FDATA
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
! Dassemble_ewaldsr.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is a dummy routine.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_ewaldsr (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                   !< counter over atoms
        integer num_neigh               !< number of neighbors

        ! forces
        type(T_forces), pointer :: pfi

! Allocate Arrays
! ===========================================================================
!       allocate (s%ewaldsr (s%natoms))

! Procedure
! ===========================================================================
! We build the ewald forces here, so we allocate and initialize
        do iatom = 1, s%natoms
          num_neigh = s%neighbors(iatom)%neighn
          pfi=>s%forces(iatom)
          allocate (pfi%ewaldsr (3, num_neigh)); pfi%ewaldsr = 0.0d0
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
        end subroutine Dassemble_ewaldsr


! ===========================================================================
! Dassemble_ewaldlr.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is a dummy routine.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_ewaldlr (s)
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
! None

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! None

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_ewaldlr


! ===========================================================================
! destroy_assemble_ewald
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is a dummy routine.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_Dassemble_ewald (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                            !< counter over atoms

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          deallocate (s%forces(iatom)%ewaldsr)
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
        end subroutine destroy_Dassemble_ewald

! End Module
! ===========================================================================
        end module M_Dassemble_ewald
