! copyright info:
!
!                             @Copyright 2023
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

! Module Description
! ============================================================================
!      This module initializes and finalizes the MPI space.
! However, it is a dummy routine.
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
! Module declaration
! ============================================================================
        module M_mpi
        use mpi

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer my_proc
        integer nprocessors

        logical iammaster
        logical iammpi

! module procedures
        contains

! ===========================================================================
! initialize_mpi
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine initializes the MPI space; however, it is a
! dummy routine.
!
! Program Declaration
! ===========================================================================
        subroutine initialize_mpi
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        my_proc = 0
        nprocessors = 1
        iammaster = .true.
        iammpi = .false.

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_mpi


! ===========================================================================
! finalize_mpi
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine finalizes the MPI space; however, it is
! a dummy routine.
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
        subroutine finalize_mpi
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end subroutine finalize_mpi

! End the module
! ===========================================================================
        end module M_mpi
