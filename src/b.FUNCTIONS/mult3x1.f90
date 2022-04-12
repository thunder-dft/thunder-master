! copyright info:
!
!                             @Copyright 2014
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

! mult3x1.f90
! Program Description
! ===========================================================================
!       This subroutine will multiply a matrix and a vector together.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
!
! Program Declaration
! ===========================================================================
        subroutine mult3x1 (a, b)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent (in), dimension (3, 3) :: a

! Output:
        real, intent (inout), dimension (3) :: b

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real, dimension (3) :: c

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! 1. column
       c(1) = a(1,1)*b(1) + a(1,2)*b(2) + a(1,3)*b(3)
       c(2) = a(2,1)*b(1) + a(2,2)*b(2) + a(2,3)*b(3)
       c(3) = a(3,1)*b(1) + a(3,2)*b(2) + a(3,3)*b(3)

! copy the result into b matrix
       b = c

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine mult3x1
