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

! mult3x3.f90
! Program Description
! ===========================================================================
!       This subroutine will multiply two 3x3 matrices.
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
        subroutine mult3x3 (a, b)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent (in), dimension (3, 3) :: a

! Output:
        real, intent (inout), dimension (3, 3) :: b

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real, dimension (3, 3) :: c

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! 1. colums
       c(1,1) = a (1,1)*b(1,1) + a(1,2)*b(2,1) + a(1,3)*b(3,1)
       c(1,2) = a (1,1)*b(1,2) + a(1,2)*b(2,2) + a(1,3)*b(3,2)
       c(1,3) = a (1,1)*b(1,3) + a(1,2)*b(2,3) + a(1,3)*b(3,3)
! 2. column
       c(2,1) = a (2,1)*b(1,1) + a(2,2)*b(2,1) + a(2,3)*b(3,1)
       c(2,2) = a (2,1)*b(1,2) + a(2,2)*b(2,2) + a(2,3)*b(3,2)
       c(2,3) = a (2,1)*b(1,3) + a(2,2)*b(2,3) + a(2,3)*b(3,3)
! 3. column
       c(3,1) = a (3,1)*b(1,1) + a(3,2)*b(2,1) + a(3,3)*b(3,1)
       c(3,2) = a (3,1)*b(1,2) + a(3,2)*b(2,2) + a(3,3)*b(3,2)
       c(3,3) = a (3,1)*b(1,3) + a(3,2)*b(2,3) + a(3,3)*b(3,3)

! copy the result into b matrix
       b = c

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine mult3x3
