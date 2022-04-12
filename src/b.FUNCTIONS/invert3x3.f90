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

! invert3x3.f90
! Program Description
! ===========================================================================
!       This subroutine will calculate the inverse of a 3x3 matrix.
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
        subroutine invert3x3 (amatrix, ainverse)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent (in), dimension (3, 3) :: amatrix

! Output:
        real, intent (out), dimension (3, 3) :: ainverse

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real determinant

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
        determinant =  &
     &    amatrix(1,1)*(amatrix(2,2)*amatrix(3,3) - amatrix(2,3)*amatrix(3,2))&
     &  + amatrix(1,2)*(amatrix(3,1)*amatrix(2,3) - amatrix(2,1)*amatrix(3,3))&
     &  + amatrix(1,3)*(amatrix(2,1)*amatrix(3,2) - amatrix(3,1)*amatrix(2,2))

! Now calculate inverse of inertia tensor
! ainv(i,j) = (-1**(i+j)) * cofactor(j,i) / det(a)
        if (abs(determinant) .gt. 1.0d-5) then
          ainverse(1,1) =                                                    &
     &     (amatrix(2,2)*amatrix(3,3) - amatrix(3,2)*amatrix(2,3))/determinant
          ainverse(2,1) =                                                    &
     &     - (amatrix(2,1)*amatrix(3,3) - amatrix(3,1)*amatrix(2,3))/determinant
          ainverse(3,1) =                                                    &
     &     (amatrix(2,1)*amatrix(3,2) - amatrix(3,1)*amatrix(2,2))/determinant
          ainverse(1,2) =                                                    &
     &     - (amatrix(1,2)*amatrix(3,3) - amatrix(3,2)*amatrix(1,3))/determinant
          ainverse(2,2) =                                                    &
     &     (amatrix(1,1)*amatrix(3,3) - amatrix(3,1)*amatrix(1,3))/determinant
          ainverse(3,2) =                                                    &
     &     - (amatrix(1,1)*amatrix(3,2) - amatrix(3,1)*amatrix(1,2))/determinant
          ainverse(1,3) =                                                    &
     &     (amatrix(1,2)*amatrix(2,3) - amatrix(1,3)*amatrix(2,2))/determinant
          ainverse(2,3) =                                                    &
     &     - (amatrix(1,1)*amatrix(2,3) - amatrix(1,3)*amatrix(2,1))/determinant
          ainverse(3,3) =                                                    &
     &     (amatrix(2,2)*amatrix(1,1) - amatrix(1,2)*amatrix(2,1))/determinant
        else
          ainverse = 0.0d0
          open (11, file = 'WARNINGS', status = 'unknown', position = 'append')
          write (11,*) ' ********* WARNING ********* '
          write (11,*) ' The determinant of the amatrix in invert3x3 is '
          write (11,*) ' equal to zero. Be careful and make sure that '
          write (11,*) ' you really want to continue. '
          write (11,*) ' If you have a dimer molecule, then you probably '
          write (11,*) ' want to set iconstraint_L = 0 '
          close (11)
          stop
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine invert3x3
