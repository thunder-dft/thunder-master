! copyright info:
!
!                             @Copyright 2013
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

! simpson.f90
! Program Description
! ===========================================================================
!       Computes integral of f(x) using simple rectangular gridding.
!
! ===========================================================================
! Code written by:
! Hong Wang and James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function simpson (n, f, dx)
        implicit none

        real simpson

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: n

        real, intent(in) :: dx
        real, intent(in), dimension (n) :: f

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint

! Procedure
! ===========================================================================
        simpson = 0.0d0
        do ipoint = 1, n
          if (ipoint .ne. 1 .or. ipoint .ne. n) then
            simpson = simpson + f(ipoint)*dx
          else
            simpson = simpson + f(ipoint)*0.5d0*dx
          end if
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end function simpson
