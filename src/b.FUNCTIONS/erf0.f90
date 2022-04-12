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

! erf0.f90
! Function Description
! ============================================================================
!       Error function, see Abramowitz and Stegun.
!
! ===========================================================================
! Code written by:
! James P. Lewis
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
        real function erf0 (x)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: x

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: a1 = 0.254829592d0
        real, parameter :: a2 = -0.284496736d0
        real, parameter :: a3 = 1.421413741d0
        real, parameter :: a4 = -1.453152027d0
        real, parameter :: a5 = 1.061405429d0
        real, parameter :: p = 0.3275911d0

! Local Variable Declaration and Description
! ===========================================================================
        real t, y, z

! Procedure
! ===========================================================================
        y = x
        if (y .lt. 0.0d0) y = -y
        if (y .gt. 1.0d-7) then
         t = 1.0d0/(1.0d0 + p*y)
         z = 1.0d0 - exp(-y*y)*(a1*t + a2*t*t + a3*(t**3) + a4*(t**4) + a5*(t**5))
         if (x .lt. 0.0d0) z = -z
        else
         z = 0.0d0
        end if
        erf0 = z

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function erf0
