! copyright info:
!
!                             @Copyright 2012
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

! Dsmoother.f90
! Function Description
! ============================================================================
!       This is the derivative of the smoother function.  The smoothing
! function is:
!
!       smoother(r,rbegin,rend).
!
! We define our final functions on a grid as:
!
!       smoother(r)*exact + (1 - smoother(r))*longrange.
!
! We also define x = ratio = (r - rbegin)/(rend - rbegin) in some equations.
!
! The following are the requirements for a sane smoother function:
!
!       smoother(rend) = 0
!       smoother(rbegin) = 1
!
! The following are nice to have:
! 1) You might like symmetry: stn(x) = 1 - stn(1 - x)
! 3) Quick to evaluate on the computer: not (cos(pi*x)+1)/2.
!
! There are different functional forms: (1-x^n)^m is the "old" method
! which obeys the requirements for n,m > 1.  It is the smoothest
! for n, m = 2 (smallest possible values).  Thus, a generic forth-order
! polynomial might be a better choice, since n, m = 2 is 1 - 2x^2 + x^4
! (no x or x^3 terms).  A forth-order polynomial is the "new" method, but the
! above conditions force it to be of the form:
!
!       1 + 0*x + (n-3)x^2 + (2 - 2n)x^3 + nx^4
!
! Where n must be within [-3,3], inclusive.
!
! For n = 1, you recover the "old" method.
! To minimize the derivative, n is set to zero.
! To make symmetric, n is set to zero
! To get second-derivatives to be zero at x = 0, n is to 3
! To get second-derivatives to be zero at x = 1, n is to -3
!
! Finally, if you want to get second-derivatives to be zero at x = 0, 1 because
! you are doing a Hessian, then you need a higher-order function.
! 1 + (-10-n)x^3 + (15+3n)x^4 + (-6-3n)x^5 + nx^6
! For symmetry, and mininizing the derivative, n = 0 is probably best.
! This is not implemented.
!
! Now a word about choosing x. Since,
! E = (short range)*f + (long range)*(1 - f), then this is well-behaved
! (for x > 0), since E is bound by (short range) and (long range). But what
! about the derivatives (i.e. forces):
! E' = (short range)'*f + (long range)'*(1 - f) + (extra term)?
! which is well-behaved, except for the extra term which makes it:
! E' = (short range)'*f + (long range)'*(1-f) + f'*(short range - long range).
! Keeping f' small is good, which is achieved by making x smaller.
! Keeping (short range - long range) small is good also, which is achieved
! by making x bigger. Clearly some sort of balance must be achieved,
! when choosing x.
!
! ============================================================================
! Code written by:
! James P. Lewis
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
        real function Dsmoother (r, rend, x)
        implicit none

        include "../include/constants.h"

! Argument Declaration and Description
! ===========================================================================
        real, intent (in) :: r
        real, intent (in) :: rend
        real, intent (in) :: x

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: scaler = 1

! Local Variable Declaration and Description
! ===========================================================================
        real ratio
        real rbegin

! Procedure
! ===========================================================================
        rbegin = x*rend
        if (r .gt. rend) then
          Dsmoother = 0.0d0
        else if (r .lt. rbegin) then
          Dsmoother = 1.0d0
        else
          ratio = (r - rbegin)/(rend - rbegin)
          Dsmoother = 2.0d0*(1.0d0 - ratio**2)*(2.0d0*ratio)/(rend - rbegin)
! two methods - if you want newer method then uncomment these lines
!         Dsmoother = 2.0d0*(scaler - 3)*ratio + 3.0d0*(2 - 2*scaler)*ratio**2 &
!     &                + 4.0d0*scaler*ratio**3
        end if

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function Dsmoother

