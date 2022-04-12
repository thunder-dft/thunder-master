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

! Ylm.f90
! Function Description
! ============================================================================
!      This function accepts the value of l and m and return the value of
! Ylm for those quantum numbers.  The value of Ylm is returned for a given
! value of x, y, z (so, for Cartesian coordinate system).
!
! --------------------------------------------------------------------
!                -------- Spherical harmonics -------
! --------------------------------------------------------------------
! s-orbital (l = 0)
!     Y0m = sqrt(1/4*Pi)
!
!     d(Y0m)/dr = 0.d0
!
! --------------------------------------------------------------------
! p-orbital (l = 1)
!     Y1m = sqrt(3/4*Pi)/r * Y1(r)
!            m=-1:   Y1(r) = y
!            m= 0:   Y1(r) = z
!            m=+1:   Y1(r) = x
!
!     d(Y1m)/dr = sqrt(2/4*Pi)d(Y1m(r)/r)/dr = 0.d0
!
! --------------------------------------------------------------------
! d-orbital (l = 2)
!     Y2m = sqrt(15/4*Pi)/r**2 * Y2m(r)
!            m=-2:   Y2(r) = x*y
!            m=-1:   Y2(r) = y*z
!            m= 0:   Y2(r) = (3z**2-r**2)/sqrt(12)
!            m=+1:   Y2(r) = x*z
!            m=+2:   Y2(r) = (x**2-y**2)/2
!
!     d(Y2m)/dr = sqrt(15/4*Pi)/r**2 * d(Y2m(r))/dr = 0.d0
!
! ============================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
!
! Program Declaration
! ===========================================================================
        real function Ylm (r, l, m)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: l, m

        real, intent(in), dimension (3) :: r   !< vector of r (point on grid)

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: lmax = 3
        double precision, parameter :: pi = 3.141592653589793238462643d0  ! pi = 4.0d0*atan(1.0d0)

! Local Variable Declaration and Description
! ===========================================================================
        real x, y, z
        real r_magnitude

        interface
          function clm (l, m)
            integer, intent (in) :: l, m
            real clm
          end function clm
        end interface

! Procedure
! ===========================================================================

! The Ylm functions
        if (l .gt. lmax) then
          write (*,*) ' Error: Wrong quantum number l !'
          stop
        end if

! Set the values of x, y, and z
        x = r(1)
        y = r(2)
        z = r(3)
        r_magnitude = sqrt(x**2 + y**2 + z**2)
        if (r_magnitude .lt. 1.0d-4) Ylm = 0.0d0

        select case ( l )
! -------------------------------------------------------------
! s-orbital
! -------------------------------------------------------------
        case (0)
           Ylm = clm(0,0)/sqrt(4.0d0*pi)

! -------------------------------------------------------------
! p-orbital
! -------------------------------------------------------------
        case (1)
! values
          if (m .eq. -1) Ylm = clm(1,-1)/sqrt(2.0d0*pi)*(y/r_magnitude)
          if (m .eq. 0) Ylm = clm(1,0)/sqrt(4.0d0*pi)*(z/r_magnitude)
          if (m .eq. 1) Ylm = clm(1,1)/sqrt(2.0d0*pi)*(x/r_magnitude)

! -------------------------------------------------------------
! d-orbital
! -------------------------------------------------------------
        case (2)
! values
           if (m .eq. -2) Ylm = clm(2,-2)/sqrt(2.0d0*pi)*(x*y/r_magnitude**2)
           if (m .eq. -1) Ylm = clm(2,-1)/sqrt(2.0d0*pi)*(y*z/r_magnitude**2)
           if (m .eq. 0) Ylm = clm(2,0)/sqrt(4.0d0*pi)*(2.0d0*z**2 - x**2 - y**2)/r_magnitude**2
           if (m .eq. 1) Ylm = clm(2,1)/sqrt(2.0d0*pi)*(x*z/r_magnitude**2)
           if (m .eq. 2) Ylm = clm(2,2)/sqrt(2.0d0*pi)*(x**2 - y**2)/r_magnitude**2
        end select

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function Ylm

