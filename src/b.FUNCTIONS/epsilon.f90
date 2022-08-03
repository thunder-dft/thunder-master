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

! epsilon.f90
! Subroutine Description
! ===========================================================================
! input: r1, r2
! output: epsilon the metric tensor
!
! note: the third column of epsilon is eta(3)
!
! spe = epsilon backwards.
!
! R1vector points toward the point O while R2vector points
!   away from the point O.
!
!                         *O    (XP,YP,ZP)
!                      *    *
!       R1VECTOR    *        *  R2VECTOR
!                *            *
!             *                *
! (X,Y,Z)  *                    *
!       *
!
!            |  ^     ^       ^     ^        ^     ^   |
!            |  X-dot-XP      X-dot-YP       X-dot-ZP  |
!            |                                         |
!            |  ^     ^       ^     ^        ^     ^   |
!      spe = |  Y-dot-XP      Y-dot-YP       Y-dot-ZP  |
!            |                                         |
!            |  ^     ^       ^     ^        ^     ^   |
!            |  Z-dot-XP      Z-dot-YP       Z-dot-ZP  |
!            |                                         |
!
! ===========================================================================
! Code written by:
! James P. Lewis
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
        subroutine epsilon_function (r1, r2, spe)
        use M_species
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: r1 (3)      !< Input vector r1
        real, intent (in) :: r2 (3)      !< Input vector r2

! Output
        real, intent (out) :: spe (3,3)  !< Output spe (epsilon backwards)

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ix, jx, kx                    ! dimension counter

        real r1_mag, r2_mag                   ! magnitudes of r1 and r2
        real yp_mag                           ! magnitude of unit yphat
        real unit

! unit vectors, zphat is unit vector of r2
        real, dimension (3) :: r1_hat
        real, dimension (3) :: xp_hat, yp_hat, zp_hat

        interface
          function a_cross_b (a, b)
            real, dimension (3) :: a_cross_b
            real, intent (in), dimension (3) :: a, b
          end function a_cross_b

          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance

          function magnitude (a)
            real magnitude
            real, intent(in), dimension (3) :: a
          end function magnitude
        end interface

! Procedure
! ===========================================================================
        r1_mag = magnitude (r1)
        r2_mag = magnitude (r2)
        if (r2_mag .lt. 1.0d-4) then
!         r2vector = (0,0,0) ----- set eps = unit matrix
          spe = 0.0d0
          do ix = 1, 3
            spe(ix,ix) = 1.0d0
          end do
          open (11, file = 'WARNINGS', status = 'unknown', position = 'append')
          write (11,*) ' ******* r2 = 0 : spe = unit matrix! '
          close (11)
          return
        end if

! zphat lies along r2vector
        zp_hat = r2/r2_mag

! yphat = zphat-cross-r1hat
        yp_mag = 0
        if (r1_mag .gt. 1.0d-4) then
          r1_hat = r1/r1_mag
          yp_hat = a_cross_b (zp_hat, r1_hat)
          yp_mag = magnitude (yp_hat)
        end if

! If zphat and r1hat are colinear or r1vector = (0,0,0)
! Find the first non-zero component of zphat:
        if (yp_mag .lt. 1.0d-4) then
          if (abs(zp_hat(1)) .gt. 1.0d-4)then
! zphat(1) not equal to zero
            yp_hat(1) = -(zp_hat(2) + zp_hat(3))/zp_hat(1)
            yp_hat(2) = 1.0d0
            yp_hat(3) = 1.0d0
            yp_mag = magnitude (yp_hat)
          else if(abs(zp_hat(2)) .gt. 1.0d-4) then
! zphat(2) not equal to zero
            yp_hat(1) = 1.0d0
            yp_hat(2) = -(zp_hat(1) + zp_hat(3))/zp_hat(2)
            yp_hat(3) = 1.0d0
            yp_mag = magnitude (yp_hat)
          else
! zphat(3) not equal to zero
            yp_hat(1) = 1.0d0
            yp_hat(2) = 1.0d0
            yp_hat(3) = -(zp_hat(1) + zp_hat(2))/zp_hat(3)
            yp_mag = magnitude (yp_hat)
          end if
        end if
        yp_hat = yp_hat/yp_mag

! find pihat
        xp_hat = a_cross_b (yp_hat, zp_hat)

! find epsilon matrix
        spe(:,1) = xp_hat(:)
        spe(:,2) = yp_hat(:)
        spe(:,3) = zp_hat(:)

! test by computing spe*spe(dagger)
        do ix = 1, 3
          do jx = 1, 3
            unit = 0.0d0
            do kx = 1, 3
              unit = unit + spe(ix,kx)*spe(jx,kx)
            end do
            if (ix .eq. jx .and. abs(unit - 1.d0) .gt. 1.0d-4) then
              open (11, file = 'WARNINGS', status = 'unknown')
              write (11,*) ' ******* error in epsilon spe*spedag = ', unit
              close (11)
            end if
            if (ix .ne. jx .and. abs(unit) .gt. 1.0d-4) then
              open (11, file = 'WARNINGS', status = 'unknown')
              write (11,*) ' ******* error in epsilon spe*spedag = ', unit
              close (11)
            end if
          end do
        end do

! format statements
! ===========================================================================
        return
      end subroutine epsilon_function
