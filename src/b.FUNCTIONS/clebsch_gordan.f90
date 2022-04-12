! copyright info:
!
!                             @Copyright 2008
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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

! clebsch_gordan.f90
! Program Description
! ===========================================================================
!       This routine calculates the Clebsch-Gordan coefficients which
! are represented by - <l1,l2;m1,m2|l,m>.
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
        function clebsch_gordan (l1, m1, l2, m2, l, m)
        implicit none
        real clebsch_gordan

! Argument Declaration and Description
! ===========================================================================
        integer l, l1, l2
        integer m, m1, m2
!
! Local Parameters and Data Declaration
! ===========================================================================
! The maximum z value is to be only 6, since the lmax value is only 6 in the
! routine which calls this function.
        integer izmax
        parameter (izmax = 50)

! Local Variable Declaration and Description
! ===========================================================================
        integer iz

        real piece1
        real piece2
        real piece3
        real part3

        interface
          function factorial (i)
            integer, intent (in) :: i
            integer factorial
          end function factorial
        end interface

! Procedure
! ===========================================================================
! Initialize the coefficient to zero.
        clebsch_gordan = 0.0d0

! First determine that the values of the l1, l2, and l satisfiy the
! triangle equation - |l1-l2| =< l >= l1 + l2.
        if (l .lt. abs(l1 - l2)) then
          return
        end if
        if (l .gt. (l1 + l2)) then
          return
        end if

! The other condition is that m = m1 + m2
        if (m .ne. (m1 + m2)) return

! The clebsch_gordon coefficient will be written as a product of three
! pieces. So clebsch_gordan = piece1*piece2*piece3
        piece1 = ((2.0d0*l + 1)*factorial(l1 + l2 - l)                      &
     &            *factorial(l1 - l2 + l)*factorial(-l1 + l2 + l))          &
     &          /factorial(l1 + l2 + l + 1)
        piece1 = sqrt(piece1)

        piece2 = factorial(l1 + m1)*factorial(l1 - m1)                      &
     &          *factorial(l2 + m2)*factorial(l2 - m2)                      &
     &          *factorial(l + m)*factorial(l - m)
        piece2 = sqrt(piece2)

        piece3 = 0.0d0
        do iz = 0, izmax
          if (((l1+l2-l-iz) .ge. -0.1) .and. ((l1-m1-iz) .ge. -0.1)         &
     &        .and. ((l2 + m2 - iz) .ge. -0.1)                              &
     &        .and. ((l - l2 + m1 + iz) .ge. -0.1)                          &
     &        .and. ((l - l1 - m2 + iz) .ge. -0.1)) then
            part3 = factorial(iz)*factorial(l1 + l2 - l - iz)               &
     &             *factorial(l1 - m1 - iz)*factorial(l2 + m2 - iz)         &
     &             *factorial(l - l2 + m1 + iz)*factorial(l - l1 - m2 + iz)
            piece3 = piece3 + (-1)**iz/part3
          end if
        end do

! Now calculate the coefficient
        clebsch_gordan = piece1*piece2*piece3

! Format Statements
! ===========================================================================

        return
        end function clebsch_gordan
