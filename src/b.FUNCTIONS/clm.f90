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

! clm.f90
! Function Description
! ============================================================================
!      This function accepts the value of l and m and return the clm (i.e. the
! Ylm coefficient) for those quantum numbers.
!
! ============================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ============================================================================
!
! Program Declaration
! ===========================================================================
        real function clm (l, m)
        implicit none


! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: l, m

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: lmax = 3

! Local Variable Declaration and Description
! ===========================================================================
! These are the factors which come from the coefficient in the Ylm
! These factors are solved as POST integral coefficients, they are not just
! Ylm coefficients.  The factor of 1/2sqrt(pi) is taken care of in two way.
! First, the 1/2 is taken into account by phifactor
! Second, the 1/sqrt(pi) is taken into account based on whether or not the
! integration over the area has a factor of pi (post)integration
        real fclm (0:lmax, -lmax:lmax)

! Procedure
! ===========================================================================
! The Ylm coefficients
        fclm(0,0)  = 1.0d0

        fclm(1,-1) = sqrt(3.0d0/2.0d0)
        fclm(1,0)  = sqrt(3.0d0)
!       fclm(1,1)  = - sqrt(3.0d0/2.0d0)
        fclm(1,1)  = sqrt(3.0d0/2.0d0)

        fclm(2,-2) = sqrt(15.0d0/8.0d0)
        fclm(2,-1) = sqrt(15.0d0/2.0d0)
        fclm(2,0)  = sqrt(5.0d0/4.0d0)
!       fclm(2,1)  = - sqrt(15.0d0/2.0d0)
        fclm(2,1)  = sqrt(15.0d0/2.0d0)
        fclm(2,2)  = sqrt(15.0d0/8.0d0)

        fclm(3,-3) = sqrt(35.0d0/16.0d0)
        fclm(3,-2) = sqrt(105.0d0/8.0d0)
        fclm(3,-1) = sqrt(21.0d0/16.0d0)
        fclm(3,0)  = sqrt(7.0d0/4.0d0)
!       fclm(3,1)  = - sqrt(21.0d0/16.0d0)
        fclm(3,1)  = sqrt(21.0d0/16.0d0)
        fclm(3,2)  = sqrt(105.0d0/8.0d0)
!       fclm(3,3)  = - sqrt(35.0d0/16.0d0)
        fclm(3,3)  = sqrt(35.0d0/16.0d0)

        clm = fclm(l,m)

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function clm

