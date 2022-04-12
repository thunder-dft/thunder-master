! copyright info:
!
!                             @Copyright 2012
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

! M_xc_1c
! Program Description
! ===========================================================================
!      This is a module containing a collection of exchange-correlation
! functionals - energies and potentials along with derivatives of both are
! calculated.
! ===========================================================================
! Code written by:
! Logan Shamberger
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
! Module Declaration
! ===========================================================================
        module M_xc_1c
          use M_precision

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! get_potxc_1c
! ===========================================================================
! Program Description
! ===========================================================================
!  This program will access the requested exchange and correlation
!  functionals from the following list:
!          1  LDA   Wigner
!          2  LDA   Hedin/Lundqvist
!          3  LDA   Ceperley/Alder Perdew/Zunger (1980)
!          4  GGA   Perdew/Wang (1991)
!          5  GGA   Becke (1988) X, Perdew (1986) C
!          6  GGA   Perdew/Burke/Ernzerhof (1996)
!          7  LDA   Zhao/Parr
!          8  LDA   Ceperley/Alder Perdew/Wang (1991)
!          9  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
!         10  GGA   Perdew/Wang (1991) X, Lee/Yang/Parr (1988) C
!         11  LSDA  Volko/Wilk/Nusair (1980)
! The numerical value above is assigned to the variable iexc output
! exchange-correlation potential. This program has been modified to account
! for the fact that the density is a sum of two densities at two different
! centers.  Also, the potential is evaluated at one point in space and the
! integrals (matrix elements) are evaluated elsewhere.
!
! input
!    iexc        xc scheme
!    r           radial coordinate
!    rho         sum of atomic density in au
!    rhop        sum of atomic density gradient (with respect to r) in au
!    rhopp       sum of second gradient (with respect to r) in au
!    rhoz        sum of atomic density gradient (with respect to z) in au
!    rhozz       sum of second gradient (with respect to z) in au
!
! output
!    vpxc        xc potential
!    newexc      xc energy
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
        subroutine get_potxc_1c (iexc, xc_fraction, r, rho, rhop, rhopp,   &
     &                              newexc, vpxc, dnuxc, dnuxcs, dexc)
        implicit none

!! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc

        real, intent (in) :: xc_fraction
        real, intent (inout) :: r        ! radial distance

        ! density and derivatives
        real, intent (inout) :: rho      ! return zero value if rho is small
        real, intent (in) :: rhop
        real, intent (in) :: rhopp

! Output
        real, intent (out) :: newexc
        real, intent (out) :: vpxc
        real, intent (out) :: dnuxc
        real, intent (out) :: dnuxcs
        real, intent (out) :: dexc

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real dummy

        real exc                     ! exchange-correlation energy
        real ex                      ! exchange potential
        real dec, dex                ! derivative of correlation and exchange

        real, dimension (2) :: xpot
        real, dimension (2) :: cpot
        real, dimension (2) :: d

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        dummy = xc_fraction
        dummy = rhop
        dummy = rhopp

! Initialize variables
        newexc = 0.0d0
        vpxc = 0.0d0
        dnuxc = 0.0d0
        dnuxcs = 0.0d0
        dexc = 0.0d0

! If r is really small, then set to manageably small number.
        if (r .lt. 1.0d-4) r = 1.0d-4

! Rho must be positive, but not too small
        if (rho .lt. 1.0d-8) then
          rho = 0.0d0
          return
        else if (rho .lt. 1.0d-5) then
          rho = 1.0d-5
        end if

! Note that we are evaluating the potential at the unpolarized limit.
! We will then treat the spin-polarization as a perturbation.
! Determine exchange-correlation potentials
! exchange (X) only
        if (iexc .eq. 11) then
          zeta = 0.0d0
          d(1) = rho*0.5*(1 + zeta)
          d(2) = rho*0.5*(1 - zeta)
          call lsdavwn (d, dex, dec, xpot, cpot, dnuxc, dnuxcs)
          newexc = dex + dec
          vpxc = xpot(2) + cpot(2)                ! Holds for zeta = 0.0d0 only

! If the improper iexc option was entered then the program will stop.
        else
          write (*,*) ' In get_potxc_1c.f90 - '
          write (*,*) ' stop: xc option not implemented', iexc
          stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
        newexc = dec + dex

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine get_potxc_1c


! lsdavwn.f90
! Program Description
! ===========================================================================
!       This routine computes the exchange and correlation potenials and
! energies for the Vosko, Wilk, Nusair LSDA functional. Each spin component
! is considered.
!
! See
!      S.H. VOSKO and L. WILK and M. NUSAIR
!      Can. J. Phys., 58, 1200 (1980)
!
! ===========================================================================
! Code written by:
! Eduardo Mendez
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
        subroutine lsdavwn (rho, ex, ec, xpot, cpot, dnuxc, dnuxcs)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in), dimension (2) :: rho

! Output
        real, intent (out) :: dnuxc
        real, intent (out) :: dnuxcs
        real, intent (out) :: ec
        real, intent (out) :: ex

        real, intent (out), dimension (2) :: cpot
        real, intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: Ap = 0.0621814d0
        real, parameter :: bp = 3.72744d0
        real, parameter :: cp = 12.9352d0
        real, parameter :: x0p = -0.10498d0

        real, parameter :: Aa = 0.033773728d0
        real, parameter :: ba = 1.13107d0
        real, parameter :: ca = 13.0045d0
        real, parameter :: x0a = -0.00475840d0

        real, parameter :: Af = 0.0310907d0
        real, parameter :: bf = 18.0578d0
        real, parameter :: cf = 7.06042d0
        real, parameter :: x0f = -0.32500d0

! Local Variable Declaration and Description
! ===========================================================================
        real density
        real densitys

        real, dimension (3) :: cdpot
        real, dimension (3) :: xdpot

! spin polarization and derivatives
        real zeta, zp1, zp2, zp1p2, zpp1, zpp2
        real x, xp, xpp
        real g, gp, gpp
        real XXp, XXf, XXa , Qp, Qf, Qa, jp, jf, ja
        real ecP, ecF, ecA, ecPp, ecFp, ecAp, ecPpp, ecFpp, ecApp
        real cte, h, hp, hpp
        real ecpx, ecpz, ecppx, ecppz, ecpxpz
        real d1ec, d2ec, dd1ec, dd2ec, d1d2ec
        real exP, exPp, exPpp
        real expd, expz, exppd, exppz, expdpz
        real d1ex, d2ex, dd1ex, dd2ex, d1d2ex

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! =========================================================================
! Initialize some parameters
        stop ! zpp1 and zpp2 are not set
        density = rho(1) + rho(2)
        densitys = rho(1) - rho(2)
        zeta = densitys/density
        if (density .le. epsilon) then
          zeta = 0.0d0
          ec = 0.0d0
          ex = 0.0d0
          cpot = 0.0d0
          xpot = 0.0d0
          cdpot = 0.0d0
          xdpot = 0.0d0
          return
        end if

! Define simple derivatives
! *************************************************************************
        zp1 = 2.0d0*rho(2)/density**2
        zp2 = -2.0d0*rho(1)/density**2
        zp1p2 = 2.0d0*zeta/density**2

        x = (3.0d0/(4.0d0*pi*density))**(1.0d0/6.0d0)
        xp = - (1.0d0/6.0d0)*x/density
        xpp = 7.0d0*x/(36.0d0*density**2)

        g = (8.0d0/9.0d0)*((1.0d0 + zeta)**(4.0d0/3.0d0)                     &
          + (1.0d0 - zeta)**(4.0d0/3.0d0) - 2.0d0)
        gp = (3.0d0/2.0d0)*((1.0d0 + zeta)**(1.0d0/3.0d0)                    &
           - (1.0d0 - zeta)**(1.0d0/3.0d0))
        gpp = (1.0d0/2.0d0)*((1.0d0 + zeta)**(-2.0d0/3.0d0)                  &
            - (1.0d0 - zeta)**(-2.0d0/3.0d0))

! Intermediate variables
        XXp = x**2.0d0 + bp*x + cp
        XXf = x**2.0d0 + bf*x + cf
        XXa = x**2.0d0 + ba*x + ca
        Qp = (4.0d0*cp - bp*bp)**0.5d0
        Qf = (4.0d0*cf - bf*bf)**0.5d0
        Qa = (4.0d0*ca - ba*ba)**0.5d0
        jp = 2.0d0*log(x - x0p) - log(XXp)                                   &
           + 2.0d0*((2.0d0*x0p + bp)/Qp)*atan(Qp/(2.0d0*x + bp))
        jf = 2.0d0*log(x - x0f) - log(XXf)                                   &
           + 2.0d0*((2.0d0*x0f + bf)/Qf)*atan(Qf/(2.0d0*x + bf))
        ja = 2.0d0*log(x - x0a) - log(XXa)                                   &
           + 2.0d0*((2.0d0*x0a + ba)/Qa)*atan(Qa/(2.0d0*x + ba))

! epsilon derivatives
        ecP = Ap*(2.0d0*log(x) - log(XXp)                                    &
                 + (2.0d0*bp/Qp)*atan(Qp/(2.0d0*x + bp))                     &
                 - (bp*x0p/(x0p*x0p + bp*x0p + cp))*jp)
        ecF = Af*(2.0d0*log(x) - log(XXf)                                    &
                 + (2.0d0*bf/Qp)*atan(Qf/(2.0d0*x + bf))                     &
                 - (bf*x0f/(x0f*x0f + bf*x0f + cf))*jp)
        ecA = Aa*(2.0d0*log(x) - log(XXa)                                    &
                 + (2.0d0*ba/Qa)*atan(Qa/(2.0d0*x + ba))                     &
                 - (ba*x0a/(x0a*x0a + ba*x0a + ca))*ja)

        ecPp = 2.0d0*Ap*cp/(XXp*x) - 2.0d0*Ap*bp*x0p/((x - x0p)*XXp)
        ecFp = 2.0d0*Af*cf/(XXf*x) - 2.0d0*Af*bf*x0f/((x - x0f)*XXf)
        ecAp = 2.0d0*Aa*ca/(XXa*x) - 2.0d0*Aa*ba*x0a/((x - x0a)*XXa)

        ecPpp = - 2.0d0*Ap*cp*(3.0d0*x**2 + 2.0d0*bp*x + cp)/(x*XXp)**2      &
               + 2.0d0*Ap*bp*x0p*((2.0d0*x + bp)*(x - x0p) + XXp)            &
                      /(XXp*(x - x0p))**2
        ecFpp = - 2.0d0*Af*cf*(3.0d0*x**2 + 2.0d0*bf*x + cf)/(x*XXf)**2      &
               + 2.0d0*Af*bf*x0f*((2.0d0*x + bf)*(x - x0f) + XXf)            &
                      /(XXf*(x - x0f))**2
        ecApp = - 2.0d0*Aa*ca*(3.0d0*x**2 + 2.0d0*ba*x + ca)/(x*XXa)**2      &
               + 2.0d0*Aa*ba*x0a*((2.0d0*x + ba)*(x - x0a) + XXa)            &
                      /(XXa*(x - x0a))**2

        cte = 4.0d0/(9.0d0*(2.0d0**(1.0d0/3.0d0) - 1.0d0))

        h = cte*((ecF - ecP)/ecA) - 1.d0
        hp = cte*((ecFp - ecPp)/ecA - (ecF - ecP)*(ecAp/ecA))
        hpp = cte*((ecFpp - ecPpp)/ecA - (ecFp - ecPp)*ecAp/ecA**2           &
                 - (ecF - ecP)*ecApp/ecA - (ecFp - ecPp)*ecAp/ecA            &
                 + (ecF - ecP)*(ecAp/ecA)**2)

! Correlation functional ( and partials to z and x ):
        if (zeta .ne. 0.0d0) then
         ec = ecP + ecA*g*(1 + h*zeta**4)
        else
         ec = ecP
        end if

        ecpx = ecPp + ecAp*g*(1 + h*zeta**4) + ecA*g*hp*zeta**4
        ecpz = ecA*gp*(1.0d0 + h*zeta**4) + ecA*g*h*4*zeta**3

        ecppx = ecPp + ecApp*g*(1.0d0 + h*zeta**4) + 2.0d0*ecAp*g*hp*zeta**4 &
              + ecA*g*hpp*zeta**4
        ecppz = ecA*gpp*(1.0d0 + h*zeta**4) + ecA*gp*h*zeta**3               &
              + ecA*g*h*12.0d0*zeta**2
        ecpxpz = ecAp*gp*(1.0d0 + h*zeta**4) + ecA*gp*hp*zeta**4             &
               + ecAp*g*h*4.0d0*zeta**3 + ecA*g*hp*4.0d0*zeta**3

! Partial derivatives VWN exchanche functional
        d1ec = xp*ecpx + zp1*ecpz
        d2ec = xp*ecpx + zp2*ecpz

! Second partial derivatives
        dd1ec = xp**2*ecpp + 2.0d0*xp*zp1*ecpxpz + xpp*ecpx                  &
              + zp1*zp1*ecppz + zpp1*ecpz
        dd2ec = xp**2*ecpp + 2.0d0*xp*zp2*ecpxpz + xpp*ecpx                  &
              + zp2*zp1*ecppz + zpp2*ecpz
        d1d2ec = xp**2*ecpp+ xp*(zp1 + zp2)*ecpxpz + xpp*ecpx                &
               + zp1*zp2*ecppz + zp1p2*ecpz

! ****************************************************************************
!
!       VNN EXCHANGE FUNCTIONAL
!
! ****************************************************************************
        exP = (-3.0d0/2.0d0)*(3.0d0*density/pi)**(1.0d0/3.0d0)
        exPp = exP/(3.0d0*density)
        exPpp = - 2.0d0*exP/(3.0d0*density)**2

        ex =(1.0d0 + 4.0d0*g/9.0d0)*exP
        expd = ex/(3.0d0*density)
        exppd = -2.0d0*ex/(9.0d0*density**2)
        expz = exP*gp
        exppz = exP*gpp
        expdpz = exPp*gp

        d1ex = expd + zp1*expz
        d2ex = expd + zp2*expz

        dd1ex = exppd + 2.0d0*zp1*expdpz + expd + zp1*zp1*exppz + zpp1*expz
        dd2ex = exppd + 2.0d0*zp2*expdpz + expd + zp2*zp2*exppz + zpp2*expz
        d1d2ex = exppd + (zp1 + zp2)*expdpz + expd + zp1*zp2*exppz + zp1p2*expz

! Functions in Rydberg units - divide by factor of 2 to get Hartree
! ****************************************************************************
        xpot(1) = 0.5d0*(density*d1ex + ex)
        xpot(2) = 0.5d0*(density*d2ex + ex)
        cpot(1) = 0.5d0*(density*d1ec + ec)
        cpot(2) = 0.5d0*(density*d2ec + ec)
        ex = 0.5d0*ex
        ec = 0.5d0*ec

        cdpot(1) = 0.5d0*dd1ec
        cdpot(2) = 0.5d0*d1d2ec
        cdpot(3) = 0.5d0*dd2ec
        xdpot(1) = 0.5d0*dd1ex
        xdpot(2) = 0.5d0*d1d2ex
        xdpot(3) = 0.5d0*dd2ex

        dnuxc = 0.25d0*density*(xdpot(1) + 2.0d0*xdpot(2) + xdpot(3))        &
              + 0.5d0*(d1ec + d2ec) + 0.5d0*(d1ex + d2ex)                   &
              + 0.25d0*density*(cdpot(1) + 2.0d0*cdpot(2) + cdpot(3))

        dnuxcs = 0.25d0*density*(xdpot(1) - 2.0d0*xdpot(2) + xdpot(3))       &
                + 0.5d0*(d1ec - d2ec) + 0.5d0*(d1ex - d2ex)                 &
                + 0.25d0*density*(cdpot(1) - 2.0d0*cdpot(2) + cdpot(3))

        dnuxcs = 0.5d0*(ecA + 4.0d0*exP/9.0d0)/density

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================

        return
        end subroutine lsdavwn

! End Module
! =============================================================================
        end module
