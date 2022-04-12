! copyright info:
!
!                             @Copyright 2013
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
! James P. Lewis
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
     &                           newexc, vpxc, dnuxc, dnuxcs, dexc)
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
        real fx
        real fxc

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

! XC Ceperley - Alder
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
        if (iexc .eq. 3) then
          call ceperley_alder (rho, ex, fx, exc, fxc, dexc, dnuxc)
          vpxc = fxc
          dex = ex
          dec = exc - ex

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


! ===========================================================================
! ceperley_alder
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine compute the ceperley-alder form of the LDA as
! parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981).  The units
! of this program are in atomic units, so the density but be changed to atomic
! units after input and the final answer converted to eV-Angstrom units.
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
        subroutine ceperley_alder (rho_in, epsx, potx, epsxc, potxc, drvexc, &
                                   dpotxc)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: rho_in

! Output
        real, intent (out) :: dpotxc
        real, intent (out) :: drvexc
        real, intent (out) :: epsx
        real, intent (out) :: epsxc
        real, intent (out) :: potx
        real, intent (out) :: potxc

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        real density
        real densityp
        real densitypp
        real dpotc
        real dpotx
        real depsc
        real ddepsc
        real rho
        real rs
        real rsl
        real sqrs

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Convert density to Angstrom units.
        rho = rho_in

        if (rho .le. epsilon) then
          epsx = 0.d0
          potx = 0.d0
          epsxc = 0.d0
          potxc = 0.d0
          dpotxc = 0.0d0
          drvexc = 0.0d0
          return
        end if

! Initialize some constants related to density.
        rs = 0.62035049d0/rho**(1.0d0/3.0d0)
        if (rho .lt. 0.23873241d0) then
          sqrs = sqrt(rs)
          density = 1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs
          epsxc = -0.4581652d0/rs - 0.1423d0/density
          potxc = epsxc - rs*(0.15273333d0/rs**2                             &
     &                  + (0.02497128d0/sqrs + 0.01581427d0)/density**2)

          densityp =  1.0529d0/(2.0d0*sqrs) + 0.3334d0
          densitypp = -0.5d0*1.0529d0/(2.0d0*rs*sqrs)
          depsc = 0.1423d0*densityp/(density*density)
          ddepsc = - 2.0d0*0.1423d0*densityp*densityp/density**3             &
     &             + 0.1423d0*densitypp/(density*density)
        else
          rsl = log(rs)
          epsxc = - 0.4581652d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs   &
     &            + 0.002d0*rs*rsl
          potxc = epsxc - rs*(0.15273333d0/rs**2 + 0.01036667d0/rs           &
     &                  - 0.003866667d0 + 0.00066667d0*(1.0d0 + rsl))
          depsc = 0.0311d0/rs - 0.0116d0 + 0.0020d0*(rsl + 1.0d0)
          ddepsc = -0.0311d0/(rs*rs) + 0.0020d0/rs
        end if

! Exchange-only energy and potential
        epsx = - 0.7385587664d0*rho**(1.0d0/3.0d0)
        potx = 4.0d0/3.0d0*epsx

! Extended hubbard additions.
        drvexc = (potxc - epsxc)/rho

! Now dpotxc; we compute dpot/dn. We use dpot/dn = 2*dexc/dn + n*d2(exc)/dn2.
! Here dexc/dn = drvexc, and
! dpot/dn = (-2/(9*n*n))*ex + 4*rs/(9*n*n)*dec + rs*rs/(9*n*n)*ddec
! Let dpotc = dpot/dn = 4*rs/(9*n*n)*dec + rs*rs/(9*n*n)*ddec
        dpotc = (4.0d0*rs/(9.0d0*rho*rho))*depsc                             &
     &         + (rs*rs/(9.0d0*rho*rho))*ddepsc
        dpotx = - (2.0d0/(9.0d0*rho*rho))*epsx
        dpotxc = 2.0d0*drvexc + rho*(dpotx + dpotc)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine ceperley_alder


! End Module
! =============================================================================
        end module
