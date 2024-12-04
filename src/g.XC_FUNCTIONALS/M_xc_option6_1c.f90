! copyright info:
!
!                             @Copyright 2022
!                           Fireball Committee
! Hong Kong Quantum AI Laboratory, Ltd. - James P. Lewis, Chair
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek
! Arizona State University - Otto F. Sankey

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! California Institute of Technology - Brandon Keith
! Czech Institute of Physics - Prokop Hapala
! Czech Institute of Physics - Vladimír Zobač
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Synfuels China Technology Co., Ltd. - Pengju Ren
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

! M_xc
! Program Description
! ===========================================================================
!      This is a module containing a collection of exchange-correlation
! functionals - energies and potentials along with derivatives of both are
! calculated.
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
! Module Declaration
! ===========================================================================
        module M_xc_1c

! /GLOBAL
        use M_precision

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! get_potxc1c.f90
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
!         11  LSDA  Vosko/Wilk/Nusair (1980)
!         12  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
!                   with exact exchange
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
!
! output
!    vpxc        xc potential
!    newexc      xc energy
!    dnuxc
!    dnuxcs
! ===========================================================================
! Code written by:
! James P. Lewis
! Hong Kong Quantum AI Laboratory
! Room 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, Hong Kong
! Office Telephone +852 9865 3077
!
! Program Declaration
! ===========================================================================
        subroutine get_potxc_1c (iexc, xc_fraction, r, rho, rhop, rhopp,   &
     &                           newexc, vpxc, dnuxc, dnuxcs, dexc)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc

        real, intent (in) :: xc_fraction
        real, intent (inout) :: r

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
        integer ix

        real dec
        real dex

! density and derivatives - spin cases
        real, dimension (2) :: d, dp, dpp

! exchange and correlation potentials - for spin
        real, dimension (2) :: xpot
        real, dimension (2) :: cpot

        ! output of easypbe
!       real expbe, vxuppbe, vxdnpbe, ecpbe, vcuppbe, vcdnpbe

! Procedure
! ===========================================================================
! Initialize to zero.
        newexc = 0.0d0
        vpxc = 0.0d0
        dnuxc = 0.0d0
        dnuxcs = 0.0d0
        dexc = 0.0d0

! If r is really small, then set to manageably small number.
        if (r .lt. 1.0d-4) r = 1.0d-4

! rho must be positive, but not too small
        if (rho .lt. 1.0d-8) then
          rho = 0.0d0
          return
        else if (rho .lt. 1.0d-5) then
          rho = 1.0d-5
        end if

! Determine exchange-correlation potentials
! Perdew-Burke-Ernzerhof GGA
! as parameterized by J.P. Perdew, K. Burke, M. Ernzerhof, Phys Rev Lett 77, 3865 (1996)
        if (iexc .eq. 6) then

! X PBE GGA by default
          ix = 5

! Inputs for easypbe
! input densities - d(1) = dup; d(2) = ddn
! derivative of input densities
          d = 0.5d0*rho
          dp = 0.5d0*rhop
          dpp = 0.5d0*rhopp
          call ggaxrad_1c (ix, r, d, dp, dpp, xpot, dex)
          call ggacrad_1c (5, r, d, dp, dpp, cpot, dec)
          vpxc = xpot(1) + cpot(1)

! If the improper iexc option was entered then the program will stop.
        else
          write (*,*) ' In get_potxc_1c.f90 - '
          write (*,*) ' You must recomplile create.x for iexc = 6 '
          write (*,*) ' Set XC = BLYP in include/OPTIONS and recompile. '
          stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
! This comment seems to be old fashioned
        newexc = dec + dex

! Format Statements
! ===========================================================================
        return
        end subroutine get_potxc_1c


! ===========================================================================
! ggaxrad_1c.f90
! ===========================================================================
! Program Description
! ===========================================================================
!      This routine calculates the Perdew-Burke-Ernzerhof exchange potential
! and energy density. Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 5    GGA-X PBE
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!
! ===========================================================================
! Code rewritten to FORTRAN 90 by:
! James P. Lewis
! Campus Box 7260
! Department of Biochemistry and Biophysics
! University of North Carolina
! Chapel Hill, NC 27599-7260
! FAX 919-966-2852
! Office telephone 919-966-4644
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ggaxrad_1c (mode, rin, rho, rhop, rhopp, xpot, xen)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp

! Output
        real, intent (out) :: xen

        real, intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ispin

        real density, densityp, densitypp
        real ex
        real fermik
        real r
        real s, u, v
        real vx

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! If r is really small, then set to manageably small number.
        r = rin
        if (rin .lt. 1.0d-4) r = 1.0d-4

! exchange GGA, loop for up & down spin
        xen = 0.0d0
        do ispin = 1, 2
          if (rho(ispin) .le. epsilon) then
            xpot(ispin) = 0.0d0
          else
            density = 2.0d0*rho(ispin)
            if (mode .eq. 5) then
              densityp = 2.0d0*rhop(ispin)
              densitypp = 2.0d0*rhopp(ispin)
              fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)

! s = abs(grad d)/(2kf*d)
! u = (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!      >>  grad(abs(grad d) has mixed derivatives ! <<
! v = (laplacian d)/(d*(2*kf)**2)
              s = abs(densityp)/(fermik*density)
              u = abs(densityp)*densitypp/(density*density*fermik**3)
              v = (densitypp + 2.0d0*densityp/r)/(density*fermik*fermik)
              call xpbe (density, s, u, v, ex, vx)
            else
              stop ' ggaxrad_1c : mode improper'
            end if
            xpot(ispin) = vx
            xen = xen + rho(ispin)*ex
          end if
        end do

! energy
        xen = xen/max(rho(1) + rho(2), epsilon)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine ggaxrad_1c


! ===========================================================================
! xpbe
! ===========================================================================
! Program Description
! ===========================================================================
!  PBE exchange for a spin-unpolarized electronic system
!
!  Gradient-corrected exchange energy based on
!   Perdew-Burke-Ernzerhof GGA
!   as parameterized by J.P. Perdew, K. Burke, M. Ernzerhof,
!   Phys Rev Lett 77, 3865 (1996)
!   Also:
!   [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
!   [b] J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!                                              {\bf 40},  3399  (1989) (E).
!  Hartree a.u.
!
!  input
!  d            density
!  s            abs(grad d)/(2kf*d)
!  u            (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!           >>  grad(abs(grad d) has mixed derivatives ! <<
!  v            (laplacian d)/(d*(2*kf)**2)
!
!  output
!  ex           exchange energy per electron
!  vx           exchange potential

!  Formulas:
!    e_x[unif] = ax*rho^(4/3)  [LDA]
!    ax = -0.75*(3/pi)^(1/3)
!    e_x[PBE] = e_x[unif]*FxPBE(s)
!    FxPBE(s) = 1 + uk - uk/(1+ul*s*s)
!    uk, ul defined after Eq. (13) of 1996 paper
!
! Martin Fuchs, FHI der MPG, Berlin, 02-1993
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
        subroutine xpbe (d, s, u, v, ex, vx)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: d                 ! density
        real, intent (in) :: s                 ! abs(grad d)/(2kf*d)
        real, intent (in) :: u                 ! grad(abs(grad d)...
        real, intent (in) :: v                 ! (laplacian d)/(d*(2*kf)**2)

! Output
        real, intent (out) :: ex               ! exchange energy
        real, intent (out) :: vx               ! exchange potential

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: ax = -0.738558766382022405884230032680836d0
        real, parameter :: um = 0.2195149727645171d0
        real, parameter :: uk = 0.8040d0
        real, parameter :: ul = um/uk

! Local Variable Declaration and Description
! ===========================================================================
        real fac
        real Fs
        real Fss
        real FxPBE
        real p0

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! LDA only
        fac = ax*d**(1.0d0/3.0d0)

! Construct PBE enhancement factor
        p0 = 1.0d0 + ul*s**2
        FxPBE = 1.0d0 + uk - uk/p0
        ex = fac*FxPBE

! find first and second derivatives of Fx w.r.t s.
! Fs = (1/s)*d FxPBE/ ds
! Fss = d(Fs)/ds
        Fs = 2.0d0*uk*ul/(p0**2)
        Fss = -4.0d0*ul*s*Fs/p0

! potential from [b](24)
        vx = fac*((4.0d0/3.0d0)*FxPBE - (u - (4.0d0/3.0d0)*s**3)*Fss - V*Fs)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

      return
      end subroutine xpbe


! ===========================================================================
! ggacrad_1c
! ===========================================================================
! Program Description
! ===========================================================================
!      This routine calculates the Perdew-Burke-Ernzerhof correlation
! potential and energy density. Spherical symmetry is used.
!
! input
!    mode = 5    GGA-C PBE
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Hong Kong Quantum AI Laboratory
! Room 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, Hong Kong
!
! Office Telephone +852 9865 3077
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine ggacrad_1c (mode, rin, rho, rhop, rhopp, cpot, cen)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp

! Output
        real, intent (out) :: cen

        real, intent (out), dimension (2) :: cpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: crs = 1.91915829267751281d0

        real, parameter :: t13 = 1.0d0/3.0d0
        real, parameter :: t23 = 2.0d0/3.0d0

! Local Variable Declaration and Description
! ===========================================================================
        real density, densityp, densitypp
        real fermik
        real g
        real gks2, gks2sq
        real h
        real r
        real rs
        real sk
        real t
        real uu, vv, ww
        real zeta, zetap

        ! energies and potentials
        real ec, vcup, vcdn, dvcup, dvcdn

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! If r is really small, then set to manageably small number.
        r = rin
        if (rin .lt. 1.0d-4) r = 1.0d-4

! LSDA
        density = rho(1) + rho(2)
        densityp = rhop(1) + rhop(2)
        densitypp = rhopp(1) + rhopp(2)
        cen = 0.0d0
        if (density .le. epsilon) then
         cen = 0.0d0
         cpot(1) = 0.0d0
         cpot(2) = 0.0d0
         return
        end if

        if (mode .eq. 5) then

! Correlation terms
          zeta = (rho(1) - rho(2))/density
          zetap = (rhop(1) - rhop(2) - zeta*densityp)/density
          fermik = (3.0d0*pi*pi*density)**t13
          rs = crs/fermik
          sk = 2.0d0*dsqrt(fermik/pi)
          g = ((1.0d0 + zeta)**t23 + (1.d0 - zeta)**t23)/2.0d0
          t = abs(densityp)/(density*(2.0d0*sk*g))
          uu = abs(densityp)*densitypp/(density*density*(2.0d0*sk*g)**3)
          vv = (densitypp + 2.0d0*densityp/r)/(density*(2.0d0*sk*g)**2)
          ww = densityp*zetap/(density*(2.0d0*sk*g)**2)

          call corpbe (.true., rs, zeta, t, uu, vv, ww, ec, vcup, vcdn, h, dvcup, dvcdn)
        else
          stop ' ggacrad_1c : mode improper'
        end if

        cen = ec + h

        cpot(1) = vcup + dvcup
        cpot(2) = vcdn + dvcdn

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

      return
      end subroutine ggacrad_1c


! ===========================================================================
! corpbe
! ===========================================================================
! Program Description
! ===========================================================================
!       Official PBE correlation code. K. Burke, May 14, 1996.
!
! Input: rs = SEITZ RADIUS=(3/4pi rho)^(1/3)
!      : zeta = RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!      : t = ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!      : uu = (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
!      : vv = (LAPLACIAN rho)/(rho * (2*KS*G)**2)
!      : ww = (GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
!      : uu, vv, ww, only needed for PBE potential
!      : lgga=flag to do gga (0=>LSD only)
!      : lpot=flag to do potential (0=>energy only)
!
! output: ec = lsd correlation energy from [a]
!      : vcup = lsd up correlation potential
!      : vcdn = lsd dn correlation potential
!      : h = nonlocal part of correlation energy per electron
!      : dvcup = nonlocal correction to vcup
!      : dvcdn = nonlocal correction to vcdn
!
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl real*8-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!
! Hartree A.U.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine corpbe (tpot, rs, zeta, t, uu, vv, ww, ec, vcm0,           &
     &                     vcp0, h, dvcm0, dvcp0)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: rs
        real, intent (in) :: t
        real, intent (in) :: uu, vv, ww
        real, intent (in) :: zeta

        logical, intent (in) :: tpot

! Output
        real, intent (out) :: ec, h
        real, intent (out) :: vcm0, dvcm0
        real, intent (out) :: vcp0, dvcp0

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: beta = 0.06672455060314922d0

        real, parameter :: t13 = 1.0d0/3.0d0
        real, parameter :: t23 = 2.0d0/3.0d0
        real, parameter :: t43 = 4.0d0/3.0d0
        real, parameter :: t89 = 8.0d0/9.0d0

        real, parameter :: tnum11 =  0.0310907D0
        real, parameter :: tnum12 =  0.21370D0
        real, parameter :: tnum13 =  7.5957D0
        real, parameter :: tnum14 =  3.5876D0
        real, parameter :: tnum15 =  1.6382D0
        real, parameter :: tnum16 =  0.49294D0
        real, parameter :: tnum21 =  0.01554535D0
        real, parameter :: tnum22 =  0.20548D0
        real, parameter :: tnum23 = 14.1189D0
        real, parameter :: tnum24 =  6.1977D0
        real, parameter :: tnum25 =  3.3662D0
        real, parameter :: tnum26 =  0.62517D0
        real, parameter :: tnum31 =  0.0168869D0
        real, parameter :: tnum32 =  0.11125D0
        real, parameter :: tnum33 = 10.357D0
        real, parameter :: tnum34 =  3.6231D0
        real, parameter :: tnum35 =  0.88026D0
        real, parameter :: tnum36 =  0.49671D0

! Local Variable Declaration and Description
! ===========================================================================
        real alfm, alfrsm, alfc
        real b, bg, bec
        real comm
        real delta
        real ecrs, eczeta, eu, eurs, ep, eprs
        real factor, factor0, factor1, factor2, factor3, factor5
        real f, fz
        real fzz
        real g, gz
        real gam
        real hb, ht, hbt, htt, hz, hzt, hrs, hrst
        real pon
        real prefactor
        real q4, q5, q8, q9
        real rtrs
        real tm16
        real z4

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize some parameters
        g = 2.0d0**(4.0d0/3.0d0) - 2.0d0
        fzz = t89/g
        gam = (1.0d0 - log(2.0d0))/pi**2
        delta = beta/gam
        tm16 = -t13/2.0d0

! Find LSD energy contributions, using [c](10) and Table I[c].
! eu = unpolarized LSD correlation energy
! eurs = d(eu)/d(rs)
! ep = fully polarized LSD correlation energy
! eprs = d(ep)/d(rs)
! alfm = -spin stiffness, [c](3).
! alfrsm = - d(alpha)/d(rs)
! f = spin-scaling factor from [c](9).
! construct ec, using [c](8)
        rtrs = sqrt(rs)
        call gcor2 (tnum11, tnum12, tnum13, tnum14, tnum15, tnum16,          &
     &              rtrs, eu, eurs)
        call gcor2 (tnum21, tnum22, tnum23, tnum24, tnum25, tnum26,          &
     &              rtrs, ep, eprs)
        call gcor2 (tnum31, tnum32, tnum33, tnum34, tnum35, tnum36,          &
     &              rtrs, alfm, alfrsm)
        alfc = -alfm
        z4 = zeta**4
        f = ((1.0d0 + zeta)**t43 + (1.0d0 - zeta)**t43 - 2.0d0)/g
        ec = eu*(1.0d0 - f*z4) + ep*f*z4 - alfm*f*(1.0d0 - z4)/fzz

! LSD potential from [c](A1)
! ecrs = d(ec)/d(rs) [c](A2)
! eczeta = d(ec)/d(zeta) [c](A3)
! fz = d(f)/d(zeta) [c](A4)
        ecrs = eurs*(1.0d0 - f*z4) + eprs*f*z4 - alfrsm*f*(1.0d0 - z4)/fzz
        fz = t43*((1.0d0 + zeta)**t13 - (1.0d0 - zeta)**t13)/g
        eczeta = 4.0d0*(zeta**3)*f*(ep - eu + alfm/fzz)                      &
     &          + fz*(z4*ep - z4*eu - (1.0d0 - z4)*alfm/fzz)
        comm = ec - rs*ecrs/3.0d0 - zeta*eczeta

        vcm0 = comm + eczeta
        vcp0 = comm - eczeta

! PBE correlation energy - add to the LSD
! G = phi(zeta), given after [a](3)
! B = A of [a](8)
        g = ((1.0d0 + zeta)**t23 +( 1.0d0 - zeta)**t23)/2.d0
        pon = -ec/(g**3*gam)
        b = delta/(EXP(pon) - 1.0d0)

        q4 = 1.D0 + b*t**2
        q5 = 1.D0 + b*t**2 + b**2*t**4
        h = g**3*(beta/delta)*log(1.0d0 + delta*q4*t**2/q5)

! Calculate potential using appendix E of [b].
        if (tpot) then
          gz = (((1.0d0 + zeta)**2 + epsilon)**tm16 - ((1.0d0 - zeta)**2 + epsilon)**tm16)/3.0d0
          factor = delta/b + 1.0d0
          bg = - 3.0d0*b**2*ec*factor/(beta*g**4)
          bec = b**2*factor/(beta*g**3)
          q8 = q5**2 + delta*q4*q5*t**2
          q9 = 1.0d0 + 2.0d0*b*t**2

          hb = - beta*g**3*b*t**6*(2.0d0 + b*t**2)/q8
          hrs = -(rs*t13)*hb*bec*ecrs
          factor0 = 2.0d0*delta - 6.0d0*b
          factor1 = q5*q9+q4*q9*q9
          hbt = 2.0d0*beta*g**3*t**4*((q4*q5*factor0 - delta*factor1)/q8)/q8
          hrst = (rs*t13)*t**2*hbt*bec*ecrs
          hz = 3.0d0*gz*h/g + hb*(bg*gz + bec*eczeta)
          ht = 2.0d0*beta*g**3*q9/q8
          hzt = 3.0d0*gz*ht/g + hbt*(bg*gz + bec*eczeta)
          factor2 = q4*q5 + b*t**2*(q4*q9 + q5)
          factor3 = 2.0d0*b*q5*q9 + delta*factor2
          htt = 4.0d0*beta*g**3*t*(2.0d0*b/q8 - (q9*factor3/q8)/q8)
          comm = h + hrs + hrst + t**2*ht/6.0d0 + 7.0d0*t**3*htt/6.0d0

          prefactor = hz - gz*t**2*ht/g
          factor5 = gz*(2.0d0*ht + t*htt)/g
          comm = comm - prefactor*zeta - uu*htt - vv*ht - ww*(hzt - factor5)

          dvcm0 = comm + prefactor
          dvcp0 = comm - prefactor
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine corpbe



! ===========================================================================
! gcor2
! ===========================================================================
! Program Description
! ===========================================================================
! Slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by Eq. (10) of
!
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
!
! Hartree A.U.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine gcor2 (a, a1, b1, b2, b3, b4, rtrs, gg, ggrs)
        implicit none


! Argument Declaration and Description
! ===========================================================================
! Input
        real a
        real a1
        real b1
        real b2
        real b3
        real b4
        real rtrs

 ! Output
        real gg
        real ggrs

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real q0
        real q1
        real q2
        real q3

! Procedure
! ===========================================================================
        q0 = -2.0d0*a*(1.0d0 + a1*rtrs**2)
        q1 = 2.0d0*a*rtrs*(b1 + rtrs*(b2 + rtrs*(b3 + b4*rtrs)))
        q2 = log(1.0d0 + 1.0d0/q1)
        gg = q0*q2
        q3 = a*(b1/rtrs + 2.0d0*b2 + rtrs*(3.0d0*b3 + 4.0d0*b4*rtrs))
        ggrs = -2.0d0*a*a1*q2 - q0*q3/(q1*(1.0d0 + q1))

! Format Statements
! ===========================================================================
        return
        end subroutine gcor2


! End Module
! =============================================================================
        end module
