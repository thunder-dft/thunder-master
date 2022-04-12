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

! M_xc
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
! Eduardo Mendez
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444

! rewritten by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
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

        real aln
        real dec
        real dex
        real drvexc
        real ecp
        real ex
        real exc
        real fx
        real fxc
        real rh
        real rs
        real x
        real zeta

        real, dimension (2) :: cpot
        real, dimension (2) :: d
        real, dimension (2) :: dp
        real, dimension (2) :: dpp
        real, dimension (2) :: xpot

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

! Rho must be positive, but not too small
        if (rho .lt. 1.0d-8) then
          rho = 0.0d0
          return
        else if (rho .lt. 1.0d-5) then
          rho = 1.0d-5
        end if

! Determine exchange-correlation potentials

! C Lee-Yang-Parr
		if (iexc .eq. 9) then

! X Becke gga by default
          ix = 2

! X Perdew-Wang gga
          d = 0.5d0*rho
          dp = 0.5d0*rhop
          dpp = 0.5d0*rhopp
          call ggaxrad_1c (ix, r, d, dp, dpp, xpot, dex)
          call ggacrad_1c (4, r, d, dp, dpp, cpot, dec)
          vpxc = xpot(1) + cpot(1)

! If the improper iexc option was entered then the program will stop.
        else
          write (*,*) ' In get_potxc_1c - '
          write (*,*) ' stop: xc option not implemented', iexc
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
!
!      This routine calculates the exchange potential and energy density.
! Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-X Becke
!    mode = 3    GGA-X Perdew
!    mode = 5    GGA-X Burke-Perdew-Ernzerhof
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

! rewritten by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
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

        real density
        real densityp
        real densitypp
        real ex
        real fermik
        real r
        real s
        real u
        real v
        real vx

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize

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
		    if (mode .eq. 2) then
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
              call xbecke (density, s, u, v, ex, vx)
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
! xbecke
! ===========================================================================
! Program Description
! ===========================================================================
!  Becke exchange for a spin-unpolarized electronic system
!
!  Gradient-corrected exchange energy based on
!     [A.D. Becke, J.Chem.Phys.96, 2155, 1992].
!  The LSDA energy functional, obtained as E{n+,n-}=(E{2n+}+E{2n-})/2,
!     and the functional derivative formula are given by
!     [J.P. Perdew , PRB 33, 8800, 1986].
!     [J.P. Perdew , PRB 34, 7406, 1986].
!  see also [G.Ortiz ...,PRB 43, 6376 (1991)] eq. (A2)
!
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
!
! Martin Fuchs, FHI der MPG, Berlin, 02-1993
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

! rewritten by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
      subroutine xbecke (d, s, u, v, ex, vx)
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
        real, parameter :: c = 0.779555417944150792d1
        real, parameter :: b = 0.42d-2
        real, parameter :: bb = -.451357747124625192d-2
        real, parameter :: ax = -.738558766382022406d0

! Local Variable Declaration and Description
! ===========================================================================
        real dd1
        real ddi
        real f
        real fac
        real fs
        real fss
        real g
        real g1
        real x
        real y0
        real y1
        real y2

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! exchange enhancement factor f
        x  = c*s
        y1 = 1.d0/sqrt(1.d0 + x*x)
        y0 = log(x + 1.d0/y1)
        y2 = -x*y1*y1*y1
        ddi= 1.d0/(1.d0 + 6.d0*b*x*y0)
        dd1= 6.d0*b*(y0 + x*y1)
        g  = 1.d0 - 0.5d0*x*dd1*ddi
        fs = -2.d0*bb*c*c*ddi
        g1 = -3.d0*b*(y0 + x*(3.d0*y1 + x*y2 - dd1*dd1*ddi/(6.0d0*b)))
        fss= fs*c*(g1 - g*dd1)*ddi
        fs = fs*g
        f  = 1.d0 - bb*x*x*ddi

! LDA only
        fac= ax*d**(1.0d0/3.0d0)

! energy
        ex = fac*f

! potential
        vx = fac*((4.0d0/3.0d0)*f - (u - (4.0d0/3.0d0)*s*s*s)*fss - v*fs)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

      return
      end subroutine xbecke


! ===========================================================================
! ggacrad_1c
! ===========================================================================
! Program Description
! ===========================================================================
!
!      This routine calculates the correlation potential and energy density.
! Spherical symmetry is used. LSDA - GGA
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-X Becke
!    mode = 3    GGA-X Perdew
!    mode = 5    GGA-X Burke-Perdew-Ernzerhof
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

! rewritten by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
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
        real, parameter :: thrd = 0.333333333333333333d0

! Local Variable Declaration and Description
! ===========================================================================
        real alfc
        real density
        real densityp
        real densityp11
        real densityp12
        real densityp22
        real densitypp
        real ec
        real ecrs
        real eczet
        real fermik
        real g
        real gsfermik
        real h
        real r
        real rs
        real sfermik
        real t
        real uu
        real vv
        real ww
        real zet
        real ztp
        real fk
        real sk

        real, dimension (2) :: dvc, vc
        real, dimension (2) :: flip

! Allocate Arrays
! ===========================================================================

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

        if (mode .eq. 4) then
         flip = rhopp + 2.0d0*rhop/r
         call corlyp_1c (.true., rho(1), rho(2), rhop(1), rhop(2), flip(1),   &
     &                   flip(2), cen, cpot(1), cpot(2))
        else
         stop 'ggacrad1c : mode improper'
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine ggacrad_1c

! ===========================================================================
! corlyp_1c
! ===========================================================================
! Program Description
! ===========================================================================
!       Lee Yang Parr correlation energy functional one-dimensional densities
! only no provisions taken against division by zero.
!
! See e.g. C. Lee et al. Phys. Rev. B 37 785 (1988)
!
! Hartree A.U.
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

! rewritten by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine corlyp_1c (tpot, dp, dm, dp1, dm1, dp2, dm2, ec, vcp0, vcm0)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: dm, dm1, dm2
        real, intent (in) :: dp, dp1, dp2

        logical, intent (in) :: tpot

! Output
        real, intent (out) :: ec
        real, intent (out) :: vcm0
        real, intent (out) :: vcp0

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: aa = 0.04918d0
        real, parameter :: bb = 0.132d0
        real, parameter :: cc = 0.2533d0
        real, parameter :: dd = 0.349d0
        real, parameter :: c5 = 4.55779986d0
        real, parameter :: c6 = 1.0d0/72.0d0
        real, parameter :: c7 = 1.0d0/18.0d0
        real, parameter :: c8 = 0.125d0
        real, parameter :: t13 = 1.0d0/3.0d0
        real, parameter :: t89 = 8.0d0/9.0d0

! Local Variable Declaration and Description
! ===========================================================================
        real c1, c2, c3, c4, c9
        real chf
        real d0xt13, d0xt53
        real d0, d1, d2
        real dmt53, dpt53
        real dxsq
        real ga
        real gafm, gafp
        real gb
        real h
        real h2
        real hf
        real hff
        real sc
        real sc2
        real scf
        real t43, t53, t83
        real yafm, yafp
        real yb, yb1, yb2
        real ybfm, ybfp
        real yy1
        real yz, yz1, yz2
        real z1, z2
        real zfm, zfp
        real zeta

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize some parameters
        c1 = -4.0d0*aa
        c2 = dd
        c3 = 2.0d0*bb
        c4 = cc
        t53 = 5.0d0*t13
        t43 = 4.0d0*t13
        t83 = 2.0d0*t43
        c9 = t43 + t89

        d0 = dp + dm
        dxsq = 1.0d0/(d0*d0)
        d1 = dp1 + dm1
        d2 = dp2 + dm2
        d0xt13 = d0**(-t13)
        d0xt53 = d0xt13*d0xt13/d0
        dpt53 = dp**t53
        dmt53 = dm**t53

! Polarization factor
        zeta = c1*(dp*dm)*dxsq

! Scaling function
        sc = 1.0d0/(1.0d0 + c2*d0xt13)
        h = c3*d0xt53*exp(-c4*d0xt13)

! kinetic energy density expansion
        ga = c5*(dp*dpt53 + dm*dmt53)
        gb = c6*(dp1*dp1 - dp*dp2 + dm1*dm1 - dm*dm2) + c7*(dp*dp2 + dm*dm2) &
           + c8*(d0*d2 - d1*d1)

! Calculate potential
        if (tpot) then
         gafp = t83*c5*dpt53
         gafm = t83*c5*dmt53

         scf = t13*c2*d0xt13/d0*sc*sc
         sc2 = scf*(d2 + 2.0d0*(scf/sc - 2.0d0*t13/d0)*d1*d1)

         chf = t13*(c4*d0xt13 - 5.0d0)/d0
         hf = chf*h
         hff = h*(chf**2 + t13*(5.0d0 - 4.0d0*t13*c4*d0xt13)*dxsq)
         h2 = (hf*d2 + hff*d1*d1)

         zfp = (c1*dm - 2.0d0*zeta*d0)*dxsq
         zfm = (c1*dp - 2.0d0*zeta*d0)*dxsq
         yz = zeta/c1
         yy1 = dp*dm1 + dm*dp1
         yz1 = (yy1 - 2.0d0*yz*d1*d0)*dxsq
         yz2 = (2.0d0*yz*d1*d1 - 2.0d0*(yz1*d1 + yz*d2)*d0                   &
             - 2.0d0*d1*yy1/d0 + (dp*dm2 + 2.0d0*dp1*dm1 + dm*dp2))*dxsq
         z1 = c1*yz1
         z2 = c1*yz2

         yafp = sc*(d0*zfp + zeta) + zeta*d0*scf
         yafm = sc*(d0*zfm + zeta) + zeta*d0*scf

         yb = sc*zeta*h
         ybfp = sc*(h*zfp + zeta*hf) + zeta*h*scf
         ybfm = sc*(h*zfm + zeta*hf) + zeta*h*scf
         yb1 = sc*(h*z1 + zeta*hf*d1) + zeta*h*scf*d1
         yb2 = (sc*hf + h*scf)*d1*z1 + h*sc*z2 + (sc*z1 + zeta*scf*d1)*hf*d1 &
             + zeta*sc*h2 + (zeta*hf*d1 + h*z1)*scf*d1 + zeta*h*sc2

! Collect contributions
         vcp0 = yafp + ybfp*(ga + gb)                                        &
              + yb*(gafp + 2.0d0*c8*(c9*dp2 + 2.0d0*dm2))                    &
              + yb1*2.0d0*c8*(c9*dp1 + 2.0d0*dm1) + yb2*c8*(t43*dp + dm)
         vcm0 = yafm + ybfm*(ga + gb)                                        &
              + yb*(gafm + 2.0d0*c8*(c9*dm2 + 2.0d0*dp2))                    &
              + yb1*2.0d0*c8*(c9*dm1 + 2.0d0*dp1) + yb2*c8*(t43*dm + dp)
        else
         vcp0 = 0.0d0
         vcm0 = 0.0d0
        endif

! Correlation energy per electron
        ec = zeta*sc*(d0 + h*(ga + gb))/d0

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine corlyp_1c

! End Module
! =============================================================================
        end module
