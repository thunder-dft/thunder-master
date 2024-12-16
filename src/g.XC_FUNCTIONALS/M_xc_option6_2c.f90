! copyright info:
!
!                             @Copyright 2010
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
        module M_xc_2c

! /GLOBAL
        use M_precision

! /XC_FUNCTIONALS
        use M_xc_1c

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! get_potxc_2c
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
        subroutine get_potxc_2c (iexc, xc_fraction, r, rho, rhop, rhopp,     &
     &                           rhoz, rhozz, rhopz, newexc, vpxc, dnuxc,    &
     &                           dnuxcs)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc

        real, intent(in) :: xc_fraction
        real, intent(inout) :: r        ! radial distance

        ! density and derivatives
        real, intent(inout) :: rho      ! return zero value if rho is small
        real, intent(in) :: rhop
        real, intent(in) :: rhopp
        real, intent(in) :: rhoz
        real, intent(in) :: rhozz
        real, intent(in) :: rhopz

! Output
        real, intent(out) :: dnuxc
        real, intent(out) :: dnuxcs
        real, intent(out) :: newexc
        real, intent(out) :: vpxc
        
! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ix

        real ex                      ! exchange potential
        real dec, dex                ! derivative of correlation and exchange
        real drvexc
        real fx

! density and derivatives - spin cases
        real, dimension (2) :: d, dp, dpp, dz, dzz, dpz

! exchange and correlation potentials - for spin
        real, dimension (2) :: xpot
        real, dimension (2) :: cpot
                
! Procedure
! ===========================================================================
! Initialize variables (only used in 3,11)
        dnuxc = 0.0d0
        dnuxcs = 0.0d0

! If r is really small, then set to manageably small number.
        if (r .lt. 1.0d-4) r = 1.0d-4

! Rho must be positive, but not too small
        if (rho .lt. 1.0d-8) then
          rho = 0.0d0
          dnuxc = 0.0d0
          newexc = 0.0d0
          vpxc = 0.0d0
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

          d = 0.5d0*rho
          dp = 0.5d0*rhop
          dpp = 0.5d0*rhopp
          dz = 0.5d0*rhoz
          dzz = 0.5d0*rhozz
          dpz = 0.5d0*rhopz
          call ggaxrad_2c (ix, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
          call ggacrad_2c (5 , r, d, dp, dpp, dz, dzz, dpz, cpot, dec)

          vpxc = xpot(1) + cpot(1)

! If the improper iexc option was entered then the program will stop.
        else
          write (*,*) ' In get_potxc_2c.f90 - '
          write (*,*) ' You must recomplile create.x for iexc = 6 '
          write (*,*) ' Set XC = PBE in include/OPTIONS and recompile. '
          stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
        newexc = dex + dec

! Format Statements
! ===========================================================================
! None

        return
        end subroutine get_potxc_2c


! ===========================================================================
! ggaxrad_2c
! ===========================================================================
! Program Description
! ===========================================================================
!      This routine calculates the Becke exchange potential and energy density.
! Spherical symmetry is used.
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992!
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
        subroutine ggaxrad_2c (mode, rin, rho, rhop, rhopp, rhoz, rhozz,     &
     &                         rhopz, xpot, xen)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp
        real, intent (in), dimension (2) :: rhoz
        real, intent (in), dimension (2) :: rhozz
        real, intent (in), dimension (2) :: rhopz

! Output
        real, intent (out) :: xen

        real, intent (out), dimension (2) :: xpot

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ispin

! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp
        real density_z, density_zz
        real density_pz

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
            if  (mode .eq. 5) then
              density_p = 2.0d0*rhop(ispin)
              density_pp = 2.0d0*rhopp(ispin)
              density_z = 2.0d0*rhoz(ispin)
              density_zz = 2.0d0*rhozz(ispin)
              density_pz = 2.0d0*rhopz(ispin)
              fermik = 2.0d0*(3.0d0*pi*pi*density)**(1.0d0/3.0d0)

! s = abs(grad d)/(2kf*d)
! u = (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!      >>  grad(abs(grad d) has mixed derivatives ! <<
! v = (laplacian d)/(d*(2*kf)**2)
              s = sqrt(density_p**2 + density_z**2)/(fermik*density)
              u = (density_p**2*density_pp                                   &
     &             + 2.0d0*density_p*density_z*density_pz &
     &             + density_z**2*density_zz)/(s*density**3*fermik**4)
              v = (density_pp + density_p/r + density_zz)/(density*fermik*fermik)

              call xpbe (density, s, u, v, ex, vx)
            else
              stop 'ggaxrad_2c : mode improper'
            end if
            xpot(ispin) = vx
            xen = xen + rho(ispin)*ex
          end if
        end do

! energy
        xen = xen/max(rho(1) + rho(2), epsilon)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine ggaxrad_2c


! ===========================================================================
! ggacrad_2c
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
        subroutine ggacrad_2c (mode, rin, rho, rhop, rhopp, rhoz, rhozz,     &
     &                         rhopz, cpot, cen)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mode

        real, intent (in) :: rin

        real, intent (in), dimension (2) :: rho
        real, intent (in), dimension (2) :: rhop
        real, intent (in), dimension (2) :: rhopp
        real, intent (in), dimension (2) :: rhoz
        real, intent (in), dimension (2) :: rhozz
        real, intent (in), dimension (2) :: rhopz

! Output
        real, intent (out) :: cen

        real, intent (out), dimension (2) :: cpot

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real density
        real r

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
        cen = 0.0d0

        if (density .le. epsilon) then
          cen = 0.0d0
          cpot(1) = 0.0d0
          cpot(2) = 0.0d0
        else
          if (mode .eq. 5) then
            call corpbe_2c (.true., r, rho(1), rho(2), rhop(1), rhop(2),     &
     &                      rhopp(1), rhopp(2), rhoz(1), rhoz(2), rhozz(1),  &
                            rhozz(2), rhopz(1), rhopz(2), cen, cpot(1), cpot(2))
          else
            stop 'ggacrad_2c : mode improper'
          end if
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine ggacrad_2c


! ===========================================================================
! corpbe_2c.f90
! ===========================================================================
! Program Description
! ===========================================================================
! This routine is a driver for the PBE subroutines, using simple inputs
! K. Burke, May 14, 1996.
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
!
! Input Variables
!
! tpot .... T  evaluate correlation energy and potential
!           F  evaluate energy only (a posteriori scheme)
! rin ....... position
! pa ...... spin up density
! pb ...... spin down density
! dpaofr .. 1st partial derivative of spin up with respect to r
! dpaofz .. 1st partial derivative of spin up with respect to z
! d2paofr . 2nd partial derivative of spin up with respect to r
! d2paofz . 2nd partial derivative of spin up with respect to z
! dpbofr .. 1st partial derivative of spin down with respect to r
! dpbofz .. 1st partial derivative of spin down with respect to z
! d2pbofr . 2nd partial derivative of spin down with respect to r
! d2pbofz . 2nd partial derivative of spin down with respect to z
!
! Output Variables
!
! ec ...... correlation energy per electron
! vp ...... correlation potential
! ===========================================================================
! Code written by:
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
        subroutine corpbe_2c (tpot, r, pa, pb, dpaofr, dpbofr, d2paofr,      &
     &                        d2pbofr, dpaofz, dpbofz, d2paofz, d2pbofz,     &
     &                        d2paofrz, d2pbofrz, ecpbe, vcuppbe, vcdnpbe)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        logical, intent (in) :: tpot

        real, intent (in) :: r

        real, intent (in) :: pa, pb
        real, intent (in) :: dpaofr, dpbofr
        real, intent (in) :: d2paofr, d2pbofr
        real, intent (in) :: dpaofz, dpbofz
        real, intent (in) :: d2paofz, d2pbofz
        real, intent (in) :: d2paofrz, d2pbofrz

! Output
        real, intent (out) :: ecpbe
        real, intent (out) :: vcuppbe, vcdnpbe

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: crs = 1.91915829267751281d0

        real, parameter :: t13 = 1.0d0/3.0d0
        real, parameter :: t23 = 2.0d0/3.0d0
        real, parameter :: t43 = 4.0d0/3.0d0
        real, parameter :: t89 = 8.0d0/9.0d0

! Local Variable Declaration and Description
! ===========================================================================
! density and derivatives
        real p             ! Sum of spin up and spin down densities
        real dpr           ! Sum of spin up and down partials wrt r
        real dpz           ! Sum of spin up and down partials wrt z
        real d2pr          ! Sum of spin 2nd partials wrt r
        real d2pz          ! Sum of spin 2nd partials wrt z
        real d2prz         ! Sum of spin 2nd partials wrt r and z

! gradient terms
        real agrad, agrup, agrdn
        real lap, uplap, dnlap
        real delgrad, delgrup, delgrdn

        real fermik
        real g
        real rs
        real sk
        real t
        real uu, vv, ww
        real zeta, zetap

! results from corpbe
        real ec, vcup, vcdn
        real h, dvcup, dvcdn

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize the inputs - take out the spin factor of 2.0d0
!       pa = 0.5d0*pa; pb = 0.5d0*pb
!       dpaofr = 0.5d0*dpaofr; dpbofr = 0.5d0*dpbofr
!       d2paofr = 0.5d0*d2paofr; d2pbofr = 0.5d0*d2pbofr
!       d2paofz = 0.5d0*d2paofz; d2pbofz = 0.5d0*d2pbofz
!       d2paofrz = 0.5d0*d2paofrz; d2pbofrz = 0.5d0*d2pbofrz

! Here the option to caldpaofrculate the potential is evaluated
        if (tpot) then
          ! total density
          p = (pa + pb)/2.0d0

          dpr = (dpaofr + dpbofr)/2.0d0
          dpz = (dpaofz + dpbofz)/2.0d0

          d2pr = (d2paofr + d2pbofr)/2.0d0
          d2pz = (d2paofz + d2pbofz)/2.0d0
          d2prz = (d2paofrz + d2pbofrz)/2.0d0

          ! agrup = |grad up|, agrdn = |grad dn|
          agrad = sqrt(dpr**2 + dpz**2)
          agrup = sqrt(dpaofr**2 + dpaofz**2)/4.0d0
          agrdn = sqrt(dpbofr**2 + dpbofz**2)/4.0d0

          uplap = d2paofr + dpaofr/r + d2paofz
          dnlap = d2pbofr + dpbofr/r + d2pbofz
          lap = (uplap + dnlap)/2.0d0

          ! might need later for spin polarization
!         delgrup = dpaofr*(dpaofr*d2paofr + dpaofz*d2paofrz)/(8.0d0*agrup)  &
!    &             + dpaofz*(dpaofr*d2paofrz + dpaofz*d2paofz)/(8.0d0*agrup)
!         delgrdn = dpbofr*(dpbofr*d2pbofz + dpbofz*d2pbofrz)/(8.0d0*agrdn)  &
!    &             + dpbofz*(dpbofr*d2pbofrz + dpbofz*d2pbofz)/(8.0d0*agrdn)
          delgrad = dpr*(dpr*d2pr + dpz*d2prz)/agrad                         &
     &             + dpz*(dpr*d2prz + dpz*d2pz)/agrad

! Input: rs = SEITZ RADIUS=(3/4pi rho)^(1/3)
!      : zeta = RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!      : t = ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!      : uu = (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
!      : vv = (LAPLACIAN rho)/(rho * (2*KS*G)**2)
!      : ww = (GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
!      : uu, vv, ww, only needed for PBE potential
! Correlation terms
          zeta = (pa - pb)/(2.0d0*p)
          zetap = (((dpaofr - dpbofr)/2.0d0) - zeta*dpr)/p
          fermik = (3.0d0*pi*pi*p)**t13
          rs = crs/fermik
          sk = 2.0d0*sqrt(fermik/pi)
          g = ((1.0d0 + zeta)**t23 + (1.d0 - zeta)**t23)/2.0d0

          t = agrad/((2.0d0*sk*g)*p)
          uu = delgrad/(p**2*(2.0d0*sk*g)**3)
          vv = lap/(p*(2.0d0*sk*g)**2)
          ww = (agrup**2 - agrdn**2 - zeta*agrad**2)/(p**2*(2.0d0*sk*g)**2)

          call corpbe (tpot, rs, zeta, t, uu, vv, ww, ec, vcup, vcdn, h, dvcup, dvcdn)

          ecpbe = ec + h
          vcuppbe = vcup + dvcup
          vcdnpbe = vcdn + dvcdn
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine corpbe_2c


! End Module
! =============================================================================
        end module
