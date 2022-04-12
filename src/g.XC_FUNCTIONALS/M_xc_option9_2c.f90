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
          use M_xc_1c
          use M_precision

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

!! Argument Declaration and Description
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

        real aln
        real exc                     ! exchange-correlation energy
        real ecp                     ! correlation potential
        real ex                      ! exchange potential
        real dec, dex                ! derivative of correlation and exchange
        real drvexc
        real fx
        real fxc
        real rs
        real x
        real zeta

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

! C Lee-Yang-Parr
        if (iexc .eq. 9) then

! X Becke gga by default
          ix = 2

          d = 0.5d0*rho
          dp = 0.5d0*rhop
          dpp = 0.5d0*rhopp
          dz = 0.5d0*rhoz
          dzz = 0.5d0*rhozz
          dpz = 0.5d0*rhopz
          call ggaxrad_2c (ix, r, d, dp, dpp, dz, dzz, dpz, xpot, dex)
          call ggacrad_2c (4 , r, d, dp, dpp, dz, dzz,      cpot, dec)

          vpxc = xpot(1) + cpot(1)


! If the improper iexc option was entered then the program will stop.
        else
          write (*,*) ' In get_potxc_2c.f90 - '
          write (*,*) ' stop: xc option not implemented', iexc
          stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy
        newexc = dec + dex

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
            if  (mode .eq. 2) then
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


                 call xbecke (density, s, u, v, ex, vx)
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
!      This routine calculates the correlation potential and energy density.
! Spherical symmetry is used. LSDA - GGA, Hartree a.u.
!
! input
!    mode = 1    LSDA
!    mode = 2    GGA-C Perdew 91
!    mode = 3    GGA-C Perdew 86
!    mode = 4    GGA-C Lee-Yang-Parr 1988
!    mode = 5    GGA-C Burke-Perdew-Ernzerhof
!    r           radius
!    rh()        spin up/down density
!    rhp()       1st derivative of rh
!    rhpp()      2nd derivative of rh
!
! output
!    zet         spin polarization
!    cpot()       - " -  correlation potential
!    cen         correlation energy density
!
! convention
!    yy(1) = spin up, yy(2) = spin down
!
! Martin Fuchs, FHI der MPG, Berlin, 07-1992
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
        subroutine ggacrad_2c (mode, rin, rho, rhop, rhopp, rhoz, rhozz,     &
     &                         cpot, cen)
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

! Output
        real, intent (out) :: cen

        real, intent (out), dimension (2) :: cpot

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: crs = 1.91915829267751281d0

! Local Variable Declaration and Description
! ===========================================================================
! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp
        real density_p11, density_p12, density_p22
        real density_pz

! inputs
        real alfc, h, rs, t, uu, vv, ww
        real zet, ztp, fk, sk, g

! answers
        real ec, ecrs, eczet, vcdn, vcup, dvcdn, dvcup


! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize

! LSDA
        density = rho(1) + rho(2)
        cen = 0.0d0

        if (density .le. epsilon) then
          cen = 0.0d0
          cpot(1) = 0.0d0
          cpot(2) = 0.0d0
        else
          if (mode .eq. 4) then
            call corlyp_2c (.true., rin, rho(1), rho(2), rhop(1), rhop(2),   &
     &                     rhopp(1), rhopp(2), rhoz(1), rhoz(2), rhozz(1),   &
                           rhozz(2), cen, cpot(1))
          else
            stop 'ggacrad : mode improper'
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


! ========================================================================
! corlyp_2c.f90
! ========================================================================
! Program Description
! ========================================================================
!
!  The subroutine corlyp calculates the correlation energy using
!  the Lee, Yang and Parr generalized gradient approximation.
!
!
! Input Variables
!
! tpot .... T  evaluate correlation energy and potential
!           F  evaluate energy only (a posteriori scheme)
! x ....... dummy
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
!
! ===========================================================================
! Code written by:
! Richard B. Evans and James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)

! Code rewritten by:
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
        subroutine corlyp_2c (tpot, r, pa, pb, dpaofr, dpbofr, d2paofr,      &
     &                        d2pbofr, dpaofz, dpbofz, d2paofz, d2pbofz,     &
     &                        ec, vp)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        logical  tpot

        real r
        real pa
        real pb
        real dpaofr
        real dpaofz
        real dpbofr
        real dpbofz
        real d2paofr
        real d2paofz
        real d2pbofr
        real d2pbofz

! Output
        real ec
        real vp

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: aa = 0.04918d0
        real, parameter :: bb = 0.132d0
        real, parameter :: cc = 0.2533d0
        real, parameter :: dd = 0.349d0
        real, parameter :: t13 = 1.0d0/3.0d0
        real, parameter :: t23 = 2.0d0/3.0d0
        real, parameter :: t43 = 4.0d0/3.0d0
        real, parameter :: t53 = 5.0d0/3.0d0
        real, parameter :: t73 = 7.0d0/3.0d0
        real, parameter :: t83 = 8.0d0/3.0d0

! Local Variable Delcaration and Descrition
! ======================================================================
! constants
        real Cf
        real dep
        real dp_one
        real ep

! density and derivatives
        real p             ! Sum of spin up and spin down densities
        real dpr           ! Sum of spin up and down partials wrt r
        real dpz           ! Sum of spin up and down partials wrt z
        real d2pr          ! Sum of spin 2nd partials wrt r
        real d2pz          ! Sum of spin 2nd partials wrt z

        real gamma
        real dgammap       ! partial derivative of gamma wrt pa
        real dgammar       ! partial derivative of gamma wrt r
        real dgammaz       ! partial derivative of gamma wrt z
        real d2gammaz      ! 2nd Partial derivative of gamma wrt z
        real d2gammar      ! 2nd Partial derivative of gamma wrt r

        real Fofp
        real dFofp         ! partial derivative of Fofp wrt pa
        real dFofpz        ! partial derivative of Fofp wrt z
        real dFofpr        ! partial derivative of Fofp wrt r
        real d2Fofpr       ! 2nd Partial derivative of F wrt r
        real d2Fofpz       ! 2nd Partial derivative of F wrt z

        real Gofp
        real dGofp         ! partial derivative of Gofp wrt pa
        real dGofpz        ! partial of Gofp wrt z
        real dGofpr        ! partial of Gofp wrt r
        real d2Gofpr       ! 2nd partial of Gofp wrt r
        real d2Gofpz       ! 2nd partial of Gofp wrt z
        real d2Gofp        ! Laplacian of Gofp

        real tW            ! Weizsacker kinetic-energy density for p
        real taW           ! Weizsacker K.E. dens. for spin up
        real tbW           ! Weizsacker K.E. dens. for spin down

! Procedure
! ===========================================================================
! The notation used for our variables is similar to the original
! notation used by Lee, Yang and Parr. (Physical Review B 37 785 (1988))
! Initialize pi

        Cf = 0.3d0*(3.0d0*pi**2)**t23

! Here the option to calculate the potential is evaluated
        if (tpot) then
          p = pa + pb

          dpr = dpaofr + dpbofr
          dpz = dpaofz + dpbofz

          d2pr = d2paofr + d2pbofr
          d2pz = d2paofz + d2pbofz

          ep = p**(-t53)*exp(-cc*p**(-t13))
          dep = cc*t13*p**(-t43) - t53/p

          dp_one = 1.0d0 + dd*p**(-t13)

          gamma = 2.0d0*(1.0d0 - (pa**2 + pb**2)/p**2)
          dgammap = - 4.0d0*pa/p**2 - 2.0d0*(gamma - 2.0d0)/p
          dgammar = - 4.0d0*(pa*dpaofr + pb*dpbofr)/p**2                     &
                   - 2.0d0*dpr*(gamma - 2.0d0)/p
          dgammaz = - 4.0d0*(pa*dpaofz + pb*dpbofz)/p**2                     &
                   - 2.0d0*dpz*(gamma - 2.0d0)/p
          d2gammar =                                                         &
           - 4.0d0*(pa*d2paofr + dpaofr**2 + pb*d2pbofr + dpbofr**2)/p**2    &
           + 8.0d0*(pa*dpaofr + pb*dpbofr)*dpr/p**3                          &
           - 2.0d0*d2pr*(gamma - 2.0d0)/p                                    &
           - 2.0d0*dpr*(p*dgammar - (gamma - 2.0d0)*dpr)/p**2
          d2gammaz =                                                         &
           - 4.0d0*(pa*d2paofz + dpaofz**2 + pb*d2pbofz + dpbofz**2)/p**2    &
           + 8.0d0*(pa*dpaofz + pb*dpbofz)*dpz/p**3                          &
           - 2.0d0*d2pz*(gamma - 2.0d0)/p                                    &
           - 2.0d0*dpz*(p*dgammaz - (gamma - 2.0d0)*dpz)/p**2

          Fofp  = gamma/dp_one
          dFofp = (dgammap + dd*t13*p**(-t43)*Fofp)/dp_one
          dFofpr = (dgammar + dd*t13*p**(-t43)*Fofp*dpr)/dp_one
          dFofpz = (dgammaz + dd*t13*p**(-t43)*Fofp*dpz)/dp_one
          d2Fofpr =                                                          &
          (dp_one*(d2gammar - dd*t13*t43*p**(-t73)*Fofp*dpr**2               &
                   + dd*t13*p**(-t43)*(dFofpr*dpr + Fofp*d2pr))              &
           + dd*t13*p**(-t43)*dpr*(dgammar + dd*t13*p**(-t43)*Fofp*dpr))     &
          /dp_one**2
         d2Fofpz =                                                           &
         (dp_one*(d2gammaz - dd*t13*t43*p**(-t73)*Fofp*dpz**2                &
                  + dd*t13*p**(-t43)*(dFofpz*dpz + Fofp*d2pz))               &
           + dd*t13*p**(-t43)*dpz*(dgammaz + dd*t13*p**(-t43)*Fofp*dpz))     &
          /dp_one**2

          Gofp  = Fofp*ep
          dGofp = dFofp*ep + Gofp*dep
          dGofpr = dFofpr*ep + Gofp*dep*dpr
          dGofpz = dFofpz*ep + Gofp*dep*dpz
          d2Gofpr =                                                          &
          (d2Fofpr + dFofpr*dpr*dep)*ep + (dGofpr*dpr + Gofp*d2pr)*dep       &
            + Gofp*dpr**2*(t53/p**2 - cc*t13*t43*p**(-t73))
          d2Gofpz =                                                          &
          (d2Fofpz + dFofpz*dpz*dep)*ep + (dGofpz*dpz + Gofp*d2pz)*dep       &
           + Gofp*dpz**2*(t53/p**2 - cc*t13*t43*p**(-t73))

          d2Gofp = d2Gofpr + dGofpr/r + d2Gofpz

          vp = - aa*(dFofp*p + Fofp)                                         &
               - 2.0d0**t53*aa*bb*Cf                                         &
                      *(dGofp*(pa**t83 + pb**t83) + t83*Gofp*pa**t53)        &
               - aa*bb/4.0d0*(p*d2Gofp + 4.0d0*(dGofpr*dpr + dGofpz*dpz)     &
                        + 4.0d0*Gofp*(d2pr + dpr/r + d2pz)                   &
                        + dGofp*(p*(d2pr + dpr/r + d2pz)                     &
                        - (dpr**2 + dpz**2)))                                &
               - aa*bb/36.0d0                                                &
                   *(3.0d0*pa*d2Gofp + 4.0d0*(dpaofr*dGofpr + dpaofz*dGofpz) &
               + 4.0d0*Gofp*(d2paofr + dpaofr/r + d2paofz)                   &
               + 3.0d0*dGofp*(pa*(d2paofr + dpaofr/r + d2paofz)              &
                              + pb*(d2pbofr + dpbofr/r + d2pbofz))           &
               + dGofp*(dpaofr**2 + dpaofz**2 + dpbofr**2 + dpbofz**2))
        else
          vp = 0.0d0
        end if

! correlation energy per electron
        tW = ((dpr**2 + dpz**2)/p - (d2pr + dpr/r + d2pz))/8.0d0
        taW = ((dpaofr**2 + dpaofz**2)/pa                                    &
              - (d2paofr + dpaofr/r + d2paofz))/8.0d0
        tbW = ((dpbofr**2 + dpbofz**2)/pb                                    &
              - (d2pbofr + dpbofr/r + d2pbofz))/8.0d0

        ec = 2.0d0**t23*Cf*(pa**t83 + pb**t83) - p*tW                        &
           + (pa*taW + pb*tbW)/9.0d0                                         &
           + (pa*(d2paofr + dpaofr/r + d2paofz)                              &
              + pb*(d2pbofr + dpbofr/r + d2pbofz))/18.0d0
        ec = - aa*gamma/dp_one*(p + 2.0d0*bb*p**(-t53)*exp(-cc*p**(-t13))*ec)
        ec = ec/p

! Format Statements
! ===========================================================================
! None

        return
        end subroutine corlyp_2c

! End Module
! =============================================================================
        end module
