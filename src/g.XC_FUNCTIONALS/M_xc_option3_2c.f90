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

! M_xc_2c
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

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! get_potxc_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This program will access the requested exchange and correlation
! functionals from the following list:
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
! Subroutine Declaration
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
        real dummy

        real exc                     ! exchange-correlation energy
        real ex                      ! exchange potential
        real dec, dex                ! derivative of correlation and exchange
        real drvexc
        real fx
        real fxc

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        dummy = xc_fraction
        dummy = rhop
        dummy = rhopp
        dummy = rhopz
        dummy = rhoz
        dummy = rhozz

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

! XC Ceperley - Alder
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
        if (iexc .eq. 3) then
          call ceperley_alder (rho, ex, fx, exc, fxc, drvexc, dnuxc)
          vpxc = fxc
          dex = ex
          dec = exc - ex

! If the improper iexc option was entered then the program will stop.
        else
          write (*,*) ' In get_potxc_2c.f90 - '
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
        end subroutine get_potxc_2c


! End Module
! =============================================================================
        end module
