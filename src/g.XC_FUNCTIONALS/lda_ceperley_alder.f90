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

! lda_ceperley_alder
! Module Description
! ===========================================================================
!>       Given the input density, this subroutine returns the LDA exchange-
!! correlation energies and potentials. This LDA is according to the
!! parameterization of Ceperley-Alder.
!!
!! "Ground state of the electron gas by a stochastic method," by D.M. Ceperley
!!  and B.J. Alder, Phys. Rev. Let. 45:566-569 (1980).
!!
!! "Self-interaction correction to density-functional approximations for
!! many-electron systems," by J.P. Perdew and A. Zunger, Phys. Rev. B.
!! 23:5048-5079 (1981).
! ===========================================================================
! Code written by:
!> @author James P. Lewis
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
        subroutine lda_ceperley_alder (rh, exc, muxc, dexc, d2exc, dmuxc,    &
     &                                 d2muxc)

        implicit none
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: rh      !< input density
        
! Output
        real, intent (out) :: dexc   !< functional derivative of exc
        real, intent (out) :: d2exc  !< 2nd functional derivative of exc
        real, intent (out) :: dmuxc  !< functional derivative of muxc
        real, intent (out) :: exc    !< energy integral of xc
        real, intent (out) :: muxc   !< potential functional of xc
        real, intent (out) :: d2muxc !< 2nd functional derivative of muxc

! Parameters and Data Declaration
! ===========================================================================
        real, parameter :: delta_rh = 1.0d-6

! Variable Declaration and Description
! ===========================================================================
        real d2nec
        real d2nex
        real d3nec
        real d3nex
        real dec
        real ddec
        real d2dec
        real den
        real dden
        real d2den
        real d3den
        real ex
        real rho_third
        real rho
        real rhx
        real rs
        real rsl
        real sqrs

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize things to zero (thero in Spanish).
        exc = 0.0d0
        muxc = 0.0d0
        dexc = 0.0d0
        dmuxc = 0.0d0
        d2muxc = 0.0d0
        d2exc = 0.0d0

        rhx = sqrt (rh**2 + delta_rh)

! Convert to a.u.
        rho = rhx*(P_abohr**3)

! Find rho^(1/3)
        rho_third = rho**(1.0e0/3.0e0)

! Effective radius
        rs = 0.62035049d0/rho_third

! Find the energy, potential, and the deerivative of the potential.
        if (rho .lt. 0.23873241d0) then
          sqrs = sqrt(rs)
          den = 1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs       ! Effective density
          exc = -0.4581652d0/rs - 0.1423d0/den
          muxc = exc - rs*(0.15273333d0/rs**2                                &
     &               + (0.02497128d0/sqrs + 0.01581427d0)/den**2)

! Stuff for dmuxc
          dden = 1.0529d0/(2.0d0*sqrs) + 0.3334d0
          d2den = (-0.5d0)*1.0529d0/(2.0d0*rs*sqrs)
          d3den = (0.75d0)*1.0529d0/(2.0d0*(rs**2)*sqrs)
          dec = 0.1423d0*dden/(den**2)
          ddec = -2.0d0*0.1423d0*dden**2/(den**3) + 0.1423d0*d2den/(den**2)
          d2dec = 6.0d0*0.1423d0*(dden*3)/(den**4)                           &
     &           - 6.0d0*0.1423d0*dden*d2den/(den**3) + 0.1423d0*d3den/(den**2)

        else
          rsl = log(rs)
          exc = -0.4581652d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs      &
     &          + 0.002d0*rs*rsl
          muxc = exc - rs*(0.15273333d0/rs**2 + 0.01036667d0/rs              &
     &               - 0.003866667d0 + 0.00066667d0*(1.0d0 + rsl))

          dec = 0.0311d0/rs - 0.0116d0 + 0.0020d0*(rsl + 1.0d0)
          ddec = -0.0311d0/(rs**2) + 0.0020d0/rs
          d2dec = 2.0d0*0.0311d0/(rs**3) - 0.0020d0/(rs**2)
        end if

! Exchange-only energy and potential
        ex = -0.7385587664d0*rho_third

! Now find dmuxc
        dexc = (muxc - exc)/rho

! Compute dmu/dn. Use d(mu)/dn = 2*dexc/dn + n*d2(exc)/dn2.
        d2nec = (4.0d0*rs/(9.0d0*rho**2))*dec + (rs**2/(9.0d0*rho**2))*ddec
        d2nex = -(2.0d0/(9.0d0*rho**2))*ex
        dmuxc = 2.0d0*dexc + rho*(d2nex + d2nec)

! Compute d2mu/dn2, using d2(mu)/dn2 = 3*d2(exc)/dn2 + n*d3(exc)/dn3
        d3nec = (-28.0d0*rs/(27.0d0*rho**3))*dec + (-4.0d0*rs**2/       &
     &   (9.0d0*rho**3))*ddec + ((rs**3)/(-27.0d0*rho**3))*d2dec
        d3nex = (10.0d0/(27.0d0*rho**3))*ex
        d2muxc = 3.0*(d2nex + d2nec) + rho*(d3nex + d3nec)
        d2exc = d2nex + d2nec

! Convert output to eV (exc and muxc) and to eV*(Angstrom**3) (dexc and dmuxc)
        exc = exc*P_Hartree
        muxc = muxc*P_Hartree
        dexc = dexc*P_Hartree*(P_abohr)**3
        d2exc = d2exc*P_Hartree*(P_abohr)**6
        dmuxc = dmuxc*P_Hartree*(P_abohr)**3
        d2muxc = d2muxc*P_Hartree*(P_abohr)**6

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine lda_ceperley_alder
