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
        subroutine get_potxc_1c (iexc, xc_fraction, r, rho, rhop, rhopp,arhop, &
     &                           newexc, vpxc, dnuxc, dnuxcs, dexc)
        implicit none

!! Argument Declaration and Description
! ===========================================================================
! Input
        integer iexc

        real, intent (in) :: xc_fraction
        real, intent (inout) :: r        ! radial distance

        ! density and derivatives
        real*8, intent (inout) :: rho      ! return zero value if rho is small
        real*8, intent (in) :: rhop
        real*8, intent (in) :: rhopp
        real*8, intent (in) :: arhop
! Output
        real, intent (out) :: newexc
        real, intent (out) :: vpxc
        real, intent (out) :: dnuxc
        real, intent (out) :: dnuxcs
        real, intent (out) :: dexc

! Local Parameters and Data Declaration
! ===========================================================================
	real*8, parameter :: thrd=1.d0/3.d0


! Local Variable Declaration and Description
! ===========================================================================
	integer ispin
	integer lpot, lgga
        
	real*8 dummy

        real*8 exc                     ! exchange-correlation energy
        real*8 ex                      ! exchange potential
        real*8 dec, dex                ! derivative of correlation and exchange
        real*8 fx
        real*8 fxc
        
        real*8 dup, ddn, d
        real*8 dpup, dpdn
        real*8 dppup, dppdn
        real*8 adpup, adpdn            !grad(abs(grad rho)) first derivative of abs of grad of rho

        real*8 agrup, agrdn, agrad
        real*8 delgrup, delgrdn, delgrad
        real*8 uplap, dnlap
        
        real*8 exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd			!Outputs
        real*8 expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91                  !Outputs
        real*8 expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe			!Outputs
        
        real*8 pi
                
! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        dummy = xc_fraction
        dummy = rhop
        dummy = rhopp
        dummy = arhop

! Initialize variables
        newexc = 0.0d0
        vpxc = 0.0d0
        dnuxc = 0.0d0
        dnuxcs = 0.0d0
        dexc = 0.0d0
	pi = 4.0d0*atan(1.0d0)

! If r is real*8ly small, then set to manageably small number.
        if (r .lt. 1.0d-4) r = 1.0d-4

! Rho must be positive, but not too small
        if (rho .lt. 1.0d-8) then
          rho = 0.0d0
          return
        else if (rho .lt. 1.0d-5) then
          rho = 1.0d-5
        end if

     
! XC Perdew-Burke-Ernzerhof GGA
! as parameterized by J.P. Perdew, K. Burke, M. Ernzerhof, Phys Rev Lett 77, 3865 (1996)
        if (iexc .eq. 6) then
        
!	PBE Parameters
	lgga = 1.0d0             !flag to do gga (0=>LSD only)
	lpot = 1.0d0             !flag to do potential (0=>energy only)

	
! Input variables are assigned to local variables for clarity
	adpup = 0.5d0*arhop  ! grad(abs(grad rho_up))
	adpdn = 0.5d0*arhop  ! grad(abs(grad rho_dn))
	dpup = 0.5d0*rhop   !grad(rho_up)
	dpdn = 0.5d0*rhop   !grad(rho_dn)
	dppup = 0.5d0*rhopp
	dppdn = 0.5d0*rhopp

	
! Input variab required by easypbe subroutine
					      ! As per the definitions of the original code
	dup = 0.5d0*rho  ! spin up density    !DUP=FUP*(ZUP**3/PI)*DEXP(-2.D0*ZUP*R)
	ddn = 0.5d0*rho  ! spin dn density    !DDN=FDN*(ZDN**3/PI)*DEXP(-2.D0*ZDN*R)
	agrup = abs(dpup)                          !agrup=2.d0*zup*dup          
	delgrup = dpup*adpup                  !8.d0*(zup**3)*dup*dup
        uplap = (dpup/r) + dppup      !4.d0*zup*dup*(zup-1.d0/r)
	
	agrdn = abs(dpdn)                          !2.d0*zdn*ddn
	delgrdn = dpdn*adpdn                  !8.d0*(zdn**3)*ddn*ddn
	dnlap = (dpdn/r) + dppdn              !4.d0*zdn*ddn*(zdn-1.d0/r)
        d = dup+ddn
        agrad = abs(rhop)                     ! agrad=|grad rho|
	delgrad = rhop*arhop                  !4.d0*agrad*(zup**2*dup+zdn**2*ddn)
	
        
    call  easypbe(dup,agrup,delgrup,uplap,ddn,agrdn,delgrdn, &		!Inputs
                dnlap,agrad,delgrad,1,1, &					!Inputs
                exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd, &			!Outputs
                expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &		!Outputs
                expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)			!Outputs
                
    newexc = expbe + ecpbe       ! exchange + corellation energy densities
	vpxc = vxuppbe + vxdnpbe + vcuppbe + vcdnpbe
!          vpxc = 0.0d0  
    newexc = real(newexc)
	vpxc = real(vpxc)
    dnuxc = real(dnuxc)
    dnuxcs = real(dnuxcs)
    dexc = real(dexc)


          
! outputs:
! exlsd = LSD exchange energy density, so that
! ExLSD = int d^3r rho(r) exlsd(r)
! vxuplsd = up LSD exchange potential
! vxdnlsd = down LSD exchange potential
! exclsd = LSD exchange-correlation energy density
! vxcuplsd = up LSD exchange-correlation potential
! vxcdnlsd = down LSD exchange-correlation potential

! expw91,vxuppw91,vxdnpw91,ecpw91,etc.=PW91 quantities

! expbe = PBE exchange energy density                             output 1
! vxuppbe = up PBE exchange potential
! vxdnpbe = down PBE exchange potential
! ecpbe = PBE exchange-correlation energy density                 output 2   -> output2 + output1 = exchange correlation energy density => final output
! vcuppbe = up  PBE exchange-correlation potential
! vcdnpbe = down  PBE exchange-correlation potential

! If the improper iexc option was entered then the program will stop.
        else
          write (*,*) ' In get_potxc_1c.f90 - '
          write (*,*) ' stop: xc option not implemented', iexc
          stop
        end if

! Calculate the exchange energy by combining the exchange and correlation
! energies and subtracting the exchange/correlation potential energy

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
 
! c----------------------------------------------------------------------
! c######################################################################
! c----------------------------------------------------------------------
      subroutine easypbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap, &
                agr,delgr,lcor,lpot, &
                exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd, &
                expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &
                expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c EASYPBE is a driver for the PBE subroutines, using simple inputs
! c K. Burke, May 14, 1996.
! c inputs: up=up density
! c	: agrup=|grad up|
! c	: delgrup=(grad up).(grad |grad up|) 
! c	: uplap=grad^2 up=Laplacian of up
! c	: dn,agrdn,delgrdn,dnlap=corresponding down quantities
! c	: agr=|grad rho|
! c	: delgr=(grad rho).(grad |grad rho|) 
! c	: lcor=flag to do correlation(=0=>don't)
! c	: lpot=flag to do potential(=0=>don't)
! c outputs: exlsd=LSD exchange energy density, so that
! c		ExLSD=int d^3r rho(r) exlsd(r)
! c	 : vxuplsd=up LSD exchange potential
! c	 : vxdnlsd=down LSD exchange potential
! c        : exclsd=LSD exchange-correlation energy density
! c	 : vxcuplsd=up LSD exchange-correlation potential
! c	 : vxcdnlsd=down LSD exchange-correlation potential
! c        : expw91,vxuppw91,vxdnpw91,ecpw91,etc.=PW91 quantities
! c        : expbe,vxuppbe,vxdnpbe,ecpbe,etc.=PBE quantities
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c needed constants:
! c pi32=3 pi**2
! c alpha=(9pi/4)**thrd
!       implicit real*8 (a-h,o-z)
      implicit none
      
      real*8, parameter :: thrd=1.d0/3.d0
      real*8, parameter :: thrd2=2.d0*thrd
      real*8, parameter :: pi32=29.608813203268075856503472999628d0
      real*8, parameter :: pi=3.1415926535897932384626433832795d0
      real*8, parameter :: alpha=1.91915829267751300662482032624669d0
      
      real*8 rho2, up, fk, s, u ,v, dn, agrdn,delgrdn,dnlap, rho
      real*8 exdnlsd, vxdnlsd, exdnpw91, vxdnpw91, exdnpbe, vxdnpbe
      real*8 exuplsd, vxuplsd, exuppw91, vxuppw91,exuppbe, vxuppbe
      real*8 exlsd,expw91, expbe, zet, g, rs, sk, twoksg, t, uu, rholap
      real*8 uplap, vv, ww, agr, eclsd, ecpbe, vcuplsd, vcuppbe,vcdnpbe
      real*8 ec,vcup,vcdn, H,DVCUP,DVCDN
      real*8 ECRS, ECZET, ALFC, ecpw91, vcuppw91, vcdnpw91
      real*8 agrup, delgrup, delgr, vcdnlsd
      integer lpot, lcor
      
      
    
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c PBE exchange
! c use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
! c do up exchange
! c fk=local Fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3) 
! c s=dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
! c u=delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
! c v=Laplacian/(rho*(2*fk)**2)_(rho=2*up)
      rho2=2.d0*up
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        
        s=2.d0*agrup/(2.d0*fk*rho2)
        u=4.d0*delgrup/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*uplap/(rho2*(2.d0*fk)**2)
        
        call exchpbe(rho2,s,u,v,0,lpot,exuplsd,vxuplsd)
        call exchpw91(rho2,s,u,v,exuppw91,vxuppw91)
        call exchpbe(rho2,s,u,v,1,lpot,exuppbe,vxuppbe)
      else
	exuplsd=0.d0
	vxuplsd=0.d0
	exuppw91=0.d0
	vxuppw91=0.d0
	exuppbe=0.d0
	vxuppbe=0.d0
      endif
! c repeat for down
      rho2=2.d0*dn
      if(rho2.gt.1d-18)then
        fk=(pi32*rho2)**thrd
        
        s=2.d0*agrdn/(2.d0*fk*rho2)
        u=4.d0*delgrdn/(rho2*rho2*(2.d0*fk)**3)
        v=2.d0*dnlap/(rho2*(2.d0*fk)**2)
        
        call exchpbe(rho2,s,u,v,0,lpot,exdnlsd,vxdnlsd)
        call exchpw91(rho2,s,u,v,exdnpw91,vxdnpw91)
        call exchpbe(rho2,s,u,v,1,lpot,exdnpbe,vxdnpbe)
      else
	exdnlsd=0.d0
	vxdnlsd=0.d0
	exdnpw91=0.d0
	vxdnpw91=0.d0
	exdnpbe=0.d0
	vxdnpbe=0.d0
      endif
10    continue 
! c construct total density and contribution to ex
      rho=up+dn
      exlsd=(exuplsd*up+exdnlsd*dn)/rho
      expw91=(exuppw91*up+exdnpw91*dn)/rho
      expbe=(exuppbe*up+exdnpbe*dn)/rho
      if(lcor.eq.0)return
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c Now do correlation
! c zet=(up-dn)/rho
! c g=phi(zeta)
! c rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk
! c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
! c twoksg=2*Ks*phi
! c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
! c uu=delgrad/(rho^2*twoksg^3)
! c rholap=Laplacian
! c vv=Laplacian/(rho*twoksg^2)
! c ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
! c ec=lsd correlation energy
! c vcup=lsd up correlation potential
! c vcdn=lsd down correlation potential
! c h=gradient correction to correlation energy
! c dvcup=gradient correction to up correlation potential
! c dvcdn=gradient correction to down correlation potential
      if(rho.lt.1.d-18)return
      zet=(up-dn)/rho
      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      fk=(pi32*rho)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      
      t=agr/(twoksg*rho)
      uu=delgr/(rho*rho*twoksg**3)
      rholap=uplap+dnlap
      vv=rholap/(rho*twoksg**2)
      ww=(agrup**2-agrdn**2-zet*agr**2)/(rho*rho*twoksg**2)
      
      call CORPBE(RS,ZET,T,UU,VV,WW,1,lpot,ec,vcup,vcdn, &
                       H,DVCUP,DVCDN)
      eclsd=ec
      ecpbe=ec+h
      vcuplsd=vcup
      vcdnlsd=vcdn
      vcuppbe=vcup+dvcup
      vcdnpbe=vcdn+dvcdn
      call CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
      call CORPW91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H,DVCUP,DVCDN)
      ecpw91=ec+h
      vcuppw91=vcup+dvcup
      vcdnpw91=vcdn+dvcdn
      return
      end subroutine easypbe
! c----------------------------------------------------------------------
! c######################################################################
! c----------------------------------------------------------------------
      SUBROUTINE EXCHPBE(rho,S,U,V,lgga,lpot,EX,VX)
! c----------------------------------------------------------------------
! C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
! c  K Burke's modification of PW91 codes, May 14, 1996
! c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! C  INPUT rho : DENSITY
! C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
! C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
! C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
! c   (for U,V, see PW86(24))
! c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
! c  input lpot:  (=0=>don't get potential and don't need U and V)
! C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c References:
! c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
! c     {\bf 40},  3399  (1989) (E).
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c Formulas:
! c   	e_x[unif]=ax*rho^(4/3)  [LDA]
! c ax = -0.75*(3/pi)^(1/3)
! c	e_x[PBE]=e_x[unif]*FxPBE(s)
! c	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! c uk, ul defined after [a](13) 
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
!       IMPLICIT real*8 (A-H,O-Z)
      implicit none
      real*8, parameter :: thrd=1.d0/3.d0
      real*8, parameter :: thrd4=4.d0/3.d0
      real*8, parameter :: pi=3.14159265358979323846264338327950d0
      real*8, parameter :: ax=-0.738558766382022405884230032680836d0
      real*8, parameter :: um=0.2195149727645171d0
      real*8, parameter :: uk=0.8040d0
      real*8, parameter :: ul=um/uk
      
      real*8 exunif, rho, ex, vx, s2,s,u,v,P0
      real*8 FxPBE, Fs, Fss
      integer lgga,  lpot 
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c construct LDA exchange energy density
      exunif = AX*rho**THRD
      if(lgga.eq.0)then
	ex=exunif
        vx=ex*thrd4
	return
      endif
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c construct PBE enhancement factor
      S2 = S*S
      P0=1.d0+ul*S2
      FxPBE = 1d0+uk-uk/P0
      EX = exunif*FxPBE
      if(lpot.eq.0)return
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! C  ENERGY DONE. NOW THE POTENTIAL:
! c  find first and second derivatives of Fx w.r.t s.
! c  Fs=(1/s)*d FxPBE/ ds
! c  Fss=d Fs/ds
      Fs=2.d0*uk*ul/(P0*P0)
      Fss=-4.d0*ul*S*Fs/P0
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c calculate potential from [b](24) 
      VX = exunif*(THRD4*FxPBE-(U-THRD4*S2*s)*FSS-V*FS)
      RETURN
      END subroutine EXCHPBE
! c----------------------------------------------------------------------
! c######################################################################
! c----------------------------------------------------------------------
      SUBROUTINE CORPBE(RS,ZET,T,UU,VV,WW,lgga,lpot,ec,vcup,vcdn, &
                       H,DVCUP,DVCDN)
! c----------------------------------------------------------------------
! c  Official PBE correlation code. K. Burke, May 14, 1996.
! C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
! C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
! C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
! C       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
! C       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
! C       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
! c       :  UU,VV,WW, only needed for PBE potential
! c       : lgga=flag to do gga (0=>LSD only)
! c       : lpot=flag to do potential (0=>energy only)
! c  output: ec=lsd correlation energy from [a]
! c        : vcup=lsd up correlation potential
! c        : vcdn=lsd dn correlation potential
! c        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
! c        : dvcup=nonlocal correction to vcup
! c        : dvcdn=nonlocal correction to vcdn
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c References:
! c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
! c     {\sl Generalized gradient approximation made simple}, sub.
! c     to Phys. Rev.Lett. May 1996.
! c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl real*8-space cutoff
! c     construction of a generalized gradient approximation:  The PW91
! c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
!       IMPLICIT real*8 (A-H,O-Z)
      implicit none
! c thrd*=various multiples of 1/3
! c numbers for use in LSD energy spin-interpolation formula, [c](9).
! c      GAM= 2^(4/3)-2
! c      FZZ=f''(0)= 8/(9*GAM)
! c numbers for construction of PBE
! c      gamma=(1-log(2))/pi^2
! c      bet=coefficient in gradient expansion for correlation, [a](4).
! c      eta=small number to stop d phi/ dzeta from blowing up at 
! c          |zeta|=1.
      real*8, parameter :: thrd=1.d0/3.d0
      real*8, parameter :: thrdm=-thrd
      real*8, parameter :: thrd2=2.d0*thrd
      real*8, parameter :: sixthm=thrdm/2.d0
      real*8, parameter :: thrd4=4.d0*thrd
      real*8, parameter :: gam=0.5198420997897463295344212145565d0
      real*8, parameter :: fzz=8.d0/(9.d0*gam)
      real*8, parameter :: gamma=0.03109069086965489503494086371273d0
      real*8, parameter :: bet=0.06672455060314922d0
      real*8, parameter :: delt=bet/gamma
      real*8, parameter :: eta=1.d-12
      
      real*8 rs,zet,t,uu,vv,ww,ec,vcup,vcdn, h,dvcup,dvcdn
      integer lgga,lpot
      real*8 rtrs, z4, eu, eurs, ep, eprs, alfm,alfrsm,alfc
      real*8  f, ecrs, fz,  comm
      
      real*8 g, g3, pon, b,b2, t2, t4, rs2,rs3
      real*8 q4, q5, g4, t6, rsthrd,gz, fac,bg,bec,q8,q9,hb 
      real*8 hrs,fact0,fact1,hbt,hrst,hz,ht,hzt,fact2,fact3
      real*8 htt, pref, fact5, eczet
      
      
      real*8 num11,num12,num13,num14,num15,num16
      real*8 num21,num22,num23,num24,num25, num26
      real*8 num31,num32,num33,num34,num35,num36
      
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c find LSD energy contributions, using [c](10) and Table I[c].
! c EU=unpolarized LSD correlation energy
! c EURS=dEU/drs
! c EP=fully polarized LSD correlation energy
! c EPRS=dEP/drs
! c ALFM=-spin stiffness, [c](3).
! c ALFRSM=-dalpha/drs
! c F=spin-scaling factor from [c](9).
! c construct ec, using [c](8)

num11=0.0310907D0;num12=0.21370D0;num13=7.5957D0;num14=3.5876D0;num15=1.6382D0;num16= 0.49294D0;
num21=0.01554535D0;num22=0.20548D0;num23=14.1189D0;num24=6.1977D0;num25=3.3662D0; num26=0.62517D0;
num31=0.0168869D0;num32=0.11125D0;num33=10.357D0;num34=3.6231D0;num35=0.88026D0;num36=0.49671D0;

      rtrs=sqrt(rs)
      CALL gcor2(num11,num12,num13,num14,num15,num16,rtrs,EU,EURS)
      CALL gcor2(num21,num22,num23,num24,num25,num26,rtRS,EP,EPRS)
      CALL gcor2(num31,num32,num33,num34,num35,num36,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c LSD potential from [c](A1)
! c ECRS = dEc/drs [c](A2)
! c ECZET=dEc/dzeta [c](A3)
! c FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
             -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if(lgga.eq.0)return
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! c PBE correlation energy
! c G=phi(zeta), given after [a](3)
! c DELT=bet/gamma
! c B=A of [a](8)
      G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(EXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      H = G3*(BET/DELT)*LOG(1.D0+DELT*Q4*T2/Q5)
      if(lpot.eq.0)return
! c----------------------------------------------------------------------
! c----------------------------------------------------------------------
! C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.D0
      GZ=(((1.d0+zet)**2+eta)**sixthm- &
     ((1.d0-zet)**2+eta)**sixthm)/3.d0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      hRST = RSTHRD*T2*hBT*BEC*ECRS
      hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.d0*BET*G3*Q9/Q8
      hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      RETURN
      END subroutine CORPBE
! c----------------------------------------------------------------------
! c######################################################################
! c----------------------------------------------------------------------
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! c slimmed down version of GCOR used in PW91 routines, to interpolate
! c LSD correlation energy, as given by (10) of
! c J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! c K. Burke, May 11, 1996.
!       IMPLICIT real*8 (A-H,O-Z)
      implicit none
      real*8 a,a1,b1,b2,b3,b4,rtrs
      real*8 q0,q1,q2,gg,q3,ggrs
      
      Q0 = -2.D0*A*(1.D0+A1*rtrs*rtrs)
      Q1 = 2.D0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = LOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.D0*B2+rtrs*(3.D0*B3+4.D0*B4*rtrs))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.d0+Q1))
      RETURN
      END subroutine GCOR2
! c----------------------------------------------------------------------
! c######################################################################
! c----------------------------------------------------------------------
      SUBROUTINE EXCHPW91(D,S,U,V,EX,VX)
! C  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
! C  INPUT D : DENSITY
! C  INPUT S:  ABS(GRAD D)/(2*KF*D)
! C  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
! C  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
! C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!       IMPLICIT real*8 (A-H,O-Z)
      implicit none
      real*8, parameter :: a1=0.19645D0
      real*8, parameter :: a2=0.27430D0
      real*8, parameter :: a3=0.15084D0
      real*8, parameter :: a4=100.d0
      real*8, parameter :: ax=-0.7385588D0
      real*8, parameter :: a=7.7956D0
      real*8, parameter :: b1=0.004d0
      real*8, parameter :: thrd=0.333333333333D0
      real*8, parameter :: thrd4=1.33333333333D0
      
      real*8 d,s,u,v
      real*8 fac,s2,s3,s4,p0,p1,p2,p3,p4,f,ex
      real*8 p5, p6, p7, fs, p8, p9, p10, p11, fss, vx
      
! c for Becke exchange, set a3=b1=0
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/SQRT(1.D0+A*A*S2)
      P1 = LOG(A*S+1.D0/P0)
      P2 = EXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
! C  LOCAL EXCHANGE OPTION
! C     EX = FAC
! C  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
! C  LOCAL EXCHANGE OPTION:
! C     VX = FAC*THRD4
      RETURN
      END subroutine EXCHPW91
! c----------------------------------------------------------------------
! c######################################################################
! c----------------------------------------------------------------------
      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
! C  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
! C  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
! C  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
! C     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
! C  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
!       IMPLICIT real*8 (A-H,O-Z)
      implicit none
      real*8, parameter :: gam=0.5198421D0
      real*8, parameter :: fzz=1.709921D0
      real*8, parameter :: thrd=0.333333333333D0
      real*8, parameter :: thrd4=1.333333333333D0
      
      real*8 rs,zet,eu, eurs, ep, eprs, alfm, alfrsm
      real*8 f, alfc, z4, ec, ecrs, fz,eczet, comm,vcup,vcdn

      real*8 num11,num12,num13,num14,num15,num16,num17
      real*8 num21,num22,num23,num24,num25,num26,num27
      real*8 num31,num32,num33,num34,num35,num36,num37
      
      num11 = 0.0310907D0; num12=0.21370D0; num13=7.5957D0; num14=3.5876D0; num15=1.6382D0; num16=0.49294D0; num17=1.00D0;
      num21 = 0.01554535D0; num22=0.20548D0; num23=14.1189D0; num24=6.1977D0; num25=3.3662D0; num26=0.62517D0; num27=1.00D0;
      num31 = 0.0168869D0;num32 =0.11125D0; num33 =10.357D0; num34 =3.6231D0; num35 =0.88026D0; num36 =0.49671D0; num37 = 1.00D0;
      
      
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      CALL GCOR(num11,num12,num13,num14,num15,num16,num17,RS,EU,EURS)
      CALL GCOR(num21,num22,num23,num24,num25,num26,num27,RS,EP,EPRS)
      CALL GCOR(num31,num32,num33,num34,num35,num36,num37,RS,ALFM,ALFRSM)
! C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
! C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
             -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END subroutine CORLSD
! c----------------------------------------------------------------------
! c######################################################################
! c----------------------------------------------------------------------
      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
! C  CALLED BY SUBROUTINE CORLSD
!       IMPLICIT real*8 (A-H,O-Z)
      implicit  none
      real*8 a,a1,b1,b2,b3,b4,rs
      real*8 p1,p,q0, rs12, rs32, rsp, q1,q2,gg,q3,ggrs       
      
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = SQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = LOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END subroutine GCOR
! c----------------------------------------------------------------------
! c######################################################################
! c----------------------------------------------------------------------
      SUBROUTINE CORpw91(RS,ZET,G,EC,ECRS,ECZET,T,UU,VV,WW,H, &
                        DVCUP,DVCDN)
! C  pw91 CORRELATION, modified by K. Burke to put all arguments 
! c  as variables in calling statement, rather than in common block
! c  May, 1996.
! C  INPUT RS: SEITZ RADIUS
! C  INPUT ZET: RELATIVE SPIN POLARIZATION
! C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
! C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
! C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
! C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
! C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
! C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
!       IMPLICIT real*8 (A-H,O-Z)
      implicit none
      real*8, parameter :: xnu = 15.75592D0
      real*8, parameter :: cc0 = 0.004235D0
      real*8, parameter :: cx = -0.001667212D0
      real*8, parameter :: alf = 0.09D0
      real*8, parameter :: c1 = 0.002568D0
      real*8, parameter :: c2 = 0.023266D0
      real*8, parameter :: c3 = 7.389D-6
      real*8, parameter :: c4 = 8.723D0
      real*8, parameter :: c5 = 0.472D0
      real*8, parameter :: c6 = 7.389D-2
      real*8, parameter :: a4 = 100.D0
      real*8, parameter :: thrdm = -0.333333333333D0
      real*8, parameter :: thrd2 = 0.666666666667D0
      
      real*8 rs,zet,g,ec,ecrs,eczet,t,uu,vv,ww
      real*8 bet, delt
      real*8 g3,g4,pon,b,b2,t2,t4,t6,rs2,rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3,h0
      real*8 h1,h,ccrs,rsthrd,r4,gz,fac,bg,bec,q8,q9,h0b,h0rs 
      real*8 fact0,fact1,h0bt,h0rst,h0z,h0t,h0zt,fact2,fact3,h0tt,h1rs,fact4,h1rst,h1z
      real*8 h1t,h1zt,h1tt,hrs,hrst,ht,htt,hz,hzt,comm,pref,fact5,dvcup,dvcdn 
            

      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(EXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = 0.663436444d0*rs
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = EXP(-R1*T2)
      H0 = G3*(BET/DELT)*LOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
! C  LOCAL CORRELATION OPTION:
! C     H = 0.0D0
! C  ENERGY DONE. NOW THE POTENTIAL:
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
! C  LOCAL CORRELATION OPTION:
! C     DVCUP = 0.0D0
! C     DVCDN = 0.0D0
      RETURN
      END subroutine CORPW91
! c----------------------------------------------------------------------

end module
