! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad Autonoma de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu

! RESTRICTED RIGHTS LEGEND
!
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.


! dpsiofr.f
! Program Description
! ===========================================================================
!       This function returns the values of the derivative d/dr psiofr(r)
! for the corresponding
! shell of the atomtype ispec.  The radial functions are normalized as:

!  int ( psiofr**2  r**2  dr ) = 1.0

! The wavefunctions for each atom must be read in by calling readpsi
! separately for each atom type.

! The value of r input to this function must be in angstrom units.
!
! ===========================================================================
! Code written from psiofr.f by
! Jose Ortega Mateo
! Dept. Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function dpsiofr (ispec, issh, r)
        implicit none
        real*8 dpsiofr

        include '../parameters.inc'
        include '../wavefunctions.inc'

! Argument Declaration and Description
! ===========================================================================
        integer ispec
        integer issh

        real*8 r

! Local Parameters and Data Declaration
! ===========================================================================
        integer norder
        parameter (norder = 5)

        integer norder2
        parameter (norder2 = norder/2)

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid
        integer max_points

        real*8 L(0:norder),mu(0:norder),Z(0:norder),alpha(0:norder)
        real*8 a(0:norder),b(0:norder),c(0:norder),d(0:norder)
        real*8 h
        integer iam
        integer i,j
        real*8 xmin
        real*8 xxp
        real*8 aaa,bbb,ccc,ddd
        real*8 dp

! Procedure
! ===========================================================================
! note : the points are equally spaced
        h=drr(issh,ispec)

        imid = int(r/h) + 1
        max_points = npoints(issh,ispec)

! Special cases
        if (r .ge. rrc(issh,ispec)) then
          dp = psi(max_points,issh,ispec) -psi(max_points-1,issh,ispec)
          dpsiofr = dp/h
          return
         else if (r .le. 0.0d0) then
          dp = psi(2,issh,ispec) -psi(1,issh,ispec)
          dpsiofr = dp/h
          return
        end if


        if (superspline) goto 1234

! Find starting point for the interpolation
        ileft = imid - norder2
        if (ileft .lt. 1) then
         ileft = 1
        else if (ileft + norder .gt. max_points) then
         ileft = max_points - norder
        end if

! Now interpolate with "natural" splines with f''(x)=0 at end points

        do i=0,norder
          a(i)=psi(i+ileft,issh,ispec)
        end do

        do i=1,norder-1
          alpha(i)=3.0d0*(a(i+1)-2*a(i)+a(i-1))/h
        end do

        L(0)=1
        mu(0)=0
        Z(0)=0
        c(0)=0
        do i=1,norder-1
          L(i)=(4.0d0-mu(i-1))*h
          mu(i)=h/L(i)
          Z(i)=(alpha(i)-h*Z(i-1))/L(i)
        end do
        L(norder)=1
        mu(norder)=0
        Z(norder)=0
        c(norder)=0

!       What curve section do we use?
        iam=imid-ileft

!       Don't need 0 to iam-1
        do j=norder-1,iam,-1
          c(j)=z(j)-mu(j)*c(j+1)
          b(j)=(a(j+1)-a(j))/h-h*(c(j+1)+2.0d0*c(j))/3.0d0
          d(j)=(c(j+1)-c(j))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp=(r-(xmin+(imid-1)*h))

!       psiofr=a(iam)+b(iam)*xxp+c(iam)*xxp**2+d(iam)*xxp**3
        dpsiofr=b(iam)+2.0d0*c(iam)*xxp+3.0d0*d(iam)*xxp**2

        return

!
! Cubic splines:  One big "super spline"
!
 1234   continue
        if (imid .gt. max_points) imid=max_points ! If we have gone off of the end
        aaa=psi_spline(1,imid,issh,ispec)
        bbb=psi_spline(2,imid,issh,ispec)
        ccc=psi_spline(3,imid,issh,ispec)
        ddd=psi_spline(4,imid,issh,ispec)

        xxp=r-(imid-1)*h

!        psiofr=aaa+bbb*xxp+ccc*xxp**2+ddd*xxp**3
        dpsiofr=bbb+2.0d0*ccc*xxp+3.0d0*ddd*xxp**2

! Format Statements
! ===========================================================================

        return
        end
