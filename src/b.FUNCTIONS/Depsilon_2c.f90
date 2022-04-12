! copyright info:
!
!                             @Copyright 2012
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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

! Depsilon_2c.f90
! Program Description
! ===========================================================================
!       This subroutine sets up deps/dr1 in the two-center molecular system
! defined by sighat=(r2-r1)/|r2-r1| and piprimehat=(r2-cross-r1)/|r2-cross-r1|.
! The eps matrix is obtained by calling epsiln(r1,rvec) where rvec=r2-r1.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Depsilon_2c (r1, r2, eps, deps)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input:
! r1    - position of atom 1
! r2    - posiiton of atom 2
! eps  - 2 center epsilon with z=r2-r1
        real, intent(in) :: r1(3), r2(3), eps (3, 3)

! Output:
! deps2 - deps/dr1
        real, intent(out) :: deps (3, 3, 3)

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer i,ii,ix

        real r2mag2,r2mag,r1mag,denom
        real crossmag,dd,dot,term,ddinv,crossinv
        real crossa(3)

! Procedure
! ===========================================================================
! Initialize
        deps = 0.0d0

! If we are doing an atom, r1=r2. Then set deps2 to zero.
        dd = sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2)
        if (dd.lt.1.0d-4)return
        ddinv = 1.0d0/dd

        r2mag2 = r2(1)**2 + r2(2)**2 + r2(3)**2
        r2mag = sqrt(r2mag2)
        r1mag = sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)
        crossa(1) = r2(2)*r1(3)-r2(3)*r1(2)
        crossa(2) = r2(3)*r1(1)-r2(1)*r1(3)
        crossa(3) = r2(1)*r1(2)-r2(2)*r1(1)
        crossmag = sqrt(crossa(1)**2+crossa(2)**2+crossa(3)**2)
        dot = 0.e0
        do i=1,3
          dot=dot+r1(i)*r2(i)
        end do
        denom=r1mag*r2mag

! check to see if atoms are colinear
        if(abs(crossmag).lt.1.0d-3)then
          return
        endif
        crossinv=1.0/crossmag
!
! now calculate deps
        do ii=1,3
         do ix=1,3
          term=xlevi(ix,ii,1)*r2(1)+   &
     &         xlevi(ix,ii,2)*r2(2)+   &
     &         xlevi(ix,ii,3)*r2(3)

          deps(ix,ii,1)=(eps(ii,1)*eps(ix,3)*ddinv)-(eps(ii,1)*  &
     &     crossinv**2)*(r2mag2*r1(ix)-dot*r2(ix))+  &
     &     ddinv*crossinv*(delk(ii,ix)*(r2mag2-dot)-  &
     &     r1(ii)*r2(ix)-r2(ii)*r2(ix)+2.e0*r1(ix)*r2(ii))

          deps(ix,ii,2)=crossinv*(term+(eps(ii,2)*crossinv)*  &
     &        (dot*r2(ix)-r2mag2*r1(ix)))

          deps(ix,ii,3)=-(delk(ii,ix)-eps(ii,3)*eps(ix,3))*ddinv
         end do
        end do
!
! Format Statements
! ===========================================================================
! None

        return
        end subroutine Depsilon_2c

