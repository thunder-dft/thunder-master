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
        subroutine Depsilon_2c (r1, r2, z, eps, deps)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent(in) :: z              ! distance between atom 1 and atom 2

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
        integer imu
        integer ix

        real crossmag
        real r1dotr2
        real r1mag2
        real r2mag2
        real sum

        real crossa(3)

! Procedure
! ===========================================================================
! Initialize
        deps = 0.0d0

! If we are doing an atom, r1=r2. Then set deps2 to zero.
        if (z .lt. 1.0d-5) return

! Calculate (r2 - r1) - cross - d
        crossa(1) = r2(2)*r1(3) - r2(3)*r1(2)
        crossa(2) = r2(3)*r1(1) - r2(1)*r1(3)
        crossa(3) = r2(1)*r1(2) - r2(2)*r1(1)
        crossmag = sqrt(crossa(1)*crossa(1) + crossa(2)*crossa(2)             &
     &                  + crossa(3)*crossa(3))

! Crossmag cannot be zero for what follows. However, this should rarely happen
! because if crossmag = 0 we probably have an "on-top" instead of a true
! three-center system.
        if (abs(crossmag) .lt. 1.0d-3) then
          open (11, file = 'WARNINGS', status = 'unknown', position = 'append')
          write (11,*) ' *********** WARNING in Depsilon_2c ************ '
          write (11,*) ' Vectors sighat and dhat dangerously colinear '
          write (11,*) ' sigma - cross - r21 = ', crossmag
          write (11,*) ' setting all 3-center deps to zero '
          close (11)
          return
        end if

        r2mag2 = r2(1)**2 + r2(2)**2 + r2(3)**2
        r1mag2 = r1(1)**2 + r1(2)**2 + r1(3)**2
        r1dotr2 = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)

! Now calculate deps
! Now calculate deps/dratm
        do ix = 1, 3
          do imu = 1, 3
            sum = xlevi(ix,imu,1)*r2(1) + xlevi(ix,imu,2)*r2(2)               &
     &           + xlevi(ix,imu,3)*r2(3)

            deps(ix,imu,1) = (eps(imu,1)*eps(ix,3)*(1.0d0/z))                 &
     &                      - (eps(imu,1)*(1.0d0/crossmag)**2)                &
     &                       *(r2mag2*r1(ix) - r1dotr2*r2(ix))                &
     &                      + (1.0d0/z)*(1.0d0/crossmag)*(delk(imu,ix)        &
     &                       *(r2mag2 - r1dotr2) - r1(imu)*r2(ix)             &
     &                         - r2(imu)*r2(ix) + 2.0d0*r1(ix)*r2(imu))

            deps(ix,imu,2) = (1.0d0/crossmag)                                 &
     &                      *(sum + (eps(imu,2)*(1.0d0/crossmag))             &
     &                      *(r1dotr2*r2(ix) - r2mag2*r1(ix)))

            deps(ix,imu,3) = -(delk(imu,ix) - eps(imu,3)*eps(ix,3))*(1.0d0/z)
          end do
        end do
!
! Format Statements
! ===========================================================================
! None

        return
        end subroutine Depsilon_2c

