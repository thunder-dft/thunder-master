! copyright info:
!
!                             @Copyright 2016
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! West Virginia University - Barry Haycock
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


! M_Fdata_3c
! Module Description
! ===========================================================================
!>       This is a module containing all of the subroutines that will read in
!! all the three-center fdata from the fdata files in the Fdata directory.  It
!! contains the following subroutines within the module:
!!
!!       read_Fdata_3c.f90 - reads in the three-center datafiles.
!!       getMEs_Fdata_3c.f90 - interpolates the three-center data and gets the
!!                             matrix elements.
!!       getDMEs_Fdata_3c   - interpolates the three-center data and gets the
!!                            derivitive matrix elements.
!!       destroy_Fdata_3c.f90 - deallocates all the three-center fdata arrays.
!!
!! For a complete list of the interactions see the files 3c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Module Declaration
! ===========================================================================
        module M_Fdata_3c
        use M_species
        !use M_Drotations

        include '../include/interactions_3c.h'

! Type Declaration
! ===========================================================================
! three-center interactions arrays

! To cut down on storage space, we actually change the storage procedure
! from previous FIREBALL code. Not all atoms have the same number of
! interactions or interaction types. Before - we would store things based
! on the maximum number of fdata points, maximum number of interactions types,
! maximum number of matrix elements - so even hydrogen-hydrogen (just ss
! and/or ss*) stored a 4x4 or an 8x8 matrix even when not needed.  This was
! quite inefficient.

! The new approach is to define some Fdata types which store the actual
! Fdata points. The smallest unit storage is called Fdata_cell_3C, containing
! all Fdata for a particular interaction/subinteraction type.
        type T_Fdata_cell_3c
          integer nME                   ! number of non-zero matrix elements
          integer nx                    ! number of data points
          integer ny                    ! number of data points

          integer, pointer :: mu_3c (:)
          integer, pointer :: nu_3c (:)
          integer, pointer :: mvalue_3c (:)

          real xmax                     ! maximum interaction range
          real ymax                     ! maximum interaction range

          ! actual Fdata points
          real, pointer :: Fdata_3c (:, :, :)
        end type T_Fdata_cell_3c

! Fdata_bundle_3c is the 3-center package of Fdata_cell_3c. It contains all the
! Fdata_cell_3c information of all interaction/subtypes for given species pair.
        type T_Fdata_bundle_3c
          integer nFdata_cell_3c                     ! number of Fdata_cell3C
          integer, dimension (P_maxtype, 0:P_maxsubtype, P_maxtheta) :: index_3c

          ! actual fdata
          type (T_Fdata_cell_3c), pointer :: Fdata_cell_3c (:)
        end type T_Fdata_bundle_3c

        type(T_Fdata_bundle_3c), pointer :: Fdata_bundle_3c (:,:,:)

! module procedures
        contains

! ===========================================================================
! read_Fdata_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in all three-center interactions. These three-center
!! interactions are contributions as described in
!!
!! "Ab initio multicenter tight-binding model for molecular-dynamics
!!  simulations and other applications in covalent systems" by O.F. Sankey
!!  and D.J. Niklewski, PRB 40:3979-3995 (1989).
!!
!! "Further developments in the local-orbital density-functional theory
!!  tight-binding method: FIREBALL" by J.P. Lewis, K.R. Glaesemann, G.A. Voth,
!!  J. Fritsch, A.A. Demkov, J. Ortega, and O.F. Sankey, PRB 64:195103 (2001).
!!
!! "Multicenter approach to the exchange-correlation interactions in ab initio
!!  tight-binding methods" by P. Jelinek, H. Wang, J.P. Lewis, O.F. Sankey,
!!  and J. Ortega, PRB 71:235101 (2005).
!
! ===========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
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
        subroutine read_Fdata_3c ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies, kspecies  !< counters for number of species
        integer icell                         !< counter for Fdata files
        integer iME                           !< counter for number of MEs
        integer itype, isubtype, itheta       !< counter over types, subtypes
        integer iindex                        !< combined type, subtype
        integer ix, iy                        !< counters over distances

        integer nME3c_max                     !< temporary value for nME

        type(T_Fdata_cell_3c), pointer :: pFdata_cell
        type(T_Fdata_bundle_3c), pointer :: pFdata_bundle

        character (len = 48) filename

! Allocate Arrays
! ===========================================================================
        allocate (Fdata_bundle_3c (nspecies, nspecies, nspecies))

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,'(A)') 'Reading three-center interactions '
        write (ilogfile,'(A)') '--------------------------------- '
        write (ilogfile,*)

        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            do kspecies = 1, nspecies
              write (ilogfile,'(4x, A18, I4, I4, I4)') '- species triplet: ',   &
      &                     species(ispecies)%nZ, species(jspecies)%nZ,      &
      &                     species(kspecies)%nZ

              ! cut some lengthy notation
              pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)

              ! open directory file
              write (filename,'("/3c.",i2.2,".",i2.2,".",i2.2,".dat")')      &
      &              species(ispecies)%nZ, species(jspecies)%nZ,             &
      &              species(kspecies)%nZ
              open (unit = 11, file = trim(Fdata_location)//trim(filename),  &
      &             status = 'old')
              read (11,*) pFdata_bundle%nFdata_cell_3c
              allocate (pFdata_bundle%Fdata_cell_3c(pFdata_bundle%nFdata_cell_3c))
              pFdata_bundle%index_3c = -1

              write (filename,'("/3c.",i2.2,".",i2.2,".",i2.2,".dir")')      &
      &              species(ispecies)%nZ, species(jspecies)%nZ,             &
      &              species(kspecies)%nZ
              open (unit = 11, file = trim(Fdata_location)//trim(filename),  &
      &             status = 'old')
              ! loop over Fdata files.
              do icell = 1, pFdata_bundle%nFdata_cell_3c
                pFdata_cell=>pFdata_bundle%Fdata_cell_3c(icell)
                read (11, *) iindex, itype, isubtype, itheta, filename,      &
      &                      pFdata_cell%nME, pFdata_cell%nx,                &
      &                      pFdata_cell%xmax, pFdata_cell%ny,               &
      &                      pFdata_cell%ymax
                nME3c_max = pFdata_cell%nME
                allocate (pFdata_cell%Fdata_3c(pFdata_cell%nx, pFdata_cell%ny,&
      &                   nME3c_max))
                allocate (pFdata_cell%mu_3c(nME3c_max))
                allocate (pFdata_cell%nu_3c(nME3c_max))
                allocate (pFdata_cell%mvalue_3c(nME3c_max))

                ! make the index file
                pFdata_bundle%index_3c(itype, isubtype, itheta) = iindex

                ! open the actual datafile
                open (12, file = trim(Fdata_location)//'/'//trim(filename))

                ! read the actual fdata
                do iy = 1, pFdata_cell%ny    ! loop over bondlengths
                  do ix = 1, pFdata_cell%nx    ! loop over 3rd atom distance
                    read (12,*) (pFdata_cell%Fdata_3c(ix, iy, iME), iME = 1, nME3c_max)
                  end do
                end do
                close (12)

                ! open mu, nu, mvalue file
                write (filename,'("/",i2.2,"_munu_3c.",i2.2,".",i2.2,".",     &
     &                                i2.2,".dat")')                          &
     &            itype, species(ispecies)%nZ, species(jspecies)%nZ,          &
     &                    species(kspecies)%nZ
                open (unit = 13, file = trim(Fdata_location)//trim(filename), &
     &                status = 'old')

                ! read the mapping (mu,nu,mvalue)
                read (13,*) (pFdata_cell%mu_3c(iindex), iindex = 1, nME3c_max)
                read (13,*) (pFdata_cell%nu_3c(iindex), iindex = 1, nME3c_max)
                read (13,*) (pFdata_cell%mvalue_3c(iindex), iindex = 1, nME3c_max)
              end do ! datafile for (ispecies, jspecies, kspecies)
              close (11)
              close (13)
            end do
          end do
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine read_Fdata_3c


! ===========================================================================
! getMEs_Fdata_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine loops over all the Legendre angles and calls the
!! interpolation to find the values at the given angles.  The function f(x) is a
!! three-center interaction function read from read_Fdata_3c. After the data
!! is interpolated, the data is multiplied by the appropriate Legendre fundtion
!! and then 'recovered' into the matrix form.

! ===========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
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
        subroutine getMEs_Fdata_3c (ispecies, jspecies, kspecies, iint, isub,&
     &                              x, y, norb_mu, norb_nu, cost, hmbox)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: ispecies, jspecies, kspecies         !< species
        integer, intent(in) :: iint, isub            !< integral type, subtype
        integer, intent(in) :: norb_mu, norb_nu      !< FIXME: What am I?

        real, intent(in) :: x, y                     !< distances between pairs
        real, intent(in) :: cost                     !< cosine bond-charge angle

! Output
        real, intent(out), dimension (norb_mu, norb_nu) :: hmbox !< FIXME

! Variable Declaration and Description
! ===========================================================================
        integer iindex, imu, inu               !< counters for building matrix
        integer itheta                         !< which angle?
        integer mvalue                         !< value of quantum number m

        real sint                              !< sin of bond-charge angle
        real cost2                             !< cost**2

        real :: p(0:P_maxtheta - 1)            !< Legendre polys p(cost)
        real, dimension (:), allocatable :: Fdata      !< F(x)

        type(T_Fdata_cell_3c), pointer :: pFdata_cell
        type(T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Set up Legendre polynomials
        cost2 = cost**2
        sint = sqrt(max(1.0d-5, 1.0d0 - cost2))
        p(0) = 1.0d0
        p(1) = cost
        p(2) = 0.5d0*(3.0d0*cost2 -1.0d0)
        p(3) = 0.5d0*cost*(5.0d0*cost2 - 3.0d0)
        p(4) = 0.125d0*(3.0d0 - cost2*(30.0d0 - 35.0d0*cost2))

        ! initialize, finding the correct fdata bundle
        pFdata_bundle => Fdata_bundle_3c(ispecies, jspecies, kspecies)
        pFdata_cell =>                                                       &
     &    pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(iint,isub,1))
        allocate (Fdata(pFdata_cell%nME))
        Fdata = 0.0d0

        ! loop over theta
        do itheta = 1, P_maxtheta
          pFdata_cell =>                                                     &
     &     pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(iint,isub,itheta))
          call addMEs_Fdata_3c (pFdata_cell, x, y, p(itheta - 1),            &
     &                          pFdata_cell%nME, Fdata)
        end do

        ! recover matrix
        hmbox = 0.0d0
        do iindex = 1, pFdata_cell%nME
          imu = pFdata_cell%mu_3c(iindex)
          inu = pFdata_cell%nu_3c(iindex)
          mvalue = pFdata_cell%mvalue_3c(iindex)
          if (mvalue .ne. 1) then
            hmbox(imu,inu) = Fdata(iindex)
          else
            hmbox(imu,inu) = Fdata(iindex)*sint
          end if
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (Fdata)

! Format Statements
! ===========================================================================
        return
        end subroutine getMEs_Fdata_3c


! addMEs_Fdata_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine uses interpolation to find the value of f(x) for any x,
!! given an array of equally spaced points for f(x).  The function f(x) is a
!! three-center interaction function read from read_Fdata_3c.
!!
!! For polynomial interpolation see Mathews and Walker, p.329

! ===========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine addMEs_Fdata_3c (pFdata_cell, x, y, ptheta, nME, result)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nME         !< number of matrix elements

        real, intent(in) :: x, y            !< distances from atom of interest
        real, intent(in) :: ptheta          !< p is the Legendre multiplier

        type(T_Fdata_cell_3c), pointer :: pFdata_cell !< FIXME

! Input and Output
        real, intent(inout), dimension (nME) :: result !< FIXME

! Variable Declaration and Description
! ===========================================================================
        integer ix, iy
        integer k
        integer nME3c_max                   !< number of matrix elements

        real dx, dy, gx, gy, gg, px, py

! data arrays for interpolation
        real, allocatable, dimension (:) :: bb0, bb1, bb2, bb3
        real, allocatable, dimension (:) :: f1m1, f1m2, f1m3, f0p3, f0p6
        real, allocatable, dimension (:) :: f1p3, f1p6, f2p1
        real, allocatable, dimension (:, :) :: g

! Local Parameters and Data Declaration
! ===========================================================================
! tolerance (may be needed to avoid roundoff error in the calling program)
! if xin > xmax but within, say, .001% of xmax then ignore
        real, parameter :: P_tiny = 1.0d-10
        real, parameter :: P_small = 1.0d0-8

! Procedure
! ===========================================================================
! Initialize
        dx = pFdata_cell%xmax/(pFdata_cell%nx - 1.0d0)
        dy = pFdata_cell%ymax/(pFdata_cell%ny - 1.0d0)

        ix = int(x/dx) + 1
        iy = int(y/dy) + 1

        if (ix .lt. 2) ix = 2
        if (iy .lt. 2) iy = 2

        if (ix .gt. pFdata_cell%nx - 2) ix = pFdata_cell%nx - 2
        if (iy .gt. pFdata_cell%ny - 2) iy = pFdata_cell%ny - 2

        px = x/dx - ix + 1.0d0
        py = y/dy - iy + 1.0d0

! ***************************************************************************
! Adaptive interpolation - estimate gradient:
        gx = (pFdata_cell%Fdata_3c(ix + 1, iy, 1)                             &
     &        - pFdata_cell%Fdata_3c(ix, iy, 1))/dx
        gy = (pFdata_cell%Fdata_3c(ix, iy + 1, 1)                             &
     &        - pFdata_cell%Fdata_3c(ix, iy, 1))/dy
        gg = gx**2 + gy**2

! We choose one of three ways to interpolate based on the size of the
! gradient (derivative). If the value of the derivative is tiny, then
! do very simple interpolation. If small, then something a little more
! complicated. If neither, then something even MORE complicated.

! METHOD 1
! Do three point linear bivariate interpolation
! Handbook of Mathematical Functions..., edited by M. Abramowitz
! and I.A. Stegun, Dover edition, Pg. 882, Eq. 25.2.65
        if (gg .lt. P_tiny) then
          result = result + ptheta                                            &
     &                     *((1.0d0 - px - py)*pFdata_cell%Fdata_3c(ix,iy,:)  &
     &                                    + px*pFdata_cell%Fdata_3c(ix+1,iy,:)&
     &                                    + py*pFdata_cell%Fdata_3c(ix,iy+1,:))
! METHOD 2
! Do quadratic bivariate interpolation (six point formula, Eq. 25.2.67)
        else if (gg .lt. P_small) then
          result = result + ptheta                                            &
     &                      *(py*(py-1)*0.5d0*pFdata_cell%Fdata_3c(ix,iy-1,:) &
     &                      + px*(px-1)*0.5d0*pFdata_cell%Fdata_3c(ix-1,iy,:) &
     &        + (1.0d0 + px*py - px**2 - py**2)*pFdata_cell%Fdata_3c(ix,iy,:) &
     &      + px*(px - 2.0d0*py + 1.0d0)*0.5d0*pFdata_cell%Fdata_3c(ix+1,iy,:)&
     &      + py*(py - 2.0d0*px + 1.0d0)*0.5d0*pFdata_cell%Fdata_3c(ix,iy+1,:)&
     &                      + px*py*pFdata_cell%Fdata_3c(ix+1,iy+1,:))
        else

! METHOD 3
! Interpolate one direction, then interpolate using these values to get
! the final value.  Use Eq. 25.2.13 (4 point intepolater).
! Total number of points used is 16
          nME3c_max = pFdata_cell%nME

          ! Allocate the arrays used in Method 3
          allocate (bb0(nME3c_max), bb1(nME3c_max), bb2(nME3c_max))
          allocate (bb3(nME3c_max), f1m1(nME3c_max), f1m2(nME3c_max))
          allocate (f1m3(nME3c_max), f0p3(nME3c_max), f0p6(nME3c_max))
          allocate (f1p3(nME3c_max), f1p6(nME3c_max), f2p1(nME3c_max))
          allocate (g(-1:2,nME3c_max))
          do k = -1, 2
            f1m1 = pFdata_cell%Fdata_3c(ix+k,iy-1,:)
            f1m2 = 2.0d0*f1m1
            f1m3 = 3.0d0*f1m1

            f0p3 = 3.0d0*pFdata_cell%Fdata_3c(ix+k,iy,:)
            f0p6 = 2.0d0*f0p3

            f1p3 = 3.0d0*pFdata_cell%Fdata_3c(ix+k,iy+1,:)
            f1p6 = 2.0d0*f1p3

            f2p1 = pFdata_cell%Fdata_3c(ix+k,iy+2,:)

            bb3 = - f1m1 + f0p3 - f1p3 + f2p1
            bb2 = f1m3 - f0p6 + f1p3
            bb1 = - f1m2 - f0p3 + f1p6 - f2p1
            bb0 = f0p6

            g(k,:) = ((bb3*py + bb2)*py + bb1)*py + bb0
          end do

          f1m1 = g(-1,:)
          f1m2 = 2.0d0*f1m1
          f1m3 = 3.0d0*f1m1

          f0p3 = 3.0d0*g(0,:)
          f0p6 = 2.0d0*f0p3

          f1p3 = 3.0d0*g(1,:)
          f1p6 = 2.0d0*f1p3

          f2p1 = g(2,:)

          bb3 = - f1m1 + f0p3 - f1p3 + f2p1
          bb2 = f1m3 - f0p6 + f1p3
          bb1 = - f1m2 - f0p3 + f1p6 - f2p1
          bb0 = f0p6
          result = result + ptheta*(((bb3*px + bb2)*px + bb1)*px + bb0)/36.0d0

          ! Deallocate the arrays used in Method 3.
          deallocate (bb0, bb1, bb2, bb3)
          deallocate (f1m1, f1m2, f1m3, f0p3, f0p6)
          deallocate (f1p3, f1p6, f2p1, g)
        end if

! A final note, if you are interested in splines, you should start with
! "Handbook on SPLINES for the User"
! by Eugene V. Shikin and Alexander I. Plis.  1995, CRC Press.  Most other
! books on slines and bivariate interpolation are nothing but proofs and
! abstract math, but this one gives the real equations you need to get
! the actual work done.

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end subroutine addMEs_Fdata_3c

! ===========================================================================
! getDMEs_Fdata_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine loops over all the Legendre angles and calls the
!! interpolation to find the values of the derivitives w.r.t.x and y at the
!! given angles.  The function f(x) is a three-center interaction function
!! read from read_Fdata_3c.  After the data is interpolated, the data is
!! multiplied by the appropriate Legendre fundtion and then 'recovered' into
!! the matrix form.
!! This subroutine is similar ot the getMES, but returns derivitives.

! ===========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
!! @author Barry J. Haycock
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
        subroutine getDMEs_Fdata_3c (ispecies, jspecies, kspecies, iint, isub,&
     &                               x, y, norb_mu, norb_nu, cost, rhat,      &
     &                               sighat, hmbox, dphmbox, dxhmbox, dyhmbox)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies, jspecies, kspecies         !< species
        integer, intent (in) :: iint, isub            !< integral type, subtype
        integer, intent (in) :: norb_mu, norb_nu

        real, intent (in) :: x, y                     !< distances between pairs
        real, intent (in) :: cost                     !< cosine bond-charge angle
        real, intent (in), dimension (3) :: rhat
        real, intent (in), dimension (3) :: sighat

! Output
        real, intent (out), dimension (norb_mu, norb_nu) :: hmbox
        real, intent (out), dimension (norb_mu, norb_nu) :: Dphmbox
        real, intent (out), dimension (norb_mu, norb_nu) :: Dxhmbox
        real, intent (out), dimension (norb_mu, norb_nu) :: Dyhmbox

! Variable Declaration and Description
! ===========================================================================
        integer iindex, imu, inu               !< counters for building matrix
        integer itheta                         !< which angle?
        integer mvalue                         !< value of quantum number m

        real sint                              !< sin of bond-charge angle
        real cost2                             !< cost**2

        real :: p(0:P_maxtheta - 1)            !< Legendre polys p(cost)
        real :: Dp(0:P_maxtheta - 1)           !< Legendre polys Dp(cost)
        real, dimension (:), allocatable :: Fdata      !< F(x)
        real, dimension (:), allocatable :: DpFdata    !< dF(x)/dTheta
        real, dimension (:), allocatable :: DxFdata    !< dF(x)/dx
        real, dimension (:), allocatable :: DyFdata    !< dF(x)/dy

        type(T_Fdata_cell_3c), pointer :: pFdata_cell
        type(T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Set up Legendre polynomials
        if (.false.) cost2 = rhat(1)
        if (.false.) cost2 = sighat(1)

        cost2 = cost**2
        sint = sqrt(max(1.0d-5, 1.0d0 - cost2))
        p(0) = 1.0d0
        p(1) = cost
        p(2) = 0.5d0*(3.0d0*cost2 -1.0d0)
        p(3) = 0.5d0*cost*(5.0d0*cost2 - 3.0d0)
        p(4) = 0.125d0*(3.0d0 - cost2*(30.0d0 - 35.0d0*cost2))

        Dp(0) = 0.0d0
        Dp(1) = 1.0d0
        Dp(2) = 3.0d0*cost
        Dp(3) = (15.0d0*cost2 - 3.0d0)/2.0d0
        Dp(4) = (35.0d0*cost*cost2 - 15.0d0*cost)/2.0d0

        ! initialize, finding the correct fdata bundle
        pFdata_bundle => Fdata_bundle_3c(ispecies, jspecies, kspecies)
        pFdata_cell =>                                                       &
     &    pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(iint,isub,1))
        allocate (Fdata(pFdata_cell%nME))
        allocate (DpFdata(pFdata_cell%nME))
        allocate (DxFdata(pFdata_cell%nME))
        allocate (DyFdata(pFdata_cell%nME))
        Fdata = 0.0d0
        DpFdata = 0.0d0
        DxFdata = 0.0d0
        DyFdata = 0.0d0

        ! loop over theta
        do itheta = 1, P_maxtheta
          pFdata_cell =>                                                     &
     &      pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(iint,isub,itheta))
          call addDMEs_Fdata_3c (pFdata_cell, x, y, p(itheta - 1),           &
     &                           dp(itheta - 1), pFdata_cell%nME, DpFdata,   &
     &                           DxFdata, DyFdata, Fdata)
        end do

        ! recover matrices
        Dphmbox = 0.0d0
        Dxhmbox = 0.0d0
        Dyhmbox = 0.0d0

        do iindex = 1, pFdata_cell%nME
          imu = pFdata_cell%mu_3c(iindex)
          inu = pFdata_cell%nu_3c(iindex)
          mvalue = pFdata_cell%mvalue_3c(iindex)
          if (mvalue .ne. 1) then
            hmbox(imu,inu) = Fdata(iindex)
            Dphmbox(imu,inu) = DpFdata(iindex)
            Dxhmbox(imu,inu) = DxFdata(iindex)
            Dyhmbox(imu,inu) = DyFdata(iindex)
          else
            hmbox(imu,inu) = Fdata(iindex)*sint
            Dphmbox(imu,inu) = DpFdata(iindex)*sint - cost*Fdata(iindex)/sint
            Dxhmbox(imu,inu) = DxFdata(iindex)*sint
            Dyhmbox(imu,inu) = DyFdata(iindex)*sint
          end if
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (Fdata)
        deallocate (DpFdata)
        deallocate (DxFdata)
        deallocate (DyFdata)


! Format Statements
! ===========================================================================
        return
        end subroutine getDMEs_Fdata_3c


! addMEs_Fdata_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine uses interpolation to find the value of f(x) for any x,   !BaZ CHANGEME!!!
!! given an array of equally spaced points for f(x).  The function f(x) is a
!! three-center interaction function read from read_Fdata_3c.
!!
!! For polynomial interpolation see Mathews and Walker, p.329

! ===========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
!! @author Barry J. Haycock
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine addDMEs_Fdata_3c (pFdata_cell, x, y, ptheta, Dptheta, &
        &                            nME, Dpresult, Dxresult, Dyresult, result)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nME         !< number of matrix elements

        real, intent(in) :: x, y            !< distances from atom of interest
        real, intent(in) :: ptheta          !< p is the Legendre multiplier
        real, intent(in) :: Dptheta         !< Dp is the Legendre multiplier

        type(T_Fdata_cell_3c), pointer :: pFdata_cell !< FIXME

! Input and Output
        real, intent(inout), dimension (nME) :: Dxresult !< FIXME
        real, intent(inout), dimension (nME) :: Dyresult !< FIXME
        real, intent(inout), dimension (nME) :: Dpresult !< FIXME
        real, intent(inout), dimension (nME) :: result !< FIXME

       real, dimension (nME) :: interim !< FIXME

! Variable Declaration and Description
! ===========================================================================
        integer ix, iy
        integer k
        integer nME3c_max                   !< number of matrix elements

        real dx, dy
        real, allocatable, dimension (:) :: gx, gy, gg
        real px, py

! data arrays for interpolation
        real, allocatable, dimension (:) :: bb0, bb1, bb2, bb3
        real, allocatable, dimension (:) :: f1m1, f1m2, f1m3, f0p3, f0p6
        real, allocatable, dimension (:) :: f1p3, f1p6, f2p1
        real, allocatable, dimension (:, :) :: g, gp

! Local Parameters and Data Declaration
! ===========================================================================
! tolerance (may be needed to avoid roundoff error in the calling program)
! if xin > xmax but within, say, .001% of xmax then ignore
        real, parameter :: P_tiny = 1.0d-10
        real, parameter :: P_small = 1.0d0-8
        nME3c_max = pFdata_cell%nME

        ! Allocate the arrays used in Method 3
        allocate (bb0(nME3c_max), bb1(nME3c_max), bb2(nME3c_max))
        allocate (bb3(nME3c_max), f1m1(nME3c_max), f1m2(nME3c_max))
        allocate (f1m3(nME3c_max), f0p3(nME3c_max), f0p6(nME3c_max))
        allocate (f1p3(nME3c_max), f1p6(nME3c_max), f2p1(nME3c_max))
        allocate (g(-1:2,nME3c_max))
        allocate (gp(-1:2,nME3c_max))

! Procedure
! ===========================================================================
! Initialize
        dx = pFdata_cell%xmax/(pFdata_cell%nx - 1.0d0)
        dy = pFdata_cell%ymax/(pFdata_cell%ny - 1.0d0)

        ix = int(x/dx) +1
        iy = int(y/dy) +1

        if (ix .lt. 2) ix = 2
        if (iy .lt. 2) iy = 2

        if (ix .gt. pFdata_cell%nx - 2) ix = pFdata_cell%nx - 2
        if (iy .gt. pFdata_cell%ny - 2) iy = pFdata_cell%ny - 2

        px = x/dx - ix + 1.0d0
        py = y/dy - iy + 1.0d0

! ***************************************************************************
! Adaptive interpolation - estimate gradient:           !Baz, think this should be an array, methinks. Yes Fdata_3c is (ix, iy, iME)
        allocate (gx (nME)); gx = 0.0d0
        allocate (gy (nME)); gy = 0.0d0
        allocate (gg (nME)); gg = 0.0d0

        gx = (pFdata_cell%Fdata_3c(ix + 1, iy, :)                             &
     &        - pFdata_cell%Fdata_3c(ix, iy, :))/dx
        gy = (pFdata_cell%Fdata_3c(ix, iy + 1, :)                             &
     &        - pFdata_cell%Fdata_3c(ix, iy, :))/dy
        gg = gx**2 + gy**2

! We choose one of three ways to interpolate based on the size of the
! gradient (derivative). If the value of the derivative is tiny, then
! do very simple interpolation. If small, then something a little more
! complicated. If neither, then something even MORE complicated.

! METHOD 1
! Do three point linear bivariate interpolation
! Handbook of Mathematical Functions..., edited by M. Abramowitz
! and I.A. Stegun, Dover edition, Pg. 882, Eq. 25.2.65
        where (gg .lt. P_tiny)
         interim = ((1.0d0 - px - py)*pFdata_cell%Fdata_3c(ix,iy,:)  &
     &                                    + px*pFdata_cell%Fdata_3c(ix+1,iy,:)&
     &                                    + py*pFdata_cell%Fdata_3c(ix,iy+1,:))
         result = result + ptheta*interim
         Dpresult = Dpresult + Dptheta*interim
         Dxresult = Dxresult + ptheta*gx
         Dyresult = Dyresult + ptheta*gy


! METHOD 2
! Do quadratic bivariate interpolation (six point formula, Eq. 25.2.67)
        elsewhere (gg .lt. P_small)
          interim = (py*(py-1)*0.5d0*pFdata_cell%Fdata_3c(ix,iy-1,:) &
     &                      + px*(px-1)*0.5d0*pFdata_cell%Fdata_3c(ix-1,iy,:) &
     &        + (1.0d0 + px*py - px**2 - py**2)*pFdata_cell%Fdata_3c(ix,iy,:) &
     &      + px*(px - 2.0d0*py + 1.0d0)*0.5d0*pFdata_cell%Fdata_3c(ix+1,iy,:)&
     &      + py*(py - 2.0d0*px + 1.0d0)*0.5d0*pFdata_cell%Fdata_3c(ix,iy+1,:)&
     &                      + px*py*pFdata_cell%Fdata_3c(ix+1,iy+1,:))

          result = ptheta*interim
          Dpresult = Dpresult + Dptheta*interim

          Dxresult = Dxresult + ptheta                                                 &
     &    *(((pFdata_cell%Fdata_3c(ix+1,iy+1,:) - pFdata_cell%Fdata_3c(ix+1,iy,:)      &
     &         - pFdata_cell%Fdata_3c(ix,iy+1,:) + pFdata_cell%Fdata_3c(ix,iy,:))*py   &
     &         + (pFdata_cell%Fdata_3c(ix-1,iy,:) + pFdata_cell%Fdata_3c(ix+1,iy,:)    &
     &         - 2.0d0*pFdata_cell%Fdata_3c(ix,iy,:))*px - 0.5d0*(pFdata_cell%Fdata_3c(ix-1,iy,:)               &
     &         - pFdata_cell%Fdata_3c(ix+1,iy,:)))/dx)

          Dyresult = Dyresult + ptheta                                                   &
           *(((pFdata_cell%Fdata_3c(ix+1,iy+1,:) - pFdata_cell%Fdata_3c(ix+1,iy,:)       &
     &      - pFdata_cell%Fdata_3c(ix,iy+1,:) + pFdata_cell%Fdata_3c(ix,iy,:))*px       &
     &            + (pFdata_cell%Fdata_3c(ix,iy-1,:) + pFdata_cell%Fdata_3c(ix,iy+1,:) &
     &            - 2.0d0*pFdata_cell%Fdata_3c(ix,iy,:))*py                              &
     &            - 0.5d0*(pFdata_cell%Fdata_3c(ix,iy-1,:) - pFdata_cell%Fdata_3c(ix,iy+1,:)))/dy)

!          dQ_Ldx = ((fun(1,1) - fun(1,0) - fun(0,1) + fun(0,0))*py           &
!     &            + (fun(-1,0) + fun(1,0) - 2.0d0*fun(0,0))*px               &
!     &            - 0.5d0*(fun(-1,0) - fun(1,0)))/hx
!          dQ_Ldy = ((fun(1,1) - fun(1,0) - fun(0,1) + fun(0,0))*px           &
!     &            + (fun(0,-1) + fun(0,1) - 2.0d0*fun(0,0))*py               &
!     &            - 0.5d0*(fun(0,-1) - fun(0,1)))/hy
        end where


! METHOD 3
! Interpolate one direction, then interpolate using these values to get
! the final value.  Use Eq. 25.2.13 (4 point intepolater).
! Total number of points used is 16
          nME3c_max = pFdata_cell%nME

          do k = -1, 2
            where (gg .gt. P_small)
            f1m1 = pFdata_cell%Fdata_3c(ix+k,iy-1,:)
            f1m2 = 2.0d0*f1m1
            f1m3 = 3.0d0*f1m1

            f0p3 = 3.0d0*pFdata_cell%Fdata_3c(ix+k,iy,:)
            f0p6 = 2.0d0*f0p3

            f1p3 = 3.0d0*pFdata_cell%Fdata_3c(ix+k,iy+1,:)
            f1p6 = 2.0d0*f1p3

            f2p1 = pFdata_cell%Fdata_3c(ix+k,iy+2,:)

            bb3 = - f1m1 + f0p3 - f1p3 + f2p1
            bb2 = f1m3 - f0p6 + f1p3
            bb1 = - f1m2 - f0p3 + f1p6 - f2p1
            bb0 = f0p6

            g(k,:) = ((bb3*py + bb2)*py + bb1)*py + bb0
            gp(k,:) = ((3*bb3*py + 2*bb2)*py + bb1)
            end where
          end do

          where (gg .gt. P_small)

          f1m1 = g(-1,:)
          f1m2 = 2.0d0*f1m1
          f1m3 = 3.0d0*f1m1

          f0p3 = 3.0d0*g(0,:)
          f0p6 = 2.0d0*f0p3

          f1p3 = 3.0d0*g(1,:)
          f1p6 = 2.0d0*f1p3

          f2p1 = g(2,:)

          bb3 = - f1m1 + f0p3 - f1p3 + f2p1
          bb2 = f1m3 - f0p6 + f1p3
          bb1 = - f1m2 - f0p3 + f1p6 - f2p1
          bb0 = f0p6
          interim = (((bb3*px + bb2)*px + bb1)*px + bb0)/36.0d0
          result = result + ptheta*interim
          Dpresult = Dpresult + Dptheta*interim

          Dxresult = Dxresult + ptheta*((3*bb3*px + 2*bb2)*px + bb1)/(36.0d0*dx)
          f1m1 = gp(-1,:)
          f1m2 = 2*f1m1
          f1m3 = 3*f1m1

          f0p3 = 3*gp(0,:)
          f0p6 = 2*f0p3

          f1p3 = 3*gp(1, :)
          f1p6 = 2*f1p3

          f2p1 = gp(2,:)

          bb3 = - f1m1 + f0p3 - f1p3 + f2p1
          bb2 = f1m3 - f0p6 + f1p3
          bb1 = - f1m2 - f0p3 + f1p6 - f2p1
          bb0 = f0p6

          Dyresult = Dyresult + ptheta*(((bb3*px + bb2)*px + bb1)*px + bb0)/(36.0d0*dy)

          end where

!           dQ_Ldx = ((3*bb3*px + 2*bb2)*px + bb1)/(36.0d0*hx)
!
!           f1m1 = gp(-1)
!           f1m2 = 2*f1m1
!           f1m3 = 3*f1m1
!
!           f0p3 = 3*gp(0)
!           f0p6 = 2*f0p3
!
!           f1p3 = 3*gp(1)
!           f1p6 = 2*f1p3
!
!           f2p1 = gp(2)
!
!           bb3 = - f1m1 + f0p3 - f1p3 + f2p1
!           bb2 = f1m3 - f0p6 + f1p3
!           bb1 = - f1m2 - f0p3 + f1p6 - f2p1
!           bb0 = f0p6
!
!           dQ_Ldy = (((bb3*px + bb2)*px + bb1)*px + bb0)/(36.0d0*hy)


 ! Deallocate the arrays used in Method 3.
          deallocate (bb0, bb1, bb2, bb3)
          deallocate (f1m1, f1m2, f1m3, f0p3, f0p6)
          deallocate (f1p3, f1p6, f2p1, g)

! A final note, if you are interested in splines, you should start with
! "Handbook on SPLINES for the User" by Eugene V. Shikin and
! Alexander I. Plis.  1995, CRC Press.  Most other books on slines and
! bivariate interpolation are nothing but proofs and abstract math, but
! this one gives the real equations you need to get the actual work done.

! Deallocate Arrays
! ===========================================================================
! None
        deallocate (gx, gy, gg)

! Format Statements
! ===========================================================================
! None

        return
        end subroutine addDMEs_Fdata_3c


! ===========================================================================
! destroy_Fdata_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the three-center interactions - these arrays
!! are read in by read2.
!
! ===========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
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
        subroutine destroy_Fdata_3c ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies, kspecies           !< counter over the species
        integer icell                        !< counter over Fdata files

! Procedure
! ===========================================================================
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            do kspecies =1, nspecies
              do icell = 1, Fdata_bundle_3c(ispecies,jspecies,kspecies)%nFdata_cell_3c
                deallocate(Fdata_bundle_3c(ispecies,jspecies,kspecies)%Fdata_cell_3c(icell)%mu_3c)
                deallocate(Fdata_bundle_3c(ispecies,jspecies,kspecies)%Fdata_cell_3c(icell)%nu_3c)
                deallocate(Fdata_bundle_3c(ispecies,jspecies,kspecies)%Fdata_cell_3c(icell)%mvalue_3c)
                deallocate(Fdata_bundle_3c(ispecies,jspecies,kspecies)%Fdata_cell_3c(icell)%Fdata_3c)
              end do
              deallocate(Fdata_bundle_3c(ispecies,jspecies,kspecies)%Fdata_cell_3c)
            end do
          end do
        end do
        deallocate (Fdata_bundle_3c)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_Fdata_3c


! End Module
! ===========================================================================
        end module M_Fdata_3c
