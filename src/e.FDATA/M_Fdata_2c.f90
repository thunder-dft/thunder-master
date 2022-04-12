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

! M_Fdata_2c
! Module Description
! ===========================================================================
!>       This is a module containing all of the subroutines that will read in
!! all the two-center fdata from the fdata files in the Fdata directory.  It
!! contains the following subroutines within the module:
!!
!!       read_Fdata_2c.f90 - reads in the two-center datafiles.
!!       getMEs_Fdata_2c.f90 - interpolates the two-center data and gets the
!!                             matrix elements.
!!       getDMEs_Fdata_2c.f90 - interpolates the two-center derivatives and
!!                              gets derivative matrix elements
!!       destroy_Fdata_2c.f90 - deallocates all the two-center fdata arrays.
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
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
        module M_Fdata_2c
        use M_species

        include '../include/interactions_2c.h'

! Type Declaration
! ===========================================================================
! two-center interactions arrays

! To cut down on storage space, we actually change the storage procedure
! from previous FIREBALL code. Not all atoms have the same number of
! interactions or interaction types. Before - we would store things based
! on the maximum number of fdata points, maximum number of interactions types,
! maximum number of matrix elements - so even hydrogen-hydrogen (just ss
! and/or ss*) stored a 4x4 or an 8x8 matrix even when not needed.  This was
! quite inefficient.

! The new approach is to define some Fdata types which store the actual
! Fdata points. The smallest unit storage is called Fdata_cell_2c, containing
! all Fdata for a particular interaction/subinteraction type.
        type T_Fdata_cell_2c
          integer nME                   !< number of non-zero matrix elements
          integer nx                    !< number of data points

          integer, pointer :: mu_2c (:)
          integer, pointer :: nu_2c (:)
          integer, pointer :: mvalue_2c (:)

          real dx                          !< distance between data points
          real xmax                        !< maximum interaction range

          ! actual Fdata points f(x)
          real, pointer :: Fdata_2c (:,:)
        end type T_Fdata_cell_2c

! Fdata_bundle_2c is the 2-center package of Fdata_cell_2c. It contains all the
! Fdata_cell_2c information of all interaction/subtypes for given species pair.
        type T_Fdata_bundle_2c
          integer nFdata_cell_2c                     !< number of Fdata_cell2C

          integer, dimension (P_maxtype, 0:P_maxsubtype) :: index_2c

          ! actual fdata
          type (T_Fdata_cell_2c), pointer :: Fdata_cell_2c (:)
        end type T_Fdata_bundle_2c

        type(T_Fdata_bundle_2c), pointer :: Fdata_bundle_2c (:,:)

! module procedures
        contains

! ===========================================================================
! read_Fdata_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in all two-center interactions. These two-center
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
        subroutine read_Fdata_2c ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies       !< counters for number of species
        integer issh                     !< counter for number of shells
        integer icell                    !< counter for Fdata files
        integer iME                      !< counter for number of MEs
        integer itype, isubtype          !< counter over types, subtypes
        integer iindex                   !< combined type, subtype to one index
        integer ix                       !< counters over distances

        integer nME2c_max                !< temporary value for nME

        type(T_Fdata_cell_2c), pointer :: pFdata_cell
        type(T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 32) filename

! Allocate Arrays
! ===========================================================================
        allocate (Fdata_bundle_2c (nspecies, nspecies))

! Procedure
! ===========================================================================
        write (ilogfile, '(A)') 'Reading Two-center interactions '
        write (ilogfile, '(A)') '------------------------------- '
        write (ilogfile, *)

        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            write (ilogfile,'(4x, A15, I4, I4)') '- species pair: ', species(ispecies)%nZ, species(jspecies)%nZ
            ! cut some lengthy notation
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

            ! get the number of interactions
            write (filename,'("/2c.",i2.2,".",i2.2,".dat")')                 &
     &        species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'old')
            read (11,*) pFdata_bundle%nFdata_cell_2c
            allocate (pFdata_bundle%Fdata_cell_2c(pFdata_bundle%nFdata_cell_2c))
            pFdata_bundle%index_2c = -1

            ! open directory file
            write (filename,'("/2c.",i2.2,".",i2.2,".dir")')                 &
     &        species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'old')
            close (11)

            ! loop over Fdata files
            do icell = 1, pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(icell)
              read (12,*) iindex, itype, isubtype, filename,                 &
      &                   pFdata_cell%nME, pFdata_cell%nx, pFdata_cell%xmax
              nME2c_max = pFdata_cell%nME
              allocate (pFdata_cell%Fdata_2c(pFdata_cell%nx,nME2c_max))
              allocate (pFdata_cell%mu_2c(nME2c_max))
              allocate (pFdata_cell%nu_2c(nME2c_max))
              allocate (pFdata_cell%mvalue_2c(nME2c_max))

              ! make the index file
              pFdata_bundle%index_2c(itype,isubtype) = iindex
              pFdata_cell%dx = pFdata_cell%xmax/(pFdata_cell%nx - 1)

              ! open the actual datafile
              open (unit = 13,                                               &
     &              file = trim(Fdata_location)//'/'//trim(filename),        &
     &              status = 'old')

              ! read the actual fdata
              do ix = 1, pFdata_cell%nx    ! loop over bondlengths
                read (13,*) (pFdata_cell%Fdata_2c(ix,iME), iME = 1, nME2c_max)
              end do
              close (13)

              ! open mu, nu, mvalue file
              write (filename,'("/", i2.2, "_munu_2c.", i2.2, ".", i2.2,     &
     &                          ".dat")')                                    &
     &          itype, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 14, file = trim(Fdata_location)//trim(filename),  &
     &          status = 'old')

              ! read the mapping - stored in mu, nu, and mvalue
              read (14,*) (pFdata_cell%mu_2c(iindex), iindex = 1, nME2c_max)
              read (14,*) (pFdata_cell%nu_2c(iindex), iindex = 1, nME2c_max)
              read (14,*) (pFdata_cell%mvalue_2c(iindex), iindex = 1, nME2c_max)
            end do ! datafile for (ispecies, jspecies)
            close (12)
            close (14)
          end do

! Get the Kleinman-Bylander coefficients
          write (filename,'("/clPP.",i2.2,".dat")') species(ispecies)%nZ
          open (unit = 14,                                                 &
     &      file = trim(Fdata_location)//'/'//trim(filename), status = 'old')
          do issh = 1, species(ispecies)%nssh_PP
            read (14,*) species(ispecies)%shell_PP(issh)%cl
          end do
          close(14)
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
        end subroutine read_Fdata_2c


! ===========================================================================
! getMEs_Fdata_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine uses interpolation to find the value of f(x) for any x,
!! given an array of equally spaced points for f(x).  The function f(x) is a
!! two-center interaction function read from read_Fdata_2c. After the data
!! is interpolated, the data is 'recovered' into the matrix form.
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
!
! Subroutine Declaration
! ===========================================================================
        subroutine getMEs_Fdata_2c (ispecies, jspecies, iint, isub, x,     &
     &                              norb_mu, norb_nu, hmbox)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: ispecies, jspecies    !< species
        integer, intent(in) :: iint, isub            !< integral type, subtype
        integer, intent(in) :: norb_mu, norb_nu      !< Index max orbital for mu and nu atoms

        real, intent(in) :: x                        !< distance between pair

! Output
        real, intent(out), dimension (norb_mu, norb_nu) :: hmbox !< Fdata indexed for mu and nu

! Local Parameters and Data Declaration
! ===========================================================================
!> @param P_tolerance (may be needed to avoid roundoff error in the calling
!! program) if xin > xmax but within, say, .001% of xmax then ignore
        real, parameter :: P_tolerance = 1.0d-5

! Variable Declaration and Description
! ===========================================================================
        integer iindex, imu, inu               !< counters for building matrix
        integer ipoint, ileft, imid, iright    !< points along the grid
        integer iprod, isum                    !< product, sum points on grid

        real pden, prod                        !< stuff for products

        real, dimension (:), allocatable :: Fdata      !< F(x)
        real, dimension (0:6) :: xx, pdenom

        type(T_Fdata_cell_2c), pointer :: pFdata_cell
        type(T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
        pFdata_bundle => Fdata_bundle_2c(ispecies,jspecies)
        pFdata_cell =>                                                       &
      &   pFdata_bundle%Fdata_cell_2c(pFdata_bundle%index_2c(iint,isub))
        allocate (Fdata(pFdata_cell%nME))

! The following should never happen.  Delete them might improve performance
        if (x .lt. 0.0d0) then
          stop ' Error in getMEs_Fdata_2c: negative x! '
        else if (x .gt. pFdata_cell%xmax*(1.0d0 + P_tolerance)) then
          write (*,*) ' interaction, subtype = ', iint, isub
          write (*,*) ' ispecies, jspecies = ', ispecies, jspecies
          write (*,*) ' x = ', x, ' xmax = ', pFdata_cell%xmax
          stop ' Error in getMEs_Fdata_2C: x too large! '
        else if (x .lt. P_tolerance) then
          ! Now, construct the matrix
          hmbox = 0.0d0
          do iindex = 1, pFdata_cell%nME
            imu = pFdata_cell%mu_2c(iindex)
            inu = pFdata_cell%nu_2c(iindex)
            hmbox(imu,inu) = pFdata_cell%Fdata_2c(1,iindex)
          end do
        end if

! now find starting and ending points for the interpolation
! note : imid is the point to the left (code assumes this)
        imid = int(x/pFdata_cell%dx) + 1
        ileft = imid - 3
        iright = imid + 3
        if (ileft .lt. 1) then
          ileft = 1
          iright = 7
        else if (iright .gt. pFdata_cell%nx) then
          ileft = pFdata_cell%nx - 6
          iright = pFdata_cell%nx
        end if

! Interpolate with polynomials of order 5
        do ipoint = ileft, iright
          xx(ipoint - ileft) = (ipoint - 1)*pFdata_cell%dx
        end do
        Fdata = 0.0d0
        do isum = ileft, iright
          prod = 1.0d0
          pdenom(isum - ileft) = 1.0d0
          do iprod = ileft, iright
            if (iprod .ne. isum) then
              pden = 1.0d0/(xx(isum - ileft) - xx(iprod - ileft))
              pdenom(isum - ileft) = pdenom(isum - ileft)*pden
              prod = prod*(x - xx(iprod - ileft))*pden
            end if
          end do
          Fdata = Fdata + pFdata_cell%Fdata_2c(isum,:)*prod
        end do

! Now, construct the matrix
        hmbox = 0.0d0
        do iindex = 1, pFdata_cell%nME
          imu = pFdata_cell%mu_2c(iindex)
          inu = pFdata_cell%nu_2c(iindex)
          hmbox(imu,inu) = Fdata(iindex)
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
        end subroutine getMEs_Fdata_2c


! ===========================================================================
! getDMEs_Fdata_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine uses interpolation to find the value of f(x) for any x,
!! given an array of equally spaced points for f(x).  The function f(x) is a
!! two-center interaction function read from read_Fdata_2c. After the data
!! is interpolated, the data is 'recovered' into the matrix form.
!!
!! For polynomial interpolation see Mathews and Walker, p.329

! ===========================================================================
! Code written by:
!> @author Barry Haycock
!! @author Ning Ma
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
        subroutine getDMEs_Fdata_2c (ispecies, jspecies, iint, isub, x,       &
     &                               norb_mu, norb_nu, hmbox, Dhmbox)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: ispecies, jspecies    !< species
        integer, intent(in) :: iint, isub            !< integral type, subtype
        integer, intent(in) :: norb_mu, norb_nu      !< Index max orbital for mu and nu atoms

        real, intent(in) :: x                        !< distance between pair

! Output
        ! Fdata indexed for mu and nu
        real, intent(out), dimension (norb_mu, norb_nu) :: hmbox

        ! Derivative of Fdata indexed for mu and nu
        real, intent(out), dimension (norb_mu, norb_nu) :: Dhmbox

! Local Parameters and Data Declaration
! ===========================================================================
!> @param P_tolerance (may be needed to avoid roundoff error in the calling
!! program) if xin > xmax but within, say, .001% of xmax then ignore
        real, parameter :: P_tolerance = 1.0d-5

! Variable Declaration and Description
! ===========================================================================
        integer iindex, imu, inu               !< counters for building matrix
        integer ipoint, ileft, imid, iright    !< points along the grid
        integer iprod, isum                    !< product, sum points on grid

        real pden, prod                        !< stuff for products

        real, dimension (:), allocatable :: Fdata     !< F(x)
        real, dimension (:), allocatable :: dFdatadx  !< dF(x)
        real, dimension (0:6) :: xx, pdenom

        type(T_Fdata_cell_2c), pointer :: pFdata_cell
        type(T_Fdata_bundle_2c), pointer :: pFdata_bundle

		integer kpoint, jpoint
		real xsumoverj, xprod

! Procedure
! ===========================================================================
        pFdata_bundle => Fdata_bundle_2c(ispecies,jspecies)
        pFdata_cell =>                                                       &
      &   pFdata_bundle%Fdata_cell_2c(pFdata_bundle%index_2c(iint,isub))
        allocate (dFdatadx(pFdata_cell%nME))
        allocate (Fdata(pFdata_cell%nME))
		dFdatadx = 0.0d0
		Fdata = 0.0d0
		xx = 0.00
		pdenom = 0.00

! The following should never happen.  Delete them might improve performance
        if (x .lt. 0.0d0) then
          stop ' Error in getDMEs_Fdata_2c: negative x! '
        else if (x .gt. pFdata_cell%xmax*(1.0d0 + P_tolerance)) then
          write (*,*) ' interaction, subtype = ', iint, isub
          write (*,*) ' ispecies, jspecies = ', ispecies, jspecies
          write (*,*) ' x = ', x, ' xmax = ', pFdata_cell%xmax
          stop ' Error in getDMEs_Fdata_2C: x too large! '
        else if (x .lt. P_tolerance) then
          ! Now, construct the matrix
          hmbox = 0.0d0
          do iindex = 1, pFdata_cell%nME
            imu = pFdata_cell%mu_2c(iindex)
            inu = pFdata_cell%nu_2c(iindex)
            hmbox(imu,inu) = pFdata_cell%Fdata_2c(1,iindex)
            Dhmbox(imu,inu) = 0.0d0
          end do          
          return
        end if

! now find starting and ending points for the interpolation
! note : imid is the point to the left (code assumes this)
        imid = int(x/pFdata_cell%dx) + 1
        ileft = imid - 3
        iright = imid + 3
        if (ileft .lt. 1) then
          ileft = 1
          iright = 7
        else if (iright .gt. pFdata_cell%nx) then
          ileft = pFdata_cell%nx - 6
          iright = pFdata_cell%nx
        end if

! Find derivatives
        do ipoint = ileft, iright
          xx(ipoint - ileft) = (ipoint - 1)*pFdata_cell%dx
        end do

        do isum = ileft, iright
          prod = 1.0d0
          pdenom(isum - ileft) = 1.0d0
          do iprod = ileft, iright
            if (iprod .ne. isum) then
              pden = 1.0d0/(xx(isum - ileft) - xx(iprod - ileft))
              pdenom(isum - ileft) = pdenom(isum - ileft)*pden
              prod = prod*(x - xx(iprod - ileft))*pden
            end if
          end do
          Fdata = Fdata + pFdata_cell%Fdata_2c(isum,:)*prod
        end do

        ! Now, construct the matrix
        hmbox = 0.0d0
        do iindex = 1, pFdata_cell%nME
          imu = pFdata_cell%mu_2c(iindex)
          inu = pFdata_cell%nu_2c(iindex)
          hmbox(imu,inu) = Fdata(iindex)
        end do

        dFdatadx = 0.0d0
        do ipoint = ileft, iright
          xsumoverj = 0.0d0
          do jpoint = ileft, iright
            if (ipoint .eq. jpoint) cycle
            xprod = 1.0d0
            do kpoint = ileft, iright
              if (kpoint .eq. jpoint) cycle
              if (kpoint .eq. ipoint) cycle
              xprod = xprod * (x - xx(kpoint - ileft))
            end do ! kpoint
            xsumoverj = xsumoverj + xprod
          end do ! jpoint

          dFdatadx = dFdatadx + pFdata_cell%Fdata_2c(ipoint,:)               &
     &                           *pdenom(ipoint - ileft)*xsumoverj
        end do ! ipoint

! Now, construct the matrix
        Dhmbox = 0.0d0
        do iindex = 1, pFdata_cell%nME
          imu = pFdata_cell%mu_2c(iindex)
          inu = pFdata_cell%nu_2c(iindex)
          Dhmbox(imu,inu) = dFdatadx(iindex)
        end do

! Deallocate Arrays
! ===========================================================================
		deallocate (dFdatadx)
		deallocate (Fdata)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine getDMEs_Fdata_2c


! ===========================================================================
! destroy_Fdata_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the two-center interactions - these arrays
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
        subroutine destroy_Fdata_2c ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies           !< counter over the species
        integer icell                        !< counter over Fdata files

! Procedure
! ===========================================================================
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            do icell = 1, Fdata_bundle_2c(ispecies,jspecies)%nFdata_cell_2c
              deallocate (Fdata_bundle_2c(ispecies,jspecies)%Fdata_cell_2c(icell)%mu_2c)
              deallocate (Fdata_bundle_2c(ispecies,jspecies)%Fdata_cell_2c(icell)%nu_2c)
              deallocate (Fdata_bundle_2c(ispecies,jspecies)%Fdata_cell_2c(icell)%mvalue_2c)
              deallocate (Fdata_bundle_2c(ispecies,jspecies)%Fdata_cell_2c(icell)%Fdata_2c)
            end do
            deallocate (Fdata_bundle_2c(ispecies,jspecies)%Fdata_cell_2c)
          end do
        end do
        deallocate (Fdata_bundle_2c)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_Fdata_2c

! End Module
! ===========================================================================
        end module M_Fdata_2c

