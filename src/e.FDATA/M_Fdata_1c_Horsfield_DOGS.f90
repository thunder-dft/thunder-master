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

! M_Fdata_1c
! Module Description
! ===========================================================================
!>       This is a module containing all of the subroutines that will read in
!! all the data from the data files in the Fdata directory.  It contains the
!! following subroutines within the module:
!!
!!       read_xc1c.f90 - read in data from one-center xc McWEDA datafiles.
!!       destroy_xc1c.f90 - destroy allocatable arrays related to the one-center
!!                          McWeda interactions.
!
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
! Module Declaration
! ===========================================================================
        module M_Fdata_1c
        use M_species

! Type Declaration
! ===========================================================================
! one-center xc data for McWEDA
        type T_vxc_1c
          real E                               !< xc-energy (shells)

          real, pointer :: dE(:)               !< 1st derivatives of energy

          real, pointer :: V(:, :)             !< Vxc potential
          real, pointer :: dV(:, :, :)         !< 1st derivatives of Vxc
        end type T_vxc_1c

        type T_gover_1c_L
          real, pointer :: f1nac(:, :) => null()   !< Vxc potential
          real, pointer :: f2nac(:, :) => null()   !< Vxc potential
        end type T_gover_1c_L

        type T_gover_1c_R
          real, pointer :: f1nac(:, :) => null()   !< Vxc potential
          real, pointer :: f2nac(:, :) => null()   !< Vxc potential
        end type T_gover_1c_R

! module variables
        type (T_vxc_1c), pointer :: vxc_1c (:) !< one-center xc for McWEDA
        type (T_gover_1c_L), pointer :: gover_1c_L (:) !< one-center xc for McWEDA
        type (T_gover_1c_R), pointer :: gover_1c_R (:) !< one-center xc for McWEDA

! module procedures
        contains

! ===========================================================================
! read_Fdata_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in the one-center (exchange-correlation)
!! interactions for McWEDA. These one-center interactions are contributions
!! as described in
!!
!! "Multicenter approach to the exchange-correlation interactions in ab initio
!!  tight-binding methods" by P. Jelinek, H. Wang, J.P. Lewis, O.F. Sankey,
!!  and J. Ortega, PRB 71:23511 (2005).
!!
!! This routine also reads in the variables which are needed to compute
!! changes of the exchange correlation for the charge transfer which is
!! contained in the one-center datafiles.
!
! ===========================================================================
! Code written by:
!> @author Daniel G. Trabada
!! @author J. Ortega
!! @author James P. Lewis\n
!! Box 6315, 209 Hodges Hall\n
!! Department of Physics\n
!! West Virginia University\n
!! Morgantown, WV 26506-6315\n
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine read_Fdata_1c ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies                  !< counter for number of species

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        integer nssh                      !< counters for number of shells
        integer issh, jssh
        logical gover_exists

        character (len=30) filename

! Allocate Arrays
! ===========================================================================
        allocate (vxc_1c (nspecies))
        allocate (gover_1c_L (nspecies))
        allocate (gover_1c_R (nspecies))

! Procedure
! ===========================================================================
! ***************************************************************************
!                 R E A D    M A T R I X    E L E M E N T S
! ***************************************************************************
        write (ilogfile, '(A)') 'Reading One-center interactions '
        write (ilogfile, '(A)') '------------------------------- '
        write (ilogfile, *)

        do ispecies = 1, nspecies
          write (ilogfile,'(4x, A15, I4)') '- species: ', species(ispecies)%nZ

          ! Open ouput file for this species pair
          write (filename, '("/vxc_1c", ".", i2.2, ".dat")') species(ispecies)%nZ
          open (11, file = trim(fdata_location)//trim(filename), status = 'old')

          nssh = species(ispecies)%nssh
          allocate (vxc_1c(ispecies)%V(nssh,nssh))

          ! 0th order
          read (11,*) vxc_1c(ispecies)%E
          do issh = 1, nssh
            read (11,*) (vxc_1c(ispecies)%V(issh,jssh), jssh = 1, nssh)
          end do
          read (11,*)

! ***************************************************************************
!                    R E A D    D E R I V A T I V E   P I E C E S
! ***************************************************************************
! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
          ideriv_min = 1
          ideriv_max = 2

          ! allocate the derivative pieces
          allocate (vxc_1c(ispecies)%dE(ideriv_min:ideriv_max))
          allocate (vxc_1c(ispecies)%dV(ideriv_min:ideriv_max,nssh,nssh))

          do ideriv = ideriv_min, ideriv_max

            ! 1st order derivative terms
            read (11,*) vxc_1c(ispecies)%dE(ideriv)
            do issh = 1, nssh
              read (11,*) (vxc_1c(ispecies)%dV(ideriv,issh,jssh), jssh = 1, nssh)
            end do
            read (11,*)
          end do  ! end loop over ideriv
          close (11)

        ! goverlap_L
          write (filename, '("/goverlap_L.",i2.2,"answer1",".dat")') species(ispecies)%nZ
          inquire (file = trim(fdata_location)//trim(filename), exist = gover_exists)
          if (gover_exists) then
            open (12, file = trim(fdata_location)//trim(filename), status = 'old')

            allocate (gover_1c_L(ispecies)%f1nac(nssh,nssh))

            do issh = 1, nssh
              read (12,*) (gover_1c_L(ispecies)%f1nac(issh,jssh), jssh = 1, nssh)
            end do
            close (12)
          end if

          write (filename, '("/goverlap_L.",i2.2,"answer2",".dat")') species(ispecies)%nZ
          inquire (file = trim(fdata_location)//trim(filename), exist = gover_exists)
          if (gover_exists) then
            open (13, file = trim(fdata_location)//trim(filename), status = 'old')

            allocate (gover_1c_L(ispecies)%f2nac(nssh,nssh))

            do issh = 1, nssh
              read (13,*) (gover_1c_L(ispecies)%f2nac(issh,jssh), jssh = 1, nssh)
            end do   
            close (13)
          end if

        ! goverlap_R
          write (filename, '("/goverlap_R.",i2.2,"answer1",".dat")') species(ispecies)%nZ
          inquire (file = trim(fdata_location)//trim(filename), exist = gover_exists)
          if (gover_exists) then
            open (14, file = trim(fdata_location)//trim(filename), status = 'old')

            allocate (gover_1c_R(ispecies)%f1nac(nssh,nssh))

            do issh = 1, nssh
              read (14,*) (gover_1c_R(ispecies)%f1nac(issh,jssh), jssh = 1, nssh)
            end do
            close (14)
          end if

          write (filename, '("/goverlap_R.",i2.2,"answer2",".dat")') species(ispecies)%nZ
          inquire (file = trim(fdata_location)//trim(filename), exist = gover_exists)
          if (gover_exists) then
            open (15, file = trim(fdata_location)//trim(filename), status = 'old')

            allocate (gover_1c_R(ispecies)%f2nac(nssh,nssh))

            do issh = 1, nssh
              read (15,*) (gover_1c_R(ispecies)%f2nac(issh,jssh), jssh = 1, nssh)
            end do  
            close (15)
          end if
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
        end subroutine read_Fdata_1c

! ===========================================================================
! get_gover_1c
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
        subroutine get_gover_1c (ispecies, iint, norb_mu, hmbox)
        implicit none
        include '../include/interactions_2c.h'     

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: ispecies              !< species
        integer, intent(in) :: iint                  !< integral type, subtype
        integer, intent(in) :: norb_mu               !< Index max orbital for mu and nu atoms

! Output
        real, intent(out), dimension (3, norb_mu, norb_mu) :: hmbox !< Fdata indexed for mu and nu

! Parameters and Data Declaration
! ===========================================================================
        real, parameter :: ur3 = 1.0d0/dsqrt(3.0d0)
        real, parameter :: ur5 = 1.0d0/dsqrt(5.0d0)
        real, parameter :: ur15 = 1.0d0/dsqrt(15.0d0)

! Variable Declaration and Description
! ===========================================================================
        integer issh, jssh
        integer l1, l2          ! Counter of subshells
        integer n1, n2          ! Counter of shells

        real f1, f2

! Procedure
! ===========================================================================
! JOM-warning : there is an overall minus (-) sign
! due to Grad wrt R (atomic position) = - Grad wrt r (electron position)
! and another overall minus (-) sign
! due to  < G mu | nu >  = - < mu | G nu >


! Loop over shells i-atom        
        n1 = 0
        do issh = 1, species(ispecies)%nssh

! Number of orbitals per the shell in x-dim
         l1 = species(ispecies)%shell(issh)%lssh
         n2 = 0

! Loop over shells i y-dimension
         do jssh = 1, species(ispecies)%nssh

! Number of orbitals per the shell in y-dim
          l2 = species(ispecies)%shell(jssh)%lssh
! The F1 and F2 radial integrals
          if (iint .eq. P_goverlapL) then
            f1 = gover_1c_L(ispecies)%f1nac(issh,jssh)
            f2 = gover_1c_L(ispecies)%f2nac(issh,jssh)
          else if (iint .eq. P_goverlapR) then 
            f1 = gover_1c_R(ispecies)%f1nac(issh,jssh)
            f2 = gover_1c_R(ispecies)%f2nac(issh,jssh)
          end if
          if (l1.eq.0 .and. l2.eq.1) then
            hmbox(1, n1 + 1, n2 + 3) = ur3*(f1 + 2.0d0*f2)
            hmbox(2, n1 + 1, n2 + 1) = ur3*(f1 + 2.0d0*f2)
            hmbox(3, n1 + 1, n2 + 2) = ur3*(f1 + 2.0d0*f2)
          else if (l1.eq.1 .and. l2.eq.0) then
            hmbox(1, n1 + 3, n2 + 1) = ur3*f1
            hmbox(2, n1 + 1, n2 + 1) = ur3*f1
            hmbox(3, n1 + 2, n2 + 1) = ur3*f1
          else if (l1.eq.1 .and. l2.eq.2) then
            hmbox(1, n1 + 1, n2 + 1) = ur5*(f1 + 3.0d0*f2)
            hmbox(1, n1 + 2, n2 + 4) = ur5*(f1 + 3.0d0*f2)
            hmbox(2, n1 + 3, n2 + 1) = ur5*(f1 + 3.0d0*f2)
            hmbox(2, n1 + 2, n2 + 2) = ur5*(f1 + 3.0d0*f2)
            hmbox(3, n1 + 1, n2 + 2) = ur5*(f1 + 3.0d0*f2)
            hmbox(3, n1 + 3, n2 + 4) = ur5*(f1 + 3.0d0*f2)
!
            hmbox(3, n1 + 2, n2 + 3) = ur15*(2.0d0*f1 + 6.0d0*f2)
!
            hmbox(1, n1 + 3, n2 + 3) = ur15*(-1.0d0*f1 - 3.0d0*f2)
            hmbox(2, n1 + 1, n2 + 3) = ur15*(-1.0d0*f1 - 3.0d0*f2)
!
            hmbox(1, n1 + 3, n2 + 5) = ur5*(1.0d0*f1 + 3.0d0*f2)
            hmbox(2, n1 + 1, n2 + 5) = ur5*(-1.0d0*f1 - 3.0d0*f2)
            else if (l1.eq.2 .and. l2.eq.1) then
            hmbox(1, n1 + 1, n2 + 1) = ur5*(f1 - f2)
            hmbox(1, n1 + 4, n2 + 2) = ur5*(f1 - f2)
            hmbox(2, n1 + 1, n2 + 3) = ur5*(f1 - f2)
            hmbox(2, n1 + 2, n2 + 2) = ur5*(f1 - f2)
            hmbox(3, n1 + 2, n2 + 1) = ur5*(f1 - f2)
            hmbox(3, n1 + 4, n2 + 3) = ur5*(f1 - f2)
!
            hmbox(3, n1 + 3, n2 + 2) = 2.0d0*ur15*(f1 - f2)
!
            hmbox(1, n1 + 3, n2 + 3) = ur15*(f2 - f1)
            hmbox(2, n1 + 3, n2 + 1) = ur15*(f2 - f1)
!
            hmbox(1, n1 + 5, n2 + 3) = ur5*(f1 - f2)
            hmbox(2, n1 + 5, n2 + 1) = ur5*(f2 - f1)
          end if
          n2 = n2 + 2*l2 + 1
         end do !do jssh = 1, species(ispecies)%nssh
         n1 = n1 + 2*l1 + 1
        end do !do issh = 1, species(ispecies)%nssh

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine get_gover_1c

! ===========================================================================
! destroy_Fdata_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the one-center (exchange-correlation)
!! interactions for McWEDA - these arrays are read in by read_vxc_1c.
!
! ===========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis\n
!! Box 6315, 209 Hodges Hall\n
!! Department of Physics\n
!! West Virginia University\n
!! Morgantown, WV 26506-6315\n
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_Fdata_1c ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies

! Procedure
! ===========================================================================
        do ispecies = 1, nspecies
          deallocate (vxc_1c(ispecies)%dE)
          deallocate (vxc_1c(ispecies)%V)
          deallocate (vxc_1c(ispecies)%dV)
          ! goverlap
          if (associated(gover_1c_L(ispecies)%f1nac)) deallocate (gover_1c_L(ispecies)%f1nac)
          if (associated(gover_1c_L(ispecies)%f2nac)) deallocate (gover_1c_L(ispecies)%f2nac)
          if (associated(gover_1c_R(ispecies)%f1nac)) deallocate (gover_1c_R(ispecies)%f1nac)
          if (associated(gover_1c_R(ispecies)%f2nac)) deallocate (gover_1c_R(ispecies)%f2nac)
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (vxc_1c)
        deallocate (gover_1c_L)
        deallocate (gover_1c_R)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_Fdata_1c


! End Module
! ===========================================================================
        end module M_Fdata_1c
