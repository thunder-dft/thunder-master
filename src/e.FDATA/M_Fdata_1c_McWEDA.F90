#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

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
!> @author Ning Ma
!> @author James P. Lewis\n
!! Box 6315, 209 Hodges Hall\n
!! Department of Physics\n
!! West Virginia University\n
!! Morgantown, WV 26506-6315\n
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
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
          real, pointer :: E(:,:)              !< xc-energy (shells)
          real, pointer :: dE(:,:,:)           !< 1st derivatives of energy

          real, pointer :: V(:, :)             !< Vxc potential
          real, pointer :: dV(:, :, :)         !< 1st derivatives of Vxc
        end type T_vxc_1c

! module variables
        type (T_vxc_1c), pointer :: vxc_1c (:) !< one-center xc for McWEDA

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
        integer ispecies                          !< counter for number of species
        integer nssh                              !< counters for number of shells
        integer issh, jssh
#ifdef DOGS
        integer kssh
#endif
        character (len=32) filename

! Allocate Arrays
! ===========================================================================
        allocate (vxc_1c (nspecies))

! Procedure
! ===========================================================================
! ***************************************************************************
!                 R E A D    M A T R I X    E L E M E N T S
! ***************************************************************************
        do ispecies = 1, nspecies

          ! Open ouput file for this species pair
          write (filename, '("/xc_1c", ".", i2.2, ".dat")')                  &
     &  	species(ispecies)%nZ
          open (11, file = trim(fdata_location)//trim(filename),             &
     &          status = 'old')

          nssh = species(ispecies)%nssh
          allocate (vxc_1c(ispecies)%E(nssh,nssh))
          allocate (vxc_1c(ispecies)%V(nssh,nssh))

          ! 0th order
          do issh = 1, nssh
            read (11,*) (vxc_1c(ispecies)%E(issh,jssh), jssh = 1, nssh)
          end do
          do issh = 1, nssh
            read (11,*) (vxc_1c(ispecies)%V(issh,jssh), jssh = 1, nssh)
          end do
        end do
        close (11)

#ifdef DOGS
! ***************************************************************************
!                       R E A D    P O T E N T I A L S
! ***************************************************************************
        do ispecies = 1, nspecies

          ! Open ouput file for this species pair
          write (filename, '("/nuxcrho_1c", ".", i2.2, ".dat")')             &
     &      species(ispecies)%nZ
          open (11, file = trim(fdata_location)//trim(filename),             &
     &            status = 'old')

          nssh = species(ispecies)%nssh
          allocate (vxc_1c(ispecies)%dV(nssh,nssh,nssh))

          ! 1st and 2nd order diagonal terms
          do kssh = 1,nssh
            do issh = 1, nssh
              read (11,*) (vxc_1c(ispecies)%dV(issh,jssh,kssh), jssh = 1, nssh)
            end do
          end do
        end do
        close (11)

! ***************************************************************************
!                       R E A D   E N E R G I E S
! ***************************************************************************
        do ispecies = 1, nspecies

          ! Open ouput file for this species pair
          write (filename, '("/excrho_1c", ".", i2.2, ".dat")')              &
     &      species(ispecies)%nZ
          open (11, file = trim(fdata_location)//trim(filename),             &
     &            status = 'old')

          nssh = species(ispecies)%nssh
          allocate (vxc_1c(ispecies)%dE(nssh,nssh,nssh))

          ! 1st and 2nd order diagonal terms
          do kssh = 1,nssh
            do issh = 1, nssh
              read (11,*) (vxc_1c(ispecies)%dE(issh,jssh,kssh), jssh = 1, nssh)
            end do
          end do
        end do
        close (11)
#endif

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
          deallocate (vxc_1c(ispecies)%E)
          deallocate (vxc_1c(ispecies)%V)
#ifdef DOGS
          deallocate (vxc_1c(ispecies)%dE)
          deallocate (vxc_1c(ispecies)%dV)
#endif
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (vxc_1c)

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
