! copyright info:
!
!                             @Copyright 2013
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

! Qmixer.f90
! Function Description
! ===========================================================================
!>      This is a dummy routine. There is no charge mixing for Harris.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Qmixer (t, iscf_iteration, sigma)
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iscf_iteration

        ! This is the pointer to the current structure
        type(T_structure), target :: t

! Output
        real, intent (inout) :: sigma

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Set sigma to something crazy in order to kick out of the scf loop.
        sigma = -999.0d0
        if (.FALSE.) t%volume = 0
        if (.FALSE.) sigma = iscf_iteration

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Qmixer
