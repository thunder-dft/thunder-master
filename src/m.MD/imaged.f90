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

! imaged.f90
! Function Description
! ===========================================================================
!       This routine will translate the atom positions during an MD
! simulation in order to keep atoms in the central cell always.
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
        subroutine imaged (t)
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        ! This is the pointer to the current structure
        type(T_structure), target :: t

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
! Transform coordinates to box centered on (0,0,0)  (note: this is shifted
! if ishiftO = 1)
        if (iimage .gt. 0 .and. icluster .eq. 0) then

!  otherwise we would always image first step
         if (mod(itime_step - nstepi, iimage) .eq. 0) then
          do iatom = 1, natoms
           from_zero = ratom(1,iatom)**2 + ratom(2,iatom)**2 + ratom(3,iatom)**2

! FIXME  - there is a 'GO TO' here that we need to eliminate
! Only need six boxes, since others are multiples
           do ibox = 1, 6
514         shifted_from_zero = (ratom(1,iatom) + xl(1,ibox))**2             &
     &                         + (ratom(2,iatom) + xl(2,ibox))**2            &
     &                         + (ratom(3,iatom) + xl(3,ibox))**2
            if (shifted_from_zero .lt.  from_zero) then
             ratom(:,iatom) = ratom(:,iatom) + xl(:,ibox)
             xdot(0,:,iatom) = xdot(0,:,iatom) + xl(:,ibox)
             ximage(:,iatom) = ximage(:,iatom) - xl(:,ibox)
             from_zero = shifted_from_zero
             go to 514

! Check again, just in case we jumped two boxes
            end if
           end do
          end do
         end if
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine imaged
