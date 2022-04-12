! copyright info:
!
!                             @Copyright 2008
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
! West Virginia University - Ning Ma, Hao Wang and Khorgolkhuu Odbadrakh
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! mpairnay.f90
! Program Description
! ===========================================================================
!       This is a neighbor routine, that we need for the nonlocal
! pseudopotential. We need this for the M_assemble_vnl_3c routine.
!
! < iatom | jatom >, given iatom, jatom and r12 = r(iatom) - r(jatom),
! then what neighbor to iatom is jatom. In other words, what is mneigh
! for a given jatom.
!
! We need this because the interaction are stored as sVNL(mu,nu,iatom,mneigh)
! and not sVNL(mu,nu,iatom,jatom). So what the heck is mneigh for jatom?
!
! This routine is a function call.
!
! Variables are iatom, jatom, rdiff(3). Other input are "constants".
! The output is mpairnay, which is the mneigh value for iatom, jatom.
!
! ===========================================================================
! Code written by:
!> @author Khorgolkhuu Odbadrakh
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================

! Program Declaration
! ===========================================================================
        subroutine mpairnay (t, iatom, jatom, r12, mpair)
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: t            !< the structure to be used

        integer, intent(in) :: iatom, jatom

        real, intent(in), dimension (3) :: r12

! Output
        integer, intent (out) :: mpair

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer imatch                   !< if we have a match then set to unity
        integer ineigh                   !< counter over neighbors
        integer jatom_match              !< jatom wanting to match
        integer mbeta                    !< cell for test atom

        real difference

        real, dimension (3) :: r1, r2, r21

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize
        mpair = 0
        r1 = t%atom(iatom)%ratom

! Loop over the neighbors of iatom.
        imatch = 0
        do ineigh = 1, t%neighbors_PP(iatom)%neighn
          jatom_match = t%neighbors_PP(iatom)%neigh_j(ineigh)

! Eliminate the obvious.
          if (jatom_match .eq. jatom) then

! OK we have a candidate. If the vector connceting them is
! rdiff(3). Then we have a match. if there is no match, then
! we should not be here. Then there is a problem.
            mbeta = t%neighbors_PP(iatom)%neigh_b(ineigh)
            r2 = t%atom(jatom)%ratom + t%xl(mbeta)%a

! Vector from 1 to 2 is r21
            r21 = r2 - r1

! Now compare rdiff to r21.
            difference = distance (r12, r21)
            ! distance between r12 (input) and r21 (calculated)
            ! If there is a match, this should be zero

! We should find only one match.
            if (difference .lt. 0.0001d0) then
              imatch = imatch + 1
              mpair = ineigh
            end if
          end if
        end do

! Sanity checks
        if (imatch .ne. 1) then
          write (*,*) ' imatch = ', imatch
          write (*,*) ' The variable imatch MUST be ONE! NO EXCEPTIONS '
          write (*,*) ' Bad imatch value in mpairnay.f90; must abort! '
          write (*,*) ' iatom, position = ', iatom, r1(:)
          write (*,*) ' jatom, position = ', jatom, r2(:)
          stop
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end subroutine mpairnay
