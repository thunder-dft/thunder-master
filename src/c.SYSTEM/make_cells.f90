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
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! ===========================================================================
! make_cells
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine makes all the periodic cells up to the maximum numbers
!! of cells which forms a 9x9x9 maxumim number of cells.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine make_cells (t)
        use M_precision
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ix, iy, iz                     !< counters for dimensions
        integer logfile                        !< writing to which unit
        integer mbeta                          !< counter over cells
        integer mbox

        real, dimension (3) :: cvec

        interface
          function a_cross_b (a, b)
            real, dimension (3) :: a_cross_b
            real, intent(in), dimension (3) :: a, b
          end function a_cross_b
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile

! Calculate the cell volume and the reciprocal lattice vectors
! First get a2 X a3.
        cvec = a_cross_b (t%lattice(2)%a, t%lattice(3)%a)

! Next find the volume of the cell.
! NOTE: t%volume actually has a sign below. At this point the sign is
! important since we form g vectors by dividing by a1 dot (a2 X a3).
! Oh, you say. what difference does it make if we change the sign of g.
! it makes no difference in principle.
        t%volume = t%lattice(1)%a(1)*cvec(1) + t%lattice(1)%a(2)*cvec(2)     &
     &            + t%lattice(1)%a(3)*cvec(3)
        t%g(1)%a = 2.0d0*pi*cvec(:)/t%volume

! Next we get a3 X a1, and g2.
        cvec = a_cross_b (t%lattice(3)%a, t%lattice(1)%a)
        t%g(2)%a = 2.0d0*pi*cvec(:)/t%volume

! Finally we get a1 X a2, and g3.
        cvec = a_cross_b (t%lattice(1)%a, t%lattice(2)%a)
        t%g(3)%a = 2.0d0*pi*cvec(:)/t%volume
        t%volume = abs(t%volume)

        write (logfile,*)
        write (logfile,*) ' Cell volume [cubic Angstroms] = ', t%volume

! If we are doing a cluster then just set all to zero and return.
        if (t%icluster .eq. 1) then
          allocate (t%xl (0:0))
          t%xl(0)%a = 0.0d0
          return
        end if

! Initially, we set mbox = P_mbox which is rather large. Later we can
! change this so that mbox can be smaller if we do not need something so
! large.
        mbox = P_mbox
        mbeta_max = (2 * P_mbox + 1)**3 - 1
        allocate (t%xl (0:mbeta_max))

! Next - fill up xl with real angstrom units.
        mbeta = 1
        t%xl(0)%a = 0.0d0
        ! mbeta = 0 corresponds to the actual unit cell.
        do iz = -mbox, mbox
          do iy = -mbox, mbox
            do ix = -mbox, mbox
              if (.not. (ix .eq. 0 .and. iy .eq. 0 .and. iz .eq. 0)) then
                t%xl(mbeta)%a = real(ix)*t%lattice(1)%a                      &
     &                         + real(iy)*t%lattice(2)%a                     &
     &                         + real(iz)*t%lattice(3)%a
                mbeta = mbeta + 1
              end if
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
        end subroutine make_cells
