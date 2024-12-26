! copyright info:
!
!                             @copyright 2001
!                           fireball committee
! brigham young university - james p. lewis, chair
! arizona state university - otto f. sankey
! university of regensburg - juergen fritsch
! universidad de madrid - jose ortega

! other contributors, past and present:
! auburn university - jian jun dong
! arizona state university - gary b. adams
! arizona state university - kevin schmidt
! arizona state university - john tomfohr
! lawrence livermore national laboratory - kurt glaesemann
! motorola, physical sciences research labs - alex demkov
! motorola, physical sciences research labs - jun wang
! ohio university - dave drabold

!
! fireball-qmd is a free (gplv3) open project.

! this program IS free software: you can redistribute it and/or modify
! it under the terms of the gnu general public license as published by
! the free software foundation, either version 3 of the license, or
! (at your option) any later version.
!
! this program IS distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose.  see the
! gnu general public license for more details.
!
! you should have received a copy of the gnu general public license
! along with this program.  if not, see <http://www.gnu.org/licenses/>.


! blacsaba.f90
! program DESCRIPTION
! ===========================================================================
!       this subroutine DOES a(ia,ja)=alpha*a(ia,ja) for a blacs distributed
! matrix.
!
! ===========================================================================
! code written by:
! kurt r. glaesemann
! henry eyring center for theoretical chemistry
! department of chemistry
! university of utah
! 315 s. 1400 e.
! salt lake city, ut 84112-0850
! fax 801-581-4353
! office telephone 801-585-1078
! ===========================================================================
!
! program DECLARATION
! ===========================================================================
subroutine blacsaba( a, ia, ja, desca, alpha,  &
&                     mycol, myrow, npcol, nprow )

! argument declaration and description
! ===========================================================================
! input
        implicit none
        integer, intent (in) :: ia, ja, mycol, myrow, npcol, nprow
        integer, intent (in) :: desca( * )
        real*8, intent (in) :: alpha

! input/output
        real*8, intent (inout) :: a( * )

! local parameters and data declaration
! ===========================================================================
        integer             lld_
        parameter          ( lld_ = 9 )

! local variable declaration and description
! ===========================================================================
        integer       iacol, iarow, iia, jja

! procedure
! ===========================================================================

!     get grid parameters.

        call infog2l( ia, ja, desca, nprow, npcol, myrow, mycol, iia, jja,  &
        &              iarow, iacol )
! do it!
        if( myrow.eq.iarow .and. mycol.eq.iacol ) then
          a( iia+(jja-1)*desca( lld_ ) ) =  &
          &     a( iia+(jja-1)*desca( lld_ ) ) * alpha
        end if
        return
end subroutine blacsaba

