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

! M_neighbors_vdW
! Module Description
! ===========================================================================
!>       This module contains all the subroutines necessary for creating the
!! neighbor mapping of the system. It contains the following subroutines
!! within this module:
!!
!!       read_neighbors_vdW - reads the neighbor information from the
!!                        NEIGHBORS file which is an output file from
!!                        a previous run
!!       find_neighbors_vdW - finds all the neighbors for the two-center and
!!                        three-center interactions
!!       writeout_neighbors - write out the neighbor mapping into a file
!!       destroy_neighbors - destroy the neighbor arrays
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
!
! Module Declaration
! ===========================================================================
        module M_neighbors_vdW
        use M_configuraciones

! Type Declaration
! ===========================================================================
        type T_neighbors_vdW
          integer neighn                         ! number of neighbors
          integer ncommon                        ! number of common neighbors

          integer, pointer :: neigh_b (:)        ! which cell is the neighbor
          integer, pointer :: neigh_j (:)        ! which atom is the neighbor
        end type T_neighbors_vdW

        type (T_neighbors_vdW), pointer :: neighbors_vdW (:)

        type node_neighbor_vdW
          type (node_neighbor_vdW), pointer :: next
          integer :: neigh_b, neigh_j
        end type node_neighbor_vdW

        ! create a neighbor linked list
        type (node_neighbor_vdW), pointer :: root, current, temp
        integer :: length

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! driver_neighbors_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is the driver for the neighbor's routines.
!
! ===========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine driver_neighbors_vdW (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! FIXME - Put these into a module later
        integer ifix_neighbors
        logical iwriteout_neighbors_vdW

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (ilogfile,*)
!        write (ilogfile,*) ' Welcome to neighbors_vdW - determine neighbor '
!        write (ilogfile,*) ' mapping for the van der Waals interactions. '

        ifix_neighbors = 0
        if (ifix_neighbors .eq. 1) then
          call read_NEIGHBORS_vdW (s)
        else
          call find_neighbors_vdW (s)
        end if

        iwriteout_neighbors_vdW = .false.
        if (iwriteout_neighbors_vdW) call writeout_neighbors_vdW (s)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine driver_neighbors_vdW


! ===========================================================================
! read_neighbors_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Sometimes we want to read in the NEIGHBORS_vdW file, thus we keep
!! the number of atoms fixed.
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
!
! Subroutine Declaration
! ===========================================================================
        subroutine read_neighbors_vdW (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom                 !< counter over atoms
        integer ineigh, jneigh               !< counter over neighbors

        integer num_neigh                    !< number of neighbors

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Open the file NEIGHBORS and read in the neighbor map.
        open (11, file = 'NEIGHBORS_vdW', status = 'old')
        read (11, *)
        do iatom = 1, s%natoms
          read (11,*) jatom, num_neigh
          neighbors_vdW(iatom)%neighn = num_neigh
          allocate (neighbors_vdW(iatom)%neigh_b(num_neigh))
          allocate (neighbors_vdW(iatom)%neigh_j(num_neigh))
          do ineigh = 1, num_neigh
            read (11,*) jneigh, neighbors_vdW(iatom)%neigh_b(ineigh),        &
     &                          neighbors_vdW(iatom)%neigh_j(ineigh)
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
        end subroutine read_neighbors_vdW


! ===========================================================================
! find_neighbors_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Finds all the neighbors and fills up the neighbor map.
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
!
! Subroutine Declaration
! ===========================================================================
        subroutine find_neighbors_vdW (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom                   !< counter over atoms
        integer in1, in2                       !< species numbers
!       integer issh                           !< counter over shells
        integer mbeta                          !< counter over shells
        integer num_neigh                      !< number of neighbors counter

        real range_max                         !< maximum range between atoms
!       real rcutoff_i, rcutoff_j              !< maximum cutoffs
        real z                                 !< distance between atom pair

        real, dimension (3) :: r1, r2          !< atom positions

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Loop over all atoms in the central basis set and find his neighbours.
        allocate (neighbors_vdW(s%natoms))

        do iatom = 1, s%natoms
          allocate(root)
          current => root
          length = 0
          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
!         rcutoff_i = -99.0d0
!         do issh = 1, species(in1)%nssh
!           rcutoff_i = max(rcutoff_i,species(in1)%shell(issh)%rcutoffA)
!         end do

! Loop over all possible neighbors
          do mbeta = 0, mbeta_max
            do jatom = 1, s%natoms
              in2 = s%atom(jatom)%imass
              r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
!             rcutoff_j = -99.0d0
!             do issh = 1, species(in2)%nssh
!               rcutoff_j = max(rcutoff_j,species(in2)%shell(issh)%rcutoffA)
!             end do

! Find the distance from (mbeta,jatom) to (0,iatom)
              z = distance (r1, r2)
!             range_max = rcutoff_i + rcutoff_j - 0.01d0
              range_max = 20.0d0
              if (z .lt. range_max) then
                if (z .lt. 0.5d0 .and. z .gt. 1.0d-5) then
                  write (ilogfile,*) ' WARNING - atoms dangerously close! '
                  write (ilogfile,*) ' iatom, jatom, distance = ',            &
     &                                iatom, jatom, z
                  write (ilogfile,*) ' Check your coordinate file, start again. '
                  stop
                end if
                !Add to the linked list
                current%neigh_j = jatom
                current%neigh_b = mbeta
                allocate (current%next)
                current => current%next
                length = length + 1
              end if
            end do ! loop over jatom (the neighbor)
          end do ! mbeta over the cells
          !transfer the linked list to an array
          allocate (neighbors_vdW(iatom)%neigh_j (length))
          allocate (neighbors_vdW(iatom)%neigh_b (length))
          current => root
          do num_neigh = 1, length
            temp => current
            neighbors_vdW(iatom)%neigh_j(num_neigh) = current%neigh_j
            neighbors_vdW(iatom)%neigh_b(num_neigh) = current%neigh_b
            current => current%next
            deallocate (temp)
          end do
          neighbors_vdW(iatom)%neighn = length
          deallocate(current)
        end do ! loop over iatom

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine find_neighbors_vdW


! ===========================================================================
! writeout_neighbors_vdW.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Write out the van der Waals neighbors for restart capabilities.
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
!
! Updated by Barry Haycock.
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine writeout_neighbors_vdW (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom                   !< the two atom-neighbors
        integer ineigh                         !< the neighbor of iatom
        integer mbeta
        integer num_neigh

        real z                                 !< distance between atoms

        real, dimension (3) :: r1, r2          !< atom positions

        logical iwriteout_neighbors

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================

        iwriteout_neighbors = .false.
        if (iwriteout_neighbors) then
          write (ilogfile,*) '  '
          write (ilogfile,*) ' Neighbors of each atom: '
          do iatom = 1, s%natoms
            r1 = s%atom(iatom)%ratom
            num_neigh =  neighbors_vdW(iatom)%neighn
            write (ilogfile,*) '  '
            write (ilogfile,100)
            write (ilogfile,101) iatom, num_neigh
            write (ilogfile,102)
            write (ilogfile,100)
            do ineigh = 1, num_neigh
              mbeta = neighbors_vdW(iatom)%neigh_b(ineigh)
              jatom = neighbors_vdW(iatom)%neigh_j(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
              z = distance (r1, r2)
              write (ilogfile,103) iatom, ineigh, mbeta, jatom, z, r2 - r1
            end do
          end do
        end if

! Open the file NEIGHBORS which contain the neighbor map information for
! restart purposes.
          open (unit = 11, file = 'NEIGHBORS_vdW', status = 'unknown')

          write (11,104) s%natoms, s%basisfile
          do iatom = 1, s%natoms
            num_neigh =  neighbors_vdW(iatom)%neighn
            write (11,*) iatom, num_neigh
            do ineigh = 1, num_neigh
              mbeta = neighbors_vdW(iatom)%neigh_b(ineigh)
              jatom = neighbors_vdW(iatom)%neigh_j(ineigh)
              write (11,*) ineigh, mbeta, jatom
            end do
          end do
          close (unit = 11)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (70('='))
101     format (2x, ' Atom: ', i4, ',', ' Number of Neighbors: ', i4)
102     format (2x, ' iatom ', ' ineigh ', ' mbeta ', ' jatom ',             &
     &              ' distance ', 10x, ' vector ')
103     format (2x, i5, 4x, i3, 5x, i3, 3x, i4, 2x, f9.4, 3f9.4)
104     format (2x, i5, 2x, a40)

        return
        end subroutine writeout_neighbors_vdW


! ===========================================================================
! destroy_neighbors_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the van der Waals
!! neighbors information.
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
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_neighbors_vdW (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom

! Procedure
! ===========================================================================
        do iatom=1, s%natoms
          deallocate (neighbors_vdW(iatom)%neigh_j)
          deallocate (neighbors_vdW(iatom)%neigh_b)
        end do

        deallocate (neighbors_vdW)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_neighbors_vdW


! End Module
! ===========================================================================
        end module M_neighbors_vdW
