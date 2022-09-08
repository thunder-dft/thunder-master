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

! M_neighbors
! Module Description
! ===========================================================================
!       This module contains all the subroutines necessary for creating the
! neighbor mapping of the system. It contains the following subroutines
! within this module:
!
!       read_NEIGHBORS - reads the neighbor information from the
!                        NEIGHBORS file which is an output file from
!                        a previous run
!       find_neigh_max - find the maximum number of neighbors
!       find_neighbors - finds all the neighbors for the two-center and
!                        three-center interactions
!       find_common_max - find maximum number of common neighbors
!       find_common_neighbors - find the actual common neighbors
!       writeout_neighbors - write out the neighbor mapping into a file
!       destroy_neighbors - destroy the neighbor arrays
!
! ===========================================================================
! Code written by:
! Barry Haycock
! James P. Lewis
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
        module M_neighbors
        use M_configuraciones

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! driver_neighbors
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This is the driver for the neighbor's routines.
!
! ===========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine driver_neighbors (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer logfile                     !< writing to which unit

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*)
        write (logfile,*) ' Welcome to neighbors - determine neighbor mapping. '
        if (ifix_neighbors .eq. 1) then
          call read_NEIGHBORS (s)
        else
          call find_neigh_max (s)
          call find_neighbors (s)
        end if
        call find_common_max (s)
        call find_common_neighbors (s)

        if (iwriteout_neighbors .eq. 1) call writeout_neighbors (s)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine driver_neighbors


! ===========================================================================
! read_NEIGHBORS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Sometimes we want to read in the neighbors file, thus we keep
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
        subroutine read_NEIGHBORS (s)
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
        integer inpfile                      !< reading from which unit
        integer num_neigh                    !< number of neighbors

        character (len = 25) :: slogfile

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        inpfile = s%inpfile

! Open the file neighbors and read in the neighbor map.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.NEIGHBORS'
        open (inpfile, file = slogfile, status = 'old')
        read (inpfile, *)
        do iatom = 1, s%natoms
          read (inpfile,*) jatom, num_neigh
          s%neighbors(iatom)%neighn = num_neigh
          allocate (s%neighbors(iatom)%neigh_b(num_neigh))
          allocate (s%neighbors(iatom)%neigh_j(num_neigh))
          do ineigh = 1, num_neigh
            read (inpfile,*) jneigh, s%neighbors(iatom)%neigh_b(ineigh),     &
     &                          s%neighbors(iatom)%neigh_j(ineigh),          &
     &                          s%neigh_self(iatom)
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
        end subroutine read_neighbors


! ===========================================================================
! find_neigh_max
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Finds all the maximum number of neighbors to atoms in the central cell.
! ===========================================================================
! Code written by:
!> @author Barry Haycock
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
        subroutine find_neigh_max (s)
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
        integer issh                           !< counter over shells
        integer mbeta                          !< counter over shells
        integer num_neigh                      !< number of neighbors counter

        real range_max                         !< maximum range between atoms
        real rcutoff_i, rcutoff_j              !< maximum cutoffs
        real z                                 !< distance between atom pair

        real, dimension (3) :: r1, r2          !< atom positions

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (s%neighbors(s%natoms))
        allocate (s%neigh_self(s%natoms))

! Procedure
! ===========================================================================
! Loop over all atoms.
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          rcutoff_i = -99.0d0
          do issh = 1, species(in1)%nssh
            rcutoff_i = max(rcutoff_i,species(in1)%shell(issh)%rcutoffA)
          end do

! Loop over all possible neighbors
          do mbeta = 0, mbeta_max
            do jatom = 1, s%natoms
              in2 = s%atom(jatom)%imass
              r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
              rcutoff_j = -99.0d0
              do issh = 1, species(in2)%nssh
                rcutoff_j = max(rcutoff_j,species(in2)%shell(issh)%rcutoffA)
              end do

! Find the distance from (mbeta,jatom) to (0,iatom)
              z = distance (r1, r2)
              range_max = rcutoff_i + rcutoff_j - 0.01d0
              if (z .lt. range_max) then
                num_neigh = num_neigh + 1
              end if
            end do ! loop over jatom (the neighbor)
          end do ! mbeta over the cells
          allocate (s%neighbors(iatom)%neigh_b(num_neigh))
          allocate (s%neighbors(iatom)%neigh_j(num_neigh))
          allocate (s%neighbors(iatom)%neigh_back(num_neigh))
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
        end subroutine find_neigh_max


! ===========================================================================
! find_neighbors
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Finds all the neighbors and fills up the neighbor map.
! ===========================================================================
! Code written by:
!> @author Barry Haycock
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
        subroutine find_neighbors (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom, katom            !< counter over atoms
        integer in1, in2                       !< species numbers
        integer ineigh, jneigh                 !< counter over neighbors
        integer issh                           !< counter over shells
        integer mbeta, jbeta                   !< counter over shells
        integer num_neigh                      !< number of neighbors counter

        real range_max                         !< maximum range between atoms
        real rcutoff_i, rcutoff_j              !< maximum cutoffs
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
! None

! Procedure
! ===========================================================================
! Loop over all atoms.
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          rcutoff_i = -99.0d0
          do issh = 1, species(in1)%nssh
            rcutoff_i = max(rcutoff_i,species(in1)%shell(issh)%rcutoffA)
          end do

! Loop over all possible neighbors
          do mbeta = 0, mbeta_max
            do jatom = 1, s%natoms
              in2 = s%atom(jatom)%imass
              r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
              rcutoff_j = -99.0d0
              do issh = 1, species(in2)%nssh
                rcutoff_j = max(rcutoff_j,species(in2)%shell(issh)%rcutoffA)
              end do

! Find the distance from (mbeta,jatom) to (0,iatom)
              z = distance (r1, r2)
              range_max = rcutoff_i + rcutoff_j - 0.01d0
              if (iatom .ne. jatom .and. z .lt. 0.5d0) then
                write (*,*) ' WARNING - atoms dangerously close! '
                write (*,*) ' iatom, jatom, distance = ', iatom, jatom, z
                write (*,*) ' iatom at: ', s%atom(iatom)%ratom
                write (*,*) ' jatom at: ', s%atom(jatom)%ratom + s%xl(mbeta)%a
                write (*,*) ' in cell mbeta = ', mbeta, ' xl = ', s%xl(mbeta)%a
                write (*,*) ' Check your coordinate file and start again. '
                write (*,*) ' Structure file is ', s%basisfile
                stop
              end if
              if (z .lt. range_max) then
                num_neigh = num_neigh + 1
                s%neighbors(iatom)%neigh_j(num_neigh) = jatom
                s%neighbors(iatom)%neigh_b(num_neigh) = mbeta
              end if
            end do ! loop over jatom (the neighbor)
          end do ! mbeta over the cells
          s%neighbors(iatom)%neighn = num_neigh
        end do ! loop over iatom

! Set up neigh_self and neigh_back.  The variable neigh_self(s%natoms) is the
! ineigh value for the "self interaction".  We will later need the
! neighbor-number of iatom with itself (neigh_self) in order to put the result
! of VNA_atom (doscentros) into VNA(mu,nu,neigh_self,iatom). The variable
! neigh_back(s%natoms) is used in the Mulliken charges since we take an average
! of two density matrix elements - iatom, ineigh pair and jatom, jneigh pair
! where jatom is the atom number of the ineigh neighbor of iatom and jneigh is
! neighbor number of iatom to jatom.

! Initialize to something ridiculous.
        s%neigh_self = -999
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            s%neighbors(iatom)%neigh_back(ineigh) = -999

            ! find neigh_self
            if (iatom .eq. jatom .and. mbeta .eq. 0) then
              if (s%neigh_self(iatom) .ne. -999) then
                write (*,*) ' We already found neigh_self, cannot find again! '
                write (*,*) ' iatom, neigh_self(iatom) = ', iatom, s%neigh_self(iatom)
                stop
              end if
              s%neigh_self(iatom) = ineigh
            end if

            ! find neigh_back
            do jneigh = 1, s%neighbors(jatom)%neighn
              jbeta = s%neighbors(jatom)%neigh_b(jneigh)
              katom = s%neighbors(jatom)%neigh_j(jneigh)
              if (iatom .ne. katom                                             &
     &             .or. abs(s%xl(mbeta)%a(1) + s%xl(jbeta)%a(1)) .gt. 1.0d-3   &
     &             .or. abs(s%xl(mbeta)%a(2) + s%xl(jbeta)%a(2)) .gt. 1.0d-3   &
     &             .or. abs(s%xl(mbeta)%a(3) + s%xl(jbeta)%a(3)) .gt. 1.0d-3) then
                  ! do nothing - there is no match here
              else
                if (s%neighbors(iatom)%neigh_back(ineigh) .ne. -999) then
                  write (*,*) ' We already found neigh_back, cannot find again! '
                  write (*,*) ' iatom, ineigh, neigh_back(iatom) = ',          &
     &                          iatom, ineigh, s%neighbors(iatom)%neigh_back(ineigh)
                  stop
                end if
                s%neighbors(iatom)%neigh_back(ineigh) = jneigh
              end if
            end do
          end do ! end loop over neighbors

          ! check neigh_back
          do ineigh = 1, s%neighbors(iatom)%neighn
            if (s%neighbors(iatom)%neigh_back(ineigh) .eq. -999) then
              write (*,*) ' Cannot find neigh_back for iatom = ', iatom
              write (*,*) ' structure = ', s%basisfile
              stop
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

        ! check neigh_self
        do iatom = 1, s%natoms
          if (s%neigh_self(iatom) .eq. -999) then
            write (*,*) ' Cannot find neigh_self for iatom = ', iatom
            write (*,*) ' structure = ', s%basisfile
            stop
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
        end subroutine find_neighbors


! ===========================================================================
! find_common_max
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Given atom ialpha (the third center), we find pairs of atoms that are
!! common neighbors of atom ialpha. At this stage we just find the number of
!! common neighbors for allocation purposes.
!!
!! Find all common neighbors of atom alpha (we call its i value ialpha).
!! atom ialpha is in the central cell, and the common neighbors are
!! at l-vector-sub-i,i and l-vector-sub-j,j
!!
!! Note:
!! You should remember that in counting the common neighbors, we do not include
!! any pairs in which one or both of the atoms is ontop atom ialpha. Thus we are
!! keeping only third party common neighbors.
!!
!! See the defined type above to see how the common neighbor information is
!! stored.
!
! ===========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine find_common_max (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha                         !< counter over atoms
        integer iatom, jatom                   !< counter over atoms
        integer ibeta, jbeta                   !< cells for iatom and jatom
        integer in1, in2                       !< species numbers
        integer ineigh, jneigh                 !< counter over neighbors
        integer issh                           !< counter over shells
        integer num_neigh                      !< number of neighbors counter

        real range_max                         !< maximum range between atoms
        real rcutoff_i, rcutoff_j              !< maximum cutoffs
        real z                                 !< distance between atom pair

        real, dimension (3) :: r1, r2          !< atom positions

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

! ===========================================================================
! Procedure
! ===========================================================================
! Loop over all atoms
        do ialpha = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

! Loop over all known neighbors of ialp. Call these atoms ineigh.
! The neighbors array was previously determined in the find_neighbors
! routine.
          do ineigh = 1, s%neighbors(ialpha)%neighn
            iatom = s%neighbors(ialpha)%neigh_j(ineigh)
            ibeta = s%neighbors(ialpha)%neigh_b(ineigh)

! Keep only third party common neighbors.
            if (.not. (iatom .eq. ialpha .and. ibeta .eq. 0)) then
              in1 = s%atom(iatom)%imass
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              rcutoff_i = -99.0d0
              do issh = 1, species(in1)%nssh
                rcutoff_i = max(rcutoff_i,species(in1)%shell(issh)%rcutoffA)
              end do

! Loop over all known neighbors of atom ialpha again, and check to see if
! atoms iatom and jatom are neighbors. If so, then we have found a pair of
! common neighbors. Of course, if iatom is jatom in all respects, that does
! not count.
              do jneigh = 1, s%neighbors(ialpha)%neighn
                jatom = s%neighbors(ialpha)%neigh_j(jneigh)
                jbeta = s%neighbors(ialpha)%neigh_b(jneigh)

! Keep only third party common neighbors.
                if (.not. (jatom .eq. ialpha .and. jbeta .eq. 0)             &
     &              .and. (ineigh .ne. jneigh)) then
                  in2 = s%atom(jatom)%imass
                  r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
                  rcutoff_j = -99.0d0
                  do issh = 1, species(in2)%nssh
                    rcutoff_j = max(rcutoff_j,species(in2)%shell(issh)%rcutoffA)
                  end do

! Find the distance from (mbeta,jatom) to (0,iatom)
                  z = distance (r1, r2)
                  range_max = rcutoff_i + rcutoff_j - 0.01d0
                  if (z .lt. range_max) then
                    num_neigh = num_neigh + 1
                  end if
                end if ! end if for jatom .ne. ialpha
              end do ! jneigh neighbor of ialpha loop
            end if ! end if for iatom .ne. ialpha
          end do ! ineigh neighbor of ialpha loop
          allocate (s%neighbors(ialpha)%iatom_common_b(num_neigh))
          allocate (s%neighbors(ialpha)%iatom_common_j(num_neigh))
          allocate (s%neighbors(ialpha)%jatom_common_b(num_neigh))
          allocate (s%neighbors(ialpha)%jatom_common_j(num_neigh))
          allocate (s%neighbors(ialpha)%neigh_common(num_neigh))
        end do  ! ialpha loop

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine find_common_max


! ===========================================================================
! find_common_neighbors
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Given atom ialpha (the third center), we find pairs of atoms that are
!! common neighbors of atom ialpha.
!!
!! Find all common neighbors of atom alpha (we call its i value ialpha).
!! atom ialpha is in the central cell, and the common neighbors are
!! at l-vector-sub-i,i and l-vector-sub-j,j
!!
!! Note:
!! You should remember that in counting the common neighbors, we do not include
!! any pairs in which one or both of the atoms is ontop atom ialp. Thus we are
!! keeping only third party common neighbors.
!!
!! See the defined type above to see how the common neighbor information is
!! stored.
!
! ===========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine find_common_neighbors (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha                         !< counter over atoms
        integer iatom, jatom, katom            !< counter over atoms
        integer ibeta, jbeta, kbeta            !< cells for iatom and jatom
        integer in1, in2                       !< species numbers
        integer ineigh, jneigh, kneigh         !< counter over neighbors
        integer issh                           !< counter over shells
        integer num_neigh                      !< number of neighbors counter

        real range_max                         !< maximum range between atoms
        real rcutoff_i, rcutoff_j              !< maximum cutoffs
        real z                                 !< distance between atom pair

        real, dimension (3) :: r1, r2, r3      !< atom positions

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

! ===========================================================================
! Procedure
! ===========================================================================
! Loop over all atoms
        do ialpha = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

! Loop over all known neighbors of ialp. Call these atoms ineigh.
! The neighbors array was previously determined in the find_neighbors
! routine.
          do ineigh = 1, s%neighbors(ialpha)%neighn
            iatom = s%neighbors(ialpha)%neigh_j(ineigh)
            ibeta = s%neighbors(ialpha)%neigh_b(ineigh)

! Keep only third party common neighbors.
            if (.not. (iatom .eq. ialpha .and. ibeta .eq. 0)) then
              in1 = s%atom(iatom)%imass
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              rcutoff_i = -99.0d0
              do issh = 1, species(in1)%nssh
                rcutoff_i = max(rcutoff_i,species(in1)%shell(issh)%rcutoffA)
              end do

! Loop over all known neighbors of atom ialpha again, and check to see if
! atoms iatom and jatom are neighbors. If so, then we have found a pair of
! common neighbors. Of course, if iatom is jatom in all respects, that does
! not count.
              do jneigh = 1, s%neighbors(ialpha)%neighn
                jatom = s%neighbors(ialpha)%neigh_j(jneigh)
                jbeta = s%neighbors(ialpha)%neigh_b(jneigh)

! Keep only third party common neighbors.
                if (.not. (jatom .eq. ialpha .and. jbeta .eq. 0)             &
     &              .and. (ineigh .ne. jneigh)) then
                  in2 = s%atom(jatom)%imass
                  r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
                  rcutoff_j = -99.0d0
                  do issh = 1, species(in2)%nssh
                    rcutoff_j = max(rcutoff_j,species(in2)%shell(issh)%rcutoffA)
                  end do

! Find the distance from (mbeta,jatom) to (0,iatom)
                  z = distance (r1, r2)
                  range_max = rcutoff_i + rcutoff_j - 0.01d0
                  if (z .lt. range_max) then
                    num_neigh = num_neigh + 1

! Put iatom in first spot
! Put jatom in second spot
                    s%neighbors(ialpha)%iatom_common_j(num_neigh) = iatom
                    s%neighbors(ialpha)%iatom_common_b(num_neigh) = ibeta
                    s%neighbors(ialpha)%jatom_common_j(num_neigh) = jatom
                    s%neighbors(ialpha)%jatom_common_b(num_neigh) = jbeta

! We also need to know for a given ialpha and (iatom,ibeta), what is the mth
! value for (jatom,jbeta) with respect to iatom. That is, jatom is the
! m'th neighbor of iatom. What is m?

! Set to a crazy value, in case loop fails.
                    s%neighbors(ialpha)%neigh_common(num_neigh) = -9999
                    do kneigh = 1, s%neighbors(iatom)%neighn
                      katom = s%neighbors(iatom)%neigh_j(kneigh)
                      kbeta = s%neighbors(iatom)%neigh_b(kneigh)
                      r3 = s%xl(kbeta)%a + s%atom(katom)%ratom - s%atom(iatom)%ratom
                      if ((abs(r3(1) - r2(1) + r1(1)) .lt. 1.0d-4) .and.      &
     &                    (abs(r3(2) - r2(2) + r1(2)) .lt. 1.0d-4) .and.      &
     &                    (abs(r3(3) - r2(3) + r1(3)) .lt. 1.0d-4)) then
                        s%neighbors(ialpha)%neigh_common(num_neigh) = kneigh
                        exit
                      end if
                    end do ! end kneigh loop

! We check to see if it really was a neighbor. (Could be very close to cutoff,
! so it was a mistake to count it.)
                    if (s%neighbors(ialpha)%neigh_common(num_neigh) .eq. -9999) &
     &                num_neigh = num_neigh - 1
                  end if
                end if ! end if for jatom .ne. ialpha
              end do ! jneigh neighbor of ialpha loop
            end if ! end if for iatom .ne. ialpha
          end do ! ineigh neighbor of ialpha loop
          s%neighbors(ialpha)%ncommon = num_neigh
        end do  ! ialpha loop

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine find_common_neighbors


! ===========================================================================
! writeout_neighbors.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Write out the neighbors for restart capabilities.
!
! ===========================================================================
! Code written by:
!> @author Barry Haycock
!! @author James P. Lewis
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
        subroutine writeout_neighbors (s)
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
        integer inpfile                        !< reading from which unit
        integer logfile                        !< writing to which unit
        integer mbeta
        integer num_neigh

        real z                                 !< distance between atoms

        real, dimension (3) :: r1, r2          !< atom positions

        character (len = 25) :: slogfile

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
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

        write (logfile,*) '  '
        write (logfile,*) ' Neighbors of each atom: '
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          num_neigh =  s%neighbors(iatom)%neighn
          write (logfile,*) '  '
          write (logfile,100)
          write (logfile,101) iatom, num_neigh
          write (logfile,102)
          write (logfile,100)
          do ineigh = 1, num_neigh
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            z = distance (r1, r2)
            write (logfile,103) iatom, ineigh, mbeta, jatom, z, r2 - r1
          end do
        end do

! Open the file neighbors which contain the neighbor map information for
! restart purposes.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.NEIGHBORS'
        open (unit = inpfile, file = slogfile, status = 'unknown')
        write (inpfile,104) s%natoms, s%basisfile
        do iatom = 1, s%natoms
          num_neigh =  s%neighbors(iatom)%neighn
          write (inpfile,*) iatom, num_neigh
          do ineigh = 1, num_neigh
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            write (inpfile,*) ineigh, mbeta, jatom, s%neigh_self(iatom)
          end do
        end do
        close (unit = inpfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (70('='))
101     format (2x, ' Atom: ', i4, ',', ' Number of neighbors: ', i4)
102     format (2x, ' iatom ', ' ineigh ', ' mbeta ', ' jatom ',             &
     &              ' distance ', 10x, ' vector ')
103     format (2x, i5, 4x, i3, 5x, i3, 3x, i4, 2x, f9.4, 3f9.4)
104     format (2x, i5, 2x, a40)

        return
        end subroutine writeout_neighbors


! ===========================================================================
! destroy_neighbors
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine deallocates the arrays containing the neighbors
! information.
!
! ===========================================================================
! Code written by:
! Barry Haycock
! James P. Lewis
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
        subroutine destroy_neighbors (s)
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
        do iatom = 1, s%natoms
          deallocate (s%neighbors(iatom)%iatom_common_b)
          deallocate (s%neighbors(iatom)%iatom_common_j)
          deallocate (s%neighbors(iatom)%jatom_common_b)
          deallocate (s%neighbors(iatom)%jatom_common_j)
          deallocate (s%neighbors(iatom)%neigh_common)

          deallocate (s%neighbors(iatom)%neigh_j)
          deallocate (s%neighbors(iatom)%neigh_b)
        end do

        deallocate (s%neighbors)
        deallocate (s%neigh_self)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_neighbors


! End Module
! ===========================================================================
        end module M_neighbors
