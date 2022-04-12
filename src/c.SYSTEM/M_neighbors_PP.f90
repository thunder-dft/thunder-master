! copyright info:
!                             @Copyright 2008
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Khorgolkhuu Odbadrakh
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

! M_neighbors_PP.f90
! Module Description
! ===========================================================================
!>       This module contains all the subroutines necessary for creating the
!! neighbor mapping of the system. It contains the following subroutines
!! within this module:
!!
!!       read_neighbors_PP - reads the neighbor information from the
!!                           NEIGHBORS_PP file which is an output file from
!!                           a previous run
!!       find_neigh_PP_max - find the maximum number of PP neighbors
!!       fine_neighbors_PP - finds all the neighbors for the pseudopotential
!!                           interactions.
!!       find_neigh_PPx_max - find the maximum number of PP neighbors
!!       fine_neighbors_PPx - finds all the neighbors for the pseudopotential
!!                            interactions.
!!       find_neigh_PPp_max - find the maximum number of PPp neighbors
!!       fine_neighbors_PPp - finds all the neighbors for the pseudopotential
!!                            interactions.
!!       find_common_PP_max - find the maximum number common PP neighbors
!!       find_common_neighbors_PP - finds all the neighbors for the
!!                                  pseudopotential interactions.
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
        module M_neighbors_PP
        use M_configuraciones

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer max_neigh

! module procedures
        contains

! ===========================================================================
! driver_neighbors_PP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is the driver for the neighbor's routines related to PP.
! It calls read_NEIGHBORS_PP, read_NEIGHBORS_PPx, find_neighbors_PP,
! find_neighbors_PPx, find_neighbors_PPp, find_common_neighbors_max,
! find_common_neighbors_PP and calculates all the neighbors.
!
! ===========================================================================
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine driver_neighbors_PP (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
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

        call find_neigh_PP_max (s)
        call find_neigh_PPx_max (s)

        if (ifix_neighbors .eq. 1) then
          call read_neighbors_PP (s)
        else
          call find_neighbors_PP (s)
          call find_neighbors_PPx (s)
        end if

        call find_neigh_PPp_max (s)
        call find_neighbors_PPp (s)

        call find_common_PP_max (s)
        call find_common_neighbors_PP (s)

        if (iwriteout_neighbors .eq. 1) call writeout_neighbors_PP (s)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine driver_neighbors_PP


! ===========================================================================
! read_neighbors_PP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Sometimes we want to read in the read_NEIGHBORS_PP file, thus we keep
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
        subroutine read_neighbors_PP (s)
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
        integer inpfile                     !< reading from which unit
        integer num_neigh                    !< number of neighbors

        character (len = 25) :: slogfile

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        inpfile = s%inpfile

! Open the file NEIGHBORS and read in the neighbor map.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.NEIGHBORS_PP'
        open (unit = inpfile, file = slogfile, status = 'old')
        read (inpfile, *)
        read (inpfile, *)

! neighbors_PP
! Loop over all the atoms in the system.
        do iatom = 1, s%natoms
          read (inpfile,*) jatom, num_neigh
          s%neighbors_PP(iatom)%neighn = num_neigh
          allocate (s%neighbors_PP(iatom)%neigh_b(num_neigh))
          allocate (s%neighbors_PP(iatom)%neigh_j(num_neigh))
          do ineigh = 1, num_neigh
            read (inpfile,*) jneigh, s%neighbors_PP(iatom)%neigh_b(ineigh),  &
     &                               s%neighbors_PP(iatom)%neigh_j(ineigh),  &
     &                               s%neighbors_PP_self(iatom)
          end do
        end do

! neighborsPPx
! Loop over all the atoms in the system.
        read (inpfile, *)
        do iatom = 1, s%natoms
          read (inpfile,*) jatom, num_neigh
          s%neighbors_PPx(jatom)%neighn = num_neigh
          allocate (s%neighbors_PPx(jatom)%neigh_b(num_neigh))
          allocate (s%neighbors_PPx(jatom)%neigh_j(num_neigh))
          do ineigh = 1, num_neigh
            read (inpfile,*) jneigh, s%neighbors_PPx(jatom)%neigh_b(ineigh), &
     &                               s%neighbors_PPx(jatom)%neigh_j(ineigh), &
     &                               s%neighbors_PPx_self(iatom)
          end do
        end do
        close (inpfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine read_neighbors_PP


! ===========================================================================
! find_neigh_PP_max
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> Only the allocation is performed at this stage!
!!
!!       Finds all the maximum number of neighbors to atoms in the central cell.
!! This routine is for the pseudopotential interactions.
!!
!!                   atom_VNL
!!                   +     +
!!                 +          +
!!               +               +
!!             +                    +
!!           +                         +
!!        atom 1                          +
!!                                       atom 2
!!
!! Generaly overlap of psi_1 and psi_2 can be zero but the PP interactions would
!! not be zero. So we have different number of pairs in case of PP interaction.
!! Pseudopotential interaction are calculated only for distance rcutoff_i
!! + rc_PP.  This is different in common interactions, where we have
!! rcutoff_i + rcutoff_j.
!!
!! In this case the first atom is the regular wavefunction and the neighbor
!! atom is the psuedopotential. The other case is covered in find_neighPPx_max.
! ===========================================================================
! Code written by:
!> @author Khorgolkhuu Odbadrakh
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
        subroutine find_neigh_PP_max (s)
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
        allocate (s%neighbors_PP(s%natoms))
        allocate (s%neighbors_PP_self(s%natoms))

! Procedure
! ===========================================================================
! Loop over all the atoms in the system.
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
              rcutoff_j = species(in2)%rcutoff_PP

! Find the distance from (mbeta,jatom) to (0,iatom)
              z = distance (r1, r2)
              range_max = rcutoff_i + rcutoff_j - 0.01d0
              if (z .lt. range_max) then
                num_neigh = num_neigh + 1
              end if
            end do ! loop over jatom (the neighbor)
          end do ! mbeta over the cells
          allocate (s%neighbors_PP(iatom)%neigh_b(num_neigh))
          allocate (s%neighbors_PP(iatom)%neigh_j(num_neigh))
          allocate (s%neighbors_PP(iatom)%map(num_neigh))
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
        end subroutine find_neigh_PP_max


! ===========================================================================
! find_neighbors_PP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> Only the allocation is performed at this stage!
!!
!!       Finds all the neighbors for the pseudopotential interactions, fills
!! up the neighbor_PP map. The distance of neighbor_PP pair must be less than
!! the sum of wave function and PP radius cutoff.
! ===========================================================================
! Code written by:
!> @author Khorgolkhuu Odbadrakh
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
        subroutine find_neighbors_PP (s)
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
        integer ineigh                         !< counter over neighbors
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
! None

! Procedure
! ===========================================================================
! First we'are going to build nPP list of neighbors. The list includes all
! neighbors of iatom with nonzero overlap <phi_i|Psi_j>, where
! phi_i .. is atomic wave function and Psi_j is pseudo wave function
! Loop over all atoms in the central basis set and find his neighbours.
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          rcutoff_i = -99.0d0
          do issh = 1, species(in1)%nssh
            rcutoff_i = max(rcutoff_i,species(in1)%shell(issh)%rcutoffA)
          end do

! Loop over all possible neighbors (VNL atoms)
          do mbeta = 0, mbeta_max
            do jatom = 1, s%natoms
              in2 = s%atom(jatom)%imass
              r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
              rcutoff_j = species(in2)%rcutoff_PP ! use PP cutoffs here

! Find the distance from (mbeta,jatom) to (0,iatom)
              z = distance (r1, r2)
              range_max = rcutoff_i + rcutoff_j - 0.01d0
              if (z .lt. range_max) then
                num_neigh = num_neigh + 1
                s%neighbors_PP(iatom)%neigh_j(num_neigh) = jatom
                s%neighbors_PP(iatom)%neigh_b(num_neigh) = mbeta
              end if
            end do ! loop over jatom (the neighbor)
          end do ! mbeta over the cells
          s%neighbors_PP(iatom)%neighn = num_neigh
        end do ! loop over iatom

! Set up neigh_self.  The variable neigh_self(s%natoms) is the ineigh value for
! the "self interaction".  Find the neighbor-number of iatom with itself
! (neigh_self) in order to put the result of vnl_atom (doscentros) into
! vnl(mu,nu,neigh_self,iatom).
        ! Initialize to something ridiculous.
        s%neighbors_PP_self = -999
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors_PP(iatom)%neighn
            mbeta = s%neighbors_PP(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)
            if (iatom .eq. jatom .and. mbeta .eq. 0) then
              s%neighbors_PP_self(iatom) = ineigh
              exit
            end if
          end do
          if (s%neighbors_PP_self(iatom) .eq. -999) then
            write (*,*) ' Cannot find neighbors_PP_self for atom ', iatom
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
        end subroutine find_neighbors_PP


! ===========================================================================
! find_neigh_PPx_max
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> Only the allocation is performed at this stage!
!!
!!       Finds all the maximum number of neighbors to atoms in the central cell.
!! This routine is for the pseudopotential interactions.
!!
!!                   atom_VNL
!!                   +     +
!!                 +          +
!!               +               +
!!             +                    +
!!           +                         +
!!        atom 1                          +
!!                                       atom 2
!!
!! Generaly overlap of psi_1 and psi_2 can be zero but the PP interactions would
!! not be zero. So we have different number of pairs in case of PP interaction.
!! Pseudopotential interaction are calculated only for distance rcutoff_i
!! + rc_PP.  This is different in common interactions, where we have
!! rcutoff_i + rcutoff_j.
!!
!! In this case the first atom is the regular wavefunction and the neighbor
!! atom is the psuedopotential. The other case is covered in find_neighPPx_max.
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
        subroutine find_neigh_PPx_max (s)
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
        allocate (s%neighbors_PPx(s%natoms))
        allocate (s%neighbors_PPx_self(s%natoms))

! Procedure
! ===========================================================================
! Loop over all atoms.
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          rcutoff_i = species(in1)%rcutoff_PP

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
          allocate (s%neighbors_PPx(iatom)%neigh_b(num_neigh))
          allocate (s%neighbors_PPx(iatom)%neigh_j(num_neigh))
          allocate (s%neighbors_PPx(iatom)%map(num_neigh))
          allocate (s%neighbors_PPx(iatom)%point(num_neigh))
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
        end subroutine find_neigh_PPx_max


! ===========================================================================
! find_neighbors_PPx
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> Only the allocation is performed at this stage!
!!
!!       Finds all the neighbors for the pseudopotential interactions, fills
!! up the neighbor_PPx map. The distance of neighbor_PP pair must be less than
!! the sum of wave function and PP radius cutoff.
! ===========================================================================
! Code written by:
!> @author Khorgolkhuu Odbadrakh
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
        subroutine find_neighbors_PPx (s)
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
        integer ineigh                         !< counter over neighbors
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
! None

! Procedure
! ===========================================================================
! First we'are going to build nPP list of neighbors. The list includes all
! neighbors of iatom with nonzero overlap <phi_i|Psi_j>, where
! phi_i .. is atomic wave function and Psi_j is pseudo wave function
! Loop over all atoms in the central basis set and find his neighbours.
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          rcutoff_i = species(in1)%rcutoff_PP

! Loop over all possible neighbors (VNL atoms)
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
                s%neighbors_PPx(iatom)%neigh_j(num_neigh) = jatom
                s%neighbors_PPx(iatom)%neigh_b(num_neigh) = mbeta
              end if
            end do ! loop over jatom (the neighbor)
          end do ! mbeta over the cells
          s%neighbors_PPx(iatom)%neighn = num_neigh
        end do ! loop over iatom

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SELF nPPx_self
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        s%neighbors_PPx_self = -999
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors_PPx(iatom)%neighn
            mbeta = s%neighbors_PPx(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PPx(iatom)%neigh_j(ineigh)
            if (iatom .eq. jatom .and. mbeta .eq. 0) then
              s%neighbors_PPx_self(iatom) = ineigh
              exit
            end if
          end do
          if (s%neighbors_PPx_self(iatom) .eq. -999) then
            write (*,*) ' Cannot find neighbors_PPx_self for atom ', iatom
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
        end subroutine find_neighbors_PPx


! ===========================================================================
! find_neigh_PPp_max (here, only allocation is performed)
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> Only the allocation is performed at this stage!
!!
!! We have two list of neighbors neighbors_PP <phi_i|Psi_j> and neighbors_PPx
!! <Psi_i|phi_j>. But what we really need is a following list of pairs
!! <phi_i|Vnl|phi_j> to be able map PP matrices to global Matrix <phi_i|V|phi_j>.
!! Let's call it neighbors_PPp, we're going to include all possible
!! pairs 1center, 2center and 3center. Also we need to perform mapping from
!! other lists to this central one.
!!
!! ===========================================================================
!! SCHEMATIC VIEW:
!!
!!           <Psi_k|        <phi_i|Vnl|phi_j> = <phi_i|Psi_k><Psi_k|phi_j>
!!           ialp
!!           * Rk (NL)
!!          / \                   neighbors_PP  ->  <phi_i|Psi_k>
!!         /   \                  neighbors_PPx ->  <Psi_k|phi_j>
!!neighbors_PP  \ neighbors_PPx
!!       /       \
!!      /         \
!!     + <phi_i|   + |phi_j>
!!     iatom       jatom
!!
!! ===========================================================================
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
!
! Subroutine Declaration
! ===========================================================================
        subroutine find_neigh_PPp_max (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ialpha                  !< counter over atoms
        integer ineigh, jneigh                 !< counter over neighbors
        integer ibeta, jbeta                   !< counter over shells
        integer num_neigh                      !< number of neighbors counter

! Allocate arrays
! ===========================================================================
        allocate (s%neighbors_PPp(s%natoms))
        allocate (s%neighbors_PPp_self(s%natoms))

! Procedure
! ===========================================================================
! Loop over all atoms.
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

! Loop over nPP ~ <phi_i|Psi_k>
          do ineigh = 1, s%neighbors_PP(iatom)%neighn
            ialpha = s%neighbors_PP(iatom)%neigh_j(ineigh)
            ibeta  = s%neighbors_PP(iatom)%neigh_b(ineigh)

! Loop over nPPx ~ <Psi_k|phi_j>
            do jneigh = 1, s%neighbors_PPx(ialpha)%neighn
              jbeta = s%neighbors_PPx(ialpha)%neigh_b(jneigh)
              num_neigh = num_neigh + 1
            end do ! do jneigh
          end do ! do ineigh

! Save number of pairs of iatom
          allocate (s%neighbors_PPp(iatom)%neigh_j(num_neigh))
          allocate (s%neighbors_PPp(iatom)%neigh_b(num_neigh))
        end do ! do iatom

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine find_neigh_PPp_max


! ===========================================================================
! find_neighbors_PPp
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> Only the allocation is performed at this stage!
!!
!! We have two list of neighbors neighbors_PP <phi_i|Psi_j> and neighbors_PPx
!! <Psi_i|phi_j>. But what we really need is a following list of pairs
!! <phi_i|Vnl|phi_j> to be able map PP matrices to global Matrix <phi_i|V|phi_j>.
!! Let's call it neighbors_PPp, we're going to include all possible
!! pairs 1center, 2center and 3center. Also we need to perform mapping from
!! other lists to this central one.
!!
!! ===========================================================================
!! SCHEMATIC VIEW:
!!
!!           <Psi_k|        <phi_i|Vnl|phi_j> = <phi_i|Psi_k><Psi_k|phi_j>
!!           ialp
!!           * Rk (NL)
!!          / \                   neighbors_PP  ->  <phi_i|Psi_k>
!!         /   \                  neighbors_PPx ->  <Psi_k|phi_j>
!!neighbors_PP  \ neighbors_PPx
!!       /       \
!!      /         \
!!     + <phi_i|   + |phi_j>
!!     iatom       jatom
!!
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
!
! Subroutine Declaration
! ===========================================================================
        subroutine find_neighbors_PPp (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom, ialpha, katom    !< counter over atoms
        integer ineigh, jneigh, kneigh         !< counter over neighbors
        integer mbeta, ibeta, jbeta, kbeta     !< counter over shells
        integer num_neigh                      !< number of neighbors counter

        real z                                 !< distance between atom pair

        real, dimension (3) :: r1, r2, r3, r   !< atom positions

        logical flag

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

! Procedure
! ===========================================================================
! Loop over all the atoms in the system.
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

! Loop over nPP ~ <phi_i|Psi_k>
          do ineigh = 1, s%neighbors_PP(iatom)%neighn
            ialpha = s%neighbors_PP(iatom)%neigh_j(ineigh)
            ibeta  = s%neighbors_PP(iatom)%neigh_b(ineigh)

! Loop over nPPx ~ <Psi_k|phi_j>
            do jneigh = 1, s%neighbors_PPx(ialpha)%neighn
              jatom = s%neighbors_PPx(ialpha)%neigh_j(jneigh)
              jbeta = s%neighbors_PPx(ialpha)%neigh_b(jneigh)

! Keeping in mind ialp is not in base unit cell, so we have to add lattice
! vector of ialp's unit cell to get real lattice vector with respect to
! base unit cell
             r3 = s%xl(jbeta)%a + s%xl(ibeta)%a

! now we have to find mbeta, which is the real lattice vector of jatom with
! respect to base unit cell.
! loop over all image unit cells
              do mbeta = 0 , mbeta_max
                r = s%xl(mbeta)%a
                z = distance (r, r3)
                if (z .lt. 1.0d-04) exit
              end do

! This is not all. Now we have to check if the jatom is not already inlcuded
! in the list, so let's search in temporary list of pairs neighPP
              flag = .false.
              do kneigh = 1, num_neigh
                 katom = s%neighbors_PPp(iatom)%neigh_j(kneigh)
                 kbeta = s%neighbors_PPp(iatom)%neigh_b(kneigh)
                 if (katom .eq. jatom .and. kbeta .eq. mbeta) flag = .true.
              end do ! do kneigh

! there isn't already included on the list
              if (.not. flag) then
                num_neigh = num_neigh + 1

! Check that # of common neighbor pairs is less than dimensions.
! Put jatom in second spot
                s%neighbors_PPp(iatom)%neigh_j(kneigh) = jatom
                s%neighbors_PPp(iatom)%neigh_b(kneigh) = mbeta
              end if ! if (flag)
            end do ! do jneigh
          end do ! do ineigh

! Save number of pairs of iatom
          s%neighbors_PPp(iatom)%neighn = num_neigh
        end do ! do iatom

! SELF nPP_self
! ===========================================================================
! What's the neighbor of the atom itself?
        s%neighbors_PPp_self = -999
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors_PPp(iatom)%neighn
            mbeta = s%neighbors_PPp(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PPp(iatom)%neigh_j(ineigh)
            if (iatom .eq. jatom .and. mbeta .eq. 0) then
              s%neighbors_PPp_self(iatom) = ineigh
              exit
            end if
          end do
          if (s%neighbors_PPp_self(iatom) .eq. -999) then
            write (*,*) ' Cannot find neighbors_PPp_self for atom ', iatom
            stop
          end if
        end do

! ===========================================================================
! MAP:  (neighbors_PP -> neighbors_PPp)
! ===========================================================================
! Now we're going to map neighbors_PP list into the global list neighbors_PPp
! Loop over atoms
        do iatom = 1, s%natoms
          num_neigh = s%neighbors_PP(iatom)%neighn

! Loop over neighbors_PPx-neighbors <Psi_i|phi_j> of iatom
          do ineigh = 1, num_neigh
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)
            jbeta = s%neighbors_PP(iatom)%neigh_b(ineigh)

! Loop over neighbors_PPp-neighbors <phi_i|phi_j> of iatom
            do kneigh = 1, s%neighbors_PPp(iatom)%neighn
              katom = s%neighbors_PPp(iatom)%neigh_j(kneigh)
              kbeta = s%neighbors_PPp(iatom)%neigh_b(kneigh)

! Test on iatom == katom
              if (katom .eq. jatom .and. kbeta .eq. jbeta) then
                s%neighbors_PP(iatom)%map(ineigh) = kneigh
                exit
              end if ! if (katom)
            end do ! do kneigh
          end do ! do ineigh
        end do ! do iatom

! ===========================================================================
!   MAP neighbors_PPx_map  (neighbors_PPx -> neighbors_PPp)
! ===========================================================================
! We need also a knowledge how place a local matrix into the global matrix.
! So, we're going to map the list neighbors_PPx to the global list
! neighbors_PPp
! Loop over atoms
        do iatom = 1, s%natoms

! Loop over nPPx-neighbors <Psi_i|phi_j> of iatom
          do ineigh = 1, s%neighbors_PPx(iatom)%neighn
            jatom = s%neighbors_PPx(iatom)%neigh_j(ineigh)
            jbeta = s%neighbors_PPx(iatom)%neigh_b(ineigh)

! Loop over neighPP-neighbors <phi_i|phi_j> of iatom
            do kneigh = 1, s%neighbors_PPp(iatom)%neighn
              katom = s%neighbors_PPp(iatom)%neigh_j(kneigh)
              kbeta = s%neighbors_PPp(iatom)%neigh_b(kneigh)

! Test on iatom == katom
              if (katom .eq. jatom .and. kbeta .eq. jbeta) then
                s%neighbors_PPx(iatom)%map(ineigh) = kneigh
                exit
              end if ! if (katom)
            end do ! do kneigh
          end do ! do ineigh
        end do ! do iatom

! ===========================================================================
!   MAP neighbors_PPx_point  (neighbors_PPx -> neighbors_PP)
! ===========================================================================
! Now we need mapping between neighbors_PP and neighbors_PPx.
! In the case of ontop left case we have <phi_i|Psi_i><Psi_i|phi_j>,
! so we will calculate it if <Psi_i|phi_j> /= 0. OK, but matrix sVNL gives
! you only <phi_j|Psi_i> overlap. So we need connection between the
! neighbors_PPx and neighbors_PP lists.

! Loop over atoms
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom

! Loop over neighbors_PPx-neighbors <Psi_i|phi_j> of iatom
          do ineigh = 1, s%neighbors_PPx(iatom)%neighn
            jatom = s%neighbors_PPx(iatom)%neigh_j(ineigh)
            jbeta = s%neighbors_PPx(iatom)%neigh_b(ineigh)

! Loop over neighbors_PP-neighbors <phi_j|Psi_k> of jatom
            do kneigh = 1, s%neighbors_PP(jatom)%neighn
              katom = s%neighbors_PP(jatom)%neigh_j(kneigh)
              kbeta = s%neighbors_PP(jatom)%neigh_b(kneigh)
              r2 = s%atom(katom)%ratom + s%xl(kbeta)%a + s%xl(jbeta)%a
              z = distance (r1, r2)
              if (z .lt. 1.0d-04) then
                s%neighbors_PPx(iatom)%point(ineigh) = kneigh
                exit
              end if ! if (distance)
            end do ! do jneigh
          end do ! do ineigh
        end do ! do iatom

! Deallocate arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine find_neighbors_PPp


! ===========================================================================
! find_common_PP_max
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>        Given atom ialpha (the third center), we find pairs of atoms that are
!! common neighbors of atom ialpha. At this stage we just find the number of
!! common neighbors for allocation purposes
!!
!! Find all common neighbors (Pseudopotential) of atom alpha
!! (we call its i value ialp). Atom alpha is in the central cell, and
!! the common neighbors are at l-vector-sub-i,i and l-vector-sub-j,j
!!
!! Note:
!! You should remember that in counting the common neighbors, we do not include
!! any pairs in which one or both of the atoms is ontop atom ialp. Thus we are
!! keeping only third party common neighbors.
!!
!! See the defined type above to see how the common neighbor information
!! is stored.
!! ****************************************************************************
!!      P P        C O M M O N      N E I G H B O R S     L I S T
!! ****************************************************************************
!! SCHEMATIC VIEW:
!!
!!         <Psi_k|        <phi_i|Vnl|phi_j> = <phi_i|Psi_k><Psi_k|phi_j>
!!          ialp
!!           * R2 (NL)
!!          / \                   nPPx  ->  <Psi_k|phi_i>
!!         /   \                  nPPx  ->  <Psi_k|phi_j>
!!   nPPx /     \ nPPx
!!       /       \
!!      /         \
!!     + <phi_i|   + |phi_j>
!!   iatom       jatom
!!    R1           R3
!! ****************************************************************************
!!
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
!===========================================================================
! Program Declaration
! ===========================================================================
        subroutine find_common_PP_max (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha, iatom, jatom           !< counter over atoms
        integer ibeta, jbeta                   !< cells for iatom and jatom
        integer ineigh, jneigh                 !< counter over neighbors
        integer num_neigh                      !< number of neighbors counter

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over all atoms.
        do ialpha = 1, s%natoms
          num_neigh = 0

! Loop over all known neighbors of ialp. Call these atoms ineigh.
          do ineigh = 1, s%neighbors_PPx(ialpha)%neighn
            iatom = s%neighbors_PPx(ialpha)%neigh_j(ineigh)
            ibeta = s%neighbors_PPx(ialpha)%neigh_b(ineigh)

! Keep only third party common neighbors.
            if (.not. (iatom .eq. ialpha .and. ibeta .eq. 0)) then

! Loop over all known neighbors of atom ialpha again, and check to see if
! atoms iatom and jatom are neighbors. If so, then we have found a pair of
! common neighbors. Of course, if iatom is jatom in all respects, that does
! not count.
              do jneigh = 1, s%neighbors_PPx(ialpha)%neighn
                jatom = s%neighbors_PPx(ialpha)%neigh_j(jneigh)
                jbeta = s%neighbors_PPx(ialpha)%neigh_b(jneigh)

! Keep only third party common neighbors.
                if (.not. (jatom .eq. ialpha .and. jbeta .eq. 0)             &
     &              .and. (ineigh .ne. jneigh))  &
                 then
                    num_neigh = num_neigh + 1
                end if ! end if for jatom .ne. ialpha
              end do ! jneigh neighbor of ialpha loop
            end if ! end if for iatom .ne. ialpha
          end do ! ineigh neighbor of ialpha loop
          allocate (s%neighbors_PP(ialpha)%iatom_common_b(num_neigh))
          allocate (s%neighbors_PP(ialpha)%iatom_common_j(num_neigh))
          allocate (s%neighbors_PP(ialpha)%jatom_common_b(num_neigh))
          allocate (s%neighbors_PP(ialpha)%jatom_common_j(num_neigh))
          allocate (s%neighbors_PP(ialpha)%neigh_common(num_neigh))
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
        end subroutine find_common_PP_max


! ===========================================================================
! find_common_neighbors_PP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> Given atom ialp, we find pairs of atoms that are common neighbors
!! of atom ialp. Note that on atom ialp is placed VNL.
!!
!! Find all common neighbors (Pseudopotential) of atom alpha
!! (we call its i value ialp). Atom alpha is in the central cell, and
!! the common neighbors are at l-vector-sub-i,i and l-vector-sub-j,j
!!
!! Note:
!! You should remember that in counting the common neighbors, we do not include
!! any pairs in which one or both of the atoms is ontop atom ialp. Thus we are
!! keeping only third party common neighbors.
!!
!! See the defined type above to see how the common neighbor information is
!! stored.
!!
!! ****************************************************************************
!!      P P        C O M M O N      N E I G H B O R S     L I S T
!! ****************************************************************************
!! SCHEMATIC VIEW:
!!
!!         <Psi_k|        <phi_i|Vnl|phi_j> = <phi_i|Psi_k><Psi_k|phi_j>
!!          ialp
!!           * R2 (NL)
!!          / \                   nPPx  ->  <Psi_k|phi_i>
!!         /   \                  nPPx  ->  <Psi_k|phi_j>
!!   nPPx /     \ nPPx
!!       /       \
!!      /         \
!!     + <phi_i|   + |phi_j>
!!   iatom       jatom
!!    R1           R3
!! ****************************************************************************
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
! Subroutine Declaration
! ===========================================================================
        subroutine find_common_neighbors_PP (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom, ialpha, katom    !< counter over atoms
        integer ibeta, jbeta, kbeta            !< counter over cells
        integer in1, in2
        integer ineigh, jneigh, kneigh         !< counter over neighbors
        integer num_neigh                      !< number of neighbors counter

        real z                                 !< distance between atom pair
        real, dimension (3) :: r1, r2, r3, r   !< atom positions

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over all atoms.
        do ialpha = 1, s%natoms
          num_neigh = 0

! 1.loop over <Psi_k|phi_i>  -> nPPx
          do ineigh = 1, s%neighbors_PPx(ialpha)%neighn
            iatom = s%neighbors_PPx(ialpha)%neigh_j(ineigh)
            ibeta = s%neighbors_PPx(ialpha)%neigh_b(ineigh)
            in1 = s%atom(iatom)%imass

! Keep only third party common neighbors. It means iatom /= ialp.
            if (.not. (iatom .eq. ialpha .and. ibeta .eq. 0)) then
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a

! Loop over all known neighbors of atom ialp again, and check to see if atoms
! iatom and jatom are neighbors. If so, then we have found a pair of common
! neighbors. Of course, if iatom is jatom in all respects, that does not count.
! 2.loop over <Psi_k|psi_j>  -> nPPx
              do jneigh = 1, s%neighbors_PPx(ialpha)%neighn
                jatom = s%neighbors_PPx(ialpha)%neigh_j(jneigh)
                jbeta = s%neighbors_PPx(ialpha)%neigh_b(jneigh)
                in2 = s%atom(jatom)%imass

! Keep only third party common neighbors. It means jatom /= ialp .and.
! jatom /= iatom
                if (.not. (jatom .eq. ialpha .and. jbeta .eq. 0)             &
       &         .and. (ineigh .ne. jneigh)) then
                r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a

! If we're here it means we've found third party common neighbors. Increase
! number of instances.
                num_neigh = num_neigh + 1
                if (num_neigh .gt. size(s%neighbors_PP(ialpha)%neigh_common)) then
                  write (*,*)                                                &
     &             ' num_neigh, size(s%neighbors_PP(ialpha)%neigh_common) ', &
     &               num_neigh, size(s%neighbors_PP(ialpha)%neigh_common)
                  stop
                end if ! if (num_neigh)

! Put iatom in first spot
! Put jatom in second spot
                s%neighbors_PP(ialpha)%iatom_common_j(num_neigh) = iatom
                s%neighbors_PP(ialpha)%iatom_common_b(num_neigh) = ibeta
                s%neighbors_PP(ialpha)%jatom_common_j(num_neigh) = jatom
                s%neighbors_PP(ialpha)%jatom_common_b(num_neigh) = jbeta

! We also need to know for a given ialpha and (iatom,ibeta), what is the m value
! for (jatom,jbeta) with respect to iatom. That is, jatom is the m'th neighbor
! of iatom. What is m?
! Set to a crazy value, in case loop fails.
                r = r2 - r1
                s%neighbors_PP(ialpha)%neigh_common(num_neigh) = -9999
                do kneigh = 1, s%neighbors_PPp(iatom)%neighn
                  katom = s%neighbors_PPp(iatom)%neigh_j(kneigh)
                  kbeta = s%neighbors_PPp(iatom)%neigh_b(kneigh)
                  r3 = s%xl(kbeta)%a + s%atom(katom)%ratom - s%atom(iatom)%ratom
                  z = distance (r, r3)
                  if (z .lt. 1.0d-4) then
                    s%neighbors_PP(ialpha)%neigh_common(num_neigh) = kneigh
                    exit
                  end if
                end do ! end kneigh loop

! We check to see if it really was a neighbor. (Could be very close to cutoff,
! so it was a mistake to count it.)
                if (s%neighbors_PP(ialpha)%neigh_common(num_neigh) .eq. -9999) &
     &           num_neigh = num_neigh - 1
                end if ! end if for jatom .ne. ialpha
              end do ! jneigh neighbor of ialpha loop
            end if ! end if for iatom .ne. ialpha
          end do ! ineigh neighbor of ialpha loop
          s%neighbors_PP(ialpha)%ncommon = num_neigh
        end do  ! ialpha loop

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================

! End Subroutine
! ===========================================================================
        return
        end subroutine find_common_neighbors_PP


! ===========================================================================
! writeout_neighbors_PP.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Write out the neighborsPP for restart capabilities.
!
! ===========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine writeout_neighbors_PP (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s         !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom                !< the two atom-neighbors
        integer ineigh                      !< the neighbor of iatom
        integer inpfile                     !< reading from which unit
        integer logfile                     !< writing to which unit
        integer mbeta
        integer num_neigh
        integer self                        !<self of neighbors_PP for writeout
        integer xself                       !<self of neighbors_PPx for writeout

        real z                              !< distance between atoms

        real, dimension (3) :: r1, r2       !< atom positions

        character (len = 25) :: slogfile

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
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
        write (logfile,*) ' Neighbors_PP of each atom: '
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          num_neigh =  s%neighbors_PP(iatom)%neighn
          write (logfile,*) '  '
          write (logfile,100)
          write (logfile,101) iatom, num_neigh
          write (logfile,102)
          write (logfile,100)
          do ineigh = 1, num_neigh
            mbeta = s%neighbors_PP(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            z = distance (r1, r2)
            write (logfile,103) iatom, ineigh, mbeta, jatom, z, r2 - r1
          end do
        end do

        write (logfile,*) '  '
        write (logfile,*) ' Neighbors_PPx of each atom: '
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          num_neigh =  s%neighbors_PPx(iatom)%neighn
          write (logfile,*) '  '
          write (logfile,100)
          write (logfile,101) iatom, num_neigh
          write (logfile,102)
          write (logfile,100)
          do ineigh = 1, num_neigh
            mbeta = s%neighbors_PPx(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PPx(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            z = distance (r1, r2)
            write (logfile,103) iatom, ineigh, mbeta, jatom, z, r2 - r1
          end do
        end do

! ****************************************************************************
!
! W R I T E    O U T    N E I G H B O R S    F I L E
! ****************************************************************************
! Open the file NEIGHBORS_PP which contain the neighbor map information for
! restart purposes.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.NEIGHBORS_PP'
        open (unit = inpfile, file = slogfile, status = 'unknown')
        write (inpfile, 104) s%natoms, s%basisfile

! neighborsPP
        write (inpfile,*) ' neighbors_PP '
        do iatom = 1, s%natoms
          self = s%neighbors_PP_self(iatom)
          num_neigh =  s%neighbors_PP(iatom)%neighn
          write (inpfile,*) iatom, num_neigh
          do ineigh = 1, num_neigh
            mbeta = s%neighbors_PP(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)
            write (inpfile,*) ineigh, mbeta, jatom, self
          end do
        end do

! neighborsPPx
        write (inpfile,*) ' neighbors_PPx '
        do iatom = 1, s%natoms
          xself = s%neighbors_PPx_self(iatom)
          num_neigh =  s%neighbors_PPx(iatom)%neighn
          write (inpfile,*) iatom, num_neigh
          do ineigh = 1, num_neigh
            mbeta = s%neighbors_PPx(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PPx(iatom)%neigh_j(ineigh)
            write (inpfile,*) ineigh, mbeta, jatom, xself
          end do
        end do
        close (unit = inpfile)

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
        end subroutine writeout_neighbors_PP


! ===========================================================================
! destroy_neighbors_PP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the neighbors
!! information.
!
! ===========================================================================
! Code written by:
!> @author Khorgolkhuu Odbadrakh
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
        subroutine destroy_neighbors_PP (s)
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
          deallocate (s%neighbors_PP(iatom)%iatom_common_b)
          deallocate (s%neighbors_PP(iatom)%iatom_common_j)
          deallocate (s%neighbors_PP(iatom)%jatom_common_b)
          deallocate (s%neighbors_PP(iatom)%jatom_common_j)
          deallocate (s%neighbors_PP(iatom)%neigh_common)

          deallocate (s%neighbors_PP(iatom)%neigh_j)
          deallocate (s%neighbors_PP(iatom)%neigh_b)
          deallocate (s%neighbors_PP(iatom)%map)

          deallocate (s%neighbors_PPx(iatom)%neigh_j)
          deallocate (s%neighbors_PPx(iatom)%neigh_b)
          deallocate (s%neighbors_PPx(iatom)%map)
          deallocate (s%neighbors_PPx(iatom)%point)

          deallocate (s%neighbors_PPp(iatom)%neigh_j)
          deallocate (s%neighbors_PPp(iatom)%neigh_b)
        end do

        deallocate (s%neighbors_PP)
        deallocate (s%neighbors_PP_self)
        deallocate (s%neighbors_PPx)
        deallocate (s%neighbors_PPx_self)
        deallocate (s%neighbors_PPp)
        deallocate (s%neighbors_PPp_self)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_neighbors_PP

! End Module
! ===========================================================================
        end module M_neighbors_PP
