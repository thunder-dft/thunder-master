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

! M_vdW
! Module Description
! ===========================================================================
!>       This module contains all the subroutines necessary for the
!! van der Waals interactions. It contains the following subroutines
!! within this module:
!!
!!       initialize_vdW - reads in the van der Waals parameter file (if exists)
!!       read_vdW - reads the specific vdW connectivity from basisfile.vdW
!!       read_neighbors_vdW - reads the neighbor information from the
!!                        NEIGHBORS file which is an output file from
!!                        a previous run
!!       find_neighbors_vdW - finds all the neighbors for the two-center and
!!                        three-center interactions
!!       writeout_neighbors - write out the neighbor mapping into a file
!!       calculate_vdW - calculate the van der Waals energies and forces
!!       destroy_neighbors - destroy the neighbor arrays
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
        module M_vdW

! /GLOBAL
        use M_precision

! /SYSTEM
        use M_configuraciones
        use M_build_forces

! Type Declaration
! ===========================================================================
        type T_neighbors_vdW
          integer neighn                         ! number of neighbors

          integer, pointer :: neigh_b (:)        ! which cell is the neighbor
          integer, pointer :: neigh_j (:)        ! which atom is the neighbor
        end type T_neighbors_vdW

        type (T_neighbors_vdW), pointer :: neighbors_vdW (:)

        type T_vdW_parameter
          real C6
          real alpha
          real R0
        end type T_vdW_parameter

        type (T_vdW_parameter), pointer :: vdW_parameter (:)

! Variable Declaration and Description
! ===========================================================================
! Later we will replace this with a molecular mapping
        integer npairs_vdW

        integer, allocatable, dimension (:) :: iatom_vdW
        integer, allocatable, dimension (:) :: jatom_vdW

        real range_vdW                       ! set range of van der Waals

! module procedures
        contains


! ===========================================================================
! initialize_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine will read in the van der Waals parameters from
! the structures.vdW file.
!
! ===========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initialize_vdW
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies                    !< counter over species

        character (len = 25) :: filename

        logical file_exists

        type(T_vdW_parameter), pointer :: pvdW

! Allocate Arrays
! ===========================================================================
        allocate (vdW_parameter(nspecies))

! Procedure
! ===========================================================================
        filename = 'structures.vdW'
        inquire (file = filename, exist = file_exists)   ! file_exists will be TRUE if the file
                                                         ! exists and FALSE otherwise
        if ( file_exists ) then
           write (ilogfile,*) ' Reading: >'//filename//'<'
        end if

        open (unit = 222, file = filename, status = 'unknown')
        read (222,*) range_vdW
        write (ilogfile,*) ' Range of van der Waals interactions = ', range_vdW
        write (ilogfile,*)
        write (ilogfile,*) ' Parameters (units of C6 are in eV-Angstrom^6): '
        write (ilogfile,*)
        write (ilogfile,201)
        write (ilogfile,200)
        do ispecies = 1, nspecies
          ! cut some lengthy notation
          pvdW=>vdW_parameter(ispecies)
          read (222,*) jspecies, pvdW%C6, pvdW%alpha, pvdW%R0
          pvdW%C6 = pvdW%C6*P_Hartree*P_abohr**6
          pvdW%alpha = pvdW%alpha*P_abohr**3
          write (ilogfile,202) jspecies, pvdW%C6, pvdW%alpha, pvdW%R0
        end do
        write (ilogfile,200)
        close (unit = 222)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
200     format (2x, 70('='))
201     format (2x, ' Species # ', 2x, ' C6 ', 6x, ' alpha ', 5x, ' R0 ')
202     format (2x, i5, 2x, 3(2x,f9.3))

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_vdW


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
        integer logfile                     !< writing to which unit

        character (len = 25) :: slogfile

        logical read_vdW_parameters

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Skip this routine if there are no van der Waals parameters
        slogfile = 'structures.vdW'
        inquire (file = slogfile, exist = read_vdw_parameters)
        if (.not. read_vdW_parameters) return

        write (logfile,*) ' Welcome to neighbors_vdW - determine neighbor '
        write (logfile,*) ' mapping for the van der Waals interactions. '
        call find_neigh_vdW_max (s)
        call find_neighbors_vdW (s)
        call writeout_neighbors_vdW (s)

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
! read_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine will read in the van der Waals parameters from
! specific basis input file file.  In this file, the user can specify
! which atoms will include van der Waals connectivity. If all atoms will
! interact with van der Waals, then no file is specified.
!
! ===========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine read_vdW_parameters (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ipair                       !< loop over van der Waals pairs
        integer inpfile                     !< reading from which unit

        character (len = 25) :: slogfile

        logical read_vdW

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        inpfile = s%inpfile

        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.vdW'
        inquire (file = slogfile, exist = read_vdW)

! Here we read in the atom pairs that should be included in calculating
! the van der Waals interactions.
        open (unit = inpfile, file = slogfile, status = 'unknown')
        read (inpfile, *) npairs_vdW
        allocate (iatom_vdW(npairs_vdW))
        allocate (jatom_vdW(npairs_vdW))
        do ipair = 1, npairs_vdW
          read (inpfile, *) iatom_vdW(ipair), jatom_vdW(ipair)
        end do
        close (unit = inpfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
200     format (2x, 70('='))
201     format (2x, ' Species # ', 2x, ' C6 ', 6x, ' alpha ', 5x, ' R0 ')
202     format (2x, i5, 2x, 3(2x,f9.3))

! End Subroutine
! ===========================================================================
        return
        end subroutine read_vdW_parameters


! ===========================================================================
! find_neigh_vdW_max
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Finds all the maximum number of neighbors to atoms in the central
!> cell - for van der Waals interactions
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
        subroutine find_neigh_vdW_max (s)
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
        integer ipair                          !< counter of pairs
        integer mbeta                          !< counter over shells
        integer num_neigh                      !< number of neighbors counter

        real z                                 !< distance between atom pair

        real, dimension (3) :: r1, r2          !< atom positions

        character (len = 25) :: slogfile

        logical found
        logical read_vdW

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (neighbors_vdW(s%natoms))

! Procedure
! ===========================================================================
! Skip this routine if there are no van der Waals parameters
        slogfile = 'structures.vdW'
        inquire (file = slogfile, exist = read_vdW)
        if (.not. read_vdW) return

        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.vdW'
        inquire (file = slogfile, exist = read_vdW)

! Loop over all the atoms
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0

          ! check to see if iatom is in iatom_vdW list
          if (read_vdW) then
            if (.not. any(iatom_vdW .eq. iatom)) cycle
          end if
          r1 = s%atom(iatom)%ratom

! Loop over all possible neighbors
          do mbeta = 0, mbeta_max
            do jatom = 1, s%natoms
              ! check to see if iatom, jatom is in iatom_vdW, jatom_vdW pair
              if (read_vdW) then
                found = .false.
                do ipair = 1, npairs_vdW
                  if ((iatom .eq. iatom_vdW(ipair)) .and.                     &
     &                (jatom .eq. jatom_vdW(ipair))) found = .true.
                end do
              end if
              if (found) then
                r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
                if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.
                else

! Find the distance from (mbeta,jatom) to (0,iatom)
                  z = distance (r1, r2)
                  if (z .lt. range_vdW) then
                    num_neigh = num_neigh + 1
                  end if
                end if
              end if
            end do ! loop over jatom (the neighbor)
          end do ! mbeta over the cells
          allocate (neighbors_vdW(iatom)%neigh_b(num_neigh))
          allocate (neighbors_vdW(iatom)%neigh_j(num_neigh))
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
        end subroutine find_neigh_vdW_max


! ===========================================================================
! find_neighbors_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Finds all the neighbors and fills up the neighbor map.
! ===========================================================================
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
        integer ipair                          !< counter of pairs
        integer mbeta                          !< counter over shells
        integer num_neigh                      !< number of neighbors counter

        real z                                 !< distance between atom pair

        real, dimension (3) :: r1, r2          !< atom positions

        character (len = 25) :: slogfile

        logical found
        logical read_vdW

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
! Skip this routine if there are no van der Waals parameters
        slogfile = 'structures.vdW'
        inquire (file = slogfile, exist = read_vdW)
        if (.not. read_vdW) return

        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.vdW'
        inquire (file = slogfile, exist = read_vdW)

! Loop over all atoms in the central basis set and find his neighbours.
        do iatom = 1, s%natoms

! The variable num_neigh counts neighbors.
          num_neigh = 0
          neighbors_vdW(iatom)%neighn = num_neigh
          
          ! check to see if iatom is in iatom_vdW list
          if (read_vdW) then
            if (.not. any(iatom_vdW .eq. iatom)) cycle
          end if
          r1 = s%atom(iatom)%ratom

! Loop over all possible neighbors
          do mbeta = 0, mbeta_max
            do jatom = 1, s%natoms
              ! check to see if iatom, jatom is in iatom_vdW, jatom_vdW pair
              if (read_vdW) then
                found = .false.
                do ipair = 1, npairs_vdW
                  if ((iatom .eq. iatom_vdW(ipair)) .and.                     &
     &                (jatom .eq. jatom_vdW(ipair))) found = .true.
                end do
              end if
              if (found) then
                r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
                if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.
                else

! Find the distance from (mbeta,jatom) to (0,iatom)
                  z = distance (r1, r2)
                  if (z .lt. range_vdW) then
                    num_neigh = num_neigh + 1
                    neighbors_vdW(iatom)%neigh_j(num_neigh) = jatom
                    neighbors_vdW(iatom)%neigh_b(num_neigh) = mbeta
                  end if
                end if
              end if
            end do ! loop over jatom (the neighbor)
          end do ! mbeta over the cells
          neighbors_vdW(iatom)%neighn = num_neigh
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
! None

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

        logical iwriteout_neighbors_vdW
        logical read_vdW
        logical read_vdW_parameters

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

! Skip this routine if there are no van der Waals parameters
        slogfile = 'structures.vdW'
        inquire (file = slogfile, exist = read_vdw_parameters)
        if (.not. read_vdW_parameters) return

        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.vdW'
        inquire (file = slogfile, exist = read_vdW)

        iwriteout_neighbors_vdW = .false.
        if (iwriteout_neighbors_vdW) then
          write (logfile,*) '  '
          write (logfile,*) ' Neighbors (van der Waal''s) of each atom: '
          do iatom = 1, s%natoms
            r1 = s%atom(iatom)%ratom
            num_neigh =  neighbors_vdW(iatom)%neighn
            write (logfile,*) '  '
            write (logfile,100)
            write (logfile,101) iatom, num_neigh
            write (logfile,102)
            write (logfile,100)
            if (num_neigh .eq. 0) cycle
            do ineigh = 1, num_neigh
              mbeta = neighbors_vdW(iatom)%neigh_b(ineigh)
              jatom = neighbors_vdW(iatom)%neigh_j(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
              z = distance (r1, r2)
              write (logfile,103) iatom, ineigh, mbeta, jatom, z, r2 - r1
            end do
          end do

! Open the file NEIGHBORS which contain the neighbor map information for
! restart purposes.
          slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
          slogfile = trim(slogfile)//'.NEIGHBORS-vdW'
          open (unit = inpfile, file = slogfile, status = 'unknown')
          write (inpfile,104) s%natoms, s%basisfile
          do iatom = 1, s%natoms
            num_neigh =  neighbors_vdW(iatom)%neighn
            write (inpfile,*) iatom, num_neigh
            if (num_neigh .eq. 0) cycle
            do ineigh = 1, num_neigh
              mbeta = neighbors_vdW(iatom)%neigh_b(ineigh)
              jatom = neighbors_vdW(iatom)%neigh_j(ineigh)
              write (inpfile,*) ineigh, mbeta, jatom
            end do
          end do
          close (unit = inpfile)
        end if

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
! calculate_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Calculate the van der Waals interactions.
! ===========================================================================
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine calculate_vdW (s, vdW)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

        real vdW                                 !< van der Waal's energy

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh             !< counter over atoms and neighbors
        integer in1, in2                  !< species numbers
        integer jatom                     !< neighbor of iatom
        integer num_neigh                 !< number of neighbors
        integer mbeta                     !< the cell containing neighbor of iatom

        real alpha                        !< for making asymtpotic to R**6
        real C6factor                     !< combined two-atom C6
        real factor, dfactor              !< factor, derivative of the factor
        real Rfactor                      !< Rfactor for scaling ranges
        real vdW_piece                    !< atom pair piece of van-der Waal's
        real z                            !< distance between atom pair

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3) :: r1, r2     !< positions of iatom and jatom
        real, dimension (3) :: sighat     !< unit vector along r2 - r1

        character (len = 25) :: slogfile

        logical read_vdW_parameters

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_vdW_parameter), pointer :: pvdWi
        type(T_vdW_parameter), pointer :: pvdWj
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Skip this routine if there are no van der Waals parameters
        slogfile = 'structures.vdW'
        inquire (file = slogfile, exist = read_vdw_parameters)
        if (.not. read_vdW_parameters) return

! Loop over all atoms - van der Waal's interactions.
        vdw = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom

          ! cut some lengthy notation
          pvdWi=>vdW_parameter(in1)
          pfi=>s%forces(iatom)

! Loop over all possible neighbors
          num_neigh =  neighbors_vdW(iatom)%neighn
          if (num_neigh .eq. 0) cycle
          do ineigh = 1, num_neigh
            mbeta = neighbors_vdW(iatom)%neigh_b(ineigh)
            jatom = neighbors_vdW(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a

            ! cut some lengthy notation
            pvdWj=>vdW_parameter(in2)
            pfj=>s%forces(jatom)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! Scale the van der Waals interactions according to these factors.
! Basically, closer to the nucleus the Rfactor will have the van der Waals
! term go to zero, but asymptotically allows the interactions to be 1/R^6
! at large distances.  See Elstner et al. J. Chem. Phys. v. 114 p. 5149 (2001)
            Rfactor = (pvdWi%R0**3 + pvdWj%R0**3)/(pvdWi%R0**2 + pvdWj%R0**2)
            C6factor = 2.0d0*pvdWi%C6*pvdWi%alpha*pvdWj%C6*pvdWj%alpha        &
     &                      /(pvdWi%C6*pvdWi%alpha**2 + pvdWj%C6*pvdWj%alpha**2)

! damping term
            alpha = -3.0d0*(z/Rfactor)**7
            factor = (1.0d0 - exp(alpha))**4

! assemble the energy term
            vdw_piece = - factor*C6factor/z**6
            vdw = vdw + vdw_piece

! assemble the forces
            dfactor = (-4.0d0*exp(alpha) + 12.0d0*(exp(alpha)**2)             &
     &                 - 12.0d0*(exp(alpha)**3) + 4.0d0*(exp(alpha)**4))*7.0d0*alpha/z
            pfi%vdW = pfi%vdW + 6.0d0*vdw_piece*eta/z + dfactor*C6factor*eta/z**6
            pfj%vdW = pfj%vdW - 6.0d0*vdw_piece*eta/z - dfactor*C6factor*eta/z**6
          end do ! loop over the neighbors
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
        end subroutine calculate_vdW


! ===========================================================================
! destroy_neighbors_vdW
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the van der Waals
!! neighbors information.
!
! ===========================================================================
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

        character (len = 25) :: slogfile

        logical read_vdW
        logical read_vdW_parameters

! Procedure
! ===========================================================================
! Skip this routine if there are no van der Waals parameters
        slogfile = 'structures.vdW'
        inquire (file = slogfile, exist = read_vdw_parameters)
        if (.not. read_vdW_parameters) return

        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.vdW'
        inquire (file = slogfile, exist = read_vdW)

        do iatom = 1, s%natoms
          ! check to see if iatom is in iatom_vdW list
          ! if not, then nothing was allocated before
          if (read_vdW) then
            if (.not. any(iatom_vdW .eq. iatom)) cycle
          end if
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
        end module M_vdW
