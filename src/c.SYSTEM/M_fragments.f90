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

! M_fragments
! Module Description
! ===========================================================================
!>       This module
!!
!!       initialize_fragments - reads in the fragments file (if exists)
!!       destroy_fragments - destroy the neighbor arrays
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
        module M_fragments

! /GLOBAL
        use M_precision

! /SYSTEM
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
! initialize_fragments
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine will read in the fragments file from
! the basisfile.fragments file.
!
! ===========================================================================
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initialize_fragments
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
! fragments
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       Calculate the
! ===========================================================================
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine fragments (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

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
! destroy_fragments
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the fragments
!! information.
!
! ===========================================================================
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_fragments (s)
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
        end module M_fragments
