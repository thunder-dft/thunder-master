! copyright info:
!
!                             @Copyright 2013
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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

! M_atomPP_ion_functions
! Module description
!==========================================================================
! This is a module containing two subroutines which will read in wavefunctions
! and interpolate and two other subroutines which will read in the neutral
! atom potentials and interpolate. The subroutines existing in this module
! are as following:
!
!       read_vPP_ion.f90 - read in the wavefunctions
!       calculate_vPP_ion.f90 - calculate pseudopotential terms
!       vPP_NLofr.f90 - interpolate the non-local part of the pseudopotential
!       vPP_shortofr.f90 - interpolate the short range term of pseudopotential

! ============================================================================
! Code written by:
! Hong Wang
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ============================================================================
! Module Declaration
! ============================================================================
        module M_atomPP_ion_functions
        use M_species
        use M_atom_functions
        use M_atomPP_functions

        implicit none

! Type Declaration
! ============================================================================
! two-center interactions arrays

! To cut down on storage space, we actually change the storage procedure
! from previous FIREBALL code. Not all atoms have the same number of
! interactions or interaction types. Before - we would store things based
! on the maximum number of Fdata points, maximum number of interactions types,
! maximum number of matrix elements - so even hydrogen-hydrogen (just ss
! and/or ss*) stored a 4x4 or an 8x8 matrix even when not needed.  This was
! quite inefficient.

! The new approach is to define some Fdata types which store the actual
! Fdata points. The smallest unit storage is called Fdata_cell_2C, containing
! all Fdata for a particular interaction/subinteraction type.
! ===========================================================================
        type(T_species_PP), pointer :: species_PP_ion (:)

! module procedures
        contains


! =============================================================================
! psiofr_ion
! =============================================================================
! Function Description
! ===========================================================================
!       This function returns the values psiofr(r) for the corresponding
! shell of the atom described in read_atoms currently in the PSITEMP Variable.
! This is done so as the Interpolator Adaptive Simpson or other, is passed a
! function that is depondent on 'r' and only 'r'.

! The radial functions are normalized so that:

!  int ( psiofr**2  r**2  dr ) = 1.0

! The value of r input to this function must be in angstrom units.
!
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function psiofr_ion (r, ispecies, issh)
        implicit none

        real psiofr_ion

! Argument Declaration and Description
! ===========================================================================
        real, intent (in) :: r                           ! position

        integer, intent (in) :: ispecies, issh           ! species and shell

! Local Parameters and Data Declaration
! ===========================================================================
        integer norder
        parameter (norder = 5)

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid
        integer mesh

        real L(0:norder), mu(0:norder), Z(0:norder), alpha(0:norder)
        real a(0:norder), b(0:norder), c(0:norder), d(0:norder)

        integer iam
        integer ix

        real h
        real xmin
        real xxp

! Procedure
! ===========================================================================
! Special cases
        if (r .ge. species(ispecies)%shell(issh)%rcutoff) then
          psiofr_ion = 0.0d0
          return
        else if (r .lt. 0.0d0) then
          psiofr_ion = wf_ion(ispecies)%shell_data(issh)%FofR(1)
          return
        end if

! note : the points are equally spaced
        h = wf_ion(ispecies)%shell_data(issh)%dr
        imid = int(r/h) + 1
        mesh = wf_ion(ispecies)%shell_data(issh)%mesh

! Find starting point for the interpolation
        ileft = imid - norder/2
        if (ileft .lt. 1) then
          ileft = 1
        else if (ileft + norder .gt. mesh) then
          ileft = mesh - norder
        end if

! Now interpolate with "natural" splines with f''(x)=0 at end points
        do ix = 0, norder
          a(ix) = wf_ion(ispecies)%shell_data(issh)%FofR(ix + ileft)
        end do

        do ix = 1, norder - 1
          alpha(ix) = 3.0d0*(a(ix+1) - 2.0d0*a(ix) + a(ix-1))/h
        end do

        L(0) = 1
        mu(0) = 0
        Z(0) = 0
        c(0) = 0
        do ix = 1, norder - 1
          L(ix) = (4.0d0 - mu(ix-1))*h
          mu(ix) = h/L(ix)
          Z(ix) = (alpha(ix) - h*Z(ix-1))/L(ix)
        end do
        L(norder) = 1
        mu(norder) = 0
        Z(norder) = 0
        c(norder) = 0

        ! What curve section do we use?
        iam = imid - ileft

        ! Don't need 0 to iam-1
        do ix = norder - 1, iam, -1
          c(ix) = z(ix)-mu(ix)*c(ix+1)
          b(ix) = (a(ix+1) - a(ix))/h - h*(c(ix+1) + 2.0d0*c(ix))/3.0d0
          d(ix) = (c(ix+1) - c(ix))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp = (r - (xmin + (imid-1)*h))

        psiofr_ion = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function psiofr_ion


! ===========================================================================
! read_vPP_ion
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine is for the ion atomic pseudopotential.
!
!        This subroutine read_vPP_ion reads the values for the ion
! pseudopotential of the different (s,p,d,f) orbitals from the data file
! Z++.pp which came from a pseudopotential generator. Note: The potential is
! already in Angstrom units on the data file. The variable ispec is a flag
! which would indicate in the multi-species MD code which species to read.

! The cutoff (rc) must be input to this subroutine in abohr units.
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine read_vPP_ion
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iexc                        ! exchange-correlation option
        integer ispecies                    ! counter over species
        integer issh                        ! counter over shells
        integer iline                       ! number of lines to skip
        integer ipoint
        integer mesh                        ! maximum number of mesh points
        integer nssh_PP                     ! number of shells in the ion

        integer, allocatable :: lssh_PP (:) ! quantum number for shell

        real rcutoff                        ! cutoff for this shell
        real Z_in                           ! atomic number

        character (len = 14) file_in        ! name of pseudo-potential file

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Open the input file
        do ispecies = 1, nspecies
          file_in = species(ispecies)%PPfile_ion
          open (unit = 88, file = trim(Fdata_location)//'/'//trim(file_in),  &
     &          status = 'old')

          write (ilogfile,*)
          write (ilogfile,*) ' *-----------------------------------------------* '
          write (ilogfile,*) ' |           Welcome to READvPP for ion          | '
          write (ilogfile,*) ' |  Reading pseudo-potential files of the atom   | '
          write (ilogfile,*) ' *-----------------------------------------------* '
          write (ilogfile, 101) file_in

! There are 14 message lines in each pseudopotential file
          do iline = 1, 14
           read (88,*)
          end do

! Read in which exchange-correlation approximation we are using.
          read (88,*) species_PP_ion(ispecies)%iexc
          if (species_PP_ion(ispecies)%iexc .eq. 12) then
            rewind (88)
            do iline = 1, 14
              read (88,*)
            end do
            read (88,*) species_PP_ion(ispecies)%iexc, species_PP_ion(ispecies)%xc_fraction
          end if
          iexc = species_PP_ion(ispecies)%iexc
          if (iexc .ne. species_PP(ispecies)%iexc) then
            write (ilogfile,*) ' The exchange-correlation options of the Z.pp and '
            write (ilogfile,*) ' the Z++.pp files do not match! '
            stop
          end if

! Read the number of shells
          read (88,*) nssh_PP

          if (nssh_PP .ne. species(ispecies)%nssh_PP) then
            write (ilogfile,*) ' nssh_PP does not match in Z.pp and Z++.pp '
            stop
          end if
          allocate (lssh_PP (nssh_PP))
          read (88,*) (lssh_PP (issh), issh = 1, nssh_PP)
          do issh = 1, nssh_PP
            if (lssh_PP(issh) .ne. species(ispecies)%shell_PP(issh)%lssh) then
              write (ilogfile,*) ' For issh = ', issh
              write (ilogfile,*) ' lssh_PP does not match in Z.pp and Z++.pp '
              stop
            end if
          end do
          deallocate (lssh_PP)

! allocate shells
          allocate (species_PP_ion(ispecies)%shell_PP(species(ispecies)%nssh_PP))

! Read Z_val
          read (88,*) Z_in
          species_PP_ion(ispecies)%Z_val = Z_in

! Read in alpha for longranged local part => -Z*e**2*erf(alpha*r)/r
          read (88,*) species_PP_ion(ispecies)%alpha

! Read cutoff radius of PP
          read (88,*) species(ispecies)%rcutoff_PP_ion

! Read in the short-range local part  - this is not needed for the crtor
          write (ilogfile,*) ' Reading short-range part of pseudopotential '
          read (88,*) mesh
          species_PP_ion(ispecies)%mesh_PP = mesh
          allocate (species_PP_ion(ispecies)%r_short(mesh))
          allocate (species_PP_ion(ispecies)%vPP_short(mesh))
          do ipoint = 1, mesh
            read (88,*) species_PP_ion(ispecies)%r_short(ipoint),            &
     &                  species_PP_ion(ispecies)%vPP_short(ipoint)
          end do
          ! Set values for dr_min and rcutoffA_max
          species_PP_ion(ispecies)%rcutoffA_max = species_PP_ion(ispecies)%r_short(mesh)
          species_PP_ion(ispecies)%dr_min = species_PP_ion(ispecies)%rcutoffA_max/real(mesh - 1)

! Read in the pseudopotential - this is not needed for the crtor
          write (ilogfile,*) ' Reading non-local part of pseudopotential '
          do issh = 1, species(ispecies)%nssh_PP
            read (88,201) mesh
            species_PP_ion(ispecies)%shell_PP(issh)%mesh_NL = mesh
            allocate (species_PP_ion(ispecies)%shell_PP(issh)%r_NL (mesh))
            allocate (species_PP_ion(ispecies)%shell_PP(issh)%vPP_NL (mesh))
            do ipoint = 1, mesh
              read (88,*) species_PP_ion(ispecies)%shell_PP(issh)%r_NL(ipoint),&
     &                    species_PP_ion(ispecies)%shell_PP(issh)%vPP_NL(ipoint)
            end do

            ! Set values for dr and rcutoff
            species_PP_ion(ispecies)%shell_PP(issh)%rcutoff_NL =             &
     &        species_PP_ion(ispecies)%shell_PP(issh)%r_NL(mesh)
            species_PP_ion(ispecies)%shell_PP(issh)%dr_NL =                  &
     &        species_PP_ion(ispecies)%shell_PP(issh)%rcutoff_NL/real(mesh - 1)
          end do

! Now read in the points for the non-local part
          write (ilogfile,*) ' Reading pseudopotential '
          do issh = 1, species(ispecies)%nssh_PP
            read (88,202) mesh
            species_PP_ion(ispecies)%shell_PP(issh)%mesh = mesh

            ! allocate some arrays
            allocate (species_PP_ion(ispecies)%shell_PP(issh)%r(mesh))
            allocate (species_PP_ion(ispecies)%shell_PP(issh)%vPP(mesh))

            do ipoint = 1, mesh
              read (88,*) species_PP_ion(ispecies)%shell_PP(issh)%r(ipoint), &
                        species_PP_ion(ispecies)%shell_PP(issh)%vPP(ipoint)
            end do
            species_PP_ion(ispecies)%shell_PP(issh)%rcutoff =                &
     &        species_PP_ion(ispecies)%shell_PP(issh)%r(mesh)
            species_PP_ion(ispecies)%shell_PP(issh)%dr =                     &
     &        species_PP_ion(ispecies)%shell_PP(issh)%r(2) -                 &
     &        species_PP_ion(ispecies)%shell_PP(issh)%r(1)
          end do
          close (unit = 88)

          do issh = 1, species(ispecies)%nssh_PP
            rcutoff = (mesh-1)*species_PP_ion(ispecies)%shell_PP(issh)%dr
          end do ! issh

          write (ilogfile,*)
          write (ilogfile,*) ' *------------- END READvPP for ion ----------------*'
          write (ilogfile,*)
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (2x, a14)
201     format (12x, i5)
202     format (12x, i5, 4x, f14.7)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_vPP_ion


! ===========================================================================
! calculate_vPP_ion
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine is for the ion atomic pseudopotential.
!
!    Subroutine will take the information read in from the specified
! pseudopotential file of the ion atom to calculate the ion potentials vc
! (the core potential) and vnl (the non-local part of the pseudopotential).
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine calculate_vPP_ion (ispecies)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer ispecies                     !< counters for number of species

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint                   ! loop over grid points
        integer issh, issh_PP            ! loop over shells
        integer lssh, lssh_PP            ! quantum number of shells
        integer mesh                     ! mesh size for the wf

        real alpha                       ! alpha of the pesudopotential
        real dr                          ! distance between grid points
        real r                           ! value of radial point
        real Z_val                       ! valence of the pseudopotential

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================

! Initialize Z_val and alpha
        Z_val = species_PP(ispecies)%Z_val

        alpha = species_PP_ion(ispecies)%alpha

! We calculate vc and vnl from the values given in the ppfile.
! UNITS NOTICE :the units that are read in the .pp file are in eV-Angstrom.
! However, we want them to be in Rydberg-abohr units.
! Iterate through each wfmesh point for the local potential (long-range)
! also the erf-potential
        mesh = wf(ispecies)%mesh_max
        dr = wf(ispecies)%dr_min
        allocate (species_PP_ion(ispecies)%vc (mesh))
        allocate (species_PP_ion(ispecies)%vnl (species(ispecies)%nssh, mesh))

! Initialize
        species_PP_ion(ispecies)%vc = 0.0d0
        species_PP_ion(ispecies)%vnl = 0.0d0

! First point on the mesh
        r = - dr
        r = r + dr
        ! core potential
        species_PP_ion(ispecies)%vc(1) = -2.0d0*Z_val*alpha*P_abohr/sqrt(pi)
        ! add to the core potential
        species_PP_ion(ispecies)%vc(1) =                                         &
     &    species_PP_ion(ispecies)%vc(1) + vPP_ion_shortofr(r*P_abohr, ispecies)/P_Hartree

        ! non-local potential
        do issh = 1, species(ispecies)%nssh
          lssh = species(ispecies)%shell(issh)%lssh
          species_PP_ion(ispecies)%vnl(issh, 1) = 0.0d0
          ! need to put the vnl into the correct shell
          do issh_PP = 1, species(ispecies)%nssh_PP
            lssh_PP = species(ispecies)%shell_PP(issh_PP)%lssh
            if (lssh_PP .eq. lssh) then
              species_PP_ion(ispecies)%vnl(issh, 1) = vPP_ion_NLofr(r*P_abohr, ispecies, issh_PP)/P_Hartree
            end if
          end do
        end do

! Loop over remaining mesh points
        do ipoint = 2, mesh
          r = r + dr
          ! core potentials
          species_PP_ion(ispecies)%vc(ipoint) = - Z_val*erf(alpha*P_abohr*r)/r

          ! add to the core potential
          species_PP_ion(ispecies)%vc(ipoint) =                                  &
     &      species_PP_ion(ispecies)%vc(ipoint) + vPP_ion_shortofr(r*P_abohr, ispecies)/P_Hartree

          ! non-local potential
          do issh = 1, species(ispecies)%nssh
            lssh = species(ispecies)%shell(issh)%lssh
            species_PP_ion(ispecies)%vnl(issh, ipoint) = 0.0d0
            ! need to put the vnl into the correct shell
            do issh_PP = 1, species(ispecies)%nssh_PP
              lssh_PP = species(ispecies)%shell_PP(issh_PP)%lssh
              if (lssh_PP .eq. lssh) then
                species_PP_ion(ispecies)%vnl(issh, ipoint) = vPP_ion_NLofr(r*P_abohr, ispecies, issh_PP)/P_Hartree
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
        end subroutine calculate_vPP_ion


! =============================================================================
! vPP_ion_NLofr
! =============================================================================
! Function Description
! ===========================================================================
!       This routine is for the ion atomic pseudopotential.
!
!       This function returns the values vPPofr(r) (the pseudopotential)
! for the corresponding shell of the atom described in read_atoms .
! This is done so as the Interpolator Adaptive Simpson or other, is passed a
! function that is depondent on 'r' and only 'r'.

! The value of r input to this function must be in angstrom units.
!
! ===========================================================================
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
!
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function vPP_ion_NLofr (r, ispecies, issh)
        implicit none

        real vPP_ion_NLofr

! Argument Declaration and Description
! ===========================================================================
        real, intent (in) :: r                           ! position

        integer, intent (in) :: ispecies, issh           ! species and shell

! Local Parameters and Data Declaration
! ===========================================================================
        integer norder
        parameter (norder = 5)

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid
        integer mesh

        real L(0:norder), mu(0:norder), Z(0:norder), alpha(0:norder)
        real a(0:norder), b(0:norder), c(0:norder), d(0:norder)

        integer iam
        integer ix

        real h
        real xmin
        real xxp

! Procedure
! ===========================================================================
! Special cases
        if (issh .ne. 0 .and.                                                &
     &      r .ge. species_PP_ion(ispecies)%shell_PP(issh)%rcutoff_NL) then
          vPP_ion_NLofr = 0.0d0
          return
        else if (r .le. 0.0d0) then
          vPP_ion_NLofr = species_PP_ion(ispecies)%shell_PP(issh)%vPP_NL(1)
          return
        end if

! note : the points are equally spaced
        h = species_PP_ion(ispecies)%shell_PP(issh)%dr
        imid = int(r/h) + 1
        mesh = species_PP_ion(ispecies)%shell_PP(issh)%mesh

! This is correct, for example, there are just not enought points to do this
! with hydrogen
        if (mesh < 6) then
          vPP_ion_NLofr = 0.0d0
          return
        end if

! Find starting point for the interpolation
        ileft = imid - norder/2
        if (ileft .lt. 1) then
          ileft = 1
        else if (ileft + norder .gt. mesh) then
          ileft = mesh - norder
        end if

! Now interpolate with "natural" splines with f''(x)=0 at end points
        do ix = 0, norder
          a(ix) = species_PP_ion(ispecies)%shell_PP(issh)%vPP_NL(ix+ileft)
        end do

        do ix = 1, norder - 1
          alpha(ix) = 3.0d0*(a(ix+1) - 2.0d0*a(ix) + a(ix-1))/h
        end do

        L(0) = 1
        mu(0) = 0
        Z(0) = 0
        c(0) = 0
        do ix = 1, norder - 1
          L(ix) = (4.0d0 - mu(ix-1))*h
          mu(ix) = h/L(ix)
          Z(ix) = (alpha(ix) - h*Z(ix-1))/L(ix)
        end do
        L(norder) = 1
        mu(norder) = 0
        Z(norder) = 0
        c(norder) = 0

!       What curve section do we use?
        iam = imid - ileft

!       Don't need 0 to iam-1
        do ix = norder - 1, iam, -1
          c(ix) = z(ix)-mu(ix)*c(ix+1)
          b(ix) = (a(ix+1) - a(ix))/h - h*(c(ix+1) + 2.0d0*c(ix))/3.0d0
          d(ix) = (c(ix+1) - c(ix))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp = (r - (xmin + (imid-1)*h))

        vPP_ion_NLofr = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function vPP_ion_NLofr


! ===========================================================================
! vPP_ion_shortofr
! ===========================================================================
! Function Description
! ===========================================================================
!       This routine is for the ion atomic pseudopotential.
!
!       This function returns the values vPP_shortofr(r) (the short range
! piece of the pseudopotential) for the corresponding shell of the atom
! described in read_atoms. This is done so as the Interpolator Adaptive
! Simpson or other, is passed a function that is depondent on 'r' and only 'r'.

! The value of r input to this function must be in angstrom units.
!
! ===========================================================================
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
!
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function vPP_ion_shortofr (r, ispecies)
        implicit none

        real vPP_ion_shortofr

! Argument Declaration and Description
! ===========================================================================
        real, intent (in) :: r                           ! position

        integer, intent (in) :: ispecies                 ! species and shell

! Local Parameters and Data Declaration
! ===========================================================================
        integer norder
        parameter (norder = 5)

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid
        integer mesh

        real L(0:norder), mu(0:norder), Z(0:norder), alpha(0:norder)
        real a(0:norder), b(0:norder), c(0:norder), d(0:norder)

        integer iam
        integer ix

        real h
        real xmin
        real xxp

! Procedure
! ===========================================================================
! Note : the points are equally spaced
        h = species_PP_ion(ispecies)%dr_min
        imid = int(r/h) + 1
        mesh = species_PP_ion(ispecies)%mesh_PP

! Special cases
        if (r .ge. species_PP_ion(ispecies)%r_short(mesh)) then
          vPP_ion_shortofr = 0.0d0
          return
        else if (r .le. 0.0d0) then
          vPP_ion_shortofr = species_PP_ion(ispecies)%vPP_short(1)
          return
        end if

! This is correct, for example, there are just not enought points to do this
! with hydrogen
        if (mesh < 6) then
          vPP_ion_shortofr = 0.0d0
          return
        end if

! Find starting point for the interpolation
        ileft = imid - norder/2
        if (ileft .lt. 1) then
          ileft = 1
        else if (ileft + norder .gt. mesh) then
          ileft = mesh - norder
        end if

! Now interpolate with "natural" splines with f''(x)=0 at end points
        do ix = 0, norder
          a(ix) = species_PP_ion(ispecies)%vPP_short(ix+ileft)
        end do

        do ix = 1, norder - 1
          alpha(ix) = 3.0d0*(a(ix+1) - 2.0d0*a(ix) + a(ix-1))/h
        end do

        L(0) = 1
        mu(0) = 0
        Z(0) = 0
        c(0) = 0
        do ix = 1, norder - 1
          L(ix) = (4.0d0 - mu(ix-1))*h
          mu(ix) = h/L(ix)
          Z(ix) = (alpha(ix) - h*Z(ix-1))/L(ix)
        end do
        L(norder) = 1
        mu(norder) = 0
        Z(norder) = 0
        c(norder) = 0

!       What curve section do we use?
        iam = imid - ileft

!       Don't need 0 to iam-1
        do ix = norder - 1, iam, -1
          c(ix) = z(ix)-mu(ix)*c(ix+1)
          b(ix) = (a(ix+1) - a(ix))/h - h*(c(ix+1) + 2.0d0*c(ix))/3.0d0
          d(ix) = (c(ix+1) - c(ix))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp = (r - (xmin + (imid-1)*h))

        vPP_ion_shortofr = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function vPP_ion_shortofr


! End Module
! =============================================================================
        end module M_atomPP_ion_functions
