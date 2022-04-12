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

! M_atomPP_functions
! Module description
!==========================================================================
! This is a module containing two subroutines which will read in wavefunctions
! and interpolate and two other subroutines which will read in the neutral
! atom potentials and interpolate. The subroutines existing in this module
! are as following:
!
!       read_vPP.f90 - read in the wavefunctions
!       vPPofr.f90 - interpolates the wavefunctions
!       calculate_vPP.f90 - calculate pseudopotential terms
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
        module M_atomPP_functions
        use M_species
        use M_atom_functions

        implicit none

! Type Declaration
! ============================================================================
! two-center interactions arrays

! To cut down on storage space, we actually change the storage procedure
! from previous FIREBALL code. Not all atoms have the same number of
! interactions or interaction types. Before - we would store things based
! on the maximum number of fdata points, maximum number of interactions types,
! maximum number of matrix elements - so even hydrogen-hydrogen (just ss
! and/or ss*) stored a 4x4 or an 8x8 matrix even when not needed.  This was
! quite inefficient.

! The new approach is to define some Fdata types which store the actual
! Fdata points. The smallest unit storage is called Fdata_cell_2C, containing
! all Fdata for a particular interaction/subinteraction type.
! ===========================================================================
        type T_shell_PP_data
          integer mesh

          real dr                           ! spacing between points
          real rcutoff                      ! rcutoff

          real, allocatable :: r (:)        ! value of r at point
          real, allocatable :: vPP (:)      ! value of the function

          integer mesh_NL

          real dr_NL                        ! spacing between points
          real rcutoff_NL                   ! rcutoff

          real, allocatable :: r_NL (:)
          real, allocatable :: vPP_NL (:)
        end type T_shell_PP_data

        type T_species_PP
          integer mesh_PP                   ! maximum mesh size for species
          integer iexc                      ! the xc used.

          real alpha
          real dr_min                       ! minimum spacing dr for species
          real rcutoffA_max                 ! maximum cutoff for this species
          real xc_fraction                  ! fraction of exact exchange added
          real Z_val                        ! valence Z number

          real, allocatable :: r_short (:)
          real, allocatable :: vPP_short (:)
          real, allocatable :: vc (:)
          real, allocatable :: vnl (:, :)

          type(T_shell_PP_data), pointer :: shell_PP (:)
        end type T_species_PP

        type(T_species_PP), pointer :: species_PP (:)

! module procedures
        contains


! ===========================================================================
! read_vPP
! ===========================================================================
! Program Description
! ===========================================================================
!        This subroutine read_vPP reads the values for the pseudopotential
! of the different (s,p,d,f) orbitals from the data file Z.pp which came
! from a pseudopotential generator. Note: The potential is already in
! Angstrom units on the data file. The variable ispec is a flag which would
! indicate in the multi-species MD code which species to read.

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
        subroutine read_vPP
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

        real rcutoff                        ! cutoff for this shell
        real Z_in                           ! atomic number

        character (len = 12) file_in        ! name of pseudo-potential file
        character (len = 12) filename       ! name of clPP file

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Open the input file
        do ispecies = 1, nspecies
          file_in = species(ispecies)%PPfile
          open (unit = 88, file = trim(Fdata_location)//'/'//trim(file_in),  &
     &          status = 'old')

          write (ilogfile,*)
          write (ilogfile,*) ' *-----------------------------------------------* '
          write (ilogfile,*) ' |               Welcome to READvPP              | '
          write (ilogfile,*) ' |  Reading pseudo-potential files of the atom   | '
          write (ilogfile,*) ' *-----------------------------------------------* '
          write (ilogfile, 101) file_in

! There are 14 message lines in each pseudopotential file
          do iline = 1, 14
           read (88,*)
          end do

! Read in which exchange-correlation approximation we are using.
          read (88,*) species_PP(ispecies)%iexc
          if (species_PP(ispecies)%iexc .eq. 12) then
            rewind (88)
            do iline = 1, 14
              read (88,*)
            end do
            read (88,*) species_PP(ispecies)%iexc, species_PP(ispecies)%xc_fraction
          end if
          iexc = species_PP(ispecies)%iexc

          write (ilogfile,*) ' We are reading in iexc from the pseudopotential file. '
          write (ilogfile,*) ' This tells us which exchange-correlation approximation '
          write (ilogfile,*) ' we are using. You have chosen iexc = ', iexc
          write (ilogfile,*)
          write (ilogfile,*) ' The available exchange-correlation options are: '
          write (ilogfile,*) ' 1  LDA Wigner'
          write (ilogfile,*) ' 2  LDA Hedin/Lundqvist'
          write (ilogfile,*) ' 3  LDA Ceperley/Alder Perdew/Zunger (1980) '
          write (ilogfile,*) ' 4  GGA Perdew/Wang (1991)'
          write (ilogfile,*) ' 5  GGA Becke (1988) X, Perdew (1986)'
          write (ilogfile,*) ' 6  GGA Perdew/Burke/Ernzerhof (1996)'
          write (ilogfile,*) ' 7  LDA Zhao/Parr'
          write (ilogfile,*) ' 8  LDA Ceperley/Alder Perdew/Wang (1991)'
          write (ilogfile,*) ' 9  GGA Becke (1988) X, Lee/Yang/Parr (1988)'
          write (ilogfile,*) ' 10 GGA Perdew/Wang (1991) X, Lee/Yang/Parr (1988)'
          write (ilogfile,*) ' 11 LDA exchange only'
          write (ilogfile,*) ' 12 GGA Becke (1988) X, Lee/Yang/Parr (1988), but '
          write (ilogfile,*) '    with mixing of exact exchange. '

! Read the number of shells
          read (88,*) species(ispecies)%nssh_PP
          allocate (species(ispecies)%shell_PP(species(ispecies)%nssh_PP))
          read (88,*) (species(ispecies)%shell_PP(issh)%lssh,                &
     &                         issh = 1, species(ispecies)%nssh_PP)
          write (ilogfile,*)
          write (ilogfile,*) ' # of pseudopotential shells = ',              &
     &                species(ispecies)%nssh_PP
          write (ilogfile,*) ' The L values = ',                             &
     &      (species(ispecies)%shell_PP(issh)%lssh,                          &
     &                                  issh = 1, species(ispecies)%nssh_PP)

! allocate shells
          allocate (species_PP(ispecies)%shell_PP(species(ispecies)%nssh_PP))

! Read Z_val
          read (88,*) Z_in
          species_PP(ispecies)%Z_val = Z_in

! Read in alpha for longranged local part => -Z*e**2*erf(alpha*r)/r
          read (88,*) species_PP(ispecies)%alpha

! Read cutoff radius of PP
          read (88,*) species(ispecies)%rcutoff_PP

! Read in the short-range local part  - this is not needed for the crtor
          write (ilogfile,*) ' Reading short-range part of pseudopotential '
          read (88,*) mesh
          species_PP(ispecies)%mesh_PP = mesh
          allocate (species_PP(ispecies)%r_short(mesh))
          allocate (species_PP(ispecies)%vPP_short(mesh))
          do ipoint = 1, mesh
            read (88,*) species_PP(ispecies)%r_short(ipoint),                &
     &                  species_PP(ispecies)%vPP_short(ipoint)
          end do
          ! Set values for dr_min and rcutoffA_max
          species_PP(ispecies)%rcutoffA_max = species_PP(ispecies)%r_short(mesh)
          species_PP(ispecies)%dr_min = species_PP(ispecies)%rcutoffA_max/real(mesh - 1)

! Read in the pseudopotential - this is not needed for the crtor
          write (ilogfile,*) ' Reading non-local part of pseudopotential '
          do issh = 1, species(ispecies)%nssh_PP
            read (88,201) mesh
            species_PP(ispecies)%shell_PP(issh)%mesh_NL = mesh
            allocate (species_PP(ispecies)%shell_PP(issh)%r_NL (mesh))
            allocate (species_PP(ispecies)%shell_PP(issh)%vPP_NL (mesh))
            do ipoint = 1, mesh
              read (88,*) species_PP(ispecies)%shell_PP(issh)%r_NL(ipoint),  &
     &                    species_PP(ispecies)%shell_PP(issh)%vPP_NL(ipoint)
            end do

            ! Set values for dr and rcutoff
            species_PP(ispecies)%shell_PP(issh)%rcutoff_NL =                 &
     &        species_PP(ispecies)%shell_PP(issh)%r_NL(mesh)
            species_PP(ispecies)%shell_PP(issh)%dr_NL =                      &
     &        species_PP(ispecies)%shell_PP(issh)%rcutoff_NL/real(mesh - 1)
          end do

! Get the Kleinman-Bylander coefficients
          write (filename,'("clPP.",i2.2,".dat")') species(ispecies)%nZ
          open (unit = 87, file = trim(Fdata_location)//'/'//trim(filename), &
     &          status = 'unknown')

! Now read in the points for the non-local part
          write (ilogfile,*) ' Reading pseudopotential '
          do issh = 1, species(ispecies)%nssh_PP
            read (88,202) mesh, species(ispecies)%shell_PP(issh)%cl
            write (87,*) species(ispecies)%shell_PP(issh)%cl
            species_PP(ispecies)%shell_PP(issh)%mesh = mesh

            ! allocate some arrays
            allocate (species_PP(ispecies)%shell_PP(issh)%r(mesh))
            allocate (species_PP(ispecies)%shell_PP(issh)%vPP(mesh))

            do ipoint = 1, mesh
              read (88,*) species_PP(ispecies)%shell_PP(issh)%r(ipoint),     &
                          species_PP(ispecies)%shell_PP(issh)%vPP(ipoint)
            end do
            species_PP(ispecies)%shell_PP(issh)%rcutoff =                    &
     &        species_PP(ispecies)%shell_PP(issh)%r(mesh)
            species_PP(ispecies)%shell_PP(issh)%dr =                         &
     &        species_PP(ispecies)%shell_PP(issh)%r(2) -                     &
     &        species_PP(ispecies)%shell_PP(issh)%r(1)
          end do
          close (unit = 87)
          close (unit = 88)

          do issh = 1, species(ispecies)%nssh_PP
            rcutoff = (mesh-1)*species_PP(ispecies)%shell_PP(issh)%dr
          end do ! issh

          if (iexc .eq. 4 .or. iexc .eq. 5 .or.                              &
              iexc .eq. 6 .or. iexc .eq. 10) then
            write (ilogfile,*) ' The exchange-correlation option that you chose '
            write (ilogfile,*) ' has not been implemented into CREATE yet. '
            write (ilogfile,*) ' Choose a different one, and restart. '
          end if
          write (ilogfile,*)
          write (ilogfile,*) ' *---------------- END READvPP -------------------*'
          write (ilogfile,*)
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (2x, a12)
201     format (12x, i5)
202     format (12x, i5, 4x, f14.7)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_vPP


! =============================================================================
! vPPofr
! =============================================================================
! Function Description
! ===========================================================================
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
        function vPPofr (r, ispecies, issh)
        implicit none

        real vPPofr

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
     &      r .ge. species_PP(ispecies)%shell_PP(issh)%rcutoff) then
          vPPofr = 0.0d0
          return
        else if (r .le. 0.0d0) then
          vPPofr = species_PP(ispecies)%shell_PP(issh)%vPP(1)
          return
        end if

! note : the points are equally spaced
        h = species_PP(ispecies)%shell_PP(issh)%dr
        imid = int(r/h) + 1
        mesh = species_PP(ispecies)%shell_PP(issh)%mesh

! This is correct, for example, when there are just not enough points to do
! an interpolation.  For instance, with hydrogen
        if (mesh < 6) then
          vPPofr = 0.0d0
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
          a(ix) = species_PP(ispecies)%shell_PP(issh)%vPP(ix+ileft)
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
          c(ix) = z(ix) - mu(ix)*c(ix+1)
          b(ix) = (a(ix+1) - a(ix))/h - h*(c(ix+1) + 2.0d0*c(ix))/3.0d0
          d(ix) = (c(ix+1) - c(ix))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp = (r - (xmin + (imid-1)*h))

        vPPofr = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function vPPofr


! ===========================================================================
! calculate_vPP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!    Subroutine will take the information read in from the specified
! pseudopotential file to calculate the potentials vc (the core potential)
! and vnl (the non-local part of the pseudopotential).
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine calculate_vPP (ispecies)
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
        alpha = species_PP(ispecies)%alpha

! We calculate vc and vnl from the values given in the ppfile.
! UNITS NOTICE :the units that are read in the .pp file are in eV-Angstrom.
! However, we want them to be in Rydberg-abohr units.
! Iterate through each wfmesh point for the local potential (long-range)
! also the erf-potential
        mesh = wf(ispecies)%mesh_max
        dr = wf(ispecies)%dr_min
        allocate (species_PP(ispecies)%vc(mesh))
        allocate (species_PP(ispecies)%vnl(species(ispecies)%nssh, mesh))

! Initialize
        species_PP(ispecies)%vc = 0.0d0
        species_PP(ispecies)%vnl = 0.0d0

! First point on the mesh
        r = - dr
        r = r + dr
        ! core potential
        species_PP(ispecies)%vc(1) = -2.0d0*Z_val*alpha*P_abohr/sqrt(pi)
        ! add to the core potential
        species_PP(ispecies)%vc(1) =                                         &
     &    species_PP(ispecies)%vc(1) + vPP_shortofr(r*P_abohr, ispecies)/P_Hartree

        ! non-local potential
        do issh = 1, species(ispecies)%nssh
          lssh = species(ispecies)%shell(issh)%lssh
          species_PP(ispecies)%vnl(issh, 1) = 0.0d0
          ! need to put the vnl into the correct shell
          do issh_PP = 1, species(ispecies)%nssh_PP
            lssh_PP = species(ispecies)%shell_PP(issh_PP)%lssh
            if (lssh_PP .eq. lssh) then
              species_PP(ispecies)%vnl(issh, 1) = vPP_NLofr(r*P_abohr, ispecies, issh_PP)/P_Hartree
            end if
          end do
        end do

! Loop over remaining mesh points
        do ipoint = 2, mesh
          r = r + dr
          ! core potentials
          species_PP(ispecies)%vc(ipoint) = - Z_val*erf(alpha*P_abohr*r)/r

          ! add to the core potential
          species_PP(ispecies)%vc(ipoint) =                                  &
     &      species_PP(ispecies)%vc(ipoint) + vPP_shortofr(r*P_abohr, ispecies)/P_Hartree

          ! non-local potential
          do issh = 1, species(ispecies)%nssh
            lssh = species(ispecies)%shell(issh)%lssh
            species_PP(ispecies)%vnl(issh, ipoint) = 0.0d0
            ! need to put the vnl into the correct shell
            do issh_PP = 1, species(ispecies)%nssh_PP
              lssh_PP = species(ispecies)%shell_PP(issh_PP)%lssh
              if (lssh_PP .eq. lssh) then
                species_PP(ispecies)%vnl(issh, ipoint) = vPP_NLofr(r*P_abohr, ispecies, issh_PP)/P_Hartree
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
        end subroutine calculate_vPP


! =============================================================================
! vPP_NLofr
! =============================================================================
! Function Description
! ===========================================================================
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
        function vPP_NLofr (r, ispecies, issh)
        implicit none

        real vPP_NLofr

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
     &      r .ge. species_PP(ispecies)%shell_PP(issh)%rcutoff_NL) then
          vPP_NLofr = 0.0d0
          return
        else if (r .le. 0.0d0) then
          vPP_NLofr = species_PP(ispecies)%shell_PP(issh)%vPP_NL(1)
          return
        end if

! note : the points are equally spaced
        h = species_PP(ispecies)%shell_PP(issh)%dr
        imid = int(r/h) + 1
        mesh = species_PP(ispecies)%shell_PP(issh)%mesh_NL

! This is correct, for example, there are just not enought points to do this
! with hydrogen
        if (mesh < 6) then
          vPP_NLofr = 0.0d0
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
          a(ix) = species_PP(ispecies)%shell_PP(issh)%vPP_NL(ix+ileft)
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

        vPP_NLofr = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function vPP_NLofr


! ===========================================================================
! vPP_shortofr
! ===========================================================================
! Function Description
! ===========================================================================
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
        function vPP_shortofr (r, ispecies)
        implicit none

        real vPP_shortofr

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
        h = species_PP(ispecies)%dr_min
        imid = int(r/h) + 1
        mesh = species_PP(ispecies)%mesh_PP

! Special cases
        if (r .ge. species_PP(ispecies)%r_short(mesh)) then
          vPP_shortofr = 0.0d0
          return
        else if (r .le. 0.0d0) then
          vPP_shortofr = species_PP(ispecies)%vPP_short(1)
          return
        end if

! This is correct, for example, there are just not enought points to do this
! with hydrogen
        if (mesh < 6) then
          vPP_shortofr = 0.0d0
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
          a(ix) = species_PP(ispecies)%vPP_short(ix+ileft)
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

        vPP_shortofr = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function vPP_shortofr


! End Module
! =============================================================================
        end module M_atomPP_functions
