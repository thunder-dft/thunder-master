! copyright info:
!
!                             @Copyright 2013
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

! M_atom_functions.f90
! Module description
!==========================================================================
! This is a module containing two subroutines which will read in wavefunctions
! and interpolate and two other subroutines which will read in the neutral
! atom potentials and interpolate. The subroutines existing in this module
! are as following:
!
!       read_wavefunctions.f90 - read in the wavefunctions
!       psiofr.f90 - interpolates the wavefunctions

!       read_napotentials.f90 - read in the potentials
!       vnaofr.f90 - interpolates the potentials

! ===========================================================================
! Code written by:
! Hong Wang
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1430 (office)
! (304) 293-5732 (FAX)
! ===========================================================================

! Module Declaration
! ===========================================================================
        module M_atom_functions
        use M_species

! Type Declaration
! ===========================================================================
        type T_data_shell
          integer mesh

          real dr
          real rcutoff
          real eigenvalue

          real, pointer :: r (:)       ! value of r at grid point
          real, pointer :: FofR (:)    ! value of the wave function

! Exact exchange potential for the single atom case - by shell
          real, pointer :: vex (:)
        end type T_data_shell

        type T_data_species
          integer mesh_max             ! maximum mesh size for species

          real dr_min                  ! minimum spacing dr for species
          real rcutoffA_max            ! maximum cutoff for this species

! Atomic energies:
          real ebs                     ! band-structure energy
          real eee                     ! electron-electron energy
          real exc                     ! exchange-correlation correction

! Charge density information
          real, pointer :: r (:)
          real, pointer :: rho (:)
          real, pointer :: sigma (:)

! Electrostatic potential for the single atom case
          real, pointer :: vee (:)

! Exchange-correlation potential for the single atom case
          real, pointer :: vxc (:)

! Neutral atom potential - shell dependencies
          real, pointer :: vna (:)

          type(T_data_shell), pointer :: shell_data (:)
        end type T_data_species

! Wavefunction information
        type(T_data_species), pointer :: wf (:)
        type(T_data_species), pointer :: wf_ion (:)

! Neutral Atom information
        type(T_data_species), pointer :: na (:)

        integer mesh                   ! number of points in file

        real dr                        ! spacing between points
        real rcutoffA                  ! cutoff for this shell (Angstrom)

        real, dimension (:), allocatable :: psitemp

! module procedures
        contains


! ===========================================================================
! read_wavefunctions
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine reads in the wavefunctions of atom.
! ===========================================================================
! Code written by:
! Hong Wang
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1430 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine read_wavefunctions
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! Data read in from the wavefunction file - compare with that in the
! create.input file read in from the user.
        integer nZ_in                       ! value of nZ from wf_wfile
        integer lssh_in                     ! angular momentum of state

        real rcutoff_in                     ! cutoff for this shell
        real rcutoff_max_in                 ! maximum cutoff for all shells
        real xnocc_in                       ! occupation number for shell

        character (len = 20) file_in        ! name of wave function file
        character (len = 15) name_in        ! input name of atom species

! Rest of the variables used:
        integer iline                       ! loop over lines in file
        integer ispecies                    ! counter over species
        integer issh                        ! counter over shells
        integer ipoint
        integer mesh_max                    ! maximum number of mesh points
        integer nssh                        ! number of states

        real rcutoff                        ! cutoff for this shell

        real r, dr                          ! value of r
        real rcutoffA_max                   ! maximum cutoff for this shell
        real xnorm                          ! value of normalization

        interface
          function simpson (n, f, dx)
            integer, intent(in) :: n
            real, intent(in) :: dx
            real, intent (in), dimension (n) :: f
            real simpson
          end function simpson
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (wf (nspecies))

! Procedure
! ============================================================================
! Loop over all the species.

        write (ilogfile,*)
        write (ilogfile,'(A)') 'Reading atomic wavefunctions '
        write (ilogfile,'(A)') '---------------------------- '

        do ispecies = 1, nspecies
          nssh = species(ispecies)%nssh
          allocate (wf(ispecies)%shell_data(nssh))
          wf(ispecies)%rcutoffA_max = -99.0d0
          wf(ispecies)%dr_min = 99.0d0

          write (ilogfile,*)
          write (ilogfile,'(A35,I3,A4,I3,A1)') '### Atomic wavefunctions for specie: ', &
     &           ispecies, ' (Z=', species(ispecies)%nZ, ')'
          write (ilogfile,*)

! Loop over the shells.
          do issh = 1, species(ispecies)%nssh

! Read in the input files
            file_in = species(ispecies)%shell(issh)%wffile
            open (unit = 12, file = trim(Fdata_location)//'/'//trim(file_in),&
     &            status = 'old', form = 'unformatted')
            read (12) file_in
            read (12) name_in, nZ_in
            read (12) lssh_in, xnocc_in
            read (12) rcutoff_in, rcutoff_max_in
            read (12) mesh
            wf(ispecies)%shell_data(issh)%mesh = mesh

!Check the input files
            if (nZ_in .ne. species(ispecies)%nZ) then
              write (*,*) ' nZ_in = ', nZ_in, ' nZ = ', species(ispecies)%nZ
              write (*,*) ' The Z number the wavefunction file, for this shell '
              write (*,*) ' does not match the cutoff radius that you put into '
              write (*,*) ' the create.inp or info.dat file. '
              write (*,*) '  '
              write (*,*) ' Please double check everything and rerun create.x'
              stop ' error in readpsi '
            end if

            rcutoff = species(ispecies)%shell(issh)%rcutoff
            if (rcutoff_in .le. (rcutoff - 1.0d-2) .or.                      &
     &          rcutoff_in .ge. (rcutoff + 1.0d-2)) then
              write (*,*) ' rcutoff_in = ', rcutoff_in, ' rcutoff = ', rcutoff
              write (*,*) ' The cutoff radius in the wavefunction file, for '
              write (*,*) ' this shell, does not match the cutoff radius that '
              write (*,*) ' you put into your create.inp or info.dat file. '
              write (*,*) '  '
              write (*,*) ' Please double check everything and rerun create.x'
              stop ' error in readpsi '
            end if

            if (xnocc_in .ne. species(ispecies)%shell(issh)%Qneutral) then
              write (*,*) ' xnocc_in = ', xnocc_in,                           &
     &                    ' xnocc = ', species(ispecies)%shell(issh)%Qneutral
              write (*,*) ' The occupation number in the wavefunction file, '
              write (*,*) ' for this shell, does not match the occupation '
              write (*,*) ' number that you put into your create.inp or   '
              write (*,*) ' info.dat file.'
              write (*,*) '  '
              write (*,*) ' Please double check everything and rerun create.x'
              stop ' error in readpsi '
            end if

            if (lssh_in .ne. species(ispecies)%shell(issh)%lssh) then
              write (*,*) ' lssh_in = ', lssh_in,                              &
     &                    ' lssh = ', species(ispecies)%shell(issh)%lssh
              write (*,*) ' The l quantum number in the wavefunction file, for '
              write (*,*) ' this shell, does not match the l quantum  number '
              write (*,*) ' that you put into your create.input file. '
              write (*,*) '  '
              write (*,*) ' Please double check everything and rerun create.x'
              stop ' error in readpsi '
            end if

! Allocate some arrays
            allocate (wf(ispecies)%shell_data(issh)%r(mesh))
            allocate (wf(ispecies)%shell_data(issh)%FofR(mesh))

! Set up calculations
            rcutoffA = species(ispecies)%shell(issh)%rcutoffA
            write (ilogfile,101) rcutoff_in, rcutoffA

            wf(ispecies)%shell_data(issh)%rcutoff = rcutoff
            wf(ispecies)%rcutoffA_max = max(wf(ispecies)%rcutoffA_max,rcutoffA)
            dr = rcutoffA/dfloat(mesh - 1)
            wf(ispecies)%shell_data(issh)%dr = dr
            wf(ispecies)%dr_min = min(wf(ispecies)%dr_min,dr)
            close (unit = 12)
          end do

! Now test the mesh and cutoffs
          mesh_max = int(wf(ispecies)%rcutoffA_max/wf(ispecies)%dr_min) + 1
          wf(ispecies)%mesh_max = mesh_max
          rcutoffA_max = dfloat(mesh_max - 1)*wf(ispecies)%dr_min
          write (ilogfile,*)
          write (ilogfile,'(4x, A)') ' Compare max cutoffs:  '
          write (ilogfile,'(4x, A31, F9.5, A2)') ' rcutoffA_max (from mesh*dr) = ', &
     &                                              rcutoffA_max, '  '
          write (ilogfile,'(4x, A31, F9.5, A2)') ' rcutoffA_max (from files)   = ', &
     &                                              wf(ispecies)%rcutoffA_max, '  '

! Loop over the shells and read in the data.
          do issh = 1, species(ispecies)%nssh
            file_in = species(ispecies)%shell(issh)%wffile
            open (unit = 12, file = trim(Fdata_location)//'/'//trim(file_in),&
     &            status = 'old', form = 'unformatted')
            do iline = 1, 4
              read (12)
            end do
            read (12) mesh
            allocate (psitemp (mesh))

! Shift the rcutoff value just a little or else the exchange-correlation
! interactions for gradients go haywire at the endpoints where rho = 0.0d0
            wf(ispecies)%rcutoffA_max = wf(ispecies)%rcutoffA_max - 1.0d-5

! Read in the points
            do ipoint = 1, mesh
              read (12) wf(ispecies)%shell_data(issh)%r(ipoint),             &
     &                  wf(ispecies)%shell_data(issh)%FofR(ipoint)
            end do
            close (unit = 12)

! Check normalization
            write (ilogfile,*)
            write (ilogfile,'(4x, A)') ' Checking normalization [NORM(l) should be 1]  '
            dr = wf(ispecies)%shell_data(issh)%dr
            r = - dr
            do ipoint = 1, mesh
              r = r + dr
              wf(ispecies)%shell_data(issh)%r(ipoint) = r
              psitemp(ipoint) = wf(ispecies)%shell_data(issh)%r(ipoint)**2   &
      &                        *wf(ispecies)%shell_data(issh)%FofR(ipoint)**2
            end do

            xnorm = simpson (mesh, psitemp, dr)
            write (ilogfile,201) issh, xnorm
            wf(ispecies)%shell_data(issh)%FofR = wf(ispecies)%shell_data(issh)%FofR/sqrt(xnorm)
            deallocate (psitemp)
          end do ! end loop over shells
        end do ! end loop over nspecies

        write (ilogfile,*)

! Format Statements
! ===========================================================================
101     format (4x, ' Rcutoff read in = ', f8.4, ' (in abohr) = ',           &
     &                                     f8.4, ' (in Angstroms)')
201     format (4x, ' NORM (shell = ', i1, ') = ', f16.12)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_wavefunctions


! =============================================================================
! psiofr
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
        function psiofr (r, rmax, ispecies, issh)
        implicit none

        real psiofr

! Argument Declaration and Description
! ===========================================================================
        real, intent (in) :: r                          ! position
        real, intent (in) :: rmax                       ! max vlaue of r

        integer, intent (in) :: ispecies, issh         ! species and shell

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
!       if (r .ge. species(ispecies)%shell(issh)%rcutoffA) then
!       if (r .ge. species(ispecies)%shell(issh)%rcutoff) then
        if (r .ge. rmax) then
          psiofr = 0.0d0
          return
        else if (r .lt. 0.0d0) then
          psiofr = wf(ispecies)%shell_data(issh)%FofR(1)
          return
        end if

! note : the points are equally spaced
        h = wf(ispecies)%shell_data(issh)%dr
        imid = int(r/h) + 1
        mesh = wf(ispecies)%shell_data(issh)%mesh

! Find starting point for the interpolation
        ileft = imid - norder/2       
        if (ileft .lt. 1) then
          ileft = 1
        else if (ileft + norder .gt. mesh) then
          ileft = mesh - norder
        end if

! Now interpolate with "natural" splines with f''(x)=0 at end points
        do ix = 0, norder
          a(ix) = wf(ispecies)%shell_data(issh)%FofR(ix + ileft)
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
          c(ix) = z(ix) - mu(ix)*c(ix+1)
          b(ix) = (a(ix+1) - a(ix))/h - h*(c(ix+1) + 2.0d0*c(ix))/3.0d0
          d(ix) = (c(ix+1) - c(ix))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp = (r - (xmin + (imid-1)*h))

        psiofr = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function psiofr


! =============================================================================
! read_napotentials
! =============================================================================
! Subroutine Description
! =============================================================================
!       This subroutine reads in the neutral atom potentials of atoms.
! Really, what we have is the true neutral atom potential in file *.na0
! and the coulomb potentials for each wave function state - s, p, d, ..., etc.
! Note: The potential is already in Angstrom units on the data file. The
! variable ispecies is a flag which would indicate in the multi-species
! MD code which species to read.

! rc must be input to this subroutine in abohr units.
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1430 (office)
! (304) 293-5732 (FAX)
!
! and Barry Haycock
! FOCAS Institute, Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine read_napotentials
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! Data read in from the wavefunction file - compare with that in the
! create.input file read in from the user.
        real rcutoff_in                     ! cutoff for this shell
        real rcutoff_max_in                 ! maximum cutoff for all shells

        character (len = 20) filename       ! name of wave function file

! Rest of the variables used:
        integer ispecies                    ! counter over species
        integer issh                        ! counter over shells
        integer ipoint
        integer mesh_max                    ! maximum number of mesh points
        integer nssh                        ! number of states

        real rcutoff
        real rcutoffA_max

! Allocate Arrays
! ===========================================================================
        allocate (na (nspecies))

! Procedure
! ===========================================================================
! Loop over all the species.
        do ispecies = 1, nspecies
          nssh = species(ispecies)%nssh
          allocate (na(ispecies)%shell_data (0:nssh))
          na(ispecies)%rcutoffA_max = -99.0d0
          na(ispecies)%dr_min = 99.0d0

          write (ilogfile,*)
          write (ilogfile,*) ' *-----------------------------------------------* '
          write (ilogfile,*) ' |              Welcome to READVNN               | '
          write (ilogfile,*) ' |  Reading (non)-neutral potential of the atoms | '
          write (ilogfile,*) ' *-----------------------------------------------* '

! Loop over the shells.
          do issh = 0, species(ispecies)%nssh    ! 0 is the pure Neutral Atom

! Open the input file
            if (issh .eq. 0) then
              filename = species(ispecies)%na0file
            else
              filename = species(ispecies)%shell(issh)%nafile
            end if
            open (unit = 12, file = trim(Fdata_location)//'/'//trim(filename),&
     &            status = 'old', form = 'unformatted')
            if (issh .eq. 0) then
              read (12) rcutoff_in
            else
              read (12) rcutoff_in, rcutoff_max_in
            end if
            read (12) mesh

! We now add the total energy of the atom to the vna data file
! (Only na, not charged parts)
            if (issh .eq. 0) then
              rcutoff =  species(ispecies)%rcutoffA_max/P_abohr
              read (12) species(ispecies)%atomicE
              write (ilogfile,*) ' Atomic total energy for species:', ispecies,&
     &                           ' is ', species(ispecies)%atomicE
            else
              rcutoff = species(ispecies)%shell(issh)%rcutoff
            end if
            write (ilogfile,*) ' issh, rcutoff = ', issh, rcutoff

! Perform some checks
            if (abs(rcutoff_in - rcutoff) .gt. 1.0d-5) then
              write (*,*) ' rcutoff_in = ', rcutoff_in, ' rcutoff = ', rcutoff
              write (*,*) ' The cutoff radius in the neutral atom file, for '
              write (*,*) ' this shell, does not match the cutoff radius '
              write (*,*) ' that you put into your create.input file. '
              write (*,*) ' Double check everything and rerun creator.'
              stop 'error in readvnn'
            end if

! Allocate some arrays
            allocate (na(ispecies)%shell_data(issh)%r(mesh))
            allocate (na(ispecies)%shell_data(issh)%FofR(mesh))

! Set up variables
            na(ispecies)%shell_data(issh)%mesh = mesh
            na(ispecies)%shell_data(issh)%rcutoff = rcutoff

            if (issh .eq. 0) then
              rcutoffA = species(ispecies)%rcutoffA_max
            else if (issh .ne. 0) then
              rcutoffA = species(ispecies)%shell(issh)%rcutoffA
            end if

            na(ispecies)%rcutoffA_max = max(na(ispecies)%rcutoffA_max, rcutoffA)
            close (12)
          end do  ! end loop over shells

! Loop over the shells.
          do issh = 0, species(ispecies)%nssh    ! 0 is the pure Neutral Atom
            dr = na(ispecies)%rcutoffA_max/dfloat(mesh - 1)
            na(ispecies)%shell_data(issh)%dr = dr
            na(ispecies)%dr_min = min(na(ispecies)%dr_min, dr)
          end do

! Now test the mesh and cutoffs
          mesh_max = int(na(ispecies)%rcutoffA_max/na(ispecies)%dr_min) + 1
          na(ispecies)%mesh_max = mesh_max
          rcutoffA_max = dfloat(mesh_max - 1)*na(ispecies)%dr_min
          write (ilogfile,*)
          write (ilogfile,*) ' rcutoffA_max (from mesh*dr) = ', rcutoffA_max
          write (ilogfile,*) ' rcutoffA_max (from files) = ', na(ispecies)%rcutoffA_max

! Loop over the shells and read in the data.
          do issh = 0, species(ispecies)%nssh
            if (issh .eq. 0) then
              filename = species(ispecies)%na0file
            else
              filename = species(ispecies)%shell(issh)%nafile
            end if
            open (unit = 12, file = trim(Fdata_location)//'/'//trim(filename),&
     &            status = 'old', form = 'unformatted')
            read (12)
            read (12) mesh
            if (issh .eq. 0) read (12)  ! line holding atomicE value

! Shift the rrc_rho value just a little or else the exchange-correlation
! interactions for gradients go haywire at the endpoints where rho = 0.0d0
            na(ispecies)%rcutoffA_max = na(ispecies)%rcutoffA_max - 1.0d-4

! Read in the points
            do ipoint = 1, mesh
              read (12) na(ispecies)%shell_data(issh)%r(ipoint),             &
                        na(ispecies)%shell_data(issh)%FofR(ipoint)
            end do
            close (unit = 12)
          end do
        end do

        write (ilogfile,*)
        write (ilogfile,*) ' *---------------- END READVNN -----------------* '
        write (ilogfile,*)

! Format Statements
! ===========================================================================
100     format (2d24.16)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_napotentials


! =============================================================================
! vnaofr
! =============================================================================
! Function Description
! ===========================================================================
!       This function returns the values vnaofr(r) (the Hartree potentials)
! for the corresponding shell of the atom described in read_atoms .
! This is done so as the Interpolator Adaptive Simpson or other, is passed a
! function that is depondent on 'r' and only 'r'.

! The value of r input to this function must be in angstrom units.
!
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Technology,
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function vnaofr (r, ispecies, issh)
        implicit none

        real vnaofr

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
! Special case
        if (issh .eq. 0 .and. r .ge. species(ispecies)%rcutoffA_max) then
          vnaofr = 0.0d0
          return
        else if (issh .ne. 0) then                                           
          if (r .gt. species(ispecies)%shell(issh)%rcutoffA) then
            vnaofr = 1.0d0/r
            return
          end if 
        else if (r .lt. 0.0d0) then
          vnaofr = na(ispecies)%shell_data(issh)%FofR(1)
          return
        end if

! note : the points are equally spaced
        h = na(ispecies)%shell_data(issh)%dr
        imid = int(r/h) + 1
        mesh = na(ispecies)%shell_data(issh)%mesh

! Find starting point for the interpolation
        ileft = imid - norder/2
        if (ileft .lt. 1) then
          ileft = 1
        else if (ileft + norder .gt. mesh) then
          ileft = mesh - norder
        end if

! Now interpolate with "natural" splines with f''(x)=0 at end points
        do ix = 0, norder
          a(ix) = na(ispecies)%shell_data(issh)%FofR(ix + ileft)
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

        vnaofr = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function vnaofr

! End Module
! ===========================================================================
        end module M_atom_functions
