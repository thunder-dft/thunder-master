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

! M_species
! Module Description
! ===========================================================================
!       This is a module containing all information related to the species
! - Z number, neutral charges, atomic energies, angular momentums, etc.
! This module contains the following routines:
!
!       read_Fdata_location.f90 - input location of the Fdata which is read
!       read_begin.f90 - reads information from the Z-begin.inp file
!       write_create.f90 - writes out the information for running create
!       read_create.f90 - reads in the Z-create.inp file
!       write_info.f90 - writes out the info.dat file
!       read_info.f90 - reads in the info.dat file
!       destroy_species.f90 - destroys all the arrays related to species storage
!
! ===========================================================================
! Code written by:
! @author Ning Ma
! @author James P. Lewis
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
        module M_species

        use M_precision
        implicit none

! Type Declaration
! ===========================================================================
! Atom shells
        type T_shell
          integer lssh                     ! angular momentum for shells

          real Qneutral                    ! neutral charge
          real Qneutral_ion                ! neutral charge for ion
          real Qin, Qout                   ! input charge
          real dQ                          ! input charge minus neutral charge

          ! for wavefunctions
          ! *******************************************************************
          real a0                          ! intial guess for begin SCF

          real rcutoff                     ! wavefunction cutoffs
          real rcutoffA                    ! wavefunction cutoffs in Angstroms

          real V0                          ! confinment potential value
          real r0                          ! confinment radial extension
          real V0_excited                  ! confinment potential value
          real r0_excited                  ! confinment radial extension

          character (len = 20) nafile        ! charged atom potential files
          character (len = 20) wffile        ! wave function file

          character (len = 20) nafile_excited ! charged atom potential files
          character (len = 20) wffile_excited ! wave function file
          ! *******************************************************************
        end type T_shell

! Atom orbitals
        type T_orbital
          integer issh
          integer l
          integer m
        end type T_orbital

! PP shells
        type T_shell_PP
          integer lssh                     ! angular momentum for PP shells

          real cl                          ! Kleinman-Bylander coefficient
        end type T_shell_PP

! Element information
        type T_species
          integer imass                    ! species number
          integer nZ                       ! atomic number
          integer nssh                     ! number of shells
          integer nssh_PP                  ! number of shells for PP
          integer norb_max                 ! maximum number of orbitals
          integer norb_PP_max              ! maximum number of orbitals for PP

          ! for begin
          integer ioptimize                ! ixc_option for this atom - check
          integer ixc_option               ! ixc_option for this atom - check
          integer nexcite                  ! uses excited states; 1 = yes

          real atomicE                     ! atomic energy
          real rcutoff_PP                  ! PP cutoff
          real rcutoff_PP_ion              ! PP cutoff for ion
          real xmass                       ! atomic mass
          real Zval                        ! valence charge
          real Zval_ion                    ! valence charge ion

          character (len = 12) PPfile        ! pseudopotential file
          character (len = 14) PPfile_ion    ! pseudopotential file
          character (len = 20) na0file       ! true neutral atom file
          character (len = 15) name          ! atomic name
          character (len = 2) symbol         ! atomic symbol

          ! for wavefunctions
          real rcutoff_max                 ! maximum cutoff
          real rcutoffA_max                ! maximum cutoff in Angstroms

          type(T_shell), pointer :: shell (:)
          type(T_shell_PP), pointer :: shell_PP (:)
          type(T_orbital), pointer :: orbital (:)
        end type T_species

! Variable Declaration and Description
! ===========================================================================
        integer nspecies                   ! number of species
        integer nspecies_Fdata             ! number of header lines to skip

        logical have_dorbitals             ! need this for rotations
        logical have_dorbitals_PP

        logical, pointer :: skip(:)        ! no match in Fdata

        character (len=128) Fdata_location ! directory location of Fdata
        character (len=70) signature

        type(T_species), pointer :: species (:) ! species for this simulation

! FIX ME - We will eliminate this ixc option - the options are resolved
! through object oriented linking of functions and modules.
! Read in the exchange-correlation (which is a number)
        integer ixc_option
! ***************************************************************************

! module procedures
        contains


! ===========================================================================
! read_Fdata_location
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in all the species information from the info.dat
!! file which is contained in the Fdata directory.
!
! ===========================================================================
! Code written by:
!> @author Ning Ma
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
        subroutine read_Fdata_location ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies                     !< counters for number of species
        logical :: file_exists

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Read in the species file.

        write (ilogfile,*) ' Reading: Fdata.inp  '

        INQUIRE(FILE="Fdata.inp", EXIST=file_exists)   ! file_exists will be TRUE if the file
                                                       ! exists and FALSE otherwise
        if ( file_exists ) then
           open (unit = 11, file = 'Fdata.inp', status = 'old')
        else
           write(*,*) 'ERROR: Could not open: "Fdata.inp"'
           call exit(1)
        end if

        read (11,*) nspecies
        allocate (species (nspecies))
        do ispecies = 1, nspecies
          read (11,*) species(ispecies)%nZ
        end do

        read (11,*) Fdata_location
        write (ilogfile,101) Fdata_location
        close (11)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (2x, ' Fdata is located at: ', a128)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_Fdata_location


! ===========================================================================
! read_begin
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine will read the begin.input file and define the
! variables assoiated with the type of species and exchange correlation.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine read_begin
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
        character*2 periodic (103)
        data periodic  / 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ',   &
             'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar',    &
             'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',    &
             'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',    &
             'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',    &
             'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce',    &
             'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',    &
             'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt',    &
             'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',    &
             'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',    &
             'Es', 'Fm', 'Md', 'No', 'Lw' /

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint                     ! counter over character point
        integer iten, ione, itwo, ithree   ! character array places
        integer ispecies                   ! counter over species
        integer issh                       ! counter over shells
        integer lssh                              ! l quantum number of shell
        integer nssh                       ! number of shells
        integer nzx_max                    ! what is the maximum Z
        integer nZ                         ! atomic number

        real rc                            ! radial cutoff
        real rcbohr                        ! cutoff in Bohr radii

        logical read_input

        character (len = 2) atomcheck
        character (len = 11) buffer        ! buffer for generating wavefunction
        character (len = 12) buffere       ! buffer for excited wavefunction
        character (len = 15) inputfile     ! input file for species info
        character (len = 3) rcchar

        character (len = 1), dimension (0:9) :: z

! Procedure
! ===========================================================================
! Initialize some symbols
        do ipoint = 0, 9
          z(ipoint) = char(48 + ipoint)
        end do

! We now read in a begin.input file. This determines the number of atoms
! and the types of atoms.
        write (ilogfile,*) ' We now read begin.inp '
        open (unit = 11, file = 'begin.inp', status = 'old')
        read (11, 101) signature

        nzx_max = 0
        do ispecies = 1, nspecies
          read (11,102) inputfile
          inquire (file = inputfile, exist = read_input)
          if (read_input) then
            open (unit = 12, file = inputfile, status = 'old')
          else
            write (ilogfile,*) ' The following input file does not exist! '
            write (ilogfile,102) inputfile
            stop ' Error in read_begin! '
          end if

          read (12,102) species(ispecies)%name
          read (12,103) species(ispecies)%symbol
          read (12,*) nZ
          if (nZ .ne. species(ispecies)%nZ) then
            stop ' inconsistency between Fdata.inp and create.inp in nZ '
          end if
          nzx_max = max(species(ispecies)%nZ,nzx_max)
          if (species(ispecies)%nZ .lt. nzx_max) then
            write (ilogfile,*) ' ispecies = ', ispecies
            write (ilogfile,*) ' Z(ispecies) .lt. Z(ispecies-1) '
            write (ilogfile,*) ' nzx(ispecies) = ', species(ispecies)%nZ
            write (ilogfile,*) ' nzx(ispecies-1) = ', species(ispecies-1)%nZ
            stop ' Must stop bad order.  Z1 < Z2 < Z3 ... violated'
          end if

! Check whether you put in the correct nz for that atom.
! Go through the periodic table and check.
          atomcheck = periodic(species(ispecies)%nZ)
          write (ilogfile,104) ispecies, species(ispecies)%nZ,               &
     &                 species(ispecies)%symbol, atomcheck
          if (species(ispecies)%symbol .ne. atomcheck) stop                  &
     &      ' Wrong nZ(nuc) for atom!!'

! The variable symbol is a three character string which has the atomic number.
! For example, 001 is hydrogen, 008 is oxygen, and 094 is plutonium.
          nZ = species(ispecies)%nZ
          iten = nZ/10
          if (nZ .lt. 10) buffer(1:3) = z(0)//z(0)//z(nZ)
          if (nZ .ge. 10 .and. nZ .le. 99) then
            iten = nZ/10
            ione = nZ - 10*iten
            buffer(1:3) = z(0)//z(iten)//z(ione)
          end if

! Read mass
          read (12,*) species(ispecies)%xmass

! Read number of valence electrons
          read (12,*) species(ispecies)%Zval

! For the pseudopotential filename
          species(ispecies)%PPfile(1:6) = 'basis/'
          species(ispecies)%PPfile(7:9) = buffer(1:3)
          species(ispecies)%PPfile(10:12) = '.pp'

          species(ispecies)%PPfile_ion(1:6) = 'basis/'
          species(ispecies)%PPfile_ion(7:9) = buffer(1:3)
          species(ispecies)%PPfile_ion(10:11) = '++'
          species(ispecies)%PPfile_ion(12:14) = '.pp'

! Read in the exchange-correlation (which is a number)
          read (12,*) species(ispecies)%ixc_option

! Now read stuff in related to the orbital information
          read (12,*) nssh
          species(ispecies)%nssh = nssh
          allocate (species(ispecies)%shell(nssh))

! Loop over the number of shells:
! Read the l quantum number (0, 1, 2, and 3 => s, p, d, and f), the occupation
! number, cutoff radius (in bohr), wavefunction, and neutral atom potential
! for each shell.
          species(ispecies)%rcutoff_max = -99.0d0
          do issh = 1, species(ispecies)%nssh
            read (12,*) species(ispecies)%shell(issh)%lssh
            read (12,*) species(ispecies)%shell(issh)%Qneutral
            read (12,*) species(ispecies)%shell(issh)%rcutoff

! Use the cutoff radius for each shell as part of the name of the data files.
! For instance if the cutoff for the 1st shell is 5.1234, then the name will
! have 512 as part of its name. The computer logic for this is a bit ugly, so
! breeze over this.
            rcbohr = species(ispecies)%shell(issh)%rcutoff
            species(ispecies)%rcutoff_max = max(rcbohr, species(ispecies)%rcutoff_max)
            if (rcbohr .ge. 10.0d0) rcbohr = rcbohr/10.0d0
            ione = rcbohr
            itwo = aint((rcbohr - float(ione))*10.0d0 + 0.10d0)
            ithree = aint(((rcbohr - float(ione))*10.0d0 - float(itwo))*10.0d0 + 0.10d0)
            rcchar = z(ione)//z(itwo)//z(ithree)       ! a 3 character string

! For the wavefunction filename
            buffer(4:4) = '_'
            buffer(5:7) = rcchar
            buffer(8:10) = '.wf'
            buffer(11:11) = z(issh)
            species(ispecies)%shell(issh)%wffile(1:11) = buffer

! For the orbital potential filename
            buffer(4:4) = '_'
            buffer(5:7) = rcchar
            buffer(8:10) = '.na'
            buffer(11:11) = z(issh)
            species(ispecies)%shell(issh)%nafile(1:11) = buffer
          end do

! For the neutral atom potential filename
          rcbohr = species(ispecies)%rcutoff_max
          if (rcbohr .ge. 10.0d0) rcbohr = rcbohr/10.0d0
          ione = rcbohr
          itwo = aint((rcbohr - float(ione))*10.0d0 + 0.10d0)
          ithree = aint(((rcbohr - float(ione))*10.0d0 - float(itwo))*10.0d0 + 0.10d0)
          rcchar = z(ione)//z(itwo)//z(ithree)       ! a 3 character string

          buffer(4:4) = '_'
          buffer(5:7) = rcchar
          buffer(8:10) = '.na'
          issh = 0
          buffer(11:11) = z(issh)
          species(ispecies)%na0file(1:11) = buffer

! Read in nexcite flag
          read (12,*) species(ispecies)%nexcite

          if (species(ispecies)%nexcite .ne. 0) then
            read (12,*) species(ispecies)%Zval_ion
            do issh = 1, species(ispecies)%nssh
              read (12,*) species(ispecies)%shell(issh)%Qneutral_ion

! Use the cutoff radius for each shell as part of the name of the data files.
! For instance if the cutoff for the 1st shell is 5.1234, then the name will
! have 512 as part of its name. The computer logic for this is a bit ugly, so
! breeze over this.
              rcbohr = species(ispecies)%shell(issh)%rcutoff
              if (rcbohr .ge. 10.0d0) rcbohr = rcbohr/10.0d0
              ione = rcbohr
              itwo = aint((rcbohr - float(ione))*10.0d0 + 0.10d0)
              ithree = aint(((rcbohr - float(ione))*10.0d0 - float(itwo))*10.0d0 + 0.10d0)

! rcchar is a 3 character string.
              rcchar = z(ione)//z(itwo)//z(ithree)
              buffere(1:3) = buffer(1:3)

! For the wavefunction filename (excited)
              buffere(4:4) = '_'
              buffere(5:7) = rcchar
              buffere(8:11) = '.ewf'
              buffere(12:12) = z(issh)
              species(ispecies)%shell(issh)%wffile_excited(1:12) = buffere

! For the orbital potential filename (excited)
              buffere(4:4) = '_'
              buffere(5:7) = rcchar
              buffere(8:11) = '.ena'
              buffere(12:12) = z(issh)
              species(ispecies)%shell(issh)%nafile_excited(1:12) = buffere
            end do
          end if

! Read in ioptimize flag
          read (12,*) species(ispecies)%ioptimize
          if (species(ispecies)%ioptimize .eq. 1) then
            do issh = 1, species(ispecies)%nssh
              read (12,*) species(ispecies)%shell(issh)%V0
              read (12,*) species(ispecies)%shell(issh)%r0
            end do
            if (species(ispecies)%nexcite .ne. 0) then
              do issh = 1, species(ispecies)%nssh
                read (12,*) species(ispecies)%shell(issh)%V0_excited
                read (12,*) species(ispecies)%shell(issh)%r0_excited
              end do
            end if
          end if

! What is a0?  It is the initial guess for the wavefunctions.
! We assume a hydrogenic atom, so psi = exp(-a0/r)
! Set the values here.
          do issh = 1, species(ispecies)%nssh
            lssh = species(ispecies)%shell(issh)%lssh
            if (lssh .eq. 0) species(ispecies)%shell(issh)%a0 = 2.0d0
            if (lssh .eq. 1) species(ispecies)%shell(issh)%a0 = 1.0d0
            if (lssh .eq. 2) species(ispecies)%shell(issh)%a0 = 0.8d0
          end do

! Figure out some stuff in Angstroms
          species(ispecies)%rcutoffA_max = -99.0d0
          do issh =1, species(ispecies)%nssh
            species(ispecies)%shell(issh)%rcutoffA =                        &
     &        species(ispecies)%shell(issh)%rcutoff*P_abohr
            rc = species(ispecies)%shell(issh)%rcutoffA
            species(ispecies)%rcutoffA_max = max(rc, species(ispecies)%rcutoffA_max)
          end do

! End loop over species
          close (unit = 12)
        end do
        close (unit = 11)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (a70)
102     format (a15)
103     format (a2)
104     format (2x, ' Species = ', i2, ' Nuclear Z = ', i3, ' atomname = ',  &
     &           a2, ' in periodic table = ', a2)

! ===========================================================================
        return
        end subroutine read_begin


! ===========================================================================
! write_create
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This routine writes out the species information for input files
! to run the program - create.x
! ===========================================================================
! Code written by:
! James P. Lewis
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
        subroutine write_create

        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ispecies                     ! counter over species
        integer issh                         ! counter over shells
        integer nssh                         ! number of shells

        character (len = 30) filename

! Procedure
! ===========================================================================
! Read from the create.inp file to determine which files to write out.
        write (ilogfile,*)
        write (ilogfile,*) ' Writing out the Z.inp files needed for CREATE. '
        write (ilogfile,*) ' Number of species = ', nspecies

        do ispecies = 1, nspecies
          if (species(ispecies)%symbol(2:2) .ne. ' ') then
            write (filename,'(a2,"-create.inp")') species(ispecies)%symbol
          else
            write (filename,'(a1,"-create.inp")') species(ispecies)%symbol
          end if
          write (ilogfile,*) ' filename = ', filename
          open (unit = 12, file = filename, status = 'unknown')

! Write out the atomic symbol and number
          write (12,101) species(ispecies)%symbol
          if (species(ispecies)%nZ .lt. 10) then
            write (12,102) species(ispecies)%nZ
          else
            write (12,202) species(ispecies)%nZ
          end if

! Write out the atomic mass
          if (species(ispecies)%xmass .lt. 10.0d0) then
            write (12,103) species(ispecies)%xmass
          else if (species(ispecies)%xmass .lt. 100.0d0) then
            write (12,203) species(ispecies)%xmass
          else
            write (12,303) species(ispecies)%xmass
          end if

! Write out the number of valence electrons
          write (12,104) species(ispecies)%Zval

! Write out the name of the pseudopotential file and the neutral atom
! potential file.
          write (12,105) species(ispecies)%PPfile(1:12)

          write (filename,'("basis/",a11)') species(ispecies)%na0file
          write (12,106) filename

! Write out the number of wavefunctions
          if (species(ispecies)%nexcite .ne. 0) then
            nssh = 2*species(ispecies)%nssh
          else
            nssh = species(ispecies)%nssh
          end if
          write (12,107) nssh

! Loop over the number of shells and write out shell specific information.
          do issh = 1, species(ispecies)%nssh
            write (12,108) species(ispecies)%shell(issh)%lssh
            write (12,109) species(ispecies)%shell(issh)%Qneutral
            write (12,110) species(ispecies)%shell(issh)%rcutoff

            write (filename,'("basis/",a11)') species(ispecies)%shell(issh)%wffile
            write (12,111) filename

            write (filename,'("basis/",a11)') species(ispecies)%shell(issh)%nafile
            write (12,112) filename
          end do
          if (species(ispecies)%nexcite .ne. 0) then
            do issh = 1, species(ispecies)%nssh
              write (12,113) species(ispecies)%shell(issh)%lssh
              write (12,114)
              write (12,115) species(ispecies)%shell(issh)%rcutoff

              write (filename,'("basis/",a12)') species(ispecies)%shell(issh)%wffile_excited
              write (12,116) filename

              write (filename,'("basis/",a12)') species(ispecies)%shell(issh)%nafile_excited
              write (12,117) filename
            end do
          end if

          write (12, 118)
          write (12, 119)
          write (12, 120)
          close (unit = 12)
        end do

! Write out the create.inp file
        write (ilogfile,*)
        write (ilogfile,*) ' Writing out the create.inp file needed for CREATE. '
        filename = 'create.inp'
        write (ilogfile,*) ' filename = ', filename
        open (unit = 12, file = filename, status = 'unknown')

        ! write the person responsible
        write (12,401) signature

        ! write the filenames
        do ispecies = 1, nspecies
          if (species(ispecies)%symbol(2:2) .ne. ' ') then
            write (filename,'(a2,"-create.inp")') species(ispecies)%symbol
          else
            write (filename,'(a1,"-create.inp")') species(ispecies)%symbol
          end if
          write (12,403) filename
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (a2, 36x, ' atomic symbol ')
102     format (i1, 37x, ' atomic number ')
103     format (f6.3, 32x, ' atomic mass ')
104     format (f5.2, 33x, ' number of valence electrons ')
105     format (a12, 26x, ' pseudopotential filename ')
106     format (a17, 21x, ' neutral atom filename ')
107     format (i1, 37x, ' number of shells ')
108     format (i1, 37x, ' l of the shell ')
109     format (f5.2, 33x, ' occupation number of the shell')
110     format (f4.2, 34x, ' cutoff ')
111     format (a17, 21x, ' wavefunction filename ')
112     format (a17, 21x, ' orbital potential filename ')
113     format (i1, 37x, ' l of the shell (excited)')
114     format ('0.0', 35x, ' occupation number of the (excited) shell')
115     format (f4.2, 34x, ' cutoff ')
116     format (a18, 20x, ' wavefunction filename (excited) ')
117     format (a18, 20x, ' orbital potential filename (excited) ')
118     format ('2', 37x, ' for xc - shell of changed charge ')
119     format ('0.5', 35x, ' dq of the changed shell for xc ')
120     format ('0.125  0.125')

202     format (i2, 36x, ' atomic number ')
203     format (f6.3, 32x, ' atomic mass ')
303     format (f7.3, 31x, ' atomic mass ')

401     format (a70)
402     format (i3)
403     format (a30)

        return
        end subroutine write_create


! ===========================================================================
! read_create
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine reads in all the species information from the create.inp
! file which contains all informations of cinput files.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine read_create
        implicit none

        include '../include/constants.h'

! Parameters and Data Declaration
! ===========================================================================
! FIXME - These are for the dq stuff which we will soon eliminate.
        integer, parameter :: P_nspecies = 6
        integer, parameter :: P_nssh = 3

        character*2 periodic (103)
        data periodic  / 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ',    &
             'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar',    &
             'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',    &
             'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',    &
             'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',    &
             'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce',    &
             'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',    &
             'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt',    &
             'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',    &
             'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',    &
             'Es', 'Fm', 'Md', 'No', 'Lw' /

! Argument Declaration and Description
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies                          ! counter over species
        integer issh                              ! counter over shells
        integer nssh                              ! number of shells
        integer nzx_max                           ! what is the maximum Z
        integer nZ                                ! atomic number

! FIX ME - Later we will get rid of this dq stuff!
        integer, dimension (P_nspecies) :: iderorb

        real add

        real, dimension (P_nspecies) :: dqorb
!        real, dimension (P_nssh, P_nspecies) :: dqint

        logical read_input

        character (len = 2) atomcheck
        character (len = 15) inputfile           ! input file for species info
        character (len = 70) what          ! short description of species info

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! We now read in a create.input file. This determines the number of atoms
! and the types of atoms.
        write (ilogfile,*) '  '
        write (ilogfile,*) ' We now read create.inp '
        open (unit = 11, file = 'create.inp', status = 'old')
        read (11, 101) signature

        nzx_max = 0
        do ispecies = 1, nspecies
          read (11,102) inputfile
          inquire (file = inputfile, exist = read_input)
          if (read_input) then
            open (unit = 12, file = inputfile, status = 'old')
          else
            write (*,*) ' The following input file does not exist! '
            write (*,102) inputfile
            stop ' Error in read_create! '
          end if

          read (12,103) species(ispecies)%symbol
          read (12,*) nZ
          if (nZ .ne. species(ispecies)%nZ) then
            stop ' inconsistency between Fdata.inp and create.inp in nZ '
          end if
          nzx_max = max(species(ispecies)%nZ,nzx_max)
          if (species(ispecies)%nZ .lt. nzx_max) then
            write (*,*) ' ispecies = ', ispecies
            write (*,*) ' Z(ispecies) .lt. Z(ispecies-1) '
            write (*,*) ' nzx(ispecies) = ', species(ispecies)%nZ
            write (*,*) ' nzx(ispecies-1) = ', species(ispecies-1)%nZ
            stop ' Must stop bad order.  Z1 < Z2 < Z3 ... violated'
          end if

! Check whether you put in the correct nz for that atom.
! Go through the periodic table and check.
          atomcheck = periodic(species(ispecies)%nZ)
          write (ilogfile,201) ispecies, species(ispecies)%nZ,               &
     &                 species(ispecies)%symbol, atomcheck
          if (species(ispecies)%symbol .ne. atomcheck) stop                  &
     &      ' wrong nz(nuc) for atom!!'

! Read masses
          read (12,*) species(ispecies)%xmass

! Read number of valence electrons
          read (12,*) species(ispecies)%Zval

! Read filename for the pseudopotential
          read (12,104) species(ispecies)%PPfile

! Read filename for the neutral atom potential
          read (12,105) species(ispecies)%na0file

! Now read stuff in related to the orbital information
          read (12,*) nssh
          species(ispecies)%nssh = nssh
          allocate (species(ispecies)%shell(nssh))

! Loop over the number of shells:
! Read the l quantum number (0, 1, 2, and 3 => s, p, d, and f), the occupation
! number, cutoff radius (in bohr), wavefunction, and neutral atom potential
! for each shell.
          do issh = 1, species(ispecies)%nssh
            read (12,*) species(ispecies)%shell(issh)%lssh
            read (12,*) species(ispecies)%shell(issh)%Qneutral
            read (12,*) species(ispecies)%shell(issh)%rcutoff
            read (12,105) species(ispecies)%shell(issh)%wffile
            read (12,105) species(ispecies)%shell(issh)%nafile
          end do

! Check that xns + xnp .lt. nzx
          add = 0.0d0
          do issh = 1, species(ispecies)%nssh
            add = add + species(ispecies)%shell(issh)%Qneutral
          end do
          if (int(add) .gt. species(ispecies)%Zval) then
            write (*,*) ' xnocc = ',                                         &
     &         (species(ispecies)%shell(issh)%Qneutral, issh = 1, species(ispecies)%nssh)
            write (*,*) ' Nuclear Z = ', species(ispecies)%Zval
            write (*,*) ' Sorry I must stop. How can the number of '
            write (*,*) ' electrons be larger than nuclear Z? '
            stop ' error in readcreate'
          end if

! Read information for the xc interactions - iderorb is the shell where
! the charge will be changed in the +-Q derivative stuff, and dqorb is
! the charge amount changed.
! FIXME - we need to add this to a type!
          iderorb(ispecies) = species(ispecies)%nssh
          dqorb(ispecies) = 0.5d0
          if (species(ispecies)%nssh .eq. 1) dqorb(ispecies) = 0.25d0
!         do issh = 1, species(ispecies)%nssh
!           dqint(issh,ispecies) = dqorb(ispecies)/species(ispecies)%nssh
!         end do
        end do ! end loop over species
        close (unit = 11)
        close (unit = 12)

! Set up the useful Angstrom array rcutoffa, and what.
        do ispecies = 1, nspecies
          do issh = 1, species(ispecies)%nssh
            species(ispecies)%shell(issh)%rcutoffA =                         &
     &      species(ispecies)%shell(issh)%rcutoff*P_abohr
          end do
          write (what,202) species(ispecies)%symbol, species(ispecies)%nZ,   &
     &          (species(ispecies)%shell(issh)%rcutoffA, issh = 1,           &
     &           species(ispecies)%nssh)
        end do

! First find the largest rc.
        do ispecies = 1, nspecies
          species(ispecies)%rcutoffA_max = -1.0d0
          do issh = 1, species(ispecies)%nssh
           species(ispecies)%rcutoffA_max =                                  &
     &       max(species(ispecies)%rcutoffA_max,                             &
     &           species(ispecies)%shell(issh)%rcutoffA)
          end do
        end do

! Format Statements
! ===========================================================================
101     format (a70)
102     format (a15)
103     format (a2)
104     format (a12)
105     format (a20)
201     format (2x, ' Species = ', i2, ' Nuclear Z = ', i3, ' atomname = ',  &
     &          a2, ' in periodic table = ', a2)
202     format (2x, a2, '(',i2,')', ': Rc''s (A) = ', 8f6.3)

! ==========================================================================
        return
        end subroutine read_create


! ===========================================================================
! write_info
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This routine writes out the info.dat file.
! ===========================================================================
! Code written by:
! James P. Lewis
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
        subroutine write_info

        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ispecies                     ! counter over species
        integer issh                         ! counter over shells
        integer nssh                         ! number of shells
        integer nssh_PP                      ! number of shells for PP

! Procedure
! ===========================================================================
! Next write to the info.dat file.
        open (unit = 12, file = trim(Fdata_location)//'/info.dat',           &
     &        status = 'unknown')
        write (12,101) signature
        write (12,*) nspecies, ' - Number of species '

        do ispecies = 1, nspecies
          write (12,100)
          write (12,301) ispecies
          write (12,302) species(ispecies)%symbol
          write (12,303) species(ispecies)%nZ
          write (12,304) species(ispecies)%xmass
          write (12,305) species(ispecies)%nssh
          nssh = species(ispecies)%nssh
          write (12,306) (species(ispecies)%shell(issh)%lssh, issh = 1, nssh)
          write (12,307) (species(ispecies)%shell(issh)%Qneutral, issh = 1, nssh)
          write (12,308) (species(ispecies)%shell(issh)%rcutoff, issh = 1, nssh)
          write (12,309) species(ispecies)%nssh_PP
          nssh_PP = species(ispecies)%nssh_PP
          write (12,310) (species(ispecies)%shell_PP(issh)%lssh, issh = 1, nssh_PP)
          write (12,311) species(ispecies)%rcutoff_PP
          write (12,312) (species(ispecies)%shell(issh)%wffile, issh = 1, nssh)
          write (12,313) species(ispecies)%na0file,                          &
     &                  (species(ispecies)%shell(issh)%nafile, issh = 1, nssh)
          write (12,314) species(ispecies)%atomicE
          write (12,100)
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (70('='))
101     format (2x, a70)
301     format (2x, i2, 6x, ' - Species Number and Element ')
302     format (2x, a2, 6x, ' - Elemental Symbol' )
303     format (1x, i3, 6x, ' - Nuclear Z ' )
304     format (1x, f7.3, 2x, ' - Atomic Mass ')
305     format (2x, i2, 6x, ' - Number of shells ')
306     format (3x, 9(i1, 2x), ' - L for each shell ' )
307     format (2x, 9(f5.2, 2x), ' - Occupation numbers ')
308     format (2x, 9(f5.2, 2x), ' - Radial cutoffs ')
309     format (2x, i2, 6x, ' - Number of shells (Pseudopotential) ' )
310     format (3x, 9(i1, 2x), ' - L for each shell ' )
311     format (2x, 9(f5.2, 2x),' - Radial cutoffs (Pseudopotential) ' )
312     format (2x, 9(a20, 2x))
313     format (2x, 9(a20, 2x))
314     format (2x, f12.5, 2x, ' - Atomic energy ')

        return
        end subroutine write_info


! ===========================================================================
! read_info
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in all the species information from the info.dat
!! file which is contained in the Fdata directory.
!
! ===========================================================================
! Code written by:
!> @author Ning Ma
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
        subroutine read_info ()
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: P_specieslines = 12  !< total lines of species

! Variable Declaration and Description
! ===========================================================================
        integer iline                        !< counter for number of lines
        integer ispecies                     !< counters for number of species
        integer iskip                        !< skip the header stuff in files

        integer issh                         !< counters for shells
        integer lqn                          !< angular momentum
        integer nssh                         !< number of shells
        integer nssh_PP                      !< number of shells
        integer norb_max                     !< maximum number of orbitals
        integer norb_PP_max                  !< maximum number of orbitals
        integer nZ                           !< atomic number
        integer imu , iorb                   !< orbital loops and counters

        character(len=25) signature
        character(len=2) symbol

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Read in species information from 'info.dat'
        write (ilogfile,*)
        write (ilogfile,*) ' Reading from ', trim(Fdata_location)//'/info.dat'
        open (11, file = trim(Fdata_location)//'/info.dat', status = 'old')
        read (11,102) signature
        if (signature/='') then
           write (ilogfile,*)
           write (ilogfile,*) ' You are using the database created by: '
           write (ilogfile,102) signature
        end if
        read (11,*) nspecies_Fdata
        allocate (skip(nspecies_Fdata))

        write (ilogfile,*)
        write (ilogfile,'(A)') 'Species Information '
        write (ilogfile,'(A)') '------------------- '
        write (ilogfile,*)

        ! loop over all species in Fdata
        do iskip = 1, nspecies_Fdata
          read (11,*)
          read (11,*)
          read (11,103) symbol
          read (11,*) nZ

          skip(iskip) = .true.
          do ispecies = 1, nspecies
            if (species(ispecies)%nZ .eq. nZ) then
              skip(iskip) = .false.
              exit
            end if
          end do
          if (skip(iskip)) then
            ! skip the species information
            do iline = 1, P_specieslines
              read (11,*)
            end do
          else
            ! actual read
            species(ispecies)%symbol = symbol
            species(ispecies)%nZ = nZ
            read (11,*) species(ispecies)%xmass
            read (11,*) nssh  ! valence shells
            allocate (species(ispecies)%shell(nssh))
            species(ispecies)%nssh = nssh
            read (11,*) (species(ispecies)%shell(issh)%lssh, issh = 1, nssh)
            read (11,*) (species(ispecies)%shell(issh)%Qneutral, issh = 1, nssh)
            read (11,*) (species(ispecies)%shell(issh)%rcutoff, issh = 1, nssh)
            read (11,*) nssh_PP  ! PP shells
            allocate (species(ispecies)%shell_PP(nssh_PP))
            species(ispecies)%nssh_PP = nssh_PP
            read (11,*)                                                      &
     &        (species(ispecies)%shell_PP(issh)%lssh, issh = 1, nssh_PP)
            read (11,*) species(ispecies)%rcutoff_PP
            read (11,112) (species(ispecies)%shell(issh)%wffile, issh = 1, nssh)
            read (11,113) species(ispecies)%na0file,                         &
     &                  (species(ispecies)%shell(issh)%nafile, issh = 1, nssh)
            read (11,*) species(ispecies)%atomicE
            read (11,*)

! Write out what we are actually using from the info.dat file.
!            write (ilogfile,100)
            write (ilogfile,302) species(ispecies)%symbol
            write (ilogfile,303) species(ispecies)%nZ
            write (ilogfile,304) species(ispecies)%nssh
            write (ilogfile,305) species(ispecies)%xmass
            write (ilogfile,306) species(ispecies)%atomicE
            write (ilogfile,307)                                             &
     &        (species(ispecies)%shell(issh)%lssh, issh = 1, nssh)
            write (ilogfile,308)                                             &
     &        (species(ispecies)%shell(issh)%Qneutral, issh = 1, nssh)
            write (ilogfile,309)                                             &
     &        (species(ispecies)%shell(issh)%rcutoff, issh = 1, nssh)
            write (ilogfile, '(4x, A)') '- Pseudopotential (PP):  '
            write (ilogfile,310) species(ispecies)%nssh_PP
            write (ilogfile,311)                                             &
     &        (species(ispecies)%shell_PP(issh)%lssh, issh = 1, nssh_PP)
            write (ilogfile,312) species(ispecies)%rcutoff_PP

!            write (ilogfile,100)
            write (ilogfile,*)

! Convert cutoffs from bohr radii to Angstroms
            do issh = 1, nssh
              species(ispecies)%shell(issh)%rcutoffA =                       &
      &         species(ispecies)%shell(issh)%rcutoff*P_abohr
            end do

! Calculate the maximum number of orbitals for this species
            norb_max = 0
            do issh = 1 , nssh
              lqn = species(ispecies)%shell(issh)%lssh
              norb_max = norb_max + 2*lqn + 1
            end do
            norb_PP_max = 0
            do issh = 1, nssh_PP
              lqn = species(ispecies)%shell_PP(issh)%lssh
              norb_PP_max = norb_PP_max + 2*lqn + 1
            end do
            species(ispecies)%norb_max = norb_max
            species(ispecies)%norb_PP_max = norb_PP_max

            allocate(species(ispecies)%orbital(norb_max))
            imu = 0
            do issh = 1 , nssh
            lqn = species(ispecies)%shell(issh)%lssh
              do iorb = 1, 2*lqn+1
                imu = imu + 1
                species(ispecies)%orbital(imu)%issh = issh
                species(ispecies)%orbital(imu)%l = lqn
                species(ispecies)%orbital(imu)%m = iorb - lqn -1
              end do
            end do
          end if  ! finish reading information for this species
        end do
        close (11)

! Set up logical - used later for rotations
        have_dorbitals = .false.
        have_dorbitals_PP = .false.
        do ispecies = 1, nspecies
          do issh = 1, species(ispecies)%nssh
            if (species(ispecies)%shell(issh)%lssh .eq. 2)                   &
     &        have_dorbitals = .true.
          end do
          do issh = 1, species(ispecies)%nssh_PP
            if (species(ispecies)%shell_PP(issh)%lssh .eq. 2)                &
     &        have_dorbitals_PP = .true.
          end do
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
102     format (2x, a70)
103     format (2x, a2)
112     format (2x, 9(a20, 2x))
113     format (2x, 9(a20, 2x))
302     format (4x, '- Element:          ', a3)
303     format (4x, '- Nuclear Z:        ', i3)
304     format (4x, '- Number of shells: ', i3)
305     format (4x, '- Atomic Mass:      ', f12.4)
306     format (4x, '- Atomic energy:    ', f12.4)
307     format (4x, '- L; quantum number for each shell:   ', 8(2x,i5) )
308     format (4x, '- Occupation numbers:                 ', 8(2x,f5.2) )
309     format (4x, '- Radial cutoffs:                     ', 8(2x,f5.2) )
310     format (8x, '- Number of shells (PP):                 ', i2)
311     format (8x, '- L; quantum number for each shell (PP): ', 8(2x,i1) )
312     format (8x, '- Radial cutoffs (PP):                   ', f12.4 )

! End Subroutine
! ===========================================================================
        return
        end subroutine read_info


! ===========================================================================
! destroy_species
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the species
!! information - these arrays were read in by read_species.
!
! ===========================================================================
! Code written by:
!> @author Ning Ma
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
        subroutine destroy_species ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies

! Procedure
! ===========================================================================
        do ispecies = 1, nspecies
          deallocate (species(ispecies)%shell)
          deallocate (species(ispecies)%shell_PP)
          deallocate (species(ispecies)%orbital)
        end do
        deallocate (species)
        deallocate (skip)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_species

! End Module
! ===========================================================================
        end module M_species
