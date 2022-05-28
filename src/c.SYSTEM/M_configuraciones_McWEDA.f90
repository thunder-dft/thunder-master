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

! M_configuraciones
! Module Description
! ===========================================================================
!       This is a module containing all coordinate information and atomic
! information such as the charges, etc.
!
! ===========================================================================
! Code written by:
! Ning Ma
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Module Declaration
! ===========================================================================
        module M_configuraciones

! /GLOBAL
        use M_precision
        use M_assemble_blocks

! /SYSTEM
        use M_species
        use M_kpoints

! Type Declaration
! ===========================================================================
! Define vector type for lattice vectors
        type T_vector
          real, dimension (3) :: a
        end type T_vector

! Each atom in the system has its own type.
        type T_atom
          integer imass                    ! species number

          real Q                           ! atom's total electronic charge
          real Q0                          ! atom's total nuclear charge

          real, dimension (3) :: ratom     ! atom positions
          real, dimension (3) :: ximage    ! central cell position

          real, dimension (3) :: vatom     ! atom velocities
          real, dimension (:, :), allocatable :: xdot ! atom accelerations

          type(T_species), pointer :: species  ! species specific information
          type(T_shell), pointer :: shell (:)
        end type T_atom

! The neighbor mapping type definition
        type T_neighbors
          integer neighn                         ! number of neighbors
          integer ncommon                        ! number of common neighbors

          integer, pointer :: neigh_b (:)        ! which cell is the neighbor
          integer, pointer :: neigh_j (:)        ! which atom is the neighbor
          integer, pointer :: neigh_back (:)     ! the backwards neighbor number

          ! common neighbor information - iatom and jatom are common neighbors
          ! to the third center
          integer, pointer :: iatom_common_b (:) ! cell of iatom
          integer, pointer :: iatom_common_j (:) ! atom number of iatom
          integer, pointer :: jatom_common_b (:) ! cell of jatom
          integer, pointer :: jatom_common_j (:) ! atom number of iatom
          ! neighbor number of jatom to iatom
          integer, pointer :: neigh_common (:)
        end type T_neighbors

        type T_neighbors_PP
          integer neighn                         ! number of neighbors
          integer ncommon                        ! number of common neighbors

          integer, pointer :: neigh_b (:)        ! which cell is the neighbor
          integer, pointer :: neigh_j (:)        ! which atom is the neighbor
          integer, pointer :: map (:)
          integer, pointer :: point (:)

          ! common neighbor information - iatom and jatom are common neighbors
          ! to the third center
          integer, pointer :: iatom_common_b (:) ! cell of iatom
          integer, pointer :: iatom_common_j (:) ! atom number of iatom
          integer, pointer :: jatom_common_b (:) ! cell of jatom
          integer, pointer :: jatom_common_j (:) ! atom number of iatom

          ! neighbor number of jatom to iatom
          integer, pointer :: neigh_common (:)
        end type T_neighbors_PP

! This is for forces, store forces in a type called forces.
        type T_forces
          real, dimension (3) :: overlap
          real, dimension (3) :: kinetic
          real, dimension (3) :: pulay

          real, dimension (3) :: vna
          ! vna forces - add these together to vna in the end
          real, allocatable :: vna_atom (:, :)  ! atom case
          real, allocatable :: vna_ontop (:, :) ! ontop case

          ! three-center force terms for Hartree interactions
          real, dimension (3) :: f3naa
          real, dimension (3) :: f3nab
          real, dimension (3) :: f3nac

          real, dimension (3) :: vnl
          ! vnl forces - add these together to vnl in the end
          real, allocatable :: vnl_atom (:, :)  ! atom case
          real, allocatable :: vnl_ontop (:, :) ! ontop case

          real, dimension (3) :: vxc
          ! vnl forces - add these together to vxc in the end
          real, allocatable :: vxc_off_site (:, :) ! atom case
          real, allocatable :: vxc_on_site (:, :)  ! atom case

          !ewald forces
          real, dimension (3) :: ewald
          real, allocatable :: ewaldsr (:, :)
          real, allocatable :: ewaldlr (:, :)

          ! three-center force terms for exchange-correlation interactions
          real, dimension (3) :: f3xca
          real, dimension (3) :: f3xcb
          real, dimension (3) :: f3xcc

          ! three-center components
          real, allocatable :: f3naMa (:, :, :)
          real, allocatable :: f3naMb (:, :, :)
          real, allocatable :: f3naXa (:, :, :)
          real, allocatable :: f3naXb (:, :, :)

          ! Short-range force (double-counting terms)
          real, dimension (3) :: usr

          ! Total forces
          real, dimension (3) :: ftot
        end type T_forces

! Each structure being calculated has its own type definition.
        type T_structure
          character (len = 25) :: basisfile

          integer icluster                    !< cluster calculation flag
          integer inpfile                     !< unit number for input file
          integer logfile                     !< which unit to write output
          integer natoms                      !< number of atoms
          integer nfragments                  !< number of fragments
          integer norbitals                   !< total number of orbitals
                                              !< (basically size of H and S)
          integer norbitals_new               !< new size of H and S
                                              !< (if there are linear dependencies)

          real volume                         !< cell volume
          real ztot

          type (T_vector), dimension (3) :: g !< reciprocal lattice vectors

          type(T_atom), pointer :: atom (:)   ! atom's information

          integer, pointer :: iblock_slot (:)

          ! kpoints
          integer nkpoints
          type (T_kpoint), pointer :: kpoints (:)

          ! lattice vector information
          type (T_vector), dimension (3) :: lattice
          type (T_vector), pointer :: xl (:)

          ! center-of-mass information
          real, dimension (3) :: gcm         ! geometric center
          real, dimension (3) :: rcm
          real, dimension (3) :: rcm_old
          real, dimension (3) :: vcm

          ! ********************************************************
          ! Neighbor mapping information
          ! ********************************************************
          integer, dimension (:), allocatable :: neigh_self
          integer, dimension (:), allocatable :: neighbors_PP_self
          integer, dimension (:), allocatable :: neighbors_PPx_self
          integer, dimension (:), allocatable :: neighbors_PPp_self

          type (T_neighbors), pointer :: neighbors (:)
          type (T_neighbors_PP), pointer :: neighbors_PP (:)
          type (T_neighbors_PP), pointer :: neighbors_PPx (:)
          type (T_neighbors_PP), pointer :: neighbors_PPp (:)
          ! ********************************************************

          ! ********************************************************
          ! Matrix Elements
          ! ********************************************************
          ! Put all the neighbor group belonging to the atom
          ! so in the end we have something like ME(mu, nu, ineigh, iatom)
          type (T_assemble_neighbors), pointer :: overlap (:)
          type (T_assemble_neighbors), pointer :: kinetic (:)

          ! Hartree interactions
          type (T_assemble_neighbors), pointer :: vna (:)

          ! dipole interactions for smoothing
          type (T_assemble_neighbors), pointer :: dipole_z (:)

          ! ewald interactions
          real, pointer :: ewald (:, :)
          type (T_assemble_neighbors), pointer :: ewaldsr (:)
          type (T_assemble_neighbors), pointer :: ewaldlr (:)

          ! Definition of rho_in and rho_local (used for McWEDA only):
          !         rho_bond (mu,nu) = < mu | rho_i | nu >
          !           if mu and nu are in the same atom "i" : onsite case
          !         rho_bond (mu,nu) = < mu | rho_i + rho_j | nu >
          !           if mu and nu are in different atoms "i" and "j" : atom case
          !         rho_in = sum of onsite rho_bond values
          type(T_assemble_neighbors), pointer :: rho_in (:)
          type(T_assemble_neighbors), pointer :: rho_bond (:)

          type (T_assemble_neighbors), pointer :: overlap_weighted (:)
          !       rho_in_weighted : this is Eq. (19): PRB 71, 235101 (2005)
          type (T_assemble_neighbors), pointer :: rho_in_weighted (:)
          !       rho_bond_weighted : this is Eqs. (22), (25): PRB 71, 235101 (2005)
          type (T_assemble_neighbors), pointer :: rho_bond_weighted (:)

          ! exchange-correlation
          type (T_assemble_neighbors), pointer :: vxc (:)

          ! piece due to pseudopotential
          type (T_assemble_neighbors), pointer :: svnl (:)
          type (T_assemble_neighbors), pointer :: vnl (:)

          ! Total Hamiltonian
          type (T_assemble_neighbors), pointer :: Hamiltonian (:)

          ! Density matrix
          type(T_assemble_neighbors), pointer :: denmat (:)
          type(T_assemble_neighbors), pointer :: denmat_old (:)
          type(T_assemble_neighbors), pointer :: capemat (:)
          type(T_assemble_neighbors), pointer :: denmat_PP (:)
          ! ********************************************************

          ! forces
          type(T_forces), pointer :: forces (:)
        end type T_structure

        ! options namelist
        integer iquench, iensemble
        integer iconstraint_rcm, iconstraint_vcm, iconstraint_L, iconstraint_KE
        integer ifix_neighbors, ifix_CHARGES
        integer nstepi, nstepf
        integer max_scf_iterations_set

        real dt
        real T_initial, T_final
        real T_want, taurelax                !< annealing temperature and relax time

        real efermi_T
        real scf_tolerance_set
        real beta_set
        real Ecut_set

        namelist /options/ nstepi, nstepf, iquench, iensemble, T_initial,    &
     &                     T_final, T_want, taurelax,                        &
     &                     iconstraint_rcm, iconstraint_vcm, iconstraint_L,  &
     &                     iconstraint_KE, ifix_neighbors, ifix_CHARGES,     &
     &                     efermi_T, dt, max_scf_iterations_set,             &
     &                     scf_tolerance_set, beta_set, Ecut_set

        ! output namelist
        integer iwriteout_ME_SandH
        integer iwriteout_density
        integer iwriteout_cdcoeffs
        integer iwriteout_charges
        integer iwriteout_populations
        integer iwriteout_energies
        integer iwriteout_forces
        integer iwriteout_neighbors
        integer iwriteout_dos
        integer iwriteout_abs
        integer iwriteout_xyz
        integer iwriteout_ewf

        namelist /output/ iwriteout_ME_SandH, iwriteout_density,             &
     &                    iwriteout_cdcoeffs, iwriteout_charges,             &
     &                    iwriteout_populations, iwriteout_energies,         &
     &                    iwriteout_forces, iwriteout_neighbors,             &
     &                    iwriteout_dos, iwriteout_abs, iwriteout_xyz,       &
     &                    iwriteout_ewf

! Parameter Declaration and Description
! ===========================================================================
! Number of shells for making cells
        integer, parameter :: P_mbox = 4                !< number of shells

! Parameters for scf initialization
        integer, parameter :: max_scf_iterations = 50
        real, parameter :: scf_tolerance = 1.0d-6

! beta = 0.08 means: Careful mixing (only 8 percent of the new charge).
! For larger max_order values, the choice of beta becomes less important.
        real, parameter :: beta = 0.08d0    ! factor for mixing old and new

! Ecut = 200.0d0
! This is the energy cutoff for determining the coarseness of the density grid.
        real, parameter :: Ecut = 200.0d0

! Variable Declaration and Description
! ===========================================================================
        type(T_structure), pointer :: structures (:)

        ! This is the pointer to the current structure
        type(T_structure), pointer :: s

        integer mbeta_max
        integer nstructures

! module procedures
        contains


! ===========================================================================
! read_parameters
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in all the initial parameters from file.inp
!! The file.inp is a structure filename provided by the user.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine read_parameters
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        character (len = 25) filename
        character (len = 20) string

        logical read_string, file_exists

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize parameters &OUTPUT
        iwriteout_ME_SandH = 0
        iwriteout_density = 0
        iwriteout_cdcoeffs = 0
        iwriteout_charges = 0
        iwriteout_energies = 0
        iwriteout_populations = 0
        iwriteout_forces = 0
        iwriteout_neighbors = 0
        iwriteout_dos = 0
        iwriteout_abs = 0
        iwriteout_xyz = 0
        iwriteout_ewf = 0

! Initialize parameters &OPTIONS
        nstepi = 1
        nstepf = 1
        iquench = 0
        T_initial = 300.0d0
        T_final = 0.0d0
        T_want = T_initial
        taurelax = 5.0d0
        efermi_T = 100.0d0
        dt = 0.25d0
        iensemble = 0
        iconstraint_rcm = 1
        iconstraint_vcm = 1
        iconstraint_L = 1
        iconstraint_KE = 1
        ifix_neighbors = 0
        ifix_CHARGES = 0
        max_scf_iterations_set = max_scf_iterations
        scf_tolerance_set = scf_tolerance
        beta_set = beta
        Ecut_set = Ecut

! Open structures.inp file and read global &OUTPUT options
        filename = 'structures.inp'
        INQUIRE(FILE=filename, EXIST=file_exists)   ! file_exists will be TRUE if the file                                                                    
                                                    ! exists and FALSE otherwise                                                                              
        if ( file_exists ) then
           write (*,*) 'Reading: >'//filename//'<'
        else
           write(*,*) 'ERROR: Could not open: ', filename
           stop
        end if

        string = '&OUTPUT'
        call read_sections (filename, string, read_string)
        if (read_string) then
          open (unit = 222, file = filename, status = 'old')
          read (222, nml = output)
          close (unit = 222)
        end if

        string = '&OPTIONS'
        call read_sections (filename, string, read_string)
        if (read_string) then
          open (unit = 222, file = filename, status = 'old')
          read (222, nml = options)
          close (unit = 222)
        end if

! State what is being written out
        write (ilogfile,*)
        write (ilogfile,'(A)') 'Writing options '
        write (ilogfile,'(A)') '=============== '
        write (ilogfile,*)

        if (iwriteout_ME_SandH .eq. 1)                                      &
     &    write (ilogfile,*) ' - Writing out the matrix elements for S and H. '
        if (iwriteout_density .eq. 1)                                       &
     &    write (ilogfile,*) ' - Writing out the density matrix. '
        if (iwriteout_cdcoeffs .eq. 1)                                      &
     &    write (ilogfile,*) ' - Writing out the charge density coefficients. '
        if (iwriteout_charges .eq. 1)                                       &
     &    write (ilogfile,*) ' - Writing out the charges. '
        if (iwriteout_energies .eq. 1)                                      &
     &    write (ilogfile,*) ' - Writing out the component of the energies. '
        if (iwriteout_forces .eq. 1)                                        &
     &    write (ilogfile,*) ' - Writing out the component of the forces. '
        if (iwriteout_neighbors .eq. 1)                                     &
     &    write (ilogfile,*) ' - Writing out the neighbor mapping. '
        if (iwriteout_dos .eq. 1)                                           &
     &    write (ilogfile,*) ' - Writing out the electronic density of states. '
        if (iwriteout_abs .eq. 1)                                           &
     &    write (ilogfile,*) ' Writing out the absorption spectra. '
        if (iwriteout_xyz .eq. 1)                                           &
     &    write (ilogfile,*) ' Writing out the xyz file. '
        if (iwriteout_ewf .eq. 1)                                           &
     &    write (ilogfile,*) ' Writing out the isosurfaces file. '

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine read_parameters


! ===========================================================================
! read_positions
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in all the initial atomic positions from file.inp
!! The file.inp is a structure filename provided by the user.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine read_positions (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over the atoms
        integer ikpoint                     !< counter of the kpoints
        integer in1                         !< species number
        integer inpfile                     !< reading from which unit
        integer ispecies                    !< counter over the species
        integer ix                          !< counter over x, y, and z
        integer logfile                     !< writing to which unit
        integer nZread                      !< Z atomic number read from file

        ! for shifting the system away from origin
        integer ishiftO
        real, dimension (3) :: shifter

        real xmass_total

        logical match

        interface
          function magnitude (a)
            real magnitude
            real, intent(in), dimension (3) :: a
          end function magnitude
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! Read in the positions file.
        write (logfile,'(A)') 'Atom coordinates from the structure input '
        write (logfile,'(A)') '----------------------------------------- '
        open (unit = inpfile, file = s%basisfile, status = 'old')
        read (inpfile,*) s%natoms, s%icluster
        allocate (s%atom(s%natoms))

! Read in the lattice vectors file.
        write (logfile,*)
        write (logfile,'(A)') 'Lattice vectors from structure input file '
        write (logfile,'(A)') '----------------------------------------- '
        do ix = 1, 3
          read (inpfile,*) s%lattice(ix)%a
          write (logfile,101) s%lattice(ix)%a
        end do

        call make_cells (s)

        ! Reading kpoints
        write (logfile,*)
        write (logfile,'(A)') 'K-point vectors from input file '
        write (logfile,'(A)') '------------------------------- '
        ! Reading number of kpoints
        read (inpfile,*) s%nkpoints
        allocate (s%kpoints (s%nkpoints))
        do ikpoint = 1, s%nkpoints
          read (inpfile,*) s%kpoints(ikpoint)%k, s%kpoints(ikpoint)%weight
          write (logfile,102) s%kpoints(ikpoint)%k, s%kpoints(ikpoint)%weight
        end do

        ! read atom coordinates
        do iatom = 1, s%natoms
          read (inpfile,*) nZread, s%atom(iatom)%ratom
          do ispecies = 1, nspecies
            match = .false.
            if (nZread .eq. species(ispecies)%nZ) then
              s%atom(iatom)%species=>species(ispecies)
              s%atom(iatom)%imass = ispecies
              match = .true.
              exit
            end if
          end do
          if (.not. match) then
            write (logfile,*) ' **************** ERROR !!! **************** '
            write (logfile,*) ' The atomic number, Z = ', nZread
            write (logfile,*) ' that is contained in ', s%basisfile, ' is '
            write (logfile,*) ' not contained in your info.dat file!'
            stop
          end if
        end do

        ! Calculate the center of mass
        s%rcm = 0.0d0
        xmass_total = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass_total = xmass_total + species(in1)%xmass
          s%rcm = s%rcm + species(in1)%xmass*s%atom(iatom)%ratom
        end do
        s%rcm = s%rcm/xmass_total
        write (logfile, *)
        write (logfile, *) ' System center of mass: '
        write (logfile, 101) s%rcm

        ! Find maximum number of orbitals
        s%norbitals = 0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          s%norbitals = s%norbitals + species(in1)%norb_max
        end do ! end do iatom

        ! Allocate Fermi occupation information
        do ikpoint = 1, s%nkpoints
          allocate (s%kpoints(ikpoint)%foccupy(s%norbitals))
          allocate (s%kpoints(ikpoint)%ioccupy(s%norbitals))
        end do

! Now write out the basis file information.
        write (logfile,*)
        write (logfile,*) ' Atom Coordinates from ', s%basisfile
        write (logfile,200)
        write (logfile,201)
        write (logfile,200)
        do iatom = 1, s%natoms
          write (logfile,202) iatom, s%atom(iatom)%species%symbol,          &
     &                               s%atom(iatom)%ratom, s%atom(iatom)%imass
        end do
        write (logfile,200)
        write (logfile,*) '  '

! Decide if we need to shift the atoms away from the origin.
        ishiftO = 0
        do iatom = 1, s%natoms
          if (magnitude(s%atom(iatom)%ratom) .lt. 1.0d-4) then
            ishiftO = 1
            write (logfile,*) ' Need to shift atoms away from origin. '
            exit
          end if
        end do

! Now shift coordinates to form computerese basis vectors.
        if (ishiftO .eq. 1) then
          shifter(1) = 4.0d0*atan(1.0d0)    ! pi
          shifter(2) = 1.0/exp(1.0d0)       ! 1/e
          shifter(3) = sqrt(2.0d0)          ! square root of 2
          do iatom = 1, s%natoms
            s%atom(iatom)%ratom = s%atom(iatom)%ratom + shifter
          end do
        end if

! Calculate the block_slot
        allocate (s%iblock_slot (s%natoms))
        s%iblock_slot(1) = 0
        do iatom = 2, s%natoms
          in1 = s%atom(iatom - 1)%imass
          s%iblock_slot(iatom) = s%iblock_slot(iatom - 1) + species(in1)%norb_max
        end do ! do iatom

! Format Statements
! ===========================================================================
101     format (2x, 3f10.4)
102     format (2x, 4f10.4)
200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 5x,   &
     &              ' x ', 8x, ' y ', 8x, ' z ', 6x, ' Species # ')
202     format (3x, i5, 7x, a2, 3(2x,f9.3), 7x, i2)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_positions


! ===========================================================================
! read_charges
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in all the charges from the CHARGES file.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine read_charges (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over the atoms
        integer in1                         !< species number
        integer inpfile                     !< reading from which unit
        integer issh, nssh                  !< counter over the shells
        integer logfile                     !< writing to which unit

        real Qin, Qneutral

        logical read_chargesfile

        character (len = 25) :: slogfile

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! We pass here the charges Qneutral from species to atoms
! Also, we initialize here provisionally the charges Qin and dqi
! If there is a charge file from a previous run, or the user added one,
! then initialize the input charges accordingly.
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.CHARGES'
        inquire (file = slogfile, exist = read_chargesfile)
        if (read_chargesfile) then
          write (logfile,*) ' We are reading from a charge file. '
          open (unit = inpfile, file = slogfile, status = 'old')
          write (logfile,200)
          write (logfile,401)
          write (logfile,200)
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            nssh = species(in1)%nssh
            s%atom(iatom)%Q0 = 0.0d0
            s%atom(iatom)%Q = 0.0d0
            allocate(s%atom(iatom)%shell(nssh))
            read (inpfile,*) (s%atom(iatom)%shell(issh)%Qin, issh = 1, nssh)
            write (logfile,402) iatom, s%atom(iatom)%species%symbol, nssh,   &
     &        (s%atom(iatom)%shell(issh)%Qin, issh = 1, nssh)
            do issh = 1, nssh
              Qneutral = species(in1)%shell(issh)%Qneutral
              Qin = s%atom(iatom)%shell(issh)%Qin
              s%atom(iatom)%shell(issh)%Qneutral = Qneutral
              s%atom(iatom)%shell(issh)%dQ = Qin - Qneutral
              s%atom(iatom)%Q0 = s%atom(iatom)%Q0 + Qneutral
              s%atom(iatom)%Q = s%atom(iatom)%Q + Qin
            end do
          end do
          write (logfile,200)
          close (unit = inpfile)
        else
          write (logfile,200)
          write (logfile,401)
          write (logfile,200)
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            nssh = species(in1)%nssh
            s%atom(iatom)%Q0 = 0.0d0
            s%atom(iatom)%Q = 0.0d0
            allocate(s%atom(iatom)%shell(nssh))
            do issh = 1, nssh
              Qneutral = species(in1)%shell(issh)%Qneutral
              s%atom(iatom)%shell(issh)%Qneutral = Qneutral
              s%atom(iatom)%shell(issh)%Qin = Qneutral
              s%atom(iatom)%shell(issh)%dQ = 0.0d0
              s%atom(iatom)%Q0 = s%atom(iatom)%Q0 + Qneutral
              s%atom(iatom)%Q = s%atom(iatom)%Q + Qneutral
            end do
            write (logfile,402) iatom, s%atom(iatom)%species%symbol, nssh, &
     &        (s%atom(iatom)%shell(issh)%Qin, issh = 1, nssh)
          end do
          write (logfile,200)
        end if

! Format Statements
! ===========================================================================
200     format (2x, 70('='))
401     format (2x, ' Atom # ', 2x, ' Type ', 2x, ' Shells ', 1x,' Charges ')
402     format (3x, i5, 6x, a2, 6x, i2, 4x, 10f9.3)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_charges


! ===========================================================================
! read_sections
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine reads in all the strings from the parameters in file.inp
!! The file.inp is a structure filename provided by the user. If a string is
!! provided by the user, then the option is read and initialized.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine read_sections (filename, string, flag)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        character (len = 25), intent (in) :: filename
        character (len = 20), intent (in) :: string

        logical, intent (inout) :: flag

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer eof
        logical read_string
        character (len = 20) line
        integer, parameter :: unit_id = 110

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize status
        flag = .false.

! Open .inp file
        open (unit = unit_id, file = filename, status = 'old')

! Loop through the file
        do
         read (unit_id, *, iostat = eof) line
! Compare two strings (no matter on upper & lower case)
         call compare_strings (line, string, read_string)
         if (read_string) then
           flag = .true.
           close (unit = unit_id)
           exit
         end if
! End of file
         if (eof < 0) then
           close (unit = unit_id)
           exit
         end if
       end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine read_sections


! ===========================================================================
! compare_strings
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine compares two text strings.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine compare_strings (string1, string2, flag)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        character (len = 20), intent(in) :: string1
        character (len = 20), intent(in) :: string2

! Output
        logical, intent (out) :: flag

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ilen

        character (len = 20) astring1
        character (len = 20) astring2

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize flag
        flag = .false.

! copy originals
        astring1 = string1
        astring2 = string2

! now shift lower case to upper case
        do ilen = 1, len(astring1)
          if (astring1(ilen:ilen) >= 'a' .and. astring1(ilen:ilen) <= 'z') then
            astring1(ilen:ilen) = achar (iachar (astring1(ilen:ilen)) - 32)
          end if
        end do
        do ilen = 1, len(astring2)
          if (astring2(ilen:ilen) >= 'a' .and. astring2(ilen:ilen) <= 'z') then
            astring2(ilen:ilen) = achar (iachar (astring2(ilen:ilen)) - 32)
          end if
        end do
        if (astring1 .eq. astring2) flag = .true.

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine compare_strings


! ===========================================================================
! destroy_positions
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the atomic positions
!! information - these arrays were read in by read_positions.
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
        subroutine destroy_positions ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ikpoint

! Procedure
! ===========================================================================
        do ikpoint = 1, s%nkpoints
          deallocate (s%kpoints(ikpoint)%eigen)
          deallocate (s%kpoints(ikpoint)%foccupy)
          deallocate (s%kpoints(ikpoint)%ioccupy)
          deallocate (s%kpoints(ikpoint)%c)
          deallocate (s%kpoints(ikpoint)%c_Lowdin)
        end do
        deallocate (s%kpoints)

        deallocate (structures)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_positions


! End Module
! ===========================================================================
        end module M_configuraciones
