! copyright info:
!
!                             @Copyright 2012
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
! ===========================================================================
! Original code written by
!> @author Cesar Gonzalez Pascual
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madird SPAIN
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
! with modification by
!> @author Pavel Jelinek
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
! dos.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine calculates the electronic density of states.
!
! ===========================================================================
        subroutine dos (t)
        use M_species
        use M_configuraciones
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom              !< counter over atoms
        integer igrid                     !< counter of energy grid points
        integer ikpoint                   !< loop over ikpoints
        integer inpfile                   !< reading from which unit
        integer in1, in2                  !< species numbers
        integer imu                       !< counter of orbitals
        integer logfile                   !< writing to which unit
        integer mmu, nnu                  !< orbital number in large nxn matrix
        integer norb_mu, norb_nu          !< number atomic orbitals

! integers for matrix inversion
        integer info                      !< error information
        integer lwork                     !< size of the working arrays
        integer, allocatable :: ipiv (:)

        integer natom_initial, natom_final !< atom starting and stopping
        integer nenergy_grid               !< size of the energy grid

        real dos_atom                      !< total density of state for atom
        real energy_min, energy_max        !< intial and final energy point
        real energy_step                   !< energy stepsize
        real eta                           !< gaussian broadening

! working arrays
        real, allocatable :: identity (:, :)

        double complex, allocatable :: greenk (:, :)
        double precision, allocatable :: Hk (:, :)

! arrays for matrix inversion
        double complex, allocatable :: work (:)

! density of states for each orbital
        real, allocatable :: dos_orbital (:)
        real, allocatable :: dos_total (:)  !< total density of states for all

        complex energy                      !< value of energy for igrid

        ! final result
        complex, allocatable :: green (:, :, :)

        character (len = 30) :: filename
        character (len = 25) :: slogfile

        logical read_dos

! Allocate Arrays
! ===========================================================================
        lwork = 1
        allocate (work(lwork))

        allocate (greenk (t%norbitals, t%norbitals))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile
        inpfile = t%inpfile

! Calculate the electronic density of states.
        write (logfile,*) 
        write (logfile,*) ' Calculating the electronic density of states. '

! Write the dos files - make the output directory
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.dos'
        call system ('mkdir '//trim(slogfile)//'')

! Determine if the dos.input files exists - if not, then exit.
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.dos.inp'
        inquire (file = slogfile, exist = read_dos)
        if (read_dos) then
! Read from input file - gives dos options
          write (logfile,*) ' Reading from dos.inp file! '
          open (unit = inpfile, file = slogfile, status = 'old')
          read (inpfile,*) nenergy_grid    ! energy grid
          read (inpfile,*) natom_initial, natom_final
          if (natom_final .gt. s%natoms) stop ' natom_final cannot be greater than natoms! '
          read (inpfile,*) energy_min, energy_max
          read (inpfile,*) eta  ! gaussian broadening
          close (inpfile)
        else
          write (logfile, *) ' DOS input file non-existent! '
          return
        end if
        energy_step = (energy_max - energy_min)/nenergy_grid

! Allocate green (energy_grid)
        allocate (green (t%norbitals, t%norbitals, nenergy_grid))

! Set the minimum energy point.
        energy = energy_min*a1 + eta*ai

! Read in the Hamiltonian (Lowdin transformed) for this kpoint
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.Hk'
        open (unit = inpfile, file = slogfile, form = 'unformatted')

! Loop over kpoints
        do ikpoint = 1, t%nkpoints
          allocate (Hk (t%norbitals, t%norbitals))
          read (inpfile) Hk

! Loop over energy grid
          do igrid = 1, nenergy_grid
            greenk = 0.0d0

! Loop over atoms
            do iatom = 1, t%natoms
              in1 = t%atom(iatom)%imass
              norb_mu = species(in1)%norb_max
              mmu = t%iblock_slot(iatom)

! Loop over atoms again
              do jatom = 1, t%natoms
                in2 = t%atom(jatom)%imass
                norb_nu = species(in2)%norb_max
                nnu = t%iblock_slot(jatom)

                allocate (identity (norb_mu, norb_nu)); identity = 0.0d0
                do imu = 1, norb_mu
                  identity (imu, imu) = 1.0d0
                end do

                if (iatom .eq. jatom) then
                  greenk(mmu + 1:mmu + norb_mu, nnu + 1: nnu + norb_nu) =      &
     &              energy*identity(1:norb_mu,1:norb_nu)                       &
     &               - Hk(mmu + 1:mmu + norb_mu, nnu + 1:nnu + norb_nu)
                else
                  greenk(mmu + 1:mmu + norb_mu, nnu + 1:nnu + norb_nu) =       &
     &               - Hk(mmu + 1:mmu + norb_mu, nnu + 1:nnu + norb_nu)
                end if
                deallocate (identity)
              end do ! end loop over jatom
            end do ! end loop over iatom

! Invert the matrix, green, send the result to green
            ! perform LU decomposition of the matrix
            allocate (ipiv (t%norbitals))
            call zgetrf (t%norbitals, t%norbitals, greenk, t%norbitals, ipiv, info)

            ! find optimal length of work
            call zgetri (t%norbitals, greenk, t%norbitals, ipiv, work, -1, info)
            lwork = work(1)
            deallocate (work)
            allocate (work (lwork))

            ! now perform actual inversion
            call zgetri (t%norbitals, greenk, t%norbitals, ipiv, work, lwork, info)
            deallocate (ipiv)
  
! Loop over atoms
            do iatom = natom_initial, natom_final
              in1 = t%atom(iatom)%imass
              norb_mu = species(in1)%norb_max
              mmu = t%iblock_slot(iatom)
  
              do jatom = natom_initial, natom_final
                in2 = t%atom(jatom)%imass
                norb_mu = species(in2)%norb_max
                nnu = t%iblock_slot(jatom)

                green(mmu + 1:mmu + norb_mu, nnu + 1: nnu + norb_nu, igrid) = &
      &           green(mmu + 1:mmu + norb_mu, nnu + 1: nnu + norb_nu, igrid) &
      &              + greenk(mmu + 1:mmu + norb_mu, nnu + 1: nnu + norb_nu)
              end do ! end loop over jatom
            end do ! end loop over iatom

            ! step in direction of energy
            energy = energy + energy_step*a1

          end do ! end loop over energy grid
        end do ! end loop over kpoints
        close (unit = inpfile)

! ===========================================================================
! ---------------------------------------------------------------------------
!               W R I T E O U T    D O S
! ---------------------------------------------------------------------------
! ===========================================================================

        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.dos'

! Loop over energy grid
        energy = energy_min
        allocate (dos_total (nenergy_grid))
        do igrid = 1, nenergy_grid
          dos_total(igrid) = 0.0d0
          do iatom = natom_initial, natom_final
            in1 = t%atom(iatom)%imass
            norb_mu = species(in1)%norb_max
            mmu = t%iblock_slot(iatom)

            write (filename, '("/",i4.4,".",i3.3,".dat")') iatom, species(in1)%nZ
            open (unit = inpfile, file = trim(slogfile)//trim(filename),     &
     &            status = 'unknown', position = 'append')

            allocate (dos_orbital (norb_mu)); dos_orbital = 0.0d0
            do imu = 1, norb_mu
              dos_orbital(imu) = -1.0d0/pi*imag(green(mmu + imu, mmu + imu, igrid))
            end do ! end loop over orbitals on iatom
            dos_atom = sum (dos_orbital)
            dos_total(igrid) = dos_total(igrid) + dos_atom

            write (inpfile,100) real(energy), dos_orbital(1:norb_mu), dos_total(igrid)
            deallocate (dos_orbital)
            close (unit = inpfile)
          end do ! end loop over iatom

          ! step in direction of energy
          energy = energy + energy_step
        end do ! end loop over energy grid

! Write out the total dos on the grid
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.dos'
        open (unit = inpfile, file = trim(slogfile)//'/total.dat', status = 'unknown')
        energy = energy_min
        do igrid = 1, nenergy_grid
          ! write to total dos file
          write (inpfile,*) real(energy), dos_total(igrid)
          energy = energy + energy_step
        end do ! end loop over energy grid
        close (unit = inpfile)

! Deallocate Arrays
! ===========================================================================
        deallocate (work)
        deallocate (green, greenk)

! Format Statements
! ===========================================================================
100     format (2x, 8f18.8)

! End Subroutine
! ===========================================================================
        return
        end subroutine dos
