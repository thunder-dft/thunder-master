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
        use M_configuraciones
        use M_species
        use M_kpoints
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
        integer iatom, jatom            !< counter over atoms
        integer igrid                   !< counter of energy grid points
        integer ikpoint                 !< loop over ikpoints
        integer in1, in2                !< species numbers
        integer imu                     !< counter of orbitals
        integer mmu, nnu                !< orbital number in large nxn matrix
        integer norb_mu, norb_nu        !< number atomic orbitals

! integers for matrix inversion
        integer info                    !< error information
        integer lwork                   !< size of the working arrays
        integer, allocatable :: ipiv (:)

        integer natom_initial, natom_final !< atom starting and stopping
        integer nenergy_grid            !< size of the energy grid

        real dos_atom                    !< total density of state for atom
        real dos_total                   !< total density of states for all
        real energy_min, energy_max      !< intial and final energy point
        real energy_step                 !< energy stepsize
        real eta                         !< gaussian broadening

        real pi                          !< parameter pi

! stuff for phase factor
        real dot

        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

! working arrays
        real, allocatable :: identity (:, :)

        double complex, allocatable :: greenk (:, :)
        double complex, allocatable :: Hk (:, :)

! arrays for matrix inversion
        double complex, allocatable :: work (:)

! density of states for each orbital
        real, allocatable :: dos_orbital (:)

        complex energy                      !< value of energy for igrid

        ! final result
        complex, allocatable :: green (:, :, :)

        interface
          function block_slot (natoms)
            integer, pointer :: block_slot (:)
            integer, intent(in) :: natoms
          end function block_slot
        end interface

        integer, pointer :: pblock_slot(:)

        character (len = 30) :: filename
        character (len = 25) :: slogfile

        logical read_dos

! Allocate Arrays
! ===========================================================================
        lwork = 1
        allocate (work(lwork))

        allocate (greenk (t%norbitals, t%norbitals))
        allocate (Hk (t%norbitals, t%norbitals))

! Procedure
! ===========================================================================
! Calculate the electronic density of states.
        write (logfile,*) ' Calculating the electronic density of states. '

! Determine if the dos.input files exists - if not, then exit.
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.dos.inp'
        inquire (file = slogfile, exist = read_dos)
        if (read_dos) then
! Read from input file - gives dos options
          write (logfile,*)
          write (logfile,*) ' Reading from dos.inp file! '
          open (unit = 11, file = slogfile, status = 'old')
          read (11,*) nenergy_grid    ! energy grid
          read (11,*) natom_initial, natom_final
          read (11,*) energy_min, energy_max
          read (11,*) eta  ! gaussian broadening
          close (11)
        else
          write (logfile, *) ' DOS input file non-existent! '
          return
        end if
        energy_step = (energy_max - energy_min)/nenergy_grid

! Allocate green (energy_grid)
        allocate (green (t%norbitals, t%norbitals, nenergy_grid))

! Initialize
        pblock_slot=>block_slot(t%natoms)

! Set the minimum energy point.
        energy = energy_min*a1 + eta*ai

! Read in the Hamiltonian (Lowdin transformed) for this kpoint
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.Hk'
        open (unit = 12, file = slogfile, form = 'unformatted')

! Loop over kpoints
        do ikpoint = 1, t%nkpoints
          read (12) Hk

! Loop over energy grid
          do igrid = 1, nenergy_grid
            greenk = 0.0d0

! Loop over atoms
            do iatom = 1, t%natoms
              in1 = t%atom(iatom)%imass
              norb_mu = species(in1)%norb_max
              mmu = pblock_slot(iatom)

! Loop over atoms again
              do jatom = 1, t%natoms
                in2 = t%atom(jatom)%imass
                norb_nu = species(in2)%norb_max
                nnu = pblock_slot(jatom)

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
              mmu = pblock_slot(iatom)

              do jatom = natom_initial, natom_final
                in2 = t%atom(jatom)%imass
                norb_mu = species(in2)%norb_max
                nnu = pblock_slot(jatom)

! Find the phase which is based on k*r
                vec = t%atom(jatom)%ratom - t%atom(iatom)%ratom
                sks = t%kpoints(ikpoint)%k
                dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)

                green(mmu + 1:mmu + norb_mu, nnu + 1: nnu + norb_nu, igrid) = &
      &           green(mmu + 1:mmu + norb_mu, nnu + 1: nnu + norb_nu, igrid) &
      &           + exp(-ai*dot)*t%kpoints(ikpoint)%weight                    &
      &              + greenk(mmu + 1:mmu + norb_mu, nnu + 1: nnu + norb_nu)
              end do ! end loop over jatom
            end do ! end loop over iatom

            ! step in direction of energy
            energy = energy + energy_step*a1

          end do ! end loop over energy grid
        end do ! end loop over kpoints
        close (unit = 12)

! ===========================================================================
! ---------------------------------------------------------------------------
!               W R I T E O U T    D O S
! ---------------------------------------------------------------------------
! ===========================================================================
! Initialize pi
        pi = 4.0d0*atan(1.0d0)

! Open the dos files
! Read in the Hamiltonian (Lowdin transformed) for this kpoint
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.dos'
        call system ('mkdir '//trim(slogfile)//'')

        open (unit = 23, file = trim(slogfile)//'/total.dat',                &
     &         status = 'unknown', position = 'append')

! Set the minimum energy point.
        energy = energy_min

! Loop over energy grid
        do igrid = 1, nenergy_grid
          dos_total = 0.0d0
          do iatom = natom_initial, natom_final
            in1 = t%atom(iatom)%imass
            norb_mu = species(in1)%norb_max
            mmu = pblock_slot(iatom)

            write (filename, '("/",i2.2,".",i2.2,".dat")') iatom, species(in1)%nZ
            open (unit = 24, file = trim(slogfile)//trim(filename),         &
     &              status = 'unknown', position = 'append')

            allocate (dos_orbital (norb_mu)); dos_orbital = 0.0d0
            allocate (charge_orbital (norb_mu)); charge_orbital = 0.0d0
            do imu = 1, norb_mu
              dos_orbital(imu) = -1.0d0/pi*imag(green(mmu + imu, mmu + imu, igrid))
            end do ! end loop over orbitals on iatom
            dos_atom = sum (dos_orbital)
            dos_total = dos_total + dos_atom

            write (24,*) real(energy), dos_orbital(1:norb_mu), dos_total
            deallocate (dos_orbital)
            close (unit = 24)
          end do ! end loop over iatom

          ! write to total dos file
          write (23,*) real(energy), dos_total

          ! step in direction of energy
          energy = energy + energy_step

        end do ! end loop over energy grid
        close (unit = 23)

! Deallocate Arrays
! ===========================================================================
        deallocate (work)
        deallocate (Hk)
        deallocate (green, greenk)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine dos
