! copyright info:
!
!                             @Copyright 2014
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


! ===========================================================================
! absorption.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine calculates the absorption probabilities.
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
!
! Declaration
! ===========================================================================
        subroutine absorption (t)
        use M_species
        use M_configuraciones
        use M_atom_functions
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter, dimension (3) :: m_dipole = reshape ([-1, 0, 1], [3])

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer iband                   !< loop over bands
        integer ihomo                   !< highest occupied molecular orbital
        integer ihomo_begin             !< which homo level reading
        integer ilevel                  !< what current energy level being calculated
        integer in1, in2                !< species numbers
        integer inpfile                 !< reading from which unit
        integer jatom                   !< neighbor of iatom
        integer logfile                 !< writing to which unit
        integer nhomo                   !< number of homo states
        integer num_neigh               !< number of neighbors
        integer mbeta                   !< the cell containing neighbor of iatom

        integer itransition             !< loop over transitions
        integer itransition_begin       !< starting point of loop (for restart)

        integer igrid                   !< loop over the energy grid
        integer nenergy_grid            !< size of the energy grid

        integer ir             ! counter over dipole directions
        integer ix, nnx        ! counter over x-direction grid points up to nnx
        integer iy, nny        ! counter over y-direction grid points up to nny
        integer iz, nnz        ! counter over z-direction grid points up to nnz

        integer l, l1, l2                         ! quantum numbers
        integer m, m1, m2
        integer n1, n2

        integer imu, inu                  ! loop over shells
        integer mmu, nnu                  ! location within the coefficient vector
        integer norb_mu, norb_nu          ! number atomic orbitals

        real dgrid                         ! grid spacing for integration
        real energy                        !< value of the nergy at a point
        real energy_min, energy_max        !< intial and final energy point
        real energy_step                   !< energy stepsize
        real eta                           !< gaussian broadening
        real homo_begin                    !< from which energy are we beginning
        real r1pin, r2pin          ! distance from centers to integration point
        real r1max, r2max          ! maximum distance of the wavefunctions

        real x, xmin, xmax
        real y, ymin, ymax
        real z, zmin, zmax

        real xtemp                 ! temporary integral sum
        real xntegral, yntegral, zntegral

        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: r1p, r2p !< shifted positions of iatom and jatom
        real, dimension (3) :: r_dipole !< dipole value at integration point

        ! Simpsons Quadrature Variables.
        real, allocatable, dimension (:) :: xmult, ymult, zmult

        ! total absorption integral result
        real xntensity                                !< result on energy grid
        real, dimension (:, :), allocatable :: answer    !< result for each transition

        character (len = 25) :: slogfile

        logical read_abs

        interface
          function Ylm (r, l, m)
            integer, intent (in) :: l, m
            real, intent (in), dimension (3) :: r
            real Ylm
          end function Ylm
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! FIXME - for now just do cluster integration
        t%icluster = 1

! Initialize logfile
        logfile = t%logfile
        inpfile = t%inpfile

! Calculate the absorption probabilities.
        write (logfile,*)
        write (logfile,*) ' Calculating the absorption probabilities. '

! Determine if the abs.inp files exists - if not, then exit.
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.abs.inp'
        inquire (file = slogfile, exist = read_abs)
        if (read_abs) then
! Read from input file - gives dos options
          write (logfile,*) ' Reading from abs.inp file! '
          open (unit = inpfile, file = slogfile, status = 'old')
          read (inpfile,*) nenergy_grid    ! energy grid
          read (inpfile,*) dgrid    ! read grid spacing for integrations
          read (inpfile,*) energy_min, energy_max  ! subtracted/added from HOMO)
          read (inpfile,*) eta  ! gaussian broadening
          close (inpfile)
        else
          write (logfile,*) ' There is no absorption input file. '
          return
        end if
        write (logfile,*)

        ! renormalize eta
        eta = eta/2.0d0

! Determine the number of transitions we will make based on the energy range.
! Loop over energy grid
        if (modulo(int(t%ztot), 2) .eq. 0) then
          ihomo = int(t%ztot)/2
        else
          ihomo = int(t%ztot)/2 + 1
        end if
        energy_min = t%kpoints(1)%eigen(ihomo) + energy_min
        energy_step = energy_max/nenergy_grid

! Find out how many "HOMO" states from which we are calculating transitions.
        nhomo = -1
        do iband = 1, ihomo
          if (t%kpoints(1)%eigen(iband) .gt. energy_min) nhomo = nhomo + 1
        end do

        ntransitions = 0
        do iband = ihomo - nhomo + 1, t%norbitals_new
          if (t%kpoints(1)%eigen(iband) .lt. energy_min + energy_max) ntransitions = ntransitions + 1
        end do
        allocate (answer (-nhomo:0, ntransitions)); answer = 0.0d0

        ! Allocate transition type and initialize imap
        write (logfile,*) ' The number of transitions included in the '
        write (logfile,*) ' absorption integration is ntransitions = ', ntransitions
        write (logfile,*)
        allocate (t%kpoints(1)%transition (-nhomo:0, ntransitions))
        t%kpoints(1)%transition (-nhomo:0, 1:ntransitions)%imap = 0
        do ilevel = -nhomo, 0
          itransition = 1
          write (logfile,*) ' We will calculate dipole transitions between the two states: '
          write (logfile,*) ' Starting from energy state = ', ihomo + ilevel, ' energy = ', t%kpoints(1)%eigen(ihomo + ilevel)
          do iband = ihomo + ilevel + 1, t%norbitals_new
            if (t%kpoints(1)%eigen(iband) .lt. energy_min + energy_max) then
              t%kpoints(1)%transition(ilevel, itransition)%imap = 0
              write (logfile,200) itransition, iband, t%kpoints(1)%eigen(iband)
              ! initialize imap
              t%kpoints(1)%transition(ilevel, itransition)%imap = iband
              itransition = itransition + 1
            end if
          end do
          write (logfile,*)
        end do

! If the .abs.strengths file exists, then this is a possible restart.
! Inquire about the file and assume that it is a restart.
! Write to logfile that this is a restart so that the user knows to remove
! the abs.strengths file.
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.abs.strengths'
        inquire (file = slogfile, exist = read_abs)
        if (read_abs) then
          write (logfile,*)
          write (logfile,*) ' Reading from abs.strengths file! '
          open (unit = inpfile, file = slogfile, status = 'old')
          do ilevel = -nhomo, 0
            read (inpfile, 101, end = 10) ihomo_begin, homo_begin
            if (ilevel .eq. ihomo_begin) then
              write (logfile,*) ' Absorption strengths previously calculated for state: ', ilevel
              itransition_begin = 1
              do iband = ihomo + ilevel + 1, t%norbitals_new
                if (t%kpoints(1)%eigen(iband) .lt. energy_min + energy_max) then
                  ! initialize imap
                  t%kpoints(1)%transition(ilevel, itransition_begin)%imap = iband
                  read (inpfile, *, end = 10) t%kpoints(1)%eigen(iband), answer(ilevel, itransition_begin)
                  write (logfile,*) t%kpoints(1)%eigen(iband), answer(ilevel, itransition_begin)
                  itransition_begin = itransition_begin + 1
                  if (itransition_begin .gt. ntransitions) exit
                end if
              end do
            else
              backspace (inpfile)
            end if
          end do
10        continue
          close (unit = inpfile)
        end if

!----------------------------------------------------------------------------
! Non Adaptive Simpson's Setup
!----------------------------------------------------------------------------
! Integration is over Cartesian coordinates.
! Determine the maximum and minimum coordinate values.
! Initialize to some very large numbers
        xmin = 999.0d0; xmax = -999.0d0
        ymin = 999.0d0; ymax = -999.0d0
        zmin = 999.0d0; zmax = -999.0d0

! Loop over the atoms in the central cell.
! If we are doing a cluster, then we simply take the maximum extent
! of the atoms wavefunctions and integrate over this maximum extent.
! If we are doing a periodic cell, then we choose the minimum
! points as the coordinate of the first atom in the cell and then
! integrate along the lattice vectors.
        if (t%icluster .eq. 1) then
          do iatom = 1, t%natoms
            r1 = t%atom(iatom)%ratom
            in1 = t%atom(iatom)%imass
            xmin = min(xmin, r1(1) - wf(in1)%rcutoffA_max)
            xmax = max(xmax, r1(1) + wf(in1)%rcutoffA_max)
            ymin = min(ymin, r1(2) - wf(in1)%rcutoffA_max)
            ymax = max(ymax, r1(2) + wf(in1)%rcutoffA_max)
            zmin = min(zmin, r1(3) - wf(in1)%rcutoffA_max)
            zmax = max(zmax, r1(3) + wf(in1)%rcutoffA_max)
          end do  ! end loop over atoms

          ! Add small buffer
          xmin = xmin - 0.1d0; xmax = xmax + 0.1d0
          ymin = ymin - 0.1d0; ymax = ymax + 0.1d0
          zmin = zmin - 0.1d0; zmax = zmax + 0.1d0

          write (logfile, *)
          write (logfile, *) ' Absorption calculation, grid box size = '
          write (logfile, *) ' xmin, xmax = ', xmin, xmax
          write (logfile, *) ' ymin, ymax = ', ymin, ymax
          write (logfile, *) ' zmin, zmax = ', zmin, zmax
        else
          write (*,*) ' Cannot integrate over periodic cell yet '
          stop
        end if

! ***************************************************************************
! SET UP SIMPSON FACTORS
! ***************************************************************************
! Strictly define what the density of the mesh should be.  Make the density
! of the number of points equivalent for all cases. Change the number of
! points to be integrated to be dependent upon the distance between the
! centers and this defined density.
        nnx = int((xmax - xmin)/dgrid)
        if (mod(nnx,2) .eq. 0) nnx = nnx + 1; allocate (xmult (nnx))
        nny = int((ymax - ymin)/dgrid)
        if (mod(nny,2) .eq. 0) nny = nny + 1; allocate (ymult (nny))
        nnz = int((zmax - zmin)/dgrid)
        if (mod(nnz,2) .eq. 0) nnz = nnz + 1; allocate (zmult (nnz))

! Set up Simpson's rule factors.
        xmult(1) = dgrid/3.0d0; xmult(nnx) = dgrid/3.0d0
        do ix = 2, nnx - 1, 2
          xmult(ix) = 4.0d0*dgrid/3.0d0
        end do
        do ix = 3, nnx - 2, 2
          xmult(ix) = 2.0d0*dgrid/3.0d0
        end do

        ymult(1) = dgrid/3.0d0; ymult(nny) = dgrid/3.0d0
        do iy = 2, nny - 1, 2
          ymult(iy) = 4.0d0*dgrid/3.0d0
        end do
        do iy = 3, nny - 2, 2
          ymult(iy) = 2.0d0*dgrid/3.0d0
        end do

        zmult(1) = dgrid/3.0d0; zmult(nnz) = dgrid/3.0d0
        do iz = 2, nnz - 1, 2
          zmult(iz) = 4.0d0*dgrid/3.0d0
        end do
        do iz = 3, nnz - 2, 2
          zmult(iz) = 2.0d0*dgrid/3.0d0
        end do

! Write the absorption profile output
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.abs.strengths'
        open (unit = inpfile, file = slogfile, status = 'unknown')

! Loop over HOMO levels
        ilevel = 0
        do ilevel = -nhomo, 0
          write (inpfile, 101) ilevel, t%kpoints(1)%eigen(ihomo + ilevel)

! Loop over the transitions. Here we are doing a sum of dipole
! transitions from HOMO to other excited states, state iband.
! If the answer was read into a file, then write to file, else calculate it.
          write (logfile, *)
!$omp parallel private (itransition, iband, ir, l, m, iz, iy, ix, xtemp)     &
!$omp          private (zntegral, yntegral, xntegral, z, y, x, r_dipole)     &
!$omp          private (iatom, r1, in1, norb_mu, num_neigh, ineigh, mbeta)   &
!$omp          private (jatom, r2, in2, norb_nu)                             &
!$omp          private (inu, nnu, n2, r2max, l2, m2, r1p, r1pin)             &
!$omp          private (imu, mmu, n1, r1max, l1, m1, r2p, r2pin)
!$omp do
          write (logfile,*) ' Starting from energy state = ', ihomo + ilevel, ' energy = ', t%kpoints(1)%eigen(ihomo + ilevel)
          do itransition = 1, ntransitions
            write (logfile, *) ' Working on itransition = ', itransition
            iband = t%kpoints(1)%transition(ilevel, itransition)%imap
            if (iband .eq. 0) then
              write (logfile, *) ' The value iband = 0 is obtained for: '
              write (logfile, *) ' ilevel = ', ilevel, ' itransition = ', itransition, '.'
              write (logfile, *) ' This band is not an allowed transition; there '
              write (logfile, *) ' is no band number associated with this transition! '
              cycle
            end if
            if (t%kpoints(1)%eigen(iband) .gt. energy_min + energy_max) then
              write (logfile, *) ' The energy eigenvalue = ', t%kpoints(1)%eigen(iband)
              write (logfile, *) ' There is very likley something wrong here - this '
              write (logfile, *) ' transition should not be calculated. '
              stop
            end if
            if (answer(ilevel, itransition) .ne. 0.0d0) then
              write (logfile, *) ' Skipping calculation, this dipole strength already calculated. '
              write (inpfile, *) t%kpoints(1)%eigen(iband), answer(ilevel, itransition)
            else

! Loop over the three dipole components
              do ir = 1, 3
                l = 1  ! angular for dipole is p-like
                m = m_dipole(ir)

! ***************************************************************************
! I N T E G R A T I O N
! ***************************************************************************
! For a cluster we integrate from the minimum points in x, y, and z to the
! maximum points in x, y, and z.
                zntegral = 0.d0
                do iz = 1, nnz
                  z = zmin + dfloat(iz-1)*dgrid
                  r_dipole(3) = z

                  yntegral = 0.d0
                  do iy = 1, nny
                    y = ymin + dfloat(iy-1)*dgrid
                    r_dipole(2) = y

                    xntegral = 0.d0
                    do ix = 1, nnx
                      x = xmin + dfloat(ix-1)*dgrid
                      r_dipole(1) = x

! Loop over the atoms in the central cell.
                      do iatom = 1, t%natoms
                        r1 = t%atom(iatom)%ratom
                        in1 = t%atom(iatom)%imass
                        norb_mu = species(in1)%norb_max
                        num_neigh = t%neighbors(iatom)%neighn

                        ! Evaluate the value of r from the center of the atom
                        r1p = r_dipole - r1
                        r1pin = sqrt(r1p(1)**2 + r1p(2)**2 + r1p(3)**2)
                        if (r1pin .gt. wf(in1)%rcutoffA_max) cycle

! Loop over the neighbors of each iatom.
                        do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
                          mbeta = t%neighbors(iatom)%neigh_b(ineigh)
                          jatom = t%neighbors(iatom)%neigh_j(ineigh)
                          r2 = t%atom(jatom)%ratom + t%xl(mbeta)%a
                          in2 = t%atom(jatom)%imass
                          norb_nu = species(in2)%norb_max

                          ! Evaluate the value of r from the center of the atom
                          r2p = r_dipole - r2
                          r2pin = sqrt(r2p(1)**2 + r2p(2)**2 + r2p(3)**2)
                          if (r2pin .gt. wf(in2)%rcutoffA_max) cycle

! Loop over the orbital states of the atom and its neighbor
                          do inu = 1, norb_nu
                            nnu = inu + t%iblock_slot(jatom)
                            n2 = species(in2)%orbital(inu)%issh
                            r2max = species(in2)%shell(n2)%rcutoffA
                            l2 = species(in2)%orbital(inu)%l
                            m2 = species(in2)%orbital(inu)%m

                            do imu = 1, norb_mu
                              mmu = imu + t%iblock_slot(iatom)
                              n1 = species(in1)%orbital(imu)%issh
                              r1max = species(in1)%shell(n1)%rcutoffA
                              l1 = species(in1)%orbital(imu)%l
                              m1 = species(in1)%orbital(imu)%m

                              xtemp = t%kpoints(1)%c(mmu,ihomo + ilevel)     &
     &                               *psiofr(r1pin, r1max, in1, n1)*Ylm (r1p, l1, m1)  &
     &                               *r_dipole(ir)*t%kpoints(1)%c(nnu,iband) &
     &                               *psiofr(r2pin, r2max, in2, n2)*Ylm (r2p, l2, m2)

                              ! add to sum
                              xntegral = xntegral + xtemp * xmult(ix)
                            end do  ! end loop over the orbitals of iatom
                          end do ! end loop over the orbitals of jatom

                        end do ! end loop over neighbors
                      end do ! end loop over atoms
                    end do  ! end Simpson's for x-dimension

                    ! add to sum
                    yntegral = yntegral + xntegral * ymult(iy)
                  end do  ! end Simpson's for y-dimension

                  ! add to sum
                  zntegral = zntegral + yntegral * zmult(iz)
                end do  ! end Simpson's for z-dimension

! End of integration
! ***************************************************************************

                ! answer to the integral
                answer(ilevel, itransition) = answer(ilevel, itransition) + zntegral**2
              end do  ! end loop over the dipole components
            end if
            iband = t%kpoints(1)%transition(ilevel, itransition)%imap
            if (t%kpoints(1)%transition(ilevel, itransition)%imap .ne. 0) then
              write (inpfile, *) t%kpoints(1)%eigen(iband), answer(ilevel, itransition)
            end if
          end do  ! end loop over transitions
!$omp end do
!$omp end parallel
        end do  ! end loop over HOMO levels
        close (unit = inpfile)

! Evaluate the intensity and divide by the normalization
! Write the absorption profile output
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.abs'
        open (unit = inpfile, file = slogfile, status = 'unknown')

        do igrid = 1, nenergy_grid
          energy = energy_step*dfloat(igrid - 1)
          xntensity = 0.0d0
          do ilevel = -nhomo, 0
            do itransition = 1, ntransitions
              iband = t%kpoints(1)%transition(ilevel, itransition)%imap
              energy_min = t%kpoints(1)%eigen(ihomo + ilevel)
              if (iband .gt. ihomo) then
                xntensity = xntensity                                        &
     &           + energy*answer(ilevel, itransition)                        &
     &            /((energy - (t%kpoints(1)%eigen(iband) - energy_min))**2 + eta**2)
              end if
            end do  ! end loop over transitions
          end do  ! end loop over HOMO levels
          write (inpfile,*) energy, xntensity
        end do  ! end loop over energy grid
        close (unit = inpfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (2x, i3, f12.6)
102     format (2x, 3f10.4)
200     format (2x, ' Transition number ', i5, ' to energy state ', i5,      &
     &               ' energy = ', f7.3)

! End Subroutine
! ===========================================================================
        return
        end subroutine absorption
