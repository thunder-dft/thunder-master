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

! M_density_matrix.f90
! Program Description
! ===========================================================================
!>       This routine calculates the density matrices and the band-structure
!! energy, as well as similar density matrices.
!
! ===========================================================================
! Original code written by
!> @author Otto F. Sankey
! with modification by
!> @author Alex A. Demkov
!! @author Jose Ortega
! ===========================================================================
! Code written by:
!> @author James P. Lewis with contributions from Amanda J Neukirch
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ============================================================================
! Module declaration
! ============================================================================
        module M_density_matrix
        use M_assemble_blocks
        use M_configuraciones
        use M_kspace

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! initialize_mdet.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine gives the initial state for molecular dynamics with
!       electronic transitions (MDET) (nonadiabatic calculation)
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initialize_mdet (s, icurrent_state)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Output
        integer, intent (out) :: icurrent_state

! Local Parameters
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                      !< counter over atoms
        integer iband, ikpoint             !< counter of band and kpoint
        integer ioccupy                    !< input occupation number
        integer in1                        !< species number
        integer inpfile                     !< reading from which unit
        integer issh
        integer itransition                !< number of transitions
        integer logfile                    !< writing to which unit
        integer nfermi                     !< Fermi level state
        integer ipop

        real qztot                          !< total number of electrons
        real foccupy

        character (len = 25) :: slogfile

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! Initialize occupations
        do ikpoint = 1, s%nkpoints
          do ioccupy = 1, s%norbitals
            s%kpoints(ikpoint)%ioccupy(ioccupy) = 0
            s%kpoints(ikpoint)%foccupy(ioccupy) = 0.0d0
          end do
        end do

! Loop over the atoms.
! Total charge - ztot
        s%ztot = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            s%ztot = s%ztot + species(in1)%shell(issh)%Qneutral
          end do
        end do

        qztot = s%ztot
        nfermi = int(qztot)/2

        do ikpoint = 1, s%nkpoints
          do iband = 1, nfermi
            s%kpoints(ikpoint)%foccupy(iband) = 1.0d0
            s%kpoints(ikpoint)%ioccupy(iband) = 1
          end do
        end do

! Read information from .mdet.input file
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.mdet.inp'
        open (unit = inpfile, file = slogfile, status = 'old')
        write (logfile,*)
        write (logfile,*) ' Reading from mdet.input file! '
        read (inpfile,*) ntransitions
        do ikpoint = 1, s%nkpoints

! Allocate transition type and initialize imap
          allocate (s%kpoints(ikpoint)%transition (ntransitions))

          do itransition = 1, ntransitions
            read (inpfile,*) iband, foccupy, ipop
! FIXME! We might need to read in foccupy and set ioccupy to 1 when foccupy
! is not 0 (since including all states no fixed relationship between ioccupy
! and foccupy)

            ! initialize imap
            s%kpoints(ikpoint)%transition(itransition)%imap = iband
            ! initialize population
            s%kpoints(ikpoint)%transition(itransition)%cna = ipop    !might be errors with type diff. between ipop and cna
            if (ipop .eq. 1) then
              icurrent_state = iband
            end if
            s%kpoints(ikpoint)%foccupy(iband) = foccupy
            if (foccupy .ge. 0.5) then
			   s%kpoints(ikpoint)%ioccupy(iband) = 1
            end if
          end do   ! end loop over transitons
          write (logfile,*) ' testing imaps reach '
          write (logfile,*) s%kpoints(ikpoint)%transition(1)%imap
        end do   ! end loop over kpoints
        close (unit = inpfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_mdet


! ===========================================================================
! denmat.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine calculates the density matix rho and stores it in the
!! structure given in M_assemble_block.f90
!
! ===========================================================================
        subroutine density_matrix (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh            !< counter over atoms and neighbors
        integer iband, ikpoint           !< counter of band and kpoint
        integer ihomo                    !< highest occupied level
        integer imu, inu
        integer in1, in2                 !< species numbers
        integer issh
        integer jatom                    !< neighbor of iatom
        integer logfile                    !< writing to which unit
        integer mmu, nnu

        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor
        integer norb_mu, norb_nu         !< size of the block for the pair

        real dot                           !< dot product between K and r
        real efermi
        real gutr                          !< real part of density matrix

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        character (len = 25) :: slogfile

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat

! Allocate Arrays
! ===========================================================================
        allocate (s%denmat (s%natoms))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Loop over the atoms.
! Total charge - ztot
        s%ztot = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            s%ztot = s%ztot + species(in1)%shell(issh)%Qneutral
          end do
        end do

        write (logfile,*)
        write (logfile,*) ' Calculating density matrix elements here. '
        write (logfile,*) ' Total number of electrons in the system is ztot = ', s%ztot

! write out eigenvalues to eigen.dat
        efermi = 0.0d0
        if (modulo(int(s%ztot), 2) .eq. 0) then
          ihomo = int(s%ztot)/2
        else
          ihomo = int(s%ztot)/2 + 1
        end if
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.eigen'
        open (unit = 22, file = slogfile, status = 'replace')
        do ikpoint = 1, s%nkpoints
          write (22,100) ikpoint, efermi,                                    &
     &      s%kpoints(ikpoint)%eigen(ihomo+1) - s%kpoints(ikpoint)%eigen(ihomo)

! write out the eigenvalues
          write (22,200) (s%kpoints(ikpoint)%eigen(imu), imu = 1, s%norbitals)
        end do
        close (unit = 22)

! ****************************************************************************
!
!                      C O M P U T E    D E N S I T I E S
! ****************************************************************************
! Write out the coefficients to a file - *.cdcoeffs
        if (iwriteout_cdcoeffs .eq. 1) then
          slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
          slogfile = trim(slogfile)//'.cdcoeffs'
!         open (unit = 22, file = slogfile, status = 'replace', form = 'unformatted')
          open (unit = 22, file = slogfile, status = 'replace')

          do ikpoint = 1, s%nkpoints
            do iband = 1, s%norbitals_new
              write (22,*) iband

! write out the coefficient
              write (22,*) (s%kpoints(ikpoint)%c(inu,iband), inu = 1, s%norbitals_new)
            end do
          end do
          close (unit = 22)
        end if

! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pdenmat=>s%denmat(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pdenmat%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (pRho_neighbors%block(norb_mu, norb_nu))
            pRho_neighbors%block = 0.0d0

! Loop over the special k points.
            do ikpoint = 1, s%nkpoints

              ! Find the phase which is based on k*r
              vec = r2 - r1
              sks = s%kpoints(ikpoint)%k
              dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
              phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight*P_spin

! Loop over all bands
! Here we assume that the coefficients are real only at this point.
              do iband = 1, s%norbitals_new
                if (s%kpoints(ikpoint)%ioccupy(iband) .ne. 0) then
                  phase = phasex*s%kpoints(ikpoint)%foccupy(iband)
                  do imu = 1, norb_mu
                    mmu = imu + s%iblock_slot(iatom)
                    step1 = phase*conjg(s%kpoints(ikpoint)%c(mmu,iband))

                    do inu = 1, norb_nu
                      nnu = inu + s%iblock_slot(jatom)
                      step2 = step1*s%kpoints(ikpoint)%c(nnu,iband)
                      gutr = real(step2)

! Finally the density matrix:
                      pRho_neighbors%block(imu,inu) =                        &
     &                  pRho_neighbors%block(imu,inu) + gutr
                    end do
                  end do
                end if

! Finish loop over bands.
              end do

! Finish loop over k-points.
            end do

! Finish loop over atoms and neighbors.
          end do
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i4, 2(2x, f10.6))
200     format (4(2x, f12.4))

! End Subroutine
! ===========================================================================
        return
        end subroutine density_matrix


! ===========================================================================
! calculate_ebs.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine calculates the band-structure energy.
!
! ===========================================================================
        subroutine calculate_ebs (s, ebs)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used
! Output
        real, intent (out) :: ebs      !< band-structure energy

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iband, ikpoint           !< counter of band and kpoint

        real ztest

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ****************************************************************************
!
!  C O M P U T E    B A N D - S T R U C T U R E    E N E R G Y
! ****************************************************************************
! Compute ebs, the band structure energy.
        ebs = 0.0d0
        ztest = 0.0d0
        do ikpoint = 1, s%nkpoints
          do iband = 1, s%norbitals_new
            if (s%kpoints(ikpoint)%ioccupy(iband) .eq. 1) then
              ebs = ebs + s%kpoints(ikpoint)%weight*P_spin*s%kpoints(ikpoint)%eigen(iband)  &
     &                                             *s%kpoints(ikpoint)%foccupy(iband)
              ztest = ztest + s%kpoints(ikpoint)%weight*P_spin*s%kpoints(ikpoint)%foccupy(iband)
             end if
          end do
        end do

! Test to make sure we get the proper number of states.
        if (abs(ztest - s%ztot) .gt. 1.0d-02) then
          write (*,*) ' *************** error *************** '
          write (*,*) ' ztest = ', ztest, ' ztot = ', s%ztot
          write (*,*) ' In calculate_ebs - ztest .ne. ztot! '
          stop
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine calculate_ebs


! ===========================================================================
! writeout_density.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine is a utility to write out the density matrix.
!
! ===========================================================================
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine writeout_density (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh          !< counters for atom, neighbor loops
        integer in1, in2               !< species numbers
        integer imu, inu               !< counters for mu, nu
        integer jatom                  !<
        integer logfile                !< writing to which unit
        integer mbeta                  !< the cell containing neighbor of iatom
        integer num_neigh              !< number of neighbors
        integer norb_mu, norb_nu       !< block size for the H and S block

        real z                         !< distance between two atoms

        real, dimension (3) :: r1, r2  !< positions for iatom, jatom
        real, dimension (3) :: sigma   !< direction along sigma bond

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*) ' '
        write (logfile,*) ' In writeout_density.f90 '
        write (logfile,*) ' Writing out pieces of the density matrix. '
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          write (logfile,*)
          write (logfile,101)
          num_neigh = s%neighbors(iatom)%neighn
          write (logfile,*) ' There are ', num_neigh,                        &
     &                      ' neighbors to atom ', iatom
          write (logfile,101)

          ! cut some lengthy notation
          pdenmat=>s%denmat(iatom)
          do ineigh = 1, num_neigh
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            in2 = s%atom(jatom)%imass
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a

            sigma = r2 - r1
            z = distance (r1, r2)
            write (logfile,*) '  '
            write (logfile,102) ineigh, jatom, mbeta, z

! Density matrix elements
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            write (logfile,103)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (logfile,104)                                            &
     &         (pRho_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do
          end do
        end do
        write (logfile,101)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (75('*'))
102     format (1x, ' Matrices connected to neighbor ', i3, ','              &
     &          ' jatom = ', i4, ', ', ' mbeta = ', i4, ', ', ' d = ', f6.3)
103     format (75('='))
104     format (9f8.3)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_density


! ===========================================================================
! destroy_denmat
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing denmat
!! information.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
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
        subroutine destroy_denmat (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                             !< counter over atoms
        integer ineigh                            !< counter over neighbors

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%denmat(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%denmat(iatom)%neighbors)
        end do
        deallocate (s%denmat)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_denmat

! End Module
! ===========================================================================
        end module M_density_matrix
