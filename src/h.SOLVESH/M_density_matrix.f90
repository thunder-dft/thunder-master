! copyright info:
!
!                             @Copyright 2016
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
!       This routine calculates the density matrices and the band-structure
! energy, as well as similar density matrices.
!
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
! ============================================================================
! Module declaration
! ============================================================================
        module M_density_matrix
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_2c
        use M_kspace
        use M_rotations

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! density_matrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine calculates the density matix rho and stores it in the
! structure given in M_assemble_block.f90
!
! ===========================================================================
        subroutine density_matrix (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh              !< counter over atoms and neighbors
        integer iband, ikpoint             !< counter of band and kpoint
        integer ihomo                      !< highest occupied level
        integer imu, inu
        integer in1, in2                   !< species numbers
        integer issh
        integer jatom                      !< neighbor of iatom
        integer logfile                    !< writing to which unit
        integer mmu, nnu

        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor
        integer norb_mu, norb_nu         !< size of the block for the pair

        real dot                         !< dot product between K and r
        real efermi
        real gutr                        !< real part of density matrix
        real qstate

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        logical read_occupy

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
        write (logfile,*) ' Total number of electrons in the system is, ztot = ', s%ztot

! Calculate the Fermi energy.
! FIX ME! For now we set qstate = 0.0d0
        qstate = 0.0d0

! Inquire here regarding the occupations file
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.OCCUPATION'
        inquire (file = slogfile, exist = read_occupy)
        if (read_occupy) then
          call read_fermie (s)
        else
          call fermie (s, qstate, efermi)
        end if
        write (logfile,*) ' Fermi Level = ', efermi

! write out eigenvalues to eigen.dat
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
          open (unit = 22, file = slogfile, status = 'replace', form = 'unformatted')

          do ikpoint = 1, s%nkpoints
            do iband = 1, s%norbitals_new
              write (22) iband

! write out the coefficient
              write (22) (s%kpoints(ikpoint)%c(inu,iband), inu = 1, s%norbitals_new)
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
! cape_matrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine calculates the cape matix rho and stores it in the
! structure given in M_assemble_block.f90
!
!       The cape matrix is a*(iband,imu)*a(iband,inu)*eigen(iband) (see Eq. 47
! in Sankey_Niklewski paper). The coeficient is called a(iband,inu) but in this
! routine it is called c(nnu,iband)
! ===========================================================================
        subroutine cape_matrix (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh              !< counter over atoms and neighbors
        integer iband, ikpoint             !< counter of band and kpoint
        integer imu, inu
        integer in1, in2                   !< species numbers
        integer jatom                      !< neighbor of iatom
        integer logfile                    !< writing to which unit
        integer mmu, nnu

        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor
        integer norb_mu, norb_nu         !< size of the block for the pair

        real dot                         !< dot product between K and r
        real gutr                        !< real part of density matrix

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pCape_neighbors
        type(T_assemble_neighbors), pointer :: pcapemat
        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat

! Allocate Arrays
! ===========================================================================
        allocate (s%capemat (s%natoms))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! ****************************************************************************
!
!                 C O M P U T E   C A P E   D E N S I T I E S
! ****************************************************************************
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pcapemat=>s%capemat(iatom)
          pdenmat=>s%denmat(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pcapemat%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pCape_neighbors=>pcapemat%neighbors(ineigh)
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (pCape_neighbors%block(norb_mu, norb_nu))
            pCape_neighbors%block = 0.0d0

! Loop over the special k points.
            do ikpoint = 1, s%nkpoints
              ! Find the phase which is based on k*r
              vec = r2 - r1
              sks = s%kpoints(ikpoint)%k
              dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
              phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight*P_spin


! Loop over all bands
              do iband = 1, s%norbitals_new
                if (s%kpoints(ikpoint)%ioccupy(iband) .ne. 0) then
                  phase = phasex*s%kpoints(ikpoint)%foccupy(iband)
                  do imu = 1, norb_mu
                    mmu = imu + s%iblock_slot(iatom)
                    step1 = phase*conjg(s%kpoints(ikpoint)%c(mmu,iband))!

                    do inu = 1, norb_nu
                      nnu = inu + s%iblock_slot(jatom)
                      step2 = step1*s%kpoints(ikpoint)%c(nnu,iband)
                      gutr = real(step2)

! Finally the density matrix:
                      pCape_neighbors%block(imu,inu) =                         &
     &                  pCape_neighbors%block(imu,inu) + s%kpoints(ikpoint)%eigen(iband)*gutr

                    end do
                  end do
                end if

! Finish loop over bands.
              end do

! Finish loop over k-points.
            end do

          end do ! Finish loop over atoms and neighbors.
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
        end subroutine cape_matrix


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
          write (*,*) ' In denmat.f - ztest .ne. ztot! '
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
! read_fermie.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>    Allows for occupations beyond just the ground state - read in from file.
!
! ===========================================================================
        subroutine read_fermie (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iband
        integer ikpoint
        integer imu
        integer inpfile                     !< reading from which unit
        integer logfile                     !< writing to which unit
        integer noccupy

        real fband

        character (len = 25) :: slogfile

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! Open the file and read information.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.OCCUPATION'
        open (unit = inpfile, file = slogfile, status = 'old')

        write (logfile,*)
        write (logfile,*) ' Reading from the OCCUPATION file! '
        do ikpoint = 1, s%nkpoints
          s%kpoints(ikpoint)%ioccupy = 0
          s%kpoints(ikpoint)%foccupy = 0.0d0
          read (inpfile,*) noccupy
          if (noccupy .gt. s%norbitals) then
            write (*,*) ' noccupy > norbitals: from OCCUPATION file. '
            stop
          end if
          do imu = 1, noccupy
            read (inpfile,*) iband, fband
            s%kpoints(ikpoint)%ioccupy(imu) = iband
            s%kpoints(ikpoint)%foccupy(imu) = fband
          end do
        end do
        close (unit = inpfile)

! End Subroutine
! ===========================================================================
        return
        end subroutine read_fermie


! ===========================================================================
! fermie.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the fermie energy and fermi occupations.
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
! Program Declaration
! ===========================================================================
        subroutine fermie (s, qstate, efermi)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        real, intent (in) :: qstate      !< FIXME

! Output
        real, intent (out) :: efermi

! Local Parameters
! ===========================================================================
        integer, parameter :: imax = 1000 !< maximum sc iterations
        integer, parameter :: nmax = 5000 !< cutoff for degeneracy check

        real, parameter :: tol = 1.0d-12

! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer imu
        integer inu
        integer iteration
        integer jkpoint

        real delta
        real emin
        real emax
        real qcharge
        real qztot
        real temperature

        real, pointer :: peigen_mu_k
        real, pointer :: peigen_nu_k

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Add in the qstate to the total charge
        qztot = s%ztot + qstate

! The subroutine fermie needs a temperature to calculate the occupations of
! the states so set temperature to some low value (1 eV = 11604 K).
        temperature = efermi_T/P_kconvert

! Find emin and emax. Also make sure degenerate eigenvalues are truly
! degenerate. However, if norbitals*nkpts*norbitals*nkpts is larger than nmax,
! then skip the degeneracy checking. Otherwise, the checking can take a while.
        if (s%norbitals**2*s%nkpoints**2 .lt. nmax) then
          emin = s%kpoints(1)%eigen(1)
          emax = s%kpoints(1)%eigen(s%norbitals)
          do ikpoint = 1, s%nkpoints
            do imu = 1, s%norbitals
              peigen_mu_k=>s%kpoints(ikpoint)%eigen(imu)
              if (peigen_mu_k .lt. emin) emin = peigen_mu_k
              if (peigen_mu_k .gt. emax) emax = peigen_mu_k
              do jkpoint = ikpoint, s%nkpoints
                do inu = imu, s%norbitals
                  peigen_nu_k=>s%kpoints(ikpoint)%eigen(inu)
                  if (abs(peigen_mu_k - peigen_nu_k) .lt. tol) then
                    s%kpoints(ikpoint)%eigen(inu) =                            &
     &                (peigen_mu_k + peigen_nu_k)/2.0d0
                    s%kpoints(ikpoint)%eigen(imu) = peigen_nu_k
                  end if
                end do
              end do
            end do
          end do
        else
          open (11, file = 'WARNINGS', status = 'unknown', position = 'append')
          write (11,*) '  '
          write (11,*) ' ************ WARNING ******* WARNING *********** '
          write (11,*) '          skipping the degeneracy checking  '
          write (11,*) '               in subroutine fermie'
          write (11,*) ' ************************************************ '
          close (11)
          emin = s%kpoints(1)%eigen(1)
          emax = s%kpoints(1)%eigen(s%norbitals)
          do ikpoint = 1, s%nkpoints
            do imu = 1, s%norbitals
              peigen_mu_k=>s%kpoints(ikpoint)%eigen(imu)
              if (peigen_mu_k .lt. emin) emin = peigen_mu_k
              if (peigen_mu_k .gt. emax) emax = peigen_mu_k
            end do
          end do
        end if

! The value of efermi must be between emin and emax
        iteration = 0
        qcharge = 0.0d0
        do while (abs(qcharge - qztot) .gt. tol .and. iteration .le. imax)
          iteration = iteration + 1

! Make a guess at efermi
          efermi = (emax + emin)/2.0d0
          qcharge = 0.0d0
          do ikpoint = 1, s%nkpoints
            do imu = 1, s%norbitals
              delta = (s%kpoints(ikpoint)%eigen(imu) - efermi)/temperature

! Skip exponential for big -/+ delta
              if (delta .gt. 10.0d0) then
                s%kpoints(ikpoint)%foccupy(imu) = 0.0d0
                s%kpoints(ikpoint)%ioccupy(imu) = 0
              else if (delta .lt. -10.0d0) then
                s%kpoints(ikpoint)%foccupy(imu) = 1.0d0
                s%kpoints(ikpoint)%ioccupy(imu) = 1
              else
                s%kpoints(ikpoint)%foccupy(imu) = 1.0d0/(1.0d0 + exp(delta))
                if (s%kpoints(ikpoint)%foccupy(imu) .gt. 1.0d-5) then
                  s%kpoints(ikpoint)%ioccupy(imu) = 1
                else
                  s%kpoints(ikpoint)%ioccupy(imu) = 0
                end if
              end if
              qcharge = qcharge                                              &
     &          + P_spin*s%kpoints(ikpoint)%foccupy(imu)                     &
     &                  *s%kpoints(ikpoint)%weight
            end do
          end do
! Narrow the range that efermi can fall into
          if (qcharge .gt. qztot) then
            emax = efermi
          else
            emin = efermi
          end if
        end do

! Print warning for going over maximum iterations.
        if (iteration .gt. imax) then
          open (11, file = 'WARNINGS', status = 'unknown', position = 'append')
          write (11,*)
          write (11,*) ' ************ WARNING ******* WARNING *********** '
          write (11,*) '        not under tolerance (toll) after ',imax
          write (11,*) '          iterations in subroutine fermie'
          write (11,*)
          write (11,*) '          qcharge = ', qcharge
          write (11,*) '          qztot = ', qztot
          write (11,*) '          emax = ', emax
          write (11,*) '          emin = ', emin
          write (11,*) ' ************************************************ '
          write (11,*)
          close (11)
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
        end subroutine fermie


! ===========================================================================
! initialize_mdet.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       Dummy routine - need to have dummy initialize_mdet
!
! ===========================================================================
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

        icurrent_state=0

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_mdet


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

        character (len = 25) :: slogfile

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

! Open the file and write information.
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.denmat'
        open (unit = 22, file = slogfile, status = 'unknown')
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          r1 = s%atom(iatom)%ratom
          write (22,*)
          write (22,101)
          num_neigh = s%neighbors(iatom)%neighn
          write (22,*) ' There are ', num_neigh, ' neighbors to atom ', iatom
          write (22,101)

          ! cut some lengthy notation
          pdenmat=>s%denmat(iatom)
          do ineigh = 1, num_neigh
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            in2 = s%atom(jatom)%imass
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a

            sigma = r2 - r1
            z = distance (r1, r2)
            write (22,*) '  '
            write (22,102) ineigh, jatom, mbeta, z

! Density matrix elements
            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            write (22,103)
            norb_mu = species(in1)%norb_max
            norb_nu = species(in2)%norb_max
            do imu = 1, norb_mu
              write (22,104) (pRho_neighbors%block(imu,inu), inu = 1, norb_nu)
            end do
          end do
        end do
        write (22,101)

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
