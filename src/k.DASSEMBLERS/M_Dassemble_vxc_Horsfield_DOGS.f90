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

! M_Dassemble_vxc
! Module Description
! ===========================================================================
!>       This is a module containing all of the  programs required
!! to assemble all of the matrix elements for the Horsfield
!! exchange-correlation. See Horsfield AP, 1997,
!! Efficient ab initio tight binding, Phys. Rev. B, Vol: 56, P
!! Pages: 6594-6602, ISSN: 1098-0121
!!
!! It contains the following subroutines within the module:
!!
!!       Dassemble_vxc_2c : for calculating two-center interactions
!!       Dassemble_vxc_3c : for calculating three-center interactions
!!
! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
! located in the Fdata directory.  This list will change depending on
! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_Dassemble_vxc

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones
        use M_neighbors
        use M_rotations
        use M_Drotations

! /FDATA
        use M_Fdata_2c
        use M_Fdata_3c

! /SOLVESH
        use M_density_matrix

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains


! ===========================================================================
! Dassemble_vxc.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>  This is the main module for assembling the Vxc matrix element
!! interactions (McWEDA) required for the forces.
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_vxc (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                      !< counter over atoms
        integer logfile                    !< writing to which unit
        integer num_neigh                  !< number of neighbors

        type(T_forces), pointer :: pfi

! Allocate Arrays
! ===========================================================================
! We build the vxc_ontop forces here, so we allocate and initialize
! The vxc_atom forces are built in M_build_forces because they are assembled
! differently than the vxc_ontop forces.
        do iatom = 1, s%natoms
          pfi=>s%forces(iatom)
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pfi%vxc_on_site (3, num_neigh)); pfi%vxc_on_site = 0.0d0
          allocate (pfi%vxc_off_site (3, num_neigh)); pfi%vxc_off_site = 0.0d0
        end do

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        write (logfile,*) ' Calculating Horsfield forces here. '

! Calculate rho_in (density) matrix elements
        write (logfile,*) ' Calling two-center vxc Dassemblers. '
        call Dassemble_vxc_2c (s)

        write (logfile,*) ' Calling three-center vxc Dassemblers. '
        call Dassemble_vxc_3c (s)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_vxc


! ===========================================================================
! Dassemble_vxc.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates derivatives of the matrix elements for the
! Hartree interactions. This is the self-consistent version of the code.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_vxc_2c (s)
        implicit none

        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer ispecies, in1, in2, in3 !< species numbers
        integer imu, inu                !< counter over orbitals
        integer jatom                   !< neighbor of iatom
        integer interaction             !< which interaction and subtype
        integer num_neigh               !< number of neighbors
        integer matom                   !< matom is the self-interaction atom
        integer mbeta                   !< the cell containing neighbor of iatom

        integer issh
        integer norb_mu, norb_nu        !< size of the block for the pair

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        real z                          !< distance between r1 and r2

! The charge transfer bit = need a +-dq on each orbital.
        real, allocatable :: dqorb (:)
        real, allocatable :: dQ (:)       !< charge on atom, i.e. ionic
        real, allocatable :: Q0 (:)       !< total neutral atom charge, i.e. Ztot
        real, allocatable :: Q (:)        !< total charge on atom

        ! quadratic fitting factors
        real, dimension (0:4) :: dQ_factor

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

! bcxcm = Hartree matrix in molecular coordinates
! bcxcx = Hartree matrix in crystal coordinates
! dbcxcm = derivative of Hartree matrix in molecular coordinates
! vdbcxcm = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcxcx = vectorized derivative of Hartree matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx
        real, dimension (:, :), allocatable :: dbcxcm
        real, dimension (:, :, :), allocatable :: vdbcxcm
        real, dimension (:, :, :), allocatable :: vdbcxcx

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_block), pointer :: pRho_neighbors_matom

        type(T_forces), pointer :: pfi

! Allocate Arrays
! ===========================================================================
        ! needed for charge transfer bits
        allocate (dqorb (nspecies))
        allocate (Q0 (s%natoms))
        allocate (Q (s%natoms))
        allocate (dQ (s%natoms))

! Procedure
! ===========================================================================
! Initialize the charge transfer bit
        do ispecies = 1, nspecies
          dqorb(ispecies) = 0.5d0
          if (species(ispecies)%nssh .eq. 1) dqorb(ispecies) = 0.25d0
        end do

! Calculate nuclear charge.
! ****************************************************************************
! We do a quadratic expansion:
!       f(q) = f(0)  +  f'(0)*q   +  1/2 f''(0)*q*q
!
! The derivatives are computed as:
!       f'(0)  = [ f(dq) - f(-dq) ] / 2dq
!       f''(0) = [ f(dq) - 2f(0) + f(-dq) ] / dq*dq
!
! We introduce linear factors:
!              linfac(0)   =  1.0
!              linfac(i)   =  -/+ * (1/2) *  q/dq        i = 1, 2
!
!    and quadratic factors:
!              quadfac(0)  = -2.0 * (1/2) * (q/dq)**2
!              quadfac(i)  =  1.0 * (1/2) * (q/dq)**2    i=1,2
!
!       f(0) = f(dq=0) ; f(1) = f(-dq) ; f(2) = f(dq)
!
! With this, f(q) is given as:
!       f(q) = sum_i  (linfac(i) + quadfac(i))*f(i)
!
! ****************************************************************************
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          Q0(iatom) = 0.0d0
          Q(iatom) = 0.0d0
          do issh = 1, species(in1)%nssh
            Q0(iatom) = Q0(iatom) + species(in1)%shell(issh)%Qneutral
            Q(iatom) = Q(iatom) + s%atom(iatom)%shell(issh)%Qin
          end do
          dQ(iatom) = Q(iatom) - Q0(iatom)
        end do

! Here we assemble only the atom cases first, then assemble the ontop cases.
! This is so that we can get the correct allocation size for the different
! blocks.  We calculate the atom cases in a separate loop.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          pdenmat=>s%denmat(iatom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)
            call Depsilon_2c (r1, r2, z, eps, deps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! CALL GetDMES AND GET VXC FOR ONTOP CASE
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.

            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in vbcxcx, the vectorized matrix
! elements; x means crytal coordinates.

! FORCES - ONTOP LEFT CASE
! ****************************************************************************
! For the vxc_ontopL case, the potential is in the first atom - left (iatom):
! dbcxcm is the "scalar" derivative of the matrix; vdbcxcm is the "vector"
! derivative of the matrix in molecular coordinates.  When we are done, we get:
! vdtx as the vector derivative of the matrix in crystal coordinates.
              interaction = P_vxc_ontop
              in3 = in2
              ideriv = 0

! bcxcm = Hartree matrix in molecular coordinates
! bcxcx = Hartree matrix in crystal coordinates
! dbcxcm = derivative of Hartree matrix in molecular coordinates
! vdbcxcm = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcxcx = vectorized derivative of Hartree matrix in crystal coordinates
              allocate (bcxcm (norb_mu, norb_nu)); bcxcm = 0.0d0
              allocate (dbcxcm (norb_mu, norb_nu)); dbcxcm = 0.0d0
              allocate (vdbcxcm (3, norb_mu, norb_nu)); vdbcxcm = 0.0d0
              allocate (vdbcxcx (3, norb_mu, norb_nu)); vdbcxcx = 0.0d0

              call getDMEs_Fdata_2c (in1, in2, interaction, ideriv, z,        &
     &                               norb_mu, norb_nu, bcxcm, dbcxcm)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                end do
              end do
              call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcxcm,     &
     &                      vdbcxcm, vdbcxcx)

! Notice the explicit negative sign, this makes it force like.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfi%vxc_off_site(:,ineigh) = pfi%vxc_off_site(:,ineigh)     &
      &             - pRho_neighbors%block(imu,inu)*vdbcxcx(:,imu,inu)
                end do
              end do

! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
              ideriv_min = 1
              ideriv_max = 4

              dQ_factor(1) = -dQ(iatom)/(2.0d0*dqorb(in1))
              dQ_factor(2) =  dQ(iatom)/(2.0d0*dqorb(in1))
              dQ_factor(3) = -dQ(jatom)/(2.0d0*dqorb(in2))
              dQ_factor(4) =  dQ(jatom)/(2.0d0*dqorb(in2))
              do ideriv = ideriv_min, ideriv_max
                call getDMEs_Fdata_2c (in1, in2, interaction, ideriv, z,      &
     &                                 norb_mu, norb_nu, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do
                call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcxcm,   &
     &                        vdbcxcm, vdbcxcx)

                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfi%vxc_off_site(:,ineigh) = pfi%vxc_off_site(:,ineigh)   &
      &               - pRho_neighbors%block(imu,inu)*dQ_factor(ideriv)*vdbcxcx(:,imu,inu)
                  end do
                end do
              end do
              deallocate (bcxcm, dbcxcm, vdbcxcm, vdbcxcx)
            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms

! FORCES - ATM CASE
! ****************************************************************************
! For the vxc_atom case, the potential is in the first atom - left (iatom):
! dbcxcm is the "scalar" derivative of the matrix; vdbcxcm is the "vector"
! derivative of the matrix in molecular coordinates.  When we are done, we get:
! vdtx as the vector derivative of the matrix in crystal coordinates.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          matom = s%neigh_self(iatom)

! cut some lengthy notation
          pdenmat=>s%denmat(iatom)
          pRho_neighbors_matom=>pdenmat%neighbors(matom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)
            ! unit vector in sigma direction.
            if (z .lt. 1.0d-05) then
              sighat(1) = 0.0d0
              sighat(2) = 0.0d0
              sighat(3) = 1.0d0
            else
              sighat = (r2 - r1)/z
            end if
            call epsilon_function (r2, sighat, eps)
            call Depsilon_2c (r1, r2, z, eps, deps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! For these interactions, there are no subtypes and ideriv = 0
            interaction = P_vxc_atom
            in3 = in1
            ideriv = 0

! Allocate block size
            norb_nu = species(in3)%norb_max
! bcxcm = Hartree matrix in molecular coordinates
! dbcxcm = derivative of Hartree matrix in molecular coordinates
! vdbcxcm = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcxcx = vectorized derivative of Hartree matrix in crystal coordinates
            allocate (bcxcm (norb_mu, norb_nu)); bcxcm = 0.0d0
            allocate (bcxcx (norb_mu, norb_nu)); bcxcx = 0.0d0
            allocate (dbcxcm (norb_mu, norb_nu)); dbcxcm = 0.0d0
            allocate (vdbcxcm (3, norb_mu, norb_nu)); vdbcxcm = 0.0d0
            allocate (vdbcxcx (3, norb_mu, norb_nu)); vdbcxcx = 0.0d0

            call getDMEs_Fdata_2c (in1, in2, interaction, ideriv, z,          &
     &                             norb_mu, norb_nu, bcxcm, dbcxcm)

! dbcxcm is the "scalar" derivative of the matrix; vdbcxcm is the "vector" of
! the matrix; vdbcxcx as the vector derivative of the matrix in crystal.

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
            do inu = 1, norb_nu
              do imu = 1, norb_mu
                if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
              end do
            end do

            call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcxcm,       &
     &                    vdbcxcm, vdbcxcx)
! Notice the explicit negative sign, this makes it force like.
            do inu = 1, norb_mu
              do imu = 1, norb_mu
                pfi%vxc_on_site(:,ineigh) = pfi%vxc_on_site(:,ineigh)         &
      &           - pRho_neighbors_matom%block(imu,inu)*vdbcxcx(:,imu,inu)
               end do
            end do

! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
            ideriv_min = 1
            ideriv_max = 4

            dQ_factor(1) = -dQ(iatom)/(2.0d0*dqorb(in1))
            dQ_factor(2) =  dQ(iatom)/(2.0d0*dqorb(in1))
            dQ_factor(3) = -dQ(jatom)/(2.0d0*dqorb(in2))
            dQ_factor(4) =  dQ(jatom)/(2.0d0*dqorb(in2))
            do ideriv = ideriv_min, ideriv_max
              call getDMEs_Fdata_2c (in1, in2, interaction, ideriv, z,        &
     &                               norb_mu, norb_nu, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  if (z .gt. 1.0d-3) vdbcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                end do
              end do
              call Drotate (in1, in3, eps, deps, norb_mu, norb_nu, bcxcm,     &
     &                      vdbcxcm, vdbcxcx)

                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfi%vxc_on_site(:,ineigh) = pfi%vxc_on_site(:,ineigh)         &
      &               - pRho_neighbors_matom%block(imu,inu)*dQ_factor(ideriv)*vdbcxcx(:,imu,inu)
                  end do
                end do
              end do

            deallocate (bcxcm, bcxcx, dbcxcm, vdbcxcm, vdbcxcx)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
        deallocate (dqorb, Q0, Q, dQ)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_vxc_2c


! ===========================================================================
! Dassemble_vxc_3c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates derivative of the neutral atom potential
! matrix interactions.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
!> @author Barry Haycock
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
        subroutine Dassemble_vxc_3c (s)
        implicit none

        include '../include/interactions_3c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha, iatom, jatom        !< the three parties involved
        integer ibeta, jbeta                !< cells for three atoms
        integer ineigh, mneigh              !< counter over neighbors
        integer ispecies, in1, in2, indna   !< species numbers
        integer interaction                 !< which interaction and subtype
        integer imu, inu, iindex            !< indexing counters

        integer issh
        integer norb_mu, norb_nu            !< size of the block for the pair

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        real z                              !< distance between r1 and r2
        real x, cost                        !< dnabc and angle

! The charge transfer bit = need a +-dq on each orbital.
        real, allocatable :: dqorb (:)
        real, allocatable :: dQ (:)        !< charge on atom, i.e. ionic
        real, allocatable :: Q0 (:)        !< total neutral atom charge, i.e. Ztot
        real, allocatable :: Q (:)         !< total charge on atom

        ! quadratic fitting factors
        real, dimension (0:6) :: dQ_factor

        real, dimension (3, 3) :: eps       !< the epsilon matrix
        real, dimension (3) :: r1, r2, rna  !< positions - iatom, jatom, ialpha
        real, dimension (3) :: r21, rnabc   !< vectors
        real, dimension (3) :: sighat       !< unit vector along r2 - r1
        real, dimension (3) :: rhat         !< unit vector along bc - rna

        real, dimension (3, 3, 3) :: depsA  !< the Depsilon matrix for the bc
        real, dimension (3, 3, 3) :: depsB  !< the Depsilon matrix for the na

        real, dimension (3) :: amt, bmt

! bcxcm = Hartree matrix in molecular coordixctes
! bcxcx = Hartree matrix in crystal coordinates
! dbcxcm = derivative of Hartree matrix in molecular coordinates
! vdbcxcm = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcxcx = vectorized derivative of Hartree matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcxcm

        real, dimension (:, :), allocatable :: dpbcxcm
        real, dimension (:, :), allocatable :: dxbcxcm
        real, dimension (:, :), allocatable :: dybcxcm

        real, dimension (:, :, :), allocatable :: f3xcMa
        real, dimension (:, :, :), allocatable :: f3xcMb

        real, dimension (:, :, :), allocatable :: f3xcXa
        real, dimension (:, :, :), allocatable :: f3xcXb
        real, dimension (:, :, :), allocatable :: f3xcXc

        type(T_Fdata_cell_3c), pointer :: pFdata_cell
        type(T_Fdata_bundle_3c), pointer :: pFdata_bundle

        ! forces
        type(T_forces), pointer :: pfalpha
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================
        ! needed for charge transfer bits
        allocate (dqorb (nspecies))
        allocate (Q0 (s%natoms))
        allocate (Q (s%natoms))
        allocate (dQ (s%natoms))

! Procedure
! ===========================================================================
! Initialize the charge transfer bit
        do ispecies = 1, nspecies
          dqorb(ispecies) = 0.5d0
          if (species(ispecies)%nssh .eq. 1) dqorb(ispecies) = 0.25d0
        end do

! Calculate nuclear charge.
! ****************************************************************************
! We do a quadratic expansion:
!       f(q) = f(0)  +  f'(0)*q   +  1/2 f''(0)*q*q
!
! The derivatives are computed as:
!       f'(0)  = [ f(dq) - f(-dq) ] / 2dq
!       f''(0) = [ f(dq) - 2f(0) + f(-dq) ] / dq*dq
!
! We introduce linear factors:
!              linfac(0)   =  1.0
!              linfac(i)   =  -/+ * (1/2) *  q/dq        i = 1, 2
!
!    and quadratic factors:
!              quadfac(0)  = -2.0 * (1/2) * (q/dq)**2
!              quadfac(i)  =  1.0 * (1/2) * (q/dq)**2    i=1,2
!
!       f(0) = f(dq=0) ; f(1) = f(-dq) ; f(2) = f(dq)
!
! With this, f(q) is given as:
!       f(q) = sum_i  (linfac(i) + quadfac(i))*f(i)
!
! ****************************************************************************
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          Q0(iatom) = 0.0d0
          Q(iatom) = 0.0d0
          do issh = 1, species(in1)%nssh
            Q0(iatom) = Q0(iatom) + species(in1)%shell(issh)%Qneutral
            Q(iatom) = Q(iatom) + s%atom(iatom)%shell(issh)%Qin
          end do
          dQ(iatom) = Q(iatom) - Q0(iatom)
        end do

! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          indna = s%atom(ialpha)%imass
          rna = s%atom(ialpha)%ratom

          ! cut some lengthy notation
          pfalpha=>s%forces(ialpha)

          ! loop over the common neighbor pairs of ialp
          do ineigh = 1, s%neighbors(ialpha)%ncommon
            mneigh = s%neighbors(ialpha)%neigh_common(ineigh)
            if (mneigh .ne. 0) then
              iatom = s%neighbors(ialpha)%iatom_common_j(ineigh)
              ibeta = s%neighbors(ialpha)%iatom_common_b(ineigh)
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              in1 = s%atom(iatom)%imass
              norb_mu = species(in1)%norb_max

              jatom = s%neighbors(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%norb_max

              ! cut lengthy notation
              pfi=>s%forces(iatom); pfj=>s%forces(jatom)

              ! density matrix
              pdenmat=>s%denmat(iatom); pRho_neighbors=>pdenmat%neighbors(mneigh)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
              r21 = r2 - r1
              z = distance (r1, r2)

              ! unit vector in sigma direction.
              if (z .lt. 1.0d-05) then
                sighat(1) = 0.0d0
                sighat(2) = 0.0d0
                sighat(3) = 1.0d0
              else
                sighat = r21/z
              end if

! ****************************************************************************
! Find rnabc = vector pointing from center of bondcharge to rna
! This gives us the distance dnabc.
              rnabc = rna - (r1 + r21/2.0d0)
              x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)

              ! unit vector in rnabc direction.
              if (x .lt. 1.0d-05) then
                rhat(1) = 0.0d0
                rhat(2) = 0.0d0
                rhat(3) = 0.0d0
              else
                rhat = rnabc/x
              end if
              cost = dot_product(sighat, rhat)
              call epsilon_function (rhat, sighat, eps)

! dera3 = depsA = deps/dratm in the 3-center system
! der13 = depsB = deps/dr1 in the 3-center system
              call Depsilon_3c (r1, r2, r21, z, rna, rnabc, eps, depsA, depsB)

! The first piece will be the force with respect to atom 3.
              if (x .gt. 1.0d-5) then
                amt = (sighat - cost*rhat)/x
              else
                amt = 0.0d0
              end if

! The second piece will be the force with respect to atom 1.
              bmt = (cost*sighat - rhat)/z

! For now we just do the neutral atom interactions.
! Charged atom interactions are assembled in assemble_ca_3c.f
! So set ideriv = 0 within this subroutine.
!
!              interaction    subtypes     index
!
!      xc3c         4          0..9(max)   1..10

! ****************************************************************************
!
! Get the D-matrices from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_xc3c
              ideriv = 0

! Allocate block arrays
              allocate (bcxcm (norb_mu,norb_nu)); bcxcm = 0.00
              allocate (dpbcxcm (norb_mu,norb_nu)); dpbcxcm = 0.00
              allocate (dxbcxcm (norb_mu,norb_nu)); dxbcxcm = 0.00
              allocate (dybcxcm (norb_mu,norb_nu)); dybcxcm = 0.00
              call getDMEs_Fdata_3c (in1, in2, indna, interaction, ideriv, x, &
     &                               z, norb_mu, norb_nu, cost, rhat, sighat, &
     &                               bcxcm, dpbcxcm, dxbcxcm, dybcxcm)

              allocate (f3xcMa(3,norb_mu,norb_nu)); f3xcMa = 0.0d0
              allocate (f3xcMb(3,norb_mu,norb_nu)); f3xcMb = 0.0d0
              allocate (f3xcXa(3,norb_mu,norb_nu)); f3xcXa = 0.0d0
              allocate (f3xcXb(3,norb_mu,norb_nu)); f3xcXb = 0.0d0
              allocate (f3xcXc(3,norb_mu,norb_nu)); f3xcXc = 0.0d0

! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.
              pFdata_bundle => Fdata_bundle_3c(in1, in2, indna)
              pFdata_cell =>                                                &
     &          pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(interaction,ideriv,1))

              ! loop over matrix elements
              do iindex = 1, pFdata_cell%nME
                imu = pFdata_cell%mu_3c(iindex)
                inu = pFdata_cell%nu_3c(iindex)

! Now recover f3xcMa which is a three-dimensional array
                f3xcMa(:,imu,inu) = rhat*dxbcxcm(imu,inu) + amt*dpbcxcm(imu,inu)

! Now recover f3xcMb which is a three-dimensional array
                f3xcMb(:,imu,inu) = - sighat*dybcxcm(imu,inu)               &
     &            + bmt*dpbcxcm(imu,inu) - f3xcMa(:,imu,inu)/2.0d0
              end do ! end loop over matrix elements

! ***************************************************************************
! Convert to Crystal Coordinates
! ***************************************************************************
! The call to rotated does the rotations to crystal coordinates of these
! force things.
!
! For example:
! Suppose we have f_A(3,mu,nu), which is d/dratm M(mu,nu) where M(mu,nu)
! is in molecular. To transform M(mu,nu) to crystal, we need Udag * M * U.
! Therefore, f_A(3,mu,nu)[CRYSTAL] = (d/dratm Udag) * M * U
!                                   + Udag * M * (d/dratm U)
!                                   + Udag * f_A * U.
!
! So, to use this baby, put in deps3c (deps/dr1, deps/dr2, deps/dratm),
! and f_A and M.
!
! NOTE: rotated works on the assumption that we are adding derivatives,
! NOT forces. So f3xcMa,... etc. MUST not yet be forcelike.
! We do the - sign for forces at the end.
! ***************************************************************************
! Force on the neutral atom with respect to atom 3 (f3xcMa).
              call Drotate (in1, in2, eps, depsA, norb_mu, norb_nu,          &
     &                      bcxcm, f3xcMa, f3xcXa)

! Force on the neutral atom with respect to atom 1 (f3xcMb).
              call Drotate (in1, in2, eps, depsB, norb_mu, norb_nu,          &
     &                      bcxcm, f3xcMb, f3xcXb)

! Make things force-like and determine f3xcXc, whcih is found from Newtons Laws:
!             f3xcXa = - f3xcXa
!             f3xcXb = - f3xcXb
              f3xcXc = - f3xcXa - f3xcXb

              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfalpha%f3xca = pfalpha%f3xca                               &
      &             - pRho_neighbors%block(imu,inu)*f3xcXa(:,imu,inu)
                  pfi%f3xcb = pfi%f3xcb                                       &
      &             - pRho_neighbors%block(imu,inu)*f3xcXb(:,imu,inu)
                  pfj%f3xcc = pfj%f3xcc                                       &
      &             - pRho_neighbors%block(imu,inu)*f3xcXc(:,imu,inu)
                end do
              end do

! ****************************************************************************
! DO CHARGE DERIVATIVE CASES HERE.
! ****************************************************************************
! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
              ideriv_min = 1
              ideriv_max = 6

              dQ_factor(1) = -dQ(iatom)/(2.0d0*dqorb(in1))
              dQ_factor(2) =  dQ(iatom)/(2.0d0*dqorb(in1))
              dQ_factor(3) = -dQ(jatom)/(2.0d0*dqorb(in2))
              dQ_factor(4) =  dQ(jatom)/(2.0d0*dqorb(in2))
              dQ_factor(5) = -dQ(ialpha)/(2.0d0*dqorb(indna))
              dQ_factor(6) =  dQ(ialpha)/(2.0d0*dqorb(indna))
              do ideriv = ideriv_min, ideriv_max
                call getDMEs_Fdata_3c (in1, in2, indna, interaction, ideriv,  &
     &                                 x, z, norb_mu, norb_nu, cost, rhat,    &
     &                                 sighat, bcxcm, dpbcxcm, dxbcxcm, dybcxcm)

! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.
                pFdata_bundle => Fdata_bundle_3c(in1, in2, indna)
                pFdata_cell =>                                                &
     &            pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(interaction,ideriv,1))

                ! loop over matrix elements
                do iindex = 1, pFdata_cell%nME
                  imu = pFdata_cell%mu_3c(iindex)
                  inu = pFdata_cell%nu_3c(iindex)

! Now recover f3xcMa which is a three-dimensional array
                  f3xcMa(:,imu,inu) = rhat*dxbcxcm(imu,inu) + amt*dpbcxcm(imu,inu)

! Now recover f3xcMb which is a three-dimensional array
                  f3xcMb(:,imu,inu) = - sighat*dybcxcm(imu,inu)               &
     &              + bmt*dpbcxcm(imu,inu) - f3xcMa(:,imu,inu)/2.0d0
                end do ! end loop over matrix elements

! ***************************************************************************
! Convert to Crystal Coordinates
! ***************************************************************************
! The call to rotated does the rotations to crystal coordinates of these
! force things.
!
! For example:
! Suppose we have f_A(3,mu,nu), which is d/dratm M(mu,nu) where M(mu,nu)
! is in molecular. To transform M(mu,nu) to crystal, we need Udag * M * U.
! Therefore, f_A(3,mu,nu)[CRYSTAL] = (d/dratm Udag) * M * U
!                                   + Udag * M * (d/dratm U)
!                                   + Udag * f_A * U.
!
! So, to use this baby, put in deps3c (deps/dr1, deps/dr2, deps/dratm),
! and f_A and M.
!
! NOTE: rotated works on the assumption that we are adding derivatives,
! NOT forces. So f3xcMa,... etc. MUST not yet be forcelike.
! We do the - sign for forces at the end.
! ***************************************************************************
! Force on the neutral atom with respect to atom 3 (f3xcMa).
                call Drotate (in1, in2, eps, depsA, norb_mu, norb_nu,         &
     &                        bcxcm, f3xcMa, f3xcXa)

! Force on the neutral atom with respect to atom 1 (f3xcMb).
                call Drotate (in1, in2, eps, depsB, norb_mu, norb_nu,         &
     &                        bcxcm, f3xcMb, f3xcXb)

! Make things force-like and determine f3xcXc, whcih is found from Newtons Laws:
!               f3xcXa = - f3xcXa
!               f3xcXb = - f3xcXb
                f3xcXc = - f3xcXa - f3xcXb

                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfalpha%f3xca = pfalpha%f3xca                               &
      &               - pRho_neighbors%block(imu,inu)*dQ_factor(ideriv)*f3xcXa(:,imu,inu)
                    pfi%f3xcb = pfi%f3xcb                                       &
      &               - pRho_neighbors%block(imu,inu)*dQ_factor(ideriv)*f3xcXb(:,imu,inu)
                    pfj%f3xcc = pfj%f3xcc                                       &
      &               - pRho_neighbors%block(imu,inu)*dQ_factor(ideriv)*f3xcXc(:,imu,inu)
                  end do
                end do
              end do
              deallocate (bcxcm, dpbcxcm, dxbcxcm, dybcxcm)
              deallocate (f3xcMa, f3xcMb)
              deallocate (f3xcXa, f3xcXb, f3xcXc)
            end if ! if (mneigh .ne. 0)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
        deallocate (dqorb, Q0, Q, dQ)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_vxc_3c

! ===========================================================================
! destroy_Dassemble_vxc
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the Dassemble_2c
!! information.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_Dassemble_vxc (s)
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

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          deallocate (s%forces(iatom)%vxc_on_site)
          deallocate (s%forces(iatom)%vxc_off_site)
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
        end subroutine destroy_Dassemble_vxc


! End Module
! ===========================================================================
        end module M_Dassemble_vxc
