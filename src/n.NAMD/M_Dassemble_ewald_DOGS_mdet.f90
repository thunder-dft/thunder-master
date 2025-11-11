! copyright info:
!
!                             @Copyright 2025
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

! M_assemble_ewald
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the interactions for the short range and long range
!! ewald interactions.
!!
!! It contains the following subroutines within the module:
!!
!!       Dassemble_ewaldsr.f90 - assemble the short-range ewald matrix
!!       Dassemble_ewaldlr.f90 - assemble the long-range ewald matrix
!!       Dassemble_ewald.f90 - builds the ewald interactions - atom pair
!!       destroy_Dassemble_ewald.f90 - destroy the arrays
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_Dassemble_ewald

! /GLOBAL
        use M_precision
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones

! /FDATA
        use M_Fdata_2c

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! Dassemble_ewaldsr.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates forces for the short range
!> ewald (Coulomb) interactions. These will be used to offset the total
!> ewald interactions which are calculated in the long-range ewald
!> subroutine.
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
        subroutine Dassemble_ewaldsr (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha, iatom, jatom, katom   !< the three parties involved
        integer ibeta, jbeta           !< cells for three atoms
        integer ineigh, mneigh         !< counter over neighbors
        integer in1, in2, indna        !< species numbers
        integer issh                   !< counter over shells

        integer num_neigh              !< number of neighbors
        integer matom                  !< matom is the self-interaction atom
        integer mbeta                  !< the cell containing neighbor of iatom

        integer imu, inu, jnu                !< counter over MEs
        integer norb_mu, norb_nu       !< size of the block for the pair

        integer mmu, nnu               !< counter over coefficients of wavefunctions
        integer iband, jband           !< counter over transitions
        integer ikpoint                !< counter over kpoints
        integer nbands                 !< number of bands        

        real distance_13, distance_23  !< distance from 3rd atom
        real z                         !< distance between r1 and r2
        real x                         !< dnabc

        real dot                        !< dot product between K and r
        real cmunu                !< density matrix elements for mdet

! The charge transfer bit = need a +-dq on each orbital.
        real, allocatable :: dQ (:)       !< charge on atom, i.e. ionic
        real, allocatable :: Q0 (:)       !< total neutral atom charge, i.e. Ztot
        real, allocatable :: Q (:)        !< total charge on atom

        real, dimension (3) :: r1, r2, rna  !< positions - iatom, jatom, ialpha
        real, dimension (3) :: r21, rnabc   !< vectors
        real, dimension (3) :: sighat   !< unit vector along r2 - r1
        real, dimension (3) :: rhatA1    !< unit vector along rna - r1
        real, dimension (3) :: rhatA2    !< unit vector along rna - r2

        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        real, dimension (:, :), allocatable :: dterm
        real, dimension (:, :, :), allocatable :: dpterm
        real, dimension (:, :), allocatable :: sterm
        real, dimension (:, :, :), allocatable :: spterm

        ! derivatives of three-center overlap-dipole ewald interactions
        real, dimension (:, :, :), allocatable :: demnplA
        real, dimension (:, :, :), allocatable :: demnplB
        real, dimension (:, :, :), allocatable :: demnplC

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pdip_neighbors
        type(T_assemble_neighbors), pointer :: pdipole_z

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_block), pointer :: pRho_neighbors_matom

        ! forces
        type(T_forces), pointer :: pfalpha
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

        ! NAC Stuff
        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Allocate Arrays
! ===========================================================================
! We build the ewald forces here, so we allocate and initialize
        do iatom = 1, s%natoms
          num_neigh = s%neighbors(iatom)%neighn
          nullify (pfi)
          pfi=>s%forces(iatom)
          allocate (s%forces(iatom)%ewaldsr (3, num_neigh)); pfi%ewaldsr = 0.0d0
        end do

        ! needed for charge transfer bits
        allocate (Q0 (s%natoms))
        allocate (Q (s%natoms))
        allocate (dQ (s%natoms))

!============================================================================
! NAC ewald short range case for two-center ewald gradient
! NAC Zhaofa Li have changed inu to jnu to match the formula in
! J. Chem. Phys. 138, 154106 (2013)
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize the charge transfer bit
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

! T W O - C E N T E R   O V E R L A P   M O N O P O L E   P I E C E
!****************************************************************************
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          nullify (poverlap, pS_neighbors)
          poverlap=>s%overlap(iatom); pS_neighbors=>poverlap%neighbors(matom)

          ! density matrix
          nullify (pdenmat, pRho_neighbors)
          pdenmat=>s%denmat(iatom)
          pRho_neighbors_matom=>pdenmat%neighbors(matom)

          ! cut some lengthy notation
          nullify (pfi); pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! SET-UP STUFF
! ***************************************************************************
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

! Get the overlap matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcnam). Divide by the distance between the centers
! and store this into the short-range ewald piece - ewaldsr.

! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction and we do nothing here.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction does not exist.
            else

              ! force on iatom, ineigh
              do inu = 1, norb_mu
                do imu = 1, norb_mu
                  pfi%ewaldsr(:,ineigh) = pfi%ewaldsr(:,ineigh)              &
      &             - P_eq2*pRho_neighbors_matom%block(imu,inu)*dQ(jatom)    &
      &                    *(pS_neighbors%block(imu,inu)/z**2)*sighat(:)
                end do
              end do

!============================================================================
! NAC derivative of ewald short range case for two-center overlap monopole
! which follows atom case
! NAC Zhaofa Li have changed inu to jnu to match the formula in
! J. Chem. Phys. 138, 154106 (2013)
! ===========================================================================
              do ikpoint = 1, s%nkpoints

                ! Cut some lengthy notation
                nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

                do iband = 1, pkpoint%nbands
                
                  ! Cut some lengthy notation
                  nullify (piband); piband=>pkpoint%transition(iband)

                  do jband = iband + 1, pkpoint%nbands

                    ! Cut some lengthy notation
                    nullify (pjband); pjband=>pkpoint%transition(jband)

                    do jnu = 1, norb_mu
                      nnu = jnu + s%iblock_slot(jatom)
                      step1 = pjband%c_mdet(nnu)
                      do imu = 1, norb_mu
                        mmu = imu + s%iblock_slot(iatom)
                        step2 = step1*conjg(piband%c_mdet(mmu))
                        cmunu = real(step2)
                        
                        piband%dij(:,iatom,jband) = piband%dij(:,iatom,jband) &
     &                    + P_eq2*cmunu*dQ(jatom)*(pS_neighbors%block(imu,jnu)/z**2)*sighat(:)
                        piband%dij(:,jatom,jband) = piband%dij(:,jatom,jband) &
     &                    - P_eq2*cmunu*dQ(iatom)*(pS_neighbors%block(imu,jnu)/z**2)*sighat(:)
                      end do ! end loop over imu
                    end do ! end loop over jnu
                  end do ! end loop over jband
                end do ! end loop over iband
              end do ! end loop over kpoints
! ===========================================================================
              
            end if  ! iatom .eq. jatom
          end do ! end loop over neighbors
          nullify (pfi)
          nullify (poverlap, pS_neighbors)
          nullify (pdenmat, pRho_neighbors)
        end do ! end loop over atoms

! T W O - C E N T E R   O V E R L A P   O N T O P   D I P O L E    P I E C E
!****************************************************************************
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          nullify (poverlap, pdipole_z)
          poverlap=>s%overlap(iatom)
          pdipole_z=>s%dipole_z(iatom)

          ! density matrix
          nullify (pdenmat); pdenmat=>s%denmat(iatom)

          ! force on iatom
          nullify (pfi); pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some more lengthy notation
            nullify (pS_neighbors, pdip_neighbors)
            pS_neighbors=>poverlap%neighbors(ineigh)
            pdip_neighbors=>pdipole_z%neighbors(ineigh)

            ! density matrix
            nullify (pRho_neighbors)
            pRho_neighbors=>pdenmat%neighbors(ineigh)

! SET-UP STUFF
! ***************************************************************************
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

! Get the overlap matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcnam). Divide by the distance between the centers
! and store this into the short-range ewald piece - ewaldsr.

! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction and we do nothing here.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction does not exist.
            else

! Find sum of charge for first atom and add to ewaldsr
              allocate (sterm (norb_mu, norb_nu))
              allocate (spterm (3, norb_mu, norb_nu))
              allocate (dterm (norb_mu, norb_nu))
              allocate (dpterm (3, norb_mu, norb_nu))

              sterm = pS_neighbors%block/(2.0d0*z)
              spterm = pS_neighbors%Dblock/(2.0d0*z)
              dterm = pdip_neighbors%block/(z**2)
              dpterm = pdip_neighbors%Dblock/(z**2)

              ! force on iatom
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfi%ewaldsr(:,ineigh) = pfi%ewaldsr(:,ineigh)               &
     &             - P_eq2*pRho_neighbors%block(imu,inu)                      &
     &                    *dQ(iatom)*(spterm(:,imu,inu) + dpterm(:,imu,inu))  &
     &             - P_eq2*pRho_neighbors%block(imu,inu)                      &
     &                    *dQ(iatom)*(sterm(imu,inu)                          &
     &                                + 2.0d0*dterm(imu,inu))*(sighat(:)/z)   &
     &             - P_eq2*pRho_neighbors%block(imu,inu)                      &
     &                    *dQ(jatom)*(spterm(:,imu,inu) - dpterm(:,imu,inu))  &
     &             - P_eq2*pRho_neighbors%block(imu,inu)                      &
     &                    *dQ(jatom)*(sterm(imu,inu)                          &
     &                                - 2.0d0*dterm(imu,inu))*(sighat(:)/z)
                end do
              end do

!============================================================================
! NAC derivative of ewald short range case for 2c overlap dipole
! which follows ontop case
! NAC Zhaofa Li have changed inu to jnu to match the formula in
! J. Chem. Phys. 138, 154106 (2013)
! ===========================================================================
                do ikpoint = 1, s%nkpoints
                                  
                  ! Cut some lengthy notation
                  nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)
  
                  ! phase for non-gamma kpoints
                  vec = r2 - r1
                  sks = s%kpoints(ikpoint)%k
                  dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
                  phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
    
                  do iband = 1, pkpoint%nbands
                  
                    ! Cut some lengthy notation
                    nullify (piband); piband=>pkpoint%transition(iband)

                    do jband = iband + 1, pkpoint%nbands

                      ! Cut some lengthy notation
                      nullify (pjband); pjband=>pkpoint%transition(jband)

                      do jnu = 1, norb_nu
                        phase = phasex
                        nnu = jnu + s%iblock_slot(jatom)
                        step1 = phase*pjband%c_mdet(nnu)
                        do imu = 1, norb_mu
                          mmu = imu + s%iblock_slot(iatom)
                          step2 = step1*conjg(piband%c_mdet(mmu))
                          cmunu = real(step2)
                          
                          piband%dij(:,iatom,jband) = piband%dij(:,iatom,jband) &
     &                      + P_eq2*cmunu*dQ(iatom)*(spterm(:,imu,jnu)       &
                                                     + dpterm(:,imu,jnu))    &
     &                      + P_eq2*cmunu*dQ(iatom)*(sterm(imu,jnu)          &
     &                                               + 2.0d0*dterm(imu,jnu))*(sighat(:)/z) &
     &                      + P_eq2*cmunu*dQ(jatom)*(spterm(:,imu,jnu)       &
     &                                               - dpterm(:,imu,jnu))    &
     &                      + P_eq2*cmunu*dQ(jatom)*(sterm(imu,jnu)          &
     &                                               - 2.0d0*dterm(imu,jnu))*(sighat(:)/z)

                          piband%dij(:,jatom,jband) = piband%dij(:,jatom,jband) &
     &                      - P_eq2*cmunu*dQ(jatom)*(spterm(:,imu,jnu)       &
     &                                               - dpterm(:,imu,jnu))    &
     &                      - P_eq2*cmunu*dQ(jatom)*(sterm(imu,jnu)          &
     &                                               - 2.0d0*dterm(imu,jnu))*(sighat(:)/z) &
     &                      - P_eq2*cmunu*dQ(iatom)*(spterm(:,imu,jnu)       &
     &                                               + dpterm(:,imu,jnu))    &
     &                      - P_eq2*cmunu*dQ(iatom)*(sterm(imu,jnu)          &
     &                                               + 2.0d0*dterm(imu,jnu))*(sighat(:)/z)
                      end do ! end loop over imu
                    end do ! end loop over jnu
                  end do ! end loop over jband
                end do ! end loop over iband
              end do ! end loop over kpoints
! =================================================================================
              
              deallocate (sterm, spterm)
              deallocate (dterm, dpterm)
            end if
            nullify (pRho_neighbors)
            nullify (pS_neighbors, pdip_neighbors)
          end do ! end loop over neighbors
          nullify (pfi)
          nullify (poverlap, pdipole_z)
          nullify (pdenmat)
        end do ! end loop over atoms

!****************************************************************************
! T H R E E - C E N T E R   O V E R L A P   A N D    D I P O L E    P I E C E
!****************************************************************************
! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          indna = s%atom(ialpha)%imass
          rna = s%atom(ialpha)%ratom

          ! cut some lengthy notation
          nullify (pfalpha)
          pfalpha=>s%forces(ialpha)

          ! loop over the common neigbor pairs of ialp
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
              nullify (pfi, pfj)
              pfi=>s%forces(iatom); pfj=>s%forces(jatom)

              ! density matrix
              nullify (pdenmat, pRho_neighbors)
              pdenmat=>s%denmat(iatom); pRho_neighbors=>pdenmat%neighbors(mneigh)

              nullify (poverlap, pS_neighbors)
              poverlap=>s%overlap(iatom); pS_neighbors=>poverlap%neighbors(mneigh)
              nullify (pdipole_z, pdip_neighbors)
              pdipole_z=>s%dipole_z(iatom); pdip_neighbors=>pdipole_z%neighbors(mneigh)

! SET-UP STUFF
! ***************************************************************************
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

! Find rnabc = vector pointing from center of bondcharge to rna
! This gives us the distance dnabc (or x value in the 2D grid).
              rnabc = rna - (r1 + r21/2.0d0)
              x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)

! Find other distances -
              distance_13 = distance (rna, r1)
              distance_23 = distance (rna, r2)

! Find the unit vector in rna-r1 direction.
              if (distance_13 .gt. 1.0d-05) then
                rhatA1(:) = (rna(:) - r1(:))/distance_13
              end if

! Find the unit vector in rna-r2 direction.
              if (distance_23 .gt. 1.0d-05) then
                rhatA2(:) = (rna(:) - r2(:))/distance_23
              end if

! Now combine the pieces together to get the correct short-range Ewald
! summation. Allocate the size of dterm and sterm according to size
! of the blocks
              allocate (sterm (norb_mu, norb_nu))
              allocate (spterm (3, norb_mu, norb_nu))
              allocate (dterm (norb_mu, norb_nu))
              allocate (dpterm (3, norb_mu, norb_nu))

              allocate (demnplA (3, norb_mu, norb_nu))
              allocate (demnplB (3, norb_mu, norb_nu))
              allocate (demnplC (3, norb_mu, norb_nu))

              sterm = pS_neighbors%block/2.0d0
              spterm = pS_neighbors%Dblock/2.0d0
              dterm = pdip_neighbors%block/z
              dpterm = pdip_neighbors%Dblock/z

! Now the derivatives
! A = 3, B = 1, C = 2 <B|A|C> This is NOT force-like.

! demnplA: The net effect is a minus sign.
! There are three minus signs in all - one belonging to the sighat term,
! one belonging to the 1/y term and a minus subtraction error.
! No minus signs are "force-like". This is not a force like derivative,
! but a regular derivative d/dr1.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  dpterm(:,imu,inu) = dpterm(:,imu,inu) + dterm(imu,inu)*sighat(:)/z
                  demnplA(:,imu,inu) =                                        &
     &             - dQ(ialpha)*(sterm(imu,inu)                               &
     &                          - dterm(imu,inu))*rhatA1(:)/distance_13**2    &
     &             - dQ(ialpha)*(sterm(imu,inu)                               &
     &                          + dterm(imu,inu))*rhatA2(:)/distance_23**2
                  demnplB(:,imu,inu) =                                        &
                   + dQ(ialpha)*(sterm(imu,inu)                               &
     &                          - dterm(imu,inu))*rhatA1(:)/distance_13**2    &
     &             + dQ(ialpha)*(spterm(:,imu,inu)                            &
     &                           - dpterm(:,imu,inu))/distance_13             &
     &             + dQ(ialpha)*(spterm(:,imu,inu)                            &
     &                           + dpterm(:,imu,inu))/distance_23
                end do
              end do
              ! terms are (-) because we subtract ewaldsr from the Hamiltonian
              demnplA = - demnplA
              demnplB = - demnplB

! By Newton's laws demnplC = - demnplA - demnplB
              demnplC = - demnplA - demnplB

              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfalpha%f3naa = pfalpha%f3naa                                &
      &             - P_eq2*pRho_neighbors%block(imu,inu)*demnplA(:,imu,inu)
                  pfi%f3nab = pfi%f3nab                                        &
      &             - P_eq2*pRho_neighbors%block(imu,inu)*demnplB(:,imu,inu)
                  pfj%f3nac = pfj%f3nac                                        &
      &             - P_eq2*pRho_neighbors%block(imu,inu)*demnplC(:,imu,inu)
                end do
              end do

!============================================================================
! NAC derivative of ewald short range case for three-center overlap dipole
! NAC Zhaofa Li have changed inu to jnu to match the formula in
! J. Chem. Phys. 138, 154106 (2013)
! ===========================================================================
              do ikpoint = 1, s%nkpoints
                                
                ! Cut some lengthy notation
                nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

                ! phase for non-gamma kpoints
                vec = r2 - r1
                sks = s%kpoints(ikpoint)%k
                dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
                phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
   
                do iband = 1, pkpoint%nbands

                  ! Cut some lengthy notation
                  nullify (piband); piband=>pkpoint%transition(iband)

                  do jband = iband + 1, pkpoint%nbands

                    ! Cut some lengthy notation
                    nullify (pjband); pjband=>pkpoint%transition(jband)

                    do jnu = 1, norb_nu
                      phase = phasex
                      nnu = jnu + s%iblock_slot(jatom)
                      step1 = phase*pjband%c_mdet(nnu)
                      do imu = 1, norb_mu
                        mmu = imu + s%iblock_slot(iatom)
                        step2 = step1*conjg(piband%c_mdet(mmu))
                        cmunu = real(step2)
                        
                        piband%dij(:,ialpha,jband) =                         &
     &                    piband%dij(:,ialpha,jband) - cmunu*demnplA(:,imu,jnu)*P_eq2
                        piband%dij(:,iatom,jband) =                          &
     &                    piband%dij(:,iatom,jband) - cmunu*demnplB(:,imu,jnu)*P_eq2
                        piband%dij(:,jatom,jband) =                          &
     &                    piband%dij(:,jatom,jband) - cmunu*demnplC(:,imu,jnu)*P_eq2
                      end do ! end loop over imu
                    end do ! end loop over jnu
                  end do ! end loop over jband
                end do ! end loop over iband
              end do ! end loop over kpoints
! ===========================================================================

              deallocate (sterm, spterm)
              deallocate (dterm, dpterm)
              deallocate (demnplA, demnplB, demnplC)
            end if
            nullify (pdenmat, pRho_neighbors)
            nullify (poverlap, pS_neighbors)
            nullify (pdipole_z, pdip_neighbors)
          end do ! end loop over neighbors
          nullify (pfalpha)
          nullify (pfi, pfj)
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
        deallocate (Q0, Q, dQ)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_ewaldsr


! ===========================================================================
! Dassemble_ewaldlr.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the total ewald long-range matrix elements.
!> First, the routine assemble_ewald should be called to get the ewald
!> interactions and this routine then assembles those into matrix elements.
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
        subroutine Dassemble_ewaldlr (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer in1, in2                !< species numbers
        integer jatom, katom            !< neighbor of iatom
        integer issh                    !< counter over shells
        integer num_neigh               !< number of neighbors
        integer mbeta                   !< the cell containing neighbor of iatom

        integer imu, inu, jnu           !< counter over MEs
        integer norb_mu, norb_nu        !< size of the block for the pair
        integer mmu, nnu                 !< counter over coefficients of wavefunctions
        integer iband, jband             !< counter over transitions
        integer ikpoint                  !< counter over kpoints
        integer nbands                   !< number of bands

        integer logfile                 !< which unit to write output

        real z                          !< distance between r1 and r2

        real dot                        !< dot product between K and r
        real cmunu                !< density matrix elements for mdet

! The charge transfer bit = need a +-dq on each orbital.
        real, allocatable :: dQ (:)       !< charge on atom, i.e. ionic
        real, allocatable :: Q0 (:)       !< total neutral atom charge, i.e. Ztot
        real, allocatable :: Q (:)        !< total charge on atom

        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        real, allocatable, dimension (:, :) :: dterm
        real, allocatable, dimension (:, :, :) :: dpterm
        real, allocatable, dimension (:, :) :: sterm
        real, allocatable, dimension (:, :, :) :: spterm
        real, allocatable, dimension (:) :: sum_ewald
        real, allocatable, dimension (:, :) :: sum_dewald

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pdip_neighbors
        type(T_assemble_neighbors), pointer :: pdipole_z

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors

        ! forces
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj
        type(T_forces), pointer :: pfk

        ! NAC Stuff
        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Allocate Arrays
! ===========================================================================
        allocate (sum_ewald (s%natoms)); sum_ewald = 0.0d0
        allocate (sum_dewald (3, s%natoms)); sum_dewald = 0.0d0

        ! needed for charge transfer bits
        allocate (Q0 (s%natoms))
        allocate (Q (s%natoms))
        allocate (dQ (s%natoms))

! ============================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        write (logfile,*) ' Computing the ewald energy and forces.'
        call Dassemble_ewald (s)

! Initialize the charge transfer bit
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

! First a preliminary quantity. ewald(i,j) = sum(L) 1/(| bi-bj +/- L|
! We need to calculate SUMS of ewald sums in which the charge is included.
! In the nutshell we are calculating:
!       Sum_(i,L) q(i)/|b(i)-b(alpha)+L|
        do iatom = 1, s%natoms
          do jatom = 1, s%natoms
            sum_ewald(iatom) = sum_ewald(iatom) + dQ(jatom)*s%ewald(iatom,jatom)
            sum_dewald(:,iatom) = sum_dewald(:,iatom) + dQ(jatom)*s%dewald(:,iatom,jatom)
          end do
        end do

! Now the meat of the calculation.  Construct ewaldlr(mu,nu,i,m) ===>
! the matrix elements of the long-range parts of the Hamiltonian.
! We make matrix elements for the Long Range Ewald according to our theory:
! ewaldlr(mu,nu,ineigh,iatom) =
! {s(mu,nu,ineigh,iatom)/2}*SUM(j_basis)(Qin(jatom) - Qneutral(jatom))
!                                        *(ewald(iatom,jatom)
!                                          + ewald(ineigh,jatom))*eq2
! The value P_eq2 makes it into the units of eV.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          nullify (poverlap, pdipole_z)
          poverlap=>s%overlap(iatom)
          pdipole_z=>s%dipole_z(iatom)

          ! density matrix
          nullify (pdenmat, pfi)
          pdenmat=>s%denmat(iatom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of the atom i.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some more lengthy notation
            nullify (pS_neighbors, pdip_neighbors)
            pS_neighbors=>poverlap%neighbors(ineigh)
            pdip_neighbors=>pdipole_z%neighbors(ineigh)

            ! density matrix - neighbors
            nullify (pRho_neighbors, pfj)
            pRho_neighbors=>pdenmat%neighbors(ineigh)
            pfj=>s%forces(jatom)

! SET-UP STUFF
! ***************************************************************************
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

! Allocate the size of dterm and sterm according to size of the blocks
            allocate (sterm (norb_mu, norb_nu))
            allocate (spterm (3, norb_mu, norb_nu))
            allocate (dterm (norb_mu, norb_nu))
            allocate (dpterm (3, norb_mu, norb_nu))

! "Charge" on each atom of the bondcharge. We split the charge S and
! dipole p to be S/2-p/d on atom 1 and S/2+p/d on atom 2. If atom 1 is
! equal to atom 2, then drop the p term.
            sterm = pS_neighbors%block/2.0d0

! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction and we do nothing here.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then
              spterm = 0.0d0
              dterm = 0.0d0
              dpterm = 0.0d0

              ! force on iatom
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfi%ewaldlr = pfi%ewaldlr                                   &
     &              - P_eq2*pRho_neighbors%block(imu,inu)*sterm(imu,inu)*sum_dewald(:,iatom)
                end do
              end do

              ! force on jatom
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfj%ewaldlr = pfj%ewaldlr                                   &
     &              - P_eq2*pRho_neighbors%block(imu,inu)*sterm(imu,inu)*sum_dewald(:,jatom)
                end do
              end do

!============================================================================
! NAC derivative of ewald long range case
! NAC Zhaofa Li have changed inu to jnu to match the formula in
! J. Chem. Phys. 138, 154106 (2013)
! ===========================================================================
              do ikpoint = 1, s%nkpoints

                ! Cut some lengthy notation
                nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

                ! phase for non-gamma kpoints
                vec = r2 - r1
                sks = s%kpoints(ikpoint)%k
                dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
                phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
   
                do iband = 1, pkpoint%nbands

                  ! Cut some lengthy notation
                  nullify (piband); piband=>pkpoint%transition(iband)

                  do jband = iband + 1, pkpoint%nbands

                    ! Cut some lengthy notation
                    nullify (pjband); pjband=>pkpoint%transition(jband)

                    do jnu = 1, norb_nu
                      phase = phasex
                      nnu = jnu + s%iblock_slot(jatom)
                      step1 = phase*pjband%c_mdet(nnu)
                      do imu = 1, norb_mu
                        mmu = imu + s%iblock_slot(iatom)
                        step2 = step1*conjg(piband%c_mdet(mmu))
                        cmunu = real(step2)
                        
                        piband%dij(:,iatom,jband) = piband%dij(:,iatom,jband) &
     &                    - P_eq2*cmunu*sterm(imu,jnu)*sum_dewald(:,iatom)
                        piband%dij(:,jatom,jband) = piband%dij(:,jatom,jband) &
     &                    - P_eq2*cmunu*sterm(imu,jnu)*sum_dewald(:,jatom)
                      end do ! end loop over imu
                    end do ! end loop over jnu
                  end do ! end loop over jband
                end do ! end loop over iband
              end do ! end loop over kpoints
! ===========================================================================

            else
              spterm = pS_neighbors%Dblock/2.0d0
              dterm = pdip_neighbors%block/z
              do imu = 1, norb_mu
                do inu = 1, norb_nu
                  dpterm(:,imu,inu) = pdip_neighbors%Dblock(:,imu,inu)/z      &
     &                               + pdip_neighbors%block(imu,inu)*sighat(:)/z**2
                end do
              end do

              ! force on iatom
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfi%ewaldlr = pfi%ewaldlr                                    &
      &              - P_eq2*pRho_neighbors%block(imu,inu)                     &
      &               *(sterm(imu,inu) - dterm(imu,inu))*sum_dewald(:,iatom)   &
      &              - P_eq2*pRho_neighbors%block(imu,inu)                     &
      &               *(spterm(:,imu,inu) - dpterm(:,imu,inu))*sum_ewald(iatom)&
      &              - P_eq2*pRho_neighbors%block(imu,inu)                     &
      &               *(sterm(imu,inu) + dterm(imu,inu))                       &
      &               *dQ(iatom)*s%dewald(:,iatom,jatom)                       &
      &              - P_eq2*pRho_neighbors%block(imu,inu)                     &
      &               *(spterm(:,imu,inu) + dpterm(:,imu,inu))*sum_ewald(jatom)
                end do
              end do

              ! force on jatom
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfj%ewaldlr = pfj%ewaldlr                                    &
      &              - P_eq2*pRho_neighbors%block(imu,inu)                     &
      &               *(sterm(imu,inu)                                         &
      &                  - dterm(imu,inu))*dQ(jatom)*s%dewald(:,jatom,iatom)   &
      &              + P_eq2*pRho_neighbors%block(imu,inu)                     &
      &               *(spterm(:,imu,inu) - dpterm(:,imu,inu))*sum_ewald(iatom)&
      &              - P_eq2*pRho_neighbors%block(imu,inu)                     &
      &               *(sterm(imu,inu) + dterm(imu,inu))*sum_dewald(:,jatom)   &
      &              + P_eq2*pRho_neighbors%block(imu,inu)                     &
      &               *(spterm(:,imu,inu) + dpterm(:,imu,inu))*sum_ewald(jatom)
                end do
              end do

!============================================================================
! NAC derivative of ewald long range case
! NAC Zhaofa Li have changed inu to jnu to match the formula in
! J. Chem. Phys. 138, 154106 (2013)
! ===========================================================================
              do ikpoint = 1, s%nkpoints

                ! Cut some lengthy notation
                nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

                ! phase for non-gamma kpoints
                vec = r2 - r1
                sks = s%kpoints(ikpoint)%k
                dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
                phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
   
                do iband = 1, pkpoint%nbands

                  ! Cut some lengthy notation
                  nullify (piband); piband=>pkpoint%transition(iband)

                  do jband = iband + 1, pkpoint%nbands

                    ! Cut some lengthy notation
                    nullify (pjband); pjband=>pkpoint%transition(jband)

                    do jnu = 1, norb_nu
                      phase = phasex
                      nnu = jnu + s%iblock_slot(jatom)
                      step1 = phase*pjband%c_mdet(nnu)
                      do imu = 1, norb_mu
                        mmu = imu + s%iblock_slot(iatom)
                        step2 = step1*conjg(piband%c_mdet(mmu))
                        cmunu = real(step2)

                        piband%dij(:,iatom,jband) = piband%dij(:,iatom,jband) &
       &                 - P_eq2*cmunu*(sterm(imu,jnu)                        &
       &                                - dterm(imu,jnu))*sum_dewald(:,iatom) &
       &                 - P_eq2*cmunu*(spterm(:,imu,jnu)                     &
       &                                - dpterm(:,imu,jnu))*sum_ewald(iatom) &
       &                 - P_eq2*cmunu*(sterm(imu,jnu)                        &
       &                                + dterm(imu,jnu))*dQ(iatom)*s%dewald(:,iatom,jatom) &
       &                 - P_eq2*cmunu*(spterm(:,imu,jnu)                     &
       &                                + dpterm(:,imu,jnu))*sum_ewald(jatom)
                        piband%dij(:,jatom,jband) = piband%dij(:,jatom,jband) &
       &                 - P_eq2*cmunu*(sterm(imu,jnu)                        &
       &                                - dterm(imu,jnu))*dQ(jatom)*s%dewald(:,jatom,iatom) &
       &                 + P_eq2*cmunu*(spterm(:,imu,jnu)                     &
       &                                - dpterm(:,imu,jnu))*sum_ewald(iatom) &
       &                 - P_eq2*cmunu*(sterm(imu,jnu)                        &
       &                                + dterm(imu,jnu))*sum_dewald(:,jatom) &
       &                 + P_eq2*cmunu*(spterm(:,imu,jnu)                     &
       &                                + dpterm(:,imu,jnu))*sum_ewald(jatom)
                      end do ! end loop over imu
                    end do ! end loop over jnu
                  end do ! end loop over jband
                end do ! end loop over iband
              end do ! end loop over kpoints
! ===========================================================================

            end if ! end if for r1 .eq. r2 case

! ****************************************************************************
! There is a correction to the force on the "other" atoms, katom, which are
! implicit in the ewald sums. Here katom is neither iatom nor jatom.
! ****************************************************************************
            do katom = 1, s%natoms

              ! force on jatom
              nullify (pfk); pfk=>s%forces(katom)

              if (katom .ne. iatom .and. katom .ne. jatom) then
                ! force on katom
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfk%ewaldlr = pfk%ewaldlr                                 &
     &                - pRho_neighbors%block(imu,inu)*dQ(katom)*P_eq2         &
     &                 *(sterm(imu,inu) - dterm(imu,inu))*s%dewald(:,katom,iatom) &
     &                - pRho_neighbors%block(imu,inu)*dQ(katom)*P_eq2         &
     &                 *(sterm(imu,inu) + dterm(imu,inu))*s%dewald(:,katom,jatom)    
                  end do
                end do

!============================================================================
! NAC derivative of ewald long range case
! NAC Zhaofa Li have changed inu to jnu to match the formula in
! J. Chem. Phys. 138, 154106 (2013)
! ===========================================================================
                do ikpoint = 1, s%nkpoints

                  ! Cut some lengthy notation
                  nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

                  ! phase for non-gamma kpoints
                  vec = r2 - r1
                  sks = s%kpoints(ikpoint)%k
                  dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
                  phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight

                  do iband = 1, pkpoint%nbands

                    ! Cut some lengthy notation
                    nullify (piband); piband=>pkpoint%transition(iband)

                    do jband = iband + 1, pkpoint%nbands

                      ! Cut some lengthy notation
                      nullify (pjband); pjband=>pkpoint%transition(jband)

                      do jnu = 1, norb_nu
                        phase = phasex
                        nnu = jnu + s%iblock_slot(jatom)
                        step1 = phase*pjband%c_mdet(nnu)
                        do imu = 1, norb_mu
                          mmu = imu + s%iblock_slot(iatom)
                          step2 = step1*conjg(piband%c_mdet(mmu))
                          cmunu = real(step2)

                          piband%dij(:,katom,jband) = piband%dij(:,katom,jband) &
     &                      - cmunu*dQ(katom)*P_eq2*(sterm(imu,jnu)          &
     &                                               - dterm(imu,jnu))*s%dewald(:,katom,iatom) &
     &                      - cmunu*dQ(katom)*P_eq2*(sterm(imu,jnu)          &
     &                                               + dterm(imu,jnu))*s%dewald(:,katom,jatom)
                        end do ! end loop over imu
                      end do ! end loop over jnu
                    end do ! end loop over jband
                  end do ! end loop over iband
                end do ! end loop over kpoints
! ===========================================================================

              end if ! end if for katom not iatom or jatom
            end do ! end loop over katom
            deallocate (sterm, spterm)
            deallocate (dterm, dpterm)
            nullify (pS_neighbors, pdip_neighbors)
            nullify (pRho_neighbors, pfj)
          end do ! end loop over neighbors
          nullify (pdenmat, pfi)
          nullify (poverlap, pdipole_z)
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
        deallocate (Q0, Q, dQ)
        deallocate (sum_ewald, sum_dewald)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_ewaldlr


! ===========================================================================
! Dassemble_ewald.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the total ewald energy. These interactions
!> are then used to buld the long-range ewald matrix elements.
!
!       This routine calculates the Ewald sum for a crystal with a given
! basis. This is specially designed for molecules with a given dipole moment.
! See Ihm, Zunger, Cohen -- Momentum Space Formalism for Total Energy of Solids
! J. Phys C v.12 (79). The ewald sum calculated here is actually gamma(ewald)
! in paper. The terms are also found in M.T. Yin and M.l. Cohen, Phys Rev. B26,
! 3259 (1982).
!
! Output:
!       ewald(natoms,natom) = ewald gamma in eV for all the basis atoms.
!       ewald = gamma1 + gamma2 + gamma3 + gamma4 - vself
!
!       gamma1 = first term of eq. 21 (Yin, Cohen paper), sum over g term.
!       gamma2 = erf term of eq. 21, sum over l term
!       gamma3 = delta(s,s') term in eq. 21.
!       gamma4 = ztot1*ztot2 term in eq. 21.
!       dewald = ewald forces
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
        subroutine Dassemble_ewald (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom           !< the three parties involved
        integer in1                    !< species numbers

        integer ig1, ig2, ig3          !< counters for the g grid
        integer ig1mx, ig2mx, ig3mx    !< maximum number of the g grid
        integer il1, il2, il3          !< counters for real-space grid
        integer il1mx, il2mx, il3mx    !< maximum number of space grid
        integer issh                   !< counter over shells

        real argument, kappa, stuff
        real factor, factorf
        real g1mag2, g2mag2, g3mag2    !< magnitude's (squared) of g's
        real gmin2, gmax               !< maximum and minimum magnitude of g's
        real gdotb                     !< dot product of k-space and real space
        real qq                        !< charges
        real r1mag2, r2mag2, r3mag2    !< magnitude's (squared) of r's
        real rmin2, rmax               !< maximum and minimum magnitude of a's
        real z

! The charge transfer bit = need a +-dq on each orbital.
        real, allocatable :: Q0 (:)       !< total neutral atom charge, i.e. Ztot
        real, allocatable :: Q (:)        !< total charge on atom

        real, dimension (3) :: a1vec, a2vec, a3vec
        real, dimension (3) :: cvec    !< resulting vector
        real, dimension (3) :: g, g1, g2, g3

        interface
          function a_cross_b (a, b)
            real, dimension (3) :: a_cross_b
            real, intent(in), dimension (3) :: a, b
          end function a_cross_b

          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance

          function magnitude (a)
            real magnitude
            real, intent(in), dimension (3) :: a
          end function magnitude
        end interface

        ! forces
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

! Allocate Arrays
! ===========================================================================
        allocate (s%ewald (s%natoms, s%natoms)); s%ewald = 0.0d0
        allocate (s%dewald (3, s%natoms, s%natoms)); s%dewald = 0.0d0

        ! needed for charge transfer bits
        allocate (Q0 (s%natoms))
        allocate (Q (s%natoms))

! Procedure
! ===========================================================================
! Calculate nuclear charge.
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          Q0(iatom) = 0.0d0
          Q(iatom) = 0.0d0
          do issh = 1, species(in1)%nssh
            Q0(iatom) = Q0(iatom) + species(in1)%shell(issh)%Qneutral
            Q(iatom) = Q(iatom) + s%atom(iatom)%shell(issh)%Qin
          end do
        end do

! Initialize ewald to zero and initialize other quantities
        s%ewald = 0.0d0
        a1vec = s%lattice(1)%a
        a2vec = s%lattice(2)%a
        a3vec = s%lattice(3)%a

! Determine the reciprical lattice vectors. Do this by Ashcroft and Mermin
! physics.
! First get a2 X a3.
        cvec = a_cross_b (a2vec, a3vec)

! Next find the volume of the cell.
! NOTE: volcel actually has a sign in the above. At this point the sign is
! important since we form g vectors by dividing by a1 dot (a2 X a3).
! Oh, you say. what difference does it make if we change the sign of g.
! it makes no difference in principle.
        s%volume = a1vec(1)*cvec(1) + a1vec(2)*cvec(2) + a1vec(3)*cvec(3)
        g1(:) = 2.0d0*pi*cvec(:)/s%volume

! Next we get a3 X a1, and g2.
        cvec = a_cross_b (a3vec, a1vec)
        g2(:) = 2.0d0*pi*cvec(:)/s%volume

! Finally we get a1 X a2, and g3.
        cvec = a_cross_b (a1vec, a2vec)
        g3(:) = 2.0d0*pi*cvec(:)/s%volume
        s%volume = abs(s%volume)

! Initialize gmax. This determines how far we sum g1, g2, and g3. See below
! why gmax = 5.0d0 is a reasonable criterion.
        gmax = 5.0d0

! Initialize rmax. This determines how far we sum a1, a2, and a3. See below
! why rmax = 5.0d0 is a reasonable criterion. The parameters a1, a2, a3 are
! the direct lattice vectors.
        rmax = 5.0d0

! Determine the magnitude of the vectors g1, g2, g3, a1, a2, a3.
        g1mag2 = magnitude(g1)**2
        g2mag2 = magnitude(g2)**2
        g3mag2 = magnitude(g3)**2

        r1mag2 = magnitude(a1vec)**2
        r2mag2 = magnitude(a2vec)**2
        r3mag2 = magnitude(a3vec)**2

! ****************************************************************************
! The parameter kappa is adjustable, chosen to make the sum's converge rapidly.
! The sum over g converges as exp (-g**2/(4*kappa*kappa)), while the
! sum over l converges as exp (-r**2*kappa**2). Lets set the arguments equal
! to determine a reasonable kappa value. We set them equal for the smallest
! g value and the smallest l value.
! First find the smallest rmag.
        rmin2 = r1mag2
        if (r2mag2 .lt. rmin2) rmin2 = r2mag2
        if (r3mag2 .lt. rmin2) rmin2 = r3mag2

! Next find the smallest gmag.
        gmin2 = g1mag2
        if (g2mag2 .lt. gmin2) gmin2 = g2mag2
        if (g3mag2 .lt. gmin2) gmin2 = g3mag2

! Now set rmin2*kappa**2 = gmin22/(4*kappa**2) and solve for kappa.
        kappa = sqrt(sqrt(gmin2/(4.0d0*rmin2)))

! ****************************************************************************
! In gamma1 we must sum over g vectors. The decay is exp(-g**2/4*kappa*kappa).
! We require the exponent for a given direction in g-space to be gmax**2.
! For instance gmax = 5.0, corresponding to an exponent of gmax**2 = 25.0 seems
! to be a reasonable choice. This gives us g = ig1mx*g1 where
! ig1mx**2 g1**2/(4*kappa**2) = gmax**2. Solve for ig1mx, and add 1.0 for
! good measure.
        ig1mx = int(gmax * sqrt(4.0d0*kappa**2/g1mag2) + 1.0d0)

! Now we do the same thing for g2 and g3
        ig2mx = int(gmax * sqrt(4.0d0*kappa**2/g2mag2) + 1.0d0)
        ig3mx = int(gmax * sqrt(4.0d0*kappa**2/g3mag2) + 1.0d0)

        if (ig1mx .le. 1) ig1mx = 2
        if (ig2mx .le. 1) ig2mx = 2
        if (ig3mx .le. 1) ig3mx = 2

! In gamma2 we must sum over l vectors. The asymptotic decay is
! exp(-kappa*kappa*r**2). We require the exponent for a given direction
! in r-space to be rmax**2. For instance rmax = 5.0, corresponding to an
! exponent of rmax**2 = 25.0 seems to be a reasonable choice. This gives us
! r = il1mx*a1 where il1mx**2 r**2 * kappa**2 = rmax**2. Solve for ir1mx, and
! add 1.0 for good measure.
        il1mx = int(rmax * sqrt(1.0d0/(kappa**2*r1mag2)) + 1.0d0)

! Now we do the same thing for r2 and r3
        il2mx = int(rmax * sqrt(1.0d0/(kappa**2*r2mag2)) + 1.0d0)
        il3mx = int(rmax * sqrt(1.0d0/(kappa**2*r3mag2)) + 1.0d0)

        if (il1mx .le. 1) il1mx = 2
        if (il2mx .le. 1) il2mx = 2
        if (il3mx .le. 1) il3mx = 2

! The real answer: now compute gamma ewald.
! ***********************************************************************
! Compute gamma1:
! ***********************************************************************
! Initialize fewald1

! Sum over g vectors.  If we are doing only a cluster, then only the gamma
! point is considered in the sum.
        if (s%icluster .eq. 1) then
          ig1mx = 0
          ig2mx = 0
          ig3mx = 0
        end if
        do ig1 = -ig1mx, ig1mx
          do ig2 = -ig2mx, ig2mx
            do ig3 = -ig3mx, ig3mx

! skip the origin
              if (.not. (ig1 .eq. 0 .and. ig2 .eq. 0 .and. ig3 .eq. 0)) then
                g(:) = ig1*g1(:) + ig2*g2(:) + ig3*g3(:)
                argument = magnitude(g)**2/(4.0d0*kappa**2)

! The variable stuff contains a number of factors, including the exponential
! which is expensive to compute. That is why its outside the iatom,jatom loop.
                stuff = 4.0d0*pi*exp(-argument)/(magnitude(g)**2*s%volume)

! Sum over s and s', the basis indices.
                do iatom = 1, s%natoms

                  ! cut some lengthy notation
                  nullify (pfi)
                  pfi=>s%forces(iatom)

                  do jatom = iatom, s%natoms

                    ! cut some lengthy notation
                    nullify (pfj)
                    pfj=>s%forces(jatom)

                    factor = 1.0d0*stuff
                    factorf = 2.0d0*stuff
                    if (jatom .eq. iatom) factor = 0.5d0*stuff
! g dot b:
                    gdotb =                                                   &
     &                g(1)*(s%atom(iatom)%ratom(1) - s%atom(jatom)%ratom(1))  &
     &              + g(2)*(s%atom(iatom)%ratom(2) - s%atom(jatom)%ratom(2))  &
     &              + g(3)*(s%atom(iatom)%ratom(3) - s%atom(jatom)%ratom(3))

                    s%ewald(iatom,jatom) = s%ewald(iatom,jatom) + factor*cos(gdotb)
                    s%ewald(jatom,iatom) = s%ewald(jatom,iatom) + factor*cos(gdotb)

                    ! forces
                    qq = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
!                   qq = (Q(iatom) - Q0(iatom))*(Q(jatom) - Q0(jatom))
                    pfi%ewald = pfi%ewald + qq*factorf*sin(gdotb)*g
                    pfj%ewald = pfj%ewald - qq*factorf*sin(gdotb)*g

                    s%dewald(:,iatom,jatom) = s%dewald(:,iatom,jatom) - factor*sin(gdotb)*g
                    s%dewald(:,jatom,iatom) = s%dewald(:,jatom,iatom) + factor*sin(gdotb)*g
                  end do
                end do
              end if
            end do
          end do
        end do

! ***********************************************************************
! Compute gamma2:
! ***********************************************************************
! Initialize fewald2

! If we are doing only a cluster, then only the central cell is considered
! in the sum.
        if (s%icluster .eq. 1) then
          il1mx = 0
          il2mx = 0
          il3mx = 0
          kappa = 0.0d0
        end if

! Now carry out the sum over the cells.
        do il1 = -il1mx, il1mx
          do il2 = -il2mx, il2mx
            do il3 = -il3mx, il3mx

! Sum over atoms iatom and atoms jatom. Note that we sum over jatom .ge. iatom
! which yields an extra factor of two for iatom .ne. jatom.
              do iatom = 1, s%natoms

                ! cut some lengthy notation
                nullify (pfi)
                pfi=>s%forces(iatom)

                do jatom = iatom, s%natoms

                  ! cut some lengthy notation
                  nullify (pfj)
                  pfj=>s%forces(jatom)

                  factor = 1.0d0
                  factorf = 2.0d0
                  if (jatom .eq. iatom) factor = 0.5d0

                  cvec = il1*a1vec + il2*a2vec + il3*a3vec
                  cvec = cvec + s%atom(iatom)%ratom - s%atom(jatom)%ratom
                  z = sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2)

! skip the infinite self term.
                  if (z .gt. 0.0001d0) then
                    argument = kappa*z

                    s%ewald(iatom,jatom) =                                    &
     &                s%ewald(iatom,jatom) + factor*erfc(argument)/z
                    s%ewald(jatom,iatom) =                                    &
     &                s%ewald(jatom,iatom) + factor*erfc(argument)/z

                    ! forces
                    qq = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
!                   qq = (Q(iatom) - Q0(iatom))*(Q(jatom) - Q0(jatom))
                    pfi%ewald =                                               &
     &                pfi%ewald + qq*factorf*cvec                             &
     &                           *(2.0d0*exp(-argument**2)*kappa*z/sqrt(pi)   &
     &                             + erfc(argument))/z**3
                    pfj%ewald =                                               &
     &                pfj%ewald - qq*factorf*cvec                             &
     &                           *(2.0d0*exp(-argument**2)*kappa*z/sqrt(pi)   &
     &                             + erfc(argument))/z**3

                    ! forces
                    s%dewald(:,iatom,jatom) = s%dewald(:,iatom,jatom)         &
     &                - factor*cvec*(2.0d0*exp(-argument**2)*kappa*z/sqrt(pi) &
     &                               + erfc(argument))/z**3
                    s%dewald(:,jatom,iatom) = s%dewald(:,jatom,iatom)         &
     &                + factor*cvec*(2.0d0*exp(-argument**2)*kappa*z/sqrt(pi) &
     &                               + erfc(argument))/z**3
                  end if
                end do
              end do
            end do
          end do
        end do

! ***********************************************************************
! Compute gamma3:
! ***********************************************************************
! This term should remain constant always - unless the charges as a function
! of r are changing. Also if the parameter kappa changes corresponding to
! the lattice vectors.
! There are no forces.
        do iatom = 1, s%natoms
          s%ewald(iatom,iatom) = s%ewald(iatom,iatom) - 2.0d0*kappa/sqrt(pi)
        end do

! ***********************************************************************
! Compute gamma4:
! ***********************************************************************
! gamma4 is zero!

! ***********************************************************************
! Combine ewald pieces
! ***********************************************************************

! Deallocate Arrays
! ===========================================================================
        deallocate (Q0, Q)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_ewald


! ===========================================================================
! destroy_Dassemble_ewald
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the assemble_2c_DOGS
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
        subroutine destroy_Dassemble_ewald (s)
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
        deallocate (s%dewald)
        do iatom = 1, s%natoms
          deallocate (s%forces(iatom)%ewaldsr)
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
        end subroutine destroy_Dassemble_ewald

! End Module
! ===========================================================================
        end module M_Dassemble_ewald
