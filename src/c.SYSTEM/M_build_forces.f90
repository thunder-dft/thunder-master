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

! M_build_forces
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the forces for the two-center interactions for
!! the Harris interactions.
!!
!! It contains the following subroutines within the module:
!!
!!      build_forces - sums all dF/dx terms from Dassemblers to get Forces.
!!                     and then finally sums all together to get an Ftot.
!!      writeout_forces - write out the forces components for each atom
!!
! ===========================================================================
        module M_build_forces

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! initialize_forces
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This subroutine initializes the force arrays.
!
! ===========================================================================
! Code written by:
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
        subroutine initialize_forces (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                   !< counter over atoms/neighbors

! Allocate Arrays
! ===========================================================================
! Forces are stored in a Type with each piece, this makes acessing them and use
! pretty easy across the game.
        allocate (s%forces (s%natoms))

! Procedure
! ===========================================================================
! Initialize forces to zero
        do iatom = 1, s%natoms
          ! band-structure interactions
          s%forces(iatom)%febs = 0.0d0
          s%forces(iatom)%kinetic = 0.0d0
          s%forces(iatom)%vna = 0.0d0
          s%forces(iatom)%vxc = 0.0d0
          s%forces(iatom)%vnl = 0.0d0
          s%forces(iatom)%ewald = 0.0d0
          s%forces(iatom)%ewaldlr = 0.0d0

          ! corrections to the force
          s%forces(iatom)%usr = 0.0d0
          s%forces(iatom)%pulay = 0.0d0

          ! three-center interactions
          ! non-local forces
          s%forces(iatom)%f3nla = 0.0d0
          s%forces(iatom)%f3nlb = 0.0d0
          s%forces(iatom)%f3nlc = 0.0d0
          ! Hartree forces
          s%forces(iatom)%f3naa = 0.0d0
          s%forces(iatom)%f3nab = 0.0d0
          s%forces(iatom)%f3nac = 0.0d0
          ! exchange-correlation forces
          s%forces(iatom)%f3xca = 0.0d0
          s%forces(iatom)%f3xcb = 0.0d0
          s%forces(iatom)%f3xcc = 0.0d0

          ! total force
          s%forces(iatom)%ftot  = 0.0d0
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_forces


! ===========================================================================
! build_forces
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This subroutine builds the total forces by adding contributions from the
!! kinetic, Vna, etc. and stores to a T_force variable called forces (:).
!
! ===========================================================================
! Code written by:
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
        subroutine build_forces (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh, matom !< counter over atoms/neighbors
        integer in1, in2             !< species number
        integer jatom, num_neigh     !< counters over neighbors
        integer logfile              !< writing to which unit
        integer mbeta                !< the cell containing neighbor of iatom
        integer norb_mu, norb_nu     !< size of the (mu, nu) block for pair
        integer ix                   !< counter over dimensions
        integer imu, inu             !< counter over MEs

        real maxf, minf, rms         !< for calculating the RMS force
        real sumT

        ! forces
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

        ! kinetic energy
        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic

        ! for overlap repulsive force
        type(T_assemble_block), pointer :: pCape_neighbors
        type(T_assemble_neighbors), pointer :: pcapemat
        type(T_assemble_block), pointer :: poverlap_neighbors
        type(T_assemble_neighbors), pointer :: poverlap

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! ***************************************************************************
!       T W O - C E N T E R   B A N D - S T R U C T U R E   F O R C E S
! ***************************************************************************

! ***************************************************************************
! KINETIC FORCES (TWO-CENTER)
! ***************************************************************************
! loop over atoms in central cell
        do iatom = 1, s%natoms
          matom = s%neigh_self(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          pfi=>s%forces(iatom)

          ! density matrix
          pdenmat=>s%denmat(iatom)

          ! interactions for each contribution
          pkinetic=>s%kinetic(iatom)

! Now loop over all neighbors ineigh of iatom.
          do ineigh = 1, num_neigh
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some lengthy notation
            pfj=>s%forces(jatom)

            ! density matrix - neighbors
            pRho_neighbors=>pdenmat%neighbors(ineigh)

            ! interactions - neighbors
            pK_neighbors=>pkinetic%neighbors(ineigh)

! ***************************************************************************
! The derivatives are tpx and, where p means derivative and x means crytal
! coordinates. The derivative is a vector in crystal
! coordinates and is stored in pK_neighbors%Dblock. The subroutine
! returns the derivative for just that one value of iatom and ineigh, and the
! result is returned in the arguement list, tpx(3,4,4).
            do ix = 1, 3
              sumT = 0.0d0
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  sumT = sumT                                                &
                   + pRho_neighbors%block(imu,inu)*pK_neighbors%Dblock(ix,imu,inu)
                end do
              end do

! Now add sum to appropriate force term. see notes "the total band structure
! The (-1.d0) makes it "force-like".
              ! direct term
              pfi%kinetic(ix) = pfi%kinetic(ix) + (-1.0d0)*sumT
              ! cross term
              pfj%kinetic(ix) = pfj%kinetic(ix) - (-1.0d0)*sumT
            end do ! do ix
          end do ! end loop over neighbors
        end do ! end loop over atoms

! ADD KINETIC CONTRIBUTIONS TO TOTAL FORCE
! ***************************************************************************
! loop over atoms in central cell
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)

! overlap repulsive contribution to total force
! ****************************************************************************
          pfi%febs = pfi%febs + pfi%kinetic
          pfi%ftot = pfi%ftot + pfi%kinetic
        end do
! ***************************************************************************
! END KINETIC FORCES (TWO-CENTER)
! ***************************************************************************

! ***************************************************************************
! ADD OTHER CONTRIBUTIONS TO GET TOTAL BAND-STRUCTURE FORCE (TWO-CENTER)
! ***************************************************************************
! loop over atoms in central cell
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)

! Loop over all neighbors of iatom and add in the neighbor-contributed forces
! ****************************************************************************
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh

! cut some lengthy notation
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            pfj => s%forces(jatom)

! vna contribution to total force
! ****************************************************************************
! Hartree forces - atom case
            pfi%vna = pfi%vna + pfi%vna_atom(:,ineigh)
            pfj%vna = pfj%vna - pfi%vna_atom(:,ineigh)

            pfi%febs = pfi%febs + pfi%vna_atom(:,ineigh)
            pfj%febs = pfj%febs - pfi%vna_atom(:,ineigh)

            pfi%ftot = pfi%ftot + pfi%vna_atom(:,ineigh)
            pfj%ftot = pfj%ftot - pfi%vna_atom(:,ineigh)

! Hartree forces - ontop terms
            pfi%vna = pfi%vna + pfi%vna_ontop(:,ineigh)
            pfj%vna = pfj%vna - pfi%vna_ontop(:,ineigh)

            pfi%febs = pfi%febs + pfi%vna_ontop(:,ineigh)
            pfj%febs = pfj%febs - pfi%vna_ontop(:,ineigh)

            pfi%ftot = pfi%ftot + pfi%vna_ontop(:,ineigh)
            pfj%ftot = pfj%ftot - pfi%vna_ontop(:,ineigh)

! vxc contribution to total force
! ****************************************************************************
! off site interactions
            pfi%vxc = pfi%vxc + pfi%vxc_off_site(:,ineigh)
            pfj%vxc = pfj%vxc - pfi%vxc_off_site(:,ineigh)

            pfi%febs = pfi%febs + pfi%vxc_off_site(:,ineigh)
            pfj%febs = pfj%febs - pfi%vxc_off_site(:,ineigh)

            pfi%ftot = pfi%ftot + pfi%vxc_off_site(:,ineigh)
            pfj%ftot = pfj%ftot - pfi%vxc_off_site(:,ineigh)

! on site interactions
            pfi%vxc = pfi%vxc + pfi%vxc_on_site(:,ineigh)
            pfj%vxc = pfj%vxc - pfi%vxc_on_site(:,ineigh)

            pfi%febs = pfi%febs + pfi%vxc_on_site(:,ineigh)
            pfj%febs = pfj%febs - pfi%vxc_on_site(:,ineigh)

            pfi%ftot = pfi%ftot + pfi%vxc_on_site(:,ineigh)
            pfj%ftot = pfj%ftot - pfi%vxc_on_site(:,ineigh)

! ewald contributions to total ewald force
! ****************************************************************************
! ewaldsr interactions
            pfi%febs = pfi%febs - pfi%ewaldsr(:,ineigh)
            pfj%febs = pfj%febs + pfi%ewaldsr(:,ineigh)

            pfi%ftot = pfi%ftot - pfi%ewaldsr(:,ineigh)
            pfj%ftot = pfj%ftot + pfi%ewaldsr(:,ineigh)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! TOTAL FORCE AFTER EWALD CORRECTION
! ***************************************************************************
! loop over atoms in central cell
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)
          pfi%febs = pfi%febs + pfi%ewaldlr
          pfi%ftot = pfi%ftot + pfi%ewaldlr
        end do

! Vnl contribution to total force
! ****************************************************************************
! Loop over all atoms iatom in the central cell.
        do iatom = 1, s%natoms

! cut some lengthy notation
          pfi=>s%forces(iatom)

! Loop over all neighbors of iatom and add in the neighbor-contributed forces
! Note - the neighbor mapping for vnl is different than neighbor mapping for
! other terms, so, we need to add in the contributions correctly.
          do ineigh = 1, s%neighbors_PP(iatom)%neighn
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)

! cut some lengthy notation
            pfj => s%forces(jatom)

! atom contribution
            pfi%vnl = pfi%vnl + pfi%vnl_atom(:,ineigh)
            pfj%vnl = pfj%vnl - pfi%vnl_atom(:,ineigh)

            pfi%febs = pfi%febs + pfi%vnl_atom(:,ineigh)
            pfj%febs = pfj%febs - pfi%vnl_atom(:,ineigh)

            pfi%ftot = pfi%ftot + pfi%vnl_atom(:,ineigh)
            pfj%ftot = pfj%ftot - pfi%vnl_atom(:,ineigh)
          end do

! ontop left contribution
          do ineigh = 1, s%neighbors_PPx(iatom)%neighn
            jatom = s%neighbors_PPx(iatom)%neigh_j(ineigh)

! cut some lengthy notation
            pfj => s%forces(jatom)

            pfi%vnl = pfi%vnl + pfi%vnl_ontop(:,ineigh)
            pfj%vnl = pfj%vnl - pfi%vnl_ontop(:,ineigh)

            pfi%febs = pfi%febs + pfi%vnl_ontop(:,ineigh)
            pfj%febs = pfj%febs - pfi%vnl_ontop(:,ineigh)

            pfi%ftot = pfi%ftot + pfi%vnl_ontop(:,ineigh)
            pfj%ftot = pfj%ftot - pfi%vnl_ontop(:,ineigh)
          end do ! end loop over neighbors
        end do  ! end loop over atoms

! ***************************************************************************
!                                   E N D
!       T W O - C E N T E R   B A N D - S T R U C T U R E   F O R C E S
! ***************************************************************************

! ***************************************************************************
! ADD OTHER CONTRIBUTIONS TO GET TOTAL BAND-STRUCTURE FORCE (THREE-CENTER)
! ***************************************************************************
! Loop over all atoms iatom in the central cell.
! Single-source loops (not dependent on neighbours)
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)

! vna three-center contribution to the total force
! ****************************************************************************
          pfi%febs = pfi%febs + pfi%f3naa + pfi%f3nab + pfi%f3nac
          pfi%ftot = pfi%ftot + pfi%f3naa + pfi%f3nab + pfi%f3nac

! vxc three-center contribution to the total force
! ****************************************************************************
          pfi%febs = pfi%febs + pfi%f3xca + pfi%f3xcb + pfi%f3xcc
          pfi%ftot = pfi%ftot + pfi%f3xca + pfi%f3xcb + pfi%f3xcc

! vnl three-center contribution to the total force
! ****************************************************************************
          pfi%febs = pfi%febs + pfi%f3nla + pfi%f3nlb + pfi%f3nlc
          pfi%ftot = pfi%ftot + pfi%f3nla + pfi%f3nlb + pfi%f3nlc
        end do ! end loop over atoms
! ***************************************************************************
!                                   E N D
!     T H R E E - C E N T E R   B A N D - S T R U C T U R E   F O R C E S
! ***************************************************************************

! ***************************************************************************
!
!            P U L A Y   C O R R E C T I O N S   (T W O - C E N T E R)
! ***************************************************************************
! loop over atoms in central cell
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          pfi=>s%forces(iatom)

          ! density matrix with eigenvalues
          pcapemat=>s%capemat(iatom)

          ! interactions for each contribution
          poverlap=>s%overlap(iatom)

! Now loop over all neighbors ineigh of iatom.
          do ineigh = 1, num_neigh
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some lengthy notation
            pfj=>s%forces(jatom)

            ! density matrix - neighbors
            pCape_neighbors=>pcapemat%neighbors(ineigh)

            ! interactions - neighbors
            poverlap_neighbors=>poverlap%neighbors(ineigh)

! The derivatives are tpx and, where p means derivative and x means crytal
! coordinates. The derivative is a vector in crystal
! coordinates and is stored in pK_neighbors%Dblock. The subroutine
! returns the derivative for just that one value of iatom and ineigh, and the
! result is returned in the arguement list, tpx(3,4,4).
            do ix = 1, 3
              sumT = 0.0d0
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  sumT = sumT                                                &
                   + pCape_neighbors%block(imu,inu)*poverlap_neighbors%Dblock(ix,imu,inu)
                end do
              end do

! Now add sum to appropriate force term. see notes "the total band structure
! The (-1.d0) makes it "force-like".
              ! direct term
              pfi%pulay(ix) = pfi%pulay(ix) + (-1.0d0)*sumT
              ! cross term
              pfj%pulay(ix) = pfj%pulay(ix) - (-1.0d0)*sumT
            end do ! do ix
          end do ! end loop over neighbors
        end do ! end loop over atoms

! TOTAL FORCE AFTER PULAY CORRECTION
! ***************************************************************************
! loop over atoms in central cell
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)

! overlap repulsive contribution to total force
! ****************************************************************************
          pfi%ftot = pfi%ftot - pfi%pulay
        end do
! ***************************************************************************
!                                  E N D
!            P U L A Y   C O R R E C T I O N S   (T W O - C E N T E R)
! ***************************************************************************

! ***************************************************************************
!
!            U S R   C O R R E C T I O N S   (T W O - C E N T E R)
! ***************************************************************************
! Loop over all atoms iatom in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)
          pfi%ftot = pfi%ftot + pfi%usr
        end do ! end loop over atoms

! Calculate the RMS force and MAX force
        rms = 0.0d0
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)
          rms = rms + pfi%ftot(1)**2 + pfi%ftot(2)**2 + pfi%ftot(3)**2
          maxf = max(maxf, maxval(pfi%ftot))
          minf = abs(minval(pfi%ftot))
          if (minf .gt. maxf) maxf = minf
        end do
        rms = sqrt(rms/(3*s%natoms))
        write (logfile,*)
        write (logfile,400) maxf, rms
        write (logfile,*)

! Format Statements
! ===========================================================================
400     format (2x, ' Cartesian Forces:  Max = ', f16.8, '    RMS = ', f16.8)

! End Subroutine
! ===========================================================================
        return
        end subroutine build_forces


! ===========================================================================
! writeout_forces
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine is a utility to write out the components of the forces.
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
        subroutine writeout_forces (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter for atom loop
        integer logfile                     !< writing to which unit

        type(T_forces), pointer :: pfi

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*)
        write (logfile,103) 'The kinetic force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_kinetic', iatom, s%atom(iatom)%species%symbol,&
     &                                            s%forces(iatom)%kinetic
        end do
        write (logfile,100)

        write (logfile,*)
        write (logfile,103) 'The Hartree (vna) two-center force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_vna_2c', iatom, s%atom(iatom)%species%symbol,  &
     &                                           s%forces(iatom)%vna
        end do
        write (logfile,100)

        write (logfile,*)
        write (logfile,103) 'The Hartree (vna) three-center force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)
          write (logfile,102) 'f_vna_3c', iatom, s%atom(iatom)%species%symbol,  &
     &      s%forces(iatom)%f3naa + s%forces(iatom)%f3nab + s%forces(iatom)%f3nac
        end do
        write (logfile,100)

        write (logfile,*)
        write (logfile,103) 'The non-local pseudopotential (vnl) two-center force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_vnl_2c', iatom, s%atom(iatom)%species%symbol,  &
     &                                           s%forces(iatom)%vnl
        end do
        write (logfile,100)

        write (logfile,*)
        write (logfile,103) 'The non-local pseudopotential (vnl) three-center force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_vnl_3c', iatom, s%atom(iatom)%species%symbol,  &
     &      s%forces(iatom)%f3nla + s%forces(iatom)%f3nlb + s%forces(iatom)%f3nlc
        end do
        write (logfile,100)

        write (logfile,100)
        write (logfile,*)
        write (logfile,103) 'The exchange correlation (vxc) two-center force: '
        write (logfile,100)  
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_vxc_2c', iatom, s%atom(iatom)%species%symbol,  &
     &                                           s%forces(iatom)%vxc
        end do
        write (logfile,100)

        write (logfile,100)
        write (logfile,*)
        write (logfile,103) 'The exchange correlation (vxc) three-center force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)
          write (logfile,102) 'f_vxc_3c', iatom, s%atom(iatom)%species%symbol,  &
     &       s%forces(iatom)%f3xca + s%forces(iatom)%f3xcb + s%forces(iatom)%f3xcc
        end do
        write (logfile,100)

        write (logfile,100)
        write (logfile,*)
        write (logfile,103) 'The long-range electrostatics (ewald) force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_ewald', iatom, s%atom(iatom)%species%symbol,  &
     &                                          s%forces(iatom)%ewaldlr
        end do
        write (logfile,100)

        write (logfile,*)
        write (logfile,103) 'The total band-structure force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_ebs', iatom, s%atom(iatom)%species%symbol,    &
     &                                        s%forces(iatom)%febs
        end do
        write (logfile,100)

        write (logfile,*)
        write (logfile,103) 'The overlap repulsive force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_pulay', iatom, s%atom(iatom)%species%symbol,&
     &                                          s%forces(iatom)%pulay
        end do
        write (logfile,100)

        write (logfile,*)
        write (logfile,103) 'The short-range (double-counting) (usr) force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_usr', iatom, s%atom(iatom)%species%symbol, &
     &                                        s%forces(iatom)%usr
        end do

        write (logfile,*)
        write (logfile,103) 'The total force: '
        write (logfile,100)
        write (logfile,101)
        write (logfile,100)
        do iatom = 1, s%natoms
          write (logfile,102) 'f_total', iatom, s%atom(iatom)%species%symbol, &
     &                                          s%forces(iatom)%ftot
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (4x, 70('='))
101     format (4x, 'Force ', 'Atom # ', 2x, ' Type ', 5x,   &
     &              ' x ', 9x, ' y ', 9x, ' z ')
102     format (4x, A,  i5, 7x, a2, 3(2x,ES10.3))
103     format (4x, A)


! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_forces


! End Module
! ===========================================================================
        end module M_build_forces
