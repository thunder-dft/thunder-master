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

! M_non_adiabatic
! Module Description
! ===========================================================================
!>       This module handles all of the routines and global variables for
! the non-adiabatic molecular-dynamics formalism.
!
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        module M_non_adiabatic
        use M_assemble_blocks
        use M_configuraciones

! Type Declaration
! ===========================================================================
! Each atom in the system has its own type.
        type T_atom_NAC
          real, dimension (3) :: ratom_old ! atom previous positions
          real, dimension (3) :: dratom    ! electronic interpolation of atom positions
          real, dimension (3) :: vatom_old ! atom previous velocities
          real, dimension (3) :: dvatom    ! electronic interpolation of atom velocities
        end type T_atom_NAC

        type(T_atom_NAC), pointer :: atom_NAC (:)   ! atom's information

! Define gradH which is based on the "forces" type
        type(T_forces), pointer :: gradH (:)

! Density for a mu by nu block - special case - only density matrix for
! states involved in the transitions
        type(T_assemble_neighbors), pointer :: denmat_MDET (:)
        real, allocatable :: djk(:,:,:,:)    !<Phij|delPhik>
        real, allocatable :: gks_atom(:,:,:)    !each atoms contribution to djk dot velocity
        real, allocatable :: gks(:,:)            ! djk dot velocity for each transition

! module procedures
        contains


! ===========================================================================
! calculate_nac
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This nacouplings subroutinge will (for now) generate random numbers for the
! non adiabatic coupling constants. This is used as a fill in until forces are
! added and more realistic NAC constants can be added.
!
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        subroutine calculate_nac (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Local Variable Declaration and Description
! ===========================================================================
!        integer itransition
!        integer iatom
!        integer ix
!        integer iorbital
!        integer jorbital
!        integer jtransition

!        real :: xrand

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
        call build_gradH (s)

! End Subroutine
! ===========================================================================
        return
        end subroutine calculate_nac


! ===========================================================================
! build_gradH
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This subroutine builds the gradient of the Hamiltonian by adding
! contributions from the kinetic, Vna, etc. and stores to a T_gradH variable
! called gradH (:). This will be used for the construction of the NAC constants
!
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        subroutine build_gradH (s)
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
        integer iatom, ineigh, matom !, ialpha    !< counter over atoms and neighbors
        integer in1, in2 !, in3           !< species number
        integer jatom, num_neigh        !< counters over neighbors
        integer norb_mu, norb_nu        !< size of the (mu, nu) block for pair
!        integer ix                      !< counter over dimensions
        integer imu, inu                !< counter over MEs
        integer mmu, nnu
        integer mbeta                   !< the cell containing iatom's neighbors
        integer itransition             !< counter of transitions
        integer jtransition
        integer iband                   !< counter of bands
        integer jband
        real, dimension (3) :: r1, r2 !, r3 !< position of atoms
!        real, dimension (3) :: gHtot !analogous to ftot in build forces
        real eigeni
        real eigenj
        real deltaE

!        integer mneigh
!        integer ibeta, jbeta

!        real sumT

!       type(T_assemble_block), pointer :: pH_neighbors
!       type(T_assemble_neighbors), pointer :: pHamiltonian
        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic
        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna
!       type(T_assemble_block), pointer :: pSR_neighbors
!       type(T_assemble_neighbors), pointer :: pewaldsr
!       type(T_assemble_block), pointer :: pLR_neighbors
!       type(T_assemble_neighbors), pointer :: pewaldlr
!       type(T_assemble_block), pointer :: pvxc_neighbors
!       type(T_assemble_neighbors), pointer :: pvxc

        ! Density matrix stuff
!        type(T_assemble_neighbors), pointer :: pdenmat
!        type(T_assemble_block), pointer :: pRho_neighbors
!        type(T_assemble_block), pointer :: pRho_neighbors_matom

        type(T_forces), pointer :: pgradH
        type(T_forces), pointer :: pgradH_iatom
        type(T_forces), pointer :: pgradH_jatom
!       type(T_forces), pointer :: gHtot

! Allocate Arrays & Initialize things
! ===========================================================================
! Forces are stored in a Type with each piece, this makes acessing them and use
! pretty easy across the game.
        allocate (gradH (s%natoms))
        allocate (djk(3,s%natoms,ntransitions,ntransitions))
        allocate (gks_atom(s%natoms,ntransitions,ntransitions))
        allocate (gks(ntransitions,ntransitions))
        djk = 0.0d0
        gks_atom = 0.0d0
        gks = 0.0d0

! Procedure
! ===========================================================================
! Loop over atoms in central cell
        do itransition = 1, ntransitions
        do jtransition = 1, ntransitions
           iband = s%kpoints(1)%transition(itransition,1)%imap
           jband = s%kpoints(1)%transition(jtransition,1)%imap
           eigeni = s%kpoints(1)%eigen(iband)
           eigenj = s%kpoints(1)%eigen(jband)
!          if (iband .ne. jband) then
!           deltaE = 1/(eigeni-eigenj)
           do iatom = 1, s%natoms
!             write (*,*) ' in loop over atoms'
             ! cut some lengthy notation
             pkinetic=>s%kinetic(iatom)
             pvna=>s%vna(iatom)
!            pewaldsr=>ewaldsr(iatom)
!            pewaldlr=>ewaldlr(iatom)
!            pvxc=>vxc(iatom)
             pgradh=>gradH(iatom)
!            pdenmat=>denmat_MDET(iatom)

             matom = s%neigh_self(iatom)
             in1 = s%atom(iatom)%imass
             norb_mu = species(in1)%norb_max
             num_neigh = s%neighbors(iatom)%neighn
!            allocate (pdenmat%neighbors(num_neigh))
             r1 = s%atom(iatom)%ratom

! allocate force terms and initialize to zero
             if (itransition .eq. 1 .and. jtransition .eq. 1) then
                allocate (pgradH%vna_atom (3, num_neigh)); pgradH%vna_atom = 0.0d0
                allocate (pgradH%vna_ontop (3, num_neigh)); pgradH%vna_ontop = 0.0d0
             end if

! Now loop over all neighbors ineigh of iatom.
             do ineigh = 1, num_neigh
                ! cut some lengthy notation
!               pH_neighbors=>pHamiltonian%neighbors(ineigh)
                pK_neighbors=>pkinetic%neighbors(ineigh)
                pvna_neighbors=>pvna%neighbors(ineigh)
!               pSR_neighbors=>pewaldsr%neighbors(ineigh)
!               pLR_neighbors=>pewaldlr%neighbors(ineigh)
!               pvxc_neighbors=>pvxc%neighbors(ineigh)
!               pRho_neighbors=>pdenmat%neighbors(ineigh)
!               pRho_neighbors_matom=>pdenmat%neighbors(matom)

               mbeta = s%neighbors(iatom)%neigh_b(ineigh)
               jatom = s%neighbors(iatom)%neigh_j(ineigh)
               in2 = s%atom(jatom)%imass
               norb_nu = species(in2)%norb_max
               r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a

! ****************************************************************************
!
!  ASSEMBLE T FORCES
! ****************************************************************************
! The derivatives are tpx and, where p means derivative and x means crytal
! coordinates. The derivative is a vector in crystal
! coordinates and is stored in pK_neighbors%Dblock. The subroutine
! returns the derivative for just that one value of iatom and ineigh, and the
! result is returned in the arguement list, tpx(3,4,4).
!           do ix = 1, 3
!             sumT = 0.0d0
!             do inu = 1, norb_nu
!               do imu = 1, norb_mu
!                 sumT = sumT + pRho_neighbors%block(imu,inu)*pK_neighbors%Dblock(ix,imu,inu)
!               end do
!             end do

! Now add sum to appropriate force term. see notes "the total band structure
! force", ofs 9/14/88. Also see notes "total force due to fr offsites"
!                                 and "total force due to fs offsites" 11/10/88
! The (-1.d0) makes it "force-like".

! Direct terms.
!             forces(iatom)%kinetic(ix) = forces(iatom)%kinetic(ix) + (-1.0d0)*sumT

! Cross terms.
!             forces(jatom)%kinetic(ix) = forces(jatom)%kinetic(ix) - (-1.0d0)*sumT
!           end do ! do ix

! ****************************************************************************
!
!  ASSEMBLE Vna ONTOP FORCES
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
               if (iatom .eq. jatom) then

! Do nothing here - special case. Interaction already calculated in atm case.

               else

! Notice the explicit negative sign, this makes it force like.

              do inu = 1, norb_nu
                nnu = inu + s%iblock_slot(iatom)
!                write (*,*) (s%kpoints(1)%c(nnu,icurrent_state))
                do imu = 1, norb_mu
                   mmu = imu + s%iblock_slot(jatom)
 !                  write (*,*) 'block_slot(jatom)', s%iblock_slot(jatom)
!                   write (*,*) s%kpoints(1)%c(mmu,iband)
                   pgradH%vna_ontop(:,ineigh) = pgradH%vna_ontop(:,ineigh)    &
     &              - real(conjg(s%kpoints(1)%c(nnu,iband))  &
     &                 *s%kpoints(1)%c(mmu,jband))*pvna_neighbors%Dblocko(:,imu,inu)*P_eq2
                end do
              end do
!            write (*,*) 'trying out ineigh stuff', pgradH%vna_ontop(:,1)
            end if


! ****************************************************************************
!
! ASSEMBLE NEUTRAL ATOM FORCE FOR ATOM CASE
! ****************************************************************************
! The vna 2 centers are: ontop (L), ontop (R), and atm.
! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
!
! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
! The derivative wrt the "atom r1" position (not the NA position) are
! stored in bcnapx.
!
! Form the "force-like" derivative of the atom terms for NA,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.

! Note that the loop below involves num_orb(in1) ONLY. Why?
! Because the potential is somewhere else (or even at iatom), but we are
! computing the vna_atom term, i.e. < phi(i) | v | phi(i) > but V=v(j) )
! interactions.

! Notice the explicit negative sign, this makes it force like.
               do inu = 1, norb_nu
                  nnu = inu + s%iblock_slot(iatom)
                  do imu = 1, norb_mu
                     mmu = imu + s%iblock_slot(jatom)
                     pgradH%vna_atom(:,ineigh) = pgradH%vna_atom(:,ineigh) -      &
     &               real(conjg(s%kpoints(1)%c(nnu,iband))     &
     &                 *s%kpoints(1)%c(mmu,jband))*pvna_neighbors%Dblock(:,imu,inu)*P_eq2
                  end do
              end do
            end do ! end loop over neighbors
         end do ! end loop over atoms
 !        end if
       end do ! end loop over itransitions
       end do  ! end loop over jtransitions
! ****************************************************************************
! Finally, sum all gradh terms into Total GradH
! ****************************************************************************
! Loop over all atoms iatom in the central cell.
      do itransition = 1, ntransitions
      do jtransition = 1, ntransitions
        iband = s%kpoints(1)%transition(itransition,1)%imap
        jband = s%kpoints(1)%transition(jtransition,1)%imap
        if (iband .ne. jband) then
        do iatom = 1, s%natoms
          num_neigh = s%neighbors(iatom)%neighn

! cut some lengthy notation
          pgradH_iatom=>gradH(iatom)
!          gHtot=>ftot

! Single-source contributions (not dependent on neighbours

! ****************************************************************************
! Kinetic contribution to Total Force
! ****************************************************************************
!          pforce_iatom%ftot = pforce_iatom%ftot + pforce_iatom%kinetic

! ****************************************************************************
! Loop over all neighbors of iatom and add in the neighbor-contributed forces
! ****************************************************************************
           do ineigh = 1, num_neigh

! cut some lengthy notation
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            pgradH_jatom => gradH(jatom)

! ****************************************************************************
! Vna contribution to Total Force
! ****************************************************************************
! There may be a confusion here about the signs.
! The variable "sum" is already force-like because it came from forces%vna_ontop, etc.
! another +/- comes from direct (+) and cross (-) terms.
! See p. 5, "the total band structure force".
! The atom terms.

! neutral atom forces - atm case
            pgradH_iatom%vna = pgradH_iatom%vna + pgradH_iatom%vna_atom(:,ineigh)
            pgradH_jatom%vna = pgradH_jatom%vna - pgradH_iatom%vna_atom(:,ineigh)

            pgradH_iatom%ftot =                                              &
     &        pgradH_iatom%ftot + pgradH_iatom%vna_atom(:,ineigh)
            pgradH_jatom%ftot =                                              &
     &        pgradH_jatom%ftot - pgradH_iatom%vna_atom(:,ineigh)

! The ontop terms.
! The factor 2.0d0 comes from ontop the bra or ontop the ket.
! neutral atom forces - ontop case
            pgradH_iatom%vna =                                               &
     &        pgradH_iatom%vna + 2.0d0*pgradH_iatom%vna_ontop(:,ineigh)
            pgradH_jatom%vna =                                               &
     &        pgradH_jatom%vna - 2.0d0*pgradH_iatom%vna_ontop(:,ineigh)

            pgradH_iatom%ftot =                                              &
     &        pgradH_iatom%ftot + 2.0d0*pgradH_iatom%vna_ontop(:,ineigh)
            pgradH_jatom%ftot =                                              &
     &        pgradH_jatom%ftot - 2.0d0*pgradH_iatom%vna_ontop(:,ineigh)
          end do ! end loop over neighbors
          djk(:,iatom,itransition,jtransition) = pgradH_jatom%ftot
!          write (*,*) 'iatom= ', iatom
!          write (*,*) 'iband and jband', iband, jband
 !         write (*,*) 'djk=', djk(:,iatom,itransition,jtransition)
 !         write (*,*) 'v=', s%atom(iatom)%vatom
          gks_atom(iatom,itransition,jtransition) = dot_product(djk(:,iatom,itransition,jtransition), s%atom(iatom)%vatom)
!          write (*,*) 'gks_atom(iatom,iband,jband)', gks_atom(iatom,itransition,jtransition)
        end do ! end loop over atoms
        end if
        end do ! end of loop over itransitions
        end do ! end of loop over jtranstions
        do itransition = 1, ntransitions
           do jtransition = 1, ntransitions
              iband = s%kpoints(1)%transition(itransition,1)%imap
              jband = s%kpoints(1)%transition(jtransition,1)%imap
              eigeni = s%kpoints(1)%eigen(iband)
              eigenj = s%kpoints(1)%eigen(jband)
              if (iband .ne. jband) then
              deltaE = 1/(eigeni-eigenj)
              do iatom = 1, s%natoms
                  gks (itransition,jtransition) = gks(itransition,jtransition) + deltaE*gks_atom(iatom,itransition,jtransition)  !WARNING NO longer tied to map iband jband *!!!!!
              end do
              end if
!              write(*,*) 'gks=', gks(itransition,jtransition)
           end do
        end do
!write (*,*) gks(1,1)

! Format Statements
! ===========================================================================
! None

! Deallocate Arrays
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine build_gradH


! ===========================================================================
! Subroutine Declaration
! ===========================================================================
!       This routine integrates, with the verlet method, the TD equations for
!       the TDSE coefficients (cna). The TD WF coefficients (c_wf), the NAC
!       constants (dnac), eigenvalues (eigen), velocity (vatom), position(ratom)
!       are interpolated in time

        subroutine evolve_ks_states (s,itime_step,icurrent_state)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.
        integer, intent(in) :: itime_step
        integer, intent(inout) :: icurrent_state
! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer iatom
        integer itransition
        integer jtransition
        integer ntransitions

        integer iorbital
        integer jorbital
        integer ix
        real dt
        real ddt
        real nddt
        real delta

! Allocate Arrays
! ===========================================================================
        do ikpoint = 1, s%nkpoints
!          ntransitions = s%kpoints(ikpoint)%ntransitions
          do itransition = 1, ntransitions
               allocate (s%kpoints(ikpoint)%transition(itransition,1)%ddjk(3))
          end do
        end do
!
! Initialize
! ===========================================================================

!       Initialize c_wf to complex, pointer :: c (:, :) (in M_kpoints)
!       added ratom_old, vatom_old to M_configuraciones, but that may not have been a good idea
!       added eigen_old to M_kpoints

        if (itime_step .eq. 1) then
         write (*,*) 'time step is', itime_step
         do iatom = 1, s%natoms
!		    s%atom(iatom)%ratom_old = s%atom(iatom)%ratom + 0.5
            atom_NAC(iatom)%vatom_old = s%atom(iatom)%vatom
!           write (*,*) s%atom(iatom)%ratom, s%atom(iatom)%ratom_old
            write (*,*) s%atom(iatom)%vatom, atom_NAC(iatom)%vatom_old
         end do

         do ikpoint = 1, s%nkpoints
!           ntransitions = s%kpoints(ikpoint)%ntransitions ! don't think I need anymore
            do itransition = 1, ntransitions
               s%kpoints(ikpoint)%transition(itransition,1)%cna_old =          &
     &                 s%kpoints(ikpoint)%transition(itransition,1)%cna
            end do
         end do

         do ikpoint = 1, s%nkpoints
            s%kpoints(ikpoint)%eigen_old = s%kpoints(ikpoint)%eigen
!            write (*,*) 'eigenvalues', s%kpoints(ikpoint)%eigen_old
         end do

        end if
! Procedure
! ===========================================================================
! ===========================================================================
! Calculate d/dt c_{ak} at different time steps in between t and t+dt
! tt(it) = t + dt/Nsteps * it . We need to interpolate the values for
! eigen_k, vatom, gks, using their values at t and t+dt. We start using
! a simple linear interpolation.
! ===========================================================================
! Interpolation stuff
        nddt = 1000
       ddt = dt / nddt

        do ikpoint = 1, s%nkpoints
          s%kpoints(ikpoint)%deigen = s%kpoints(ikpoint)%eigen - s%kpoints(ikpoint)%eigen_old
        end do

        do iatom = 1, s%natoms
           atom_NAC(iatom)%dvatom = s%atom(iatom)%vatom - atom_NAC(iatom)%vatom_old
           atom_NAC(iatom)%dratom = s%atom(iatom)%ratom - atom_NAC(iatom)%ratom_old
        end do

        do ikpoint = 1, s%nkpoints
!           ntransitions = s%kpoints(ikpoint)%ntransitions
!          do itransition = 1, ntransitions
              do iorbital = 1, s%norbitals
              if (iorbital .eq. icurrent_state) then
                 do jtransition = 1, ntransitions
                 do jorbital = 1, s%norbitals
                 if (jorbital .eq. s%kpoints(ikpoint)%transition(jtransition,1)%imap) then
                    do ix = 1, 3
                       s%kpoints(ikpoint)%transition(jtransition,1)%ddjk(ix)  =  &
                            & s%kpoints(ikpoint)%transition(jtransition,1)%djk(ix) - &
                            & s%kpoints(ikpoint)%transition(jtransition,1)%djk_old(ix)
                       write (*,*) s%kpoints(ikpoint)%transition(jtransition,1)%ddjk(ix)
                    end do
                 end if
                 end do
                 end do
              end if
              end do
!           end do
        end do

!         do ikpoint = 1, s%nkpoints
!            ntransitions = s%kpoints(ikpoint)%ntransitions ! don't think I need anymore
!            do itransition = 1, ntransitions
!               s%kpoints(ikpoint)%transition(itransition,1)%cna_old =          &
!     &                 s%kpoints(ikpoint)%transition(itransition,1)%cna
!            end do
!         end do
        delta = 0.0d0
!       call dcdt_nac (s, dc_na)
! End Subroutine
! ===========================================================================
        return
        end subroutine evolve_ks_states
! ===========================================================================
! End Module
! ===========================================================================
        end module M_non_adiabatic
