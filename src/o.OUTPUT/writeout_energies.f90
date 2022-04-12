! copyright info:
!
!                             @Copyright 2013
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
! Ohio University - Dave Drabold
! University of Texas at Austin - Alex Demkov
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

! writeout_energies.f90
! Subroutine Description
! ===========================================================================
!       This routine prints out the energy components.
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine writeout_energies (t, ebs, uii_uee, uxcdcc)
        use M_assemble_blocks
        use M_species
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: t             ! the structure to be used

        real, intent (in) :: ebs                   ! band-structure energy
        real, intent (in) :: uii_uee, uxcdcc       ! short-range energies

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh              ! counter over atoms and neighbors
        integer in1, in2                   ! species numbers
        integer imu, inu                   ! counters for mu, nu
        integer inpfile                    ! reading from which unit
        integer jatom                      ! neighbor of iatom
        integer jneigh                     ! neighbor counter for vnl
        integer katom, kbeta               ! neighbor atom and cell number
        integer mbeta                      ! the cell containing neighbor of iatom
        integer num_neigh                  ! number of neighbors
        integer num_neighPPp               ! number of neighbors for vnl elements
        integer norb_mu, norb_nu           ! block size for the ME blocks

! Energies
        real atomic_energy                     ! total atomic energy
        real etot                              ! total energy
        real etot_per_atom                     ! total energy per atom

! Energy Components
        real comp_EBS
        real comp_KE
        real comp_VNA
        real comp_VXC
        real comp_VNL
        real comp_EWDSR
        real comp_EWDLR

        character (len = 25) :: slogfile

        logical header

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat

        type(T_assemble_block), pointer :: pH_neighbors
        type(T_assemble_neighbors), pointer :: pHamiltonian
        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic
        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna
        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc
        type(T_assemble_block), pointer :: pvnl_neighbors
        type(T_assemble_neighbors), pointer :: pvnl
        type(T_assemble_block), pointer :: pSR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldsr
        type(T_assemble_block), pointer :: pLR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldlr

! Procedure
! ===========================================================================
! Initialize logfile
        inpfile = t%inpfile

! Writing out the energy pieces
        write (t%logfile, *)
        write (t%logfile, '(A)') 'Total Energy'
        write (t%logfile, '(A)') '------------'
        write (t%logfile, *)
        write (t%logfile, 102) ebs
        write (t%logfile, 103) uii_uee
        write (t%logfile, 105) uxcdcc

        ! Evaluate total energy
        etot = ebs + uii_uee + uxcdcc
        write (t%logfile, 107) etot

        ! Total energy per atom
        etot_per_atom = etot/s%natoms
        write (t%logfile, 108) etot_per_atom

        ! Cohesive Energy
        atomic_energy = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          atomic_energy = atomic_energy + species(in1)%atomicE
        end do
        write (s%logfile, 109) atomic_energy
        write (s%logfile, 110) etot - atomic_energy
        write (s%logfile, *)
        write (s%logfile, 111) (etot - atomic_energy)/s%natoms
        write (s%logfile, *)

! Writeout the components of the energies if requested.
        if (iwriteout_energies .eq. 1) then

! We use comp to mean component. So comp_KE means kinetic energy component
! to the total energy.  First initialize everything to zero.
          comp_EBS = 0.0d0
          comp_KE = 0.0d0
          comp_VNA = 0.0d0
          comp_VXC = 0.0d0
          comp_VNL = 0.0d0
          comp_EWDSR = 0.0d0
          comp_EWDLR = 0.0d0

! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
! Loop over the atoms in the central cell.
          do iatom = 1, t%natoms
            ! cut some lengthy notation
            pdenmat=>t%denmat(iatom)
            pHamiltonian=>t%Hamiltonian(iatom)
            pkinetic=>t%kinetic(iatom)
            pvna=>t%vna(iatom)
            pvxc=>t%vxc(iatom)
            pvnl=>s%vnl(iatom)
            pewaldsr=>s%ewaldsr(iatom)
            pewaldlr=>s%ewaldlr(iatom)

            in1 = t%atom(iatom)%imass
            norb_mu = species(in1)%norb_max
            num_neigh = t%neighbors(iatom)%neighn

! Loop over the neighbors of each iatom.
            do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
              ! cut some more lengthy notation
              pRho_neighbors=>pdenmat%neighbors(ineigh)
              pH_neighbors=>pHamiltonian%neighbors(ineigh)
              pK_neighbors=>pkinetic%neighbors(ineigh)
              pvna_neighbors=>pvna%neighbors(ineigh)
              pvxc_neighbors=>pvxc%neighbors(ineigh)
              pSR_neighbors=>pewaldsr%neighbors(ineigh)
              pLR_neighbors=>pewaldlr%neighbors(ineigh)

              jatom = t%neighbors(iatom)%neigh_j(ineigh)
              mbeta = t%neighbors(iatom)%neigh_b(ineigh)
              in2 = t%atom(jatom)%imass

! Allocate the block size
              norb_nu = species(in2)%norb_max
              do imu = 1, norb_mu
                do inu = 1, norb_nu

                  ! band structure energy
                  comp_EBS = comp_EBS + pH_neighbors%block(imu,inu)*pRho_neighbors%block(imu,inu)

                  ! kinetic energy
                  comp_KE = comp_KE + pK_neighbors%block(imu,inu)*pRho_neighbors%blocko(imu,inu)

                  ! Hartree energy
                  comp_VNA = comp_VNA + pvna_neighbors%block(imu,inu)*pRho_neighbors%block(imu,inu)

                  ! exchange-correlation energy
                  comp_VXC = comp_VXC + pvxc_neighbors%block(imu,inu)*pRho_neighbors%block(imu,inu)

                  ! short-range Ewald energy
                  comp_EWDSR = comp_EWDSR + pSR_neighbors%block(imu,inu)*pRho_neighbors%block(imu,inu)

                  ! long-range Ewald energy
                  comp_EWDLR = comp_EWDLR + pLR_neighbors%block(imu,inu)*pRho_neighbors%block(imu,inu)
                end do
              end do

! The vnl terms have a different neighbor mapping.  So, for the
! pseudo-potential part - we need to do a comparison - find the piece
! of the vnl matrix element that belongs to iatom, ineigh.
              num_neighPPp = t%neighbors_PPp(iatom)%neighn
              do jneigh = 1, num_neighPPp
                katom = t%neighbors_PPp(iatom)%neigh_j(jneigh)
                kbeta = t%neighbors_PPp(iatom)%neigh_b(jneigh)
                if (katom .eq. jatom .and. kbeta .eq. mbeta) then
                  ! cut some more lengthy notation
                  pvnl_neighbors=>pvnl%neighbors(jneigh)

                  norb_mu = species(in1)%norb_max
                  norb_nu = species(in2)%norb_max
                  do imu = 1, norb_mu
                    do inu = 1, norb_nu
                      ! non-local pseudopotential energy
                      comp_VNL = comp_VNL + pvnl_neighbors%block(imu,inu)*pRho_neighbors%block(imu,inu)
                      comp_EBS = comp_EBS + pvnl_neighbors%block(imu,inu)*pRho_neighbors%block(imu,inu)
                    end do
                  end do
                end if
              end do

! Finish loop over atoms and neighbors.
            end do
          end do

! Writout the energies to a .ENERGIES file
          slogfile = t%basisfile(:len(trim(t%basisfile))-4)
          slogfile = trim(slogfile)//'.ENERGIES'
          inquire (file = slogfile, exist = header)
          open (unit = inpfile, file = slogfile, status = 'unknown', position = 'append')
          ! write out the heading of components if the file did not exist before
          if (.not. header) write (inpfile,201)
          write (inpfile,202) comp_EBS, comp_KE, comp_VNA, comp_VXC,         &
     &                        comp_VNL, comp_EWDSR, comp_EWDLR, uii_uee, uxcdcc
          close (unit = inpfile)
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
102     format (2x, '           ebs = ', f15.6)
103     format (2x, '     uii - uee = ', f15.6)
105     format (2x, '        uxcdcc = ', f15.6)
107     format (2x, '          ETOT = ', f15.6)
108     format (2x, '     Etot/atom = ', f15.6)
109     format (2x, ' Atomic Energy = ', f15.6)
110     format (2x, '     CohesiveE = ', f15.6)
111     format (2x, ' Cohesive Energy per atom  = ', f15.6)

201     format (3x, ' comp_EBS ', 5x, ' comp_KE ', 6x, ' comp_VNA ', 5x,     &
     &              ' comp_VXC ', 5x, ' comp_VNL ', 5x, ' comp_EWDSR ', 3x,  &
     &              ' comp_EWDLR ', 4x, ' uii_uee ', 6x, ' uxcdcc ')
202     format (9f15.6)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_energies
