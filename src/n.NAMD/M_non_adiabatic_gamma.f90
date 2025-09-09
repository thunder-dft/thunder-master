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

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones

! Type Declaration
! ===========================================================================
! Each atom in the system has its own type.




! module procedures
        contains

! ===========================================================================
! build_gHmatrix
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This subroutine builds the Hamiltonian by adding matrices of kinetic, Vna,
! ... etc. involved in the nonadiabatic dynamics and stores to a T_neighbors
! type variable called gHmatrix(:).
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine build_gHmatrix (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer in1, in2                !< species number
        integer jatom, num_neigh        !< counters over neighbors
        integer norb_mu, norb_nu        !< size of the (mu, nu) block for pair

        type(T_assemble_block), pointer :: H_neighbors
        type(T_assemble_neighbors), pointer :: gHamiltonian
        type(T_assemble_block), pointer :: K_neighbors
        type(T_assemble_neighbors), pointer :: gkinetic
        type(T_assemble_block), pointer :: vna_neighbors
        type(T_assemble_neighbors), pointer :: gvna
        type(T_assemble_block), pointer :: SR_neighbors
        type(T_assemble_neighbors), pointer :: gewaldsr
        type(T_assemble_block), pointer :: LR_neighbors
        type(T_assemble_neighbors), pointer :: gewaldlr
        type(T_assemble_block), pointer :: vxc_neighbors
        type(T_assemble_neighbors), pointer :: gvxc

! Allocate Arrays
! ===========================================================================
        allocate (s%gHamiltonian(s%natoms))

! Procedure
! ===========================================================================
! loop over atoms in central cell
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          nullify (gHamiltonian)
          gHamiltonian=>s%gHamiltonian(iatom)

          nullify (gkinetic, gvna, gewaldsr, gewaldlr, gvxc)
          gkinetic=>s%gkinetic(iatom)
          gvna=>s%gvna(iatom)
          gewaldsr=>s%gewaldsr(iatom)
          gewaldlr=>s%gewaldlr(iatom)
          gvxc=>s%gvxc(iatom)

          in1 = s%atom(iatom)%imass
! ===========================================================================
! ===========================================================================
          norb_mu = species(in1)%norb_max
! ===========================================================================
! ===========================================================================
          num_neigh = s%neighbors(iatom)%neighn
          allocate(s%gHamiltonian(iatom)%neighbors(num_neigh))

! Now loop over all neighbors ineigh of iatom.
          do ineigh = 1, num_neigh
            ! cut some lengthy notation
            nullify (H_neighbors)
            H_neighbors=>gHamiltonian%neighbors(ineigh)

            nullify (K_neighbors, vna_neighbors)
            nullify (SR_neighbors, LR_neighbors, vxc_neighbors)
            K_neighbors=>gkinetic%neighbors(ineigh)
            vna_neighbors=>gvna%neighbors(ineigh)
            SR_neighbors=>gewaldsr%neighbors(ineigh)
            LR_neighbors=>gewaldlr%neighbors(ineigh)
            vxc_neighbors=>gvxc%neighbors(ineigh)

            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
! ===========================================================================
! ===========================================================================
            norb_nu = species(in2)%norb_max
! ===========================================================================
! ===========================================================================
            allocate(s%gHamiltonian(iatom)%neighbors(ineigh)%Dblock(3, norb_mu, norb_nu))
            H_neighbors%Dblock = K_neighbors%Dblock + vna_neighbors%Dblock     &
        &                       + vna_neighbors%Dblocko + vxc_neighbors%Dblock &
        &                       - SR_neighbors%Dblock + LR_neighbors%Dblock
            nullify (H_neighbors)
            nullify (K_neighbors, vna_neighbors)
            nullify (SR_neighbors, LR_neighbors, vxc_neighbors)
          end do
          nullify (gHamiltonian)
          nullify (gkinetic, gvna, gewaldsr, gewaldlr, gvxc)
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine build_gHmatrix

! ===========================================================================
! End Module
! ===========================================================================
        end module M_non_adiabatic
