! copyright info:
!
!                             @Copyright 2008
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

! M_assemble_vxc
! Module Description
! ===========================================================================
!>       This is a module containing all of the  programs required
!! to assemble all of the matrix elements for exchange-correlation -
!! see PRB 71, 235101 (2005).
!!
!! It contains the following subroutines within the module:
!!
!!       assemble_vxc : XC-main  program. Driver
!!
!! Definition of rho_in and rho_local:
!!       rho_bond (mu,nu) = < mu | rho_i | nu >
!!         if mu and nu are in the same atom "i" : onsite case
!!       rho_bond (mu,nu) = < mu | rho_i + rho_j | nu >
!!         if mu and nu are in different atoms "i" and "j" : atom case
!!       rho_in = sum of onsite rho_bond values
!!
!!
!!           vxc_bond (mu,nu) = < mu | V_xc (rho_i) | nu >
!!                            if mu and nu are in the same atom "i"; or
!!           vxc_bond (mu,nu) = < mu | V_xc(rho_i + rho_j) | nu >
!!                            if mu and nu are in different atoms "i" and "j"
!
! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
! located in the Fdata directory.  This list will change depending on
! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_assemble_vxc
        use M_assemble_rho_McWEDA
        use M_assemble_2c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_1c
        use M_Fdata_2c
        use M_Fdata_3c
        use M_neighbors
        use M_rotations

! Type Declaration
! ===========================================================================
! Output
        type(T_assemble_neighbors), pointer :: vxc (:)

! module procedures
        contains

! ===========================================================================
! assemble_vxc.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>  This is the main module for assembling the Vxc matrix element
!! interactions (McWEDA).  Subroutines from M_assemble_rho_McWEDA_rho
!! are used. The results are stored in vxc (potential).
! ===========================================================================
! Code written by:
!> @author Daniel G. Trabada
!! @author Jose Ortega (JOM)
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_vxc (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh            !< counter over atoms and neighbors
        integer in1, in2                 !< species numbers
        integer jatom                    !< neighbor of iatom
        integer num_neigh                !< number of neighbors

        integer norb_mu, norb_nu         !< size of the block for the pair

        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc

! Allocate Arrays
! ===========================================================================
        allocate (s%vxc(s%natoms))

! Procedure
! ===========================================================================
! Calculate rho_in (density) matrix elements
        write (s%logfile,*) ' Calling rho (density) input assemblers. '
        call assemble_rho_2c (s)
        call assemble_rho_3c (s)

! calculate average_rho matrix elements
! See PRB 71, 235101 (2005), Eqs. (19), (22) and (25)
        call assemble_rho_weighted_2c (s)
        call assemble_rho_weighted_3c (s)

! calculate  XC-potential matrix elements
! See PRB 71, 235101 (2005), Eqs. (16), (21) and (24)
        write (s%logfile,*) ' Calling vxc assemblers. '

! (3) Sum the 3-contributions in Eq. (16):
!  vxc = vxc_bond + vxc_SN - vxc_SN_bond
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pvxc=>s%vxc(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pvxc%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pvxc_neighbors=>pvxc%neighbors(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            allocate (pvxc_neighbors%block(norb_mu, norb_nu))

! ecuacion (16) PRB 71, 235101 (2005)
!  vxc = vxc_bond + vxc_SN - vxc_SN_bond
   !         pvxc_neighbors%block = vxc_SN(iatom)%neighbors(ineigh)%block     &
   !  &                            + vxc_bond(iatom)%neighbors(ineigh)%block  &
   !  &                            - vxc_SN_bond(iatom)%neighbors(ineigh)%block
          end do
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
        end subroutine assemble_vxc


! ===========================================================================
! assemble_vxc_SN
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates/assembles the generalized Sankey-Niklewski
!! exchange-correlation potential matrix elements for
!!
!!               rho_in (input density) --> (vxc_SN)
!!       and
!!               rho_bond (local or "atomic" (_at) density) --> vxc_SN_bond
!!
!! Definition of rho_in and rho_local:
!!       rho_bond (mu,nu) = < mu | rho_i | nu >
!!         if mu and nu are in the same atom "i" : onsite case
!!       rho_bond (mu,nu) = < mu | rho_i + rho_j | nu >
!!         if mu and nu are in different atoms "i" and "j" : atom case
!!       rho_in = sum of onsite rho_bond values
!
! ===========================================================================
! Code written by:
!> @author Daniel G. Trabada
!! @author Jose Ortega (JOM)
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================


!
! Program Declaration
! ===========================================================================
! assemble_vxc_bond.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine assembles the XC-potential matrix elements for vxc_bond
!!
!!    vxc_bond (mu,nu) = < mu | V_xc (rho_i) | nu >
!!                   if mu and nu are in the same atom "i"; or
!!    vxc_bond (mu,nu) = < mu | V_xc(rho_i + rho_j) | nu >
!!                  if mu and nu are in different atoms "i" and "j"
!!
!!    (first term on the right in Eqs. (16), (21) and (24)
!!    see PRB 71, 235101 (2005)
!

! ===========================================================================
! destroy_assemble_vxc_McWEDA
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the vxc
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
        subroutine destroy_assemble_vxc_McWEDA (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                            !< counter over atoms
        integer ineigh                           !< counter over neighbors

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%vxc(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%vxc(iatom)%neighbors)
        end do
        deallocate (s%vxc)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_assemble_vxc_McWEDA

! End Module
! ===========================================================================
        end module M_assemble_vxc

