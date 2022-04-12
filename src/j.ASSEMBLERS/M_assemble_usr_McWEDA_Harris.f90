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

! M_assemble_usr
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble the short range interactions related to the Hartree energies.
!! It contains the following subroutines within the module:
!!
!!       assemble_uee - assemble the double counting Hartree energy
!!       assemble_uxc - assemble the double counting xc correction
!!
! ===========================================================================
        module M_assemble_usr
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_1c
        use M_Fdata_2c

! module procedures
        contains

! ===========================================================================
! assemble_uee.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates all the double counting corrections to the
!! Hartree interactions and is stored in uee.
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
        subroutine assemble_uee (s, uii_uee)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s  !< type structure

! Output
        real uii_uee           !< double counting correction for coulomb part

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh             !< counter over atoms and neighbors
        integer in1, in2, in3             !< species numbers
        integer jatom                     !< neighbor of iatom
        integer interaction, isubtype     !< which interaction and subtype
        integer logfile                   !< writing to which unit
        integer num_neigh                 !< number of neighbors
        integer mbeta                     !< the cell containing iatom's neighbor
        integer issh, jssh                !< counting over orbitals
        integer norb_mu, norb_nu          !< size of the block for the pair

        real u0_tot
        real uee_self_tot
        real z                           !< distance between atom pairs
        real Zi, Zj

        real, allocatable :: coulomb (:, :)
        real, allocatable :: Q0 (:)
        real, dimension (3) :: r1, r2
        real, allocatable :: uee_self (:)

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_neighbors), pointer :: u0 (:)

! Allocate Arrays
! ===========================================================================
        allocate (u0 (s%natoms))             !< neutral atom uee
        allocate (uee_self (s%natoms))       !< self-interaction for uee
        allocate (Q0 (s%natoms))             !< neutral atom charge, i.e. ionic

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        write (logfile,*)
        write (logfile,*) ' Welcome to assemble_usr.f! '

! Initialize arrays
        u0_tot = 0.0d0
        uee_self_tot = 0.0d0
        ! Calculate nuclear charge.
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          Q0(iatom) = 0.0d0
          do issh = 1, species(in1)%nssh
            Q0(iatom) = Q0(iatom) + species(in1)%shell(issh)%Qneutral
          end do
        end do

! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%nssh
          num_neigh = s%neighbors(iatom)%neighn
          allocate (u0(iatom)%neighbors(num_neigh))
          Zi = Q0(iatom)

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%nssh
            Zj = Q0(jatom)

            ! distance between the two atoms
            z = distance (r1, r2)

! GET COULOMB INTERACTIONS
! ****************************************************************************
! Now find the three coulomb integrals need to evaluate the neutral
! atom/neutral atom hartree interaction.
! For these Harris interactions, there are no subtypes and isorp = 0
            isubtype = 0
            interaction = P_coulomb
            in3 = in2

            allocate (coulomb (norb_mu, norb_nu))
            call getMEs_Fdata_2c (in1, in2, interaction, isubtype, z,         &
     &                             norb_mu, norb_nu, coulomb)

! Actually, now we calculate not only the neutral atom contribution,
! but also the short-ranged contribution due to the transfer of charge
! between the atoms:
!
! (Eii - Eee)neut - SUM(short range)(n(i) + dn(i))*dn(j)*J(i,j),
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! SPECIAL CASE: SELF-INTERACTION
              uee_self(iatom) = 0.0d0
              do issh = 1, species(in1)%nssh
                do jssh = 1, species(in1)%nssh
                  uee_self(iatom) = uee_self(iatom)                           &
     &              + (P_eq2/2.0d0)*species(in1)%shell(issh)%Qneutral         &
     &                *species(in1)%shell(jssh)%Qneutral*coulomb(issh,jssh)
                end do
              end do

! put the half and the units in:
              u0(iatom)%neighbors(ineigh)%E = 0.0d0
            else

! BONAFIDE TWO ATOM CASE
! Compute u0
              u0(iatom)%neighbors(ineigh)%E = 0.0d0
              do issh = 1, species(in1)%nssh
                do jssh = 1, species(in2)%nssh
                  u0(iatom)%neighbors(ineigh)%E =                             &
     &              u0(iatom)%neighbors(ineigh)%E                             &
     &              + species(in1)%shell(issh)%Qneutral                       &
     &                *species(in2)%shell(jssh)%Qneutral*coulomb(issh,jssh)
                end do
              end do
              u0(iatom)%neighbors(ineigh)%E = (P_eq2/2.0d0)*(Zi*Zj/z          &
     &          - u0(iatom)%neighbors(ineigh)%E)
            end if
            deallocate (coulomb)

! Compute the total cell value of uii-uee; uii-uee = sum u0(i,m) - sum uee00(i)
            u0_tot = u0_tot + u0(iatom)%neighbors(ineigh)%E
          end do ! end loop over neighbors
          uee_self_tot = uee_self_tot + uee_self(iatom)
        end do ! end loop over atoms

        ! Final double-counting correction for coulomb part
        uii_uee = u0_tot - uee_self_tot

! Deallocate Arrays
! ===========================================================================
        deallocate (Q0, uee_self)
        do iatom = 1, s%natoms
          deallocate (u0(iatom)%neighbors)
        end do
        deallocate (u0)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_uee


!===========================================================================
! assemble_uxc
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the exchange-correlation double-counting
!! energy for McWEDA (Harris).
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
        subroutine assemble_uxc (s, uxcdcc)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Output
        real uxcdcc           !< double counting correction for coulomb part

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                    !< counter over atoms
        integer in1                      !< species numbers
        integer matom                    !< matom is the self-interaction atom

        integer imu
        integer norb_mu                  !< size of the block for the pair

! Loops
        integer issh, n1, l1, ind1

! Ceperley-Adler
        real arhoin
        real arhobond
        real ovlap
        real rhoin
        real rhobond
        real dexc_in, dexc_bond
        real d2exc_in, d2exc_bond
        real dmuxc_in, dmuxc_bond
        real exc_in, exc_bond
        real muxc_in, muxc_bond
        real d2muxc_in, d2muxc_bond

        real e_xc_sn, e_xc_bond_sn, e_xc_bond
        real e_vxc_sn, e_vxc_bond_sn, e_vxc_bond
        real uxcdc_sn, uxcdc_bond_sn, uxcdc_bond
        real q_mu

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize
        e_xc_sn = 0.d0
        e_xc_bond_sn = 0.0d0
        e_xc_bond = 0.d0
        e_vxc_sn = 0.d0
        e_vxc_bond_sn = 0.0d0
        e_vxc_bond = 0.d0

! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          matom = s%neigh_self(iatom)
! only "diagonal" terms for the exchange-correlation energy

! Loop over shells i-atom
          n1 = 0
          do issh = 1, species(in1)%nssh
! Calculate exc_bond : xc-energy from first term on the right in Eq. (16)
!    see PRB 71, 235101 (2005)
! e_vxc_bond : energy already included in the band-structure through vxc_bond
            e_xc_bond = e_xc_bond + vxc_1c(in1)%E(issh,issh)*s%atom(iatom)%shell(issh)%Qneutral
            e_vxc_bond = e_vxc_bond + vxc_1c(in1)%V(issh,issh)*s%atom(iatom)%shell(issh)%Qneutral
            l1 = species(in1)%shell(issh)%lssh
            n1 = n1 + l1 + 1
            q_mu = s%atom(iatom)%shell(issh)%Qneutral / (2*l1 +1)
            arhoin = s%rho_in_weighted(iatom)%neighbors(matom)%block(issh,issh)
            call lda_ceperley_alder (arhoin, exc_in, muxc_in, dexc_in, d2exc_in, dmuxc_in, d2muxc_in)
            arhobond = s%rho_bond_weighted(iatom)%neighbors(matom)%block(issh,issh)
            call lda_ceperley_alder (arhobond, exc_bond, muxc_bond, dexc_bond, d2exc_bond, dmuxc_bond, d2muxc_bond)

! Set the xc-submatrices
! loop over orbitals in the iatom-shell
            do ind1 = -l1, l1
              imu = n1 + ind1
              ovlap = s%overlap(iatom)%neighbors(matom)%block(imu,imu)
              rhoin = s%rho_in(iatom)%neighbors(matom)%block(imu,imu)
              rhobond = s%rho_bond(iatom)%neighbors(matom)%block(imu,imu)
! calculate SN part
! exc_sn : xc-energy from second term on the right in Eq. (16): PRB 71, 235101 (2005)
! e_vxc_sn : energy already included in the band-structure through vxc_sn
              e_xc_sn = e_xc_sn + q_mu*(exc_in*ovlap + dexc_in*(rhoin - arhoin*ovlap))
              e_vxc_sn = e_vxc_sn + q_mu*(muxc_in*ovlap + dmuxc_in*(rhoin - arhoin*ovlap))
! calculate SN-AT part ("atomic" correction)
! exc_sn_bond : xc-energy from third term on the right in Eq. (16): PRB 71, 235101 (2005)
! e_vxc_bond_sn : energy already included in the band-structure through vxc_sn_bond
              e_xc_bond_sn = e_xc_bond_sn + q_mu*(exc_bond*ovlap + dexc_bond*(rhobond - arhobond*ovlap))
              e_vxc_bond_sn = e_vxc_bond_sn + q_mu*(muxc_bond*ovlap + dmuxc_bond*(rhobond - arhobond*ovlap))
            end do !do ind1 = -l1, l1
            n1 = n1 + l1
          end do !do issh = 1, nssh(in1)
        end do ! end loop over atoms

! The different contributions then are:
        uxcdc_sn = e_xc_sn - e_vxc_sn
        uxcdc_bond_sn = e_xc_bond_sn - e_vxc_bond_sn
        uxcdc_bond = e_xc_bond - e_vxc_bond

! The total double-counting XC term ( Eq. (16): PRB 71, 235101 (2005) )
        uxcdcc = uxcdc_sn + uxcdc_bond - uxcdc_bond_sn

! Deallocate Arrays
! ===========================================================================
! None
 
! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_uxc

! End Module
! ===========================================================================
        end module M_assemble_usr

