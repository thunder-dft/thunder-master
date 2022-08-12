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

! M_Dassemble_usr
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
        module M_Dassemble_usr

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones

! /FDATA
        use M_Fdata_1c
        use M_Fdata_2c

! module procedures
        contains

! ===========================================================================
! Dassemble_uee.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates all the double counting corrections to the
!! Hartree interactions and is stored in uee.
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
        subroutine Dassemble_uee (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s  !< type structure

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

        real z                            !< distance between atom pairs
        real Zi, Zj

        real, dimension (:, :), allocatable :: coulomb
        real, dimension (:, :), allocatable :: dcoulomb
        real, dimension (:, :, :), allocatable :: vdcoulomb

        real, allocatable :: Q0 (:)

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2     !< positions of iatom and jatom
        real, dimension (3) :: sighat     !< unit vector along r2 - r1

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

! Allocate Arrays
! ===========================================================================
       allocate (Q0 (s%natoms))             !< neutral atom charge, i.e. ionic

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Initialize arrays
        ! Calculate nuclear charge.
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          Q0(iatom) = 0.0d0
          do issh = 1, species(in1)%nssh
            Q0(iatom) = Q0(iatom) + species(in1)%shell(issh)%Qneutral
          end do

          ! cut some lengthy notation
          pfi=>s%forces(iatom); pfi%usr = 0.0d0
        end do

! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%nssh

          Zi = Q0(iatom)
          r1 = s%atom(iatom)%ratom

          ! cut some lengthy notation
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%nssh

            Zj = Q0(jatom)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a

            ! cut some lengthy notation
            pfj=>s%forces(jatom)

            ! distance between the two atoms
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
            call Depsilon_2c (r1, r2, eps, deps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! GET COULOMB INTERACTIONS
! ****************************************************************************
! Now find the three coulomb integrals need to evaluate the neutral
! atom/neutral atom hartree interaction.
! For these Harris interactions, there are no subtypes and isorp = 0
            isubtype = 0
            interaction = P_coulomb
            in3 = in2
            allocate (coulomb (norb_mu, norb_nu))
            allocate (dcoulomb (norb_mu, norb_nu))
            allocate (vdcoulomb (3, norb_mu, norb_nu))
            call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,        &
     &                             norb_mu, norb_nu, coulomb, dcoulomb)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
            do issh = 1, norb_mu
              do jssh = 1, norb_nu
                if (z .gt. 1.0d-3) vdcoulomb(:,issh,jssh) = - eta(:)*dcoulomb(issh,jssh)
              end do
            end do

! Actually, now we calculate not only the neutral atom contribution,
! but also the short-ranged contribution due to the transfer of charge
! between the atoms:
!
! (Eii - Eee)neut - SUM(short range)(n(i) + dn(i))*dn(j)*J(i,j),
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! SPECIAL CASE: SELF-INTERACTION - NO FORCE

            else

! BONAFIDE TWO ATOM CASE
! Compute u0
!             u0(iatom)%neighbors(ineigh)%E = 0.0d0
              do issh = 1, species(in1)%nssh
                do jssh = 1, species(in2)%nssh
                  pfi%usr = pfi%usr                                           &
      &             + (P_eq2/2.0d0)*species(in1)%shell(issh)%Qneutral         &
      &              *species(in2)%shell(jssh)%Qneutral*vdcoulomb(:,issh,jssh)
                  pfj%usr = pfj%usr                                           &
      &             - (P_eq2/2.0d0)*species(in1)%shell(issh)%Qneutral         &
      &              *species(in2)%shell(jssh)%Qneutral*vdcoulomb(:,issh,jssh)
                end do
              end do
              pfi%usr = pfi%usr - (P_eq2/2.0d0)*eta(:)*(Zi*Zj/z**2)
              pfj%usr = pfj%usr + (P_eq2/2.0d0)*eta(:)*(Zi*Zj/z**2)
            end if
            deallocate (coulomb, dcoulomb, vdcoulomb)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
        deallocate (Q0)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_uee


!===========================================================================
! Dassemble_uxc
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the exchange-correlation double-counting
!! energy for Horsfield (Harris).
!
! We use a finite difference approach to change the densities and then
! find the corresponding changes in the exchange-correlation potential.
! This a bit clumsy sometimes, and probably can use a rethinking to improve.
! However, it actually works quite well for many molecular systems.
!
! We compute neutral cases for ideriv = 1. For other ideriv's we have the
! following KEY:
!
! Case 1 (KEY=1), neutral neutral corresponds to (00)
! KEY = 1,2,3,4,5 for ideriv = 1,2,3,4,5
!
! For the one-center case we only use ideriv = 1,2,3 as there is only one atom.
!
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
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
        subroutine Dassemble_uxc (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None
        
! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh, matom     !< counter over atoms/neighbors
                                         !< matom is the self-interaction atom
        integer in1, in2                  !< species numbers
        integer jatom                     !< neighbor of iatom
        integer interaction, isubtype     !< which interaction and subtype
        integer num_neigh                 !< number of neighbors
        integer mbeta                     !< the cell containing iatom's neighbor
        integer norb_mu, norb_nu          !< size of the block for the pair

        real z                           !< distance between atom pairs

        ! results
        real, dimension (0:0) :: uxc
        real, dimension (0:0) :: duxc

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2   !< positions of iatom and jatom
        real, dimension (3) :: sighat   !< unit vector along r2 - r1

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom

          ! cut some lengthy notation
          pfi=>s%forces(iatom)

! Now loop over all neighbors ineigh of iatom and put the force within the neighbors.
          do ineigh = 1, num_neigh
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some lengthy notation
            pfj=>s%forces(jatom)

            ! distance between the two atoms
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
            call Depsilon_2c (r1, r2, eps, deps)

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! GET DOUBLE-COUNTING EXCHANGE-CORRELATION FORCES
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! SPECIAL CASE: SELF-INTERACTION - NO FORCE
            else

! Now find the value of the integrals needed to evaluate the neutral
! atom/neutral atom exchange-correlation double-counting interaction.
! For these Harris interactions, there are no subtypes and ideriv = 0
              isubtype = 0
              interaction = P_uxc

              ! we set norb_mu = 1 and norb_nu = 1 because this interaction
              ! is not a matrix, but rather just a energy based on distances
              norb_mu = 1
              norb_nu = 1

              call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,      &
     &                               norb_mu, norb_nu, uxc, duxc)

! BONAFIDE TWO ATOM CASE
              pfi%usr = pfi%usr + eta(:)*(duxc(0)/2.0d0)
              pfj%usr = pfj%usr - eta(:)*(duxc(0)/2.0d0)
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None
 
! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_uxc

! End Module
! ===========================================================================
        end module M_Dassemble_usr
