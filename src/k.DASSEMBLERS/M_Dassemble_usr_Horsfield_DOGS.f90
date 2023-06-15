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

        real, dimension (3) :: dcorksr

        real, dimension (:, :), allocatable :: coulomb
        real, dimension (:, :), allocatable :: dcoulomb
        real, dimension (:, :, :), allocatable :: vdcoulomb

        real, allocatable :: Q0 (:)
        real, allocatable :: Q(:)         !< total charge on atom

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
        allocate (Q (s%natoms))              !< total input charge on atom

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        write (logfile,*)
        write (logfile,*) ' Welcome to Dassemble_usr.f! '

! Initialize arrays
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
            call Depsilon_2c (r1, r2, z, eps, deps)

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
                   pfi%usr = pfi%usr                                          &
      &              + (P_eq2/2.0d0)*s%atom(iatom)%shell(issh)%Qin            &
      &                *s%atom(jatom)%shell(jssh)%Qin*vdcoulomb(:,issh,jssh)
                   pfj%usr = pfj%usr                                          &
      &              - (P_eq2/2.0d0)*s%atom(iatom)%shell(issh)%Qin            &
      &                *s%atom(jatom)%shell(jssh)%Qin*vdcoulomb(:,issh,jssh)
                end do
              end do
              pfi%usr = pfi%usr - (P_eq2/2.0d0)*eta(:)*(Zi*Zj/z**2)
              pfj%usr = pfj%usr + (P_eq2/2.0d0)*eta(:)*(Zi*Zj/z**2)

              ! force due dcorksr
              dcorksr(:) = - (P_eq2/2.0d0)*eta(:)*(Zi*Zj - Q(iatom)*Q(jatom))/z**2
              pfi%usr = pfi%usr - dcorksr
              pfj%usr = pfj%usr + dcorksr
            end if
            deallocate (coulomb, dcoulomb, vdcoulomb)
          end do ! end loop over neighbors

          ! add in ewald contributions
          pfi%usr = pfi%usr - (P_eq2/2.0d0)*pfi%ewald
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
        deallocate (Q, Q0)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_uee


!===========================================================================
! assemble_uxc
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the exchange-correlation double-counting
!! energy for McWEDA (Harris).
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
!> @author Daniel G. Trabada
!! @author Jose Ortega (JOM)
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
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
        integer iatom, ineigh             !< counter over atoms and neighbors
        integer ispecies, in1, in2        !< species numbers
        integer jatom                     !< neighbor of iatom
        integer interaction, isubtype     !< which interaction and subtype
        integer num_neigh                 !< number of neighbors
        integer mbeta                     !< the cell containing iatom's neighbor

        integer issh                      !< counting over orbitals
        integer norb_mu, norb_nu          !< size of the block for the pair

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

! Some parameters for the derivative parts:
        integer, parameter, dimension (0:4) :: jsign = (/0, -1, +1, -1, +1/)
        integer, parameter, dimension (0:4) :: kspecies_key = (/1, 1, 1, 2, 2/)

        real z                           !< distance between atom pairs

! The charge transfer bit = need a +-dq on each orbital.
        real, allocatable :: dqorb (:)
        real, allocatable :: dQ (:)    !< charge on atom, i.e. ionic
        real, allocatable :: Q0 (:)    !< total neutral atom charge, i.e. Ztot
        real, allocatable :: Q(:)      !< total charge on atom

        ! results
        real duxcdc_bond
        real, dimension (0:4) :: uxc
        real, dimension (0:4) :: duxc

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
        allocate (dqorb (nspecies))
        allocate (dQ (s%natoms))
        allocate (Q0 (s%natoms))
        allocate (Q (s%natoms))

! Procedure
! ===========================================================================
! Initialize the charge transfer bit
        do ispecies = 1, nspecies
          dqorb(ispecies) = 0.5d0
          if (species(ispecies)%nssh .eq. 1) dqorb(ispecies) = 0.25d0
        end do

! Calculate nuclear charge.
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
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          pfi=>s%forces(iatom)

! Calculate exc_bond : energy from first one center piece of exchange-
! correlation already included in the band-structure through vxc_1c
! No force due to this one-center term
!         uxcdc_bond = uxcdc_bond + vxc_1c(in1)%E

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
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
            call Depsilon_2c (r1, r2, z, eps, deps)

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

! GET DOUBLE-COUNTING EXCHANGE-CORRELATION INTERACTIONS
! ****************************************************************************
! Now find the value of the integrals needed to evaluate the neutral
! atom/neutral atom exchange-correlation double-counting interaction.
!
! Now we should interpolate - the fast way is:
!
!   e(dqi,dqj) = exc(0,0) + dqi*(exc(1,0) - exc(0,0))/dQ
!                         + dqj*(exc(0,1) - exc(0,0))/dQ
!
! The good way is to use a three point Lagrange interpolation along the axis.
!
! Lagrange: f(x) = f(1)*L1(x) + f(2)*L2(x) + f(3)*L3(x)
!
! L1(x) = (x - x2)/(x1 - x2)*(x - x3)/(x1 - x3)
! L2(x) = (x - x3)/(x2 - x1)*(x - x3)/(x2 - x3)
! L3(x) = (x - x1)/(x3 - x1)*(x - x2)/(x3 - x2)
!
! in our case:
!
! L1(dq) = (1/2)*dq*(dq - 1)
! L2(dq) = -(dq + 1)*(dq - 1)
! L3(dq) = (1/2)*dq*(dq + 1)
!
! The interpolation does not depend on the (qi,qj) quadrant:
!
!  f(dqi,dqj) = f(0,dqj) + f(dqi,0) - f(0,0)
! ****************************************************************************
              isubtype = 0
              interaction = P_uxc

              ! we set norb_mu = 1 and norb_nu = 1 because this interaction
              ! is not a matrix, but rather just a energy based on distances
              norb_mu = 1
              norb_nu = 1
              call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,      &
     &                               norb_mu, norb_nu, uxc(0), duxc(0))

              ! Add in the neutral atom piece
              pfi%usr = pfi%usr + eta(:)*(duxc(0)/2.0d0)
              pfj%usr = pfj%usr - eta(:)*(duxc(0)/2.0d0)

! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
              ideriv_min = 1
              ideriv_max = 4

              ! Now loop over the different ideriv subtypes
              do ideriv = ideriv_min, ideriv_max
                isubtype = ideriv
                call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,    &
     &                                 norb_mu, norb_nu, uxc(isubtype),       &
     &                                 duxc(isubtype))
              end do

              ! add the uxc from interpolation to the total
              ! First, we look at the sign of the charge and "pick" the
              ! correct differences to add accordingly
              duxcdc_bond = 0.0d0
              if (dQ(iatom) .gt. 0.0d0 .and. dQ(jatom) .gt. 0.0d0) then
                duxcdc_bond = dQ(iatom)*(duxc(2) - duxc(0))/dqorb(in1)        &
     &                        + dQ(jatom)*(duxc(4) - duxc(0))/dqorb(in2)
              else if (dQ(iatom) .lt. 0.0d0 .and. dQ(jatom) .lt. 0.0d0) then
                duxcdc_bond = dQ(iatom)*(duxc(0) - duxc(1))/dqorb(in1)        &
     &                        + dQ(jatom)*(duxc(0) - duxc(3))/dqorb(in2)
              else if (dQ(iatom) .gt. 0.0d0 .and. dQ(jatom) .lt. 0.0d0) then
                duxcdc_bond = Q(iatom)*(duxc(2) - duxc(0))/dqorb(in1)         &
     &                        + dQ(jatom)*(duxc(0) - duxc(3))/dqorb(in2)
              else if (dQ(iatom) .lt. 0.0d0 .and. dQ(jatom) .gt. 0.0d0) then
                duxcdc_bond = dQ(iatom)*(duxc(0) - duxc(1))/dqorb(in1)        &
     &                        + dQ(jatom)*(duxc(4) - duxc(0))/dqorb(in2)
              end if

              ! Add in the charged atom piece
              pfi%usr = pfi%usr + eta(:)*(duxcdc_bond/2.0d0)
              pfj%usr = pfj%usr - eta(:)*(duxcdc_bond/2.0d0)
            end if ! end if for r1 .eq. r2 case
          end do ! end loop over neighbors
        end do ! end loop over atoms


! Deallocate Arrays
! ===========================================================================
        deallocate (dqorb, dQ, Q, Q0)

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
