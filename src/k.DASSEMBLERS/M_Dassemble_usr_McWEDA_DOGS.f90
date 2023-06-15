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

! /DASSEMBLERS
        use M_Dassemble_rho_McWEDA 
        use M_Dassemble_2c 

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
        integer in1                      !< species numbers
        integer jatom, num_neigh         !< counters over neighbors
        integer mbeta                    !< the cell containing neighbor of iatom

        integer imu                      !< counter over orbitals
        integer issh                     !< counter over shells
        integer n1, l1, m1               !< quantum numbers n, l, and m
        integer norb_mu                  !< size of the block for the pair

        real prho_in                     !< temporary storage
        real prho_in_shell               !< temporary storage
        real prho_bond                   !< temporary storage
        real prho_bond_shell             !< temporary storage

        ! Ceperley-Alder terms
        real dexc_in, dexc_bond          !< 1st derivative of xc energy
        real d2exc_in, d2exc_bond        !< 2nd derivative of xc energy
        real dmuxc_in, dmuxc_bond        !< 1st derivative of xc potential
        real exc_in, exc_bond            !< xc energy
        real muxc_in, muxc_bond          !< xc potential
        real d2muxc_in, d2muxc_bond      !< 2nd derivative of xc potential

        real q_mu                        !< charges
        real z                           !< distance between atom pairs

        ! results
        real, dimension (3) :: de_xc_sn      ! this terms are the forces
        real, dimension (3) :: de_xc_bond_sn
        real, dimension (3) :: de_vxc_sn
        real, dimension (3) :: de_vxc_bond_sn

        ! vector derivitives of rho_in and rho_bond interactions
        real, dimension (3) :: Dprho_in
        real, dimension (3) :: Dprho_in_shell
        real, dimension (3) :: Dprho_bond
        real, dimension (3) :: Dprho_bond_shell

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

!           Find forces for uxcdcc
!           uxcdcc = uxcdc_bond + uxcdc_sn - uxcdc_bond_sn

! This term has zero force because it is one-center energy
!           uxcdc_bond = e_xc_bond - e_vxc_bond

!           forces for the two energy terms:
!           uxcdc_sn = e_xc_sn - e_vxc_sn
!           uxcdc_bond_sn = e_xc_bond_sn - e_vxc_bond_sn
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! SPECIAL CASE: SELF-INTERACTION - NO FORCE

            else
! SPECIAL LOOP: we want to minimize the number of calls to lda-function
! we only need to call lda_ceperley-adler for each pair of shells
! but then we need to calculate the (mu,nu)-block of matrix elements

! Initialize
              de_xc_sn = 0.d0; de_vxc_sn = 0.0d0
              de_xc_bond_sn = 0.0d0; de_vxc_bond_sn = 0.0d0

! Loop over shells i-atom
              n1 = 0
              do issh = 1, species(in1)%nssh

! n1 : counter used to determine orbitals imu
                l1 = species(in1)%shell(issh)%lssh
                n1 = n1 + l1 + 1
                q_mu = s%atom(iatom)%shell(issh)%Qin/(2*l1+1)

! Call lda-function for rho_in
                prho_in_shell =                                               &
     &            s%rho_in_weighted(iatom)%neighbors(matom)%block(issh,issh)
                Dprho_in_shell =                                              &
                  s%rho_in_weighted(iatom)%neighbors(ineigh)%Dblock(:,issh,issh)
                  call lda_ceperley_alder (prho_in_shell, exc_in, muxc_in,    &
     &                                     dexc_in, d2exc_in, dmuxc_in, d2muxc_in)
                  
                prho_bond_shell =                                             &
     &            s%rho_bond_weighted(iatom)%neighbors(matom)%block(issh,issh)
                Dprho_bond_shell = Dprho_in_shell
                call lda_ceperley_alder (prho_bond_shell, exc_bond,           &
     &                                   muxc_bond, dexc_bond, d2exc_bond,    &
     &                                   dmuxc_bond, d2muxc_bond)

! Calculate vxc_SN and vxc_SN_bond for (mu,nu)-block
! loop over orbitals in the iatom-shell (imu)
                do m1 = -l1, l1
                  imu = n1 + m1
                  prho_in = s%rho_in(iatom)%neighbors(matom)%block(imu,imu)
                  Dprho_in = s%rho_in(iatom)%neighbors(ineigh)%Dblock(:,imu,imu)

                  prho_bond = s%rho_bond(iatom)%neighbors(matom)%block(imu,imu)
                  Dprho_bond = Dprho_in

! calculate SN part forces
! exc_sn : xc-energy from second term on the right in Eq. (16): PRB 71, 235101 (2005)
! e_vxc_sn : energy already included in the band-structure through vxc_sn
                  de_xc_sn = de_xc_sn + q_mu*(dexc_in*Dprho_in_shell             &
     &              + Dprho_in_shell*d2exc_in*(prho_in - prho_in_shell)          &
     &              + dexc_in*(Dprho_in - Dprho_in_shell))
                  de_vxc_sn = de_vxc_sn + q_mu*(dmuxc_in*Dprho_in_shell          &
     &              + Dprho_in_shell*d2muxc_in*(prho_in - prho_in_shell)         &
     &              + dmuxc_in*(Dprho_in - Dprho_in_shell))

! calculate SN-AT part ("atomic" correction) forces
! exc_sn_bond : xc-energy from third term on the right in Eq. (16): PRB 71, 235101 (2005)
! e_vxc_bond_sn : energy already included in the band-structure through vxc_sn_bond
                 de_xc_bond_sn = de_xc_bond_sn                                   &
     &             + q_mu*(dexc_bond*Dprho_bond_shell                            &
     &             + Dprho_bond_shell*d2exc_bond*(prho_bond - prho_bond_shell)   &
     &             + dexc_bond*(Dprho_bond - Dprho_bond_shell))
                 de_vxc_bond_sn = de_vxc_bond_sn                                 &
     &             + q_mu*(dmuxc_bond*Dprho_bond_shell                           &
     &             + Dprho_bond_shell*d2muxc_bond*(prho_bond - prho_bond_shell)  &
     &             + dmuxc_bond*(Dprho_bond - Dprho_bond_shell))
                end do ! end loop m1 = -l1, l1
                n1 = n1 + l1
              end do  ! do issh = 1, nssh(in1)
!             Find forces for uxcdcc
!             uxcdcc = uxcdc_bond + uxcdc_sn - uxcdc_bond_sn

! This term has zero force because it is one-center energy
!             uxcdc_bond = e_xc_bond - e_vxc_bond

!             forces for the two energy terms:
!             uxcdc_sn = e_xc_sn - e_vxc_sn
!             uxcdc_bond_sn = e_xc_bond_sn - e_vxc_bond_sn

              pfi%usr = pfi%usr - ((de_xc_sn - de_vxc_sn) + (de_xc_bond_sn - de_vxc_bond_sn))
              pfj%usr = pfj%usr + ((de_xc_sn - de_vxc_sn) - (de_xc_bond_sn - de_vxc_bond_sn))
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
