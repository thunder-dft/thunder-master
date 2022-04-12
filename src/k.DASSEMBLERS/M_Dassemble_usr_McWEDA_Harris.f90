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
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_1c
        use M_Fdata_2c
        use M_Dassemble_rho_McWEDA 
        use M_Dassemble_2c 

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

        real z                           !< distance between atom pairs
        real Zi, Zj

        real, dimension (:, :), allocatable :: coulomb
        real, dimension (:, :), allocatable :: dcoulomb
        real, dimension (:, :, :), allocatable :: vdcoulomb

        real, allocatable :: Q0 (:)

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
        end do

! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)

          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%nssh
          Zi = Q0(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%nssh
            Zj = Q0(jatom)

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
            call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,       &
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
      &              + (P_eq2/2.0d0)*species(in1)%shell(issh)%Qneutral         &
      &                *species(in2)%shell(jssh)%Qneutral*vdcoulomb(:,issh,jssh)
                   pfj%usr = pfj%usr                                           &
      &              - (P_eq2/2.0d0)*species(in1)%shell(issh)%Qneutral         &
      &                *species(in2)%shell(jssh)%Qneutral*vdcoulomb(:,issh,jssh)
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
        integer iatom                     !< counter over atoms and neighbors
!        integer ineigh, mneigh            !< counter over atoms and neighbors
        integer in1!                     !< species numbers
!        integer in2, in3
        integer jatom                    !< jatom is the neighbor of iatom ialpha is the third atom for 3c
!        integer ialpha
!       integer logfile                  !< writing to which unit
        integer num_neigh                 !< number of neighbors
        integer matom                      !< matom is the self-interaction atom
 !       integer mbeta                     !< the cell containing neighbor of iatom

        integer imu                      !< counter over orbitals
 !       integer inu                      !< counter over orbitals
        integer issh                     !< counter over shells
 !       integer jssh                     !< counter over shells
        integer n1, l1, m1               !< quantum numbers n, l, and m
 !       integer n2, l2, m2               !< quantum numbers n, l, and m
        integer norb_mu, norb_nu         !< size of the block for the pair

        real prho_in_shell               !< temporary storage
        real prho_bond_shell             !< temporary storage
!        real poverlap                    !< temporary stroage
        real prho_in                     !< temporary storage
        real prho_bond                   !< temporary storage
        real dexc_in, dexc_bond          !< 1st derivative of xc energy
        real d2exc_in, d2exc_bond        !< 2nd derivative of xc energy
        real dmuxc_in, dmuxc_bond        !< 1st derivative of xc potential
        real exc_in, exc_bond            !< xc energy
        real muxc_in, muxc_bond          !< xc potential
        real d2muxc_in, d2muxc_bond      !< 2nd derivative of xc potential
        real q_mu                        !< charge on shell
        
        real, dimension (3) :: e_xc_sn      ! this terms are now forces
        real, dimension (3) :: e_xc_bond_sn
        real, dimension (3) :: e_vxc_sn
        real, dimension (3) :: e_vxc_bond_sn
        
        ! vector derivitives of rho_in and rho_bond interactions
        real, dimension (3) :: Dprho_in_weighted
        real, dimension (3) :: Dprho_in
        real, dimension (3) :: Dprho_bond_weighted
        real, dimension (3) :: Dprho_bond
!        real, dimension (3) :: Dpoverlap
        
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize
        e_xc_sn = 0.d0
        e_xc_bond_sn = 0.0d0
        e_vxc_sn = 0.d0
        e_vxc_bond_sn = 0.0d0
        
! Loop over the atoms in the central cell.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)
          matom = s%neigh_self(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          jatom = s%neighbors(iatom)%neigh_j(matom)
          pfj=>s%forces(jatom)
          
! Allocate block size
          norb_nu = species(in1)%norb_max

! SPECIAL LOOP: we want to minimize the number of calls to lda-function
! we only need to call lda_ceperley-adler for each pair of shells
! but then we need to calculate the (mu,nu)-block of matrix elements

! Loop over shells i-atom
          n1 = 0
          do issh = 1, species(in1)%nssh

! n1 : counter used to determine orbitals imu
            l1 = species(in1)%shell(issh)%lssh
            n1 = n1 + l1 + 1
            q_mu = s%atom(iatom)%shell(issh)%Qneutral / (2*l1 +1)
! Call lda-function for rho_in
            prho_in_shell =                                                  &
     &        s%rho_in_weighted(iatom)%neighbors(matom)%block(issh,issh)
              call lda_ceperley_alder (prho_in_shell, exc_in, muxc_in,       &
     &                                 dexc_in, d2exc_in, dmuxc_in, d2muxc_in)
                  
            prho_bond_shell =                                                &
     &        s%rho_bond_weighted(iatom)%neighbors(matom)%block(issh,issh)
              call lda_ceperley_alder (prho_bond_shell, exc_bond,            &
     &                                 muxc_bond, dexc_bond, d2exc_bond,     &
     &                                 dmuxc_bond, d2muxc_bond)

! Calculate vxc_SN and vxc_SN_bond for (mu,nu)-block
! loop over orbitals in the iatom-shell (imu)
            do m1 = -l1, l1
              imu = n1 + m1
              prho_in = s%rho_in(iatom)%neighbors(matom)%block(imu,imu)
              Dprho_in_weighted = s%rho_in_weighted(iatom)%neighbors(matom)%Dblock(:,issh,issh)
              Dprho_in = s%rho_in(iatom)%neighbors(matom)%Dblock(:,imu,imu)
              prho_bond = s%rho_bond(iatom)%neighbors(matom)%block(imu,imu)
              Dprho_bond_weighted(:) = s%rho_bond_weighted(iatom)%neighbors(matom)%Dblock(:,issh,issh)
              Dprho_bond(:) = s%rho_bond(iatom)%neighbors(matom)%Dblock(:,imu,imu)
! calculate GSN for rho_in
              e_vxc_sn = q_mu*(dmuxc_in*Dprho_in_weighted                    &
     &          + Dprho_in_weighted*d2muxc_in*(prho_in -prho_in_shell)       &
     &          + dmuxc_in*(Dprho_in - Dprho_in_weighted)) 
     
              e_xc_sn = q_mu*(dexc_in*Dprho_in_weighted                      &
     &          + Dprho_in_weighted*d2exc_in*(prho_in -prho_in_shell)        &
     &          + dexc_in*(Dprho_in - Dprho_in_weighted))        
              
              e_vxc_bond_sn = q_mu*(dmuxc_bond*Dprho_bond_weighted           &
     &          + Dprho_bond_weighted*d2muxc_bond*(prho_bond -prho_bond_shell) &
     &          + dmuxc_bond*(Dprho_bond - Dprho_bond_weighted)) 
     
              e_xc_bond_sn = q_mu*(dexc_bond*Dprho_bond_weighted             &
     &          + Dprho_bond_weighted*d2exc_bond*(prho_bond -prho_bond_shell)  &
     &          + dexc_bond*(Dprho_bond - Dprho_bond_weighted))
             
              pfi%usr = pfi%usr - (e_xc_sn - e_vxc_sn - (e_xc_bond_sn - e_vxc_bond_sn))
              pfj%usr = pfj%usr - (e_xc_sn - e_vxc_sn - (e_xc_bond_sn - e_vxc_bond_sn))
            end do
            n1 = n1 + l1
          end do  ! do issh = 1, nssh(in1)
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
