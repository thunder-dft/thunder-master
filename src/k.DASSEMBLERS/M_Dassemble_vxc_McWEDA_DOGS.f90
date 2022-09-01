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

! M_Dassemble_vxc
! Module Description
! ===========================================================================
!>       This is a module containing all of the  programs required
!! to assemble all the deriatives of the matrix elements for exchange-correlation -
!! see PRB 71, 235101 (2005).
!!
!! It contains the following subroutines within the module:
!!
!!       Dassemble_vxc : XC-main  program. Driver
!!       Dassemble_vxc_SN : assembles the generalized Sankey-Niklewski both for
!!
!!               rho_in (input density) --> (vxc_SN)
!!           and
!!               rho_bond (local or "atomic" (_at) density) --> vxc_SN_bond
!!
!! Definition of rho_in and rho_local:
!!       rho_bond (mu,nu) = < mu | rho_i | nu >
!!         if mu and nu are in the same atom "i" : onsite case
!!       rho_bond (mu,nu) = < mu | rho_i + rho_j | nu >
!!         if mu and nu are in different atoms "i" and "j" : atom case
!!       rho_in = sum of onsite rho_bond values
!!
!!       Dassemble_vxc_bond : assembles the derivative XC-potential matrix elements
!!
!!           d (vxc_bond (mu,nu)/dR =  d (< mu | V_xc(rho_i + rho_j) | nu >)/dR
!!
!!      Dassemble_vxc_3c : assmebles the derivatives for the three center interactions
!!                         of the XC-potential matrix elements
!
! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
! located in the Fdata directory.  This list will change depending on
! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_Dassemble_vxc

! /SYSTEM
        use M_assemble_blocks
        use M_configuraciones

! /ASSEMBLERS
        use M_assemble_vxc

! /DASSEMBLERS
        use M_Dassemble_rho_McWEDA
        use M_Dassemble_2c
        use M_Dassemble_vxc_3c

! Type Declaration
! ===========================================================================

! module procedures
        contains

! ===========================================================================
! Dassemble_vxc.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>  This is the main module for assembling the derivatires of the Vxc matrix elements
!! interactions (McWEDA). The results are stored in vxc (potential) Dblok part.
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
        subroutine Dassemble_vxc (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                      !< counter over atoms
        integer logfile                    !< writing to which unit
        integer num_neigh                  !< number of neighbors

        type(T_forces), pointer :: pfi

! Allocate Arrays
! ===========================================================================
! We build the vna_ontop forces here, so we allocate and initialize
! The vna_atom forces are built in M_build_forces because they are assembled
! differently than the vna_ontop forces.
        do iatom = 1, s%natoms
          pfi=>s%forces(iatom)
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pfi%vxc_on_site (3, num_neigh)); pfi%vxc_on_site = 0.0d0
          allocate (pfi%vxc_off_site (3, num_neigh)); pfi%vxc_off_site = 0.0d0
        end do

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Calculate rho_in (density) matrix elements
        write (logfile,*) ' Calling rho (density) input assemblers. '

! calculate  derivatires for average_rho matrix elements
! See Review OSLXC method
        call Dassemble_rho_2c (s)
        call Dassemble_rho_weighted_2c (s)

! calculate derivatives of XC-potential matrix elements
        write (logfile,*) ' Calling vxc Dassemblers. '

! Calculate the derivativies of the matrix elements
        ! this term has the 2 center SNXC part and the OLSXC part too
        call Dassemble_vxc_SN (s)

        ! this term has the d(<mu_i|vxc(rho_ij)|nu_j>)/dR part
        call Dassemble_vxc_bond (s)

! Three-center forces for exchange-correlation
        call Dassemble_vxc_3c (s)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_vxc


! ===========================================================================
! Dassemble_vxc_SN
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates/assembles the SNXC 2 center part, the OLSXC 2
!!       center part, and also the d (< mu_i| V_xc (rho) | nu_i>)/dR called
!!       "on-site" in the notes.
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
        subroutine Dassemble_vxc_SN (s)
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
        integer in1, in2, in3            !< species numbers
        integer jatom                    !< neighbor of iatom
        integer interaction, isorp       !< which interaction and subtype
        integer logfile                  !< writing to which unit
        integer num_neigh                !< number of neighbors
        integer matom                    !< matom is the self-interaction atom
        integer mbeta                    !< the cell containing neighbor of iatom

        integer imu, inu                 !< counter over orbitals
        integer issh, jssh               !< counter over shells
        integer n1, n2, l1, l2, m1, m2   !< quantum numbers n, l, and m
        integer norb_mu, norb_nu         !< size of the block for the pair

        integer nssh_i                   !< size of the block for the pair

        ! inputs for xc functional
        real prho_in_shell               !< temporary storage
        real prho_bond_shell             !< temporary storage
        real poverlap                    !< temporary storage
        real prho_in                     !< temporary storage
        real prho_bond                   !< temporary storage
        real dexc_in, dexc_bond          !< 1st derivative of xc energy
        real d2exc_in, d2exc_bond        !< 2nd derivative of xc energy
        real dmuxc_in, dmuxc_bond        !< 1st derivative of xc potential
        real exc_in, exc_bond            !< xc energy
        real muxc_in, muxc_bond          !< xc potential_
        real d2muxc_in, d2muxc_bond      !< 2nd derivative of xc potential

        real z                           !< distance between r1 and r2
        real Qin

        real, dimension (3) :: eta       !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sighat    !< unit vector along r2 - r1

        ! vector derivatives of rho pieces
        real, dimension (3) :: Dprho_in_shell
        real, dimension (3) :: Dprho_in
        real, dimension (3) :: Dprho_bond_shell
        real, dimension (3) :: Dprho_bond
        real, dimension (3) :: Dpoverlap

! bcxcm = density matrix in molecular coordinates
! bcxcx = density matrix in crystal coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of density matrix in molecular coordinates
! vdbcxcx = vectorized derivative of density matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: dbcxcm
        real, dimension (:, :, :), allocatable :: vdbcxcm

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_block), pointer :: pRho_neighbors_matom

        type(T_forces), pointer :: pfi

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Procedure
! ===========================================================================
! ***************************************************************************
!
!                       O F F  -  S I T E    P I E C E(off diagonal terms)
!
! ***************************************************************************
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          pdenmat=>s%denmat(iatom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)

! If r1 .eq. r2, then this is a case of a self-interaction or "on-site" term;
! therefore, we do not calculate here.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in "on-site" case.
            else

! SPECIAL LOOP: we want to minimize the number of calls to lda-function
! we only need to call lda_ceperley-adler for each pair of shells
! but then we need to calculate the (mu,nu)-block of matrix elements
              n1 = 0
              do issh = 1, species(in1)%nssh
! n1 : counter used to determine orbitals imu
                l1 = species(in1)%shell(issh)%lssh
                n1 = n1 + l1 + 1
! Loop over shells ineigh
                n2 = 0
                do jssh = 1, species(in2)%nssh
! n2 : counter used to determine orbitals inu
                   l2 = species(in2)%shell(jssh)%lssh
                   n2 = n2 + l2 + 1
! Call lda-function with rho_in_weighted to get the coefficients for vxc expancion
                   prho_in_shell =                                           &
                    s%rho_in_weighted(iatom)%neighbors(ineigh)%block(issh,jssh)
                   Dprho_in_shell=                                           &
     &              s%rho_in_weighted(iatom)%neighbors(ineigh)%Dblock(:,issh,jssh)
                   call lda_ceperley_alder (prho_in_shell, exc_in, muxc_in,  &
     &                                      dexc_in, d2exc_in, dmuxc_in, d2muxc_in)

! Call lda-function with rho_bon_weighted to get the coefficients for vxc expancion
                   prho_bond_shell =                                         &
     &              s%rho_bond_weighted(iatom)%neighbors(ineigh)%block(issh,jssh)
                   Dprho_bond_shell=                               &
     &              s%rho_bond_weighted(iatom)%neighbors(ineigh)%Dblock(:,issh,jssh)     
                   call lda_ceperley_alder (prho_bond_shell, exc_bond, muxc_bond, &
     &                                      dexc_bond, d2exc_bond, dmuxc_bond, d2muxc_bond)

! Calculate forces from vxc_SN and vxc_SN_bond for (mu,nu)-block
! loop over orbitals in the iatom-shell (imu)
                   do m1 = -l1, l1
                     imu = n1 + m1
! loop over orbitals in the ineigh-shell (inu)
                     do m2 = -l2, l2
                       inu = n2 + m2
! terms needed to build up the forces for the two-center parts
                       poverlap =                                              &
     &                   s%overlap(iatom)%neighbors(ineigh)%block(imu,inu)
                       Dpoverlap =                                             &
     &                   s%overlap(iatom)%neighbors(ineigh)%Dblock(:,imu,inu)

                       prho_in =                                               &
     &                   s%rho_in(iatom)%neighbors(ineigh)%block(imu,inu)
                       Dprho_in =                                              &
     &                   s%rho_in(iatom)%neighbors(ineigh)%Dblock(:,imu,inu)
                       
                       prho_bond =                                             &
     &                   s%rho_bond(iatom)%neighbors(ineigh)%block(imu,inu)
                       Dprho_bond =                                            &
     &                   s%rho_bond(iatom)%neighbors(ineigh)%Dblock(:,imu,inu)

! calculate GSN force based on rho_in
                       pfi%vxc_off_site(:,ineigh) = pfi%vxc_off_site(:,ineigh)&
     &                   - pRho_neighbors%block(imu,inu)                      &
     &                    *(Dpoverlap*(muxc_in - dmuxc_in*prho_in_shell)      &
     &                      + Dprho_in_shell*d2muxc_in                        &
     &                       *(prho_in - prho_in_shell*poverlap)              &
     &                      + dmuxc_in*Dprho_in)

! calculate GSN for rho_bond ("atomic" correction)
! Use "+" here because the energy contribution for the bond-part is "-"
                       pfi%vxc_off_site(:,ineigh) = pfi%vxc_off_site(:,ineigh)&
     &                   + pRho_neighbors%block(imu,inu)                      &
     &                    *(Dpoverlap*(muxc_bond - dmuxc_bond*prho_bond_shell)&
     &                      + Dprho_bond_shell*d2muxc_bond                    &
     &                       *(prho_bond - prho_bond_shell*poverlap)          &
     &                      + dmuxc_bond*Dprho_bond)
                     end do !** m2 = -l2, l2
                   end do !** m1 = -l1, l1
                end do !** jssh = 1, species(in2)%nssh
              end do !** issh = 1, species(in1)%nssh
            end if !** differenciate between on site and off site
          end do !** over the neighbors
        end do !** over the atoms

! ***************************************************************************
!
!                       O N  -  S I T E    P I E C E (diagonal terms)
!
! ***************************************************************************
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          nssh_i = species(in1)%nssh
          matom = s%neigh_self(iatom)

          ! cut some lengthy notation
          pdenmat=>s%denmat(iatom)
          pRho_neighbors_matom=>pdenmat%neighbors(matom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Calculate the distance between the two centers.
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

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! If r1 .eq. r2, then this is a case of a self-interaction or "on-site" term;
! therefore, we do not calculate here.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in "on-site" case.
            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). No rotations here
              interaction = P_rhoS_atom
              in3 = in1

! bcxcm = density matrix in molecular coordinates
! dbcxcm = derivative of density matrix in molecular coordinates
! vdbcxcm = vectorized derivative of denstiy matrix in molecular coordinates
              allocate (bcxcm (nssh_i, nssh_i)); bcxcm = 0.0d0
              allocate (dbcxcm (nssh_i, nssh_i)); dbcxcm = 0.0d0
              allocate (vdbcxcm (3, nssh_i, nssh_i)); vdbcxcm = 0.0d0

              do isorp = 1, species(in2)%nssh
                Qin = s%atom(jatom)%shell(isorp)%Qin
                call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,       &
     &                                 nssh_i, nssh_i, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, nssh_i ! norb_nu
                  do imu = 1, nssh_i ! norb_mu
                    if (z .gt. 1.0d-3) then
                      vdbcxcm(:,imu,inu) = vdbcxcm(:,imu,inu)                 &
     &                                    - eta(:)*dbcxcm(imu,inu)*Qin
                    end if
                  end do
                end do
              end do

! SPECIAL LOOP: we want to minimize the number of calls to lda-function
! we only need to call lda_ceperley-adler for each pair of shells
! but then we need to calculate the (mu,nu)-block of matrix elements
! Loop over shells i-atom
              n1 = 0
              do issh = 1, species(in1)%nssh
! n1 : counter used to determine orbitals imu
                l1 = species(in1)%shell(issh)%lssh
                n1 = n1 + l1 + 1

! Call lda-function for rho_in
                prho_in_shell =                                                &
     &           s%rho_in_weighted(iatom)%neighbors(matom)%block(issh,issh)
                Dprho_in_shell = vdbcxcm(:,issh,issh)
                call lda_ceperley_alder (prho_in_shell, exc_in, muxc_in,       &
     &                                   dexc_in, d2exc_in, dmuxc_in, d2muxc_in)

!               prho_bond_shell =                                              &
!    &           s%rho_bond_weighted(iatom)%neighbors(matom)%block(issh,issh)
!               Dprho_bond_shell=                                              &
!    &           s%rho_bond_weighted(iatom)%neighbors(matom)%Dblock(:,issh,issh)
!               call lda_ceperley_alder (prho_bond_shell, exc_bond, muxc_bond, &
!    &                                   dexc_bond, d2exc_bond, dmuxc_bond,    &
!    &                                   d2muxc_bond)

! Calculate vxc_SN and vxc_SN_bond for (mu,nu)-block
! loop over orbitals in the iatom-shell (imu)
                do m1 = -l1, l1
                  imu = n1 + m1
                  prho_in = s%rho_in(iatom)%neighbors(matom)%block(imu,imu)
                  Dprho_in = s%rho_in(iatom)%neighbors(matom)%Dblock(:,imu,imu)
     
!                 prho_bond = s%rho_bond(iatom)%neighbors(matom)%block(imu,imu)
!                 Dprho_bond = s%rho_bond(iatom)%neighbors(matom)%Dblock(:,imu,imu)

! calculate GSN for rho_in
                  pfi%vxc_on_site(:,ineigh) = pfi%vxc_on_site(:,ineigh)        &
     &              - pRho_neighbors_matom%block(imu,imu)                      &
     &               *(dmuxc_in*Dprho_in_shell                                 &
     &                 + Dprho_in_shell*d2muxc_in*(prho_in - prho_in_shell))

! calculate GSN for rho_bond ("atomic" correction)
! Use "+" here because the energy contribution for the bond-part is "-"
!                 pfi%vxc_on_site(:,ineigh)= pfi%vxc_on_site(:,ineigh)         &
!    &              + pRho_neighbors_matom%block(imu,imu)                      &
!    &               *(dmuxc_bond*Dprho_bond_shell                             &
!    &                 + Dprho_bond_shell*d2muxc_bond*(prho_bond - prho_bond_shell))
                end do ! m1 = -l1, l1

                n1 = n1 + l1
              end do  ! do issh = 1, nssh(in1)
              deallocate (bcxcm, dbcxcm, vdbcxcm)
            end if  ! end if iatom .eq. jatom
          end do  ! loop over the neighbors
        end do  ! loop over the atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_vxc_SN


! ===========================================================================
! Dassemble_vxc_bond.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> Dassemble_vxc_bond : assembles the derivative XC-potential matrix elements
!!
!!         d (vxc_bond (mu,nu)/dR =  d (< mu | V_xc(rho_i + rho_j) | nu >)/dR
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
        subroutine Dassemble_vxc_bond (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh             !< counter over atoms and neighbors
        integer in1, in2, in3             !< species numbers
        integer jatom                     !< neighbor of iatom
        integer interaction, isorp        !< which interaction and subtype
        integer num_neigh                 !< number of neighbors
        integer mbeta                     !< the cell containing iatom's neighbor

        integer imu, inu
        integer norb_mu, norb_nu          !< size of the block for the pair

        real z                            !< distance between r1 and r2
        real dQ                           !< charge in orbital

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3, 3) :: eps     !< the epsilon matrix
        real, dimension (3) :: r1, r2     !< positions of iatom and jatom
        real, dimension (3) :: sighat     !< unit vector along r2 - r1

        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: dbcxcm
        real, dimension (:, :, :), allocatable :: vdcxcm
        real, dimension (:, :, :), allocatable :: vdbcxcx

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors

        type(T_forces), pointer :: pfi

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          pdenmat=>s%denmat(iatom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some more lengthy notation
            pRho_neighbors=>pdenmat%neighbors(ineigh)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
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
! call epsilon (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! CALL DOSCENTROS AND GET rho_in FOR ONTOP CASE (i.e. OFF-SITE matrix elements)
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - this is the on-site case

            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcnam). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in vbcnax, the vectorized matrix
! elements; x means crytal coordinates.

! FORCES - VXC ONTOP CASE
! ****************************************************************************
! For the vxc_ontop case, the potential is in the first atom - left (iatom):
! dbcnam is the "scalar" derivative of the matrix; vdbcnam is the "vector"
! derivative of the matrix in molecular coordinates.  When we are done, we get:
! vdtx as the vector derivative of the matrix in crystal coordinates.
              interaction = P_vxc_ontop
              in3 = in2

! Allocate array blocks
              allocate (bcxcm (norb_mu, norb_nu)); bcxcm = 0.0d0
              allocate (dbcxcm (norb_mu, norb_nu)); dbcxcm = 0.0d0
              allocate (vdcxcm (3, norb_mu, norb_nu)); vdcxcm = 0.0d0
              allocate (vdbcxcx (3, norb_mu, norb_nu)); vdbcxcx = 0.0d0

! Neutral atom case
              isorp = 0
              call getDMEs_Fdata_2c (in1, in3, interaction, isorp, z,          &
     &                               norb_mu, norb_nu, bcxcm, dbcxcm)

! dtm is the "scalar" derivative of the matrix; vdtm is the "vector" derivative
! of the matrix in molecular coordinates.  When we are done, we get: vdtx as
! the vector derivative of the matrix in crystal coordinates.

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  if (z .gt. 1.0d-3) vdcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                end do
              end do

! Drotate then puts the vectors in coordinates alone the bond-charge.
              call Drotate (in1, in3, eps, deps, norb_mu, norb_nu,            &
     &                      bcxcm, vdcxcm, vdbcxcx)
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfi%vxc_off_site(:,ineigh) = pfi%vxc_off_site(:,ineigh)     &
     &              - pRho_neighbors%block(imu,inu)*vdbcxcx(:,imu,inu)
                end do
              end do

! FORCES - DNUXC ONTOP LEFT AND RIGHT CASES
! ****************************************************************************
! For the dnuxc_ontop case, the potential is in the first atom - left (iatom):
! dbcnam is the "scalar" derivative of the matrix; vdbcnam is the "vector"
! derivative of the matrix in molecular coordinates.  When we are done, we get:
! vdtx as the vector derivative of the matrix in crystal coordinates.

! Charged atom cases
              interaction = P_dnuxc_ontopL
              in3 = in2

              do isorp = 1, species(in1)%nssh
                call getDMEs_Fdata_2c (in1, in3, interaction, isorp, z,        &
     &                                 norb_mu, norb_nu, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

! Drotate then puts the vectors in coordinates alone the bond-charge.
                call Drotate (in1, in3, eps, deps, norb_mu, norb_nu,           &
     &                        bcxcm, vdcxcm, vdbcxcx)

                dQ = s%atom(iatom)%shell(isorp)%dQ
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfi%vxc_off_site(:,ineigh) = pfi%vxc_off_site(:,ineigh)    &
     &                - pRho_neighbors%block(imu,inu)*vdbcxcx(:,imu,inu)*dQ
                  end do
                end do
              end do  ! end loop over isorp

! Charged atom cases
              interaction = P_dnuxc_ontopR
              in3 = in2

              do isorp = 1, species(in2)%nssh
                call getDMEs_Fdata_2c (in1, in3, interaction, isorp, z,        &
     &                                 norb_mu, norb_nu, bcxcm, dbcxcm)

! Note the minus sign. d/dr1 = - eta * d/dd.
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    if (z .gt. 1.0d-3) vdcxcm(:,imu,inu) = - eta(:)*dbcxcm(imu,inu)
                  end do
                end do

! Drotate then puts the vectors in coordinates alone the bond-charge.
                call Drotate (in1, in3, eps, deps, norb_mu, norb_nu,           &
     &                        bcxcm, vdcxcm, vdbcxcx)

                dQ = s%atom(jatom)%shell(isorp)%dQ
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfi%vxc_off_site(:,ineigh) = pfi%vxc_off_site(:,ineigh)    &
     &                - pRho_neighbors%block(imu,inu)*vdbcxcx(:,imu,inu)*dQ
                  end do
                end do
              end do  ! end loop over isorp
              deallocate (bcxcm, dbcxcm, vdcxcm, vdbcxcx)
            end if ! different atoms condition
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
        end subroutine Dassemble_vxc_bond


! End Module
! ===========================================================================
        end module M_Dassemble_vxc
