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
!!       assemble_vxc_SN : assembles the generalized Sankey-Niklewski both for
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
!!       assemble_vxc_bond : assembles the XC-potential matrix elements
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
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_1c
        use M_Fdata_2c
        use M_Fdata_3c
        use M_assemble_rho_McWEDA
        use M_rotations

! Type Declaration
! ===========================================================================
! two-center interactions
! Calculate SN-exchange-correlation potential using input from
! M_assemble_McWEDA_rho
        type(T_assemble_neighbors), pointer :: vxc_bond (:)
        type(T_assemble_neighbors), pointer :: vxc_SN (:)
        type(T_assemble_neighbors), pointer :: vxc_SN_bond (:)

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
        integer logfile                    !< writing to which unit
        integer num_neigh                !< number of neighbors

        integer norb_mu, norb_nu         !< size of the block for the pair

        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc

! Allocate Arrays
! ===========================================================================
        allocate (s%vxc(s%natoms))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Calculate rho_in (density) matrix elements
        write (logfile,*) ' Calling rho (density) input assemblers. '
        call assemble_rho_2c (s)
        call assemble_rho_3c (s)

! calculate average_rho matrix elements
! See PRB 71, 235101 (2005), Eqs. (19), (22) and (25)
        call assemble_rho_weighted_2c (s)
        call assemble_rho_weighted_3c (s)

! calculate  XC-potential matrix elements
! See PRB 71, 235101 (2005), Eqs. (16), (21) and (24)
        write (logfile,*) ' Calling vxc assemblers. '

! (1) calculate Sankey-Niklewski XC-potential matrix elements
! See PRB 71, 235101 (2005), Eq. (20),
! for the second and third terms on the right of Eq. (16)
        call assemble_vxc_SN (s)

! (2) calculate the "atomic" (first term on the right of Eq. (16))
! matrix elements
        call assemble_vxc_bond (s)

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
            pvxc_neighbors%block = vxc_SN(iatom)%neighbors(ineigh)%block     &
     &                            + vxc_bond(iatom)%neighbors(ineigh)%block  &
     &                            - vxc_SN_bond(iatom)%neighbors(ineigh)%block
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
        subroutine assemble_vxc_SN (s)
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

        integer imu, inu                 !< counter over orbitals
        integer issh, jssh               !< counter over shells
        integer n1, n2, l1, l2, m1, m2   !< quantum numbers n, l, and m
        integer norb_mu, norb_nu         !< size of the block for the pair

! Ceperley-Adler
        real prho_in_shell               !< temporary storage
        real prho_bond_shell             !< temporary storage
        real poverlap                    !< temporary stroage
        real prho_in                     !< temporary storage
        real prho_bond                   !< temporary storage
        real dexc_in, dexc_bond          !< 1st derivative of xc energy
        real d2exc_in, d2exc_bond        !< 2nd derivative of xc energy
        real dmuxc_in, dmuxc_bond        !< 1st derivative of xc potential
        real exc_in, exc_bond            !< xc energy
        real muxc_in, muxc_bond          !< xc potential_
        real d2muxc_in, d2muxc_bond      !< 2nd derivative of xc potential

        type(T_assemble_block), pointer :: pvxc_SN_neighbors
        type(T_assemble_neighbors), pointer :: pvxc_SN

        type(T_assemble_block), pointer :: pvxc_SN_bond_neighbors
        type(T_assemble_neighbors), pointer :: pvxc_SN_bond

! Allocate Arrays
! ===========================================================================
        allocate (vxc_SN (s%natoms))
        allocate (vxc_SN_bond (s%natoms))

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pvxc_SN=>vxc_SN(iatom)
          pvxc_SN_bond=>vxc_SN_bond(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pvxc_SN%neighbors(num_neigh))
          allocate (pvxc_SN_bond%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pvxc_SN_neighbors=>pvxc_SN%neighbors(ineigh)
            pvxc_SN_bond_neighbors=>pvxc_SN_bond%neighbors(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pvxc_SN_neighbors%block(norb_mu, norb_nu))
            pvxc_SN_neighbors%block = 0.0d0
            allocate (pvxc_SN_bond_neighbors%block(norb_mu, norb_nu))
            pvxc_SN_bond_neighbors%block = 0.0d0

! SPECIAL LOOP: we want to minimize the number of calls to lda-function
! we only need to call lda_ceperley-adler for each pair of shells
! but then we need to calculate the (mu,nu)-block of matrix elements

! Loop over shells i-atom
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

! Call lda-function for rho_in
                prho_in_shell =                                              &
     &           s%rho_in_weighted(iatom)%neighbors(ineigh)%block(issh,jssh)
                call lda_ceperley_alder (prho_in_shell, exc_in, muxc_in,     &
     &                                   dexc_in, d2exc_in, dmuxc_in, d2muxc_in)

! Call lda-function for rho_bond
                prho_bond_shell =                                            &
     &           s%rho_bond_weighted(iatom)%neighbors(ineigh)%block(issh,jssh)
                call lda_ceperley_alder (prho_bond_shell, exc_bond,          &
     &                                   muxc_bond, dexc_bond, d2exc_bond,   &
     &                                   dmuxc_bond, d2muxc_bond)

! Calculate vxc_SN and vxc_SN_bond for (mu,nu)-block
! loop over orbitals in the iatom-shell (imu)
                do m1 = -l1, l1
                  imu = n1 + m1
! loop over orbitals in the ineigh-shell (inu)
                  do m2 = -l2, l2
                    inu = n2 + m2
                    poverlap = s%overlap(iatom)%neighbors(ineigh)%block(imu,inu)
                    prho_in = s%rho_in(iatom)%neighbors(ineigh)%block(imu,inu)
                    prho_bond = s%rho_bond(iatom)%neighbors(ineigh)%block(imu,inu)

! calculate GSN for rho_in
                    pvxc_SN_neighbors%block(imu,inu) = muxc_in*poverlap      &
     &                + dmuxc_in*(prho_in - prho_in_shell*poverlap)

! calculate GSN for rho_bond ("atomic" correction)
                    pvxc_SN_bond_neighbors%block(imu,inu) =                  &
     &                + muxc_bond*poverlap                                   &
     &                + dmuxc_bond*(prho_bond - prho_bond_shell*poverlap)

                  end do !do m2 = -l2, l2
                end do !do m1 = -l1, l1

                n2 = n2 + l2
              end do  ! do jssh = 1, nssh(in1)

              n1 = n1 + l1
            end do  ! do issh = 1, nssh(in1)

! end of SPECIAL LOOP
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_vxc_SN


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
! Code written by:
!> @author Daniel G. Trabada
!! @author Jose Ortega (JOM)
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_vxc_bond (s)
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
        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor

        integer imu, inu
        integer issh, jssh, kssh
        integer norb_mu, norb_nu         !< size of the block for the pair

        real z                           !< distance between r1 and r2
        real dQ

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sighat    !< unit vector along r2 - r1

        real, dimension (:, :), allocatable :: bcxcm
        real, dimension (:, :), allocatable :: bcxcx

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pvxc_bond_neighbors
        type(T_assemble_neighbors), pointer :: pvxc_bond

! Allocate Arrays
! ===========================================================================
        allocate (vxc_bond(s%natoms))

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pvxc_bond=>vxc_bond(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pvxc_bond%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pvxc_bond_neighbors=>pvxc_bond%neighbors(ineigh)
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pvxc_bond_neighbors%block(norb_mu, norb_nu))
            pvxc_bond_neighbors%block = 0.0d0

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

! CALL DOSCENTROS AND GET rho_in FOR ONTOP CASE (i.e. OFF-SITE matrix elements)
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! We calculate here the "on-site" term: mu and nu are in the same atom "i"
!    vxc_bond (mu,nu) = < mu | V_xc (rho_i) | nu >

! transform from l-block to mu-block and get the matrix element from the data files
! Fdata_1c (onecenter:uncentro)

! this loop defining dQ maybe outside this subroutine
              do inu = 1, norb_nu
                jssh = species(in1)%orbital(inu)%issh
                do imu = 1, norb_mu
                  issh = species(in1)%orbital(imu)%issh
                  if (species(in1)%orbital(imu)%l                            &
     &                .eq. species(in1)%orbital(inu)%l .and.                 &
     &                species(in1)%orbital(imu)%m                            &
     &                .eq. species(in1)%orbital(inu)%m) then
                    pvxc_bond_neighbors%block(imu,inu) = vxc_1c(in1)%V(issh,jssh)
                    do kssh = 1, species(in1)%nssh
                      dQ = s%atom(iatom)%shell(kssh)%dQ
                      pvxc_bond_neighbors%block(imu,inu) =                   &
     &                  pvxc_bond_neighbors%block(imu,inu)                   &
     &                    + vxc_1c(in1)%dV(issh,jssh,kssh)*dQ
                    end do
                  end if
                end do
              end do

            else

! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcxcm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in bcxcx, where x means crytal
! coordinates.
              interaction = P_vxc_ontop
              in3 = in2

! Allocate array blocks
              allocate (bcxcm (norb_mu, norb_nu))
              allocate (bcxcx (norb_mu, norb_nu))

! Neutral atom case
              isorp = 0
              call getMEs_Fdata_2c (in1, in3, interaction, isorp, z,         &
     &                              norb_mu, norb_nu, bcxcm)
              call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)
              pvxc_bond_neighbors%block = bcxcx

! Charged atom cases
              interaction = P_dnuxc_ontopL
              in3 = in2
              do isorp = 1, species(in1)%nssh
                call getMEs_Fdata_2c (in1, in3, interaction, isorp, z,       &
     &                                norb_mu, norb_nu, bcxcm)
                call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)
                dQ = s%atom(iatom)%shell(isorp)%dQ
                pvxc_bond_neighbors%block = pvxc_bond_neighbors%block + bcxcx*dQ
              end do

              interaction = P_dnuxc_ontopR
              in3 = in2
              do isorp = 1, species(in2)%nssh
                call getMEs_Fdata_2c (in1, in3, interaction, isorp, z,       &
     &                                norb_mu, norb_nu, bcxcm)
                call rotate (in1, in3, eps, norb_mu, norb_nu, bcxcm, bcxcx)
                dQ = s%atom(jatom)%shell(isorp)%dQ
                pvxc_bond_neighbors%block = pvxc_bond_neighbors%block + bcxcx*dQ
              end do
              deallocate (bcxcm)
              deallocate (bcxcx)

            end if ! end if for r1 .eq. r2 case

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
        end subroutine assemble_vxc_bond


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
            deallocate (vxc_SN(iatom)%neighbors(ineigh)%block)
            deallocate (vxc_SN_bond(iatom)%neighbors(ineigh)%block)
            deallocate (vxc_bond(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%vxc(iatom)%neighbors)
          deallocate (vxc_SN(iatom)%neighbors)
          deallocate (vxc_SN_bond(iatom)%neighbors)
          deallocate (vxc_bond(iatom)%neighbors)
        end do
        deallocate (s%vxc)
        deallocate (vxc_SN)
        deallocate (vxc_SN_bond)
        deallocate (vxc_bond)

        call destroy_assemble_rho (s)

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

