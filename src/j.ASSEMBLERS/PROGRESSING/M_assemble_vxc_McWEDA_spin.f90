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
!>       This is a module containing all of the assembler programs required
!! to assemble all of the matrix elements for exchange-correlation
!! It contains the following subroutines within the module:
!!       assemble_vxc_SN - assembles the Generalized-Sankey-Niklewski
!!       ...
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================

! OLD STUFF!!!!!!!!!! Need to be deleted when we're sure that we'll not use it!!!

if (.false.) then !! for not using this
        module M_assemble_vxc
        use M_assemble_McWEDA_rho
        use M_assemble_McWEDA_weight
        use M_assemble_2c
        use M_assemble_blocks
        use M_cells
        use M_configuraciones
        use M_Fdata_2c
        use M_Fdata_3c
        use M_neighbors
        use M_rotations
        use M_species

! Type Declaration
! ===========================================================================
! Calculate SN-exchange-correlation potential using inpurt from M_assemble_McWEDA_rho
! and M_assemble_McWEDA_weight
! Put all the neighbor group belonging to the atom
! so in the end we have something like overlap(mu, nu, ineigh, iatom)
!        type(T_assemble_neighbors), pointer :: vxc_sn (:) ! This should be out?
        ! Maybe this definition should be done in another place
        type T_assemble_spin
          type(T_assemble_neighbors), pointer :: vxc_sn_spin (:)
          type(T_assemble_neighbors), pointer :: rho_in_spin (:)
          type(T_assemble_neighbors), pointer :: rho_at_spin (:)
          type(T_assemble_neighbors), pointer :: w_rho_in_spin (:)
          type(T_assemble_neighbors), pointer :: w_rho_at_spin (:)
        end type T_assemble_spin

        type(T_assemble_spin), pointer :: xc_spin_sn(:)

! module procedures
        contains

! ===========================================================================
! assemble_vxc_SN_spin.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates exchange-correlation potential
!!       matrix interactions.
!
! ===========================================================================
! Code written by:
!> @author Enrique Abad Gonzalez
! Dpto. de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! Phone: +34 91 497 86 48
! Fax: +34 91 497 49 50
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_vxc_spin_SN (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

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
        integer interaction, isubtype    !< which interaction and subtype
        integer num_neigh                !< number of neighbors
        integer matom                    !< matom is the self-interaction atom
        integer mbeta                    !< the cell containing neighbor of iatom

        integer imu, inu
        integer norb_mu, norb_nu         !< size of the block for the pair




! Ceperley-Adler
        real, dimension(2) :: arho
        real ovlap
        real, dimension(2) :: rhoin
        real, dimension(2) :: dexc
        real, dimension(3) :: d2exc
        real, dimension(4) :: dmuxc
        real exc
        real, dimension(2) :: muxc
        real, dimension(6) :: d2muxc

        type(T_assemble_block), dimension(2), pointer :: pvxc_sn_neighbors
        type(T_assemble_neighbors), dimension(2), pointer :: pvxc_sn


! Allocate Arrays
! ===========================================================================
!        allocate (s%vxc_sn(s%natoms))
        allocate (xc_spin_sn(2))

! Procedure
! ===========================================================================

      write (logfile,*) ' Calling  rho assemblers: '
      do ispin = 1, 2
        s%rho_in => xc_spin_sn(ispin)%rho_in_spin
        s%rho_at => xc_spin_sn(ispin)%rho_at_spin
        call assemble_rho_2c(Qin)
        call assemble_rho_3c(Qin)
      end do

      write (logfile,*) ' Calling  w_rho assemblers: '
      do ispin = 1, 2
        w_rho_in => xc_spin_sn(ispin)%w_rho_in_spin
        w_rho_at => xc_spin_sn(ispin)%w_rho_at_spin
        call assemble_w_rho_2c(Qin)
        call assemble_w_rho_3c(Qin)
      end do

      do ispin = 1, 2
        av_rho_at_shell => xc_spin_sn(ispin)%av_rho_at_shell_spin
        av_rho_in_shell => xc_spin_sn(ispin)%av_rho_in_shell_spin
        call assemble_av_rho()
      end do


    write (logfile,*) ' Calling  vxc assemblers: '








! Loop over the atoms in the central cell.
!      pxc_spin_sn=> xc_spin_sn(ispin)
      allocate(xc_spin_sn(ispin)%vxc_sn(s%natoms))
        do iatom = 1, s%natoms
          ! cut some lengthy notation
!          pvxc_sn=>pvxc_spin_sn%vxc_sn(iatom)
          do ispin = 1, 2
            pvxc_sn(ispin)=>xc_spin_sn(ispin)%vxc_sn(iatom)
          end do
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = neighbors(iatom)%neighn
          do ispin = 1, 2
            allocate (pvxc_sn(ispin)%neighbors(num_neigh))
          end do

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            do ispin = 1,2
              pvxc_sn_neighbors(ispin)=>pvxc_sn(ispin)%neighbors(ineigh)
            end do
            jatom = neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
 write (logfile,*) ' Matrices vxc SN',iatom, ineigh, jatom, neighbors(iatom)%neigh_b(ineigh)
 write (logfile,*)'            '
! Allocate block size
            norb_nu = species(in2)%norb_max
            do ispin = 1, 2
              allocate (pvxc_sn_neighbors(ispin)%block(norb_mu, norb_nu))
              pvxc_sn_neighbors(ispin)%block = 0.0d0
            end do






! Loop over shells i-atom
            n1 = 0
            do issh = 1, species(in1)%nssh

! Number of orbitals per the shell in x-dim
              l1 = species(in1)%shell(issh)%lssh
              n1 = n1 + l1 + 1

! Loop over shells ineigh
              n2 = 0
              do jssh = 1, species(in2)%nssh
                l2 = species(in2)%shell(jssh)%lssh
! Number of orbitals per the shell
                n2 = n2 + l2 + 1
                do ispin = 1, 2
                  arho(ispin)=xc_spin_sn(ispin)%av_rho_in_shell_spin(iatom)%neighbors(ineigh)%block(issh,jssh)
                end do
                call lsdavwn (arho, ex, ec, xpot, cpot, dxpot, dcpot)

! Set the XC-submatrices
! loop over orbitals in the iatom-shell
                do ind1 = -l1, l1
                  imu = n1 + ind1
! loop over orbitals in the ineigh-shell
                  do ind2 = -l2, l2
                    inu = n2 + ind2

                    ovlap = s%overlap(iatom)%neighbors(ineigh)%block(imu,inu)
                    do ispin = 1, 2
                      rhoin(ispin) = xc_spin_sn(ispin)%rho_in_spin(iatom)%neighbors(ineigh)%block(imu,inu)
                    end do
                    do ispin = 1, 2
                      pvxc_sn_neighbors(ispin)%block(imu,inu) = (xpot(ispin)+cpot(ispin))*ovlap + &
         & (dxpot(ispin,1)+dcpot(ispin,1))*(rhoin(1)-arho(1)*ovlap) +(dxpot(ispin,2)+dcpot(ispin,2))*(rhoin(2)-arho(2)*ovlap)
                     end do


                  end do !do ind2 = -l2, l2
                end do !do ind1 = -l1, l1
          n2 = n2 + l2
         end do !do jssh = 1, nssh(in1)
         n1 = n1 + l1
        end do !do issh = 1, nssh(in1)
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
        end subroutine assemble_vxc_sn




        ! lda_ceperley_alder
! Module Description
! ===========================================================================
!>       Given the input density, this subroutine returns the LDA exchange-
!! correlation energies and potentials. This LDA is according to the
!! parameterization of Ceperley-Alder.
!!
!! "Ground state of the electron gas by a stochastic method," by D.M. Ceperley
!!  and B.J. Alder, Phys. Rev. Let. 45:566-569 (1980).
!!
!! "Self-interaction correction to density-functional approximations for
!! many-electron systems," by J.P. Perdew and A. Zunger, Phys. Rev. B.
!! 23:5048-5079 (1981).
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
        subroutine lda_ceperley_alder (rh, exc, muxc, dexc, d2exc, dmuxc,    &
     &                                 d2muxc)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: rh      !< FIXME

! Output
        real, intent (out) :: dexc   !< FIXME
        real, intent (out) :: d2exc  !< FIXME
        real, intent (out) :: dmuxc  !< FIXME
        real, intent (out) :: exc    !< FIXME
        real, intent (out) :: muxc   !< FIXME
        real, intent (out) :: d2muxc !< FIXME

! Parameters and Data Declaration
! ===========================================================================
        real, parameter :: eps = 1.0d-3 !< FIXME

! Variable Declaration and Description
! ===========================================================================
        real d2nec
        real d2nex
        real d3nec
        real d3nex
        real dec
        real ddec
        real d2dec
        real den
        real dden
        real d2den
        real d3den
        real ex
        real hartree
        real rho_third
        real rho
        real rs
        real rsl
        real sqrs

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize things to zero (thero in Spanish).
        exc = 0.0d0
        muxc = 0.0d0
        dexc = 0.0d0
        dmuxc = 0.0d0
        d2muxc = 0.0d0
        d2exc = 0.0d0

        if (rh .le. eps) return

! Convert to a.u.
        rho = rh*(P_abohr**3)

! Find rho^(1/3)
        rho_third = rho**(1.0e0/3.0e0)

! Effective radius
        rs = 0.62035049d0/rho_third

! Find the energy, potential, and the deerivative of the potential.
        if (rho .lt. 0.23873241d0) then
          sqrs = sqrt(rs)
          den = 1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs       ! Effective density
          exc = -0.4581652d0/rs - 0.1423d0/den
          muxc = exc - rs*(0.15273333d0/rs**2                                &
     &               + (0.02497128d0/sqrs + 0.01581427d0)/den**2)

! Stuff for dmuxc
          dden = 1.0529d0/(2.0d0*sqrs) + 0.3334d0
          d2den = (-0.5d0)*1.0529d0/(2.0d0*rs*sqrs)
          d3den = (0.75d0)*1.0529d0/(2.0d0*rs**2*sqrs)
          dec = 0.1423d0*dden/(den**2)
          ddec = -2.0d0*0.1423d0*dden**2/(den**3) + 0.1423d0*d2den/(den**2)
          d2dec = 6.0d0*0.1423d0*(dden*3)/(den**4)                           &
     &           - 6.0d0*0.1423d0*dden*d2den/(den**3) + 0.1423d0*d3den/(den**2)

        else
          rsl = log(rs)
          exc = -0.4581652d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs      &
     &          + 0.002d0*rs*rsl
          muxc = exc - rs*(0.15273333d0/rs**2 + 0.01036667d0/rs              &
     &               - 0.003866667d0 + 0.00066667d0*(1.0d0 + rsl))

          dec = 0.0311d0/rs - 0.0116d0 + 0.0020d0*(rsl + 1.0d0)
          ddec = -0.0311d0/(rs**2) + 0.0020d0/rs
          d2dec = 2.0d0*0.0311d0/(rs**3) - 0.0020d0/(rs**2)
        end if

! Exchange-only energy and potential
        ex = -0.7385587664d0*rho_third

! Now find dmuxc
        dexc = (muxc - exc)/rho

! Compute dmu/dn. Use d(mu)/dn = 2*dexc/dn + n*d2(exc)/dn2.
        d2nec = (4.0d0*rs/(9.0d0*rho**2))*dec + (rs**2/(9.0d0*rho**2))*ddec
        d2nex = -(2.0d0/(9.0d0*rho**2))*ex
        dmuxc = 2.0d0*dexc + rho*(d2nex + d2nec)

! Compute d2mu/dn2, using d2(mu)/dn2 = 3*d2(exc)/dn2 + n*d3(exc)/dn3
        d3nec = (-28.0d0*rs/(27.0d0*rho**3))*dec + (-4.0d0*rs**2/       &
     &   (9.0d0*rho**3))*ddec + (rs**3/(-27.0d0*rho**3))*d2dec
        d3nex = (10.0d0/(27.0d0*rho**3))*ex
        d2muxc = 3.0*(d2nex + d2nec) + rho*(d3nex + d3nec)
        d2exc = d2nex + d2nec

! Convert output to eV (exc and muxc) and to eV*(Angstrom**3) (dexc and dmuxc)
        hartree = P_eq2/P_abohr
        exc = exc*hartree
        muxc = muxc*hartree
        dexc = dexc*hartree*(P_abohr)**3
        d2exc = d2exc*hartree*(P_abohr)**6
        dmuxc = dmuxc*hartree*(P_abohr)**3
        d2muxc = d2muxc*hartree*(P_abohr)**6

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine lda_ceperley_alder

end if        ! end if false, for not using old stuff


        module M_assemble_vxc
        use M_assemble_rho_McWEDA_spin
        use M_assemble_2c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_1c
        use M_Fdata_2c
        use M_Fdata_3c
        use M_neighbors
        use M_rotations
        use M_species


! Type Declaration
! ===========================================================================
! two-center interactions
! Calculate SN-exchange-correlation potential using input from
! M_assemble_McWEDA_rho

! Output


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
! modifications for spin-dependent formalism by: Enrique Abad
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

		integer ispin
        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc

! Allocate Arrays
! ===========================================================================
        do ispin = 1, 2
          allocate (s%spinstuff(ispin)%vxc(s%natoms))
        end do

! Procedure
! ===========================================================================
! Calculate rho_in (density) matrix elements
        write (logfile,*) ' Calling rho (density) input assemblers: '
        do ispin = 1, 2
          call assemble_rho_2c (s,ispin)
          call assemble_rho_3c (s,ispin)
          s%spinstuff(ispin)%rho_in = rho_in
          s%spinstuff(ispin)%rho_bond = rho_bond
        end do

! calculate average_rho matrix elements
! See PRB 71, 235101 (2005), Eqs. (19), (22) and (25)
        call assemble_S_weighted (s)
!        This is not necessary for the SN, isn't it?
!        do ispin = 1, 2
!          call assemble_rho_weighted_2c (s,ispin)
!          call assemble_rho_weighted_3c (s,ispin)
!          call assemble_rho_average (s,ispin)
!        end do

! calculate  XC-potential matrix elements
! See PRB 71, 235101 (2005), Eqs. (16), (21) and (24)
        write (logfile,*) ' Calling  vxc assemblers: '

! (1) calculate Sankey-Niklewski XC-potential matrix elements
! See PRB 71, 235101 (2005), Eq. (20),
! for the second and third terms on the right of Eq. (16)
        call assemble_vxc_SN (s)

! (2) calculate the "atomic" (first term on the right of Eq. (16))
! matrix elements
!        call assemble_vxc_bond (s)

! (3) Sum the 3-contributions in Eq. (16):
!  vxc = vxc_bond + vxc_SN - vxc_SN_bond
        do iatom = 1, s%natoms
         do ispin = 1, 2
          ! cut some lengthy notation
          pvxc=>s%spinstuff(ispin)%vxc(iatom)
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = neighbors(iatom)%neighn
          allocate (pvxc%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            pvxc_neighbors=>pvxc%neighbors(ineigh)
            jatom = neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            allocate (pvxc_neighbors%block(norb_mu, norb_nu))

! ecuacion (16) PRB 71, 235101 (2005)
!  vxc = vxc_bond + vxc_SN - vxc_SN_bond
!            pvxc_neighbors%block = s%spinstuff(ispin)%vxc_SN(iatom)%neighbors(ineigh)%block     &
!     &                            + s%spinstuff(ispin)%vxc_bond(iatom)%neighbors(ineigh)%block  &
!     &                            - s%spinstuff(ispin)%vxc_SN_bond(iatom)%neighbors(ineigh)%block
            pvxc_neighbors%block = s%spinstuff(ispin)%vxc_SN(iatom)%neighbors(ineigh)%block
          end do
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
! modifications for spin-dependent formalism by: Enrique Abad
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
        integer ispin                    !< the spin index

! Ceperley-Adler
        real, dimension(2) :: prho_in_shell        !< temporary storage
        real, dimension(2) :: prho_bond_shell      !< temporary storage
        real, dimension(2) :: prho_in              !< temporary storage
        real, dimension(2) :: prho_bond            !< temporary storage
        real poverlap                    !< temporary stroage
!        real prho_in                     !< temporary storage
!        real prho_bond                   !< temporary storage
        real ex_in, ec_in                !< exchange and correlation energy
        real, dimension(2) :: xpot_in, cpot_in ! xc potential
        real, dimension(2,2) :: dxpot_in, dcpot_in  ! derivative of the xc potential

        type(T_assemble_block), pointer :: pvxc_SN_neighbors_up
        type(T_assemble_block), pointer :: pvxc_SN_neighbors_down
        type(T_assemble_neighbors), pointer :: pvxc_SN_up
        type(T_assemble_neighbors), pointer :: pvxc_SN_down

        type(T_assemble_block), dimension(:), pointer :: pvxc_SN_bond_neighbors
        type(T_assemble_neighbors), dimension(:), pointer :: pvxc_SN_bond

! Allocate Arrays
! ===========================================================================
		do ispin = 1, 2
          allocate (s%spinstuff(ispin)%vxc_SN (s%natoms))
          allocate (s%spinstuff(ispin)%vxc_SN_bond (s%natoms))
        end do


! Procedure
! ===========================================================================

! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = neighbors(iatom)%neighn
           pvxc_SN_up=>s%spinstuff(1)%vxc_SN(iatom)
           pvxc_SN_down=>s%spinstuff(2)%vxc_SN(iatom)
!           pvxc_SN_bond(:)=>spinstuff(:)%vxc_SN_bond(iatom)
           allocate (pvxc_SN_up%neighbors(num_neigh))
           allocate (pvxc_SN_down%neighbors(num_neigh))
!           allocate (pvxc_SN_bond(ispin)%neighbors(num_neigh))


! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            jatom = neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
             pvxc_SN_neighbors_up=>pvxc_SN_up%neighbors(ineigh)
             pvxc_SN_neighbors_down=>pvxc_SN_up%neighbors(ineigh)
!             pvxc_SN_bond_neighbors(ispin)=>pvxc_SN_bond(ispin)%neighbors(ineigh)
             allocate (pvxc_SN_neighbors_up%block(norb_mu, norb_nu))
             allocate (pvxc_SN_neighbors_down%block(norb_mu, norb_nu))
             pvxc_SN_neighbors_up%block = 0.0d0
             pvxc_SN_neighbors_down%block = 0.0d0
!             allocate (pvxc_SN_bond_neighbors%block(norb_mu, norb_nu))
!             pvxc_SN_bond_neighbors(ispin)%block = 0.0d0


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

! Call lsdavwn-function for rho_in
! We only need to call it once for the two spins (that's the reason why the code has not a loop in spins)
                 prho_in_shell(1) =                                              &
     &            s%spinstuff(1)%rho_in_shell(iatom)%neighbors(ineigh)%block(issh,jssh)
                 prho_in_shell(2) =                                              &
     &            s%spinstuff(2)%rho_in_shell(iatom)%neighbors(ineigh)%block(issh,jssh)
                call lsdavwn (prho_in_shell, ex_in, ec_in, xpot_in, cpot_in, dxpot_in, dcpot_in)

! Call lda-function for rho_bond
!                prho_bond_shell =                                            &
!     &           rho_bond_shell(iatom)%neighbors(ineigh)%block(issh,jssh)
!                call lda_ceperley_alder (prho_bond_shell, exc_bond,          &
!     &                                   muxc_bond, dexc_bond, d2exc_bond,   &
!     &                                   dmuxc_bond, d2muxc_bond)
!                call lsdavwn (prho_in_shell, ex, ec, xpot, cpot, dxpot, dcpot)

! Calculate vxc_SN and vxc_SN_bond for (mu,nu)-block
! loop over orbitals in the iatom-shell (imu)
                do m1 = -l1, l1
                  imu = n1 + m1
! loop over orbitals in the ineigh-shell (inu)
                  do m2 = -l2, l2
                    inu = n2 + m2

                    poverlap = overlap(iatom)%neighbors(ineigh)%block(imu,inu)
                    do ispin = 1, 2
                    prho_in(ispin) = s%spinstuff(ispin)%rho_in(iatom)%neighbors(ineigh)%block(imu,inu)
!                    prho_bond(ispin) = spinstuff(ispin)%rho_bond(iatom)%neighbors(ineigh)%block(imu,inu)
					end do

! calculate GSN for rho_in
					ispin = 1
                    pvxc_SN_neighbors_up%block(imu,inu) = (xpot_in(ispin)+cpot_in(ispin))*poverlap +    &
     &        (dxpot_in(ispin,1)+dcpot_in(ispin,1))*(prho_in(1)-prho_in_shell(1)*poverlap) +         &
     &         (dxpot_in(ispin,2)+dcpot_in(ispin,2))*(prho_in(2)-prho_in_shell(2)*poverlap)
			        ispin = 2
                    pvxc_SN_neighbors_down%block(imu,inu) = (xpot_in(ispin)+cpot_in(ispin))*poverlap +    &
     &        (dxpot_in(ispin,1)+dcpot_in(ispin,1))*(prho_in(1)-prho_in_shell(1)*poverlap) +         &
     &         (dxpot_in(ispin,2)+dcpot_in(ispin,2))*(prho_in(2)-prho_in_shell(2)*poverlap)



! calculate GSN for rho_bond ("atomic" correction)
!                    pvxc_SN_bond_neighbors%block(imu,inu) =                  &
!     &                + muxc_bond*poverlap                                   &
!     &                + dmuxc_bond*(prho_bond - prho_bond_shell*poverlap)

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
! Adapted for spin-polarized case by Enrique Abad
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

        do ispin = 1, 2
          do iatom = 1, s%natoms
            do ineigh = 1, neighbors(iatom)%neighn
              deallocate (s%spinstuff(ispin)%vxc(iatom)%neighbors(ineigh)%block)
              deallocate (s%spinstuff(ispin)%vxc_SN(iatom)%neighbors(ineigh)%block)
              deallocate (s%spinstuff(ispin)%vxc_SN_bond(iatom)%neighbors(ineigh)%block)
              deallocate (s%spinstuff(ispin)%vxc_bond(iatom)%neighbors(ineigh)%block)
            end do
            deallocate (s%spinstuff(ispin)%vxc(iatom)%neighbors)
            deallocate (s%spinstuff(ispin)%vxc_SN(iatom)%neighbors)
            deallocate (s%spinstuff(ispin)%vxc_SN_bond(iatom)%neighbors)
            deallocate (s%spinstuff(ispin)%vxc_bond(iatom)%neighbors)
          end do
          deallocate (s%spinstuff(ispin)%vxc)
          deallocate (s%spinstuff(ispin)%vxc_SN)
          deallocate (s%spinstuff(ispin)%vxc_SN_bond)
          deallocate (s%spinstuff(ispin)%vxc_bond)
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
        end subroutine destroy_assemble_vxc_McWEDA




! lsdavwn.f90
! Program Description
! ===========================================================================
!>       This routine computes the exchange and correlation potenials and
!! energies for the Vosko, Wilk, Nusair LSDA functional. Each spin component
!! is considered.
!!
!! See
!!      S.H. VOSKO and L. WILK and M. NUSAIR
!!      Can. J. Phys., 58, 1200 (1980)
!
! ===========================================================================
! Code written by:
!> @author Eduardo Mendez
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! Adapted to Lightning by:
!> @author Enrique Abad
! Departamento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! Phone: +34 91 497 8648
! Fax: +34 497 4950
! ===========================================================================
!
! Program declaration
! ===========================================================================
        subroutine lsdavwn (rh, ex, ec, xpot, cpot, dxpot, dcpot)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        real*8, intent (in), dimension (2) :: rh        !< FIXME

! Output
!        real*8, intent (out) :: dnuxc
!        real*8, intent (out) :: dnuxcs
        real*8, intent (out), dimension (2,2) :: dxpot  !< FIXME
        real*8, intent (out), dimension (2,2) :: dcpot  !< FIXME
        real*8, intent (out) :: ec                      !< FIXME
        real*8, intent (out) :: ex                      !< FIXME

        real*8, intent (out), dimension (2) :: cpot     !< FIXME
        real*8, intent (out), dimension (2) :: xpot     !< FIXME

! Local Parameters and Data Declaration
! ===========================================================================
        real*8, parameter :: Ap = 0.0621814d0  !< FIXME: What are all these magic numbers?
        real*8, parameter :: bp = 3.72744d0
        real*8, parameter :: cp = 12.9352d0
        real*8, parameter :: x0p = -0.10498d0

        real*8, parameter :: Aa = 0.033773728d0
        real*8, parameter :: ba = 1.13107d0
        real*8, parameter :: ca = 13.0045d0
        real*8, parameter :: x0a = -0.00475840d0

        real*8, parameter :: Af = 0.01554535d0
        real*8, parameter :: bf = 7.06042d0
        real*8, parameter :: cf = 18.0578d0
        real*8, parameter :: x0f = -0.32500d0

        real*8, parameter :: epsilon = 1.0d-10
        real*8, parameter :: pi = 3.141592654

! Local Variable Declaration and Description
! ===========================================================================
        real*8 density
        real*8 densitys
        real*8, dimension (2) :: rho

        real*8, dimension (3) :: cdpot
        real*8, dimension (3) :: xdpot

! spin polarization and derivatives
        real*8 zeta, zp1, zp2, zp1p2, zpp1, zpp2
        real*8 x, xp, xpp
        real*8 g, gp, gpp
        real*8 XXp, XXf, XXa , Qp, Qf, Qa, jp, jf, ja
        real*8 ecP, ecF, ecA, ecPp, ecFp, ecAp, ecPpp, ecFpp, ecApp
        real*8 cte, h, hp, hpp
        real*8 ecpx, ecpz, ecppx, ecppz, ecpxpz
        real*8 d1ec, d2ec, dd1ec, dd2ec, d1d2ec
        real*8 exP, exPp, exPpp
        real*8 expd, expz, exppd, exppz, expdpz
        real*8 d1ex, d2ex, dd1ex, dd2ex, d1d2ex
        real*8 F, Fs
! This should be out?
        real*8 dnuxc
        real*8 dnuxcs

! Allocate Arrays
! ===========================================================================

! Procedure
! =========================================================================
! Another change of units
        rho(1) = rh(1)*(P_abohr**3)
        rho(2) = rh(2)*(P_abohr**3)
! Initialize some parameters
        density = rho(1) + rho(2)
        densitys = rho(1) - rho(2)
        zeta = densitys/density
        if (density .le. epsilon) then
         zeta = 0.0d0
         ec = 0.0d0
         ex = 0.0d0
         cpot = 0.0d0
         xpot = 0.0d0
         cdpot = 0.0d0
         xdpot = 0.0d0
         return
        end if

! Define simple derivatives
! *************************************************************************
        zp1 = 2.0d0*rho(2)/density**2
        zp2 = -2.0d0*rho(1)/density**2
        zp1p2 = 2.0d0*(rho(1)*rho(1)+rho(2)*rho(2))/density**4
        zpp1 = -4.0*(rho(1)*rho(2)-rho(2)*rho(2))/density**4
        zpp2 = 4.0*(rho(1)*rho(1)+rho(1)*rho(2))/density**4

        x = (3.0d0/(4.0d0*pi*density))**(1.0d0/6.0d0)
        xp = - (1.0d0/6.0d0)*x/density
        xpp = 7.0d0*x/(36.0d0*density**2)

        g = (9.0d0/8.0d0)*((1.0d0 + zeta)**(4.0d0/3.0d0)                     &
     &     + (1.0d0 - zeta)**(4.0d0/3.0d0) - 2.0d0)
        gp = (3.0d0/2.0d0)*((1.0d0 + zeta)**(1.0d0/3.0d0)                    &
     &      - (1.0d0 - zeta)**(1.0d0/3.0d0))
        gpp = (1.0d0/2.0d0)*((1.0d0 + zeta)**(-2.0d0/3.0d0)                  &
     &       - (1.0d0 - zeta)**(-2.0d0/3.0d0))

! Intermediate variables
        XXp = x**2.0d0 + bp*x + cp
        XXf = x**2.0d0 + bf*x + cf
        XXa = x**2.0d0 + ba*x + ca
        Qp = (4.0d0*cp - bp**2)**0.5d0
        Qf = (4.0d0*cf - bf**2)**0.5d0
        Qa = (4.0d0*ca - ba**2)**0.5d0
        jp = 2.0d0*log(x - x0p) - log(XXp)                                   &
     &      + 2.0d0*((2.0d0*x0p + bp)/Qp)*atan(Qp/(2.0d0*x + bp))
        jf = 2.0d0*log(x - x0f) - log(XXf)                                   &
     &      + 2.0d0*((2.0d0*x0f + bf)/Qf)*atan(Qf/(2.0d0*x + bf))
        ja = 2.0d0*log(x - x0a) - log(XXa)                                   &
     &      + 2.0d0*((2.0d0*x0a + ba)/Qa)*atan(Qa/(2.0d0*x + ba))

! epsilon derivatives
        ecP = Ap*(2.0d0*log(x) - log(XXp)                                    &
     &            + (2.0d0*bp/Qp)*atan(Qp/(2.0d0*x + bp))                    &
     &            - (bp*x0p/(x0p*x0p + bp*x0p + cp))*jp)
        ecF = Af*(2.0d0*log(x) - log(XXf)                                    &
     &            + (2.0d0*bf/Qp)*atan(Qf/(2.0d0*x + bf))                    &
     &            - (bf*x0f/(x0f*x0f + bf*x0f + cf))*jp)
        ecA = Aa*(2.0d0*log(x) - log(XXa)                                    &
     &            + (2.0d0*ba/Qa)*atan(Qa/(2.0d0*x + ba))                    &
     &            - (ba*x0a/(x0a*x0a + ba*x0a + ca))*ja)

        ecPp = 2.0d0*Ap*cp/(XXp*x) - 2.0d0*Ap*bp*x0p/((x - x0p)*XXp)
        ecFp = 2.0d0*Af*cf/(XXf*x) - 2.0d0*Af*bf*x0f/((x - x0f)*XXf)
        ecAp = 2.0d0*Aa*ca/(XXa*x) - 2.0d0*Aa*ba*x0a/((x - x0a)*XXa)

        ecPpp = - 2.0d0*Ap*cp*(3.0d0*x**2 + 2.0d0*bp*x + cp)/(x*XXp)**2      &
     &          + 2.0d0*Ap*bp*x0p*((2.0d0*x + bp)*(x - x0p) + XXp)           &
     &                 /(XXp*(x - x0p))**2
        ecFpp = - 2.0d0*Af*cf*(3.0d0*x**2 + 2.0d0*bf*x + cf)/(x*XXf)**2      &
     &          + 2.0d0*Af*bf*x0f*((2.0d0*x + bf)*(x - x0f) + XXf)           &
     &                 /(XXf*(x - x0f))**2
        ecApp = - 2.0d0*Aa*ca*(3.0d0*x**2 + 2.0d0*ba*x + ca)/(x*XXa)**2      &
     &          + 2.0d0*Aa*ba*x0a*((2.0d0*x + ba)*(x - x0a) + XXa)           &
     &                 /(XXa*(x - x0a))**2

        cte = 4.0d0/(9.0d0*(2.0d0**(1.0d0/3.0d0) - 1.0d0))

        h = cte*((ecF - ecP)/ecA) - 1.d0
        hp = cte*((ecFp - ecPp)/ecA - (ecF - ecP)*(ecAp/ecA))
        hpp = cte*((ecFpp - ecPpp)/ecA - (ecFp - ecPp)*ecAp/ecA**2           &
     &            - (ecF - ecP)*ecApp/ecA - (ecFp - ecPp)*ecAp/ecA           &
     &            + (ecF - ecP)*(ecAp/ecA)**2)


! Correlation functional ( and partials to z and x ):
        if (zeta .ne. 0.0d0) then
         ec = ecP + ecA*g*(1 + h*zeta**4)
        else
         ec = ecP
        end if

        ecpx = ecPp + ecAp*g*(1 + h*zeta**4) + ecA*g*hp*zeta**4
        ecpz = ecA*gp*(1.0d0 + h*zeta**4) + ecA*g*h*4*zeta**3

        ecppx = ecPp + ecApp*g*(1.0d0 + h*zeta**4) + 2.0d0*ecAp*g*hp*zeta**4 &
     &         + ecA*g*hpp*zeta**4
        ecppz = ecA*gpp*(1.0d0 + h*zeta**4) + ecA*gp*h*zeta**3               &
     &         + ecA*g*h*12.0d0*zeta**2
        ecpxpz = ecAp*gp*(1.0d0 + h*zeta**4) + ecA*gp*hp*zeta**4             &
     &          + ecAp*g*h*4.0d0*zeta**3 + ecA*g*hp*4.0d0*zeta**3

! Partial derivatives VWN exchanche functional
        d1ec = xp*ecpx + zp1*ecpz
        d2ec = xp*ecpx + zp2*ecpz

! Second partial derivatives
        dd1ec = xp**2*ecpp + 2.0d0*xp*zp1*ecpxpz + xpp*ecpx                  &
     &         + zp1*zp1*ecppz + zpp1*ecpz
        dd2ec = xp**2*ecpp + 2.0d0*xp*zp2*ecpxpz + xpp*ecpx                  &
     &         + zp2*zp1*ecppz + zpp2*ecpz
        d1d2ec = xp**2*ecpp+ xp*(zp1 + zp2)*ecpxpz + xpp*ecpx                &
     &          + zp1*zp2*ecppz + zp1p2*ecpz

! ****************************************************************************
!
!       VNN EXCHANGE FUNCTIONAL
!
! ****************************************************************************
        exP = (-3.0d0/2.0d0)*(3.0d0*density/pi)**(1.0d0/3.0d0)
        exPp = exP/(3.0d0*density)
        exPpp = - 2.0d0*exP/(3.0d0*density)**2

        ex =(1.0d0 + 4.0d0*g/9.0d0)*exP
        expd = ex/(3.0d0*density)
        exppd = -2.0d0*ex/(9.0d0*density**2)
        expz = exP*gp
        exppz = exP*gpp
        expdpz = exPp*gp

        d1ex = expd + zp1*expz
        d2ex = expd + zp2*expz

        dd1ex = exppd + 2.0d0*zp1*expdpz + expd + zp1*zp1*exppz + zpp1*expz
        dd2ex = exppd + 2.0d0*zp2*expdpz + expd + zp2*zp2*exppz + zpp2*expz
        d1d2ex = exppd + (zp1 + zp2)*expdpz + expd + zp1*zp2*exppz + zp1p2*expz

! Functions in Rydberg units - divide by factor of 2 to get Hartree
! after that multiply by a factor of 27.2 to get eV
! ****************************************************************************
        xpot(1) = 0.5d0*(density*d1ex + ex)*(P_hartree)
        xpot(2) = 0.5d0*(density*d2ex + ex)*(P_hartree)
        cpot(1) = 0.5d0*(density*d1ec + ec)*(P_hartree)
        cpot(2) = 0.5d0*(density*d2ec + ec)*(P_hartree)
        ex = 0.5d0*ex/(P_hartree)
        ec = 0.5d0*ec/(P_hartree)

        cdpot(1) = 0.5d0*dd1ec*(P_hartree)
        cdpot(2) = 0.5d0*d1d2ec*(P_hartree)
        cdpot(3) = 0.5d0*dd2ec*(P_hartree)
        xdpot(1) = 0.5d0*dd1ex*(P_hartree)
        xdpot(2) = 0.5d0*d1d2ex*(P_hartree)
        xdpot(3) = 0.5d0*dd2ex*(P_hartree)

! This are the 2nd derivatives of the energy respect to n and n_s. In theory we don't need them
        dnuxc = 0.25d0*density*(xdpot(1) + 2.0d0*xdpot(2) + xdpot(3))        &
     &         + 0.5d0*(d1ec + d2ec) + 0.5d0*(d1ex + d2ex)                   &
     &         + 0.25d0*density*(cdpot(1) + 2.0d0*cdpot(2) + cdpot(3))

        dnuxcs = 0.25d0*density*(xdpot(1) - 2.0d0*xdpot(2) + xdpot(3))       &
     &           + 0.5d0*(d1ec - d2ec) + 0.5d0*(d1ex - d2ex)                 &
     &           + 0.25d0*density*(cdpot(1) - 2.0d0*cdpot(2) + cdpot(3))

        dnuxcs = 0.5d0*(ecA + 4.0d0*exP/9.0d0)/density
! What we need is the derivatives of the potential respect to n_up and n_down
! dxcpot(sigma,sigmaprime) = dV_xc_sigma/dn_sigmaprime

  dxpot(1,1) = 0.5d0*(density*dd1ex + 2*d1ex)*(P_hartree)
  dcpot(1,1) = 0.5d0*(density*dd1ec + 2*d1ec)*(P_hartree)
  dxpot(1,2) = 0.5d0*(density*d1d2ex + 2*d2ex)*(P_hartree)
  dcpot(1,2) = 0.5d0*(density*d1d2ec + 2*d2ec)*(P_hartree)
  dxpot(2,1) = 0.5d0*(density*d1d2ex + 2*d1ex)*(P_hartree)
  dcpot(2,1) = 0.5d0*(density*d1d2ec + 2*d1ec)*(P_hartree)
  dxpot(2,2) = 0.5d0*(density*dd2ex + 2*d2ex)*(P_hartree)
  dcpot(2,2) = 0.5d0*(density*dd2ec + 2*d2ec)*(P_hartree)




! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine lsdavwn

     end module M_assemble_vxc

