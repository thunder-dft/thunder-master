! copyright info:
!
!                             @Copyright 2014
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
        use M_configuraciones
        use M_grid
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
        integer ipoint                    !< counter over grid points
        integer in1, in2, in3             !< species numbers
        integer jatom                     !< neighbor of iatom
        integer interaction, isubtype     !< which interaction and subtype
        integer logfile                   !< writing to which unit
        integer num_neigh                 !< number of neighbors
        integer mbeta                     !< the cell containing iatom's neighbor
        integer issh, jssh                !< counting over orbitals
        integer norb_mu, norb_nu          !< size of the block for the pair

        real density                      !< value of density at grid point
        real u0_tot

        real uhdcc                        !< double counting corrections - grid
        real dvhdn, dvhn0

        real uee_self_tot
        real z                            !< distance between atom pairs
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

! Now calculate the double-counting term for the Poisson solution.
! Hartree term
        uhdcc = 0.0d0
        dvhn0 = 0.0d0
        dvhdn = 0.0d0

! Loop over points in the grid
        do ipoint = 0, nrm - 1
          density = rhoG0(ipoint) + drhoG(ipoint)
          uhdcc = uhdcc + vcaG(ipoint)*density

          ! functional derivative pieces integrated
          dvhn0 = dvhn0 + vcaG(ipoint)*rhoG0(ipoint)
          dvhdn = dvhdn + vcaG(ipoint)*drhoG(ipoint)
        end do

        ! multiply by integration factor
        dvhn0 = -1.0d0*dvhn0*dvolume
        dvhdn = -0.5d0*dvhdn*dvolume
        uhdcc = dvhn0 + dvhdn

        ! Final double-counting correction for coulomb part
        uii_uee = u0_tot - uee_self_tot + uhdcc

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
! assebmle_uxc
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the exchange-correlation double-counting
!! energy for density functional theory on the grid.
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_uxc (s, uxcdcc)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

        !< double counting correction for coulomb part
        real, intent (out) :: uxcdcc

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ipoint                  !< counter over grid point

! Loops
        real density                     !< input density

        ! parameters for lda_ceperley_alder
        real exc, muxc, dexc, d2exc, dmuxc, d2muxc

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize
        uxcdcc = 0.d0

! Loop over points in the grid
        do ipoint = 0, nrm - 1
          density = drhoG(ipoint) + rhoG0(ipoint)

! calculate the exchange-correlation potential
          call lda_ceperley_alder (density, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
          uxcdcc = uxcdcc + (exc - muxc)*density
        end do

        ! multiply by integration factor
        uxcdcc = uxcdcc*dvolume

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

