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

! M_assemble_PP_3c.f90
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the matrix elements for the three-center interactions
!! related to the pseudopotential interactions.
!! It contains the following subroutines within the module:
!!
!! Dassemble_vnl_3a - assemble total pseudopotential forces for three center
!
! ===========================================================================
        module M_Dassemble_PP_3c

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones
        use M_Drotations_PP

! /FDATA
        use M_Fdata_3c

! /ASSEMBLERS
        use M_assemble_PP_3c

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! Dassemble_vnl_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions.
! These are true three-center in that iatom .ne. jatom .ne. katom.
! The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
! A third party term is when we consider the NA (or etc.) to be at the origin
! and take the matrix element between a pair of neighbors, neither of which is
! the NA (or etc.), and the pair is truly a pair, and not an atom.
! the results are stored in: f3na(ix,i), f3xc(ix,i), etc.
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms.!
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
! Subroutine Declaration
! ===========================================================================
        subroutine Dassemble_vnl_3c (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_3c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha, iatom, jatom         !< counter over atoms and neighbors
        integer in1, in2, indna              !< species numbers
        integer num_neigh                    !< number of neighbors
        integer ibeta, jbeta
        integer ineigh, mneigh, m31, m32
        integer ncc

        integer imu, inu
        integer norb_mu, norb_nu         !< size of the block for the pair

        real, pointer :: cl_value (:)

        real, dimension (3) :: r1, r2, rna, r31, r32

        interface
          function cl(itype)
            real, pointer :: cl (:)
            integer, intent(in) :: itype
          end function cl
        end interface

! interaction arrays
        type(T_assemble_block), pointer :: psvnl1_neighbors
        type(T_assemble_block), pointer :: psvnl2_neighbors

        type(T_assemble_neighbors), pointer :: psvnl1
        type(T_assemble_neighbors), pointer :: psvnl2

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors

        ! forces
        real, dimension (:, :, :), allocatable :: f3nlXa
        real, dimension (:, :, :), allocatable :: f3nlXb
        real, dimension (:, :, :), allocatable :: f3nlXc

        type(T_forces), pointer :: pfalpha
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ASSEMBLE VNL ATM CASE  <phi_i|Psi_j><Psi_j|phi_i>
        do ialpha = 1, s%natoms
          rna = s%atom(ialpha)%ratom
          indna = s%atom(ialpha)%imass
          num_neigh = s%neighbors_PP(ialpha)%ncommon

          ! cut some lengthy notation
          pfalpha=>s%forces(ialpha)

! For the non-local potential pieces call find the coefficients corresponding
! to indna.
          ! memory is allocated inside function
          cl_value => cl(indna)

! Loop over the neighbors of each ialp.
! Now look at all common neighbor pairs which were figured out in main.
! The common neighbors were figured out in common_neighbors.f90
          do ineigh = 1, num_neigh
            mneigh = s%neighbors_PP(ialpha)%neigh_common(ineigh)

! The second atom (jatom) is the mneigh'th neighbor of iatom.
            if (mneigh .ne. 0) then
              iatom = s%neighbors_PP(ialpha)%iatom_common_j(ineigh)
              ibeta = s%neighbors_PP(ialpha)%iatom_common_b(ineigh)
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              in1 = s%atom(iatom)%imass
              norb_mu = species(in1)%norb_max

              jatom = s%neighbors_PP(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors_PP(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%norb_max

! PERFORM ACTUAL CALCULATIONS
! ASSEMBLE NON-LOCAL PSEUDOPOTENTIAL PIECE.
! Here 1=iatom, 2=jatom, 3=ialp:  <1 | V(3) | 2>.
              r31 = rna - r1
              r32 = rna - r2

! Find m value for iatom, ialp pair, and for jatom, ialp pair.
              call mpairnay (s, iatom, ialpha, r31, m31)
              call mpairnay (s, jatom, ialpha, r32, m32)

              psvnl1 => s%svnl(iatom)
              psvnl1_neighbors => psvnl1%neighbors(m31)

              psvnl2 => s%svnl(jatom)
              psvnl2_neighbors => psvnl2%neighbors(m32)

              ! cut lengthy notation - forces
              pfi=>s%forces(iatom); pfj=>s%forces(jatom)

              ! density matrix
              pdenmat=>s%denmat_PP(iatom); pRho_neighbors=>pdenmat%neighbors(mneigh)

              allocate (f3nlXa(3,norb_mu,norb_nu)); f3nlXa = 0.0d0
              allocate (f3nlXb(3,norb_mu,norb_nu)); f3nlXb = 0.0d0
              allocate (f3nlXc(3,norb_mu,norb_nu)); f3nlXc = 0.0d0

! A "general" matrix element is:
! <iatom|VNL(3)|jatom> = <iatom|V(3)><V(3)|jatom> =
! SUM_ncc cl(ncc)*sVNL(mu,ncc,iatom,m31)*sVNL(nu,ncc,jatom,m32)
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  do ncc = 1, species(indna)%norb_PP_max
                    f3nlXb(:,imu,inu) = f3nlXb(:,imu,inu)                      &
                      + cl_value(ncc)*psvnl1_neighbors%Dblock(:,imu,ncc)       &
                                     *psvnl2_neighbors%block(inu,ncc)
                    f3nlXc(:,imu,inu) = f3nlXc(:,imu,inu)                      &
                      + cl_value(ncc)*psvnl1_neighbors%block(imu,ncc)          &
                                     *psvnl2_neighbors%Dblock(:,inu,ncc)
                  end do ! do ncc
                end do ! do imu
              end do ! do inu

! Make things force-like and determine f3nlXc, which is found from Newtons Laws:
!             f3nlXb = - f3nlXb
!             f3nlXc = - f3nlXc
              f3nlXa = - f3nlXb - f3nlXc

              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfalpha%f3nla = pfalpha%f3nla                              &
      &             - pRho_neighbors%block(imu,inu)*f3nlXa(:,imu,inu)
                  pfi%f3nlb = pfi%f3nlb                                      &
      &             - pRho_neighbors%block(imu,inu)*f3nlXb(:,imu,inu)
                  pfj%f3nlc = pfj%f3nlc                                      &
      &             - pRho_neighbors%block(imu,inu)*f3nlXc(:,imu,inu)
                end do
              end do

! Deallocate Arrays
              deallocate (f3nlXa, f3nlXb, f3nlXc)
            end if  ! end if mneigh .ne. 0
          end do ! end loop over neighbors
          deallocate (cl_value)
        end do ! end loop over ialpha

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Dassemble_vnl_3c

! End Module
! ===========================================================================
        end module M_Dassemble_PP_3c
