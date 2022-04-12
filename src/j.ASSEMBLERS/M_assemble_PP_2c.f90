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
! Dublin Institute of Technology - Khorgolkhuu Odbadrakh
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! M_assemble_PP_2c
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the matrix elements for the two-center interactions
!! related to the pseudopotential interactions.
!! It contains the following subroutines within the module:
!!
!!       assemble_svnl.f90 - assemble separable pseudopotential pieces
!!       assemble_vnl.f90 - assemble total pseudopotential pieces
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_assemble_PP_2c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_2c

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! Subroutine: assemble_svnl()
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates separable non-local pseudo-potential matrix
!! interactions.
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
        subroutine assemble_svnl (s)
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
        integer in1, in2                 !< species numbers
        integer jatom                    !< neighbor of iatom
        integer interaction, isubtype    !< which interaction and subtype
        integer num_neigh                !< number of neighbors
        integer mbeta                    !< cell containing neighbor of iatom

        integer norb_mu, norb_nu         !< size of the block for the pair

        real z                            !< distance between r1 and r2

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sighat    !< unit vector along r2 - r1

        real, dimension (:, :), allocatable :: svnlm
        real, dimension (:, :), allocatable :: svnlx

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: psvnl_neighbors
        type(T_assemble_neighbors), pointer :: psvnl

! Allocate Arrays
! ===========================================================================
        allocate (s%svnl(s%natoms))

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          psvnl=>s%svnl(iatom)

          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors_PP(iatom)%neighn
          allocate (psvnl%neighbors(num_neigh))

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            psvnl_neighbors=>psvnl%neighbors(ineigh)
            mbeta = s%neighbors_PP(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_PP_max
            allocate (psvnl_neighbors%block(norb_mu, norb_nu))
            psvnl_neighbors%block = 0.0d0

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

! CALL DOSCENTROSPP AND GET SVNL
! ****************************************************************************
! For the vna_ontopL case, the potential is in the first atom - left (iatom):
            isubtype = 0
            interaction = P_vnl

! Allocate block size
            allocate (svnlm (norb_mu, norb_nu))
            allocate (svnlx (norb_mu, norb_nu))
            call getMEs_Fdata_2c (in1, in2, interaction, isubtype, z,        &
     &                            norb_mu, norb_nu, svnlm)
            call rotate_PP (in1, in2, eps, norb_mu, norb_nu, svnlm, svnlx)

            psvnl_neighbors%block = svnlx

            deallocate (svnlm, svnlx)
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
        end subroutine assemble_svnl


! ===========================================================================
! assemble_vnl_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine takes all the separable pseudopotential interactions and
!! combines them together to build the Hamiltonian matrix elements.
!
! ===========================================================================
! Code written by:
!> @author Khorgolkhuu Odbadrakh
!! @author James P. Lewis
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
        subroutine assemble_vnl_2c (s)
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
        integer in1, in2                 !< species numbers
        integer jatom                    !< neighbor of iatom
        integer matom                    !< matom is the self-interaction atom
        integer mbeta                    !< cell containing neighbor of iatom
        integer mneigh_self, jneigh
        integer kneigh, num_neigh

        integer imu, inu
        integer ncc                      !< counter over pseudo-orbitals
        integer norb_mu, norb_nu         !< size of the block for the pair

        real, pointer :: cl_value (:)
        real, dimension (:, :), allocatable :: PPx

        type(T_assemble_block), pointer :: pvnl_neighbors
        type(T_assemble_block), pointer :: psvnl_neighbors
        type(T_assemble_block), pointer :: psvnl1_neighbors
        type(T_assemble_block), pointer :: psvnl2_neighbors

        type(T_assemble_neighbors), pointer :: pvnl
        type(T_assemble_neighbors), pointer :: psvnl
        type(T_assemble_neighbors), pointer :: psvnl1
        type(T_assemble_neighbors), pointer :: psvnl2

        interface
          function cl(ispecies)
            real, pointer :: cl (:)
            integer, intent(in) :: ispecies
          end function cl
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (s%vnl(s%natoms))
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          pvnl=>s%vnl(iatom)
          num_neigh = s%neighbors_PPp(iatom)%neighn
          allocate (pvnl%neighbors(num_neigh))

          do ineigh = 1, num_neigh   !  <==== loop over i's neighbors
            jatom = s%neighbors_PPp(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some lengthy notation
            pvnl_neighbors=>pvnl%neighbors(ineigh)
            allocate (pvnl_neighbors%block(norb_mu, norb_nu))
            pvnl_neighbors%block = 0.0d0
          end do
        end do

! Procedure
! ===========================================================================
! ASSEMBLE VNL ATM CASE  <phi_i|Psi_j><Psi_j|phi_i>
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pvnl=>s%vnl(iatom)
          psvnl=>s%svnl(iatom)

          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          matom = s%neighbors_PPp_self(iatom)
          num_neigh = s%neighbors_PP(iatom)%neighn

          ! cut some lengthy notation
          pvnl_neighbors=>pvnl%neighbors(matom)

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh        !  <==== loop over iatom's neighbors
            mbeta = s%neighbors_PP(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass

            ! cut some lengthy notation
            psvnl_neighbors=>psvnl%neighbors(ineigh)

! Get the coefficients
! We now loop though all shells, and create cl for each orbital.  For example,
! sp^3 has two shells; cl(1) = cl_PP(0) and cl(2) = cl(3) = cl(4) = cl_PP(1).
            ! memory is allocated inside function
            cl_value => cl(in2)

! Now we combine and sum:
! in1 twice because it is an atom case.
            allocate (PPx (norb_mu, norb_mu)); PPx = 0.0d0
            do inu = 1, norb_mu
              do imu = 1, norb_mu
                do ncc = 1, species(in2)%norb_PP_max
                  PPx(imu,inu) = PPx(imu,inu)                                &
     &             + cl_value(ncc)*psvnl_neighbors%block(imu,ncc)            &
     &                            *psvnl_neighbors%block(inu,ncc)
                end do
              end do
            end do

! Final (not nearly as final!) assembly of vnl - the energy piece.
            pvnl_neighbors%block = pvnl_neighbors%block + PPx

            deallocate (PPx)
            deallocate (cl_value)
          end do ! end do ineigh
        end do ! end do iatom

! ===========================================================================
! ASSEMBLE VNL ONTOP LEFT CASE   <phi_i|Psi_i><Psi_i|phi_j>
! ===========================================================================
! Loop over iatom
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pvnl=>s%vnl(iatom)
          psvnl1=>s%svnl(iatom)

          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors_PPx(iatom)%neighn

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh        !  <==== loop over i's neighbors
            mbeta = s%neighbors_PPx(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PPx(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some lengthy notation
            psvnl2=>s%svnl(jatom)

! if r1 .ne. r2, then this is a case where the potential is located
! at one of the sites of a wavefunction (ontop case).
! Case 1. PP is iatom.  <i | VNL(i) |j>.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then ! do nothing
              if (s%neighbors_PPx_self(iatom) .ne. ineigh) then
                write (*,*) ' neighbors_PPx_self(iatom) .ne. ineigh',        &
     &                       s%neighbors_PPx_self(iatom), ineigh
                stop
              end if ! if(neighPP_self)
            else

              ! memory is allocated inside function
              cl_value => cl(in1)

              mneigh_self = s%neighbors_PP_self(iatom)
              jneigh = s%neighbors_PPx(iatom)%point(ineigh)

              ! cut some lengthy notation
              psvnl1_neighbors=>psvnl1%neighbors(mneigh_self)
              psvnl2_neighbors=>psvnl2%neighbors(jneigh)

! <phi_i|Psi_i>  ->  nPP(mneigh_self,iatom)
! Now we combine and sum:
              allocate (PPx (norb_mu, norb_nu)); PPx = 0.0d0
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  do ncc = 1, species(in1)%norb_PP_max
                    PPx(imu,inu) = PPx(imu,inu)                              &
     &                + cl_value(ncc)*psvnl1_neighbors%block(imu,ncc)        &
     &                       *psvnl2_neighbors%block(inu,ncc)
                  end do ! do ncc
                end do ! do imu
              end do ! do inu

! Mapping to the global matrix
              kneigh = s%neighbors_PPx(iatom)%map(ineigh)
              pvnl_neighbors=>pvnl%neighbors(kneigh)

! Assemble the global matrix
              pvnl_neighbors%block = pvnl_neighbors%block + PPx

              deallocate (PPx)
              deallocate (cl_value)
            end if
          end do ! do ineigh
        end do ! do iatom

! ===========================================================================
! ASSEMBLE VNL ONTOP RIGHT CASE   <phi_i|Psi_j><Psi_j|phi_j>
! ===========================================================================
! Loop over iatom
        do iatom = 1,s%natoms
          ! cut some lengthy notation
          pvnl => s%vnl(iatom)
          psvnl1 => s%svnl(iatom)

          in1 = s%atom(iatom)%imass
          norb_mu =  species(in1)%norb_max
          num_neigh = s%neighbors_PP(iatom)%neighn

! Loop over the neighbors of each iatom.
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors_PP(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some lengthy notation
            psvnl2=>s%svnl(jatom)

! if r1 .ne. r2, then this is a case where the potential is located
! at one of the sites of a wavefunction (ontop case).
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.
! sanity check
              if (s%neighbors_PP_self(iatom) .ne. ineigh) then
                write (*,*) 'neighbors_PP_self(iatom) .ne. ineigh',          &
     &                       s%neighbors_PP_self(iatom), ineigh
                stop
              end if ! if (neighPP_self)
              else ! if (iatom .eq. jatom)

              !memory is allocated inside function
              cl_value => cl(in2)

! Now the second case. <i | V(j) | j>.
! Looking for <phi_j|Psi_j>, what is jneigh of jatom itself in the nPPx list
              mneigh_self = s%neighbors_PP_self(jatom)

              ! cut some lengthy notation
              psvnl1_neighbors => psvnl1%neighbors(ineigh)
              psvnl2_neighbors => psvnl2%neighbors(mneigh_self)

! Now we combine and sum:
              allocate (PPx(norb_mu, norb_nu)); PPx = 0.0d0
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  do ncc = 1, species(in2)%norb_PP_max
                    PPx(imu,inu) = PPx(imu,inu)                              &
     &              + cl_value(ncc)*psvnl1_neighbors%block(imu,ncc)          &
     &                             *psvnl2_neighbors%block(inu,ncc)
                  end do
                end do
              end do

! Mapping to the global matrix
              kneigh = s%neighbors_PP(iatom)%map(ineigh)
              pvnl_neighbors => pvnl%neighbors(kneigh)

! Assemble the global matrix
              pvnl_neighbors%block = pvnl_neighbors%block + PPx

              deallocate(PPx)
              deallocate (cl_value)
            end if ! if(iatom .eq. jatom)

! End loop over iatom and its neighbors - jatom.
          end do ! do ineigh
        end do ! do iatom


! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_vnl_2c


! ===========================================================================
! destroy_assemblePP_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the assemblePP_2c
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
        subroutine destroy_assemble_PP_2c (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                             !< counter over atoms
        integer ineigh                            !< counter over neighbors

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors_PP(iatom)%neighn
            deallocate (s%svnl(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%svnl(iatom)%neighbors)
        end do
        deallocate (s%svnl)

        do iatom = 1, s%natoms
          do ineigh=1, s%neighbors_PPp(iatom)%neighn
            deallocate (s%vnl(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%vnl(iatom)%neighbors)
        end do
        deallocate (s%vnl)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_assemble_PP_2c


! End Module
! ===========================================================================
        end module M_assemble_PP_2c
