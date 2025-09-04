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

! M_Dassemble_PP_2c
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
        module M_Dassemble_PP_2c

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones
        use M_Drotations_PP

! /FDATA
        use M_Fdata_2c

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! denmat_PP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This subroutine calculates the density matix rho and stores it in the
!! structure given in M_assemble_block.f90
!
! ===========================================================================
        subroutine densityPP_matrix (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh              !< counter over atoms and neighbors
        integer iband, ikpoint             !< counter of band and kpoint
        integer imu, inu
        integer in1, in2                   !< species numbers
        integer jatom                      !< neighbor of iatom
        integer mmu, nnu

        integer num_neigh                !< number of neighbors
        integer mbeta                    !< the cell containing iatom's neighbor
        integer norb_mu, norb_nu         !< size of the block for the pair

        real dot                         !< dot product between K and r
        real gutr                        !< real part of density matrix

        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom

        complex phase, phasex            !< phase between K and r
        complex step1, step2

        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_neighbors), pointer :: pdenmat

! Allocate Arrays
! ===========================================================================
        allocate (s%denmat_PP (s%natoms))

! Procedure
! ===========================================================================
! ****************************************************************************
!
!                      C O M P U T E    D E N S I T I E S
! ****************************************************************************
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          nullify (pdenmat)
          pdenmat=>s%denmat_PP(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors_PPp(iatom)%neighn
          allocate (s%denmat_PP(iatom)%neighbors(num_neigh))
          do ineigh = 1, num_neigh        !  <==== loop over iatom's neighbors
            mbeta = s%neighbors_PPp(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PPp(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

            ! cut some more lengthy notation
            nullify (pRho_neighbors)
            pRho_neighbors=>pdenmat%neighbors(ineigh)

! Allocate the block size
            norb_nu = species(in2)%norb_max
            allocate (s%denmat_PP(iatom)%neighbors(ineigh)%block(norb_mu, norb_nu))
            pRho_neighbors%block = 0.0d0

! Loop over the special k points.
            do ikpoint = 1, s%nkpoints

              ! Find the phase which is based on k*r
              vec = r2 - r1
              sks = s%kpoints(ikpoint)%k
              dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
              phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight*P_spin

! Loop over all bands
! Here we assume that the coefficients are real only at this point.
              do iband = 1, s%norbitals_new
                if (s%kpoints(ikpoint)%ioccupy(iband) .ne. 0) then
                  phase = phasex*s%kpoints(ikpoint)%foccupy(iband)
                  do imu = 1, norb_mu
                    mmu = imu + s%iblock_slot(iatom)
                    step1 = phase*conjg(s%kpoints(ikpoint)%c(mmu,iband))

                    do inu = 1, norb_nu
                      nnu = inu + s%iblock_slot(jatom)
                      step2 = step1*s%kpoints(ikpoint)%c(nnu,iband)
                      gutr = real(step2)

! Finally the density matrix:
                      pRho_neighbors%block(imu,inu) =                        &
     &                  pRho_neighbors%block(imu,inu) + gutr
                    end do
                  end do
                end if

! Finish loop over bands.
              end do

! Finish loop over k-points.
            end do

! Finish loop over atoms and neighbors.
            nullify (pRho_neighbors)
          end do
          nullify (pdenmat)
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i4, 2(2x, f10.6))
200     format (4(2x, f12.4))

! End Subroutine
! ===========================================================================
        return
        end subroutine densityPP_matrix


! ===========================================================================
! Dassemble_svnl
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the derivatives of the separable non-local
!>       pseudo-potential matrix interactions.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine Dassemble_svnl (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
! Input
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

        integer imu, inu
        integer norb_mu, norb_nu         !< size of the block for the pair

        real z                           !< distance between r1 and r2

        real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sighat    !< unit vector along r2 - r1

! sm = overlap matrix in molecular coordinates
! sx = overlap matrix in crystal coordinates
! dsm = derivative of overlap matrix in molecular coordinates
! vdsm = vectorized derivative of overlap matrix in molecular coordinates
! vdsx = vectorized derivative of overlap matrix in crystal coordinates
        real, dimension (:, :), allocatable :: svnlm
        real, dimension (:, :), allocatable :: svnlx
        real, dimension (:, :), allocatable :: dsvnlm
        real, dimension (:, :, :), allocatable :: vdsvnlm
        real, dimension (:, :, :), allocatable :: vdsvnlx

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
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          nullify (psvnl)
          psvnl=>s%svnl(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors_PP(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            ! cut some more lengthy notation
            nullify (psvnl_neighbors)
            psvnl_neighbors=>psvnl%neighbors(ineigh)
            mbeta = s%neighbors_PP(iatom)%neigh_b(ineigh)
            jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_PP_max
            allocate (s%svnl(iatom)%neighbors(ineigh)%Dblock(3, norb_mu, norb_nu))
            psvnl_neighbors%Dblock = 0.0d0

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
            call Depsilon_2c (r1, r2, z, eps, deps)

! CALL DOSCENTROSPP AND GET SVNL
! ****************************************************************************
! For the vna_ontopL case, the potential is in the first atom - left (iatom):
            isubtype = 0
            interaction = P_vnl

! Allocate block size
! svnlm = matrix in molecular coordinates
! svnlx = matrix in crystal coordinates
! dsvnlm = derivative of matrix in molecular coordinates
! vdsvnlm = vectorized derivative of matrix in molecular coordinates
! vdsvnlx = vectorized derivative of matrix in crystal coordinates
            allocate (svnlm (norb_mu, norb_nu)); svnlm = 0.0d0
            allocate (svnlx (norb_mu, norb_nu)); svnlx = 0.0d0
            allocate (dsvnlm (norb_mu, norb_nu)); dsvnlm = 0.0d0
            allocate (vdsvnlm (3, norb_mu, norb_nu)); vdsvnlm = 0.0d0
            allocate (vdsvnlx (3, norb_mu, norb_nu)); vdsvnlx = 0.0d0
            call getDMEs_Fdata_2c (in1, in2, interaction, isubtype, z,       &
     &                             norb_mu, norb_nu, svnlm, dsvnlm)

! Apply epsilon, the direction of the bondcharge.
! ****************************************************************************
!
! FORCES
! ****************************************************************************
! dsm is the "scalar" derivative of the matrix; vstm is the "vector" derivative
! of the matrix in molecular coordinates.  When we are done, we get: vdsx as
! the vector derivative of the matrix in crystal coordinates.

! As long as epsilon is called with sighat in the second "spot" as
! call epsilon_function (R1, sighat, spe), then eps(ix,3) = eta(ix).
            eta(:) = eps(:,3)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
            do inu = 1, norb_nu
              do imu = 1, norb_mu
                if (z .gt. 1.0d-3) vdsvnlm(:,imu,inu) = - eta(:)*dsvnlm(imu,inu)
              end do
            end do

! Drotate then puts the vectors in coordinates along the bond-charge.
            call Drotate_PP (in1, in2, eps, deps, norb_mu, norb_nu, svnlm, vdsvnlm, vdsvnlx) 

! Store the derivitive, rotate vector matrix.
            psvnl_neighbors%Dblock = vdsvnlx
            deallocate (svnlm, svnlx, dsvnlm, vdsvnlm, vdsvnlx)
            nullify (psvnl_neighbors)
          end do ! end loop over neighbors
          nullify (psvnl)
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
        end subroutine Dassemble_svnl


! ===========================================================================
! Dassemble_vnl
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine takes all the derivatives of the separable
!> pseudopotential interactions and combines them together to build the
!> Hamiltonian matrix elements.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine Dassemble_vnl_2c (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
! Input
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
        real, dimension (:, :, :), allocatable :: PPx

        ! pointers for matrix elements and derivatives
        type(T_assemble_block), pointer :: psvnl_neighbors
        type(T_assemble_block), pointer :: psvnl1_neighbors
        type(T_assemble_block), pointer :: psvnl2_neighbors

        type(T_assemble_neighbors), pointer :: psvnl
        type(T_assemble_neighbors), pointer :: psvnl1
        type(T_assemble_neighbors), pointer :: psvnl2

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_block), pointer :: pRho_neighbors_matom

        type(T_forces), pointer :: pfi
        
        interface
          function cl(ispecies)
            real, pointer :: cl (:)
            integer, intent(in) :: ispecies
          end function cl
        end interface

! Allocate Arrays
! ===========================================================================
        do iatom = 1, s%natoms
          pfi=>s%forces(iatom)
          num_neigh = s%neighbors_PPp(iatom)%neighn
          allocate (pfi%vnl_atom (3, num_neigh)); pfi%vnl_atom = 0.0d0
          allocate (pfi%vnl_ontop (3, num_neigh)); pfi%vnl_ontop = 0.0d0
        end do

! Procedure
! ===========================================================================
! ASSEMBLE VNL ATM CASE  <phi_i|Psi_j><Psi_j|phi_i>
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max ! do not use PP here
          matom = s%neighbors_PPp_self(iatom)

          ! cut some lengthy notation
          psvnl=>s%svnl(iatom)
          pdenmat=>s%denmat_PP(iatom); pRho_neighbors_matom=>pdenmat%neighbors(matom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors_PP(iatom)%neighn
          do ineigh = 1, num_neigh        !  <==== loop over iatom's neighbors
            if (ineigh .eq. s%neighbors_PP_self(iatom)) cycle
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
            allocate (PPx (3, norb_mu, norb_mu)); PPx = 0.0d0
            do inu = 1, norb_mu
              do imu = 1, norb_mu
                do ncc = 1, species(in2)%norb_PP_max
                  PPx(:,imu,inu) = PPx(:,imu,inu)                            &
     &             + cl_value(ncc)*(psvnl_neighbors%Dblock(:,imu,ncc)        &
     &                              *psvnl_neighbors%block(inu,ncc)          &
     &                              + psvnl_neighbors%block(imu,ncc)         &
     &                               *psvnl_neighbors%Dblock(:,inu,ncc))
                end do
              end do
            end do

! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
! The derivative wrt the "atom r1" position (not the NA position) are
! stored in bcnapx.
!
! Form the "force-like" derivative of the atom terms for NA,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.

! Note that the loop below involves num_orb(in1) ONLY. Why?
! Because the potential is somewhere else (or even at iatom), but we are
! computing the vna_atom term, i.e. < phi(i) | v | phi(i) > but V=v(j) )
! interactions.

! Notice the explicit negative sign, this makes it force like.
            do inu = 1, norb_mu
              do imu = 1, norb_mu
                pfi%vnl_atom(:,ineigh) = pfi%vnl_atom(:,ineigh)              &
      &            - pRho_neighbors_matom%block(imu,inu)*PPx(:,imu,inu)
              end do
            end do
            deallocate (PPx)
            deallocate (cl_value)
          end do ! end do ineigh
        end do ! end do iatom

! ===========================================================================
! ASSEMBLE VNL ONTOP LEFT CASE   <phi_i|Psi_i><Psi_i|phi_j>
! ===========================================================================
! Loop over iatom
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

          ! cut some lengthy notation
          psvnl1=>s%svnl(iatom)
          pdenmat=>s%denmat_PP(iatom)
          pfi=>s%forces(iatom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors_PPx(iatom)%neighn
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
                write (*,*) 'neighbors_PPx_self(iatom) .ne. ineigh',         &
      &                    s%neighbors_PPx_self(iatom), ineigh
                 stop
              end if ! if (neighPP_self)
            else

              ! memory is allocated inside function
              cl_value => cl(in1)

              mneigh_self = s%neighbors_PP_self(iatom)
              jneigh = s%neighbors_PPx(iatom)%point(ineigh)

              ! cut lengthy notation
              psvnl1_neighbors=>psvnl1%neighbors(mneigh_self)
              psvnl2_neighbors=>psvnl2%neighbors(jneigh)

! <phi_i|Psi_i>  ->  nPP(mneigh_self,iatom)
! Note - spVNL are always the derivative of the orbital and NOT the potential.
! Repeat. sVNL(i,m) = <i|V(m)> and spVNL = d/dri <i|VNL(m)>.
! What we need in the next line is
! d/drm <i|VNL(m)> = -d/dr1<i|VNL(m)> = -spVNL(i,m).
! Now we combine and sum:
              allocate (PPx (3, norb_mu, norb_nu)); PPx = 0.0d0
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  do ncc = 1, species(in1)%norb_PP_max
                    PPx(:,imu,inu) = PPx(:,imu,inu)                          &
                      - cl_value(ncc)*psvnl1_neighbors%block(imu,ncc)        &
     &                               *psvnl2_neighbors%Dblock(:,inu,ncc)
                  end do ! do ncc
                end do ! do imu
              end do ! do inu

! Mapping to the global matrix
              kneigh = s%neighbors_PPx(iatom)%map(ineigh)
              pRho_neighbors=>pdenmat%neighbors(kneigh)

! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
! The derivative wrt the "atom r1" position (not the NA position) are
! stored in bcnapx.
!
! Form the "force-like" derivative of the atom terms for NA,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.

! Note that the loop below involves num_orb(in1) ONLY. Why?
! Because the potential is somewhere else (or even at iatom), but we are
! computing the vna_atom term, i.e. < phi(i) | v | phi(i) > but V=v(j) )
! interactions.

! Notice the explicit negative sign, this makes it force like.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfi%vnl_ontop(:,ineigh) = pfi%vnl_ontop(:,ineigh)          &
     &              - pRho_neighbors%block(imu,inu)*PPx(:,imu,inu)
                end do
              end do
              deallocate (PPx)
              deallocate (cl_value)
            end if
          end do ! do ineigh
        end do ! do iatom

! ===========================================================================
! ASSEMBLE VNL ONTOP RIGHT CASE   <phi_i|Psi_j><Psi_j|phi_j>
! ===========================================================================
! Loop over iatom
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          psvnl1 => s%svnl(iatom)
          pdenmat=>s%denmat_PP(iatom)
          pfi=>s%forces(iatom)

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
                write (*,*) ' neighbors_PP_self(iatom) .ne. ineigh',         &
     &                       s%neighbors_PP_self(iatom), ineigh
                stop
              end if ! if (neighPP_self)
            else ! if (iatom .eq. jatom)

              ! memory is allocated inside function
              cl_value => cl(in2)

              mneigh_self = s%neighbors_PP_self(jatom)

              ! cut lengthy notation
              psvnl1_neighbors => psvnl1%neighbors(ineigh)
              psvnl2_neighbors => psvnl2%neighbors(mneigh_self)

! Note - spVNL are always the derivative of the orbital and NOT the potential.
! Repeat. sVNL(i,m) = <i|V(m)> and spVNL = d/dri <i|VNL(m)>.
! What we need in the next line is
! d/drm <i|VNL(m)> = -d/dr1<i|VNL(m)> = -spVNL(i,m).
! Now we combine and sum:
              allocate (PPx(3, norb_mu, norb_nu)); PPx = 0.0d0
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  do ncc = 1, species(in2)%norb_PP_max
                    PPx(:,imu,inu) = PPx(:,imu,inu)                            &
      &              - cl_value(ncc)*psvnl1_neighbors%Dblock(:,imu,ncc)        &
      &                             *psvnl2_neighbors%block(inu,ncc)
                  end do
                end do
              end do

! Mapping to the global matrix
              kneigh = s%neighbors_PP(iatom)%map(ineigh)
              pRho_neighbors=>pdenmat%neighbors(kneigh)

! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
! The derivative wrt the "atom r1" position (not the NA position) are
! stored in bcnapx.
!
! Form the "force-like" derivative of the atom terms for NA,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.

! Note that the loop below involves num_orb(in1) ONLY. Why?
! Because the potential is somewhere else (or even at iatom), but we are
! computing the vna_atom term, i.e. < phi(i) | v | phi(i) > but V=v(j) )
! interactions.

! Notice the explicit negative sign, this makes it force like.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
!                 pfi%vnl_ontop(:,ineigh) = pfi%vnl_ontop(:,ineigh)         &
!    &             + pRho_neighbors%block(imu,inu)*PPx(:,imu,inu)
                end do
              end do
              deallocate (PPx)
              deallocate (cl_value)
            end if  ! if (iatom .eq. jatom)
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
        end subroutine Dassemble_vnl_2c

! ===========================================================================
! destroy_Dassemble_PP_2c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the Dassemble_2c
!! information.
!
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
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_Dassemble_PP_2c (s)
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
            deallocate (s%svnl(iatom)%neighbors(ineigh)%Dblock)
          end do
        end do

        do iatom = 1, s%natoms
          deallocate (s%forces(iatom)%vnl_atom)
          deallocate (s%forces(iatom)%vnl_ontop)
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
        end subroutine destroy_Dassemble_PP_2c

! End Module
! ===========================================================================
        end module M_Dassemble_PP_2c
