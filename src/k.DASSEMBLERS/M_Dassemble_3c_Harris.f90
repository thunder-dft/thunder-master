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

! Dassemble_3c.f90
! Program Description
! ===========================================================================
!!      This routine assembles all of the true three-center interactions.
!! These are true three-center in that iatom .ne. jatom .ne. katom.
!! The matrix elements look like <psi(1)|V(3)|psi(2)>.

!! It contains the following subroutines within the module:
!!
!!       Dassemble_vna_3c.f90 - assemble neutral atom potential matrix derivatives
!!
!! For a complete list of the interactions see the files 3c.Z1.Z2.Z3.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================

! Module Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions.
! These are true three-center in that iatom .ne. jatom .ne. katom.
! The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
!! A third party term is when we consider the NA (or etc.) to be at the
!! origin and take the matrix element between a pair of neighbors, neither
!! of which is the NA (or etc.), and the pair is truly a pair, and not an
!! atom.
!
!! The generic form of the matrix element is
!! v = <mu(r-(r1+ratm)) | v(r-ratm) | nu(r-(r2+ratm))>.
!! We want the derivative with respect to ratm keeping the position of
!! orbitals mu, nu fixed (constant r1 + ratm).
!
!! The derivative f3 looks like for the neutral atom case:
!! f3na = - sum (all neighbors of atom alpha at (li,bi) but bi.ne.balph)
!!    * sum (all neighbors m of (li,bi), and not having b value balph)
!!    * rho(mu,nu,i,m)* deriv wrt balpha <i! v(balph) !j>.
!! Note the minus sign to make it "force-like".
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms. Breifly fa=-d/dratm, fb=-d/d(bi), fc=-d/d(bj)
! where matrix elements are < psi(r-bi) ! v(r-ratm) ! psi(r-bj)>
!
! ===========================================================================
        module M_Dassemble_3c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_3c
        use M_neighbors
        use M_rotations
        use M_Drotations
        use M_assemble_3c
        use M_density_matrix

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! Dassemble_vna_3c.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates derivative of the neutral atom potential
! matrix interactions.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
!> @author Barry Haycock
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
        subroutine Dassemble_vna_3c (s)

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
        integer ialpha, iatom, jatom     !< the three parties involved
        integer ibeta, jbeta             !< cells for three atoms
        integer ineigh, mneigh           !< counter over neighbors
        integer in1, in2, in3            !< species numbers
        integer interaction, isorp       !< which interaction and subtype
        integer imu, inu, iindex         !< indexing counters
        integer norb_mu, norb_nu         !< size of the block for the pair

        real z                           !< distance between r1 and r2
        real x, cost                     !< dnabc and angle

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3) :: r1, r2, rna  !< positions - iatom, jatom, ialpha
        real, dimension (3) :: r12, r21  !< vectors
        real, dimension (3) :: sighat    !< unit vector along r2 - r1
        real, dimension (3) :: rhat      !< unit vector along bc - rna

        real, dimension (3, 3, 3) :: depsA  !< the Depsilon matrix for the bc
        real, dimension (3, 3, 3) :: depsB  !< the Depsilon matrix for the na

        real, dimension (3) :: amt, bmt

! bcnam = Hartree matrix in molecular coordinates
! bcnax = Hartree matrix in crystal coordinates
! dbcnam = derivative of Hartree matrix in molecular coordinates
! vdbcnam = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcnax = vectorized derivative of Hartree matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcnam

        real, dimension (:, :), allocatable :: dpbcnam
        real, dimension (:, :), allocatable :: dxbcnam
        real, dimension (:, :), allocatable :: dybcnam

        real, dimension (:, :, :), allocatable :: f3naMa
        real, dimension (:, :, :), allocatable :: f3naMb

        real, dimension (:, :, :), allocatable :: f3naXa
        real, dimension (:, :, :), allocatable :: f3naXb
        real, dimension (:, :, :), allocatable :: f3naXc

        type(T_Fdata_cell_3c), pointer :: pFdata_cell
        type(T_Fdata_bundle_3c), pointer :: pFdata_bundle

        ! forces
        type(T_forces), pointer :: pfalpha
        type(T_forces), pointer :: pfi
        type(T_forces), pointer :: pfj

        ! density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors

        interface
          function distance (a, b)
            real distance
            real, intent (in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          in3 = s%atom(ialpha)%imass
          rna = s%atom(ialpha)%ratom

          ! cut some lengthy notation
          pfalpha=>s%forces(ialpha)

          ! loop over the common neighbor pairs of ialp
          do ineigh = 1, s%neighbors(ialpha)%ncommon
            mneigh = s%neighbors(ialpha)%neigh_common(ineigh)
            if (mneigh .ne. 0) then
              iatom = s%neighbors(ialpha)%iatom_common_j(ineigh)
              ibeta = s%neighbors(ialpha)%iatom_common_b(ineigh)
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              in1 = s%atom(iatom)%imass
              norb_mu = species(in1)%norb_max

              jatom = s%neighbors(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%norb_max

              ! cut lengthy notation
              pfi=>s%forces(iatom); pfj=>s%forces(jatom)

              ! density matrix
              pdenmat=>s%denmat(iatom); pRho_neighbors=>pdenmat%neighbors(mneigh)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
              r21 = r2 - r1
              z = distance (r1, r2)

              ! unit vector in sigma direction.
              if (z .lt. 1.0d-05) then
                sighat(1) = 0.0d0
                sighat(2) = 0.0d0
                sighat(3) = 1.0d0
              else
                sighat = (r2 - r1)/z
              end if

! ****************************************************************************
! Find rnabc = vector pointing from center of bondcharge to rna
! This gives us the distance dnabc.
              r12 = 0.5d0*(r1 + r2)
              x = distance (r12, rna)

              ! unit vector in rnabc direction.
              if (x .lt. 1.0d-05) then
                rhat(1) = 0.0d0
                rhat(2) = 0.0d0
                rhat(3) = 0.0d0
              else
                rhat = (rna - 0.5d0*(r1 + r2))/x
              end if
              cost = dot_product(sighat, rhat)
              call epsilon_function (rhat, sighat, eps)

! dera3 = depsA = deps/dratm in the 3-center system
! der13 = depsB = deps/dr1 in the 3-center system
              call Depsilon_3c (r1, r2, r21, z, rna, rhat, eps, depsA, depsB)

! For now we just do the neutral atom interactions.
! Charged atom interactions are assembled in assemble_ca_3c.f
! So set isorp = 0 within this subroutine.
!
!              interaction    subtypes     index
!
!      bcna         1           0..9(max)   1..10

! ****************************************************************************
!
! Get the D-matrices from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_bcna
              isorp = 0

! Allocate block arrays
              allocate (bcnam (norb_mu, norb_nu)); bcnam = 0.00
              allocate (dpbcnam (norb_mu, norb_nu)); dpbcnam = 0.00
              allocate (dxbcnam (norb_mu, norb_nu)); dxbcnam = 0.00
              allocate (dybcnam (norb_mu, norb_nu)); dybcnam = 0.00
              call getDMEs_Fdata_3c (in1, in2, in3, interaction, isorp, x,     &
     &                               z, norb_mu, norb_nu, cost, rhat, sighat,  &
     &                               bcnam, dpbcnam, dxbcnam, dybcnam)

              allocate (f3naMa(3, norb_mu, norb_nu)); f3naMa = 0.0d0
              allocate (f3naMb(3, norb_mu, norb_nu)); f3naMb = 0.0d0
              allocate (f3naXa(3, norb_mu, norb_nu)); f3naXa = 0.0d0
              allocate (f3naXb(3, norb_mu, norb_nu)); f3naXb = 0.0d0
              allocate (f3naXc(3, norb_mu, norb_nu)); f3naXc = 0.0d0

! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.

! The first piece will be the force with respect to atom 3.
              if (x .gt. 1.0d-5) then
                amt = (sighat - cost*rhat)/x
              else
                amt = 0.0d0
              end if

! The second piece will be the force with respect to atom 1.
              bmt = (cost*sighat - rhat)/z

              pFdata_bundle => Fdata_bundle_3c(in1, in2, in3)
              pFdata_cell =>                                                &
     &          pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(interaction,isorp,1))

              ! loop over matrix elements
              do iindex = 1, pFdata_cell%nME
                imu = pFdata_cell%mu_3c(iindex)
                inu = pFdata_cell%nu_3c(iindex)

! Now recover f3naMa which is a three-dimensional array
                f3naMa(:,imu,inu) = rhat*dxbcnam(imu,inu) + amt*dpbcnam(imu,inu)

! Now recover f3naMb which is a three-dimensional array
                f3naMb(:,imu,inu) = - sighat*dybcnam(imu,inu)               &
     &           + bmt*dpbcnam(imu,inu) - f3naMa(:,imu,inu)/2.0d0
              end do ! end loop over matrix elements

! ***************************************************************************
! Convert to Crystal Coordinates
! ***************************************************************************
! The call to rotated does the rotations to crystal coordinates of these
! force things.
!
! For example:
! Suppose we have f_A(3,mu,nu), which is d/dratm M(mu,nu) where M(mu,nu)
! is in molecular. To transform M(mu,nu) to crystal, we need Udag * M * U.
! Therefore, f_A(3,mu,nu)[CRYSTAL] = (d/dratm Udag) * M * U
!                                   + Udag * M * (d/dratm U)
!                                   + Udag * f_A * U.
!
! So, to use this baby, put in deps3c (deps/dr1, deps/dr2, deps/dratm),
! and f_A and M.
!
! NOTE: rotated works on the assumption that we are adding derivatives,
! NOT forces. So f3naMa,... etc. MUST not yet be forcelike.
! We do the - sign for forces at the end.
! ***************************************************************************
! Force on the neutral atom with respect to atom 3 (f3naMa).
              call Drotate (in1, in2, eps, depsA, norb_mu, norb_nu,          &
     &                      bcnam, f3naMa, f3naXa)

! Force on the neutral atom with respect to atom 1 (f3naMb).
              call Drotate (in1, in2, eps, depsB, norb_mu, norb_nu,          &
     &                      bcnam, f3naMb, f3naXb)

! Make things force-like and determine f3naXc, whcih is found from Newtons Laws:
              f3naXa = - f3naXa
              f3naXb = - f3naXb
              f3naXc = - f3naXa - f3naXb

              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfalpha%f3naa = pfalpha%f3naa                              &
      &             - pRho_neighbors%block(imu,inu)*f3naXa(:,imu,inu)
                  pfi%f3nab = pfi%f3nab                                      &
      &             - pRho_neighbors%block(imu,inu)*f3naXb(:,imu,inu)
                  pfj%f3nac = pfj%f3nac                                      &
      &             - pRho_neighbors%block(imu,inu)*f3naXc(:,imu,inu)
                end do
              end do

              deallocate (bcnam, dpbcnam, dxbcnam, dybcnam)
              deallocate (f3naMa, f3naMb)
              deallocate (f3naXa, f3naXb, f3naXc)
            end if ! if (mneigh .ne. 0)
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
        end subroutine Dassemble_vna_3c

! End Module
! ===========================================================================
        end module M_Dassemble_3c
