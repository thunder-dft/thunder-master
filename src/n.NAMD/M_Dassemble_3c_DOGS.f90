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

! /GLOBAL
        use M_assemble_blocks

! /SYSTEM
        use M_configuraciones
        use M_neighbors
        use M_rotations
        use M_Drotations

! /FDATA
        use M_Fdata_3c

! /ASSEMBLERS
        use M_assemble_3c

! /SOLVESH
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
        integer in1, in2, indna          !< species numbers
        integer interaction, isorp       !< which interaction and subtype
        integer imu, inu, iindex         !< indexing counters
        integer issh
        integer norb_mu, norb_nu         !< size of the block for the pair

        real distance_13, distance_23    !< distance from 3rd atom

        ! smoothing variables
        real rcutoff1_min, rcutoff2_min, rcutoff3_min  !< cutoffs for smoothing
        real rend1, rend2                !< smoothing end points
        real stinky, stinky13, stinky23  !< smoothing values
        real xsmooth                     !< for smoothing function

        ! smoothing derivative variables
        real dstinky13, dstinky23
        real, dimension (3) :: dstnA     !< smoothing value
        real, dimension (3) :: dstnB     !< smoothing value
        real, dimension (3) :: dstnC     !< smoothing value

        real dQ                          !< isorp charge
        real z                           !< distance between r1 and r2
        real x, cost                     !< dnabc and angle

        real, dimension (3, 3) :: eps    !< the epsilon matrix
        real, dimension (3, 3, 3) :: depsA  !< the Depsilon matrix for the bc
        real, dimension (3, 3, 3) :: depsB  !< the Depsilon matrix for the na

        real, dimension (3) :: r1, r2, rna  !< positions - iatom, jatom, ialpha
        real, dimension (3) :: r21, rnabc   !< vectors
        real, dimension (3) :: sighat    !< unit vector along r2 - r1
        real, dimension (3) :: rhat      !< unit vector along bc - rna
        real, dimension (3) :: rhatA1    !< unit vector along rna - r1
        real, dimension (3) :: rhatA2    !< unit vector along rna - r2

        real, dimension (3) :: amt, bmt

        real, dimension (:, :), allocatable :: dterm
        real, dimension (:, :, :), allocatable :: dpterm
        real, dimension (:, :), allocatable :: sterm
        real, dimension (:, :, :), allocatable :: spterm
        real, dimension (:, :), allocatable :: emnpl

        ! derivatives of three-center overlap-dipole ewald interactions
        real, dimension (:, :, :), allocatable :: demnplA
        real, dimension (:, :, :), allocatable :: demnplB
        real, dimension (:, :, :), allocatable :: demnplC

! bcnam = Hartree matrix in molecular coordinates
! bcnax = Hartree matrix in crystal coordinates
! dbcnam = derivative of Hartree matrix in molecular coordinates
! vdbcnam = vectorized derivative of Hartree matrix in molecular coordinates
! vdbcnax = vectorized derivative of Hartree matrix in crystal coordinates
        real, dimension (:, :), allocatable :: bcnam
        real, dimension (:, :), allocatable :: bcnax

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

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pdip_neighbors
        type(T_assemble_neighbors), pointer :: pdipole_z

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

        interface
          function smoother (z, rend, x)
            real smoother
            real, intent(in) :: z, rend, x
          end function smoother
        end interface

        interface
           real function Dsmoother (r, rend, x)
             real, intent (in) :: r
             real, intent (in) :: rend
             real, intent (in) :: x
           end function Dsmoother
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
!****************************************************************************
! Neutral atom piece:
! For now we just do the neutral atom interactions.
! So set isorp = 0 within this subroutine.
!
!              interaction    subtypes     index
!
!      bcna         1           0..9(max)   1..10
! ****************************************************************************
!
!****************************************************************************
! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          indna = s%atom(ialpha)%imass
          rna = s%atom(ialpha)%ratom

          ! cut some lengthy notation
          nullify (pfalpha)
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
              nullify (pfi, pfj)
              pfi=>s%forces(iatom); pfj=>s%forces(jatom)

              ! density matrix
              nullify (pdenmat, pRho_neighbors)
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
                sighat = r21/z
              end if

! ****************************************************************************
! Find rnabc = vector pointing from center of bondcharge to rna
! This gives us the distance dnabc.
              rnabc = rna - (r1 + r21/2.0d0)
              x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)

              ! unit vector in rnabc direction.
              if (x .lt. 1.0d-05) then
                rhat(1) = 0.0d0
                rhat(2) = 0.0d0
                rhat(3) = 0.0d0
              else
                rhat = rnabc/x
              end if
              cost = dot_product(sighat, rhat)
              call epsilon_function (rhat, sighat, eps)

! dera3 = depsA = deps/dratm in the 3-center system
! der13 = depsB = deps/dr1 in the 3-center system
              call Depsilon_3c (r1, r2, r21, z, rna, rnabc, eps, depsA, depsB)

! The first piece will be the force with respect to atom 3.
              if (x .gt. 1.0d-5) then
                amt = (sighat - cost*rhat)/x
              else
                amt = 0.0d0
              end if

! The second piece will be the force with respect to atom 1.
              bmt = (cost*sighat - rhat)/z

! Get the D-matrices from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_bcna
              isorp = 0

! Allocate block arrays
              allocate (bcnam (norb_mu,norb_nu)); bcnam = 0.00
              allocate (dpbcnam (norb_mu,norb_nu)); dpbcnam = 0.00
              allocate (dxbcnam (norb_mu,norb_nu)); dxbcnam = 0.00
              allocate (dybcnam (norb_mu,norb_nu)); dybcnam = 0.00
              call getDMEs_Fdata_3c (in1, in2, indna, interaction, isorp, x,  &
     &                               z, norb_mu, norb_nu, cost, rhat, sighat, &
     &                               bcnam, dpbcnam, dxbcnam, dybcnam)

! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.
              nullify (pFdata_bundle, pFdata_cell)
              pFdata_bundle => Fdata_bundle_3c(in1, in2, indna)
              pFdata_cell =>                                                  &
     &          pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(interaction,isorp,1))

              allocate (f3naMa(3,norb_mu,norb_nu)); f3naMa = 0.0d0
              allocate (f3naMb(3,norb_mu,norb_nu)); f3naMb = 0.0d0
              allocate (f3naXa(3,norb_mu,norb_nu)); f3naXa = 0.0d0
              allocate (f3naXb(3,norb_mu,norb_nu)); f3naXb = 0.0d0
              allocate (f3naXc(3,norb_mu,norb_nu)); f3naXc = 0.0d0
              ! loop over matrix elements
              do iindex = 1, pFdata_cell%nME
                imu = pFdata_cell%mu_3c(iindex)
                inu = pFdata_cell%nu_3c(iindex)

! Now recover f3naMa which is a three-dimensional array
                f3naMa(:,imu,inu) = rhat*dxbcnam(imu,inu) + amt*dpbcnam(imu,inu)


! Now recover f3naMb which is a three-dimensional array
                f3naMb(:,imu,inu) = - sighat*dybcnam(imu,inu)                 &
     &            + bmt*dpbcnam(imu,inu) - f3naMa(:,imu,inu)/2.0d0
              end do ! end loop over matrix elements
              nullify (pFdata_bundle, pFdata_cell)

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
              call Drotate (in1, in2, eps, depsA, norb_mu, norb_nu,           &
     &                      bcnam, f3naMa, f3naXa)

! Force on the neutral atom with respect to atom 1 (f3naMb).
              call Drotate (in1, in2, eps, depsB, norb_mu, norb_nu,           &
     &                      bcnam, f3naMb, f3naXb)

! Make things force-like and determine f3naXc, whcih is found from Newtons Laws:
!             f3naXa = - f3naXa
!             f3naXb = - f3naXb
              f3naXc = - f3naXa - f3naXb

              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  pfalpha%f3naa = pfalpha%f3naa                               &
      &             - pRho_neighbors%block(imu,inu)*f3naXa(:,imu,inu)*P_eq2
                  pfi%f3nab = pfi%f3nab                                       &
      &             - pRho_neighbors%block(imu,inu)*f3naXb(:,imu,inu)*P_eq2
                  pfj%f3nac = pfj%f3nac                                       &
      &             - pRho_neighbors%block(imu,inu)*f3naXc(:,imu,inu)*P_eq2
                end do
              end do

              deallocate (bcnam, dpbcnam, dxbcnam, dybcnam)
              deallocate (f3naMa, f3naMb)
              deallocate (f3naXa, f3naXb, f3naXc)
              nullify (pfi, pfj)
              nullify (pdenmat, pRho_neighbors)
            end if ! if (mneigh .ne. 0)
          end do ! end loop over neighbors
          nullify (pfalpha)
        end do ! end loop over atoms

! ****************************************************************************
! DO CHARGE CASES HERE.
! Charged atom interactions are assembled here.
! So set isorp = 1, nssh within this subroutine.
!
!              interaction    subtypes     index
!
!      bcna         1           0..9(max)   1..10

! ****************************************************************************
! There is a lot here because of the smoothing nonsense.
! ****************************************************************************
! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          indna = s%atom(ialpha)%imass
          rna = s%atom(ialpha)%ratom

          ! cut some lengthy notation
          nullify (pfalpha)
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
              nullify (pfi, pfj)
              pfi=>s%forces(iatom); pfj=>s%forces(jatom)

              ! density matrix
              nullify (pdenmat, pRho_neighbors)
              pdenmat=>s%denmat(iatom); pRho_neighbors=>pdenmat%neighbors(mneigh)

              ! interactions for smoothing part
              nullify (poverlap, pS_neighbors)
              poverlap=>s%overlap(iatom); pS_neighbors=>poverlap%neighbors(mneigh)
              nullify (pdipole_z, pdip_neighbors)
              pdipole_z=>s%dipole_z(iatom); pdip_neighbors=>pdipole_z%neighbors(mneigh)

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
              rnabc = rna - (r1 + r21/2.0d0)
              x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)

! Find other distances -
              distance_13 = distance (rna, r1)
              distance_23 = distance (rna, r2)

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
              call Depsilon_3c (r1, r2, r21, z, rna, rnabc, eps, depsA, depsB)

! Find the unit vector in rna-r1 direction.
              if (distance_13 .gt. 1.0d-05) then
                rhatA1(:) = (rna(:) - r1(:))/distance_13
              end if

! Find the unit vector in rna-r2 direction.
              if (distance_23 .gt. 1.0d-05) then
                rhatA2(:) = (rna(:) - r2(:))/distance_23
              end if

! The first piece will be the force with respect to atom 3.
              if (x .gt. 1.0d-5) then
                amt = (sighat - cost*rhat)/x
              else
                amt = 0.0d0
              end if

! The second piece will be the force with respect to atom 1.
              bmt = (cost*sighat - rhat)/z

! ****************************************************************************
! BEGIN
! Monopole and Dipole interactions for "smoothing"
! ****************************************************************************
! Find the smoothing quantity - here we calculate the long-range effective
! monopole.  This term is included so that we obtain no discontinuities when
! atoms leave or enter the rcutoff_1 + rcutoff_2 range criteria.
! Therefore, "close" two-center interactions are exact, while more distant
! two-center integrals go to effective monopoles.  The monopoles are effective
! in the sense that the two atoms in the matrix element, each has a different
! charge.  Since they are separated, this gives a monopole contribution at long
! range.

! The smoothing function is found by calling smoother(r,rbegin,rend).
! We define our final matrix element answer as
! smoother(r)*exact_piece + (1 - smoother(r))*longrange.  The distance r is the
! distance of the third center from the "effective" center of the bondcharge.
! The effective center of the bondcharge is (d + rc1 - rc2)/2 from r1 except in
! weird cases (see below). The distance rbegin is the distance at which we
! include only exact answers and do not smooth. The distance rend is the
! distance past which smooth(r) is zero, so that the result is long-range only.
! We skipped self-interaction terms.
              allocate (sterm (norb_mu, norb_nu)); sterm = 0.0d0
              allocate (spterm (3, norb_mu, norb_nu)); spterm = 0.0d0

              allocate (dterm (norb_mu, norb_nu)); dterm = 0.0d0
              allocate (dpterm (3, norb_mu, norb_nu)); dpterm = 0.0d0

              allocate (emnpl (norb_mu, norb_nu)); emnpl = 0.0d0
              allocate (DemnplA (3, norb_mu, norb_nu)); DemnplA = 0.0d0
              allocate (DemnplB (3, norb_mu, norb_nu)); DemnplB = 0.0d0
              allocate (DemnplC (3, norb_mu, norb_nu)); DemnplC = 0.0d0

              xsmooth = 0.8d0  ! parameter for smoothing
              rcutoff1_min = 99.0d0
              do issh = 1, species(in1)%nssh
                rcutoff1_min = min(rcutoff1_min, species(in1)%shell(issh)%rcutoffA)
              end do

              rcutoff2_min = 99.0d0
              do issh = 1, species(in2)%nssh
                rcutoff2_min = min(rcutoff2_min, species(in2)%shell(issh)%rcutoffA)
              end do

              rcutoff3_min = 99.0d0
              do issh = 1, species(indna)%nssh
                rcutoff3_min = min(rcutoff3_min, species(indna)%shell(issh)%rcutoffA)
              end do

              rend1 = rcutoff1_min + rcutoff3_min
              stinky13 = smoother (distance_13, rend1, xsmooth)
              dstinky13 = Dsmoother (distance_13, rend1, xsmooth)

              rend2 = rcutoff2_min + rcutoff3_min
              stinky23 = smoother (distance_23, rend2, xsmooth)
              dstinky23 = Dsmoother (distance_23, rend2, xsmooth)

              stinky = stinky13*stinky23

              dstnB = - dstinky13*stinky23*rhatA1
              dstnC = - stinky13*dstinky23*rhatA2
              dstnA = - dstnB - dstnC  ! using Newton's third law

! Set value for emnpl
              sterm = pS_neighbors%block/2.0d0
              spterm = pS_neighbors%Dblock/2.0d0
              dterm = pdip_neighbors%block/z
              dpterm = pdip_neighbors%Dblock/z
              emnpl = (sterm - dterm)/distance_13 + (sterm + dterm)/distance_23

! Now the derivatives
! A = 3, B = 1, C = 2 <B|A|C> This is NOT force-like.

! demnplA: The net effect is a minus sign.
! There are three minus signs in all - one belonging to the sighat term,
! one belonging to the 1/y term and a minus subtraction error.
! No minus signs are "force-like". This is not a force like derivative,
! but a regular derivative d/dr1.
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  dpterm(:,imu,inu) = dpterm(:,imu,inu) + dterm(imu,inu)*sighat(:)/z
                  demnplA(:,imu,inu) =                                        &
     &              - (sterm(imu,inu)                                         &
     &                 - dterm(imu,inu))*rhatA1(:)/distance_13**2             &
     &              - (sterm(imu,inu)                                         &
     &                 + dterm(imu,inu))*rhatA2(:)/distance_23**2
                  demnplB(:,imu,inu) =                                        &
                    + (sterm(imu,inu)                                         &
     &                 - dterm(imu,inu))*rhatA1(:)/distance_13**2             &
     &              + (spterm(:,imu,inu)                                      &
     &                  - dpterm(:,imu,inu))/distance_13                      &
     &              + (spterm(:,imu,inu)                                      &
     &                 + dpterm(:,imu,inu))/distance_23
                end do
              end do
              demnplC = - demnplA - demnplB
! ****************************************************************************
! END "smoothing" preparations - now build the forces.
! ****************************************************************************

! Get the D-matrices from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
              interaction = P_bcna
              do isorp = 1, species(indna)%nssh

! Allocate block arrays
                allocate (bcnam (norb_mu,norb_nu)); bcnam = 0.00
                allocate (bcnax (norb_mu,norb_nu)); bcnax = 0.00
                allocate (dpbcnam (norb_mu,norb_nu)); dpbcnam = 0.00
                allocate (dxbcnam (norb_mu,norb_nu)); dxbcnam = 0.00
                allocate (dybcnam (norb_mu,norb_nu)); dybcnam = 0.00
                call getDMEs_Fdata_3c (in1, in2, indna, interaction, isorp, x,&
     &                                 z, norb_mu, norb_nu, cost, rhat,       &
     &                                 sighat, bcnam, dpbcnam, dxbcnam,       &
     &                                 dybcnam)

                ! Rotate into crystal coordinates
                call rotate (in1, in2, eps, norb_mu, norb_nu, bcnam, bcnax)

                allocate (f3naMa(3,norb_mu,norb_nu)); f3naMa = 0.0d0
                allocate (f3naMb(3,norb_mu,norb_nu)); f3naMb = 0.0d0
                allocate (f3naXa(3,norb_mu,norb_nu)); f3naXa = 0.0d0
                allocate (f3naXb(3,norb_mu,norb_nu)); f3naXb = 0.0d0
                allocate (f3naXc(3,norb_mu,norb_nu)); f3naXc = 0.0d0

! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.
                nullify (pFdata_bundle, pFdata_cell)
                pFdata_bundle => Fdata_bundle_3c(in1, in2, indna)
                pFdata_cell =>                                                &
     &            pFdata_bundle%Fdata_cell_3c(pFdata_bundle%index_3c(interaction,isorp,1))

                ! loop over matrix elements
                do iindex = 1, pFdata_cell%nME
                  imu = pFdata_cell%mu_3c(iindex)
                  inu = pFdata_cell%nu_3c(iindex)

! Now recover f3naMa which is a three-dimensional array
                  f3naMa(:,imu,inu) = rhat*dxbcnam(imu,inu) + amt*dpbcnam(imu,inu)

! Now recover f3naMb which is a three-dimensional array
                  f3naMb(:,imu,inu) = - sighat*dybcnam(imu,inu)               &
     &              + bmt*dpbcnam(imu,inu) - f3naMa(:,imu,inu)/2.0d0
                end do ! end loop over matrix elements
                nullify (pFdata_bundle, pFdata_cell)

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
                call Drotate (in1, in2, eps, depsA, norb_mu, norb_nu,         &
     &                        bcnam, f3naMa, f3naXa)

! Force on the neutral atom with respect to atom 1 (f3naMb).
                call Drotate (in1, in2, eps, depsB, norb_mu, norb_nu,         &
     &                        bcnam, f3naMb, f3naXb)

! Make things force-like and determine f3naXc, whcih is found from Newtons Laws:
!               f3naXa = - f3naXa
!               f3naXb = - f3naXb
                f3naXc = - f3naXa - f3naXb

                dQ = s%atom(ialpha)%shell(isorp)%dQ
                do inu = 1, norb_nu
                  do imu = 1, norb_mu
                    pfalpha%f3naa = pfalpha%f3naa                             &
      &               - P_eq2*dQ*pRho_neighbors%block(imu,inu)                &
      &                *(stinky*f3naXa(:,imu,inu) + dstnA(:)*bcnax(imu,inu)   &
      &                  - dstnA(:)*emnpl(imu,inu) + (1.0d0 - stinky)*demnplA(:,imu,inu))
                    pfi%f3nab = pfi%f3nab                                     &
      &               - P_eq2*dQ*pRho_neighbors%block(imu,inu)                &
      &                *(stinky*f3naXb(:,imu,inu) + dstnB(:)*bcnax(imu,inu)   &
      &                  - dstnB(:)*emnpl(imu,inu) + (1.0d0 - stinky)*demnplB(:,imu,inu))
                    pfj%f3nac = pfj%f3nac                                     &
      &               - P_eq2*dQ*pRho_neighbors%block(imu,inu)                &
      &                *(stinky*f3naXc(:,imu,inu) + dstnB(:)*bcnax(imu,inu)   &
      &                  - dstnC(:)*emnpl(imu,inu) + (1.0d0 - stinky)*demnplC(:,imu,inu))
                  end do
                end do
                deallocate (bcnam, bcnax, dpbcnam, dxbcnam, dybcnam)
                deallocate (f3naMa, f3naMb)
                deallocate (f3naXa, f3naXb, f3naXc)
              end do ! end loop over isorp
              nullify (pfi, pfj)
              nullify (pdenmat, pRho_neighbors)
              nullify (poverlap, pS_neighbors)
              nullify (pdipole_z, pdip_neighbors)
            end if ! if (mneigh .ne. 0)
            deallocate (sterm, spterm, dterm, dpterm, emnpl)
            deallocate (DemnplA, DemnplB, DemnplC)
          end do ! end loop over neighbors
          nullify (pfalpha)
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
