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

! M_Drotations
! Module Description
! ===========================================================================
!>      This module determines the (two-center) matrix elements (mu,nu) in a 
!! matrix form.  First, we read the data from Fdata_2c, but we get the data
!! in a one-dimensional array of non-zero matrix elements. Next, we need to 
!! interpolate the data for the input distance between the two species. After
!! this, we then recover the data to get the data in matrix form.  The matrix
!! is initially in molecular coordinates; thus, we need to rotate it to crystal
!! coordinates. It contains the following subroutines within the module:
!! 
!!       Drotate.f90 - rotates the matrix derivatives to crystal coordinates
!!       Dchooser.f90 - choose left or right derivative matrices for rotation
!!       Dtwister.f90 - prepare the derivative matrices for rotation
!!       Dmatrix.f90 - make the derivative matrix
!
! ===========================================================================
! =========================================================================== 
        module M_Drotations
        use M_rotations
        use M_species

! Type Declaration
! =========================================================================== 
! None

! module procedures
        contains
        
! ===========================================================================
! Drotate
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine rotates a matrix derivative from molecular coordinates
! to crystal coordinates.
!
! The variable matm is the matrix-box in molecular coordinates
! The variable dmatm is a matrix which is the derivative of the matrix-box
! in molecular coordinates.
! The output is dmatx which is the derivative of the matrix-box in crystal
! coordinates.
!
! In the molecular coordinates we have atoms 1 and 2 (bondcharge) along the
! z-axis; the third atom is in the XZ-plane.
! The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            1
!
!   P-shell :           py   pz   px
!                       -1   0    1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1    0    1      2
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Drotate (in1, in2, eps, Deps, norb_mu, norb_nu, matm,     &
     &                      Dmatm, Dmatx)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: in1, in2
        integer, intent(in) :: norb_mu, norb_nu             !< size of block

        real, intent(in), dimension (3, 3, 3) :: Deps
        real, intent(in), dimension (3, 3) :: eps

! Note that all dmatm arrays are vectors. That means the corresponding
! scalar derivatives were multiplied by -eta, which is dD/dr1, after the
! interpolation.
        real, intent(in), dimension (3, norb_mu, norb_nu)  :: Dmatm
        real, intent(in), dimension (norb_mu, norb_nu) :: matm

! Output
        real, intent(out), dimension (3, norb_mu, norb_nu) :: Dmatx

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        real pmat (3, 3)       ! p-states matrix part
        real dmat (5, 5)       ! d-states matrix part

        real Dpmat (3, 3, 3)   ! derivative of p-states matrix part
        real Ddmat (3, 5, 5)   ! derivative of d-states matrix part

! Procedure
! ===========================================================================
! Set rotational matrices and their respective derivatives.
        call twister (eps, dmat, pmat)
        call Dtwister (eps, Deps, Ddmat, Dpmat)

! Now we have to do the cases A, B and C, which  are added in makeDmat
! page 1 10/1/98 notes.
! term A: dAleft/dr*Aright*matrix
! term B: Aleft*dAright/dr*matrix
! term C Aleft*Aright*dmatrix/dr
        call Dmatrix (in1, in2, norb_mu, norb_nu, matm, dmatm, dmat, pmat,  &
     &                Ddmat, Dpmat, Dmatx)

! Format Statements
! ===========================================================================
! None
        return
        end subroutine Drotate


! ===========================================================================
! Dtwister
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine prepares the derivatives of the D matrices for a given
! geometry of the matrix element.  This routine needs to know the variable
! epsilon.  This epsilon is either eps2 or eps3, depending on what matrix
! element we are dealing with (e.g. T is a 2C, while Vna is a 3C case).
!
! Input:
! The variable eps is a 3x3 output of the subroutine epsilon.
! The variable Deps is a 3x3X3 output of either Deps2cent, or Deps3cent.

! Output:
! The variable ddmat is a 3x5x5 matrix  derivative rotating d-orbitals.
! The variable dpmat is a 3x3x3 matrix  derivative rotating p-orbitals.
!
! Here is the famous Ortega convention:
! In the molecular coordinates we have atoms 1 and 2 (bondcharge) along the
! z-axis; the third atom is in the XZ-plane.
!
! The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dtwister (eps, Deps, Ddmat, Dpmat)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent(in), dimension (3, 3) :: eps
        real, intent(in), dimension (3, 3, 3) :: Deps

! Output
        real, intent(out), dimension (3, 5, 5) :: Ddmat
        real, intent(out), dimension (3, 3, 3) :: Dpmat

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer ix, jx, kx

        real aterm12, aterm32, aterm33, aterm13, aterm11

        real amat_term

! Procedure
! ===========================================================================
! Nothing for S orbitals.
! Set the dP/drk matrices: eps: x,y,z; pmat: y,z,x
        Dpmat(:,1,1) = Deps(:,2,2)
        Dpmat(:,1,2) = Deps(:,2,3)
        Dpmat(:,1,3) = Deps(:,2,1)

        Dpmat(:,2,1) = Deps(:,3,2)
        Dpmat(:,2,2) = Deps(:,3,3)
        Dpmat(:,2,3) = Deps(:,3,1)

        Dpmat(:,3,1) = Deps(:,1,2)
        Dpmat(:,3,2) = Deps(:,1,3)
        Dpmat(:,3,3) = Deps(:,1,1)

! ***************************************************************************
! Set the dD/dr matrices according to the general formula
! (see p.3 notes 9/29/98):
!
! dD(mu|m)/dr_k ~ dLAMBDA(mu)_[a,b]/d_r = SUM_(i,j)
!
! a(mu)_[i,j]((Deps(i,a)/dr_k)*eps(j,b) + eps(i,a)*(Deps(j,b)/dr_k) )
!
! and m defines [a,b], e.g. m = -2 requires [1,2] (see p. 3 notes 9/29/98),
! Note the formula is "mu-independent" for a given m!
! ***************************************************************************
        do imu = 1, 5
          do kx = 1, 3
            aterm12 = 0.0d0
            aterm32 = 0.0d0
            aterm33 = 0.0d0
            aterm13 = 0.0d0
            aterm11 = 0.0d0
            do ix = 1, 3
              do jx = 1, 3
                if (amat(jx,ix,imu) .ne. 0.0d0) then
                  amat_term = amat(jx,ix,imu)
                  aterm12 = aterm12 + amat_term*(Deps(kx,ix,1)*eps(jx,2) + eps(ix,1)*Deps(kx,jx,2))
                  aterm32 = aterm32 + amat_term*(Deps(kx,ix,3)*eps(jx,2) + eps(ix,3)*Deps(kx,jx,2))
                  aterm33 = aterm33 + amat_term*(Deps(kx,ix,3)*eps(jx,3) + eps(ix,3)*Deps(kx,jx,3))
                  aterm13 = aterm13 + amat_term*(Deps(kx,ix,1)*eps(jx,3) + eps(ix,1)*Deps(kx,jx,3))
                  aterm11 = aterm11 + amat_term*(Deps(kx,ix,1)*eps(jx,1) + eps(ix,1)*Deps(kx,jx,1))
                end if
              end do
            end do

! Creating derivatives of D(mu|-2) ==> D(mu|1), LAMBDA_[12]
            Ddmat(kx,imu,1) = 2.0d0*aterm12

! Creating derivatives of D(mu|-1) ==> D(mu|2)  , LAMBDA_[32]
            Ddmat(kx,imu,2) = 2.0d0*aterm32

! Creating derivatives of D(mu|0)  ==> D(mu|3), LAMBDA_[33]
            Ddmat(kx,imu,3) = sqrt(3.0d0)*aterm33

! Creating derivatives of D(mu|1)  ==> D(mu|4), LAMBDA_[13]
            Ddmat(kx,imu,4) = 2.0d0*aterm13

! Creating derivatives of D(mu|2)  ==> D(mu|5), 2*LAMBDA_[11]+LAMBDA_[33]
            Ddmat(kx,imu,5) = 2.0d0*aterm11 + aterm33
          end do
        end do

! Format Statements
! ===========================================================================
! None
        return
        end subroutine Dtwister


! ===========================================================================
! Dchooser
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine obtains the derivatives of the left and right matrices
! for the equation L*M*R.
!
! The variable ddmat is a 3x5x5 derivative matrix rotating d-orbitals.
! The variable dpmat is a 3x3x3 derivative matrix rotating p-orbitals.
! The variable dmatrix is a derivative of a rotation matrix for l-value.
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dchooser (l, Ddmat, Dpmat, Dmatrix)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l

        real, intent (in), dimension(3, 5, 5) :: Ddmat
        real, intent (in), dimension(3, 3, 3) :: Dpmat

! Output
        real, intent (out), dimension(3, 5, 5) :: Dmatrix

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize to zero
        Dmatrix = 0.0d0

! Choose pmatrix or dmatrix.
        if (l .eq. 0) then
        else if (l .eq. 1) then
          Dmatrix(:,1:3,1:3) = Dpmat
        else if (l .eq. 2) then
          Dmatrix = Ddmat
        end if

! Format Statements
! ===========================================================================
! None
        return
        end subroutine Dchooser


! ===========================================================================
! Dmatrix
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine assembles components of the matrix derivative.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dmatrix (in1, in2, norb_mu, norb_nu, matm, Dmatm, dmat,   &
     &                      pmat, Ddmat, Dpmat, term)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: in1, in2
        integer, intent(in) :: norb_mu, norb_nu             !< size of block

! Note that all dmatm arrays are vectors. That means the corresponding
! scalar derivatives were multiplied by -eta, which is dD/dr1, after the
! interpolation.
        real, intent(in), dimension (3, norb_mu, norb_nu)  :: Dmatm
        real, intent(in), dimension (norb_mu, norb_nu) :: matm

        real, intent(in), dimension (3, 3) :: pmat       ! p-states matrix part
        real, intent(in), dimension (5, 5) :: dmat       ! d-states matrix part

        ! derivative of p-states matrix part
        real, intent(in), dimension (3, 3, 3) :: Dpmat

        ! derivative of d-states matrix part
        real, intent(in), dimension (3, 5, 5) ::  Ddmat

! Output
        !print *, norb_mu, norb_nu, "norb_mu norb_nu for term" 
        real, intent(out), dimension (3, norb_mu, norb_nu) :: term
         
! Local Variable Declaration and Description
! ===========================================================================
        integer issh, jssh
        integer ix

        integer k1, k2
        integer l1, l2
        integer m1, m2
        integer n1, n2

        real left (5, 5)
        real right (5, 5)

        real Dleft (3, 5, 5)
        real Dright (3, 5, 5)

! Procedure
! ===========================================================================
        n1 = 0
        do issh = 1, species(in1)%nssh
          l1 =  species(in1)%shell(issh)%lssh
          call chooser (l1, dmat, pmat, left)
          call Dchooser (l1, Ddmat, Dpmat, Dleft)

          n2 = 0
          do jssh = 1, species(in2)%nssh
            l2 =  species(in2)%shell(jssh)%lssh
            call chooser (l2, dmat, pmat, right)
            call Dchooser (l2, Ddmat, Dpmat, Dright)

! Initialize
            do k2 = 1, 2*l2 + 1
              do k1 = 1, 2*l1 + 1
                term(:,n1+k1,n2+k2) = 0.0d0
              end do
            end do

! Calculate
            do m2 = 1, 2*l2 + 1
              do k2 = 1, 2*l2 + 1
                do m1 = 1, 2*l1 + 1
                  do k1 = 1, 2*l1 + 1
                    do ix = 1, 3

! A, B, C terms in order
                      term(ix,n1+k1,n2+k2) = term(ix,n1+k1,n2+k2)            &
     &                  + Dleft(ix,k1,m1)*right(k2,m2)*matm(n1+m1,n2+m2)     &
     &                  + left(k1,m1)*Dright(ix,k2,m2)*matm(n1+m1,n2+m2)     &
     &                  + left(k1,m1)*right(k2,m2)*Dmatm(ix,n1+m1,n2+m2)
                    end do
                  end do
                end do
              end do
            end do

            n2 = n2 + 2*l2 + 1
          end do
          n1 = n1 + 2*l1 + 1
        end do

! Format Statements
! ===========================================================================
! None
        return
        end subroutine Dmatrix

! End Module
! ===========================================================================
        end module M_Drotations
