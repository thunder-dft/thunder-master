! copyright info:
!
!                             @Copyright 2012
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

! M_Drotations_PP
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
        module M_Drotations_PP
        use M_rotations
        use M_Drotations
        use M_species

! Type Declaration
! =========================================================================== 

! module procedures
        contains
        
! ===========================================================================
! Drotate_PP
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>      This subroutine handles rotations and derivatives of rotations 
!! specific for two-center pseudopotential interactions. 
!!
!!       This routine rotates a matrix from molecular to crystal coordinates.
!!
!! The variable eps is a 3x3 output of the subroutine epsilon.
!! The variable dmat is a 5x5 matrix rotating d-orbitals.
!! The variable pmat is a 3x3 matrix rotating p-orbitals.
!!
!! Here is the famous Ortega convention:
!! In the molecular coordinates we have atoms 1 and 2 (bondcharge) along the
!! z-axis; the third atom is in the XZ-plane.
!!
!! The labelling of the orbitals is as follows:
!!
!!   S-shell :                s
!!                            1
!!
!!   P-shell :           py   pz   px
!!                       1    2    3
!!
!!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!!                  1     2   3     4      5
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
        subroutine Drotate_PP (in1, in2, eps, Deps, norb_mu, norb_nu, matm,&
     &                            Dmatm, Dmatx)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: in1, in2
        integer, intent(in) :: norb_mu, norb_nu             !< size of block

        real, intent(in), dimension (3, 3, 3) :: Deps
        real, intent(in), dimension (3, 3) :: eps

! Note that all dmatm ar-            write (42,*) iatom, jatom, tm, dtm
!rays are vectors. That means the corresponding
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
        call Dmatrix_PP (in1, in2, norb_mu, norb_nu, matm, dmatm,  dmat,     &
     &                   pmat, Ddmat, Dpmat, Dmatx)

! Format Statements
! ===========================================================================
! None
        return
        end subroutine Drotate_PP


! ===========================================================================
! Dmatrix_PP
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine assembles components of the matrix derivative for the
! pseudopotential terms.
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
        subroutine Dmatrix_PP (in1, in2, norb_mu, norb_nu, matm, Dmatm, dmat,&
     &                         pmat, Ddmat, Dpmat, term)
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
          do jssh = 1, species(in2)%nssh_PP
            l2 =  species(in2)%shell_PP(jssh)%lssh
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
        end subroutine Dmatrix_PP

! End Module
! ===========================================================================
        end module M_Drotations_PP
