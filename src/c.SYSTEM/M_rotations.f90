! copyright info:
!
!                             @Copyright 2013
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

! M_rotations
! Module Description
! ===========================================================================
!      This module determines the (two-center) matrix elements (mu,nu) in a
! matrix form.  First, we read the data from Fdata_2c, but we get the data
! in a one-dimensional array of non-zero matrix elements. Next, we need to
! interpolate the data for the input distance between the two species. After
! this, we then recover the data to get the data in matrix form.  The matrix
! is initially in molecular coordinates; thus, we need to rotate it to crystal
! coordinates. It contains the following subroutines within the module:
!
!       rotate.f90 - rotates the matrix to crystal coordinates.
!       chooser.f90 - choose left or right matrices for the rotation - L*A*R
!       twister.f90 - prepare the matrices for the rotation
!
! ===========================================================================
! =========================================================================== 
        module M_rotations
        use M_species

! Type Declaration
! =========================================================================== 

! module procedures
        contains
        

! ===========================================================================
! rotate.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine rotates a matrix from molecular to crystal coordinates.
!
! The variable eps is a 3x3 output of the subroutine epsilon.
! The variable dmat is a 5x5 matrix rotating d-orbitals.
! The variable pmat is a 3x3 matrix rotating p-orbitals.
!
! Here is the famous Ortega convention:
! In the molecular coordinates we have atoms 1 and 2 (bondcharge) along the
! z-axis; the third atom is in the XZ-plane.
!
! The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            1
!
!   P-shell :           py   pz   px
!                       1    2    3
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                  1     2   3     4      5
! ===========================================================================
! Code written by:
! James P. Lewis
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
        subroutine rotate (in1, in3, eps, norb_mu, norb_nu, mmatrix, xmatrix)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input:
        integer, intent(in) :: in1, in3                     ! species pair
        integer, intent(in) :: norb_mu, norb_nu             ! size of block
 
        real, intent(in), dimension (3, 3) :: eps           ! epsilon matrix
        real, intent(in), dimension (norb_mu, norb_nu) :: mmatrix
 
! Output:
        real, intent(out), dimension (norb_mu, norb_nu) :: xmatrix
 
! Local Parameters and Data Declaration
! ===========================================================================
! None
 
! Local Variable Declaration and Description
! ===========================================================================
        integer issh
        integer jssh
        integer k1, k2
        integer n1, l1, m1
        integer n2, l2, m2
 
        real, dimension (5, 5) :: dmat
        real, dimension (5, 5) :: left
        real, dimension (3, 3) :: pmat
        real, dimension (5, 5) :: right 

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Set up the matrices to rotate p and d orbitals.
        call twister (eps, dmat, pmat)

        xmatrix = 0.0d0
        n1 = 0
        do issh = 1, species(in1)%nssh
          l1 = species(in1)%shell(issh)%lssh
          call chooser (l1, dmat, pmat, left)
 
          n2 = 0
          do jssh = 1, species(in3)%nssh
            l2 = species(in3)%shell(jssh)%lssh
            call chooser (l2, dmat, pmat, right)
 
! Compute
            do m2 = 1, 2*l2 + 1
              do m1 = 1, 2*l1 + 1
                do k2 = 1, 2*l2 + 1
                  do k1 = 1, 2*l1 + 1
                    xmatrix(n1+k1,n2+k2) = xmatrix(n1+k1,n2+k2)              &
     &                         + left(k1,m1)*mmatrix(n1+m1,n2+m2)*right(k2,m2)
                  end do
                end do
              end do
            end do
            n2 = n2 + 2*l2 + 1
          end do
          n1 = n1 + 2*l1 + 1
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
        end subroutine rotate
       

! ===========================================================================
! twister.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine prepares the D matrices for a given geometry of the
!! matrix element. The set of vectors stored in eps(3,3) are needed.
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
!!                            0
!!
!!   P-shell :           py   pz   px
!!                       -1   0    +1
!!
!!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!!                 -2    -1   0    +1     +2
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
! Program Declaration
! ===========================================================================
        subroutine twister (eps, dmat, pmat)
        implicit none

        include '../include/constants.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent(in), dimension (3, 3) :: eps    ! epsilon matrix
 
! Output
        real, intent(out), dimension (5, 5) :: dmat  ! d-orbital type matrix
        real, intent(out), dimension (3, 3) :: pmat  ! p-orbital type matrix
 
! Local Parameters and Data Declaration
! ===========================================================================
! None
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer jx
        integer ix
        
        real amat_term
        real xlambda11
        real xlambda12
        real xlambda13
        real xlambda32
        real xlambda33 

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Nothing for S orbitals
! Set the P matrices: eps - x,y,z => pmat - y,z,x
        pmat(1,1) = eps(2,2) 
        pmat(1,2) = eps(2,3) 
        pmat(1,3) = eps(2,1)
  
        pmat(2,1) = eps(3,2) 
        pmat(2,2) = eps(3,3)
        pmat(2,3) = eps(3,1)
         
        pmat(3,1) = eps(1,2)
        pmat(3,2) = eps(1,3)
        pmat(3,3) = eps(1,1)

        if (.not. have_dorbitals) return

! Set the lambda matrices:
        do imu = 1, 5
          xlambda11 = 0.0d0
          xlambda12 = 0.0d0
          xlambda13 = 0.0d0
          xlambda32 = 0.0d0
          xlambda33 = 0.0d0
          do jx = 1, 3
            do ix = 1, 3
              if (amat(ix,jx,imu) .ne. 0.0d0) then
                amat_term = amat(ix,jx,imu)
                xlambda11 = xlambda11 + amat_term*eps(jx,1)*eps(ix,1)
                xlambda12 = xlambda12 + amat_term*eps(jx,1)*eps(ix,2)
                xlambda13 = xlambda13 + amat_term*eps(jx,1)*eps(ix,3)
                xlambda32 = xlambda32 + amat_term*eps(jx,3)*eps(ix,2)
                xlambda33 = xlambda33 + amat_term*eps(jx,3)*eps(ix,3)
              end if
            end do
          end do
         
! Set the D matrices:
          dmat(imu,1) = 2.0d0*xlambda12
          dmat(imu,2) = 2.0d0*xlambda32
          dmat(imu,3) = sqrt(3.0d0)*xlambda33
          dmat(imu,4) = 2.0d0*xlambda13
          dmat(imu,5) = 2.0d0*xlambda11 + xlambda33
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
        end subroutine twister


! ===========================================================================
! chooser
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine prepares the L (left) and R (right) matrices for the
! rotation of L*Matrix*R.
!
! ===========================================================================
! Original code from 
! Alex Demkov.
! Code rewritten by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
! 
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine chooser (lqn, dmat, pmat, rmatrix)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input:
        integer, intent(in) :: lqn  !< l quantum number 
 
        real, intent(in), dimension (5, 5) :: dmat  !< d-orbital type matrix 
        real, intent(in), dimension (3, 3) :: pmat  !< p-orbital type matrix
 
! Output:
        real, intent(out), dimension (5, 5) :: rmatrix  !< rotation output 
 
! Local Parameters and Data Declaration
! ===========================================================================
! None
 
! Local Variable Declaration and Description
! ===========================================================================
! None 
 
! Procedure
! ===========================================================================
! Initialize the rotation matrix
        rmatrix = 0.0d0
        if (lqn .eq. 0) then
          rmatrix(1,1) = 1.0d0
        else if (lqn .eq. 1) then
          rmatrix(1:3,1:3) = pmat
        else if (lqn .eq. 2) then
          rmatrix = dmat
        end if
 
! Format Statements
! ===========================================================================
! None

! End Subroutine
! =========================================================================== 
        return
        end subroutine chooser


! End Module
! ===========================================================================
        end module M_rotations
