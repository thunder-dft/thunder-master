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

! rotate_PP.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>      This subroutine handles rotations specific for two-center 
!! pseudopotential interactions. 
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
        subroutine rotate_PP (in1, in2, eps, norb_mu, norb_nu, mmatrix,      &
     &                        xmatrix)
        use M_rotations
        use M_species
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input:
        integer, intent(in) :: in1, in2                     !< species pair
        integer, intent(in) :: norb_mu, norb_nu             !< size of block
 
        real, intent(in), dimension (3, 3) :: eps           !< epsilon matrix
        real, intent(in), dimension (norb_mu, norb_nu) :: mmatrix !< matrix in molecular coordinates
 
! Output:
        real, intent(out), dimension (norb_mu, norb_nu) :: xmatrix  !< matrix in crystal coordinates 
 
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
          do jssh = 1, species(in2)%nssh_PP
            l2 = species(in2)%shell_PP(jssh)%lssh
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
        end subroutine rotate_PP
