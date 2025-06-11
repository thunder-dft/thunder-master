! copyright info:
!
!                             @Copyright 2024
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
! Computer Network Information Center, Chinese Academy of Sciences    &
!     & University of Chinese Academy of Sciences - Runfeng Jin
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman
!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
!
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
!
! Runfeng Jin
! Computer Network Information Center, Chinese Academy of Sciences
! Beijing, China
!
! E-mail: jsfaraway@gmail.com
! M_Dassemble_ewald_DOGS_kokkos_interface
! ===========================================================================
! Program Description
! ===========================================================================
!>       This is a file for interacing with Kokkos C.
!!
!! It contains the following modules:
!!
!!       kokkos_data - the structure for transfering data between C and FORTRAN
!!       kokkos_interface - the C function interface called by FORTRAN code
!!
!
! ===========================================================================
! Module Declaration
! ============================================================================

module kokkos_data
        use :: precision_kinds
        use :: iso_c_binding
        implicit none
        type :: Gamma_data
                integer(cik) :: ig1mx
                integer(cik) :: ig2mx
                integer(cik) :: ig3mx
                real(crk) :: kappa
                type(c_ptr) :: v_c_ratom
                type(c_ptr) :: v_c_ewald
                type(c_ptr) :: v_c_dewald
                type(c_ptr) :: v_c_atom_ewald
                type(c_ptr) :: v_c_g123
                type(c_ptr) :: v_c_Qall
                real(crk) :: volume
                integer(cik) :: natoms
                integer(cik) :: icluster
                type(c_ptr) :: v_c_a123vec
                integer(cik), DIMENSION(3) :: c_il123mx
        end type Gamma_data

end module kokkos_data


module kokkos_interface
        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env
        use :: flcl_mod
        implicit none
        interface
          subroutine f_print_ratom( ratom ) &
          & bind(c, name='c_print_ratom')
                  use, intrinsic :: iso_c_binding
                  use :: flcl_mod
                  type(c_ptr), intent(in) :: ratom
          end subroutine f_print_ratom

          subroutine fgamma(gamma_data_f) &
          & bind(c, name='c_gamma')
                  use :: flcl_mod
                  use :: kokkos_data, only: Gamma_data
                  implicit none
                  type(Gamma_data), intent(inout) :: gamma_data_f
          end subroutine fgamma
        end interface

contains

subroutine gamma1(ig1mx, ig2mx, ig3mx, kappa,  v_c_ratom, v_c_ewald, v_c_dewald,  v_c_atom_ewald, v_c_g123, v_c_Qall, volume, natoms, icluster, v_c_a123vec, c_il123mx)
        use :: kokkos_data, only: Gamma_data
        use :: precision_kinds
        implicit none
        integer, intent(in) :: ig1mx, ig2mx, ig3mx, natoms, icluster
        real, intent(in) :: kappa, volume
        type(gamma_view2d), intent(in) :: v_c_ratom
        type(gamma_view2d), intent(inout) :: v_c_ewald
        type(gamma_view3d), intent(inout) :: v_c_dewald
        type(gamma_view2d), intent(inout) :: v_c_atom_ewald
        type(gamma_view2d), intent(in) :: v_c_g123
        type(gamma_view2d), intent(in) :: v_c_Qall
        type(gamma_view2d), intent(in) :: v_c_a123vec
        type(Gamma_data) :: gamma_data_f
        integer(cik), DIMENSION(3)  :: c_il123mx

        gamma_data_f%kappa = kappa
        gamma_data_f%volume = volume
        gamma_data_f%ig1mx = ig1mx
        gamma_data_f%ig2mx = ig2mx
        gamma_data_f%ig3mx = ig3mx
        gamma_data_f%v_c_ratom = v_c_ratom%view%ptr()
        gamma_data_f%v_c_ewald = v_c_ewald%view%ptr()
        gamma_data_f%v_c_dewald = v_c_dewald%view%ptr()
        gamma_data_f%v_c_atom_ewald = v_c_atom_ewald%view%ptr()
        gamma_data_f%v_c_g123 = v_c_g123%view%ptr()
        gamma_data_f%v_c_Qall = v_c_Qall%view%ptr()
        gamma_data_f%natoms = natoms
        gamma_data_f%icluster = icluster
        gamma_data_f%v_c_a123vec = v_c_a123vec%view%ptr()
        gamma_data_f%c_il123mx = c_il123mx

        call fgamma(gamma_data_f)

end subroutine gamma1

subroutine print_ratom( ratom )
        use :: precision_kinds
        implicit none
        type(gamma_view2d), intent(in) :: ratom

        call f_print_ratom(ratom%view%ptr())

end subroutine print_ratom

end module kokkos_interface
