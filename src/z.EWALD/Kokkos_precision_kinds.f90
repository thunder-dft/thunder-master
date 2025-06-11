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
! precision_kinds
! ===========================================================================
! Program Description
! ===========================================================================
!>       This is a module used for manage the precision used in GPU for 
!!  Ewald computation.
!!
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
! ===========================================================================
! Module Declaration
! ============================================================================

! #define SINGLE_PRECISION
! #ifndef SINGLE_PRECISION
 ! use iso_fortran_env, only: rk => REAL64
 ! use  iso_c_binding, only: crk => C_DOUBLE, cik => C_INT, c_ptr
 ! use flcl_mod : gamma_view2d => dualview_r64_2d_t, gamma_view3d => dualview_r64_3d_t
! #else
! use iso_fortran_env, only: rk => REAL32
! use  iso_c_binding, only: crk => C_FLOAT, cik => C_INT, c_ptr
! use flcl_mod, only : gamma_view2d => dualview_r32_2d_t, gamma_view3d => dualview_r32_3d_t
! ! #endif

!single precision
! module precision_kinds
!         use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64, &
!                 input_unit, output_unit, error_unit
!         use iso_c_binding, only: c_int, c_float, c_double, c_ptr
!         use :: flcl_mod
!         implicit none
!         integer, parameter :: rk = real32
!         integer, parameter :: crk = c_float
!         integer, parameter :: cik = c_int
!         type gamma_view2d
!                 type(dualview_r32_2d_t) :: view
!         end type gamma_view2d
!         type gamma_view3d
!                 type(dualview_r32_3d_t) :: view
!         end type gamma_view3d
! end module precision_kinds

! double precision
module precision_kinds
        use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64, &
                input_unit, output_unit, error_unit
        use iso_c_binding, only: c_int, c_float, c_double, c_ptr
        use :: flcl_mod
        implicit none
        integer, parameter :: rk = real64
        integer, parameter :: crk = c_double
        integer, parameter :: cik = c_int
        type gamma_view2d
                type(dualview_r64_2d_t) :: view
        end type gamma_view2d
        type gamma_view3d
                type(dualview_r64_3d_t) :: view
        end type gamma_view3d
end module precision_kinds
