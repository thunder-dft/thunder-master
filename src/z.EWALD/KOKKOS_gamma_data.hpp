// copyright info:
//
//                             @Copyright 2024
//                           Fireball Committee
// Hong Kong Quantum AI Laboratory, Ltd. - James P. Lewis, Chair
// Universidad de Madrid - Jose Ortega
// Academy of Sciences of the Czech Republic - Pavel Jelinek
// Arizona State University - Otto F. Sankey
//
// Previous and/or current contributors:
// Auburn University - Jian Jun Dong
// California Institute of Technology - Brandon Keith
// Czech Institute of Physics - Prokop Hapala
// Czech Institute of Physics - Vladimír Zobač
// Dublin Institute of Technology - Barry Haycock
// Pacific Northwest National Laboratory - Kurt Glaesemann
// University of Texas at Austin - Alex Demkov
// Ohio University - Dave Drabold
// Synfuels China Technology Co., Ltd. - Pengju Ren
// Washington University - Pete Fedders
// West Virginia University - Ning Ma and Hao Wang
// Computer Network Information Center, Chinese Academy of Sciences    &
//     & University of Chinese Academy of Sciences - Runfeng Jin
// also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
//      and Spencer Shellman
//
// RESTRICTED RIGHTS LEGEND
// Use, duplication, or disclosure of this software and its documentation
// by the Government is subject to restrictions as set forth in subdivision
// { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
// clause at 52.227-7013.
// ===========================================================================
// Code written by:
// James P. Lewis
// Unit 909 of Building 17W
// 17 Science Park West Avenue
// Pak Shek Kok, New Territories 999077
// Hong Kong
//
// Phone: +852 6612 9539 (mobile)
//
// Runfeng Jin
// Computer Network Information Center, Chinese Academy of Sciences
// Beijing, China
//
// E-mail: jsfaraway@gmail.com
// ===========================================================================
//
// M_Dassemble_ewald_DOGS_kokkos_interface
// ===========================================================================
// Program Description
// ===========================================================================
//>       This file contains the c Structure to receive the data from Fortran.
//
//
// ===========================================================================
#ifndef KOKKOS_GAMMA_DATA_HPP
#define KOKKOS_GAMMA_DATA_HPP

#include "Kokkos_Core.hpp"
#include "flcl-cxx.hpp"
namespace kokkos_gamma_data
{
/* manage mix precision */
// #define SINGLE_PRECISION
#ifdef SINGLE_PRECISION
 using gamma_view2d = flcl::dualview_r32_2d_t;
 using gamma_view3d = flcl::dualview_r32_3d_t;
 using gamma_real = float;
#else
 using gamma_view2d = flcl::dualview_r64_2d_t;
 using gamma_view3d = flcl::dualview_r64_3d_t;
 using gamma_real = double;
#endif

 struct gamma_data
 {
  int ig1mx;
  int ig2mx;
  int ig3mx;
  gamma_real kappa;
  gamma_view2d *v_c_ratom;
  gamma_view2d *v_c_ewald;
  gamma_view3d *v_c_dewald;
  gamma_view2d *v_c_atom_ewald;
  gamma_view2d *v_c_g123;
  gamma_view2d *v_c_Qall;
  gamma_real volume;
  int natoms;
  int icluster;
  gamma_view2d *v_c_a123vec;
  int il123mx[3];
 };
 // extern gamma_data gamma_data_f;

}

#endif // KOKKOS_GAMMA_DATA_HPP