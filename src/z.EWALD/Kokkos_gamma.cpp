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
//>       This file impelment the initial kokkos version of ewald and output function.
//
// It contains the following function:
//          gamma1_c : Initial Kokkos version of Ewald.
//          gamma_res: Print out the key result of Ewald.
//
// ===========================================================================
#include "Kokkos_gamma.hpp"

namespace kokkos_gamma
{

 void gamma1_c(int ig1mx, int ig2mx, int ig3mx, float kappa, gamma_view2d v_c_ratom, gamma_view2d v_c_ewald, gamma_view3d v_c_dewald, gamma_view2d v_c_atom_ewald, gamma_view2d v_c_g123, gamma_view2d v_c_Qall, float volume, int natoms)
 {
  const double PI = 3.141592653589793238462643;

  for (int ig1 = -ig1mx; ig1 <= ig1mx; ++ig1)
   for (int ig2 = -ig2mx; ig2 <= ig2mx; ++ig2)
    for (int ig3 = -ig3mx; ig3 <= ig3mx; ++ig3)
    {
     // Skip the zero vector
     if (!(ig1 == 0 && ig2 == 0 && ig3 == 0))
     {
      float g[3];
      g[0] = ig1 * v_c_g123.h_view(0, 0) + ig2 * v_c_g123.h_view(0, 1) + ig3 * v_c_g123.h_view(0, 2);
      g[1] = ig1 * v_c_g123.h_view(1, 0) + ig2 * v_c_g123.h_view(1, 1) + ig3 * v_c_g123.h_view(1, 2);
      g[2] = ig1 * v_c_g123.h_view(2, 0) + ig2 * v_c_g123.h_view(2, 1) + ig3 * v_c_g123.h_view(2, 2);

      float magnitude_g = Kokkos::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
      float argument = magnitude_g * magnitude_g / (4.0 * kappa * kappa);
      double stuff = 4.0 * PI * Kokkos::exp(-argument) / (magnitude_g * magnitude_g * volume);

      for (int iatom = 0; iatom < natoms; ++iatom)
      {
       for (int jatom = iatom; jatom < natoms; ++jatom)
       {
        double gdotb = g[0] * (v_c_ratom.h_view(iatom, 0) - v_c_ratom.h_view(jatom, 0)) +
                       g[1] * (v_c_ratom.h_view(iatom, 1) - v_c_ratom.h_view(jatom, 1)) +
                       g[2] * (v_c_ratom.h_view(iatom, 2) - v_c_ratom.h_view(jatom, 2));

        double factor = (jatom == iatom) ? 0.5 * stuff : stuff;
        double qq = (v_c_Qall.h_view(iatom, 1) * v_c_Qall.h_view(jatom, 1) - v_c_Qall.h_view(iatom, 0) * v_c_Qall.h_view(jatom, 0));

        // Atomic updates need atomic operations to prevent race conditions
        Kokkos::atomic_add(&v_c_ewald.h_view(iatom, jatom), factor * Kokkos::cos(gdotb));
        Kokkos::atomic_add(&v_c_ewald.h_view(jatom, iatom), factor * Kokkos::cos(gdotb));

        for (int dim = 0; dim < 3; ++dim)
        {
         Kokkos::atomic_add(&v_c_atom_ewald.h_view(iatom, dim), qq * 2.0 * stuff * Kokkos::sin(gdotb) * g[dim]);
         Kokkos::atomic_add(&v_c_atom_ewald.h_view(jatom, dim), -qq * 2.0 * stuff * Kokkos::sin(gdotb) * g[dim]);
         Kokkos::atomic_add(&v_c_dewald.h_view(dim, iatom, jatom), -factor * Kokkos::sin(gdotb) * g[dim]);
         Kokkos::atomic_add(&v_c_dewald.h_view(dim, jatom, iatom), factor * Kokkos::sin(gdotb) * g[dim]);
        }
       }
      }
     }
    }

 } // gamma1_c

 void gamma_res(gamma_view2d v_c_ewald, gamma_view2d v_c_atom_ewald, gamma_view3d v_c_dewald, int natoms)
 {
  // Create three files, print out v_c_ewald, v_c_dewald, v_c_atom_ewald into different file
  FILE *fp_ewald = fopen("v_c_ewald.txt", "w");
  FILE *fp_dewald = fopen("v_c_dewald.txt", "w");
  FILE *fp_atom_ewald = fopen("v_c_atom_ewald.txt", "w");
  for (int i = 0; i < natoms; i++)
  {
   fprintf(fp_atom_ewald, "%15.13f %15.13f %15.13f\n", v_c_atom_ewald.h_view(i, 0), v_c_atom_ewald.h_view(i, 1), v_c_atom_ewald.h_view(i, 2));
   for (int j = 0; j < natoms; j++)
   {
    fprintf(fp_ewald, "%15.13f\n", v_c_ewald.h_view(i, j));
    fprintf(fp_dewald, "%15.13f %15.13f %15.13f\n", v_c_dewald.h_view(0, i, j), v_c_dewald.h_view(1, i, j), v_c_dewald.h_view(2, i, j));
   }
  }
 }

} // namespace kokkos_gamma