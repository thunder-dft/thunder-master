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
//>       This file contains the current kokkos version of ewald and output function.
//      This version seperate the diagonal and off-diagonal element.
//
// It contains the following function:
//          FunctorGamma1_diag :  Diagonal part of Gamma1 and Gamma3
//          FunctorGamma1_offdiag : Offdiagonal part of Gamma1
//          FunctorGamma2_diag : Diagonal part of Gamma2
//          FunctorGamma2_offdiag : Offdiagonal part of Gamma2
//
// ===========================================================================
#include "Kokkos_Core.hpp"
#include "flcl-cxx.hpp"
#include <cmath>

#include "Kokkos_gamma_data.hpp"

namespace kokkos_gamma
{
  using kokkos_gamma_data::gamma_view2d;
  using kokkos_gamma_data::gamma_view3d;
  using kokkos_gamma_data::gamma_data;
  using kokkos_gamma_data::gamma_real;

  class FunctorGamma2_diag
  {
  public:
    FunctorGamma2_diag(gamma_data gamma_data_f)
        : ig1mx(gamma_data_f.ig1mx), ig2mx(gamma_data_f.ig2mx), ig3mx(gamma_data_f.ig3mx),
          kappa(gamma_data_f.kappa),
          v_c_ratom(*(gamma_data_f.v_c_ratom)),
          v_c_ewald(*(gamma_data_f.v_c_ewald)),
          v_c_dewald(*(gamma_data_f.v_c_dewald)),
          v_c_atom_ewald(*(gamma_data_f.v_c_atom_ewald)),
          v_c_g123(*(gamma_data_f.v_c_g123)),
          v_c_Qall(*(gamma_data_f.v_c_Qall)),
          volume(gamma_data_f.volume),
          natoms(gamma_data_f.natoms),
          icluster(gamma_data_f.icluster),
          v_c_a123vec(*(gamma_data_f.v_c_a123vec)),
          il1mx(gamma_data_f.il123mx[0]),
          il2mx(gamma_data_f.il123mx[1]),
          il3mx(gamma_data_f.il123mx[2])
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(const int &iorder) const
    {
      // int iatom = iorder;
      // int jatom = iorder;
      // v_c_ewald.d_view(iorder, iorder) = 0;
      double tmp = 0.0;
      for (int il1 = -il1mx; il1 <= il1mx; ++il1)
        for (int il2 = -il2mx; il2 <= il2mx; ++il2)
          for (int il3 = -il3mx; il3 <= il3mx; ++il3)
          {
            double cvec[3];

            cvec[0] = il1 * v_c_a123vec.d_view(0, 0) + il2 * v_c_a123vec.d_view(1, 0) + il3 * v_c_a123vec.d_view(2, 0);
            cvec[1] = il1 * v_c_a123vec.d_view(0, 1) + il2 * v_c_a123vec.d_view(1, 1) + il3 * v_c_a123vec.d_view(2, 1);
            cvec[2] = il1 * v_c_a123vec.d_view(0, 2) + il2 * v_c_a123vec.d_view(1, 2) + il3 * v_c_a123vec.d_view(2, 2);

            double z = Kokkos::sqrt(cvec[0] * cvec[0] + cvec[1] * cvec[1] + cvec[2] * cvec[2]);
            if (z > 0.0001)
            {
              tmp += Kokkos::erfc(kappa * z) / z;
            }
          }
      v_c_ewald.d_view(iorder, iorder) += tmp;
    } // gamma2 diagnal elements
    int ig1mx, il1mx;
    int ig2mx, ig3mx, il2mx, il3mx;

  private:
    gamma_real kappa;
    const gamma_real volume;
    gamma_view2d v_c_ratom, v_c_ewald;
    gamma_view2d v_c_atom_ewald; //(natoms, xyz)
    gamma_view2d v_c_g123, v_c_Qall;
    gamma_view3d v_c_dewald; //(3, natoms, natoms)
    const double PI = 3.141592653589793238462643;
    const int natoms;
    int icluster;
    gamma_view2d v_c_a123vec;
  };
  /**********Gamma 2 offdiag*************/
  class FunctorGamma2_offdiag
  {
  public:
    FunctorGamma2_offdiag(gamma_data gamma_data_f)
        : ig1mx(gamma_data_f.ig1mx), ig2mx(gamma_data_f.ig2mx), ig3mx(gamma_data_f.ig3mx),
          kappa(gamma_data_f.kappa),
          v_c_ratom(*(gamma_data_f.v_c_ratom)),
          v_c_ewald(*(gamma_data_f.v_c_ewald)),
          v_c_dewald(*(gamma_data_f.v_c_dewald)),
          v_c_atom_ewald(*(gamma_data_f.v_c_atom_ewald)),
          v_c_g123(*(gamma_data_f.v_c_g123)),
          v_c_Qall(*(gamma_data_f.v_c_Qall)),
          volume(gamma_data_f.volume),
          natoms(gamma_data_f.natoms),
          icluster(gamma_data_f.icluster),
          v_c_a123vec(*(gamma_data_f.v_c_a123vec)),
          il1mx(gamma_data_f.il123mx[0]),
          il2mx(gamma_data_f.il123mx[1]),
          il3mx(gamma_data_f.il123mx[2])
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int &iorder) const
    {
      // only iterate off-diag elements
      int iatom = natoms - 2 - std::floor(Kokkos::sqrt(-8 * iorder + 4 * natoms * (natoms - 1) - 7) / 2.0 - 0.5);
      int jatom = iorder + iatom + 1 - natoms * (natoms - 1) / 2 + (natoms - iatom) * ((natoms - iatom) - 1) / 2;
      // initialize memory
      // v_c_atom_ewald.d_view(iatom, 0) = v_c_atom_ewald.d_view(iatom, 1) = v_c_atom_ewald.d_view(iatom, 2) = 0.0;
      // v_c_atom_ewald.d_view(jatom, 0) = v_c_atom_ewald.d_view(jatom, 1) = v_c_atom_ewald.d_view(jatom, 2) = 0.0;
      // v_c_ewald.d_view(iatom, jatom) = v_c_ewald.d_view(jatom, iatom) = 0.0;
      // v_c_dewald.d_view(0, jatom, iatom) = v_c_dewald.d_view(1, jatom, iatom) = v_c_dewald.d_view(2, jatom, iatom) = 0.0;
      // v_c_dewald.d_view(0, iatom, jatom) = v_c_dewald.d_view(1, iatom, jatom) = v_c_dewald.d_view(2, iatom, jatom) = 0.0;

      double diff1 = v_c_ratom.d_view(iatom, 0) - v_c_ratom.d_view(jatom, 0);
      double diff2 = v_c_ratom.d_view(iatom, 1) - v_c_ratom.d_view(jatom, 1);
      double diff3 = v_c_ratom.d_view(iatom, 2) - v_c_ratom.d_view(jatom, 2);
      double tmp_ewald = 0.0;
      double tmp_dewald1 = 0.0, tmp_dewald2 = 0.0, tmp_dewald3 = 0.0;
      double tmp_atom_ewald1 = 0.0, tmp_atom_ewald2 = 0.0, tmp_atom_ewald3 = 0.0;
      for (int il1 = -il1mx; il1 <= il1mx; ++il1)
        for (int il2 = -il2mx; il2 <= il2mx; ++il2)
          for (int il3 = -il3mx; il3 <= il3mx; ++il3)
          {
            // double factor =  1.0;
            double factorf = 2.0;
            double cvec[3];

            cvec[0] = il1 * v_c_a123vec.d_view(0, 0) + il2 * v_c_a123vec.d_view(1, 0) + il3 * v_c_a123vec.d_view(2, 0) + diff1;
            cvec[1] = il1 * v_c_a123vec.d_view(0, 1) + il2 * v_c_a123vec.d_view(1, 1) + il3 * v_c_a123vec.d_view(2, 1) + diff2;
            cvec[2] = il1 * v_c_a123vec.d_view(0, 2) + il2 * v_c_a123vec.d_view(1, 2) + il3 * v_c_a123vec.d_view(2, 2) + diff3;
            double z = Kokkos::sqrt(cvec[0] * cvec[0] + cvec[1] * cvec[1] + cvec[2] * cvec[2]);
            if (z > 0.0001)
            {
              double argument = kappa * z;
              // double erfc_arg = Kokkos::erfc(argument) / z;
              double force_factor = (2.0 * Kokkos::exp(-argument * argument) * kappa * z / Kokkos::sqrt(PI) + Kokkos::erfc(argument));
              // Kokkos::atomic_add(&v_c_ewald.d_view(iatom, jatom), erfc_arg);
              // Kokkos::atomic_add(&v_c_ewald.d_view(jatom, iatom), erfc_arg);
              // v_c_ewald.d_view(iatom, jatom) += erfc_arg;
              // v_c_ewald.d_view(jatom, iatom) += erfc_arg;
              tmp_ewald += Kokkos::erfc(argument) / z;

              double qq = v_c_Qall.d_view(iatom, 1) * v_c_Qall.d_view(jatom, 1) - v_c_Qall.d_view(iatom, 0) * v_c_Qall.d_view(jatom, 0);
              tmp_dewald1 += cvec[0] * force_factor / (z * z * z);
              tmp_dewald2 += cvec[1] * force_factor / (z * z * z);
              tmp_dewald3 += cvec[2] * force_factor / (z * z * z);
              tmp_atom_ewald1 += qq * factorf * cvec[0] * force_factor / (z * z * z);
              tmp_atom_ewald2 += qq * factorf * cvec[1] * force_factor / (z * z * z);
              tmp_atom_ewald3 += qq * factorf * cvec[2] * force_factor / (z * z * z);

              // for (int d = 0; d < 3; ++d)
              // {
              //   double atomewald_contribution = qq * factorf * cvec[d] * force_factor / (z * z * z);
              //   double dewald_contribution = cvec[d] * force_factor / (z * z * z);
              //   Kokkos::atomic_add(&v_c_atom_ewald.d_view(iatom, d), atomewald_contribution);
              //   Kokkos::atomic_add(&v_c_atom_ewald.d_view(jatom, d), -atomewald_contribution);
              //   // Kokkos::atomic_add(&v_c_dewald.d_view(d, iatom, jatom), -dewald_contribution);
              //   // Kokkos::atomic_add(&v_c_dewald.d_view(d, jatom, iatom), dewald_contribution);
              //   v_c_dewald.d_view(d, iatom, jatom) -= dewald_contribution;
              //   v_c_dewald.d_view(d, jatom, iatom) += dewald_contribution;
              // }
            }
          }

      v_c_ewald.d_view(iatom, jatom) += tmp_ewald;
      v_c_ewald.d_view(jatom, iatom) += tmp_ewald;
      v_c_dewald.d_view(0, iatom, jatom) -= tmp_dewald1;
      v_c_dewald.d_view(1, iatom, jatom) -= tmp_dewald2;
      v_c_dewald.d_view(2, iatom, jatom) -= tmp_dewald3;
      v_c_dewald.d_view(0, jatom, iatom) += tmp_dewald1;
      v_c_dewald.d_view(1, jatom, iatom) += tmp_dewald2;
      v_c_dewald.d_view(2, jatom, iatom) += tmp_dewald3;
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(iatom, 0), tmp_atom_ewald1);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(iatom, 1), tmp_atom_ewald2);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(iatom, 2), tmp_atom_ewald3);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(jatom, 0), -tmp_atom_ewald1);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(jatom, 1), -tmp_atom_ewald2);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(jatom, 2), -tmp_atom_ewald3);
    }

    int ig1mx, il1mx;
    int ig2mx, ig3mx, il2mx, il3mx;

  private:
    gamma_real kappa;
    const gamma_real volume;
    gamma_view2d v_c_ratom, v_c_ewald;
    gamma_view2d v_c_atom_ewald; //(natoms, xyz)
    gamma_view2d v_c_g123, v_c_Qall;
    gamma_view3d v_c_dewald; //(3, natoms, natoms)
    const double PI = 3.141592653589793238462643;
    const int natoms;
    int icluster;
    gamma_view2d v_c_a123vec;
  }; // functor gamma2 off-diag

  /**********Gamma 1 diag*************/
  class FunctorGamma1_diag
  {
  public:
    FunctorGamma1_diag(gamma_data gamma_data_f)
        : ig1mx(gamma_data_f.ig1mx), ig2mx(gamma_data_f.ig2mx), ig3mx(gamma_data_f.ig3mx),
          kappa(gamma_data_f.kappa),
          v_c_ratom(*(gamma_data_f.v_c_ratom)),
          v_c_ewald(*(gamma_data_f.v_c_ewald)),
          v_c_dewald(*(gamma_data_f.v_c_dewald)),
          v_c_atom_ewald(*(gamma_data_f.v_c_atom_ewald)),
          v_c_g123(*(gamma_data_f.v_c_g123)),
          v_c_Qall(*(gamma_data_f.v_c_Qall)),
          volume(gamma_data_f.volume),
          natoms(gamma_data_f.natoms),
          icluster(gamma_data_f.icluster),
          v_c_a123vec(*(gamma_data_f.v_c_a123vec)),
          il1mx(gamma_data_f.il123mx[0]),
          il2mx(gamma_data_f.il123mx[1]),
          il3mx(gamma_data_f.il123mx[2])
    {
      // print out il123mx and ig123mx
      // printf("c il123mx[0] is %d\n", il1mx);
      // printf("c il123mx[1] is %d\n", il2mx);
      // printf("c il123mx[2] is %d\n", il3mx);
      // printf("c ig1mx, ig2mx, ig3mx is %d, %d, %d\n", ig1mx, ig2mx, ig3mx);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int &iorder) const
    {
      double tmp = 0.0;
      for (int ig1 = -ig1mx; ig1 <= ig1mx; ++ig1)
        for (int ig2 = -ig2mx; ig2 <= ig2mx; ++ig2)
          for (int ig3 = -ig3mx; ig3 <= ig3mx; ++ig3)
          {
            // Skip the zero vector
            if (!(ig1 == 0 && ig2 == 0 && ig3 == 0))
            {
              double g[3];
              g[0] = ig1 * v_c_g123.d_view(0, 0) + ig2 * v_c_g123.d_view(0, 1) + ig3 * v_c_g123.d_view(0, 2);
              g[1] = ig1 * v_c_g123.d_view(1, 0) + ig2 * v_c_g123.d_view(1, 1) + ig3 * v_c_g123.d_view(1, 2);
              g[2] = ig1 * v_c_g123.d_view(2, 0) + ig2 * v_c_g123.d_view(2, 1) + ig3 * v_c_g123.d_view(2, 2);

              double magnitude_g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
              // double magnitude_g = Kokkos::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
              double argument = magnitude_g2 / (4.0 * kappa * kappa);
              double stuff = 4.0 * PI * Kokkos::exp(-argument) / (magnitude_g2 * volume);

              // Atomic updates need atomic operations to prevent race conditions
              // v_c_ewald.d_view(iorder, iorder) += stuff;
              tmp += stuff;
              // Kokkos::atomic_add(&v_c_ewald.d_view(iorder, iorder), stuff);
              // Kokkos::atomic_add(&v_c_ewald.d_view(jatom, iatom), factor );
            }
          }
      /*  For gamma 3*/
      v_c_ewald.d_view(iorder, iorder) += tmp - 2.0 * kappa / PI_sqrt;
      // v_c_ewald.d_view(iorder, iorder) -= 2.0 * kappa / Kokkos::sqrt(PI);
    }

    int ig1mx, il1mx;
    int ig2mx, ig3mx, il2mx, il3mx;

  private:
    gamma_real kappa;
    const gamma_real volume;
    gamma_view2d v_c_ratom, v_c_ewald;
    gamma_view2d v_c_atom_ewald; //(natoms, xyz)
    gamma_view2d v_c_g123, v_c_Qall;
    gamma_view3d v_c_dewald; //(3, natoms, natoms)
    const double PI = 3.141592653589793238462643;
    const double PI_sqrt = 1.77245385090551588;
    const int natoms;
    int icluster;
    gamma_view2d v_c_a123vec;
  }; // kokkos gamma1_diag

  /**********Gamma 1 offdiag*************/
  class FunctorGamma1_offdiag
  {
  public:
    FunctorGamma1_offdiag(gamma_data gamma_data_f)
        : ig1mx(gamma_data_f.ig1mx), ig2mx(gamma_data_f.ig2mx), ig3mx(gamma_data_f.ig3mx),
          kappa(gamma_data_f.kappa),
          v_c_ratom(*(gamma_data_f.v_c_ratom)),
          v_c_ewald(*(gamma_data_f.v_c_ewald)),
          v_c_dewald(*(gamma_data_f.v_c_dewald)),
          v_c_atom_ewald(*(gamma_data_f.v_c_atom_ewald)),
          v_c_g123(*(gamma_data_f.v_c_g123)),
          v_c_Qall(*(gamma_data_f.v_c_Qall)),
          volume(gamma_data_f.volume),
          natoms(gamma_data_f.natoms),
          icluster(gamma_data_f.icluster),
          v_c_a123vec(*(gamma_data_f.v_c_a123vec)),
          il1mx(gamma_data_f.il123mx[0]),
          il2mx(gamma_data_f.il123mx[1]),
          il3mx(gamma_data_f.il123mx[2])
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int &iorder) const
    {
      // only iterate off-diag elements
      int iatom = natoms - 2 - std::floor(Kokkos::sqrt(-8 * iorder + 4 * natoms * (natoms - 1) - 7) / 2.0 - 0.5);
      int jatom = iorder + iatom + 1 - natoms * (natoms - 1) / 2 + (natoms - iatom) * ((natoms - iatom) - 1) / 2;
      double tmp_ewald = 0.0;
      double tmp_dewald1 = 0.0, tmp_dewald2 = 0.0, tmp_dewald3 = 0.0;
      double tmp_atom_ewald1 = 0.0, tmp_atom_ewald2 = 0.0, tmp_atom_ewald3 = 0.0;

      for (int ig1 = -ig1mx; ig1 <= ig1mx; ++ig1)
        for (int ig2 = -ig2mx; ig2 <= ig2mx; ++ig2)
          for (int ig3 = -ig3mx; ig3 <= ig3mx; ++ig3)
          {
            // Skip the zero vector
            if (!(ig1 == 0 && ig2 == 0 && ig3 == 0))
            {
              double g[3];
              g[0] = ig1 * v_c_g123.d_view(0, 0) + ig2 * v_c_g123.d_view(0, 1) + ig3 * v_c_g123.d_view(0, 2);
              g[1] = ig1 * v_c_g123.d_view(1, 0) + ig2 * v_c_g123.d_view(1, 1) + ig3 * v_c_g123.d_view(1, 2);
              g[2] = ig1 * v_c_g123.d_view(2, 0) + ig2 * v_c_g123.d_view(2, 1) + ig3 * v_c_g123.d_view(2, 2);

              double magnitude_g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
              double argument = magnitude_g2 / (4.0 * kappa * kappa); // 这里可以优化，magnitude_g不求平方
              double stuff = 4.0 * PI * Kokkos::exp(-argument) / (magnitude_g2* volume);
              /* Above is same for the same igxxx */
              // for (int iatom = 0; iatom < natoms; ++iatom)
              // for (int jatom = iatom; jatom < natoms; ++jatom)
              double gdotb = g[0] * (v_c_ratom.d_view(iatom, 0) - v_c_ratom.d_view(jatom, 0)) +
                             g[1] * (v_c_ratom.d_view(iatom, 1) - v_c_ratom.d_view(jatom, 1)) +
                             g[2] * (v_c_ratom.d_view(iatom, 2) - v_c_ratom.d_view(jatom, 2));

              double factor = stuff * Kokkos::sin(gdotb);
              double qq = 2.0 * factor * (v_c_Qall.d_view(iatom, 1) * v_c_Qall.d_view(jatom, 1) - v_c_Qall.d_view(iatom, 0) * v_c_Qall.d_view(jatom, 0));

              // Atomic updates need atomic operations to prevent race conditions
              // v_c_ewald.d_view(iatom, jatom) += factor * Kokkos::cos(gdotb);
              // v_c_ewald.d_view(jatom, iatom) += factor * Kokkos::cos(gdotb);
              tmp_ewald += stuff * Kokkos::cos(gdotb);
              tmp_dewald1 += factor * g[0];
              tmp_dewald2 += factor * g[1];
              tmp_dewald3 += factor * g[2];
              tmp_atom_ewald1 += qq * g[0];
              tmp_atom_ewald2 += qq * g[1];
              tmp_atom_ewald3 += qq * g[2];

              // for (int dim = 0; dim < 3; ++dim)
              // {
              //   Kokkos::atomic_add(&v_c_atom_ewald.d_view(iatom, dim), qq * 2.0 * stuff * Kokkos::sin(gdotb) * g[dim]);
              //   Kokkos::atomic_add(&v_c_atom_ewald.d_view(jatom, dim), -qq * 2.0 * stuff * Kokkos::sin(gdotb) * g[dim]);
              //   v_c_dewald.d_view(dim, iatom, jatom) -= factor * Kokkos::sin(gdotb) * g[dim];
              //   v_c_dewald.d_view(dim, jatom, iatom) += factor * Kokkos::sin(gdotb) * g[dim];
              //   // Kokkos::atomic_add(&v_c_dewald.d_view(dim, iatom, jatom), -factor * Kokkos::sin(gdotb) * g[dim]);
              //   // Kokkos::atomic_add(&v_c_dewald.d_view(dim, jatom, iatom), factor * Kokkos::sin(gdotb) * g[dim]);
              // }
            }
          }
      v_c_ewald.d_view(iatom, jatom) += tmp_ewald;
      v_c_ewald.d_view(jatom, iatom) += tmp_ewald;
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(iatom, 0), tmp_atom_ewald1);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(iatom, 1), tmp_atom_ewald2);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(iatom, 2), tmp_atom_ewald3);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(jatom, 0), -tmp_atom_ewald1);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(jatom, 1), -tmp_atom_ewald2);
      Kokkos::atomic_add(&v_c_atom_ewald.d_view(jatom, 2), -tmp_atom_ewald3);
      v_c_dewald.d_view(0, iatom, jatom) -= tmp_dewald1;
      v_c_dewald.d_view(1, iatom, jatom) -= tmp_dewald2;
      v_c_dewald.d_view(2, iatom, jatom) -= tmp_dewald3;
      v_c_dewald.d_view(0, jatom, iatom) += tmp_dewald1;
      v_c_dewald.d_view(1, jatom, iatom) += tmp_dewald2;
      v_c_dewald.d_view(2, jatom, iatom) += tmp_dewald3;

    } // opertaor

    int ig1mx, il1mx;
    int ig2mx, ig3mx, il2mx, il3mx;

  private:
    gamma_real kappa;
    const gamma_real volume;
    gamma_view2d v_c_ratom, v_c_ewald;
    gamma_view2d v_c_atom_ewald; //(natoms, xyz)
    gamma_view2d v_c_g123, v_c_Qall;
    gamma_view3d v_c_dewald; //(3, natoms, natoms)
    const double PI = 3.141592653589793238462643;
    const int natoms;
    int icluster;
    gamma_view2d v_c_a123vec;
  }; // kokkos gamma1_offdiag

} // CLASS kokkos_gamma