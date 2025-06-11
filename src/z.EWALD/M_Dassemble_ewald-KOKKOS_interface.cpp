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
//>       This file contains the interface with the FORTRAN.
//
// It contains the following modules:
//
//       c_print_ratom - Print out the data
//       gamma1_cluster - Set the ig*mx value 
//       gamma2_cluster - Set the il*mx value
//       c_gamma - Receive the data from fortran and call GPU to compute.
//
// ===========================================================================
#include <Kokkos_Core.hpp>
#include "flcl-cxx.hpp"
#include <chrono>

#include "KOKKOS_gamma.hpp"
#include "KOKKOS_gamma_s.hpp"

// using gamma_view2d = flcl::gamma_view2d;

extern "C"
{

  void c_print_ratom(kokkos_gamma_data::gamma_view2d **v_ratom)
  {
    using kokkos_gamma_data::gamma_view2d;
    using kokkos_gamma_data::gamma_view3d;
    using flcl::view_from_ndarray;

    gamma_view2d ratom = **v_ratom;
    int natoms = ratom.h_view.extent(0);
    printf("dimension 0 is %d\n", natoms);
    // print all ratom(:,3) into a file c_ratom.txt, with each line has 3 columns
    FILE *fp = fopen("c_ratom.txt", "w");
    for (int i = 0; i < natoms; i++)
    {
      fprintf(fp, "%10.6f %10.6f %10.6f\n", ratom.h_view(i, 0), ratom.h_view(i, 1), ratom.h_view(i, 2));
    }

    // ratom.template modify<typename view_type::host_mirror_space>();
    // y.template sync<typename view_type::execution_space>();

    // Kokkos::parallel_for( "axpy", y.extent(0), KOKKOS_LAMBDA( const size_t idx)
    // {
    //   y.d_view(idx) += d_alpha * x.d_view(idx);
    // });

    // y.template modify<typename view_type::execution_space>();
    // y.template sync<typename view_type::host_mirror_space>();

    return;
  }

  void gamma1_cluster(kokkos_gamma_data::gamma_data *gamma_data_f)
  {
    if (gamma_data_f->icluster == 1)
    {
      gamma_data_f->ig1mx = 0;
      gamma_data_f->ig2mx = 0;
      gamma_data_f->ig3mx = 0;
    }
  }

  void gamma2_cluster(kokkos_gamma_data::gamma_data *gamma_data_f)
  {
    if (gamma_data_f->icluster == 1)
    {
      gamma_data_f->il123mx[0] = 0;
      gamma_data_f->il123mx[1] = 0;
      gamma_data_f->il123mx[2] = 0;
      gamma_data_f->kappa = 0.0;
    }
  }

  void c_gamma(kokkos_gamma_data::gamma_data *gamma_data_f)
  {
    using kokkos_gamma_data::gamma_view2d;
    using kokkos_gamma_data::gamma_view3d;
    using Kokkos::MDRangePolicy;
    using Kokkos::parallel_for;
    using Kokkos::RangePolicy;
    using kokkos_gamma_data::gamma_data;

    int natoms = gamma_data_f->natoms;

    // transfer data HOST -> DEVICE
    auto start = std::chrono::high_resolution_clock::now();
    gamma_data_f->v_c_ratom->template modify<typename gamma_view2d::host_mirror_space>();
    gamma_data_f->v_c_ewald->template modify<typename gamma_view2d::host_mirror_space>();
    gamma_data_f->v_c_dewald->template modify<typename gamma_view3d::host_mirror_space>();
    gamma_data_f->v_c_atom_ewald->template modify<typename gamma_view2d::host_mirror_space>();
    gamma_data_f->v_c_g123->template modify<typename gamma_view2d::host_mirror_space>();
    gamma_data_f->v_c_Qall->template modify<typename gamma_view2d::host_mirror_space>();
    gamma_data_f->v_c_a123vec->template modify<typename gamma_view2d::host_mirror_space>();

    gamma_data_f->v_c_ratom->template sync<typename gamma_view2d::execution_space>();
    gamma_data_f->v_c_ewald->template sync<typename gamma_view2d::execution_space>();
    gamma_data_f->v_c_dewald->template sync<typename gamma_view3d::execution_space>();
    gamma_data_f->v_c_atom_ewald->template sync<typename gamma_view2d::execution_space>();
    gamma_data_f->v_c_g123->template sync<typename gamma_view2d::execution_space>();
    gamma_data_f->v_c_Qall->template sync<typename gamma_view2d::execution_space>();
    gamma_data_f->v_c_a123vec->template sync<typename gamma_view2d::execution_space>();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    // std::cout << "KOKKOS Time taken by host->device data: " << diff.count() << " seconds" << std::endl;

    // test the time of gamma constructor
    start = std::chrono::high_resolution_clock::now();
    // kokkos_gamma::FunctorGamma gamma(*gamma_data_f);

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    // std::cout << "KOKKOS Time taken by  construct functor : " << diff.count() << " seconds" << std::endl;

    /* gamma1 */
    gamma1_cluster(gamma_data_f);
    kokkos_gamma::FunctorGamma1_diag gamma1_diag(*gamma_data_f);
    kokkos_gamma::FunctorGamma1_offdiag gamma1_offdiag(*gamma_data_f);

    start = std::chrono::high_resolution_clock::now();
    parallel_for("gamma", RangePolicy<Kokkos::Cuda, Kokkos::IndexType<int>, Kokkos::LaunchBounds<256, 4>>(0, natoms), gamma1_diag);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    // std::cout << "KOKKOS Time taken by gamma 1 diag: " << diff.count() << " seconds" << std::endl;

    int num_upper = natoms * (natoms + 1) / 2 - natoms;
    start = std::chrono::high_resolution_clock::now();
    parallel_for("gamma", RangePolicy<Kokkos::Cuda, Kokkos::IndexType<int>, Kokkos::LaunchBounds<256, 4>>(0, num_upper), gamma1_offdiag);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    // std::cout << "KOKKOS Time taken by gamma 1 offdiag: " << diff.count() << " seconds" << std::endl;

    /* gamma2 diag */
    gamma2_cluster(gamma_data_f);
    kokkos_gamma::FunctorGamma2_diag gamma2_diag(*gamma_data_f);
    kokkos_gamma::FunctorGamma2_offdiag gamma2_offdiag(*gamma_data_f);

    start = std::chrono::high_resolution_clock::now();
    parallel_for("KOKKOS gamma2_DIAG", RangePolicy<Kokkos::Cuda, Kokkos::IndexType<int>, Kokkos::LaunchBounds<256, 4>>(0, natoms), gamma2_diag);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    // std::cout << "KOKKOS Time taken by gamma 2 diag: " << diff.count() << " seconds" << std::endl;

    /* gamma2 off-diag */
    start = std::chrono::high_resolution_clock::now();
    parallel_for("gamma2_offDiag", RangePolicy<Kokkos::Cuda, Kokkos::IndexType<int>, Kokkos::LaunchBounds<256, 4>>(0, num_upper), gamma2_offdiag);
    Kokkos::fence();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    // std::cout << "KOKKOS Time taken by gamma 2 off-diag: " << diff.count() << " seconds" << std::endl;

    gamma_data_f->v_c_ewald->template modify<typename gamma_view2d::execution_space>();
    gamma_data_f->v_c_dewald->template modify<typename gamma_view3d::execution_space>();
    gamma_data_f->v_c_atom_ewald->template modify<typename gamma_view2d::execution_space>();
 
    gamma_data_f->v_c_ewald->template sync<typename gamma_view2d::host_mirror_space>();
    gamma_data_f->v_c_dewald->template sync<typename gamma_view3d::host_mirror_space>();
    gamma_data_f->v_c_atom_ewald->template sync<typename gamma_view2d::host_mirror_space>();

    // gamma1_c( *ig1mx,  *ig2mx,  *ig3mx,  *kappa,  **v_c_ratom,  **v_c_ewald,  **v_c_dewald,  **v_c_atom_ewald,  **v_c_g123,  **v_c_Qall,  *volume,  natoms);
    /*** print out c gamma result  ***/
    // kokkos_gamma::gamma_res(*((*gamma_data_f).v_c_ewald), *((*gamma_data_f).v_c_atom_ewald), *((*gamma_data_f).v_c_dewald), (*gamma_data_f).natoms);

    return;
  }
} // extern c