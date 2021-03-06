! copyright info:
!
!                             @Copyright 2022
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
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! M_kpoints
! Module Description
! ===========================================================================
!>       This is a module which will either read in the kpoints from a file
!! provided by the user - necessary for doing periodic boundary conditions.
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
! Module Declaration
! ===========================================================================
        module M_kpoints

! Type Declaration
! ===========================================================================
! If we are interested in transitions, such as for the evolving the
! coefficients in order to calculate non-adiabatic couplings or for calculating
! absorption. We consider only a range of states near the valence and conduction
! band edges. The number of states that we consider is defined by ntransitions.
        integer ntransitions
        type T_transition
          integer imap                 ! to which band does this transition map

          ! these are the time-dependent Schrodinger coefficients - these will
          ! be evolved in time according to Verlet or Runge-Kutta as desired
          ! these correspond to the non-adiabatic coupled states
          complex cna_old
          complex cna

          ! these are the coefficients of the adiabatic eigenstates - these will
          ! be evolved in time similarly to cna and cna_old
          complex, allocatable :: c_wf (:)

          ! this is the non-adiabatic coupling belonging to the transition state
          complex, allocatable :: djk (:)
          complex, allocatable :: djk_old (:)
          complex, allocatable :: ddjk (:)
        end type

        type T_kpoint
          real weight                          ! weight of kpoint

          integer, pointer :: ioccupy (:)    ! integer occupation number

          real, pointer :: eigen (:)          ! eigenvalues for k
          real, pointer :: eigen_old (:)      ! previous eigenvlues for k
          real, pointer :: deigen (:)         ! interpolted eigen values
          real, pointer :: foccupy (:)        ! occupation real value for k

          real, dimension (3) :: k           ! kpoint vector

          ! the real and imaginary components of the eigenvectors
          complex, pointer :: c (:, :)

          ! the real and imaginary components of the Lowdin transformed eigenvectors
          complex, pointer :: c_Lowdin (:, :)

! If we are interested in transitions, such as for the evolving the
! coefficients in order to calculate non-adiabatic couplings or for calculating
! absorption. We consider only a range of states near the valence and conduction
! band edges. The number of states that we consider is defined by ntransitions.
          type (T_transition), pointer :: transition (:)
        end type

! End Module
! ===========================================================================
        end module M_kpoints
