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

! M_integrals_2c.f90
! Program Description
! ============================================================================
!      This is a module calculating the integrals of two centers
!
! ============================================================================
! Code written by:
! Hong Wang
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
!
! Changes by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! ============================================================================
! Module Declaration
! ============================================================================
       module M_integrals_1c

! /GLOBAL
       use M_precision

! /SYSTEM
       use M_atom_functions
       use M_species

       include '../include/interactions_2c.h'

! Type Declaration
! ============================================================================
! Two-center interactions arrays

! To cut down on storage space, we actually change the storage procedure
! from previous FIREBALL code. Not all atoms have the same number of
! interactions or interaction types. Before - we would store things based
! on the maximum number of fdata points, maximum number of interactions types,
! maximum number of matrix elements - so even hydrogen-hydrogen (just ss
! and/or ss*) stored a 4x4 or an 8x8 matrix even when not needed.  This was
! quite inefficient.

! The new approach is to define some Fdata types which store the actual
! Fdata points. The smallest unit storage is called Fdata_cell_2C, containing
! all Fdata for a particular interaction/subinteraction type.
!       type T_Fdata_cell_2c
!         integer nME                   ! number of non-zero matrix elements

!         integer, pointer :: mu_2c (:)
!         integer, pointer :: nu_2c (:)
!         integer, pointer :: mvalue_2c (:)

!         integer, pointer :: N_mu (:)
!         integer, pointer :: N_nu (:)
!         integer, pointer :: L_mu (:)
!         integer, pointer :: L_nu (:)
!         integer, pointer :: M_mu (:)
!          integer, pointer :: M_nu (:)

!         real, pointer :: dx (:)                ! distance between data points
!         real, pointer :: xmax (:)              ! maximum interaction range

          ! actual Fdata points f(x)
!         real, pointer :: fofx (:)
!        end type T_Fdata_cell_2c

! Fdata_bundle_2c is the 2-center package of Fdata_cell_2c. It contains all the
! Fdata_cell_2c information of all interaction/subtypes for given species pair.
!       type T_Fdata_bundle_2c
!         integer index_2c_overlapS                  ! location of overlapS
!         integer nFdata_cell_2c                     ! number of Fdata_cell2C

          ! actual fdata
!         type (T_Fdata_cell_2c), pointer :: Fdata_cell_2c (:)
!       end type T_Fdata_bundle_2c

!       type(T_Fdata_bundle_2c), pointer :: Fdata_bundle_2c (:,:)

        real, parameter :: tolerance = 1d0-6  ! tolerance for Adaptive Simpsons

! module procedures
        contains


! ===========================================================================
! evaluate_integral_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      Here we integrate the two-center integral - first along the z-direction
! (zint) which is the direction along the bond and is determined by the
! distance between the two centers. Second, we integrate in the "r" direction
! (rint) which is the distance from the first center. Thus, we perform an
! integral in cylindrical coordinates.
!
! ===========================================================================
        subroutine evaluate_integral_1c (itype, ispecies, jspecies, isorp,   &
     &                                   ideriv, rcutoff1, rcutoff2, d,      &
     &                                   nz, nrho, rint, phunction, zmin,    &
     &                                   zmax, rhomin, rhomax, Fdata)
        implicit none

! Auguments Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itype
        integer, intent (in) :: ispecies, jspecies      ! two species
        integer, intent (in) :: isorp
        integer, intent (in) :: ideriv
        integer, intent (in) :: nz             ! number of intergration points
        integer, intent (in) :: nrho

        real, intent (in) ::  d
        real, intent (in) :: rcutoff1, rcutoff2    ! cutoffs for wavefunctions
        real, intent (in) :: zmin, zmax, rhomin, rhomax  !limits of integration

        interface
          function rint (itype, ispecies, jspecies, isorp, d, rho, z1, z2, &
     &                     ideriv, index_2c)
            integer, intent(in) :: ispecies, jspecies, isorp, itype, ideriv, &
     &                             index_2c
            real, intent(in) :: rho, z1, z2, d
            real rint
          end function rint
        end interface

        interface
          function phunction (l1, m1, l2, m2)
            integer, intent (in) :: l1, m1, l2, m2
            real phunction
          end function phunction
        end interface

! Output
        real, pointer :: Fdata (:)

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration adn Description
! ===========================================================================
        integer iz, nnz        ! counting over z-direction grid points up to nnz
        integer lmu, mmu, nmu  ! mu quantum numbers
        integer lnu, mnu, nnu  ! nu quantum numbers
        integer nME2c_max      ! number of matrix elements
        integer index_2c       ! counter 1, ..., nME2c_max

        real z1
        real dz
        real temp

        ! Simpsons Quadrature Variables.
        real, allocatable, dimension (:) :: zmult

        real xntegral          ! the answer

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)
        nME2c_max = pFdata_cell%nME

! Loop over the matrix elements:
        do index_2c = 1, nME2c_max
          nmu = pFdata_cell%N_mu(index_2c)
          lmu = pFdata_cell%L_mu(index_2c)
          mmu = pFdata_cell%M_mu(index_2c)
          nnu = pFdata_cell%N_nu(index_2c)
          lnu = pFdata_cell%L_nu(index_2c)
          mnu = pFdata_cell%M_nu(index_2c)

! Initialize the sum to zero
          xntegral = 0.0d0

! Integration is over z (z-axis points from atom 1 to atom 2) and rho (rho is
! radial distance from z-axis).

!----------------------------------------------------------------------------
! Non Adaptive Simpson's Setup
!----------------------------------------------------------------------------
! Strictly define what the density of the mesh should be.  Make the density
! of the number of points equivalent for all cases. Change the number of
! points to be integrated to be dependent upon the distance between the
! centers and this defined density.
          dz = ((rcutoff1 + rcutoff2)/2.0d0)/float(nz)
          nnz = int((zmax - zmin)/dz)
          if (mod(nnz,2) .eq. 0) nnz = nnz + 1; allocate (zmult (nnz))

! Set up Simpson's rule factors. First for the z integration and then for
! the rho integration.
          zmult(1) = dz/3.0d0; zmult(nnz) = dz/3.0d0
          do iz = 2, nnz - 1, 2
            zmult(iz) = 4.0d0*dz/3.0d0
          end do
          do iz = 3, nnz - 2, 2
            zmult(iz) = 2.0d0*dz/3.0d0
          end do

!----------------------------------------------------------------------------
! Actually use the normal Simpson. Set up z integration here, and rho
! integration the zint function.
!----------------------------------------------------------------------------
          xntegral = 0.00
          do iz = 1, nnz
            z1 = zmin + float(iz-1)*dz
            temp = zint (itype, ispecies, jspecies, isorp, z1, ideriv,       &
     &                   index_2c, d, nrho, rint, rcutoff1, rcutoff2,        &
     &                   rhomin, rhomax)
            temp = temp * zmult(iz)
            xntegral = xntegral + temp
          end do

! This factor (phifactor) comes from a result of multiplying the two Ylm's
! together, after the factor of pi is multiplied out after the phi
! integration.
          xntegral = phunction(lmu,mmu,lnu,mnu)*xntegral
          Fdata(index_2c) = xntegral

! deallocate
          deallocate (zmult)
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine evaluate_integral_1c


! End Module
! =============================================================================
        end module M_integrals_2c
