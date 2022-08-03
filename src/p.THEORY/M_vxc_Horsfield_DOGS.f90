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

! M_vxc_Harris
! Program Description
! ===========================================================================
!      This is a module calculating the integrals of two centers for the
! the exchange-correlation interactions.
!
! We use a finite difference approach to change the densities and then
! find the corresponding changes in the exchange-correlation potential.
! This a bit clumsy sometimes, and probably can use a rethinking to improve.
! However, it actually works quite well for many molecular systems.
!
! We compute neutral cases for ideriv = 1. For other ideriv's we have the
! following KEY:
!
! Case 1 (KEY=1), neutral neutral corresponds to (00)
! KEY = 1,2,3,4,5 for ideriv = 1,2,3,4,5
!
! For the one-center case we only use ideriv = 1,2,3 as there is only one atom.
!
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Module Declaration
! ===========================================================================
        module M_vxc_DOGS

! /GLOBAL
        use M_precision

! /SYSTEM
        use M_atom_functions
        use M_atomPP_functions
        use M_atomPP_ion_functions
        use M_species
        use M_integrals_2c

! /CREATE
        use M_vxc_Harris
        use M_xc_1c
        use M_xc_2c

! Type Declaration
! ===========================================================================
! The charge transfer bit = need a +-dq on each orbital.
        real, allocatable :: dqorb (:)

! Some parameters for the derivative parts:
        integer, parameter, dimension (0:4) :: jsign = (/0, -1, +1, -1, +1/)
        integer, parameter, dimension (0:4) :: kspecies_key = (/1, 1, 1, 2, 2/)

! module procedures
        contains

! ===========================================================================
! initialize_vxc_DOGS
! ===========================================================================
! Program Description
! ===========================================================================
!       We need to determine how many interactions belong to each nspecies
! bundle pair. This routine just counts how many total interactions contribute
! to that bundle.  Something like overlap is obviously only 1 interaction
! added, but something like vna needs number of interactions based on the
! number of shells.
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_vxc_DOGS
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species
        do ispecies = 1, nspecies
          pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)

! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
          ideriv_min = 1
          ideriv_max = 2

! For the once center case we only do +- dq changes in the density.
          do ideriv = ideriv_min, ideriv_max
            ! add one for vxc_1c
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
          end do

          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
            ideriv_min = 1
            ideriv_max = 4

! For the once center case we only do +- dq changes in the density.
            do ideriv = ideriv_min, ideriv_max

              ! add one for uxc
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1

              ! add one for vxc_ontop
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1

              ! add one for vxc_atom
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            end do ! ideriv
          end do ! jspecies
        end do ! ispecies

! Initialize the charge transfer bit
        allocate (dqorb (nspecies))
        do ispecies = 1, nspecies
          dqorb(ispecies) = 0.5d0
          if (species(ispecies)%nssh .eq. 1) dqorb(ispecies) = 0.25d0
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
        end subroutine initialize_vxc_DOGS


! ===========================================================================
! vxc_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calls the subroutines required to calculate the vna_ontop
! (both left/right cases) and atom cases.  Also, the over-counting correction
! to the exchange-correlation energy is calculated here.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine vxc_DOGS
        implicit none

! Parameters and Data Declaration
! ===========================================================================
! None

! ===========================================================================
! Argument Declaration and Description
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*) '        E X C H A N G E   C O R R E L A T I O N          '
        write (ilogfile,*) '        C H A R G E D   I N T E R A C T I O N S          '
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*)

        write (ilogfile,*) ' Calling one-center case. '
        call vxc_1c_DOGS

        write (ilogfile,*)
        write (ilogfile,*) ' Calling correction case. '
        call uxc_DOGS

        write (ilogfile,*)
        write (ilogfile,*) ' Calling ontop case. '
        call vxc_ontop_DOGS

        write (ilogfile,*)
        write (ilogfile,*) ' Calling atom case. '
        call vxc_atom_DOGS

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ==========================================================================
        return
        end subroutine vxc_DOGS


! ===========================================================================
! vxc_1c_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the one-center
! matrix elements of the form <psi1|V(1)|psi2>.  The potential V(1) and
! wavefunctions, psi1 and psi2, are located at the same atom site.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
!
! Also, the energy piece is calculated which is needed for the correction
! terms in the USR part of the total energy.
!
! We calculate   (n1+n2)*(exc(1+2)-muxc(1+2)) - n1*(exc(1)-xcmu(1))
!                                             - n2*(exc(2)-xcmu(2))
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine vxc_1c_DOGS
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies                    !< counters for number of species
        integer issh, jssh
        integer iexc                        ! which type of exchange-correlation
        integer index_1c, nME1c_max         !< basically the number of non-zero

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        integer nFdata_cell_1c              !< indexing of interactions

        real dqint                          !< delta of charge in shell
        real rhomin, rhomax
        real xc_fraction

! Simpsons' Variables.
        integer irho, nnrho               ! integrate from irho to nnrho

        real xntegral
        real rho                          ! variable of integration

        real drho                         ! interval between rho points
        real temprho

        ! Simpsons Quadrature Variables.
        real, allocatable, dimension (:) :: rhomult

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 25) filename

! Procedure
! ============================================================================
! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
        ideriv_min = 1
        ideriv_max = 2

! Loop over species
        do ispecies = 1, nspecies

! For the once center case we only do +- dq changes in the density.
          do ideriv = ideriv_min, ideriv_max

            ! change the charge state to affect the density
            dqint = dqorb(ispecies)/species(ispecies)%nssh
            do issh = 1, species(ispecies)%nssh
              species(ispecies)%shell(issh)%Qneutral_ion =                    &
     &          species(ispecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
            end do

            ! cut some lengthy notation with pointers
            pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_1c = pFdata_bundle%nFdata_cell_2c
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_1c)

            call make_munuS (nFdata_cell_1c, ispecies, ispecies)
            nME1c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME1c_max))

            ! Open ouput file for this species pair
            write (filename, '("/vxc_1c", ".", i2.2, ".dat")')                &
     &        species(ispecies)%nZ
            open (unit = 11, file = trim(Fdata_location)//trim(filename),     &
     &            status = 'unknown', position = 'append')

            ! Set integration limits
            rhomin = 0.0d0
            rhomax = species(ispecies)%rcutoffA_max

! Loop over grid
            write (ilogfile,100) species(ispecies)%nZ

! ***************************************************************************
! exc_1c piece:
! We do not have shells here, just a single energy.
! ***************************************************************************
! Find which exchange-correlation we are calcuating:
            iexc = species_PP(ispecies)%iexc
            xc_fraction = species_PP(ispecies)%xc_fraction

! Non Adaptive Simpson's Setup
! ---------------------------------------------------------------------------
! Strictly define what the density of the mesh should be.  Make the density of
! the number of points equivalent for all cases. Change the number of points
! to be integrated to be dependent upon the distance between the centers and
! this defined density.
            drho = rhomax/float(nrho_rho_store)
            nnrho = int((rhomax - rhomin)/drho)
            if (mod(nnrho,2) .eq. 0) nnrho = nnrho + 1; allocate (rhomult(nnrho))

            rhomult(1) = drho/3.0d0; rhomult(nnrho) = drho/3.0d0
            do irho = 2, nnrho - 1, 2
              rhomult(irho) = 4.0d0*drho/3.0d0
            end do
            do irho = 3, nnrho - 2, 2
              rhomult(irho) = 2.0d0*drho/3.0d0
            end do

! Actually use the normal Simpson. Set up rho integration
!----------------------------------------------------------------------------
            xntegral = 0.00
            do irho = 1, nnrho
              rho = rhomin + float(irho - 1)*drho
              temprho = 4.0d0*pi*dexc_1c (iexc, xc_fraction, ispecies, rho)*rho**2
              temprho = temprho * rhomult(irho)
              xntegral = xntegral + temprho
            end do

            ! Write out details.
            index_1c = 1
            write (11,*) xntegral

! ***************************************************************************
! vxc_1c piece:
! We do not have shells here, just a single energy.
! ***************************************************************************
! Loop over the matrix elements:
            do index_1c = 1, nME1c_max
              xntegral = 0.00
              do irho = 1, nnrho
                rho = rhomin + float(irho - 1)*drho
                temprho = rint_vxc_1c (nFdata_cell_1c, ispecies, rho, index_1c)
                temprho = temprho * rhomult(irho)
                xntegral = xntegral + temprho
              end do
              pFdata_cell%fofx(index_1c) = xntegral
            end do

            ! Write out details.
            index_1c = 1
            do issh = 1, species(ispecies)%nssh
              write (11,*) (pFdata_cell%fofx(jssh),                           &
     &                      jssh = index_1c, index_1c + species(ispecies)%nssh - 1)
              index_1c = index_1c + species(ispecies)%nssh
            end do
            write (11,*)

            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c - 1

            ! Deallocate rhomult so we can recalculate for the next species
            deallocate (rhomult)
          end do ! derivatives for charge transfer
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Evaluating vxc_1c (ideriv) integrals for nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vxc_1c_DOGS


! ===========================================================================
! uxc_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the exchange-correlation correction term used
! in the short range (usr) interactions.
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine uxc_DOGS
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer issh
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp                       !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        real dmax                           !< max distance between two centers
        real dqint                          !< delta of charge in shell
        real drr                            !< distance between mesh points
        real d                              !< distance between the two centers
        real rcutoff1, rcutoff2             !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 30) filename
        character (len = 25) interactions

        logical skip

! Procedure
! ============================================================================
! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
        ideriv_min = 1
        ideriv_max = 4

! We do not need isorp for Horsfield exchange-correlation, so set isorp = 0
        isorp = 0

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

! For the once center case we only do +- dq changes in the density.
            do ideriv = ideriv_min, ideriv_max

              ! Set the values of Qneutral_ion to the original Qneutral
              do issh = 1, species(ispecies)%nssh
                species(ispecies)%shell(issh)%Qneutral_ion =                  &
     &            species(ispecies)%shell(issh)%Qneutral
              end do
              do issh = 1, species(jspecies)%nssh
                species(jspecies)%shell(issh)%Qneutral_ion =                  &
     &            species(jspecies)%shell(issh)%Qneutral
              end do

              ! change the charge state to affect the density
              if (kspecies_key(ideriv) .eq. 1) then
                dqint = dqorb(ispecies)/species(ispecies)%nssh
                do issh = 1, species(ispecies)%nssh
                  species(ispecies)%shell(issh)%Qneutral_ion =                &
     &              species(ispecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
                end do
              else
                dqint = dqorb(jspecies)/species(jspecies)%nssh
                do issh = 1, species(jspecies)%nssh
                  species(jspecies)%shell(issh)%Qneutral_ion =                &
     &              species(jspecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
                end do
              end if

! ****************************************************************************
! Calling rho2c_store for new densities.
! Because we change the charge states for the densities, then we need to
! recalculate the density for the ispecies, jspecies pair and store again.
! ****************************************************************************
              call rho_2c_store (ispecies, jspecies)

              ! cut some lengthy notation with pointers
              pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
              nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

              ! This correction term does not have matrix elements;
              ! set nME2c_max = 1
              call make_munu_uxc (nFdata_cell_2c, ispecies, jspecies)
              nME2c_max = pFdata_cell%nME
              allocate (pFdata_cell%fofx(nME2c_max))

              ! Open ouput file for this species pair
              write (filename, '("/uxc_", i2.2,".",i2.2,".",i2.2,".dat")')    &
     &               ideriv, species(ispecies)%nZ, species(jspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
              if (skip) cycle
              open (unit = 11, file = trim(Fdata_location)//trim(filename),   &
     &              status = 'unknown')

              ! Set up grid loop control constants
              rcutoff1 = species(ispecies)%rcutoffA_max
              rcutoff2 = species(jspecies)%rcutoffA_max
              dmax = wf(ispecies)%rcutoffA_max + na(jspecies)%rcutoffA_max
              drr = dmax/float(ndd_uxc - 1)
              d = -drr

              ! Set final integration limit
              rhomax = min(rcutoff1, rcutoff2)

              ! open directory file
              write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')            &
     &          species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown', position = 'append')
              write (13,100) pFdata_bundle%nFdata_cell_2c, P_uxc, ideriv,     &
     &                       filename(2:30), pFdata_cell%nME, ndd_uxc, dmax
              close (unit = 13)

              ! Open mu, nu, mvalue file and write out values.
              write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")')&
     &               P_uxc, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
     &              status = 'unknown', position = 'append')

              ! write the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%mvalue_2c(index_2c),                  &
     &                      index_2c = 1, nME2c_max)

! Loop over grid
              write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ, ideriv
              do igrid = 1, ndd_uxc
                d = d + drr

                ! Set integration limits
                zmin = max(-rcutoff1, d - rcutoff2)
                zmax = min(rcutoff1, d + rcutoff2)

                call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies,&
     &                                     isorp, ideriv, rcutoff1, rcutoff2, &
     &                                     d, nz_uxc, nrho_uxc, rint_uxc,     &
     &                                     phifactor, zmin, zmax,             &
     &                                     rhomin, rhomax, pFdata_cell%fofx)

                ! Write out details.
                write (11,*) (pFdata_cell%fofx(index_2c),                     &
     &                                         index_2c = 1, nME2c_max)
              end do ! igrid
              write (11,*)

            end do ! ideriv
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating uxc integrals for nZ = ', i3,                &
     &              ' and nZ = ', i3, ', ideriv = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine uxc_DOGS


! ===========================================================================
! vxc_ontop_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the left site of one of the orbitals.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine vxc_ontop_DOGS
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer issh
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp                       !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        real dmax                           !< max distance between two centers
        real dqint                          !< delta of charge in shell
        real drr                            !< distance between mesh points
        real d                              !< distance between the two centers
        real rcutoff1, rcutoff2             !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 30) filename
        character (len = 25) interactions

        logical skip

! Procedure
! ============================================================================
! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
        ideriv_min = 1
        ideriv_max = 4

! We do not need isorp for Horsfield exchange-correlation, so set isorp = 0
        isorp = 0

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

! For the once center case we only do +- dq changes in the density.
            do ideriv = ideriv_min, ideriv_max

              ! Set the values of Qneutral_ion to the original Qneutral
              do issh = 1, species(ispecies)%nssh
                species(ispecies)%shell(issh)%Qneutral_ion =                  &
     &            species(ispecies)%shell(issh)%Qneutral
              end do
              do issh = 1, species(jspecies)%nssh
                species(jspecies)%shell(issh)%Qneutral_ion =                  &
     &            species(jspecies)%shell(issh)%Qneutral
              end do

              ! change the charge state to affect the density
              if (kspecies_key(ideriv) .eq. 1) then
                dqint = dqorb(ispecies)/species(ispecies)%nssh
                do issh = 1, species(ispecies)%nssh
                  species(ispecies)%shell(issh)%Qneutral_ion =                &
     &              species(ispecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
                end do
              else
                dqint = dqorb(jspecies)/species(jspecies)%nssh
                do issh = 1, species(jspecies)%nssh
                  species(jspecies)%shell(issh)%Qneutral_ion =                &
     &              species(jspecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
                end do
              end if

! ****************************************************************************
! Calling rho2c_store for new densities.
! Because we change the charge states for the densities, then we need to
! recalculate the density for the ispecies, jspecies pair and store again.
! ****************************************************************************
              call rho_2c_store (ispecies, jspecies)

              ! cut some lengthy notation with pointers
              pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
              nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

              call make_munu (nFdata_cell_2c, ispecies, jspecies)
              nME2c_max = pFdata_cell%nME
              allocate (pFdata_cell%fofx(nME2c_max))

              ! Open ouput file for this species pair
              write (filename, '("/vxc_ontop_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &  	         ideriv, species(ispecies)%nZ, species(jspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
              if (skip) cycle
              open (unit = 11, file = trim(Fdata_location)//trim(filename),   &
     &            status = 'unknown')

              ! Set up grid loop control constants
              rcutoff1 = species(ispecies)%rcutoffA_max
              rcutoff2 = species(jspecies)%rcutoffA_max
              dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
              drr = dmax/float(ndd_vxc - 1)
              d = -drr

              ! Set integration limits
	          rhomin = 0.0d0
    	      rhomax = min(rcutoff1, rcutoff2)

              ! open directory file
              write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')            &
     &          species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown', position = 'append')
              write (13,100) pFdata_bundle%nFdata_cell_2c, P_vxc_ontop, ideriv,&
     &                       filename(2:30), pFdata_cell%nME, ndd_vxc, dmax
              close (unit = 13)

              ! Open mu, nu, mvalue file and write out values.
              write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")')&
     &               P_vxc_ontop, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
     &              status = 'unknown', position = 'append')

              ! write the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%mvalue_2c(index_2c),                  &
     &                      index_2c = 1, nME2c_max)

! Loop over grid
              write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ, ideriv
              do igrid = 1, ndd_vxc
                d = d + drr

                ! Set integration limits
                zmin = max(-rcutoff1, d - rcutoff2)
                zmax = min(rcutoff1, d + rcutoff2)

                call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies,&
     &                                     isorp, ideriv, rcutoff1, rcutoff2, &
     &                                     d, nz_vxc, nrho_vxc,               &
     &                                     rint_vxc_ontop, phifactor, zmin,   &
     &                                      zmax, rhomin, rhomax, pFdata_cell%fofx)

                ! Write out details.
                write (11,*) (pFdata_cell%fofx(index_2c),                     &
     &                                         index_2c = 1, nME2c_max)
              end do ! igrid
              write (11,*)
            end do ! ideriv
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating vxc ontop integrals for nZ = ', i3,          &
     &              ' and nZ = ', i3, ', ideriv = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vxc_ontop_DOGS


! ===========================================================================
! vxc_atom_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the left site of one of the orbitals.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine vxc_atom_DOGS
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer issh
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp                       !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        real dmax                           !< max distance between two centers
        real dqint                          !< delta of charge in shell
        real drr                            !< distance between mesh points
        real d                              !< distance between the two centers
        real rcutoff1, rcutoff2             !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 30) filename
        character (len = 25) interactions

        logical skip

! Procedure
! ============================================================================
! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
        ideriv_min = 1
        ideriv_max = 4

! We are doing only Harris here, so set isorp = 0
        isorp = 0

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

! For the once center case we only do +- dq changes in the density.
            do ideriv = ideriv_min, ideriv_max

              ! Set the values of Qneutral_ion to the original Qneutral
              do issh = 1, species(ispecies)%nssh
                species(ispecies)%shell(issh)%Qneutral_ion =                  &
     &            species(ispecies)%shell(issh)%Qneutral
              end do
              do issh = 1, species(jspecies)%nssh
                species(jspecies)%shell(issh)%Qneutral_ion =                  &
     &            species(jspecies)%shell(issh)%Qneutral
              end do

              ! change the charge state to affect the density
              if (kspecies_key(ideriv) .eq. 1) then
                dqint = dqorb(ispecies)/species(ispecies)%nssh
                do issh = 1, species(ispecies)%nssh
                  species(ispecies)%shell(issh)%Qneutral_ion =                &
     &              species(ispecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
                end do
              else
                dqint = dqorb(jspecies)/species(jspecies)%nssh
                do issh = 1, species(jspecies)%nssh
                  species(jspecies)%shell(issh)%Qneutral_ion =                &
     &              species(jspecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
                end do
              end if

! ****************************************************************************
! Calling rho2c_store for new densities.
! Because we change the charge states for the densities, then we need to
! recalculate the density for the ispecies, jspecies pair and store again.
! ****************************************************************************
              call rho_2c_store (ispecies, jspecies)

              ! cut some lengthy notation with pointers
              pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
              nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

              call make_munu_atom (nFdata_cell_2c, ispecies, jspecies)
              nME2c_max = pFdata_cell%nME
              allocate (pFdata_cell%fofx(nME2c_max))

              ! Open ouput file for this species pair
              write (filename, '("/vxc_atom_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &               ideriv, species(ispecies)%nZ, species(jspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
              if (skip) cycle
              open (unit = 11, file = trim(Fdata_location)//trim(filename),   &
     &              status = 'unknown')

              ! Set up grid loop control constants
              rcutoff1 = species(ispecies)%rcutoffA_max
              rcutoff2 = species(jspecies)%rcutoffA_max
              dmax = wf(ispecies)%rcutoffA_max + na(jspecies)%rcutoffA_max
              drr = dmax/float(ndd_vxc - 1)
              d = -drr

              ! Set final integration limit
              rhomax = min(rcutoff1, rcutoff2)

              ! open directory file
              write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')            &
     &          species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown', position = 'append')
              write (13,100) pFdata_bundle%nFdata_cell_2c, P_vxc_atom, ideriv,&
     &                       filename(2:30), pFdata_cell%nME, ndd_vxc, dmax
              close (unit = 13)

              ! Open mu, nu, mvalue file and write out values.
              write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")')&
     &               P_vxc_atom, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 12, file = trim(Fdata_location)//trim(filename),   &
     &              status = 'unknown', position = 'append')

              ! write the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%mvalue_2c(index_2c),                  &
     &                      index_2c = 1, nME2c_max)

! Loop over grid
              write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ, ideriv
              do igrid = 1, ndd_vxc
                d = d + drr

                ! Set integration limits
                zmin = max(-rcutoff1, d - rcutoff2)
                zmax = min(rcutoff1, d + rcutoff2)

                call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies,&
     &                                     isorp, ideriv, rcutoff1, rcutoff2, &
     &                                     d, nz_vxc, nrho_vxc,               &
     &                                     rint_vxc_atom, phifactor, zmin,    &
     &                                     zmax, rhomin, rhomax, pFdata_cell%fofx)

                ! Write out details.
                write (11,*) (pFdata_cell%fofx(index_2c),                       &
     &                                         index_2c = 1, nME2c_max)
              end do ! igrid
              write (11,*)
            end do ! ideriv
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating vxc atom integrals for nZ = ', i3,           &
     &              ' and nZ = ', i3, ', ideriv = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vxc_atom_DOGS

! End Module
! =============================================================================
      end module
