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
! M_case2.f90
! Program Description
! ===========================================================================
!      This is a module calculating the integrals of two centers for the
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
! Module Declaration
! ===========================================================================
        module M_vna_DOGS

! /GLOBAL
        use M_precision

! /SYSTEM
        use M_atom_functions
        use M_species
        use M_integrals_2c

! /CREATE
        use M_vna_HARRIS

! /MPI
        use M_MPI

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains


! ===========================================================================
! initialize_vna_DOGS
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
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_vna_DOGS
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
        integer isorp                       !< counter for shells

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species for
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

! For vna_ontopL_DOGS
            do isorp = 0, species(ispecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            end do

! For vna_ontopR_DOGS
            do isorp = 0, species(jspecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            end do

! For vna_atom_DOGS
            do isorp = 0, species(jspecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            end do
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_vna_DOGS


! ===========================================================================
! vna_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calls the subroutines required to calculate the vna_ontop
! (both left/right cases) and atom cases.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine vna_DOGS
        implicit none

! Parameters and Data Declaration
! ===========================================================================
! None

! Argument Declaration and Description
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
! begin iammaster
        if (my_proc .eq. 0) then
          write (ilogfile,*)
          write (ilogfile,*) ' ******************************************************* '
          write (ilogfile,*) '        C H A R G E D   A T O M    H A R T R E E         '
          write (ilogfile,*) '                  I N T E R A C T I O N S                '
          write (ilogfile,*) ' ******************************************************* '
          write (ilogfile,*)
          write (ilogfile,*) ' Calling Ontop Left case. '
        end if
        call vna_ontopL_DOGS

        if (my_proc .eq. 0) write (ilogfile,*)
        if (my_proc .eq. 0) write (ilogfile,*) ' Calling Ontop Right case. '
        call vna_ontopR_DOGS

        if (my_proc .eq. 0) write (ilogfile,*)
        if (my_proc .eq. 0) write (ilogfile,*) ' Calling Atom case. '
        call vna_atom_DOGS

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ==========================================================================
        return
        end subroutine vna_DOGS


! ===========================================================================
! vna_ontopL_DOGS
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
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine vna_ontopL_DOGS
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
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        real dmax                           !< max distance between two centers
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
! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_rho - 1)

            ! Set final integration limit
            rhomax = min(rcutoff1, rcutoff2)

! Loop over shells
            do isorp = 1, species(ispecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
              nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

              call make_munu (nFdata_cell_2c, ispecies, jspecies)
              nME2c_max = pFdata_cell%nME
              allocate (pFdata_cell%fofx(nME2c_max))

! begin iammaster
              if (my_proc .eq. 0) then
                write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ, isorp
              end if

              ! Open ouput file for this species pair
              write (filename, '("/vna_ontopL_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &          isorp, species(ispecies)%nZ, species(jspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
! Skipping causes issues in parallel if resetting bundle size
              if (skip) cycle
!             if (skip) then
!               pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c - 1
!               pFdata_bundle%nFdata_cell_2c =                                &
!                 pFdata_bundle%nFdata_cell_2c + species(ispecies)%nssh
!               cycle
!             end if
              open (unit = 11, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown')

              ! open directory file
              write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')           &
     &          species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown', position = 'append')
              write (13,100) pFdata_bundle%nFdata_cell_2c, P_vna_ontopL,     &
     &                       isorp, filename(2:30), pFdata_cell%nME, ndd_vna,&
     &                       dmax
              close (unit = 13)

              ! Open mu, nu, mvalue file and write out values.
              write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")')&
     &               P_vna_ontopL, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 12, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown', position = 'append')

             ! write the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%mvalue_2c(index_2c),                 &
     &                      index_2c = 1, nME2c_max)
              close (unit = 12)

! Loop over grid
              d = -drr
              do igrid = 1, ndd_vna
                d = d + drr

                ! Set integration limits
                zmin = max(-rcutoff1, d - rcutoff2)
                zmax = min(rcutoff1, d + rcutoff2)

                call evaluate_integral_2c (nFdata_cell_2c, ispecies,        &
     &                                     jspecies, isorp, ideriv, rcutoff1,&
     &                                     rcutoff2, d, nz_vna, nrho_vna,    &
     &                                     rint_vna_ontopL, phifactor,       &
     &                                     zmin, zmax, rhomin, rhomax,       &
     &                                     pFdata_cell%fofx)
                ! Write out details.
                write (11,*) (pFdata_cell%fofx(index_2c),                    &
	 &                                         index_2c = 1, nME2c_max)
              end do ! igrid
              write (11,*)
            end do ! isorp
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating vna ontopL integrals for nZ = ', i3,        &
     &              ' and nZ = ', i3, ', isorp = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vna_ontopL_DOGS


! ===========================================================================
! vna_ontopR_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(2)|psi2>.  Thus V(2) is located at
! the right site of the orbitals.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine vna_ontopR_DOGS
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
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        real dmax                           !< max distance between two centers
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
! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_rho - 1)

            ! Set final integration limit
            rhomax = min(rcutoff1, rcutoff2)

! Loop over shells
            do isorp = 1, species(jspecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
              nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

              call make_munu (nFdata_cell_2c, ispecies, jspecies)
              nME2c_max = pFdata_cell%nME
              allocate (pFdata_cell%fofx(nME2c_max))

! begin iammaster
              if (my_proc .eq. 0) then
                write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ, isorp
              end if

              ! Open ouput file for this species pair
              write (filename, '("/vna_ontopR_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &               isorp, species(ispecies)%nZ, species(jspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
! Skipping causes issues in parallel if resetting bundle size
              if (skip) cycle
!             if (skip) then
!               pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c - 1
!               pFdata_bundle%nFdata_cell_2c =                                &
!                 pFdata_bundle%nFdata_cell_2c + species(jspecies)%nssh
!               cycle
!             end if
              open (unit = 11, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown')

              ! open directory file
              write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')           &
     &          species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown', position = 'append')
              write (13,100) pFdata_bundle%nFdata_cell_2c, P_vna_ontopR,     &
     &                       isorp, filename(2:30), pFdata_cell%nME, ndd_vna,&
     &                       dmax
              close (unit = 13)

              ! Open mu, nu, mvalue file and write out values.
              write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")')&
     &               P_vna_ontopR, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 12, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown', position = 'append')

              ! write the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%mvalue_2c(index_2c),                 &
     &                      index_2c = 1, nME2c_max)
              close (unit = 12)

! Loop over grid
              d = -drr
              do igrid = 1, ndd_vna
                d = d + drr

                ! Set integration limits
                zmin = max(-rcutoff1, d - rcutoff2)
                zmax = min(rcutoff1, d + rcutoff2)

                call evaluate_integral_2c (nFdata_cell_2c, ispecies,         &
     &                                     jspecies, isorp, ideriv, rcutoff1,&
     &                                     rcutoff2, d, nz_vna, nrho_vna,    &
     &                                     rint_vna_ontopR, phifactor, zmin, &
     &                                     zmax, rhomin, rhomax,             &
     &                                     pFdata_cell%fofx)
                ! Write out details
                write (11,*) (pFdata_cell%fofx(index_2c),                    &
     &                                         index_2c = 1, nME2c_max)
              end do !igrid
              write (11,*)
            end do ! isorp
          end do ! end loop over jspecies
        end do ! end loop over ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating vna ontopR integrals for nZ = ', i3,        &
     &              ' and nZ = ', i3, ', isorp = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vna_ontopR_DOGS


! ===========================================================================
! vna_atom
! ===========================================================================
! SUbroutine Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V|psi2>; the two wavefunctions psi1, and
! psi2 are located on one site and the potential V is located on the other
! site.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! This subroutine then writes the results to file.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! Subroutine Declaration
! ===========================================================================
        subroutine vna_atom_DOGS
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
        integer igrid                       !< number of grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_2c              !< indexing of interactions

        real dmax                           !< max distance between two centers
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
! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + na(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_vna - 1)

            ! Set final integration limit
            rhomax = rcutoff1
            zmin = - rcutoff1
            zmax = rcutoff1

! Loop over shells
            do isorp = 1, species(jspecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
              nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

              call make_munu_atom (nFdata_cell_2c, ispecies, jspecies)
              nME2c_max = pFdata_cell%nME
              allocate (pFdata_cell%fofx(nME2c_max))

! begin iammaster
              if (my_proc .eq. 0) then
                write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ, isorp
              end if

              ! Open ouput file for this species pair
              write (filename, '("/vna_atom_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &          isorp, species(ispecies)%nZ, species(jspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
! Skipping causes issues in parallel if resetting bundle size
              if (skip) cycle
!             if (skip) then
!               pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c - 1
!               pFdata_bundle%nFdata_cell_2c =                                &
!                 pFdata_bundle%nFdata_cell_2c + species(jspecies)%nssh
!               cycle
!             end if
              open (unit = 11, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown')

              ! open directory file
              write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')           &
     &          species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown', position = 'append')
              write (13,100) pFdata_bundle%nFdata_cell_2c, P_vna_atom, isorp,&
     &                       filename(2:30), pFdata_cell%nME, ndd_vna, dmax
              close (unit = 13)

              ! Open mu, nu, mvalue file and write out values.
              write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")')&
     &               P_vna_atom, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 12, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown', position = 'append')

              ! write the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%mvalue_2c(index_2c),                 &
     &                      index_2c = 1, nME2c_max)
              close (unit = 12)

! Loop over grid - first do d = 0.0d0 case, then loop over other distances.
              d = -drr

              ! d = 0.0d0 case - no smoothing function here.
              d = d + drr
              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, nz_vna, nrho_vna, rint_vna_atom, &
     &                                   phifactor, zmin, zmax, rhomin,      &
     &                                   rhomax, pFdata_cell%fofx)
              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                      &
	 &                                       index_2c = 1, nME2c_max)

              do igrid = 2, ndd_vna
                d = d + drr
                call evaluate_integral_2c (nFdata_cell_2c, ispecies,         &
     &                                     jspecies, isorp, ideriv, rcutoff1,&
     &                                     rcutoff2, d, nz_vna, nrho_vna,    &
     &                                     rint_vna_atom, phifactor,         &
     &                                     zmin, zmax, rhomin, rhomax,       &
     &                                     pFdata_cell%fofx)
                ! Write out details.
                write (11,*) (pFdata_cell%fofx(index_2c),                    &
	 &                                         index_2c = 1, nME2c_max)
              end do !igrid
              write (11,*)
            end do ! isorp
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating vna atom integrals for nZ = ', i3,          &
     &              ' and nZ = ', i3, ', isorp = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vna_atom_DOGS


! End Module
! =============================================================================
        end module
