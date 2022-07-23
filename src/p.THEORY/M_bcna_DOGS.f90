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

! Module Declaration
! ===========================================================================
        module M_bcna_DOGS

! /GLOBAL
        use M_precision

! /SYSTEM
        use M_species
        use M_atom_functions
        use M_integrals_3c
        use M_bcna_Harris             ! need for phiint_bcna function

! Type Declaration
! ============================================================================
! three-center interactions arrays

! To cut down on storage space, we actually change the storage procedure
! from previous FIREBALL code. Not all atoms have the same number of
! interactions or interaction types. Before - we would store things based
! on the maximum number of fdata points, maximum number of interactions types,
! maximum number of matrix elements - so even hydrogen-hydrogen (just ss
! and/or ss*) stored a 4x4 or an 8x8 matrix even when not needed.  This was
! quite inefficient.

! The new approach is to define some Fdata types which store the actual
! Fdata points. The smallest unit storage is called Fdata_cell_3C, containing
! all Fdata for a particular interaction/subinteraction type.

! module procedures
        contains

! ===========================================================================
! initialize_bcna_DOGS
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
        subroutine initialize_bcna_DOGS
        implicit none

        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies, kspecies !< counters for number of species
        integer isorp                        !< counter for shells
        integer itheta                       !< counter for theta

        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species for
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

! For bcna_DOGS
            do kspecies = 1, nspecies
              pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
              do isorp = 1, species(kspecies)%nssh
                do itheta = 1, P_ntheta
                  pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1
                end do
              end do
            end do ! kspecies
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
        end subroutine initialize_bcna_DOGS


! ===========================================================================
! bcna_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general three-center
! matrix elements of the form <psi1|V(1)|psi2> for the spherical density.
! ===========================================================================
        subroutine bcna_DOGS
        implicit none

        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ibcba, inaba, itheta    ! looping counters
        integer index_3c, nME3c_max     ! different mu, nu types
        integer iounit                  ! file for writing
        integer isorp                   ! loop over shells
        integer ispmin, ispmax
        integer ispecies, jspecies, kspecies  ! species numbers
        integer nFdata_cell_3c          !< indexing of interactions

        real dbc, dna                  ! distances between centers
        real dbcx, dnax, distance_bc
        real rcutoff1, rcutoff2, rcutoff3   ! atom cutoffs

        real, dimension (P_ntheta) :: ctheta
        real, dimension (P_ntheta) :: ctheta_weights

        real, pointer :: qpl (:, :, :)

        character (len = 30) filename
        character (len = 30) interactions

        logical skip

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*) '        B O N D   C H A R G E D   A T O M   D O G S      '
        write (ilogfile,*) '                  (B C C A) M A T R I X                  '
        write (ilogfile,*) '                  I N T E R A C T I O N S                '
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*)

! Initialize the Legendre coefficients
        call gleg (ctheta, ctheta_weights, P_ntheta)

! Loop over species along the bond charge
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            do kspecies = 1, nspecies
              pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
              pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1
              nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c
              pFdata_cell=>pFdata_bundle%Fdata_cell_3c(nFdata_cell_3c)

              call make_munu_3c (nFdata_cell_3c, ispecies, jspecies, kspecies)
              nME3c_max = pFdata_cell%nME

              ispmin = 1
              ispmax = species(kspecies)%nssh
              allocate (qpl(P_ntheta, nME3c_max, ispmin:(ispmax - ispmin + 1)))
              qpl = 0.0d0

              ! Test ouput file for this species triplet
              itheta = 1
              isorp = 1
              write (filename, '("/", "bcna_", i2.2, "_", i2.2, ".", i2.2,     &
     &                                       ".", i2.2, ".", i2.2, ".dat")')   &
     &               itheta, isorp, species(ispecies)%nZ,                      &
     &                              species(jspecies)%nZ, species(kspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
              if (skip) cycle

              ! Set up grid loop control constants
              rcutoff1 = species(ispecies)%rcutoffA_max
              rcutoff2 = species(jspecies)%rcutoffA_max
              rcutoff3 = species(kspecies)%rcutoffA_max
              dbc = rcutoff1 + rcutoff2
              dna = rcutoff3 + max(rcutoff1,rcutoff2)

              pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c - 1
              do isorp = ispmin, ispmax
                do itheta = 1, P_ntheta
                  pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1

                  write (filename, '("/", "bcna_", i2.2, "_", i2.2, ".", i2.2, &
     &                                    ".", i2.2, ".", i2.2, ".dat")')      &
     &              itheta, isorp, species(ispecies)%nZ, species(jspecies)%nZ, &
     &                             species(kspecies)%nZ

                  ! open directory file
                  write (interactions,                                         &
     &                   '("/3c.",i2.2,".",i2.2,".",i2.2,".dir")')             &
     &              species(ispecies)%nZ, species(jspecies)%nZ,                &
     &              species(kspecies)%nZ
                  open (unit = 13,                                             &
     &                  file = trim(Fdata_location)//trim(interactions),       &
     &                  status = 'unknown', position = 'append')
                  write (13,100) pFdata_bundle%nFdata_cell_3c, P_bcna, isorp,  &
     &                           itheta, filename(2:30), pFdata_cell%nME,      &
     &                           nna_bcna, dna, nbc_bcna, dbc
                  close (unit = 13)

                  ! Open mu, nu, mvalue file and write out values.
                  write (filename, '("/",i2.2, "_munu_3c.",                    &
     &                                   i2.2,".",i2.2,".",i2.2,".dat")')      &
     &              P_bcna, species(ispecies)%nZ, species(jspecies)%nZ,        &
     &                       species(kspecies)%nZ
                  open (unit = 12, file = trim(Fdata_location)//trim(filename),&
     &                  status = 'unknown', position = 'append')

                  ! Write out the mapping - stored in mu, nu, and mvalue
                  write (12,*) (pFdata_cell%mu_3c(index_3c),                   &
     &                                            index_3c = 1, nME3c_max)
                  write (12,*) (pFdata_cell%nu_3c(index_3c),                   &
     &                                                index_3c = 1, nME3c_max)
                  write (12,*) (pFdata_cell%mvalue_3c(index_3c),               &
     &                                                    index_3c = 1, nME3c_max)
                end do
              end do

              write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ, &
     &                             species(kspecies)%nZ

! Open all the output files.
              iounit = 12
              do isorp = ispmin, ispmax
                do itheta = 1, P_ntheta
                  iounit = iounit + 1
                  write (filename, '("/", "bcna_", i2.2, "_", i2.2, ".",       &
     &                                i2.2, ".", i2.2, ".", i2.2, ".dat")')    &
     &                   itheta, isorp, species(ispecies)%nZ,                  &
     &                   species(jspecies)%nZ, species(kspecies)%nZ

! Read the mapping - stored in mu, nu, and mvalue
                  open (unit = (iounit),                                       &
     &                  file = trim(Fdata_location)//trim(filename),           &
     &                  status = 'unknown')
                end do
              end do  ! end isorp loop

! ----------------------------------------------------------------------------
! Begin the big loops over dbc and dna.
! ----------------------------------------------------------------------------
! Loop over all bondcharge distances.
              do ibcba = 1, nbc_rho
                dbcx = float(ibcba - 1)*dbc/float(nbc_bcna - 1)

! for all bondcharges-- we set b=dbcx/2.
                distance_bc = dbcx/2.0d0

! Loop over all neutral atom distances.
! The distance is measured from the bondcharge center (b=dbcx/2)
                do inaba = 1, nna_bcna
                  dnax = float(inaba - 1)*dna/float(nna_bcna - 1)
                  call evaluate_integral_3c (nFdata_cell_3c, ispecies,       &
     &                                       jspecies, kspecies, ispmin,     &
     &                                       ispmax, ctheta, ctheta_weights, &
     &                                       dbcx, dnax, nnr_bcna,           &
     &                                       nntheta_bcna, psiofr, &
     &                                       phiint_bcna, qpl)

! ----------------------------------------------------------------------------
! qpl's are the answer
! ------------------------------------- ---------------------------------------
! Write the qpl coefficients into the data files each combination in1, in2, in3,
! itheta(=1,ntheta_max), isorp gives an individual file.  The values for the
! different non-zero matrix elements of a given combination are written out
! after the index loop.
! ----------------------------------------------------------------------------
                  iounit = 12
                  do isorp = ispmin, ispmax
                    do itheta = 1, P_ntheta
                      iounit = iounit + 1
                      write (iounit,*)                                       &
     &                  (qpl(itheta,index_3c,isorp), index_3c = 1, nME3c_max)
                    end do
                  end do  ! end isporp loop
                end do  ! end of the dna loop
              end do  ! the end of the dbc loop

              ! close files
              iounit = 12
              do isorp = ispmin, ispmax
                do itheta = 1, P_ntheta
                  iounit = iounit + 1
                  close (unit = iounit)
                end do
                deallocate (qpl)
              end do  ! end isporp loop
            end do  ! end loop over kspecies
          end do  ! end loop over jspecies
        end do  ! end loop over ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i2, 1x, i2, 1x, i2, 1x, i3, 1x, a29, 1x, i3,             &
     &          1x, i4, 1x, f9.6, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating bcna integrals for nZ = ', i3,              &
     &              ' and nZ = ', i3, ', potential on nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine bcna_DOGS

! End Module
! =============================================================================
        end module
