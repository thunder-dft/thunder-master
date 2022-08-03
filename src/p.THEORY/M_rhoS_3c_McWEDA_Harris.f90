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
        module M_rhoS_3c_Harris

! /GLOBAL
        use M_precision

! /SYSTEM
        use M_species
        use M_atom_functions
        use M_integrals_3c

        ! need for overlapS to weight the densities
        use M_integrals_2c
        use M_rhoS_2c_Harris

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
! Fdata points. The smallest unit storage is called Fdata_cell_2C, containing
! all Fdata for a particular interaction/subinteraction type.

! module procedures
        contains

! ===========================================================================
! initialize_rho_3c_Harris
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
        subroutine initialize_rhoS_3c_HARRIS
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
        integer ispecies, jspecies, kspecies ! counters for number of species
        integer isorp                        ! counter for shells
        integer itheta                       ! counter for theta

        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species for
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

! For rho_3c_Harris
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
        end subroutine initialize_rhoS_3c_HARRIS


! ===========================================================================
! rhoS_3c_Harris
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general three-center
! matrix elements of the form <psi1|V(1)|psi2> for the spherical density.
! ===========================================================================
        subroutine rhoS_3c_Harris
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
        integer index_2c_overlapS       !< where overlapS starts
        integer iounit                  ! file for writing
        integer isorp, ideriv           ! loop over shells
        integer ispmin, ispmax
        integer ispecies, jspecies, kspecies  ! species numbers
        integer nFdata_cell_3c          !< indexing of interactions

        real dbc, dna                  ! distances between centers
        real dbcx, dnax, distance_bc
        real rcutoff1, rcutoff2, rcutoff3   ! atom cutoffs

        ! integration limits for overlapS
        real zmin, zmax
        real rhomin, rhomax

        real, dimension (P_ntheta) :: ctheta
        real, dimension (P_ntheta) :: ctheta_weights

        real, allocatable :: qpl (:, :, :)

        character (len = 30) filename
        character (len = 30) interactions

        logical skip

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! for overlapS
        type (T_Fdata_cell_2c), pointer :: pFdata_overlapS_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle_2c

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*) '      S P H E R I C A L   D E N S I T Y    M A T R I X   '
        write (ilogfile,*) '                 I N T E R A C T I O N S                 '
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*)

! Initialize the Legendre coefficients
        call gleg (ctheta, ctheta_weights, P_ntheta)

! Loop over species along the bond charge
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

! Set up stuff for overlapS interactions - we need to later divide the density by
! this spherical overlap in order to weight the density for McWEDA.
            pFdata_bundle_2c=>Fdata_bundle_2c(ispecies, jspecies)
            index_2c_overlapS = pFdata_bundle_2c%index_2c_overlapS
            pFdata_overlapS_cell=>pFdata_bundle_2c%Fdata_cell_2c(index_2c_overlapS)

            do kspecies = 1, nspecies
              pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
              pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1
              nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c
              pFdata_cell=>pFdata_bundle%Fdata_cell_3c(nFdata_cell_3c)

              call make_munuS_3c (nFdata_cell_3c, ispecies, jspecies, kspecies)
              nME3c_max = pFdata_cell%nME

              ispmin = 1
              ispmax = species(kspecies)%nssh
              if (allocated (qpl)) deallocate (qpl)
              allocate (qpl (P_ntheta, nME3c_max, ispmin:(ispmax - ispmin + 1)))
              qpl = 0.0d0

              ! Test ouput file for this species triplet
              itheta = 1
              isorp = 1
              write (filename, '("/", "rhoS_3c_", i2.2, "_", i2.2, ".", i2.2, &
     &                                       ".", i2.2, ".", i2.2, ".dat")')  &
     &          itheta, isorp, species(ispecies)%nZ,                          &
     &                         species(jspecies)%nZ, species(kspecies)%nZ
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

                  write (filename, '("/", "rhoS_3c_", i2.2, "_", i2.2, ".",   &
     &                               i2.2, ".", i2.2, ".", i2.2, ".dat")')    &
     &              itheta, isorp, species(ispecies)%nZ,                      &
     &                             species(jspecies)%nZ, species(kspecies)%nZ

                  ! open directory file
                  write (interactions,                                        &
     &                   '("/3c.",i2.2,".",i2.2,".",i2.2,".dir")')            &
     &              species(ispecies)%nZ, species(jspecies)%nZ,               &
     &              species(kspecies)%nZ
                  open (unit = 13,                                            &
     &                  file = trim(Fdata_location)//trim(interactions),      &
     &                  status = 'unknown', position = 'append')
                  write (13,100) pFdata_bundle%nFdata_cell_3c, P_rhoS_3c,     &
     &                           isorp, itheta, filename(2:30),               &
     &                           pFdata_cell%nME, nna_rho, dna, nbc_rho, dbc
                  close (unit = 13)

                  ! Open mu, nu, mvalue file and write out values.
                  write (filename, '("/",i2.2, "_munu_3c.",                   &
     &                                   i2.2,".",i2.2,".",i2.2,".dat")')     &
     &               P_rhoS_3c, species(ispecies)%nZ, species(jspecies)%nZ,   &
     &                          species(kspecies)%nZ
                   open (unit = 12, file = trim(Fdata_location)//trim(filename),&
     &                   status = 'unknown', position = 'append')

                   ! Write out the mapping - stored in mu, nu, and mvalue
                   write (12,*) (pFdata_cell%mu_3c(index_3c),                 &
     &                                             index_3c = 1, nME3c_max)
                   write (12,*) (pFdata_cell%nu_3c(index_3c),                 &
     &                                             index_3c = 1, nME3c_max)
                   write (12,*) (pFdata_cell%mvalue_3c(index_3c),             &
     &                                                 index_3c = 1, nME3c_max)
                end do
              end do

              write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ,&
     &                             species(kspecies)%nZ

! Open all the output files.
              iounit = 12
              do isorp = ispmin, ispmax
                do itheta = 1, P_ntheta
                  iounit = iounit + 1
                  write (filename, '("/", "rhoS_3c_", i2.2, "_", i2.2, ".",   &
     &                               i2.2, ".", i2.2, ".", i2.2, ".dat")')    &
     &                  itheta, isorp, species(ispecies)%nZ,                  &
     &                  species(jspecies)%nZ, species(kspecies)%nZ
                  open (unit = (iounit),                                      &
     &                  file = trim(Fdata_location)//trim(filename),          &
     &                  status = 'unknown')
                end do
              end do

! ----------------------------------------------------------------------------
! Begin the big loops over dbc and dna.
! ----------------------------------------------------------------------------
! Loop over all bondcharge distances.
              do ibcba = 1, nbc_rho
                dbcx = float(ibcba - 1)*dbc/float(nbc_rho - 1)

! for all bondcharges-- we set b=dbcx/2.
                distance_bc = dbcx/2.0d0

! Evaluate overlapS for this distance of dbcx
                ! Set integration limits for overlapS
                zmin = max(-rcutoff1, dbcx - rcutoff2)
                zmax = min(rcutoff1, dbcx + rcutoff2)
                rhomin = 0.0d0
                rhomax = min(rcutoff1, rcutoff2)

! Set ideriv to zero here
                ideriv = 0
                call evaluate_integral_2c (index_2c_overlapS, ispecies,       &
     &                                     jspecies, isorp, ideriv, rcutoff1, &
     &                                     rcutoff2, dbcx, nz_rho, nrho_rho,  &
     &                                     rint_overlapS, nopi, zmin, zmax,   &
     &                                     rhomin, rhomax,                    &
     &                                     pFdata_overlapS_cell%fofx)

! Loop over all neutral atom distances.
! The distance is measured from the bondcharge center (b=dbcx/2)
                do inaba = 1, nna_rho
                  dnax = float(inaba - 1)*dna/float(nna_rho - 1)

! Loop over the shells
                  call evaluate_integral_3c (nFdata_cell_3c, ispecies,        &
     &                                       jspecies, kspecies, ispmin,      &
     &                                       ispmax, ctheta, ctheta_weights,  &
     &                                       dbcx, dnax, nnr_rho,             &
     &                                       nntheta_rho, psiSofr,            &
     &                                       phiint_rhoS, qpl)

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

! Now, divide by the spherical overlap so to get the weighted spherical density.
                      qpl(itheta,:,isorp) = qpl(itheta,:,isorp)               &
     &                                      /(pFdata_overlapS_cell%fofx + 1.0d-6)
                      iounit = iounit + 1

                      write (iounit,*)                                        &
     &                  (qpl(itheta,index_3c,isorp), index_3c = 1, nME3c_max)
                    end do
                  end do
                end do   ! end of the dna loop
              end do  ! the end of the dbc loop

              ! close files
              iounit = 12
              do isorp = ispmin, ispmax
                do itheta = 1, P_ntheta
                  iounit = iounit + 1
                  close (unit = iounit)
                end do
              end do
              deallocate (qpl)

            end do  ! end loop over kspecies
          end do  ! end loop over jspecies
        end do  ! end loop over ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i2, 1x, i2, 1x, i2, 1x, i3, 1x, a29, 1x, i3,              &
     &          1x, i4, 1x, f9.6, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating rhoS_3c integrals for nZ = ', i3,            &
     &              ' and nZ = ', i3, ', density on nZ = ', i3)


! End Subroutine
! ===========================================================================
        return
        end subroutine rhoS_3c_Harris


! ===========================================================================
! phiint_rhoS
! ===========================================================================
! Program Description
! ===========================================================================
! This program contains the actual final phi function for integration purposes.
! This is the function in center coordinates - r1 and r2 (the variable from
! the center's of the atoms), and rna, which is distance from the bond-center.
!
! The result is written into avgVmat which is are the integral stored
! in a mu, nu form - single dimension array form.
! ===========================================================================
        subroutine phiint_rhoS (itype, ispecies, jspecies, kspecies, ispmin,  &
                                   ispmax, r, ds, zr, r1, r2, rna, avgVmat)
        implicit none

        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies, jspecies, kspecies  ! the species
        integer, intent (in) :: ispmin, ispmax, itype    ! which type doing

        ! integration coordinates
        real, intent (in) :: r, ds, zr
        real, intent (in) :: r1, r2, rna(3)

! Output
        real, intent (out) :: avgVmat (ispmin:, :)

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer index_3c                ! counter for matrix location - mu, nu
        integer iphi                    ! integration over theta
        integer isorp                   ! loop over shells
        integer nME3c_max               ! number of matrix elements

        integer, allocatable :: mleft (:)  ! m quantum numbers
        integer, allocatable :: mright (:)

        real dummy
        real dphi                       ! interval between theta points
        real xr, yr, r3                 ! coordinate points
        real rmax                       ! maximum value of r along grid
        real phi                        ! value of phi at point of integration
        real prod
        real vpot                       ! value of potential

        real phifactor (-3:3)           ! phifactors for integration
        real, allocatable :: phimult (:)

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Allocate Arrays
! ===========================================================================
        allocate (phimult (nnphi_rho))

! Procedure
! ===========================================================================
        dummy = r1
        dummy = r2

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(itype)
        nME3c_max = pFdata_cell%nME

! allocate the quantum numbers and set values
        allocate (mleft (nME3c_max))
        allocate (mright (nME3c_max))
        mleft = pFdata_cell%M_mu(:)
        mright = pFdata_cell%M_nu(:)

! Initialize
        avgVmat = 0.0d0

! Set up integration factors
        dphi = pi/float(nnphi_rho - 1)
        phimult(1) = dphi*41.0d0/140.0d0
        do iphi = 2, nnphi_rho - 1
          if (mod(iphi,6) .eq. 2) phimult(iphi) = dphi*216.0d0/140.0d0
          if (mod(iphi,6) .eq. 3) phimult(iphi) = dphi*27.0d0/140.0d0
          if (mod(iphi,6) .eq. 4) phimult(iphi) = dphi*272.0d0/140.0d0
          if (mod(iphi,6) .eq. 5) phimult(iphi) = dphi*27.0d0/140.0d0
          if (mod(iphi,6) .eq. 0) phimult(iphi) = dphi*216.0d0/140.0d0
          if (mod(iphi,6) .eq. 1) phimult(iphi) = dphi*82.0d0/140.0d0
        end do
        phimult(nnphi_rho)= dphi*41.0d0/140.0d0

! ***************************************************************************
!              Do integral over phi:
! ***************************************************************************
! The phi factors depend only on m.
!
! Note: We order the p-orbitals here x,y,z (or pi,pi',sig), NOT z,x,y.
! Note that px, and xz now are +1. And so on!
        do iphi = 1, nnphi_rho
          phi = float(iphi - 1)*dphi

! set up the phifactors
          phifactor(0) = 1.0d0
          phifactor(1) = cos(phi)   ! p-states
          phifactor(-1) = sin(phi)
          phifactor(2) = cos(phi)**2 - sin(phi)**2
          phifactor(-2) = cos(phi)*sin(phi)

! Find coordinate x and y:
          xr = r*ds*cos(phi)
          yr = r*ds*sin(phi)
          r3 = sqrt((xr - rna(1))**2 + (yr - rna(2))**2 + (zr - rna(3))**2)

          do isorp = ispmin, ispmax
            rmax = species(kspecies)%shell(isorp)%rcutoffA
            vpot = (psiofr(r3, rmax, kspecies, isorp)**2)/(4.0d0*pi)

            prod = vpot*phimult(iphi)
            do index_3c = 1, nME3c_max
              avgVmat(isorp,index_3c) = avgVmat(isorp,index_3c)               &
        &       + prod*phifactor(mleft(index_3c))*phifactor(mright(index_3c))
            end do
          end do ! isorp
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (mleft, mright)
        deallocate (phimult)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine phiint_rhoS


! ===========================================================================
! psiSofr
! ===========================================================================
! Function Declaration
! ===========================================================================
        function psiSofr (r, rmax, ispecies, n) ! psiS change
        implicit none

        real psiSofr

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies, n

        real, intent (in) :: r
        real, intent (in) :: rmax

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ==========================================================================
! None

! Procedure
! ===========================================================================
        psiSofr = sqrt(psiofr(r, rmax, ispecies, n)**2)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function psiSofr

! End Module
! ===========================================================================
        end module
