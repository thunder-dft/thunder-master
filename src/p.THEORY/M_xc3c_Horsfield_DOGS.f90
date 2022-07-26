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
        module M_xc3c_DOGS

! /GLOBAL
        use M_precision

! /SYSTEM
        use M_species
        use M_atom_functions
        use M_integrals_3c

! /THEORY
        use M_vxc_Harris

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

! The charge transfer bit = need a +-dq on each orbital.
        real, allocatable :: dqorb (:)

! Some parameters for the derivative parts:
        integer, parameter, dimension (0:6) :: jsign = (/0, -1, +1, -1, +1, -1, +1/)
        integer, parameter, dimension (0:6) :: kspecies_key = (/1, 1, 1, 2, 2, 3, 3/)

! module procedures
        contains

! ===========================================================================
! initialize_xc3c_DOGS
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
        subroutine initialize_xc3c_DOGS
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
        integer itheta                       !< counter for theta

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max

        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
        ideriv_min = 1
        ideriv_max = 6

! Loop over species for
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

! For xc3c_DOGS; ideriv = ideriv_min, ideriv_max, so loop over ideriv
            do kspecies = 1, nspecies
              do ideriv = ideriv_min, ideriv_max
                pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
                do itheta = 1, P_ntheta
                  pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1
                end do
              end do ! ideriv
            end do ! kspecies
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
        end subroutine initialize_xc3c_DOGS


! ===========================================================================
! xc3c_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general three-center
! matrix elements of the form <psi1|V(1)|psi2> for the spherical density.
! ===========================================================================
        subroutine xc3c_DOGS
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
!       integer isuperloop              ! counter over species**3 - parallel
        integer ispecies, jspecies, kspecies  ! species numbers
!       integer itemp                   ! used to find species values
        integer nFdata_cell_3c          !< indexing of interactions

        ! different derivative cases
        integer ideriv, ideriv_min, ideriv_max
! MPI
!       integer my_proc, nproc
!       logical iammaster, iammpi

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
        write (ilogfile,*) '           E X C H A N G E   C O R R E L A T I O N       '
        write (ilogfile,*) '    T H R E E - C E N T E R    I N T E R A C T I O N S   '
        write (ilogfile,*) '                 (X C 3 C) M A T R I X                   '
        write (ilogfile,*) '                I N T E R A C T I O N S                  '
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*)

! Initialize the Legendre coefficients
        call gleg (ctheta, ctheta_weights, P_ntheta)

! Loop over the different derivative types. We already did ideriv=0 case.
! Set ideriv_min = 1 and ideriv_max = 2 for one center case.
        ideriv_min = 1
        ideriv_max = 6

! Initialize MPI
!       call initialize_MPI (iammaster, iammpi, my_proc, nproc)

! Loop over species along the bond charge - do in a super loop
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            do kspecies = 1, nspecies
!             do isuperloop = 1, nspecies*nspecies*nspecies

! Establish which species are on this processor
!               if (mod(isuperloop, nproc) .eq. my_proc) then
!               itemp = isuperloop
!               kspecies = 1 + int((itemp - 1)/(nspecies*nspecies))
!               itemp = itemp - (kspecies - 1)*(nspecies*nspecies)
!               jspecies = 1 + int((itemp - 1)/(nspecies))
!               itemp   = itemp - (jspecies - 1)*(nspecies)
!               ispecies = itemp

! Each ideriv evaluates +- dq changes in the density.
              do ideriv = ideriv_min, ideriv_max

                ! cut some lengthy notation with pointers
                pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
                pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1
                nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c
                pFdata_cell=>pFdata_bundle%Fdata_cell_3c(nFdata_cell_3c)

                call make_munu_3c (nFdata_cell_3c, ispecies, jspecies, kspecies)
                nME3c_max = pFdata_cell%nME

                allocate (qpl(P_ntheta, nME3c_max, ideriv_min:(ideriv_max - ideriv_min + 1)))
                qpl = 0.0d0

                ! Test output file for this species triplet
                itheta = 1
                write (filename, '("/", "xc3c_", i2.2, "_", i2.2, ".", i2.2,  &
     &                                         ".", i2.2, ".", i2.2, ".dat")')&
     &                 itheta, ideriv, species(ispecies)%nZ,                  &
     &                                 species(jspecies)%nZ, species(kspecies)%nZ
                inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
                if (skip) cycle

                ! Set up grid loop control constants
                rcutoff1 = species(ispecies)%rcutoffA_max
                rcutoff2 = species(jspecies)%rcutoffA_max
                rcutoff3 = species(kspecies)%rcutoffA_max
                dbc = rcutoff1 + rcutoff2
                dna = rcutoff3 + max(rcutoff1,rcutoff2)

                pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c - 1
                do itheta = 1, P_ntheta
                  pFdata_bundle%nFdata_cell_3c = pFdata_bundle%nFdata_cell_3c + 1

                  write (filename, '("/", "xc3c_", i2.2, "_", i2.2, ".", i2.2,&
     &                                        ".", i2.2, ".", i2.2, ".dat")') &
     &              itheta, ideriv, species(ispecies)%nZ, species(jspecies)%nZ,&
     &                              species(kspecies)%nZ

                  ! open directory file
                  write (interactions,                                        &
     &                   '("/3c.",i2.2,".",i2.2,".",i2.2,".dir")')            &
     &              species(ispecies)%nZ, species(jspecies)%nZ,               &
     &              species(kspecies)%nZ
                  open (unit = 13,                                            &
     &                  file = trim(Fdata_location)//trim(interactions),      &
     &                  status = 'unknown', position = 'append')
                  write (13,100) pFdata_bundle%nFdata_cell_3c, P_xc3c, ideriv,&
     &                           itheta, filename(2:30), pFdata_cell%nME,     &
     &                           nna_xc3c, dna, nbc_xc3c, dbc
                  close (unit = 13)

                  ! Open mu, nu, mvalue file and write out values.
                  write (filename, '("/",i2.2, "_munu_3c.",                   &
     &                                   i2.2,".",i2.2,".",i2.2,".dat")')     &
     &               P_xc3c, species(ispecies)%nZ, species(jspecies)%nZ,      &
     &                     species(kspecies)%nZ
                  open (unit = 12, file = trim(Fdata_location)//trim(filename),&
     &                  status = 'unknown', position = 'append')

                  ! Write out the mapping - stored in mu, nu, and mvalue
                  write (12,*) (pFdata_cell%mu_3c(index_3c),                  &
     &                                            index_3c = 1, nME3c_max)
                  write (12,*) (pFdata_cell%nu_3c(index_3c),                  &
     &                                                index_3c = 1, nME3c_max)
                  write (12,*) (pFdata_cell%mvalue_3c(index_3c),              &
     &                                                    index_3c = 1, nME3c_max)
                end do
              end do ! end loop over ideriv

              write (ilogfile,200) species(ispecies)%nZ, species(jspecies)%nZ,&
     &                             species(kspecies)%nZ

! Open all the output files.
              iounit = 12
              do ideriv = ideriv_min, ideriv_max
                do itheta = 1, P_ntheta
                  iounit = iounit + 1
                  write (filename, '("/", "xc3c_", i2.2, "_", i2.2, ".", i2.2,&
     &                               ".", i2.2, ".", i2.2, ".dat")')          &
     &                   itheta, ideriv, species(ispecies)%nZ,                &
     &                   species(jspecies)%nZ, species(kspecies)%nZ

! Write out the data...
                  open (unit = (iounit),                                      &
     &                  file = trim(Fdata_location)//trim(filename),          &
     &                  status = 'unknown')
                end do
              end do ! end loop over ideriv

! ----------------------------------------------------------------------------
! Begin the big loops over dbc and dna.
! ----------------------------------------------------------------------------
! Loop over all bondcharge distances.
             do ibcba = 1, nbc_xc3c
                dbcx = float(ibcba - 1)*dbc/float(nbc_xc3c - 1)

! for all bondcharges-- we set b=dbcx/2.
                distance_bc = dbcx/2.0d0

! Loop over all neutral atom distances.
! The distance is measured from the bondcharge center (b=dbcx/2)
                do inaba = 1, nna_xc3c
                  dnax = float(inaba - 1)*dna/float(nna_xc3c - 1)
                  call evaluate_integral_3c (nFdata_cell_3c, ispecies,        &
     &                                       jspecies, kspecies, ideriv_min,  &
     &                                       ideriv_max, ctheta,              &
     &                                       ctheta_weights, dbcx, dnax,      &
     &                                       nnr_xc3c, nntheta_xc3c, psiofr,  &
     &                                       phiint_xc3c, qpl)

! ----------------------------------------------------------------------------
! qpl's are the answer
! ------------------------------------- ---------------------------------------
! Write the qpl coefficients into the data files each combination in1, in2, in3,
! itheta(=1,ntheta_max), ideriv gives an individual file.  The values for the
! different non-zero matrix elements of a given combination are written out
! after the index loop.
! ----------------------------------------------------------------------------
                  iounit = 12
                  do ideriv = ideriv_min, ideriv_max
                    do itheta = 1, P_ntheta
                      iounit = iounit + 1
                      write (iounit,*)                                        &
     &                  (qpl(itheta,index_3c,ideriv), index_3c = 1, nME3c_max)
                    end do
                  end do  ! end loop over ideriv
                end do  ! end of the dna loop
              end do  ! the end of the dbc loop

              ! close files
              iounit = 12
              do ideriv = ideriv_min, ideriv_max
                do itheta = 1, P_ntheta
                  iounit = iounit + 1
                  close (unit = iounit)
                end do
              end do  ! end loop over ideriv

!         end if ! MPI which node end if
!       end do ! end loop over isuperloop

              deallocate (qpl)
            end do  ! end loop over kspecies
          end do  ! end loop over jspecies
        end do  ! end loop over ispecies

! Finalize MPI
!       call finalize_MPI

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i2, 1x, i2, 1x, i2, 1x, i3, 1x, a29, 1x, i3,             &
     &          1x, i4, 1x, f9.6, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating xc3c integrals for nZ = ', i3,              &
     &              ' and nZ = ', i3, ', potential on nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine xc3c_DOGS


! ===========================================================================
! phiint_xc3c
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
        subroutine phiint_xc3c (itype, ispecies, jspecies, kspecies,          &
     &                          ideriv_min, ideriv_max, r, ds, zr, r1, r2,    &
     &                          rna, avgVmat)
        implicit none

        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies, jspecies, kspecies    ! the species
        integer, intent (in) :: itype                      ! which type doing

        ! different derivative cases
        integer, intent (in) :: ideriv_min, ideriv_max

        ! integration coordinates
        real, intent (in) :: r, ds, zr
        real, intent (in) :: r1, r2, rna(3)

! Output
        real, intent (out) :: avgVmat(ideriv_min:, :)

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer index_3c                ! counter for matrix location - mu, nu
        integer iphi                    ! integration over theta
        integer ideriv, issh            ! loop over derivatives and shells
        integer nME3c_max               ! number of matrix elements

        integer, allocatable :: mleft (:)  ! m quantum numbers
        integer, allocatable :: mright (:)

        real dphi                       ! interval between theta points
        real xr, yr, r3                 ! coordinate points
        real phi                        ! value of phi at point of integration
        real prod
        real vpot                       ! value of potential

        real dqint                      ! charge transfer bit for +-dq

        real phifactor (-3:3)           ! phifactors for integration
        real, allocatable :: phimult (:)

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Allocate Arrays
! ===========================================================================
        allocate (phimult (nnphi_xc3c))

! Procedure
! ===========================================================================
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
        dphi = pi/float(nnphi_xc3c - 1)
        phimult(1) = dphi*41.0d0/140.0d0
        do iphi = 2, nnphi_xc3c - 1
          if (mod(iphi,6) .eq. 2) phimult(iphi) = dphi*216.0d0/140.0d0
          if (mod(iphi,6) .eq. 3) phimult(iphi) = dphi*27.0d0/140.0d0
          if (mod(iphi,6) .eq. 4) phimult(iphi) = dphi*272.0d0/140.0d0
          if (mod(iphi,6) .eq. 5) phimult(iphi) = dphi*27.0d0/140.0d0
          if (mod(iphi,6) .eq. 0) phimult(iphi) = dphi*216.0d0/140.0d0
          if (mod(iphi,6) .eq. 1) phimult(iphi) = dphi*82.0d0/140.0d0
        end do
        phimult(nnphi_xc3c)= dphi*41.0d0/140.0d0

! ***************************************************************************
!              Do integral over phi:
! ***************************************************************************
! The phi factors depend only on m.
!
! Note: We order the p-orbitals here x,y,z (or pi,pi',sig), NOT z,x,y.
! Note that px, and xz now are +1. And so on!
        do iphi = 1, nnphi_xc3c
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

          do ideriv = ideriv_min, ideriv_max

            ! Set the values of Qneutral_ion to the original Qneutral
            do issh = 1, species(ispecies)%nssh
              species(ispecies)%shell(issh)%Qneutral_ion =                    &
     &          species(ispecies)%shell(issh)%Qneutral
            end do
            do issh = 1, species(jspecies)%nssh
              species(jspecies)%shell(issh)%Qneutral_ion =                    &
     &          species(jspecies)%shell(issh)%Qneutral
            end do
            do issh = 1, species(kspecies)%nssh
              species(kspecies)%shell(issh)%Qneutral_ion =                    &
     &          species(kspecies)%shell(issh)%Qneutral
            end do

            ! Now change the charge state to affect the density
            if (kspecies_key(ideriv) .eq. 1) then
              dqint = dqorb(ispecies)/species(ispecies)%nssh
              do issh = 1, species(ispecies)%nssh
                species(ispecies)%shell(issh)%Qneutral_ion =                  &
     &            species(ispecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
              end do
            else if (kspecies_key(ideriv) .eq. 2) then
              dqint = dqorb(jspecies)/species(jspecies)%nssh
              do issh = 1, species(jspecies)%nssh
                species(jspecies)%shell(issh)%Qneutral_ion =                  &
     &            species(jspecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
              end do
            else
              dqint = dqorb(kspecies)/species(kspecies)%nssh
              do issh = 1, species(kspecies)%nssh
                species(kspecies)%shell(issh)%Qneutral_ion =                  &
     &            species(kspecies)%shell(issh)%Qneutral + jsign(ideriv)*dqint
              end do
            end if

            vpot = dvxc_3c (itype, ispecies, jspecies, kspecies, r1, r2, r3)

            prod = vpot*phimult(iphi)
            do index_3c = 1, nME3c_max
              avgVmat(ideriv,index_3c) = avgVmat(ideriv,index_3c)             &
        &       + prod*phifactor(mleft(index_3c))*phifactor(mright(index_3c))
            end do
          end do ! ideriv
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
        end subroutine phiint_xc3c


! ===========================================================================
! dvxc_3c
! ===========================================================================
! Program Description
! ===========================================================================
!       This subroutine computes vxc(n1+n2+n3) - vxc(n1+n2) with the
! densities ni of atom in_i at the distance ri from their centers.  Due to
! the small contribution of the three center case to the overall energy,
! only the lda level of theory will be used for this calculation.
!
! On output:
!     dvxc3c = vxc(n1+n2+n3) - vxc(n1+n2)
!
! ===========================================================================
! Original code from Juergen Fritsch

! Code rewritten by:
! Richard B. Evans
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
        real function dvxc_3c (itype, ispecies, jspecies, kspecies, r1, r2, r3)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itype                        ! which type doing
        integer, intent (in) :: ispecies, jspecies, kspecies ! the species

        real, intent (in) :: r1, r2, r3            ! locations of three centers

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer iexc                    !< which flavor of exchange-correlation

        real drho                       !< distance between mesh points
        real rin                        !< value of r in Bohr radii

! Value of density and corresponding derivatives at the point r, z
! ....for one center piece
        real density, density_p, density_pp

! Value of density and corresponding derivatives at the point r, z
! Exchange-correlation potential and energies for two-center density
        real density_2c
        real density_2c_p, density_2c_pp
        real dnuxc_2c, dnuxcs_2c
        real exc_2c, vxc_2c, dexc_2c

! Value of density and corresponding derivatives at the point r, z
! Exchange-correlation potential and energies for three-center density
        real density_3c
        real density_3c_p, density_3c_pp
        real dnuxc_3c, dnuxcs_3c
        real exc_3c, vxc_3c, dexc_3c

        real xc_fraction                ! fraction of exact exchange

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
        idummy = itype

! This exchange-correlation routine deals only with three-center interactions,
! therefore, we do not do exact exchange for the three-center terms and
! fraction should be initiallized to 1.0d0.
        xc_fraction = 1.0d0

! We always calculate LDA (iexc = 3) for the three-center interactions because
! doing gradient corrections for three-centers is a bug-a-boo and does not
! yield significant improvements.
        iexc = 3

! Establish drho for this one-center case.
        drho = min(species(ispecies)%rcutoffA_max,                            &
                   species(jspecies)%rcutoffA_max, species(kspecies)%rcutoffA_max)
        drho = drho/dfloat(nrho_rho_store)

! Find the total density of all three centers - located at r, r2, and r3
! Evaluate the density of each center individually and then sum to get
! the total density.

! One-center piece: vxc[n1(r)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call rho_1c (ispecies, r1, drho, density, density_p, density_pp)
        density_2c = density
        density_2c_p = density_p
        density_2c_pp = density_pp

! One-center piece: vxc[n2(r2)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call rho_1c (jspecies, r2, drho, density, density_p, density_pp)
        density_2c = density_2c + density
        density_2c_p = density_2c_p + density_p
        density_2c_pp = density_2c_pp + density_pp

! One-center piece: vxc[n3(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call rho_1c (kspecies, r3, drho, density, density_p, density_pp)
        density_3c = density_2c + density
        density_3c_p = density_2c_p + density_p
        density_3c_pp = density_2c_pp + density_pp

! Three-center-piece: vxc_3c[n1(r1) + n2(r2) + n3(r3)]
! Compute the exchange correlation potential for the three-center case
! ***************************************************************************
! The total three-center density is the sum of the three.
! Note that rin, densp, and denspp are not used in the LDA limits.
        rin = r1/P_abohr
        density_3c = density_3c*P_abohr**3
        density_3c_p = density_3c_p*P_abohr**4
        density_3c_pp = density_3c_pp*P_abohr**5
        call get_potxc_1c (iexc, xc_fraction, rin, density_3c, density_3c_p,  &
     &                     density_3c_pp, exc_3c, vxc_3c, dnuxc_3c, dnuxcs_3c,&
     &                     dexc_3c)

! Two-center-piece: vxc_2c[n1(r1) + n2(r2)]
! Compute the exchange correlation potential for the three-center case
! ***************************************************************************
! The total two-center density is the sum of the two - dens1 + dens2.
! Note that rin, densp, and denspp are not used in the LDA limits.
        rin = r1/P_abohr
        density_2c = density_2c*P_abohr**3
        density_2c_p = density_2c_p*P_abohr**4
        density_2c_pp = density_2c_pp*P_abohr**5
        call get_potxc_1c (iexc, xc_fraction, rin, density_2c, density_2c_p,  &
     &                     density_2c_pp, exc_2c, vxc_2c, dnuxc_2c, dnuxcs_2c,&
     &                     dexc_2c)

! Answers are in Hartrees convert to eV.
        dvxc_3c = P_Hartree*(vxc_3c - vxc_2c)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end function dvxc_3c

! End Module
! =============================================================================
        end module
