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
        module M_vxc_Harris

! /GLOBAL
        use M_precision

! /SYSTEM
        use M_atom_functions
        use M_atomPP_functions
        use M_atomPP_ion_functions
        use M_species
        use M_integrals_2c

! /CREATE
        use M_xc_1c
        use M_xc_2c

! /MPI
        use M_MPI

! Type Declaration
! ===========================================================================
! This type contains the two-center density and derivatives on rho, z grid
        type T_rho_2c_store
          integer nnrho                 ! number of grid points along rho
          integer nnz                   ! number of grid points along z

          real drho                     ! spacing between rho points
          real dz                       ! spacing between z points
          real rhomin                   ! minimum rho for that density
          real rhomax                   ! maximum rho for that density
          real zmin                     ! minimum z for that density
          real zmax                     ! maximum z for that density

          real, pointer :: rho (:, :)   ! the density - function of r and z
          real, pointer :: rhop (:, :)  ! derivative with respect to rho
          real, pointer :: rhopp (:, :) ! second derivative with respect to rho
          real, pointer :: rhoz (:, :)  ! derivative with respect to z
          real, pointer :: rhozz (:, :) ! second derivative with respect to z
          real, pointer :: rhopz (:, :) ! mixed derivative - rho and z
        end type T_rho_2c_store

! the grid that contains the two-center densities
        type T_rho_2c_bundle

          ! the density on grid of d (distance between two centers)
          type (T_rho_2c_store), pointer :: rho_2c_store (:)
        end type T_rho_2c_bundle

! stored two-center densities into bundles based on ispecies and jspecies pair
        type(T_rho_2c_bundle), pointer :: rho_2c_bundle (:, :)

! module procedures
        contains

! ===========================================================================
! initialize_vxc_Harris
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
        subroutine initialize_vxc_Harris
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

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species
        do ispecies = 1, nspecies
          pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)

          ! add one for vxc_1c
          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1

          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

            ! add one for uxc
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1

            ! add one for vxc_ontop
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1

            ! add one for vxc_atom
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
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
        end subroutine initialize_vxc_Harris


! ===========================================================================
! vxc_Harris
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calls the subroutines required to calculate the vna_ontop
! (both left/right cases) and atom cases.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine vxc_Harris
        implicit none

! Parameters and Data Declaration
! ===========================================================================
! None

! ===========================================================================
! Argument Declaration and Description
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer issh                        !< counters for number of shells

! Procedure
! ===========================================================================
! begin iammaster
        if (my_proc .eq. 0) then
          write (ilogfile,*)
          write (ilogfile,*) ' ******************************************************* '
          write (ilogfile,*) '        E X C H A N G E   C O R R E L A T I O N          '
          write (ilogfile,*) '                  I N T E R A C T I O N S                '
          write (ilogfile,*) ' ******************************************************* '
          write (ilogfile,*)
        end if

! We are doing Harris interactions here, so we set Qneutral_ion to the
! non-charged atom before calculating the densities.  The Qneutral is used
! throughout, so we set it here for convenience.
! Loop over species
        do ispecies = 1, nspecies
          do issh = 1, species(ispecies)%nssh
            species(ispecies)%shell(issh)%Qneutral_ion =                      &
     &        species(ispecies)%shell(issh)%Qneutral
         end do
        end do

        if (my_proc .eq. 0) then
          write (ilogfile,*) ' Calling one-center case. '
          call vxc_1c_Harris
        end if

! begin iammaster
        if (my_proc .eq. 0) write (ilogfile,*)
        if (my_proc .eq. 0) then
          write (ilogfile,*) ' Building the two center density on grid '
        end if

        allocate (rho_2c_bundle (nspecies, nspecies))
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            call rho_2c_store (ispecies, jspecies)
          end do
        end do

        if (my_proc .eq. 0) write (ilogfile,*)
        if (my_proc .eq. 0) write (ilogfile,*) ' Calling correction case. '
        call uxc_Harris

        if (my_proc .eq. 0) write (ilogfile,*)
        if (my_proc .eq. 0) write (ilogfile,*) ' Calling ontop case. '
        call vxc_ontop_Harris

        if (my_proc .eq. 0) write (ilogfile,*)
        if (my_proc .eq. 0) write (ilogfile,*) ' Calling atom case. '
        call vxc_atom_Harris

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ==========================================================================
        return
        end subroutine vxc_Harris


! ===========================================================================
! rho_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This routine calculates and stores the combined density of a single species
! as a function of r, z and d. This is important because later the gradients
! with respect to these variables are calculated.
!
! On output: rho:      The density as a function of species type,
!                      r, z and d.  Output is placed in common block
!                      density located in wavefunctions.inc
!            rhop:     derivative with respect to rho
!            rhopp:    second derivative with respect to rho
!
! ====================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ====================================================================
!
! Program Declaration
! ====================================================================
        subroutine rho_1c (ispecies, r, dr, rho, rhop, rhopp)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispecies  !< the species number for the density

        real, intent (in) :: r            !< the radial coordinate
        real, intent (in) :: dr           !< the grid spacing

! Output
        ! density value and derivatives
        real, intent (out) :: rho, rhop, rhopp

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer issh                        ! index for looping over shells
        integer iexc                        ! which type of exchange-correlation

        real density                        ! the one-center density
        real density_pdr                    ! one-center density at r + dr
        real density_mdr                    ! one-center density at r - dr
        real density_p2dr                   ! one-center density at r + 2dr
        real xnocc                          ! the occupation charge
        real rin                            ! r + 2*dr for call to psiofr
        real rmax                           ! maximum value of r along grid

! Local Parameters and Data Declaration
! ====================================================================
! None

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Set iexc
        iexc = species_PP(ispecies)%iexc

! Here the density is computed for r, r+dr and r-dr
        density = 0.0d0
        do issh = 1, species(ispecies)%nssh
          xnocc = species(ispecies)%shell(issh)%Qneutral_ion
          rmax = species(ispecies)%shell(issh)%rcutoffA
          density = density + xnocc*psiofr(r, rmax, ispecies, issh)**2
        end do
        rho = density/(4.0d0*4.0d0*atan(1.0d0))

! **************************************************************************
! Now calculate the derivatives
! Only calculate the derivatives if doing GGA exchange-correlation.
        rhop = 0.0d0
        rhopp = 0.0d0
!       rhopap = 0.0d0
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 .or.                &
     &      iexc .eq. 9 .or. iexc .eq. 10) then

          density_pdr = 0.0d0
          density_mdr = 0.0d0
          do issh = 1, species(ispecies)%nssh
            xnocc = species(ispecies)%shell(issh)%Qneutral_ion
            rmax = species(ispecies)%shell(issh)%rcutoffA
            density_pdr = density_pdr + xnocc*psiofr(r + dr,rmax,ispecies,issh)**2
            density_mdr = density_mdr + xnocc*psiofr(r - dr,rmax,ispecies,issh)**2
          end do

! Here the first and second derivatives of the density is computed.
          if ((r - dr) .gt. 1.0d-5) then
            rhop = (density_pdr - density_mdr)/(2.0d0*dr)
            rhopp = (density_pdr - 2.0d0*density + density_mdr)/dr**2
!           rhopap = (abs(density_pdr - density)                             &
!    &                - abs(density - density_mdr))/dr**2
          else

! At the endpoint do a forward difference. First, we need the point at r+2dr.
            density_p2dr = 0.0d0
            do issh = 1, species(ispecies)%nssh
              xnocc = species(ispecies)%shell(issh)%Qneutral_ion
              rmax = species(ispecies)%shell(issh)%rcutoffA
              rin = r + 2.0d0*dr
              density_p2dr =                                                 &
     &          density_p2dr + xnocc*psiofr(rin, rmax, ispecies, issh)**2
            end do

            rhop = (density_pdr - density)/dr
            rhopp = (density_p2dr - 2.0d0*density_pdr                        &
     &                            + density)/dr**2
!           rhopap = (abs(density - density_p2dr)                            &
!                     - abs(density_pdr - density_p2dr))/2.0*dr**2
          end if
        end if
        rhop = rhop/(4.0d0*4.0d0*atan(1.0d0))
        rhopp = rhopp/(4.0d0*4.0d0*atan(1.0d0))
!       rhopap = rhopap/(4.0d0*4.0d0*atan(1.0d0))

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ============================================================================
! None

        return
        end subroutine rho_1c


! ===========================================================================
! vxc_1c_Harris
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
        subroutine vxc_1c_Harris
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
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_1c              !< indexing of interactions

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
! We are doing only Harris here, so set ideriv and isorp = 0
        ideriv = 0
        isorp = 0

! Loop over species
        do ispecies = 1, nspecies
          pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
          nFdata_cell_1c = pFdata_bundle%nFdata_cell_2c
          pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_1c)

          call make_munuS (nFdata_cell_1c, ispecies, ispecies)
          nME1c_max = pFdata_cell%nME
          allocate (pFdata_cell%fofx(nME1c_max))

          ! Open ouput file for this species pair
          write (filename, '("/vxc_1c", ".", i2.2, ".dat")')                  &
     &      species(ispecies)%nZ
          open (unit = 11, file = trim(Fdata_location)//trim(filename),       &
     &          status = 'unknown')

          ! Set integration limits
          rhomin = 0.0d0
          rhomax = species(ispecies)%rcutoffA_max

! Loop over grid
          if (my_proc .eq. 0) write (ilogfile,100) species(ispecies)%nZ

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
            write (11,*) (pFdata_cell%fofx(jssh),                             &
     &                    jssh = index_1c, index_1c + species(ispecies)%nssh - 1)
            index_1c = index_1c + species(ispecies)%nssh
          end do
          write (11,*)

          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c - 1

          ! Deallocate rhomult so we can recalculate for the next species
          deallocate (rhomult)
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Evaluating vxc_1c integrals for nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vxc_1c_Harris


! ===========================================================================
! dexc_1c
! ===========================================================================
! Program Description
! ===========================================================================
! This function computes only exc(n1) with the densities ni of atom in_i.
!
! On input:  in1: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  exc: exc[n1(r1)]  (the one center exc)

! We calculate exc(n1).
!
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
! Program Declaration
! ===========================================================================
        real function dexc_1c (iexc, xc_fraction, ispecies, r)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc      !< which exchange-correlation
        integer, intent (in) :: ispecies  !< two centers species

        real, intent (in) :: xc_fraction  !< fraction of exact exchange
        real, intent (in) :: r            !< integration coordinate

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real drho                           !< distance between mesh points
        real rin                            !< value of r in Bohr radii

! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp

! Output from calling get_potxc_1c
        real exc
        real vxc
        real dnuxc
        real dnuxcs
        real dexc

! Procedure
! ===========================================================================
! Establish drho for this one-center case.
        drho = species(ispecies)%rcutoffA_max/float(nrho_rho_store)

! One-center piece: exc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call rho_1c (ispecies, r, drho, density, density_p, density_pp)

        rin = r/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5

        exc = 0.0d0
        vxc = 0.0d0
        call get_potxc_1c (iexc, xc_fraction, rin, density, density_p,        &
     &                     density_pp, exc, vxc, dnuxc, dnuxcs, dexc)

! Answers are in Hartrees convert to eV.
        dexc_1c = P_hartree*(density/P_abohr**3)*(exc - vxc)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end function dexc_1c


! ===========================================================================
! rint_vxc_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This is the radial integral function for the one-center exchange-
! correlation interactions.
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
        real function rint_vxc_1c (itype, ispecies, rho, index_2c)
        implicit none

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies               ! two centers species
        integer, intent (in) :: itype                  ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent(in) :: rho

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iexc                     ! which exchange-correlation option
        integer l1, n1                   ! quantum  numbers
        integer l2, n2

        real psi1val, psi2val
        real rmax                        ! maximum value of r along grid

        real vofr
        real xc_fraction

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Find which exchange-correlation we are calcuating:
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Pick up "passed" data
        l1 =  pFdata_cell%L_mu(index_2c)
        n1 =  pFdata_cell%N_mu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        n2 =  pFdata_cell%N_nu(index_2c)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n1)%rcutoffA
        psi1val = psiofr (rho, rmax, ispecies, n1)

! Find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n2)%rcutoffA
        psi2val = psiofr (rho, rmax, ispecies, n2)

! One-center piece: exc[n1(r1)]
        if (l1 .eq. l2) then
          vofr = dvxc_1c (iexc, xc_fraction, ispecies, rho)
        else
          vofr = 0.0d0
        end if

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_vxc_1c = psi1val*vofr*psi2val*rho**2

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_vxc_1c


! ===========================================================================
! dvxc_1c
! ===========================================================================
! Program Description
! ===========================================================================
! This function computes only vxc(n1) with the densities ni of atom in_i.
!
! On input:  in1: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  vxc: vxc[n1(r1)]  (the one center vxc)

! We calculate vxc(n1) here.
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
! Program Declaration
! ===========================================================================
        real function dvxc_1c (iexc, xc_fraction, ispecies, r1)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc      !< which exchange-correlation
        integer, intent (in) :: ispecies  !< two centers species

        real, intent (in) :: xc_fraction  !< fraction of exact exchange
        real, intent (in) :: r1           !< location of the atom center

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real drho                           !< distance between mesh points
        real rin                            !< value of r in Bohr radii

! Value of density and corresponding derivatives at the point r, z
! Exchange-correlation potential and energies for three-center density
        real density
        real density_p, density_pp

! Output from calling get_potxc_1c
        real dnuxc, dnuxcs
        real exc, vxc, dexc

! Procedure
! ===========================================================================
! Establish drho for this one-center case.
        drho = species(ispecies)%rcutoffA_max/dfloat(nrho_rho_store)

! One-center piece: vxc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call rho_1c (ispecies, r1, drho, density, density_p, density_pp)

        rin = r1/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5

        vxc = 0.0d0
        call get_potxc_1c (iexc, xc_fraction, rin, density, density_p,        &
     &                     density_pp, exc, vxc, dnuxc, dnuxcs, dexc)

! Answers are in Hartrees convert to eV.
        dvxc_1c = P_hartree*vxc

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end function dvxc_1c


! ===========================================================================
! rho_2c_store
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This routine calculates and stores the combined density of two species
! as a function of r, z and d. This is important because later the gradients
! with respect to these variables are calculated.
!
! On output: rho:      The density as a function of species type,
!                      r, z and d.  Output is placed in common block
!                      density located in wavefunctions.inc
!            rhop:     derivative with respect to rho
!            rhopp:    second derivative with respect to rho
!            rhoz:     derivative with respect to z
!            rhozz:    second derivative with respect to rho
!            rhopz:    mixed derivative with respect to rho and z
!
! ====================================================================
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
! Program Declaration
! ====================================================================
        subroutine rho_2c_store (ispecies, jspecies)
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies          !< which pair

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer issh, jssh                  ! index for looping over shells
        integer iexc                        ! which type of exchange-correlation
        integer igrid                       !< counter for grid points
        integer iz, irho                    !< counter for grid points
        integer nnz, nnrho                  !< number of intergration points

        real d                              !< distance between the two centers
        real dmax                           !< max distance between two centers
        real drr, dz, drho                  !< distance between mesh points
        real rcutoff1, rcutoff2             !< cutoffs for the two centers
        real r1, r2, z1, z2                 !< value of r and z for two centers
        real xnocc                          !< the occupation charge

        real density                        !< the two-center density

        real zmin, zmax
        real rho, rhomin, rhomax

        real rmax                           ! maximum value of r along grid

        type (T_rho_2c_bundle), pointer :: prho_bundle
        type (T_rho_2c_store), pointer :: prho_2c

! Local Parameters and Data Declaration
! ====================================================================
! None

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Set iexc
        iexc = species_PP(ispecies)%iexc

! Loop over species
        prho_bundle=>rho_2c_bundle(ispecies, jspecies)

        ! Set up grid loop control constants
        rcutoff1 = species(ispecies)%rcutoffA_max
        rcutoff2 = species(jspecies)%rcutoffA_max

        dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
        drr = dmax/float(ndd_vxc - 1)
        d = -drr

        ! Set integration limits
        rhomin = 0.0d0
        rhomax = max(rcutoff1, rcutoff2)

! Loop over grid
!       write (ilogfile,100) species(ispecies)%nZ, species(jspecies)%nZ
        allocate (prho_bundle%rho_2c_store(ndd_vxc))
        do igrid = 1, ndd_vxc
          prho_2c => prho_bundle%rho_2c_store(igrid)
          d = d + drr

          ! Set integration limits
          zmin = min(-rcutoff1, d - rcutoff2)
          zmax = max(rcutoff1, d + rcutoff2)

          dz = ((rcutoff1 + rcutoff2)/2.0d0)/float(nz_rho_store)
          nnz = int((zmax - zmin)/dz)
          if (mod(nnz, 2) .eq. 0) nnz = nnz + 1

          drho = ((rcutoff1 + rcutoff2)/2.0d0)/float(nrho_rho_store)
          nnrho = int((rhomax - rhomin)/drho)
          if (mod(nnrho,2) .eq. 0) nnrho = nnrho + 1

! Allocate all of the rho_2c arrays and assign variables
          allocate (prho_2c%rho (nnrho, nnz))

          prho_2c%nnrho = nnrho
          prho_2c%drho = drho

          prho_2c%nnz = nnz
          prho_2c%dz = dz

          prho_2c%rhomin = rhomin
          prho_2c%rhomax = rhomax

          prho_2c%zmin = zmin
          prho_2c%zmax = zmax

! **************************************************************************
! Here we loop over z and r computing the sum of the densities
! for species in1 and in2 at each value of d, z and r.
          do iz = 1, nnz
            z1 = zmin + (iz - 1)*dz
            z2 = z1 - d

            do irho = 1, nnrho
              rho = rhomin + (irho - 1)*drho
              r1 = sqrt(z1**2 + rho**2)
              r2 = sqrt(z2**2 + rho**2)

! Evaluate the density at this grid point.
              density = 0.0d0
              do issh = 1, species(ispecies)%nssh
                xnocc = species(ispecies)%shell(issh)%Qneutral_ion
                rmax = species(ispecies)%shell(issh)%rcutoffA
                density = density + xnocc*psiofr(r1, rmax, ispecies, issh)**2
              end do

              do jssh = 1, species(jspecies)%nssh
                xnocc = species(jspecies)%shell(jssh)%Qneutral_ion
                rmax = species(jspecies)%shell(jssh)%rcutoffA
                density = density + xnocc*psiofr(r2, rmax, jspecies, jssh)**2
              end do
              prho_2c%rho(irho, iz) = density/(4.0d0*4.0d0*atan(1.0d0))
            end do
          end do

! **************************************************************************
! Now calculate the derivatives
! Only calculate the derivatives if doing GGA exchange-correlation.
          if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 .or.              &
     &        iexc .eq. 9 .or. iexc .eq. 10) then
            allocate (prho_2c%rhop (nnrho, nnz))
            allocate (prho_2c%rhopp (nnrho, nnz))
            allocate (prho_2c%rhoz (nnrho, nnz))
            allocate (prho_2c%rhozz (nnrho, nnz))
            allocate (prho_2c%rhopz (nnrho, nnz))

! Calculate rhop_2c and rhopp_2c.
            do iz = 1, nnz
              do irho = 2, nnrho - 1

! First derivative:
                prho_2c%rhop(irho, iz) =                                      &
     &            (prho_2c%rho(irho + 1, iz) - prho_2c%rho(irho - 1, iz))/(2.0d0*drho)

! Second derivative:
                prho_2c%rhopp(irho, iz) =                                     &
     &            (prho_2c%rho(irho + 1, iz) - 2.0d0*prho_2c%rho(irho, iz)    &
     &             + prho_2c%rho(irho - 1, iz))/(drho**2)
              end do

! Find endpoint values for the derivatives calculated above.
              prho_2c%rhop(1, iz) =                                           &
     &          2.0d0*prho_2c%rhop(2, iz) - prho_2c%rhop(3, iz)

              prho_2c%rhop(nnrho, iz) =                                       &
     &          2.0d0*prho_2c%rhop(nnrho - 1, iz) - prho_2c%rhop(nnrho - 2, iz)

              prho_2c%rhopp(1, iz) =                                          &
     &          2.0d0*prho_2c%rhopp(2, iz) - prho_2c%rhopp(3, iz)

              prho_2c%rhopp(nnrho, iz) =                                      &
     &          2.0d0*prho_2c%rhopp(nnrho - 1, iz) - prho_2c%rhopp(nnrho - 2, iz)
            end do

! Calculate rhoz_2c and rhozz_2c.
            do irho = 1, nnrho
              do iz = 2, nnz - 1

! First derivative:
                prho_2c%rhoz(irho, iz) =                                      &
     &            (prho_2c%rho(irho, iz + 1) - prho_2c%rho(irho, iz - 1))/(2.0d0*dz)

! Second derivative:
                prho_2c%rhozz(irho, iz) =                                     &
     &            (prho_2c%rho(irho, iz + 1) - 2.0d0*prho_2c%rho(irho, iz)    &
     &             + prho_2c%rho(irho, iz - 1))/(dz**2)
              end do

! Find endpoint values for the derivatives calculated above.
              prho_2c%rhoz(irho, 1) =                                         &
     &          2.0d0*prho_2c%rhoz(irho, 2) - prho_2c%rhoz(irho, 3)

              prho_2c%rhoz(irho, nnz) =                                       &
     &          2.0d0*prho_2c%rhoz(irho, nnz - 1) - prho_2c%rhoz(irho, nnz - 2)

              prho_2c%rhozz(irho, 1) =                                        &
     &          2.0d0*prho_2c%rhozz(irho, 2) - prho_2c%rhozz(irho, 3)

              prho_2c%rhozz(irho, nnz) =                                      &
     &          2.0d0*prho_2c%rhozz(irho, nnz - 1) - prho_2c%rhozz(irho, nnz - 2)
            end do

! **************************************************************************
! Now calculate the cross term derivatives - rhopz2c.
            do irho = 1, nnrho
              do iz = 2, nnz - 1
                prho_2c%rhopz(irho, iz) =                                     &
     &            (prho_2c%rhop(irho, iz + 1) - prho_2c%rhop(irho, iz - 1))/(2.0d0*dz)
              end do

! Now calculate the derivatives at the endpoints.
              prho_2c%rhopz(irho, 1) =                                        &
     &          2.0d0*prho_2c%rhopz(irho,2) - prho_2c%rhopz(irho, 3)

              prho_2c%rhopz(irho, nnz) =                                      &
     &          2.0d0*prho_2c%rhopz(irho, nnz - 1) - prho_2c%rhopz(irho, nnz - 2)
            end do
          end if
        end do ! end loop over the grid (distance between centers)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ============================================================================
100     format (2x, ' Evaluating rho_2c_store arrays for nZ = ', i3,         &
     &              ' and nZ = ', i3)

        return
        end subroutine rho_2c_store


! ===========================================================================
! uxc_Harris
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
        subroutine uxc_Harris
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
! We are doing only Harris here, so set ideriv = 0 and isorp = 0
        ideriv = 0
        isorp = 0

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

            ! This correction term does not have matrix elements;
            ! set nME2c_max = 1
            call make_munu_uxc (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

! begin iammaster
            if (my_proc .eq. 0)                                              &
     &        write (ilogfile,200) ideriv, species(ispecies)%nZ, species(jspecies)%nZ

            ! Open ouput file for this species pair
            write (filename, '("/uxc_", i2.2,".",i2.2,".",i2.2,".dat")')     &
     &             ideriv, species(ispecies)%nZ, species(jspecies)%nZ
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown')

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + na(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_uxc - 1)
            d = -drr

            ! Set final integration limit
            rhomax = min(rcutoff1, rcutoff2)

            ! open directory file
            write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')             &
     &        species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &            status = 'unknown', position = 'append')
            write (13,100) pFdata_bundle%nFdata_cell_2c, P_uxc, ideriv,      &
     &                     filename(2:30), pFdata_cell%nME, ndd_uxc, dmax
            close (unit = 13)

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")') &
     &             P_uxc, species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                    index_2c = 1, nME2c_max)
            close (unit = 12)

! Loop over grid
            do igrid = 1, ndd_uxc
              d = d + drr

              ! Set integration limits
              zmin = max(-rcutoff1, d - rcutoff2)
              zmax = min(rcutoff1, d + rcutoff2)

              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, nz_uxc, nrho_uxc, rint_uxc,      &
     &                                   phifactor, zmin, zmax,              &
     &                                   rhomin, rhomax, pFdata_cell%fofx)

              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                      &
     &                                       index_2c = 1, nME2c_max)
            end do ! igrid
            write (11,*)
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating uxc integrals for ideriv = ', i3,           &
     &              ' nZ = ', i3, ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine uxc_Harris


! ===========================================================================
! rint_uxc
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
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
! Program Declaration
! ===========================================================================
        real function rint_uxc (itype, ispecies, jspecies, isorp, d,         &
     &                          rho, z1, z2, ideriv, index_2c)
        implicit none

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent (in) ::  d           !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer iexc                     ! which exchange-correlation option
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1, r2                      ! location of the atom centers

        real vofr
        real xc_fraction

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
        idummy = isorp
        dummy = d

! Find which exchange-correlation we are calcuating:
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

!Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)
        r2 = sqrt(z2**2 + rho**2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
! No wavefunctions used here.
        psi1val = 1.0d0
        psi2val = 1.0d0
        vofr = dexc_2c (iexc, xc_fraction, ispecies, jspecies, rho, z1, d, r1, r2)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_uxc = psi1val*vofr*psi2val*rho*(4.0d0*pi)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_uxc


! ===========================================================================
! dexc_2c
! ===========================================================================
! Program Description
! ===========================================================================
! This function computes the exchange-correlation correction case of
! exc(n1 + n2) - vxc(n1 + n2).
!
! On input:   r, z: geometry information for the charge gradient
!
! On output:  dexc_2c: exc[n1(r1) + n2(r2)] - vxc[n1(r1) + n2(r2)]
!
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
! Program Declaration
! ===========================================================================
        real function dexc_2c (iexc, xc_fraction, ispecies, jspecies,        &
     &                         r, z, d, r1, r2)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc      !< which exchange-correlation
        integer, intent (in) :: ispecies, jspecies      !< two centers species

        real, intent (in) :: d            !< distance between the two centers
        real, intent (in) :: xc_fraction  !< fraction of exact exchange
        real, intent (in) :: r, z         !< evaluation variables
        real, intent (in) :: r1, r2       !< location of atom centers

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer igrid                       !< number of grid points
        integer nnz, nnrho                  !< number of intergration points

        real dmax                           !< max distance between two centers
        real drr, drho, dz                  !< distance between mesh points
        real rin                            !< value of r in Bohr radii

        real zmin, zmax
        real rhomin, rhomax

! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp

        real density_z, density_zz
        real density_pz

! Output from calling get_potxc_1c and potxc_2c
        real exc_1c, vxc_1c, dnuxc_1c, dnuxcs_1c, dexc_1c
        real exc_2c, vxc_2c, dnuxc_2c, dnuxcs_2c

        type (T_rho_2c_bundle), pointer :: prho_bundle
        type (T_rho_2c_store), pointer :: prho_2c

! Procedure
! ======================================================================
! Find which grid point we are located:
        dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
        drr = dmax/float(ndd_uxc - 1)
        igrid = int(d/drr + 0.6) + 1

! Set up pointer for which grid point:
        prho_bundle=>rho_2c_bundle(ispecies, jspecies)
        prho_2c=>prho_bundle%rho_2c_store(igrid)

! Set up grid loop control constants
        rhomin = prho_2c%rhomin
        rhomax = prho_2c%rhomax
        nnrho = prho_2c%nnrho
        drho = prho_2c%drho

        zmin = prho_2c%zmin
        zmax = prho_2c%zmax
        nnz = prho_2c%nnz
        dz = prho_2c%dz

! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                                 nnrho, nnz, prho_2c%rho, density)

! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 .or.               &
            iexc .eq. 9 .or. iexc .eq. 10) then
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
     &                             dz, nnrho, nnz, prho_2c%rhop, density_p)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhopp, density_pp)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhoz, density_z)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhozz, density_zz)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhopz, density_pz)
        else
          density_p = 0.0d0
          density_pp = 0.0d0

          density_z = 0.0d0
          density_zz = 0.0d0
          density_pz = 0.0d0
        end if

! Convert to atomic units
        rin = r/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5

        density_z = density_z*P_abohr**4
        density_zz = density_zz*P_abohr**5
        density_pz = density_pz*P_abohr**5

! Two-center piece: vxc[n1 + n2(r,z)]
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).

! Here energy and potential due to exchange and correlation are calculated.
        call get_potxc_2c (iexc, xc_fraction, rin, density, density_p,       &
     &                     density_pp, density_z, density_zz, density_pz,    &
     &                     exc_2c, vxc_2c, dnuxc_2c, dnuxcs_2c)

! Answers are in Hartrees convert to eV.
        dexc_2c = P_Hartree*(density/P_abohr**3)*(exc_2c - vxc_2c)

! ***************************************************************************
! Subtract off the one-center pieces.
! ***************************************************************************
! We cannot just call dvxc_1c and dexc_1c here because we need to multiply
! by the densities of the atom subtracted.  Just recalculate dvxc_1c and
! dexc_1c for simplicity.
! ***************************************************************************

! One-center piece: vxc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Establish drho for this one-center case.
        drho = species(ispecies)%rcutoffA_max/dfloat(nrho_rho_store)

! Evaluate the density for the one-center - in1
        call rho_1c (ispecies, r1, drho, density, density_p, density_pp)

        rin = r1/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5

        exc_1c = 0.0d0
        vxc_1c = 0.0d0
        call get_potxc_1c (iexc, xc_fraction, rin, density, density_p,       &
     &                     density_pp, exc_1c, vxc_1c, dnuxc_1c, dnuxcs_1c, dexc_1c)

! Answers are in Hartrees convert to eV.
        dexc_2c = dexc_2c - P_Hartree*(density/P_abohr**3)*(exc_1c - vxc_1c)

! One-center piece: vxc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Establish drho for this one-center case.
        drho = species(jspecies)%rcutoffA_max/dfloat(nrho_rho_store)

! Evaluate the density for the one-center - in1
        call rho_1c (jspecies, r2, drho, density, density_p, density_pp)

        rin = r2/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5

        exc_1c = 0.0d0
        vxc_1c = 0.0d0
        call get_potxc_1c (iexc, xc_fraction, rin, density, density_p,       &
     &                     density_pp, exc_1c, vxc_1c, dnuxc_1c, dnuxcs_1c, dexc_1c)

! Answers are in Hartrees convert to eV.
        dexc_2c = dexc_2c - P_Hartree*(density/P_abohr**3)*(exc_1c - vxc_1c)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end function dexc_2c


! ===========================================================================
! vxc_ontop_Harris
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
        subroutine vxc_ontop_Harris
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
! We are doing only Harris here, so set ideriv = 0 and isorp = 0
        ideriv = 0
        isorp = 0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

            call make_munu (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

! begin iammaster
            if (my_proc .eq. 0)                                              &
     &        write (ilogfile,200) ideriv, species(ispecies)%nZ, species(jspecies)%nZ

            ! Open ouput file for this species pair
            write (filename, '("/vxc_ontop_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &  	       ideriv, species(ispecies)%nZ, species(jspecies)%nZ
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
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
            write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')             &
     &        species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &            status = 'unknown', position = 'append')
            write (13,100) pFdata_bundle%nFdata_cell_2c, P_vxc_ontop, ideriv,&
     &                     filename(2:30), pFdata_cell%nME, ndd_vxc, dmax
            close (unit = 13)

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")') &
     &             P_vxc_ontop, species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                    index_2c = 1, nME2c_max)
            close (unit = 12)

! Loop over grid
            do igrid = 1, ndd_vxc
              d = d + drr

              ! Set integration limits
              zmin = max(-rcutoff1, d - rcutoff2)
              zmax = min(rcutoff1, d + rcutoff2)

              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, nz_vxc, nrho_vxc,                &
     &                                   rint_vxc_ontop, phifactor, zmin,    &
     &                                   zmax, rhomin, rhomax, pFdata_cell%fofx)

              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                      &
     &                                       index_2c = 1, nME2c_max)
            end do ! igrid
            write (11,*)
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating vxc ontop integrals for ideriv = ', i3,     &
     &              ' nZ = ', i3, ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vxc_ontop_Harris


! ===========================================================================
! rint_vxc_ontop
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
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
! Program Declaration
! ===========================================================================
        real function rint_vxc_ontop (itype, ispecies, jspecies, isorp, d,   &
     &                                rho, z1, z2, ideriv, index_2c)
        implicit none

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent (in) ::  d           !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer iexc                     ! which exchange-correlation option
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1, r2                      ! location of the atom centers
        real rmax                        ! maximum value of r along grid

        real vofr
        real xc_fraction

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
        idummy = isorp
        dummy = d

! Find which exchange-correlation we are calcuating:
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

! Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)
        r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n1)%rcutoffA
        psi1val = psiofr (r1, rmax, ispecies, n1)

! Find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        rmax = species(jspecies)%shell(n2)%rcutoffA
        psi2val = psiofr (r2, rmax, jspecies, n2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
        psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
        psi2val = rescaled_psi (l2, m2, rho, r2, z2, psi2val)
        vofr = vxc_2c (iexc, xc_fraction, ispecies, jspecies, rho, z1, d)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_vxc_ontop = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_vxc_ontop


! ===========================================================================
! vxc_2c
! ===========================================================================
! Program Description
! ===========================================================================
! This function computes vxc(n1 + n2) with the densities ni of atom
! in_i at the distance r_i from their centers.
!
! On input:   r, z: geometry information for the charge gradient
!
! On output:  vxc: vxc[n1(r1) + n2(r2)]

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
! Program Declaration
! ===========================================================================
        real function vxc_2c (iexc, xc_fraction, ispecies, jspecies, r, z, d)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc      !< which exchange-correlation
        integer, intent (in) :: ispecies, jspecies      !< two centers species

        real, intent (in) :: d            !< distance between the two centers
        real, intent (in) :: xc_fraction  !< fraction of exact exchange
        real, intent (in) :: r, z         !< evaluation variables

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer igrid                       !< number of grid points
        integer nnz, nnrho                  !< number of intergration points

        real dmax                           !< max distance between two centers
        real drr, drho, dz                  !< distance between mesh points
        real rin                            !< value of r in Bohr radii

        real zmin, zmax
        real rhomin, rhomax

! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp

        real density_z, density_zz
        real density_pz

! Output from calling get_potxc_1c
        real dnuxc, dnuxcs
        real exc, vxc

        type (T_rho_2c_bundle), pointer :: prho_bundle
        type (T_rho_2c_store), pointer :: prho_2c

! Procedure
! ======================================================================
! Find which grid point we are located:
        dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
        drr = dmax/float(ndd_vxc - 1)
        igrid = int(d/drr + 0.6) + 1

! Set up pointer for which grid point:
        prho_bundle=>rho_2c_bundle(ispecies, jspecies)
        prho_2c=>prho_bundle%rho_2c_store(igrid)

! Set up grid loop control constants
        rhomin = prho_2c%rhomin
        rhomax = prho_2c%rhomax
        nnrho = prho_2c%nnrho
        drho = prho_2c%drho

        zmin = prho_2c%zmin
        zmax = prho_2c%zmax
        nnz = prho_2c%nnz
        dz = prho_2c%dz

! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                                 nnrho, nnz, prho_2c%rho, density)

! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 .or.               &
            iexc .eq. 9 .or. iexc .eq. 10) then
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
     &                             dz, nnrho, nnz, prho_2c%rhop, density_p)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhopp, density_pp)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhoz, density_z)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhozz, density_zz)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhopz, density_pz)
        else
          density_p = 0.0d0
          density_pp = 0.0d0

          density_z = 0.0d0
          density_zz = 0.0d0
          density_pz = 0.0d0
        end if

! Convert to atomic units
        rin = r/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5

        density_z = density_z*P_abohr**4
        density_zz = density_zz*P_abohr**5
        density_pz = density_pz*P_abohr**5

! Two-center piece: vxc[n1 + n2(r,z)]
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).

! Here energy and potential due to exchange and correlation are calculated.
        call get_potxc_2c (iexc, xc_fraction, rin, density, density_p,       &
     &                     density_pp, density_z, density_zz, density_pz,    &
     &                     exc, vxc, dnuxc, dnuxcs)

! Answers are in Hartrees convert to eV.
        vxc_2c = P_Hartree*vxc

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end function vxc_2c


! ===========================================================================
! vxc_atom_Harris
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
        subroutine vxc_atom_Harris
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
! We are doing only Harris here, so set ideriv = 0 and isorp = 0
        ideriv = 0
        isorp = 0

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

            call make_munu_atom (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

! begin iammaster
            if (my_proc .eq. 0)                                              &
     &        write (ilogfile,200) ideriv, species(ispecies)%nZ, species(jspecies)%nZ

            ! Open ouput file for this species pair
            write (filename, '("/vxc_atom_", i2.2,".",i2.2,".",i2.2,".dat")')&
     &             ideriv, species(ispecies)%nZ, species(jspecies)%nZ
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown')

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + na(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_vxc - 1)
            d = -drr

            ! Set final integration limit
            rhomax = min(rcutoff1, rcutoff2)

            ! open directory file
            write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')             &
     &        species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &            status = 'unknown', position = 'append')
            write (13,100) pFdata_bundle%nFdata_cell_2c, P_vxc_atom, ideriv, &
     &                     filename(2:30), pFdata_cell%nME, ndd_vxc, dmax
            close (unit = 13)

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")') &
     &             P_vxc_atom, species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                    index_2c = 1, nME2c_max)
            close (unit = 12)

! Loop over grid
            do igrid = 1, ndd_vxc
              d = d + drr

              ! Set integration limits
              zmin = max(-rcutoff1, d - rcutoff2)
              zmax = min(rcutoff1, d + rcutoff2)

              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies, &
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, nz_vxc, nrho_vxc,                &
     &                                   rint_vxc_atom, phifactor, zmin,     &
     &                                   zmax, rhomin, rhomax, pFdata_cell%fofx)

              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                      &
     &                                       index_2c = 1, nME2c_max)
            end do ! igrid
            write (11,*)
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating vxc atom integrals for ideriv = ', i3,      &
                    ' nZ = ', i3, ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine vxc_atom_Harris


! ===========================================================================
! rint_vxc_atom
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
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
! Program Declaration
! ===========================================================================
        real function rint_vxc_atom (itype, ispecies, jspecies, isorp, d,    &
     &                               rho, z1, z2, ideriv, index_2c)
        implicit none

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent (in) ::  d           !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer iexc                     ! which exchange-correlation option
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1, r2                      ! location of the atom centers
        real rmax                        ! maximum value of r along grid

        real vofr
        real xc_fraction

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
        idummy = isorp
        dummy = d

! Find which exchange-correlation we are calcuating:
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

!Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)
        r2 = sqrt(z2**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n1)%rcutoffA
        psi1val = psiofr (r1, rmax, ispecies, n1)

!find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n2)%rcutoffA
        psi2val = psiofr (r1, rmax, ispecies, n2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
        psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
        psi2val = rescaled_psi (l2, m2, rho, r1, z1, psi2val)
        vofr = dvxc_2c (iexc, xc_fraction, ispecies, jspecies, rho, z1, d, r1)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_vxc_atom = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_vxc_atom


! ===========================================================================
! dvxc_2c
! ===========================================================================
! Program Description
! ===========================================================================
! This function computes the atom case of vxc(n1 + n2) with the densities
! and hence the potential at one atom and the wavefunctions at a neighboring
! atom in_i at the distance r_i from their centers.
!
! On input:   r, z: geometry information for the charge gradient
!
! On output:  vxc: vxc[n1(r1) + n2(r2)]
!
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
! Program Declaration
! ===========================================================================
        real function dvxc_2c (iexc, xc_fraction, ispecies, jspecies,        &
     &                         r, z, d, r1)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc      !< which exchange-correlation
        integer, intent (in) :: ispecies, jspecies      !< two centers species

        real, intent (in) :: d            !< distance between the two centers
        real, intent (in) :: xc_fraction  !< fraction of exact exchange
        real, intent (in) :: r, z         !< evaluation variables
        real, intent (in) :: r1           !< location of the atom center

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer igrid                       !< number of grid points
        integer nnz, nnrho                  !< number of intergration points

        real dmax                           !< max distance between two centers
        real drr, drho, dz                  !< distance between mesh points
        real rin                            !< value of r in Bohr radii

        real zmin, zmax
        real rhomin, rhomax

! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp

        real density_z, density_zz
        real density_pz

! Output from calling get_potxc_2c
        real dnuxc_2c, dnuxcs_2c
        real exc_2c, vxc_2c

        type (T_rho_2c_bundle), pointer :: prho_bundle
        type (T_rho_2c_store), pointer :: prho_2c

! Procedure
! ======================================================================
! Find which grid point we are located:
        dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
        drr = dmax/float(ndd_vxc - 1)
        igrid = int(d/drr + 0.6) + 1

! Set up pointer for which grid point:
        prho_bundle=>rho_2c_bundle(ispecies, jspecies)
        prho_2c=>prho_bundle%rho_2c_store(igrid)

! Set up grid loop control constants
        rhomin = prho_2c%rhomin
        rhomax = prho_2c%rhomax
        nnrho = prho_2c%nnrho
        drho = prho_2c%drho

        zmin = prho_2c%zmin
        zmax = prho_2c%zmax
        nnz = prho_2c%nnz
        dz = prho_2c%dz

! Interpolate the density and gradients of the density at the given
! point (r, z).
        call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax, dz, &
                                 nnrho, nnz, prho_2c%rho, density)

! Only interpolate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 .or.               &
            iexc .eq. 9 .or. iexc .eq. 10) then
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
     &                             dz, nnrho, nnz, prho_2c%rhop, density_p)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhopp, density_pp)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhoz, density_z)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhozz, density_zz)
          call interpolate_rho_2c (r, rhomin, rhomax, drho, z, zmin, zmax,   &
                                   dz, nnrho, nnz, prho_2c%rhopz, density_pz)
        else
          density_p = 0.0d0
          density_pp = 0.0d0
          density_z = 0.0d0
          density_zz = 0.0d0
          density_pz = 0.0d0
        end if

! Convert to atomic units
        rin = r/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5
        density_z = density_z*P_abohr**4
        density_zz = density_zz*P_abohr**5
        density_pz = density_pz*P_abohr**5

! Two-center piece: vxc[n1 + n2(r,z)]
! ***************************************************************************
! Interpolate the density and gradients of the density at the given
! point (r, z).

! Here energy and potential due to exchange and correlation are calculated.
        call get_potxc_2c (iexc, xc_fraction, rin, density, density_p,       &
     &                     density_pp, density_z, density_zz, density_pz,    &
     &                     exc_2c, vxc_2c, dnuxc_2c, dnuxcs_2c)

! Answers are in Hartrees convert to eV.
        dvxc_2c = P_Hartree*vxc_2c - dvxc_1c (iexc, xc_fraction, ispecies, r1)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end function dvxc_2c


! ===========================================================================
! interpolate_rho_2c
! ===========================================================================
! Program Description
! ===========================================================================
!       This routine is a two-dimensional interpolater on a 4x4 sub-grid for
! the density and corresponding derivatives.
!
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
! Program Declaration
! ===========================================================================
        subroutine interpolate_rho_2c (rho, rhomin, rhomax, drho, z, zmin,  &
                                       zmax, dz, nnrho, nnz, frho, answer)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nnrho, nnz

        real, intent (in) :: drho, dz
        real, intent (in) :: rho, rhomin, rhomax
        real, intent (in) :: z, zmin, zmax

! This is the function being interpolated
        real, intent (in) :: frho (nnrho, nnz)

! Output
        real answer

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer imidrho
        integer imidz
        integer j
        integer jj
        integer k

        real gradrho, gradz, gradtest    ! test gradient to simply interpolation

        real e36t
        real f0p3
        real f0p6
        real f1m2
        real f1m3
        real f1p3
        real f1p6
        real ftp
        real prho
        real prod
        real pz
        real tp

        real b (0:5)
        real bb (0:5, -2:3)
        real g (-2:3)
        real fun (-1:2, -1:2)

! Procedure
! ===========================================================================
! Check and make sure that the point (r, z) is within the limits of the
! stored density.
        if (rho .gt. (rhomax + 1.0d-5) .or. rho .lt. rhomin) then
         write (*,*) ' What the heck is going on in interpolate_rho_2c !'
         write (*,*) ' ************* error !!! ************* '
         write (*,*) ' Input rho = ', rho
         write (*,*) ' Max. data = ', rhomax, ' Min. data = ', rhomin
        end if

        if (z .gt. (zmax + 1.0d-5) .or. z .lt. (zmin - 1.0d-5)) then
         write (*,*) ' What the heck is going on in interpolate_rho_2c !'
         write (*,*) ' ************* error !!! ************* '
         write (*,*) ' Input z = ', z
         write (*,*) ' Max. data = ', zmax, ' Min. data = ', zmin
        end if

! Much of the code in this routine is a dodge to avoid interpolating, even
! though the mission of this subprogram is to interpolate. The dodge can be
! approached at two levels: first, if the gradient is locally very small, we
! can just return. This is CRUCIAL for multirc problems, but seems to help for
! monorc.

! Set-up everything for interpolation.
        imidrho = int((rho - rhomin)/drho) + 1
        imidz = int((z - zmin)/dz) + 1

        if (imidrho .lt. 2) imidrho = 2
        if (imidrho .gt. nnrho - 2) imidrho = nnrho - 2
        if (imidz .lt. 2) imidz = 2
        if (imidz .gt. nnz - 2) imidz = nnz - 2

        prho = (rho - rhomin)/drho - float(imidrho - 1)
        pz = (z - zmin)/dz - float(imidz - 1)

        fun(-1,-1) = frho(imidrho - 1, imidz - 1)
        fun(-1, 0) = frho(imidrho - 1, imidz)
        fun(-1, 1) = frho(imidrho - 1, imidz + 1)
        fun(-1, 2) = frho(imidrho - 1, imidz + 2)

        fun(0,-1) = frho(imidrho, imidz - 1)
        fun(0, 0) = frho(imidrho, imidz)
        fun(0, 1) = frho(imidrho, imidz + 1)
        fun(0, 2) = frho(imidrho, imidz + 2)

        fun(1,-1) = frho(imidrho + 1, imidz - 1)
        fun(1, 0) = frho(imidrho + 1, imidz)
        fun(1, 1) = frho(imidrho + 1, imidz + 1)
        fun(1, 2) = frho(imidrho + 1, imidz + 2)

        fun(2,-1) = frho(imidrho + 2, imidz - 1)
        fun(2, 0) = frho(imidrho + 2, imidz)
        fun(2, 1) = frho(imidrho + 2, imidz + 1)
        fun(2, 2) = frho(imidrho + 2, imidz + 2)

! **************************************************************************
! If the gradient is small, then do quadratic quick bivariate interpolation.
        gradrho = (fun(0,0) - fun(1,0))/drho
        gradz = (fun(0,0) - fun(0,1))/dz

        gradtest = abs(gradrho) + abs(gradz)

! Form the criterion for a quick interpolation. Empirically, I find that
! gradtest < 1.0d-05 is adequate. If you dont want this option, change
! 1.0d-05 in the next line to 0.0d0!
        if (gradtest .lt. 1.0d-05) then
          answer = (1.0d0 - prho - pz)*fun(0,0) + prho*fun(1,0) + pz*fun(0,1)
          return
        end if

! **************************************************************************
! Phase III. All else fails. Interpolate carefully. Original pfed
! interpolator with minimal multiplies.
        e36t = 1.0d0/36.0d0

        do k = - 1, 2
          f1m2 = fun(k,-1) + fun(k,-1)
          f1m3 = f1m2 + fun(k,-1)

          f0p3 = fun(k,0) + fun(k,0) + fun(k,0)
          f0p6 = f0p3 + f0p3

          f1p3 = fun(k,1) + fun(k,1) + fun(k,1)
          f1p6 = f1p3 + f1p3

          bb(3,k) = - fun(k,-1) + f0p3 - f1p3 + fun(k,2)
          bb(2,k) = f1m3 - f0p6 + f1p3
          bb(1,k) = - f1m2 - f0p3 + f1p6 - fun(k,2)

          tp = fun(k,0)
          tp = tp + tp + tp
          bb(0,k) = tp + tp

          prod = bb(3,k)*pz
          do j = 1, 2
            jj = 3 - j
            ftp = bb(jj,k)
            prod = (prod + ftp)*pz
          end do
          g(k) = prod + bb(0,k)
        end do

        f1m2 = g(-1) + g(-1)
        f1m3 = f1m2 + g(-1)

        f0p3 = g(0) + g(0) + g(0)
        f0p6 = f0p3 + f0p3

        f1p3 = g(1) + g(1) + g(1)
        f1p6 = f1p3 + f1p3

        b(3) = -g(-1) + f0p3 - f1p3 + g(2)
        b(2) = f1m3 - f0p6 + f1p3
        b(1) = -f1m2 - f0p3 + f1p6 - g(2)
        tp = g(0) + g(0) + g(0)
        b(0) = tp + tp

        prod = b(3)*prho
        do j = 1, 2
          jj = 3 - j
          ftp = b(jj)
          prod = (prod + ftp)*prho
        end do
        prod = prod + b(0)

! Final answer
        answer = e36t*prod

! Format Statements
! ===========================================================================
! None

        return
        end subroutine interpolate_rho_2c


! End Module
! =============================================================================
      end module
