! copyright info:
!
!                             @Copyright 2016
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
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

! M_vxc_DOGS
! Program Description
! ===========================================================================
!      This is a module calculating the integrals of two centers for the
! the exchange-correlation interactions.
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute, Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
! Module Declaration
! ===========================================================================
        module M_vxc_DOGS
        use M_precision
        use M_atom_functions
        use M_species
        use M_integrals_2c
        use M_vnl
        use M_xc_1c
        use M_xc_2c
        use M_vxc_Harris

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
        integer isorp                       !< counter for shells

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

            ! add one for nuxc_1c and exc_1c
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1

            ! add one each for dnuxc_ontopL and dnuxc_ontopR - loop over shells
            do isorp = 1, species(ispecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            end do

            do isorp = 1, species(jspecies)%nssh
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
        end subroutine initialize_vxc_DOGS


! ===========================================================================
! vxc_DOGS
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calls the subroutines required to calculate the vna_ontop
! (both left/right cases) and atom cases.
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine vxc_DOGS
        implicit none


! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer logfile                     !< writing to which unit

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = 21

        write (logfile,*)
        write (logfile,*) ' ******************************************************* '
        write (logfile,*) '        E X C H A N G E   C O R R E L A T I O N          '
        write (logfile,*) '        C H A R G E D   I N T E R A C T I O N S          '
        write (logfile,*) ' ******************************************************* '
        write (logfile,*)

        write (logfile,*) ' Calling one-center case. '
        call nuxc_1c

        write (logfile,*)
        call exc_1c

        write (logfile,*)
        write (logfile,*) ' Calling two-center dnuxc_ontopL case. '
        call dnuxc_ontopL

        write (logfile,*)
        write (logfile,*) ' Calling two-center dnuxc_ontopR case. '
        call dnuxc_ontopR

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
! nuxc_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the one-center
! matrix elements of the form <psi1|dV(1)|psi2>.  The potential dV(1) and
! wavefunctions, psi1 and psi2, are located at the same atom site.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
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
        subroutine nuxc_1c
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
        integer index_1c, nME1c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer logfile                     !< writing to which unit
        integer nFdata_cell_1c              !< indexing of interactions

        real d                              !< distance between the two centers
        real rcutoff1                       !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 25) filename

! Procedure
! ============================================================================
! Initialize logfile
        logfile = 21

! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

! Loop over species
        do ispecies = 1, nspecies
          pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
          nFdata_cell_1c = pFdata_bundle%nFdata_cell_2c
          pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_1c)

          ! Open ouput file for this species pair
          write (filename, '("/nuxcrho_1c", ".", i2.2, ".dat")')             &
     &      species(ispecies)%nZ
          open (unit = 11, file = trim(Fdata_location)//trim(filename),      &
     &          status = 'unknown')

          write (logfile,100) species(ispecies)%nZ

          do isorp = 1, species(ispecies)%nssh
            call make_munuS (nFdata_cell_1c, ispecies, ispecies)
            nME1c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME1c_max))

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max

            ! Set integration limits
            rhomin = 0.0d0
            rhomax = rcutoff1

! Loop over grid
            d = 0.0d0

            ! Set integration limits
            zmin = -rcutoff1
            zmax = rcutoff1

            call evaluate_integral_2c (nFdata_cell_1c, ispecies, ispecies,   &
     &                                 isorp, ideriv, rcutoff1, rcutoff1, d, &
     &                                 nz_vxc, nrho_vxc, rint_nuxc_1c,       &
     &                                 phifactor, zmin, zmax, rhomin, rhomax,&
     &                                 pFdata_cell%fofx)

            ! Write out details.
            index_1c = 1
            do issh = 1, species(ispecies)%nssh
              write (11,*) (pFdata_cell%fofx(jssh),                          &
     &                      jssh = index_1c, index_1c + species(ispecies)%nssh - 1)
              index_1c = index_1c + species(ispecies)%nssh
            end do

          end do ! isorp

          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c - 1
        end do ! ispecies


! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Evaluating nuxc_1c integrals for nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine nuxc_1c


! ===========================================================================
! rint_nuxc_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This is the radial integral function for the one-center exchange-
! correlation interactions.
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
        real function rint_nuxc_1c (itype, ispecies, jspecies, isorp, d, rho,&
     &                                 z1, z2, ideriv, index_2c)
        implicit none

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent (in) :: d            !< distance between the two centers
        real, intent (in) :: rho
        real, intent (in) :: z1, z2

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
        real r1
        real rmax                        ! maximum value of r along grid

        real vofr
        real xc_fraction

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
        dummy = d
        if (.FALSE.) dummy = jspecies

! Find which exchange-correlation we are calcuating:
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

! Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)+0.0*z2

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
        rmax = species(ispecies)%shell(n2)%rcutoffA
        psi2val = psiofr (r1, rmax, ispecies, n2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
        psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
        psi2val = rescaled_psi (l2, m2, rho, r1, z1, psi2val)

! One-center piece: exc[n1(r1)]
        vofr = dnuxc_1c (iexc, xc_fraction, ispecies, r1)
        rmax = species(ispecies)%shell(isorp)%rcutoffA
        vofr = vofr*(psiofr(r1, rmax, ispecies, isorp))**2/(4.0d0*pi)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_nuxc_1c = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_nuxc_1c


! ===========================================================================
! dnuxc_1c
! ===========================================================================
! Program Description
! ===========================================================================
! This function computes only vxc(n1) with the densities ni of atom in_i.
!
! On input:  in1: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  dnuxc: nuxc[n1(r1)]  (the one center nuxc)

! We calculate vxc(n1). We compute neutral atoms only.
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
!
! Program Declaration
! ===========================================================================
        real function dnuxc_1c (iexc, xc_fraction, ispecies, r1)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc      !< which exchange-correlation
        integer, intent (in) :: ispecies  !< two centers species

        real, intent (in) :: xc_fraction  !< fraction of exact exchange
        real, intent (in) :: r1           !< evaluation variables

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real drho                           !< distance between mesh points
        real rin                            !< value of r in Bohr radii
        real rcutoff1                       !< cutoffs for the two centers

! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp

        real exc
        real vxc
        real dnuxc
        real dnuxcs
        real dexc

! Procedure
! ===========================================================================
! Establish drho for this one-center case.
        rcutoff1 = species(ispecies)%rcutoffA_max
        drho = rcutoff1/dfloat(nrho_rho_store)

! One-center piece: vxc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call rho_1c (ispecies, r1, drho, density, density_p, density_pp)

        rin = r1/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5
        call get_potxc_1c (iexc, xc_fraction, rin, density, density_p,       &
     &                     density_pp, exc, vxc, dnuxc, dnuxcs, dexc)

! Answers are in Hartrees convert to eV.
        dnuxc_1c = P_hartree*dnuxc*P_abohr**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end function dnuxc_1c


! ===========================================================================
! exc_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the one-center
! exchange-correlation interactions - dE.  The potential dE(1) and
! wavefunctions, psi1 and psi2, are located at the same atom site.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
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
        subroutine exc_1c
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
        integer index_1c, nME1c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer logfile                     !< writing to which unit
        integer nFdata_cell_1c              !< indexing of interactions

        real d                              !< distance between the two centers
        real rcutoff1                       !< cutoffs for the two centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 25) filename

! Procedure
! ============================================================================
! Initialize logfile
        logfile = 21

! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

! Loop over species
        do ispecies = 1, nspecies
          pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
          nFdata_cell_1c = pFdata_bundle%nFdata_cell_2c
          pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_1c)

          ! Open ouput file for this species pair
          write (filename, '("/excrho_1c", ".", i2.2, ".dat")')              &
     &      species(ispecies)%nZ
          open (unit = 11, file = trim(Fdata_location)//trim(filename),      &
     &          status = 'unknown')

          write (logfile,100) species(ispecies)%nZ

          do isorp = 1, species(ispecies)%nssh
            call make_munuS (nFdata_cell_1c, ispecies, ispecies)
            nME1c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME1c_max))

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max

            ! Set integration limits
            rhomin = 0.0d0
            rhomax = rcutoff1

! Loop over grid
            d = 0.0d0

            ! Set integration limits
            zmin = -rcutoff1
            zmax = rcutoff1

            call evaluate_integral_2c (nFdata_cell_1c, ispecies, ispecies,   &
     &                                 isorp, ideriv, rcutoff1, rcutoff1, d, &
     &                                 nz_vxc, nrho_vxc, rint_dexc_1c,       &
     &                                 phifactor, zmin, zmax, rhomin, rhomax,&
     &                                 pFdata_cell%fofx)

            ! Write out details.
            index_1c = 1
            do issh = 1, species(ispecies)%nssh
              write (11,*) (pFdata_cell%fofx(jssh),                          &
     &                      jssh = index_1c, index_1c + species(ispecies)%nssh - 1)
              index_1c = index_1c + species(ispecies)%nssh
            end do

          end do ! isorp

          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c - 1
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Evaluating exc_1c integrals for nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine exc_1c


! ===========================================================================
! rint_dexc_1c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This is the radial integral function for the one-center exchange-
! correlation interactions.
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
        real function rint_dexc_1c (itype, ispecies, jspecies, isorp, d, rho,&
     &                                 z1, z2, ideriv, index_2c)
        implicit none

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent (in) :: d               !< distance between the two centers
        real, intent (in) :: rho
        real, intent (in) :: z1, z2

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
        real r1
        real rmax                         ! maximum value of r along grid

        real vofr
        real xc_fraction

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
        dummy = d
        if (.false.) dummy = jspecies

! Find which exchange-correlation we are calcuating:
        iexc = species_PP(ispecies)%iexc
        xc_fraction = species_PP(ispecies)%xc_fraction

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)

! Pick up "passed" data
        n1 =  pFdata_cell%N_mu(index_2c)
        l1 =  pFdata_cell%L_mu(index_2c)
        m1 =  pFdata_cell%M_mu(index_2c)

        n2 =  pFdata_cell%N_nu(index_2c)
        l2 =  pFdata_cell%L_nu(index_2c)
        m2 =  pFdata_cell%M_nu(index_2c)

! Set parameters for actual function
        r1 = sqrt(z1**2 + rho**2)+0.0*z2

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
        rmax = species(ispecies)%shell(n2)%rcutoffA
        psi2val = psiofr (r1, rmax, ispecies, n2)

! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
        psi1val = rescaled_psi (l1, m1, rho, r1, z1, psi1val)
        psi2val = rescaled_psi (l2, m2, rho, r1, z1, psi2val)

! One-center piece: exc[n1(r1)]
        vofr = ddexc_1c (iexc, xc_fraction, ispecies, r1)
        rmax = species(ispecies)%shell(isorp)%rcutoffA
        vofr = vofr*(psiofr(r1, rmax, ispecies, isorp))**2/(4.0d0*pi)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_dexc_1c = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_dexc_1c


! ===========================================================================
! ddexc_1c
! ===========================================================================
! Program Description
! ===========================================================================
! This function computes only vxc(n1) with the densities ni of atom in_i.
!
! On input:  in1: atomic indices
!            r, z: geometry information for the charge gradient
!            ix: switch for the derivatives
!
! On output:  dnuxc: nuxc[n1(r1)]  (the one center nuxc)

! We calculate vxc(n1). We compute neutral atoms only.
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
!
! Program Declaration
! ===========================================================================
        real function ddexc_1c (iexc, xc_fraction, ispecies, r1)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc      !< which exchange-correlation
        integer, intent (in) :: ispecies  !< two centers species

        real, intent (in) :: xc_fraction  !< fraction of exact exchange
        real, intent (in) :: r1           !< evaluation variables

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        real drho                           !< distance between mesh points
        real rin                            !< value of r in Bohr radii
        real rcutoff1                       !< cutoffs for the two centers

! Value of density and corresponding derivatives at the point r, z
        real density
        real density_p, density_pp

        real exc
        real vxc
        real dnuxc
        real dnuxcs
        real dexc

! Procedure
! ===========================================================================
! Establish drho for this one-center case.
        rcutoff1 = species(ispecies)%rcutoffA_max
        drho = rcutoff1/dfloat(nrho_rho_store)

! One-center piece: vxc[n1(r1)]
! Compute the exchange correlation potential for the one-center case
! ***************************************************************************
! Evaluate the density for the one-center - in1
        call rho_1c (ispecies, r1, drho, density, density_p, density_pp)

        rin = r1/P_abohr
        density = density*P_abohr**3
        density_p = density_p*P_abohr**4
        density_pp = density_pp*P_abohr**5
        call get_potxc_1c (iexc, xc_fraction, rin, density, density_p,       &
     &                     density_pp, exc, vxc, dnuxc, dnuxcs, dexc)

! Answers are in Hartrees convert to eV.
        ddexc_1c = P_hartree*dexc*P_abohr**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end function ddexc_1c


! ===========================================================================
! dnuxc_ontopL
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
        subroutine dnuxc_ontopL
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
        integer logfile                     !< writing to which unit
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
! Initialize logfile
        logfile = 21

! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999
        if (.false.) ideriv = jspecies


! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_vxc - 1)

            ! Set integration limits
            rhomin = 0.0d0
            rhomax = min(rcutoff1, rcutoff2)

! Loop over shells
            do isorp = 1, species(ispecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
              nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

              call make_munu (nFdata_cell_2c, ispecies, jspecies)
              nME2c_max = pFdata_cell%nME
              allocate (pFdata_cell%fofx(nME2c_max))

              ! Open ouput file for this species pair
              write (filename, '("/dnuxc_ontopL_", i2.2, ".", i2.2, ".",     &
     &                           i2.2, ".dat")')                             &
     &          isorp, species(ispecies)%nZ, species(jspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
              if (skip) cycle
              open (unit = 11, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown')

              ! open directory file
              write (interactions, '("/2c.", i2.2, ".", i2.2, ".dir")')      &
     &          species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown', position = 'append')
              write (13,100) pFdata_bundle%nFdata_cell_2c, P_dnuxc_ontopL,   &
     &          isorp, filename(2:30), pFdata_cell%nME, ndd_vxc, dmax
              close (unit = 13)

              ! Open mu, nu, mvalue file and write out values.
              write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")')&
     &               P_dnuxc_ontopL, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 12, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown', position = 'append')

              ! write the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                      index_2c = 1, nME2c_max)

! Loop over grid
              write (logfile,200) species(ispecies)%nZ, species(jspecies)%nZ, isorp
              d = -drr
              do igrid = 1, ndd_vxc
                d = d + drr

                ! Set integration limits
                zmin = min(-rcutoff1, d - rcutoff2)
                zmax = max(rcutoff1, d + rcutoff2)

                call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies,&
     &                                     isorp, ideriv, rcutoff1, rcutoff2, &
     &                                     d, nz_vxc, nrho_vxc,               &
     &                                     rint_dnuxc_ontopL, phifactor, zmin,&
     &                                     zmax, rhomin, rhomax, pFdata_cell%fofx)

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
200     format (2x, ' Evaluating dnuxc ontopL integrals for nZ = ', i3,      &
     &              ' and nZ = ', i3, ', isorp = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine dnuxc_ontopL


! ===========================================================================
! rint_dnuxc_ontopL
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
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
!
! Program Declaration
! ===========================================================================
        real function rint_dnuxc_ontopL (itype, ispecies, jspecies, isorp, d,&
     &                                   rho, z1, z2, ideriv, index_2c)
        implicit none

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent (in) :: d            !< distance between the two centers
        real, intent (in) :: rho
        real, intent (in) :: z1, z2

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
        real r1, r2
        real rmax                        ! maximum value of r along grid

        real vofr
        real xc_fraction

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
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

! Multiply by the "left" density to get dnuxc*rho
        vofr = dnuxc (iexc, xc_fraction, ispecies, jspecies, rho, z1, d)
        rmax = species(ispecies)%shell(isorp)%rcutoffA
        vofr = vofr*(psiofr(r1, rmax, ispecies, isorp))**2/(4.0d0*pi)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_dnuxc_ontopL = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_dnuxc_ontopL


! ===========================================================================
! dnuxc_ontopR
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the right site of one of the orbitals.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
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
        subroutine dnuxc_ontopR
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
        integer logfile                     !< writing to which unit
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
! Initialize logfile
        logfile = 21

! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_vxc - 1)

            ! Set integration limits
            rhomin = 0.0d0
            rhomax = min(rcutoff1, rcutoff2)

! Loop over shells
            do isorp = 1, species(jspecies)%nssh
              pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
              nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
              pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

              call make_munu (nFdata_cell_2c, ispecies, jspecies)
              nME2c_max = pFdata_cell%nME

              allocate (pFdata_cell%fofx(nME2c_max))

              ! Open ouput file for this species pair
              write (filename, '("/dnuxc_ontopR_", i2.2, ".", i2.2, ".",     &
     &                           i2.2, ".dat")')                             &
     &          isorp, species(ispecies)%nZ, species(jspecies)%nZ
              inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
              if (skip) cycle
              open (unit = 11, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown')

              ! open directory file
              write (interactions, '("/2c.", i2.2, ".", i2.2, ".dir")')      &
     &          species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &              status = 'unknown', position = 'append')
              write (13,100) pFdata_bundle%nFdata_cell_2c, P_dnuxc_ontopR,   &
     &          isorp, filename(2:30), pFdata_cell%nME, ndd_vxc, dmax
              close (unit = 13)

              ! Open mu, nu, mvalue file and write out values.
              write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")')&
     &               P_dnuxc_ontopR, species(ispecies)%nZ, species(jspecies)%nZ
              open (unit = 12, file = trim(Fdata_location)//trim(filename),  &
     &              status = 'unknown', position = 'append')

              ! write the mapping - stored in mu, nu, and mvalue
              write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
              write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                      index_2c = 1, nME2c_max)

! Loop over grid
              write (logfile,200) species(ispecies)%nZ, species(jspecies)%nZ, isorp
              d = -drr
              do igrid = 1, ndd_vxc
                d = d + drr

                ! Set integration limits
                zmin = min(-rcutoff1, d - rcutoff2)
                zmax = max(rcutoff1, d + rcutoff2)

                call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies,&
     &                                     isorp, ideriv, rcutoff1, rcutoff2, &
     &                                     d, nz_vxc, nrho_vxc,               &
     &                                     rint_dnuxc_ontopR, phifactor, zmin,&
     &                                     zmax, rhomin, rhomax, pFdata_cell%fofx)

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
200     format (2x, ' Evaluating dnuxc ontopR integrals for nZ = ', i3,      &
     &              ' and nZ = ', i3, ', isorp = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine dnuxc_ontopR


! ===========================================================================
! rint_dnuxc_ontopR
! ===========================================================================
! Program Description
! ===========================================================================
! The dnuxc_ontopR part of the twocenter_overlap.
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
!
! Program Declaration
! ===========================================================================
        real function rint_dnuxc_ontopR (itype, ispecies, jspecies, isorp, d,&
     &                                   rho, z1, z2, ideriv, index_2c)
        implicit none

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent (in) :: d            !< distance between the two centers
        real, intent (in) :: rho
        real, intent (in) :: z1, z2

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
        real r1, r2
        real rmax                        ! maximum values of r along grid

        real vofr
        real xc_fraction

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================
! Initialize some dummy variables for warning removal
        idummy = ideriv
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

! Multiply by the "left" density to get dnuxc*rho
        vofr = dnuxc (iexc, xc_fraction, ispecies, jspecies, rho, z1, d)
        rmax = species(jspecies)%shell(isorp)%rcutoffA
        vofr = vofr*(psiofr(r2, rmax, jspecies, isorp))**2/(4.0d0*pi)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_dnuxc_ontopR = psi1val*vofr*psi2val*rho

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_dnuxc_ontopR


! ===========================================================================
! dnuxc
! ===========================================================================
! Program Description
! ===========================================================================
! This function computes vxc(n1 + n2) with the densities ni of atom
! in_i at the distance r_i from their centers.
!
! On input:   r, z: geometry information for the charge gradient
!
! On output:  dnuxc: dnuxc[n1(r1) + n2(r2)]

! We calculate dnuxc(n1+n2). We compute neutral atoms here only.
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
!
! Program Declaration
! ===========================================================================
        real function dnuxc (iexc, xc_fraction, ispecies, jspecies, r, z, d)
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
        prho_2c => prho_bundle%rho_2c_store(igrid)

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
        dnuxc = P_hartree*dnuxc_2c*P_abohr**3

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end function dnuxc


! End Module
! =============================================================================
        end module
