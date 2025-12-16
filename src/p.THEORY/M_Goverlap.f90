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

! M_Goverlap_1c.f90
! Program Description
! ===========================================================================
!      This is a module calculating the integrals for calculating the overlap
! of a matrix element and the derivative of a matrix element
! - these are two center interactions.
!
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
!
! ===========================================================================
! Module Declaration
! ===========================================================================
        module M_Goverlap

! /GLOBAL
        use M_precision

! /SYSTEM
        use M_species
        use M_integrals_2c

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! initialize_overlap
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
!> @author Amanda Neukirch
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
        subroutine initialize_goverlap
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies          !< counters for number of species

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species
        do ispecies = 1, nspecies
          pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
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
        end subroutine initialize_goverlap


! ===========================================================================
! Goverlap
! ===========================================================================
! Program Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the site of one of the orbitals.  The potential V(1) is something like Vxc
! for the exchange correlation potential, Vna for the neutral atom potential,
! or 1 for the overlap term.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! This subroutine then writes the results to file.
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        subroutine goverlap_L
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
        integer nFdata_cell_1c              !< indexing of interactions

        real d                              !< distance between the two centers
        real rcutoff1              !< cutoffs for one centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 30) filename

! Procedure
! ============================================================================
        write (ilogfile,*)
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*) '        G O V E R L A P   I N T E R A C T I O N S        '
        write (ilogfile,*) ' ******************************************************* '

! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

! We are doing only Harris here, so set isorp = 0
        isorp = 0

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
          nFdata_cell_1c = pFdata_bundle%nFdata_cell_2c
          pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_1c)

          call make_munuS (nFdata_cell_1c, ispecies, ispecies)
          nME1c_max = pFdata_cell%nME

          allocate (pFdata_cell%fofx(nME1c_max))

          ! Open ouput file for this species pair and 1st = answer1 in old code
          write (filename, '("/goverlap_L.",i2.2,"answer1",".dat")')                   &
     &           species(ispecies)%nZ
          open (unit = 11, file = trim(Fdata_location)//trim(filename),      &
     &          status = 'unknown')

          ! Open ouput file for this species pair and 2nd= answer2 in old code
          write (filename, '("/goverlap_L.",i2.2,"answer2",".dat")')                   &
     &           species(ispecies)%nZ
          open (unit = 12, file = trim(Fdata_location)//trim(filename),      &
     &          status = 'unknown')

          ! Set up grid loop control constants
          rcutoff1 = species(ispecies)%rcutoffA_max

          ! Set final integration limit
          rhomin = 0.0d0
          rhomax = rcutoff1

! Loop over grid
          write (ilogfile,100) species(ispecies)%nZ
          d = 0.0d0

          ! Set integration limits
          zmin = -rcutoff1
          zmax = rcutoff1

          call evaluate_integral_2c (nFdata_cell_1c, ispecies, ispecies,     &
     &                               isorp, ideriv, rcutoff1, rcutoff1,      &
     &                               d, nz_overlap, nrho_overlap,            &
     &                               rint_goverlap_answer1_L, nopi, zmin, zmax,      &
     &                               rhomin, rhomax, pFdata_cell%fofx)

          index_1c = 1
          do issh = 1, species(ispecies)%nssh
            write (11,*) (pFdata_cell%fofx(jssh),                           &
     &                    jssh = index_1c, index_1c + species(ispecies)%nssh - 1)
            index_1c = index_1c + species(ispecies)%nssh
          end do

          call evaluate_integral_2c (nFdata_cell_1c, ispecies, ispecies,     &
     &                               isorp, ideriv, rcutoff1, rcutoff1,      &
     &                               d, nz_overlap, nrho_overlap,            &
     &                               rint_goverlap_answer2_L, nopi, zmin, zmax,      &
     &                               rhomin, rhomax, pFdata_cell%fofx)

          index_1c = 1
          do issh = 1, species(ispecies)%nssh
            write (12,*) (pFdata_cell%fofx(jssh),                           &
     &                    jssh = index_1c, index_1c + species(ispecies)%nssh - 1)
            index_1c = index_1c + species(ispecies)%nssh
          end do
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating goverlap_R integrals for nZ = ', i3,        &
     &              ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine goverlap_L

! ===========================================================================
! Goverlap
! ===========================================================================
! Program Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the site of one of the orbitals.  The potential V(1) is something like Vxc
! for the exchange correlation potential, Vna for the neutral atom potential,
! or 1 for the overlap term.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
! This subroutine then writes the results to file.
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        subroutine goverlap_R
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
        integer ispecies           !< counters for number of species
        integer issh, jssh
        integer index_1c, nME1c_max         !< basically the number of non-zero
        integer isorp, ideriv               !< the number of different types
        integer nFdata_cell_1c              !< indexing of interactions

        real d                              !< distance between the two centers
        real rcutoff1              !< cutoffs for one centers

        real zmin, zmax
        real rhomin, rhomax

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 30) filename

! Procedure
! ============================================================================
        write (ilogfile,*)
        write (ilogfile,*) ' ******************************************************* '
        write (ilogfile,*) '        G O V E R L A P   I N T E R A C T I O N S        '
        write (ilogfile,*) ' ******************************************************* '

! Assign values to the unrequired variables for this specific interaction.
        ideriv = 999

! We are doing only Harris here, so set isorp = 0
        isorp = 0

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          pFdata_bundle=>Fdata_bundle_2c(ispecies, ispecies)
          pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
          nFdata_cell_1c = pFdata_bundle%nFdata_cell_2c
          pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_1c)

          call make_munuS (nFdata_cell_1c, ispecies, ispecies)
          nME1c_max = pFdata_cell%nME

          allocate (pFdata_cell%fofx(nME1c_max))

          ! Open ouput file for this species pair and 1st = answer1 in old code
          write (filename, '("/goverlap_R.",i2.2,"answer1",".dat")')                   &
     &           species(ispecies)%nZ
          open (unit = 13, file = trim(Fdata_location)//trim(filename),      &
     &          status = 'unknown')

          ! Open ouput file for this species pair and 2nd= answer2 in old code
          write (filename, '("/goverlap_R.",i2.2,"answer2",".dat")')                   &
     &           species(ispecies)%nZ
          open (unit = 14, file = trim(Fdata_location)//trim(filename),      &
     &          status = 'unknown')

          ! Set up grid loop control constants
          rcutoff1 = species(ispecies)%rcutoffA_max

          ! Set final integration limit
          rhomin = 0.0d0
          rhomax = rcutoff1

! Loop over grid
          write (ilogfile,100) species(ispecies)%nZ
          d = 0.0d0

          ! Set integration limits
          zmin = -rcutoff1
          zmax = rcutoff1

          call evaluate_integral_2c (nFdata_cell_1c, ispecies, ispecies,     &
     &                               isorp, ideriv, rcutoff1, rcutoff1,      &
     &                               d, nz_overlap, nrho_overlap,            &
     &                               rint_goverlap_answer1_R, nopi, zmin, zmax,      &
     &                               rhomin, rhomax, pFdata_cell%fofx)

          index_1c = 1
          do issh = 1, species(ispecies)%nssh
            write (13,*) (pFdata_cell%fofx(jssh),                           &
     &                    jssh = index_1c, index_1c + species(ispecies)%nssh - 1)
            index_1c = index_1c + species(ispecies)%nssh
          end do

          call evaluate_integral_2c (nFdata_cell_1c, ispecies, ispecies,     &
     &                               isorp, ideriv, rcutoff1, rcutoff1,      &
     &                               d, nz_overlap, nrho_overlap,            &
     &                               rint_goverlap_answer2_R, nopi, zmin, zmax,      &
     &                               rhomin, rhomax, pFdata_cell%fofx)

          index_1c = 1
          do issh = 1, species(ispecies)%nssh
            write (14,*) (pFdata_cell%fofx(jssh),                           &
     &                    jssh = index_1c, index_1c + species(ispecies)%nssh - 1)
            index_1c = index_1c + species(ispecies)%nssh
          end do
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating goverlap_L integrals for nZ = ', i3,        &
     &              ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine goverlap_R

! ===========================================================================
! rint_overlap_answer1_L
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! Here we take the gradient on the left wavefunction.
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        function rint_goverlap_answer1_L (itype, ispecies, jspecies, isorp, d, rho,  &
     &                              z1, z2, ideriv, index_2c)
        implicit none

        real rint_goverlap_answer1_L

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent(in) :: d             !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1
        real rmax                        ! maximum value of r

        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================

! Initialize some dummy variables for warning removal
        idummy = isorp
        idummy = ideriv
        dummy = d
        if (.FALSE.) dummy = z2

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
! GAF dummy modification to use z2
        r1 = sqrt(z1**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n1)%rcutoffA
        psi1val = dpsiofr (r1, rmax, ispecies, n1)

! Find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n2)%rcutoffA
        psi2val = psiofr (r1, rmax, ispecies, n2)
        vofr = 1.0d0

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_goverlap_answer1_L = psi1val*vofr*psi2val*rho


! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_goverlap_answer1_L

! ===========================================================================
! rint_overlap_answer2_L
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! Here we take the gradient on the left wavefunction.
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        function rint_goverlap_answer2_L (itype, ispecies, jspecies, isorp, d, rho,  &
     &                              z1, z2, ideriv, index_2c)
        implicit none

        real rint_goverlap_answer2_L

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent(in) :: d             !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1
        real rmax                        ! maximum value of r

        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================

! Initialize some dummy variables for warning removal
        idummy = isorp
        idummy = ideriv
        dummy = d
        if (.FALSE.) dummy = z2

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
! GAF: Dummy modification to use z2
        r1 = sqrt(z1**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson

! Find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n2)%rcutoffA
        psi2val = psiofr (r1, rmax, ispecies, n2)
        vofr = 1.0d0

! Actual function (Ylm's are calculated above and multilplied after integration
!        rint_goverlap_L = psi1val*vofr*psi2val*rho

! Redefine and add new term to rint_goverlap_L
        rmax = species(ispecies)%shell(n1)%rcutoffA
        psi1val = psiofr (r1, rmax, ispecies, n1)
        vofr = 2.0d0

!       int_goverlap_L = rint_goverlap_L + psi1val*vofr*psi2val/(pi)
        rint_goverlap_answer2_L = psi1val*vofr*psi2val/(pi)

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_goverlap_answer2_L


! ===========================================================================
! ===========================================================================
! rint_overlap_answer1_R
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! Here we take the gradient on the left wavefunction.
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        function rint_goverlap_answer1_R (itype, ispecies, jspecies, isorp, d, rho,  &
     &                              z1, z2, ideriv, index_2c)
        implicit none

        real rint_goverlap_answer1_R

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent(in) :: d             !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1
        real rmax                        ! maximum value of r

        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================

! Initialize some dummy variables for warning removal
        idummy = isorp
        idummy = ideriv
        dummy = d
        if (.FALSE.) dummy = z2

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
        psi2val = dpsiofr (r1, rmax, ispecies, n2)
        vofr = 1.0d0

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_goverlap_answer1_R = psi1val*vofr*psi2val*rho


! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_goverlap_answer1_R


! ===========================================================================
! rint_overlap_answer2_R
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! Here we take the gradient on the left wavefunction.
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        function rint_goverlap_answer2_R (itype, ispecies, jspecies, isorp, d, rho,  &
     &                              z1, z2, ideriv, index_2c)
        implicit none

        real rint_goverlap_answer2_R

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element

        real, intent(in) :: d             !< distance between the two centers
        real, intent(in) :: rho
        real, intent(in) :: z1, z2

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val, psi2val
        real r1
        real rmax                        ! maximum value of r

        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================

! Initialize some dummy variables for warning removal
        idummy = isorp
        idummy = ideriv
        dummy = d
        if (.FALSE.) dummy = z2

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
! GAF: Dummy modification to use z2
        r1 = sqrt(z1**2 + rho**2)

! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.

! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! Find psi1 value at point r1 as based on above and rho called for by
! adaptive_simpson

! Find psi2 value at point r2 as based on above and rho  called for by
! adaptive_simpson
        rmax = species(ispecies)%shell(n2)%rcutoffA
        psi2val = psiofr (r1, rmax, ispecies, n2)
        vofr = 1.0d0

! Actual function (Ylm's are calculated above and multilplied after integration
!        rint_goverlap_L = psi1val*vofr*psi2val*rho

! Redefine and add new term to rint_goverlap_L
        rmax = species(ispecies)%shell(n1)%rcutoffA
        psi1val = psiofr (r1, rmax, ispecies, n1)
        vofr = 2.0d0

!       rint_goverlap_L = rint_goverlap_L + psi1val*vofr*psi2val/(pi)
        rint_goverlap_answer2_R = psi1val*vofr*psi2val/(pi)

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_goverlap_answer2_R


! ===========================================================================
! dpsiofr
! ===========================================================================
! Program Description
! ===========================================================================
!      This function returns the values of the derivative d/dr psiofr(r)
! for the corresponding shell of the atomtype ispec.  The radial functions
! are normalized as:

!  int ( psiofr**2  r**2  dr ) = 1.0

! The wavefunctions for each atom must be read in by calling readpsi
! separately for each atom type.

! The value of r input to this function must be in angstrom units.
! ===========================================================================
! Code written by:
!> @author Amanda Neukirch
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
        function dpsiofr (r, rmax, ispecies, issh)
        implicit none

        real dpsiofr

! Argument Declaration and Description
! ===========================================================================
        real, intent (in) :: r                          ! position
        real, intent (in) :: rmax                       ! max vlaue of r

        integer, intent (in) :: ispecies, issh         ! species and shell

! Local Parameters and Data Declaration
! ===========================================================================
        integer norder
        parameter (norder = 5)

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft, imid
        integer mesh

        real L(0:norder), mu(0:norder), Z(0:norder), alpha(0:norder)
        real a(0:norder), b(0:norder), c(0:norder), d(0:norder)

        integer iam
        integer ix

        real dp
        real h
        real xmin
        real xxp

! Procedure
! ===========================================================================
! note : the points are equally spaced
        h = wf(ispecies)%shell_data(issh)%dr
        imid = int(r/h) + 1
        mesh = wf(ispecies)%shell_data(issh)%mesh

! Special cases
        if (r .ge. rmax) then
          dp = wf(ispecies)%shell_data(issh)%FofR(mesh) - wf(ispecies)%shell_data(issh)%FofR(mesh-1)
          dpsiofr = dp/h
          return
         else if (r .le. 0.0d0) then
          dp = wf(ispecies)%shell_data(issh)%FofR(2) - wf(ispecies)%shell_data(issh)%FofR(1)
          dpsiofr = dp/h
          return
        end if

! Find starting point for the interpolation
        ileft = imid - norder/2
        if (ileft .lt. 1) then
          ileft = 1
        else if (ileft + norder .gt. mesh) then
          ileft = mesh - norder
        end if

! Now interpolate with "natural" splines with f''(x)=0 at end points
        do ix = 0, norder
          a(ix) = wf(ispecies)%shell_data(issh)%FofR(ix + ileft)
        end do

        do ix = 1, norder - 1
          alpha(ix) = 3.0d0*(a(ix+1) - 2.0d0*a(ix) + a(ix-1))/h
        end do

        L(0) = 1
        mu(0) = 0
        Z(0) = 0
        c(0) = 0
        do ix = 1, norder - 1
          L(ix) = (4.0d0 - mu(ix-1))*h
          mu(ix) = h/L(ix)
          Z(ix) = (alpha(ix) - h*Z(ix-1))/L(ix)
        end do
        L(norder) = 1
        mu(norder) = 0
        Z(norder) = 0
        c(norder) = 0

        ! What curve section do we use?
        iam = imid - ileft

        ! Don't need 0 to iam-1
        do ix = norder - 1, iam, -1
          c(ix) = z(ix) - mu(ix)*c(ix+1)
          b(ix) = (a(ix+1) - a(ix))/h - h*(c(ix+1) + 2.0d0*c(ix))/3.0d0
          d(ix) = (c(ix+1) - c(ix))/(3.0d0*h)
        end do

        xmin = 0.0d0
        xxp = (r - (xmin + (imid-1)*h))

!       psiofr = a(iam) + b(iam)*xxp + c(iam)*xxp**2 + d(iam)*xxp**3
        dpsiofr = b(iam) + 2.0d0*c(iam)*xxp + 3.0d0*d(iam)*xxp**2

! End Function
! =============================================================================
        return
        end function dpsiofr


! ===========================================================================
! End Module
! =============================================================================
        end module M_Goverlap
