! copyright info:
!
!                             @Copyright 2012
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

! M_Coulomb.f90
! Program Description
! ============================================================================
!      This is a module calculating the integrals of two centers for the
! Coulomb - Harris interactions only.
!
! ============================================================================
! Code written by:
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ============================================================================
! Module Declaration
! ============================================================================
        module M_Coulomb
        use M_species
        use M_atom_functions
        use M_integrals_2c

        implicit none

! Type Declaration
! ============================================================================
! Type Declaration
! ============================================================================
! None

! module procedures
        contains


! ===========================================================================
! initialize_Coulomb
! ===========================================================================
! Program Description
! ===========================================================================
!       We need to determine how many interactions belong to each nspecies
! bundle pair. This routine just counts how many total interactions contribute
! from the Coulomb interactions to that bundle.  Something like overlap
! only adds 1 interaction, but something like vna needs number of interactions
! based on the number of shells.
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_Coulomb
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
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
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
        end subroutine initialize_Coulomb


! ===========================================================================
! Coulomb
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
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
!
! Subroutine Declaration
! ===========================================================================
        subroutine Coulomb
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

        write (logfile,*)
        write (logfile,*) ' ******************************************************* '
        write (logfile,*) '          C O U L O M B   I N T E R A C T I O N S        '
        write (logfile,*) ' ******************************************************* '

! Assign values to the unrequired variables for this specific interaction.
        isorp = 0
        ideriv = 0

        ! Set initial integration limit
        rhomin = 0.0d0

! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

            call make_munuS (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

            ! Open ouput file for this species pair
            write (filename, '("/coulomb",".",i2.2,".",i2.2,".dat")')        &
     &              species(ispecies)%nZ, species(jspecies)%nZ
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown')

            ! Set up grid loop control constants
            rcutoff1 = species(ispecies)%rcutoffA_max
            rcutoff2 = species(jspecies)%rcutoffA_max
            dmax = na(ispecies)%rcutoffA_max + na(jspecies)%rcutoffA_max
            drr = dmax/float(ndd_coulomb - 1)
            d = -drr

            ! Set integration limits
            zmin = -rcutoff1
            zmax = rcutoff1
            rhomax = rcutoff1

            ! open directory file
            write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')             &
     &        species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &            status = 'unknown', position = 'append')
            write (13,100) pFdata_bundle%nFdata_cell_2c, P_coulomb, isorp,   &
     &                     filename(2:30), pFdata_cell%nME, ndd_coulomb, dmax
            close (unit = 13)

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")') &
     &             P_coulomb, species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                     index_2c = 1, nME2c_max)

! Loop over grid
            write (logfile,200) species(ispecies)%nZ, species(jspecies)%nZ
            do igrid = 1, ndd_coulomb
              d = d + drr
              call evaluate_integral_2c (nFdata_cell_2c, ispecies, jspecies,&
     &                                   isorp, ideriv, rcutoff1, rcutoff2,  &
     &                                   d, nz_coulomb, nrho_coulomb,        &
     &                                   rint_Coulomb, twopi, zmin, zmax,    &
     &                                   rhomin, rhomax, pFdata_cell%fofx)
              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                     &
     &                                       index_2c = 1, nME2c_max)
            end do !igrid
            write (11,*)
          end do ! end loop over jspecies
        end do ! end loop over ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating Coulomb integrals for nZ = ', i3,          &
     &              ' and nZ = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine Coulomb


! ===========================================================================
! rint_Coulomb
! ===========================================================================
! Program Description
! ===========================================================================
! The rho part of the twocenter_overlap with adaptive simpsons routine.
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
!
! Program Declaration
! ===========================================================================
        real function rint_Coulomb (itype, ispecies, jspecies, isorp, d,     &
     &                               rho, z1, z2, ideriv, index_2c)
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
        integer l1, m1, n1               ! quantum  numbers
        integer l2, m2, n2

        real dummy
        real psi1val
        real r1, r2
        real rmax                         ! maximum value of r along grid

        real vofr

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ===========================================================================

! Initialize some dummy variables for warning removal
        idummy = isorp
        idummy = ideriv
        dummy = d

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
        vofr = vnaofr (r2, jspecies, n2)

! Actual function (Ylm's are calculated above and multilplied after integration
        rint_Coulomb = psi1val**2*vofr*rho/(4.0d0*pi)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! =============================================================================
        return
        end function rint_Coulomb

! End Module
! =============================================================================
        end module
