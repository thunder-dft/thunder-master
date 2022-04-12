! copyright info:
!
!                             @Copyright 2009
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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
       module M_integrals_2c
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
        type T_Fdata_cell_2c
          integer nME                   ! number of non-zero matrix elements

          integer, pointer :: mu_2c (:)
          integer, pointer :: nu_2c (:)
          integer, pointer :: mvalue_2c (:)

          integer, pointer :: N_mu (:)
          integer, pointer :: N_nu (:)
          integer, pointer :: L_mu (:)
          integer, pointer :: L_nu (:)
          integer, pointer :: M_mu (:)
          integer, pointer :: M_nu (:)

          real, pointer :: dx (:)                ! distance between data points
          real, pointer :: xmax (:)              ! maximum interaction range

          ! actual Fdata points f(x)
          real, pointer :: fofx (:)
        end type T_Fdata_cell_2c

! Fdata_bundle_2c is the 2-center package of Fdata_cell_2c. It contains all the
! Fdata_cell_2c information of all interaction/subtypes for given species pair.
        type T_Fdata_bundle_2c
          integer index_2c_overlapS                  ! location of overlapS
          integer nFdata_cell_2c                     ! number of Fdata_cell2C

          ! actual fdata
          type (T_Fdata_cell_2c), pointer :: Fdata_cell_2c (:)
        end type T_Fdata_bundle_2c

        type(T_Fdata_bundle_2c), pointer :: Fdata_bundle_2c (:,:)

        real, parameter :: tolerance = 1d0-6  ! tolerance for Adaptive Simpsons

! module procedures
        contains

! ============================================================================
! size_Fdata_2c
! ============================================================================
! Subroutine Description
! ============================================================================
!       This subroutine calculates the basic size of Fdata_2c.In another word,
! it allocates the value of mu and nu which decide the size of Fdata_2c overlap
! matrix elements, or vna matrix elements, or....
! ============================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ============================================================================
        subroutine size_Fdata_2c
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies        ! counters over species
        integer logfile                     !< writing to which unit

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = 21

        write (logfile, *)
        write (logfile, *) ' Sizing two-center integrals: '
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies

            ! cut some lengthy notation
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            write (logfile,100) ispecies, jspecies, pFdata_bundle%nFdata_cell_2c
            allocate (pFdata_bundle%Fdata_cell_2c(pFdata_bundle%nFdata_cell_2c))

! Set this back to zero and then start counting as interactions are computed.
            pFdata_bundle%nFdata_cell_2c = 0
          end do
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' For ispecies = ', i3, ' and jspecies = ', i3,          &
     &              ' the bundle size = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine size_Fdata_2c


! ============================================================================
! make_munu
! ============================================================================
! Subroutine Description
! ============================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)):
!
! num_orb (in1): number of orbitals in atom-type in1
! mu (index,in1,in2): the mu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2
! nu (index,in1,in2): the nu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2

! (on the BOX ( num_orb(in1) x num_orb(in2)))
!
! Atoms 1 and 2 (bondcharge) are along the z-axis; the third atom is in the
! xz-plane. The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
! Subroutine Declaration
! ============================================================================
        subroutine make_munu (itype, ispecies, jspecies)
        implicit none

! Auguments Declaration and Description
! ============================================================================
! Input
        integer, intent(in) :: itype
        integer, intent(in) :: ispecies, jspecies    ! which pair for mu, nu

! Parameters and Data Declaration
! ============================================================================
! None

! Local Variable Declaration and Description
! ============================================================================
        integer index_2c                ! counter for matrix location - mu, nu
        integer issh, jssh              ! index for looping over shells
        integer mvalue

        integer n1, l1                  ! left quantum numbers
        integer n2, l2                  ! right quantum numbers
        integer nME2c_max

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================
! Loop over the pairs of species.  For each species pair, establish what the
! quantum number values for the orbital mu (the left orbital) and nu (the
! right orbital).
          ! cut some lengthy notation
          pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
          nME2c_max = 0

! First, find the maximum number of matrix elements for each species pair.
          do issh = 1, species(ispecies)%nssh
            l1 = species(ispecies)%shell(issh)%lssh
            do jssh = 1, species(jspecies)%nssh
              l2 = species(jspecies)%shell(jssh)%lssh

              do mvalue = -min(l1,l2), min(l1,l2)
                nME2c_max = nME2c_max + 1
              end do
            end do
          end do

! Now allocate the sizes for mu_2c, nu_2c, and the quantum numbers NLM
! for each mu and nu pair.
          pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)
          pFdata_cell%nME = nME2c_max

          allocate (pFdata_cell%mu_2c(nME2c_max))
          allocate (pFdata_cell%nu_2c(nME2c_max))
          allocate (pFdata_cell%mvalue_2c(nME2c_max))

          allocate (pFdata_cell%N_mu(nME2c_max))
          allocate (pFdata_cell%L_mu(nME2c_max))
          allocate (pFdata_cell%M_mu(nME2c_max))

          allocate (pFdata_cell%N_nu(nME2c_max))
          allocate (pFdata_cell%L_nu(nME2c_max))
          allocate (pFdata_cell%M_nu(nME2c_max))

! Set the values for NLM of each mu, nu pair.
          index_2c = 0
          n1 = 0
          do issh = 1, species(ispecies)%nssh
            l1 = species(ispecies)%shell(issh)%lssh
            n1 = n1 + l1 + 1
            n2 = 0
            do jssh = 1, species(jspecies)%nssh
              l2 = species(jspecies)%shell(jssh)%lssh
              n2 = n2 + l2 + 1
              do mvalue = -min(l1,l2), min(l1,l2)
                index_2c = index_2c + 1
                pFdata_cell%mu_2c(index_2c) = n1 + mvalue
                pFdata_cell%nu_2c(index_2c) = n2 + mvalue
                pFdata_cell%mvalue_2c(index_2c) = 0

                pFdata_cell%N_mu(index_2c) = issh
                pFdata_cell%L_mu(index_2c) = l1
                pFdata_cell%M_mu(index_2c) = mvalue

                pFdata_cell%N_nu(index_2c) = jssh
                pFdata_cell%L_nu(index_2c) = l2
                pFdata_cell%M_nu(index_2c) = mvalue
              end do
              n2 = n2 + l2
            end do
            n1 = n1 + l1
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
        end subroutine make_munu


! ============================================================================
! make_munu_atom
! ============================================================================
! Subroutine Description
! ============================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)):
!
! num_orb (in1): number of orbitals in atom-type in1
! mu (index,in1,in2): the mu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2
! nu (index,in1,in2): the nu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2

! (on the BOX ( num_orb(in1) x num_orb(in2)))
!
! Atoms 1 and 2 (bondcharge) are along the z-axis; the third atom is in the
! xz-plane. The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
!
! For the "atom" case, the number of non-zero matrix elements is dependent
! only on one atom. Both wavefunctions are located at the same site. This
! does that by "acknowledging" the jspecies call, but completely ignoring
! it, except for addressing the N_nu, L_nu, M_nu's
!
! Subroutine Declaration
! ============================================================================
        subroutine make_munu_atom (itype, ispecies, jspecies)
        implicit none

! Auguments Declaration and Description
! ============================================================================
! Input
        integer, intent(in) :: itype
        integer, intent(in) :: ispecies, jspecies    ! which pair for mu, nu

! Parameters and Data Declaration
! ============================================================================
! None

! Local Variable Declaration adn Description
! ============================================================================
        integer index_2c                ! counter for matrix location - mu, nu
        integer issh, jssh              ! index for looping over shells
        integer mvalue

        integer n1, l1                  ! left quantum numbers
        integer n2, l2                  ! right quantum numbers
        integer nME2c_max

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================
! Loop over the pairs of species.  For each species pair, establish what the
! quantum number values for the orbital mu (the left orbital) and nu (the
! right orbital).
        ! cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        nME2c_max = 0

! First, find the maximum number of matrix elements for each species pair.
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          do jssh = 1, species(ispecies)%nssh
            l2 = species(ispecies)%shell(jssh)%lssh
            do mvalue = -min(l1,l2), min(l1,l2)
              nME2c_max = nME2c_max + 1
            end do
          end do
        end do

! Now allocate the sizes for mu_2c, nu_2c, and the quantum numbers NLM
! for each mu and nu pair.
! Overlap Interactions:
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)
        pFdata_cell%nME = nME2c_max

        allocate (pFdata_cell%mu_2c(nME2c_max))
        allocate (pFdata_cell%nu_2c(nME2c_max))
        allocate (pFdata_cell%mvalue_2c(nME2c_max))

        allocate (pFdata_cell%N_mu(nME2c_max))
        allocate (pFdata_cell%L_mu(nME2c_max))
        allocate (pFdata_cell%M_mu(nME2c_max))

        allocate (pFdata_cell%N_nu(nME2c_max))
        allocate (pFdata_cell%L_nu(nME2c_max))
        allocate (pFdata_cell%M_nu(nME2c_max))

! Set the values for NLM of each mu, nu pair.
        index_2c = 0
        n1 = 0
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, species(ispecies)%nssh
            l2 = species(ispecies)%shell(jssh)%lssh
            n2 = n2 + l2 + 1

            do mvalue = -min(l1,l2), min(l1,l2)
              index_2c = index_2c + 1
              pFdata_cell%mu_2c(index_2c) = n1 + mvalue
              pFdata_cell%nu_2c(index_2c) = n2 + mvalue
              pFdata_cell%mvalue_2c(index_2c) = 0

              pFdata_cell%N_mu(index_2c) = issh
              pFdata_cell%L_mu(index_2c) = l1
              pFdata_cell%M_mu(index_2c) = mvalue

              pFdata_cell%N_nu(index_2c) = jssh
              pFdata_cell%L_nu(index_2c) = l2
              pFdata_cell%M_nu(index_2c) = mvalue
            end do
            n2 = n2 + l2
          end do
          n1 = n1 + l1
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
        end subroutine make_munu_atom


! ============================================================================
! make_munuS
! ============================================================================
! Subroutine Description
! ============================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)):
!
! num_orb (in1): number of orbitals in atom-type in1
! mu (index,in1,in2): the mu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2
! nu (index,in1,in2): the nu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2

! (on the BOX ( num_orb(in1) x num_orb(in2)))
!
! Atoms 1 and 2 (bondcharge) are along the z-axis; the third atom is in the
! xz-plane. The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
! Subroutine Declaration
! ============================================================================
        subroutine make_munuS (itype, ispecies, jspecies)
        implicit none

! Auguments Declaration and Description
! ============================================================================
! Input
        integer, intent(in) :: itype
        integer, intent(in) :: ispecies, jspecies    ! which pair for mu, nu

! Parameters and Data Declaration
! ============================================================================
! None

! Local Variable Declaration adn Description
! ============================================================================
        integer index_2c                ! counter for matrix location - mu, nu
        integer issh, jssh              ! index for looping over shells
        integer mvalue

        integer n1, l1                  ! left quantum numbers
        integer n2, l2                  ! right quantum numbers
        integer nME2c_max

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================
! Loop over the pairs of species.  For each species pair, establish what the
! quantum number values for the orbital mu (the left orbital) and nu (the
! right orbital).
        ! cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        nME2c_max = 0

! First, find the maximum number of matrix elements for each species pair.
        do issh = 1, species(ispecies)%nssh
          do jssh = 1, species(jspecies)%nssh
            nME2c_max = nME2c_max + 1
          end do
        end do

! Now allocate the sizes for mu_2c, nu_2c, and the quantum numbers NLM
! for each mu and nu pair.
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)
        pFdata_cell%nME = nME2c_max

        allocate (pFdata_cell%mu_2c(nME2c_max))
        allocate (pFdata_cell%nu_2c(nME2c_max))
        allocate (pFdata_cell%mvalue_2c(nME2c_max))

        allocate (pFdata_cell%N_mu(nME2c_max))
        allocate (pFdata_cell%L_mu(nME2c_max))
        allocate (pFdata_cell%M_mu(nME2c_max))

        allocate (pFdata_cell%N_nu(nME2c_max))
        allocate (pFdata_cell%L_nu(nME2c_max))
        allocate (pFdata_cell%M_nu(nME2c_max))

! Set the values for NLM of each mu, nu pair.
        index_2c = 0
        n1 = 0
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          n1 = n1 + 1
          n2 = 0
          do jssh = 1, species(jspecies)%nssh
            l2 = species(jspecies)%shell(jssh)%lssh
            n2 = n2 + 1

            mvalue = 0
            index_2c = index_2c + 1
            pFdata_cell%mu_2c(index_2c) = n1 + mvalue
            pFdata_cell%nu_2c(index_2c) = n2 + mvalue
            pFdata_cell%mvalue_2c(index_2c) = 0

            pFdata_cell%N_mu(index_2c) = issh
            pFdata_cell%L_mu(index_2c) = l1
            pFdata_cell%M_mu(index_2c) = mvalue

            pFdata_cell%N_nu(index_2c) = jssh
            pFdata_cell%L_nu(index_2c) = l2
            pFdata_cell%M_nu(index_2c) = mvalue
          end do
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
        end subroutine make_munuS


! ============================================================================
! make_munuS_atom
! ============================================================================
! Subroutine Description
! ============================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)):
!
! num_orb (in1): number of orbitals in atom-type in1
! mu (index,in1,in2): the mu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2
! nu (index,in1,in2): the nu-position for each matrix-element (index) between
!                     atom-type in1 and atom-type in2

! (on the BOX ( num_orb(in1) x num_orb(in2)))
!
! Atoms 1 and 2 (bondcharge) are along the z-axis; the third atom is in the
! xz-plane. The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy   yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
!
! For the "atom" case, the number of non-zero matrix elements
! is dependent only on one atom. Both wavefunctions are located at
! the same site. This does that by "acknowledging" the jspecies call,
! But comopletely ignoring it, except for addressing the N_Nu L_Nu M_Nu's
!
! Subroutine Declaration
! ============================================================================
        subroutine make_munuS_atom (itype, ispecies, jspecies)
        implicit none

! Auguments Declaration and Description
! ============================================================================
! Input
        integer, intent(in) :: itype
        integer, intent(in) :: ispecies, jspecies    ! which pair for mu, nu

! Parameters and Data Declaration
! ============================================================================
! None

! Local Variable Declaration adn Description
! ============================================================================
        integer index_2c                ! counter for matrix location - mu, nu
        integer issh, jssh              ! index for looping over shells
        integer mvalue

        integer n1, l1                  ! left quantum numbers
        integer n2, l2                  ! right quantum numbers
        integer nME2c_max

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================
! Loop over the pairs of species.  For each species pair, establish what the
! quantum number values for the orbital mu (the left orbital) and nu (the
! right orbital).
        ! cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
        nME2c_max = 0

! First, find the maximum number of matrix elements for each species pair.
        do issh = 1, species(ispecies)%nssh
          do jssh = 1, species(ispecies)%nssh
            nME2c_max = nME2c_max + 1
          end do
        end do

! Now allocate the sizes for mu_2c, nu_2c, and the quantum numbers NLM
! for each mu and nu pair.
! Overlap Interactions:
        pFdata_cell=>pFdata_bundle%Fdata_cell_2c(itype)
        pFdata_cell%nME = nME2c_max

        allocate (pFdata_cell%mu_2c(nME2c_max))
        allocate (pFdata_cell%nu_2c(nME2c_max))
        allocate (pFdata_cell%mvalue_2c(nME2c_max))

        allocate (pFdata_cell%N_mu(nME2c_max))
        allocate (pFdata_cell%L_mu(nME2c_max))
        allocate (pFdata_cell%M_mu(nME2c_max))

        allocate (pFdata_cell%N_nu(nME2c_max))
        allocate (pFdata_cell%L_nu(nME2c_max))
        allocate (pFdata_cell%M_nu(nME2c_max))

! Set the values for NLM of each mu, nu pair.
        index_2c = 0
        n1 = 0
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          n1 = n1 + 1
          n2 = 0
          do jssh = 1, species(ispecies)%nssh
            l2 = species(ispecies)%shell(jssh)%lssh
            n2 = n2 + 1

            mvalue = 0
            index_2c = index_2c + 1
            pFdata_cell%mu_2c(index_2c) = n1 + mvalue
            pFdata_cell%nu_2c(index_2c) = n2 + mvalue
            pFdata_cell%mvalue_2c(index_2c) = 0

            pFdata_cell%N_mu(index_2c) = issh
            pFdata_cell%L_mu(index_2c) = l1
            pFdata_cell%M_mu(index_2c) = mvalue

            pFdata_cell%N_nu(index_2c) = jssh
            pFdata_cell%L_nu(index_2c) = l2
            pFdata_cell%M_nu(index_2c) = mvalue
          end do
        end do

! End Subroutine
! =============================================================================
        return
        end subroutine make_munuS_atom


! ===========================================================================
! evaluate_integral_2c
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
        subroutine evaluate_integral_2c (itype, ispecies, jspecies, isorp,   &
     &                                     ideriv, rcutoff1, rcutoff2, d,    &
     &                                     nz, nrho, rint, phunction, zmin,  &
     &                                     zmax, rhomin, rhomax, Fdata)
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
          dz = ((rcutoff1 + rcutoff2)/2.0d0)/dfloat(nz)
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
            z1 = zmin + dfloat(iz-1)*dz
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
        end subroutine evaluate_integral_2c


! ===========================================================================
! zint
! ===========================================================================
! Program Description
! ===========================================================================
! This program contains, for the overlap integral, the actual function.
! It is written ins such a way as the function is called by a single variable
! only, this means that the Adaptive_Simpson subroutine can integrate it.
! All other required variables are global.
! This is the function in z, which contains the function in rho (as explained
! in twocenter_overlap) adaptive simpson is called within zint to integrate over
! rho (rint) in situ.
!
! ===========================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute,
! Dublin Institute of Techology,
! Dublin 2,
! Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
!
! Program Declaration
! ===========================================================================
        function zint (itype, ispecies, jspecies, isorp, z, ideriv, index_2c,&
      &                d, nrho, rint, rcutoff1, rcutoff2, rhomin, rhomax)
        implicit none

        real zint

! Argument Declaration
! ===========================================================================
        integer, intent (in) :: ispecies, jspecies     ! two centers species
        integer, intent (in) :: itype, isorp, ideriv   ! which interaction
        integer, intent (in) :: index_2c               ! which matrix element
        integer, intent (in) :: nrho          ! number of integration points

        real, intent (in) :: d                ! distance between two centers
        real, intent(in) :: z
        real, intent(in) :: rcutoff1, rcutoff2
        real, intent(in) :: rhomin, rhomax

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration adn Description
! ===========================================================================
! Simpsons' Variables.
       integer irho, nnrho               ! integrate from irho to nnrho

       real xntegral
       real rho, z1, z2                  ! variables of integration

       real drho                         ! interval between rho points
       real temprho

       ! Simpsons Quadrature Variables.
       real, allocatable, dimension (:) :: rhomult

       interface
         function rint (itype, ispecies, jspecies, isorp, d, rho, z1, z2,    &
     &                  ideriv, index_2c)
           real, intent (in) :: rho, z1, z2, d
           integer, intent (in) :: ispecies, jspecies
           integer, intent (in) :: isorp, itype, ideriv
           integer, intent (in) :: index_2c
           real rint
         end function rint
       end interface

! Procedure
! ===========================================================================
! Set parameters for integration
        z1 = z
        z2 = z1 - d

! ---------------------------------------------------------------------------
! Non Adaptive Simpson's Setup
! ---------------------------------------------------------------------------
! Strictly define what the density of the mesh should be.  Make the density of
! the number of points equivalent for all cases. Change the number of points
! to be integrated to be dependent upon the distance between the centers and
! this defined density.
        drho = max(rcutoff1,rcutoff2)/dfloat(nrho)
        nnrho = int((rhomax - rhomin)/drho)
        if (mod(nnrho,2) .eq. 0) nnrho = nnrho + 1; allocate (rhomult(nnrho))

        rhomult(1) = drho/3.0d0; rhomult(nnrho) = drho/3.0d0
        do irho = 2, nnrho - 1, 2
          rhomult(irho) = 4.0d0*drho/3.0d0
        end do
        do irho = 3, nnrho - 2, 2
          rhomult(irho) = 2.0d0*drho/3.0d0
        end do

!----------------------------------------------------------------------------
! Non Adaptive Simpson's Setup
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Actually use the normal Simpson. Set uprho integration
!----------------------------------------------------------------------------
        xntegral = 0.00
        do irho = 1, nnrho
          rho = rhomin + dfloat(irho - 1)*drho
          temprho = rint (itype, ispecies, jspecies, isorp, d, rho, z1, z2,  &
     &                    ideriv, index_2c)
          temprho = temprho * rhomult(irho)
          xntegral = xntegral + temprho
        end do
        zint = xntegral

! Deallocate Arrays
! ===========================================================================
        deallocate (rhomult)

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function zint


! ===========================================================================
! rescaled_psi
! ===========================================================================
! Function Declaration
! ===========================================================================
        function rescaled_psi (l, m, rho, r, z, psi)
        implicit none

        real rescaled_psi

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l
        integer, intent (in) :: m

        real, intent (in) :: r
        real, intent (in) :: rho
        real, intent (in) :: z
        real, intent (in) :: psi

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
!                       >>> THE MAGIC FORMULA <<<

! The i'th matrix element XY_Z - where X,Y = s, p, d, or f and Z = sigma, pi,
! delta, or phi - is:
!
!      faktor(i)*integral[X_Z(1)*Y_Z(2)*R1*R2*rho*drho*dz],
!
! where the 1, 2 mean the orbital is on atom 1,2. The faktor(i) are given
! below (e.g., faktor(5) is for the matrix element pp_pi and is 3/4)
! and the X_Z and Y_Z are:
!
!          -------------------------------------------------------
!          |                Magic formula table                  |
!          -------------------------------------------------------
!          |             s(sigma)  = 1                           |
!          |             p_sigma   = z/r                         |
!          |             p_pi      = rho/r                       |
!          |             d_sigma   = (2*z**2-rho**2)/r**2        |
!          |             d_pi      = rho*z/r**2                  |
!          |             d_delta   = rho**2/r**2                 |
!          |             f_sigma   = z*(2*z**2-3*rho**2)/r**3    |
!          |             f_pi      = rho*(4*z**2-rho**2)/r**3    |
!          |             f_delta   = z*rho**2/r**3               |
!          |             f_phi     = rho**3/r**3                 |
!          -------------------------------------------------------
! *************************************************************************
! Add magic factors based on what type of orbital is involved in the
! integration.
        rescaled_psi = 0

! s-states:
        if (l .eq. 0) then
          rescaled_psi = psi
          return
        end if

! Quick returns
        if (psi .eq. 0.0d0) then
          rescaled_psi = 0.0d0
          return
        end if
        if (r .le. 1.0d-5)then
          rescaled_psi = 0.0d0
          return
        end if

! p-states:
        if (l .eq. 1) then
          if (abs(m) .eq. 1) rescaled_psi = psi*rho/r
          if (m .eq. 0) rescaled_psi = psi*z/r
          return
        end if

! d-states:
        if (l .eq. 2) then
          if (abs(m) .eq. 2) rescaled_psi = psi*rho**2/r**2
          if (abs(m) .eq. 1) rescaled_psi = psi*rho*z/r**2
          if (m .eq. 0) rescaled_psi = psi*(2.0d0*z**2 - rho**2)/r**2
          return
        end if

! f states:
        if (l .eq. 3) then
          if (abs(m) .eq. 3) rescaled_psi = psi*rho**3/r**3
          if (abs(m) .eq. 2) rescaled_psi = psi*rho**2*z/r**3
          if (abs(m) .eq. 1) rescaled_psi = psi*rho*(4.0d0*z**2 - rho**2)/r**3
          if (m .eq. 0) rescaled_psi = psi*z*(2.0d0*z**2 - 3.0d0*rho**2)/r**3
          return
        end if

! These terms for l > 4 are only used in the case of exact exchange.
        if (l .eq. 4) then
          if (abs(m) .eq. 4) rescaled_psi = psi*rho**4/r**4
          if (abs(m) .eq. 3) rescaled_psi = psi*z*rho**3/r**4
          if (abs(m) .eq. 2) rescaled_psi = psi*(6.0d0*z**2 - rho**2)*rho**2/r**4
          if (abs(m) .eq. 1)                                                 &
     &     rescaled_psi = psi*z*(4.0d0*z**2 - 3.0d0*rho**2)*rho/r**4
         if (m .eq. 0)                                                       &
     &     rescaled_psi = psi*(35.0d0*z**4/r**4                              &
     &                   + (3.0d0*rho**2 - 27.0d0*z**2)/r**2)
         return
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function rescaled_psi


! ===========================================================================
! phunction
! ===========================================================================
! Function Declaration
! ===========================================================================
        function phifactor (l1, m1, l2, m2)
        implicit none

        real phifactor

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l1, m1, l2, m2

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        interface
          function clm (l, m)
            integer, intent (in) :: l, m
            real clm
          end function clm
        end interface

! Procedure
! ===========================================================================
! Function to accept value of l and m for two states and return the value of
! the product of the Ylm for those quantum numbers divided by two, known in
! FIREBALL as the phifactor.
        phifactor = clm(l1,m1)*clm(l2,m2)/2.0d0

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function phifactor


! ===========================================================================
! twopi
! ===========================================================================
! Function Declaration
! ===========================================================================
        function twopi (l1, m1, l2, m2)
        implicit none

        real twopi

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l1, m1, l2, m2

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy

! Procedure
! ===========================================================================
! Set up dummy values.
        idummy = l1
        idummy = m1
        idummy = l2
        idummy = m2

! Here the phifactor is just a constant.
        twopi = 2.0d0*pi

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function twopi

! ===========================================================================
! nopi
! ===========================================================================
! Function Declaration
! ===========================================================================
        function nopi (l1, m1, l2, m2)
        implicit none

        real nopi

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l1, m1, l2, m2

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer idummy

! Procedure
! ===========================================================================
! Set up dummy values.
        idummy = l1
        idummy = m1
        idummy = l2
        idummy = m2

! Here the phifactor is just a constant, and that is one.
        nopi = 1.0d0/2.0d0

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function nopi

! End Module
! =============================================================================
        end module M_integrals_2c
