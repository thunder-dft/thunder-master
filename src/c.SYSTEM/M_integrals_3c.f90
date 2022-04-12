! copyright info:
!
!                             @Copyright 2013
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

! M_integrals_3c.f90
! Program Description
! ============================================================================
!      This is a module calculating the integrals of three centers.
!
! ============================================================================
! Code written by:
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
       module M_integrals_3c
       use M_atom_functions
       use M_species

       include '../include/interactions_3c.h'

! Type Declaration
! ============================================================================
! Three-center interactions arrays

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
        type T_Fdata_cell_3c
          integer nME                   ! number of non-zero matrix elements
          integer nx                    ! number of data points

          integer, pointer :: mu_3c (:)
          integer, pointer :: nu_3c (:)
          integer, pointer :: mvalue_3c (:)

          integer, pointer :: N_mu (:)
          integer, pointer :: N_nu (:)
          integer, pointer :: L_mu (:)
          integer, pointer :: L_nu (:)
          integer, pointer :: M_mu (:)
          integer, pointer :: M_nu (:)

          real, pointer :: dx (:)                ! distance between data points
          real, pointer :: xmax (:)              ! maximum interaction range

          ! actual Fdata points f(x)
          real, pointer :: Fdata_3c (:,:)
        end type T_Fdata_cell_3c

! Fdata_bundle_3c is the 3-center package of Fdata_cell_3c. It contains all the
! Fdata_cell_3c information of all interaction/subtypes for given species pair.
        type T_Fdata_bundle_3c
          integer nFdata_cell_3c                     ! number of Fdata_cell3C

          ! actual fdata
          type (T_Fdata_cell_3c), pointer :: Fdata_cell_3c (:)
        end type T_Fdata_bundle_3c

        type(T_Fdata_bundle_3c), pointer :: Fdata_bundle_3c (:,:,:)

        integer, parameter :: nthetamax = 5       ! Gauss-Legendre expansion

! module procedures
        contains

! ============================================================================
! size_Fdata_3c
! ============================================================================
! Subroutine Description
! ============================================================================
! This subroutine calculates the basic size of Fdata_3c.  In another words, it
! allocates the value of mu and nu which decide the size of all the Fdata_3c
! matrix elements.
! ============================================================================
! Code written by:
! Barry Haycock
! FOCAS Institute
! Dublin Institute of Techology
! Dublin 2, Ireland
! bhaycock@dit.ie
! (+353) 1 402 7960
! ===========================================================================
!
! Subroutine Declaration
! ============================================================================
        subroutine size_Fdata_3c
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies, kspecies        ! counters over species
        integer logfile                             !< writing to which unit

        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = 21

        write (logfile, *)
        write (logfile, *) ' Sizing three-center integrals: '
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            do kspecies = 1, nspecies
              ! cut some lengthy notation
              pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
              write (logfile,100) ispecies, jspecies, kspecies, pFdata_bundle%nFdata_cell_3c
              allocate (pFdata_bundle%Fdata_cell_3c(pFdata_bundle%nFdata_cell_3c))

! Set this back to zero and then start counting as interactions are computed.
              pFdata_bundle%nFdata_cell_3c = 0
           end do
         end do
       end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' For ispecies = ', i3, ', jspecies = ', i3,             &
     &              ', and kspecies = ', i3, ' the bundle size = ', i3)

! End Subroutine
! ===========================================================================
        return
        end subroutine size_Fdata_3c


! ============================================================================
! make_munu_3c
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
        subroutine make_munu_3c (itype, ispecies, jspecies, kspecies)
        implicit none

        include '../include/constants.h'
        include '../include/gridsizes.h'


! Auguments Declaration and Description
! ============================================================================
! Input
        integer, intent(in) :: itype
        integer, intent(in) :: ispecies, jspecies, kspecies

! Parameters and Data Declaration
! ============================================================================
! None

! Local Variable Declaration adn Description
! ============================================================================
        integer index_3c                ! counter for matrix location - mu, nu
        integer issh, jssh              ! index for looping over shells
        integer mvalue
        integer nME3c_max

        integer n1, l1                  ! left quantum numbers
        integer n2, l2                  ! right quantum numbers

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================
! Loop over the pairs of species.  For each species pair, establish what the
! quantum number values for the orbital mu (the left orbital) and nu (the
! right orbital).

! cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)

! First, find the maximum number of matrix elements for each species triplet.
! Actually, its dependent on the ispecies-jspecies pair.
        nME3c_max = 0
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          do jssh = 1, species(jspecies)%nssh
            l2 = species(jspecies)%shell(jssh)%lssh
            do mvalue = -min(l1,l2), min(l1,l2)
              nME3c_max = nME3c_max + 1
            end do
          end do
        end do

! Thats the end of the Delta-m = 0 stuff. Next is Delta-m = 1
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          do jssh = 1, species(jspecies)%nssh
            l2 = species(jspecies)%shell(jssh)%lssh

            if (l1 .eq. 0 .and. l2 .ne. 0) nME3c_max = nME3c_max + 1

            if (l1 .eq. 1) then
              nME3c_max = nME3c_max +1
              if (l2 .ne. 0) nME3c_max = nME3c_max + 1
              if (l2 .eq. 2) nME3c_max = nME3c_max + 2
            end if

            if (l1 .eq. 2) then
              nME3c_max = nME3c_max + 1
              if (l2 .ne. 0) nME3c_max = nME3c_max + 3
              if (l2 .eq. 2) nME3c_max = nME3c_max + 2
            end if
          end do
        end do

! That is the end of the Delta-m = 1 stuff. Next is Delta-m = 2
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          do jssh = 1, species(jspecies)%nssh
            l2 = species(jspecies)%shell(jssh)%lssh
            if (l1 .eq. 2) nME3c_max = nME3c_max + 1
            if (l2 .eq. 2) nME3c_max = nME3c_max + 1
          end do
        end do

! Now allocate the sizes for mu_3c, nu_3c, and the quantum numbers NLM
! for each mu and nu pair.
        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(itype)
        pFdata_cell%nME = nME3c_max

        allocate (pFdata_cell%mu_3c(nME3c_max))
        allocate (pFdata_cell%nu_3c(nME3c_max))
        allocate (pFdata_cell%mvalue_3c(nME3c_max))

        allocate (pFdata_cell%N_mu(nME3c_max))
        allocate (pFdata_cell%L_mu(nME3c_max))
        allocate (pFdata_cell%M_mu(nME3c_max))

        allocate (pFdata_cell%N_nu(nME3c_max))
        allocate (pFdata_cell%L_nu(nME3c_max))
        allocate (pFdata_cell%M_nu(nME3c_max))

! Set the values for NLM of each mu, nu pair.
        index_3c = 0
        n1 = 0
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, species(jspecies)%nssh
            l2 = species(jspecies)%shell(jssh)%lssh
            n2 = n2 + l2 + 1
            do mvalue = -min(l1,l2), min(l1,l2)
              index_3c = index_3c + 1
              pFdata_cell%mu_3c(index_3c) = n1 + mvalue
              pFdata_cell%nu_3c(index_3c) = n2 + mvalue
              pFdata_cell%mvalue_3c(index_3c) = 0

              pFdata_cell%N_mu(index_3c) = issh
              pFdata_cell%L_mu(index_3c) = l1
              pFdata_cell%M_mu(index_3c) = mvalue

              pFdata_cell%N_nu(index_3c) = jssh
              pFdata_cell%L_nu(index_3c) = l2
              pFdata_cell%M_nu(index_3c) = mvalue
            end do
            n2 = n2 + l2
          end do
          n1 = n1 + l1
        end do

! That is the end of the Delta-m = 0 stuff. Next is Delta-m = 1
        n1 = 0
        mvalue = 1
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, species(jspecies)%nssh
            l2 = species(jspecies)%shell(jssh)%lssh
            n2 = n2 + l2 + 1

            if (l1 .eq. 0 .and. l2 .ne. 0) then
              index_3c = index_3c +1
              pFdata_cell%mu_3c(index_3c) = n1 + 0
              pFdata_cell%nu_3c(index_3c) = n2 + 1
              pFdata_cell%mvalue_3c(index_3c) = 1

              pFdata_cell%N_mu(index_3c) = issh
              pFdata_cell%L_mu(index_3c) = l1
              pFdata_cell%M_mu(index_3c) = 0

              pFdata_cell%N_nu(index_3c) = jssh
              pFdata_cell%L_nu(index_3c) = l2
              pFdata_cell%M_nu(index_3c) = +1
            end if

            if (l1 .eq. 1) then
              index_3c = index_3c +1
              pFdata_cell%mu_3c(index_3c) = n1 + 1
              pFdata_cell%nu_3c(index_3c) = n2 + 0
              pFdata_cell%mvalue_3c(index_3c) = 1

              pFdata_cell%N_mu(index_3c) = issh
              pFdata_cell%L_mu(index_3c) = l1
              pFdata_cell%M_mu(index_3c) = +1

              pFdata_cell%N_nu(index_3c) = jssh
              pFdata_cell%L_nu(index_3c) = l2
              pFdata_cell%M_nu(index_3c) = 0

              if (l2 .ne. 0) then
                index_3c = index_3c +1
                pFdata_cell%mu_3c(index_3c) = n1 + 0
                pFdata_cell%nu_3c(index_3c) = n2 + 1
                pFdata_cell%mvalue_3c(index_3c) = 1

                pFdata_cell%N_mu(index_3c) = issh
                pFdata_cell%L_mu(index_3c) = l1
                pFdata_cell%M_mu(index_3c) = 0

                pFdata_cell%N_nu(index_3c) = jssh
                pFdata_cell%L_nu(index_3c) = l2
                pFdata_cell%M_nu(index_3c) = +1
              end if

              if (l2 .eq. 2) then
                index_3c = index_3c +1
                pFdata_cell%mu_3c(index_3c) = n1 + 1
                pFdata_cell%nu_3c(index_3c) = n2 + 2
                pFdata_cell%mvalue_3c(index_3c) = 1

                pFdata_cell%N_mu(index_3c) = issh
                pFdata_cell%L_mu(index_3c) = l1
                pFdata_cell%M_mu(index_3c) = +1

                pFdata_cell%N_nu(index_3c) = jssh
                pFdata_cell%L_nu(index_3c) = l2
                pFdata_cell%M_nu(index_3c) = +2

                index_3c = index_3c +1
                pFdata_cell%mu_3c(index_3c) = n1 - 1
                pFdata_cell%nu_3c(index_3c) = n2 - 2
                pFdata_cell%mvalue_3c(index_3c) = 1

                pFdata_cell%N_mu(index_3c) = issh
                pFdata_cell%L_mu(index_3c) = l1
                pFdata_cell%M_mu(index_3c) = -1

                pFdata_cell%N_nu(index_3c) = jssh
                pFdata_cell%L_nu(index_3c) = l2
                pFdata_cell%M_nu(index_3c) = -2
              end if
            end if

            if (l1 .eq. 2) then
              index_3c = index_3c +1
              pFdata_cell%mu_3c(index_3c) = n1 + 1
              pFdata_cell%nu_3c(index_3c) = n2 + 0
              pFdata_cell%mvalue_3c(index_3c) = 1

              pFdata_cell%N_mu(index_3c) = issh
              pFdata_cell%L_mu(index_3c) = l1
              pFdata_cell%M_mu(index_3c) = +1

              pFdata_cell%N_nu(index_3c) = jssh
              pFdata_cell%L_nu(index_3c) = l2
              pFdata_cell%M_nu(index_3c) = 0

              if (l2 .ne. 0) then
                index_3c = index_3c +1
                pFdata_cell%mu_3c(index_3c) = n1 + 0
                pFdata_cell%nu_3c(index_3c) = n2 + 1
                pFdata_cell%mvalue_3c(index_3c) = 1

                pFdata_cell%N_mu(index_3c) = issh
                pFdata_cell%L_mu(index_3c) = l1
                pFdata_cell%M_mu(index_3c) = 0

                pFdata_cell%N_nu(index_3c) = jssh
                pFdata_cell%L_nu(index_3c) = l2
                pFdata_cell%M_nu(index_3c) = +1

                index_3c = index_3c +1
                pFdata_cell%mu_3c(index_3c) = n1 - 2
                pFdata_cell%nu_3c(index_3c) = n2 - 1
                pFdata_cell%mvalue_3c(index_3c) = 1

                pFdata_cell%N_mu(index_3c) = issh
                pFdata_cell%L_mu(index_3c) = l1
                pFdata_cell%M_mu(index_3c) = -2

                pFdata_cell%N_nu(index_3c) = jssh
                pFdata_cell%L_nu(index_3c) = l2
                pFdata_cell%M_nu(index_3c) = -1

                index_3c = index_3c +1
                pFdata_cell%mu_3c(index_3c) = n1 + 2
                pFdata_cell%nu_3c(index_3c) = n2 + 1
                pFdata_cell%mvalue_3c(index_3c) = 1

                pFdata_cell%N_mu(index_3c) = issh
                pFdata_cell%L_mu(index_3c) = l1
                pFdata_cell%M_mu(index_3c) = +2

                pFdata_cell%N_nu(index_3c) = jssh
                pFdata_cell%L_nu(index_3c) = l2
                pFdata_cell%M_nu(index_3c) = +1
              end if

              if (l2 .eq. 2) then
                index_3c = index_3c +1
                pFdata_cell%mu_3c(index_3c) = n1 + 1
                pFdata_cell%nu_3c(index_3c) = n2 + 2
                pFdata_cell%mvalue_3c(index_3c) = 1

                pFdata_cell%N_mu(index_3c) = issh
                pFdata_cell%L_mu(index_3c) = l1
                pFdata_cell%M_mu(index_3c) = +1

                pFdata_cell%N_nu(index_3c) = jssh
                pFdata_cell%L_nu(index_3c) = l2
                pFdata_cell%M_nu(index_3c) = +2

                index_3c = index_3c +1
                pFdata_cell%mu_3c(index_3c) = n1 - 1
                pFdata_cell%nu_3c(index_3c) = n2 - 2
                pFdata_cell%mvalue_3c(index_3c) = 1

                pFdata_cell%N_mu(index_3c) = issh
                pFdata_cell%L_mu(index_3c) = l1
                pFdata_cell%M_mu(index_3c) = -1

                pFdata_cell%N_nu(index_3c) = jssh
                pFdata_cell%L_nu(index_3c) = l2
                pFdata_cell%M_nu(index_3c) = -2
              end if
            end if
            n2 = n2 + l2
          end do
          n1 = n1 + l1
        end do

! Thats the end of the Delta-m = 1 stuff. Next is Delta-m = 2
        n1 = 0
        mvalue = 2
        do issh = 1, species(ispecies)%nssh
          l1 = species(ispecies)%shell(issh)%lssh
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, species(jspecies)%nssh
            l2 = species(jspecies)%shell(jssh)%lssh
            n2 = n2 + l2 + 1

            if (l1 .eq. 2) then
              index_3c = index_3c +1
              pFdata_cell%mu_3c(index_3c) = n1 + 2
              pFdata_cell%nu_3c(index_3c) = n2 + 0
              pFdata_cell%mvalue_3c(index_3c) = 2

              pFdata_cell%N_mu(index_3c) = issh
              pFdata_cell%L_mu(index_3c) = l1
              pFdata_cell%M_mu(index_3c) = +2

              pFdata_cell%N_nu(index_3c) = jssh
              pFdata_cell%L_nu(index_3c) = l2
              pFdata_cell%M_nu(index_3c) = 0
            end if

            if (l2 .eq. 2) then
              index_3c = index_3c +1
              pFdata_cell%mu_3c(index_3c) = n1 + 0
              pFdata_cell%nu_3c(index_3c) = n2 + 2
              pFdata_cell%mvalue_3c(index_3c) = 2

              pFdata_cell%N_mu(index_3c) = issh
              pFdata_cell%L_mu(index_3c) = l1
              pFdata_cell%M_mu(index_3c) = 0

              pFdata_cell%N_nu(index_3c) = jssh
              pFdata_cell%L_nu(index_3c) = l2
              pFdata_cell%M_nu(index_3c) = +2
            end if

            n2 = n2 + l2
          end do
          n1 = n1 + l1
        end do
        nME3c_max = index_3c
        pFdata_cell%nME = nME3c_max

! End Subroutine
! =============================================================================
        return
        end subroutine make_munu_3c


! ============================================================================
! make_munuS_3c
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

! This is for the spherical 3C case, we ony look at the s shell (l=0), this
! is a modified version of the munu above for this case.

! Subroutine Declaration
! ============================================================================
        subroutine make_munuS_3c (itype, ispecies, jspecies, kspecies)
        implicit none

! Auguments Declaration and Description
! ============================================================================
! Input
        integer, intent(in) :: itype
        integer, intent(in) :: ispecies, jspecies, kspecies

! Parameters and Data Declaration
! ============================================================================
! None

! Local Variable Declaration adn Description
! ============================================================================
        integer index_3c                ! counter for matrix location - mu, nu
        integer issh, jssh              ! index for looping over shells
        integer mvalue

        integer nME3c_max

        integer n1, l1                  ! left quantum numbers
        integer n2, l2                  ! right quantum numbers

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Allocate Arrays
! ============================================================================
! None

! Procedure
! ============================================================================
! Loop over the pairs of species.  For each species pair, establish what the
! quantum number values for the orbital mu (the left orbital) and nu (the
! right orbital).

! cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
        nME3c_max = 0

! First, find the maximum number of matrix elements for each species triplet.
! Actually, its dependent on the ispecies-jspecies pair.
        do issh = 1, species(ispecies)%nssh
          l1 = 0
          do jssh = 1, species(jspecies)%nssh
            l2 = 0
            do mvalue = -min(l1,l2), min(l1,l2)
              nME3c_max = nME3c_max + 1
            end do
          end do
        end do

! Now allocate the sizes for mu_3c, nu_3c, and the quantum numbers NLM
! for each mu and nu pair.
        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(itype)
        pFdata_cell%nME = nME3c_max

        allocate (pFdata_cell%mu_3c(nME3c_max))
        allocate (pFdata_cell%nu_3c(nME3c_max))
        allocate (pFdata_cell%mvalue_3c(nME3c_max))

        allocate (pFdata_cell%N_mu(nME3c_max))
        allocate (pFdata_cell%L_mu(nME3c_max))
        allocate (pFdata_cell%M_mu(nME3c_max))

        allocate (pFdata_cell%N_nu(nME3c_max))
        allocate (pFdata_cell%L_nu(nME3c_max))
        allocate (pFdata_cell%M_nu(nME3c_max))

! Set the values for NLM of each mu, nu pair.
! Now we do delta-m = 0 stuff.
        index_3c = 0
        n1 = 0
        do issh = 1, species(ispecies)%nssh
          l1 = 0
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, species(jspecies)%nssh
            l2 = 0
            n2 = n2 + l2 + 1

            mvalue = 0
            index_3c = index_3c + 1
            pFdata_cell%mu_3c(index_3c) = n1 + mvalue
            pFdata_cell%nu_3c(index_3c) = n2 + mvalue
            pFdata_cell%mvalue_3c(index_3c) = 0

            pFdata_cell%N_mu(index_3c) = issh
            pFdata_cell%L_mu(index_3c) = l1
            pFdata_cell%M_mu(index_3c) = mvalue

            pFdata_cell%N_nu(index_3c) = jssh
            pFdata_cell%L_nu(index_3c) = l2
            pFdata_cell%M_nu(index_3c) = mvalue

            n2 = n2 + l2
          end do
          n1 = n1 + l1
        end do

! Thats the end of the delta-m = 0 stuff. Now Delta-m = 1
        n1 = 0
        mvalue = 1
        do issh = 1, species(ispecies)%nssh
          l1 = 0
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, species(jspecies)%nssh
            l2 = 0
            n2 = n2 + l2 + 1
          end do
        end do

! Thats the end of the Delta-m = 1 stuff. Next is Delta-m = 2
        n1 = 0
        mvalue = 2
        do issh = 1, species(ispecies)%nssh
          l1 = 0
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, species(jspecies)%nssh
            l2 = 0
            n2 = n2 + l2 + 1
          end do
        end do

! End Subroutine
! =============================================================================
        return
        end subroutine make_munuS_3c


! ============================================================================
! gleg
! ============================================================================
! Subroutine Description
! ============================================================================
!       Calculates the roots of the Legendre polynomials.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
        subroutine gleg (ctheta, ctheta_weights, ntheta_max)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ntheta_max

! Output
        real, intent (out), dimension (ntheta_max) :: ctheta
        real, intent (out), dimension (ntheta_max) :: ctheta_weights

! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iterm
        integer itheta
        integer jterm
        integer nw
        integer ndown, nup
        integer node

        real pl1
        real pl2
        real pl3
        real plp
        real pdown, pup
        real xx

! Procedure
! ===========================================================================
! Calculate the positions and weights for an n point Gauss-Legendre
! integration formula. The positions are returned in x and the weights in w.
! The integration interval is from -1 to 1.
        ctheta = 0.0d0
        ctheta_weights = 0.0d0

! Use symmetry -- only need 1/2 the points
! Find the nodes of P sub n.
        xx = 1.0d0
        nw = (ntheta_max + 1)/2
        pl2 = 2.0d0
        plp = 1.0d0

        ndown = nw
        nup = 0
        do iterm = 1, nw

! The variables pup and pdown are upper and lower bounds for the nodes.
          pup = xx
          pdown = 0.0d0
          do jterm = 1, 50

! Use either binary chop or Newton-Raphson to get to node. Use the property
! of the recursion relation that the number of sign changes tells how many
! nodes are located between xx and 1.
            if (nup .eq. iterm - 1 .and. ndown .eq. iterm) then
              xx = xx - pl2/plp
              if (xx .lt. pdown .or. xx .gt. pup) xx = 0.5d0*(pup + pdown)
            else
              xx = 0.5d0*(pup + pdown)
            end if

! Calculate Legendre polynomial from recursion relation and record the number
! of sign changes.
            pl1 = 1.0d0
            pl2 = xx
            node = 0
            do itheta = 2, ntheta_max
              pl3 = (2.0d0*real(itheta) - 1.0d0)*xx*pl2                      &
      &             - (real(itheta) - 1.0d0)*pl1
              if (sign(1.00,pl3) .ne. sign(1.0,pl2)) node = node + 1
              pl1 = pl2
              pl2 = pl3/real(itheta)
            end do

! Calculate the derivative
            plp = ntheta_max*(pl1 - xx*pl2)/(1.0d0 - xx*xx)
            if (node .ge. iterm) then
              pdown = xx
              ndown = node
            else
              pup = xx
              nup = node
            end if
            if (abs(pl2) .lt. 1.0d-10) exit
          end do
          if (abs(pl2) .gt. 1.0d-10)                                         &
     &     write (*,*) ' Warning no convergence in gleg after', jterm,       &
     &                 ' iterations'
          ctheta(iterm) = xx

! Gauss-Legendre weights
          ctheta_weights(iterm) = 2.0d0/((1.0d0 - xx*xx)*plp*plp)
        end do

! Fill in other half of points and weights from symmetry
        do iterm = 1, nw
          ctheta(ntheta_max - iterm + 1) = ctheta(iterm)
          ctheta(iterm) = - ctheta(iterm)
          ctheta_weights(ntheta_max - iterm + 1) = ctheta_weights(iterm)
        end do
        if (int(nw/2)*2 .ne. nw) ctheta(nw) = 0.0d0 ! Symmetry

! Format Statements
! ===========================================================================
        return
        end subroutine gleg


! ===========================================================================
! evaluate_integral_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!      Here we integrate the three-center integral for all non-zero orbital
! combinations. First, we integrate in the "r" direction which is taken from
! the first atom, second, phi, and third theta.
!
!      The integrals are a Legendre polynomial expansion in terms of
! multipoles. We take only the first 5 terms in the expansion. Hence, there
! are 5 special angles that are used so that these 5 terms solve the exact
! solution with their coefficients.
! ===========================================================================
        subroutine evaluate_integral_3c (itype, ispecies, jspecies, kspecies,&
     &                                   ispmin, ispmax, ctheta,             &
     &                                   ctheta_weights, dbcx, dnax, nnr,    &
     &                                   nntheta, psifunction, phiint, qpl)
        implicit none

        include '../include/gridsizes.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itype
        integer, intent (in) :: ispecies, jspecies, kspecies
        integer, intent (in) :: ispmin, ispmax

        integer, intent (in) :: nnr, nntheta  ! mesh for integration

        real, intent (in) :: dbcx, dnax

        real, dimension (P_ntheta), intent(in) :: ctheta
        real, dimension (P_ntheta), intent(in) :: ctheta_weights

! Output
        real, intent (out) :: qpl (:, :, ispmin:)

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer index_3c, nME3c_max     ! different mu, nu types
        integer isorp                   ! loop over shells
        integer itheta, jtheta          ! loop over special angles
        integer l1, l2                  ! l quantum numbers

        integer, allocatable :: nabs (:)
        integer, allocatable :: m1 (:), m2 (:)   ! m quantum numbers

        real cost, sint
        real pl, plm, plmm              ! Legendre polynomials
        real rcutoff1, rcutoff2, rcutoff3   ! atom cutoffs
        real rna (3)                    ! location of "potential"
        real temp

        real, dimension (P_ntheta) :: answer

! Final results after sin(theta) or cos(theta) factor.
        real, allocatable :: ggstore (:, :, :)

        real, allocatable :: gmat (:, :)

        real, allocatable :: prod (:)
        real, allocatable :: znormMU (:), znormNU (:)

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

        interface
          function psifunction (r, rmax, ispecies, n) ! psiofr
            integer, intent (in) :: n, ispecies
            real, intent (in) :: r
            real, intent (in) :: rmax
            real psifunction
          end function psifunction
        end interface

        interface
          subroutine phiint (itype, ispecies, jspecies, kspecies, ispmin,    &
                             ispmax, r, ds, zr, rna, avgVmat)
            integer, intent (in) :: ispecies, jspecies, kspecies,           &
     &                              ispmin, ispmax, itype
            real, intent (in) :: ds, r, zr
            real, intent (in) :: rna (3)
            real, intent (out) :: avgVmat(ispmin:, :)
          end subroutine phiint
        end interface

! Procedure
! ===========================================================================
! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(itype)
        nME3c_max = pFdata_cell%nME

! allocate the quantum numbers and set values
        allocate (znormMu (nME3c_max))
        allocate (znormNu (nME3c_max))

! Initialize some variables
        rcutoff1 = species(ispecies)%rcutoffA_max
        rcutoff2 = species(jspecies)%rcutoffA_max
        rcutoff3 = species(kspecies)%rcutoffA_max

! ***************************************************************************
! Add in integration factors.
! ***************************************************************************
! Here is the correct list:
!   m:     -2       -1         0        1         2
!          xy       yz      3z^2-r^2   xz       x^2-y^2
! sq15 * (  1        1        1/sq12    1       1/sq4  )
!
! ***************************************************************************
! Perform a test on nabs here
        allocate (m1 (nME3c_max))
        allocate (m2 (nME3c_max))
        m1 = pFdata_cell%M_mu(:)
        m2 = pFdata_cell%M_nu(:)

        do index_3c = 1, nME3c_max
          l1 = pFdata_cell%L_mu(index_3c)
          l2 = pFdata_cell%L_nu(index_3c)
          if (l1 .eq. 0) znormMu(index_3c) = 1.0d0
          if (l2 .eq. 0) znormNu(index_3c) = 1.0d0
          if (l1 .eq. 1) znormMu(index_3c) = sqrt(3.0d0)
          if (l2 .eq. 1) znormNu(index_3c) = sqrt(3.0d0)
          if (l1 .eq. 2) znormMu(index_3c) = sqrt(15.0d0)
          if (l2 .eq. 2) znormNu(index_3c) = sqrt(15.0d0)

! Working on d for mu
          if (l1 .eq. 2)then
            if (m1(index_3c) .eq. 0)                                         &
     &        znormMu(index_3c) = znormMu(index_3c)/sqrt(12.0d0)
            if (m1(index_3c) .eq. 2)                                         &
     &        znormMu(index_3c) = znormMu(index_3c)/sqrt(4.0d0)
          end if

! Working on d for nu
          if (l2 .eq. 2) then
            if (m2(index_3c) .eq. 0)                                         &
     &        znormNu(index_3c) = znormNu(index_3c)/sqrt(12.0d0)
            if (m2(index_3c) .eq. 2)                                         &
     &        znormNu(index_3c) = znormNu(index_3c)/sqrt(4.0d0)
          end if
        end do ! end loop over index_3c

        allocate (nabs(nME3c_max))
        nabs = abs(m1 - m2)
        if (any(nabs .gt. 2)) stop 'WRONG NABS IN BCNA!!!!!'

! ***************************************************************************
!              Do integral over r:
! ***************************************************************************
! ---------------------------------------------------------------------------
! Since threecenter_integral internaly loops over ispmin to ispmax.
! The thetas are roots of P_(ntheta_max + 1)
! ---------------------------------------------------------------------------
        allocate (gmat(ispmin:(ispmax - ispmin + 1), nMe3c_max)); gmat = 0.0d0
        allocate (ggstore(ispmin:(ispmax - ispmin + 1), nMe3c_max, P_ntheta))
        ggstore = 0.0d0
        do itheta = 1, P_ntheta
          cost = ctheta(itheta)
          sint = sqrt(1 - cost*cost)

          rna(1) = sint*dnax
          rna(2) = 0.0d0
          rna(3) = cost*dnax

          gmat = 0.00
          call Rint (itype, ispecies, jspecies, kspecies, ispmin, ispmax,    &
                     dbcx, rna, nnr, nntheta, psifunction, phiint, gmat)

! Finally,the normalization factors for the different orbitals. For instance,
! a p orbital is sqrt(3) * x/r * r(r). The sqrt(3) factor (and ones like it)
! are now included. Also averagetheta*averagephi = 1/(2 pi) is added.
          allocate (prod (nME3c_max)); prod = 0.0d0
          prod = znormMu*znormNu*0.5d0/(4.0d0*atan(1.0d0))
          do isorp = ispmin, ispmax
            gmat(isorp,:) = gmat(isorp,:)*prod
          end do

! ***************************************************************************
! SUMMARY
! ***************************************************************************
! We have computed gmat(isorp,mu,nu). mu, and nu are 1 to 9 for sp^3d^5.

! We are in molecular coordinates with z along sigma, x along pi, and
! y along pi'. Many matrix elements og gmat are zero. We have computed them
! anyway, and indeed find they are zero. Just to avoid at a later time,
! any trouble with roundoffs, I will now set those that are supposed to be
! zero, identically to zero.

! Here is the matrix, and the zero's are indicated.
!
!          \ s   x   y   z    3z^2-r^2    x^2-y^2    xz    xy      yz
!
!  s                 0                                     0        0
!
!  x                 0                                     0        0
!
!  y         0   0       0        0          0        0
!
!  z                 0                                     0        0
!
! 3z^2-r^2           0                                     0        0
!
! x^2-y^2            0                                     0        0
!
!  xz                0                                     0        0
!
!  xy        0   0       0        0          0        0
!
!  yz        0   0       0        0          0        0
!
! ***************************************************************************
! ---------------------------------------------------------------------------
! Correct integrals as either type A or B
! ---------------------------------------------------------------------------
          do index_3c = 1, nME3c_max
            if (nabs(index_3c) .eq. 1) then
              do isorp = ispmin, ispmax
!               type B: V=sin(theta)*Sum(l)* P*Q (nabs=1)
                if (sint .lt. 0.001d0) stop 'sin theta is zero!!!'
                ggstore(isorp,index_3c,itheta) = gmat(isorp,index_3c)/sint
              end do
            else
!             type A: V=Sum(l) P*Q. (nabs=0,2):do nothing
              do isorp = ispmin, ispmax
                ggstore(isorp,index_3c,itheta) = gmat(isorp,index_3c)
              end do
            end if
          end do
          deallocate (prod)
        end do ! end loop over P_ntheta (above)
        deallocate (gmat)


! ***************************************************************************
! F I N A L   P R O D U C T
! ***************************************************************************
! Now all qpl coefficients are computed for a fixed dbc, dna pair, for all
! potentials  isorp = ispmin, ispmax and for all itheta(...as). The results
! for different matrix elements of one combination of isorp and itheta
! are written sequentially into one file.

! Begin the loop over all possible potentials on the third site (in3).
        do isorp = ispmin, ispmax

! ----------------------------------------------------------------------------
!    Gauss-Legendre integration. Since the Gauss-Legendre integration has to
!    be done separately for each kind of potential and for each kind of orbital
!    combination, we need to loop here once more over the number of
!    non-vanishing three center integrals
! ----------------------------------------------------------------------------
          do index_3c = 1, nME3c_max

! Looping over the thetas, which are the roots of a Pl
            answer = 0.0d0
            do itheta = 1, P_ntheta
              plmm = 1.0d0
              plm = ctheta(itheta)
              temp = ctheta_weights(itheta)*ggstore(isorp,index_3c,itheta)
              if (abs(temp) .lt. 1.0d-10) temp = 0.0d0
              answer(1) = answer(1) + plmm*temp
              answer(2) = answer(2) + plm*temp
              do jtheta = 3, P_ntheta
                pl = (plm*ctheta(itheta)                                     &
     &              *(2.0d0*jtheta - 3.0d0) - (jtheta - 2.0d0)*plmm)         &
     &              /(jtheta - 1.0d0)
                answer(jtheta) = answer(jtheta) + pl*temp
                plmm = plm
                plm = pl
              end do
            end do

! Normalize the coefficient, and write them out to qpl's.
            do itheta = 1, P_ntheta
              qpl(itheta,index_3c,isorp) =                                   &
     &          answer(itheta)*(2.0d0*itheta - 1.d0)*0.5d0
            end do
          end do  ! end of the GL loop
        end do  ! end of the isorp loop

! Deallocate Arrays
! ===========================================================================
        deallocate (ggstore)
        deallocate (m1, m2)
        deallocate (nabs)
        deallocate (znormMu, znormNu)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine evaluate_integral_3c


! ===========================================================================
! Rint
! ===========================================================================
! Program Description
! ===========================================================================
! This program contains the actual R function for integration purposes.
! This is the function in rna (the distance between the two atom centers)
! and dbcx, which is the bond-center distance to the potential.
!
! The result is written into gmat which is are the R integrals stored
! in a mu, nu form.
! ===========================================================================
        subroutine Rint (itype, ispecies, jspecies, kspecies, ispmin, ispmax,&
     &                   dbcx, rna, nnr, nntheta, psifunction, phiint, gmat)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispmin, ispmax, itype
        integer, intent (in) :: ispecies, jspecies, kspecies
        integer, intent (in) :: nnr, nntheta    ! mesh for integration

        real, intent (in) :: dbcx, rna(3)

! Output
        real, intent (out) :: gmat (ispmin:, :)

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ir                      ! integration over theta
        integer isorp                   ! counter over shells
        integer nME3c_max               ! number of matrix elements

        real prod
        real r                          ! radial integration point
        real rcutoff1, rcutoff2         ! atom cutoffs
        real rmin, rmax                 ! minimum and maximum integration

        real, allocatable :: rmult (:)
        real, allocatable :: thetamat (:, :)

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

        interface
          function psifunction (r, rmax, ispecies, n)  ! psiofr
            integer, intent (in) :: n, ispecies
            real, intent (in) :: r
            real, intent (in) :: rmax
            real psifunction
          end function psifunction
        end interface

        interface
          subroutine phiint (itype, ispecies, jspecies, kspecies, ispmin,    &
                             ispmax, r, ds, zr, rna, avgVmat)
            integer, intent (in) :: ispecies, jspecies, kspecies,           &
     &                              ispmin, ispmax, itype
            real, intent (in) :: ds, r, zr
            real, intent (in) :: rna (3)
            real, intent (out) :: avgVmat(ispmin:, :)
          end subroutine phiint
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (rmult (nnr))

! Procedure
! ===========================================================================
! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(itype)
        nME3c_max = pFdata_cell%nME

! Allocate arrays
        allocate (thetamat(ispmin:(ispmax - ispmin + 1), nME3c_max)); thetamat = 0.0d0

! Get the cutoffs for each atom.
        rcutoff1 = species(ispecies)%rcutoffA_max
        rcutoff2 = species(jspecies)%rcutoffA_max

        rmin = 0.0d0
        rmax = 0.5d0*(rcutoff1 + rcutoff2)

! Set up integration factors
        dr = (rmax - rmin)/dfloat(nnr - 1)
        rmult(1) = dr*41.0D0/140.0d0
        do ir = 2, nnr - 1
          if (mod(ir,6) .eq. 2) rmult(ir) = dr*216.0d0/140.0d0
          if (mod(ir,6) .eq. 3) rmult(ir) = dr*27.0d0/140.0d0
          if (mod(ir,6) .eq. 4) rmult(ir) = dr*272.0d0/140.0d0
          if (mod(ir,6) .eq. 5) rmult(ir) = dr*27.0d0/140.0d0
          if (mod(ir,6) .eq. 0) rmult(ir) = dr*216.0d0/140.0d0
          if (mod(ir,6) .eq. 1) rmult(ir) = dr*82.0d0/140.0d0
        end do
        rmult(nnr)= dr*41.0D0/140.0d0

! ***************************************************************************
!              Do integral over r:
! ***************************************************************************
        do ir = 1, nnr
          r = rmin + dfloat(ir - 1)*dr

          thetamat = 0.0d0
          call thetaint (itype, ispecies, jspecies, kspecies, ispmin, ispmax,&
     &                   r, dbcx, rna, nntheta, psifunction, phiint, thetamat)

          prod = rmult(ir)*r*r
          do isorp = ispmin, ispmax
            gmat(isorp,:) = gmat(isorp,:) + prod*thetamat(isorp,:)
          end do
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (thetamat)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Rint


! ===========================================================================
! thetaint
! ===========================================================================
! Program Description
! ===========================================================================
! This program contains the actual theta function for integration purposes.
! This is the function in r (the variable from the center of the first atom)
! and dbcx, which is the bond-center distance.
!
! The result is written into thetamat which is are the theta integrals stored
! in a mu, nu form.
! ===========================================================================
        subroutine thetaint (itype, ispecies, jspecies, kspecies, ispmin,  &
     &                       ispmax, r, dbcx, rna, nntheta, psifunction,   &
     &                       phiint, thetamat)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ispmin, ispmax, itype    ! which type doing
        integer, intent (in) :: ispecies, jspecies, kspecies  ! the species
        integer, intent (in) :: nntheta  ! mesh for integration

        real, intent (in) :: r, dbcx, rna(3)   ! the integration coordinates

        interface
          function psifunction (r, rmax, ispecies, n)  ! psi  change
            integer, intent (in) :: n, ispecies
            real, intent (in) :: r
            real, intent (in) :: rmax
            real psifunction
          end function psifunction
        end interface

        interface
          subroutine phiint (itype, ispecies, jspecies, kspecies, ispmin,    &
                             ispmax, r, ds, zr, rna, avgVmat)
            integer, intent (in) :: ispecies, jspecies, kspecies,           &
     &                              ispmin, ispmax, itype
            real, intent (in) :: ds, r, zr
            real, intent (in) :: rna (3)
            real, intent (out) :: avgVmat(ispmin:, :)
          end subroutine phiint
        end interface

! Output
        real, intent (out) :: thetamat (ispmin:, :)

! Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer isorp                   ! counter over shells
        integer itheta                  ! integration over theta
        integer lmu, mmu, nmu           ! mu quantum numbers
        integer lnu, mnu, nnu           ! nu quantum numbers
        integer nME3c_max               ! number of matrix elements
        integer index_3c                ! counter 1, ..., nME3c_max3c_max

        real ds                         ! dcos(theta), dsin(theta)
        real dtheta                     ! interval between theta points
        real prod                       ! final product for factors
        real r1, r2                     ! radial distances for each center
        real rcutoff1, rcutoff2         ! atom cutoffs
        real stuffmunu
        real theta                      ! value of theta at point of integration
        real thetaL, thetaR             ! theta values for mu, nu wavefunctions
        real zr
        real rmax1, rmax2               ! max values of r

        real, allocatable :: avgVmat (:, :)
        real, allocatable :: thetamult (:)

        type (T_Fdata_cell_3c), pointer :: pFdata_cell
        type (T_Fdata_bundle_3c), pointer :: pFdata_bundle

! Allocate Arrays
! ===========================================================================
        allocate (thetamult (nntheta))

! Procedure
! ===========================================================================
        zr = 0.0d0

! Cut some lengthy notation
        pFdata_bundle=>Fdata_bundle_3c(ispecies, jspecies, kspecies)
        pFdata_cell=>pFdata_bundle%Fdata_cell_3c(itype)
        nME3c_max = pFdata_cell%nME

! Get the cutoffs for each atom.
        rcutoff1 = species(ispecies)%rcutoffA_max
        rcutoff2 = species(jspecies)%rcutoffA_max

! Deallocate arrays if already allocated
        allocate (avgVmat(ispmin:(ispmax - ispmin + 1), nME3c_max)); avgVmat = 0.0d0

! Set up integration factors
        dtheta = 4.0d0*atan(1.0d0)/dfloat(nntheta - 1)
        thetamult(1) = dtheta*41.0D0/140.0d0
        do itheta = 2, nntheta - 1
          if (mod(itheta,6) .eq. 2) thetamult(itheta) = dtheta*216.0d0/140.0d0
          if (mod(itheta,6) .eq. 3) thetamult(itheta) = dtheta*27.0d0/140.0d0
          if (mod(itheta,6) .eq. 4) thetamult(itheta) = dtheta*272.0d0/140.0d0
          if (mod(itheta,6) .eq. 5) thetamult(itheta) = dtheta*27.0d0/140.0d0
          if (mod(itheta,6) .eq. 0) thetamult(itheta) = dtheta*216.0d0/140.0d0
          if (mod(itheta,6) .eq. 1) thetamult(itheta) = dtheta*82.0d0/140.0d0
        end do
        thetamult(nntheta) = dtheta*41.0d0/140.0d0

! ***************************************************************************
!              Do integral over theta:
! ***************************************************************************
        do itheta = 1, nntheta
          theta = dfloat(itheta - 1)*dtheta
          ds = sin(theta)

          r1 = r**2 + 0.25d0*(dbcx**2) + r*dbcx*cos(theta)
          if (r1 .le. 0.0d0) r1 = 0.0d0
          r1 = sqrt(r1)

          r2 = r**2 + 0.25d0*(dbcx**2) - r*dbcx*cos(theta)
          if (r2 .le. 0.0d0) r2 = 0.0d0
          r2 = sqrt(r2)

          ! go to next angle if this theta causes integration point outside
          if (r1 .gt. rcutoff1 .or. r2 .gt. rcutoff2) cycle
          zr = r*cos(theta)

          thetaL = 0.0d0
          if (abs(r1 - (zr + 0.5d0*dbcx)) .le. 1.0d-10 .or.                  &
     &        abs(r1 + (zr + 0.5d0*dbcx)) .le. 1.0d-10) then
            thetaL = 0.0d0
          else if (r1 .gt. 1.0d-4) then
            thetaL = acos((zr + 0.5d0*dbcx)/r1)
          end if

          thetaR = 0.0d0
          if (abs(r2 - (zr - 0.5d0*dbcx)) .le. 1.0d-10 .or.                  &
     &        abs(r2 + (zr - 0.5d0*dbcx)) .le. 1.0d-10) then
            thetaR = 0.0d0
          else if (r2 .gt. 1.0d-4) then
            thetaR = acos((zr - 0.5d0*dbcx)/r2)
          end if

! Do integral over phi:
          call phiint (itype, ispecies, jspecies, kspecies, ispmin, ispmax,  &
                       r, ds, zr, rna, avgVmat)

          prod = thetamult(itheta)*sin(theta)
          do index_3c = 1, nME3c_max
            nmu = pFdata_cell%N_mu(index_3c)
            lmu = pFdata_cell%L_mu(index_3c)
            mmu = pFdata_cell%M_mu(index_3c)
            nnu = pFdata_cell%N_nu(index_3c)
            lnu = pFdata_cell%L_nu(index_3c)
            mnu = pFdata_cell%M_nu(index_3c)
            rmax1 = species(ispecies)%shell(nmu)%rcutoffA
            rmax2 = species(jspecies)%shell(nnu)%rcutoffA

            stuffmunu = prod*thetafactor(lmu, mmu, thetaL)                   &
     &                      *psifunction(r1, rmax1, ispecies, nmu)           &
     &                      *thetafactor(lnu, mnu, thetaR)                   &
     &                      *psifunction(r2, rmax2, jspecies, nnu)

            do isorp = ispmin, ispmax
              thetamat(isorp,index_3c) =                                     &
     &          thetamat(isorp,index_3c) + avgVmat(isorp,index_3c)*stuffmunu
            end do
          end do ! index_3c
        end do ! itheta

! Deallocate Arrays
! ===========================================================================
        deallocate (thetamult)
        deallocate (avgVmat)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine thetaint


! ===========================================================================
! thetafactor
! ===========================================================================
! Function Declaration
! ===========================================================================
        function thetafactor (l, m, theta)
        implicit none

        real thetafactor

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: l, m

        real, intent (in) :: theta

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Theta factors for A and B.
! Use the (l,m) notation. For example 3z^2-1 becomes (2,0),
! px becomes (1,1), and so on.

! s-states
        if (l .eq. 0) then
          thetafactor = 1.0d0
          return
        end if

! p-states:
        if (l .eq. 1) then
          if (abs(m) .eq. 1) thetafactor = sin(theta)
          if (m .eq. 0) thetafactor = cos(theta)
          return
        end if

! d-states:
        if (l .eq. 2) then
          if (abs(m) .eq. 2) thetafactor = sin(theta)**2
          if (abs(m) .eq. 1) thetafactor = sin(theta)*cos(theta)
          if (m .eq. 0) thetafactor = 3.0d0*cos(theta)**2 - 1.0d0
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
        end function thetafactor

! End Module
! =============================================================================
        end module
