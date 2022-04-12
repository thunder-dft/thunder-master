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

! M_assemble_ewald
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the interactions for the short range and long range
!! ewald interactions.
!!
!! It contains the following subroutines within the module:
!!
!!       assemble_ewaldsr.f90 - assemble the short-range ewald matrix
!!       assemble_ewaldlr.f90 - assemble the long-range ewald matrix
!!       assemble_ewald.f90 - builds the ewald interactions - atom pair
!!
!! For a complete list of the interactions see the files 2c.Z1.Z2.dir now
!! located in the Fdata directory.  This list will change depending on
!! the datafiles included there. This list is an output from running create.x
! ===========================================================================
        module M_assemble_ewald
        use M_precision
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_2c

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains


! ===========================================================================
! assemble_ewaldsr.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates matrix elements for the short range
!> ewald (Coulomb) interactions. These will be used to offset the total
!> ewald interactions which are calculated in the long-range ewald
!> subroutine.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
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
        subroutine assemble_ewaldsr (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha, iatom, jatom   !< the three parties involved
        integer ibeta, jbeta           !< cells for three atoms
        integer ineigh, mneigh         !< counter over neighbors
        integer in1, in2, in3          !< species numbers
        integer isorp                  !< which interaction and subtype

        integer num_neigh              !< number of neighbors
        integer matom                  !< matom is the self-interaction atom
        integer mbeta                  !< the cell containing neighbor of iatom

        integer norb_mu, norb_nu       !< size of the block for the pair

        real distance_13, distance_23  !< distance from 3rd atom
        real dQ_sum                    !< net charge on atoms
        real z                         !< distance between r1 and r2
        real x                         !< dnabc

        real, dimension (3) :: r1, r2, r3, r12!< positions

        real, dimension (:, :), allocatable :: dterm
        real, dimension (:, :), allocatable :: sterm

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pdip_neighbors
        type(T_assemble_neighbors), pointer :: pdipole_z

        type(T_assemble_block), pointer :: pSR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldsr

! Allocate Arrays
! ===========================================================================
        allocate (s%ewaldsr (s%natoms))
        do iatom = 1, s%natoms
          pewaldsr=>s%ewaldsr(iatom)
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pewaldsr%neighbors(num_neigh))
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          do ineigh = 1, num_neigh   !  <==== loop over i's neighbors
            pSR_neighbors=>pewaldsr%neighbors(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max
            allocate (pSR_neighbors%block(norb_mu, norb_nu))
            pSR_neighbors%block = 0.0d0
          end do
        end do

! Procedure
! ===========================================================================
!
! T W O - C E N T E R   O V E R L A P   M O N O P O L E   P I E C E
!****************************************************************************
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          poverlap=>s%overlap(iatom)
          pewaldsr=>s%ewaldsr(iatom)

          matom = s%neigh_self(iatom)
          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass

          ! cut some more lengthy notation
          pS_neighbors=>poverlap%neighbors(matom)
          pSR_neighbors=>pewaldsr%neighbors(matom)

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)

! Get the overlap matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcnam). Divide by the distance between the centers
! and store this into the short-range ewald piece - ewaldsr.

! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction and we do nothing here.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in ontop case.
            else
              dQ_sum = 0.0d0 
              do isorp = 1, species(in2)%nssh
                dQ_sum = dQ_sum + s%atom(jatom)%shell(isorp)%dQ
              end do
              pSR_neighbors%block = pSR_neighbors%block                      &
     &          + dQ_sum*(pS_neighbors%block/z)*P_eq2
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

!
! T W O - C E N T E R   O V E R L A P   O N T O P   D I P O L E    P I E C E
!****************************************************************************
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          poverlap=>s%overlap(iatom)
          pdipole_z=>s%dipole_z(iatom)
          pewaldsr=>s%ewaldsr(iatom)

          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max

! Loop over the neighbors of each iatom.
          num_neigh = s%neighbors(iatom)%neighn
          do ineigh = 1, num_neigh  ! <==== loop over i's neighbors
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some more lengthy notation
            pS_neighbors=>poverlap%neighbors(ineigh)
            pdip_neighbors=>pdipole_z%neighbors(ineigh)
            pSR_neighbors=>pewaldsr%neighbors(ineigh)

! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)

! Get the overlap matrix from the data files - which is the matrix in molecular
! coordinates (stored in bcnam). Divide by the distance between the centers
! and store this into the short-range ewald piece - ewaldsr.

! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction and we do nothing here.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in ontop case.
            else

! Find sum of charge for first atom and add to ewaldsr
              allocate (sterm (norb_mu, norb_nu))
              allocate (dterm (norb_mu, norb_nu))

              dQ_sum = 0.0d0
              do isorp = 1, species(in1)%nssh
                dQ_sum = dQ_sum + s%atom(iatom)%shell(isorp)%dQ
              end do
              sterm = dQ_sum*pS_neighbors%block/(2.0d0*z)
              dterm = dQ_sum*pdip_neighbors%block/(z**2)

              pSR_neighbors%block = pSR_neighbors%block                      &
     &          + (sterm + dterm)*P_eq2

! Find sum of charge for second atom and add to ewaldsr
              dQ_sum = 0.0d0
              do isorp = 1, species(in2)%nssh
                dQ_sum = dQ_sum + s%atom(jatom)%shell(isorp)%dQ
              end do
              sterm = dQ_sum*pS_neighbors%block/(2.0d0*z)
              dterm = dQ_sum*pdip_neighbors%block/(z**2)

              pSR_neighbors%block = pSR_neighbors%block                      &
     &          + (sterm - dterm)*P_eq2

              deallocate (sterm)
              deallocate (dterm)
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

!****************************************************************************
! T H R E E - C E N T E R   O V E R L A P   A N D    D I P O L E    P I E C E
!****************************************************************************
! Loop over the atoms in the central cell.
        do ialpha = 1, s%natoms
          in3 = s%atom(ialpha)%imass
          r3 = s%atom(ialpha)%ratom

          ! loop over the common neigbor pairs of ialp
          do ineigh = 1, s%neighbors(ialpha)%ncommon
            mneigh = s%neighbors(ialpha)%neigh_common(ineigh)
            if (mneigh .ne. 0) then
              iatom = s%neighbors(ialpha)%iatom_common_j(ineigh)
              ibeta = s%neighbors(ialpha)%iatom_common_b(ineigh)
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              in1 = s%atom(iatom)%imass
              norb_mu = species(in1)%norb_max

              ! cut some lengthy notation
              pS_neighbors=>s%overlap(iatom)%neighbors(mneigh)
              pdip_neighbors=>s%dipole_z(iatom)%neighbors(mneigh)
              pSR_neighbors=>s%ewaldsr(iatom)%neighbors(mneigh)

              jatom = s%neighbors(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%norb_max

! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
              z = distance (r1, r2)

! Find rnabc = vector pointing from center of bondcharge to r3
! This gives us the distance dnabc (or x value in the 2D grid).
              r12 = 0.5d0*(r1 + r2)
              x = distance (r12, r3)

! Find other distances -
              distance_13 = distance (r3, r1)
              distance_23 = distance (r3, r2)

! Now combine the pieces together to get the correct short-range Ewald
! summation. Allocate the size of dterm and sterm according to size
! of the blocks
              allocate (sterm (norb_mu, norb_nu))
              allocate (dterm (norb_mu, norb_nu))

              dQ_sum = 0.0d0 
              do isorp = 1, species(in3)%nssh
                dQ_sum = dQ_sum + s%atom(ialpha)%shell(isorp)%dQ
              end do

              sterm = dQ_sum*pS_neighbors%block/2.0d0
              dterm = dQ_sum*pdip_neighbors%block/z
              pSR_neighbors%block =                                          &
     &          pSR_neighbors%block + ((sterm - dterm)/distance_13           &
     &                               + (sterm + dterm)/distance_23)*P_eq2
              deallocate (sterm)
              deallocate (dterm)
            end if
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_ewaldsr


! ===========================================================================
! assemble_ewaldlr.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the total ewald long-range matrix elements.
!> First, the routine assemble_ewald should be called to get the ewald
!> interactions and this routine then assembles those into matrix elements.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
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
        subroutine assemble_ewaldlr (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh           !< counter over atoms and neighbors
        integer in1, in2                !< species numbers
        integer jatom                   !< neighbor of iatom
        integer issh                    !< counter over shells
        integer num_neigh               !< number of neighbors
        integer mbeta                   !< the cell containing neighbor of iatom

        integer norb_mu, norb_nu        !< size of the block for the pair

        real dQ_sum                     !< net charge on atom
        real z                          !< distance between r1 and r2

        real, dimension (3) :: r1, r2   !< positions of iatom and jatom

        real, allocatable, dimension (:, :) :: dterm
        real, allocatable, dimension (:, :) :: sterm
        real, allocatable, dimension (:) :: sum_ewald

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap
        type(T_assemble_block), pointer :: pdip_neighbors
        type(T_assemble_neighbors), pointer :: pdipole_z

        type(T_assemble_block), pointer :: pLR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldlr

! Allocate Arrays
! ===========================================================================
        allocate (s%ewaldlr (s%natoms))
        allocate (sum_ewald (s%natoms))

! Procedure
! ===========================================================================
        call assemble_ewald (s)

! First a preliminary quantity. ewald(i,j) = sum(L) 1/(| bi-bj +/- L|
! We need to calculate SUMS of ewald sums in which the charge is included.
! In the nutshell we are calculating:
!              Sum_(i,L) q(i)/|b(i)-b(alpha)+L|
        sum_ewald = 0.0d0
        do iatom = 1, s%natoms
          do jatom = 1, s%natoms
            in2 = s%atom(jatom)%imass
            dQ_sum = 0.0d0
            do issh = 1, species(in2)%nssh
              dQ_sum = dQ_sum + s%atom(jatom)%shell(issh)%dQ
            end do
            sum_ewald(iatom) = sum_ewald(iatom) + dQ_sum*s%ewald(iatom,jatom)
          end do
        end do

! Now the meat of the calculation.  Construct ewaldlr(mu,nu,i,m) ===>
! the matrix elements of the long-range parts of the Hamiltonian.
! We make matrix elements for the Long Range Ewald according to our theory:
! ewaldlr(mu,nu,ineigh,iatom) =
! {s(mu,nu,ineigh,iatom)/2}*SUM(j_basis)(Qin(jatom) - Qneutral(jatom))
!                                        *(ewald(iatom,jatom)
!                                          + ewald(ineigh,jatom))*eq2
! The value P_eq2 makes it into the units of eV.
! Loop over the atoms in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          poverlap=>s%overlap(iatom)
          pdipole_z=>s%dipole_z(iatom)
          pewaldlr=>s%ewaldlr(iatom)

          r1 = s%atom(iatom)%ratom
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn
          allocate (pewaldlr%neighbors(num_neigh))

! Loop over the neighbors of the atom i.
          do ineigh = 1, num_neigh
            ! cut some more lengthy notation
            pS_neighbors=>poverlap%neighbors(ineigh)
            pdip_neighbors=>pdipole_z%neighbors(ineigh)
            pLR_neighbors=>pewaldlr%neighbors(ineigh)

            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
            in2 = s%atom(jatom)%imass

! Allocate block size
            norb_nu = species(in2)%norb_max
            allocate (pLR_neighbors%block(norb_mu, norb_nu))
            pLR_neighbors%block = 0.0d0

! SET-UP STUFF
! ***************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
            z = distance (r1, r2)

! Allocate the size of dterm and sterm according to size of the blocks
            allocate (sterm (norb_mu, norb_nu))
            allocate (dterm (norb_mu, norb_nu))

! "Charge" on each atom of the bondcharge. We split the charge S and
! dipole p to be S/2-p/d on atom 1 and S/2+p/d on atom 2. If atom 1 is
! equal to atom 2, then drop the p term.
            sterm = pS_neighbors%block/2.0d0

! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction and we do nothing here.
            if (iatom .eq. jatom .and. mbeta .eq. 0) then
              dterm = 0.0d0
            else
              dterm = pdip_neighbors%block/z
            end if

            pLR_neighbors%block =                                            &
     &        pLR_neighbors%block + (sterm - dterm)*sum_ewald(iatom)*P_eq2   &
     &                            + (sterm + dterm)*sum_ewald(jatom)*P_eq2
            deallocate (sterm)
            deallocate (dterm)
          end do ! end loop over neighbors
        end do ! end loop over atoms

! Deallocate Arrays
! ===========================================================================
        deallocate (sum_ewald)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_ewaldlr


! ===========================================================================
! assemble_ewald.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the total ewald energy. These interactions
!> are then used to buld the long-range ewald matrix elements.
!
!       This routine calculates the Ewald sum for a crystal with a given
! basis. This is specially designed for molecules with a given dipole moment.
! See Ihm, Zunger, Cohen -- Momentum Space Formalism for Total Energy of Solids
! J. Phys C v.12 (79). The ewald sum calculated here is actually gamma(ewald)
! in paper. The terms are also found in M.T. Yin and M.l. Cohen, Phys Rev. B26,
! 3259 (1982).
!
! Output:
!       ewald(natoms,natom) = ewald gamma in eV for all the basis atoms.
!       ewald = gamma1 + gamma2 + gamma3 + gamma4 - vself
!
!       gamma1 = first term of eq. 21 (Yin, Cohen paper), sum over g term.
!       gamma2 = erf term of eq. 21, sum over l term
!       gamma3 = delta(s,s') term in eq. 21.
!       gamma4 = ztot1*ztot2 term in eq. 21.
!       dewald = ewald forces
! ===========================================================================
! Code written by:
!> @author James P. Lewis
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
        subroutine assemble_ewald (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_2c.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom           !< the three parties involved
        integer ig1, ig2, ig3          !< counters for the g grid
        integer ig1mx, ig2mx, ig3mx    !< maximum number of the g grid
        integer il1, il2, il3          !< counters for real-space grid
        integer il1mx, il2mx, il3mx    !< maximum number of space grid
        integer logfile                !< which unit to write output

        real argument, kappa, stuff
        real factor
        real g1mag2, g2mag2, g3mag2    !< magnitude's (squared) of g's
        real gmin2, gmax               !< maximum and minimum magnitude of g's
        real gdotb                     !< dot product of k-space and real space
        real r1mag2, r2mag2, r3mag2    !< magnitude's (squared) of r's
        real rmin2, rmax               !< maximum and minimum magnitude of a's
        real z

        real, dimension (3) :: a1vec, a2vec, a3vec
        real, dimension (3) :: cvec    !< resulting vector
        real, dimension (3) :: g, g1, g2, g3

        interface
          function a_cross_b (a, b)
            real, dimension (3) :: a_cross_b
            real, intent(in), dimension (3) :: a, b
          end function a_cross_b

          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance

          function magnitude (a)
            real magnitude
            real, intent(in), dimension (3) :: a
          end function magnitude
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (s%ewald (s%natoms, s%natoms))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        write (logfile,*) ' Computing the ewald energy.'

! Initialize ewald to zero and initialize other quantities
        s%ewald = 0.0d0
        a1vec = s%lattice(1)%a
        a2vec = s%lattice(2)%a
        a3vec = s%lattice(3)%a

! Determine the reciprical lattice vectors. Do this by Ashcroft and Mermin
! physics.
! First get a2 X a3.
        cvec = a_cross_b (a2vec, a3vec)

! Next find the volume of the cell.
! NOTE: volcel actually has a sign in the above. At this point the sign is
! important since we form g vectors by dividing by a1 dot (a2 X a3).
! Oh, you say. what difference does it make if we change the sign of g.
! it makes no difference in principle.
        s%volume = a1vec(1)*cvec(1) + a1vec(2)*cvec(2) + a1vec(3)*cvec(3)
        g1(:) = 2.0d0*pi*cvec(:)/s%volume

! Next we get a3 X a1, and g2.
        cvec = a_cross_b (a3vec, a1vec)
        g2(:) = 2.0d0*pi*cvec(:)/s%volume

! Finally we get a1 X a2, and g3.
        cvec = a_cross_b (a1vec, a2vec)
        g3(:) = 2.0d0*pi*cvec(:)/s%volume
        s%volume = abs(s%volume)

! Initialize gmax. This determines how far we sum g1, g2, and g3. See below
! why gmax = 5.0d0 is a reasonable criterion.
        gmax = 5.0d0

! Initialize rmax. This determines how far we sum a1, a2, and a3. See below
! why rmax = 5.0d0 is a reasonable criterion. The parameters a1, a2, a3 are
! the direct lattice vectors.
        rmax = 5.0d0

! Determine the magnitude of the vectors g1, g2, g3, a1, a2, a3.
        g1mag2 = magnitude(g1)**2
        g2mag2 = magnitude(g2)**2
        g3mag2 = magnitude(g3)**2

        r1mag2 = magnitude(a1vec)**2
        r2mag2 = magnitude(a2vec)**2
        r3mag2 = magnitude(a3vec)**2

! ****************************************************************************
! The parameter kappa is adjustable, chosen to make the sum's converge rapidly.
! The sum over g converges as exp (-g**2/(4*kappa*kappa)), while the
! sum over l converges as exp (-r**2*kappa**2). Lets set the arguments equal
! to determine a reasonable kappa value. We set them equal for the smallest
! g value and the smallest l value.
! First find the smallest rmag.
        rmin2 = r1mag2
        if (r2mag2 .lt. rmin2) rmin2 = r2mag2
        if (r3mag2 .lt. rmin2) rmin2 = r3mag2

! Next find the smallest gmag.
        gmin2 = g1mag2
        if (g2mag2 .lt. gmin2) gmin2 = g2mag2
        if (g3mag2 .lt. gmin2) gmin2 = g3mag2

! Now set rmin2*kappa**2 = gmin22/(4*kappa**2) and solve for kappa.
        kappa = sqrt(sqrt(gmin2/(4.0d0*rmin2)))

! ****************************************************************************
! In gamma1 we must sum over g vectors. The decay is exp(-g**2/4*kappa*kappa).
! We require the exponent for a given direction in g-space to be gmax**2.
! For instance gmax = 5.0, corresponding to an exponent of gmax**2 = 25.0 seems
! to be a reasonable choice. This gives us g = ig1mx*g1 where
! ig1mx**2 g1**2/(4*kappa**2) = gmax**2. Solve for ig1mx, and add 1.0 for
! good measure.
        ig1mx = int(gmax * sqrt(4.0d0*kappa**2/g1mag2) + 1.0d0)

! Now we do the same thing for g2 and g3
        ig2mx = int(gmax * sqrt(4.0d0*kappa**2/g2mag2) + 1.0d0)
        ig3mx = int(gmax * sqrt(4.0d0*kappa**2/g3mag2) + 1.0d0)

        if (ig1mx .le. 1) ig1mx = 2
        if (ig2mx .le. 1) ig2mx = 2
        if (ig3mx .le. 1) ig3mx = 2

! In gamma2 we must sum over l vectors. The asymptotic decay is
! exp(-kappa*kappa*r**2). We require the exponent for a given direction
! in r-space to be rmax**2. For instance rmax = 5.0, corresponding to an
! exponent of rmax**2 = 25.0 seems to be a reasonable choice. This gives us
! r = il1mx*a1 where il1mx**2 r**2 * kappa**2 = rmax**2. Solve for ir1mx, and
! add 1.0 for good measure.
        il1mx = int(rmax * sqrt(1.0d0/(kappa**2*r1mag2)) + 1.0d0)

! Now we do the same thing for r2 and r3
        il2mx = int(rmax * sqrt(1.0d0/(kappa**2*r2mag2)) + 1.0d0)
        il3mx = int(rmax * sqrt(1.0d0/(kappa**2*r3mag2)) + 1.0d0)

        if (il1mx .le. 1) il1mx = 2
        if (il2mx .le. 1) il2mx = 2
        if (il3mx .le. 1) il3mx = 2

! The real answer: now compute gamma ewald.
! ***********************************************************************
! Compute gamma1:
! ***********************************************************************
! Initialize fewald1

! Sum over g vectors.  If we are doing only a cluster, then only the gamma
! point is considered in the sum.
        if (s%icluster .eq. 1) then
          ig1mx = 0
          ig2mx = 0
          ig3mx = 0
        end if
        do ig1 = -ig1mx, ig1mx
          do ig2 = -ig2mx, ig2mx
            do ig3 = -ig3mx, ig3mx

! skip the origin
              if (.not. (ig1 .eq. 0 .and. ig2 .eq. 0 .and. ig3 .eq. 0)) then
                g(:) = ig1*g1(:) + ig2*g2(:) + ig3*g3(:)
                argument = magnitude(g)**2/(4.0d0*kappa**2)

! The variable stuff contains a number of factors, including the exponential
! which is expensive to compute. That is why its outside the iatom,jatom loop.
                stuff = 4.0d0*pi*exp(-argument)/(magnitude(g)**2*s%volume)

! Sum over s and s', the basis indices.
                do iatom = 1, s%natoms
                  do jatom = iatom, s%natoms
                    factor = 1.0d0*stuff
                    if (jatom .eq. iatom) factor = 0.5d0*stuff
! g dot b:
                    gdotb =                                                   &
     &                g(1)*(s%atom(iatom)%ratom(1) - s%atom(jatom)%ratom(1))  &
     &              + g(2)*(s%atom(iatom)%ratom(2) - s%atom(jatom)%ratom(2))  &
     &              + g(3)*(s%atom(iatom)%ratom(3) - s%atom(jatom)%ratom(3))

                    s%ewald(iatom,jatom) = s%ewald(iatom,jatom) + factor*cos(gdotb)
                    s%ewald(jatom,iatom) = s%ewald(jatom,iatom) + factor*cos(gdotb)
                  end do
                end do
              end if
            end do
          end do
        end do

! ***********************************************************************
! Compute gamma2:
! ***********************************************************************
! Initialize fewald2

! If we are doing only a cluster, then only the central cell is considered
! in the sum.
        if (s%icluster .eq. 1) then
          il1mx = 0
          il2mx = 0
          il3mx = 0
          kappa = 0.0d0
        end if

! Now carry out the sum over the cells.
        do il1 = -il1mx, il1mx
          do il2 = -il2mx, il2mx
            do il3 = -il3mx, il3mx

! Sum over atoms iatom and atoms jatom. Note that we sum over jatom .ge. iatom
! which yields an extra factor of two for iatom .ne. jatom.
              do iatom = 1, s%natoms
                do jatom = iatom, s%natoms
                  factor = 1.0d0
                  if (jatom .eq. iatom) factor = 0.5d0

                  cvec = il1*a1vec + il2*a2vec + il3*a3vec
		  		  cvec = cvec + s%atom(iatom)%ratom
                  z = distance (s%atom(jatom)%ratom, cvec)

! skip the infinite self term.
                  if (z .gt. 0.0001d0) then
                    argument = kappa*z

                    s%ewald(iatom,jatom) =                                   &
     &                s%ewald(iatom,jatom) + factor*erfc(argument)/z
                    s%ewald(jatom,iatom) =                                   &
     &                s%ewald(jatom,iatom) + factor*erfc(argument)/z
                  end if
                end do
              end do
            end do
          end do
        end do

! ***********************************************************************
! Compute gamma3:
! ***********************************************************************
! This term should remain constant always - unless the charges as a function
! of r are changing. Also if the parameter kappa changes corresponding to
! the lattice vectors.
! There are no forces.
        do iatom = 1, s%natoms
          s%ewald(iatom,iatom) = s%ewald(iatom,iatom) - 2.0d0*kappa/sqrt(pi)
        end do

! ***********************************************************************
! Compute gamma4:
! ***********************************************************************
! gamma4 is zero!

! ***********************************************************************
! Combine ewald pieces
! ***********************************************************************

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_ewald


! ===========================================================================
! destroy_assemble_ewaldsr
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the assemble_2c_DOGS
!! information.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
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
        subroutine destroy_assemble_ewald (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                             !< counter over atoms
        integer ineigh                            !< counter over neighbors

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          do ineigh = 1, s%neighbors(iatom)%neighn
            deallocate (s%ewaldsr(iatom)%neighbors(ineigh)%block)
            deallocate (s%ewaldlr(iatom)%neighbors(ineigh)%block)
          end do
          deallocate (s%ewaldsr(iatom)%neighbors)
          deallocate (s%ewaldlr(iatom)%neighbors)
        end do
        deallocate (s%ewald)
        deallocate (s%ewaldsr)
        deallocate (s%ewaldlr)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_assemble_ewald

! End Module
! ===========================================================================
        end module M_assemble_ewald
