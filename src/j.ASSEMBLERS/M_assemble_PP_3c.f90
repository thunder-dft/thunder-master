! copyright info:
!
!                             @Copyright 2016
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
! West Virginia University - Khorgolkhuu Odbadrakh
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman
!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! M_assemble_PP_3c.f90
! Module Description
! ===========================================================================
!>       This is a module containing all of the assembler programs required
!! to assemble all of the matrix elements for the three-center interactions
!! related to the pseudopotential interactions.
!! It contains the following subroutines within the module:
!!
!! assemble_PP_3 - assemble total pseudopotential pieces for three center
!
! ===========================================================================
        module M_assemble_PP_3c
        use M_assemble_blocks
        use M_configuraciones
        use M_Fdata_3c
        use M_rotations

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! assemble_vnl_3c
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions.
! These are true three-center in that iatom .ne. jatom .ne. katom.
! The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
! A third party term is when we consider the NA (or etc.) to be at the origin
! and take the matrix element between a pair of neighbors, neither of which is
! the NA (or etc.), and the pair is truly a pair, and not an atom.
! the results are stored in: f3na(ix,i), f3xc(ix,i), etc.
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms.!
! ===========================================================================
! Code written by:
!> @author Khorgolkhuu Odbadrakh
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine assemble_vnl_3c (s)
        implicit none

        include '../include/constants.h'
        include '../include/interactions_3c.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ialpha, iatom, jatom         !< counter over atoms and neighbors
        integer in1, in2, indna              !< species numbers
        integer num_neigh                    !< number of neighbors
        integer ibeta, jbeta
        integer ineigh, mneigh, m31, m32
        integer ncc

        integer imu, inu
        integer norb_mu, norb_nu         !< size of the block for the pair

        real, pointer :: cl_value (:)
        real, allocatable :: bcnlx (:,:)

        real, dimension (3) :: r1, r2, rna, r31, r32

        interface
          function cl(itype)
            real, pointer :: cl (:)
            integer, intent(in) :: itype
          end function cl
        end interface

! interaction arrays
        type(T_assemble_block), pointer :: pvnl_neighbors
        type(T_assemble_block), pointer :: psvnl1_neighbors
        type(T_assemble_block), pointer :: psvnl2_neighbors

        type(T_assemble_neighbors), pointer :: pvnl
        type(T_assemble_neighbors), pointer :: psvnl1
        type(T_assemble_neighbors), pointer :: psvnl2

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ASSEMBLE VNL ATM CASE  <phi_i|Psi_j><Psi_j|phi_i>
        do ialpha = 1, s%natoms
          rna = s%atom(ialpha)%ratom
          indna = s%atom(ialpha)%imass
          num_neigh = s%neighbors_PP(ialpha)%ncommon

! For the non-local potential pieces call find the coefficients corresponding
! to indna.
          !memory is allocated inside function
          cl_value => cl(indna)

! Loop over the neighbors of each ialp.
! Now look at all common neighbor pairs which were figured out in main.
! The common neighbors were figured out in common_neighbors.f90
          do ineigh = 1, num_neigh
            mneigh = s%neighbors_PP(ialpha)%neigh_common(ineigh)

! The second atom (jatom) is the mneigh'th neighbor of iatom.
            if (mneigh .ne. 0) then
              iatom = s%neighbors_PP(ialpha)%iatom_common_j(ineigh)
              ibeta = s%neighbors_PP(ialpha)%iatom_common_b(ineigh)
              r1 = s%atom(iatom)%ratom + s%xl(ibeta)%a
              in1 = s%atom(iatom)%imass
              norb_mu = species(in1)%norb_max

              jatom = s%neighbors_PP(ialpha)%jatom_common_j(ineigh)
              jbeta = s%neighbors_PP(ialpha)%jatom_common_b(ineigh)
              r2 = s%atom(jatom)%ratom + s%xl(jbeta)%a
              in2 = s%atom(jatom)%imass
              norb_nu = species(in2)%norb_max

! PERFORM ACTUAL CALCULATIONS
! ASSEMBLE NON-LOCAL PSEUDOPOTENTIAL PIECE.
! Here 1=iatom, 2=jatom, 3=ialp:  <1 | V(3) | 2>.
              r31 = rna - r1
              r32 = rna - r2

! Find m value for iatom, ialp pair, and for jatom, ialp pair.
              call mpairnay (s, iatom, ialpha, r31, m31)
              call mpairnay (s, jatom, ialpha, r32, m32)

              ! cut some lenghty notation
              pvnl => s%vnl(iatom)
              pvnl_neighbors => pvnl%neighbors(mneigh)

              psvnl1 => s%svnl(iatom)
              psvnl1_neighbors => psvnl1%neighbors(m31)

              psvnl2 => s%svnl(jatom)
              psvnl2_neighbors => psvnl2%neighbors(m32)

! A "general" matrix element is:
! <iatom|VNL(3)|jatom> = <iatom|V(3)><V(3)|jatom> =
! SUM_ncc cl(ncc)*sVNL(mu,ncc,iatom,m31)*sVNL(nu,ncc,jatom,m32)
              allocate (bcnlx (norb_mu,norb_nu)); bcnlx = 0.0d0
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  do ncc = 1, species(indna)%norb_PP_max
                    bcnlx(imu,inu) = bcnlx(imu,inu)                          &
     &               + cl_value(ncc)*psvnl1_neighbors%block(imu,ncc)         &
     &                              *psvnl2_neighbors%block(inu,ncc)
                  end do ! do ncc
                end do ! do imu
              end do ! do inu

! Add this piece for iatom, jatom, and katom into the total (bcnlx ===> vnl)
!             pvnl_neighbors%block = pvnl_neighbors%block + bcnlx

! Deallocate Arrays
              deallocate (bcnlx)
            end if
          end do ! end loop over neighbors
          deallocate (cl_value)
        end do ! end loop over ialpha

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine assemble_vnl_3c


! ===========================================================================
! mpairnay
! ===========================================================================
! Program Description
! ===========================================================================
!       This is a neighbor routine, that we need for the nonlocal
! pseudopotential. We need this for the M_assemble_vnl_3c routine.
!
! < iatom | jatom >, given iatom, jatom and r12 = r(iatom) - r(jatom),
! then what neighbor to iatom is jatom. In other words, what is mneigh
! for a given jatom.
!
! We need this because the interaction are stored as sVNL(mu,nu,iatom,mneigh)
! and not sVNL(mu,nu,iatom,jatom). So what the heck is mneigh for jatom?
!
! This routine is a function call.
!
! Variables are iatom, jatom, rdiff(3). Other input are "constants".
! The output is mpairnay, which is the mneigh value for iatom, jatom.
!
! ===========================================================================
! Code written by:
!> @author Khorgolkhuu Odbadrakh
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================

! Program Declaration
! ===========================================================================
        subroutine mpairnay (t, iatom, jatom, r12, mpair)
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: t            !< the structure to be used

        integer, intent(in) :: iatom, jatom

        real, intent(in), dimension (3) :: r12

! Output
        integer, intent (out) :: mpair

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer imatch                   !< if we have a match then set to unity
        integer ineigh                   !< counter over neighbors
        integer jatom_match              !< jatom wanting to match
        integer mbeta                    !< cell for test atom

        real difference

        real, dimension (3) :: r1, r2, r21

        interface
          function distance (a, b)
            real distance
            real, intent(in), dimension (3) :: a, b
          end function distance
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize
        mpair = 0
        r1 = t%atom(iatom)%ratom

! Loop over the neighbors of iatom.
        imatch = 0
        do ineigh = 1, t%neighbors_PP(iatom)%neighn
          jatom_match = t%neighbors_PP(iatom)%neigh_j(ineigh)

! Eliminate the obvious.
          if (jatom_match .eq. jatom) then

! OK we have a candidate. If the vector connceting them is
! rdiff(3). Then we have a match. if there is no match, then
! we should not be here. Then there is a problem.
            mbeta = t%neighbors_PP(iatom)%neigh_b(ineigh)
            r2 = t%atom(jatom)%ratom + t%xl(mbeta)%a

! Vector from 1 to 2 is r21
            r21 = r2 - r1

! Now compare rdiff to r21.
            difference = distance (r12, r21)
            ! distance between r12 (input) and r21 (calculated)
            ! If there is a match, this should be zero

! We should find only one match.
            if (difference .lt. 0.0001d0) then
              imatch = imatch + 1
              mpair = ineigh
            end if
          end if
        end do

! Sanity checks
        if (imatch .ne. 1) then
          write (*,*) ' imatch = ', imatch
          write (*,*) ' The variable imatch MUST be ONE! NO EXCEPTIONS '
          write (*,*) ' Bad imatch value in mpairnay.f90; must abort! '
          write (*,*) ' iatom, position = ', iatom, r1(:)
          write (*,*) ' jatom, position = ', jatom, r2(:)
          stop
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
        end subroutine mpairnay


! End Module
! ===========================================================================
        end module M_assemble_PP_3c
