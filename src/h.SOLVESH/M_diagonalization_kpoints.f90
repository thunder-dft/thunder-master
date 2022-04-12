! copyright info:
!                             @Copyright 2013
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel JelinekmNZxbnmb

!
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

! M_diagonalization
! ===========================================================================
! Program Description
! ===========================================================================
!       This is a version of matrix diagonalization for the Gamma kpoint.
! The set of routines here use the blas library.
!      It contains the following subroutines within the module:
!
!      diagonalization_initialize - initialize the n x n matrices
!      diagonalize_H - diagonalize the generalized eigenvalue equation of the
!                      Hamiltonian in non-orthogonal basis
!
! ============================================================================
! Module declaration
! ============================================================================
        module M_diagonalization
        use M_configuraciones

! Type declarations for Hamiltonian matrix in k-space
! =========================================================================

! Define eigenvalues (lam) as 1 dimensional array
        double precision, pointer :: eigen (:)

! define matrices
        double complex, pointer :: Smatrix (:,:)
        double complex, pointer :: Hmatrix (:,:)

! define parameter for linear dependence criteria
        double precision, parameter :: overtol = 1.0d-4

! this value detemines if the eigenvalue should be accounted for:
! lam<overtol then it is considered zero and otherwise if lam>overtol
! ===========================================================================
! module procedures
        contains


! ===========================================================================
! diagonalization_initialize
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine initializes real dummy matrices used in kspace.
! The dimensions of these dummy matrices must be initialized to
! norbitals by norbitals.
!
! ===========================================================================
! Code written by:
! Prokop Hapala
! Department of Thin Films
! Institute of Physics
! Czech Academy of Sciences
! Prague, Czech Republic

! with modifications by:
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
        subroutine diagonalization_initialize (s, iscf_iteration)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           ! the structure to be used

        integer, intent (in) :: iscf_iteration   ! which scf iteration?

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
         if (iscf_iteration .eq. 1) then
           allocate (eigen(s%norbitals))

           allocate (Smatrix (s%norbitals, s%norbitals))
           allocate (Hmatrix (s%norbitals, s%norbitals))
         end if
         Smatrix = 0.0d0
         Hmatrix = 0.0d0

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine diagonalization_initialize


! ===========================================================================
! diagonalize_H
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine performs diagonalizaton on the generalized eigenvalue
! non-orthogonal Hamiltonian equation in k-space: H*psi = E*S*psi. After
! diagonalization, the eigenvalues are stored in the structure
! kpoints(ikpoint)%eigen(:)

! ===========================================================================
! Code written by:
! Prokop Hapala
! Department of Thin Films
! Institute of Physics
! Czech Academy of Sciences
! Prague, Czech Republic

! with modifications by:
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
        subroutine diagonalize_H (s, ikpoint)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: ikpoint           !< which kpoint

        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! ===========================================================================
! Local variables declarations
        integer info                   ! error information
        integer lrwork, lwork          ! size of the working arrays

        double precision, pointer :: rwork (:)        ! working vector

        double complex, allocatable :: work(:)

! Allocate Arrays
! ===========================================================================
        lwork = 1
        allocate (work(lwork))
        lrwork = 3*s%norbitals - 2
        allocate (rwork(lrwork))

! Procedure
! ===========================================================================
! Initialize some constants
        s%norbitals_new = s%norbitals  ! we need to set this for forces later

! Now, DIAGONALIZE THE HAMILTONIAN in the non-orthogonal basis set
! ****************************************************************************
! Eigenvectors are needed to calculate the charges and for forces!
        call zhegv (1, 'V', 'U', s%norbitals, Hmatrix, s%norbitals, Smatrix,         &
     &               s%norbitals, eigen, work, -1, rwork , info)
        ! first find optimal length of work
        lwork = work(1)
        deallocate (work)
        allocate (work(lwork))

        call zhegv (1, 'V', 'U', s%norbitals, Hmatrix, s%norbitals, Smatrix,         &
     &               s%norbitals, eigen, work, lwork, rwork, info)

! INFORMATION FOR THE LOWDIN CHARGES
! After calling diagonalization_H_Lowdin what we have remaining is Hmatrix
! which are the eigenvectors of the transformed H matrix and Smatrix which
! is at this point the eigenstates of the untransformed H matrix (or the
! eigenstates in the atomic orbital basis).
        s%kpoints(ikpoint)%c = Hmatrix

! Deallocate Arrays
! ===========================================================================
        deallocate (rwork)
        deallocate (work)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine diagonalize_H


! ===========================================================================
! phase
! ===========================================================================
! Function Description
! ===========================================================================
!       Returns the phase of k*r.
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
        function phase (dot)
        implicit none

        complex phase

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent(in) :: dot    !> the dot product of k*r

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        phase = cmplx(cos(dot),sin(dot))

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end function phase

! End the module
! ===========================================================================
        end module M_diagonalization

