#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

! copyright info:
!                             @Copyright 2016
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
!
! M_diagonalization
! ===========================================================================
! Program Description
! ===========================================================================
!>       This is a version of matrix diagonalization for the Gamma kpoint.
!  The set of routines here use the blas library.
!       It contains the following subroutines within the module:
! 
!       diagonalization_initialize - initialize the n x n matrices
!       diagonalize_S - diagonalizes the overlap matrix
!       diagonalize_H_Lowdin - perform Lowdin transformation of Hamiltonian
!                              and then diagonalize the transformed Hamiltonian
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
! ============================================================================
! Module declaration
! ============================================================================
        module M_diagonalization
        use M_configuraciones

! Type declarations for Hamiltonian matrix in k-space
! =========================================================================

! Define eigenvalues as 1 dimensional array
        double precision, allocatable :: eigen (:)

! define matrices
#ifdef GAMMA
        double precision, allocatable :: Smatrix (:, :)
        double precision, allocatable, save :: S12matrix (:, :)
        double precision, allocatable :: Hmatrix (:, :)
#else
        double complex, allocatable :: Smatrix (:, :)
        double complex, allocatable, save :: S12matrix (:, :)
        double complex, allocatable :: Hmatrix (:, :)
#endif

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
!>       This subroutine initializes real dummy matrices used in kspace.
! The dimensions of these dummy matrices must be initialized to
! norbitals by norbitals.
!
! ===========================================================================
! Code written by:
!> @author Kh. Odbadrakh
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
        integer, intent (in) :: iscf_iteration   !< which scf iteration?

        type(T_structure), target :: s            !< the structure to be used

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
! diagonalize_S
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is diagonalization subroutine for real matrix using blas
! libraries....
!
! ===========================================================================
! Code written by:
!> @author Kh. Odbadrakh
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
        subroutine diagonalize_S (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! The lam value detemines if the eigenvalue should be accounted for:
! If lam < overtol then it is considered zero and otherwise if lam > overtol
! For zero lam, then there is a linear dependency problem and we need to
! only consider a submatrix.
        integer info                        !< error information
        integer lrwork                      !< size of the working array
        integer logfile                     !< writing to which unit
        integer mineig                      !< minimum non-zero eigenvalue
        integer imu, jmu                    !< counters over eigenstates

#ifdef GAMMA
        double precision, allocatable :: rwork (:)   ! working vector
#else
        integer lwork
        double complex, allocatable :: rwork (:)   ! working vector
        double complex, allocatable :: work (:)
#endif

! Allocate Arrays
! ===========================================================================
#ifdef GAMMA
        lrwork = 1
#else
        lwork = 1
        allocate (work(lwork))
        lrwork= 3*s%norbitals - 2
#endif
        allocate (rwork(lrwork))

! Procedure
! ===========================================================================
! DIAGONALIZE THE OVERLAP MATRIX
! ***************************************************************************
! Initialize logfile
        logfile = s%logfile
        write (logfile,*)
        write (logfile,*) ' Call diagonalizer '

#ifdef GAMMA
        call dsyev ('V', 'U', s%norbitals, Smatrix, s%norbitals, eigen,      &
     &               rwork, -1, info)
        ! first find optimal length of rwork
        lrwork = rwork(1)
        deallocate (rwork)
        allocate(rwork(lrwork))
        call dsyev ('V', 'U', s%norbitals, Smatrix, s%norbitals, eigen,      &
     &               rwork, lrwork, info)

#else
        call zheev ('V', 'U', s%norbitals, Smatrix, s%norbitals, eigen, work,  &
     &              -1, rwork , info)
        !first find optimal length of work
        lwork = work(1)
        deallocate (work)
        allocate (work (lwork))
        call zheev ('V', 'U', s%norbitals, Smatrix, s%norbitals, eigen, work,  &
     &              lwork, rwork , info)
#endif

! NOTE: After calling dsyev, Smatrix now becomes the eigenvectors of the
! diagonalized Smatrix!

! Error check in diagonalization process
        if (info .ne. 0) then
          write (*,*) '  '
          write (*,*) ' Diagonalization not successful, info = ', info
          if (info .lt. 0) then
            write (*,*) ' The ', info, '-th argument had an illegal value'
          else
! LAPACK style errors
            write (*,*) ' It failed to converge.'
            write (*,*) info, ' off-diagonal elements of an intermediate'
            write (*,*) ' tridiagonal form did not converge to zero. '
          end if
          stop
        end if

! ***************************************************************************
! CHECK THE LINEAR DEPENDENCE
! ***************************************************************************
! Fix the linear dependence problem. References: Szabo and Ostlund, Modern
! Quantum Chem. McGraw Hill 1989 p. 142; Szabo and Ostlund, Modern Quantum
! Chem. Dover 1996 p. 145. A tolerance for a small overlap eigenvalue is
! set by overtol.

! Determine the smallest active eigenvector
          mineig = 0
          do imu = 1, s%norbitals
            if (eigen(imu) .lt. overtol) mineig = imu
          end do

          mineig = mineig + 1
          s%norbitals_new = s%norbitals + 1 - mineig
          if (s%norbitals_new .ne. s%norbitals) then
            write (logfile,*) '  '
            write (logfile,*) ' WARNING. ### ### ### '
            write (logfile,*) ' Linear dependence encountered in eigenvectors. '
            write (logfile,*) ' Eigenvalue is very small. '
            write (logfile,*) s%norbitals - s%norbitals_new, ' vectors removed. '
            do imu = mineig, s%norbitals
              jmu = imu - mineig + 1
              Smatrix(:,jmu) = Smatrix(:,imu)
              eigen(jmu) = eigen(imu)
            end do
          end if

! Deallocate Arrays
! ===========================================================================
        deallocate (rwork)
#ifdef GAMMA
        !No more deallocations
#else
        deallocate (work)
#endif
! None
! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine diagonalize_S


! ===========================================================================
! diagonalize_H_Lowdin
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This subroutine performs Lowdin transformation on the Hamiltonian
!! matix H in k space: (S^-1/2)*H*(S^-1/2) and then diagonalizes it, and
!! stores the eigenvalues in the structure kpoints(ikpoint)%eigen(:)


! ===========================================================================
! Code written by:
!> @author Kh. Odbadrakh
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
        subroutine diagonalize_H_Lowdin (s, iscf_iteration, ikpoint)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

        integer, intent (in) :: iscf_iteration   !< which scf iteration?
        integer, intent (in) :: ikpoint

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer info                   ! error information
        integer lrwork                 ! size of the working array
        integer imu                    ! counters over eigenstates

        real sqlami        ! square root of overlap eigenvalues

#ifdef GAMMA
        double precision, allocatable :: rwork (:)        ! working vector
        double precision a0r, a1r
#else
        integer lwork
        double complex, allocatable :: work(:)
#endif

        ! checks to see if structure has changed
        type(T_structure), pointer, save :: current

        character (len = 25) :: slogfile

! Allocate Arrays
! ===========================================================================
#ifdef GAMMA
        lrwork = 1
#else
        lwork = 1
        allocate (work(lwork))
        lrwork = 3*s%norbitals - 2
#endif
        allocate (rwork(lrwork))

! Procedure
! ===========================================================================
! Check to see if the structure has changed
        if (.not. associated (current, target=s)) then
          current => s
        else if (.not. associated(current, target=s) .and. allocated (S12matrix)) then
          deallocate (S12matrix)
          current => s
        end if

        if (iscf_iteration .eq. 1) then
          allocate (S12matrix (s%norbitals, s%norbitals)); S12matrix = 0.0d0
        end if

! Initialize some constants
        a0r = 0.0d0
        a1r = 1.0d0
        s%norbitals_new = size(eigen,1)

! ****************************************************************************
! CALCULATE (S^-1/2) --> lam12
! ****************************************************************************
! In a diagonal representation (Udagger*S*U = s, s is a diagonal matrix)
! We just take the inverse of the square roots of the eigenvalues to get
! s^-1/2. Then we 'undiagonalize' the s^-1/2 matrix back to get
! S^-1/2 = U*s^-1/2*Udagger.
! Note: We do S^-1/4 here, because the sqlami contribution get squared
! after it is combined with overlap.
! NOTE: Smatrix here is really the eigenstates now of S^-1/4

! Only do this part if this is the first scf iteration
        if (iscf_iteration .eq. 1) then
          do imu = 1, s%norbitals_new
            sqlami = eigen(imu)**(-0.25d0)
            Smatrix(:,imu) = Smatrix(:,imu)*sqlami
          end do
#ifdef GAMMA 
         call dgemm ('N', 'C', s%norbitals, s%norbitals, s%norbitals_new,  &
     &                a1r, Smatrix, s%norbitals, Smatrix, s%norbitals, a0r,   &
     &                S12matrix, s%norbitals)
#else
          call zgemm ('N', 'C', s%norbitals, s%norbitals, s%norbitals_new,   &
     &                 a1, Smatrix, s%norbitals, Smatrix, s%norbitals, a0,   &
     &                 S12matrix, s%norbitals)
#endif
        end if

! NOTE: S12matrix here NOW contains the matrix S^-1/2

! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
! lam12 = S^-1/2 in AO basis
! zz = Unused dummy matrix
! yy = Hamiltonian eigenvectors in the MO basis
! Set M=H*(S^-.5)
! We do not need Smatrix any more at this point, so use it as a dummy
        Smatrix = 0.0d0

#ifdef GAMMA
        call dsymm ('R', 'U', s%norbitals, s%norbitals, a1r, S12matrix,       &
     &               s%norbitals, Hmatrix, s%norbitals, a0r, Smatrix,         &
     &               s%norbitals)

! Set Z=(S^-.5)*M
        call dsymm ('L', 'U', s%norbitals, s%norbitals, a1r, S12matrix,       &
     &               s%norbitals, Smatrix, s%norbitals, a0r, Hmatrix,         &
     &               s%norbitals)
#else
        call zhemm ( 'R', 'U', s%norbitals, s%norbitals, a1, S12matrix,      &
     &               s%norbitals, Hmatrix, s%norbitals, a0, Smatrix,         &
     &               s%norbitals )

! Set Z=(S^-.5)*M
        call zhemm ( 'L', 'U', s%norbitals, s%norbitals, a1, S12matrix,      &
     &               s%norbitals, Smatrix, s%norbitals, a0, Hmatrix,         &
     &               s%norbitals)
#endif

        if (iwriteout_dos .eq. 1) then
          slogfile = s%basisfile(:len(trim(s%basisfile))-4)
          slogfile = trim(slogfile)//'.Hk'
          open (unit = 22, file = slogfile, status = 'replace', form = 'unformatted')
          write (22) Hmatrix
          close (unit = 22)
        end if

! Now, DIAGONALIZE THE HAMILTONIAN in the orthogonal basis set
! ****************************************************************************
! Eigenvectors are needed to calculate the charges and for forces!

#ifdef GAMMA
        call dsyev ('V', 'U', s%norbitals, Hmatrix, s%norbitals, eigen,      &
     &               rwork, -1, info)
        ! first find optimal length of rwork
        lrwork = rwork(1)
        deallocate (rwork)
        allocate (rwork(lrwork))
        call dsyev ('V', 'U', s%norbitals, Hmatrix, s%norbitals, eigen,      &
     &               rwork, lrwork, info)
#else
        call zheev ('V', 'U', s%norbitals, Hmatrix, s%norbitals, eigen,      &
     &              work, -1, rwork , info)
        ! first find optimal length of work
        lwork = work(1)
        deallocate (work)
        allocate (work (lwork))
        call zheev ('V', 'U', s%norbitals, Hmatrix, s%norbitals, eigen,      &
     &              work, lwork, rwork , info)
#endif

! INFORMATION FOR THE LOWDIN CHARGES
! ****************************************************************************
! S12matrix = S^-1/2 in AO basis
! We do not need Smatrix any more at this point, so use it as a dummy
        Smatrix = 0.0d0

#ifdef GAMMA
        call dsymm ('L', 'U', s%norbitals, s%norbitals, a1r, S12matrix,       &
     &               s%norbitals, Hmatrix, s%norbitals, a0r, Smatrix,         &
     &               s%norbitals)
#else
        call zhemm ('L', 'U', s%norbitals, s%norbitals, a1, S12matrix,       &
     &               s%norbitals, Hmatrix, s%norbitals, a0, Smatrix,         &
     &               s%norbitals)
#endif

! NOTE: After multiplication Smatrix = S^-1/2 * yy
! We did a symmetric orthogonalization followed by a diagonalization
! of the Hamiltonian in this "MO" basis set. This yields a net
! canonical diagonalization with matrix bbnk.

! INFORMATION FOR THE LOWDIN CHARGES
! After calling diagonalization_H_Lowdin what we have remaining is Hmatrix
! which are the eigenvectors of the transformed H matrix and Smatrix which
! is at this point the eigenstates of the untransformed H matrix (or the
! eigenstates in the atomic orbital basis).
#ifdef GAMMA
        do imu = 1, s%norbitals_new
          s%kpoints(ikpoint)%c_Lowdin(:,imu) = Hmatrix(:,imu)
          s%kpoints(ikpoint)%c(:,imu) = Smatrix(:,imu)
        end do
#else
        s%kpoints(ikpoint)%c_Lowdin = Hmatrix
        s%kpoints(ikpoint)%c = Smatrix
#endif

! Deallocate Arrays
! ===========================================================================
        deallocate (rwork)
#ifdef GAMMA
        ! No more deallocations
#else
        deallocate (work)
#endif

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine diagonalize_H_Lowdin


! phase.f90
! Function Description
! ===========================================================================
!>       Returns the phase of k*r.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        function phase (dot)
        implicit none

#ifdef GAMMA
        real phase
#else
        complex phase
#endif

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
#ifdef GAMMA
        phase = 1.0d0
        if(.FALSE.) phase=dot
#else
        phase = cmplx(cos(dot),sin(dot))
#endif

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
