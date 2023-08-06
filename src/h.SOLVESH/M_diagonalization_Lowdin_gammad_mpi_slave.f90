! copyright info:
!
!                             @Copyright 2023
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
!
! M_diagonalization
! ===========================================================================
! Program Description
! ===========================================================================
!>       This is a version of matrix diagonalization for the Gamma kpoint.
!! The set of routines here use the blas library.
!!
!!      It contains the following subroutines within the module:
!!
!!      diagonalization_initialize - initialize the n x n matrices
!!      diagonalize_S - diagonalizes the overlap matrix
!!      diagonalize_H_Lowdin - perform Lowdin transformation of Hamiltonian
!!                             and then diagonalize the transformed Hamiltonian
!!
!!      This version uses mpi for parallel diagonalization.
!
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
! Module Declaration
! ============================================================================
        module M_diagonalization_slave

! /SYSTEM
        use M_configuraciones

! /MPI
        use M_mpi

! Type declarations for Hamiltonian matrix in k-space
! =========================================================================
! Define eigenvalues as 1 dimensional array
!       double precision, allocatable :: eigen (:)

! define matrices
!       double precision, allocatable :: Smatrix (:, :)
!       double precision, allocatable :: Hmatrix (:, :)

! define parameter for linear dependence criteria
!       double precision, parameter :: overtol = 1.0d-4

! MPI Stuff:
        integer, dimension (3) :: mybuffer

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
! Program Declaration
! ===========================================================================
        subroutine diagonalization_initialize_slave (s, iscf_iteration, ikpoint)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: iscf_iteration   !< which scf iteration?
        integer, intent (in) :: ikpoint          !< which kpoint

        type(T_structure), target :: s           !< the structure to be used

        type(T_kpoint), pointer :: pkpoint       !< point to specific kpoint

! Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: desc_length = 10

! Variable Declaration and Description
! ===========================================================================
! MPI Broadcasting information
        integer ierror
        integer npcol, nprow, nb

        integer ic, ir, mq0, np0
        integer mycol, myrow
        integer icontext, info

        integer, external :: numroc

        integer, dimension (desc_length) :: desc_x
        integer, dimension (desc_length) :: desc_y
        integer, dimension (desc_length) :: desc_z

        real, dimension (:, :), allocatable :: xxxx
        real, dimension (:, :), allocatable :: yyyy
        real, dimension (:, :), allocatable :: zzzz

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize BLACS
        call MPI_BCAST (mybuffer, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        nprow = mybuffer(1)                ! number of processors per row
        npcol = mybuffer(2)                ! number of processors per column
        nb = mybuffer(3)

        call blacs_pinfo (my_proc, nprocessors)
        call blacs_get (0, 0, icontext)
        call blacs_gridinit (icontext, 'R', nprow, npcol)
        call blacs_gridinfo (Icontext, nprow, npcol, myrow, mycol)

        if (myrow .eq. -1) then
          call MPI_FINALIZE (ierror)
          stop
        end if

! Reduce memory requirements for parallel diagonalization
        ir = max(1, numroc(s%norbitals, nb, myrow, 0, nprow))
        ic = max(1, numroc(s%norbitals, nb, mycol, 0, npcol))
        np0 = numroc(s%norbitals, nb, 0, 0, nprow)
        mq0 = numroc(s%norbitals, nb, 0, 0, npcol)

        if (iscf_iteration .eq. 1) then
          ! cut some lengthy notation
          pkpoint=>s%kpoints(ikpoint)
          allocate (pkpoint%S12matrix (1:ir, 1:ic))
          pkpoint%S12matrix = 0.0d0

          ! allocate working arrays
          allocate (xxxx (1:ir, 1:ic))
          allocate (yyyy (1:ir, 1:ic))
          allocate (zzzz (1:ir, 1:ic))
        end if

        call descinit (desc_x, s%norbitals, s%norbitals, nb, nb, 0, 0,        &
     &                 icontext, s%norbitals, info)
        call descinit (desc_y, s%norbitals, s%norbitals, nb, nb, 0, 0,        &
     &                 icontext, s%norbitals, info)
        call descinit (desc_z, s%norbitals, s%norbitals, nb, nb, 0, 0,        &
     &                 icontext, s%norbitals, info)

! End subroutine
! ===========================================================================
        return
        end subroutine diagonalization_initialize_slave


! ===========================================================================
! diagonalize_S_slave
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is diagonalization subroutine for real matrix using blas
! libraries....
!
! Program Declaration
! ===========================================================================
        subroutine diagonalize_S_slave (s)
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
        integer info                    !< error information
        integer liwork                  !< size of the integer working array
        integer lrwork                  !< size of the working array
        integer logfile                 !< writing to which unit
        integer mineig                  !< minimum non-zero eigenvalue
        integer imu, jmu                !< counters over eigenstates

        integer, allocatable :: iwork (:)           ! integer working vector
        double precision, allocatable :: rwork (:)  ! real working vector

! Allocate Arrays
! ===========================================================================
        lrwork = 1
        liwork = 1
        allocate (rwork(lrwork))
        allocate (iwork(liwork))

! Procedure
! ===========================================================================
! DIAGONALIZE THE OVERLAP MATRIX
! ***************************************************************************
        ! first find optimal length of rwork
!       call dsyevd ('V', 'U', s%norbitals, Smatrix, s%norbitals, eigen,      &
!    &               rwork, -1, iwork, -1, info)
!       lrwork = rwork(1)
!       liwork = iwork(1)
!       deallocate (rwork, iwork)
!       allocate (rwork(lrwork)); allocate (iwork(liwork))
!       call dsyevd ('V', 'U', s%norbitals, Smatrix, s%norbitals, eigen,      &
!    &               rwork, lrwork, iwork, liwork, info)
! NOTE: After calling dsyev, Smatrix now becomes the eigenvectors of the
! diagonalized Smatrix!

! Error check in diagonalization process
!       if (info .ne. 0) then
!         write (*,*) '  '
!         write (*,*) ' Diagonalization not successful, info = ', info
!         if (info .lt. 0) then
!           write (*,*) ' The ', info, '-th argument had an illegal value'
!          else
! LAPACK style errors
!           write (*,*) ' It failed to converge.'
!           write (*,*) info, ' off-diagonal elements of an intermediate'
!           write (*,*) ' tridiagonal form did not converge to zero. '
!         end if
!         stop
!       end if

! ***************************************************************************
! CHECK THE LINEAR DEPENDENCE
! ***************************************************************************
! Fix the linear dependence problem. References: Szabo and Ostlund, Modern
! Quantum Chem. McGraw Hill 1989 p. 142; Szabo and Ostlund, Modern Quantum
! Chem. Dover 1996 p. 145. A tolerance for a small overlap eigenvalue is
! set by overtol.

! Determine the smallest active eigenvector
!       mineig = 0
!       do imu = 1, s%norbitals
!         if (eigen(imu) .lt. overtol) mineig = imu
!       end do

!       mineig = mineig + 1
!       s%norbitals_new = s%norbitals + 1 - mineig
!       if (s%norbitals_new .ne. s%norbitals) then
!         write (logfile,*) '  '
!         write (logfile,*) ' WARNING. ### ### ### '
!         write (logfile,*) ' Linear dependence encountered in eigenvectors. '
!         write (logfile,*) ' Eigenvalue is very small. '
!         write (logfile,*) s%norbitals - s%norbitals_new, ' vectors removed.'
!         do imu = mineig, s%norbitals
!           jmu = imu - mineig + 1
!           Smatrix(:,jmu) = Smatrix(:,imu)
!           eigen(jmu) = eigen(imu)
!         end do
!       end if

! Deallocate Arrays
! ===========================================================================
        deallocate (iwork, rwork)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine diagonalize_S_slave


! End the module
! ===========================================================================
        end module M_diagonalization_slave

