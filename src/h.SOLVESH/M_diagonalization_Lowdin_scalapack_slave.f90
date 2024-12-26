! copyright info:
!
!                             @Copyright 2024
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
! Computer Network Information Center, Chinese Academy of Sciences    &
!     & University of Chinese Academy of Sciences - Runfeng Jin
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
!>       This is a version of matrix diagonalization for the Gamma kpoint using   ScaLAPACK library.
!! The set of routines here use the blas library.
!!
!!      It contains the following subroutines within the module:
!!
!!      diagonalization_initialize - initialize the n x n matrices
!!      diagonalize_S - diagonalizes the overlap matrix using ScaLAPACK
!!      diagonalize_H_Lowdin - perform Lowdin transformation of Hamiltonian
!!                             and then diagonalize the transformed Hamiltonian using ScaLAPACK
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
!
! Runfeng Jin
! Computer Network Information Center, Chinese Academy of Sciences
! Beijing, China
!
! E-mail: jsfaraway@gmail.com
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
        integer, parameter, private :: desc_length = 10
        ! Define eigenvalues as 1 dimensional array
        double precision, allocatable, private :: eigen(:)
        ! MPI information
        integer, private:: npcol, nprow, nb

        integer, private :: ic, ir, mq0, np0
        integer, private:: mycol, myrow, icontext

! working matrices - for parallel divide and conquer

        integer, dimension(desc_length) :: desc_x
        integer, dimension(desc_length) :: desc_y
        integer, dimension(desc_length) :: desc_z

        double precision, allocatable, private :: xxxx(:, :)
        double precision, allocatable, private :: yyyy(:, :)
        double precision, allocatable, private :: zzzz(:, :)
        double precision, allocatable, private :: trival_full(:, :)
        double precision, private :: time_begin, time_end

! define parameter for linear dependence criteria
        double precision, parameter :: overtol = 1.0d-4

! MPI Stuff:
        integer, dimension(3), private :: mybuffer

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
        subroutine diagonalization_initialize_slave(s, iscf_iteration, ikpoint)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: iscf_iteration   !< which scf iteration?
        integer, intent(in) :: ikpoint          !< which kpoint

        type(T_structure), target :: s           !< the structure to be used

        type(T_kpoint), pointer :: pkpoint       !< point to specific kpoint

! Parameters and Data Declaration
! ===========================================================================
!       integer, parameter :: desc_length = 10

! Variable Declaration and Description
! ===========================================================================
! MPI Broadcasting information
        integer ierror
        integer info

        integer, external :: numroc

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize BLACS
        ! find the optimal process distribute
        do npcol = NINT(SQRT(REAL(nprocs, kind=8))), 2, -1
          if (mod(nprocs, npcol) == 0) exit
        end do
        nprow = nprocs/npcol;
        nb = 64

        call blacs_pinfo(my_proc, nprocs)
        call blacs_get(0, 0, icontext)
        call blacs_gridinit(icontext, 'R', nprow, npcol)
        call blacs_gridinfo(Icontext, nprow, npcol, myrow, mycol)

        if (myrow .eq. -1) then
          call MPI_FINALIZE(ierror)
          stop
        end if

! Reduce memory requirements for parallel diagonalization
        ! print rank id and ir ic np0 mq0
        ir = max(1, numroc(s%norbitals, nb, myrow, 0, nprow))
        ic = max(1, numroc(s%norbitals, nb, mycol, 0, npcol))
        np0 = numroc(s%norbitals, nb, 0, 0, nprow)
        mq0 = numroc(s%norbitals, nb, 0, 0, npcol)
        ! write(*, *) 'rank id, ir, ic, np0, mq0: ', my_proc, ir, ic, np0, mq0

        if (iscf_iteration .eq. 1) then
          ! if has been allocated, remove these
          if (allocated(xxxx)) deallocate(xxxx)
          if (allocated(yyyy)) deallocate(yyyy)
          if (allocated(zzzz)) deallocate(zzzz)
          if (allocated(trival_full)) deallocate(trival_full)
          ! cut some lengthy notation
          pkpoint => s%kpoints(ikpoint)
        !   if (associated(pkpoint%S12matrix)) deallocate(pkpoint%S12matrix)
          ! have allocate this in the main function, so here can directly deallocate it.
          ! This due to associated cannot always work well
          deallocate(pkpoint%S12matrix)
          allocate (pkpoint%S12matrix(1:ir, 1:ic))
          pkpoint%S12matrix = 0.0d0

          ! allocate working arrays
          allocate (xxxx(1:ir, 1:ic))
          allocate (yyyy(1:ir, 1:ic))
          allocate (zzzz(1:ir, 1:ic))
          allocate (trival_full(1:s%norbitals, 1:s%norbitals))
        end if

        call descinit(desc_x, s%norbitals, s%norbitals, nb, nb, 0, 0,        &
        &                 icontext, ir, info)
        call descinit(desc_y, s%norbitals, s%norbitals, nb, nb, 0, 0,        &
        &                 icontext,ir, info)
        call descinit(desc_z, s%norbitals, s%norbitals, nb, nb, 0, 0,        &
        &                 icontext, ir, info)
        if (iscf_iteration .eq. 1 .and. ikpoint .eq. 1) then
          if(allocated(eigen)) deallocate(eigen)
          allocate (eigen(s%norbitals))
        end if
! End subroutine
! ===========================================================================
        return
        end subroutine diagonalization_initialize_slave

! ===========================================================================
! diagonalize_S_slave
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This is diagonalization subroutine for real matrix using ScaLAPACK
! libraries....
!
! Program Declaration
! ===========================================================================
        subroutine diagonalize_S_slave(s)
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
        integer matrixfile

        integer, allocatable :: iwork(:)           ! integer working vector
        double precision, allocatable :: rwork(:)  ! real working vector

! write out the matrix
        INTEGER :: i, j, unit
        CHARACTER(LEN=20) :: format_string

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
        ! acquire the block-cylic Overlap matrix S from master in yyyy, xxxx is trival
        ! call pclaputter(yyyy, desc_y, xxxx, s%norbitals)
        call pclaputter(yyyy, desc_y, xxxx, s%norbitals)

        call pdsyevd('V', 'U', s%norbitals, yyyy, 1, 1, desc_y, eigen,       &
        &                xxxx, 1, 1, desc_x, rwork, -1, iwork, -1, info)
        lrwork = rwork(1)
        liwork = iwork(1)
        deallocate (rwork, iwork)
        allocate (rwork(lrwork)); allocate (iwork(liwork))
        ! compute the eigenvalue(in eigen) and eigenvector(in xxxx, stored block-cyclic
        ! as S in yyyy)
        call cpu_time(time_begin)
        call pdsyevd('V', 'U', s%norbitals, yyyy, 1, 1, desc_y, eigen,       &
        &                xxxx, 1, 1, desc_x, rwork, lrwork, iwork, liwork, info)
        call cpu_time(time_end)

! NOTE: After calling pdsyevd, xxxx now becomes the eigenvectors of the
! diagonalized Smatrix!

! Error check in diagonalization process
        if (info .ne. 0) then
          write (*, *) '  '
          write (*, *) ' Diagonalization not successful, info = ', info
          if (info .lt. 0) then
            write (*, *) ' The ', info, '-th argument had an illegal value'
          else
! LAPACK style errors
            write (*, *) ' It failed to converge.'
            write (*, *) info, ' off-diagonal elements of an intermediate'
            write (*, *) ' tridiagonal form did not converge to zero. '
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
          eigen(1:mineig - 1) = 0.0d0
        end if

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



! ===========================================================================
! diagonalize_H_Lowdin_slave
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This subroutine performs Lowdin transformation on the Hamiltonian
!! matix H in k space: (S^-1/2)*H*(S^-1/2) and then diagonalizes it, and
!! stores the eigenvalues in the structure kpoints(ikpoint)%eigen(:)
!
! Program Declaration
! ===========================================================================
        subroutine diagonalize_H_Lowdin_slave(s, iscf_iteration, ikpoint)
        implicit none
        ! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

        integer, intent(in) :: iscf_iteration   !< which scf iteration?
        integer, intent(in) :: ikpoint

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer info                   ! error information
        integer lrwork                 ! size of the working real array
        integer liwork                 ! size of the working integer array
        integer imu, inu                    ! counters over eigenstates

        double precision a0, a1
        double precision sqlami        ! square root of overlap eigenvalues

        integer, allocatable :: iwork(:)           ! integer working vector
        double precision, allocatable :: rwork(:)  ! real working vector

        character(len=25) :: slogfile

        type(T_kpoint), pointer :: pkpoint   !< point to current kpoint

! Allocate Arrays
! ===========================================================================
        lrwork = 1
        liwork = 1
        allocate (rwork(lrwork))
        allocate (iwork(liwork))
        ! Procedure
! ===========================================================================
! Initialize some constants
        a0 = 0.0d0
        a1 = 1.0d0
        s%norbitals_new = size(eigen, 1)

! Cut some lengthy notation
        pkpoint => s%kpoints(ikpoint)
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
          do imu = 1, s%norbitals
            if (eigen(imu) .lt. overtol) then
              sqlami = 0.0d0
            else
              sqlami = eigen(imu)**(-0.25d0)
            end if
            do inu = 1, s%norbitals
              call blacsaba(xxxx, inu, imu, desc_x, sqlami, mycol, myrow,      &
              &                    npcol, nprow)
            end do
          end do
          call pdgemm('N', 'C', s%norbitals, s%norbitals, s%norbitals, a1, xxxx,  &
          &             1, 1, desc_x, xxxx, 1, 1, desc_x, a0, yyyy, 1, 1, desc_y)
          ! Now put S^-1/2 into s(k)^-1/2, this will be remembered for the duration of
          ! the scf cycle.
          pkpoint%S12matrix(1:ir, 1:ic) = yyyy(1:ir, 1:ic)
        else ! not the first iteration
          ! restore S^-1/2 from the structure s
          yyyy(1:ir, 1:ic) = pkpoint%S12matrix(1:ir, 1:ic)
        end if ! end of if (iscf_iteration .eq. 1)
! NOTE: yyyy here NOW contains the matrix S^-1/2
! xxxx = nothing
! yyyy = S^-1/2 in AO basis
! zzzz = nothing

! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
! put Hmatrix into block-cyclic in xxxx
        call pclaputter(xxxx, desc_x, yyyy, s%norbitals)

! zzzz = Unused dummy matrix
! yyyy = S^-1/2 in AO basis
! xxxx = Hamiltonian eigenvectors in the AO basis, in block-cyclic
! Set M=H*(S^-.5)
! We do not need Smatrix any more at this point, so use it as a dummy

        call pdsymm('R', 'U', s%norbitals, s%norbitals, a1, yyyy, 1, 1,    &
        &                desc_y, xxxx, 1, 1, desc_x, a0, zzzz, 1, 1, desc_z)

! Set Z=(S^-.5)*M
        call pdsymm('L', 'U', s%norbitals, s%norbitals, 1.0d0, yyyy, 1, 1,    &
        &                desc_y, zzzz, 1, 1, desc_z, 0.0d0, xxxx, 1, 1, desc_x)

! xxxx = H in the MO basis，
! yyyy = S^-1/2 in AO basis
! zzzz = nothing

! Now, DIAGONALIZE THE HAMILTONIAN in the orthogonal basis set
! ****************************************************************************
! Eigenvectors are needed to calculate the charges and for forces!
! first find optimal length of rwork

        call pdsyev('V', 'U', s%norbitals, xxxx, 1, 1, desc_x, eigen, zzzz,   &
        &               1, 1, desc_z, rwork, -1, info)
        lrwork = rwork(1)
        deallocate (rwork, iwork)
        allocate (rwork(lrwork)); ! allocate (iwork(liwork))

        call pdsyev('V', 'U', s%norbitals, xxxx, 1, 1, desc_x, eigen, zzzz,   &
        &                1, 1, desc_z, rwork, lrwork, info)

        ! if (info .ne. 0) call diag_error(info, 1)
        ! xxxx = H in the MO basis，
        ! yyyy = S^-1/2 in AO basis
        ! zzzz = eigenvector of transformed H matrix

! INFORMATION FOR THE LOWDIN CHARGES
! ****************************************************************************
! Save the eigenvectors of the transformed H matrix
! put the eigenvector into Smatrix by MPI
        call pclagetter(zzzz, desc_z, trival_full, s%norbitals)

! Save the eigenvectors of untransformed H matrix

        call pdsymm('L', 'U', s%norbitals, s%norbitals, 1.0d0, yyyy, 1, 1,     &
        &               desc_y, zzzz, 1, 1, desc_z, 0.0d0, xxxx, 1, 1, desc_x)

        call pclagetter(xxxx, desc_x, trival_full, s%norbitals)

        deallocate (rwork)
        call blacs_gridexit (icontext)

        end subroutine diagonalize_H_Lowdin_slave

! End the module
! ===========================================================================
        end module M_diagonalization_slave
