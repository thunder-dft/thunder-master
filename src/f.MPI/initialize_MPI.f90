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

! initialize_MPI.f90
! Function Description
! ============================================================================
!      This function initializes the MPI space.
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
        subroutine initialize_MPI (iammaster, iammpi, my_proc, nproc)
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Output
        integer, intent (out) :: my_proc
        integer, intent (out) :: nproc

        logical, intent (out) :: iammaster
        logical, intent (out) :: iammpi

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi

! Procedure
! ===========================================================================
        call MPI_init (ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write(*,*) 'mpi error MPI_Init'
          stop
        end if
        call MPI_comm_rank (MPI_COMM_WORLD, my_proc, ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write (*,*) ' mpi error MPI_Comm_rank'
          stop
        end if
        call MPI_comm_size (MPI_COMM_WORLD, nproc, ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write (*,*) 'mpi error MPI_Comm_size'
          stop
        end if
        iammaster = .false.
        if (my_proc .eq. 0) iammaster = .true.
        iammpi = .true.

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end subroutine initialize_MPI
