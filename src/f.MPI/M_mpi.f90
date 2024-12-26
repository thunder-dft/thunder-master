! copyright info:
!
!                             @Copyright 2023
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

! Module Description
! ============================================================================
!      This module initializes and finalizes the MPI space.
! However, it is a dummy routine.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
! Module declaration
! ============================================================================
        module M_mpi
        use mpi

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer my_proc
        integer nprocs

! module procedures
        contains

! ===========================================================================
! initialize_mpi
! ===========================================================================
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
        subroutine initialize_mpi
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi

! Procedure
! ===========================================================================
        call MPI_INIT (ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write(*,*) 'mpi error MPI_Init'
          stop
        end if

        call MPI_COMM_RANK (MPI_COMM_WORLD, my_proc, ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write (*,*) ' mpi error MPI_COMM_RANK'
          stop
        end if

        call MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ierr_mpi)
        if (ierr_mpi .ne. 0) then
          write (*,*) 'mpi error MPI_COMM_SIZE'
          stop
        end if

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end subroutine initialize_mpi


! ===========================================================================
! finalize_mpi
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine finalizes the MPI space.
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine finalize_mpi
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi

! Procedure
! ===========================================================================
        call MPI_barrier (MPI_COMM_WORLD, ierr_mpi)
        call MPI_finalize (ierr_mpi)

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end subroutine finalize_mpi


! ===========================================================================
! awake_slave
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine is used to wake up the slave processes and stop them.
! ===========================================================================
! Code written by:
! Runfeng Jin
!
! jsfaraway@gmail.com 
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine awake_slave
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ierr_mpi

! Procedure
! ===========================================================================
        if (nprocs .gt. 0) then
            if(my_proc .gt. 0) then
                stop "Slave process is not allowed to call this subroutine"
            endif
            call MPI_Bcast (-1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
        endif

! Format Statements
! ===========================================================================
! None

! End Function
! ===========================================================================
        return
        end subroutine awake_slave
        

! End the module
! ===========================================================================
        end module M_mpi
