! copyright info:
!
!                             @Copyright 2008
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

! cl.f90
! Program Description
! ===========================================================================
!>       This routine returns the Kleinman Bylander cl values for atom itype.
!! The "raw" date is read in vnl.z1.z2.dat by  program read2c. We include up to
!! 5 non-local (L values) of the pseudopotential. Usually, you will have 2
!! (L = 0, 1) and sometimes 3 (L = 2).
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
        function cl (itype)
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, pointer :: cl (:)             !< cl values
        integer, intent (in) :: itype       !< atom of interest

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer imu                  !< counting over L values
        integer iorb                 !< counting over orbitals
        integer issh
        integer Lvalue, Lmax         !< angular quntum numbers
        integer norb_PP_max
        
! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize to zero
        norb_PP_max = species(itype)%norb_PP_max
        allocate (cl(norb_PP_max))
        iorb = 0

! We now loop though all shells, and create cl for each orbital.  For example,
! sp^3 has two shells; cl(1) = cl_PP(0) and cl(2) = cl(3) = cl(4) = cl_PP(1).
        do issh = 1, species(itype)%nssh_PP
          Lvalue = species(itype)%shell_PP(issh)%lssh
          Lmax = (2*Lvalue + 1)
          do imu = 1, Lmax
            iorb = iorb + 1
            cl(iorb) = species(itype)%shell_PP(issh)%cl
          end do
        end do
        
! Sanity check.
        if (iorb .ne. species(itype)%norb_PP_max) then
          write (*,*) ' itype = ', itype
          write (*,*) ' index of orbitals for pseudopotential = ', iorb
          write (*,*) ' Program has norb_PP_max = ',                         &
     &                        species(itype)%norb_PP_max
          write (*,*) ' cl: index and norb_PP_max DO NOT agree. '
          stop
        end if

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

        return
        end


