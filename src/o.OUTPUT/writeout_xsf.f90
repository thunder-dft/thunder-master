! copyright info:
!
!                             @Copyright 2014
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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

! writeout_xsf.f90
! Program Description
! ===========================================================================
!       The subroutine writes out data in xsf-format file used by xcrysden
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine writeout_xsf (t, xsf, xsfname, message)
        use M_configuraciones
        use M_species
        use M_grid
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: t           !< the structure to be used
        real, dimension (:), pointer :: xsf
        character (len=25), intent (in) :: xsfname
        character (len=25), intent (in) :: message

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over atoms
        integer i, j, k                     !< counters ofer grid points
        integer i0, j0, k0                  !< grid index points
        integer in1                         !< species numbers
        integer index0                      !< index point
        integer inpfile                     !< reading from which unit
        integer logfile                     !< writing to which unit
        !real, dimension (:,:), allocatable :: ratom2g

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile
        inpfile = t%inpfile

! write out grid array information into *.xsf file
! (format of xcrysden visual code - for details see www.xcrysden.org)
        write (logfile, *) '  Write out xsf file ', xsfname
        open (unit = inpfile, file = xsfname, status = 'unknown' )
        !allocate (ratom2g (3, t%natoms))


! print the list of atoms
        if (t%icluster .eq. 1) then
          write (inpfile,*) 'ATOMS'
          do iatom = 1,t%natoms
            in1 = t%atom(iatom)%imass
            write (inpfile,'(i2,3f14.8)') species(in1)%nz, ratom2g(:,iatom)
          end do
        else
          write (inpfile,*) 'CRYSTAL'
          write (inpfile,*) 'PRIMVEC'
          write (inpfile,*) t%lattice(1)%a
          write (inpfile,*) t%lattice(2)%a
          write (inpfile,*) t%lattice(3)%a

          write (inpfile,*) 'PRIMCOORD'
          write (inpfile,*) t%natoms, 1
          do iatom = 1, t%natoms
            in1 = t%atom(iatom)%imass
            write (inpfile, 100) species(in1)%nz, ratom2g(:,iatom)
          end do
        end if

        write (inpfile,*)
        write (inpfile,*) 'BEGIN_BLOCK_DATAGRID_3D'
        write (inpfile,*) message
        write (inpfile,*) 'DATAGRID_3D_DENSITY'
        write (inpfile,*) irm1 + 1, irm2 + 1, irm3 + 1
        write (inpfile,*) 0.0d0, 0.0d0, 0.0d0

        ! print lattice vector
        write (inpfile,*) xsfgrid(1)%a
        write (inpfile,*) xsfgrid(2)%a
        write (inpfile,*) xsfgrid(3)%a

! print values of the grid point
        do k = 0, irm3
          if (k .eq. irm3) then
            k0 = 0
          else
            k0 = k
          end if
          do j = 0, irm2
            if (j .eq. irm2) then
              j0 = 0
            else
              j0 = j
            end if
            do i = 0, irm1
              if (i .eq. irm1) then
                i0 = 0
              else
                i0 = i
              end if

! mapping index within the regular mesh
              index0 = i0 + irm1*j0 + irm1*irm2*k0
              write (inpfile,200) xsf(index0)
            end do ! do i
          end do ! do j
        end do ! do k
        write (inpfile,*) 'END_DATAGRID_3D'
        write (inpfile,*) 'END_BLOCK_DATAGRID_3D'

! close file
        close (inpfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (i2, 3f14.8)
200     format (e14.6)

        return
        end subroutine writeout_xsf
