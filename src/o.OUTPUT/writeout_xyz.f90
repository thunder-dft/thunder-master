! copyright info:
!
!                             @Copyright 2013
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
! Ohio University - Dave Drabold
! University of Texas at Austin - Alex Demkov
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

! writeout_xyz.f90
! Subroutine Description
! ===========================================================================
!       This routine writes out the xyz file.
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
        subroutine writeout_xyz (t, ebs, uii_uee, uxcdcc)
        use M_species
        use M_configuraciones
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: t             ! the structure to be used

        real, intent (in) :: ebs                   ! band-structure energy
        real, intent (in) :: uii_uee, uxcdcc       ! short-range energies

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                      ! counter over atoms and neighbors
        integer in1                        ! species numbers
        integer inpfile                    ! reading from which unit

        real etot                          ! total energy

        character (len = 25) :: slogfile

! Procedure
! ===========================================================================
! Initialize logfile
        inpfile = t%inpfile

        ! Evaluate total energy
        etot = ebs + uii_uee + uxcdcc

! Open the xyz output file for this structure
! If the .xyz file exists, then this is a possible restart.
! Inquire about the file and assume that it is a restart.
! Write to logfile that this is a restart so that the user knows to remove
! the .xyz file if need be.
        slogfile = t%basisfile(:len(trim(t%basisfile))-4)
        slogfile = trim(slogfile)//'.xyz'
        open (unit = inpfile, file = slogfile, status = 'unknown', position = 'append')

! Loop over the atoms in the central cell.
        write (inpfile,*) t%natoms
        write (inpfile,*) etot
        do iatom = 1, t%natoms
          in1 = t%atom(iatom)%imass
          write (inpfile,101) species(in1)%symbol, t%atom(iatom)%ratom

! Finish loop over atoms.
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (2x, a2, 3(2x,f12.6))

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_xyz
