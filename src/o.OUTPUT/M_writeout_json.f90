! copyright info:
!
!                             @Copyright 2022
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

! M_MD_setting
! Module Description
! ===========================================================================
!>       This is a module containing all molecular-dynamics algorithms for
! doing the predictor-corrector Gear algorithm.
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
!
! Module Declaration
! ===========================================================================
        module M_writeout_file
        use M_configuraciones

! Type Declaration
! ===========================================================================
! None

! Parameter Declaration and Description
! ===========================================================================
! None


! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains    


! ===========================================================================
! openfile
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine openfile (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        character (len = 25) :: slogfile
        character (len = 25) :: sjsonfile

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        sjsonfile = trim(slogfile)//'.json'
        open (unit = s%jsonfile, file = sjsonfile, status = 'replace')


! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine openfile        


! ===========================================================================
! writeout_file_head
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_head (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
        write (s%jsonfile,'(A)') '{"fireball":['


! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_head        

! ===========================================================================
! writeout_file_step_head
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_step_head (s, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
        write (s%jsonfile,'(A)') '{'
        write (s%jsonfile,'(6x, A, I5)') '"nstep":', itime_step

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_step_head       

! ===========================================================================
! writeout_file_cell
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_cell (s, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer idi
        integer idj

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items
        write(s%jsonfile,'(6x, A)') ','

        write (s%jsonfile,'(6x, A)') '"cell":['

        do idi = 1, 3
           write(s%jsonfile,'(6x, A)', advance='no') '['

           do idj = 1, 3
              write(s%jsonfile,'(3x, F15.6)', advance='no') s%lattice(idi)%a(idj)
              if (idj < 3) write(s%jsonfile,'(A)', advance='no') ','
           end do

           write(s%jsonfile,'(A)', advance='no') ']'
           if (idi < 3) write(s%jsonfile,'(A)', advance='no') ','
           write(s%jsonfile,*)
        end do

        write(s%jsonfile,'(6x, A)') ']'


! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_cell       

! ===========================================================================
! writeout_file_atoms_number
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_atoms_number (s, itime_step)
        
        use M_species

        implicit none
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                     !< counter over atoms and neighbors
        integer in1

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items
        write(s%jsonfile,'(6x, A)') ','     

        write (s%jsonfile,'(6x, A)') '"numbers":['
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          write (s%jsonfile,'(16x, i3)', advance='no') species(in1)%nZ
          if (iatom < s%natoms) write(s%jsonfile,'(A)', advance='no') ','
          write(s%jsonfile,*)                            
        end do
        write(s%jsonfile,'(6x, A)') ']'


! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_atoms_number               

! ===========================================================================
! writeout_file_position
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_position (s, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                     !< counter over atoms
        integer idx                        !< counter over spatial dimension

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items
        write(s%jsonfile,'(6x, A)') ','            
        
        write (s%jsonfile,'(6x, A)') '"positions":['
        do iatom = 1, s%natoms
          write(s%jsonfile,'(6x, A)', advance='no') '['
          do idx = 1, 3
            write(s%jsonfile,'(3x, F15.6)', advance='no') &            
     &            s%atom(iatom)%ratom(idx)     
            if (idx < 3) write(s%jsonfile,'(A)', advance='no') ','
          end do
          write(s%jsonfile,'(A)', advance='no') ']'
          if (iatom < s%natoms) write(s%jsonfile,'(A)', advance='no') ','
          write(s%jsonfile,*)     
        end do
        write(s%jsonfile,'(6x, A)') ']'


! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_position                       

! ===========================================================================
! writeout_file_charges
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_charges (s, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                     !< counter over atoms and neighbors
        integer in1
        integer nssh                      !< number of shells
        integer issh                      !< counter over shells

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items
        write(s%jsonfile,'(6x, A)') ','            
        
        write (s%jsonfile,'(6x, A)') '"charges":['
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          nssh = species(in1)%nssh
          write(s%jsonfile,'(6x, A)', advance='no') '['   

          do issh = 1, nssh
            write(s%jsonfile,'(3x, F15.6)', advance='no') &            
     &              s%atom(iatom)%shell(issh)%Qin     
            if (issh < nssh) write(s%jsonfile,'(A)', advance='no') ','
          end do

          write(s%jsonfile,'(A)', advance='no') ']'
          if (iatom < s%natoms) write(s%jsonfile,'(A)', advance='no') ','
          write(s%jsonfile,*)    

        end do
        write(s%jsonfile,'(6x, A)') ']'          


! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_charges      

! ===========================================================================
! writeout_file_efermi
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_efermi (s, efermi, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        real, intent(in) :: efermi
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items
        write(s%jsonfile,'(6x, A)') ','              
        write (s%jsonfile,'(6x, A, F15.6)') '"fermi":', efermi

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_efermi


! ===========================================================================
! writeout_file_etot
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_etot (s, etot, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        real, intent(in) :: etot
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items
        write(s%jsonfile,'(6x, A)') ','              
        write (s%jsonfile,'(6x, A, F15.6)') '"energy":', etot

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_etot   


! ===========================================================================
! writeout_file_force
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_force (s, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                     !< counter over atoms and neighbors
        integer idx                       !< counter over spatial dimension

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items
        write(s%jsonfile,'(6x, A)') ','                     
        
        write (s%jsonfile,'(6x, A)') '"forces":['
        do iatom = 1, s%natoms
          write(s%jsonfile,'(6x, A)', advance='no') '['
          do idx = 1, 3
            write(s%jsonfile,'(3x, F15.6)', advance='no') &            
     &            s%forces(iatom)%ftot(idx)     
            if (idx < 3) write(s%jsonfile,'(A)', advance='no') ','
          end do
          write(s%jsonfile,'(A)', advance='no') ']'
          if (iatom < s%natoms) write(s%jsonfile,'(A)', advance='no') ','
          write(s%jsonfile,*)     
        end do
        write(s%jsonfile,'(6x, A)') ']'            

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_force   


! ===========================================================================
! writeout_file_rms
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_rms (s, rms, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        real, intent(in) :: rms
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items    
        write(s%jsonfile,'(6x, A)') ','  
        write (s%jsonfile,'(6x, A, F15.6)') '"RMS":', rms          

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_rms   


! ===========================================================================
! writeout_file_dij_nac
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_dij_nac (s, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                     !< counter over atoms and neighbors
        integer in1
        integer ikpoint
        integer iband
        integer jband
        integer idx

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
! Add comma to speparate above items
        write(s%jsonfile,'(6x, A)') ','                             
        
        write(s%jsonfile,'(6x, A)') '"dij": ['

        do ikpoint = 1, s%nkpoints
          write(s%jsonfile,'(6x, A)') '['
        
          do iband = 1, s%kpoints(ikpoint)%nbands - 1
            write(s%jsonfile,'(6x, A)') '['
        
            do jband = iband + 1, s%kpoints(ikpoint)%nbands
              write(s%jsonfile,'(6x, A)') '['

              do iatom = 1, s%natoms
                write(s%jsonfile,'(6x, A)', advance='no') '['
                do idx = 1, 3
                  write(s%jsonfile,'(3x, F15.6)', advance='no') &            
     &                         s%kpoints(ikpoint)%transition(iband)%dij(idx,iatom,jband)     
                  if (idx < 3) write(s%jsonfile,'(A)', advance='no') ','
                end do
                write(s%jsonfile,'(A)', advance='no') ']'
                if (iatom < s%natoms) write(s%jsonfile,'(A)', advance='no') ','
                write(s%jsonfile,*)     
              end do
        
              write(s%jsonfile,'(6x, A)', advance='no') ']'
              if (jband < s%kpoints(ikpoint)%nbands) write(s%jsonfile,'(A)', advance='no') ','
              write(s%jsonfile,*)
            end do
        
            write(s%jsonfile,'(6x, A)', advance='no') ']'
            if (iband < s%kpoints(ikpoint)%nbands - 1) write(s%jsonfile,'(A)', advance='no') ','
            write(s%jsonfile,*)
          end do
        
          write(s%jsonfile,'(6x, A)', advance='no') ']'
          if (ikpoint < s%nkpoints) write(s%jsonfile,'(A)', advance='no') ','
          write(s%jsonfile,*)
        
        end do
        
        write(s%jsonfile,'(A)') '      ]'

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_dij_nac 

! ===========================================================================
! writeout_file_step_end
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_step_end (s, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used
        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
        write(s%jsonfile,'(A)', advance='no') '}'
        if (itime_step < nstepf) write(s%jsonfile,'(A)', advance='no') ','

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_step_end      

! ===========================================================================
! writeout_file_end 
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the following constraints - center-of-mass
! coordinates, center-of-mass velocities, kinetic energy rescaling of
! velocities, and angular momentum.
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
        subroutine writeout_file_end (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
!

! Procedure
! ===========================================================================
        write (s%jsonfile,'(A)') ']}'
        close(s%jsonfile)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_file_end        

 ! End Module
! ===========================================================================
        end module M_writeout_file
       