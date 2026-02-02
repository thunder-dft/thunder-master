! copyright info:
!
!                             @Copyright 2025
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

! M_dynamics
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
        module M_dynamics
        use M_configuraciones

! Type Declaration
! ===========================================================================
! None

! Parameter Declaration and Description
! ===========================================================================
! Gear parameters - found by running set_gear

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! set_constraints
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
        subroutine set_constraints (s)
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
        integer iatom                       !< counter over the atoms
        integer in1                         !< species number
        integer inpfile                     !< reading from which unit
        integer ix                          !< counter over x, y, and z
        integer logfile                     !< writing to which unit

! random number generation
        integer :: seed_size, i
        integer, allocatable :: seed(:)
        integer :: clock(8)
        
        real xmass_total
        real vscale                         !< random number initial velocity

        character (len = 25) :: slogfile

        logical velocity

! Allocate Arrays
! ===========================================================================
        allocate (s%md)
     
! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! random number generation
        call random_seed (size=seed_size)
        allocate (seed(seed_size))

        call date_and_time (values=clock)
        seed = clock(8) + 37 * (/ (i-1, i=1,seed_size) /)

        call random_seed (put=seed)
        deallocate (seed)

! If a file called VELOCITIES exist, then read from the file.
! Read velocities from a velocities file. Note: if this is done, then it will
! wipe out the velocities originally initialized from a random temperature
! distribution.
        write (logfile,*)
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.VATOM'
        inquire (file = slogfile, exist = velocity)
        if (velocity) then
          write (logfile,*) ' We are reading from a velocity file. '
          open (unit = inpfile, file = slogfile, status = 'old')
          do iatom = 1, s%natoms
            read (inpfile,*) s%atom(iatom)%vatom
          end do
          close (unit = inpfile)

        ! Initialize atom velocities if not set by user.
        else
          write (logfile,*) ' We are setting random velocities. '
          call random_number (vscale)
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            do ix = 1, 3
              vscale = sqrt(-2.0d0*log(vscale))
              s%atom(iatom)%vatom(ix) = vscale*sqrt(P_fovermp*T_initial/   &
       &                                            (P_kconvert*species(in1)%xmass))
              call random_number (vscale)
              s%atom(iatom)%vatom(ix) =  s%atom(iatom)%vatom(ix)*cos(2.0d0*pi*vscale)
            end do
          end do
        end if

! Calculate the center-of-mass position.
        s%rcm = 0.0d0
        xmass_total = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass_total = xmass_total + species(in1)%xmass
          s%rcm = s%rcm + species(in1)%xmass*s%atom(iatom)%ratom
        end do
        s%rcm = s%rcm/xmass_total
        if (ishiftO .eq. 1) s%rcm = s%rcm - shifter
        write (logfile,100) s%rcm
        s%rcm_old = s%rcm

! Constraint #1
! ----------------
! Shift new ratom so they are measured from the center of mass.
! in this way, rcmmol = 0.
        write (logfile,*)
        if (iconstraint_rcm .eq. 1) then
          write (logfile,*) ' Constraining the positions about the center-of-mass. '
          do iatom = 1, s%natoms
            s%atom(iatom)%ratom = s%atom(iatom)%ratom - s%rcm
          end do
        else
          write (logfile,*) ' No constraining the positions about the center-of-mass. '
        end if

! Constraint #2
! -----------------
! Now adjust the velocities to get velocity of vcm = 0
        if (iconstraint_vcm .eq. 1) then
          write (logfile,*)
          write (logfile,*) ' Constraining the velocities about the center-of-mass. '
          s%vcm = 0.0d0
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            s%vcm = s%vcm + species(in1)%xmass*s%atom(iatom)%vatom
          end do
          s%vcm = s%vcm/xmass_total
          write (logfile,101) s%vcm

          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom - s%vcm
          end do

          ! recalculate vcm
          s%vcm = 0.0d0
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            s%vcm = s%vcm + species(in1)%xmass*s%atom(iatom)%vatom
          end do
          s%vcm = s%vcm/xmass_total
          write (logfile,102) s%vcm
        else
          write (logfile,*) ' No constraining the velocities about the center-of-mass. '
        end if

! Constraint #3
! -----------------
! Finally rescale the velocities to get the average temp = temperature_want
! tkinetic = average kinetic energy per particle in ev.
        if (iconstraint_KE .eq. 1) then
          write (logfile,*)
          write (logfile,*) ' Rescaling the velocities based on T_intial. '
          s%md%tkinetic = 0.0d0
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            s%md%tkinetic = s%md%tkinetic                                    &
     &        + (0.5d0/P_fovermp)*species(in1)%xmass                         &
     &         *(s%atom(iatom)%vatom(1)**2 + s%atom(iatom)%vatom(2)**2       &
     &                                     + s%atom(iatom)%vatom(3)**2)
          end do
          s%md%T_instantaneous = (2.0d0/3.0d0)*(s%md%tkinetic/s%natoms)*P_kconvert

! The temperature we now have (3/2 kb * T_instantaneous = tkinetic )
          write (logfile,*) ' T_initial = ', T_initial
          write (logfile,*) ' T_instantaneous, before rescaling = ', s%md%T_instantaneous
          if (s%md%T_instantaneous .gt. 0.0d0) then
            vscale = sqrt(T_initial/s%md%T_instantaneous)
          else
            vscale = 0.0d0
          end if
          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom*vscale
          end do

          ! check final temperature
          s%md%tkinetic = 0.0d0
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            s%md%tkinetic = s%md%tkinetic                                              &
     &        + (0.5d0/P_fovermp)*species(in1)%xmass                         &
     &         *(s%atom(iatom)%vatom(1)**2 + s%atom(iatom)%vatom(2)**2       &
     &                                     + s%atom(iatom)%vatom(3)**2)
          end do
          s%md%T_instantaneous = (2.0d0/3.0d0)*(s%md%tkinetic/s%natoms)*P_kconvert
          s%md%T_previous = s%md%T_instantaneous
          write (logfile,*) ' T_instantaneous, after rescaling = ', s%md%T_instantaneous
        else
          write (logfile,*) ' No rescaling the velocities based on T_intial. '
        end if

! Constraint #4
! ----------------
        if (iconstraint_L .eq. 1) then
          write (logfile,*)
          write (logfile,*) ' Constraining the angular momentum. '
          call zero_ang_mom (s)
        else
          write (logfile,*) ' No constraining the angular momentum. '
        end if

! Writeout the velocities
        write (logfile,*)
        write (logfile,*) ' Atom Velocities: '
        write (logfile,200)
        write (logfile,201)
        write (logfile,200)
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          write (logfile,202) iatom, species(in1)%symbol, s%atom(iatom)%vatom, in1
        end do

! Format Statements
! ===========================================================================
100     format (2x, ' Calculated position of the Center-of-Mass: ', 3(2x,f7.3))
101     format (2x, ' initial vcm = ', 3d16.7)
102     format (2x, '   final vcm = ', 3d16.7)

200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 6x,   &
     &              ' x ', 9x, ' y ', 9x, ' z ', 6x, ' Species # ')
202     format (3x, i5, 7x, a2, 3(2x,f10.5), 7x, i2)

! End Subroutine
! ===========================================================================
        return
        end subroutine set_constraints

! ===========================================================================
! zero_lin_mom
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       Contains subroutine zero_lin_mom which takes a set of random
! velocity and adjusts them to get total center mass momentum = 0.

! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine zero_lin_mom (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over the atoms
        integer in1                         !< species number
        integer logfile                     !< writing to which unit

        real xmass_total

! Procedure
! ===========================================================================
! Calculate the center-of-mass position.
        s%rcm = 0.0d0
        xmass_total = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass_total = xmass_total + species(in1)%xmass
          s%rcm = s%rcm + species(in1)%xmass*s%atom(iatom)%ratom
        end do
        s%rcm = s%rcm/xmass_total
        if (ishiftO .eq. 1) s%rcm = s%rcm - shifter
        write (logfile,100) s%rcm
        s%rcm_old = s%rcm

! Linear momentum in position  
        write (logfile,*)
        if (iconstraint_rcm .eq. 1) then
          write (logfile,*) ' Constraining the positions about the center-of-mass. '
          do iatom = 1, s%natoms
            s%atom(iatom)%ratom = s%atom(iatom)%ratom - s%rcm
          end do
        else
          write (logfile,*) ' No constraining the positions about the center-of-mass. '
        end if

! Linear momentum in velocity        
          s%vcm = 0.0d0
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            s%vcm = s%vcm + species(in1)%xmass*s%atom(iatom)%vatom
          end do
          s%vcm = s%vcm/xmass_total
          write (logfile,101) s%vcm

          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom - s%vcm
          end do

          ! recalculate vcm
          s%vcm = 0.0d0
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            s%vcm = s%vcm + species(in1)%xmass*s%atom(iatom)%vatom
          end do
          s%vcm = s%vcm/xmass_total
          write (logfile,102) s%vcm

! Format Statements
! ===========================================================================
100     format (2x, ' Calculated position of the Center-of-Mass: ', 3(2x,f7.3))
101     format (2x, ' initial vcm = ', 3d16.7)
102     format (2x, '   final vcm = ', 3d16.7)

! End Subroutine
! ===========================================================================
        return
        end subroutine zero_lin_mom

! ===========================================================================
! zero_ang_mom
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       Contains subroutine zero_ang_mom which takes a set of random
! velocity and adjusts them to get total angular momentum = 0.

! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine zero_ang_mom (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over the atoms
        integer in1                         !< species number
        integer ix                          !< counter over x, y, and z
        integer logfile                     !< writing to which unit

        real xmass

        real, dimension (3) :: crossa
        real, dimension (3, 3) :: xinertia

        real, dimension (3, 3) :: xinvert
        real, dimension (3) :: xlcm
        real, dimension (3) :: wvec

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        xlcm = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass
          xlcm(1) = xlcm(1) + xmass*(s%atom(iatom)%ratom(2)*s%atom(iatom)%vatom(3) - &
     &                               s%atom(iatom)%ratom(3)*s%atom(iatom)%vatom(2))
          xlcm(2) = xlcm(2) + xmass*(s%atom(iatom)%ratom(3)*s%atom(iatom)%vatom(1) - &
     &                               s%atom(iatom)%ratom(1)*s%atom(iatom)%vatom(3))
          xlcm(3) = xlcm(3) + xmass*(s%atom(iatom)%ratom(1)*s%atom(iatom)%vatom(2) - &
     &                               s%atom(iatom)%ratom(2)*s%atom(iatom)%vatom(1))
        end do
        write (logfile,100) xlcm

! Calculate the inertia tensor
! I(i,j) = sumoverk m(k) * ( r**2(k) delk(i,j)  -  ri(k) * rj(k) )
        xinertia = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass
          xinertia(1,1) = xinertia(1,1)                                     &
     &     + xmass*(s%atom(iatom)%ratom(2)**2 + s%atom(iatom)%ratom(3)**2)
          xinertia(1,2) = xinertia(1,2)                                     &
     &     - xmass*(s%atom(iatom)%ratom(1)*s%atom(iatom)%ratom(2))
          xinertia(1,3) = xinertia(1,3)                                     &
     &     - xmass*(s%atom(iatom)%ratom(1)*s%atom(iatom)%ratom(3))
          xinertia(2,1) = xinertia(2,1)                                     &
     &     - xmass*(s%atom(iatom)%ratom(2)*s%atom(iatom)%ratom(1))
          xinertia(2,2) = xinertia(2,2)                                     &
     &     + xmass*(s%atom(iatom)%ratom(1)**2 + s%atom(iatom)%ratom(3)**2)
          xinertia(2,3) = xinertia(2,3)                                     &

     &     - xmass*(s%atom(iatom)%ratom(2)*s%atom(iatom)%ratom(3))
          xinertia(3,1) = xinertia(3,1)                                     &
     &     - xmass*(s%atom(iatom)%ratom(3)*s%atom(iatom)%ratom(1))
          xinertia(3,2) = xinertia(3,2)                                     &
     &     - xmass*(s%atom(iatom)%ratom(3)*s%atom(iatom)%ratom(2))
          xinertia(3,3) = xinertia(3,3)                                     &
     &     + xmass*(s%atom(iatom)%ratom(1)**2 + s%atom(iatom)%ratom(2)**2)
        end do

! Here the inertia tensor has units of mp*A**2. This is okay
! because the units will work out in the end.
        call invert3x3 (xinertia, xinvert)

! L = I - dot - omega  so   omega = I(inverse) - dot - L
        wvec = 0.0d0
        do ix = 1, 3
          wvec(:) = wvec(:) + xinvert(:,ix)*xlcm(ix)
        end do

        do iatom = 1, s%natoms
          crossa(1) = wvec(2)*s%atom(iatom)%ratom(3) - wvec(3)*s%atom(iatom)%ratom(2)
          crossa(2) = wvec(3)*s%atom(iatom)%ratom(1) - wvec(1)*s%atom(iatom)%ratom(3)
          crossa(3) = wvec(1)*s%atom(iatom)%ratom(2) - wvec(2)*s%atom(iatom)%ratom(1)
          s%atom(iatom)%vatom = s%atom(iatom)%vatom - crossa
        end do

        xlcm = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass
          xlcm(1) = xlcm(1) + xmass*(s%atom(iatom)%ratom(2)*s%atom(iatom)%vatom(3) -  &
     &                               s%atom(iatom)%ratom(3)*s%atom(iatom)%vatom(2))
          xlcm(2) = xlcm(2) + xmass*(s%atom(iatom)%ratom(3)*s%atom(iatom)%vatom(1) -  &
     &                               s%atom(iatom)%ratom(1)*s%atom(iatom)%vatom(3))
          xlcm(3) = xlcm(3) + xmass*(s%atom(iatom)%ratom(1)*s%atom(iatom)%vatom(2) -  &
     &                               s%atom(iatom)%ratom(2)*s%atom(iatom)%vatom(1))
        end do
        write (logfile,101) xlcm

! Format Statements
! ===========================================================================
100     format (2x, ' initial Lcm = ', 3d16.7)
101     format (2x, '   final Lcm = ', 3d16.7)

! End Subroutine
! ===========================================================================
        return
        end subroutine zero_ang_mom


! ===========================================================================
! set_gear
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine is fake
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine set_gear ()
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Fake Subroutine

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine set_gear


! ===========================================================================
! md
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the driver for molecular-dynamics simulations.
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine md (s, itime_step)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

        integer itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
        integer logfile                      !< writing to which unit

        real vscale                          !< velocity rescaling

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*)
        write (logfile,*) ' Nuclear Evolution '
        
! velocity - Verlet
        if (itime_step .gt. 1) then
          write (s%logfile,*) ' Second half velocity update '
          if ( iensemble .eq. 0 ) then
            call vverlet_velocity (s)   
          else if ( iensemble .eq. 1 ) then  
            call langevin_velocity (s)
          end if
          call writeout_momentum (s)
          write(s%logfile, *) ' T_instantaneous ', s%md%T_instantaneous
        end if
        write (s%logfile,*) ' First half velocity update '
        write (s%logfile,*) ' And coordinate update '   
        if ( iensemble .eq. 0 ) then
          write (s%logfile, *) ' NVE NOW '
          call vverlet_velocity (s)
          call vverlet_coordinate (s)              
        else if ( iensemble .eq. 1 ) then
          write (s%logfile, *) ' Langevin NOW '
          do iatom = 1, s%natoms
            allocate (s%atom(iatom)%langevin)
          end do
          call langevin (s)
          call langevin_velocity (s)
          call langevin_coordinate (s)             
        end if

! Write out coordinates and velocity
        call writeout_coordinate_velocity (s)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine md

! ===========================================================================
! vverlet_coordinate
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the velocity verlet to update coordinate
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine vverlet_coordinate(s)

        implicit none
  
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer :: iatom	        ! counter of atoms

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================  
        do iatom = 1, s%natoms     
          s%atom(iatom)%ratom = s%atom(iatom)%ratom                          &
    &                          + s%atom(iatom)%vatom*dt
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return     
        end subroutine vverlet_coordinate


! ===========================================================================
! vverlet_velocity
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the velocity verlet to update velocity
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine vverlet_velocity(s)

        implicit none
  
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used


! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer :: iatom	        ! counter of atoms
        integer :: in1	                ! counter of species

        real :: xmass	                ! mass of atom

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================  
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass          
          s%atom(iatom)%vatom = s%atom(iatom)%vatom                          &
    &                          + 0.5d0*(P_fovermp*s%forces(iatom)%ftot/xmass)*dt
        end do

        ! Update kinetic energy and temperature after updating velocity              
        s%md%tkinetic = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass            
          s%md%tkinetic = s%md%tkinetic                                      &
   &                     + (0.5d0/P_fovermp)*species(in1)%xmass              &
   &                      *(s%atom(iatom)%vatom(1)**2                        &
   &                        + s%atom(iatom)%vatom(2)**2                      &
   &                        + s%atom(iatom)%vatom(3)**2)
        end do
        s%md%T_instantaneous = (2.0d0/3.0d0)*(s%md%tkinetic/s%natoms)*P_kconvert   
        s%md%T_previous = s%md%T_instantaneous

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return     
        end subroutine vverlet_velocity

! ===========================================================================
! langevin
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the velocity Verlet to update velocity
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine langevin (s)

        implicit none
  
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
        real, parameter :: friction = 0.01d0 / 0.009822694788464065d0

! Variable Declaration and Description
! ===========================================================================
        integer :: iatom	            ! counter of atoms
        integer :: in1	                ! counter of species
        integer :: ix	                ! counter of spatial dimension

        real :: xmass	                ! mass of atom
        real :: c1, c2, c3, c4, c5	    ! parameters for Langevin
        real :: sigma
        real :: u1, u2

        real, dimension (3) :: xi, eta
        real, dimension (3) :: ratom_random, vatom_random

        type(T_md), pointer :: pmd

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================  
        do iatom = 1, s%natoms

          ! Gaussian random numbers
          do ix = 1, 3
            call random_number(u1)
            call random_number(u2)
            if (u1 .lt. 1.0d-12) u1 = 1.0d-12   ! avoid log(0)
            xi(ix) = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * acos(-1.0d0) * u2)
            call random_number(u1)              
            call random_number(u2)
            if (u1 .lt. 1.0d-12) u1 = 1.0d-12   ! avoid log(0)
            eta(ix) = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * acos(-1.0d0) * u2)     
          end do

          ! Temperature in energy units      
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass                    
          sigma = sqrt(2.0d0*friction * P_fovermp*T_want/(P_kconvert*xmass))

          ! Coefficients (ASE Langevin)
          c1 = dt / 2.0d0 - dt * dt * friction / 8.0d0
          c2 = dt * friction / 2.0d0 - dt * dt * friction * friction / 8.0d0
          c3 = sqrt(dt) * sigma / 2.0d0 - dt**1.5d0 * friction * sigma / 8.0d0
          c5 = dt**1.5d0 * sigma / (2.0d0 * sqrt(3.0d0))
          c4 = friction / 2.0d0 * c5

          ! Random displacement & velocity
          ratom_random = c5 * eta
          vatom_random = c3 * xi - c4 * eta

          ! Cut off length
          nullify (pmd); pmd => s%atom(iatom)%langevin

          ! store variables
          pmd%c1 = c1; pmd%c2 = c2
          pmd%ratom_random = ratom_random
          pmd%vatom_random = vatom_random
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return     
        end subroutine langevin

! ===========================================================================
! langevin_coordinate
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the velocity verlet to update coordinate
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine langevin_coordinate(s)

        implicit none
  
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer :: iatom	        ! counter of atoms

        type(T_md), pointer :: pmd

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================  
        do iatom = 1, s%natoms  

          ! Cut some lengthy notation
          nullify (pmd); pmd => s%atom(iatom)%langevin   

          s%atom(iatom)%ratom = s%atom(iatom)%ratom                          &
&                              + dt*s%atom(iatom)%vatom                      &
&                              + pmd%ratom_random
        end do

! Format Statements
! ===========================================================================
! None


! End Subroutine
! ===========================================================================
        return     
        end subroutine langevin_coordinate

! ===========================================================================
! langevin_velocity
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the velocity verlet to update velocity for langevin
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine langevin_velocity (s)

        implicit none
  
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer :: iatom	        ! counter of atoms
        integer :: in1	                ! counter of species

        real :: xmass	                ! mass of atom

        type(T_md), pointer :: pmd
! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================  
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass

          ! Cut some lengthy notation
          nullify (pmd); pmd => s%atom(iatom)%langevin

          s%atom(iatom)%vatom = s%atom(iatom)%vatom  &
&                              + pmd%c1*P_fovermp*s%forces(iatom)%ftot/xmass &
&                              - pmd%c2*s%atom(iatom)%vatom                  &
&                              + pmd%vatom_random
        end do

        ! Update kinetic energy and temperature after updating velocity              
        s%md%tkinetic = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass            
          s%md%tkinetic = s%md%tkinetic                                      &
   &                 + (0.5d0/P_fovermp)*species(in1)%xmass                  &
   &                  *(s%atom(iatom)%vatom(1)**2 +  s%atom(iatom)%vatom(2)**2 &
   &                                              +  s%atom(iatom)%vatom(3)**2)
        end do
        s%md%T_instantaneous = (2.0d0/3.0d0)*(s%md%tkinetic/s%natoms)*P_kconvert   
        s%md%T_previous = s%md%T_instantaneous

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return     
        end subroutine langevin_velocity

! ===========================================================================
! berendson
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the berendson
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine berendson (s)

        implicit none
  
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used


! Parameters and Data Declaration
! ===========================================================================
        real, parameter :: tau = 100.0d0 * 0.009822694788464065

! Variable Declaration and Description
! ===========================================================================
        integer :: iatom	              ! counter of atoms

        real :: T_scale

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================  
        s%md%T_instantaneous = (2.0d0/3.0d0)*(s%md%tkinetic/s%natoms)*P_kconvert
        
        ! sqrt(1 + (T0/T - 1) * dt/tau)
        T_scale = sqrt(1.0d0 + (s%md%T_instantaneous/s%md%T_previous - 1.0d0) &
                              * dt/tau)
        
        ! limit scaling factor
        if (T_scale .gt. 1.1d0) T_scale = 1.1d0
        if (T_scale .lt. 0.9d0) T_scale = 0.9d0
        do iatom = 1, s%natoms
          s%atom(iatom)%vatom = s%atom(iatom)%vatom * T_scale
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return     
        end subroutine berendson


! ===========================================================================
! writeout_coordinate_velocity
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the driver for molecular-dynamics simulations.
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine writeout_coordinate_velocity (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
        integer logfile                      !< writing to which unit

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Now write out the basis file information.
        write (logfile,*)
        write (logfile,*) ' Atom Coordinates: '
        write (logfile,200)
        write (logfile,201)
        write (logfile,200)
        do iatom = 1, s%natoms
          if (ishiftO .eq. 1) then
            write (logfile,202) iatom, s%atom(iatom)%species%symbol,         &
     &                                 s%atom(iatom)%ratom - shifter,        &
     &                                 s%atom(iatom)%imass
          else
            write (logfile,202) iatom, s%atom(iatom)%species%symbol,         &
     &                                 s%atom(iatom)%ratom, s%atom(iatom)%imass
          end if
        end do
        write (logfile,200)

! Writeout the velocities
        write (logfile,*)
        write (logfile,*) ' Atom Velocities: '
        write (logfile,200)
        write (logfile,201)
        write (logfile,200)
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          write (logfile,202) iatom, species(in1)%symbol, s%atom(iatom)%vatom, in1
        end do

! Format Statements
! ===========================================================================
200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 6x,   &
     &              ' x ', 9x, ' y ', 9x, ' z ', 6x, ' Species # ')
202     format (3x, i5, 7x, a2, 3(2x,f10.5), 7x, i2)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_coordinate_velocity


! ===========================================================================
! writeout_momentum
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the subroutin writeout angular and linear momenta.
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
        subroutine writeout_momentum (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
        integer logfile

        real xmass
        real xmass_total

        real, dimension (3) :: xlcm

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Calculate the center of mass position and velocity.
        xmass_total = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass_total = xmass_total + species(in1)%xmass
          s%rcm = s%rcm + species(in1)%xmass*s%atom(iatom)%ratom
          s%vcm = s%vcm + species(in1)%xmass*s%atom(iatom)%vatom
        end do
        s%rcm = s%rcm/xmass_total
        if (ishiftO .eq. 1) s%rcm = s%rcm - shifter
        s%vcm = s%vcm/xmass_total

        write (logfile, 100) s%rcm
        write (logfile, 101) s%vcm

! Calculate the center of mass angular momentum.
        xlcm = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass
          xlcm(1) = xlcm(1) + xmass*(s%atom(iatom)%ratom(2)*s%atom(iatom)%vatom(3) - &
     &                               s%atom(iatom)%ratom(3)*s%atom(iatom)%vatom(2))
          xlcm(2) = xlcm(2) + xmass*(s%atom(iatom)%ratom(3)*s%atom(iatom)%vatom(1) - &
     &                               s%atom(iatom)%ratom(1)*s%atom(iatom)%vatom(3))
          xlcm(3) = xlcm(3) + xmass*(s%atom(iatom)%ratom(1)*s%atom(iatom)%vatom(2) - &
     &                               s%atom(iatom)%ratom(2)*s%atom(iatom)%vatom(1))
        end do
        write (logfile, 102) xlcm

! Deallocate Arrays
! ===========================================================================
!

! Format Statements
! ===========================================================================
100     format (2x, '         center of mass position = ', 3d12.4)
101     format (2x, '         center of mass velocity = ', 3d12.4)
102     format (2x, ' center of mass angular momentum = ', 3d12.4)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_momentum
        
        
! End Module
! ===========================================================================
        end module M_dynamics
