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
! Synfuels China Technology Co., Ltd. - Zhaofa Li
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
        module M_dynamics_settings
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

! Initialize the center-of-mass position, velocity and angular momentum.
        call writeout_momentum (s) 
! Initialize the kinetic energy and the instantaneous temperature.
        call writeout_kin_T (s)
        s%md%T_previous = s%md%T_instantaneous

! Constraint #1
! ----------------
! Shift new ratom so they are measured from the center of mass.
! in this way, rcmmol = 0.
        write (logfile,*)
        if (iconstraint_rcm .eq. 1) then
          write (logfile,*) ' Constraining the positions about the center-of-mass. '
          do iatom = 1, s%natoms
            s%atom(iatom)%ratom = s%atom(iatom)%ratom - s%md%rcm
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
          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom - s%md%vcm
          end do
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
! Recalculate the kinetic energy and the instantaneous temperature after rescaling.
          call writeout_kin_T (s)
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
! Recalculate the center-of-mass position, velocity and angular momentum after all constraints are applied.
        call writeout_momentum (s)

! Writeout the coordinates and velocities
        call writeout_coodinate_velocity (s)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine set_constraints


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
        real, dimension (3) :: wvec

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

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
          wvec(:) = wvec(:) + xinvert(:,ix)*s%md%xlcm(ix)
        end do

        do iatom = 1, s%natoms
          crossa(1) = wvec(2)*s%atom(iatom)%ratom(3) - wvec(3)*s%atom(iatom)%ratom(2)
          crossa(2) = wvec(3)*s%atom(iatom)%ratom(1) - wvec(1)*s%atom(iatom)%ratom(3)
          crossa(3) = wvec(1)*s%atom(iatom)%ratom(2) - wvec(2)*s%atom(iatom)%ratom(1)
          s%atom(iatom)%vatom = s%atom(iatom)%vatom - crossa
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine zero_ang_mom


! ===========================================================================
! writeout_kin_T
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
        subroutine writeout_kin_T (s)
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

! Calculate the kinetic energy and the instantaneous temperature
        s%md%tkinetic = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          s%md%tkinetic = s%md%tkinetic                                                 &
     &       + (0.5d0/P_fovermp)*species(in1)%xmass                            &
     &        *dot_product(s%atom(iatom)%vatom, s%atom(iatom)%vatom)
        end do
        s%md%T_instantaneous = (2.0d0/3.0d0)*(s%md%tkinetic/s%natoms)*P_kconvert
        write (logfile, 101) s%md%tkinetic
        write (logfile, 102) s%md%T_instantaneous

! Format Statements
! ===========================================================================
101     format (2x, '    Kinetics = ', 3d16.7)
102     format (2x, ' Temperature = ', 3d16.7)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_kin_T

! ===========================================================================
! writeout_coodinate_velocity
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
        subroutine writeout_coodinate_velocity (s)
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
        end subroutine writeout_coodinate_velocity


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

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Calculate the center of mass position and velocity.
        s%md%rcm = 0.0d0; s%md%vcm = 0.0d0; s%md%xlcm = 0.0d0
        xmass_total = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass_total = xmass_total + species(in1)%xmass

          s%md%rcm = s%md%rcm + species(in1)%xmass*s%atom(iatom)%ratom
          s%md%vcm = s%md%vcm + species(in1)%xmass*s%atom(iatom)%vatom

          s%md%xlcm(1) = s%md%xlcm(1) + xmass*(s%atom(iatom)%ratom(2)*s%atom(iatom)%vatom(3) - &
     &                               s%atom(iatom)%ratom(3)*s%atom(iatom)%vatom(2))
          s%md%xlcm(2) = s%md%xlcm(2) + xmass*(s%atom(iatom)%ratom(3)*s%atom(iatom)%vatom(1) - &
     &                               s%atom(iatom)%ratom(1)*s%atom(iatom)%vatom(3))
          s%md%xlcm(3) = s%md%xlcm(3) + xmass*(s%atom(iatom)%ratom(1)*s%atom(iatom)%vatom(2) - &
     &                               s%atom(iatom)%ratom(2)*s%atom(iatom)%vatom(1))          
        end do
        s%md%rcm = s%md%rcm/xmass_total
        if (ishiftO .eq. 1) s%md%rcm = s%md%rcm - shifter
        s%md%vcm = s%md%vcm/xmass_total

        write (logfile, 100) s%md%rcm
        write (logfile, 101) s%md%vcm
        write (logfile, 102) s%md%xlcm

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


! destroy_MD_setting
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing MD settings
!>       information.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Subroutine Declaration
! ===========================================================================
        subroutine destroy_settings (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        deallocate (s%md)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_settings
        
! End Module
! ===========================================================================
        end module M_dynamics_settings
