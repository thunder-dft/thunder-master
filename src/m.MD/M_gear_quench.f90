! copyright info:
!
!                             @Copyright 2012
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

! M_dynamics
! Module Description
! ===========================================================================
!>       This is a module containing all molecular-dynamics algorithms for
! doing the predictor-corrector Gear algorithm.
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
        integer, parameter :: ngear = 5

        real, dimension (0:7) :: cfactor

! Variable Declaration and Description
! ===========================================================================
        real tkinetic                        !< kinetic energy
        real T_average                       !< average temperature
        real T_instantaneous                 !< instantaneous temperature
        real T_previous                      !< previous temperature

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
!> @author James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
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

        real xmass_total
        real tkinetic                       !< kinetic energy
        real vscale                         !< random number initial velocity

        character (len = 25) :: slogfile

        logical velocity

! Allocate Arrays
! ===========================================================================
        do iatom = 1, s%natoms
          allocate (s%atom(iatom)%xdot(0:ngear,3)); s%atom(iatom)%xdot = 0.0d0
        end do

! Procedure
! ===========================================================================

! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! If a file called VELOCITIES exist, then read from the file.
! Read velocities from a velocities file. Note: if this is done, then it will
! wipe out the velocities originally initialized from a random temperature
! distribution.
        write (logfile,*)
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.vatom'
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
            s%atom(iatom)%vatom = 0.0d0
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
        write (logfile,100) s%rcm
        s%rcm_old = s%rcm

! Constraint #1
! ----------------
! Shift new ratom so they are measured from the center of mass.
! in this way, rcmmol = 0.
        if (iconstraint_rcm .eq. 1) then
          write (logfile,*) ' Constraining the positions about the center-of-mass. '
          do iatom = 1, s%natoms
            s%atom(iatom)%ratom = s%atom(iatom)%ratom - s%rcm
          end do
        end if

! Constraint #4
! ----------------
        if (iconstraint_L .eq. 1) then
          write (logfile,*) ' Constraining angular momentum. '
          call zero_ang_mom (s)
        end if

! Constraint #2
! -----------------
! Now adjust the velocities to get velocity of vcm = 0
        if (iconstraint_vcm .eq. 1) then
          write (logfile,*) ' Constraining the velocities about the center-of-mass. '
          s%vcm = 0.0d0
          do iatom = 1, s%natoms
            s%vcm = s%vcm + species(in1)%xmass*s%atom(iatom)%vatom
          end do
          s%vcm = s%vcm/xmass_total
          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom - s%vcm
          end do
        end if

! Constraint #3
! -----------------
! Finally rescale the velocities to get the average temp = temperature_want
! tkinetic = average kinetic energy per particle in ev.
        if (iconstraint_KE .eq. 1) then
          write (logfile,*) ' Rescaling the velocities based on desired temperature. '
          tkinetic = 0.0d0
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            tkinetic = tkinetic                                              &
     &        + 0.5d0*species(in1)%xmass/P_fovermp                           &
     &               *(s%atom(iatom)%vatom(1)**2 + s%atom(iatom)%vatom(2)**2 &
     &                                           + s%atom(iatom)%vatom(3)**2)
          end do
          tkinetic = tkinetic/(s%natoms - s%nfragments)

! The temperature we now have (3/2 kb * T_instantaneous = tkinetic )
          T_instantaneous = (2.0d0/3.0d0)*tkinetic*P_kconvert
          if (T_instantaneous .gt. 0.0d0) then
            vscale = sqrt(T_initial/T_instantaneous)
          else
            vscale = 0.0d0
          end if
          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom*vscale
          end do
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
200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 6x,   &
     &              ' x ', 9x, ' y ', 9x, ' z ', 6x, ' Species # ')
202     format (3x, i5, 7x, a2, 3(2x,f10.5), 7x, i2)

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
        write (logfile,200) xlcm

! Format Statements
! ===========================================================================
100     format (2x, ' initial Lcm = ', 3d16.7)
200     format (2x, '   final Lcm = ', 3d16.7)

! End Subroutine
! ===========================================================================
        return
        end subroutine zero_ang_mom


! ===========================================================================
! set_gear
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine establishes the gear algorithm parameters.
!
! Gear proceedures use the Gear algorithm to move the atoms
!
! Gear algorithm --- use a Taylor series expansion up to fifth-order
!                   in the derivatives.
!
!                   Prediction step:
!                     x(t+dt)=x(t)+dx(t)*dt+(1/2)*ddx(t)*dt**2+.....
!                    dx(t+dt)=dx(t)+ddx(t)*dt+(1/2)*dddx(t)*dt**2+....
!                   ddx(t+dt)=ddx(t)+dddx(t)*dt+(1/2)*ddddx(t)*dt**2+...
!                  dddx(t+dt)=dddx(t)+ddddx(t)*dt+(1/2)*dddddx(t)*dt**2
!                 ddddx(t+dt)=ddddx(t)+dddddx(t)*dt
!                dddddx(t+dt)=dddddx(t)
!
!                    (here dx(t) means first derivative of x wrt t,
!                       ddx(t) means second derivative of x wrt t,etc.)
!
!                   Correction step:
!                      First calculate the force at the predicted
!                        position x(t+dt) from above.  From the force
!                        calculate the acceleration a.  Then take
!                        the difference between the predicted and
!                        actual accelerations:
!
!                                D = a - ddx(t+dt)
!
!                      Now correct all the derivatives:
!                        xnew(t+dt) = x(t+dt)      + c0*D*(0.5*dt**2)
!                       dxnew(t+dt) = dx(t+dt)     + c1*D*(0.5*dt)
!                      ddxnew(t+dt) = ddx(t+dt)    + c2*D
!                     dddxnew(t+dt) = dddx(t+dt)   + c3*D*(3/dt)
!                    ddddxnew(t+dt) = ddddx(t+dt)  + c4*D*(12/dt**2)
!                   dddddxnew(t+dt) = dddddx(t+dt) + c5*D*(60/dt**3)
!
!                     References:
!                        C.W. Gear, ANL Report No. ANL-7126, 1966.
!                        Young Hee Lee, Ph.D. Dissertation, Kent
!                            State University, August 1986.
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
! Gear coeficients for second-order equation (thus cfactor(2) = 1.0d0 always)
        cfactor = 0.0d0
        if (ngear .eq. 2) then   ! velocity verlet style
          cfactor(0) = 0.0d0
          cfactor(1) = 1.0d0
          cfactor(2) = 1.0d0
        else if (ngear .eq. 3) then
          cfactor(0) = 1.0d0/6.0d0
          cfactor(1) = 5.0d0/6.0d0
          cfactor(2) = 1.0d0
          cfactor(3) = 1.0d0/3.0d0
        else if (ngear .eq. 4) then
          cfactor(0) = 19.0d0/120.0d0
          cfactor(1) = 3.0d0/4.0d0
          cfactor(2) = 1.0d0
          cfactor(3) = 1.0d0/2.0d0
          cfactor(4) = 1.0d0/12.0d0
        else if (ngear .eq. 5) then
          cfactor(0) = 3.0d0/20.0d0
          cfactor(1) = 251.0d0/360.0d0
          cfactor(2) = 1.0d0
          cfactor(3) = 11.0d0/18.0d0
          cfactor(4) = 1.0d0/6.0d0
          cfactor(5) = 1.0d0/60.0d0
        else if (ngear .eq. 6) then
          cfactor(0) = 863.0d0/6048.0d0
          cfactor(1) = 665.0d0/1008.0d0
          cfactor(2) = 1.0d0
          cfactor(3) = 25.0d0/36.0d0
          cfactor(4) = 35.0d0/144.0d0
          cfactor(5) = 1.0d0/24.0d0
          cfactor(6) = 1.0d0/360.0d0
        else if (ngear .eq. 7) then
          cfactor(0) = 1925.0d0/14112.0d0
          cfactor(1) = 19087.0d0/30240.0d0
          cfactor(2) = 1.0d0
          cfactor(3) = 137.0d0/180.0d0
          cfactor(4) = 5.0d0/16.0d0
          cfactor(5) = 17.0d0/240.0d0
          cfactor(6) = 1.0d0/120.0d0
          cfactor(7) = 1.0d0/2520.0d0
        end if

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
        integer logfile                     !< writing to which unit

        real vscale                          !< velocity rescaling
        real xmass

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! Corrector
        write (logfile,*)
        write (logfile,*) ' Predictor-Corrector: correct the positions. '
        call corrector (s, itime_step)
        do iatom = 1, s%natoms
          s%atom(iatom)%ratom = s%atom(iatom)%xdot(0,:)
          s%atom(iatom)%vatom = s%atom(iatom)%xdot(1,:)
        end do

! Calculate the kinetic energy and the instantaneous temperature
        tkinetic = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass
          tkinetic = tkinetic + (0.5d0/P_fovermp)*xmass                      &
     &      *(s%atom(iatom)%vatom(1)**2 + s%atom(iatom)%vatom(2)**2          &
     &        + s%atom(iatom)%vatom(3)**2)
        end do
        T_instantaneous = (2.0d0/3.0d0)*tkinetic*P_kconvert/s%natoms
        write (logfile,*) ' T_instantaneous = ', T_instantaneous

! Rescale the temperature if iensemble = 1 (constant temperature MD)
        if (iensemble .eq. 1 .and. .not. T_instantaneous .le. 0) then
          vscale = sqrt(T_want/T_instantaneous)
          write (logfile,*) ' Constant temperature ensemble (iensemble = 1): scaling velocities. '
          write (logfile,*) ' Scaling = ', vscale, ' It should be near 1.0!'
          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom*vscale
            s%atom(iatom)%xdot(1,:) = s%atom(iatom)%xdot(1,:)*vscale
          end do
        end if
        T_average = ((itime_step - 1)*T_average + T_instantaneous)/itime_step

! Predictor
        write (logfile,*) ' Predictor-Corrector: predict the positions. '
        call predictor (s, itime_step)

! Scale velocities again after predictor step
        if (iensemble .eq. 1) then
          tkinetic = 0.0d0
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            xmass = species(in1)%xmass
            tkinetic = tkinetic + (0.5d0/P_fovermp)*xmass                      &
     &        *(s%atom(iatom)%vatom(1)**2 + s%atom(iatom)%vatom(2)**2          &
     &          + s%atom(iatom)%vatom(3)**2)
          end do
          T_instantaneous = (2.0d0/3.0d0)*tkinetic*P_kconvert/s%natoms
          vscale = sqrt(T_want/T_instantaneous)
          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom*vscale
            s%atom(iatom)%xdot(1,:) = s%atom(iatom)%xdot(1,:)*vscale
          end do
          if (vscale .gt. 1.3d0 .or. vscale .lt. 0.7d0) then
            write (logfile,*)
            write (logfile,*) ' Warning: significant velocity scaling needed after '
            write (logfile,*) ' predictor. Scaling factor = ', vscale
            write (logfile,*)
          end if
        end if

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine md


! ===========================================================================
! corrector
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the corrector part for the positions and time derivatives.
! We use a 5th order Gear algorithm here.
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
        subroutine corrector (s, itime_step)
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
        integer iorder

        real dtfactor
        real xmass

        real, dimension (3) :: acceleration
        real, dimension (:, :), allocatable :: difference

        interface
          function factorial (i)
            integer, intent (in) :: i
            integer factorial
          end function factorial
        end interface

! Allocate Arrays
! ===========================================================================
        allocate (difference (3, s%natoms))

! Procedure
! ===========================================================================
! First calculate the accelerations at the predicted points and the
! differences between the predicted and actual accelerations
        if (itime_step .eq. 1) then
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            xmass = species(in1)%xmass
            s%atom(iatom)%xdot(0,:) = s%atom(iatom)%ratom
            s%atom(iatom)%xdot(1,:) = s%atom(iatom)%vatom
            s%atom(iatom)%xdot(2,:) = P_fovermp*s%forces(iatom)%ftot/xmass
          end do
        else
          do iatom = 1, s%natoms
            in1 = s%atom(iatom)%imass
            xmass = species(in1)%xmass
            acceleration = P_fovermp*s%forces(iatom)%ftot/xmass
            difference(:,iatom) = acceleration - s%atom(iatom)%xdot(2,:)
          end do
        end if

! Gear (often fifth-order)
        do iatom = 1, s%natoms
          do iorder = 0, ngear
            if (iorder .eq. 2) then
              dtfactor = 1.0d0
            else
              dtfactor = dt**(2-iorder)
            end if
              s%atom(iatom)%xdot(iorder,:) = s%atom(iatom)%xdot(iorder,:)    &
     &         + cfactor(iorder)*difference(:,iatom)*(factorial(iorder)/2.0d0)*dtfactor
          end do
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine corrector


! ===========================================================================
! predictor
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the predictor part for the positions and time derivatives.
! We use a 5th order Gear algorithm here.
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
        subroutine predictor (s, itime_step)
        implicit none

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
        integer ifactor, iorder, isum       !< for Gear loop
        integer ix                          !< loop over spatial coordinates
        integer kquench                     !< query if quench or not
        integer logfile                     !< writing to which unit

        real dot                             !< calculate power
        real rfactor                         !< annealing factor
        real xmass
        real xmass_total

        real, dimension (3) :: rcm, vcm, xlcm

        interface
          function factorial (i)
            integer, intent (in) :: i
            integer factorial
          end function factorial
        end interface

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
          rcm = rcm + species(in1)%xmass*s%atom(iatom)%ratom
          vcm = vcm + species(in1)%xmass*s%atom(iatom)%vatom
        end do
        rcm = rcm/xmass_total
        vcm = vcm/xmass_total

        write (logfile, *)
        write (logfile, 100) rcm
        write (logfile, 101) vcm

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

! Move the atoms
! For the gear algorithm update both the positions and the velocities after
! each prediction and correction
! Actual prediction step - usually 5th order Gear algorithm.
        do iatom = 1, s%natoms
          do iorder = 0, ngear - 1
            do isum = iorder + 1, ngear
              ifactor = isum - iorder
              s%atom(iatom)%xdot(iorder,:) = s%atom(iatom)%xdot(iorder,:)    &
     &         + s%atom(iatom)%xdot(isum,:)*(dt**ifactor)/factorial(ifactor)
            end do
          end do
        end do

! Project out some forces if using FRAGMENTS
!       if (numfrags .ne. 0) then
!         call fixfrags (xmass, natoms)
!       end if
        do iatom = 1, s%natoms
          s%atom(iatom)%ratom = s%atom(iatom)%xdot(0,:)
          s%atom(iatom)%vatom = s%atom(iatom)%xdot(1,:)
        end do

        if (iquench .eq. 0) then
          write (logfile, *)
          write (logfile, 200) T_instantaneous, T_previous
          T_previous = T_instantaneous
          return
        end if

! Completely quench the velocities if necessary.  Quench the velocities on
! every n`th step if iquench = +n or whenever the instantaneous temperature
! (T_instantaneous) is lower than the instantaneous temperature on the previous ! step (T_previous) if iquench = -1.
        kquench = 0
        if (iquench .gt. 0 .and. mod(itime_step,iquench) .eq. 0) kquench = 1
        if (iquench .eq. -1 .and. T_instantaneous .lt. T_previous) kquench = 1

        if (kquench .eq. 1) then
          write (logfile, *)
          write (logfile, 201)  T_instantaneous, T_previous
          do iatom = 1, s%natoms
            s%atom(iatom)%xdot(1,:) = 0.0d0
            s%atom(iatom)%vatom = 0.0d0
          end do
          T_instantaneous = 0.0d0  ! quenched temperature
        else
          write (logfile, *)
          write (logfile, 202) T_instantaneous, T_previous
        end if

! To really "anneal" a cell we need to do something like constant temperature
! molecular dynamics.  The simplest way to proceed is to rescale velocities at
! each time step get some desired temperature.  If you set iquench = -2, then
! this is accomplished by the  following snippet of code. The input parameters
! are the anneal temperature ("T_want") and taurelax, which specifies how
! rapidly the cell is forced to T_want.  The variable taurelax is essentially
! the relaxation time to make the cell recover T_want. See the code!
        if (iquench .eq. -2) then
          write (logfile, *)
          write (logfile,*) ' Annealing: quench rate = ', taurelax
          write (logfile,*) ' Annealing temperature to reach = ', T_want
          rfactor = sqrt((1.0d0 + dt/taurelax*(T_want/T_instantaneous - 1.0d0)))
          do iatom = 1, s%natoms
            s%atom(iatom)%xdot(1,:) = s%atom(iatom)%xdot(1,:)*rfactor
            s%atom(iatom)%vatom = s%atom(iatom)%xdot(1,:)
          end do
        end if

! Coordinate power quench
! Quench velocities one coordinate at a time.
        if (iquench .eq. -3) then
          do iatom = 1, s%natoms
            do ix = 1, 3
              dot = s%forces(iatom)%ftot(ix)*s%atom(iatom)%vatom(ix)
              if (dot .lt. 0.0d0) then
                s%atom(iatom)%vatom(ix) = 0.0d0
                s%atom(iatom)%xdot(1,ix) = 0.0d0
              end if
            end do
          end do
        end if

! Format Statements
! ===========================================================================
100     format (2x, '         center of mass position = ', 3d12.4)
101     format (2x, '         center of mass velocity = ', 3d12.4)
102     format (2x, ' center of mass angular momentum = ', 3d12.4)
200     format (2x, ' Free Dynamics: T_instantaneous =  ', f12.4,      &
      &         ' T_previous = ', f12.4)
201     format (2x, ' Quenching !! T_instantaneous =  ', f12.4,      &
      &         ' T_previous = ', f12.4)
202     format (2x, ' No quenching. T_instantaneous =  ', f12.4,      &
      &         ' T_previous = ', f12.4)

! End Subroutine
! ===========================================================================
        return
        end subroutine predictor


! End Module
! ===========================================================================
        end module M_dynamics
