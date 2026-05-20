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
        module M_dynamics_integrator
        use M_configuraciones
        use M_dynamics_settings

! Type Declaration
! ===========================================================================
! None

! Parameter Declaration and Description
! ===========================================================================
        

! Variable Declaration and Description
! ===========================================================================
! None

! Allocate Arrays
! ===========================================================================
! None


! module procedures
        contains

! ===========================================================================
! initialize_dynamics
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
        subroutine initialize_dynamics (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                    !< counter of atoms

! Allocate Arrays
! ===========================================================================
        s%md%ngear = 5
        do iatom = 1, s%natoms
          allocate (s%atom(iatom)%xdot(0:s%md%ngear,3)); s%atom(iatom)%xdot = 0.0d0
        end do

! Procedure
! ===========================================================================
! Gear coeficients for second-order equation (thus cfactor(2) = 1.0d0 always)
        s%md%cfactor = 0.0d0
        if (s%md%ngear .eq. 2) then   ! velocity verlet style
          s%md%cfactor(0) = 0.0d0
          s%md%cfactor(1) = 1.0d0
          s%md%cfactor(2) = 1.0d0
        else if (s%md%ngear .eq. 3) then
          s%md%cfactor(0) = 1.0d0/6.0d0
          s%md%cfactor(1) = 5.0d0/6.0d0
          s%md%cfactor(2) = 1.0d0
          s%md%cfactor(3) = 1.0d0/3.0d0
        else if (s%md%ngear .eq. 4) then
          s%md%cfactor(0) = 19.0d0/120.0d0
          s%md%cfactor(1) = 3.0d0/4.0d0
          s%md%cfactor(2) = 1.0d0
          s%md%cfactor(3) = 1.0d0/2.0d0
          s%md%cfactor(4) = 1.0d0/12.0d0
        else if (s%md%ngear .eq. 5) then
          s%md%cfactor(0) = 3.0d0/20.0d0
          s%md%cfactor(1) = 251.0d0/360.0d0
          s%md%cfactor(2) = 1.0d0
          s%md%cfactor(3) = 11.0d0/18.0d0
          s%md%cfactor(4) = 1.0d0/6.0d0
          s%md%cfactor(5) = 1.0d0/60.0d0
        else if (s%md%ngear .eq. 6) then
          s%md%cfactor(0) = 863.0d0/6048.0d0
          s%md%cfactor(1) = 665.0d0/1008.0d0
          s%md%cfactor(2) = 1.0d0
          s%md%cfactor(3) = 25.0d0/36.0d0
          s%md%cfactor(4) = 35.0d0/144.0d0
          s%md%cfactor(5) = 1.0d0/24.0d0
          s%md%cfactor(6) = 1.0d0/360.0d0
        else if (s%md%ngear .eq. 7) then
          s%md%cfactor(0) = 1925.0d0/14112.0d0
          s%md%cfactor(1) = 19087.0d0/30240.0d0
          s%md%cfactor(2) = 1.0d0
          s%md%cfactor(3) = 137.0d0/180.0d0
          s%md%cfactor(4) = 5.0d0/16.0d0
          s%md%cfactor(5) = 17.0d0/240.0d0
          s%md%cfactor(6) = 1.0d0/120.0d0
          s%md%cfactor(7) = 1.0d0/2520.0d0

        end if

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_dynamics


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
        
! Corrector
        write (logfile,*)
        write (logfile,*) ' Predictor-Corrector: correct the positions. '

        call corrector (s, itime_step)
        do iatom = 1, s%natoms
          s%atom(iatom)%ratom = s%atom(iatom)%xdot(0,:)
          s%atom(iatom)%vatom = s%atom(iatom)%xdot(1,:)
        end do

! Calculate the kinetic energy and the instantaneous temperature
        call writeout_kin_T (s)

! Rescale the temperature if iensemble = 1 (constant temperature MD)
        if (iensemble .eq. 1 .and. .not. s%md%T_instantaneous .le. 0) then
          vscale = sqrt(T_want/s%md%T_instantaneous)
          write (logfile,*) ' Constant temperature ensemble (iensemble = 1): scaling velocities. '
          write (logfile,*) ' Scaling = ', vscale, ' It should be near 1.0!'
          do iatom = 1, s%natoms
            s%atom(iatom)%vatom = s%atom(iatom)%vatom*vscale
             s%atom(iatom)%xdot(1,:) = s%atom(iatom)%xdot(1,:)*vscale
          end do
        end if
        s%md%T_average = ((itime_step - 1)*s%md%T_average + s%md%T_instantaneous)/itime_step
        write (logfile,*) ' Average Temperature = ', s%md%T_average

! Predictor
        write (logfile,*)
        write (logfile,*) ' Predictor-Corrector: predict the positions. '
        call predictor (s, itime_step)

! Scale velocities again after predictor step
        if (iensemble .eq. 1) then
          call writeout_kin_T (s)
          vscale = sqrt(T_want/s%md%T_instantaneous)
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

! Write out coodinates and velocity
        call writeout_coodinate_velocity (s)

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
        integer logfile

        real dtfactor
        real xmass
        real xmass_total

        real, dimension (3) :: xlcm

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
! Initialize logfile
        logfile = s%logfile

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
        end if
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          xmass = species(in1)%xmass
          acceleration = P_fovermp*s%forces(iatom)%ftot/xmass
          difference(:,iatom) = acceleration - s%atom(iatom)%xdot(2,:)
        end do

! Gear (often fifth-order)
        do iorder = 0, s%md%ngear
          dtfactor = dt**(2-iorder)
          do iatom = 1, s%natoms
            s%atom(iatom)%xdot(iorder,:) = s%atom(iatom)%xdot(iorder,:)        &
     &       + s%md%cfactor(iorder)*difference(:,iatom)*(factorial(iorder)/2.0d0)*dtfactor
          end do
        end do

! Calculate the center of mass position and velocity and angular momentum.
        call writeout_momentum (s)

! Deallocate Arrays
! ===========================================================================
        deallocate (difference)

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

        real, dimension (3) :: xlcm

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

! Move the atoms
! For the gear algorithm update both the positions and the velocities after
! each prediction and correction
! Actual prediction step - usually 5th order Gear algorithm.
        do iorder = 0, s%md%ngear - 1
          do isum = iorder + 1, s%md%ngear
            ifactor = isum - iorder
            do iatom = 1, s%natoms
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

! Calculate the center of mass position, velocity and angular momentum.
        call writeout_momentum (s)

! Completely quench the velocities if necessary.  Quench the velocities on
! every n`th step if iquench = +n or whenever the instantaneous temperature
! (T_instantaneous) is lower than the instantaneous temperature on the previous
! step (T_previous) if iquench = -1.
        if (iquench .eq. 0) then
          write (logfile, *)
          write (logfile, 200) s%md%T_instantaneous, s%md%T_previous
          s%md%T_previous = s%md%T_instantaneous
          return
        end if

        kquench = 0
        if (iquench .gt. 0 .and. mod(itime_step,iquench) .eq. 0) kquench = 1
        if (iquench .eq. -1 .and. s%md%T_instantaneous .lt. s%md%T_previous) kquench = 1

        if (kquench .eq. 1) then
          write (logfile, *)
          write (logfile, 201)  s%md%T_instantaneous, s%md%T_previous
          do iatom = 1, s%natoms
            s%atom(iatom)%xdot(1,:) = 0.0d0
            s%atom(iatom)%vatom = 0.0d0
          end do
          s%md%T_instantaneous = 0.0d0  ! quenched temperature
        else
          write (logfile, *)
          write (logfile, 202) s%md%T_instantaneous, s%md%T_previous
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
          rfactor = sqrt((1.0d0 + dt/taurelax*(T_want/s%md%T_instantaneous - 1.0d0)))
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
        s%md%T_previous = s%md%T_instantaneous

! Format Statements
! ===========================================================================
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

! destroy_dynamics
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
        subroutine destroy_dynamics (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                    !< counter of atoms

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          deallocate (s%atom(iatom)%xdot)
        end do

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine destroy_dynamics

! End Module
! ===========================================================================
        end module M_dynamics
