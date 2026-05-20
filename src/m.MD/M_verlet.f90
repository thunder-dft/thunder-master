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
! None

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! initialize_dynamics
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
! None

! Allocate Arrays
! ===========================================================================
        integer iatom                    !< counter of atoms


! Procedure
! ===========================================================================
        if (iensemble .eq. 1) then
          do iatom = 1, s%natoms
            allocate (s%atom(iatom)%langevin)
          end do
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

        integer, intent(in) :: itime_step

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer logfile                      !< writing to which unit

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
          else if ( iensemble .eq. 2 ) then  
            call vverlet_velocity (s)
            call nose_hoover (s)
          end if
        end if

        call writeout_kin_T (s)
        call writeout_momentum (s)
        call writeout_coodinate_velocity (s)

        write (s%logfile,*) ' First half velocity update '
        write (s%logfile,*) ' And coordinate update '   
        if ( iensemble .eq. 0 ) then
          write (s%logfile, *) ' NVE NOW '
          call vverlet_velocity (s)
          call vverlet_coordinate (s)              
        else if ( iensemble .eq. 1 ) then
          write (s%logfile, *) ' Langevin NOW '
          call langevin (s)
          call langevin_velocity (s)
          call langevin_coordinate (s)       
        else if ( iensemble .eq. 2 ) then
          write (s%logfile, *) ' Nose-Hoover NOW '
          if (itime_step .eq. 1) then
            call initialize_nose_hoover (s)
          end if
          call nose_hoover (s)
          call vverlet_velocity (s)
          call vverlet_coordinate (s)   
        end if    
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

          ! Cut off length
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

          ! Cut off length
          nullify (pmd); pmd => s%atom(iatom)%langevin
          s%atom(iatom)%vatom = s%atom(iatom)%vatom  &
&                              + pmd%c1*P_fovermp*s%forces(iatom)%ftot/xmass &
&                              - pmd%c2*s%atom(iatom)%vatom                  &
&                              + pmd%vatom_random
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return     
        end subroutine langevin_velocity

! ===========================================================================
! initialize_nose_hoover
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the subroutine to initialize nose-hoover
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
        subroutine initialize_nose_hoover (s)

        implicit none
  
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! Boltzmann's constant
        real, parameter :: kb = 8.617343693082104d-5

! Variable Declaration and Description
! ===========================================================================
        integer inos
        real tau, pi, omega
! Gaussian
        real num1, num2
        real x1, x2, w, rnum1, rnum2
        real sv
        logical gaus_stored

! Procedure
! ===========================================================================
        ! total chain number is four
        s%md%nnos = 4
        s%md%kT = kb*T_want        

        allocate (s%md%xi (s%md%nnos))
        allocate (s%md%v_xi (s%md%nnos))
        allocate (s%md%Q_i (s%md%nnos))

        ! default 500 fs relaxation time
        tau = 500 
        pi = 4.0 * atan(1.0)        
        omega = 2.0*pi/tau 
        s%md%xi = 0.0
        do inos = 1, s%md%nnos
          !Q is eV*fs^2=M*L^2
          if (inos .eq. 1) s%md%Q_i(inos) = s%natoms*3.0*s%md%kT/omega**2
          s%md%Q_i(inos) = s%md%kT/omega**2
        end do

        ! initialize the thermostat velocities with a Gaussian distribution
        gaus_stored = .false.

        do inos = 1, s%md%nnos

           ! ---- Gaussian generator inline ----
           if (gaus_stored) then
              num1 = num2
              gaus_stored = .false.
           else
              do
                 call random_number(rnum1)
                 call random_number(rnum2)

                 x1 = 2.0*rnum1 - 1.0
                 x2 = 2.0*rnum2 - 1.0
                 w  = x1*x1 + x2*x2

                 if (w > 0.0 .and. w < 1.0) exit
              end do

              w = sqrt((-2.0*log(w))/w)

              num1 = x1*w
              num2 = x2*w
              gaus_stored = .true.
           end if
           ! ---- End Gaussian ----
           sv = sqrt(s%md%kT/s%md%Q_i(inos))
           s%md%v_xi(inos) = sv*num1

        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        end subroutine initialize_nose_hoover


! ===========================================================================
! nose_hoover
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the subroutine of nose-hoover
! This routine does the nose-hoover part of the integration from t=0 to t=dt/2
! get the total kinetic energy
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
        subroutine nose_hoover(s)

        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: nresn = 1 ! Orders of yoshida integrators
        integer, parameter :: nyosh = 5 ! Orders of yoshida integrators
        
! Variable Declaration and Description
! ===========================================================================
! it is taken from (G.J. Martyna et.al. Mol Phys 1996 87(5) 1117-1157)

        integer iatom	                                ! counter of atoms
        integer iresn, iyosh, inos

        real akin, scale, G_i, aa
        real www                                        ! Coeffients of yoshida integrators
        real, dimension(5) :: wdti2, wdti4, wdti8       ! Coeffients of yoshida integrators

! Procedure
! ===========================================================================
        if (nyosh .eq. 1) then

         wdti2(1) = dt / ( 2.d0 * nresn )
         wdti4(1) = wdti2(1) * 0.5d0
         wdti8(1) = wdti4(1) * 0.5d0

        else if (nyosh .eq. 3) then

         www = 1.d0 / ( 2.d0 - 2.d0**(1.d0/3.d0) )

         wdti2(1) = www * dt / ( 2.d0 * nresn )
         wdti2(2) = ( 1.d0 - 2.d0 * www ) * dt / ( 2.d0 * nresn )
         wdti2(3) = wdti2(1)

         wdti4(1) = wdti2(1) * 0.5d0
         wdti4(2) = wdti2(2) * 0.5d0
         wdti4(3) = wdti4(1)

         wdti8(1) = wdti4(1) * 0.5d0
         wdti8(2) = wdti4(2) * 0.5d0
         wdti8(3) = wdti8(1)

        else if (nyosh .eq. 5) then

         www = 1.d0 / ( 4.d0 - 4.d0**(1.d0/3.d0) )

        ! the extra factor of 2 is because we are doing dt/2.0         
         wdti2(1) = www * dt / ( 2.0*nresn )
         wdti2(2) = wdti2(1)
         wdti2(3) = ( 1.d0 - 4.d0 * www ) * dt / ( 2.0*nresn )
         wdti2(4) = wdti2(1)
         wdti2(5) = wdti2(1)

         wdti4(1) = wdti2(1) * 0.5d0
         wdti4(2) = wdti4(1)
         wdti4(3) = wdti2(3) * 0.5d0
         wdti4(4) = wdti4(1)
         wdti4(5) = wdti4(1)

         wdti8(1) = wdti4(1) * 0.5d0
         wdti8(2) = wdti8(1)
         wdti8(3) = wdti4(3) * 0.5d0
         wdti8(4) = wdti8(1)
         wdti8(5) = wdti8(1)

        end if

        scale = 1.0
        akin = s%md%tkinetic

        ! start the multiple time step procedure
        do iresn = 1, nresn
          do iyosh = 1, nyosh
 
        !   update the thermostat velocities
            do inos = s%md%nnos, 1, -1
              if (inos .lt. s%md%nnos) then
                aa = exp(-wdti8(iyosh)*s%md%v_xi(inos+1))
                s%md%v_xi(inos) = s%md%v_xi(inos)*aa
              end if
  
              if (inos .eq. 1) then
                G_i = (2*akin - s%natoms*3.0*s%md%kT)/s%md%Q_i(inos)
              else if (inos .gt. 1) then
                G_i = (s%md%Q_i(inos-1)*s%md%v_xi(inos-1)**2 &
                  - s%md%kT)/s%md%Q_i(inos)
              end if
  
              s%md%v_xi(inos) = s%md%v_xi(inos) + wdti4(iyosh)*G_i
  
              if (inos .lt. s%md%nnos) s%md%v_xi(inos) = s%md%v_xi(inos)*aa
            end do  
 
        !   update the particle velocities
            aa = exp(-wdti2(iyosh)*s%md%v_xi(1))
            scale = scale*aa
            akin = akin*aa*aa
 
        !   update the thermostat positions
            do inos = 1, s%md%nnos
             s%md%xi(inos) = s%md%xi(inos) + s%md%v_xi(inos)*wdti2(iyosh)
            end do
 
        !  update the thermostat velocities
            do inos = 1, s%md%nnos
              if (inos .lt. s%md%nnos) then
                aa = exp(-wdti8(iyosh)*s%md%v_xi(inos+1))
                s%md%v_xi(inos) = s%md%v_xi(inos)*aa
              end if
  
              if (inos .eq. 1) then
                G_i = (2*akin - s%natoms*3.0*s%md%kT)/s%md%Q_i(1)
              else if (inos .gt. 1) then
                G_i = (s%md%Q_i(inos-1)*s%md%v_xi(inos-1)**2 &
                  - s%md%kT)/s%md%Q_i(inos)
              end if
  
              s%md%v_xi(inos) = s%md%v_xi(inos) + wdti4(iyosh)*G_i
  
              if (inos .lt. s%md%nnos) s%md%v_xi(inos) = s%md%v_xi(inos)*aa
            end do
 
          end do
        end do

        !update the particle velocities
        do iatom = 1, s%natoms
          s%atom(iatom)%vatom = s%atom(iatom)%vatom*scale
        end do

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        end subroutine nose_hoover

! ===========================================================================
! update nose-hoover
! ===========================================================================
! Subroutine Description
! ===========================================================================
! This is the subroutine to update nose-hoover
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
        subroutine update_nose_hoover (s)

        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input  
        type(T_structure), target :: s           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! Boltzmann's constant
        real, parameter :: kb = 8.617343693082104d-5

! Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
          s%md%kT = kb*T_want
         s%md%Q_i = s%md%Q_i*T_want/T_want
        s%md%v_xi = s%md%v_xi*sqrt(T_want/T_want)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        end subroutine update_nose_hoover


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
        ! sqrt(1 + (T0/T - 1) * dt/tau)
        T_scale = sqrt(1.0d0 + (T_want/s%md%T_instantaneous - 1.0d0) &
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
        if (iensemble .eq. 1) then
          do iatom = 1, s%natoms
            deallocate (s%atom(iatom)%langevin)
          end do
        else if (iensemble .eq. 2) then
          deallocate (s%md%xi)
          deallocate (s%md%v_xi)
          deallocate (s%md%Q_i)
        end if

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
