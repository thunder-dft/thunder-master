! Various constants
       real, parameter :: P_abohr = 0.529177d0   !< Bohr radii to Angstrom
       real, parameter :: P_Debye = 0.208194d0   !< dipole units e-Angstrom
       real, parameter :: P_eq2 = 14.39975d0     !< conversion to get eV
       real, parameter :: P_Hartree = P_eq2/P_abohr  !< Hartree energy
       real, parameter :: P_kconvert = 11604.49558d0 !< conversion to kT
       real, parameter :: P_fovermp = 0.009648957597d0

! For now - spin independent - so set parameter
       real, parameter :: P_spin = 2.0d0

! Kronecker delta
       real, parameter, dimension (3, 3) ::    delk =                        &
     &   reshape ([1.0d0, 0.0d0, 0.0d0,                                      &
     &             0.0d0, 1.0d0, 0.0d0,                                      &
     &             0.0d0, 0.0d0, 1.0d0], [3,3])

! Levi-civita parameter
       real, parameter, dimension (3, 3, 3) :: xlevi =                       &
     &   reshape ([0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, -1.0d0, 0.0d0,  1.0d0,&
     &             0.0d0, 0.0d0, 0.0d0,  1.0d0, 0.0d0,  0.0d0, 0.0d0, -1.0d0,&
     &             0.0d0, 0.0d0, 0.0d0, -1.0d0, 0.0d0,  1.0d0, 0.0d0,  0.0d0,&
     &             0.0d0, 0.0d0, 0.0d0], [3,3,3])

! The A matrix - for rotations
       real, parameter, dimension (3, 3, 5) :: amat =                        &
     &   reshape ([0.0d0, 0.5d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,   &
     &             0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0,   &
     &             0.5d0, 0.0d0, -1.0d0/(2.0d0*sqrt(3.0d0)), 0.0d0, 0.0d0,   &
     &             0.0d0, -1.0d0/(2.0d0*sqrt(3.0d0)), 0.0d0, 0.0d0, 0.0d0,   &
     &             2.0d0/(2.0d0*sqrt(3.0d0)), 0.0d0, 0.0d0, 0.5d0, 0.0d0,    &
     &             0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0,   &
     &             0.0d0, -0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], [3,3,5])

! Tolerances
       real, parameter :: xc_overtol = 1.0d-6  !< tolerance  for the definition of average rho

! Complex numbers
       complex, parameter :: a0 = (0.0d0, 0.0d0)
       complex, parameter :: a1 = (1.0d0, 0.0d0)
       complex, parameter :: ai = (0.0d0, 1.0d0)
