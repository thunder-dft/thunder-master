! Definition of three-center interaction types.
        integer, parameter :: P_maxtype = 4     !< current max interactions
        integer, parameter :: P_maxsubtype = 4  !< current max subtypes
        integer, parameter :: P_maxtheta = 5    !< current max angles

        integer, parameter :: P_rho_3c = 1      ! density for McWEDA
        integer, parameter :: P_rhoS_3c = 2     ! spherical density for McWEDA
        integer, parameter :: P_bcna = 3        ! Neutral Atom (bcna)
        integer, parameter :: P_xc3c = 4        ! Horsfield (xc3c)
