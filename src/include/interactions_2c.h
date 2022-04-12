! Definition of two-center interaction types
        integer, parameter :: P_maxtype = 25       ! current max interactions
        integer, parameter :: P_maxsubtype = 9     ! current max subtypes

        integer, parameter :: P_overlap = 1        ! overlap
        integer, parameter :: P_kinetic = 2        ! kinetic
        integer, parameter :: P_vna_ontopL = 3     ! vna ontop - left
        integer, parameter :: P_vna_ontopR = 4     ! vna ontop - right
        integer, parameter :: P_vna_atom = 5       ! vna atom
        integer, parameter :: P_vnl = 6            ! non-local pp
        integer, parameter :: P_vxc_ontop = 7      ! xc ontop
        integer, parameter :: P_vxc_atom = 8       ! xc atom
        integer, parameter :: P_xc_correction = 9  ! xc correction
        integer, parameter :: P_dipole_z = 10      ! z-dipole
        integer, parameter :: P_dipole_y = 11      ! y-dipole
        integer, parameter :: P_dipole_x = 12      ! x-dipole
        integer, parameter :: P_coulomb = 13       ! coulomb
        integer, parameter :: P_eh = 14            ! extended-hubbard
        integer, parameter :: P_rho_ontopL = 15    ! density_ontopl
        integer, parameter :: P_rho_ontopR = 16    ! density_ontopr
        integer, parameter :: P_rho_atom = 17      ! density_atom
        integer, parameter :: P_overlapS = 18      ! spherical overlap
        integer, parameter :: P_rhoS_ontopL = 19   ! spherical density
        integer, parameter :: P_rhoS_ontopR = 20   ! spherical density
        integer, parameter :: P_rhoS_atom = 21     ! spherical density
        integer, parameter :: P_dnuxc_ontopL = 22  ! dnuxc ontopl
        integer, parameter :: P_dnuxc_ontopR = 23  ! dnuxc ontopr
        integer, parameter :: P_goverlapL = 24     ! overlap gradient left
        integer, parameter :: P_goverlapR = 25     ! overlap gradient right

