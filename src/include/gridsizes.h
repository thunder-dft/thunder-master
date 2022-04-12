! ***************************************************************************
! T W O - C E N T E R    P A R A M E T E R S
! ***************************************************************************
! Must be accurate and odd numbers as they determine othogonality

! OVERLAP
        integer, parameter :: ndd_overlap = 107
        integer, parameter :: nz_overlap = 106
        integer, parameter :: nrho_overlap = 106

! KINETIC
        integer, parameter :: nqke = 400
        integer, parameter :: nrke = 240
        integer, parameter :: ndd_ke = 107

        real, parameter :: ecutke = 40.0d3

! HARTREE (VNA)
        integer, parameter :: ndd_vna = 107
        integer, parameter :: nz_vna = 106
        integer, parameter :: nrho_vna = 106

! COULOMB
        integer, parameter :: ndd_coulomb = 107
        integer, parameter :: nz_coulomb = 106
        integer, parameter :: nrho_coulomb = 106

! RHO (Density for McWeda)
        integer, parameter :: ndd_rho = 107
        integer, parameter :: nz_rho = 106
        integer, parameter :: nrho_rho = 106

! RHO_STORE (Density for LDA and GGA - Vxc)
        integer, parameter :: nz_rho_store = 67
        integer, parameter :: nrho_rho_store = 67

! DIPOLE_Z
        integer, parameter :: ndd_dipole_z = 107
        integer, parameter :: nz_dipole_z = 106
        integer, parameter :: nrho_dipole_z = 106

! VNL
        integer, parameter :: ndd_vnl = 107
        integer, parameter :: nz_vnl = 106
        integer, parameter :: nrho_vnl = 106

! VXC
        integer, parameter :: ndd_vxc = 127
        integer, parameter :: nz_vxc = 116
        integer, parameter :: nrho_vxc = 116


! ***************************************************************************
! T H R E E - C E N T E R    P A R A M E T E R S
! ***************************************************************************
! Must be accurate and odd numbers as they determine othogonality
        integer, parameter :: P_ntheta = 5

! RHO (Density for McWeda)
        integer, parameter :: nbc_rho = 29
        integer, parameter :: nna_rho = 29

        ! number of mesh integration points
        integer, parameter :: nnr_rho = 63
        integer, parameter :: nntheta_rho = 63
        integer, parameter :: nnphi_rho = 31

! BCNA
        integer, parameter :: nbc_bcna = 29
        integer, parameter :: nna_bcna = 29

        ! number of mesh integration points
        integer, parameter :: nnr_bcna = 63
        integer, parameter :: nntheta_bcna = 63
        integer, parameter :: nnphi_bcna = 31

! XC3C
        integer, parameter :: nbc_xc = 29
        integer, parameter :: nna_xc = 29

        ! number of mesh integration points
        integer, parameter :: nnr_xc = 63
        integer, parameter :: nntheta_xc = 63
        integer, parameter :: nnphi_xc = 31
