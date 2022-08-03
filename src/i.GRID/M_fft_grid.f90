! copyright info:
!
!                             @Copyright 2014
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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

! M_fft_grid.f90
! Module Description
! ===========================================================================
!       This is a module containing all subroutines related to solving the
! Laplace equation for the density functional theory grid code. The following
! subroutines are called here:
!
!       Laplace_fft.f90 - solve the Laplace equation via Fourier Transform
!
! ===========================================================================
! Code written by:
! Prokop Hapala
! Pavel Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
!
! James P. Lewis
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
        module M_fft_grid
        use M_species
        use M_configuraciones
        use M_grid

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer, parameter :: fftw_estimate = 64

! module procedures
        contains


! ===========================================================================
! Laplace_fft
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calculates the Hartree Potential via the Poisson equation
! in reciprocal space. The subroutine utilizes the FFTW library.
!
! Given:
! d^2 Vh(r) = -4*Pi*n(r)   ; n(r)=sum_k(n(k)*exp(ikr))
! d^2 Vh(k)*exp(ikr) = i^2*k^2*Vh(k)*exp(ikr)= -1*k^2*Vh(k)*exp(ikr)
!
! final formulae for k-conponent of Vh:  Vh(k) = 4*Pi*n(k)/k^2
! ===========================================================================
! Code written by:
! Prokop Hapala
! Pavel Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
!
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ============================================================================
! Subroutine Declaration
! ===========================================================================
        subroutine Laplace_fft (t)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: t           !< the structure to be used

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer i, j, k
        integer index1                      !< counting
        integer kx, ky, kz
        integer logfile                     !< writing to which unit

        real dmax
        real factor

        real, dimension (3) :: kvector

! variables used by FFTW
        integer*8 :: plan

        double precision, allocatable, dimension(:,:,:) :: Vin
        double complex, allocatable, dimension(:,:,:) :: Gout
        real, target, allocatable, dimension (:) :: result

! For the visualization part
!        integer iatom

        character (len=25) xsfname
        character (len=25) message

        real, dimension (:), pointer   :: pmat
        interface
          subroutine writeout_xsf (t, xsf, xsfname, message)
            use M_configuraciones
            use M_species
            implicit none
            type(T_structure), target :: t           !< the structure to be used
            real, dimension (:), pointer :: xsf
            character (len=25), intent(in) :: xsfname
            character (len=25), intent(in) :: message
! Inconsistency Interface declaration
! Following implemented subroutine
!           real, pointer :: xsf (:)
          end subroutine writeout_xsf
        end interface

! Allocate Arrays
! ===========================================================================
! FFTW can save half of the space for the real data.
! Hence, this is why allocation of four3D uses just irm1 / 2.
        allocate (Vin (irm1, irm2, irm3))
        allocate (Gout (irm1/2 + 1, irm2, irm3))
        allocate (result (0:nrm - 1))

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = t%logfile
        write (logfile,*)
        write (logfile,*) ' Assemble Kohn-Sham potentials - Vxc(G) and Vna(G) '

! Save the previous result
        result = vcaG

! Get data from 1D input to 3D: drhoG, rm are from the 'grid' module
        do k = 1, irm3
          do j = 1, irm2
            do i = 1, irm1
              index1 = i + irm1*(j - 1) + irm1*irm2*(k - 1) - 1
! (to remind the current units are 1/Ang^3)
              Vin(i,j,k) = drhoG(index1)
            end do
          end do
        end do

! Forward FFT  n(r) -> n(k)
! If needed, there are more precise FFTW modes than FFTW_ESTIMATE. But so far,
! the 'ESTIMATE' has proved to be good enough for all purposes (even better).
        call dfftw_plan_dft_r2c_3d (plan, irm1, irm2, irm3, Vin, Gout, fftw_estimate)
        call dfftw_execute (plan)
        call dfftw_destroy_plan (plan)

! Comparing coefficients in Fourier series of source term and solution
! FFT saves 'positive' frequencies in the first half of the arrays,
! the 'negative' ones in the second half.
        do i = 1, (irm1 / 2 + 1)
          if (i <= (irm1 / 2 + 1)) then
            kx = i - 1
          else
            kx = irm1 - i + 1
          end if

          do j = 1, irm2
            if (j <= (irm2 / 2 + 1)) then
              ky = j - 1
            else
              ky = irm2 - j + 1
            end if

            do k = 1, irm3
              if (k <= (irm3 / 2 + 1)) then
                kz = k-1
              else
                kz = irm3 - k + 1
              end if

              ! variable transform
              kvector(1) = s%g(1)%a(1)*kx + s%g(1)%a(2)*ky + s%g(1)%a(3)*kz
              kvector(2) = s%g(2)%a(1)*kx + s%g(2)%a(2)*ky + s%g(2)%a(3)*kz
              kvector(3) = s%g(3)%a(1)*kx + s%g(3)%a(2)*ky + s%g(3)%a(3)*kz

! we divide both part of the complex number by the real number
              factor = 4.0d0*pi/(kvector(1)*kvector(1)                       &
     &                           + kvector(2)*kvector(2)                     &
     &                           + kvector(3)*kvector(3) + 1.0d-08)
              Gout(i,j,k) = factor*Gout(i,j,k)
            end do
          end do
        end do
        Gout(1,1,1) = 0.0d0 ! zero Fourier term

! Backward FFT Vh(k) -> Vh(r)
        call dfftw_plan_dft_c2r_3d (plan, irm1, irm2, irm3, Gout, Vin, fftw_estimate)
        call dfftw_execute (plan)
        call dfftw_destroy_plan (plan)

! Convert units of the potential from a.u. to eV and divide by nrm as for
! FFT renormalization
        factor = P_eq2/nrm

! Transform back to 1D, returning potential
        do k = 1, irm3
          do j = 1, irm2
            do i = 1, irm1
              index1 = i + (j - 1)*irm1 + (k - 1)*irm1*irm2 - 1
! renormalize
             vcaG(index1) = Vin(i,j,k)*factor
            end do
          end do
        end do

! evaluate residual of vcaG
        dmax = 0.0d0
        do index1 = 0, nrm - 1
          dmax = max(dmax, abs(vcaG(index1) - result(index1)))
        end do
        write (logfile,*) ' residual vcaG = ', dmax, dmax*dvolume

! Visualization in XCrysDen
! ===========================================================================
! Write out hartree potential piece resulting from drho
        pmat => vcaG
        xsfname = 'dVHartree.xsf'
        message = '3D-plot of Hartree Potential'
        call writeout_xsf (t, pmat, xsfname, message)

! write out Hartree potential into xsf file
        result = vcaG + vnaG
        pmat => result
        xsfname = 'VHartree.xsf'
        message = '3D-plot of Hartree Potential'
        call writeout_xsf (t, pmat, xsfname, message)

! Deallocate Arrays
! ===========================================================================
        deallocate (Vin)
        deallocate (Gout)
        deallocate (result)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine Laplace_fft

! End Module
! ===========================================================================
        end module M_fft_grid
