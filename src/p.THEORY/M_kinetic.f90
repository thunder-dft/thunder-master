! copyright info:
!
!                             @Copyright 2016
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
! M_case2.f90
! Program Description
! ============================================================================
!      This is a module calculating the integrals of two centers for the
! kinetic energy.
!
! ============================================================================
! Code written by:
! Barry Haycock
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ============================================================================
! Module Declaration
! ============================================================================
        module M_kinetic
        use M_atom_functions
        use M_species
        use M_integrals_2c

        implicit none

! Type Declaration
! ===========================================================================


! module procedures
        contains

! ===========================================================================
! initialize_kinetic
! ===========================================================================
! Program Description
! ===========================================================================
!       We need to determine how many interactions belong to each nspecies
! bundle pair. This routine just counts how many total interactions contribute
! to that bundle.  Something like overlap is obviously only 1 interaction
! added, but something like vna needs number of interactions based on the
! number of shells.
!
! ===========================================================================
! Code written by:
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
! Subroutine Declaration
! ===========================================================================
        subroutine initialize_kinetic
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species

        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

! Procedure
! ============================================================================
! Loop over species
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
          end do ! jspecies
        end do ! ispecies

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_kinetic


! kinetic.f90
! Program Description
! ===========================================================================
!       Calculates kinetic energy matrix elements. Specifically, this routine
! will create kinetic energy matrix elements for two atomic species designated
! by in1 and in2, which gets passed in the call list. Appropriately named
! output files are created. Essentially, we calculate this kinetic energy in
! k-space for easier manipulation so that there is not need to take
! derivatives.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1430 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine kinetic
        implicit none

        include "../include/gridsizes.h"

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
        integer kmax
        parameter (kmax = 6)           !< maximum number of surviving bessels

! FIX ME! This is a clumsy way to do this....
        integer nsh_max
        parameter (nsh_max = 3)        !< maximum number of shells

! eV= ((h/2pi)**2)/m/pi (in eV times angstroms squared)
        real eV
        parameter (eV = 2.4255d0)      !< conversion factor

! Local Variable Declaration and Description
! ===========================================================================
        integer ispecies, jspecies          !< counters for number of species
        integer isorp
        integer idgrid, iqgrid, irgrid      !< counter for grid points
        integer index_2c, nME2c_max         !< basically the number of non-zero
        integer logfile                     !< writing to which unit
        integer nFdata_cell_2c              !< indexing of interactions

        integer l1, m1, n1                  !< quantum  numbers
        integer l2, m2, n2

        integer imu
        integer issh                        !< counter to loop over shells
        integer ktype
        integer lssh
        integer nqq, nrr                    ! number of grid points

        real ang_integral                   !< angular integral
        real d                              !< distance between centers
        real factor                         !< Simpson's integration factor
        real q, qmax, dq                    !< kinetic energy integration
        real r, rmin, dr, dmax              !< radial integration

        real rmax                           !< maximum value of r along grid

        real, allocatable :: angular (:, :)
        real answer (0:kmax)

! Stuff for radial bessel integration
        type T_Bessel_species
          ! we allocate this according to number of shells and grid points
          real, pointer :: FofR (:, :)      ! value of the function
        end type T_Bessel_species

        type(T_Bessel_species), pointer :: Bessel (:)

        type (T_Fdata_cell_2c), pointer :: pFdata_cell
        type (T_Fdata_bundle_2c), pointer :: pFdata_bundle

        character (len = 30) filename
        character (len = 25) interactions

        logical skip

! FOR TESTING PURPOSES ONLY - REMOVE COMMENTS TO TEST
! The 5 in the esplit is to split q into 5 ranges for comparison purposes.
! This is similarly done for xnormq and qsplit below.
        integer isplit
        integer nsplit
        integer nqtop

        real xnqtop
        real xtra

        real esplit (5)
        real qsplit (5)

        real, allocatable :: sumq (:)
        real, allocatable :: xnormq (:, :)

        logical testing

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = 21

! Initialize some variables
        isorp = 0
        qmax = sqrt(2.0d0*ecutke/7.62d0)

        ! Set up integration parameters
        nqq = 2*nqke + 1
        nrr = 2*nrke + 1
        rmin = 0.0d0
        dq = qmax/real(nqq - 1)

! ***************************************************************************
!
! F O U R I E R   T R A N S F O R M   O F   R A D I A L   C O M P O N E N T
! ***************************************************************************
        allocate (Bessel(nspecies))
        do ispecies = 1, nspecies
          q = - dq
          allocate (Bessel(ispecies)%FofR(species(ispecies)%nssh, nqq))
          Bessel(ispecies)%FofR = 0.0d0

! Use Simpson's rule
! Now do integral over r for this fixed q value. We are integrating bessel
! functions here. The r integral goes from 0 to rcutoff. Do this for both
! atoms.
          do iqgrid = 1, nqq
            q = q + dq
            dr = (wf(ispecies)%rcutoffA_max - 0.0d0)/real(nrr - 1)
            r = -dr
            do irgrid = 1, nrr
              r = r + dr
              factor = 4.0d0/3.0d0
              if (mod(irgrid,2) .eq. 1) factor = 2.0d0/3.0d0
              if (irgrid .eq. 1 .or. irgrid .eq. nrr) factor = 1.0d0/3.0d0
              do issh = 1, species(ispecies)%nssh
                lssh = species(ispecies)%shell(issh)%lssh
                rmax = species(ispecies)%shell(issh)%rcutoffA
                Bessel(ispecies)%FofR(issh, iqgrid) =                        &
     &            Bessel(ispecies)%FofR(issh, iqgrid)                        &
     &            + factor*jl2(lssh,q*r)*psiofr(r, rmax, ispecies, issh)*r*r*dr
              end do ! end loop over shells
            end do ! end loop over irgrid
          end do ! end loop over iqgrid

! FOR TESTING PURPOSES ONLY - REMOVE COMMENTS TO TEST
! Test the normalization - we first normalize in q space. Then, in q-space we
! get a factor (4.0d0*pi)/((2.0d0*pi)**3) = 2.0d0/pi. Afterwards, normalize
! in r-space.
          testing = .false.
          if (testing) then
            allocate (sumq(species(ispecies)%nssh))
            allocate (xnormq(0:species(ispecies)%nssh, 5))
            xtra = 2.0d0/pi
            nsplit = 5
            do isplit = 1, nsplit
              xnqtop = real(nqq)*sqrt(real(isplit)/real(nsplit))
              nqtop = int(xnqtop)
              sumq = 0.0d0
              q = - dq
              do iqgrid = 1, nqtop
                q = q + dq
                factor = (4.0d0/3.0d0)
                if (mod(iqgrid,2) .eq. 1) factor = (2.0d0/3.0d0)
                if (iqgrid .eq. 1 .or. iqgrid .eq. nqq) factor = (1.0d0/3.0d0)
                do issh = 1, species(ispecies)%nssh
                  sumq(issh) = sumq(issh)                                    &
     &              + factor*q*q*dq*Bessel(ispecies)%FofR(issh,iqgrid)**2*xtra
                end do
              end do

              qsplit(isplit) = q
              esplit(isplit) = ecutke*((q/qmax)**2)
              do issh = 1, species(ispecies)%nssh
                xnormq(issh,isplit) = sumq(issh)
              end do
            end do

            write (logfile,*)
            write (logfile,*) ' Normalization in q space for ispecies = ', ispecies
            do issh = 1, species(ispecies)%nssh
              write (logfile,*)
              write (logfile,*) ' Shell: ', issh
              write (logfile,*)
              write (logfile,301) issh
              do isplit = 1, 5
                write (logfile,302) qsplit(isplit), esplit(isplit), xnormq(issh,isplit)
              end do
            end do
            write (logfile,*)
            deallocate (sumq)
            deallocate (xnormq)
          end if
        end do ! end loop over species

! ***************************************************************************
! E N D   F O U R I E R   T R A N S F O R M
!
! ***************************************************************************

! ***************************************************************************
!
! C A L C U L A T E   K I N E T I C   E N E R G Y
! ***************************************************************************
        write (logfile,*)
        write (logfile,*) ' ******************************************************* '
        write (logfile,*) '          K I N E T I C   I N T E R A C T I O N S        '
        write (logfile,*) ' ******************************************************* '

! We are ready to go
        do ispecies = 1, nspecies
          do jspecies = 1, nspecies
            pFdata_bundle=>Fdata_bundle_2c(ispecies, jspecies)
            pFdata_bundle%nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c + 1
            nFdata_cell_2c = pFdata_bundle%nFdata_cell_2c
            pFdata_cell=>pFdata_bundle%Fdata_cell_2c(nFdata_cell_2c)

! Open necessary file
            call make_munu (nFdata_cell_2c, ispecies, jspecies)
            nME2c_max = pFdata_cell%nME
            allocate (pFdata_cell%fofx(nME2c_max))

            ! Open ouput file for this species pair
            write (filename, '("/kinetic.",i2.2,".",i2.2,".dat")')           &
     &             species(ispecies)%nZ, species(jspecies)%nZ
            inquire (file = trim(Fdata_location)//trim(filename), exist = skip)
            if (skip) cycle
            open (unit = 11, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown')

! This section caculates and stores R0(q), R1(q), R2(q), R3(q).
! RL(q)=int {jL(qr) * RL(r) * r**2} dr
! Remember that psi_nlm (r) = Y_lm(Omega_r) R_nl(r), and
!               psi_nlm (q) = 4 pi (-i)^l Y_lm(Omega_q) R_nl(q)

! ***************************************************************************
!
! E S T A B L I S H   A N G U L A R   I N T E G R A T I O N   F A C T O R S
! ***************************************************************************
            allocate (angular (0:kmax, nME2c_max))
            do index_2c = 1, nME2c_max
              n1 =  pFdata_cell%N_mu(index_2c)
              l1 =  pFdata_cell%L_mu(index_2c)
              m1 =  pFdata_cell%M_mu(index_2c)

              n2 =  pFdata_cell%N_nu(index_2c)
              l2 =  pFdata_cell%L_nu(index_2c)
              m2 =  pFdata_cell%M_nu(index_2c)

! The variable angular(k) is the integral of Ylm*Pk[cos(theta)]*Yl'm'.
! We only consider terms up to k = 6. This is for the following reason
! The selection rules are such that |l1-k|<=l2<=l1+k.  Therefore, if we
! are only doing up to f-orbitals here, l1 = 0, 1, 2, 3 and l2 = 0, 1, 2, 3.
! For these values of l1 and l2 k is bounded from above by 6.
! See my notes or check this out for yourself if you don't believe me!
! kmax is dimensioned in parameters above.
              call Pintegral (l1, m1, l2, m2, kmax, answer)

              angular(0,index_2c) = answer(0)
              angular(1:kmax,index_2c) = 0.0d0
              do imu = 1, kmax
                if (mod(l1-l2-imu,2) .eq. 0)                                 &
     &            angular(imu,index_2c) =                                    &
     &              answer(imu)*(2*imu+1)*((-1)**((l1-l2-imu)/2))
              end do
            end do   ! end loop for matrix elements
! ***************************************************************************

! Now calculate the kinetic energy for all distances.
! The following section calculates the kinetic energy matrix elements
! for all the possible values of d
! Some preliminaries. Set up simpson rule factors and eV.
! Do a convergence test for d = 0.0
! ***************************************************************************
            write (logfile,200) species(ispecies)%nZ, species(jspecies)%nZ
            dr = (wf(ispecies)%rcutoffA_max                                  &
     &            + wf(jspecies)%rcutoffA_max)/real(ndd_ke - 1)
            dmax = (wf(ispecies)%rcutoffA_max + wf(jspecies)%rcutoffA_max)
            d = - dr

            ! open directory file
            write (interactions,'("/2c.",i2.2,".",i2.2,".dir")')             &
     &        species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 13, file = trim(Fdata_location)//trim(interactions),&
     &            status = 'unknown', position = 'append')
            write (13,100) pFdata_bundle%nFdata_cell_2c, P_kinetic, isorp,   &
     &                     filename(2:30), pFdata_cell%nME, ndd_ke, dmax
            close (unit = 13)

            ! Open mu, nu, mvalue file and write out values.
            write (filename, '("/",i2.2, "_munu_2c.",i2.2,".",i2.2,".dat")') &
     &             P_kinetic, species(ispecies)%nZ, species(jspecies)%nZ
            open (unit = 12, file = trim(Fdata_location)//trim(filename),    &
     &            status = 'unknown', position = 'append')

            ! write the mapping - stored in mu, nu, and mvalue
            write (12,*) (pFdata_cell%mu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%nu_2c(index_2c), index_2c = 1, nME2c_max)
            write (12,*) (pFdata_cell%mvalue_2c(index_2c),                   &
     &                    index_2c = 1, nME2c_max)

            do idgrid = 1, ndd_ke
              d = d + dr

! Initialize
              pFdata_cell%fofx = 0.0d0

! Now the integral over q loop.
              q = - dq
              do iqgrid = 1, nqq
                q = q + dq

! Simpson rule for integration.
                factor = (4.0d0/3.0d0)*eV
                if (mod(iqgrid,2) .eq. 1) factor = (2.0d0/3.0d0)*eV
                if (iqgrid .eq. 1 .or. iqgrid .eq. nqq) factor = (1.0d0/3.0d0)*eV

                do index_2c = 1, nME2c_max
                  n1 = pFdata_cell%N_mu(index_2c)
                  n2 = pFdata_cell%N_nu(index_2c)

! We only consider terms up to k = 6. This is for the following reason.
! The selection rules are such that |l1-k|<=l2<=l1+k.  Therefore, if we
! are only doing up to f-orbitals here, l1 = 0, 1, 2, 3 and l2 = 0, 1, 2, 3.
! The maximum value of k is 6!
                  ang_integral = 0.0d0
                  do ktype = 0, 6
                    ang_integral =                                           &
     &                ang_integral + angular(ktype,index_2c)*jl2(ktype,q*d)
                  end do
                  pFdata_cell%fofx(index_2c) = pFdata_cell%fofx(index_2c)    &
     &              + factor*Bessel(ispecies)%FofR(n1,iqgrid)                &
     &                      *Bessel(jspecies)%FofR(n2,iqgrid)                &
     &                      *ang_integral*q**4*dq
                end do
              end do

              ! Write out details.
              write (11,*) (pFdata_cell%fofx(index_2c),                      &
     &                                       index_2c = 1, nME2c_max)
            end do ! idgrid
            write (11,*)

            deallocate (angular)
          end do ! ispecies
        end do ! jspecies

! Format Statements
! ===========================================================================
100     format (2x, i3, 1x, i3, 1x, i3, 1x, a29, 1x, i3, 1x, i4, 1x, f9.6)
200     format (2x, ' Evaluating kinetic integrals for nZ = ', i3,           &
     &              ' and nZ = ', i3)
301     format (8x, 'q', 12x, 'E', 5x, 'Int[psi (l = ', i1,',q)**2]')
302     format (3x, f9.3, 1x, f11.2, 1x, f12.7)

        return
        end subroutine kinetic

! ===========================================================================
! Pintegral.f
! ===========================================================================
! Program Description
! ===========================================================================
! This subroutine calculates the integral of two spherical harmonics of
! l,m and lp,mp multiplied by Legendre Polynomials of order from
! 0 to kmax.
!   integral (Ylm*Pn(cos(theta))Ylpmp)
!
! This will handle polynomials up to kmax = 10, in accordance
! with the subroutine integral which calculates the integral of two
! spherical harmonics multiplied by cos(theta) to the kth power.
!
! This will use the recursion relation for the Legendre Polynomials
!  (l+1)Pl+1(x) = (2l+1)*x*Pl(x) - l*Pl-1(x)
!   where x = cos(theta)
!
!  input:
!     l,m = l and m value for first spherical harmonic
!     lp,mp = l and m value for second spherical harmonic
!     kmax = highest order Legendre Polynomial to calculate
!           will calculate polynomials from 0 to kmax
!     Pint(i) = array of values for the integral of spherical harmonics
!          times Leg. Polynomial of order i
!
! ===========================================================================
! Code written by:
! Kirk VanOpdorp
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-5909
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Pintegral (l, m, lp, mp, kmax, Pint)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer l
        integer m
        integer lp
        integer mp
        integer kmax

! Output
        real Pint (0:kmax)

! Local Parameters and Data Declaration
! ===========================================================================
        integer nmax, imax, jmax
        parameter (nmax = 10, imax = 10, jmax=10)

! Local Variable Declaration and Description
! ===========================================================================
        integer i, j

        real arg

        real pn (0:imax,0:jmax)
        real theta (0:nmax)

        complex sum

! Procedure
! ===========================================================================
!----------------------------------------------------------------------
! Get the cosine integrals
!----------------------------------------------------------------------
        do i = 0, kmax
          call Tintegral (l, m, lp, mp, i, sum)
          theta(i) = dble(sum)
        end do

!----------------------------------------------------------------------
! Get the Legendre Polynomials by the recurrence relation
!----------------------------------------------------------------------
        do i = 0, imax
          do j = 0, jmax
            pn(i,j) = 0.0d0
          end do
        end do

        pn(0,0) = 1.0d0
        pn(1,1) = 1.0d0

        do i = 2, kmax
          do j = 0,i
            arg = dble(i)
            pn(i,j+1) = pn(i,j+1) + (2.0d0*arg - 1.0d0)*pn(i-1,j)/arg
            pn(i,j) = pn(i,j) - (arg - 1.0d0)*pn(i-2,j)/arg
          end do
        end do

!---------------------------------------------------------------------
! Calculate the integral of Legendre Polynomial with spherical
! harmonics.
!---------------------------------------------------------------------
        do i = 0, kmax
          Pint(i) = 0.0d0
          do j = 0, i
            Pint(i) = Pint(i) + pn(i,j)*theta(j)
          end do
        end do

! Format Statements
! ===========================================================================
! None

        return
        end subroutine Pintegral

! ===========================================================================
! Tintegral.f
! ===========================================================================
! Program Description
! ===========================================================================
! This subroutine will calculate the Tintegral of a product of two
! different spherical harmonics and cos(theta) to the power of k
!
!          integral ( Ylm*[cos(theta)]**kYl'm' )
!
! This routine uses this formula to expand the cos(theta)*Yl'm'
!
!      cos(theta)*Ylm = c1(l,m)*Yl+1,m + c2(l,m)*Yl-1,m
!
! where c1(l,m) = ((l-m+1)*(l+m+1)/(2*l+1)/(2*l+3))^1/2
!  and  c2(l,m) = ((l-m)(l+m)/(2l-1)/(2l+3))^1/2
!
! This is expanded using a tree algorithm as explained in the subroutine
!  tree.
!
! Once the (cos(theta))^k*Ylm is expanded, orthogonality relationships
!   are used to calculate the value of the integral.
!
! Right now the maximum k that can be done is k = 10
!   The number of array elements needed for any k is: 2**(k+1) - 1
!   so the cl and iyl arrays are dimensioned as cl(0:2046) and iyl(0:2046)
! =========================================================================

!
! ===========================================================================
! Code rewritten by:
! Kirk VanOpdorp
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-5909
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Tintegral (l1, m1, l2, m2, k, sum)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        integer l1
        integer m1
        integer l2
        integer m2
        integer k

! Output:
        complex sum             ! value of the integral

! Local Parameters and Data Declaration
! ===========================================================================
        integer kmax
        integer maxnode
        parameter (kmax = 10, maxnode = 2046)

! Local Variable Declaration and Description
! ===========================================================================
        integer nnodes
        integer istart
        integer iend
        integer i
        integer lnew

        integer iyl (0:maxnode)

        complex cl (0:maxnode)

! Procedure
! ===========================================================================
        if (k .gt. kmax) then
          write (*,*) ' The value of k',k,' given to routine integral '
          write (*,*) ' is too large, k must be less than or equal to ', kmax
          stop ' k > kmax in Tintegral.f '
        else if (k .lt. 0) then
          write (*,*) ' The value of k ', k, ' given to routine '
          write (*,*) ' integral is less than 0, k must be a positive '
          write (*,*) ' number for this subroutine. '
          stop ' k < 0 in Tintegral.f '
        end if

        nnodes = 2**(k + 1) - 1
        if (k .gt. 0) then
          call tree (l2, m2, nnodes, cl, iyl)
          sum = 0
          istart = 2**k - 1
          iend = 2*istart
          do i = istart, iend
            lnew = iyl(i)
            sum = sum + cl(i)*delk(l1,lnew)
          end do
          sum = sum*delk(m1,m2)
        else
          sum = delk(l1,l2)*delk(m1,m2)
        end if

! Format Statements
! ===========================================================================
! None

        return
        end subroutine Tintegral

! Kronecker delta for l or m
! ===========================================================================
        function delk(i,j)

        integer i, j
        real delk

        if (i .eq. j) then
         delk = 1.0d0
        else
         delk = 0.0d0
        end if

        return
        end function delk

! ===========================================================================
! bessel.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       This is the spherical bessel function.
!
! ===========================================================================
! Original code written by Otto F. Sankey

! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        real function jl2 (l, x)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer l

        real x

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        if (l .gt. 6) then
         write (*,*) ' You are trying to calculate a spherical bessel '
         write (*,*) ' function for L = 7 or higher.  Currently, terms '
         write (*,*) ' for L > 6 are not considered in this function. '
         stop 'error in bessel'
        end if

! L = 0
! j0 = sinx/x
        if (l .eq. 0) then
         if (x .gt. 1.0d-4) then
          jl2 = sin(x)/x
         else
          jl2 = (1.0d0 - x**2/6.0d0 + x**4/120.0d0)
         end if
        end if

! L = 1
! j1 = sinx/x**2 - cosx/x
        if (l .eq. 1) then
         if (x .gt. 1.0d-1) then
          jl2 = sin(x)/(x**2) - cos(x)/x
         else
          jl2 = 1.0d0 - x**2/10.0d0 + x**4/280.0d0 - x**6/15120.0d0 &
     &                + x**8/1330560.0d0
          jl2 = (x/3.0d0)*jl2
         end if
        end if

! L = 2
! j2 = (3/x**3 - 1/x)*sinx - 3*cosx/x**2
        if (l .eq. 2) then
         if (x .gt. 5.0d-1) then
          jl2 = (3.0d0/x**3 - 1.0d0/x)*sin(x) - (3.0d0/x**2)*cos(x)
         else
          jl2 = 1.0d0 - x**2/14.0d0 + x**4/504.0d0 - x**6/33264.0d0 &
     &                + x**8/3459456.0d0
          jl2 = (x**2/15.0d0)*jl2
         end if
        end if

! L = 3
        if (l .eq. 3) then
         if (x .gt. 5.0d-1) then
           jl2 = (15.0d0/x**4 - 6.0d0/x**2)*sin(x) &
     &         - (15.0d0/x**3 - 1.0d0/x)*cos(x)
         else
           jl2 = 1.0d0 - x**2/18.0d0 + x**4/792.0d0 - x**6/61776.0d0 &
     &                + x**8/7413120.0d0
           jl2 = (x**3/105.0d0)*jl2
         end if
        end if

! L = 4
        if (l .eq. 4) then
         if (x .gt. 5.0d-1) then
           jl2 = (105.0d0/x**5 - 45.0d0/x**3 + 1.0d0/x)*sin(x) &
     &         - (105.0d0/x**4 - 10.0d0/x**2)*cos(x)
         else
           jl2 = 1.0d0 - x**2/22.0d0 + x**4/1144.0d0 - x**6/102960.0d0 &
     &                + x**8/14002560.0d0
           jl2 = (x**4/945.0d0)*jl2
         end if
        end if

! L = 5
        if (l .eq. 5) then
         if (x .gt. 5.0d-1) then
          jl2 = (945.0d0/x**6 - 420.0d0/x**4 + 15.0d0/x**2)*sin(x) &
     &         - (945.0d0/x**5 - 105.0d0/x**3 + 1.0d0/x)*cos(x)
         else
          jl2 = 1.0d0 - x**2/26.0d0 + x**4/1560.0d0 - x**6/159120.0d0
          jl2 = (x**5/10395.0d0)*jl2
         end if
        end if

! L = 6
        if (l .eq. 6) then
         if (x .gt. 5.0d-1) then
          jl2 = (10395.0d0/x**7 - 4725.0d0/x**5 + 210.0d0/x**3 &
     &                                          - 1.0d0/x)*sin(x) &
     &         - (10395.0d0/x**6 - 1260.0d0/x**4 + 21.0d0/x**2)*cos(x)
         else
          jl2 = 1.0d0 - x**2/30.0d0 + x**4/2040.0d0 - x**6/232560.0d0
          jl2 = (x**6/135135.0d0)*jl2
         end if
        end if

! Format Statements
! ===========================================================================
! None

        return
        end function jl2

! ===========================================================================
! tree.f
! ===========================================================================
! Program Description
! ===========================================================================
! This subroutine will expand the product of cos(theta)^k*Ylm
! where Ylm is a spherical harmonic for l and m.
! This is done by the formula:
!      cos(theta)*Yl,m = c1(l,m)*Yl+1,m + c2(l,m)*Yl-1,m
!  If one does this k times, a tree is made with each Yl,m getting two
!  offspring until you reach the desired k.
!
!                      k=0     Yl,m
!                             /    \
!             k=1  c1(l,m)Yl+1,m   c2(l,m)*Yl-1,m
!
! Each keeps branching doubly until the desired n.
!
! Note: As the branching continues, the coefficient of each child
!       is multiplied by the coefficient of the parent on down to
!       the end.
!
! Some useful formulas: The array starts at 0. The array index for the
! left child of indek i is 2*i+1. The array index for the right child
! of index i is 2*i+2.
! The array index for last parent: (n-2)/2.
! The array goes from 0 to n, where n is the number of nodes-1
!   nnodes = 2**(k+1) - 1
!
! The array index at the beginning of a row n, is given by:
!       2**n - 1
! The array index at the end of a row n, is given by:
!       2*(2**n - 1)
!
! input: l,m = l and m value for spherical harmonic being expanded
!          n = number of array elements.
! output: cl = array of values for coefficient
!        iyl = array of l-values for spherical harmonics being expanded
!
! Note: The arrays cl() and iyl() must be dimensioned in the calling
!       routine to be from zero to n or greater. Or we start overwriting
!       array values and get a bunch of garbage.
!
! ===========================================================================
! Code rewritten by:
! Kirk VanOpdorp
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-5909
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine tree (l, m, n, cl, iyl)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        integer l
        integer m
        integer n               !number of nodes

! Output:
        integer iyl (0:n)       !array of l-values

        complex cl(0:n)      !array of coefficients

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ifinal          !final parent
        integer i
        integer ileft           !left child
        integer iright          !right child
        integer lnow

! Procedure
! ===========================================================================
        cl(0) = 1
        iyl(0) = l

        ifinal = (n - 2)/2

        do i = 0, ifinal
         ileft = 2*i+1
         iright= 2*i+2
         lnow = iyl(i)
         cl(ileft) = c1(lnow,m)*cl(i)
         cl(iright) = c2(lnow,m)*cl(i)
         iyl(ileft) = iyl(i) + 1
         iyl(iright) = iyl(i) - 1
        end do

        return
        end subroutine tree

! Functions for calculating coefficients in angular integrals.
! ===========================================================================
        complex function c1(l, m)
        implicit none

! Input
        integer l, m

        integer i, j
        complex arg

        i = (l - m + 1)*(l + m + 1)
        j = (2*l + 1)*(2*l + 3)
        arg = dcmplx(dble(i)/dble(j))
        c1 = sqrt(arg)

        return
        end function c1

! ===========================================================================
        complex function c2(l,m)
        implicit none

! Input
        integer l, m

        integer i, j
        complex arg

        i = (l - m)*(l + m)
        j = (2*l - 1)*(2*l + 1)
        arg = dcmplx(dble(i)/dble(j))
        c2 = sqrt(arg)

        return
        end function c2

! ===========================================================================
        end module M_kinetic
