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

! M_nonadiabatic_mdet.f90
! Program Description
! ===========================================================================
!>       This routine calculates the density matrices pieces that are related
!> building the nonadiabatic coupling vectors for caclulating the pieces of
!> the time dependent Schrodinger equation for nonadiabatic molecular dynamics.
!
! ===========================================================================
! Code written by:
! James P. Lewis (with Zhaofa Li at Synfuels China Technology)
! Unit 909 of Building 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ============================================================================
! Module declaration
! ============================================================================
        module M_nonadiabatic_mdet

! /SYSTEM
        use M_configuraciones

! Type Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! initialize_mdet.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine gives the initial state for molecular dynamics with
! electronic transitions (MDET) (nonadiabatic calculation)
! ===========================================================================
! Code written by:
!> @author James P. Lewis (with Zhaofa Li at Synfuels China Technology)
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initialize_mdet (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                      !< counter over atoms
        integer iband, ikpoint             !< counter of band and kpoint
        integer iband_in                   !< counter over transitions
        integer in1                        !< species number
        integer inpfile                    !< reading from which unit
        integer issh
        integer logfile                    !< writing to which unit
        integer nfermi                     !< Fermi level state
        integer ipop

        real qztot                         !< total number of electrons
        real foccupy

        character (len = 25) :: slogfile

        type(T_transition), pointer :: piband
        type(T_kpoint), pointer :: pkpoint

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile
        inpfile = s%inpfile

! Loop over the atoms.
! Total charge - ztot
        s%ztot = 0.0d0
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          do issh = 1, species(in1)%nssh
            s%ztot = s%ztot + species(in1)%shell(issh)%Qneutral
          end do
        end do

! Initialize the fermi occupations
        qztot = s%ztot
        nfermi = int(qztot)/2
        do ikpoint = 1, s%nkpoints
          do iband = 1, nfermi
            s%kpoints(ikpoint)%foccupy(iband) = 1.0d0
            s%kpoints(ikpoint)%ioccupy(iband) = 1
          end do
        end do

! Read information from .mdet.input file
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.mdet.inp'
        open (unit = inpfile, file = slogfile, status = 'old')
        write (logfile,*)
        write (logfile,*) ' Reading from mdet.input file! '
        do ikpoint = 1, s%nkpoints       

          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
        
! Allocate transition type and initialize imap
          read (inpfile,*) pkpoint%nbands
          allocate (pkpoint%transition(pkpoint%nbands))
          do iband = 1, pkpoint%nbands
          
            ! cut some more lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)
            allocate (pkpoint%eigen_old (s%norbitals)); pkpoint%eigen_old = 0.0d0
            read (inpfile,*) iband_in, foccupy, ipop
! FIXME! We might need to read in foccupy and set ioccupy to 1 when foccupy
! is not 0 (since including all states no fixed relationship between ioccupy
! and foccupy)

            ! initialize imap
            piband%imap = iband_in
            ! initialize mdet occupations
            pkpoint%foccupy(piband%imap) = foccupy

            ! NAC initialize the dij, dij_old, and c_mdet
            allocate (piband%c_mdet(s%norbitals)); piband%c_mdet = cmplx(0.0d0, 0.0d0)
            allocate (piband%dij(3, s%natoms, pkpoint%nbands))
            piband%dij = 0.0d0
            allocate (piband%dij_old(3, s%natoms, pkpoint%nbands))
            piband%dij_old = 0.0d0

            ! MDET initialize the c_na
            allocate (piband%c_na(pkpoint%nbands))
            piband%c_na = cmplx(0.0d0, 0.0d0)
            piband%c_na(iband) = cmplx(1.0d0, 0.0d0)  
  
            if (foccupy .ge. 0.5d0) pkpoint%ioccupy(piband%imap) = 1
            write (logfile,*) ' testing imaps reach '
            write (logfile,*) piband%imap
          end do   ! end loop over bands
        end do   ! end loop over kpoints
        close (unit = inpfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine initialize_mdet


! ===========================================================================
! density_matrix_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This subroutine calculates the density matix rho and stores it in the
! structure given in M_assemble_block.f90 for molecular dynamics with electronic
! transition (MDET)
!
! ===========================================================================
        subroutine density_matrix_mdet (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iband, ikpoint           !< counter of band and kpoint
        integer logfile                  !< writing to which unit

        character (len = 25) :: slogfile

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

! ****************************************************************************
!
!       C O M P U T E    N O N A D I A B A T I C   C O E F F I C I E N T S
! ****************************************************************************
        write (logfile,*)
        write (logfile,*) ' Calculating density matrix elements for '
        write (logfile,*) ' nonadiabatic coupling vectors based on transitions '
        write (logfile,*) ' from iband to jband. '

! Loop over the special k points.
        do ikpoint = 1, s%nkpoints
          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)

! Loop over all bands
          do iband = 1, pkpoint%nbands
            ! cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            ! set the coefficient for the iband transition
            piband%c_mdet = pkpoint%c(:, piband%imap)

          end do   ! end loop over bands
        end do   ! end loop over kpoints
        
! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, i4, 2(2x, f10.6))
200     format (4(2x, f12.4))

! End Subroutine
! ===========================================================================
        return
        end subroutine density_matrix_mdet


! ===========================================================================
! writeout_density_mdet.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine is a utility to write out the density matrix.
!
! ===========================================================================
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
! Program Declaration
! ===========================================================================
        subroutine writeout_density_mdet (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint                !< counter of band and kpoint
        integer iband                  !< counter of band
        integer inu                    !< counters for mu, nu

        character (len = 25) :: slogfile

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Formatted file needed for Multimwfn
! Only for gamma
        slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
        slogfile = trim(slogfile)//'.cdcoeffs-mdet-mwfn'
        open (unit = 22, file = slogfile, status = 'replace')
        do ikpoint = 1, s%nkpoints
          write (22,"('Kpoint=',i10)") ikpoint

          ! cut some lengthy notation
          nullify (pkpoint); pkpoint => s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands
            ! cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            write (22,*)
            write (22,"('Index=',i10)") piband%imap
            write (22,"('Type=',i2)") 0
            write (22,"('Energy=',1PE16.8)") 0
            write (22,"('Occ=',f12.8)") 0
            write (22,"('Sym= ?')")
            write (22,"('$Coeff')")

! write out the coefficient
            write (22,"(5(1PE16.8))") (real(piband%c_mdet(inu)), inu = 1, s%norbitals)
          end do   ! end loop over bands
          write (22,*)
        end do   ! end loop over kpoints
        close (unit = 22)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
101     format (75('*'))
102     format (1x, ' Matrices connected to neighbor ', i3, ',',               &
     &          ' jatom = ', i4, ', ', ' mbeta = ', i4, ', ', ' d = ', f6.3)
103     format (75('='))
104     format (9f8.3)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_density_mdet


! ===========================================================================
! build_dij_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine build nonadiabatic coupliongs dij.
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
        subroutine build_dij_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
        real, parameter :: tolnac = 0.0001d0

! Variable Declaration and Description
! ===========================================================================
        integer iatom              !< counter over atoms and neighbors
        integer ikpoint                    !< counter of band and kpoint
        integer iband, jband                     !< counter of transitions

        real eigen_i, eigen_j             !< eigen values in band i and j
        real diff

        real, dimension (:, :), allocatable :: temp

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands

            ! Cut some lengthy notation
            nullify (piband); piband=>pkpoint%transition(iband)

            ! set the iband eigenvalue
            eigen_i = pkpoint%eigen(piband%imap)

            do jband = iband + 1, pkpoint%nbands

              ! Cut some lengthy notation
              nullify (pjband); pjband=>pkpoint%transition(jband)

              ! set the jband eigenvalue
              eigen_j = pkpoint%eigen(pjband%imap)

              diff = abs(eigen_i - eigen_j )
              if (diff .lt. tolnac) then
                write (s%logfile,*) ' TWO EIGENVALUES VERY CLOSE'
                write (s%logfile,*) ' iband', piband%imap, eigen_i
                write (s%logfile,*) ' jband', pjband%imap, eigen_j
                write (s%logfile,*) ' The nonadiabatic coupling is'
                write (s%logfile,*) ' NOT CALCULATED '
              else 

                piband%dij(:, :, jband) = piband%dij(:,:,jband)/(eigen_i - eigen_j)

                ! NAC force anti-symmetry for NAC
                allocate (temp (3, s%natoms)); temp = 0.0d0
                temp = - piband%dij(:, :, jband)
                pjband%dij(:, :, iband) = temp
                deallocate (temp)

              end if ! end if check for degeneracy
            end do ! end loop over jband
          end do ! end loop over iband
        end do  ! end loop over kpoints

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine build_dij_mdet


! ===========================================================================
! writeout_dij_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine writes out nonadiabatic couplings dij.
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
        subroutine writeout_dij_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom              !< counter over atoms and neighbors
        integer ikpoint                    !< counter of band and kpoint
        integer iband, jband                     !< counter of transitions
        integer logfile                     !< writing to which unit

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
        logfile = s%logfile
        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          do iband = 1, pkpoint%nbands

            ! Cut some lengthy notation
            nullify (piband); piband=>pkpoint%transition(iband)

            do jband = iband + 1, pkpoint%nbands

              ! Cut some lengthy notation
              nullify (pjband); pjband=>pkpoint%transition(jband)

              write (logfile,*)
              write (logfile,103) ' The nonadibaitc couplings: '
              write (logfile,*) piband%imap, pjband%imap
              write (logfile,100)
              write (logfile,101)
              write (logfile,100)
              do iatom = 1, s%natoms
                write (logfile,102) 'dij', iatom,                            &
     &                                     s%atom(iatom)%species%symbol,     &
        &                                  piband%dij(:, iatom, jband)
              end do 
              write (logfile,100)
              write (logfile,*)
              write (logfile,*)
              write (logfile,103) ' The nonadiabatic couplings: '
              write (logfile,*) pjband%imap, piband%imap
              write (logfile,100)
              write (logfile,101)
              write (logfile,100)
              do iatom = 1, s%natoms
                write (logfile,102) 'dij', iatom,                            &
        &                                  s%atom(iatom)%species%symbol,     &
        &                                  pjband%dij(:, iatom, iband)
              end do 
              write (logfile,100)
              write (logfile,*)
            end do ! end loop over jband
          end do ! end loop over iband
        end do  ! end loop over kpoints

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (4x, 70('='))
101     format (4x, 'dij ', 'Atom # ', 2x, ' Type ', 5x,   &
     &              ' x ', 9x, ' y ', 9x, ' z ')
102     format (4x, A,  i5, 7x, a2, 3(2x,ES10.3))
103     format (4x, A)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_dij_mdet


! ===========================================================================
! evolve_ks_states_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This routine integrates the TD equations for the coefficients c_na of the
!> TD-wfs: phi(i) = \Sum_ij c_na(i,j)*psi(j)
!> use Runge-Kutta 4th order
!>
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
        subroutine evolve_ks_states_mdet (s, itime_step)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used

        integer itime_step        

! Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: nddt = 100          !< number of electron steps

! Variable Declaration and Description
! ===========================================================================
        integer iatom              !< counter over atoms and neighbors
        integer ikpoint            !< counter of band and kpoint
        integer iband, jband       !< counter of transitions
        integer kband              !< counter of transitions
        integer iteration          !< counter of electron step
        integer ix                 !< counter of spatial dimension

        real ddt                   !< length of electron step    
        real step                  !< length of rk4 step    

        complex, dimension(:), allocatable :: dc_na, dc_tmp, c_tmp

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
        if (itime_step .eq. 1) call save_stuff_mdet (s)

        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          allocate (dc_na (pkpoint%nbands))
          allocate (dc_tmp (pkpoint%nbands))
          allocate (c_tmp (pkpoint%nbands))

          ddt = dt/nddt
          do iteration = 1, nddt
            do iband = 1, pkpoint%nbands

              ! Cut some lengthy notation
              nullify (piband); piband => pkpoint%transition(iband)

! step 1      
! Interpolation       
              c_tmp = piband%c_na
              step = float((iteration - 1.0d0)/nddt)
              call couplings  (s, ikpoint, pkpoint%nbands, step, c_tmp, dc_tmp)
              dc_na = dc_tmp/6.0d0
! step 2      
! Interpolation 
              c_tmp = piband%c_na + dc_tmp*ddt*0.5d0
              step = float((iteration - 0.5d0)/nddt)
              call couplings (s, ikpoint, pkpoint%nbands, step, c_tmp, dc_tmp)
              dc_na = dc_na + dc_tmp/3.0d0     
! step 3      
! Interpolation 
              c_tmp = piband%c_na + dc_tmp*ddt*0.5d0
              step = float((iteration - 0.5d0)/nddt) 
              call couplings (s, ikpoint, pkpoint%nbands, step, c_tmp, dc_tmp)
              dc_na = dc_na + dc_tmp/3.0d0       
! step 4      
! Interpolation 
              c_tmp = piband%c_na + dc_tmp*ddt
              step = float(iteration/nddt)
              call couplings (s, ikpoint, pkpoint%nbands, step, c_tmp, dc_tmp)
              dc_na = dc_na + dc_tmp/6.0d0 

! Integrate coefficients c_na
              piband%c_na = piband%c_na + dc_na*ddt    

            end do ! end loop over iband
          end do ! end time loop
          deallocate (c_tmp, dc_tmp, dc_na)
        end do ! end loop over kpoints
  
        call save_stuff_mdet (s)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine evolve_ks_states_mdet


! ===========================================================================
! save_stuff_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine saves some variables needed for the time
!> integration of the TD-wfs
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
        subroutine save_stuff_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                           !< counter over atoms
        integer ikpoint                         !< counter of kpoint
        integer iband, jband                    !< counter of transitions
        integer logfile                         !< writing to which unit

! NAC stuff
        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
! Save all the velocities
        do iatom = 1, s%natoms
          s%atom(iatom)%vatom_old = s%atom(iatom)%vatom
        end do

        do ikpoint = 1, s%nkpoints

          ! Cut some lengthy notation
          nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)

          pkpoint%eigen_old = pkpoint%eigen

          do iband = 1, pkpoint%nbands

            ! Cut some lengthy notation
            nullify (piband); piband => pkpoint%transition(iband)

            piband%dij_old = piband%dij
          end do ! end loop over iband
        end do ! end loop over kpoints

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine save_stuff_mdet


! ===========================================================================
! couplings
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine gives the derivative wrt time for the coefficients
!>       c_wf of the TD-wfs
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
        subroutine couplings (s, ikpoint, nbands, step, c_tmp, dc_tmp)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        integer, intent(in) :: ikpoint, nbands   !< number of bands

        real, intent(in) :: step                 !< step of rk4

        !< coefficients of wavefunctions
        complex, intent(in), dimension(nbands) :: c_tmp

! Output
        complex, intent(out), dimension (nbands) :: dc_tmp   !< derivative of coefficients

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom              !< counter over atoms
        integer jband, kband       !< counter of transitions
        integer ix                 !< counter over spatial dimension

        real djkov                 !< dot product of dij and veolocity
        real deigen, eigen_tmp     !< step and temporary eigen value in rk4
        real dv, v_tmp             !< step and temporary velocity in rk4
        real ddij, dij_tmp         !< step and temporary NAC in rk4

! NAC stuff
        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
        nullify (pkpoint); pkpoint => s%kpoints(ikpoint)

        ! Initialize dc_tmp
        dc_tmp = cmplx(0.0d0, 0.0d0)

        do jband = 1, pkpoint%nbands

          ! Cut off length
          nullify (pjband); pjband => pkpoint%transition(jband)

          deigen = pkpoint%eigen(pjband%imap) -  pkpoint%eigen_old(pjband%imap)
          eigen_tmp = pkpoint%eigen_old(pjband%imap) + deigen*step

          do kband = 1, pkpoint%nbands

            djkov = 0.0d0
            do iatom = 1, s%natoms
              do ix = 1, 3
                dv = s%atom(iatom)%vatom(ix) - s%atom(iatom)%vatom_old(ix)
                v_tmp = s%atom(iatom)%vatom_old(ix) + dv*step

                ddij = pjband%dij(ix,iatom,kband) - pjband%dij_old(ix,iatom,kband)
                dij_tmp = pjband%dij_old(ix,iatom,kband) + ddij*step

                ! dij_tmp index (jband, kband)
                djkov = djkov + v_tmp*dij_tmp
              end do ! end loop over ix
            end do ! end loop over iatom

            ! c_tmp index (iband, kband) and  dc_tmp index (iband, jband)
            dc_tmp(jband) = dc_tmp(jband) - djkov*c_tmp(kband)
          end do ! end loop over kband

          ! eigen_tmp index (pjband%imap), c_tmp index (iband, jband)
          ! and dc_tmp index (iband, jband)
          dc_tmp(jband) =                                                    &
     &      dc_tmp(jband) - (cmplx(0.0d0, 1.0d0)/hbar)*eigen_tmp*c_tmp(jband)
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
        end subroutine couplings


! ===========================================================================
! fewest_switches_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This routine determines the hoppings between Kohn-Sham state.
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
        subroutine fewest_switches_mdet (s, itime_step)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

        integer itime_step                       !< counter of nuclear step   

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                           !< counter over atoms and neighbors
        integer ikpoint                         !< counter of band and kpoint
        integer iband, jband, nbands            !< counter of transitions
        integer ix                              !< counter of spatial dimension
        integer iswitch                         !< switch band

        real djiov                              !< inner product of NAC and velocity
        real xrand                              !< random number

        !< probability of switches
        real, dimension(:), allocatable :: probability

        complex aii, aji, bji                   !< coefficients of FSSH algorithm

! NAC stuff
        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
! Initialize seed for random_number
        call random_seed

! Calculate by fewest switches algorithm gives hopping probabilities from
! the current state we follow the possible transitions associated with states
        if (s%nkpoints .gt. 1) stop ' Cannot execute FSSH for non-gamma '

        ! Cut some lengthy notation
        nullify (pkpoint); pkpoint=>s%kpoints(ikpoint)
        allocate (probability (pkpoint%nbands))

        do iband = 1, pkpoint%nbands
 
            ! Cut some lengthy notation
            nullify (piband); piband=>pkpoint%transition(iband)

            ! Random numbers for Monte-Carlo
            call random_number(xrand)
            write (s%logfile, *) 'random', xrand

            aii = real(conjg(piband%c_na(iband))*piband%c_na(iband))

            do jband = 1, pkpoint%nbands
              ! Cut some lengthy notation
              nullify (pjband); pjband=>pkpoint%transition(jband)
   
              aji = piband%c_na(jband)*conjg(piband%c_na(iband))
              djiov = 0.0d0
              do iatom = 1, s%natoms
                do ix = 1, 3
                  djiov = djiov + s%atom(iatom)%vatom(ix)*pjband%dij(ix,iatom,iband)
                end do
              end do
   
              bji = -2.0d0*real(conjg(aji)*djiov)
! JOM-warning: may be later we can "improve" this by using eq(29) in JCP 101 4657 (1994)
! JOM-info : probability of the iband ---> jband transition
              probability(jband) = bji*dt/aii
              write (s%logfile,*) ' probability = ', piband%imap, pjband%imap, probability(jband)
              if (probability(jband) .lt. 0.0d0) probability(jband) = 0.0d0
            end do ! do jband = 1, nbands

! Monte-Carlo calculation for possible transitions      
! JOM-warning : we should also allow transitions to states that are not fully occupied
! (ioccupy_na = 0, 1) [from states that are occupied ioccupy_na = 1, 2 ].
! Use iocc for this (fix later)
! iocc (jband) = ioccupy_na (jband, ikpoint)
            call mc_switch (s, iband, pkpoint%nbands, xrand, probability, iswitch)
            if (iswitch .ne. 0) then
              nullify(pjband); pjband => pkpoint%transition(iswitch)
              write(s%logfile, *) 'SWITCH!!' , piband%imap, '--->', pjband%imap
! Perform transition current band ---> switch band
              call transition (s, iband, iswitch, itime_step)
              return  ! we can only have one switch
            end if

          end do ! end loop over ibands
          deallocate (probability)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format ( 1i4, 4f8.4)
102     format ( 3i4, 5f8.4)

! End Subroutine
! ===========================================================================
        return
        end subroutine fewest_switches_mdet


! ===========================================================================
! mc_switch
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This routine determines the hoppings switches between Kohn-Sham states
!> based on a Monte-Carlo approach (see J.C. Tully, JCP 93, 1061 (1990).
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
        subroutine mc_switch (s, iband, nbands, xrand, probability, iswitch)
        implicit none

        include '../include/constants.h'        

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        integer, intent(in) :: iband
        integer, intent(in) :: nbands

        real, intent(in) ::  xrand

        real, dimension (nbands), intent(in) :: probability

! Output
        integer, intent(out) :: iswitch

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer :: iband, jband
        real :: aux

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

! Procedure
! ===========================================================================
! Check that sum of probabilities is smaller than 1
! This only works for gamma kpoint
        if (s%nkpoints .gt. 1) stop ' Cannot execute FSSH for non-gamma '
        ikpoint = 1

        ! Cut some lengthy notion
        nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
        nullify (piband); piband => pkpoint%transition(iband)

        aux = 0.0d0
        do jband = 1, nbands

          ! cut some lenghty notation
          nullify (pjband); pjband => pkpoint%transition(jband)

          ! Consider only allowed transitions
          if (pkpoint%ioccupy(piband%imap) .gt. 0) then     
            if (pkpoint%ioccupy(pjband%imap) .lt. 2) then     
              aux = aux + prob(jband)
            end if
          end if
        end do ! end loop over jband

        if (aux .gt. 1.0d0) then
          write (s%logfile, *) ' Sum of probabilities greater than 1 in mc_switches.f90 '
          write (s%logfile, *) ' total probabilty', aux, ' for state', piband%imap
          do jband = 1, nbands
            write (s%logfile, *) ' probability: ',                           &
     &                             pkpoint%transition(jband)%imap, prob(jband)
          end do
        end if

        switch = 0
        aux = 0.0d0
        do jband = 1, nbands

          ! cut some lengthy notation
          nullify (pjband); pjband => pkpoint%transition(jband)

          ! Consider only allowed transitions
          if (pkpoint%ioccupy(piband%imap) .gt. 0) then     
            if (pkpoint%ioccupy(pjband%imap) .lt. 2) then    
              aux = aux + prob(jband)
              if (aux .gt. xrand) then
                iswitch = jband
                do iband = 1, nbands
                  write (s%logfile, *) ' probability ',                      &
     &                                   pkpoint%transition(iband)%imap, probability(iband)
                end do
                exit
              end if
            end if
          end if
        end do
!----------------------------------------------------------
! If iswitch = 0, no switch (hopping) between states
! Otherwise, switch from current state (piband%imap in
! fewest_switches subroutine) to state "iswitch"
!----------------------------------------------------------

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine mc_switch


! ===========================================================================
! transition
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>     This routine performs the electronic transition
!>     current band --> switch band
!>     and rescales the velocities after the transition 
!>     to conserve total energy
!>     The velocities are re-scaled along the direction of the nonadiabatic
!>     coupling (see JCP 101, 4657 (1994) )
!
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
        subroutine transition (s, iband, iswitch, itime_step)
        implicit none

        include '../include/constants.h'        

! Argument Declaration and Description
! ===========================================================================
! Input
        type(T_structure), target :: s           !< the structure to be used.

        integer, intent(in) :: iband             !< current band
        integer, intent(in) :: iswitch           !< switch band
        integer, intent(in) :: itime_step        !< counter of nuclear step

! Parameters and Data Declaration
! ===========================================================================
        real, parameter :: tolaa = 1.0d-08

! Variable Declaration and Description
! ===========================================================================
        integer iatom                            !< counter of atomn
        integer ix                               !< counter of spatial dimension

        real aa, bb, cc, alfa                    !< coefficients of quadratic equation
        real ejump                               !< energy difference in eV
        real energy, tkinetic                    !< potential and kinetics energy
        real etot, etot_before, etot_after       !< For other theories

!       NAC stuff
        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband
        type(T_transition), pointer :: pjband

        integer itheory

! Procedure
! ===========================================================================
! Switch from current band ---> switch band
! This only works for gamma kpoint
        if (s%nkpoints .gt. 1) stop ' Cannot execute FSSH for non-gamma '
        ikpoint = 1

        ! cut some lengthy notation
        nullify (pkpoint); pkpoint => s%kpoints(ikpoint)
        nullify (piband); piband => pkpoint%transition(iband)
        nullify (pjband); pjband => pkpoint%transition(iswitch)

        pkpoint%ioccupy(piband%imap) = pkpoint%ioccupy(piband%imap) - 1
        pkpoint%foccupy(piband%imap) = pkpoint%foccupy(piband%imap) - 0.50d0
        pkpoint%ioccupy(pjband%imap) = pkpoint%ioccupy(pjband%imap) + 1
        pkpoint%foccupy(pjband%imap) = pkpoint%foccupy(pjband%imap) + 0.50d0

! Calculate energy jump
        itheory = 0
        if (itheory .eq. 0) then
          ejump = pkpoint%eigen(pjband%imap) - pkpoint%eigen(piband%imap)
        else
          etot_before = etot
          call scf_loop(itime_step)
          call getenergy(itime_step)
          etot_after = etot
          ejump = etot_after - etot_before
          write (s%logfile,*) 'ETOT-BEFORE = ', etot_before
          write (s%logfile,*) 'ETOT-AFTER  = ', etot_after
        end if

! Transform energy shift from eV to atomic units ( amu*(angs/fs)**2 )
        write (s%logfile,*) 'ejump (eV) = ', ejump
        energy = ejump*P_fovermp

        write (s%logfile,*) 'energy (dynamical units) = ', energy

! ===========================================================================
! Find out if transition current band --> switch band is accesible 
! (i.e. if there is enough kinetic energy) 
        aa = 0.0d0
        bb = 0.0d0
        do iatom = 1, s%natoms
          do ix = 1, 3
            aa = aa + 0.50d0*piband%dij(ix,iatom,iswitch)                    &
        &                   *piband%dij(ix,iatom,iswtich)/s%atom(iatom)%imass
            bb = bb + s%atom(iatom)%vatom(ix)*piband%dij(ix,iatom,iswitch)
          end do
        end do
        write (s%logfile,*) ' aa, 4*energy*aa =', aa, 4.0d0*aa*energy
        write (s%logfile,*) ' bb, bb**2 = ', bb, bb**2
        cc = bb**2 - 4.0d0*aa*energy
        write (s%logfile,*) ' cc = ', cc

        if (aa .gt. tolaa) then
! The transition is accepted
          if (cc .ge. 0.0d0) then
            write (s%logfile,*) ' transition accepted '
            if (bb .ge. 0.0d0) then
              alfa = (bb - sqrt(cc))/(2.0d0*aa)
            else
              alfa = (bb + sqrt(cc))/(2.0d0*aa)
            end if
          else !(cc .ge. 0.0d0)

! The transition is NOT accepted
            write (s%logfile,*) 'transition NOT accepted '
 
! Revert transition written by vlada
            pkpoint%ioccupy(piband%imap) = pkpoint%ioccupy(piband%imap) + 1
            pkpoint%foccupy(piband%imap) = pkpoint%foccupy(piband%imap) + 0.50d0
            pkpoint%ioccupy(pjband%imap) = pkpoint%ioccupy(pjband%imap) - 1
            pkpoint%foccupy(pjband%imap) = pkpoint%foccupy(pjband%imap) - 0.50d0

! Define alfa for re-scaling velocities
! JOM-info : two possibilities: a) Do nothing; b) Reflection
! a) Do nothing
 
! b) Reflection (see JCP 101, 4657 (1994), pag. 4664)
            alfa = bb/aa 
        
          end if !(cc .ge. 0.0d0)
        else  !(aa .gt. tolaa)
          alfa = 0.0d0
        end if  !(aa .gt. tolaa)
        write (s%logfile, *) ' alfa = ', alfa
          
        write (s%logfile, *) ' Velocity before Rescaling'
        do iatom = 1, s%natoms
          write (s%logfile, *) (s%atom(iatom)%vatom(ix), ix = 1,3)
        end do
 
        tkinetic = 0.0d0
        do iatom = 1, s%natoms
          tkinetic = tkinetic                                                 &
     &              + (0.5d0/P_fovermp)*s%atom(iatom)%imass                   &
     &               *(s%atom(iatom)%vatom(1)**2 +  s%atom(iatom)%vatom(2)**2 &
     &                                           +  s%atom(iatom)%vatom(3)**2)
        end do
        write (s%logfile, *) 'KINETIC = ', tkinetic
 
! RESCALING VELOCITIES
        do iatom = 1, s%natoms
          do ix = 1, 3
            s%atom(iatom)%vatom(ix) = s%atom(iatom)%vatom(ix)                 &
     &                               - alfa*piband%dij(ix,iatom,iswitch)/s%atom(iatom)%imass
          end do
        end do

        write (s%logfile, *) ' Velocity After Rescaling'
        do iatom = 1, s%natoms
          write (s%logfile, *) (s%atom(iatom)%vatom(ix), ix = 1, 3)
        end do
 
        tkinetic = 0.0d0
        do iatom = 1, s%natoms
         tkinetic = tkinetic                                                  &
     &              + (0.5d0/P_fovermp)*s%atom(iatom)%imass                   &
     &               *(s%atom(iatom)%vatom(1)**2 +  s%atom(iatom)%vatom(2)**2 &
     &                                           +  s%atom(iatom)%vatom(3)**2)
        end do
        write(s%logfile, *) ' KINETIC = ', tkinetic
 
! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine transition


! destroy_denmat_mdet
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing denmat_mdet
!! information.
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
        subroutine destroy_denmat_mdet (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom              !< counter over atoms and neighbors
        integer ikpoint                    !< counter of band and kpoint
        integer iband                     !< counter of transitions

        type(T_kpoint), pointer :: pkpoint
        type(T_transition), pointer :: piband

! Procedure
! ===========================================================================
        ! destroy the density matrix pieces - forces are already evaluated
        do ikpoint = 1, s%nkpoints
          do iband = 1, s%kpoints(ikpoint)%nbands
            deallocate (s%kpoints(ikpoint)%transition(iband)%c_mdet)
            deallocate (s%kpoints(ikpoint)%transition(iband)%c_na)
            deallocate (s%kpoints(ikpoint)%transition(iband)%dij)
            deallocate (s%kpoints(ikpoint)%transition(iband)%dij_old)
          end do
          deallocate (s%kpoints(ikpoint)%transition)
          deallocate (s%kpoints(ikpoint)%eigen_old)
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
        end subroutine destroy_denmat_mdet

! End Module
! ===========================================================================
        end module M_density_matrix_mdet
