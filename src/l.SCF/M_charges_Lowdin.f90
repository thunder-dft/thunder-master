! copyright info:
!
!                             @Copyright 2013
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

! M_charges_Lowdin
! Module Description
! ===========================================================================
!       This is a module containing all coordinate information and atomic
! information such as the charges, etc.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Box 6315, 135 Willey St.
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-5141 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Module Declaration
! ===========================================================================
        module M_charges
        use M_species
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
! calculate_charges.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine calculates the Lowdin population analysis for then
! determining the charges for the atoms.
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
        subroutine calculate_charges (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                      !< counter over the atoms
        integer ikpoint                    !< counter over kpoints
        integer iorbital                   !< counter over orbitals
        integer imu                        !< another counting to find spot
        integer in1                        !< species number for iatom
        integer inpfile                    !< reading from which unit
        integer issh, jssh                 !< counter over shells
        integer nssh                       !< number of shells
        integer mmu                        !< spot in the array - block_slot

        real aux1, aux2, aux3              !< temporary multipliers

        character (len = 25) :: slogfile

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        inpfile = s%inpfile

! Loop over the atoms.
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          s%atom(iatom)%Q = 0.0d0
          do issh = 1, species(in1)%nssh
            s%atom(iatom)%shell(issh)%Qout = 0.0d0
          end do

! Loop over the special k points.
          do ikpoint = 1, s%nkpoints
            aux1 = s%kpoints(ikpoint)%weight*P_spin
            do iorbital = 1, s%norbitals
              if (s%kpoints(ikpoint)%ioccupy(iorbital) .eq. 1) then
                aux2 = aux1*s%kpoints(ikpoint)%foccupy(iorbital)

! Finally the imu loop.
                imu = 0
                do issh = 1, species(in1)%nssh
                  do jssh = 1, 2*species(in1)%shell(issh)%lssh + 1
                    imu = imu + 1
                    mmu = imu + s%iblock_slot(iatom)
                    aux3 = aux2*cabs(s%kpoints(ikpoint)%c_Lowdin(mmu,iorbital))**2
                    s%atom(iatom)%shell(issh)%Qout =                           &
     &                s%atom(iatom)%shell(issh)%Qout + aux3
                    s%atom(iatom)%Q = s%atom(iatom)%Q + aux3
                  end do
                end do
              end if
            end do  ! end loop over orbitals
          end do  ! end loop over kpoints

! End loop over atoms
        end do

! Writout the charges to a .CHARGES file
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.CHARGES'
        open (unit = inpfile, file = slogfile, status = 'unknown')
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          nssh = species(in1)%nssh
          write (inpfile,*) (s%atom(iatom)%shell(issh)%Qin, issh = 1, nssh)
        end do
        close (unit = inpfile)

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
        return
        end subroutine calculate_charges


! ===========================================================================
! writeout_charges.f90
! ===========================================================================
! Subroutine Description
! ===========================================================================
!       This routine writes out the charges for each atom - both total and
! charge per shell.
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
        subroutine writeout_charges (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                       !< counter over the atoms
        integer in1                         !< species number for iatom
        integer issh                        !< counter over shells
        integer logfile                     !< writing to which unit

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! Initialize logfile
        logfile = s%logfile

        write (logfile,*)
        write (logfile,*) ' LOWDIN CHARGES (by shell): '
        write (logfile,500)
        write (logfile,501)
        write (logfile,500)
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          write (logfile,502) iatom, species(in1)%symbol, species(in1)%nssh, &
     &      (s%atom(iatom)%shell(issh)%Qout, issh = 1, species(in1)%nssh)
        end do

        write (logfile,500)
        write (logfile,*)
        write (logfile,*) ' Total Lowdin charges for each atom: '
        write (logfile,500)
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          write (logfile,503) iatom, species(in1)%symbol, s%atom(iatom)%Q
        end do
        write (logfile,500)
        write (logfile,*) '  '

! Format Statements
! ===========================================================================
500     format (70('='))
501     format (2x, ' Atom # ', 2x, ' Type ', 2x, ' Shells ', 1x,' Output Charges ')
502     format (3x, i5, 7x, a2, 5x, i2, 4x, 8(1x, f5.2))
503     format (3x, i5, 7x, a2, 4x, f10.4)

! End Subroutine
! ===========================================================================
        return
        end subroutine writeout_charges


! ===========================================================================
! calculate_populations
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine calculates the Lowdin population analysis and writes
! these out for each band.  This population analysis can be used to calculate
! things like W (number of accessible atoms), etc.
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
        subroutine calculate_populations (s)
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh              !< counter over atoms and neighbors
        integer iband                      !< counter over the bands
        integer ikpoint                    !< counter over kpoints
        integer imu, inu                   !< counter over shells
        integer in1, in2                   !< species number for iatom
        integer inpfile                    !< reading from which unit
        integer jatom                      !< neighbor of iatom
        integer mmu, nnu                   !< spot in the array - block_slot
        integer num_neigh                  !< number of neighbors
        integer mbeta                      !< cell containing neighbor of iatom
        integer norb_mu, norb_nu           !< size of the block for the pair

        real dot                           !< dot product between K and r
        real gutr
        real step1, step2
        real W                             !< localization parameter

        real, dimension (3) :: r1, r2    !< positions of iatom and jatom
        real, dimension (3) :: sks       !< k point value
        real, dimension (3) :: vec

        real, dimension (:, :, :), allocatable :: pcharge

        complex phase, phasex            !< phase between K and r

        character (len = 25) :: slogfile

        type(T_assemble_block), pointer :: pS_neighbors
        type(T_assemble_neighbors), pointer :: poverlap

! Allocate Arrays
! ===========================================================================
        allocate (pcharge (s%nkpoints, s%norbitals_new, s%natoms)); pcharge = 0.0d0

! Procedure
! ===========================================================================
! Initialize logfile
        inpfile = s%inpfile

! Writout the populations to a .pop file
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.pop'
        open (unit = inpfile, file = slogfile, status = 'unknown', form = 'unformatted')

! Loop over the special k points.
        write (inpfile) s%nkpoints, s%norbitals_new, s%natoms
        do ikpoint = 1, s%nkpoints

! Loop over the bands - we do this first and then we loop over atoms.
          do iband = 1, s%norbitals_new

! Loop over the atoms and the neighbors
            do iatom = 1, s%natoms
              ! cut some lengthy notation
              poverlap=>s%overlap(iatom)

              r1 = s%atom(iatom)%ratom
              in1 = s%atom(iatom)%imass
              norb_mu = species(in1)%norb_max
              num_neigh = s%neighbors(iatom)%neighn

              ! Initialize pcharge
              pcharge(ikpoint, iband, iatom) = 0.0d0

              do ineigh = 1, num_neigh
                ! cut some more lengthy notation
                pS_neighbors=>poverlap%neighbors(ineigh)

                mbeta = s%neighbors(iatom)%neigh_b(ineigh)
                jatom = s%neighbors(iatom)%neigh_j(ineigh)
                r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
                in2 = s%atom(jatom)%imass
                norb_nu = species(in2)%norb_max

                ! Find the phase which is based on k*r
                vec = r2 - r1
                sks = s%kpoints(ikpoint)%k
                dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
                phasex = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
                phase = phasex*s%kpoints(ikpoint)%foccupy(iband)

! Finally the imu loop.
                do imu = 1, norb_mu
                  mmu = imu + s%iblock_slot(iatom)
                  step1 = phasex*conjg(s%kpoints(ikpoint)%c(mmu,iband))

                  do inu = 1, norb_nu
                    nnu = inu + s%iblock_slot(jatom)
                    step2 = step1*s%kpoints(ikpoint)%c(nnu,iband)
                    gutr = real(step2)
                    pcharge(ikpoint, iband, iatom) =                                &
     &                pcharge(ikpoint, iband, iatom) + gutr*pS_neighbors%block(imu,inu)
                  end do
                end do
              end do  ! end loop over neighbors
              write (inpfile) pcharge(ikpoint, iband, iatom)
            end do  ! end loop over atoms
          end do  ! end loop over bands
        end do  ! end loop over kpoints
        close (unit = inpfile)

! Calculate W - number of accessible atoms in each band
        slogfile = s%basisfile(:len(trim(s%basisfile))-4)
        slogfile = trim(slogfile)//'.W'
        open (unit = inpfile, file = slogfile, status = 'unknown')
        do ikpoint = 1, s%nkpoints
          do iband = 1, s%norbitals_new
            W = 0.0d0
            do iatom = 1, s%natoms
              W = W - pcharge(ikpoint, iband, iatom)*log(abs(pcharge(ikpoint, iband, iatom)))
            end do  ! end loop over atoms
            W = exp(W)
            write (inpfile, *) s%kpoints(ikpoint)%eigen(iband), W
          end do  ! end loop over bands
        end do  ! end loop over kpoints
        close (unit = inpfile)

! Format Statements
! ===========================================================================
! None

! Deallocate Arrays
! ===========================================================================
        deallocate (pcharge)

! End Subroutine
! ===========================================================================
        return
        end subroutine calculate_populations


! ===========================================================================
! destroy_charges
! ===========================================================================
! Subroutine Description
! ===========================================================================
!>       This routine deallocates the arrays containing the assemble_2c
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
        subroutine destroy_charges (s)
        implicit none

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s           !< the structure to be used.

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iatom                             !< counter over atoms

! Procedure
! ===========================================================================
        do iatom = 1, s%natoms
          deallocate (s%atom(iatom)%shell)
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
        end subroutine destroy_charges

! End Module
! ===========================================================================
        end module M_charges
