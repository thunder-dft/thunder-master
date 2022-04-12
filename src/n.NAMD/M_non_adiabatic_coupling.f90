module M_non_adiabatic_coupling

  use M_configuraciones
  use M_assemble_blocks  ! Using this to be able to use block and dblock structures

  ! Implicit none is a good idea, always
  implicit none

  type T_NAC_vars
     integer ::  nac_inpfile = 123
     integer :: ntransitions
     ! Matrix elements involving the Gradients
     integer :: nddt 
     real, allocatable :: gover(:,:,:,:,:)
     real, allocatable :: gover1c(:,:,:)
     real, allocatable :: gh_2c(:,:,:,:,:)
     real, allocatable :: gh_atm(:,:,:,:,:)
     real, allocatable :: gh_3c(:,:,:,:,:,:)
     real, allocatable :: gh_pp_otr(:,:,:,:,:)
     real, allocatable :: gh_pp_otl(:,:,:,:,:)
     real, allocatable :: gh_pp_atm(:,:,:,:,:)
     real, allocatable :: gh_pp_3c(:,:,:,:,:,:)
     real, allocatable :: gks (:, :, :, :)
     real, allocatable :: gks_old (:, :, :, :)
     real, allocatable :: gks_1(:, :, :, :)
     real, allocatable :: gks_0(:, :, :, :)

     ! numeric derivative (g*v)
     real, allocatable :: dnac (:,:)
     real, allocatable :: dnac_old (:,:)
     ! maps of lists of evolving KS-states:
     integer, allocatable :: map_ks (:)
     integer, allocatable :: map_proj (:)
     integer, allocatable :: iocc(:)
     complex, allocatable :: c_na(:,:,:)

     real, allocatable :: ratom_old(:,:)
     real, allocatable :: vatom_old(:,:)
     real, allocatable :: eigen_old(:,:)
     real, allocatable ::eigen_0(:,:) ! these eigens can be moved to density_MDET
     real, allocatable :: eigen_1(:,:)

 
     real :: tolnac = 0.0001d0

     ! Include all the variables that are exclusive for NAC
     ! Arrays should be allocatable in priciple, except for
     ! those you know their size in advance

     !             integer :: size
     !             real, allocatable :: m(:,:)

  end type T_NAC_vars

contains

  subroutine NAC_initialize(n, natoms, numorb_max, neigh_max, neighPP_max, ztot, nkpoints)

    !Arguments
    type(T_NAC_vars), intent(inout) :: n
    integer, intent(in) :: natoms, numorb_max, neigh_max, neighPP_max, nkpoints
    real, intent(in) :: ztot

    integer istate
    integer stage
    integer iband

    open (unit = n%nac_inpfile, file = 'mdet.inp', status = 'old')
    ! Note : ntransitions is equal to nele in the old code
    read (n%nac_inpfile,*) n%ntransitions
    stage = 1
    ! allocating map_ks, map_proj and iocc, these are in the nonadiabatic.f90 module in the old code, i have to create them here later
    allocate(n%map_ks(n%ntransitions))
    allocate(n%map_proj(n%ntransitions))
    allocate(n%iocc(n%ntransitions))
    ! Reading the transitions from mdet.inp file

    do istate = 1, n%ntransitions
       read (n%nac_inpfile,*) iband, n%iocc(istate)
       n%map_ks(istate) = iband
    end do
    close(n%nac_inpfile)

    ! Allocating variables gks, dnac, don't know what they do yet, these variables are in old nonadiabatic.f90 module
    allocate (n%gks(3, natoms,n%ntransitions,n%ntransitions))
    allocate (n%gks_old(3, natoms,n%ntransitions,n%ntransitions))
    allocate (n%dnac(n%ntransitions, n%ntransitions))
    allocate (n%dnac_old(n%ntransitions, n%ntransitions))
   
 ! Need allocation for imdet = 2, deal with it later

    ! we need to work with foccupy and ioccupy, they are in module density.f90, i will come back to them when I understand the reason of part

    ! Allocations
    allocate (n%gover (3, numorb_max, numorb_max, neigh_max, natoms))
    allocate (n%gover1c (3, numorb_max, numorb_max))
    allocate (n%gh_2c    (3, numorb_max, numorb_max, neigh_max, natoms))
    allocate (n%gh_atm   (3, numorb_max, numorb_max, neigh_max, natoms))
    allocate (n%gh_3c    (3, natoms, numorb_max, numorb_max, neigh_max, natoms))
    allocate (n%gh_pp_otr (3, numorb_max, numorb_max, neighPP_max, natoms))
    allocate (n%gh_pp_otl (3, numorb_max, numorb_max, neighPP_max, natoms))
    allocate (n%gh_pp_atm (3, numorb_max, numorb_max, neighPP_max, natoms))
    allocate (n%gh_pp_3c (3, natoms, numorb_max, numorb_max, neighPP_max**2, natoms))

    allocate(n%c_na(n%ntransitions,n%ntransitions,nkpoints))
    allocate(n%ratom_old(3,natoms))
    allocate(n%vatom_old(3,natoms))
    allocate(eigen_old(numorb_max,nkpoints))
    allocate(eigen_1(n%ntransitions,nkpoints))
    allocate(eigen_0(n%ntransitions,nkpoints))
    call NAC_io(n, stage)

  end subroutine NAC_initialize

  subroutine NAC_finalize(n)
    implicit none
    type(T_NAC_vars), target :: n

    deallocate(n%gover)
    deallocate(n%gover1c)
    deallocate(n%gh_2c)
    deallocate(n%gh_atm)
    deallocate(n%gh_3c)
    deallocate(n%gh_pp_otr)
    deallocate(n%gh_pp_otl)
    deallocate(n%gh_pp_atm)
    deallocate(n%gh_pp_3c)
    deallocate(n%gks)
    deallocate(n%gks_old)
    deallocate(n%dnac)
    deallocate(n%dnac_old)
    deallocate(n%map_ks)
    deallocate(n%map_proj)
    deallocate(n%iocc)
    deallocate(n%c_na)
    deallocate(n%eigen_1)
    deallocate(n%eigen_0)
    deallocate(n%eigen_old)

  end subroutine NAC_finalize


  subroutine NAC_io(n, stage)

    implicit none
    type(T_NAC_vars), intent(inout) :: n
    integer :: stage

    if (stage == 1) then
       write(*,*) 'ilogfile', ilogfile

       write (ilogfile,*)
       write (ilogfile,*) ' Reading: mdet.inp '
       write (ilogfile,*) ' Number of transitions', n%ntransitions

    else if (stage == 'degenerate bands') then
      write(ilogfile,*)'TWO EIGENVALUES VERY CLOSE'
      write(ilogfile,*)'band', iband, eigen_k(map_ks(iband),ikpoint)
      write(ilogfile,*)'band', jband, eigen_k(map_ks(jband),ikpoint)
      write(ilogfile,*)'The nonadiabatic coupling is'
      write(ilogfile,*)'NOT CALCULATED'
    end if
    !        if (stage == 2) then
    !            write (ilogfile,*)
    !            write (ilogfile,*)'Call SCF_LOOP'
    !        end if
    !        if (stage == 3) then
    !            write (ilogfile,*)
    !            write (ilogfile,8)'Call getenery'





    !         write(*,*) 'The value of the elemnt (', i, i,') is:', this%m(i,i)


    ! Here you print on screen or ilogfile information about your MD

  end subroutine NAC_io


  subroutine NAC_fileio
    ! Here you should write info in those files  with massive data
    ! They should be written on HDF5, I will teach you this later
    ! Need to use the HDF5

  end subroutine NAC_fileio

  subroutine NACs(nac_vars, den_vars, s)
    type(T_NAC_vars), intent(inout) :: nac_vars
    type(T_density_MDET), intent(in) :: den_vars
    type(T_structure), intent(in) :: s


    integer iband
    integer jband
    integer ikpoint
    integer iatom, jatom, katom
    integer ineigh
    integer natoms
    integer mbeta
    integer imu, inu
    integer in1, in2
    integer mmu, nnu



    character (len = 20) :: stage

    real diff
    real :: vec(3)
    real dot
    real cmunu

    complex pahse
    complex step1, step2

    gks = 0.0d0
    do iband =1, nac_vars%ntransitions
      do jband = 1, nac_vars%ntransitions
        do ikpoints =1 , den_vars%nkpoints !Loop over special kpoints
          if (iband .ne. jband) then
! Check for degeneracy
            diff = abs(den_vars%eigen_k(nac_vars%map_ks(iband),ikpoint) - den_vars%eigen_k(nac_vars%map_ks(jband),ikpoint) )
            if (diff .lt. nac_vars%tolnac) then
              stage == 'degenerate bands'
              call NAC_io(stage)
            else
! ===========================================================================
!                 LOPP-2c
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
! ===========================================================================
            do iatom = 1 , s%natoms
              in1 = s%atom(iatom)%imass
              do ineigh = 1 , s%neighbors(iatom)%neighn
                mbeta = s%neighbors(iatom)%neigh_b(ineigh)
                jatom = s%neighbors(iatom)%neigh_j(ineigh)
                in2 = s%neighbors(jatom)%imass
                vec =  s%atom(jatom)%ratom - s%atom(jatom)%ratom + s%xl(mbeta)%a
! TWO THINGS : 1.CREATE A FUNCTION FOR DOT PRODUCT
!              2. SPECIAL_K ARE NOT IN THE NEW FIREBALL (ACCORDING TO ARTURO)
!                 FIND A WAY CALCULATE THEM, THEY ARE IN GETKPOINTS.F90 IN THE OLD CODE
                dot = DOT(special_k(ikpoint),vec)
! NOTE from old code :! JOM : I guess we do not need now any phase (icluster.eq.1) but we may
! need it later; keep it just in case, buy without foccupy
!              phase = phasex*foccupy(map_ks(iband),ikpoint)
                phase = cmplx(cos(dot),sin(dot))* s%kpoints(ikpoint)%weight
                norb_mu = species(in1)%norb_max
                norb_nu = species(in2)%norb_max
                do imu = 1,norb_mu
! CAN NOT FIND deglect(iatom) in the new code It's easy to calculate it, it's in initbasics.f90 in the old code
                  mmu = imu + degelec(iatom)
                  step1 = phase*den_vars%bbnkre(mmu, nac_vars%map_ks(iband),ikpoint)
                  do inu = 1, norb_nu
                    nnu = inu + degelec(jatom)
                    step2 = step1*den_vars%bbnkre(nnu, nac_vars%map_ks(iband),ikpoint)
                    cmunu = real(step2)
                    nac_vars%gks(:,iatom,iband,jband) = nac_vars%gks(:,iatom,iband,jband) +    &
                    &   cmunu*( nac_vars%gover(:,imu,inu,ineigh,iatom)*den_var%eigen_k(nac_var%map_ks(jband),ikpoint)    &
                    &           - nac_vars%gh_2c(:,imu,inu,ineigh,iatom) )

                    nac_vars%gks(:,jatom,iband,jband) = nac_vars%gks(:,jatom,iband,jband) +    &
                    &   cmunu*( - nac_vars%gover(:,imu,inu,ineigh,iatom)*den_vars%eigen_k(nac_vars%map_ks(iband),ikpoint)  &
                    &           + nac_vars%gh_2c(:,imu,inu,ineigh,iatom) )
! ===========================================================================
!                LOOP to add 3-c contributions
! ===========================================================================
                    do katom = 1 , s%natoms
                      nac_vars%gks(:,katom,iband,jband) = nac_vars%gks(:,katom,iband,jband) -   &
                      &   cmunu*nac_vars%gh_3c(:,katom,imu,inu,ineigh,iatom)
                    end do ! end loop on katom
                  end do ! end loop on inu
                end do ! end loop on imu
! ===========================================================================
!               Special case : atom-case
! ===========================================================================
                do imu =1, norb_mu
                                   mmu = imu + degelec(iatom)
!              step1 = bbnkre(mmu,map_ks(iband),ikpoint)*spin
                 step1 = bbnkre(mmu,map_ks(iband),ikpoint)
                 do inu = 1, norb_nu
                  nnu = inu + degelec(iatom)
                  step2 = step1*den_vars%bbnkre(nnu,nac_vars%map_ks(jband),ikpoint)
! JOM : careful with this once we include periodicity
                  cmunu = real(step2)
! Finally the expressions.........
                  nac_vars%gks(:,iatom,iband,jband) = nac_vars%gks(:,iatom,iband,jband) +    &
                  &   cmunu*( - nac_vars%gh_atm(:,imu,inu,ineigh,iatom) )
!
                  nac_vars%gks(:,jatom,iband,jband) = nac_vars%gks(:,jatom,iband,jband) +    &
                  &   cmunu*(  nac_vars%gh_atm(:,imu,inu,ineigh,iatom) )
                end do ! end loop imu
              end do ! end loop on ineigh
            end do ! end loop on iatom
! ===========================================================================
!                 PP-neighbors-2C    
! ===========================================================================
! Loop over all atoms
            do iatom 1 , s%natoms
              in1 = s%atom(iatom)%imass
! ===========================================================================
! ONtop-Left case
! ===========================================================================
! Loop over neighbors of each iatom for ontop-Left case
              do ineigh = 1 , s%neighborsPP(iatom)%neighn
                mbeta = s%neighborsPP(iatom)%neigh_b(ineigh)
                jatom = s%neighborsPP(iatom)%neigh_j(ineigh)
                in2 = s%neighorsPP(jatom)%imass
                vec = s%atom(jatom)%ratom - s%atom(jatom)%ratom + s%xl(mbeta)%a
                dot = DOT(special_k(ikpoint),vec))
                phase = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%wight
                norb_mu = species(in1)%norb_max
                norb_nu = species(in2)%norb_max
                do imu = 1, norb_mu
                  mmu = imu + degelec(iatom)
                  step1 = phase*den_vars%bbnkre(mmu,nac_vars%map_ks(iband),ikpoint)
                  do inu = 1, norb_nu
                    nnu = inu + degelec(jatom)
                    step2 = astep1*den_vars%bbnkre(nnu,nac_vars%map_ks(iband),ikpoint)
                    cmunu = real(step2)
                    nac_vars%gks(:,iatom,iband,jband) = nac_vars%gks(:,iatom,iband,jband)+ &
                    & cmunu*(- nac_vars%gh_pp_otl(:,imu,inu,ineigh,iatom))
                    nac_vars%gks(:,iatom,iband,jband) = nac_vars%gks(:,iatom,iband,jband)+ &
                    & cmunu*(  nac_vars%gh_pp_otl(:,imu,inu,ineigh,iatom))
                  end do ! end loop on inu
                end do ! end loop on imu
              end do ! end loop over ineigh  
! ===========================================================================
! ONtop-Right case
! ===========================================================================
! Loop over neighbors of each iatom for ontop-Right case
              do ineigh = 1, s%neighborsPP(iatom)%neighn
                mbeta = s%neighborsPP(iatom)%neigh_b(ineigh)
                jatom = s%neighnorsPP(iatom)%neigh_j(ineigh)
                in2 = s%neighorsPP(jatom)%imass
                vec = s%atom(jatom)%ratom - s%atom(jatom)%ratom + s%xl(mbeta)%a
                dot = DOT(special_k(ikpoint),vec)
                phase = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
                norb_mu = species(in1)%norb_max
                norb_nu = special(in2)%norb_max
                do imu = 1 , norb_mu
                  mmu = imu + degelec(iatom)
                  step1 = phase*den_vars%bbnkre(mmu,nac_vars%map_ks(iband),ikpoint)
                  do inu = 1, norb_nu
                    nnu = inu + degelec(jatom)
                    step2 = step1*den_vars%bbnkre(nnu, nac_vars%map_ks(jband),ikpoint)
                    cmunu = real(step2)
                    ! Finally the expressions.........
                    nac_vars%gks(:,iatom,iband,jband) = nac_vars%gks(:,iatom,iband,jband) +    &
                    &   cmunu*( - nac_vars%gh_pp_otr(:,imu,inu,ineigh,iatom) )
!
                    nac_vars%gks(:,jatom,iband,jband) = nac_vars%gks(:,jatom,iband,jband) +    &
                    &   cmunu*(  nac_vars%gh_pp_otr(:,imu,inu,ineigh,iatom) )
                  end do ! end loop over inu
                end do ! end lopp over imu
              end do ! end loop over ineigh      
! ===========================================================================
! ATOM case
! ===========================================================================
            do ineigh = 1,  s%neighborsPP(iatom)%neighn
              mbeta = s%neighborsPP(iatom)%neigh_b(ineigh)
              jatom = s%neighnorsPP(iatom)%neigh_j(ineigh)
              in2 = s%neighorsPP(jatom)%imass
              vec = s%atom(jatom)%ratom - s%atom(jatom)%ratom + s%xl(mbeta)%a
              dot = DOT(special_k(ikpoint), vec)
              phase = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
              norb_mu = species(in1)%norb_max
              norb_nu = special(in2)%norb_max
              do imu = 1, norb_mu
                mmu = imu + degelec(iatom)
                step1 = phase*den_vars%bbnkre(mmu,nac_vars%map_ks(iband),ikpoint)
                do inu = 1, norb_nu
                  nnu = inu + degelec(jatom)
                  step2 = step1*den_vars%bbnkre(nnu, nac_vars%map_ks(jband),ikpoint)
                  cmunu = real(step2)
                  nac_vars%gks(:,iatom,iband,jband) = nac_vars%gks(:,iatom,iband,jband) +    &
                  &   cmunu*( - nac_vars%gh_pp_atm(:,imu,inu,ineigh,iatom) )
                  nac_vars%gks(:,jatom,iband,jband) = nac_vars%gks(:,jatom,iband,jband) +    &
                  &   cmunu*(  nac_vars%gh_pp_atm(:,imu,inu,ineigh,iatom) )
                end do ! end loop inu
              end do ! end loop imu
            end do ! end loop ineigh
            end do ! end loop over iatom 

! ===========================================================================
! PP-neighbors-3C
! ===========================================================================
! Loop over all atoms iatom in the unit cell
            do iatom = 1, s%natoms
              in1 = s%atom(iatom)%imass
              do ineigh = 1, s%neighborsPP(iatom)%neighn
                mbeta = s%neighborsPP(iatom)%neigh_b(ineigh)
                jatom = s%neighnorsPP(iatom)%neigh_j(ineigh)
                in2 = s%neighorsPP(jatom)%imass
                vec = s%atom(jatom)%ratom - s%atom(jatom)%ratom + s%xl(mbeta)%a
                dot = DOT(special_k(ikpoint), vec)
                phase = cmplx(cos(dot),sin(dot))*s%kpoints(ikpoint)%weight
                norb_mu = species(in1)%norb_max
                norb_nu = special(in2)%norb_max
                do imu = 1, norb_mu
                  mmu = imu + degelec(iatom)
                  step1 = phase*den_vars%bbnkre(mmu,nac_vars%map_ks(iband),ikpoint)
                  do inu = 1, norb_nu
                    nnu = inu + degelec(jatom)
                    step2 = step1*den_vars%bbnkre(nnu, nac_vars%map_ks(jband),ikpoint)
                    cmunu = real(step2)
                    do katom = 1, s%natom
                      nac_vars%gks(:,katom,iband,jband) = nac_vars%gks(:,katom,iband,jband) +    &
                      &   cmunu*( - nac_vars%gh_pp_3c(:,katom,imu,inu,ineigh,iatom) )
                    end do ! end loop over katom
                  end do ! end loop over inu
                end do ! end loop over imu
              end do ! end loop over ineigh
            end do ! end loop over iatom
            do katom = 1, s%natoms
              nac_vars%gks(:,katom,iband,jband) = nac_vars%gks(:,katom,iband,jband) /        &
     &    (den_vars%eigen_k(nac_vars%map_ks(iband),ikpoint) - den_vars%eigen_k(nac_vars%map_ks(jband),ikpoint) )
            end do ! end loop on katom
            end if ! end if for tolnac
          end if ! end if iband .ne jband
        end do ! end loop pn ikpoints
      end do ! end loop on jband
    end do ! end loop on iband
  stage = 'print dij in files'
  ! ask Guillermo what kind of file do we need to print these d_ij
  end subroutine NACs

  subroutine delta_t_ks(nac_vars, s , species)
  ! you can get nssh fro species(ispecies)%nssh
    type(T_NAC_vars), intent(inout) :: nac_vars
    type(T_density_MDET), intent(in) :: den_vars
    type(T_structure), intent(in) :: s
    type(T_species), intent(in) :: species(:)
    do iatom = 1, s%natoms
      r1 = nac_vars%ratom_old
      rcutoff_i = 0.0d0
      in1 = s%atom(iatom)%imass
      do imu = 1, species(in1)%nssh
        if (species(in1)%shell(imu)%rcutoff .gt. rcutoff_i) rcutoff_i = species(in1)%shell(imu)%rcutoff )
      end do ! end loop on imu
      do jatom = 1,natoms
        r2 = s%atom(jatom)%ratom
        in2 = s%atom(iatom)%imass
        z = distance (r1, r2)
        rcutoff_j = 0.0d0
        do imu =1, species(in2)%nssh
          if (species(in2)%shell(imu)%rcutoff .gt. rcutoff_j) rcutoff_j = species(in2)%shell(imu)%rcutoff )
        end do ! end loop over imu
        range = (rcutoff_i + rcutoff_j - 0.01d0)**2
        range = sqrt(range) !? why do they square it and then take the squre root ?!
        if (z .gt. range) then
          s = 0.0d0
        else 
          if (z .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
          else
            sighat = (r2 - r1)/z
          end if ! if on z
          call epsilon_function (r2,sighat,eps)
          call Depsilon_2c (r1,r2,eps,deps)
          isorb = 0
          interaction = P_overlap
          in3 = in2
          norb_mu = species(in1)%norb_max
          norb_nu = species(in2)%norb_max
          !! in the old call they call doscentros here ask arturo about this,
          !I'm folowing M_assemblr_2c_Harris, all we need is to get sx which
          !means we want crystal coordinates so there has to be a rotation after
          !getDME....
          allocate (sm (norb_mu, norb_nu)); sm = 0.0d0
          allocate (sx (norb_mu, norb_nu)); sx = 0.0d0
          call getDMEs_Fdata_2c (in1, in2, interaction, isorp, z,          &
   &                             norb_mu, norb_nu, sm, dsm)
          call rotate (in1, in3, eps, norb_mu, norb_nu, sm, sx)
          ! find out what's degelec(iatom) if it's not nessecary you can replace
          ! the following loops just using s = sx
          do inu = 1, norb_nu
            jnu = inu + degelec(jatom)
              do imu = 1, norb_mu
                jmu = imu + degelec(iatom)
                s(jmu,jnu) = sx(imu,inu)
              end do ! end loop over imu
            end do ! end loop on inu
          end if ! end if on z.gt. range
        end do ! end loop on  jatom
    end do ! end loop over iatom
! ===========================================================================
! Calculate overlap between Kohn-Sham states at different
! time steps
! Non-adiabatic term: dot pruduct sum

    sumb = 0.0d0
    do ikpoint = 1, den_vars%nkpoints
      do ij = 1, nac_vars%ntransitions
        do ik = 1, nac_vars%ntransitions
          do imu = 1, s%norbitals
            do inu = 1, s%norbitals
              sumb(ik,ij) = sumb(ik,ij) +                               &
     &        den_vars%bbnkre_old(imu,nac_vars%map_ks(ik),ikpoint)*den_vars%bbnkre(inu,nac_vars%map_ks(ij),ikpoint)*s(imu,inu)
            end do ! end loop over inu
          end do ! end loop over imu
        end do ! end loop over ik
      end do ! end loop over ij
    end do ! end loop over ikpoints

    do ij = 1, nac_vars%ntransitions
      do ik = 1, nac_vars%ntransitions
        nac_vars%dnac (ik,ij) = (sumb(ik,ij) - sumb(ij,ik))/(2.0d0*dt)
!         write (552,*) ik, ij, dnac(ik,ij), sumb(ik,ij)/dt, sumb(ij,ik)/dt
!         write(*,301) ij, ik, sumb(ij,ik)
      end do ! end loop over ik
    end do ! end loop over ij

    if (.true.) then ! I don't know what should be true here?
 ! this part calculates the V.d_ij again why?


! Calculate non-adiabatic dot-product sum using non-adiabatic couplings
! gks
! ===========================================================================
!       ddt = dt / nddt
!       deig = eigen_k - eigen_old
    !    dvatom = vatom - vatom_old

    !    dgks = gks - gks_old
!write(*,*) "DEBUG1"
! ===========================================================================
!        delta = 1.0d0
!        delta = 0.0d0
         delta = 0.5d0
!write(*,*) "DEBUG2"
! Interpolation
!        v = vatom_old + dvatom*delta
! JOM-test
         v = (ratom - ratom_old)/dt
    !    g = gks_old + dgks*delta
!write(*,*) "DEBUG3"
! JOM-test
!         g = gks_old
!        v = vatom_old
        g = gks
!        v = vatom
!        eig = eigen_old + deig*delta
! Non-adiabatic term: dot pruduct sum
         suma = 0.0d0
         do ik = 1, nele
          do ij = 1, nele
           do iatom = 1, natoms
            do ix = 1, 3
          suma(ik,ij) = suma(ik,ij) + v(ix,iatom)*g(ix,iatom,ik,ij)
            end do
           end do
          end do
         end do
!write(*,*) "DEBUG4"
! Compare both ways to calculate non-adaiabatic contribution
!write(552,*) '---------------------------------'
           do ik = 1, nele
            do ij = 1, nele
            diff = suma(ik,ij)-dnac(ik,ij)
            write(210,300)ik,ij,suma(ik,ij),dnac(ik,ij)
            end do
           end do
!        write(*,300)3,2,suma(3,2),dnac(3,2)
! Write c_na
!         do ia = 1, norbitals
!          do ik = 1, norbitals
!         write(*,100)ia,ik,cabs(c_na(ia,ik,1)),c_na(ia,ik,1)
!         end do
!        end do
end if ! end if false     


  end subroutine delta_t_ks
  
  subroutine evolve_ks_state(nac_vars,den_vars,s)
  type(T_NAC_vars), intent(inout) :: nac_vars
  type(T_density_MDET), intent(inout) :: den_vars
  type(T_structure), intent(in) :: s
  
  do iele = 1, nac_vars%ntransitions
    iband = nac_vars%map_ks(iele)
    do ikpoint = 1, den_vars%nkpoints
      eigen_1(iele,ikpoint) = eigen_k(iband,ikpoint)
      eigen_0(iele,ikpoint) = eigen_old(iband,ikpoint)
    end do ! end loop on ikpoint
  end do ! end loop on iele
  ddt = dt / nac_vars%nddt
  deig = (eigen_1 - eigen_0)/nac_vars%nddt
  dvatom = (vatom - vatom_old)/nac_vars%nddt
  dgks = (gks - gks_old)/nac_vars%nddt
  ddnac = (dnac - dnac_old)/nac_vars%nddt
! Time Loop
  do it = 1 , nddt
! step 1
! Interpolation
    v = vatom_old + dvatom*(it - 1)
    g = gks_old + dgks*(it - 1)
    eig = eigen_0 + deig*(it - 1)
    nonac = dnac + ddnac*(it - 1)
    c_aux = c_na
! Calculate derivative d/dt c_na(t)
    call dcdt_nac (v, g, nonac, eig, c_aux, dc_na, s%natoms, nac_vars%ntransitions,den_vars%nkpoints, s%norbitals)
    dc_aux = dc_na/6.0d0
! step 2
! Interpolation
    v = vatom_old + dvatom*(it - 0.5d0)
    g = gks_old + dgks*(it - 0.5d0)
    eig = eigen_0 + deig*(it - 0.5d0)
    nonac = dnac + ddnac*(it - 0.5d0)
    c_aux = c_na + dc_na*ddt*0.5d0
! Calculate derivative d/dt c_na(t)
    call dcdt_nac (v,g,nonac,eig,c_aux,dc_na,s%natoms,nac_vars%ntransitions,den_vars%nkpoints,s%norbitals)
    dc_aux = dc_aux + dc_na/3.0d0
!
! step 3
! Interpolation
    v = vatom_old + dvatom*(it - 0.5d0)
    g = gks_old + dgks*(it - 0.5d0)
    eig = eigen_0 + deig*(it - 0.5d0)
    nonac = dnac + ddnac*(it - 0.5d0)
    c_aux = c_na + dc_na*ddt*0.5d0
! Calculate derivative d/dt c_na(t)
    call dcdt_nac (v,g,nonac,eig,c_aux,dc_na,s%natoms,nac_vars%ntransitions,den_vars%nkpoints,s%norbitals)
!         call dcdt_nac (nonac,eig,c_aux,dc_na,natoms,nele,nkpoints)
! JOM-test
!       write(*,*) 'tercera llamada a dcdt, it =',it
!        write(*,*)'c_na',c_na(2,2,1),c_aux(2,2,1),dc_na(2,2,1)
    dc_aux = dc_aux + dc_na/3.0d0
! step 4
! Interpolation
    v = vatom_old + dvatom*(it)
    g = gks_old + dgks*(it)
    eig = eigen_0 + deig*(it)
    nonac = dnac + ddnac*(it - 0.5d0)
    c_aux = c_na + dc_na*ddt
! Calculate derivative d/dt c_na(t)
    call dcdt_nac (v,g,nonac,eig,c_aux,dc_na,s%natoms,nac_vars%ntransitions,den_vars%nkpoints,s%norbitals)
!         call dcdt_nac (nonac,eig,c_aux,dc_na,natoms,nele,nkpoints)
! JOM-test
!       write(*,*) 'cuarta llamada a dcdt, it =',it
!        write(*,*)'c_na',c_na(2,2,1),c_aux(2,2,1),dc_na(2,2,1)
    dc_aux = dc_aux + dc_na/6.0d0
!
!
! Integrate coefficients c_na
    nac_vars%c_na = nac_vars%c_na + dc_aux*ddt
! JOM-test
!       write(*,*) 'hacemos c_na = c_na + dc_aux*dt, it =',it
!        write(*,*)'c_na',c_na(2,2,1),c_aux(2,2,1),dc_na(2,2,1)

  end do ! end loop on it
  nac_vars%stage = 'C_na'
  call NAC_io(nac_vars%stage)
!  do iele = 1 , nac_vars%ntransitions
!    do jele = 1 , nac_vars%ntransitions
!   end do ! end loop on jele
!  end do ! end loop iele

  end subroutine evolve_ks_state  
  subroutine FSSH(nac_vars,den_vars,s, itime_step)
  type(T_NAC_vars), intent(inout) :: nac_vars
  type(T_density_MDET), intent(in) :: den_vars
  type(T_structure), intent(in) :: s
! Procedure
! ===========================================================================
    aim = cmplx(0.0d0, 1.0d0)
    a0 = cmplx(0.0d0, 0.0d0)
    a1 = cmplx(1.0d0, 0.0d0)
!----------------------------------------------------------
! Initialize seed for random_number
    call random_seed
!----------------------------------------------------------

!JOM-info : I assume that nkpoints = 1
!----------------------------------------------------------
!    if (nkpoints .gt. 1) then
!      write(*,*)'nkpoints=',nkpoints
!      write(*,*)'nkpoints is greater then 1'
!      write(*,*)'in subroutine fewest_switches'
!      write(*,*)'not ready, must stop'
!      stop
!    end if
!----------------------------------------------------------
! Calculate hopping probabilities for fewest switches
! map_ks(iele) gives back the corresponding adiabatic KS state
! we follow the possible transitions associated with states
! iele=1,nele
    do ikpoint = 1, den_vars%nkpoints
      do ij = 1, nac_vars%ntransitions
!----------------------------------------------------------
! Random numbers for Monte-Carlo
        call random_number(xrand)
!----------------------------------------------------------
! JOM-test
       nac_vars%stage = 'Random no'
       call NACio(nac_vars%stage)
!      read (2121,*) xrand   ! debug vlad
!----------------------------------------------------------
       ajj = real(conjg(nac_vars%c_na(ij,ij,den_vars%ikpoint))*c_na(ij,ij,den_vars%ikpoint))
       do ik = 1, nac_vars%ntransitions
         akj = c_na(ij,ik,den_vars%ikpoint)*conjg(c_na(ij,ij,den_vars%ikpoint))
         bkj = -2.0d0*real(conjg(akj)*nac_vars%dnac(ik,ij))
! JOM-warning: may be later we can "imporve" this by using eq(29) in JCP
! 101 4657 (1994)
!----------------------------------------------------------
!JOM-info : probability of the j ---> k transition
         prob(ik) = bkj*dt/ajj
         nac_vars%stage = 'probability'
         call NACio(nac_vars%stage)
!         write(*,*)'prob',ij,ik,prob(ik)
         if (prob(ik) .lt. 0.0d0) then
           prob(ik) = 0.0d0
         end if
!----------------------------------------------------------
! JOM-test
!          if(prob(ik) .gt. 0.0001) then
!          write(*,*)'prob',ij,ik,prob(ik)
!          write(*,*)'akk', real(conjg(c_na(ia,ik,ikpoint))*c_na(ia,ik,ikpoint))
!          write(*,*)'ajj',ajj
!          write(*,*)'akj',akj
!          write(*,*)'bkj,dt',bkj,dt
!          end if
!----------------------------------------------------------
        end do ! do ik = 1, nele
!----------------------------------------------------------
!----------------------------------------------------------
! Monte-Carlo calculation for possible transitions
! JOM-warning : we should also allow transitions to states that are not
! fully occupied (ioccupy_na = 0, 1) [from states that are occupied
! ioccupy_na = 1, 2 ]. Use iocc for this (fix later)
!         iocc (ik) = ioccupy_na (ik, ikpoint)
!----------------------------------------------------------
        call mc_switch (xrand, nele, prob, ij, ikpoint, iswitch)
!----------------------------------------------------------
! JOM-test
!       write(*,*)'nele',nele,iele,ia,ij
!       write(*,*)'xrand',xrand
!       write(*,*)'iswitch',iswitch
!----------------------------------------------------------
        if (iswitch .ne. 0) then
          nac_vars%stage = 'switch'
          call NACio(nac_vars%stage)
!----------------------------------------------------------
! perform transition ij ---> iswitch
          call transition (itime_step, ij, iswitch, ikpoint)
!----------------------------------------------------------
          return  ! we can only have one switch
        end if

      end do !nele
    end do !kpoints
! ===========================================================================
  end FSSH


!          subroutine NAC_do_something1(this)

  ! Working routines, they act on the variables defined NAC_vars and
  ! operate changes
  !            type(NAC_vars), intent(inout) :: this
  !            integer :: i

  !            do i=1, this%size

  !               this%m(i,i) = i**2

  !            end do


  !         end subroutine NAC_do_something1


  !        subroutine NAC_do_something2

  !        end subroutine NAC_do_something2

  subroutine find_neigh_max_NAC(s,neigh_max, neighPP_max) !< this should not be here but we need it for allocating gh arrays, we can change it later

  type(T_structure), intent(in) :: s
  integer neigh_max,neighPP_max ,iatom

  neigh_max = -99
  do iatom = 1, s%natoms
    neigh_max = max(neigh_max,size(s%neighbors(iatom)%neigh_j))
  end do
  do iatom = 1, s%natoms
    neighPP_max = max(neigh_max, size(s%neighbors_PP(iatom)%neigh_j))
  end do

  end subroutine find_neigh_max_NAC


  subroutine NAC_normalization(this, nkpoints)
    type(T_NAC_vars), intent(inout) :: this

    complex :: a0 = cmplx(0.0d0,0.0d0)
    complex :: a1 = cmplx(1.0d0,0.0d0)
    complex :: cnorm, caux
    integer :: ikpoint, iele, jele, nkpoints
    real :: norm

    this%c_na = a0
    do ikpoint = 1, nkpoints
       do iele = 1, this%ntransitions
          this%c_na (iele, iele, ikpoint) = a1
          !       write(*,*)'c_na',iband,c_na (iband, iband, ikpoint)
       end do
    end do

    ! check Normalization !
    do ikpoint = 1, nkpoints
       do iele = 1, this%ntransitions
          cnorm = a0
          do jele = 1, this%ntransitions
             caux = this%c_na(iele,jele,ikpoint)
             cnorm = cnorm + caux*conjg(caux)
          end do
          norm = cabs (cnorm)
          write(*,*)'Norm of initial states',iele, norm
       end do
    end do

  end subroutine NAC_normalization


end module M_non_adiabatic_coupling
