recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
#ifdef RT
  use rt_hydro_commons
  use SED_module
  use UV_module
  use coolrates_module, only: update_coolrates_tables
  use rt_cooling_module, only: update_UVrates
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar
  logical::ok_defrag
  logical,save::first_step=.true.

  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  if(levelmin.lt.nlevelmax .and..not. static)then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then
              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)
              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
#ifdef SOLVERmhd
                 do ivar=1,nvar+3
#else
                 do ivar=1,nvar
#endif
                    call make_virtual_fine_dp(uold(1,ivar),i)
#ifdef SOLVERmhd
                 end do
#else
                 end do
#endif
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
#ifdef RT
              if(rt)then
                 do ivar=1,nrtvar
                    call make_virtual_fine_dp(rtuold(1,ivar),i)
                 end do
                 if(simple_boundary)call rt_make_boundary_hydro(i)
              end if
#endif
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if
           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0.and.first_step)then
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if

  !-----------------
  ! Update sink cloud particle properties
  !-----------------
  if(sink)call update_cloud(ilevel)

  !-----------------
  ! Particle leakage
  !-----------------
  if(pic)call make_tree_fine(ilevel)
  
  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout).or.output_now.EQV..true.)then
        if(.not.ok_defrag)then
           call defrag
        endif

        call dump_all

        ! Run the clumpfinder, (produce output, don't keep arrays alive on output)
        if(clumpfind .and. ndim==3) call clump_finder(.true.,.false.)
        ! Dump lightcone
        if(lightcone) call output_cone()

        if (output_now.EQV..true.) then
          output_now=.false.
        endif

     endif
     ! Important can't be done in sink routines because it must be done after dump all
     if(sink)acc_rate=0.
  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(imov.le.imovout)then ! ifort returns error for next statement if looking
                             ! beyond what is allocated as an array amovout/tmovout
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
           call output_frame()
        endif
     endif
  end if

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  if(ilevel==levelmin)then
     !----------------------------------------------------
     ! Kinetic feedback from giant molecular clouds
     !----------------------------------------------------
     if(hydro.and.star.and.eta_sn>0.and.f_w>0)call kinetic_feedback
  endif

  !$acc enter data copyin(active,reception,emission,boundary) if(ilevel==levelmin.or.icount>1)
  !$acc update device(son,nbor,father,xg,cpu_map)             if(ilevel==levelmin.or.icount>1)
  call acc_communication_buffers(ilevel)
#if defined(_OPENACC)
  if(ilevel==levelmin.or.icount>1) then
  do i=ilevel,nlevelmax
     do ivar=1,nvar
        call update_globalvar_dp_to_device(uold(1,ivar),i)
     end do
     
     call update_globalvar_dp_to_device(phi,i)

     do idim=1,ndim
        call update_globalvar_dp_to_device(f(1,idim),i)
     end do
  end do
  end if
#endif

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)
     call rho_fine(ilevel,icount)
  endif

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic)then
     ! Remove particles to finer levels
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end if

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then

     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif     

     ! Compute gravitational potential
     if(ilevel>levelmin)then
        if(ilevel .ge. cg_levelmin) then
           print*,"If you want to use phi_fine_cg you must port it in GPU!"
           print*,"Entering here with ACCyes, then WRONG RESULTS!!"
           call phi_fine_cg(ilevel,icount)
        else
           call multigrid_fine(ilevel,icount)
        end if
     else
        call multigrid_fine(levelmin,icount)
     end if
     !when there is no old potential...
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     ! Synchronize remaining particles for gravity
     if(pic)then
#if defined(_OPENACC)
        do idim=1,ndim
           call update_globalvar_dp_to_host(f(1,idim),ilevel)
        end do
#endif
        ! This subroutine needs f for the particles
        call synchro_fine(ilevel)
     end if    

     if(hydro)then

        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

        ! Update boundaries
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp_acc(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)

        ! Compute Bondi-Hoyle accretion parameters
        if(sink)call collect_acczone_avg(ilevel)

     end if
  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#endif

  !----------------------
  ! Compute new time step
  !----------------------
#if defined(_OPENACC)
  ! The courant_fine subroutine (from newdt_fine) needs uold & f
  do ivar=1,nvar
     call update_globalvar_dp_to_host(uold(1,ivar),ilevel)
  end do
#endif

  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
  if(hydro)call set_unew(ilevel)

#ifdef RT
  ! Set rtunew equal to rtuold
  if(rt)call rt_set_unew(ilevel)
#endif

#if defined(_OPENACC)
  ! Copies that later must be eliminated
  call update_globalvar_dp_to_host(phi,ilevel)
  call update_globalvar_dp_to_host(phi_old,ilevel)
#endif

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           !$acc exit data delete(active,reception,emission,boundary)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
        if(sink)call update_sink(ilevel)
     end if
  else
     call update_time(ilevel)
     if(sink)call update_sink(ilevel)
  end if

  ! Thermal feedback from stars
  if(hydro.and.star.and.eta_sn>0)call thermal_feedback(ilevel)

#ifdef RT
  ! Add stellar radiation sources
  !if(rt.and.rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))
#endif

  ! Density threshold or Bondi accretion onto sink particle
  if(sink)then
     call grow_sink(ilevel,.false.)
  end if

  !-----------
  ! Hydro step
  !-----------
  if(hydro)then

     ! Hyperbolic solver
#if defined(_OPENACC) || defined(USE_ACC_VERSION)
     call godunov_fine_acc(ilevel)
#else
     call godunov_fine(ilevel)
#endif
     ! Reverse update boundaries
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_reverse_dp_acc(unew(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(pressure_fix)then
        call make_virtual_reverse_dp_acc(enew(1),ilevel)
        call make_virtual_reverse_dp_acc(divu(1),ilevel)
     endif

     ! Set uold equal to unew
     call set_uold(ilevel)

     ! ! Density threshold or Bondi accretion onto sink particle
     ! if(sink)then
     !    !this is a trick to temporarily solve the issue with sink accretion 
     !    !from ghost zones. Only an option for simulations without dark matter.
     !    if (.not. cosmo)then
     !       call make_tree_fine(ilevel)
     !       call virtual_tree_fine(ilevel)
     !       ! assuming all sink cloud parts sit on levelmax 
     !       ! it's better to compute the accretion_rate based on
     !       ! the updated values
     !       call collect_acczone_avg(ilevel)
     !    end if
     !    call grow_sink(ilevel,.false.)
     ! end if

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step 
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

#if defined(_OPENACC)
     ! copy to CPU for the non-ported parts (only the specified level)
     call copy_hydro_to_host(ilevel)
#endif

     ! Restriction operator
     call upload_fine(ilevel)
  endif

#ifdef RT
  ! here were do the main RT calls
  if(rt .and. rt_advect) then
     call rt_step(ilevel)
  else
     ! Still need a chemistry call if RT is defined but not
     ! actually doing radiative transfer (i.e. rt==false):
     if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
  endif
  ! Regular updates and book-keeping:
  if(ilevel==levelmin) then
     if(cosmo) call update_rt_c
     if(cosmo .and. haardt_madau) call update_UVrates(aexp)
     if(cosmo .and. rt_isDiffuseUVsrc) call update_UVsrc
     if(cosmo) call update_coolrates_tables(dble(aexp))
     if(ilevel==levelmin) call output_rt_stats
  endif
#else
  if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
#endif

  !---------------
  ! Move particles
  !---------------
  if(pic)then
     call move_fine(ilevel) ! Only remaining particles
  end if

  !-------------------------------
  ! Source term in leaf cells only
  !-------------------------------


  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
  if(hydro.and.star)call star_formation(ilevel)

  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if(hydro)then
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif
#ifdef RT
  ! Here was the RT boundary stuff
#endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
 if(hydro)then
     if(eta_mag>0d0.and.ilevel==levelmin)then
        call diffusion
     endif
  end if
#endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
  if(.not.static) call flag_fine(ilevel,icount)


  !----------------------------
  ! Merge finer level particles
  !----------------------------
  if(pic)call merge_tree_fine(ilevel)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin)then
     call rad_step(dtnew(ilevel))
  endif
#endif

  if(sink)then
     !-------------------------------
     ! Update coarser level sink velocity
     !-------------------------------
     if(ilevel>levelmin)then
        vsold(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel-1)
        if(nsubcycle(ilevel-1)==1)vsnew(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel)
        if(icount==2)vsnew(1:nsink,1:ndim,ilevel-1)= &
             (vsold(1:nsink,1:ndim,ilevel)*dtold(ilevel)+vsnew(1:nsink,1:ndim,ilevel)*dtnew(ilevel))/ &
             (dtold(ilevel)+dtnew(ilevel))
     end if
     !---------------
     ! Sink production
     !---------------
     if(ilevel==levelmin)call create_sink
  end if

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step

#ifdef RT
!*************************************************************************
subroutine rt_step(ilevel)

!  Radiative transfer and chemistry step. Either do one step on ilevel,
!  with radiation field updates in coarser level neighbours, or, if
!  rt_nsubsteps>1, do many substeps in ilevel only, using Dirichlet
!  boundary conditions for the level boundaries. 
!-------------------------------------------------------------------------
  use amr_parameters, only: dp
  use amr_commons,    only: levelmin, dtnew, myid
  use rt_parameters, only: rt_isDiffuseUVsrc
  use rt_cooling_module, only: update_UVrates
  use rt_hydro_commons
  use UV_module
  use SED_module,     only: star_RT_feedback
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ilevel
  real(dp) :: dt_hydro, t_left, dt_rt
  integer  :: i_substep, ivar
!-------------------------------------------------------------------------
  dt_hydro = dtnew(ilevel)                   ! Store hydro timestep length
  t_left = dt_hydro

  i_substep = 0
  do while (t_left > 0)                      !                RT sub-cycle
     i_substep = i_substep + 1
     call get_rt_courant_coarse(dt_rt)
     ! Temporarily change timestep length to rt step:
     dtnew(ilevel) = MIN(t_left, dt_rt/2.0**(ilevel-levelmin))

     !if (myid==1) write(*,900) dt_hydro, dtnew(ilevel), i_substep, ilevel    
     if (i_substep > 1) call rt_set_unew(ilevel)

     if(rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))

     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
     do ivar=1,nrtvar
        call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
     end do

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

     if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)

     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)

     t_left = t_left - dtnew(ilevel)
  end do                                   !          End RT subcycle loop
  dtnew(ilevel) = dt_hydro                 ! Restore hydro timestep length

  ! Restriction operator to update coarser level split cells
  call rt_upload_fine(ilevel)

  if (myid==1 .and. rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep

900 format (' dt_hydro=', 1pe12.3, ' dt_rt=', 1pe12.3, ' i_sub=', I5, ' level=', I5)
901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')

end subroutine rt_step
#endif


