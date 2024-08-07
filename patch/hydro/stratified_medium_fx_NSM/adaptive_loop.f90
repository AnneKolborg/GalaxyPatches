subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,idim,ivar,info
  real(kind=8)::tt1,tt2,wallsec,dumpsec
  real(kind=4)::real_mem,real_mem_tot
  real(kind=8),save::tstart=0.0

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
  if (tstart.eq.0.0)then
      tstart=MPI_WTIME()
  endif
#endif

  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(nrestart==0)call init_refine    ! Build initial AMR grid
  if(cooling)call set_table(dble(aexp))  ! Initialize cooling look up table
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again

#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  if(myid==1)write(*,*)'Starting time integration' 

  do ! Main time loop

#ifndef WITHOUTMPI
     tt1=MPI_WTIME(info)
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax)then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),ilevel)
                 end do
              end if
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     ! Call base level
     call amr_step(levelmin,1)

     if(levelmin.lt.nlevelmax)then
        ! Hydro book-keeping
        if(hydro)then
           do ilevel=levelmin-1,1,-1
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),ilevel)
                 end do
              end if
           end do
        end if
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

#ifndef WITHOUTMPI
     tt2=MPI_WTIME(info)
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(myid==1)then
           write(*,*)'Time elapsed since last coarse step:',tt2-tt1
           call writemem(real_mem_tot)
        endif
        if(walltime_hrs.gt.0d0) then
            wallsec = walltime_hrs*3600.     ! Convert from hours to seconds
            dumpsec = minutes_dump*60.       ! Convert minutes before end to seconds
            if(wallsec-dumpsec.lt.tt2-tstart) then
               output_now=.true.
               if(myid==1) write(*,*) 'Dumping snapshot before walltime runs out'
               ! Now set walltime to a negative number so we don't keep printing outputs
               walltime_hrs = -1d0
            endif
         endif
     endif
#endif

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
