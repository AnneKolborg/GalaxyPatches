!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine dump_all
   !-----------------------------------------------------------------------------------
   ! SUBROUTINE: dump_all
   !
   ! DESCRIPTION:
   !   This subroutine handles the creation and storage of various output files related
   !   to the state of the system being simulated. It creates directories and files based 
   !   on various conditions and writes data into them. Outputs include headers, amr, 
   !   hydro, rt, parts, poisson, radiation, gadget format, and timers.
   !
   ! MODULE DEPENDENCIES:
   !   - amr_commons     : Common AMR (Adaptive Mesh Refinement) variables and methods.
   !   - pm_commons      : Common PM (Particle Mesh) variables and methods.
   !   - hydro_commons   : Common hydrodynamic variables and methods.
   !   - cooling_module  : Cooling-related variables and methods.
   !
   ! GLOBALS MODIFIED:
   !   - None directly by this subroutine; dependent subroutines called may modify globals.
   !
   ! INPUT/OUTPUT ARGUMENTS:
   !   - No explicit arguments; uses module variables and parameters.
   !
   ! NOTES:
   !   1. If the subroutine is compiled with MPI support (WITHOUTMPI not defined), MPI
   !      functionalities are invoked to ensure synchronization across processes.
   !   2. Data is written to directories structured with 'checkpoint_' naming conventions
   !      and various descriptive subdirectories.
   !   3. The subroutine checks the status of certain flags (e.g., hydro) to determine which
   !      modules' data needs to be written to the output.
   !   4. Verbose outputs can be triggered to notify the user about the progress of file 
   !      backups.
   !
   ! EXAMPLES OF USAGE:
   !   call dump_all
   !-----------------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
#if ! defined (WITHOUTMPI) || defined (NOSYSTEM)
  integer::info
#endif
  character::nml_char
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::filename,filename_desc,filedir
  integer::ierr

  if(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0)return
  if(nstep_coarse==0.and.nrestart>0)return
  if(verbose)write(*,*)'Entering dump_all'

  call write_screen
  call title(ifout,nchar)
  ifout=ifout+1
  if(t>=tout(iout).or.aexp>=aout(iout))iout=iout+1
  output_done=.true.

  if(IOGROUPSIZEREP>0)call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)

  if(ndim>1)then

     if(IOGROUPSIZEREP>0) then
        filedir='checkpoint_'//TRIM(nchar)//'/checkpoint_'//TRIM(ncharcpu)//'/'
     else
        filedir='checkpoint_'//TRIM(nchar)//'/'
     endif

     call create_output_dirs(filedir)
   
     if(myid==1.and.print_when_io) write(*,*)'Start backup header'
     ! Output header: must be called by each process !
     filename=TRIM(filedir)//'header_'//TRIM(nchar)//'.txt'
     call output_header(filename)
#ifndef WITHOUTMPI
     if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(myid==1.and.print_when_io) write(*,*)'End backup header'

     if(myid==1.and.print_when_io) write(*,*)'Start backup info etc.'
     ! Only master process
     if(myid==1)then
        filename=TRIM(filedir)//'info_'//TRIM(nchar)//'.txt'
        call output_info(filename)
        filename=TRIM(filedir)//'makefile.txt'
        call output_makefile(filename)
        filename=TRIM(filedir)//'patches.txt'
        call output_patch(filename)
        if(cooling .and. .not. neq_chem)then
           filename=TRIM(filedir)//'cooling_'//TRIM(nchar)//'.out'
           call output_cool(filename)
        end if
        if(sink)then
           filename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.info'
           call output_sink(filename)
           filename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.csv'
           call output_sink_csv(filename)
        endif
        ! Copy namelist file to output directory
        filename=TRIM(filedir)//'namelist.txt'
        OPEN(10, FILE=namelist_file, ACCESS="STREAM", ACTION="READ")
        OPEN(11, FILE=filename,      ACCESS="STREAM", ACTION="WRITE")
        DO
           READ(10, IOSTAT=IERR)nml_char
           IF (IERR.NE.0) EXIT
           WRITE(11)nml_char
        END DO
        CLOSE(11)
        CLOSE(10)
        ! Copy compilation details to output directory
        filename=TRIM(filedir)//'compilation.txt'
        OPEN(UNIT=11, FILE=filename, FORM='formatted')
        write(11,'(" compile date = ",A)')TRIM(builddate)
        write(11,'(" patch dir    = ",A)')TRIM(patchdir)
        write(11,'(" remote repo  = ",A)')TRIM(gitrepo)
        write(11,'(" local branch = ",A)')TRIM(gitbranch)
        write(11,'(" last commit  = ",A)')TRIM(githash)
        CLOSE(11)
     endif
#ifndef WITHOUTMPI
     if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(myid==1.and.print_when_io) write(*,*)'End backup info etc.'

     if(myid==1.and.print_when_io) write(*,*)'Start backup amr'
     filename=TRIM(filedir)//'amr_'//TRIM(nchar)//'.out'
     call backup_amr(filename)
#ifndef WITHOUTMPI
     if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(myid==1.and.print_when_io) write(*,*)'End backup amr'

     if(hydro)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup hydro'
        filename=TRIM(filedir)//'hydro_'//TRIM(nchar)//'.out'
        filename_desc = trim(filedir)//'hydro_file_descriptor.txt'
        call backup_hydro(filename, filename_desc)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup hydro'
     end if

#ifdef RT
     if(rt.or.neq_chem)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup rt'
        filename=TRIM(filedir)//'rt_'//TRIM(nchar)//'.out'
        filename_desc = trim(filedir) // 'rt_file_descriptor.txt'
        call rt_backup_hydro(filename, filename_desc)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup rt'
     endif
#endif

     if(pic)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup part'
        filename=trim(filedir)//'part_'//trim(nchar)//'.out'
        filename_desc=TRIM(filedir)//'part_file_descriptor.txt'
        call backup_part(filename, filename_desc)
        if(sink)then
           filename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.out'
           call backup_sink(filename)
        end if
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup part'
     end if

     if(poisson)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup poisson'
        filename=TRIM(filedir)//'grav_'//TRIM(nchar)//'.out'
        call backup_poisson(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup poisson'
     end if
#ifdef ATON
     if(aton)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup rad'
        filename=TRIM(filedir)//'rad_'//TRIM(nchar)//'.out'
        call backup_radiation(filename)
        filename=TRIM(filedir)//'radgpu_'//TRIM(nchar)//'.out'
        call store_radiation(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup rad'
     end if
#endif
     if (gadget_output) then
        if(myid==1.and.print_when_io) write(*,*)'Start backup gadget format'
        filename=TRIM(filedir)//'gsnapshot_'//TRIM(nchar)
        call savegadget(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup gadget format'
     end if

     if(myid==1.and.print_when_io) write(*,*)'Start timer'
     ! Output timer: must be called by each process !
     filename=TRIM(filedir)//'timer_'//TRIM(nchar)//'.txt'
     call output_timer(.true., filename)
#ifndef WITHOUTMPI
     if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(myid==1.and.print_when_io) write(*,*)'End output timer'

  end if

end subroutine dump_all
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine dump_plotfile
   !-----------------------------------------------------------------------------------
   ! SUBROUTINE: dump_plotfile
   !
   ! DESCRIPTION:
   !   This subroutine is responsible for dumping necessary data to plot files. The routine
   !   addresses various conditions for dumping, including MPI parallelization, the current
   !   step in coarse-grained simulations, output directory structures, and specific output
   !   components like headers, info, patches, cooling data, and timers. Depending on the 
   !   configurations and condition checks, the subroutine writes to the relevant files. 
   !
   ! MODULE DEPENDENCIES:
   !   - amr_commons     : Common functionalities for Adaptive Mesh Refinement (AMR).
   !   - pm_commons      : Common functionalities related to the particle-mesh algorithm.
   !   - hydro_commons   : Common functionalities related to hydrodynamics.
   !   - cooling_module  : Module containing routines for cooling.
   !
   ! GLOBALS MODIFIED:
   !   - None directly by this subroutine; dependent subroutines called may modify globals.
   !
   ! INPUT/OUTPUT ARGUMENTS:
   !   - No explicit arguments; uses module variables and parameters.
   !
   ! NOTES:
   !   1. If the subroutine is compiled with MPI support (WITHOUTMPI not defined), MPI
   !      functionalities are invoked to ensure synchronization across processes.
   !   2. Depending on conditions, data is written to directories structured with 
   !      'plotfile_' naming conventions and various descriptive subdirectories.
   !   3. The subroutine checks the status of certain flags (e.g., hydro) to determine which
   !      modules' data needs to be written to the output.
   !   4. Verbose outputs can be triggered to notify the user about the progress of file 
   !      backups.
   !
   ! EXAMPLES OF USAGE:
   !   call dump_plotfile
   !-----------------------------------------------------------------------------------
   use amr_commons
   use pm_commons
   use hydro_commons
   use cooling_module
   implicit none
 #ifndef WITHOUTMPI
   include 'mpif.h'
 #endif
 #if ! defined (WITHOUTMPI) || defined (NOSYSTEM)
   integer::info
 #endif
   character::nml_char
   character(LEN=5)::nchar,ncharcpu
   character(LEN=80)::filename,filename_desc,filedir
   integer::ierr
 
   if(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0)return
   if(nstep_coarse==0.and.nrestart>0)return
   if(verbose)write(*,*)'Entering dump_plot_file'
 
   call write_screen
   call title(ifout,nchar)
   ifout=ifout+1
   if(t>=tout(iout).or.aexp>=aout(iout))iout=iout+1
   output_done=.true.
 
   if(IOGROUPSIZEREP>0)call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
 
   if(ndim>1)then
 
      if(IOGROUPSIZEREP>0) then
         filedir='plotfile_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
      else
         filedir='plotfile_'//TRIM(nchar)//'/'
      endif
 
      call create_output_dirs(filedir)
    
      if(myid==1.and.print_when_io) write(*,*)'Start backup header'
      ! Output header: must be called by each process !
      filename=TRIM(filedir)//'header_'//TRIM(nchar)//'.txt'
      call output_header(filename)
 #ifndef WITHOUTMPI
      if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
 #endif
      if(myid==1.and.print_when_io) write(*,*)'End backup header'
 
      if(myid==1.and.print_when_io) write(*,*)'Start backup info etc.'
      ! Only master process
      if(myid==1)then
         filename=TRIM(filedir)//'info_'//TRIM(nchar)//'.txt'
         call output_info(filename)
         filename=TRIM(filedir)//'makefile.txt'
         call output_makefile(filename)
         filename=TRIM(filedir)//'patches.txt'
         call output_patch(filename)
         if(cooling .and. .not. neq_chem)then
            filename=TRIM(filedir)//'cooling_'//TRIM(nchar)//'.out'
            call output_cool(filename)
         end if
         ! Copy namelist file to output directory
         filename=TRIM(filedir)//'namelist.txt'
         OPEN(10, FILE=namelist_file, ACCESS="STREAM", ACTION="READ")
         OPEN(11, FILE=filename,      ACCESS="STREAM", ACTION="WRITE")
         DO
            READ(10, IOSTAT=IERR)nml_char
            IF (IERR.NE.0) EXIT
            WRITE(11)nml_char
         END DO
         CLOSE(11)
         CLOSE(10)
         ! Copy compilation details to output directory
         filename=TRIM(filedir)//'compilation.txt'
         OPEN(UNIT=11, FILE=filename, FORM='formatted')
         write(11,'(" compile date = ",A)')TRIM(builddate)
         write(11,'(" patch dir    = ",A)')TRIM(patchdir)
         write(11,'(" remote repo  = ",A)')TRIM(gitrepo)
         write(11,'(" local branch = ",A)')TRIM(gitbranch)
         write(11,'(" last commit  = ",A)')TRIM(githash)
         CLOSE(11)
      endif
 #ifndef WITHOUTMPI
      if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
 #endif
      if(myid==1.and.print_when_io) write(*,*)'End backup info etc.'
 
      if(myid==1.and.print_when_io) write(*,*)'Start backup amr'
      filename=TRIM(filedir)//'amr_'//TRIM(nchar)//'.out'
      call backup_amr(filename)
 #ifndef WITHOUTMPI
      if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
 #endif
      if(myid==1.and.print_when_io) write(*,*)'End backup amr'
 
      if(hydro)then
         if(myid==1.and.print_when_io) write(*,*)'Start backup hydro'
         filename=TRIM(filedir)//'hydro_'//TRIM(nchar)//'.out'
         filename_desc = trim(filedir)//'hydro_file_descriptor.txt'
         call backup_hydro(filename, filename_desc)
 #ifndef WITHOUTMPI
         if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
 #endif
         if(myid==1.and.print_when_io) write(*,*)'End backup hydro'
      end if
 
      if(myid==1.and.print_when_io) write(*,*)'Start timer'
      ! Output timer: must be called by each process !
      filename=TRIM(filedir)//'timer_'//TRIM(nchar)//'.txt'
      call output_timer(.true., filename)
 #ifndef WITHOUTMPI
      if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
 #endif
      if(myid==1.and.print_when_io) write(*,*)'End output timer'
 
   end if
 
 end subroutine dump_plotfile
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
 
subroutine backup_amr(filename)
   !-----------------------------------------------------------------------------------
   ! SUBROUTINE: backup_amr
   !
   ! DESCRIPTION:
   !   This subroutine manages the backup process of the Adaptive Mesh Refinement (AMR)
   !   data. It converts units, calculates local constants, and writes various AMR grid
   !   data to a file, including time variables, level variables, boundaries, coarse level,
   !   and fine levels. MPI synchronization mechanisms are employed for parallel 
   !   execution (when compiled with MPI support).
   !
   ! MODULE DEPENDENCIES:
   !   - amr_commons    : Common AMR (Adaptive Mesh Refinement) variables and methods.
   !   - pm_commons     : Common PM (Particle Mesh) variables and methods.
   !   - hydro_commons  : Common hydrodynamic variables and methods.
   !
   ! GLOBALS MODIFIED:
   !   - Various AMR-related global variables can be modified by this subroutine.
   !
   ! INPUT/OUTPUT ARGUMENTS:
   !   - filename (character, input)  : The name of the file to which the AMR data
   !                                   will be backed up.
   !
   ! NOTES:
   !   1. If the subroutine is compiled with MPI support (WITHOUTMPI not defined), MPI
   !      functionalities are invoked to ensure synchronization across processes.
   !   2. The output file structure includes units conversion, local constants, and 
   !      various AMR grid data.
   !   3. The subroutine uses an efficient token-passing mechanism for IO in parallel 
   !      MPI environments to prevent simultaneous IO operations.
   !   4. Depending on configurations, the subroutine may interact with boundary 
   !      conditions and the overall domain decomposition strategy.
   !   5. Verbose outputs can be triggered to notify the user about entering the backup
   !      phase.
   !
   ! EXAMPLES OF USAGE:
   !   call backup_amr('my_amr_data.dat')
   !-----------------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::dummy_io,info2
#endif
  character(LEN=80)::filename

  integer::nx_loc,ny_loc,nz_loc,ilun
  integer::ilevel,ibound,ncache,istart,i,igrid,idim,ind,iskip
  integer,allocatable,dimension(:)::ind_grid,iig
  real(dp),allocatable,dimension(:)::xdp
  real(dp),dimension(1:3)::skip_loc
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale
  integer,parameter::tag=1120

  if(verbose)write(*,*)'Entering backup_amr'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  !-----------------------------------
  ! Output amr grid in file
  !-----------------------------------
  ilun=myid+10
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

   ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

  open(unit=ilun,file=fileloc,form='unformatted')
  ! Write grid variables
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)nx,ny,nz
  write(ilun)nlevelmax
  write(ilun)ngridmax
  write(ilun)nboundary
  write(ilun)ngrid_current
  write(ilun)boxlen
  ! Write time variables
  write(ilun)noutput,iout,ifout
  write(ilun)tout(1:noutput)
  write(ilun)aout(1:noutput)
  write(ilun)t
  write(ilun)dtold(1:nlevelmax)
  write(ilun)dtnew(1:nlevelmax)
  write(ilun)nstep,nstep_coarse
  write(ilun)einit,mass_tot_0,rho_tot
  write(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
  write(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
  write(ilun)mass_sph
  ! Write levels variables
  write(ilun)headl(1:ncpu,1:nlevelmax)
  write(ilun)taill(1:ncpu,1:nlevelmax)
  write(ilun)numbl(1:ncpu,1:nlevelmax)
  write(ilun)numbtot(1:10,1:nlevelmax)
  ! Read boundary linked list
  if(simple_boundary)then
     write(ilun)headb(1:nboundary,1:nlevelmax)
     write(ilun)tailb(1:nboundary,1:nlevelmax)
     write(ilun)numbb(1:nboundary,1:nlevelmax)
  end if
  ! Write free memory
  write(ilun)headf,tailf,numbf,used_mem,used_mem_tot
  ! Write cpu boundaries
  write(ilun)ordering
  if(ordering=='bisection') then
     write(ilun)bisec_wall(1:nbinodes)
     write(ilun)bisec_next(1:nbinodes,1:2)
     write(ilun)bisec_indx(1:nbinodes)
     write(ilun)bisec_cpubox_min(1:ncpu,1:ndim)
     write(ilun)bisec_cpubox_max(1:ncpu,1:ndim)
  else
     write(ilun)bound_key(0:ndomain)
  endif

  ! Write coarse level
  write(ilun)son(1:ncoarse)
  write(ilun)flag1(1:ncoarse)
  write(ilun)cpu_map(1:ncoarse)
  ! Write fine levels
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache),iig(1:ncache))
           ! Write grid index
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           write(ilun)ind_grid
           ! Write next index
           do i=1,ncache
              iig(i)=next(ind_grid(i))
           end do
           write(ilun)iig
           ! Write prev index
           do i=1,ncache
              iig(i)=prev(ind_grid(i))
           end do
           write(ilun)iig
           ! Write grid center
           do idim=1,ndim
              do i=1,ncache
                 xdp(i)=xg(ind_grid(i),idim)
              end do
              write(ilun)xdp
           end do
           ! Write father index
           do i=1,ncache
              iig(i)=father(ind_grid(i))
           end do
           write(ilun)iig
           ! Write nbor index
           do ind=1,twondim
              do i=1,ncache
                 iig(i)=nbor(ind_grid(i),ind)
              end do
              write(ilun)iig
           end do
           ! Write son index
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=son(ind_grid(i)+iskip)
              end do
              write(ilun)iig
           end do
           ! Write cpu map
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=cpu_map(ind_grid(i)+iskip)
              end do
              write(ilun)iig
           end do
           ! Write refinement map
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=flag1(ind_grid(i)+iskip)
              end do
              write(ilun)iig
           end do
           deallocate(xdp,iig,ind_grid)
        end if
     end do
  end do
  close(ilun)

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

end subroutine backup_amr
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_info(filename)
   !==============================================================================
   ! Subroutine output_info
   !
   ! Description:
   !   This subroutine is designed to write essential information about the simulation
   !   setup to a file, such as the run parameters, physical parameters, and
   !   the domain decomposition ordering.
   !
   ! Parameters:
   !   - filename: The name of the file to which the information will be written.
   !
   ! Modules Used:
   !   - amr_commons: Common data structures and constants related to adaptive mesh refinement.
   !   - hydro_commons: Common data structures and constants for hydrodynamics calculations.
   !   - pm_commons: Common data structures and constants for particle mesh calculations.
   !   - mpif.h: MPI (Message Passing Interface) library header. Used for parallel computing.
   !
   ! Local Variables:
   !   - nx_loc, ny_loc, nz_loc: Local constants determining the resolution along each axis.
   !   - ilun: Unit number associated with the file operation.
   !   - icpu, idom: Iterative counters.
   !   - ierr: Error status from file operations.
   !   - scale, scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v: Scaling factors.
   !   - fileloc: File location obtained by trimming filename.
   !
   ! Notes:
   !   - The subroutine starts by calculating conversion factors from user units to cgs units.
   !   - It then writes a number of parameters and physical constants to the file for further reference.
   !   - If there's an error during file operations, the subroutine can abort the whole program.
   !   - MPI is conditionally used, if it is enabled (controlled by the WITHOUTMPI flag).
   !
   ! Example usage:
   !   call output_info('simulation_info.dat')
   !
   !==============================================================================
   subroutine output_info(filename)
      ...
      end subroutine output_info
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom,ierr
  real(dp)::scale
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_info'

  ilun=11

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  scale=boxlen/dble(nx_loc)

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted',iostat=ierr)
  if(ierr .ne. 0)then
     write(*,*) 'Error - Could not write ',fileloc
#ifndef WITHOUTMPI
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#else
     stop
#endif
  endif

  ! Write run parameters
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("ndim        =",I11)')ndim
  write(ilun,'("levelmin    =",I11)')levelmin
  write(ilun,'("levelmax    =",I11)')nlevelmax
  write(ilun,'("ngridmax    =",I11)')ngridmax
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("boxlen      =",E23.15)')scale
  write(ilun,'("time        =",E23.15)')t
  write(ilun,'("aexp        =",E23.15)')aexp
  write(ilun,'("H0          =",E23.15)')h0
  write(ilun,'("omega_m     =",E23.15)')omega_m
  write(ilun,'("omega_l     =",E23.15)')omega_l
  write(ilun,'("omega_k     =",E23.15)')omega_k
  write(ilun,'("omega_b     =",E23.15)')omega_b
  write(ilun,'("unit_l      =",E23.15)')scale_l
  write(ilun,'("unit_d      =",E23.15)')scale_d
  write(ilun,'("unit_t      =",E23.15)')scale_t
  write(ilun,*)

  ! Write ordering information
  write(ilun,'("ordering type=",A80)')ordering
  if(ordering=='bisection') then
     do icpu=1,ncpu
        ! write 2*ndim floats for cpu bound box
        write(ilun,'(E23.15)')bisec_cpubox_min(icpu,:),bisec_cpubox_max(icpu,:)
        ! write 1 float for cpu load
        write(ilun,'(E23.15)')dble(bisec_cpu_load(icpu))
     end do
  else
     write(ilun,'("   DOMAIN   ind_min                 ind_max")')
     do idom=1,ndomain
        write(ilun,'(I8,1X,E23.15,1X,E23.15)')idom,bound_key(idom-1),bound_key(idom)
     end do
  endif

  close(ilun)

end subroutine output_info
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_header(filename)
   !==============================================================================
   ! Subroutine output_header
   !
   ! Description:
   !   This subroutine writes a header to the given file containing information 
   !   about the particle families and other particle-related properties in the 
   !   simulation. The header includes a count of particles for each family type, 
   !   as well as a listing of particle fields that are present.
   !
   ! Parameters:
   !   - filename: The name of the file to which the header will be written.
   !
   ! Modules Used:
   !   - amr_commons: Common data structures and constants related to adaptive mesh refinement.
   !   - hydro_commons: Common data structures and constants for hydrodynamics calculations.
   !   - pm_commons: Common data structures and constants for particle mesh calculations.
   !   - mpif.h: MPI (Message Passing Interface) library header. Used for parallel computing.
   !
   ! Local Variables:
   !   - ilun: Unit number associated with the file operation.
   !   - fileloc: File location obtained by trimming filename.
   !   - npart_family_loc: Local count of particles in each family.
   !   - npart_family: Global count of particles in each family (across all processes).
   !   - npart_all_loc: Local count of all particles.
   !   - npart_all: Global count of all particles (across all processes).
   !   - ifam, ipart: Iterative counters.
   !
   ! Notes:
   !   - The subroutine calculates both local (to the current process) and global (across all processes)
   !     particle counts for each family type.
   !   - MPI is conditionally used to aggregate particle counts across processes, if it is enabled 
   !     (controlled by the WITHOUTMPI flag).
   !   - Only process with ID "1" writes to the output file, ensuring a single header is written.
   !
   ! Example usage:
   !   call output_header('particles_header.dat')
   !
   !==============================================================================
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
  integer::info
#endif
  character(LEN=80)::filename

  integer::ilun
  character(LEN=80)::fileloc
#ifdef LONGINT
  integer(i8b)::npart_family_loc(-5:5), npart_family(-5:5), npart_all_loc, npart_all
#else
  integer::npart_family_loc(-5:5), npart_family(-5:5), npart_all_loc, npart_all
#endif
  integer :: ifam, ipart

  if(verbose)write(*,*)'Entering output_header'
  if(myid==1)then
     ! Open file
     fileloc=TRIM(filename)
     open(newunit=ilun,file=fileloc,form='formatted')
  end if

  ! Compute total number of particles
  ! Count number of particles
  npart_family_loc = 0; npart_all_loc = 0
  do ipart = 1, npartmax
     ! Only used particles have a levelp > 0
     if (levelp(ipart) > 0) then
        npart_all_loc = npart_all_loc + 1
        do ifam = -5, 5
           if (typep(ipart)%family == ifam) then
              npart_family_loc(ifam) = npart_family_loc(ifam) + 1
           end if
        end do
     end if
  end do

#ifndef WITHOUTMPI
#ifdef LONGINT
  call MPI_ALLREDUCE(npart_family_loc,npart_family,11,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(npart_all_loc,npart_all,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#else
  call MPI_ALLREDUCE(npart_family_loc,npart_family,11,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(npart_all_loc,npart_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#else
  npart_family = npart_family_loc
  npart_all = npart_all_loc
#endif

  if (myid == 1) then
     write(ilun, '(a1,a12,a10)') '#', 'Family', 'Count'
     do ifam = -NFAMILIES, NFAMILIES
        write(ilun, '(a13, i10)') &
             trim(particle_family_keys(ifam)), npart_family(ifam)
     end do
     write(ilun, '(a13, i10)') &
          'undefined', npart_all - sum(npart_family)
  end if

  if (myid == 1) then
     ! Keep track of what particle fields are present
     write(ilun,*)'Particle fields'
     write(ilun,'(a)',advance='no')'pos vel mass iord level family tag '
#ifdef OUTPUT_PARTICLE_POTENTIAL
     write(ilun,'(a)',advance='no')'phi '
#endif
     if(star.or.sink) then
        write(ilun,'(a)',advance='no')'tform '
        if(metal) then
           write(ilun,'(a)',advance='no')'metal '
        endif
     endif
     close(ilun)

  endif

end subroutine output_header
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine savegadget(filename)
   !==============================================================================
   ! Subroutine savegadget
   !
   ! Description:
   !   This subroutine saves particle data in the Gadget file format. The Gadget 
   !   file format is commonly used in cosmological simulations. This subroutine 
   !   prepares the header, as well as the positions, velocities, and IDs of the 
   !   particles, and then invokes the gadgetwritefile subroutine to write 
   !   these to the specified file.
   !
   ! Parameters:
   !   - filename: The name of the file to which the particle data should be written.
   !
   ! Modules Used:
   !   - amr_commons: Common data structures and constants related to adaptive mesh refinement.
   !   - hydro_commons: Common data structures and constants for hydrodynamics.
   !   - pm_commons: Common data structures and constants related to particle mesh operations.
   !   - gadgetreadfilemod: Module to handle reading and writing in Gadget file format.
   !   - mpif.h: MPI (Message Passing Interface) library header. Used for parallel computing.
   !
   ! Local Variables:
   !   - header: Gadget file header containing metadata for the particles.
   !   - pos: Array containing particle positions.
   !   - vel: Array containing particle velocities.
   !   - ids: Array containing particle IDs.
   !   - gadgetvfact: Factor used in calculating the particle velocities.
   !   - npart_tot: Total number of particles across all processes.
   !   - RHOcrit: Critical density of the universe.
   !   ... (other intermediate and utility variables)
   !
   ! Notes:
   !   - Particles in this routine are described using their positions (xp), 
   !     velocities (vp), and IDs (idp). The npart and npartmax variables seem 
   !     to represent the current number and the maximum number of particles, respectively.
   !   - Depending on the LONGINT flag, 64-bit integers might be used for some operations.
   !   - The subroutine handles parallel operations using MPI, but can also work 
   !     in a serial mode when the WITHOUTMPI flag is defined.
   !   - The "RHOcrit" parameter is the critical density of the universe and is used 
   !     to compute the mass of the particles for the Gadget file.
   !
   ! Example usage:
   !   call savegadget('/path/to/output/file')
   !
   !==============================================================================
  use amr_commons
  use hydro_commons
  use pm_commons
  use gadgetreadfilemod
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  integer(i8b)::npart_loc
#endif
  character(LEN=80)::filename
  TYPE(gadgetheadertype)::header
  real,allocatable,dimension(:,:)::pos,vel
  integer(i8b),allocatable,dimension(:)::ids
  integer::i,idim,ipart
  real(dp)::gadgetvfact
  integer(i8b)::npart_tot
  real(dp),parameter::RHOcrit=2.7755d11

#ifndef WITHOUTMPI
  npart_loc=npart
#ifndef LONGINT
  call MPI_ALLREDUCE(npart_loc,npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  call MPI_ALLREDUCE(npart_loc,npart_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#else
  npart_tot=npart
#endif

  allocate(pos(ndim, npart), vel(ndim, npart), ids(npart))
  gadgetvfact = 100.0 * boxlen_ini / aexp / SQRT(aexp)

  header%npart = 0
  header%npart(2) = npart
  header%mass = 0
  header%mass(2) = omega_m*RHOcrit*(boxlen_ini)**3/npart_tot/1.d10
  header%time = aexp
  header%redshift = 1.d0/aexp-1.d0
  header%flag_sfr = 0
  header%nparttotal = 0
#ifndef LONGINT
  header%nparttotal(2) = npart_tot
#else
  header%nparttotal(2) = MOD(npart_tot,4294967296_i8b)
#endif
  header%flag_cooling = 0
  header%numfiles = ncpu
  header%boxsize = boxlen_ini
  header%omega0 = omega_m
  header%omegalambda = omega_l
  header%hubbleparam = h0/100.0
  header%flag_stellarage = 0
  header%flag_metals = 0
  header%totalhighword = 0
#ifndef LONGINT
  header%totalhighword(2) = 0
#else
  header%totalhighword(2) = npart_tot/4294967296_i8b
#endif
  header%flag_entropy_instead_u = 0
  header%flag_doubleprecision = 0
  header%flag_ic_info = 0
  header%lpt_scalingfactor = 0
  header%unused = ' '

  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           if (ipart .gt. npart) then
                write(*,*) myid, "Ipart=",ipart, "exceeds", npart
                call clean_stop
           endif
           pos(idim, ipart)=real(xp(i,idim) * boxlen_ini , kind=4)
           vel(idim, ipart)=real(vp(i,idim) * gadgetvfact , kind=4)
           if (idim.eq.1) ids(ipart) = idp(i)
        end if
     end do
  end do

  call gadgetwritefile(filename, myid-1, header, pos, vel, ids)
  deallocate(pos, vel, ids)

end subroutine savegadget
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine create_output_dirs(filedir)
   !==============================================================================
   ! Subroutine create_output_dirs
   !
   ! Description:
   !   This subroutine creates the specified output directory and, if necessary,
   !   its parent directories. The directory creation can be achieved either 
   !   using a system command or platform-specific Fortran calls depending on 
   !   the build flags (e.g., NOSYSTEM, WITHOUTMPI).
   !
   ! Parameters:
   !   - filedir: The directory path that needs to be created.
   !
   ! Modules Used:
   !   - amr_commons: Common data structures and constants related to adaptive mesh refinement.
   !   - mpif.h: MPI (Message Passing Interface) library header. Used for parallel computing.
   !
   ! Local Variables:
   !   - ierr: Error code returned from system operations.
   !   - filecmd: System command to create the directory.
   !   - filedirini: Initial directory segment for NOSYSTEM use.
   !   - info: Error information for MPI calls.
   !
   ! Notes:
   !   - If the NOSYSTEM flag is defined at compile time, the subroutine uses platform-specific
   !     Fortran calls (e.g., PXFMKDIR) to make directories.
   !   - If the WITHOUTMPI flag isn't defined, MPI calls are used to synchronize operations 
   !     across different processes.
   !   - The "withoutmkdir" flag controls whether the directory creation should be skipped.
   !   - If the system command to create directories fails, an error message is displayed 
   !     and the program can be terminated (conditionally based on the WITHOUTMPI flag).
   !
   ! Example usage:
   !   call create_output_dirs('/path/to/output/dir/')
   !
   !==============================================================================
  use amr_commons
  implicit none
  character(LEN=80), intent(in):: filedir
  integer :: ierr
  character(LEN=80)::filecmd
#ifdef NOSYSTEM
  character(LEN=80)::filedirini
#endif
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer :: info
#endif


  filecmd='mkdir -p '//TRIM(filedir)

  if (.not.withoutmkdir) then
#ifdef NOSYSTEM
    filedirini = filedir(1:13)
    call PXFMKDIR(TRIM(filedirini),LEN(TRIM(filedirini)),O'755',info)
    call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'755',info)
#else
    call EXECUTE_COMMAND_LINE(filecmd,exitstat=ierr,wait=.true.)
    if(ierr.ne.0 .and. ierr.ne.127)then
      write(*,*) 'Error - Could not create ',trim(filedir),' error code=',ierr
#ifndef WITHOUTMPI
      call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
      stop
#endif
    endif
#endif
  endif


#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif


end subroutine create_output_dirs
