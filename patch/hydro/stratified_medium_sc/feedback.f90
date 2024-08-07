!################################################################
!################################################################
!################################################################
!################################################################
subroutine thermal_feedback(ilevel,icount)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel,icount
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::t0,scale,dx_min,vsn,rdebris,ethermal
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  real(dp),dimension(1:3)::skip_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  if(icount==2)return

  ! Gather star particles only.

#if NDIM==3
  ! Loop over cpus
  do icpu=1,ncpu
     
     !  headl is linked list
     ! head list in level

     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids 
     !
     !  in level (level is passed as argument to function)
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid) !  linked list of pointers to particles in linked list
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).gt.0.and.tp(ipart).ne.0)then !  idp is particle ID; tp is birth epoch
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather star particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                !  feedbk routine called for a number of particles equal to
                !       nvector.
                !
                ! nvector is size of vector sweeps
 
                 call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

#endif

111 format('   Entering thermal_feedback for level ',I2)

end subroutine thermal_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each debris particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! uold.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(kind=8)::RandNum
  real(dp)::SN_BOOST,mstar,dx_min,vol_min
  real(dp)::xxx,mmm,t0,ESN,mejecta,zloss
  real(dp)::dx,dx_loc,scale,vol_loc,birth_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !  function is called for a specific level. It is called for each level,
  ! recurviely in amr_step
 
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif

  ! Compute stochastic boost to account for target GMC mass
  SN_BOOST=MAX(mass_gmc*2d33/(scale_d*scale_l**3)/mstar,1d0)

  ! Massive star lifetime from Myr to code units
  t0=10.*1d6*(365.*24.*3600.)/scale_t

  ! Type II supernova specific energy from cgs to code units
  ESN=1d51/(10.*2d33)/scale_v**2

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=2.0D0.or.x(j,idim)>=4.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in sn2'
     write(*,*)ilevel,ng,np
  end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Compute individual time steps
  ! WARNING: the time step is always the coarser level time step
  ! since we do not have feedback for icount=2
  if(ilevel==levelmin)then
     do j=1,np
        if(ok(j))then
           dteff(j)=dtold(levelp(ind_part(j)))
        end if
     end do
  else
     do j=1,np
        if(ok(j))then
           dteff(j)=dtold(levelp(ind_part(j))-1)
        end if
     end do
  endif

  ! Reset ejected mass, metallicity, thermal energy
  do j=1,np
     if(ok(j))then
        mloss(j)=0d0
        mzloss(j)=0d0
        ethermal(j)=0d0
     endif
  end do

  ! Compute stellar mass loss and thermal feedback due to supernovae
  !
  !  f_w is SN mass loading factor
  if(f_w==0)then
     do j=1,np
        
       if(ok(j))then
           birth_time=tp(ind_part(j))
           ! Make sure that we don't count feedback twice
           if(birth_time.lt.(t-t0).and.birth_time.ge.(t-t0-dteff(j)))then
              ! Stellar mass loss
              mejecta=eta_sn*mp(ind_part(j))
              mloss(j)=mloss(j)+mejecta/vol_loc
              ! Thermal energy
              ethermal(j)=ethermal(j)+mejecta*ESN/vol_loc
              ! Metallicity
              if(metal)then
                 zloss=yield+(1d0-yield)*zp(ind_part(j))
                 mzloss(j)=mzloss(j)+mejecta*zloss/vol_loc
              endif
              ! Reduce star particle mass
              mp(ind_part(j))=mp(ind_part(j))-mejecta
              ! Boost SNII energy and depopulate accordingly
              if(SN_BOOST>1d0)then
                 call ranf(localseed,RandNum)
                 if(RandNum<1d0/SN_BOOST)then
                    mloss(j)=SN_BOOST*mloss(j)
                    mzloss(j)=SN_BOOST*mzloss(j)
                    ethermal(j)=SN_BOOST*ethermal(j)
                 else
                    mloss(j)=0d0
                    mzloss(j)=0d0
                    ethermal(j)=0d0
                 endif
              endif
           endif
        end if
     end do
  endif

  ! Update hydro variables due to feedback
  
  do j=1,np
     if(ok(j))then
        ! Specific kinetic energy of the star
        ekinetic(j)=0.5*(vp(ind_part(j),1)**2 &
             &          +vp(ind_part(j),2)**2 &
             &          +vp(ind_part(j),3)**2)
        ! Update hydro variable in NGP cell

        uold(indp(j),1)=uold(indp(j),1)+mloss(j) !  mass
        uold(indp(j),2)=uold(indp(j),2)+mloss(j)*vp(ind_part(j),1) !  momentum
        uold(indp(j),3)=uold(indp(j),3)+mloss(j)*vp(ind_part(j),2)
        uold(indp(j),4)=uold(indp(j),4)+mloss(j)*vp(ind_part(j),3)
        uold(indp(j),5)=uold(indp(j),5)+mloss(j)*ekinetic(j)+ethermal(j) !  energy
        
     endif
  end do

  ! Add metals
  if(metal)then
     do j=1,np
        if(ok(j))then
           uold(indp(j),imetal)=uold(indp(j),imetal)+mzloss(j)
        endif
     end do
  endif

  ! Add delayed cooling switch variable
  if(delayed_cooling)then
     do j=1,np
        if(ok(j))then
           uold(indp(j),idelay)=uold(indp(j),idelay)+mloss(j)
        endif
     end do
  endif

#endif
  
end subroutine feedbk
!################################################################
!################################################################
!################################################################
!################################################################
subroutine blast_wave_feedback(ilevel,icount)
  use pm_commons
  use amr_commons
  use hydro_commons
  use random
  use poisson_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
!  this routine is used in this patch to inject thermal
! and kinetic energy when a SNR explodes.
  integer::ind,ix,iy,iz,ngrid,iskip,ilev
  integer::i,j,nx_loc,igrid,ncache,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,drr,dtprob
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  logical ,dimension(1:nvector)::ok
  real(dp)::RandNum
  integer::isnr,nrandom_local,n_expl
  integer,dimension(:),allocatable::num_snr_local
  real(dp),dimension(:),allocatable::xsnr,ysnr,zsnr
  real(dp),dimension(:),allocatable::rho_average,rho_average_all,weight,weight_all
  real(dp)::PoissMean,t_s,boost,dcgs,ddyn,vri,vxi,vyi,vzi,mu,a1,a2,z0
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  real(dp)::Ei,Vi,Pi_,M_ej,rho_ej,pi
  real(dp),dimension(1:nvector,1:3)::pos
  integer,dimension(1:nvector)::cell_index,cell_levl,cc

  real(dp)::alpha,rcool,rrise,rnod,rbreak

  integer::exp_type, funit, istat
  character(80)::filename

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  if(icount==2)return

!! Open the file to write the SNR data to
  if(myid==1)then
   filename = 'explosions.txt'
   open(funit, file=trim(filename), position='APPEND', status = 'OLD', iostat=istat)

   if(istat .NE.0)then !File didn't exsits so create it
      open(funit, file=trim(filename), position='APPEND', status = 'NEW')
   endif
  endif

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  pi=3.14159265359

  ! Maximum radius of the ejecta
  rmax=rc/aexp !  Check the units of the ejecta radius. Not sure about aexp, no problem for cosmo=.false..
  rmax2=rmax*rmax

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  ! Check if it is the right time for a SN to explode.
  nsnr=int(sn_rate*tout(noutput)) ! Total number of SNR in the evolution of the sim 
  if(.not.init_marker)then
     allocate(marker(1:nsnr))
     marker=0
     do i=1,nblast_current
        marker(i)=1
     end do
     if(verbose .and. myid==1) write(*,*)'nblast_current =',nblast_current
     if(verbose .and. myid==1) write(*,*)'nblast_current_internal (before update) =',nblast_current_internal
     nblast_current_internal=nblast_current
     if(verbose .and. myid==1) write(*,*)'nblast_current_internal (after update) =',nblast_current_internal
     init_marker=.true.
  end if

  ! Loop to identify the SN centers
  nrandom_local=nrandom/ncpu
  allocate(xsnr(1:nrandom))
  allocate(ysnr(1:nrandom))
  allocate(zsnr(1:nrandom))
  allocate(rho_average(1:nrandom))
  allocate(rho_average_all(1:nrandom))
  allocate(weight(1:nrandom))
  allocate(weight_all(1:nrandom))
  allocate(num_snr_local(1:nrandom))
  xsnr=0.0d0
  ysnr=0.0d0
  zsnr=0.0d0
  num_snr_local=0
  isnr=1+(myid-1)*nrandom_local
  ! Loop over levels
  do ilev=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities)                                                                                
     dx=0.5D0**ilev
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position                                                                                           
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do
  
     ! Loop over grids                                                                                                                                  
     ncache=active(ilev)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilev)%igrid(igrid+i-1)
        end do
        
        ! Loop over cells                                                                                                                               
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells                                                                                                                            
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position                                                                                                                
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
             
                 ! Compute expected no. of SN per cell
                 boost=100.
                 dcgs=uold(ind_cell(i),1)*scale_d
                 ! Density for star formation timescale
                 ! Parameters of the model - convert to cgs
                 a1=gravity_params(1) ! disc component factor
                 a2=gravity_params(2) ! halo component factor
                 z0=gravity_params(3) ! scale height of the disc
                 a1=a1*3.08d21/(1d6*365.*24.*3600.)**2 ! cm/s2
                 a2=a2/(1d6*365.*24.*3600.)**2 ! s-2
                 z0=z0*3.08d21 ! cm
                 ddyn= 1.0/(4.0*3.14159*6.66d-8)*(a2+a1/sqrt((z*scale_l)**2+z0**2)-a1*(z*scale_l)**2/((z*scale_l)**2+z0**2)**1.5)
                 ! Star formation time
                 t_s=sqrt(3.0*3.14159/16.0/6.66d-8/ddyn)/scale_t*boost               
                 !  the star formation time is set by the external potential generating the force.
                 PoissMean=dtnew(ilevel)/t_s*dcgs*vol_loc*(scale_l**3)/2.0d35
                 ! Compute Poisson realisation
                 call poissdev(localseed,PoissMean,n_expl)
                 num_snr_local(isnr)=n_expl
                 !if(num_snr_local(isnr)>0)write(*,*)PoissMean,num_snr_local(isnr)
                 if(num_snr_local(isnr)>0)then
                    xsnr(isnr)=x+dx_loc/2
                    ysnr(isnr)=y+dx_loc/2
                    zsnr(isnr)=z+dx_loc/2
                    isnr=isnr+1
                 end if
              endif
           end do
        end do
        ! End loop over cells                                                                                                                           
     end do
     ! End loop over grids                                                                                                                              
  end do
  ! End loop over levels    
  
  ! Share the locations of the explosions with the other cores
#ifndef WITHOUTMPI
  x_expl=0.0d0
  y_expl=0.0d0
  z_expl=0.0d0
  random_expl=0
  call MPI_ALLREDUCE(xsnr,x_expl,nrandom,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)  
  call MPI_ALLREDUCE(ysnr,y_expl,nrandom,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(zsnr,z_expl,nrandom,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(num_snr_local,random_expl,nrandom,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  do i=1,nrandom
     x_expl(i)=xsnr(i)
     y_expl(i)=ysnr(i)
     z_expl(i)=zsnr(i)
     random_expl(i)=num_snr_local(i)
  end do
#endif

  n_expl=0
  do i=1,nrandom
     if(random_expl(i)>0)then
        n_expl=n_expl+1
     end if
  end do

  if (n_expl > 0) then
     ! Compute average density
     rho_average=0.0d0
     weight=1.0d-20
     ! Loop over levels
     do ilev=levelmin,nlevelmax
        ! Computing local volume (important for averaging hydro quantities)
        dx=0.5D0**ilev
        dx_loc=dx*scale
        vol_loc=dx_loc**ndim
        ! Cells center position relative to grid center position
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do

        ! Loop over grids
        ncache=active(ilev)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(ilev)%igrid(igrid+i-1)
           end do
           
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              
              
              ! Flag leaf cells
              do i=1,ngrid
                 ok(i)=son(ind_cell(i))==0
              end do
              
              do i=1,ngrid
                 if(ok(i))then
                    ! Get gas cell position
                    x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                    y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                    z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                    
                    ! Loop over SNR
                    do isnr=1,nrandom
                       ! Check if the cell lies within the ejecta radius
                       dxx=x-x_expl(isnr)
                       if(dxx> 0.5*scale)then
                          dxx=dxx-scale
                       endif
                       if(dxx<-0.5*scale)then
                          dxx=dxx+scale
                       endif
                       dyy=y-y_expl(isnr)
                       if(dyy> 0.5*scale)then
                          dyy=dyy-scale
                       endif
                       if(dyy<-0.5*scale)then
                          dyy=dyy+scale
                       endif
                       dzz=z-z_expl(isnr)
                       if(dzz> 0.5*scale)then
                          dzz=dzz-scale
                       endif
                       if(dzz<-0.5*scale)then
                          dzz=dzz+scale
                       endif
                       drr=dxx*dxx+dyy*dyy+dzz*dzz
                       
                       if(drr.lt.rmax2.and.random_expl(isnr)>0)then
                          nblast_current_internal=nblast_current_internal+1
                          marker(nblast_current_internal)=1
                          rho_average(isnr)=rho_average(isnr)+uold(ind_cell(i),1)*vol_loc
                          weight(isnr)=weight(isnr)+vol_loc
                       endif
                    end do
                    ! End loop over SNR
                 end if
              end do
           end do
           ! End loop over cells
        end do
        ! End loop over grids
     end do
     ! End loop over levels

#ifndef WITHOUTMPI
     rho_average_all=0.0d0
     weight_all=0.0d0
     call MPI_ALLREDUCE(rho_average,rho_average_all,nrandom,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(weight,weight_all,nrandom,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     do i=1,nrandom
        rho_average(i)=rho_average_all(i)
        weight(i)=weight_all(i)
        rho_average(i)=rho_average(i)/weight(i)
        if(random_expl(i)>0.and.myid==1)then
           rbreak=4.001*3.08d18/scale_l*(rho_average(i)/scale_nH/100.)**(-0.429)
           rnod=0.969*3.08d18/scale_l*(rho_average(i)/scale_nH/100.)**(-0.330)
           !rbreak=7.97*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.458)                                                                                      
           !rnod=2.438*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.349)                                                                                       
           !Momentum                                                                                                                                           
           if (rc<rbreak) then
              P_SN_rad=3.479830d42*(rc/rnod)**1.5
           else
              P_SN_rad=3.479830d42*(rbreak/rnod)**1.5
           end if
           !Thermal energy                                                                                                                                     
           alpha=-7.83*(rho_average(i)/scale_nH/100.)**(0.024)
           rcool=3.021*3.08d18/scale_l*(rho_average(i)/scale_nH/100.)**(-0.421)
           rrise=5.472*3.08d18/scale_l*(rho_average(i)/scale_nH/100.)**(-0.403)
           !alpha=-11.3*(rhocool/scale_nH/100.)**(0.072)                                                                                                       
           !rcool=6.35*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.418)                                                                                       
           !rrise=9.19*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.439)                                                                                       
           if (rc<rcool) then
              E_SN_th=7.1d50
           else
              if (rcool <= rc .and. rc < rrise) then
                 E_SN_th=7.1d50*(rc/rcool)**alpha
              else
                 E_SN_th=7.1d50*(rrise/rcool)**alpha
              end if
           end if
           if (myid == 1) write(funit,'(A, I7.5)') 'Output number is', ifout
           if (myid == 1) write(funit,'(A, I10.1)') 'Explosion number:', nblast_current_internal
           if (myid == 1) write(funit, '(A, ES13.5, ES13.5, ES13.5)') 'Position is:',x_expl(i), y_expl(i), z_expl(i)
           if(myid==1)write(funit,'(A, ES12.5)')'Average density ', rho_average(i)
           if(myid==1)write(funit,'(A, ES12.5, ES12.5)')'E_th, P_rad (cgs) =', E_SN_th, P_SN_rad

           nblast_current_internal = nblast_current_internal + 1

         end if
     end do
#else
     do i=1,nrandom
        rho_average(i)=rho_average(i)/weight(i)
        if(random_expl(i)>0.and.myid==1)then
            !update the SNR counter 
            nblast_current_internal = nblast_current_internal + 1
           rbreak=4.001*3.08d18/scale_l*(rho_average(i)/scale_nH/100.)**(-0.429)
           rnod=0.969*3.08d18/scale_l*(rho_average(i)/scale_nH/100.)**(-0.330)
           !rbreak=7.97*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.458)
           !rnod=2.438*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.349)
           !Momentum
           if (rc<rbreak) then
              P_SN_rad=3.479830d42*(rc/rnod)**1.5
           else
              P_SN_rad=3.479830d42*(rbreak/rnod)**1.5
           end if
           !Thermal energy
           alpha=-7.83*(rho_average(i)/scale_nH/100.)**(0.024)
           rcool=3.021*3.08d18/scale_l*(rho_average(i)/scale_nH/100.)**(-0.421)
           rrise=5.472*3.08d18/scale_l*(rho_average(i)/scale_nH/100.)**(-0.403)
           !alpha=-11.3*(rhocool/scale_nH/100.)**(0.072)
           !rcool=6.35*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.418)
           !rrise=9.19*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.439)
           if (rc<rcool) then
              E_SN_th=7.1d50
           else 
              if (rcool <= rc .and. rc < rrise) then
                 E_SN_th=7.1d50*(rc/rcool)**alpha
              else
                 E_SN_th=7.1d50*(rrise/rcool)**alpha
              end if
           end if
        end if
     end do
#endif

     ! Inject energy/momentum with subgrid model
     ! Loop over levels
     do ilev=levelmin,nlevelmax
        ! Computing local volume (important for averaging hydro quantities)
        dx=0.5D0**ilev
        dx_loc=dx*scale
        vol_loc=dx_loc**ndim
        ! Cells center position relative to grid center position
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do
        
        ! Loop over grids
        ncache=active(ilev)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=active(ilev)%igrid(igrid+i-1)
           end do
           
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              
              ! Flag leaf cells
              do i=1,ngrid
                 ok(i)=son(ind_cell(i))==0
              end do
              
              do i=1,ngrid
                 if(ok(i))then
                    ! Get gas cell position
                    x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                    y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                    z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                    
                    ! Loop over SNR
                    do isnr=1,nrandom
                       ! Check if the cell lies within the ejecta radius
                       dxx=x-x_expl(isnr)
                       if(dxx> 0.5*scale)then
                          dxx=dxx-scale
                       endif
                       if(dxx<-0.5*scale)then
                          dxx=dxx+scale
                       endif
                       dyy=y-y_expl(isnr)
                       if(dyy> 0.5*scale)then
                          dyy=dyy-scale
                       endif
                       if(dyy<-0.5*scale)then
                          dyy=dyy+scale
                       endif
                       dzz=z-z_expl(isnr)
                       if(dzz> 0.5*scale)then
                          dzz=dzz-scale
                       endif
                       if(dzz<-0.5*scale)then
                          dzz=dzz+scale
                       endif
                       drr=dxx*dxx+dyy*dyy+dzz*dzz
                       
                       if(drr.lt.rmax2.and.random_expl(isnr)>0)then                 
                          rbreak=4.001*3.08d18/scale_l*(rho_average(isnr)/scale_nH/100.)**(-0.429)
                          rnod=0.969*3.08d18/scale_l*(rho_average(isnr)/scale_nH/100.)**(-0.330)
                          !rbreak=7.97*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.458)
                          !rnod=2.438*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.349)
                          !Momentum
                          if (rc<rbreak) then
                             P_SN_rad=3.479830d42*(rc/rnod)**1.5
                          else
                             P_SN_rad=3.479830d42*(rbreak/rnod)**1.5
                          end if
                          !Thermal energy
                          alpha=-7.83*(rho_average(isnr)/scale_nH/100.)**(0.024)
                          rcool=3.021*3.08d18/scale_l*(rho_average(isnr)/scale_nH/100.)**(-0.421)
                          rrise=5.472*3.08d18/scale_l*(rho_average(isnr)/scale_nH/100.)**(-0.403)
                          !alpha=-11.3*(rhocool/scale_nH/100.)**(0.072)
                          !rcool=6.35*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.418)
                          !rrise=9.19*3.08d18/scale_l*(rhocool/scale_nH/100.)**(-0.439)
                          if (rc<rcool) then
                             E_SN_th=7.1d50
                          else
                             if (rcool <= rc .and. rc < rrise) then
                                E_SN_th=7.1d50*(rc/rcool)**alpha
                             else
                                E_SN_th=7.1d50*(rrise/rcool)**alpha
                             end if
                          end if
                          
                          ! Update the gas properties
                          
                          ! Thermal and kinetic energy have been alreadyconverted to code units
                          ! in condinit. They are global variables. They should not be converted here.
                          
                          ! Initial pressure in the injection radius, in code units.
                          ! Use for injection in thermal form.
                          Vi = (4.0*pi*rc**3.0)/3.0 ! code units
                          Pi_ = (2.0/3.0)*(E_SN_th/Vi/(scale_v**2.0*scale_d*scale_l**3.0))
                          
                          ! Density of SN ejecta, assuming uniformly distributed within rc, as
                          ! in Thornton et al. (1998).
                          M_ej = 5.97*1d33 ! 3 Msun, cgs
                          
                          ! Mass from cgs to code units.
                          M_ej = M_ej/(scale_d*scale_l**3.0)
                          rho_ej = M_ej/Vi
                          
                          ! Radial velocity of the gas within rc
                          vri = 3.0*P_SN_rad/(4.0*pi*(rho_ej+rho_average(isnr))*rc**3.0)/(scale_v*scale_d*scale_l**3.0)
                          
                          ! As in courant_file.f90.
                          mu = 0.6
                          
                          vxi = dxx*vri/sqrt(drr)
                          vyi = dyy*vri/sqrt(drr)
                          vzi = dzz*vri/sqrt(drr)
                          
                          ! This is a what is done when the blast-wave is simply added
                          uold(ind_cell(i),1)=rho_ej+rho_average(isnr) 
                          uold(ind_cell(i),2)=uold(ind_cell(i),2)+(rho_ej+rho_average(isnr))*vxi
                          uold(ind_cell(i),3)=uold(ind_cell(i),3)+(rho_ej+rho_average(isnr))*vyi
                          uold(ind_cell(i),4)=uold(ind_cell(i),4)+(rho_ej+rho_average(isnr))*vzi
                          uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5*(rho_ej+rho_average(isnr))*(vxi**2+vyi**2+vzi**2)+Pi_/(gamma-1.0d0)
                          if(metal)uold(ind_cell(i),6)=uold(ind_cell(i),6)+yield*rho_ej ! metal yield
                          
                          if (myid == 1) write(funit,'(A, I10.1)') 'Explosion number:',nblast_current_internal
                       end if
                    end do
                 end if
              end do
           end do
           ! End loop over cells
        end do
        ! End loop over grids
     end do
     ! End loop over levels
  end if

  deallocate(xsnr)
  deallocate(ysnr)
  deallocate(zsnr)
  deallocate(rho_average)
  deallocate(rho_average_all)
  deallocate(weight)
  deallocate(weight_all)
  deallocate(num_snr_local)
  
111 format('   Entering blast_wave_feedback for level ',I2)
  
end subroutine blast_wave_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,sSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! implement this using exploding GMC particles. 
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::ip,icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,info,iSN,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::scale,dx_min,vol_min,nISM,nCOM,d0,mstar
  integer::nx_loc
  integer,dimension(:),allocatable::ind_part,ind_grid
  logical,dimension(:),allocatable::ok_free
  integer ,dimension(:),allocatable::indSN
  real(dp),dimension(:),allocatable::mSN,sSN,ZSN,m_gas,vol_gas,ekBlast
  real(dp),dimension(:,:),allocatable::xSN,vSN,u_gas,dq

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering make_sn'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Initial star particle mass
  mstar=n_star/(scale_nH*aexp**3)*vol_min

  ! Lifetime of Giant Molecular Clouds from Myr to code units
  t0=10.*(1d6*365.*24.*3600.)/scale_t

  !------------------------------------------------------
  ! Gather GMC particles eligible for disruption
  !------------------------------------------------------
  nSN_loc=0
  ! Loop over levels
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count old enough GMC particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).le.0.and. tp(ipart).lt.(t-t0))then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
            end do
        endif
        nSN_loc=nSN_loc+npart2   ! Add SNe to the total
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  ! End loop over levels
  nSN_icpu=0
  nSN_icpu(myid)=nSN_loc
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all
#endif

  nSN_tot=sum(nSN_icpu(1:ncpu))

  if (nSN_tot .eq. 0) return
  
  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of GMC to explode=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN(1:nSN_tot,1:3),vSN(1:nSN_tot,1:3))
  allocate(mSN(1:nSN_tot),sSN(1:nSN_tot),ZSN(1:nSN_tot))
  xSN=0.;vSN=0.;mSN=0.;sSN=0.;ZSN=0.
  ! Allocate arrays for particles index and parent grid
  if(nSN_loc>0)then
     allocate(ind_part(1:nSN_loc),ind_grid(1:nSN_loc),ok_free(1:nSN_loc))
  endif

  !------------------------------------------------------
  ! Store position and mass of the GMC into the SN array
  !------------------------------------------------------
  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif
  ! Loop over levels
  ip=0
  do icpu=1,ncpu
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).le.0.and. tp(ipart).lt.(t-t0))then
                 iSN=iSN+1
                 xSN(iSN,1)=xp(ipart,1)
                 xSN(iSN,2)=xp(ipart,2)
                 xSN(iSN,3)=xp(ipart,3)
                 vSN(iSN,1)=vp(ipart,1)
                 vSN(iSN,2)=vp(ipart,2)
                 vSN(iSN,3)=vp(ipart,3)
                 mSN(iSN)=mp(ipart)
                 sSN(iSN)=dble(-idp(ipart))*mstar
                 if(metal)ZSN(iSN)=zp(ipart)
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_part(ip)=ipart
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        igrid=next(igrid)   ! Go to next grid
     end do
  end do 
  ! End loop over levels

  ! Remove GMC particle
  IF(nSN_loc>0)then
     ok_free=.true.
     call remove_list(ind_part,ind_grid,ok_free,nSN_loc)
     call add_free_cond(ind_part,ok_free,nSN_loc)
     deallocate(ind_part,ind_grid,ok_free)
  endif

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:3),vSN_all(1:nSN_tot,1:3),mSN_all(1:nSN_tot),sSN_all(1:nSN_tot),ZSN_all(1:nSN_tot))
  call MPI_ALLREDUCE(xSN,xSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vSN,vSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sSN,sSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZSN,ZSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  xSN=xSN_all
  vSN=vSN_all
  mSN=mSN_all
  sSN=sSN_all
  ZSN=ZSN_all
  deallocate(xSN_all,vSN_all,mSN_all,sSN_all,ZSN_all)
#endif

  nSN=nSN_tot
  allocate(m_gas(1:nSN),u_gas(1:nSN,1:3),vol_gas(1:nSN),dq(1:nSN,1:3),ekBlast(1:nSN))
  allocate(indSN(1:nSN))

  ! Compute the grid discretization effects
  call average_SN(xSN,vol_gas,dq,ekBlast,indSN,nSN)

  ! Modify hydro quantities to account for a Sedov blast wave
  call Sedov_blast(xSN,vSN,mSN,sSN,ZSN,indSN,vol_gas,dq,ekBlast,nSN)

  deallocate(xSN,vSN,mSN,sSN,ZSN,indSN,m_gas,u_gas,vol_gas,dq,ekBlast)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

end subroutine kinetic_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_SN(xSN,vol_gas,dq,ekBlast,ind_blast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,j,iSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,d,u,v,w,ek,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  real(dp),dimension(1:nSN)::mSN,m_gas,vol_gas,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,u_gas,dq,u2Blast
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN)::m_gas_all,vol_gas_all,ekBlast_all
  real(dp),dimension(1:nSN,1:3)::u_gas_all,dq_all,u2Blast_all
#endif
  logical ,dimension(1:nvector),save::ok

  if(nSN==0)return
  if(verbose)write(*,*)'Entering average_SN'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0.0;dq=0.0;u2Blast=0.0;ekBlast=0.0;ind_blast=-1

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                    if(dr_SN.lt.rmax2)then
                       vol_gas(iSN)=vol_gas(iSN)+vol_loc
                       ! Take account for grid effects on the conservation of the
                       ! normalized linear momentum
                       u=dxx/rmax
                       v=dyy/rmax
                       w=dzz/rmax
                       ! Add the local normalized linear momentum to the total linear
                       ! momentum of the blast wave (should be zero with no grid effect)
                       dq(iSN,1)=dq(iSN,1)+u*vol_loc
                       dq(iSN,2)=dq(iSN,2)+v*vol_loc
                       dq(iSN,3)=dq(iSN,3)+w*vol_loc
                       u2Blast(iSN,1)=u2Blast(iSN,1)+u*u*vol_loc
                       u2Blast(iSN,2)=u2Blast(iSN,2)+v*v*vol_loc
                       u2Blast(iSN,3)=u2Blast(iSN,3)+w*w*vol_loc
                    endif
                    if(dr_cell.le.dx_loc/2.0)then
                       ind_blast(iSN)=ind_cell(i)
                       ekBlast  (iSN)=vol_loc
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(vol_gas,vol_gas_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dq     ,dq_all     ,nSN*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(u2Blast,u2Blast_all,nSN*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekBlast,ekBlast_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas=vol_gas_all
  dq     =dq_all
  u2Blast=u2Blast_all
  ekBlast=ekBlast_all
#endif
  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        dq(iSN,1)=dq(iSN,1)/vol_gas(iSN)
        dq(iSN,2)=dq(iSN,2)/vol_gas(iSN)
        dq(iSN,3)=dq(iSN,3)/vol_gas(iSN)
        u2Blast(iSN,1)=u2Blast(iSN,1)/vol_gas(iSN)
        u2Blast(iSN,2)=u2Blast(iSN,2)/vol_gas(iSN)
        u2Blast(iSN,3)=u2Blast(iSN,3)/vol_gas(iSN)
        u2=u2Blast(iSN,1)-dq(iSN,1)**2
        v2=u2Blast(iSN,2)-dq(iSN,2)**2
        w2=u2Blast(iSN,3)-dq(iSN,3)**2
        ekBlast(iSN)=max(0.5d0*(u2+v2+w2),0.0d0)
     endif
  end do

  if(verbose)write(*,*)'Exiting average_SN'

end subroutine average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine Sedov_blast(xSN,vSN,mSN,sSN,ZSN,indSN,vol_gas,dq,ekBlast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,j,iSN,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,d,u,v,w,ek,u_r,ESN
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,sSN,ZSN,m_gas,p_gas,d_gas,d_metal,vol_gas,uSedov,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,u_gas,dq
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector),save::ok

  if(nSN==0)return
  if(verbose)write(*,*)'Entering Sedov_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax
  
  ! Supernova specific energy from cgs to code units
  ESN=(1d51/(10d0*2d33))/scale_v**2

  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        d_gas(iSN)=mSN(iSN)/vol_gas(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/vol_gas(iSN)
        if(ekBlast(iSN)==0d0)then
           p_gas(iSN)=eta_sn*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=0d0
        else
           p_gas(iSN)=(1d0-f_ek)*eta_sn*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=sqrt(f_ek*eta_sn*sSN(iSN)*ESN/mSN(iSN)/ekBlast(iSN))
        endif
     else
        d_gas(iSN)=mSN(iSN)/ekBlast(iSN)
        p_gas(iSN)=eta_sn*sSN(iSN)*ESN/ekBlast(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/ekBlast(iSN)
     endif
  end do

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    if(dr_SN.lt.rmax2)then
                       ! Compute the mass density in the cell
                       uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas(iSN)
                       ! Compute the metal density in the cell
                       if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_metal(iSN)
                       ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=uSedov(iSN)*(dxx/rmax-dq(iSN,1))+vSN(iSN,1)
                       v=uSedov(iSN)*(dyy/rmax-dq(iSN,2))+vSN(iSN,2)
                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))+vSN(iSN,3)
                       ! Add each momentum component of the blast wave to the gas
                       uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas(iSN)*u
                       uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas(iSN)*v
                       uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas(iSN)*w
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5*d_gas(iSN)*(u*u+v*v+w*w)+p_gas(iSN)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do iSN=1,nSN
     if(vol_gas(iSN)==0d0)then
        u=vSN(iSN,1)
        v=vSN(iSN,2)
        w=vSN(iSN,3)
        if(indSN(iSN)>0)then
           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas(iSN)
           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas(iSN)*u
           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas(iSN)*v
           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas(iSN)*w
           uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas(iSN)*0.5*(u*u+v*v+w*w)+p_gas(iSN)
           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_metal(iSN)
        endif
     endif
  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine Sedov_blast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
