!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine SN_feedback(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info2,dummy_io
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ivar
  integer::ig,ip,npart1,npart2,icpu,ilun=0,idim
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  character(LEN=80)::filename,filedir,fileloc,filedirini
  character(LEN=5)::nchar,ncharcpu
  logical::file_exist
  integer,parameter::tag=1120

  if(sf_log_properties) then
     call title(ifout-1,nchar)
     if(IOGROUPSIZEREP>0) then
        filedirini='output_'//TRIM(nchar)//'/'
        filedir='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
     else
        filedir='output_'//TRIM(nchar)//'/'
     endif
     filename=TRIM(filedir)//'stars_'//TRIM(nchar)//'.out'
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

     inquire(file=fileloc,exist=file_exist)
     if(.not.file_exist) then
        open(ilun, file=fileloc, form='formatted')
        write(ilun,'(A24)',advance='no') '# event id  ilevel  mp  '
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'xp',idim,'  '
        enddo
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'vp',idim,'  '
        enddo
        do ivar=1,nvar
           if(ivar.ge.10) then
              write(ilun,'(A1,I2,A2)',advance='no') 'u',ivar,'  '
           else
              write(ilun,'(A1,I1,A2)',advance='no') 'u',ivar,'  '
           endif
        enddo
        write(ilun,'(A5)',advance='no') 'tag  '
        write(ilun,'(A1)') ' '
     else
        open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
     endif
  endif


  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather star particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if ( is_star(typep(ipart)) ) then
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
              if ( is_star(typep(ipart)) ) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 if (.not.adaptive_feedback)then
                    call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 else
                    call adp_feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel) ! Davide Martizzi: this is the new routine
                 end if
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
     if(ip>0)then ! Davide Martizzi: added 2 options
        if (.not.adaptive_feedback)then
           call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
        else
           call adp_feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel) ! Davide Martizzi: this is the new routine
        end if
     end if
  end do
  ! End loop over cpus

  if(sf_log_properties) close(ilun)

111 format('   Entering thermal_feedback for level ',I2)

end subroutine SN_feedback
#endif
!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
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
  ! This routine is called by subroutine feedback. Each stellar particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! unew.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ivar,ilun=0
  real(kind=8)::RandNum
  real(dp)::SN_BOOST,mstar,dx_min,vol_min
  real(dp)::t0,ESN,mejecta,zloss,e,uvar
  real(dp)::ERAD,RAD_BOOST,tauIR,msne_min,mstar_max,eta_sn2
  real(dp)::delta_x,tau_factor,rad_factor
  real(dp)::dx,dx_loc,scale,birth_time,current_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
#if NENER>0
  integer::irad
#endif

  if(sf_log_properties) ilun=myid+10
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif
  msne_min=mass_sne_min*2d33/(scale_d*scale_l**3)
  mstar_max=mass_star_max*2d33/(scale_d*scale_l**3)

  ! Compute stochastic boost to account for target GMC mass
  SN_BOOST=MAX(mass_gmc*2d33/(scale_d*scale_l**3)/mstar,1d0)

  ! Massive star lifetime from Myr to code units
  ! Davide Martizzi: slowly ramp up the t0 from 0 to t_sne
  if(use_proper_time)then
     t0=t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2)
     current_time=texp
     !if(current_time.lt.t0) t0=t0/2
     !t0 = MIN(current_time,t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2))
  else
     t0=t_sne*1d6*(365.*24.*3600.)/scale_t
     current_time=t
     !if(current_time.lt.t0) t0=t0/2
     !t0 = MIN(current_time,t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2))
  endif

  ! Type II supernova specific energy from cgs to code units
  ESN=1d51/(10.*2d33)/scale_v**2 ! Davide Martizzi: This is the energy in the original model

  ! Life time radiation specific energy from cgs to code units
  ERAD=1d53/(10.*2d33)/scale_v**2 ! Davide Martizzi: This is the energy in the original model

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

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
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
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
        vol_loc(j)=vol_loc(j)*2**ndim ! ilevel-1 cell volume
     end if
  end do

  ! Compute individual time steps
  do j=1,np
     dteff(j)=dtnew(levelp(ind_part(j)))
  end do

  if(use_proper_time)then
     do j=1,np
        dteff(j)=dteff(j)*aexp**2
     end do
  endif

  ! Reset ejected mass, metallicity, thermal energy
  do j=1,np
     mloss(j)=0d0
     mzloss(j)=0d0
     ethermal(j)=0d0
  end do

  ! Compute stellar mass loss and thermal feedback due to supernovae
  if(f_w==0)then
     do j=1,np
        if(current_time.lt.t0)then ! Davide Martizzi
           birth_time=tp(ind_part(j))-t0
        else
           birth_time=tp(ind_part(j))
        end if
        ! Make sure that we don't count feedback twice
        if(birth_time.lt.(current_time-t0).and.birth_time.ge.(current_time-t0-dteff(j)))then
           eta_sn2   = eta_sn
           if(sf_imf)then
              if(mp(ind_part(j)).le.mstar_max)then
                 if(mp(ind_part(j)).ge.msne_min) eta_sn2 = eta_ssn
                 if(mp(ind_part(j)).lt.msne_min) eta_sn2 = 0.0
              endif
           endif
           ! Stellar mass loss
           mejecta=eta_sn2*mp(ind_part(j))
           mloss(j)=mloss(j)+mejecta/vol_loc(j)
           ! Thermal energy
           ethermal(j)=ethermal(j)+mejecta*ESN/vol_loc(j)
           ! Metallicity
           if(metal)then
              zloss=yield+(1d0-yield)*zp(ind_part(j))
              mzloss(j)=mzloss(j)+mejecta*zloss/vol_loc(j)
           endif
           ! Reduce star particle mass ! Davide Martizzi
           if(current_time.ge.t0) mp(ind_part(j))=mp(ind_part(j))-mejecta
           !if(use_proper_time)then           
           !   if(t0.eq.t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2)) mp(ind_part(j))=mp(ind_part(j))-mejecta
           !else
           !   if(t0.eq.t_sne*1d6*(365.*24.*3600.)/(scale_t)) mp(ind_part(j))=mp(ind_part(j))-mejecta
           !end if

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
           if(sf_log_properties) then
              write(ilun,'(I10)',advance='no') 1
              write(ilun,'(2I10,E24.12)',advance='no') idp(ind_part(j)),ilevel,mp(ind_part(j))
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') xp(ind_part(j),idim)
              enddo
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') vp(ind_part(j),idim)
              enddo
              write(ilun,'(E24.12)',advance='no') unew(indp(j),1)
              do ivar=2,nvar
                 if(ivar.eq.ndim+2)then
                    e=0.0d0
                    do idim=1,ndim
                       e=e+0.5*unew(ind_cell(i),idim+1)**2/max(unew(ind_cell(i),1),smallr)
                    enddo
#if NENER>0
                    do irad=0,nener-1
                       e=e+unew(ind_cell(i),inener+irad)
                    enddo
#endif
#ifdef SOLVERmhd
                    do idim=1,ndim
                       e=e+0.125d0*(unew(ind_cell(i),idim+ndim+2)+unew(ind_cell(i),idim+nvar))**2
                    enddo
#endif
                    ! Temperature
                    uvar=(gamma-1.0)*(unew(ind_cell(i),ndim+2)-e)*scale_T2
                 else
                    uvar=unew(indp(j),ivar)
                 endif
                 write(ilun,'(E24.12)',advance='no') uvar/unew(indp(j),1)
              enddo
              write(ilun,'(I10)',advance='no') typep(ind_part(i))%tag
              write(ilun,'(A1)') ' '
           endif
        endif
     end do
  endif

  ! Update hydro variables due to feedback

  ! For IR radiation trapping,
  ! we use a fixed length to estimate the column density of gas
  delta_x=200.*3d18
  if(metal)then
     tau_factor=kappa_IR*delta_x*scale_d/0.02
  else
     tau_factor=kappa_IR*delta_x*scale_d*z_ave
  endif
  rad_factor=ERAD/ESN

  do j=1,np

     ! Infrared photon trapping boost
     if(metal)then
        tauIR=tau_factor*max(uold(indp(j),imetal),smallr)
     else
        tauIR=tau_factor*max(uold(indp(j),1),smallr)
     endif
     if(uold(indp(j),1)*scale_nH > 10.)then
        RAD_BOOST=rad_factor*(1d0-exp(-tauIR))
     else
        RAD_BOOST=0.0
     endif

     ! Specific kinetic energy of the star
     ekinetic(j)=0.5*(vp(ind_part(j),1)**2 &
          &          +vp(ind_part(j),2)**2 &
          &          +vp(ind_part(j),3)**2)

     ! Update hydro variable in NGP cell
     unew(indp(j),1)=unew(indp(j),1)+mloss(j)
     unew(indp(j),2)=unew(indp(j),2)+mloss(j)*vp(ind_part(j),1)
     unew(indp(j),3)=unew(indp(j),3)+mloss(j)*vp(ind_part(j),2)
     unew(indp(j),4)=unew(indp(j),4)+mloss(j)*vp(ind_part(j),3)
     unew(indp(j),5)=unew(indp(j),5)+mloss(j)*ekinetic(j)+ &
          & ethermal(j)*(1d0+RAD_BOOST)
  end do

  ! Add metals
  if(metal)then
     do j=1,np
        unew(indp(j),imetal)=unew(indp(j),imetal)+mzloss(j)
     end do
  endif

  ! Add delayed cooling switch variable
  if(delayed_cooling)then
     do j=1,np
        unew(indp(j),idelay)=unew(indp(j),idelay)+mloss(j)
     end do
  endif

end subroutine feedbk
#endif
!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine adp_feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each stellar particle
  ! dumps mass, momentum and energy in the surrounding 9+1 cells
  !-----------------------------------------------------------------------
  integer::i,j,k,l,idim,nx_loc,ivar,ilun=0
  integer::indcell_inj,indgrid_inj
  real(kind=8)::RandNum
  real(dp)::mstar,dx_min,vol_min
  real(dp)::t0,ESN,mejecta,zloss,e,uvar
  real(dp)::ERAD,RAD_BOOST,tauIR,msne_min,mstar_max,eta_sn2
  real(dp)::delta_x,tau_factor,rad_factor
  real(dp)::dx,dx_loc,scale,birth_time,current_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::rr
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:ndim)::kfather,lson
  real(dp),dimension(1:ndim)::x_inj,drvec
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic
  real(dp),dimension(1:nvector),save::pradial,dteff,vol_loc
  real(dp),dimension(1:nvector),save::rho_FB,met_FB,vol_FB,R_inj_FB
  real(dp),dimension(1:nvector,1:ndim),save::x
  real(dp),dimension(1:nvector,1:ndim),save::x_FB
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::flag_SN
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
#if NENER>0
  integer::irad
#endif
  ! Parameters of subgrid model
  real(dp)::alpha,rcool,rrise,rnod,rbreak,met_floor

  if(sf_log_properties) ilun=myid+10
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif
  msne_min=mass_sne_min*2d33/(scale_d*scale_l**3)
  mstar_max=mass_star_max*2d33/(scale_d*scale_l**3)

  ! Compute individual time steps
  do j=1,np
     dteff(j)=dtnew(levelp(ind_part(j)))
  end do

  if(use_proper_time)then
     do j=1,np
        dteff(j)=dteff(j)*aexp**2
     end do
  endif

  ! Massive star lifetime from Myr to code units
  ! Davide Martizzi: slowly ramp up the t0 from 0 to t_sne
  if(use_proper_time)then
     t0=t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2)
     current_time=texp
     !if(current_time.lt.t0) t0=t0/2
     !t0 = MIN(current_time,t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2))
  else
     t0=t_sne*1d6*(365.*24.*3600.)/scale_t
     current_time=t
     !if(current_time.lt.t0) t0=t0/2
     !t0 = MIN(current_time,t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2))
  endif

  ! Type II supernova specific energy from cgs to code units
  ESN=1d51/(10.*2d33)/scale_v**2 ! Davide Martizzi: This is the energy in the original model, now not used

  ! Life time radiation specific energy from cgs to code units
  ERAD=1d53/(10.*2d33)/scale_v**2 ! Davide Martizzi: This is the energy in the original model

  ! Davide Martizzi: Reset local ISM properties and
  ! flag particles that should inject SNe
  do j=1,np
     rho_FB(j)=0.0d0
     met_FB(j)=0.0d0
     vol_FB(j)=1.0d-20
     R_inj_FB(j)=0.0d0
     flag_SN(j)=0
     do idim=1,ndim
        x_FB(j,idim)=0.0d0
        x_inj(idim)=0.0d0
     end do
     if(current_time.lt.t0)then ! Davide Martizzi
        birth_time=tp(ind_part(j))-t0
     else
        birth_time=tp(ind_part(j))
     end if
     if(birth_time.lt.(current_time-t0).and.birth_time.ge.(current_time-t0-dteff(j))) flag_SN(j)=1
  end do

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

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
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
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
        vol_loc(j)=vol_loc(j)*2**ndim ! ilevel-1 cell volume
     end if
  end do

  ! Davide Martizzi: Compute ISM properties for all particles
  if(f_w==0)then
     do j=1,np
        if(ok(j).and.flag_SN(j).gt.0)then
           ! Feedback injection radius and position
           R_inj_FB(j) = 2.05*dx*SQRT(dble(ndim))
           do idim=1,ndim
              x_FB(j,idim) = xp(ind_part(j),idim)/scale+skip_loc(idim)
           end do
           ! Loop over 27 neighbor father cells 
           do k=1,3**ndim
              ! Grid index associated to father cell
              indgrid_inj=son(nbors_father_cells(ind_grid_part(j),k))
              if (indgrid_inj.gt.0) then 
                 ! If neighbor father cell has sons
                 ! include neighbor cell at the same level
                 kfather(3) = (k-1)/9
                 kfather(2) = (k-1-9*kfather(3))/3
                 kfather(1) = (k-1-9*kfather(3)-3*kfather(2))                 
                 do l=1,2**ndim ! loop on all cells in the neighbour father grid
                    lson(3) = (l-1)/4
                    lson(2) = (l-1-4*lson(3))/2
                    lson(1) = (l-1-4*lson(3)-2*lson(2))
                    do idim=1,ndim
                       ! Grid center
                       x_inj(idim) = xg(indgrid_inj,idim)
                       ! Add offset with respect to grid center
                       x_inj(idim) = x_inj(idim)+(dble(lson(idim))-0.5)*dx
                    end do
                    ! Only add if radius is within R_inj_FB
                    rr = SQRT((x_inj(1)-x_FB(j,1))**2+(x_inj(2)-x_FB(j,2))**2+(x_inj(3)-x_FB(j,3))**2)
                    if (rr.lt.R_inj_FB(j)) then
                       indcell_inj = ncoarse+(l-1)*ngridmax+indgrid_inj
                       vol_FB(j) = vol_FB(j) + dx_loc**ndim
                       rho_FB(j) = rho_FB(j) + unew(indcell_inj,1)*dx_loc**ndim
                       if(metal) met_FB(j) = met_FB(j) + unew(indcell_inj,imetal)*dx_loc**ndim
                    end if
                 end do
              else
                 ! If neighbor father cell does not have sons
                 ! include neighbor father cell at lower level
                 kfather(3) = (k-1)/9
                 kfather(2) = (k-1-9*kfather(3))/3
                 kfather(1) = (k-1-9*kfather(3)-3*kfather(2))
                 do idim=1,ndim
                    x_inj(idim) = xg(ind_grid(ind_grid_part(j)),idim)-2.0D0*dx+2*dx*dble(kfather(idim))
                 end do
                 ! Only add if radius is within 2*R_inj_FB
                 rr = SQRT((x_inj(1)-x_FB(j,1))**2+(x_inj(2)-x_FB(j,2))**2+(x_inj(3)-x_FB(j,3))**2)
                 if (rr.lt.R_inj_FB(j)) then
                    indcell_inj = nbors_father_cells(ind_grid_part(j),k)
                    vol_FB(j) = vol_FB(j) + dx_loc**ndim*2**ndim
                    rho_FB(j) = rho_FB(j) + unew(indcell_inj,1)*dx_loc**ndim*2**ndim
                    if(metal) met_FB(j) = met_FB(j) + unew(indcell_inj,imetal)*dx_loc**ndim*2**ndim
                 end if
              end if
           end do
        end if
     end do
     
     ! Normalize ISM properties 
     do j=1,np
        rho_FB(j) = rho_FB(j)/vol_FB(j)
        met_FB(j) = met_FB(j)/vol_FB(j)
     end do

  end if
  
  ! Reset ejected mass, metallicity, thermal energy, radial momentum
  do j=1,np
     mloss(j)=0d0
     mzloss(j)=0d0
     ethermal(j)=0d0
     pradial(j)=0.0d0
  end do

  ! For IR radiation trapping,
  ! we use a fixed length to estimate the column density of gas
  delta_x=200.*3d18
  if(metal)then
     tau_factor=kappa_IR*delta_x*scale_d/0.02
  else
     tau_factor=kappa_IR*delta_x*scale_d*z_ave
  endif
  rad_factor=ERAD/ESN

  ! Compute stellar mass loss and feedback due to supernovae
  if(f_w==0)then
     do j=1,np
        if(ok(j).and.flag_SN(j).gt.0)then

           ! Infrared photon trapping boost
           if(metal)then
              tauIR=tau_factor*max(uold(indp(j),imetal),smallr)
           else
              tauIR=tau_factor*max(uold(indp(j),1),smallr)
           endif
           if(uold(indp(j),1)*scale_nH > 10.)then
              RAD_BOOST=rad_factor*(1d0-exp(-tauIR))
           else
              RAD_BOOST=0.0
           endif
           
           ! Make sure that we don't count SN feedback twice
           eta_sn2   = eta_sn
           if(sf_imf)then
              if(mp(ind_part(j)).le.mstar_max)then
                 if(mp(ind_part(j)).ge.msne_min) eta_sn2 = eta_ssn
                 if(mp(ind_part(j)).lt.msne_min) eta_sn2 = 0.0
              endif
           endif
           
           ! Stellar mass loss
           mejecta=eta_sn2*mp(ind_part(j))
           mloss(j)=mloss(j)+mejecta/vol_FB(j)
           
           ! Metallicity
           if(metal)then
              zloss=yield+(1d0-yield)*zp(ind_part(j))
              mzloss(j)=mzloss(j)+mejecta*zloss/vol_FB(j)
           endif
           
           ! Reduce star particle mass ! Davide Martizzi
           if(current_time.ge.t0) mp(ind_part(j))=mp(ind_part(j))-mejecta

           ! Parameters of SN subgrid model
           met_floor = MAX(2.0d-4,met_FB(j)/rho_FB(j)) ! the min metallicity for the model is 1d-2 Z_sun
           rbreak=4.001*3.08d18/scale_l*(rho_FB(j)*scale_nH/100.)**(-0.429)*(met_floor/0.02)**(-0.077)
           rnod=0.969*3.08d18/scale_l*(rho_FB(j)*scale_nH/100.)**(-0.330)*(met_floor/0.02)**(0.046)
           !rbreak=7.97*3.08d18/scale_l*(rho_FB(j)*scale_nH/100.)**(-0.458)*(met_floor/0.02)**(-0.058)
           !rnod=2.438*3.08d18/scale_l*(rho_FB(j)*scale_nH/100.)**(-0.349)*(met_floor/0.02)**(0.021)
           alpha=-7.83*(rho_FB(j)*scale_nH/100.)**(0.024)*(met_floor/0.02)**(0.050)
           rcool=3.021*3.08d18/scale_l*(rho_FB(j)*scale_nH/100.)**(-0.421)*(met_floor/0.02)**(-0.082)
           rrise=5.472*3.08d18/scale_l*(rho_FB(j)*scale_nH/100.)**(-0.403)*(met_floor/0.02)**(-0.074)
           !alpha=-11.3*(rho_FB(j)*scale_nH/100.)**(0.072)*(met_floor/0.02)**(0.070)
           !rcool=6.35*3.08d18/scale_l*(rho_FB(j)*scale_nH/100.)**(-0.418)*(met_floor/0.02)**(-0.050)
           !rrise=9.19*3.08d18/scale_l*(rho_FB(j)*scale_nH/100.)**(-0.439)*(met_floor/0.02)**(-0.067)           
           
           ! Radial momentum per SN
           if (R_inj_FB(j)<rbreak/scale) then
              pradial(j)=3.479830d42*(R_inj_FB(j)*scale/rnod)**1.5
           else
              pradial(j)=3.479830d42*(rbreak/rnod)**1.5
           end if
           ! Normalize by ejecta mass scaling and by volume
           pradial(j) = pradial(j)*(mejecta*scale_d*scale_l**3/3.0d0/2.0d33)**0.5
           pradial(j) = pradial(j)/(scale_d*scale_l**3*scale_v)/vol_FB(j)
              
           ! Thermal energy per SN 
           if (R_inj_FB(j)<rcool/scale) then
              ethermal(j)=7.1d50
           else
              if (rcool/scale <= R_inj_FB(j) .and. R_inj_FB(j) < rrise/scale) then
                 ethermal(j)=7.1d50*(R_inj_FB(j)*scale/rcool)**alpha
              else
                 ethermal(j)=7.1d50*(rrise/rcool)**alpha
              end if
           end if
           ! Normalize by ejecta mass and by volume
           ethermal(j)=ethermal(j)*(mejecta/3.0d0/2.0d33)/vol_FB(j)/scale_v**2
           
           ! Kinetic energy of the star plus SN feedback
           ekinetic(j)=0.5*mloss(j)*(vp(ind_part(j),1)**2 &
                &          +vp(ind_part(j),2)**2 &
                &          +vp(ind_part(j),3)**2) &
                &          +0.5*pradial(j)**2/(rho_FB(j) + mloss(j))
           
           if(sf_log_properties) then
              write(ilun,'(I10)',advance='no') 1
              write(ilun,'(2I10,E24.12)',advance='no') idp(ind_part(j)),ilevel,mp(ind_part(j))
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') xp(ind_part(j),idim)
              enddo
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') vp(ind_part(j),idim)
              enddo
              write(ilun,'(E24.12)',advance='no') unew(indp(j),1)
              do ivar=2,nvar
                 if(ivar.eq.ndim+2)then
                    e=0.0d0
                    do idim=1,ndim
                       e=e+0.5*unew(ind_cell(i),idim+1)**2/max(unew(ind_cell(i),1),smallr)
                    enddo
#if NENER>0
                    do irad=0,nener-1
                       e=e+unew(ind_cell(i),inener+irad)
                    enddo
#endif
#ifdef SOLVERmhd
                    do idim=1,ndim
                       e=e+0.125d0*(unew(ind_cell(i),idim+ndim+2)+unew(ind_cell(i),idim+nvar))**2
                    enddo
#endif
                    ! Temperature
                    uvar=(gamma-1.0)*(unew(ind_cell(i),ndim+2)-e)*scale_T2
                 else
                    uvar=unew(indp(j),ivar)
                 endif
                 write(ilun,'(E24.12)',advance='no') uvar/unew(indp(j),1)
              enddo
              write(ilun,'(I10)',advance='no') typep(ind_part(i))%tag
              write(ilun,'(A1)') ' '
           endif
        end if
        
        ! Feedback loop!
        ! Loop over 27 neighbor father cells
        do k=1,3**ndim
           ! Grid index associated to father cell
           indgrid_inj=son(nbors_father_cells(ind_grid_part(j),k))
           if (indgrid_inj.gt.0) then
              ! If neighbor father cell has sons
              ! include neighbor cell at the same level
              kfather(3) = (k-1)/9
              kfather(2) = (k-1-9*kfather(3))/3
              kfather(1) = (k-1-9*kfather(3)-3*kfather(2))
              do l=1,2**ndim ! loop on all cells in the neighbour father grid
                 lson(3) = (l-1)/4
                 lson(2) = (l-1-4*lson(3))/2
                 lson(1) = (l-1-4*lson(3)-2*lson(2))
                 do idim=1,ndim
                    ! Grid center
                    x_inj(idim) = xg(indgrid_inj,idim)
                    ! Add offset with respect to grid center
                    x_inj(idim) = x_inj(idim)+(dble(lson(idim))-0.5)*dx
                 end do
                 ! Only add if radius is within R_inj_FB
                 rr = SQRT((x_inj(1)-x_FB(j,1))**2+(x_inj(2)-x_FB(j,2))**2+(x_inj(3)-x_FB(j,3))**2)
                 if (rr.lt.R_inj_FB(j)) then
                    indcell_inj=ncoarse+(l-1)*ngridmax+indgrid_inj
                    drvec(1) = (x_inj(1)-x_FB(j,1))/(1.0d-20+rr)
                    drvec(2) = (x_inj(2)-x_FB(j,2))/(1.0d-20+rr)
                    drvec(3) = (x_inj(3)-x_FB(j,3))/(1.0d-20+rr)
                    ! Update hydro variables
                    !unew(indcell_inj,1) = rho_FB(j) + mloss(j)
                    unew(indcell_inj,1) = unew(indcell_inj,1) + mloss(j)
                    unew(indcell_inj,2) = unew(indcell_inj,2) + pradial(j)*drvec(1) + mloss(j)*vp(ind_part(j),1)
                    unew(indcell_inj,3) = unew(indcell_inj,3) + pradial(j)*drvec(2) + mloss(j)*vp(ind_part(j),2)
                    unew(indcell_inj,4) = unew(indcell_inj,4) + pradial(j)*drvec(3) + mloss(j)*vp(ind_part(j),3)
                    unew(indcell_inj,5) = unew(indcell_inj,5) + ekinetic(j) + ethermal(j)*(1d0+RAD_BOOST) 
                    !write(*,*)'internal', rho_FB(j), mloss(j), pradial(j), drvec(1), drvec(2), drvec(3), &
                    !     & vp(ind_part(j),1), vp(ind_part(j),2), vp(ind_part(j),3), ekinetic(j), ethermal(j)
                    ! Add metals
                    !if(metal)unew(indcell_inj,imetal)=met_FB(j)+mzloss(j)
                    if(metal)unew(indcell_inj,imetal)=unew(indcell_inj,imetal)+mzloss(j)
                    !if(met_FB(j).lt.0.0d0)write(*,*)'BRO!',rho_FB(j),met_FB(j),mzloss(j),zp(ind_part(j))
                    ! Add delayed cooling switch variable
                    if(delayed_cooling)unew(indcell_inj,idelay)=unew(indcell_inj,idelay)+mloss(j)
                 end if
              end do
           else
              ! If neighbor father cell does not have sons
              ! include neighbor father cell at lower level
              kfather(3) = (k-1)/9
              kfather(2) = (k-1-9*kfather(3))/3
              kfather(1) = (k-1-9*kfather(3)-3*kfather(2))
              do idim=1,ndim
                 x_inj(idim) = xg(ind_grid(ind_grid_part(j)),idim)-2.0D0*dx+2*dx*dble(kfather(idim))
              end do
              ! Only add if radius is within 2*R_inj_FB
              rr = SQRT((x_inj(1)-x_FB(j,1))**2+(x_inj(2)-x_FB(j,2))**2+(x_inj(3)-x_FB(j,3))**2)
              if (rr.lt.R_inj_FB(j)) then
                 indcell_inj= nbors_father_cells(ind_grid_part(j),k)
                 drvec(1) = (x_inj(1)-x_FB(j,1))/(1.0d-20+rr)
                 drvec(2) = (x_inj(2)-x_FB(j,2))/(1.0d-20+rr)
                 drvec(3) = (x_inj(3)-x_FB(j,3))/(1.0d-20+rr)
                 ! Update hydro variables
                 !unew(indcell_inj,1) = rho_FB(j) + mloss(j)
                 unew(indcell_inj,1) = unew(indcell_inj,1) + mloss(j)
                 unew(indcell_inj,2) = unew(indcell_inj,2) + pradial(j)*drvec(1) + mloss(j)*vp(ind_part(j),1)
                 unew(indcell_inj,3) = unew(indcell_inj,3) + pradial(j)*drvec(2) + mloss(j)*vp(ind_part(j),2)
                 unew(indcell_inj,4) = unew(indcell_inj,4) + pradial(j)*drvec(3) + mloss(j)*vp(ind_part(j),3)
                 unew(indcell_inj,5) = unew(indcell_inj,5) + ekinetic(j) + ethermal(j)*(1d0+RAD_BOOST) 
                 !write(*,*)'external', rho_FB(j), mloss(j), pradial(j), drvec(1), drvec(2), drvec(3), &
                 !     & vp(ind_part(j),1), vp(ind_part(j),2), vp(ind_part(j),3), ekinetic(j), ethermal(j)
                 ! Add metals
                 !if(metal)unew(indcell_inj,imetal)=met_FB(j)+mzloss(j)
                 if(metal)unew(indcell_inj,imetal)=unew(indcell_inj,imetal)+mzloss(j)
                 !if(met_FB(j).lt.0.0d0)write(*,*)'BROBRO!',rho_FB(j),met_FB(j),mzloss(j),zp(ind_part(j))
                 ! Add delayed cooling switch variable
                 if(delayed_cooling)unew(indcell_inj,idelay)=unew(indcell_inj,idelay)+mloss(j)
              end if
           end if
        end do
        ! End Feedback loop!
     end do
  endif

end subroutine adp_feedbk
#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,sSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! imolement this using exploding GMC particles.
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::ip,icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,iSN,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::current_time
  real(dp)::scale,dx_min,vol_min,mstar
  integer::nx_loc
  integer,dimension(:),allocatable::ind_part,ind_grid
  logical,dimension(:),allocatable::ok_free
  integer,dimension(:),allocatable::indSN
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

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif


  ! Lifetime of Giant Molecular Clouds from Myr to code units
  ! Massive star lifetime from Myr to code units
  if(use_proper_time)then
     t0=t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2)
     current_time=texp
  else
     t0=t_sne*1d6*(365.*24.*3600.)/scale_t
     current_time=t
  endif

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
              if ( is_debris(typep(ipart)) .and. tp(ipart).lt.(current_time-t0) ) then
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
              if ( is_debris(typep(ipart)) .and. tp(ipart).lt.(current_time-t0) ) then
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
  if(nSN_loc>0)then
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
  integer::info
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,iSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,u,v,w,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  real(dp),dimension(1:nSN)::vol_gas,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,dq,u2Blast
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN)::vol_gas_all,ekBlast_all
  real(dp),dimension(1:nSN,1:3)::dq_all,u2Blast_all
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
  integer::ilevel,iSN,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,u,v,w,ESN,mstar,eta_sn2,msne_min,mstar_max
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,vol_min
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,sSN,ZSN,p_gas,d_gas,d_metal,vol_gas,uSedov,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,dq
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
  vol_min=dx_min**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif
  msne_min=mass_sne_min*2d33/(scale_d*scale_l**3)
  mstar_max=mass_star_max*2d33/(scale_d*scale_l**3)
  ! Supernova specific energy from cgs to code units
  ESN=(1d51/(10d0*2d33))/scale_v**2

  do iSN=1,nSN
     eta_sn2    = eta_sn
     if(sf_imf)then
        if(mSN(iSN).le.mstar_max)then
           if(mSN(iSN).ge.msne_min) eta_sn2 = eta_ssn
           if(mSN(iSN).lt.msne_min) eta_sn2 = 0.0
        endif
     endif
     if(vol_gas(iSN)>0d0)then
        d_gas(iSN)=mSN(iSN)/vol_gas(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/vol_gas(iSN)
        if(ekBlast(iSN)==0d0)then
           p_gas(iSN)=eta_sn2*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=0d0
        else
           p_gas(iSN)=(1d0-f_ek)*eta_sn2*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=sqrt(f_ek*eta_sn2*sSN(iSN)*ESN/mSN(iSN)/ekBlast(iSN))
        endif
     else
        d_gas(iSN)=mSN(iSN)/ekBlast(iSN)
        p_gas(iSN)=eta_sn2*sSN(iSN)*ESN/ekBlast(iSN)
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
