subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
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
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do

        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)+fana(ind_leaf(i),idim) ! Davide Martizzi
              end do
           end do
        end if

        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do

        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,ndim+2)*vol
        end do

        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,ndim+2)*vol
        end do
        do ivar=1,ndim
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/max(uu(i,1),smallr)*vol
           end do
        end do
#if NENER>0
        do ivar=1,nener
           do i=1,nleaf
              eint_loc=eint_loc-uu(i,ndim+2+ivar)*vol
           end do
        end do
#endif

        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  dt_all=dt_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine

subroutine check_cons(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Check mass and energy conservation
  !----------------------------------------------------------------------
  integer::i,ivar,ind,ncache,igrid,iskip
  integer::nleaf,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dx,vol
  real(kind=8)::mass_loc,ekin_loc,eint_loc
  real(kind=8)::mass_all,ekin_all,eint_all
  real(dp),dimension(1:nvector,1:nvar),save::uu

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0

  ! Mesh spacing at that level
  dx=0.5D0**ilevel*boxlen
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
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
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do

        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do

        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,ndim+2)*vol
        end do

        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,ndim+2)*vol
        end do
        do ivar=1,ndim
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/max(uu(i,1),smallr)*vol
           end do
        end do
#if NENER>0
        do ivar=1,nener
           do i=1,nleaf
              eint_loc=eint_loc-uu(i,ndim+2+ivar)*vol
           end do
        end do
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all

111 format('   Entering check_cons for level ',I2)

end subroutine check_cons

! Davide Martizzi: Introduced this subroutine to include ceilings
subroutine ceiling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use hydro_parameters, only: gamma, smallr, smallc, cceil
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Enforce a temperature and velocity ceiling
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(dp)::ekin,eth,vv

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
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
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables and enforce ceiling
        do i=1,nleaf
           do ivar=1,nvar
              if(ISNAN(uold(ind_leaf(i),ivar)).and.ivar.le.ndim+2)then
                 uold(ind_leaf(i),1) = smallr
                 uold(ind_leaf(i),2) = 0.0d0
                 uold(ind_leaf(i),3) = 0.0d0
                 uold(ind_leaf(i),4) = 0.0d0
                 uold(ind_leaf(i),ndim+2) = smallr*smallc**2/gamma/(gamma-1)
                 if(metal)uold(ind_leaf(i),imetal) = smallr**2
              end if
              if (ivar.eq.1)then
                 if (metal.and.uold(ind_leaf(i),1).lt.smallr) then
                    uold(ind_leaf(i),imetal) = uold(ind_leaf(i),imetal)*smallr/uold(ind_leaf(i),1)
                 end if
                 uold(ind_leaf(i),1)=max(uold(ind_leaf(i),1),smallr)
              end if
              if (ivar.eq.2.or.ivar.eq.3.or.ivar.eq.4)then
                 vv = ABS(uold(ind_leaf(i),ivar)/uold(ind_leaf(i),1))
                 if (vv.gt.3.333*cceil) then
                    uold(ind_leaf(i),ndim+2)=uold(ind_leaf(i),ndim+2)-0.5*uold(ind_leaf(i),ivar)**2/uold(ind_leaf(i),1)
                    uold(ind_leaf(i),ivar)=uold(ind_leaf(i),ivar)/ABS(uold(ind_leaf(i),ivar))*uold(ind_leaf(i),1)*3.333*cceil
                    uold(ind_leaf(i),ndim+2)=uold(ind_leaf(i),ndim+2)+0.5*uold(ind_leaf(i),1)*(3.333*cceil)**2
                 end if
              end if
              if (ivar.eq.ndim+2)then
                 ekin = 0.5*(uold(ind_leaf(i),2)**2+uold(ind_leaf(i),3)**2+uold(ind_leaf(i),4)**2)/uold(ind_leaf(i),1)
                 eth = uold(ind_leaf(i),ivar)-ekin
                 if (eth .lt. uold(ind_leaf(i),1)*smallc**2/gamma/(gamma-1)) eth = uold(ind_leaf(i),1)*smallc**2/gamma/(gamma-1)
                 if (eth .gt. uold(ind_leaf(i),1)*cceil**2/gamma/(gamma-1)) eth = uold(ind_leaf(i),1)*cceil**2/gamma/(gamma-1)
                 uold(ind_leaf(i),ivar) = eth+ekin
              end if
              if (metal.and.ivar.eq.imetal)then
                 uold(ind_leaf(i),ivar) = max(uold(ind_leaf(i),ivar),smallr**2)
              end if
           end do
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering ceiling_fine for level ',I2)

end subroutine ceiling_fine
