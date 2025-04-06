!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use ddadi_class,       only: ddadi
   use vdscalar_class,    only: vdscalar
   use viscoelastic_class, only: viscoelastic
   use hypre_uns_class,   only: hypre_uns
   use hypre_str_class,   only: hypre_str
   use lowmach_class,     only: lowmach
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use string, only: str_medium
   implicit none
   private

   !> Get an LPT solver, a lowmach solver, and corresponding time tracker, plus a couple of linear solvers
   type(hypre_uns),     public :: ps
   type(hypre_str),     public :: vs
   type(viscoelastic), public :: ve
   type(ddadi),         public :: vfs, ves
   type(lowmach),       public :: fs
   type(vdscalar),      public :: vf        ! polymer Volume fraction scalar solver
   type(timetracker),   public :: time

   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,scfile,Difffile

   public :: simulation_init,simulation_run,simulation_final

   !> Work arrays and fluid properties
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,resRHO
   real(WP), dimension(:,:,:), allocatable :: SdiffU,SdiffV,SdiffW      ! solvent diffusion velocity
   real(WP), dimension(:,:,:), allocatable :: PdiffU,PdiffV,PdiffW      ! polymer diffusion velocity
   real(WP), dimension(:,:,:), allocatable :: SdiffUi,SdiffVi,SdiffWi
   real(WP), dimension(:,:,:), allocatable :: PdiffUi,PdiffVi,PdiffWi
   real(WP) :: visc, rho
   real(WP), dimension(:,:,:), allocatable :: resVF,  vfTmp     ! volume fraction residual
   logical, dimension(:,:,:), allocatable :: vfbqflag

   real(WP), dimension(:,:,:), allocatable :: relaxTime
   real(WP), dimension(:,:,:,:),   allocatable :: resSC,SCtmp,SR
   real(WP), dimension(:,:,:,:),   allocatable :: stress
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: gradUxx

   real(WP) :: diff0, ad
   real(WP) :: inflowVelocity

   real(WP), dimension(:,:,:), allocatable:: MixViscosity
   real(WP) :: maxMixViscosity
   real(WP) :: maxRelaxTime, minRelaxTime

   real(WP) :: threshold = 1.0E-5_WP

   real(WP) :: SMALL = 1.0E-12_WP


   !> Check for stabilization
   logical :: stabilization

   real(WP) :: maxSdiffU, maxSdiffV, maxSdiffW, minSdiffU, minSdiffV, minSdiffW


contains

   subroutine computeDiffU
      integer :: i,j,k
      real(WP) :: fx, fy, fz
      real(WP) :: phi2x, phi2y, phi2z     ! interpolation of volume fraction at the cell face

      SdiffU = 0.0_WP; SdiffV = 0.0_WP; SdiffW = 0.0_WP
      PdiffU = 0.0_WP; PdiffV = 0.0_WP; PdiffW = 0.0_WP

      call vf%metric_reset()
      call cfg%sync(vf%SC)

      do k=cfg%kmin_,cfg%kmax_+1
         do j=cfg%jmin_,cfg%jmax_+1
            do i=cfg%imin_,cfg%imax_+1
               ! fx, fy, fz are the diffusion flux of polymer D*grad(u)
               fx = sum(vf%itp_x(:,i,j,k)*vf%diff(i-1:i,j,k))*sum(vf%grdsc_x(:,i,j,k)*vf%SC(i-1:i,j,k))
               fy = sum(vf%itp_y(:,i,j,k)*vf%diff(i,j-1:j,k))*sum(vf%grdsc_y(:,i,j,k)*vf%SC(i,j-1:j,k))
               fz = sum(vf%itp_z(:,i,j,k)*vf%diff(i,j,k-1:k))*sum(vf%grdsc_z(:,i,j,k)*vf%SC(i,j,k-1:k))

               ! phi2x, phi2y, phi2z are the interpolation of volume fraction at the cell face
               phi2x = sum(vf%itp_x(:,i,j,k)*vf%SC(i-1:i,j,k))
               phi2y = sum(vf%itp_y(:,i,j,k)*vf%SC(i,j-1:j,k))
               phi2z = sum(vf%itp_z(:,i,j,k)*vf%SC(i,j,k-1:k))

               if (phi2x .gt. threshold) PdiffU(i,j,k) = -fx/phi2x
               if (phi2y .gt. threshold) PdiffV(i,j,k) = -fy/phi2y
               if (phi2z .gt. threshold) PdiffW(i,j,k) = -fz/phi2z

               if ((1.0_WP-phi2x) .gt. threshold) SdiffU(i,j,k) = fx/((1.0_WP-phi2x)+SMALL)!*(1.0_WP-0.5_WP*(1.0_WP+erf((phi2x-criticalPVF)/0.03_WP)))
               if ((1.0_WP-phi2y) .gt. threshold) SdiffV(i,j,k) = fy/((1.0_WP-phi2y)+SMALL)!*(1.0_WP-0.5_WP*(1.0_WP+erf((phi2y-criticalPVF)/0.03_WP)))
               if ((1.0_WP-phi2z) .gt. threshold) SdiffW(i,j,k) = fz/((1.0_WP-phi2z)+SMALL)!*(1.0_WP-0.5_WP*(1.0_WP+erf((phi2z-criticalPVF)/0.03_WP)))

            end do
         end do
      end do

      call cfg%sync(PdiffU); call cfg%sync(PdiffV); call cfg%sync(PdiffW)
      call cfg%sync(SdiffU); call cfg%sync(SdiffV); call cfg%sync(SdiffW)

      ! interpolate the diffusion velocity to the cell center, only for output
      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_-1
               PdiffUi(i,j,k) = sum(vf%itp_x(:,i,j,k)*PdiffU(i:i+1,j,k))
               SdiffUi(i,j,k) = sum(vf%itp_x(:,i,j,k)*SdiffU(i:i+1,j,k))
            end do
         end do
      end do
      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_-1
            do i=cfg%imino_,cfg%imaxo_
               PdiffVi(i,j,k) = sum(vf%itp_y(:,i,j,k)*PdiffV(i,j:j+1,k))
               SdiffVi(i,j,k) = sum(vf%itp_y(:,i,j,k)*SdiffV(i,j:j+1,k))
            end do
         end do
      end do
      do k=cfg%kmino_,cfg%kmaxo_-1
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
               PdiffWi(i,j,k) = sum(vf%itp_z(:,i,j,k)*PdiffW(i,j,k:k+1))
               SdiffWi(i,j,k) = sum(vf%itp_z(:,i,j,k)*SdiffW(i,j,k:k+1))
            end do
         end do
      end do
      ! Add last layer in each direction
      !if ((.not.cfg%xper).and.(cfg%iproc.eq.cfg%npx)) PdiffUi(cfg%imaxo,:,:)=PdiffU(cfg%imaxo,:,:) ; SdiffUi(cfg%imaxo,:,:)=SdiffU(cfg%imaxo,:,:)
      !if ((.not.cfg%yper).and.(cfg%jproc.eq.cfg%npy)) PdiffVi(:,cfg%jmaxo,:)=PdiffV(:,cfg%jmaxo,:) ; SdiffVi(:,cfg%jmaxo,:)=SdiffV(:,cfg%jmaxo,:)
      !if ((.not.cfg%zper).and.(cfg%kproc.eq.cfg%npz)) PdiffWi(:,:,cfg%kmaxo)=PdiffW(:,:,cfg%kmaxo) ; SdiffWi(:,:,cfg%kmaxo)=SdiffW(:,:,cfg%kmaxo)

      call cfg%sync(PdiffUi); call cfg%sync(SdiffUi)
      call cfg%sync(PdiffVi); call cfg%sync(SdiffVi)
      call cfg%sync(PdiffWi); call cfg%sync(SdiffWi)

   end subroutine computeDiffU

   subroutine applyExtraGradU
      ! add extra solvent diffusion velocity to calculate the gradU for viscoelastic solver
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: dudy,dudz,dvdx,dvdz,dwdx,dwdy
      real(WP) :: phi2x, phi2y, phi2z

      ! Compute dudx, dvdy, and dwdz first
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_

               phi2x = sum(vf%itp_x(:,i,j,k)*vf%SC(i:i+1,j,k))
               phi2y = sum(vf%itp_y(:,i,j,k)*vf%SC(i,j:j+1,k))
               phi2z = sum(vf%itp_z(:,i,j,k)*vf%SC(i,j,k:k+1))

               gradU(1,1,i,j,k) = gradU(1,1,i,j,k) + sum(fs%grdu_x(:,i,j,k)*SdiffU(i:i+1,j,k))!*(1.0_WP-0.5*(1+erf((phi2x-criticalPVF)/0.05_WP)))
               gradU(2,2,i,j,k) = gradU(2,2,i,j,k) + sum(fs%grdv_y(:,i,j,k)*SdiffV(i,j:j+1,k))!*(1.0_WP-0.5*(1+erf((phi2y-criticalPVF)/0.05_WP)))
               gradU(3,3,i,j,k) = gradU(3,3,i,j,k) + sum(fs%grdw_z(:,i,j,k)*SdiffW(i,j,k:k+1))!*(1.0_WP-0.5*(1+erf((phi2z-criticalPVF)/0.05_WP)))

            end do
         end do
      end do

      ! Allocate velocity gradient components
      allocate(dudy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(dudz(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(dvdx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(dvdz(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(dwdx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(dwdy(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

      ! Calculate components of the velocity gradient at their natural locations with an extra cell for interpolation
      do k=cfg%kmin_,cfg%kmax_+1
         do j=cfg%jmin_,cfg%jmax_+1
            do i=cfg%imin_,cfg%imax_+1
               phi2x = sum(vf%itp_x(:,i,j,k)*vf%SC(i-1:i,j,k))
               phi2y = sum(vf%itp_y(:,i,j,k)*vf%SC(i,j-1:j,k))
               phi2z = sum(vf%itp_z(:,i,j,k)*vf%SC(i,j,k-1:k))
               dudy(i,j,k)=sum(fs%grdu_y(:,i,j,k)*SdiffU(i,j-1:j,k))!*(1.0_WP-0.5*(1+erf((phi2y-criticalPVF)/0.05_WP)))
               dudz(i,j,k)=sum(fs%grdu_z(:,i,j,k)*SdiffU(i,j,k-1:k))!*(1.0_WP-0.5*(1+erf((phi2z-criticalPVF)/0.05_WP)))
               dvdx(i,j,k)=sum(fs%grdv_x(:,i,j,k)*SdiffV(i-1:i,j,k))!*(1.0_WP-0.5*(1+erf((phi2x-criticalPVF)/0.05_WP)))
               dvdz(i,j,k)=sum(fs%grdv_z(:,i,j,k)*SdiffV(i,j,k-1:k))!*(1.0_WP-0.5*(1+erf((phi2z-criticalPVF)/0.05_WP)))
               dwdx(i,j,k)=sum(fs%grdw_x(:,i,j,k)*SdiffW(i-1:i,j,k))!*(1.0_WP-0.5*(1+erf((phi2x-criticalPVF)/0.05_WP)))
               dwdy(i,j,k)=sum(fs%grdw_y(:,i,j,k)*SdiffW(i,j-1:j,k))!*(1.0_WP-0.5*(1+erf((phi2y-criticalPVF)/0.05_WP)))
            end do
         end do
      end do

      ! Interpolate off-diagonal components of the velocity gradient to the cell center
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               gradu(2,1,i,j,k) = gradu(2,1,i,j,k) + 0.25_WP*sum(dudy(i:i+1,j:j+1,k))
               gradu(3,1,i,j,k) = gradu(3,1,i,j,k) + 0.25_WP*sum(dudz(i:i+1,j,k:k+1))
               gradu(1,2,i,j,k) = gradu(1,2,i,j,k) + 0.25_WP*sum(dvdx(i:i+1,j:j+1,k))
               gradu(3,2,i,j,k) = gradu(3,2,i,j,k) + 0.25_WP*sum(dvdz(i,j:j+1,k:k+1))
               gradu(1,3,i,j,k) = gradu(1,3,i,j,k) + 0.25_WP*sum(dwdx(i:i+1,j,k:k+1))
               gradu(2,3,i,j,k) = gradu(2,3,i,j,k) + 0.25_WP*sum(dwdy(i,j:j+1,k:k+1))
            end do
         end do
      end do

      ! Apply a Neumann condition in non-periodic directions
      if (.not.cfg%xper) then
         if (cfg%iproc.eq.1)            gradu(:,:,cfg%imin-1,:,:)=gradu(:,:,cfg%imin,:,:)
         if (cfg%iproc.eq.cfg%npx) gradu(:,:,cfg%imax+1,:,:)=gradu(:,:,cfg%imax,:,:)
      end if
      if (.not.cfg%yper) then
         if (cfg%jproc.eq.1)            gradu(:,:,:,cfg%jmin-1,:)=gradu(:,:,:,cfg%jmin,:)
         if (cfg%jproc.eq.cfg%npy) gradu(:,:,:,cfg%jmax+1,:)=gradu(:,:,:,cfg%jmax,:)
      end if
      if (.not.cfg%zper) then
         if (cfg%kproc.eq.1)            gradu(:,:,:,:,cfg%kmin-1)=gradu(:,:,:,:,cfg%kmin)
         if (cfg%kproc.eq.cfg%npz) gradu(:,:,:,:,cfg%kmax+1)=gradu(:,:,:,:,cfg%kmax)
      end if

      ! Ensure zero in walls
      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
               if (fs%mask(i,j,k).eq.1) gradu(:,:,i,j,k)=0.0_WP
            end do
         end do
      end do

      gradUxx = gradU(1,1,:,:,:)

      ! Sync it
      call cfg%sync(gradu)
      call cfg%sync(gradUxx)

      ! Deallocate velocity gradient storage
      deallocate(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
   end subroutine applyExtraGradU

   subroutine getDiffutionCoeff
      integer :: i,j,k
      real(WP) :: diff, phi1
      vf%diff = 0.0_WP
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               phi1 = 1.0_WP-vf%SC(i,j,k)    ! phi1 is solvent volume fraction
               diff = diff0*exp(ad*phi1)
               vf%diff(i,j,k) = diff
            end do
         end do
      end do
      call cfg%sync(vf%diff)
   end subroutine getDiffutionCoeff

   subroutine getRelaxTime
      integer :: i,j,k
      real(WP) :: phi2
      relaxTime = ve%trelax
      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
               phi2 = vf%SC(i,j,k)  ! phi2 is polymer volume fraction
               relaxTime(i,j,k) = ve%trelax*(1.0+phi2*100_WP)
            end do
         end do
      end do
      call cfg%sync(relaxTime)
   end subroutine getRelaxTime

   subroutine getMixViscosity
      integer :: i,j,k
      real(WP) :: phi2, tau, tauRatio
      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
               phi2 = vf%SC(i,j,k)  ! phi2 is polymer volume fraction
               tau = relaxTime(i,j,k)
               tauRatio = (tau-minRelaxTime)/(maxRelaxTime-minRelaxTime)
               tauRatio = max(0.0_WP,min(1.0_WP,tauRatio))
               MixViscosity(i,j,k) = maxMixViscosity*phi2*tauRatio
            end do
         end do
      end do
      call cfg%sync(MixViscosity)
   end subroutine getMixViscosity

   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin .and. j.gt.cfg%jmin) isIn=.true.
   end function left_of_domain

   !> Function that localizes the left (x-) of the domain for scalar
   function left_of_domainsc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function left_of_domainsc

   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain

   !> Function that localizes the bottom (y-) of the domain
   function bottom_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function bottom_of_domain

   !> Function that localizes the bottom (y-) of the domain for scalar
   function bottom_of_domainsc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function bottom_of_domainsc

   !> Function that localizes the top (y+) of the domain
   function top_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function top_of_domain

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      use string, only: str_medium
      implicit none

      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=1
      end block initialize_timetracker

      ! Create a low Mach flow solver with bconds
      create_flow_solver: block
         use hypre_uns_class, only: gmres_amg
         use hypre_str_class, only: pcg_pfmg
         use lowmach_class,   only: dirichlet,clipped_neumann, neumann, slip
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Define bc
         call fs%add_bcond(name='inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='outflow',type=clipped_neumann,   locator=right_of_domain,face='x',dir=+1,canCorrect=.true.)
         call fs%add_bcond(name='topfreebd',type=slip      ,locator=top_of_domain,face='y',dir=+1,canCorrect=.false.)
         call fs%add_bcond(name='bottom',type=dirichlet      ,locator=bottom_of_domainsc,face='y',dir=-1,canCorrect=.false.)
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',rho); fs%rho=rho
         ! Configure pressure solver
         ps=hypre_uns(cfg=cfg,name='Pressure',method=gmres_amg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver

      ! Create a volume fraction scalar solver
      create_vf_solver: block
         use vdscalar_class, only: dirichlet,neumann,upwind,bquick
         ! Create volume fraction scalar solver
         vf=vdscalar(cfg=cfg,scheme=bquick,name='Volume fraction')
         ! Define bc
         call vf%add_bcond(name='inflow',type=dirichlet      ,locator=left_of_domainsc,dir='x-')
         call vf%add_bcond(name='outflow',type=neumann      ,locator=right_of_domain,dir='x+')
         call vf%add_bcond(name='bottom',type=dirichlet      ,locator=bottom_of_domainsc,dir='y-')
         ! Configure implicit scalar solver
         vfs = ddadi(cfg=cfg,name='Volume fraction',nst=13)
         ! Setup the solver
         call vf%setup(implicit_solver=vfs)
      end block create_vf_solver

      ! Create a viscoelastic solver
      create_viscoelastic: block
         use multiscalar_class,   only: bquick, upwind, neumann, dirichlet
         use viscoelastic_class,  only: oldroydb, fenep
         integer :: i,j,k
         ! Create FENE model solver
         call ve%init(cfg=cfg,model=fenep,scheme=upwind,name='FENE')
         ! Maximum extensibility of polymer chain
         call param_read('Maximum polymer extensibility',ve%Lmax)
         ! Relaxation time for polymer
         call param_read('Polymer relaxation time',ve%trelax)
         ! Polymer viscosity at zero strain rate
         call param_read('Polymer viscosity',ve%visc_p)
         ! Configure implicit scalar solver
         ves=ddadi(cfg=cfg,name='scalar',nst=13)
         ! Setup the solver
         call ve%setup(implicit_solver=ves)
         ! call ve%setup()
      end block create_viscoelastic

      prepareMixViscosity: block
         use param, only: param_read
         ! Read the maximum viscosity
         call param_read('Max Mixture viscosity',maxMixViscosity)

         ! calculate the relaxation time
         minRelaxTime = ve%trelax
         maxRelaxTime = ve%trelax*101.0_WP

      end block prepareMixViscosity


      ! Allocate work arrays
      allocate_work_arrays: block
         !< FLow solver
         allocate(resRHO  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi      (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resVF   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(vfbqflag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(vfTmp   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(SdiffU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SdiffV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SdiffW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(PdiffU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(PdiffV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(PdiffW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(SdiffUi (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)) ; SdiffUi = 0.0_WP
         allocate(SdiffVi (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)) ; SdiffVi = 0.0_WP
         allocate(SdiffWi (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)) ; SdiffWi = 0.0_WP
         allocate(PdiffUi (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)) ; PdiffUi = 0.0_WP
         allocate(PdiffVi (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)) ; PdiffVi = 0.0_WP
         allocate(PdiffWi (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)) ; PdiffWi = 0.0_WP

         allocate(relaxTime(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(resSC (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
         allocate(SCtmp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
         allocate(SR    (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(stress(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
         allocate(gradU (1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(gradUxx(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(MixViscosity(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

      end block allocate_work_arrays

      ! Initialize our velocity field
      initialize_velocity: block
         use lowmach_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Form momentum
         call fs%rho_multiply()


         ! Apply all other boundary conditions
         ! Read the shear velocity
         call param_read('inflow velocity',inflowVelocity)
         call fs%apply_bcond(time%t,time%dt)
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k) = inflowVelocity ; fs%rhoU(i,j,k) = rho*inflowVelocity
            fs%V(i-1,j,k) = 0.0_WP ; fs%rhoV(i-1,j,k) = rho*0.0_WP
         end do
         call fs%get_bcond('bottom',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k) = 0.0_WP ; fs%rhoU(i,j,k) = rho*0.0_WP
            fs%V(i,j,k) = 0.0_WP ; fs%rhoV(i,j,k) = rho*0.0_WP
         end do


         call fs%interp_vel(Ui,Vi,Wi)
         ! Here it should be changed to use mvdscalar
         resRHO=0.0_WP
         call fs%get_div(drhodt=resRHO)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block initialize_velocity

      initialize_volume_fraction: block
         use param, only: param_read
         use vdscalar_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: i,j,k,n

         vf%SC=0.0_WP

         ! read diffusivity
         vf%rho = rho
         call param_read('diffusion coefficient 0',diff0)
         call param_read('scaling factor',ad)
         call getDiffutionCoeff()

         call vf%apply_bcond(time%t,time%dt)
         ! Apply all other boundary conditions
         call vf%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            vf%SC(i,j,k) = 0.0_WP
         end do
         call vf%get_bcond('bottom',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            vf%SC(i,j,k) = 0.5_WP
         end do 

      end block initialize_volume_fraction

      ! Initialize viscoelastic solver
      initialize_viscoelastic: block
         use param, only: param_read
         use multiscalar_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Check first if we use stabilization
         call param_read('Stabilization',stabilization,default=.true.)

         if (stabilization) then
            !> Allocate storage fo eigenvalues and vectors
            allocate(ve%eigenval    (1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); ve%eigenval=0.0_WP
            allocate(ve%eigenvec(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); ve%eigenvec=0.0_WP
            allocate(ve%eigenval_log(1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); ve%eigenval_log=0.0_WP
            !> Allocate storage for reconstructured C
            allocate(ve%SCrec   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6)); ve%SCrec=0.0_WP
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     ve%SCrec(i,j,k,1)=1.0_WP  !< Cxx
                     ve%SCrec(i,j,k,4)=1.0_WP  !< Cyy
                     ve%SCrec(i,j,k,6)=1.0_WP  !< Czz
                  end do
               end do
            end do
            ! ! Apply all other boundary conditions
            ! call ve%get_bcond('rightwall',mybc)
            ! do n=1,mybc%itr%no_
            !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            !    ve%SCrec(i,j,k,1)=1.0_WP  !< Cxx
            !    ve%SCrec(i,j,k,4)=1.0_WP  !< Cyy
            !    ve%SCrec(i,j,k,6)=1.0_WP  !< Czz
            ! end do
            ! call ve%apply_bcond(time%t,time%dt)
            ! Get eigenvalues and eigenvectors
            call ve%get_eigensystem()
         end if

      end block initialize_viscoelastic

      ! Add Ensight output
      create_ensight: block
         integer :: i, nsc
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='flat_plate')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('PolymerDiffVelo',PdiffUi,PdiffVi,PdiffWi)
         call ens_out%add_vector('SolventDiffVelo',SdiffUi,SdiffVi,SdiffWi)
         call ens_out%add_scalar('pressure',fs%P)
         !call ens_out%add_scalar('diffusivity',vf%diff)
         call ens_out%add_scalar('MixViscosity',MixViscosity)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('volumefraction',vf%SC)
         call ens_out%add_scalar('gradUxx',gradUxx)
         call ens_out%add_scalar('relaxTime',relaxTime)
         ! call ens_out%add_vector("eigenvalues",ve%eigenval(1,:,:,:),ve%eigenval(2,:,:,:),ve%eigenval(3,:,:,:))
         ! call ens_out%add_vector("logeigenvalues",ve%eigenval_log(1,:,:,:),ve%eigenval_log(2,:,:,:),ve%eigenval_log(3,:,:,:))
         do nsc=1,ve%nscalar
            call ens_out%add_scalar(trim(ve%SCname(nsc)),ve%SCrec(:,:,:,nsc))
            ! call ens_out%add_scalar('ln'//trim(ve%SCname(nsc)),ve%SC(:,:,:,nsc))
            ! call ens_out%add_scalar('sctmp'//trim(ve%SCname(nsc)),SCtmp(:,:,:,nsc))
         end do
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      ! Create monitor filea
      create_monitor: block
         integer :: nsc
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         if (stabilization) then
            call ve%get_max_reconstructed()
         end if
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create scalar monitor
         scfile=monitor(ve%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         if (stabilization) then
            do nsc=1,ve%nscalar
               call scfile%add_column(ve%SCrecmin(nsc),trim(ve%SCname(nsc))//'_min')
               call scfile%add_column(ve%SCrecmax(nsc),trim(ve%SCname(nsc))//'_max')
            end do
         end if
         call scfile%write()
         ! Create Diffusion monitor
         Difffile=monitor(vf%cfg%amRoot,'Diffusion')
         call Difffile%add_column(time%n,'Timestep number')
         call Difffile%add_column(time%t,'Time')
         call Difffile%write()
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      use mathtools, only: twoPi
      use parallel, only: parallel_time

      implicit none
      real(WP) :: cfl
      integer :: isc,i,j,k,n

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old density, velocity, and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         vf%rhoold=vf%rho
         vf%SCold = vf%SC

         ! Remember old scalars
         ve%SCold=ve%SC

         ! ======================= conformation tensor solver ===========================

         ! Calculate grad(U)
         call fs%get_gradu(gradu)
         call applyExtraGradU()

         ! Transport our liquid conformation tensor using log conformation
         advance_scalar: block
            use, intrinsic :: ieee_arithmetic
            integer :: ierr
            integer :: i,j,k,nsc
            ! Add streching source term for constitutive model
            call ve%get_CgradU_log(gradU,SCtmp); resSC=SCtmp
            ve%SC=ve%SC+time%dt*resSC
            ve%SCold=ve%SC
            ! Explicit calculation of dSC/dt from scalar equation
            PdiffU = PdiffU + fs%Uold
            PdiffV = PdiffV + fs%Vold
            PdiffW = PdiffW + fs%Wold
            call ve%get_drhoSCdt(resSC,PdiffU,PdiffV,PdiffW)
            call ve%get_drhoSCdt(drhoSCdt=resSC,rhoU=fs%U,rhoV=fs%V,rhoW=fs%W)
            ve%SC = ve%SCold + time%dt*resSC
         end block advance_scalar


         ! Add in relaxation forcing and reconstruct C
         if (stabilization) then
            ! Get eigenvalues and eigenvectors from lnC
            call ve%get_eigensystem()
            ! Reconstruct conformation tensor from eigenvalues and eigenvectors
            call ve%reconstruct_conformation()
            ! Add in relaxtion source from semi-anlaytical integration
            call getRelaxTime()
            call ve%get_relax_analytical_varyRelxTime(relaxTime,time%dt)
            ! call ve%get_relax_analytical(time%dt)
            ! Reconstruct lnC for next time step
            !> get eigenvalues and eigenvectors based on reconstructed C
            call ve%get_eigensystem_SCrec()
            !> Reconstruct lnC from eigenvalues and eigenvectors
            call ve%reconstruct_log_conformation()
         end if


         ! ==============================================================================

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Update density based on particle volume fraction and multi-vd
            fs%rho=0.5_WP*(fs%rho+fs%rhoold)

            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

            ! ==================== Volume Fraction Solver ==================

            call getDiffutionCoeff()

            call vf%metric_reset()

            ! Buid mid-time volume fraction
            vf%SC=0.5_WP*(vf%SC+vf%SCold)

            ! Assembly of explicit residual
            call vf%get_drhoSCdt(resVF,fs%rhoU,fs%rhoV,fs%rhoW)
            resVF = time%dt*resVF - (2.0_WP*vf%rho*vf%SC - (vf%rho+vf%rhoold)*vf%SCold)

            !< Get temperary solution for bquick
            vfTmp = 2.0_WP*vf%SC - vf%SCold + resVF/vf%rho

            vfbqflag = .false.
            do k=vf%cfg%kmino_,vf%cfg%kmaxo_
               do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                  do i=vf%cfg%imino_,vf%cfg%imaxo_
                     if ((vfTmp(i,j,k).lt.0.0_WP).or.(vfTmp(i,j,k).gt.1.0_WP)) then
                        vfbqflag(i,j,k) = .true.
                     end if
                  end do
               end do
            end do
            ! all bquick to avoid oscillations
            vfbqflag = .true.

            ! Adjust metrics
            call vf%metric_adjust(vfTmp, vfbqflag)

            ! re-Assemble explicit residual
            call vf%get_drhoSCdt(resVF,fs%rhoU,fs%rhoV,fs%rhoW)
            resVF = time%dt*resVF - (2.0_WP*vf%rho*vf%SC - (vf%rho+vf%rhoold)*vf%SCold)

            ! solve implicitlly
            call vf%solve_implicit(time%dt,resVF,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Apply these residuals
            vf%SC = 2.0_WP*vf%SC - vf%SCold + resVF

            ! Apply boundary conditions
            call vf%apply_bcond(time%t,time%dt)
            vf_bc : block
               use vdscalar_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               ! Apply all other boundary conditions
               call vf%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  vf%SC(i,j,k) = 0.0_WP
               end do
               call vf%get_bcond('bottom',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  vf%SC(i,j,k) = 0.5_WP
               end do
            end block vf_bc

            ! ===================== Velocity Solver =======================

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)


            ! ! ! Add polymer stress term
            polymer_stress: block
               use viscoelastic_class, only: fenep,lptt,eptt, oldroydb
               integer :: i,j,k,n
               real(WP), dimension(:,:,:), allocatable :: Txy,Tyz,Tzx
               real(WP) :: coeff
               real(WP) :: mixvisc

               ! Get the viscosity of the mixture
               call getMixViscosity()

               ! Allocate work arrays
               allocate(Txy   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(Tyz   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(Tzx   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               ! Build liquid stress tensor
               select case (ve%model)
                case (fenep, oldroydb)
                  call ve%get_relax_varyRelaxTime(relaxTime,stress,time%dt)
                  do k=cfg%kmino_,cfg%kmaxo_
                     do j=cfg%jmino_,cfg%jmaxo_
                        do i=cfg%imino_,cfg%imaxo_
                           if (ve%mask(i,j,k).ne.0) cycle                !< Skip non-solved cells
                           mixvisc = MixViscosity(i,j,k)
                           do n=1,6
                              stress(i,j,k,n)=-mixvisc*stress(i,j,k,n)
                           end do
                        end do
                     end do
                  end do
                  ! do n=1,6
                  !    stress(:,:,:,n)=-ve%visc_p*stress(:,:,:,n)
                  ! end do
                case (eptt,lptt)
                  stress=0.0_WP
                  coeff=ve%visc_p/(ve%trelax*(1-ve%affinecoeff))
                  do n=1,6
                     do k=cfg%kmino_,cfg%kmaxo_
                        do j=cfg%jmino_,cfg%jmaxo_
                           do i=cfg%imino_,cfg%imaxo_
                              if (ve%mask(i,j,k).ne.0) cycle                !< Skip non-solved cells
                              stress(i,j,k,1)=coeff*(ve%SC(i,j,k,1)-1.0_WP) !> xx tensor component
                              stress(i,j,k,2)=coeff*(ve%SC(i,j,k,2)-0.0_WP) !> xy tensor component
                              stress(i,j,k,3)=coeff*(ve%SC(i,j,k,3)-0.0_WP) !> xz tensor component
                              stress(i,j,k,4)=coeff*(ve%SC(i,j,k,4)-1.0_WP) !> yy tensor component
                              stress(i,j,k,5)=coeff*(ve%SC(i,j,k,5)-0.0_WP) !> yz tensor component
                              stress(i,j,k,6)=coeff*(ve%SC(i,j,k,6)-1.0_WP) !> zz tensor component
                           end do
                        end do
                     end do
                  end do
               end select
               ! Interpolate tensor components to cell edges
               do k=cfg%kmin_,cfg%kmax_+1
                  do j=cfg%jmin_,cfg%jmax_+1
                     do i=cfg%imin_,cfg%imax_+1
                        Txy(i,j,k)=sum(fs%itp_xy(:,:,i,j,k)*stress(i-1:i,j-1:j,k,2))
                        Tyz(i,j,k)=sum(fs%itp_yz(:,:,i,j,k)*stress(i,j-1:j,k-1:k,5))
                        Tzx(i,j,k)=sum(fs%itp_xz(:,:,i,j,k)*stress(i-1:i,j,k-1:k,3))
                     end do
                  end do
               end do
               ! Add divergence of stress to residual
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+(sum(fs%divu_x(:,i,j,k)*stress(i-1:i,j,k,1))&
                        &                                                 +sum(fs%divu_y(:,i,j,k)*Txy(i,j:j+1,k))     &
                        &                                                 +sum(fs%divu_z(:,i,j,k)*Tzx(i,j,k:k+1)))*time%dt
                        if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+(sum(fs%divv_x(:,i,j,k)*Txy(i:i+1,j,k))     &
                        &                                                 +sum(fs%divv_y(:,i,j,k)*stress(i,j-1:j,k,4))&
                        &                                                 +sum(fs%divv_z(:,i,j,k)*Tyz(i,j,k:k+1)))*time%dt
                        if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+(sum(fs%divw_x(:,i,j,k)*Tzx(i:i+1,j,k))     &
                        &                                                 +sum(fs%divw_y(:,i,j,k)*Tyz(i,j:j+1,k))     &
                        &                                                 +sum(fs%divw_z(:,i,j,k)*stress(i,j,k-1:k,6)))*time%dt
                     end do
                  end do
               end do
               ! Clean up
               deallocate(Txy,Tyz,Tzx)
            end block polymer_stress


            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions and update momentum
            call fs%rho_multiply()

            call fs%apply_bcond(time%t,time%dt)
            fs_bc : block
               use lowmach_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               ! Apply all other boundary conditions
               call fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k) = inflowVelocity ; fs%rhoU(i,j,k) = rho*inflowVelocity
                  fs%V(i-1,j,k) = 0.0_WP ; fs%rhoV(i-1,j,k) = rho*0.0_WP
               end do
               call fs%get_bcond('bottom',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k) = 0.0_WP ; fs%rhoU(i,j,k) = rho*0.0_WP
                  fs%V(i,j,k) = 0.0_WP ; fs%rhoV(i,j,k) = rho*0.0_WP
               end do
            end block fs_bc

            ! Solve Poisson equation
            call fs%correct_mfr(drhodt=resRHO)
            call fs%get_div(drhodt=resRHO)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide()

            ! ==============================================================

            ! Increment sub-iteration counter

            time%it=time%it+1

         end do

         ! ==================== Post Process ==================================

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(drhodt=resRHO)

         call computeDiffU()

         ! Output to ensight
         if (ens_evt%occurs()) then
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         if (stabilization) then
            call ve%get_max_reconstructed()
         end if
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         call Difffile%write()
      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker

      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)

   end subroutine simulation_final

end module simulation

