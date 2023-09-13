module dust_update_mod
! module to calculate moment dust field updates

use constants_mod, only: KAPPA, CP_AIR, RDGAS, GRAV, PI, RADIAN, seconds_per_day  

use fms_mod, only: error_mesg, FATAL,                      &
                   open_namelist_file, check_nml_error,                &
                   mpp_pe, mpp_root_pe, close_file,                    &
                   write_version_number, stdlog,                       &
                   uppercase, read_data, write_data, field_size

use fms2_io_mod,            only:  file_exists

use   mpp_domains_mod, only: domain2d

use        fms_io_mod, only: register_restart_field, restart_file_type, &
                             save_restart, restore_state, get_mosaic_tile_file

use time_manager_mod, only: time_type, get_time
use diag_manager_mod, only: register_diag_field, send_data,  diag_axis_init

use astronomy_mod,  only  :   mars_calender
#ifndef RELEASE
use tagging_method_mod,  only  : tagging_main
#endif

use mars_surface_mod,  only:  sfc_roughness, sfc_topo 
use initracer_mod
use aerosol_util_mod, only:  dust_map_scale, Reff_backgd, Reff_stress, Reff_dd

implicit none
private

!------------------- interfaces ---------------------------------------

public ::  dust_update_init,dust_update,dust_update_end,&
           sfc_dust_mass,tcol_mass,tcol_core,dust_surf_ini


! sfc_dust_mass is a globally-dimensioned array for evolving surface dust_mom (moment scheme)
!-------------------- namelist -----------------------------------------

logical,public :: firstcall=.true.              !  first call for dust_update_init
logical :: interact=.false.                     !  fully interactive dust lifting
logical :: stress_lift=.true.                  !  stress lifting
logical :: Background=.true.                    !  background dust lifting
logical :: opac_from_aerosol=.false.            !  background dust lifting from aerosol routine
logical :: convective_lift=.true.              !  dust devil lifting
logical :: DDA=.true.                          !  newman dust devil lifting with sensible heat flux threshold
logical :: DTH=.false.                          !  newman dust devil lifting with tangential wind threshold
logical :: devils_fv3=.false.                   !  fv3 style dust devil lifting
logical :: inject_pbl_stress=.false.            !  inject stress lifted dust into pbl
logical :: inject_pbl_dd=.true.                !  inject dust devil lifted dust into pbl
real    ::  alfa = 0.0065                         !  scaling factor for stress lifting
real    ::  threshold_stress=0.022              !  threshold for stress lifting
real    ::  thres_shflx=0.                      !  threshold for sensible heat flux
real    ::  ws_co2_ice_thresh_N=1.              !  threshold for shutting off stress lifting over north co2 ice
real    ::  ws_co2_ice_thresh_S=1.              !  threshold for shutting off stress lifting over south co2 ice
real    ::  ws_h2o_ice_thresh=1.                !  threshold for shutting off stress lifting over water ice
real    ::  dd_co2_ice_thresh_N=1.              !  threshold for shutting off dust devil lifting over north co2 ice
real    ::  dd_co2_ice_thresh_S=1.              !  threshold for shutting off dust devil lifting over south co2 ice
real    ::  dd_h2o_ice_thresh=1.                !  threshold for shutting off dust devil lifting over water ice
real    ::  optd_thresh=100.                    !  threshold opacity for shutting off all lifting
integer    ::  stress_scheme=1                  !  stress lifting scheme select. 1 and default: legacy formulation; 2: John; 3: John with a3=1; 4: White 1979; 5: Shao 1993
real    ::  assim_t_thresh = 120.0              !  threshold surface temperature for shutting off assimilated dust
logical ::  no_assim_over_caps = .true.        !  shut off assimilated dust over co2 ice caps
real    ::  conv_delt_threshold = 22.           !  threshold tsfc-tatm difference for dust devils
real    ::  alpha_dth = 0.5                     !  scaling for dth lifting
real    ::  alpha_dda = 1.4e-10                     !  scaling for dda lifting
real    ::  source_ddfv3 = 0.3                  !  rough conversion from opacity to mixing ratio
!real    ::  Qext_dt = 2.938                     !  for Background dust to be injected
real    ::  Dlift = 6.e-6                       !  Diameter of lifted dust particles for dust devils lifting DTH.
real    ::  reff_ddfv3=2.5e-6                   !  Particle radius for dust devils mode fv3 
real    ::  tauscale=1.                         !  factor for the dust scenario opacities
real    ::  injscale=1.                         !  Scaling the source with the dust scenario
real    ::  ltscale=0.                          !  Scaling the source with the dust scenario local time
real    ::  sinkscale=1.                        !  Scaling the sink with the dust scenario
real    ::  dust_surf_ini=30.                   !  Initial reservoir of dust for cold cases
real    ::  limres_xfac=20.                     !  Factor for limited reservoir stress threshold calculation
real    ::  limres_dec=1.e-11                       !  Linear decrease of threshold for refilling dust into the system
real    ::  delta_thres=0.2                      !  critical increase of wind stress threshold from which dust devil lifting is cut
integer, dimension(2)   :: kfix = (/ 0, 0 /)    !  levelsfrom which dust in injected in the atmosphere from the surface
logical ::  inject_pbl_bd=.true.                !  Inject dust in PBL for Background mode
logical ::  sink_bd=.false.                      !  Sink dust in entire column for Background mode
logical ::  limited_reservoir=.false.           !  Compute flux according to reservoir on surface
integer ::  dual_mode=0                         !  Run with free dust mode but with opacity constraint from dust scenario :  0 : F / 1 : mode1 / 2 : mode2
logical ::  segment=.false.                     !  segment dust scenario instead of interpolation
logical ::  bottom_sink=.false.                 !  Remove dust from the bottom (near surface)
real    ::  dual_scale=1.0                      !  factor for scaling opacity in dual mode
real    ::  dtaulim=-999.                       !  critical Dtau to allow injection
integer    ::  dgdm_type=1                      !  Type of daily global dust map (1=VIS, 2=IR)
integer ::  injec_type=1                        !  Calculation method for injection factor: 1 = constant factor, 2 = grid avg
logical ::  lifting_nosed=.false.               !  Flag to turn off sedimentation if dust was lifted
integer :: ndust_mass_rst= 0                    !  restart dust mass number
namelist /dust_update_nml/ ws_co2_ice_thresh_N,ws_co2_ice_thresh_S,threshold_stress, alfa, dd_co2_ice_thresh_N,dd_co2_ice_thresh_S, &
                            interact,stress_lift,Background,optd_thresh,stress_scheme, ws_h2o_ice_thresh,dd_h2o_ice_thresh, &
                            assim_t_thresh,no_assim_over_caps,convective_lift,  &
                            conv_delt_threshold,alpha_dth,alpha_dda,devils_fv3,DDA,DTH, &
                            source_ddfv3,inject_pbl_dd,inject_pbl_stress, &
                            !Qext_dt,
                            Dlift,reff_ddfv3,tauscale,inject_pbl_bd,dust_surf_ini,injscale,sinkscale,sink_bd,limres_xfac,limited_reservoir, &
                            limres_dec,thres_shflx,delta_thres,dual_mode,dual_scale,segment,dtaulim,ltscale,kfix,bottom_sink,dgdm_type, &
                            injec_type,lifting_nosed,opac_from_aerosol


!-------------------- Other -----------------------------------------
logical ::  mcpu0
real, dimension(:,:,:), allocatable :: tcol_mass,tcol_core
real, dimension(:,:,:),   allocatable, save  ::  sfc_dust_mass
integer, dimension(:),  allocatable  ::   id_dmcol,  id_dust_mass_src, id_dust_mass_accum,id_dust_mass_sink,id_dmtau,id_dncol,id_rn,id_dmtau2
integer, dimension(:),  allocatable  ::   id_dust_mass_src_ws,id_dust_mass_src_dd,id_dust_mass_src_bg
integer, dimension(:),  allocatable  ::   id_dccol, id_rm2d, id_inject, id_taufrac2d
integer  ::  id_threshold_stress,id_src_theo_dd,id_src_theo_ws,id_tauscenario,id_dtauscenario
integer  ::  id_tauscenario_ini,id_taucurrent,id_ltfrac,id_tau_old,id_tauadd
real, parameter   :: missing_value = -1.e10
character(len=11) :: mod_name = 'dust_source'
logical :: module_is_initialized = .false.

real, dimension(:),     allocatable :: areo_dust_scenario
real, dimension(:,:,:), allocatable :: tau_scenario_ls

logical :: custom_sources = .false.
real, dimension(:),     allocatable :: areo_sources
real, dimension(:,:,:), allocatable :: sources_scenario_ls

logical :: ltfrac = .false.
real, dimension(:),     allocatable :: areo_ltfrac
integer, dimension(:),     allocatable :: bins_ltfrac
real, dimension(:,:,:,:), allocatable :: ltfrac_scenario_ls
integer  ::  bins_length
real, dimension(:,:), allocatable :: tau_ltfrac_old
integer, dimension(:,:), allocatable :: bin_ltfrac_old

logical :: changeradius = .false.
real, dimension(:),     allocatable :: areo_radius
integer, dimension(:),     allocatable :: bins_radius
real, dimension(:,:,:,:), allocatable :: radius_scenario_ls
integer  ::  radbins_length

logical :: inidust = .false.
real, dimension(:,:), allocatable :: inidust_scenario

!--- for restart file
type(restart_file_type), pointer, save :: Dst_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()   !needed for tile restarts
logical                                :: in_different_file = .false.
integer,  dimension(:),  allocatable   ::   dust_mass_indx2

! ******************* NOTES *******************

! The way and options to rad transfer : to be improved
! opt_dst and opt_cld could merge in one subroutine
! nnb must be taken into account
! Add limited surface dust reservoirs
! Get taudust VIS and IR in output
! Put flag to check that field table is ok with options : otherwise big bugs in micro
! qextref(:) = 2.0 : used to compute tau_new (taudst_mom). This is an estimate but not the real opacity which is computed in the Rad/opt_dst/ routine. Used for all options here.
! if dust is active (flag active_dust=true) then the dust profile from the dust_nma micro is used for the radiative transfer

contains

!********************************************************************
!********************************************************************
!********************************************************************

subroutine dust_update (   is, js, lon, lat, dt, Time, &
                              p_half, p_full, p_pbl, k_pbl, tsfc,  snow, &
                              frost, stress, taudust, taudust_fix, t, tdt, r, rdt, shflx, &
                              rdt_dst, source_mom, lifting_dust ) 

use field_manager_mod, only  : MODEL_ATMOS, find_field_index      

integer, intent(in)  :: is, js
real,    intent(in)  :: dt
type(time_type), intent(in)             :: Time
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:)     :: lat
real, intent(in),    dimension(:,:,:)   :: p_half, p_full
real, intent(in),    dimension(:,:)     :: p_pbl
integer, intent(in),    dimension(:,:)  :: k_pbl
real, intent(in),    dimension(:,:)     :: tsfc
real, intent(in),    dimension(:,:)     :: snow
real, intent(in),    dimension(:,:)     :: frost
real, intent(in),    dimension(:,:)     :: stress    
real, intent(in),    dimension(:,:,:)   :: t,tdt
real, intent(in),    dimension(:,:,:)   :: taudust
real, intent(in),    dimension(:,:,:)   :: taudust_fix
real, intent(in),    dimension(:,:,:,:) :: r,rdt
real, intent(in),    dimension(:,:)     :: shflx  
real, intent(out), dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_dst
real, intent(out),    dimension(:,:)     :: source_mom  
logical, intent(out), dimension(size(r,1),size(r,2),size(r,4)) :: lifting_dust

! -----Local  variables ----------------------
integer  :: ie, je, id, jd, kd, i, j, k, l, nt,ndx
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rnew,rini
real, dimension(size(r,1),size(r,2),ndust_mass) :: dmcol,source_mass,Nadd,Madd,Nrem,Mrem,dmtau,dccol,dncol
real, dimension(size(r,1),size(r,2),ndust_mass) :: flux,tau_new,source_dd,source_ws,source_bg,fluxrem,tau_new2,injec,sink
real, dimension(size(t,1),size(t,2),size(t,3)) :: delp, rho,tnew
real, dimension(size(t,1),size(t,2),size(t,3),ndust_mass) :: Rn
real, dimension(size(r,1),size(r,2)) :: rhokd, tair, threshold, tau_current,fluxtot
real, dimension(size(r,1),size(r,2)) :: vtan_t,eta,eta_H,dp_lift,fluxfv3,fluxdda,fluxdth,&
                                        vtan,ratio_sinkN,ratio_sinkM,Mtot,Ntot
real, dimension(size(r,1),size(r,2)) :: mass1l,b,masspbl,masspblfix,Mrem_col,Nrem_col,Mtmp,Ntmp,Mtot_col,Ntot_col
real, dimension(size(r,1),size(r,2)) :: flux_theo_dd,flux_theo_ws 
real, dimension(size(r,1),size(r,2)) :: cst2
real, dimension(size(r,1),size(r,2))         :: tau_scenario,sources_scenario   ! optical depth from dust scenario dust_cycle.nc
real, dimension(size(r,1),size(r,2),bins_length-1)         :: ltfrac_scenario   ! scenario for dust injection depending on local time
real, dimension(size(r,1),size(r,2),radbins_length)         :: radius_scenario   ! scenario for lifted dust particle radius
real, dimension(size(r,1),size(r,2))         :: ltfrac_curr,localtime,dtau,dtau_current  !fraction injection local time, array of current local time
real, dimension(size(r,1),size(r,2))         :: tau_add   ! difference between optical depth from dust scenario dust_cycle.nc and current one
integer, dimension(size(r,1),size(r,2))         :: bin_curr   ! current array of bin (for ltfrac option)
real, dimension(size(r,1),size(r,2),ndust_mass)         :: Rm2D,Rss2D,tau_frac2D  ! radius and taufrac of lifted paticle size 2D 
real, dimension(size(r,3)) :: qextref 
integer :: days, seconds, iaa, iaam 
real    :: areolat, fac1, delareo,secs,sols,fjd,r_orbit,declin         
real :: dev,dens
real :: cst
real :: Rm,mass_part
real :: xsi,gama    !,p_pbl
real tauref_ddfv3
real, dimension(size(r,1),size(r,2)) :: tautmp,injscale2D
real pbar_ddfv3
real rhodust_ddfv3 
real qext_ddfv3
real reffec_ddfv3     !  Effective particle radius
real source_scaling
real ddfv3_source
logical :: used
integer :: nma_dst,nnb_dst

mcpu0 = (mpp_pe() == mpp_root_pe())
id= size(t,1); jd= size(t,2); kd= size(t,3) 
ie= is + id - 1 
je= js + jd - 1 

!if (mpp_pe() == mpp_root_pe()) print*,

!------------------- Initialisation ----------------------

qextref(:) = 2.
flux(:,:,:)=0. 
fluxrem(:,:,:)=0. 
fluxdda(:,:)=0. 
fluxdth(:,:)=0. 
fluxfv3(:,:)=0. 
tau_new(:,:,:) = 0.
tau_new2(:,:,:) = 0.
dmcol(:,:,:) = 0.
dncol(:,:,:) = 0.
dccol(:,:,:) = 0.
flux_theo_dd(:,:)=0.
flux_theo_ws(:,:)=0.
source_mom(:,:)=0.
source_dd(:,:,:)=0.
source_ws(:,:,:)=0.
source_bg(:,:,:)=0.
lifting_dust(:,:,:)=.false.

!-------------- Get time info and areolat-----------------
call get_time(Time, seconds, days)
secs= days * seconds_per_day + seconds
sols = secs / seconds_per_day
days = sols
fjd = sols - days    !sol=0 at lon=180
call mars_calender( days, fjd, r_orbit, declin, areolat )
areolat = areolat * RADIAN
if( areolat < 0.0 )  areolat = areolat + 360.0
areolat = modulo(areolat,360.)

!---------------- Pressure each layer --------------------
do k=1,kd
    delp(:,:,k)= p_half(:,:,k+1)-p_half(:,:,k)
enddo 
!---------------- Particles properties -------------------
dens    = dpden_dt  ! density dust 
dev     = dev_dt    ! Standard deviation of the dust distribution
cst     = .75 / (pi*dens) * exp( -4.5*dev**2. )   ! fixed distribution

!------------- Mass each layer and PBL -------------------
masspbl(:,:)=0.0
masspblfix(:,:)=0.0
do i=is,ie
    do j=js,je
        do k= k_pbl(i,j), kd            ! from namelist
            masspbl(i,j)= masspbl(i,j) + delp(i,j,k)/grav
        enddo
        do k= kd-kfix(2), kd-kfix(1)            ! from namelist
            masspblfix(i,j)= masspblfix(i,j) + delp(i,j,k)/grav
        enddo
    enddo
enddo

!------------------ Update tracers -----------------------
rini(:,:,:,:)= 0. 
rnew(:,:,:,:)= 0. 
rdt_dst(:,:,:,:)= 0. 
rnew(is:ie,js:je,:,:)=r(is:ie,js:je,:,:)+rdt(is:ie,js:je,:,:)*dt
rini(is:ie,js:je,:,:)=r(is:ie,js:je,:,:)+rdt(is:ie,js:je,:,:)*dt
do nt= 1, ndust_mass
    ndx= dust_mass_indx(nt)
    rnew(:,:,:,ndx)=max(rnew(:,:,:,ndx), 1.e-31 )
    rnew(:,:,:,ndx+1)=max(rnew(:,:,:,ndx+1), 1.e-31 )
    rini(:,:,:,ndx)=max(rini(:,:,:,ndx), 1.e-31 )
    rini(:,:,:,ndx+1)=max(rini(:,:,:,ndx+1), 1.e-31 )
enddo

! Scale current opacity to a reference pressure of 6.1 mbar
if (opac_from_aerosol) then
    tau_current(:,:) = taudust(:,:,dgdm_type)
else
    tau_current(:,:) = taudust(:,:,dgdm_type) / p_full(:,:,kd) * 610.
endif
if( id_taucurrent > 0 )  used = send_data( id_taucurrent, tau_current(:,:), Time, is, js )
#ifndef RELEASE
!------------------- Tagging methods --------------------
nma_dst= find_field_index( MODEL_ATMOS, 'dst_mass_mom' )
nnb_dst= find_field_index( MODEL_ATMOS, 'dst_num_mom' )
do nt= 1, ndust_mass
    ndx= dust_mass_indx(nt)
    call tagging_main( (/ "envir" /) ,lat,lon,ndx,field3d=rnew(:,:,:,ndx),taudust_curr=tau_current(:,:),dustfield=rnew(:,:,:,nma_dst), pres=p_half(:,:,:))
    call tagging_main( (/ "envir" /) ,lat,lon,ndx+1,field3d=rnew(:,:,:,ndx+1),taudust_curr=tau_current(:,:),dustfield=rnew(:,:,:,nnb_dst), pres=p_half(:,:,:))
    call tagging_main( (/ "activatestorm" /) ,lat,lon,ndx,field3d=rnew(:,:,:,ndx),taudust_curr=tau_current(:,:),dustfield=rnew(:,:,:,nma_dst),pres=p_half(:,:,:),solrun=sols,stress=stress)
    call tagging_main( (/ "activatestorm" /) ,lat,lon,ndx+1,field3d=rnew(:,:,:,ndx+1),taudust_curr=tau_current(:,:),dustfield=rnew(:,:,:,nnb_dst),pres=p_half(:,:,:),solrun=sols,stress=stress)
enddo
#endif

!------------------ Update temperature ------------------
tnew(:,:,:)=t(:,:,:)+tdt(:,:,:)*dt
!------------------ Update density ------------------
do k=1,kd
    rho(:,:,k) = p_full(:,:,k) / ( rdgas*tnew(:,:,k) )
enddo 

!------------------ Update Other Data -------------------

! Get tau scenario
if (background.or.custom_sources.or.dual_mode.gt.0) then
    if (.not.opac_from_aerosol) then
        ! get index corresponding to current ls : iaa,iaam
        do iaa= 1, size(areo_dust_scenario)
            if( areo_dust_scenario(iaa) > areolat )  then
                exit
            else
                cycle
            endif
        enddo
        iaam= iaa - 1
        if( iaam < 1 )  iaam= iaam + size(areo_dust_scenario)
        ! get delta of ls 
        delareo= areo_dust_scenario(iaa)-areo_dust_scenario(iaam)
        if( delareo  < 0. )       delareo = delareo + 360.
        fac1= (areo_dust_scenario(iaa)- areolat)/delareo
        tau_scenario(:,:)= fac1*tau_scenario_ls(is:ie,js:je,iaam) + (1.0-fac1)*tau_scenario_ls(is:ie,js:je,iaa)
        if (segment) tau_scenario(:,:) = tau_scenario_ls(is:ie,js:je,iaa)
        ! difference tau between the two current sols
        dtau(:,:)=tau_scenario_ls(is:ie,js:je,iaa)-tau_scenario_ls(is:ie,js:je,iaam)
    else
        tau_scenario(:,:) = taudust_fix(:,:,dgdm_type)
        dtau(:,:)=0.0
    endif
    ! output
    if( id_tauscenario_ini > 0 )  used = send_data( id_tauscenario_ini, tau_scenario(:,:), time, is, js )
endif

if (dual_mode.eq.1.or.dual_mode.eq.21) then
    dtau_current(:,:)=dual_scale*tau_scenario(:,:)-tau_current(:,:)  
endif

!****** data for Local Time depending lifting *******
if (ltfrac) then
    !!! get ltfrac scenario, interpolated at current ls
    ! get index corresponding to current ls : iaa,iaam
    do iaa= 1, size(areo_ltfrac)
        if( areo_ltfrac(iaa) > areolat )  then
            exit
        else
            cycle
        endif
    enddo
    iaam= iaa - 1
    if( iaam < 1 )  iaam= iaam + size(areo_ltfrac)
    ! get delta of ls 
    delareo= areo_ltfrac(iaa)-areo_ltfrac(iaam)
    if( delareo  < 0. )       delareo = delareo + 360.
    ! interpolate current ls with dust scenario and get optical depth scenario
    fac1= (areo_ltfrac(iaa)- areolat)/delareo
    ltfrac_scenario(:,:,:)= fac1*ltfrac_scenario_ls(is:ie,js:je,iaam,:) + (1.0-fac1)*ltfrac_scenario_ls(is:ie,js:je,iaa,:)

    !!! get current local time array / current bins array / current fraction array
    ! get current local time
    localtime(:,:)=modulo( modulo(sols+0.5,1.)*24.-(180.-modulo(lon(:,:)*180./pi,360.))*12./180., 24.)
    ! get current fractions and bins
    do i=1,bins_length-1
        where(localtime(:,:).lt.bins_ltfrac(i+1).and.localtime(:,:).ge.bins_ltfrac(i))
            ltfrac_curr(:,:)=ltfrac_scenario(:,:,i)
            bin_curr(:,:)=i
        end where
    enddo

    !!! Compute current tau
    if (opac_from_aerosol) then
        where (dtau(:,:).gt.0.)
            tau_scenario(:,:)= ltfrac_curr(:,:)*dtau(:,:)+tau_ltfrac_old(:,:)
            tau_scenario(:,:)=min(tau_scenario(:,:),taudust_fix(:,:,dgdm_type))
            tau_scenario(:,:)=max(tau_scenario(:,:),taudust_fix(:,:,dgdm_type))
        end where
        where (dtau(:,:).le.0.)
            tau_scenario(:,:)=taudust_fix(:,:,dgdm_type)
        end where
    else
        where (dtau(:,:).gt.0.) 
            tau_scenario(:,:)= ltfrac_curr(:,:)*dtau(:,:)+tau_ltfrac_old(:,:)
            tau_scenario(:,:)=min(tau_scenario(:,:),tau_scenario_ls(is:ie,js:je,iaa))
            tau_scenario(:,:)=max(tau_scenario(:,:),tau_scenario_ls(is:ie,js:je,iaam))
        end where
        where (dtau(:,:).le.0.) 
            tau_scenario(:,:)=tau_scenario_ls(is:ie,js:je,iaa)
        end where
    end if

    !!! check if changes of bins : if change, then tau_old and bin updated  !! should be before
    where(bin_curr(:,:).ne.bin_ltfrac_old(:,:))
        bin_ltfrac_old(:,:)=bin_curr(:,:)
        tau_ltfrac_old(:,:)=tau_scenario(:,:)
    end where
    if( id_ltfrac > 0 )  used = send_data( id_ltfrac, ltfrac_curr(:,:), time, is, js )
    if( id_tau_old > 0 )  used = send_data( id_tau_old, tau_ltfrac_old(:,:), time, is, js )
endif

!************ Custom sources ****************
if (custom_sources) then
    ! get index corresponding to current ls : iaa,iaam
    do iaa= 1, size(areo_sources)
        if( areo_sources(iaa) > areolat )  then
            exit
        else
            cycle
        endif
    enddo
    iaam= iaa - 1
    if( iaam < 1 )  iaam= iaam + size(areo_sources)
    ! get delta of ls 
    delareo= areo_sources(iaa)-areo_sources(iaam)
    if( delareo  < 0. )       delareo = delareo + 360.
    ! interpolate current ls with dust scenario and get optical depth scenario
    fac1= (areo_sources(iaa)- areolat)/delareo
    sources_scenario(:,:)= fac1*sources_scenario_ls(is:ie,js:je,iaam) + (1.0-fac1)*sources_scenario_ls(is:ie,js:je,iaa)
endif

!****** data for time-dependant lifted particle size *******
if (changeradius) then
    !!! get radius_scenario, interpolated at current ls
    ! for each bin (different dust type) 
    do iaa= 1, size(areo_radius)
        if( areo_radius(iaa) > areolat )  then
            exit
        else
            cycle
        endif
    enddo
    iaam= iaa - 1
    if( iaam < 1 )  iaam= iaam + size(areo_radius)
    ! get delta of ls 
    delareo= areo_radius(iaa)-areo_radius(iaam)
    if( delareo  < 0. )       delareo = delareo + 360.
    ! interpolate current ls with dust scenario and get optical depth scenario
    fac1= (areo_radius(iaa)- areolat)/delareo
    radius_scenario(:,:,:)= fac1*radius_scenario_ls(is:ie,js:je,iaam,:) + (1.0-fac1)*radius_scenario_ls(is:ie,js:je,iaa,:)
endif

!---------------------- Tests --------------------------
if (interact.and.Background.and.dual_mode.eq.0) then
    print*, 'STOP in dust_update : interact and Background are not compatible'
    STOP
endif

!********************************************************************
!********************************************************************
!                    1/   Interactive dust
!********************************************************************
!********************************************************************
if (interact) then 
     
!**************************
!**************************
! 1.1/ Dust devils lifting 
!**************************
!**************************
    if ( convective_lift ) then

        xsi = rdgas / cp_air +1.
        gama=0.5
        rm = reff_dd * exp(-0.5*dev_dt**2.)
        mass_part = 4./3.* pi * dens * rm**3.   

        !***********************
        ! get the injection flux
        !***********************
        if (dth) then ! newman et al 2002 sec 3.2
            vtan(:,:) = 0.0    
            b(:,:)=0.0 
            ! get density first layer : rhokd 
            tair(:,:)= tnew(:,:,kd)
            rhokd(:,:)=  p_full(:,:,kd)/( rdgas*tair(:,:) ) 
            ! get threshold tangential wind speed
            vtan_t(:,:) =  sqrt( 1. + 15. / (dens*grav*dlift) ) *sqrt( dens*grav*dlift/rhokd(:,:) )
            ! get efficiency eta and eta_h
            b(is:ie,js:je) = (p_full(is:ie,js:je,kd)**(xsi) - p_pbl(is:ie,js:je)**(xsi)) &
                        / ( (p_full(is:ie,js:je,kd)-p_pbl(is:ie,js:je)) * (xsi) * p_full(is:ie,js:je,kd)**xsi )
            eta(:,:)   = 1. - b(:,:)
            eta_h(:,:) = (tair(:,:)-tsfc(:,:)) / tsfc(:,:)
            ! get efficiency delta p
            dp_lift(:,:) = p_full(:,:,kd) * ( 1. - exp( (gama*eta(:,:))/(gama*eta(:,:)-1.)*(eta_h(:,:)/(xsi-1)) ) )
            ! get current vtan
            vtan(:,:)=dp_lift(:,:)/rhokd(:,:)

            where (vtan.ge.vtan_t)
                fluxdth(:,:) = alpha_dth * (rhokd(:,:)*vtan(:,:)**2. - 15.) / grav
            end where

        endif ! dth

        !***********************
        if (dda) then ! newman et al 2002 sec 3.2

            b(:,:)=0.0 
            do i=is,ie
                do j=js,je
                    if (shflx(i,j) .gt. thres_shflx .and. k_pbl(i,j) .lt. kd ) then  
                        b(i,j) = (p_full(i,j,kd)**(xsi) - p_pbl(i,j)**(xsi)) / ( (p_full(i,j,kd)-p_pbl(i,j)) * (xsi) * p_full(i,j,kd)**(xsi) )
                        fluxdda(i,j) = alpha_dda * (1.-b(i,j)) * shflx(i,j) 
                    endif
                enddo
            enddo

        endif ! dda
      
        !***********************
        if (devils_fv3) then ! version fv3

            ! source_ddfv3 given in namelist. rough conversion from opacity to mixing ratio
            tauref_ddfv3=1.
            pbar_ddfv3= 6.e2
            rhodust_ddfv3 = 2.50e3
            qext_ddfv3= 2.5
            reffec_ddfv3= 1.0e-6     !  effective particle radius
            source_scaling= tauref_ddfv3*(grav/pbar_ddfv3)*1.33*(rhodust_ddfv3/qext_ddfv3)*reffec_ddfv3

            ddfv3_source= source_ddfv3 * source_scaling ! equivalent column mass

            ! balance source with sink;  near-surface sedimentation varies roughly as 
            ! (rdust_core(nt)/reffec) *1.2e-2/200.0  
            fluxfv3(:,:)=ddfv3_source*( (reff_ddfv3/reffec_ddfv3)**2 )*1.2e-2/200.0

            ! require surface convective instability for dust injection
            where ( tsfc - tnew(:,:,kd) < conv_delt_threshold ) ! given in nml
                fluxfv3(:,:)= 0.0
            end where

        endif ! devils_fv3

        !***********************
        ! sum of the different dust devils source flux
        !***********************
        fluxtot(:,:)=fluxfv3(:,:)+fluxdda(:,:)+fluxdth(:,:)       !flux : kg.m-2-s-1
        !***********************
        ! cut flux options
        !***********************
        where( snow > dd_co2_ice_thresh_n .and. lat > 0.0 ) 
            fluxtot(:,:) = 0.0
        end where
        where( snow > dd_co2_ice_thresh_s .and. lat < 0.0 ) 
            fluxtot(:,:) = 0.0
        end where
        where( frost > dd_h2o_ice_thresh ) 
            fluxtot(:,:) = 0.0
        end where
        if (dual_mode.eq.1.or.dual_mode.eq.21) then
            where( dtau_current < 0. ) 
                flux_theo_dd(:,:)=fluxtot(:,:)
                fluxtot(:,:) = 0.0
            end where
        endif

        !***********************
        ! Limited reservoirs : depends on threshold for wind stress
        !***********************
        if (limited_reservoir) then
            where( threshold_stress*(exp(2.*limres_xfac*(dust_surf_ini-sfc_dust_mass(:,:,1)))-1.).gt.delta_thres)
                fluxtot(:,:) = 0.0
            end where
        endif

        !***********************
        ! FINAL FLUX AND SOURCES
        !***********************
        do nt= 1, ndust_mass
            flux(:,:,nt) = max(fluxtot(:,:)*dt,0.)

            !***********************
            ! Tagging methods 
            !***********************
            ndx= dust_mass_indx(nt)
#ifndef RELEASE
            call tagging_main( (/ "geosource" , "loctime", "antigeosource" , "solsource", "cutsource" , "inidust" , "custom_sources" /) ,lat,lon,ndx,field2D=flux(:,:,nt),solrun=sols, sources=sources_scenario*dt,inidust=inidust_scenario,flag='dd',stress=stress,source=flux(:,:,nt)/dt,taudust_curr=tau_current(:,:))
#endif
            source_dd(:,:,nt)=flux(:,:,nt)

            !***********************
            ! Injection near-surface or PBL
            !***********************
            Nadd(:,:,nt)    = 0. 
            Madd(:,:,nt)    = 0.
            if (inject_pbl_dd) then
                where (flux(:,:,nt) .gt. 0.)
                    madd(:,:,nt) = flux(:,:,nt)/masspbl(:,:)  
                    nadd(:,:,nt) = madd(:,:,nt) / mass_part  
                end where
                do i=is,ie
                    do j=js,je
                        do k= k_pbl(i,j), kd            ! from namelist
                            rnew(i,j,k,ndx)=rnew(i,j,k,ndx) + madd(i,j,nt)
                            rnew(i,j,k,ndx+1)=rnew(i,j,k,ndx+1) + nadd(i,j,nt)
                        enddo
                    enddo
                enddo
            else
                where (flux(:,:,nt) .gt. 0.)
                    madd(:,:,nt) = flux(:,:,nt)/masspblfix(:,:)   !mass1l(:,:)  ! (kgdust m-2 / kgair * m-2) = kgdust/kgair
                    nadd(:,:,nt) = madd(:,:,nt) / mass_part  
                end where
                do k= kd-kfix(2),kd-kfix(1)    ! from namelist
                    rnew(:,:,k,ndx)=rnew(:,:,k,ndx) + madd(:,:,nt)
                    rnew(:,:,k,ndx+1)=rnew(:,:,k,ndx+1) + nadd(:,:,nt)
                enddo
            endif 
            ! rm dust from ground
            sfc_dust_mass(is:ie,js:je,nt)=sfc_dust_mass(is:ie,js:je,nt)-source_dd(:,:,nt)  

        enddo  !! ndust mass

    endif ! convective dust source (dust devils)
 
!**************************
!**************************
! 1.2/ Wind stress lifting
!**************************
!**************************
    ! Inject dust according to stress. Different schemes possible to calculate flux
    if (stress_lift) then

        ! will need near surface density :
        tair(:,:)= tnew(:,:,kd)
        rhokd(:,:)=  p_full(:,:,kd)/( rdgas*tair(:,:) ) 
        ! will need these constants :
        Rm = Reff_stress * exp(-0.5*dev_dt**2.)
        mass_part = 4./3.* pi * dens * Rm**3.  

        !***********************
        !! compute flux from wind stress
        !***********************

        !----- Compute threshold ----
        if (limited_reservoir) then
            threshold(:,:)=threshold_stress*exp(2.*limres_xfac*(dust_surf_ini-sfc_dust_mass(:,:,1)))
            ! Refill System with decrease of threshold
            threshold(:,:)=threshold(:,:)-limres_dec*dt
        else
            threshold(:,:)=threshold_stress
        endif

        !--------- get flux --------- + limit flux according to ice or opacity threshold
        fluxtot(:,:)=0.
        call flux_stress( stress, threshold, alfa, rhokd, lat, snow,frost, tau_current(:,:), fluxtot(:,:))

        if (dual_mode.eq.1.or.dual_mode.eq.21) then
            where( dtau_current < 0. ) 
                flux_theo_ws(:,:)=fluxtot(:,:)
                fluxtot(:,:) = 0.0
            end where
        endif

        !***********************
        ! FINAL FLUX AND SOURCES
        !***********************
        do nt= 1, ndust_mass
            flux(:,:,nt) = max(fluxtot(:,:)*dt,0.)  !flux : kg.m-2   

            !***********************
            ! tagging methods 
            !***********************
            ndx= dust_mass_indx(nt)
#ifndef RELEASE
            call tagging_main( (/ "geosource" , "loctime", "antigeosource" , "solsource", "cutsource" , "inidust" , "custom_sources" /) ,lat,lon,ndx,field2D=flux(:,:,nt),solrun=sols, sources=sources_scenario*dt,inidust=inidust_scenario, flag="ws",stress=stress,source=flux(:,:,nt)/dt,taudust_curr=tau_current(:,:))
#endif
            source_ws(:,:,nt)=flux(:,:,nt)

            !***********************
            ! injection near-surface or pbl
            !***********************
            nadd(:,:,nt)    = 0. 
            madd(:,:,nt)    = 0.

            if (inject_pbl_stress) then
                where (flux(:,:,nt) .gt. 0.)
                    madd(:,:,nt) = flux(:,:,nt)/masspbl(:,:)  
                    nadd(:,:,nt) = madd(:,:,nt) / mass_part  
                end where
                do i=is,ie
                    do j=js,je
                        do k= k_pbl(i,j), kd            ! from namelist
                            rnew(i,j,k,ndx)=rnew(i,j,k,ndx) + madd(i,j,nt)
                            rnew(i,j,k,ndx+1)=rnew(i,j,k,ndx+1) + nadd(i,j,nt)
                        enddo
                    enddo
                enddo
            else
                where (flux(:,:,nt) .gt. 0.)
                    madd(:,:,nt) = flux(:,:,nt)/masspblfix(:,:)     !mass1l(:,:)  ! kgdust/kgair
                    nadd(:,:,nt) = madd(:,:,nt) / mass_part  
                end where
                do k= kd-kfix(2), kd-kfix(1)            ! from namelist
                    rnew(:,:,k,ndx)=rnew(:,:,k,ndx) + madd(:,:,nt)
                    rnew(:,:,k,ndx+1)=rnew(:,:,k,ndx+1) + nadd(:,:,nt)
                enddo
            endif
            ! rm dust from ground
            sfc_dust_mass(:,:,nt)=sfc_dust_mass(:,:,nt)-source_ws(:,:,nt)  

        enddo !! ndust_mass

    endif ! stress_lift

    !***************************************
    ! Calculation opacity tau_new and column 
    !***************************************
    do nt= 1, ndust_mass
        !---- for dust mass -----
        ndx= dust_mass_indx(nt)
        tau_new(:,:,nt)=0.
        dmcol(:,:,nt)=0.
        dncol(:,:,nt)=0.
        do k= 1, kd
            cst2(:,:)=delp(:,:,k)/grav*qextref(k)  
            if (nt.eq.1) then
                rn(:,:,k,nt) = ( rnew(:,:,k,ndx)/( rnew(:,:,k,ndx+1)+1.e-32) *cst )**(athird) * exp( dev**2. )
            else
                rn(:,:,k,nt) = rn(:,:,k,1)
            endif
            tau_new(:,:,nt) = tau_new(:,:,nt) + pi * rn(:,:,k,nt)**2 * rnew(:,:,k,ndx+1) * cst2(:,:)
            dmcol(:,:,nt) = dmcol(:,:,nt) + rnew(:,:,k,ndx) * ((p_half(:,:,k+1)-p_half(:,:,k))/grav)
            dncol(:,:,nt) = dncol(:,:,nt) + rnew(:,:,k,ndx+1) * ((p_half(:,:,k+1)-p_half(:,:,k))/grav)
        enddo

        !---- for core mass -----
        ndx= dust_cor_indx(nt)
        dccol(:,:,nt) = 0.
        do k= 1, kd
            dccol(:,:,nt) = dccol(:,:,nt) + rnew(:,:,k,ndx) * ((p_half(:,:,k+1)-p_half(:,:,k))/grav)
        enddo

    enddo !! ndust_mass

endif ! end interact

!********************************************************************
!********************************************************************
!                   2/  Assimilated dust
!********************************************************************
!********************************************************************
if (Background) then
    ! Initialisations
    flux(:,:,:)=0.
    rn(:,:,:,:)=0.
    do nt = 1,ndust_mass
        ndx= dust_mass_indx(nt)
        if (changeradius) then
            !rm2d(:,:,nt)=radius_scenario(:,:,1) * exp( -2.5 * dev**2. )
            rm2d(:,:,nt)=reff_mom(ndx) * exp( -2.5 * dev**2. )
            tau_frac2D(:,:,nt)=radius_scenario(:,:,nt)
        else
            rm2d(:,:,nt)=reff_mom(ndx) * exp( -2.5 * dev**2. ) 
        endif
    enddo

    rss2d(:,:,:) = rm2d(:,:,:) * exp( dev**2. )

    nadd(:,:,:)    = 0. 
    madd(:,:,:)    = 0.
    nrem(:,:,:)    = 0. 
    mrem(:,:,:)    = 0.
    injec(:,:,:)=0.
    sink(:,:,:)=0.

    ! compute tau_add
    tau_add(:,:) = tau_scenario(:,:)*tauscale - tau_current(:,:)  
    injscale2D(:,:)=injscale

    tautmp=tau_scenario(:,:)*tauscale/tau_current(:,:)
    select case (injec_type)
    
    case(2)
        injscale=SUM(tautmp,tautmp>1.)/(MAX(1.,real(count(tautmp>1.))))
        injscale=MAX(1.,injscale)
        injscale2D=injscale
    case(3)
        where (tautmp .GT. 1.)
            injscale2D=tautmp
        elsewhere
            injscale2D=injscale
        end where 
    case(4)
        where (tautmp .GT. 1.)
            injscale2D=exp(tautmp-1.)
        elsewhere
            injscale2D=injscale
        end where
    case(5)
        where (tautmp .GT. 1.)
            injscale2D=tautmp**2
        elsewhere
            injscale2D=injscale
        end where
    case(6)
        where (tautmp .GT. 1.)
            injscale2D=tautmp**3
        elsewhere
            injscale2D=injscale
        end where
    case(7)
        where (tautmp .GT. 1.)
            injscale2D=tautmp**4
        elsewhere
            injscale2D=injscale
        end where
    case default
        injscale2D=injscale
    end select

    if (ltfrac) then
        tau_add(:,:)=(ltfrac_curr(:,:)*ltscale+injscale2D)/injscale2D*tau_add(:,:)
        where (ltfrac_curr(:,:).le.0.)
            tau_add(:,:)=0.
        end where
    endif

    ! shut down injections for some values of dtau (optional, not used by default)
    where (dtau(:,:).le.dtaulim)
        tau_add(:,:)=0.
    end where
  
    !*******************************************************************************
    ! eliminate dust updates where observations are unreliable, or ice, or others...
    !*******************************************************************************
    where( tsfc < assim_t_thresh )  
        tau_add(:,:) = 0.
    end where
    if( no_assim_over_caps )  then
        where( snow > 0.9 )
            tau_add(:,:) = 0.
        end where
    endif
    if( id_tauadd > 0 )  used = send_data( id_tauadd, tau_add(:,:), Time, is, js )

    !***********************
    ! injection near-surface or pbl
    !***********************
    do nt=1, ndust_mass
        ndx= dust_mass_indx(nt)

        if (changeradius) then
          where ( tau_add .gt. 0. )
            injec(:,:,nt)   = tau_frac2D(:,:,nt)*4./3.*injscale2D*tau_add(:,:)*dens*rm2d(:,:,nt)**3. &
                    * exp( 4.5*dev**2. )*dt / (qext_dt(ndx,dgdm_type) * rss2d(:,:,nt)**2. * seconds_per_day) !! kg m-2 per timestep 
          end where  
        else
          where ( tau_add .gt. 0. )
            injec(:,:,nt)   = tau_frac(ndx)*4./3.*injscale2D*tau_add(:,:)*dens*rm2d(:,:,nt)**3. &
                    * exp( 4.5*dev**2. )*dt / (qext_dt(ndx,dgdm_type) * rss2d(:,:,nt)**2. * seconds_per_day) !! kg m-2 per timestep 
          end where  
        endif
    
        if (sink_bd) then
            sink(:,:,nt)   = tau_frac(ndx)*sinkscale*(tau_scenario(:,:)*tauscale - tau_current(:,:)) &
                     / (qext_dt(ndx,dgdm_type) * pi * rss2d(:,:,nt)**2. ) / ((p_half(:,:,kd+1)-p_half(:,:,1))/grav) !! kg-1 per sol 
            sink(:,:,nt)   = min(sink(:,:,nt),0.0)
        endif
        if( id_inject(nt) > 0 )  used = send_data( id_inject(nt), injec(:,:,nt), time, is, js )
    enddo

    !***********************
    ! final flux and sources
    !***********************
    do nt= 1, ndust_mass
        ndx= dust_mass_indx(nt)
        flux(:,:,nt)=injec(:,:,nt)

        !***********************
        ! tagging methods 
        !***********************
#ifndef RELEASE
        call tagging_main( (/ "geosource" , "loctime", "antigeosource" , "solsource",  &
                    "cutsource" , "custom_sources" , "inidust" /) ,lat,lon,ndx, &
                    field2d=flux(:,:,nt),solrun=sols,sources=sources_scenario*dt, &
                    inidust=inidust_scenario,flag='bg',stress=stress, source=flux(:,:,nt)/dt) 
#endif
        source_bg(:,:,nt)=flux(:,:,nt)

        !***********************
        ! compute mass of dust to be injected
        !***********************
        if (inject_pbl_bd) then
            madd(:,:,nt)=flux(:,:,nt)/masspbl(:,:)
            nadd(:,:,nt)=madd(:,:,nt)/(4./3.*pi*dens*rm2d(:,:,nt)**3.*exp(4.5*dev**2.))
            do i=is,ie
                do j=js,je
                    do k= k_pbl(i,j), kd            ! from namelist
                        rnew(i,j,k,ndx)=rnew(i,j,k,ndx) + madd(i,j,nt)
                        rnew(i,j,k,ndx+1)=rnew(i,j,k,ndx+1) + nadd(i,j,nt)
                    enddo
                enddo
            enddo 
        else
            madd(:,:,nt)=flux(:,:,nt)/masspblfix(:,:)   !mass1l(:,:)
            nadd(:,:,nt)=madd(:,:,nt)/(4./3.*pi*dens*rm2d(:,:,nt)**3.*exp(4.5*dev**2.))
            do k= kd-kfix(2), kd-kfix(1)            ! from namelist
                rnew(:,:,k,ndx)=rnew(:,:,k,ndx) + madd(:,:,nt)
                rnew(:,:,k,ndx+1)=rnew(:,:,k,ndx+1) + nadd(:,:,nt)
            enddo
        endif
        ! Remove dust from ground
        sfc_dust_mass(is:ie,js:je,nt)=sfc_dust_mass(is:ie,js:je,nt)-source_bg(:,:,nt) 

        !***********************
        ! removing dust by sinks
        !***********************
        if (sink_bd) then
        ratio_sinkm(:,:)=0.
        ratio_sinkn(:,:)=0.
            !! first tracer : amount of dust to be removed : nrem, mrem
            if (ndx .le. nt_nontag) then
                nrem(:,:,nt)   = -sink(:,:,nt) ! put n positive
                mrem(:,:,nt)   = 4./3.*pi * nrem(:,:,nt) * dens * rm2d(:,:,nt)**3. * exp( 4.5*dev**2. )  !! kgdust/kgair per sol
                nrem(:,:,nt)   = nrem(:,:,nt) / seconds_per_day * dt   ! per timestep
                mrem(:,:,nt)   = mrem(:,:,nt) / seconds_per_day * dt   

                !! total current mass and number 
                ntot(:,:)=0.
                mtot(:,:)=0.
                do k= 1, kd      
                    mtot(:,:)=mtot(:,:)+rnew(:,:,k,ndx)*(p_half(:,:,k+1)-p_half(:,:,k))
                    ntot(:,:)=ntot(:,:)+rnew(:,:,k,ndx+1)*(p_half(:,:,k+1)-p_half(:,:,k))
                enddo
                mtot(:,:)=mtot(:,:)/(p_half(:,:,kd+1)-p_half(:,:,1)) !total kg/kg
                ntot(:,:)=ntot(:,:)/(p_half(:,:,kd+1)-p_half(:,:,1))

                !! mass mixing ratio to be removed: security
                mrem(:,:,nt)=min(mrem(:,:,nt),mtot(:,:)+1.e-32)
                nrem(:,:,nt)=min(nrem(:,:,nt),ntot(:,:)+1.e-32)

                !! if removing from bottom
                if (bottom_sink) then
                    mtot_col(:,:)=mtot(:,:)*(p_half(:,:,kd+1)-p_half(:,:,1))/grav
                    mrem_col(:,:)=mrem(:,:,nt)*(p_half(:,:,kd+1)-p_half(:,:,1))/grav
                    ntot_col(:,:)=ntot(:,:)*(p_half(:,:,kd+1)-p_half(:,:,1))/grav
                    nrem_col(:,:)=nrem(:,:,nt)*(p_half(:,:,kd+1)-p_half(:,:,1))/grav
                    mtmp(:,:)=0.
                    ntmp(:,:)=0.
                    do i=is,ie
                        do j=js,je

                            k=kd
                            do while (Mtmp(i,j).lt.Mrem_col(i,j)) 
                                if (Mtmp(i,j)+rnew(i,j,k,ndx)*(p_half(i,j,k+1)-p_half(i,j,k))/grav.lt.Mrem_col(i,j)) then
                                    rnew(i,j,k,ndx)=0.
                                    Mtmp(i,j)=Mtmp(i,j)+rnew(i,j,k,ndx)*(p_half(i,j,k+1)-p_half(i,j,k))/grav
                                else
                                    rnew(i,j,k,ndx)=rnew(i,j,k,ndx)-(Mrem_col(i,j)-Mtmp(i,j))*grav/(p_half(i,j,k+1)-p_half(i,j,k))
                                    Mtmp(i,j)=Mrem_col(i,j)
                                endif
                                k=k-1
                            enddo
                            k=kd
                            do while (Ntmp(i,j).lt.Nrem_col(i,j)) 
                                if (Ntmp(i,j)+rnew(i,j,k,ndx+1)*(p_half(i,j,k+1)-p_half(i,j,k))/grav.lt.Nrem_col(i,j)) then
                                    rnew(i,j,k,ndx+1)=0.
                                    Ntmp(i,j)=Ntmp(i,j)+rnew(i,j,k,ndx+1)*(p_half(i,j,k+1)-p_half(i,j,k))/grav
                                else
                                    rnew(i,j,k,ndx+1)=rnew(i,j,k,ndx+1)-(Nrem_col(i,j)-Ntmp(i,j))*grav/(p_half(i,j,k+1)-p_half(i,j,k))
                                    Ntmp(i,j)=Nrem_col(i,j)
                                endif
                                k=k-1
                            enddo

                        enddo
                    enddo

                else           

                    !! pourcentage to be removed in each cell
                    ratio_sinkm(is:ie,js:je)=mrem(is:ie,js:je,nt)/mtot(is:ie,js:je)   ! per timestep
                    ratio_sinkn(is:ie,js:je)=nrem(is:ie,js:je,nt)/ntot(is:ie,js:je)   ! per timestep
                    !! remove for all tracers
                    do k= 1, kd     
                        rnew(:,:,k,ndx)=rnew(:,:,k,ndx)*(1.-ratio_sinkm(:,:))
                        rnew(:,:,k,ndx+1)=rnew(:,:,k,ndx+1)*(1.-ratio_sinkn(:,:))
                    enddo

                endif  ! if bottom sinks

                    !! update dust budget
                    fluxrem(is:ie,js:je,nt)  = min(mrem(is:ie,js:je,nt),mtot(i,j)) &
                                * ((p_half(is:ie,js:je,kd+1)-p_half(is:ie,js:je,1))/grav) 

            endif ! nt=1

            ! add dust from ground
            sfc_dust_mass(is:ie,js:je,nt)=sfc_dust_mass(is:ie,js:je,nt)+fluxrem(:,:,nt)  
        endif ! sink bd

    enddo ! ndust mass

endif ! End of background test

!********************************************************************
!********************************************************************
!                    3/ Update dust tendencies
!********************************************************************
!********************************************************************
do nt= 1, ndust_mass
    ndx= dust_mass_indx(nt)
    rdt_dst(:,:,:,ndx)= ( rnew(:,:,:,ndx) - rini(:,:,:,ndx) )/dt
    rdt_dst(:,:,:,ndx+1)= ( rnew(:,:,:,ndx+1) - rini(:,:,:,ndx+1) )/dt
    if(lifting_nosed) then
        where(rdt_dst(:,:,kd,ndx) .gt. 0.) 
            lifting_dust(:,:,ndx)=.true.
            lifting_dust(:,:,ndx+1)=.true.
        end where
    end if
enddo

!********************************************************************
!********************************************************************
!                          4/ Outputs
!********************************************************************
!********************************************************************
!*********************************
! Get new opacity and column mass 
!*********************************
do nt= 1, ndust_mass
    ndx= dust_mass_indx(nt)
    !---- dust mass opacity and column
    tau_new(:,:,nt) = 0.
    tau_new2(:,:,nt) = 0.
    dmcol(:,:,nt) = 0.
    dncol(:,:,nt) = 0.
    do k= 1, kd
        cst2(:,:)=delp(:,:,k)/grav*qextref(k)   
        rn(:,:,k,nt) = ( rnew(:,:,k,ndx)/( rnew(:,:,k,ndx+1)+1.e-32) *cst )**(athird) * exp( dev**2. )
        tau_new(:,:,nt) = tau_new(:,:,nt) + pi * rn(:,:,k,nt)**2. * rnew(:,:,k,ndx+1) * cst2(:,:)
        tau_new2(:,:,nt) = tau_new2(:,:,nt) + 3./4.*rnew(:,:,k,ndx)*2./(dpden_dt*rn(:,:,k,nt)*grav)*(p_half(:,:,k+1)-p_half(:,:,k))
        dmcol(:,:,nt) = dmcol(:,:,nt) + rnew(:,:,k,ndx) * ((p_half(:,:,k+1)-p_half(:,:,k))/grav)
        dncol(:,:,nt) = dncol(:,:,nt) + rnew(:,:,k,ndx+1) * ((p_half(:,:,k+1)-p_half(:,:,k))/grav)
    enddo

    !---- core mass column
    ndx= dust_cor_indx(nt)
    dccol(:,:,nt) = 0.
    do k= 1, kd
        dccol(:,:,nt) = dccol(:,:,nt) + rnew(:,:,k,ndx) * ((p_half(:,:,k+1)-p_half(:,:,k))/grav)
    enddo

    !---- Output sfc reservoir / sources / sinks
    if (Background.and.dual_mode.eq.0) then
        source_mass(:,:,nt)=source_bg(:,:,nt)/dt    !! / s
    elseif (Background.and.dual_mode.ge.2) then
        source_mass(:,:,nt)=(source_bg(:,:,nt)+source_dd(:,:,nt)+source_ws(:,:,nt))/dt   !! /s
    else
        source_mass(:,:,nt)=(source_dd(:,:,nt)+source_ws(:,:,nt))/dt   !! /s
    endif

    if( id_dust_mass_accum(nt) > 0 )  used = send_data( id_dust_mass_accum(nt),  &
                             sfc_dust_mass(is:ie,js:je,nt), Time, is, js )
    if( id_dust_mass_src(nt) > 0 )  used = send_data( id_dust_mass_src(nt), source_mass(:,:,nt), Time, is, js )
    if( id_dust_mass_src_ws(nt) > 0 )  used = send_data( id_dust_mass_src_ws(nt), source_ws(:,:,nt)/dt, Time, is, js )
    if( id_dust_mass_src_dd(nt) > 0 )  used = send_data( id_dust_mass_src_dd(nt), source_dd(:,:,nt)/dt, Time, is, js )
    if( id_dust_mass_src_bg(nt) > 0 )  used = send_data( id_dust_mass_src_bg(nt), source_bg(:,:,nt)/dt, Time, is, js )
    if( id_dust_mass_sink(nt) > 0 )  used = send_data( id_dust_mass_sink(nt), fluxrem(:,:,nt)/dt, Time, is, js )

    !---- Output colums mass / opacity 
    if( id_dmcol(nt) > 0 )  used = send_data( id_dmcol(nt), dmcol(:,:,nt), Time, is, js )
    if( id_dccol(nt) > 0 )  used = send_data( id_dccol(nt), dccol(:,:,nt), Time, is, js )
    if( id_dmtau(nt) > 0 )  used = send_data( id_dmtau(nt), tau_new(:,:,nt), Time, is, js )
    if( id_dmtau2(nt) > 0 )  used = send_data( id_dmtau2(nt), tau_new2(:,:,nt), Time, is, js )

    !---- Useful for diagnostics
    tcol_core(:,:,nt)=dccol(:,:,nt)
    tcol_mass(:,:,nt)=dmcol(:,:,nt)

    if( id_rm2d(nt) > 0 )  used = send_data( id_rm2d(nt), Rm2D(:,:,nt), Time, is, js )
    if( id_taufrac2d(nt) > 0 )  used = send_data( id_taufrac2d(nt), tau_frac2D(:,:,nt), Time, is, js )

enddo ! ndust_mass

! If to be injected for dust bins
source_mom(:,:)=source_mass(:,:,1)

if (dual_mode.eq.1) then
    if( id_src_theo_dd > 0 )  used = send_data( id_src_theo_dd, flux_theo_dd(:,:), Time, is, js )   !/s
    if( id_src_theo_ws > 0 )  used = send_data( id_src_theo_ws, flux_theo_ws(:,:), Time, is, js )   !/s
endif

if( id_threshold_stress > 0 )  used = send_data( id_threshold_stress, threshold(:,:), Time, is, js )
if( id_tauscenario > 0 )  used = send_data( id_tauscenario, tau_scenario(:,:), Time, is, js )
if( id_dtauscenario > 0 )  used = send_data( id_dtauscenario, dtau(:,:), Time, is, js )

return
end

!********************************************************************
!#######################################################################
!********************************************************************

subroutine flux_stress( stress, threshold, alfa, rhobot, latij,snow,frost, tau, flux) 
!
!   dust flux from stress lifting
!

use constants_mod, only: grav
real,    intent(in), dimension(:,:)     ::  threshold
real,    intent(in)     ::  alfa
real,    intent(in), dimension(:,:)     ::  stress
real,    intent(in), dimension(:,:)     ::  rhobot       ! density layer above surface
real,    intent(in), dimension(:,:)     ::  latij
real,    intent(in), dimension(:,:)     ::  snow
real,    intent(in), dimension(:,:)     ::  frost
real,    intent(in), dimension(:,:)     ::  tau
real,    intent(out), dimension(:,:)    ::  flux

flux(:,:)= 0.0
select case (stress_scheme) ! from namelist
    case(1)   ! Mel lifting
        flux(:,:)=alfa*2.43E-3*(stress(:,:)**2)*(stress(:,:)-threshold(:,:))/threshold(:,:)
    case(2)   ! John lifting
        flux(:,:)=alfa/grav*(1-(threshold(:,:)/stress(:,:))**0.5)*(1+(threshold(:,:)/stress(:,:))**0.5)**2*(stress(:,:)**3/rhobot(:,:))**0.5  
    case(3)   ! John lifting with a3=1
        flux(:,:)=alfa/grav*(stress(:,:)**3/rhobot(:,:))**0.5  
    case(4)   ! White 1979
        flux(:,:)=alfa*2.61/grav*(1-(threshold(:,:)/stress(:,:))**0.5)*(1+(threshold(:,:)/stress(:,:))**0.5)**2*(stress(:,:)**3/rhobot(:,:))**0.5  
    case(5)   ! Shao 1993
        flux(:,:)=alfa*2.3E-3*(stress(:,:)**2)*(stress(:,:)-threshold(:,:))/threshold(:,:)
    case default  ! Mel lifting
        flux(:,:)=alfa*2.43E-3*(stress(:,:)**2)*(stress(:,:)-threshold(:,:))/threshold(:,:)
end select

flux(:,:)= max( 0.0, flux(:,:) )

!!! lifting only if stress higher than threshold (namelist)
where( stress <  threshold )
    flux(:,:) =0.
end where

!!! no lifting over co2 ice surfaces  north and south
where( snow > ws_co2_ice_thresh_n .and.latij > 0.0 ) 
    flux(:,:) = 0.0
end where
where( snow > ws_co2_ice_thresh_s .and.latij < 0.0 ) 
    flux(:,:) = 0.0
end where
where( frost > ws_h2o_ice_thresh ) 
    flux(:,:) = 0.0
end where

!!! use this to limit lifting when dust optical depth is already very high
where( tau > optd_thresh ) 
    flux(:,:) = 0.0
end where

end subroutine flux_stress

!********************************************************************
!#######################################################################
!********************************************************************

subroutine dust_update_init( nlon, mlat, lonb, latb, lon, lat, axes, Time, phys_domain )
! 
!   initialize dust update
!

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 
use tracer_manager_mod, only: query_method, get_tracer_index,  &
                                 get_number_tracers, get_tracer_names

use horiz_interp_mod,  only: horiz_interp, horiz_interp_init, horiz_interp_new,  &
                             horiz_interp_end,  horiz_interp_del, horiz_interp_type

!! Arguments
integer, intent(in)                   :: nlon, mlat
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon, lat
integer, intent(in)                   :: axes(4)
type(time_type), intent(in) :: Time
type(domain2d),      intent(inout) :: phys_domain

!! Local
integer  unit, io, ierr 
integer  id, jd, km, i, j, k, is, js, ie, je, nt, ndx,n
character (len=128) :: filename, fieldname, tracer_name, tname
character (len=2) :: aString
integer   ::   id_inpt, jd_inpt, areo_length, fld_dims(4), iaa, iaa2

real,  dimension(:,:,:), allocatable  :: source_inpt
real,  dimension(:,:), allocatable  :: sfc_source_inpt
real, dimension(:),  allocatable  ::   lon_inpt,  lat_inpt
real, dimension(:),  allocatable  ::   lonb_inpt, latb_inpt
real,  dimension(size(lat,1),size(lat,2))  :: localtime

type(horiz_interp_type)  ::  Interp
integer  nma_dst
integer :: days, seconds, iaam 
real    :: areolat, fac1, delareo,secs,sols,fjd,r_orbit,declin         

mcpu0 = (mpp_pe() == mpp_root_pe())

is= 1
js= 1
id= size(lon,1)
jd= size(lat,2)
ie= is + id - 1
je= js + jd - 1

! *********************************************************
!     ----- read namelist /dust_update_nml/   -----
! *********************************************************

if (firstcall) then 
    if (file_exists('input.nml')) then
        unit = open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
            read  (unit, nml=dust_update_nml, iostat=io, end=10)
            ierr = check_nml_error (io, 'dust_update_nml')
        enddo
10     call close_file (unit)
    endif
    firstcall=.false.
endif

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=dust_update_nml)

! *********************************************************
!     ----- register diagnostic fields -----
! *********************************************************

allocate ( id_dust_mass_accum(ndust_mass) )
allocate ( id_rn(ndust_mass) )
allocate ( id_dmcol(ndust_mass) )
allocate ( id_dncol(ndust_mass) )
allocate ( id_dccol(ndust_mass) )
allocate ( id_dmtau(ndust_mass) )
allocate ( id_dmtau2(ndust_mass) )
allocate ( id_dust_mass_src(ndust_mass) )
allocate ( id_dust_mass_src_ws(ndust_mass) )
allocate ( id_dust_mass_src_dd(ndust_mass) )
allocate ( id_dust_mass_src_bg(ndust_mass) )
allocate ( id_dust_mass_sink(ndust_mass) )
allocate ( id_rm2d(ndust_mass) )
allocate ( id_taufrac2d(ndust_mass) )
allocate ( id_inject(ndust_mass) )

!! Preparing output
do nt= 1, ndust_mass
     
    !Get name tracer
    ndx= dust_mass_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)

    !nma_dst= find_field_index( MODEL_ATMOS, 'dst_mass_mom' )
    !call get_tracer_names(MODEL_ATMOS, nma_dst, tracer_name)

    tname= trim(tracer_name) // '_surf'  
    id_dust_mass_accum(nt) = register_diag_field ( mod_name, trim(tname),  &
         axes(1:2), Time, 'Surface dust mass', 'kg/m/m', &
         missing_value=missing_value)

    tname= trim(tracer_name) // '_tau'
    id_dmtau(nt) = register_diag_field ( mod_name, trim(tname), &
                  (/axes(1:2)/),  Time, &
                'Dust mass column opacity VIS', '',  missing_value=missing_value )

    tname= trim(tracer_name) // '_tau2'
    id_dmtau2(nt) = register_diag_field ( mod_name, trim(tname), &
                  (/axes(1:2)/),  Time, &
                'Dust mass column opacity IR', '',  missing_value=missing_value )

    tname= trim(tracer_name) // '_col'
    id_dmcol(nt) = register_diag_field ( mod_name, trim(tname), &
                  (/axes(1:2)/),  Time, &
                'Dust mass column', 'kg/m/m',  missing_value=missing_value )

    tname= trim(tracer_name) // '_ncol'
    id_dncol(nt) = register_diag_field ( mod_name, trim(tname), &
                  (/axes(1:2)/),  Time, &
                'Dust number column', 'num/m/m',  missing_value=missing_value )

    tname=  trim(tracer_name) // '_source'
    id_dust_mass_src(nt) = register_diag_field ( mod_name, trim(tname),  &
                                   (/axes(1:2)/), Time,           &
                                  'Dust mass Source', 'kg/m/m/s',        &
                                   missing_value=missing_value )

    tname=  trim(tracer_name) // '_sourcews'
    id_dust_mass_src_ws(nt) = register_diag_field ( mod_name, trim(tname),  &
                                   (/axes(1:2)/), Time,           &
                                  'Dust mass Source WS', 'kg/m/m/s',        &
                                   missing_value=missing_value )
    tname=  trim(tracer_name) // '_sourcedd'
    id_dust_mass_src_dd(nt) = register_diag_field ( mod_name, trim(tname),  &
                                   (/axes(1:2)/), Time,           &
                                  'Dust mass Source DD', 'kg/m/m/s',        &
                                   missing_value=missing_value )
    tname=  trim(tracer_name) // '_sourcebg'
    id_dust_mass_src_bg(nt) = register_diag_field ( mod_name, trim(tname),  &
                                   (/axes(1:2)/), Time,           &
                                  'Dust mass Source BG', 'kg/m/m/s',        &
                                   missing_value=missing_value )

    tname=  trim(tracer_name) // '_sink'
    id_dust_mass_sink(nt) = register_diag_field ( mod_name, trim(tname),  &
                                   (/axes(1:2)/), Time,           &
                                  'Dust mass Sink', 'kg/m/m/s',        &
                                   missing_value=missing_value )

    tname= trim(tracer_name) // '_rad'
    id_rn(nt) = register_diag_field ( mod_name, trim(tname), &
                  (/axes(1:3)/),  Time, &
                'Dust particle radius for opacity', '',  missing_value=missing_value )

    !!! dust effective radius for background lifting
    tname= trim(tracer_name) // '_Rm2D'
    id_rm2d(nt) = register_diag_field ( mod_name, trim(tname), axes(1:2), Time, &
                      'dust effective radius for background lifting', 'm',  &
                       missing_value=missing_value     )

    !!! dust fraction for each mode
    tname= trim(tracer_name) // '_taufrac2d'
    id_taufrac2d(nt) = register_diag_field ( mod_name, trim(tname), axes(1:2), Time, &
                      'tau fraction per mode', '',  &
                       missing_value=missing_value     )
                       
    !!! injection for background lifting per dust mass tracer
    tname= trim(tracer_name) // '_inject'
    id_inject(nt) = register_diag_field ( mod_name, trim(tname), axes(1:2), Time, &
                  ' bacground lifting rate', ' ',  &
                   missing_value=missing_value   ) 
                   
    !!! Core mass
    ndx= dust_cor_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)

    tname= trim(tracer_name) // '_col'
    id_dccol(nt) = register_diag_field ( mod_name, trim(tname), &
                  (/axes(1:2)/),  Time, &
                'Dust core column mass', 'kg/m/m',  missing_value=missing_value )
enddo

id_threshold_stress = register_diag_field ( mod_name, 'threshold_stress', axes(1:2), Time, &
                  ' threshold stress', 'N m-2',  &
                   missing_value=missing_value     )
id_tauscenario = register_diag_field ( mod_name, 'tau_scenario', axes(1:2), Time, &
                  ' dust scenario opacity for closest previous', ' ',  &
                   missing_value=missing_value     )
id_dtauscenario = register_diag_field ( mod_name, 'dtau_scenario', axes(1:2), Time, &
                  'dust scenario tau difference between closest sols', ' ',  &
                   missing_value=missing_value     )
id_tauscenario_ini = register_diag_field ( mod_name, 'tau_sce_ini', axes(1:2), Time, &
                  ' interpolated dust scenario tau', ' ',  &
                   missing_value=missing_value     )
id_taucurrent = register_diag_field ( mod_name, 'tau_current', axes(1:2), Time, &
                  ' current model dust tau', ' ',  &
                   missing_value=missing_value     )
id_ltfrac = register_diag_field ( mod_name, 'ltfrac', axes(1:2), Time, &
                  ' local time lifting fraction', ' ',  &
                   missing_value=missing_value     )
id_tau_old = register_diag_field ( mod_name, 'tau_old', axes(1:2), Time, &
                  ' previous time step dust scenario tau', ' ',  &
                   missing_value=missing_value     )
id_tauadd = register_diag_field ( mod_name, 'tauadd', axes(1:2), Time, &
                  ' dust tau to be added after modifications', ' ',  &
                   missing_value=missing_value     )

if (dual_mode.eq.1) then
    id_src_theo_dd = register_diag_field ( mod_name, 'src_theo_dd', axes(1:2), Time, &
                  ' dual mode dust devil lifting', 'kg m-2 s-1',  &
                   missing_value=missing_value     )
    id_src_theo_ws = register_diag_field ( mod_name, 'src_theo_ws', axes(1:2), Time, &
                  ' dual mode wind stress lifting', 'kg m-2 s-1',  &
                   missing_value=missing_value     )
endif

! *********************************************************
!     -----  Read BINS / Areo / lon / lat file ------- e.g local times for injection
! *********************************************************

filename= 'INPUT/localtime_frac.nc' 
if( file_exists( trim( filename ) ) ) then 
    ltfrac=.true.
endif

filename= 'INPUT/sources.nc' 
if( file_exists( trim( filename ) ) ) then 
    custom_sources=.true.
endif

filename= 'INPUT/inidust.nc' 
if( file_exists( trim( filename ) ) ) then 
    inidust=.true.
endif

filename= 'INPUT/radius_scenario.nc' 
if( file_exists( trim( filename ) ) ) then 
    changeradius=.true.
endif


! *********************************************************
! *********************************************************
!          ----- BACKGROUND DUST MODE -----
! *********************************************************
! *********************************************************
  
if( Background.or.custom_sources.or.dual_mode.gt.0 ) then 

   ! *********************************************************
   !          ----- Read dust scenario  ----- area / lon / lat
   ! *********************************************************

    filename= 'INPUT/dust_cycle.nc' 
    if( file_exists( trim( filename ) ) ) then 

        call field_size( trim(filename), 'areo', fld_dims )
        areo_length= fld_dims(1)
        call field_size( trim(filename), 'lon', fld_dims )
        id_inpt= fld_dims(1)
        call field_size( trim(filename), 'lat', fld_dims )
        jd_inpt= fld_dims(1)

        allocate (  areo_dust_scenario(areo_length)  ) 
        call read_data( trim(filename), 'areo', areo_dust_scenario, no_domain=.true. )

        if(mcpu0) print *, 'Have read tes areo data file in dust_update: ' 

#ifndef CUBE_CORE
        if( nlon==id_inpt  .and. mlat==jd_inpt ) then
            if(mcpu0) write(*,*) '   --->  Dust cycle field sizes match:  no interpolation necessary '  

            allocate (  source_inpt(id_inpt,jd_inpt,areo_length)  ) 
            allocate (  tau_scenario_ls(is:ie,js:je,  areo_length)  ) 

            call read_data( trim(filename), 'taufill', source_inpt )

            if(mcpu0) print *, 'Have read tes taufill data file in dust_update (cube)' 

            tau_scenario_ls(is:ie,js:je,:) = source_inpt(is:ie,js:je,:) 

            ! IR VS VIS opacities
            if (dgdm_type.eq.2) then
                tau_scenario_ls=tau_scenario_ls/2.75
            else
                tau_scenario_ls=tau_scenario_ls*dust_map_scale/2.75
            endif

            deallocate ( source_inpt  )
            if(mcpu0)  print *, 'Successful setup for dust cycle in dust_update:', is, js, je, areo_length 

        else
#endif

            allocate( lat_inpt (jd_inpt) )
            allocate( latb_inpt(jd_inpt+1) )
            allocate( lon_inpt (id_inpt) )
            allocate( lonb_inpt(id_inpt+1) )

            call read_data( trim(filename), 'lat',  lat_inpt,  no_domain=.true. )
            call read_data( trim(filename), 'latb', latb_inpt, no_domain=.true. )
            call read_data( trim(filename), 'lon',  lon_inpt,  no_domain=.true. )
            call read_data( trim(filename), 'lonb', lonb_inpt, no_domain=.true. )
            latb_inpt(:)= latb_inpt(:)/RADIAN
            lonb_inpt(:)= lonb_inpt(:)/RADIAN

            allocate (  source_inpt  (id_inpt,jd_inpt,areo_length)  ) 

            call read_data( trim(filename), 'taufill', source_inpt, no_domain=.true. )
            if(mcpu0) print *, 'Have read tes taufill data file in dust_update: ' 

            allocate (  tau_scenario_ls(    is:ie,js:je,areo_length)  ) 

            ! Carry out horizontal interpolation 
            call horiz_interp_init

            call horiz_interp_new( Interp,lonb_inpt,latb_inpt,lon ,lat ,interp_method= 'bilinear' )

            do iaa= 1, areo_length
                call horiz_interp( interp, source_inpt (:,:,iaa), tau_scenario_ls(is:ie,js:je,iaa) )
            enddo 
            call horiz_interp_del( Interp )

            ! IR VS VIS opacities
            if (dgdm_type.eq.2) then
                tau_scenario_ls=tau_scenario_ls/2.75
            else
                tau_scenario_ls=tau_scenario_ls*dust_map_scale/2.75
            endif

            ! Check no negative values in dust scenario
            if (any(tau_scenario_ls(:,:,:) < 0.)) then 
                Print *, 'Problem in dust scenario dust_cycle.nc : negative values' 
                stop
            endif

            if(mcpu0)  print *, 'Successful setup for dust cycle in dust_update:', is, js, je, areo_length 

            deallocate ( source_inpt )
            deallocate ( lat_inpt, latb_inpt, lon_inpt, lonb_inpt )
#ifndef CUBE_CORE
        endif        !   horizontal interpolation option 
#endif  

    else      !       else use default values 
        print *, '############  WARNING:  missing dust cycle input fields '
    endif

endif  !! Background mode
! *********************************************************
! *********************************************************
! Read Scenario for local-time dependent lifting  
  
if( ltfrac ) then 
    filename= 'INPUT/localtime_frac.nc' 
    call field_size( trim(filename), 'bins', fld_dims )
    bins_length= fld_dims(1)
    call field_size( trim(filename), 'areo', fld_dims )
    areo_length= fld_dims(1)
    call field_size( trim(filename), 'lon', fld_dims )
    id_inpt= fld_dims(1)
    call field_size( trim(filename), 'lat', fld_dims )
    jd_inpt= fld_dims(1)

    allocate (  areo_ltfrac(areo_length)  ) 
    call read_data( trim(filename), 'areo', areo_ltfrac, no_domain=.true. )
    if(mcpu0) print *, 'Have read ltfrac areo data file in dust_update: ' 

    allocate (  bins_ltfrac(bins_length)  ) 
    call read_data( trim(filename), 'bins', bins_ltfrac, no_domain=.true. )


    allocate( lat_inpt (jd_inpt) )
    allocate( latb_inpt(jd_inpt+1) )
    allocate( lon_inpt (id_inpt) )
    allocate( lonb_inpt(id_inpt+1) )

    call read_data( trim(filename), 'lat',  lat_inpt,  no_domain=.true. )
    call read_data( trim(filename), 'latb', latb_inpt, no_domain=.true. )
    call read_data( trim(filename), 'lon',  lon_inpt,  no_domain=.true. )
    call read_data( trim(filename), 'lonb', lonb_inpt, no_domain=.true. )
    latb_inpt(:)= latb_inpt(:)/RADIAN
    lonb_inpt(:)= lonb_inpt(:)/RADIAN

    ! Keep in memory value of reference tau and reference bin
    allocate (  tau_ltfrac_old(is:ie,js:je)  )
    allocate (  bin_ltfrac_old(is:ie,js:je)  )

    !! Need first time to get first Ls and first tau
    call get_time(Time, seconds, days)
    secs= days * seconds_per_day + seconds
    sols = secs / seconds_per_day
    days = sols
    fjd = sols - days    !sol=0 at lon=180
    call mars_calender( days, fjd, r_orbit, declin, areolat )

    areolat = areolat * RADIAN
    if( areolat < 0.0 )  areolat = areolat + 360.0
    areolat = modulo(areolat,360.)

    do iaa= 1, size(areo_ltfrac)
        if( areo_ltfrac(iaa) > areolat )  then
            exit
        else
            cycle
        endif
    enddo
    ! current tau is tau_scenario at the end of the day (ensure all dust needed is injected
    tau_ltfrac_old(is:ie,js:je)= tau_scenario_ls(is:ie,js:je,iaa) 
    ! current bin depends on current local time
    localtime(:,:)=modulo( modulo(sols+0.5,1.)*24.-(180.-modulo(lon(:,:)*180./pi,360.))*12./180., 24.)
    do i=1,bins_length-1
        where(localtime(:,:).lt.bins_ltfrac(i+1).and.localtime(:,:).ge.bins_ltfrac(i))
            bin_ltfrac_old(:,:)=i
        end where
    enddo

    allocate (  ltfrac_scenario_ls(is:ie,js:je,areo_length,bins_length-1)  )

    ! Carry out horizontal interpolation 
    call horiz_interp_init
    call horiz_interp_new( Interp,lonb_inpt,latb_inpt,lon ,lat ,interp_method= 'bilinear' )

    do iaa2=1,bins_length-1

        ! read source
        allocate (  source_inpt(id_inpt,jd_inpt,areo_length)  ) 
        write(astring,'(i2.2)') bins_ltfrac(iaa2)
        call read_data( trim(filename), "ltfrac"//trim(astring), source_inpt, no_domain=.true. )
        if(mcpu0) print *, 'have read ltfrac data file in dust_update: ' 

        ! interpolate
        do iaa= 1, areo_length
            call horiz_interp( interp, source_inpt (:,:,iaa), ltfrac_scenario_ls(is:ie,js:je,iaa,iaa2) )
        enddo 

        deallocate ( source_inpt )

    enddo

    call horiz_interp_del( interp )
    deallocate ( lat_inpt, latb_inpt, lon_inpt, lonb_inpt )
endif

! *********************************************************
! *********************************************************
! Read Scenario for iparticle size changes over time  
if( changeradius ) then 
    filename= 'INPUT/radius_scenario.nc' 
    call field_size( trim(filename), 'bins', fld_dims )
    radbins_length= fld_dims(1)
    call field_size( trim(filename), 'areo', fld_dims )
    areo_length= fld_dims(1)
    call field_size( trim(filename), 'lon', fld_dims )
    id_inpt= fld_dims(1)
    call field_size( trim(filename), 'lat', fld_dims )
    jd_inpt= fld_dims(1)

    allocate (  areo_radius(areo_length)  ) 
    call read_data( trim(filename), 'areo', areo_radius, no_domain=.true. )
    if(mcpu0) print *, 'Have read radius_scenario areo data file in dust_update: ' 

    allocate (  bins_radius(radbins_length)  ) 
    call read_data( trim(filename), 'bins', bins_radius, no_domain=.true. )

    allocate( lat_inpt (jd_inpt) )
    allocate( latb_inpt(jd_inpt+1) )
    allocate( lon_inpt (id_inpt) )
    allocate( lonb_inpt(id_inpt+1) )

    call read_data( trim(filename), 'lat',  lat_inpt,  no_domain=.true. )
    call read_data( trim(filename), 'latb', latb_inpt, no_domain=.true. )
    call read_data( trim(filename), 'lon',  lon_inpt,  no_domain=.true. )
    call read_data( trim(filename), 'lonb', lonb_inpt, no_domain=.true. )
    latb_inpt(:)= latb_inpt(:)/RADIAN
    lonb_inpt(:)= lonb_inpt(:)/RADIAN

    !! need first time to get first ls and first tau
    call get_time(time, seconds, days)
    secs= days * seconds_per_day + seconds
    sols = secs / seconds_per_day
    days = sols
    fjd = sols - days    !sol=0 at lon=180
    call mars_calender( days, fjd, r_orbit, declin, areolat )

    areolat = areolat * radian
    if( areolat < 0.0 )  areolat = areolat + 360.0
    areolat = modulo(areolat,360.)

    do iaa= 1, size(areo_radius)
        if( areo_radius(iaa) > areolat )  then
            exit
        else
            cycle
        endif
    enddo
       
    allocate (  radius_scenario_ls(is:ie,js:je,areo_length,radbins_length)  )

    ! Carry out horizontal interpolation 
    call horiz_interp_init
    call horiz_interp_new( Interp,lonb_inpt,latb_inpt,lon ,lat ,interp_method= 'bilinear' )

    do iaa2=1,radbins_length

        ! read source
        allocate (  source_inpt(id_inpt,jd_inpt,areo_length)  ) 
        write(astring,'(i2.2)') bins_radius(iaa2)
        call read_data( trim(filename), "rad"//trim(astring), source_inpt, no_domain=.true. )
        if(mcpu0) print *, 'have read radius_scenario data file in dust_update: ' 

        ! interpolate
        do iaa= 1, areo_length
            call horiz_interp( interp, source_inpt (:,:,iaa), radius_scenario_ls(is:ie,js:je,iaa,iaa2) )
        enddo 

        deallocate ( source_inpt )
    enddo

    call horiz_interp_del( interp )
    deallocate ( lat_inpt, latb_inpt, lon_inpt, lonb_inpt )
endif

! *********************************************************
!     -----  Read custom source  ----- Areo / lon / lat
! *********************************************************

if( custom_sources ) then 
    filename= 'INPUT/sources.nc' 
    call field_size( trim(filename), 'areo', fld_dims )
    areo_length= fld_dims(1)
    call field_size( trim(filename), 'lon', fld_dims )
    id_inpt= fld_dims(1)
    call field_size( trim(filename), 'lat', fld_dims )
    jd_inpt= fld_dims(1)

    allocate (  areo_sources(areo_length)  ) 
    call read_data( trim(filename), 'areo', areo_sources, no_domain=.true. )

    if(mcpu0) print *, 'Have read sources areo data file in dust_update: ' 

    allocate( lat_inpt (jd_inpt) )
    allocate( latb_inpt(jd_inpt+1) )
    allocate( lon_inpt (id_inpt) )
    allocate( lonb_inpt(id_inpt+1) )

    call read_data( trim(filename), 'lat',  lat_inpt,  no_domain=.true. )
    call read_data( trim(filename), 'latb', latb_inpt, no_domain=.true. )
    call read_data( trim(filename), 'lon',  lon_inpt,  no_domain=.true. )
    call read_data( trim(filename), 'lonb', lonb_inpt, no_domain=.true. )
    latb_inpt(:)= latb_inpt(:)/RADIAN
    lonb_inpt(:)= lonb_inpt(:)/RADIAN

    allocate (  source_inpt  (id_inpt,jd_inpt,areo_length)  ) 

    call read_data( trim(filename), 'data', source_inpt, no_domain=.true. )
    if(mcpu0) print *, 'Have read sources data file in dust_update: ' 

    allocate (  sources_scenario_ls(is:ie,js:je,areo_length)  ) 

    ! Carry out horizontal interpolation 
    call horiz_interp_init

    call horiz_interp_new( Interp,lonb_inpt,latb_inpt,lon ,lat ,interp_method= 'bilinear' )


    do iaa= 1, areo_length
        call horiz_interp( interp, source_inpt (:,:,iaa), sources_scenario_ls(is:ie,js:je,iaa) )
    enddo 
    call horiz_interp_del( interp )

    deallocate ( source_inpt )
    deallocate ( lat_inpt, latb_inpt, lon_inpt, lonb_inpt )

endif

! *********************************************************
! -----  Read initial dust surface dust distribution  ----- lon / lat 
! *********************************************************

if( inidust ) then 
    filename= 'INPUT/inidust.nc' 
    call field_size( trim(filename), 'lon', fld_dims )
    id_inpt= fld_dims(1)
    call field_size( trim(filename), 'lat', fld_dims )
    jd_inpt= fld_dims(1)

    if(mcpu0) print *, 'Have read inidust areo data file in dust_update: ' 

    allocate( lat_inpt (jd_inpt) )
    allocate( latb_inpt(jd_inpt+1) )
    allocate( lon_inpt (id_inpt) )
    allocate( lonb_inpt(id_inpt+1) )

    call read_data( trim(filename), 'lat',  lat_inpt,  no_domain=.true. )
    call read_data( trim(filename), 'latb', latb_inpt, no_domain=.true. )
    call read_data( trim(filename), 'lon',  lon_inpt,  no_domain=.true. )
    call read_data( trim(filename), 'lonb', lonb_inpt, no_domain=.true. )
    latb_inpt(:)= latb_inpt(:)/RADIAN
    lonb_inpt(:)= lonb_inpt(:)/RADIAN

    allocate (  sfc_source_inpt  (id_inpt,jd_inpt)  ) 

    call read_data( trim(filename), 'distrib', sfc_source_inpt, no_domain=.true. )
    if(mcpu0) print *, 'Have read inidust taufill data file in dust_update: ' 

    allocate (  inidust_scenario(is:ie,js:je)  ) 

    ! Carry out horizontal interpolation 
    call horiz_interp_init

    call horiz_interp_new( Interp,lonb_inpt,latb_inpt,lon ,lat ,interp_method= 'bilinear' )

    call horiz_interp( Interp, sfc_source_inpt (:,:), inidust_scenario(is:ie,js:je) )
    call horiz_interp_del( Interp )

    deallocate ( sfc_source_inpt )
    deallocate ( lat_inpt, latb_inpt, lon_inpt, lonb_inpt )
endif


! *********************************************************
! *********************************************************
! ------------Read surface dust restart file if present (created from restart)--------------
! *********************************************************
! *********************************************************

if(mcpu0)  print *,  'Allocating space for surface dust mass fields '

allocate (  sfc_dust_mass (is:ie,js:je,ndust_mass)  )
allocate (  tcol_mass (is:ie,js:je,ndust_mass)  )
allocate (  tcol_core (is:ie,js:je,ndust_mass)  )

filename= 'INPUT/soil_accum2.res.nc'  

call dstup_register_restart(is,ie,js,je,'soil_accum2.res.nc',phys_domain)

if( file_exists( trim( filename ) ) ) then
    
    call restore_state(Dst_restart)
    if (in_different_file) call restore_state(Til_restart)  
    if(mcpu0) print *, 'Have read soil dust accm2 restart file:',ndust_mass_rst,'Bins'

    do nt= 1, ndust_mass
        if(mcpu0)  print *, 'surface dust mass',  nt, js, sfc_dust_mass(:,:,nt)
    enddo

else    !    Default case (eg cold case)
    sfc_dust_mass(is:ie,js:je,:)= dust_surf_ini ! in namelist  
endif

#ifndef RELEASE
! Tagging method for initial sources of dust
do nt= 1, ndust_mass
    ndx= dust_mass_indx(nt)
    call tagging_main( (/ "geosource" , "antigeosource" , "cutsource" , "inidust" /) , &
                lat,lon,ndx,field2d=sfc_dust_mass(:,:,nt),inidust=inidust_scenario,flag="in")
enddo    ! nt dust mass
#endif

module_is_initialized  = .true.

end subroutine dust_update_init

! =====================================================================
! =====================================================================


subroutine dust_update_end
! 
! Write soil dust accumulations to netCDF file


implicit none

call save_restart(Dst_restart)
if(in_different_file) call save_restart(Til_restart)

end subroutine dust_update_end

! =====================================================================
! =====================================================================

subroutine read_sfc_dust_mass_rst(  is, ie, js, je, km )
!
! Read surface dust accumulations from netCDF file:   
! fname= 'INPUT/soil_accum.res.nc'        field_name= 'sfc_dust'
!

implicit none

integer,             intent(in) :: is, ie, js, je
integer,             intent(out) :: km

!  ---------- local variables --------------------
character(len=128)   :: fname, filename
integer              :: i, j, k, im, jm, fld_dims(4), nmass
real, allocatable :: w3d(:,:,:)

fname= 'INPUT/soil_accum2.res.nc'

call field_size( trim(fname), 'sfc_dust_mass', fld_dims )

im= fld_dims(1);   jm= fld_dims(2);  km= fld_dims(3)
if(mcpu0)  print *, 'Surface dust mass restart dims: ',  im, jm, km, is, ie, js, je

allocate ( w3d (is:ie,js:je,km) )
nmass= MIN( km,ndust_mass )

call read_data(trim(fname), 'sfc_dust_mass', w3d )

! copy the first nmass into the sfc_dust_mass array 
do k=1, nmass
    sfc_dust_mass(is:ie,js:je,k) = w3d(is:ie,js:je,k)
enddo

deallocate( w3d )

end  subroutine read_sfc_dust_mass_rst


!--------------------------------------------------------
!--------------------------------------------------------

subroutine dstup_register_restart(is,ie,js,je,fname,phys_domain)
! register restart field to be written to restart file.
integer,                          intent(in) :: is,ie,js,je
character(len=*),                 intent(in) :: fname
character(len=64)                            :: fname2
type(domain2d),                intent(inout) :: phys_domain
integer :: id_restart, ndx, n
character (len=128) :: tracer_name

call get_mosaic_tile_file(fname, fname2, is_no_domain=.false., domain=phys_domain )

!default restart file: read/write definition
allocate(Dst_restart)
if(trim(fname2) == trim(fname)) then
    Til_restart => Dst_restart
    in_different_file = .false.
else
    in_different_file = .true.
    allocate(Til_restart)
endif

do n=1,ndust_mass
    ndx= dust_mass_indx(n)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)
    id_restart = register_restart_field(Til_restart, fname, trim(tracer_name), sfc_dust_mass(is:ie,js:je,n), domain=phys_domain,mandatory=.false.)
end do

id_restart = register_restart_field(Dst_restart, fname, 'dgdm_type', dgdm_type, no_domain = .true.,mandatory=.false.)



end subroutine dstup_register_restart

!--------------------------------------------------------
!--------------------------------------------------------
  

end module dust_update_mod





