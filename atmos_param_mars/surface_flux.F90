module  surface_flux_mod
!
!  calculate surface drag, heat and momentum fluxes
!

use   fms_mod, only: FATAL, close_file, mpp_pe, mpp_root_pe, write_version_number

use   fms_mod, only: check_nml_error, open_namelist_file, stdlog

use fms2_io_mod,            only:  file_exists
#ifdef fv3_turb
use   monin_obukhov_mod, only: mo_drag, mo_profile
#endif
use   time_manager_mod,   only: time_type

use   constants_mod, only: cp_air, pi, stefan, rdgas, rvgas, grav, kappa

implicit none
private

!  public interface 
public  surface_flux_2d


character(len=*), parameter :: version = '$Id: surface_flux.F90,v 1.1.2.1.2.1 2011/11/22 20:56:59 rjw Exp $'
character(len=*), parameter :: tagname = '$Name: mars_feb2012_rjw $'

logical :: do_init = .true.

real, parameter :: d622   = rdgas/rvgas
real, parameter :: d378   = 1.-d622
real, parameter :: gcp    = grav/cp_air
real            :: d608   = d378/d622
! d608 set to zero at initialization if the use of 
! virtual temperatures is turned off in namelist


logical :: use_virtual_temp      = .false. 
logical :: do_mo_drag            = .true. 
logical :: alt_gustiness         = .false.
logical :: old_dtaudv            = .false.
logical :: decouple_surface      = .false.
real    :: gust_const            =  1.0
real    :: gust_min              =  0.1
real    :: cd_drag_cnst          =  0.003
real    :: sfc_heat_flx_amp      = 1.0

namelist /surface_flux_nml/                       &
                         use_virtual_temp,     &
                         alt_gustiness,        &
                         gust_const,           &
                         gust_min,             &
                         old_dtaudv,           &
                         do_mo_drag,           &
                         cd_drag_cnst,         &
                         decouple_surface,     &
                         sfc_heat_flx_amp


contains


subroutine surface_flux_1d ( lon, lat,                        &
     t_atm,     q_atm_in,   u_atm,  v_atm, p_atm,  z_atm,     &
     p_surf,    t_surf,     q_surf,                           &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,   &
     flux_t, flux_q, flux_r, flux_u, flux_v,                  &
     cd_m,      cd_t,       cd_q,                             &
     w_atm,     u_star,     b_star,     q_star,               &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,            &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,           &
     is, js, dt, Time )
! 
! calculate 1d surface fluxes

!  arguments 
integer, intent(in)                      :: is, js
real, intent(in)                      :: dt
type(time_type), intent(in)                :: Time
real, intent(in),    dimension(:)     :: lon
real, intent(in),    dimension(:)     :: lat


real, intent(in),  dimension(:) ::                         &
    t_atm,     q_atm_in,   u_atm,     v_atm,              &
    p_atm,     z_atm,                                     &
    p_surf,    t_surf,                                    &
    rough_mom, rough_heat, rough_moist, rough_scale, gust

real, intent(out), dimension(:) ::                         &
    flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
    cd_m,      cd_t,       cd_q,                          &
    dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
    dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
    w_atm,     u_star,     b_star,    q_star

real, intent(inout), dimension(:) :: q_surf


!  local constants 
! temperature increment and its reciprocal value for comp. of derivatives
real, parameter:: del_temp=0.1, del_temp_inv=1.0/del_temp

!  local vars 

real, dimension(size(t_atm(:))) ::                       &
    thv_atm,  th_atm,   tv_atm,    thv_surf,            &
    e_sat,    e_sat1,   q_sat,     q_sat1,    p_ratio,  &
    t_surf0,  t_surf1,  u_dif,     v_dif,               &
    rho_drag, drag_t,    drag_m,   drag_q,    rho,      &
    q_atm,    q_surf0,  dw_atmdu,  dw_atmdv,  w_gust,   &
    u_surf,   v_surf


if (do_init) call surface_flux_init

!---- use local value of surf temp ----

t_surf0 = t_surf

t_surf1 = t_surf0 + del_temp


! surface mixing ratio at saturation
! surface specific humidity at saturation

q_sat=  0.401 * 6.11e2 * exp( 22.5*( 1.0 - 273.16/t_surf0 ) ) / p_surf 
q_sat1= 0.401 * 6.11e2 * exp( 22.5*( 1.0 - 273.16/t_surf1 ) ) / p_surf 


! initilaize surface air humidity according to surface type
q_surf0 = q_surf ! land calculates it

! check for negative atmospheric humidities
where( q_atm_in < 0.0) q_atm = 0.0

! generate information needed by monin_obukhov
p_ratio = (p_surf/p_atm)**kappa

!!!     tv_atm  = t_atm  * (1.0 + d608*q_atm)     ! virtual temperature
tv_atm  = t_atm  * (1.0 )     ! virtual temperature
th_atm  = t_atm  * p_ratio                ! potential T, using p_surf as refernce
thv_atm = tv_atm * p_ratio                ! virt. potential T, using p_surf as reference 
!!!     thv_surf= t_surf0 * (1.0 + d608*q_surf0 ) ! surface virtual (potential) T
thv_surf= t_surf0 * (1.0  ) ! surface virtual (potential) T
!     thv_surf= t_surf0                        ! surface virtual (potential) T -- just for testing tun off the q_surf

u_surf= 0.0
v_surf= 0.0
u_dif = u_surf - u_atm                    ! velocity components relative to surface
v_dif = v_surf - v_atm


if(alt_gustiness) then
w_atm = max( sqrt(u_dif**2 + v_dif**2), gust_const)
else
w_gust = max(gust,gust_min) ! minimum gustiness
w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + w_gust*w_gust)
endif

! derivatives of surface wind w.r.t. atm. wind components
dw_atmdu = u_dif/w_atm
dw_atmdv = v_dif/w_atm

#ifdef fv3_turb
if( do_mo_drag ) then  !  monin-obukhov similarity theory 
call mo_drag (thv_atm, thv_surf, z_atm, rough_mom, rough_heat,      &
            rough_moist, w_atm, cd_m, cd_t, cd_q, u_star, b_star )
else
#endif
cd_m=  cd_drag_cnst
cd_t=  cd_drag_cnst
cd_q=  cd_drag_cnst
u_star= 0.0
b_star= 0.0
!endif

! surface layer drag coefficients
drag_t = cd_t * w_atm
drag_q = cd_q * w_atm
drag_m = cd_m * w_atm

! density
rho = p_atm / (rdgas * tv_atm)  

! sensible heat flux
rho_drag = cp_air * drag_t * rho * sfc_heat_flx_amp
flux_t = rho_drag * (t_surf0 - th_atm)  ! flux of sensible heat (W/m**2)
dhdt_surf =  rho_drag                   ! d(sensible heat flux)/d(surface temperature)
dhdt_atm  = -rho_drag*p_ratio           ! d(sensible heat flux)/d(atmos temperature)

! evaporation
rho_drag  =  drag_q * rho
flux_q    =  rho_drag * (q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))

dedq_surf = rho_drag
dedt_surf = 0

dedq_atm  = -rho_drag   ! d(latent heat flux)/d(atmospheric mixing ratio)

q_star = flux_q / (u_star * rho)             ! moisture scale
q_surf = q_atm + flux_q / (rho*cd_q*w_atm)   ! surface specific humidity


if( decouple_surface ) then   !  Eliminate sensible heat transfer with the surface
flux_t = 0.0  
dhdt_atm  = 0.0
flux_q= 0.0
dedq_atm= 0.0
endif 


! upward long wave radiation
!       This really needs to incorporate the influence of non-unity emissivity
flux_r    =   stefan*t_surf**4               ! (W/m**2)
drdt_surf = 4*stefan*t_surf**3               ! d(upward longwave)/d(surface temperature)

! stresses
rho_drag   = drag_m * rho
flux_u     = rho_drag * u_dif   ! zonal      component of stress (Nt/m**2)
flux_v     = rho_drag * v_dif   ! meridional component of stress 


!        calculate d(stress component)/d(atmos wind component)
dtaudu_atm = 0.0
dtaudv_atm = 0.0
if (old_dtaudv) then
    dtaudv_atm = -rho_drag
    dtaudu_atm = -rho_drag
else
    dtaudu_atm = -cd_m*rho*(dw_atmdu*u_dif + w_atm)
    dtaudv_atm = -cd_m*rho*(dw_atmdv*v_dif + w_atm)
endif

end subroutine surface_flux_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine surface_flux_2d ( lon, lat,                           &
     t_atm,     q_atm_in,   u_atm,  v_atm, p_atm, z_atm,         &
     p_surf,    t_surf,     q_surf,                              &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,      &
     flux_t, flux_q, flux_r, flux_u, flux_v,                     &
     cd_m,      cd_t,       cd_q,                                &
     w_atm,     u_star,     b_star,     q_star,                  &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,               &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,              &
     is, js, dt, Time )
!
! calculate 2d surface fluxes by looping over 1d subroutine

!  arguments 
integer, intent(in)                        :: is, js
real, intent(in)                        :: dt
type(time_type), intent(in)                  :: Time
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:)     :: lat


real, intent(in),  dimension(:,:) ::                        &
    t_atm,     q_atm_in,   u_atm,     v_atm,               &
    p_atm,     z_atm,                                      &
    p_surf,    t_surf,                                     &
    rough_mom, rough_heat, rough_moist, rough_scale, gust

real, intent(out), dimension(:,:) ::                       &
    flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
    cd_m,      cd_t,       cd_q,                          &
    dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
    dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
    w_atm,     u_star,     b_star,    q_star

real, intent(inout), dimension(:,:) :: q_surf


!  local vars 
integer :: j

do j = 1, size(t_atm,2)
    call surface_flux_1d ( lon(:,j), lat(:,j),                              &
        t_atm(:,j),     q_atm_in(:,j),   u_atm(:,j),  v_atm(:,j),  p_atm(:,j),  z_atm(:,j),  &
        p_surf(:,j),    t_surf(:,j),     q_surf(:,j),                                        &
        rough_mom(:,j), rough_heat(:,j), rough_moist(:,j), rough_scale(:,j), gust(:,j),      &
        flux_t(:,j),    flux_q(:,j),     flux_r(:,j),    flux_u(:,j),    flux_v(:,j),        &
        cd_m(:,j),      cd_t(:,j),       cd_q(:,j),                                          &
        w_atm(:,j),     u_star(:,j),     b_star(:,j),     q_star(:,j),                       &
        dhdt_surf(:,j), dedt_surf(:,j),  dedq_surf(:,j),  drdt_surf(:,j),                    &
        dhdt_atm(:,j),  dedq_atm(:,j),   dtaudu_atm(:,j), dtaudv_atm(:,j),                    &
        is, js, dt, Time  )
end do



end subroutine surface_flux_2d



subroutine surface_flux_init
! 
!  Initialization of the surface flux module--reads the nml.     
!

!  local vars 
  integer :: unit, ierr, io

! read namelist
if ( file_exists('input.nml')) then
    unit = open_namelist_file ()
    ierr=1; 
    do while (ierr /= 0)
        read  (unit, nml=surface_flux_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'surface_flux_nml')
    enddo
10   call close_file (unit)
endif

! write version number
call write_version_number(version, tagname)

if ( mpp_pe() == mpp_root_pe() )  write (stdlog(), nml=surface_flux_nml)

if(.not. use_virtual_temp) d608 = 0.0

do_init = .false.

end subroutine surface_flux_init




end module surface_flux_mod

