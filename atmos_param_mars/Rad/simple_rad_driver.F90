module simple_rad_driver_mod
! module for simple radiation
use constants_mod, only: KAPPA, CP_AIR, GRAV, STEFAN, PI, RADIAN,      &
                         diffac, seconds_per_day,  RAD_TO_DEG 
 

use fms_mod, only: error_mesg, FATAL,                      &
                   open_namelist_file, check_nml_error,                &
                   mpp_pe, mpp_root_pe, close_file,                    &
                   write_version_number, stdlog,                       &
                   uppercase, read_data, write_data, field_size

use fms2_io_mod,            only:  file_exists
use time_manager_mod, only: time_type, get_time
use diag_manager_mod, only: register_diag_field, send_data

use radiation_util_mod, only : scatter_calc

use astronomy_mod,  only  :   astronomy_init,  mars_calender, sol88, zenith, &
                              solar_constant


implicit none
private

!------------------- interfaces ---------------------------------------

public :: simple_rad_driver,                                  &
          simple_rad_driver_init, simple_rad_driver_end

!-------------------- namelist -----------------------------------------

logical :: use_forget_swheat= .false.
integer :: rad_calc_intv = 1800 

logical ::  do_diurnal_avg_rad = .false. 

logical :: use_newton_damping = .false.
real    :: tdamp = 1.0    !     sols 

integer :: fixed_calender_date = -1

logical :: fixed_dust = .true.
logical :: do_lat_vary_dust = .false.
real    :: optical_depth_inpt = 0.2  

real    :: qext_cnst= 2.5         
real    :: sscat_cnst= 0.92       
real    :: gfac_cnst= 0.55        
real    :: conrath= 0.001         

real    :: tauir= 1.0 
real    :: alfir= 0.0 
integer :: kwave= -1

real    ::  radscale = 1.0
real    ::  rad_amp  = 0.0 
real    ::  rad_vscale = 2.0
logical ::  modulate_heating= .true. 

real    :: reference_press = 610.0 
    
logical :: no_atmos_dn_ir_flx= .false.    !   Remove downward atmos IR flux at surface

real    :: day_zero = 0.0



namelist /simple_rad_driver_nml/ use_forget_swheat,                         &
                                 rad_calc_intv,                             &
                                 do_diurnal_avg_rad,                        &
                                 use_newton_damping, fixed_calender_date,   &
                                 optical_depth_inpt, sscat_cnst, gfac_cnst, &
                                 conrath, alfir, tauir,                     &
                                 do_lat_vary_dust,                          &
                                 no_atmos_dn_ir_flx, reference_press,       &
                                 kwave, day_zero, tdamp, radscale,          &
                                 rad_amp, rad_vscale, modulate_heating  

real, dimension(:,:,:),   allocatable  ::  heatrate 
real, dimension(:,:  ),   allocatable  ::  sfc_sw_flux, sfc_lw_flux  


real, parameter :: gcp= GRAV / CP_AIR

logical :: module_is_initialized = .false.

integer :: id_insol, id_swheat, id_lwheat, id_lwdust,                  &
           id_ir_flx, id_opac, id_tdt_rad, id_solar_flx, id_vis_od,    &
           id_vis_od_dust, id_cosz

real, parameter   :: missing_value = -1.e10
character(len=16), parameter :: model='simple_rad'

character(len=128) :: version='$Id: simple_rad_driver.F90,v 1.1.2.1.4.1 2012/02/23 20:04:33 rjw Exp $'
character(len=128) :: tagname='$Name: mars_feb2012_rjw $'

logical ::  mcpu0

integer ::  radiation_counter 


contains

!=====================================================================
subroutine simple_rad_driver ( is, js, lon, lat, dt, Time,                &
                               p_half, p_full, tsfc, albedo, sfc_emiss,   &
                               t, r, tdt, rdt, swfsfc, lwfsfc, cosz )             
!  driver subroutine for simple radiation
!=======================================================
integer, intent(in)  :: is, js
real,    intent(in)  :: dt
type(time_type), intent(in)             :: Time
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:)     :: lat
real, intent(in),    dimension(:,:,:)   :: p_half, p_full
real, intent(in),    dimension(:,:)     :: tsfc
real, intent(in),    dimension(:,:)     :: albedo
real, intent(in),    dimension(:,:)     :: sfc_emiss
real, intent(in),    dimension(:,:,:)   :: t
real, intent(in),    dimension(:,:,:,:) :: r
real, intent(inout), dimension(:,:,:,:) :: rdt
real, intent(inout), dimension(:,:,:)   :: tdt
real, intent(out),   dimension(:,:)     :: swfsfc, lwfsfc  
real, intent(out),   dimension(:,:)     :: cosz  

! Local Variables

real, dimension(size(t,1),size(t,2)) :: psfc, hang, frac, radin, coszro 

real, dimension(size(t,1),size(t,2),size(t,3)) :: hsw, heatra
real, dimension(size(t,1),size(t,2))           :: sfc_ir_flx, trans, outflx

real, dimension(size(t,1),size(t,2),size(t,3)) :: delp, tau, opac
real, dimension(size(t,1),size(t,2),size(t,3)) :: sscat, gfac
real, dimension(size(t,1),size(t,2))           :: vis_od

real, dimension(size(t,1),size(t,2),size(p_half,3)) :: tcol

logical :: used
real    :: pref 
real :: fjd, r_orbit, declin, areolat, slag, dhr                            
real :: rsolar, fact 
real :: secs, sols                                                          

integer :: days, seconds
integer :: ie, je, id, jd, kd, ntrace
integer :: i, j, k, idjd

real :: ramp, tlen, xbar, psr, psi, psr2, psi2, rnorm   
real, dimension(size(t,1),size(t,2)) :: sinkx, coskx, sin2kx, cos2kx
real, dimension(size(t,1)) :: ttsave

real, dimension(size(t,1),size(t,2)) :: plog, damp



mcpu0 = (mpp_pe() == mpp_root_pe()) 

id= size(t,1); jd= size(t,2); kd= size(t,3);  ntrace= size(r,4)

ie= is + id - 1 
je= js + jd - 1 

idjd= id*jd

!         First get model time and calculate Mars' seasonal date  
!                  
call get_time(Time, seconds, days)

secs= days * seconds_per_day + seconds
sols = secs / seconds_per_day
days = sols
fjd = sols - days                                                     

if ( fixed_calender_date >= 0 ) then
    call mars_calender( fixed_calender_date, fjd, r_orbit, declin, areolat )
else
    call mars_calender( days, fjd, r_orbit, declin, areolat )
endif 


rsolar= solar_constant / r_orbit**2

areolat = areolat * RADIAN
if( areolat < 0.0 )  areolat = areolat + 360.0
areolat = modulo(areolat,360.)

!           Next get hang, zenith angle, and frac,  given declin 
call sol88( declin, lat, hang, coszro, frac )  


if( do_diurnal_avg_rad  ) then
    !          We already have the diurnal averaged zenith angle  
    cosz(:,:)= coszro(:,:)                            
else
    slag= 0.0
    dhr= 0.5 * (rad_calc_intv/seconds_per_day)*pi                                     
    call zenith( fjd, declin, slag, lon, lat, hang, dhr, cosz, frac )  
endif

if (id_cosz > 0)  used = send_data ( id_cosz, cosz, Time, is, js)

!           Now calculate solar heating 

radin(:,:)= rsolar*frac(:,:)*cosz(:,:)

do k=1,kd
    delp(:,:,k)= p_half(:,:,k+1)-p_half(:,:,k)
enddo 
psfc(:,:)= p_half(:,:,kd+1)

rdt(:,:,:,:)= 0.0

pref= reference_press

call dust_distrib_fixed( is, js, areolat, lon, lat, optical_depth_inpt, &
                                     pref, p_full, p_half, delp, tau )

opac(:,:,:)= tau(:,:,:) / delp(:,:,:)
if (id_opac > 0) used = send_data ( id_opac, opac,  Time, is, js )


vis_od= 0.0
do k= 1, kd
    vis_od(:,:)= vis_od(:,:) + tau(:,:,k)
enddo
if (id_vis_od > 0)   used = send_data ( id_vis_od, vis_od, Time, is, js )


tcol(:,:,1:kd)= t(:,:,:)
tcol(:,:,kd+1)= tsfc(:,:)

!         Aerosol optical properties:  single scattering albedo and asymmetry factor
sscat(:,:,:)= sscat_cnst
gfac(:,:,:) = gfac_cnst

call swheat( p_full, p_half, radin, albedo, cosz,          &
                tau, sscat, gfac, hsw, trans, outflx )


call lwheat( p_full, p_half, tcol, tsfc, sfc_emiss, heatra, sfc_ir_flx )    


!            Combine lwave and shortwave fluxes 

tlen= 5.0   
if( sols .lt. (tlen) ) then 
    ramp= 0.5*( 1.0 - cos( (sols)*pi/tlen )  )
else
    ramp= 1.0
end if

hsw=   ramp*hsw
trans= ramp*trans


#ifdef SPEC_CORE

if( modulate_heating ) then
    hang(:,:)=  rad_amp * gcp*(1.0/(1.0+rad_vscale))*trans(:,:)/psfc(:,:)
else
    hang(:,:)=  rad_amp * gcp*(1.0/(1.0+rad_vscale))*trans(:,:)/reference_press
endif 

do k= 1, kd 
    hsw(:,:,k)= hang(:,:)*(p_full(:,:,k)/psfc(:,:))**rad_vscale  + hsw(:,:,k)
enddo 

hsw= radscale * hsw

if( kwave == 0 ) then 
    rnorm= 1.0/float(id)
    do j= 1, jd
        do k= 1, kd
            xbar= sum(  hsw(:,j,k) )*rnorm
            hsw(:,j,k)= hsw(:,j,k) - xbar
        enddo
    enddo 
elseif ( kwave > 0) then 
    rnorm= 2.0/float(id)
    coskx=  cos(    lon(:,:)*kwave )     
    sinkx=  sin(    lon(:,:)*kwave )
    cos2kx= cos( 2.*lon(:,:)*kwave )     
    sin2kx= sin( 2.*lon(:,:)*kwave )
    do j= 1, jd
        do k= 1, kd
            ttsave= hsw(:,j,k)
            psr=  sum( ttsave(:)*coskx (:,j) )*rnorm
            psi=  sum( ttsave(:)*sinkx (:,j) )*rnorm
            psr2= sum( ttsave(:)*cos2kx(:,j) )*rnorm
            psi2= sum( ttsave(:)*sin2kx(:,j) )*rnorm

            hsw(:,j,k)= psr*coskx(:,j) + psi*sinkx(:,j) + psr2*cos2kx(:,j) + psi2*sin2kx(:,j)
        enddo
    enddo 

else
!            use full heating field
endif 


if( use_newton_damping  )  then 
!             add newtonian damping 

    rnorm= 1.0/float(id)
    do k= 1, kd
        plog(:,:)= -log( p_full(:,:,k)/psfc(:,:) )
        damp(:,:)= 2.0 -  1.7*0.5*(1.0+tanh( (plog-5.0)/1.0 ) )
        do j= 1, jd
            xbar= sum( t(:,j,k) )*rnorm
            heatra(:,j,k)= -( t(:,j,k) - xbar )
        enddo
        heatra(:,:,k)= heatra(:,:,k)/( damp(:,:)*88775.0 )
    enddo 

endif 




#endif 


if (id_swheat > 0)  used = send_data ( id_swheat, hsw,    Time, is, js)
if (id_lwheat > 0)  used = send_data ( id_lwheat, heatra, Time, is, js)

heatra = heatra + hsw 
swfsfc = trans
lwfsfc = sfc_ir_flx

if( no_atmos_dn_ir_flx )  lwfsfc= 0.0

!         Update global arrays which can save radiation over multiple time steps

sfc_sw_flux(is:ie,js:je)=   swfsfc
sfc_lw_flux(is:ie,js:je)=   lwfsfc
heatrate   (is:ie,js:je,:)= heatra

tdt = heatra

if (id_ir_flx     > 0)  used = send_data ( id_ir_flx,    lwfsfc, Time, is, js)
if (id_solar_flx  > 0)  used = send_data ( id_solar_flx, swfsfc, Time, is, js)
if (id_tdt_rad    > 0)  used = send_data ( id_tdt_rad,   tdt,    Time, is, js)


end subroutine simple_rad_driver
 

!=======================================================================
!=======================================================================




subroutine simple_rad_driver_init ( nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )   
!  initialize simple radiation
!=======================================================================


integer, intent(in)                    :: nlon, mlat, nlevels
real,    intent(in),  dimension(:,:)   :: lonb, latb
real,    intent(in),  dimension(:,:)   :: lon , lat
integer, intent(in) :: axes(4)
type(time_type), intent(in) :: Time

integer :: unit, io, ierr
integer :: i, j, ie, je, id, jd, is, js 

character (len=128) :: filename, fieldname


id= size(lonb,1)-1; jd= size(latb,2)-1

is= 1;   js= 1;
ie= is + id - 1 
je= js + jd - 1 

mcpu0 = (mpp_pe() == mpp_root_pe()) 

radiation_counter= 0

!     ----- read namelist -----

if (file_exists('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
        read  (unit, nml=simple_rad_driver_nml, iostat=io, end=20)
        ierr = check_nml_error (io, 'simple_rad_driver_nml')
    enddo
20    call close_file (unit)
endif

!     ----- write version info and namelist to log file -----

call write_version_number (version,tagname)

if (mpp_pe() == mpp_root_pe())    write (stdlog(),nml=simple_rad_driver_nml)


allocate ( heatrate(is:ie,js:je,nlevels)  )
allocate ( sfc_sw_flux(is:ie,js:je)  )
allocate ( sfc_lw_flux(is:ie,js:je)  )

!       Get heating rates and surface fluxes from radiation restart file 
filename= 'INPUT/radiation.res.nc' 
if (file_exists( trim( filename ) ) ) then
    if(mcpu0) print *,'Reading  restart file:   ',  trim( filename )
    call read_data( trim( filename ), 'heatrate'   , heatrate)
    call read_data( trim( filename ), 'sfc_sw_flux', sfc_sw_flux)
    call read_data( trim( filename ), 'sfc_lw_flux', sfc_lw_flux)
else
    heatrate= 0.0
    sfc_sw_flux= 0.0
    sfc_lw_flux= 0.0
endif


if(mcpu0) print *, 'Initializing heating  ',  is, ie, js, je 


!     ----- register diagnostic fields -----

id_insol = register_diag_field ( model, 'insol', axes(1:2),          &
                                       Time, 'insolation', 'W/m2', &
                                      missing_value=missing_value )

id_ir_flx = register_diag_field ( model, 'sfcirflx', axes(1:2),      &
                                 Time, 'downward IR flux', 'W/m2', &
                                       missing_value=missing_value )

id_cosz = register_diag_field ( model, 'cosz', axes(1:2),          &
                                Time, 'zenith angle', 'degrees', &
                                      missing_value=missing_value )

id_solar_flx = register_diag_field ( model, 'sfcswflx', axes(1:2),   &
                                 Time, 'downward SW flux', 'W/m2', &
                                       missing_value=missing_value )

id_swheat = register_diag_field ( model, 'swheat', axes(1:3),        &
                                       Time, 'sw heating', 'W/m2', &
                                      missing_value=missing_value )

id_lwheat = register_diag_field ( model, 'lwheat', axes(1:3),        &
                                       Time, 'lw heating', 'W/m2', &
                                      missing_value=missing_value )

id_lwdust = register_diag_field ( model, 'lwdust', axes(1:3),       &
                                  Time, 'lw dust heating', 'W/m2', &
                                      missing_value=missing_value )

id_opac = register_diag_field ( model, 'opac', axes(1:3),            &
                               Time, 'dust opacity', 'normalized', &
                                      missing_value=missing_value )
id_vis_od = register_diag_field ( model, 'vis_od', axes(1:2),      &
                               Time, 'visible opacity', 'total', &
                                      missing_value=missing_value )

id_vis_od_dust = register_diag_field ( model, 'vis_od_dust', axes(1:2), &
                               Time, 'Dust visible opacity', 'total', &
                                      missing_value=missing_value )

id_tdt_rad = register_diag_field ( model, 'tdt_rad', axes(1:3),      &
                    Time, 'radiative temperature tendency', 'K/s', &
                                      missing_value=missing_value )


module_is_initialized  = .true.

end subroutine simple_rad_driver_init

!=======================================================================
!=======================================================================

subroutine swheat( pp, pph, radin, albedo, coszro,        &
                   tau, sscat, gfac, hsw, trans, outflx )
!  calculate visible heating
!=======================================================================

real, intent(in), dimension(:,:)     :: radin   ! W/m/m
real, intent(in), dimension(:,:)     :: albedo  ! albedo
real, intent(in), dimension(:,:)     :: coszro  ! cosine of zenith angle

real, intent(in), dimension(:,:,:)   :: pp   ! kd 
real, intent(in), dimension(:,:,:)   :: pph  ! kp
real, intent(in), dimension(:,:,:)   :: tau, sscat, gfac

real, intent(out), dimension(:,:,:)  :: hsw    ! solar heating rate
real, intent(out), dimension(:,:)    :: trans  ! solar flux at surface
real, intent(out), dimension(:,:)    :: outflx ! upward visible flux at top of atmosphere 


!   LOCAL:

real, parameter :: co2heat0 = 1.3/88775.0
real, parameter :: p0nonlte = 7.5e-3

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: delp, delpi, heat, fscat

integer ::  kd, km, kp, id, jd, i, j, k, ichan 

id= size(pp,1); jd= size(pp,2); kd= size(pp,3)

!    ========calculate heating in the nir:  co2 contributions===========

if ( use_forget_swheat ) then  
    do k=1,kd
        hsw(:,:,k) = co2heat0*sqrt( coszro(:,:)*700.0/pp(:,:,k) )      &
                                      / (1.0 + p0nonlte/pp(:,:,k) )
    enddo
else
    hsw(:,:,:) = 0.0
endif 

!    =========calculate aerosol heating at solar wavelengths===

do k=1,kd
    delp (:,:,k)= pph(:,:,k+1)-pph(:,:,k)
    delpi(:,:,k)= 1.0 / ( pph(:,:,k+1)-pph(:,:,k) )
enddo

!        get aerosol optical properties distribution   
fscat(:,:,:)=  gfac(:,:,:) * gfac(:,:,:)

call scatter_calc( id*jd, kd, tau, sscat, gfac, fscat,      &
                coszro, albedo, heat, trans, outflx )

do k=1,kd
    hsw(:,:,k)= hsw(:,:,k) + gcp*radin(:,:)*heat(:,:,k)*delpi(:,:,k)
enddo

trans(:,:)= radin(:,:)*trans(:,:)

end subroutine swheat




!=======================================================================

subroutine lwheat( pp, pph, tcol, tsfc, sfc_emiss, heatra, flx_sfc )
!=======================================================================
!   Uses the Hourdin (1993) scheme to calculate CO2 heating/cooling rates
!    in the 15 micron band
!
!   Also includes a calculation of the influence of dust in the 
!   15 micron region


real, intent(in), dimension(:,:,:)   :: tcol  ! kd+1 values
real, intent(in), dimension(:,:,:)   :: pp    ! kd 
real, intent(in), dimension(:,:,:)   :: pph   ! kp

real, intent(in), dimension(:,:)     :: tsfc      ! surface temperature
real, intent(in), dimension(:,:)     :: sfc_emiss

real, intent(out), dimension(:,:,:)  :: heatra    ! heating rates 
real, intent(out), dimension(:,:)    :: flx_sfc   ! downward fluxes at surface
 

!   LOCAL:
real, dimension(size(pp,1),size(pp,2),size(pp,3))    :: delp, delpi

real, dimension(size(pph,1),size(pph,2),size(pph,3)) :: fup, fdn, net_flux
real, dimension(size(pph,1),size(pph,2),size(pph,3)) :: planck

real, dimension(size(pp,1),size(pp,2))               :: planck_sfc

real, dimension(size(pp,1),size(pp,2))               :: vfac, psfc

real, dimension(size(pph,1),size(pph,2),size(pph,3)) :: tau
real, dimension(size(pp,1),size(pp,2),size(pp,3))    :: dtrans 

real    :: tnot, pnot
integer :: id, jd, kd, kp, i, j, k 


id= size(pp,1); jd= size(pp,2); kd= size(pp,3)
kp= kd+1

do k=1,kd
    delp (:,:,k)= pph(:,:,k+1)-pph(:,:,k)
    delpi(:,:,k)= 1.0 / ( pph(:,:,k+1)-pph(:,:,k) )
enddo


!      planck function   (at layers) 
planck     = stefan*tcol*tcol*tcol*tcol
planck_sfc = stefan*tsfc*tsfc*tsfc*tsfc


pnot= reference_press

do k= 1, kp
    tau(:,:,k)= tauir * ( pph(:,:,k) / pnot )**(alfir+1.0)
enddo

do k = 1, kd
    dtrans(:,:,k) = exp( -(tau(:,:,k+1)-tau(:,:,k)) )
enddo

fup(:,:,kd+1) = planck_sfc(:,:)
do k = kd, 1, -1
    fup(:,:,k) = fup(:,:,k+1)*dtrans(:,:,k) + planck(:,:,k)*(1.0 - dtrans(:,:,k))
enddo

fdn(:,:,1) = 0.0
do k = 1, kd
    fdn(:,:,k+1) = fdn(:,:,k)*dtrans(:,:,k) + planck(:,:,k)*(1.0 - dtrans(:,:,k))
enddo

net_flux= fup - fdn

flx_sfc(:,:)= fdn(:,:,kp)

do k=1,kd
    heatra(:,:,k)= (net_flux(:,:,k+1)-net_flux(:,:,k))*delpi(:,:,k)*gcp 
enddo 

!           incorporate newtonian damping on zonal mean temp
!          i've increased damping timescale from 5.e4 to 2.e5
if( use_newton_damping  ) then 
    pnot= 600.e2 * 0.9447012000e-05
    tnot= 140.0
    do k= 1, 4
        vfac(:,:)= ( ( log(pp(:,:,k))/log(pnot) )**6  )/5.0e4

        heatra(:,:,k)= heatra(:,:,k) - vfac(:,:)*(tcol(:,:,k)-tnot)
    enddo
endif


return
end subroutine lwheat 

!   =================================================================
!   =================================================================


 subroutine dust_distrib_fixed( is, js, areolat, lon, lat, optical_depth, &
                                           pref, pp, pph, delp, tau )
! calculate fixed dust opacity
!   =================================================================

integer, intent(in)                 ::  is, js
real, intent(in)                    ::  areolat
real, intent(in), dimension(:,:)    ::  lon    !  in radians 
real, intent(in), dimension(:,:)    ::  lat    !  in radians 
real, intent(in)                    ::  optical_depth
real, intent(in)                    ::  pref
real, intent(in),  dimension(:,:,:) ::  pp   ! kd 
real, intent(in),  dimension(:,:,:) ::  pph  ! kp
real, intent(in),  dimension(:,:,:) ::  delp ! kd 
real, intent(out), dimension(:,:,:) ::  tau  ! opacity 


! evaluate optical depth between"flux" pressure levels
! Assume pressure increases with index number.
!
!    The parameter vscale (so-called Conrath paramter)
!         serves to pick at a height above which
!         the dust mixing ratio decreases stongly.

!    For vscale=0.01, this is at about 35 km.

!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt 
real, dimension(size(lat,1),size(lat,2))         :: rnorm, optical_depth2d
real, dimension(size(lat,1),size(lat,2))         :: zmax 

real     ::   a, zmaxx, fac1, delareo 
integer  ::  ie, je, id, jd, i, j, k, kd, iaa, iaam

kd= size(pp,3)
jd= size(pp,2)
id= size(pp,1)

ie= is + id - 1;   je= js + jd - 1

optical_depth2d(:,:)= optical_depth
a= conrath
zmaxx= 70.0
zmax(:,:)= 35.0 

if( do_lat_vary_dust ) then
    optical_depth2d(:,:)= optical_depth * ( 0.3 + 0.7*cos( lat(:,:) ) )
    zmax(:,:)= zmaxx - 40.0*sin( lat(:,:) )**2
endif


do k= 1, kd
    tau(:,:,k)=exp( a*(1.0 -  (pref/pp(:,:,k))**(zmaxx/zmax(:,:))   ) )   
    where( pp(:,:,k) > pref )
        tau(:,:,k)= 1.0
    end where
enddo

!         evaluate accumulated opacity at half-levels
opt(:,:,1)= 0.0
do k= 1, kd
    opt(:,:,k+1)= opt(:,:,k) + tau(:,:,k)*delp(:,:,k)
enddo

rnorm(:,:)= optical_depth2d(:,:) * pph(:,:,kd+1) / ( pref *opt(:,:,kd+1) )
do k= 1, kd+1
    opt(:,:,k)= opt(:,:,k)*rnorm(:,:)
enddo

!              Optical depth of individual layers
do k= 1, kd
    tau(:,:,k)= opt(:,:,k+1)-opt(:,:,k)
enddo



return
end subroutine dust_distrib_fixed

!   =================================================================
!   =================================================================


subroutine simple_rad_driver_end
!  release simple radiation memory
!   =================================================================

character (len=128) :: filename

mcpu0=  (mpp_pe() == mpp_root_pe()) 

deallocate ( heatrate  )
deallocate ( sfc_sw_flux )
deallocate ( sfc_lw_flux )   


return
end subroutine simple_rad_driver_end

!#######################################################################

end module simple_rad_driver_mod
