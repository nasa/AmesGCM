module aerosol_mod
!
!
!  Module to calculate binned aerosol opacities either from map or interactively from tracers
!

use constants_mod,      only: KAPPA, CP_AIR, RDGAS, GRAV, PI, RADIAN, RAD_TO_DEG

use fms_mod,            only: error_mesg, FATAL,                      &
                           open_namelist_file, check_nml_error,                &
                           mpp_pe, mpp_root_pe, close_file,                    &
                           write_version_number, stdlog,                       &
                           uppercase, read_data, write_data, field_size

use fms2_io_mod,            only:  file_exists
use time_manager_mod,   only: time_type, get_time
use diag_manager_mod,   only: register_diag_field, send_data,  diag_axis_init

use field_manager_mod,  only: MODEL_ATMOS, parse, find_field_index
use tracer_manager_mod, only: query_method, get_tracer_index,  &
                            get_number_tracers, get_tracer_names

use aerosol_util_mod, only: dust_map_scale_bin


implicit none

!------------------- interfaces ---------------------------------------

public ::      aerosol_init, aerosol_end,                      &
               aerosol_optics_init,                            &
               ndust_bins, nice_bins, nice_moms, aerosol_bins, &
               dust_indx, ice_bin_indx, ice_mom_indx,          &
               reff_dust, reff_ice,  radiative_active_inpt,    &
               dust_distrib_fix, dust_distrib_bin,             &
               do_cloud_rad
              
           
integer :: ndust_bins       !    number of dust aerosol bins 
integer :: nice_bins        !    number of ice aerosol bins 
integer :: nice_moms        !    number of ice moments  
integer :: aerosol_bins     !    number of aerosol bins 


real,     dimension(:),  allocatable  ::   reff_dust, reff_ice 
logical,  dimension(:),  allocatable  ::   radiative_active_inpt 
integer,  dimension(:),  allocatable  ::   dust_indx, ice_bin_indx, ice_mom_indx

!-------------------- namelist -----------------------------------------


logical :: fixed_dust = .true.              ! dust opacity from fixed map and analytical vertical dist
real    :: optical_depth_inpt = 0.0         ! background dust opacity

real    :: qext_cnst= 2.5                   ! dust extinction
real    :: sscat_cnst= 0.90                 ! dust single scattering
real    :: gfac_cnst= 0.65                  ! dust asymmetry parameter
real    :: conrath= 0.003         

integer :: conrath_type  = 1                ! 0 for orig conrath, 1 for new with zmax !SJG

real    :: optical_depth_pulse= 0.0         ! opacity for dust pulse
real    :: optical_depth_pulse_var= 0.0
real    :: pulse_width=   1.0               ! pulse depth in log(sigma) space
real    :: pulse_center= 30.0               ! pressure level for pulse center [mb]

logical :: do_interactive_dust_rad= .false. ! dust opacity from advected bin tracers
logical :: do_assimilated_dust = .false.    ! dust opacity from bin tracers and optional map
logical :: do_inpt_dust_cycle = .false.     ! use input dust map
logical :: do_inpt_dust_cycle_ktop = .false.! use advected tracer for dust top
integer :: inpt_dust_cycle_dust_index  = 1  ! bin tracer number for dust cycle
real    :: scale_inpt_dust_column = 1.0     ! scaling factor for dust opacity
integer :: dust_cycle_scheme = 0            ! 1 = scale dust column to zmax from map. 2 = no zmax scaling
logical :: do_lat_vary_dust = .false.       ! do LMD style latitudinal variation
integer :: irad_tracer = 1                  ! tracer number that is radiatively active
integer :: nrad_tracers = 1                 ! number of tracers that are active. WARNING should be consistent with clouds_physics_nml
integer :: dust_rad_size= 15                ! dust particle radius *10 mm

logical :: do_cloud_rad = .false.           ! bin cloud radiation
integer :: rad_ice_scheme = 6               ! different latitude distributions for ice
integer :: special_fixed_dust_case = 0
real    :: scale_rad_ice= 0.3               ! scaling factor for total ice opacity
real    :: scale_rad_ice_np= 0.1            ! scaling factor for north ice opacity
real    :: scale_rad_ice_sp= 0.1            ! scaling factor for south ice opacity
real    :: sscat_sw_ice= 0.995               ! visible single scattering albedo for ice
real    :: gfac_sw_ice= 0.80                ! visible asymmetry parameter for ice

logical :: do_spectrum_ir = .false.          ! split dust opacity to channels

logical :: reduce_low_level_dust = .false.  ! exponentially reduce low level dust opacity

logical :: add_storm_opacity = .false.      ! add dust storm
real    :: storm_sigx= 60.0                 ! storm longitude width
real    :: storm_sigy= 20.0                 ! storm latitude width
real    :: storm_lonc= 10.0                 ! storm longitude center
real    :: storm_latc= 15.0                 ! storm latitude center
real    :: add_storm_amp = 0.0              ! storm center opacity
real    :: storm_areo = 180.0               ! storm center timing


namelist /aerosol_nml/ fixed_dust, optical_depth_inpt, conrath, conrath_type, &
                        do_inpt_dust_cycle, &
                        inpt_dust_cycle_dust_index,                                  &
                        do_lat_vary_dust, do_inpt_dust_cycle_ktop,                   &
                        dust_cycle_scheme, do_interactive_dust_rad,                  &
                        special_fixed_dust_case,                                     &
                        do_assimilated_dust, optical_depth_pulse,                    &
                        pulse_width, pulse_center,                                   &
                        add_storm_opacity, storm_sigx, storm_sigy,                   &
                        storm_lonc, storm_latc,                                      &
                        add_storm_amp,  storm_areo, optical_depth_pulse_var




namelist /aerosol_optics_nml/ qext_cnst, sscat_cnst, gfac_cnst,             &
                            irad_tracer,                                     &
                            dust_rad_size,                                   &
                            sscat_sw_ice,  gfac_sw_ice, nrad_tracers,        &
                            do_cloud_rad,  scale_rad_ice,                    &
                            scale_rad_ice_np, scale_rad_ice_sp,              &
                            rad_ice_scheme,                                  &
                            do_spectrum_ir, scale_inpt_dust_column,          &
                            reduce_low_level_dust



real, parameter :: qex_ref= 2.5                             ! reference dust extinction

!     Fields for specified dust cycle
real, dimension(:),     allocatable :: areo_cycle           ! dust cycle map Ls
real, dimension(:,:,:), allocatable :: zmax_cycle           ! dust cycle map zmax
real, dimension(:,:,:), allocatable :: tau_cycle            ! dust cycle map opacity
real, dimension(:,:,:), allocatable :: taufill_cycle        ! dust cycle filled map opacity


!     The following are specific to the gfdl radiation scheme
!     and are those used in radiation_driver_mod
integer, parameter :: nchan_lw = 3                          ! number of IR channels
integer, parameter :: nchan_sw = 1                          ! number of visible channels
integer, parameter :: nradbands= nchan_sw + nchan_lw + 1    ! number of radiation bands

real, dimension(nchan_lw+1) :: qex_lw_dust, sscat_lw_dust, gfac_lw_dust
real, dimension(nchan_lw+1) :: qex_lw_ice,  sscat_lw_ice,  gfac_lw_ice

real, parameter   :: missing_value = -1.e10
character(len=7) :: mod_name = 'aerosol'

logical :: module_is_initialized = .false.

character(len=128) :: version='$Id: aerosol.F90,v 1.1.2.2 2011/11/22 18:33:07 rjw Exp $'
character(len=128) :: tagname='$Name: mars_feb2012_rjw $'

logical ::  mcpu0


contains

!#######################################################################




 subroutine aerosol_init ( nlon, mlat, lonb, latb, lon, lat, axes, Time )
!
!   Initialize dust cycle aand FV3 aerosol optical properties
!
!

integer, intent(in)                   :: nlon, mlat
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon, lat

integer, intent(in) :: axes(4)
type(time_type), intent(in) :: Time

character (len=128) :: filename, fieldname, tracer_name

integer ::  unit, io, ierr, n, nt
integer ::  ntrace, ntprog, ntdiag, ntfam, ntdust, ntice, ntmom

character(len=128) :: scheme, params

real    ::   value
real    ::   values(3)

!     ----- read namelist /aerosol_nml/-----

if (file_exists('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1; do while (ierr /= 0)
        read  (unit, nml=aerosol_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'aerosol_nml')
    enddo
10  call close_file (unit)
endif

mcpu0 = (mpp_pe() == mpp_root_pe())

call write_version_number (version,tagname)
if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=aerosol_nml)


call get_number_tracers (MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfam )

aerosol_bins= 0
ndust_bins= 0
nice_bins= 0
nice_moms= 0

do n = 1, ntrace
    if (query_method('aerosol', MODEL_ATMOS, n, scheme, params)) then
        if(mcpu0)   print *, 'Aerosol Tracer', n, trim(scheme), trim(params)
        aerosol_bins= aerosol_bins + 1

        if( trim(scheme)== 'ice_bin' ) nice_bins= nice_bins + 1
        if( trim(scheme)== 'dust_bin' ) ndust_bins= ndust_bins + 1
        if( trim(scheme)== 'ice_mom' ) nice_moms= nice_moms + 1
    endif
enddo


if ( ndust_bins > 0 ) then
    if(mcpu0)  print *, 'Allocating dust bins:  ', ndust_bins
    allocate ( dust_indx(ndust_bins) )
    allocate ( radiative_active_inpt(ndust_bins) )
    allocate ( reff_dust(ndust_bins) )

    reff_dust(:)= 2.5
    radiative_active_inpt(:)= .false.
endif

if ( nice_bins  > 0 ) then
    if(mcpu0)   print *, 'Allocating ice bins  :  ', nice_bins
    allocate ( ice_bin_indx(nice_bins) )
    allocate ( reff_ice(nice_bins) )

    reff_ice(:)= 5.0
endif

if ( nice_moms  > 0 ) allocate ( ice_mom_indx(nice_moms) )


ntdust= 0
ntice= 0
ntmom= 0
do n = 1, ntrace
    if (query_method('aerosol', MODEL_ATMOS, n, scheme, params)) then
        if(mcpu0)  print *, 'Aerosol Tracer', n, trim(scheme), trim(params)

        call get_tracer_names(MODEL_ATMOS, n, tracer_name)

        if(mcpu0) print *, 'Field', n, trim(tracer_name), '  ', trim(scheme), &
                                   '  ', trim(params)

        if( trim(scheme)== 'ice_bin' ) then
            ntice= ntice + 1
            ice_bin_indx(ntice)= n
            if (parse(params,'radius',value) == 1) then
                reff_ice(ntice)= value
            endif
        endif

        if( trim(scheme)== 'dust_bin' ) then
            ntdust= ntdust + 1
            dust_indx(ntdust)= n
            if(mcpu0)  print *, 'Setting dust index ',  trim(params)
            if (parse(params,'radius',value) == 1) then
                reff_dust(ntdust)= value
            endif
        endif

        if( trim(scheme)== 'ice_mom' ) then
            ntmom= ntmom + 1
            ice_mom_indx(ntmom)= n
        endif
    endif
enddo



if(mcpu0) print *, '---Aerosol Summary:  tracer_index, radius, Radiatively_active ------'

do n= 1, ndust_bins
    if(mcpu0) print *, 'Dust Bin:  ', dust_indx(n), reff_dust(n), radiative_active_inpt(n)
    reff_dust(n)= reff_dust(n)*1.E-6   !   convert from microns to meters
enddo

do n= 1, nice_bins
    if(mcpu0) print *, 'Ice Bin:  ', ice_bin_indx(n), reff_ice(n)
    reff_ice(n)= reff_ice(n)*1.E-6     !   convert from microns to meters
enddo

call aerosol_optics_init( nlon, mlat, lonb, latb, lon, lat, axes, Time )

module_is_initialized  = .true.

end subroutine aerosol_init


!   =========================================================================
!   =========================================================================




subroutine aerosol_optics_init ( nlon, mlat, lonb, latb, lon, lat, axes, Time )
!
!
!    routine for initializing the radiation and aerosol tables
!
!

use horiz_interp_mod,  only: horiz_interp, horiz_interp_init, horiz_interp_new,  &
                          horiz_interp_end,  horiz_interp_del, horiz_interp_type

integer,     intent(in)                    :: nlon, mlat
real,        intent(in),   dimension(:,:)  :: lonb, latb
real,        intent(in),   dimension(:,:)  :: lon
real,        intent(in),   dimension(:,:)  :: lat
integer,     intent(in)                    :: axes(4)
type(time_type), intent(in)                :: Time


!     local variables

integer :: unit, io, ierr
integer :: ie, je, id, jd, k, iaa, is, js

character (len=128) :: filename, fieldname

type(horiz_interp_type)  ::  Interp

real,  dimension(:),     allocatable  :: areoz_inpt
real,  dimension(:,:,:), allocatable  :: zmax_inpt
real,  dimension(:,:,:), allocatable  :: tau2d_inpt
real,  dimension(:,:,:), allocatable  :: taufill_inpt

real, dimension(:),  allocatable  ::   lon_inpt,  lat_inpt
real, dimension(:),  allocatable  ::   lonb_inpt, latb_inpt

real, dimension(2)   ::   lonbz  =  (/ 0.0, 6.2832 /)

integer   ::   id_inpt, jd_inpt, areo_length, fld_dims(4)


!

id= size(lonb,1)-1; jd= size(latb,2)-1

is= 1;   js= 1;
ie= is + id - 1
je= js + jd - 1

!     ----- read namelist -----

mcpu0 = (mpp_pe() == mpp_root_pe())

if( mcpu0 )  print *, 'Aerosol dimensions  :',  id, jd

if (file_exists('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
        read  (unit, nml=aerosol_optics_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'aerosol_optics_nml')
    enddo
10    call close_file (unit)
endif


if( do_inpt_dust_cycle_ktop .AND. dust_cycle_scheme < 1 )  &
call error_mesg ('aerosol_optics_init',                &
    'do_inpt_dust_cycle_ktop and dust_cycle_scheme are wrong', FATAL)


!     ----- write version info and namelist to log file -----

call write_version_number (version,tagname)
if (mpp_pe() == mpp_root_pe())  write (stdlog(),nml=aerosol_optics_nml)


!   ------------Read input dust cycle data--------------

if( do_inpt_dust_cycle  ) then

    if(mcpu0)  print *,  'Reading input dust cycle data from aerosol_init  '

    filename= 'INPUT/dust_cycle.nc'
    if( file_exists( trim( filename ) ) ) then

        call field_size( trim(filename), 'areo', fld_dims )
        areo_length= fld_dims(1)
        call field_size( trim(filename), 'lon', fld_dims )
        id_inpt= fld_dims(1)
        call field_size( trim(filename), 'lat', fld_dims )
        jd_inpt= fld_dims(1)

        if(mcpu0)  print *,  'Reading input dust cycle data from aerosol_init:   ', jd_inpt, id_inpt

        allocate (  areoz_inpt(areo_length)  )
        call read_data( trim(filename), 'areo', areoz_inpt, no_domain=.true. )
        if(mcpu0) print *, 'Have read tes areo data file in aerosol.F: ', areoz_inpt(1:3)

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

        allocate (  zmax_inpt   (id_inpt,jd_inpt,areo_length)  )
        allocate (  tau2d_inpt  (id_inpt,jd_inpt,areo_length)  )
        allocate (  taufill_inpt(id_inpt,jd_inpt,areo_length)  )


        if(mcpu0) print *, 'Prepare read tes tau data from aerosol_init: '
        call read_data( trim(filename), 'tau', tau2d_inpt, no_domain=.true. )
        if(mcpu0) print *, 'Have read tes tau data file from aerosol_init: '

        call read_data( trim(filename), 'zmax', zmax_inpt, no_domain=.true. )
        if(mcpu0) print *, 'Have read tes zmax data file from aerosol_init: '

        call read_data( trim(filename), 'taufill', taufill_inpt, no_domain=.true. )
        if(mcpu0) print *, 'Have read tes taufill data file: from aerosol_init '

        allocate (  areo_cycle(               areo_length)  )
        allocate (  zmax_cycle   (is:ie,js:je,areo_length)  )
        allocate (  tau_cycle    (is:ie,js:je,areo_length)  )
        allocate (  taufill_cycle(is:ie,js:je,areo_length)  )


!                            Carry out horizontal interpolation
!       Use (lon,lat) --> bilinear and   (lonb,latb) --> conservative

        call horiz_interp_init
        call horiz_interp_new( Interp,lonb_inpt,latb_inpt,lon ,lat ,interp_method= 'bilinear' )

        do iaa= 1, areo_length
            call horiz_interp( Interp,   tau2d_inpt(:,:,iaa),     tau_cycle(is:ie,js:je,iaa) )
            call horiz_interp( Interp, taufill_inpt(:,:,iaa), taufill_cycle(is:ie,js:je,iaa) )
            call horiz_interp( Interp,    zmax_inpt(:,:,iaa),    zmax_cycle(is:ie,js:je,iaa) )
        enddo

        call horiz_interp_del( Interp )


        if(mcpu0)  print *, 'Maxval tau_cycle, zmax   ',  maxval( taufill_cycle ), maxval( zmax_cycle )
        if(mcpu0)  print *, 'Minval tau_cycle, zmax   ',  minval( taufill_cycle ), minval( zmax_cycle )

!        ! IR VS VIS opacities
!        if (dgdm_type.eq.2) then
!            taufill_cycle=taufill_cycle/2.75
!        else
            taufill_cycle=taufill_cycle*dust_map_scale_bin/2.75
!        endif

        areo_cycle(:)=  areoz_inpt(:)

        if(mcpu0)  print *, 'Successful setup for dust cycle:', is, js, je, areo_length

        deallocate ( areoz_inpt, zmax_inpt, tau2d_inpt, taufill_inpt  )
        deallocate ( lat_inpt, latb_inpt, lon_inpt, lonb_inpt )

    endif

endif                      !            End of dust cycle input


module_is_initialized  = .true.

end subroutine aerosol_optics_init


!   =========================================================================
!   =========================================================================




subroutine dust_distrib_simple( is, js, areolat, lon, lat, optical_depth, &
                                                 pref, pp, pph, delp, tau )
!
!  find the dust opacity distribution
!

integer, intent(in)                 ::  is, js
real, intent(in)                    ::  areolat
real, intent(in), dimension(:,:)    ::  lon    !  in radians
real, intent(in), dimension(:,:)    ::  lat    !  in radians
real, intent(in)                    ::  optical_depth
real, intent(in)                    ::  pref
real, intent(in), dimension(:,:,:)  ::  pp   ! kd
real, intent(in), dimension(:,:,:)  ::  pph  ! kp
real, intent(in), dimension(:,:,:)  ::  delp   ! kd
real, intent(out), dimension(:,:,:) ::  tau  ! opacity

!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt
real, dimension(size(lat,1),size(lat,2))         :: rnorm, optical_depth2d
real, dimension(size(lat,1),size(lat,2))         :: zmax , optx

real   ::   a, zmaxx, fac1, delareo
integer  ::  ie, je, id, jd, i, j, k, kd, iaa, iaam

kd= size(pp,3)
jd= size(pp,2)
id= size(pp,1)

ie= is + id - 1;   je= js + jd - 1

optical_depth2d(:,:)= optical_depth * cos( lat(:,:) )**2
a= conrath
zmaxx= 70.0
zmax(:,:)= 35.0


if( do_inpt_dust_cycle ) then
    call season_table( areolat, is, ie, js, je, optical_depth2d, zmax  )
    zmax= min( zmax, 60.0 )

    if ( add_storm_opacity ) then
        call add_storm( areolat, lon, lat, optx );
        optical_depth2d= optical_depth2d + optx
    endif

    do k= 1, kd
        if(conrath_type .eq. 0) then
            tau(:,:,k)=exp( a*(1.0 -  (pref/pp(:,:,k))   ) )
        else
            ! lmd-type vertical dust distribution
            tau(:,:,k)=exp( 0.007*(1.0 -  (700.0/pp(:,:,k))**(zmaxx/zmax(:,:))   ) )
        end if
    enddo

elseif( do_lat_vary_dust ) then
    optical_depth2d(:,:)= optical_depth * ( 0.3 + 0.7*cos( lat(:,:) ) )
    zmax(:,:)= zmaxx - 40.0*sin( lat(:,:) )**2

    do k= 1, kd
        tau(:,:,k)=exp( a*(1.0 -  (pref/pp(:,:,k))**(zmaxx/zmax(:,:))   ) )  !
    enddo
else
    do k= 1, kd
        tau(:,:,k)=exp( a*(1.0 -  (pref/pp(:,:,k))   ) )  !
    enddo
endif         !           -------------- dust cycle look-up ----------



if( reduce_low_level_dust ) then
    fac1= 1.0/float(kd)
    do k= 1, kd
        tau(:,:,k)= tau(:,:,k)*exp( -3.0*( fac1*k )**6 )
    enddo
endif

!         evaluate accumulated opacity at half-levels
opt(:,:,1)= 0.0
do k= 1, kd
    opt(:,:,k+1)= opt(:,:,k) + tau(:,:,k)*delp(:,:,k)
enddo

rnorm(:,:) = 0.
where (opt(:,:,kd+1) .gt. 0.)  rnorm(:,:)= optical_depth2d(:,:) * pph(:,:,kd+1) / ( pref *opt(:,:,kd+1) )
do k= 1, kd+1
    opt(:,:,k)= opt(:,:,k)*rnorm(:,:)
enddo

!              optical depth of individual layers
do k= 1, kd
    tau(:,:,k)= opt(:,:,k+1)-opt(:,:,k)
enddo


return
end subroutine dust_distrib_simple

!   =========================================================================
!   =========================================================================


 subroutine add_storm( areolat, lon, lat, optical_depth )
!
!   add dust storm opacity centered in time and space
!

real, intent(in)                                     ::  areolat
real, intent(in), dimension(:,:)          ::  lon    !  in radians
real, intent(in), dimension(:,:)          ::  lat    !  in radians
real, intent(out),dimension(:,:)        ::  optical_depth

integer  ::  id, jd

real, dimension(size(lat,1),size(lat,2))              :: anom, aa, aa2

real ::  slope,  tpi, rt2, sigma_x, sigma_y, theta, a, b, c,  lonc, lonc2, latc

jd= size(lon,2)
id= size(lon,1)

tpi= 2.0*pi

sigma_x= storm_sigx*pi/180
sigma_y= storm_sigy*pi/180

theta= slope*pi/180
lonc=  storm_lonc*pi/180
latc=  storm_latc*pi/180

if( storm_lonc   < 180 )  then
    lonc2= lonc + tpi
else
    lonc2= lonc - tpi
endif


rt2= sqrt(2.0)

a= ( cos(theta)  / ( rt2*sigma_x )  )**2  +  (sin(theta)/(rt2*sigma_y) )**2;

b=  -sin(theta)  / ( 2*sigma_x ) **2  +  sin(2*theta)/(2*sigma_y)**2;

aa(:,:) =    a*(lon(:,:)-lonc  )**2 + 2.0*b*(lon(:,:)-lonc  )*(lat(:,:)-latc) + c*(lat(:,:)-latc)**2
aa2(:,:) =    a*(lon(:,:)-lonc2)**2 + 2.0*b*(lon(:,:)-lonc2)*(lat(:,:)-latc) + c*(lat(:,:)-latc)**2

anom(:,:)= exp( -aa(:,:) )    + exp( -aa2(:,:) )


optical_depth= add_storm_amp * anom *  0.5*( 1 + tanh( (areolat-storm_areo)/1.0 )  )


return
end subroutine add_storm


!   =========================================================================
!   =========================================================================


subroutine dust_distrib_pulse( is, js, areolat, lon, lat, optical_depth, &
                                             pref, pp, pph, delp, tau )
!
!   add pulse of dust to atmosphere
!
integer, intent(in)                 ::  is, js
real, intent(in)                    ::  areolat
real, intent(in), dimension(:,:)    ::  lon    !  in radians
real, intent(in), dimension(:,:)    ::  lat    !  in radians
real, intent(in)                    ::  optical_depth
real, intent(in)                    ::  pref
real, intent(in), dimension(:,:,:)  ::  pp   ! kd
real, intent(in), dimension(:,:,:)  ::  pph  ! kp
real, intent(in), dimension(:,:,:)  ::  delp   ! kd
real, intent(out), dimension(:,:,:) ::  tau  ! opacity


! evaluate optical depth between"flux" pressure levels
! Assume pressure increases with index number.
!
!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt
real, dimension(size(lat,1),size(lat,2))              :: rnorm, optical_depth2d, sig

real   ::   a,  pwid, pulse_ctr, snot
integer  ::  ie, je, id, jd, i, j, k, kd


kd= size(pp,3)
jd= size(pp,2)
id= size(pp,1)

!  ie= is + id - 1;   je= js + jd - 1

optical_depth2d(:,:)= optical_depth * cos( lat(:,:) ) **2

snot= log( pulse_center/pref )

pwid= pulse_width

mcpu0 = (mpp_pe() == mpp_root_pe())

do k= 1, kd
    sig(:,:)= log(  pp(:,:,k)/pref  )
    tau(:,:,k)= exp( -( (sig(:,:) - snot)/pwid )**2   )
enddo

!  if( mcpu0 ) print *, 'have calculated initial tau in dust-pulse '

opt(:,:,1)= 0.0
do k= 1, kd
    opt(:,:,k+1)= opt(:,:,k) + tau(:,:,k)*delp(:,:,k)
enddo

!    since this is an elevated dust layer, dont normalize with surface pressure
rnorm(:,:) = 0.
where (opt(:,:,kd+1) .gt. 0.)  rnorm(:,:)= optical_depth2d(:,:)  / ( opt(:,:,kd+1) )
do k= 1, kd+1
    opt(:,:,k)= opt(:,:,k)*rnorm(:,:)
enddo

!              optical depth of individual layers
do k= 1, kd
    tau(:,:,k)= opt(:,:,k+1)-opt(:,:,k)
enddo

return
end subroutine dust_distrib_pulse

!   =========================================================================
!   =========================================================================


subroutine dust_distrib_storm( is, js, areolat, lon, lat, optical_depth, &
                                                 pref, pp, pph, delp, tau )
!
!  add a storm to the dust field
!

integer, intent(in)                 ::  is, js
real, intent(in)                    ::  areolat
real, intent(in), dimension(:,:)    ::  lon    !  in radians
real, intent(in), dimension(:,:)    ::  lat    !  in radians
real, intent(in)                    ::  optical_depth
real, intent(in)                    ::  pref
real, intent(in), dimension(:,:,:)  ::  pp   ! kd
real, intent(in), dimension(:,:,:)  ::  pph  ! kp
real, intent(in), dimension(:,:,:)  ::  delp   ! kd
real, intent(out), dimension(:,:,:) ::  tau  ! opacity

! evaluate optical depth between"flux" pressure levels
! Assume pressure increases with index number.
!
!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt
real, dimension(size(lat,1),size(lat,2))              :: rnorm, optical_depth2d, sig

real   ::   a,  pwid, pulse_ctr, snot
integer  ::  ie, je, id, jd, i, j, k, kd


kd= size(pp,3)
jd= size(pp,2)
id= size(pp,1)

!  ie= is + id - 1;   je= js + jd - 1

optical_depth2d(:,:)= optical_depth * cos( lat(:,:) ) **2

snot= log( pulse_center/pref )

pwid= pulse_width

mcpu0 = (mpp_pe() == mpp_root_pe())

do k= 1, kd
    sig(:,:)= log(  pp(:,:,k)/pref  )
    tau(:,:,k)= exp( -( (sig(:,:) - snot)/pwid )**2   )
enddo

!  if( mcpu0 ) print *, 'have calculated initial tau in dust-pulse '

opt(:,:,1)= 0.0
do k= 1, kd
    opt(:,:,k+1)= opt(:,:,k) + tau(:,:,k)*delp(:,:,k)
enddo

!    since this is an elevated dust layer, dont normalize with surface pressure
rnorm(:,:) = 0.
where (opt(:,:,kd+1) .gt. 0.)  rnorm(:,:)= optical_depth2d(:,:)  / ( opt(:,:,kd+1) )
do k= 1, kd+1
    opt(:,:,k)= opt(:,:,k)*rnorm(:,:)
enddo

!              optical depth of individual layers
do k= 1, kd
    tau(:,:,k)= opt(:,:,k+1)-opt(:,:,k)
enddo



return
end subroutine dust_distrib_storm


!   =========================================================================
!   =========================================================================

subroutine dust_distrib_bin( is, js, areolat, lon, lat, pref,                &
                                 pp, pph, delp,  aerosolmr,               &
                                 tau_spec,  sscat_spec, gfac_spec,        &
                                 dustref                     )
!
!   call binned dust distrib for ames radiation only.
!   should eventually be merged with dust_distrib
!
integer, intent(in)                   ::  is, js
real, intent(in)                      ::  areolat
real, intent(in), dimension(:,:)      ::  lon    !  in radians
real, intent(in), dimension(:,:)      ::  lat    !  in radians
real, intent(in)                      ::  pref
real, intent(in), dimension(:,:,:)    ::  pp
real, intent(in), dimension(:,:,:)    ::  pph
real, intent(in), dimension(:,:,:)    ::  delp          !
real, intent(in), dimension(:,:,:,:)  ::  aerosolmr   !

real, intent(out), dimension(:,:,:,:) ::  tau_spec   !
real, intent(out), dimension(:,:,:,:) ::  sscat_spec !
real, intent(out), dimension(:,:,:,:) ::  gfac_spec  !

real, intent(out), dimension(:,:,:)   ::  dustref   !

!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt
real, dimension(size(lat,1),size(lat,2))         :: rnorm, optical_depth2d
real, dimension(size(lat,1),size(lat,2))         :: zmax, ystruc

!      dimension radwt(ntrace), sscatwt(ntrace), gfacwt(ntrace)

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: dustrp, aerorp, cloudvis

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: tau, sscat, gfac

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: odfac, ssfac, ggfac

real, dimension(size(aerosolmr,4)  )    :: rdust_core
real, dimension(size(aerosolmr,4),1)    :: qex_sw, sscat_sw, gfac_sw

real, dimension(size(pp,1),size(pp,2),size(pp,3),size(aerosolmr,4)) :: taucc


real     ::   a, zmaxx, optical_depth
integer  ::  ie, je, id, jd, i, j, k, kd, nrad_trace, nt, ndx

integer  ::  iaa, iaam, nchan, lwchan, ird
real     ::  fac1, fact,  tau_thresh, qex_ref, qex_ice, reff

integer  ::   rad_tracers, nice

real     ::   rhoice=   0.917E3
real     ::   rhodust = 2.500E3

kd= size(pp,3)
jd= size(pp,2)
id= size(pp,1)

mcpu0 = (mpp_pe() == mpp_root_pe())

rdust_core(:)= 1.6E-6
qex_sw  (:,1)= qext_cnst
sscat_sw(:,1)= sscat_cnst
gfac_sw (:,1)= gfac_cnst

qex_ref= 2.5      !   for visible channel

tau  (:,:,:) = 0.0
sscat(:,:,:) = 0.0
gfac (:,:,:) = 0.0
!         fscat(:,:) = 0.0

tau_spec  (:,:,:,:)= 0.0
sscat_spec(:,:,:,:)= 0.0
gfac_spec (:,:,:,:)= 0.0

optical_depth= optical_depth_inpt


if( do_interactive_dust_rad )  then

    if( optical_depth .gt. 0.0 ) then    !  ---------------- background dust

        call dust_distrib_simple( is, js, areolat, lon, lat,  &
                      optical_depth,  pref, pp, pph, delp, tau )


!          Form weighted optical depth, single-scattering albedo etc

        tau_spec  (:,:,:,1)= tau(:,:,:)
        sscat_spec(:,:,:,1)= tau(:,:,:) * sscat_cnst
        gfac_spec (:,:,:,1)= tau(:,:,:) * sscat_cnst * gfac_cnst

    endif   !   --------------------  end background dust

!                Solar channel :   Dust only radiation:
!
!   ------------Begin loop over tracers  and accumulate dust --------------

    tau  (:,:,:)= tau_spec(:,:,:,1)
    sscat(:,:,:)= sscat_spec(:,:,:,1)
    gfac (:,:,:)= gfac_spec (:,:,:,1)



!    Accumulate contributions to opacity, sscat etc  from individual dustbins

    do  nt= 1, nrad_tracers

        if ( nrad_tracers > 1 )  then
            ndx= dust_indx(nt)
        else
            ndx= dust_indx(irad_tracer)
        endif

        !             note that qex_sw is currently independent of radius

        fact= 3.0 / ( 4.0*grav*rhodust*reff_dust(ndx) )

        aerorp(:,:,:) = aerosolmr(:,:,:,ndx) * fact * delp(:,:,:)

        odfac(:,:,:)     = qex_sw(nt,1) ! dust core only

        ssfac(:,:,:)= odfac(:,:,:)*sscat_sw(nt,1)
        ggfac(:,:,:)= ssfac(:,:,:)*gfac_sw (nt,1)

        tau(:,:,:)   = tau  (:,:,:) + aerorp(:,:,:)*odfac(:,:,:)
        sscat(:,:,:) = sscat(:,:,:) + aerorp(:,:,:)*ssfac(:,:,:)
        gfac(:,:,:)  = gfac (:,:,:) + aerorp(:,:,:)*ggfac(:,:,:)

    enddo


!             Normalization of sscat and gfac
    tau_thresh = 1.0e-12

    where( tau > tau_thresh )
        tau_spec  (:,:,:,1) = tau  (:,:,:)
        gfac_spec (:,:,:,1) = gfac (:,:,:) / sscat(:,:,:)
        sscat_spec(:,:,:,1) = sscat(:,:,:) / tau  (:,:,:)
    end where


    dustref(:,:,:)= tau_spec(:,:,:,1)


else   ! do assimilated dust if not interactive

    dustref(:,:,:)= 0.0


    if( optical_depth .gt. 0.0 ) then    !  ---------------- background dust
        call dust_distrib_simple( is, js, areolat, lon, lat,             &
                              optical_depth, pref, pp, pph, delp, tau )
        dustref(:,:,:)= tau(:,:,:)
    endif

!              Add in a "pulse" of dust
    if( optical_depth_pulse .gt. 0.0 ) then    !  ---------------- background dust
        call  dust_distrib_pulse( is, js, areolat, lon, lat,            &
                        optical_depth_pulse, pref, pp, pph, delp, tau )
        dustref(:,:,:)= dustref(:,:,:) + tau(:,:,:)
    endif


!        Passing tracers  ndx,ndx+1,ndx+2....ndx+nrad_tracers-1
    ndx= dust_indx(inpt_dust_cycle_dust_index)

    call dust_opacity_assim( is, js, areolat, lon, lat,                 &
                    aerosolmr(:,:,:,ndx:ndx+nrad_tracers-1),     &
                    pref, pp, pph, delp, tau,                    &
                    taucc(:,:,:,1:nrad_tracers)   )


    tau(:,:,:)= tau(:,:,:) * scale_inpt_dust_column
    taucc=      taucc      * scale_inpt_dust_column

!              Add in the background  and "pulse" dust contributions
    tau(:,:,:)= dustref(:,:,:) + tau(:,:,:)

    tau_spec  (:,:,:,1)= tau(:,:,:)
    sscat_spec(:,:,:,1)= sscat_cnst
    gfac_spec (:,:,:,1)= gfac_cnst

    dustref(:,:,:)= tau(:,:,:)

end if            !    --------------fixed dust -----------


return
end subroutine dust_distrib_bin

!   =========================================================================
!   =========================================================================

subroutine dust_distrib_fix( is, js, areolat, lon, lat, pref,                &
                                 pp, pph, delp,  aerosolmr,               &
                                 tau_spec,  sscat_spec, gfac_spec,        &
                                 dustref                )
!
!  fixed dust distrib for ames radiation only
!

integer, intent(in)                   ::  is, js
real, intent(in)                      ::  areolat
real, intent(in), dimension(:,:)      ::  lon    !  in radians
real, intent(in), dimension(:,:)      ::  lat    !  in radians
real, intent(in)                      ::  pref
real, intent(in), dimension(:,:,:)    ::  pp
real, intent(in), dimension(:,:,:)    ::  pph
real, intent(in), dimension(:,:,:)    ::  delp          !
real, intent(in), dimension(:,:,:,:)  ::  aerosolmr   !

real, intent(out), dimension(:,:,:,:) ::  tau_spec   !
real, intent(out), dimension(:,:,:,:) ::  sscat_spec !
real, intent(out), dimension(:,:,:,:) ::  gfac_spec  !

real, intent(out), dimension(:,:,:)   ::  dustref   !

!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt
real, dimension(size(lat,1),size(lat,2))         :: rnorm, optical_depth2d
real, dimension(size(lat,1),size(lat,2))         :: zmax, ystruc

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: dustrp, aerorp, cloudvis

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: tau, sscat, gfac

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: odfac, ssfac, ggfac

real, dimension(size(aerosolmr,4)  )    :: rdust_core
real, dimension(size(aerosolmr,4),1)    :: qex_sw, sscat_sw, gfac_sw

real, dimension(size(pp,1),size(pp,2),size(pp,3),size(aerosolmr,4)) :: taucc


real     ::   a, zmaxx, optical_depth
integer  ::  ie, je, id, jd, i, j, k, kd, nrad_trace, nt, ndx

integer  ::  iaa, iaam, nchan, lwchan, ird
real     ::  fac1, fact,  tau_thresh, qex_ref, qex_ice, reff

integer  ::   rad_tracers, nice

real     ::   rhoice=   0.917E3
real     ::   rhodust = 2.500E3

kd= size(pp,3)
jd= size(pp,2)
id= size(pp,1)

mcpu0 = (mpp_pe() == mpp_root_pe())

rdust_core(:)= 1.6E-6
qex_sw  (:,1)= qext_cnst
sscat_sw(:,1)= sscat_cnst
gfac_sw (:,1)= gfac_cnst

qex_ref= 2.5      !   for visible channel

tau  (:,:,:) = 0.0
sscat(:,:,:) = 0.0
gfac (:,:,:) = 0.0

tau_spec  (:,:,:,:)= 0.0
sscat_spec(:,:,:,:)= 0.0
gfac_spec (:,:,:,:)= 0.0

optical_depth= optical_depth_inpt

dustref(:,:,:)= 0.0

!              Add in a "pulse" of dust
if( optical_depth_pulse .gt. 0.0 ) then    !  ---------------- background dust
    call  dust_distrib_pulse( is, js, areolat, lon, lat, &
                optical_depth_pulse, pref, pp, pph, delp, tau )
    dustref(:,:,:)= tau(:,:,:)
endif

!            Formulate fixed dust distribution
if( do_inpt_dust_cycle_ktop ) then
    ndx= dust_indx(inpt_dust_cycle_dust_index)

    call dust_distrib_fixed2( is, js, areolat, lon, lat, optical_depth, &
                             pref, pp, pph, delp, tau, aerosolmr(:,:,:,ndx) )

else
    call dust_distrib_fixed2( is, js, areolat, lon, lat, optical_depth,  &
                                pref, pp, pph, delp, tau )
endif



!   only use the rescaling option for dust read in from input table. Note that the IR is rescaled
tau= scale_inpt_dust_column * tau


!              Add in a "pulse" of dust
tau(:,:,:)= dustref(:,:,:) + tau(:,:,:)


tau_spec  (:,:,:,1)= tau(:,:,:)
sscat_spec(:,:,:,1)= sscat_cnst
gfac_spec (:,:,:,1)= gfac_cnst

dustref(:,:,:)= tau(:,:,:)


!            Infrared Channels
do lwchan=1,nchan_lw+1
    nchan= nchan_sw + lwchan
    tau_spec  (:,:,:,nchan)= tau_spec(:,:,:,1)*qex_lw_dust(lwchan)/qex_ref
    sscat_spec(:,:,:,nchan)= sscat_lw_dust(lwchan)
    gfac_spec (:,:,:,nchan)= gfac_lw_dust(lwchan)
enddo


return
end subroutine dust_distrib_fix


!=======================================================================


subroutine cld_distrib_fix( is, js, areolat, lon, lat, pref, pp, pph, delp, &
                               aerosolmr, cldice )



! subroutine cld_distrib_fix( is, js, areolat, lon, lat, pref,                &
!                                 pp, pph, delp,  aerosolmr,               &
!                                 tau_spec,  sscat_spec, gfac_spec,        &
!                                 dustref                     )
!
!   call binned dust distrib for ames radiation only.
!   should eventually be merged with dust_distrib
!
integer, intent(in)                   ::  is, js
real, intent(in)                      ::  areolat
real, intent(in), dimension(:,:)      ::  lon    !  in radians
real, intent(in), dimension(:,:)      ::  lat    !  in radians
real, intent(in)                      ::  pref
real, intent(in), dimension(:,:,:)    ::  pp
real, intent(in), dimension(:,:,:)    ::  pph
real, intent(in), dimension(:,:,:)    ::  delp          !
real, intent(in), dimension(:,:,:,:)  ::  aerosolmr   !
real, intent(out), dimension(:,:,:)    ::  cldice   !

!!!  real, intent(out), dimension(:,:,:,:) ::  tau_spec   !
!!  real, intent(out), dimension(:,:,:,:) ::  sscat_spec !
!!  real, intent(out), dimension(:,:,:,:) ::  gfac_spec  !

!!  real, intent(out), dimension(:,:,:)   ::  dustref   !

!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt
real, dimension(size(lat,1),size(lat,2))         :: rnorm, optical_depth2d
real, dimension(size(lat,1),size(lat,2))         :: zmax, ystruc

!      dimension radwt(ntrace), sscatwt(ntrace), gfacwt(ntrace)

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: dustrp, aerorp, cloudvis

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: tau, sscat, gfac

real, dimension(size(pp,1),size(pp,2),size(pp,3)) :: odfac, ssfac, ggfac

real, dimension(size(aerosolmr,4)  )    :: rdust_core
real, dimension(size(aerosolmr,4),1)    :: qex_sw, sscat_sw, gfac_sw

real, dimension(size(pp,1),size(pp,2),size(pp,3),size(aerosolmr,4)) :: taucc


real     ::   a, zmaxx, optical_depth
integer  ::  ie, je, id, jd, i, j, k, kd, nrad_trace, nt, ndx

integer  ::  iaa, iaam, nchan, lwchan, ird
real     ::  fac1, fact,  tau_thresh, qex_ref, qex_ice, reff

integer  ::   rad_tracers, nice

real     ::   rhoice=   0.917E3
real     ::   rhodust = 2.500E3

kd= size(pp,3)
jd= size(pp,2)
id= size(pp,1)

mcpu0 = (mpp_pe() == mpp_root_pe())


nice= ice_bin_indx(1)

cloudvis(:,:,:)= aerosolmr(:,:,:,nice)

#ifdef SKIP2

       reff=  4.0 * 1.e-6    ! effective cloud radius
       fact= scale_rad_ice * 3.0 / ( 4.0*grav*rhoice*reff )
       nchan= 1
       qex_ice= 2.5            ! effective sw cloud extinction
       cloudvis(:,:,:)= aerosolmr(:,:,:,nice)*fact*qex_ice*delp(:,:,:)

#endif SKIP2



       if( rad_ice_scheme == 2 ) then
          ystruc(:,:)= 0.5*( 1.0-tanh( abs(RAD_TO_DEG*lat(:,:)-70.0)/20.0) )
          WHERE( lat > 0.0 )
                ystruc(:,:)= MAX( ystruc(:,:), scale_rad_ice_np )
          ELSEWHERE
                ystruc(:,:)= MAX( ystruc(:,:), scale_rad_ice_sp )
          END WHERE
          DO k= 1, kd
             cloudvis(:,:,k)= cloudvis(:,:,k)*ystruc(:,:)
          ENDDO
       endif

       if( rad_ice_scheme == 3 ) then
          ystruc(:,:)= 0.5*( 1.0-tanh( ( abs(RAD_TO_DEG*lat(:,:))-50.0)/7.0) )
          WHERE( lat > 0.0 )
                ystruc(:,:)= MAX( ystruc(:,:), scale_rad_ice_np )
          ELSEWHERE
                ystruc(:,:)= MAX( ystruc(:,:), scale_rad_ice_sp )
          END WHERE

          DO k= 1, kd
             cloudvis(:,:,k)= cloudvis(:,:,k)*ystruc(:,:)
          ENDDO
       endif

       if( rad_ice_scheme == 4 ) then
          fact= 1.0/(1.0-tanh( -200.0/150.0 ) )
          DO k= 1, kd
             cloudvis(:,:,k)= fact*cloudvis(:,:,k)*(1.0-tanh((pp(:,:,k)-200.0)/200.0) )
          ENDDO
       endif


       if( rad_ice_scheme == 4 ) then
          fact= 1.0/(1.0-tanh( -200.0/150.0 ) )
          ystruc(:,:)= 0.5*( 1.0-tanh( ( abs(RAD_TO_DEG*lat(:,:))-50.0)/7.0) )
          DO k= 1, kd
             cloudvis(:,:,k)= ystruc(:,:)*cloudvis(:,:,k)*(1.0-tanh((pp(:,:,k)-200.0)/150.0) )
          ENDDO
       endif


       if( rad_ice_scheme == 5 ) then
          ystruc(:,:)= 0.5*( 1.0-tanh( ( abs(RAD_TO_DEG*lat(:,:))-45.0)/7.0) )
          fact= 1.0/(1.0-tanh( -300.0/100.0 ) )
          DO k= 1, kd
             cloudvis(:,:,k)= fact*cloudvis(:,:,k)*ystruc(:,:)*  &
                       (1.0-tanh((pp(:,:,k)-300.0)/100.0) )
          ENDDO
       endif

       if( rad_ice_scheme == 6 ) then
          ystruc(:,:)= 0.5*( 1.0-tanh( ( abs(RAD_TO_DEG*lat(:,:))-40.0)/10.0) )
          WHERE( lat > 0.0 )
                ystruc(:,:)= MAX( ystruc(:,:), scale_rad_ice_np )
          ELSEWHERE
                ystruc(:,:)= MAX( ystruc(:,:), scale_rad_ice_sp )
          END WHERE

          DO k= 1, kd
             cloudvis(:,:,k)= cloudvis(:,:,k)*ystruc(:,:)
          ENDDO
       endif


       if( rad_ice_scheme == 7 ) then
          ystruc(:,:)= 0.5*( 1.0-tanh( ( abs(RAD_TO_DEG*lat(:,:))-45.0)/7.0) )
          fact= 1.0/(1.0-tanh( -500.0/100.0 ) )
          DO k= 1, kd
             cloudvis(:,:,k)= fact*cloudvis(:,:,k)*ystruc(:,:)*  &
                       (1.0-tanh((pp(:,:,k)-500.0)/100.0) )
          ENDDO
       endif

  cldice= cloudvis



#ifdef SKIP2
       WHERE( cloudvis(:,:,:) .gt. tau_spec(:,:,:,1) )
             tau_spec  (:,:,:,1) = cloudvis(:,:,:)
             sscat_spec(:,:,:,1)= sscat_sw_ice
             gfac_spec (:,:,:,1)= gfac_sw_ice
       END WHERE

       cldref(:,:,:)= cloudvis(:,:,:)

       DO lwchan=1,nchan_lw+1
           nchan= nchan_sw + lwchan
           aerorp(:,:,:)= cloudvis(:,:,:)*qex_lw_ice(lwchan)/qex_ref
           WHERE( aerorp(:,:,:) .gt.  tau_spec(:,:,:,nchan)  )
                tau_spec  (:,:,:,nchan) = aerorp(:,:,:)
                sscat_spec(:,:,:,nchan)= sscat_lw_ice(lwchan)
                gfac_spec (:,:,:,nchan)= gfac_lw_ice (lwchan)
           END WHERE
       ENDDO
#endif SKIP2


return
end subroutine cld_distrib_fix



!=======================================================================
!=======================================================================

 subroutine dust_distrib_fixed2( is, js, areolat, lon, lat, optical_depth, &
                                      pref, pp, pph, delp, tau, aerosol )
!
!   calculate fixed dust opacity distribution
!

integer, intent(in)                 ::  is, js
real, intent(in)                    ::  areolat
real, intent(in), dimension(:,:)    ::  lon    !  in radians
real, intent(in), dimension(:,:)    ::  lat    !  in radians
real, intent(in)                    ::  optical_depth
real, intent(in)                    ::  pref
real, intent(in), dimension(:,:,:)  ::  pp   ! kd
real, intent(in), dimension(:,:,:)  ::  pph  ! kp
real, intent(in), dimension(:,:,:)  ::  delp   ! kd
real, intent(out), dimension(:,:,:) ::  tau  ! opacity

real, optional, intent(in), dimension(:,:,:)  ::  aerosol


! evaluate optical depth between"flux" pressure levels
! Assume pressure increases with index number.
!
!    The parameter vscale (so-called Conrath paramter)
!         serves to pick at a height above which
!         the dust mixing ratio decreases stongly.

!    For vscale=0.01, this is at about 35 km.

!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt
real, dimension(size(pp,1),size(pp,2),size(pp,3))     :: tau2
real, dimension(size(lat,1),size(lat,2))              :: rnorm, optical_depth2d
real, dimension(size(lat,1),size(lat,2))              :: zmax, colav, optx, ystruc

real     ::   a, zmaxx
integer  ::  ie, je, id, jd, i, j, k, kd, kzz

kd= size(pp,3);   jd= size(pp,2);   id= size(pp,1)

ie= is + id - 1;   je= js + jd - 1

a= conrath
zmaxx= 70.0
zmax(:,:)= 35.0


if( .not. present(aerosol)  .and. dust_cycle_scheme > 0 ) &
    call error_mesg ('dust_distrib_fixed2',  'no aerosol array is present  ', FATAL)



!  Establish spatial structure

if( do_inpt_dust_cycle ) then
    call season_table( areolat, is, ie, js, je, optical_depth2d, zmax  )

    if ( add_storm_opacity ) then
        call add_storm( areolat, lon, lat, optx );
        optical_depth2d= optical_depth2d + optx
    endif

    if ( optical_depth > 0.0 ) then
        optical_depth2d= optical_depth2d + optical_depth
    endif

elseif (do_lat_vary_dust )  then
    optical_depth2d(:,:)= optical_depth * ( 0.3 + 0.7*cos( lat(:,:) ) )
    zmax(:,:)= zmaxx - 40.0*sin( lat(:,:) )**2

else
    optical_depth2d(:,:)= optical_depth
    zmax(:,:)= zmaxx

endif         !           -------------- dust cycle look-up ----------


!  Establish vertical structure

if( dust_cycle_scheme < 1 ) then  !------ default dust_cycle_scheme = 0
    do k= 1, kd
        if(conrath_type .eq. 0) then
            tau(:,:,k)=exp( a*(1.0 -  (pref/pp(:,:,k))   ) )
        else
! lmd-type vertical dust distribution
            tau(:,:,k)=exp( 0.007*(1.0 -  (700.0/pp(:,:,k))**(zmaxx/zmax(:,:))   ) )
        end if
    enddo

endif


!   If present(aerosol), then use aerosol to establish height of dust
!          distribution

if( present(aerosol) ) then

    if( dust_cycle_scheme == 1 ) then  !------ dust_cycle_scheme = 1
!           Dust depth using zmax based on passive tracer field

        colav(:,:)= 0.0
        do k= 1, kd
            colav(:,:)= colav(:,:) + aerosol(:,:,k)
        enddo
        colav= colav / (1.0*kd)

        !           formulate dust cloud top (in km)   (10 km atmos scale height)
        do i= 1, id
            do j= 1, jd
                do k= 1, kd
                    if( aerosol(i,j,k) > 0.5*colav(i,j) )  then
                        kzz=k
                        exit
                    endif
                enddo
                zmax(i,j)= 10.0 * log( pref / pph(i,j,kzz) )
            enddo
        enddo

        zmax= max( 10.0, zmax )
        do k= 1, kd
            tau(:,:,k)=exp( a*(1.0 -  (pref/pp(:,:,k))**(zmaxx/zmax(:,:))   ) )
        enddo
    endif


    if( dust_cycle_scheme > 1 ) then  !------ dust_cycle_scheme = 2
!           Dust profile based on passive tracer profile
        tau(:,:,1:kd)= aerosol(:,:,1:kd)

        do k= kd-4, kd
            tau(:,:,k)= aerosol(:,:,kd-4)
        enddo

!      eliminate opacity due to dust extending to the top layers
        do k= 1, 3
            tau(:,:,k)= 0.0
        enddo

    endif

endif

!               ---- Normalize ------

!         Evaluate accumulated opacity at half-levels
opt(:,:,1)= 0.0
do k= 1, kd
    opt(:,:,k+1)= opt(:,:,k) + tau(:,:,k)*delp(:,:,k)
enddo

!               now, scale the column optical depth


rnorm = 0.
where (opt(:,:,kd+1) .gt. 0.) rnorm(:,:)= optical_depth2d(:,:) * pph(:,:,kd+1) / ( pref *opt(:,:,kd+1) )
do k= 1, kd+1
    opt(:,:,k)= opt(:,:,k)*rnorm(:,:)
enddo

!       Optionally add in a "pulse" of dust:  In this case, it is a fixed fraction of the total column

if( optical_depth_pulse_var .gt. 0.0 ) then    !  ---------------- background dust

    rnorm = 0.0
    where (opt(:,:,kd+1) .gt. 0.) rnorm(:,:)= pph(:,:,kd+1) / ( pref *opt(:,:,kd+1) )

    do k= 1, kd
        tau(:,:,k)= tau(:,:,k) * rnorm(:,:)
    enddo

    !       now create a normalized "pulse" vertical column

    call  dust_distrib_pulse2( is, js, areolat, lon, lat,            &
                1.0, pref, pp, pph, delp, tau2 )
    !

    !         assume simple latitude structure for the pulse
    ystruc(:,:)=  cos( lat(:,:) ) **2

    !           add "pulse" distribution
    do k= 1, kd
        tau(:,:,k)= tau(:,:,k) +  optical_depth_pulse_var * ystruc(:,:)*tau2(:,:,k)
    enddo


    !    normalize again
    opt(:,:,1)= 0.0
    do k= 1, kd
        opt(:,:,k+1)= opt(:,:,k) + tau(:,:,k)*delp(:,:,k)
    enddo

    rnorm= optical_depth2d(:,:) / opt(:,:,kd+1)

    do k= 1, kd+1
        opt(:,:,k)= opt(:,:,k)*rnorm(:,:)*pph(:,:,kd+1)/pref
    enddo


else

    rnorm = 0.
    where (opt(:,:,kd+1) .gt. 0.) rnorm(:,:)= optical_depth2d(:,:) * pph(:,:,kd+1) / ( pref *opt(:,:,kd+1) )
    do k= 1, kd+1
        opt(:,:,k)= opt(:,:,k)*rnorm(:,:)
    enddo


endif        !     -----------------end of pulse calculation --------


!              optical depth of individual layers
do k= 1, kd
    tau(:,:,k)= opt(:,:,k+1)-opt(:,:,k)
enddo


return
end subroutine dust_distrib_fixed2

!=======================================================================
!=======================================================================


subroutine dust_distrib_pulse2( is, js, areolat, lon, lat, optical_depth, &
                                     pref, pp, pph, delp, tau )
!
!  helper routine for pulse
!

integer, intent(in)                 ::  is, js
real, intent(in)                    ::  areolat
real, intent(in), dimension(:,:)    ::  lon    !  in radians
real, intent(in), dimension(:,:)    ::  lat    !  in radians
real, intent(in)                    ::  optical_depth
real, intent(in)                    ::  pref
real, intent(in), dimension(:,:,:)  ::  pp   ! kd
real, intent(in), dimension(:,:,:)  ::  pph  ! kp
real, intent(in), dimension(:,:,:)  ::  delp   ! kd
real, intent(out), dimension(:,:,:) ::  tau  ! opacity


! evaluate optical depth between"flux" pressure levels
! Assume pressure increases with index number.
!
!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt
real, dimension(size(lat,1),size(lat,2))              :: rnorm, optical_depth2d, sig

real   ::   a,  pwid, pulse_ctr, snot
integer  ::  ie, je, id, jd, i, j, k, kd


kd= size(pp,3)
jd= size(pp,2)
id= size(pp,1)

!  ie= is + id - 1;   je= js + jd - 1

!    optical_depth2d(:,:)= optical_depth * cos( lat(:,:) ) **2

optical_depth2d(:,:)= optical_depth


snot= log( pulse_center/pref )

pwid= pulse_width

mcpu0 = (mpp_pe() == mpp_root_pe())

do k= 1, kd
    sig(:,:)= log(  pp(:,:,k)/pref  )
    tau(:,:,k)= exp( -( (sig(:,:) - snot)/pwid )**2   )
enddo

!  if( mcpu0 ) print *, 'have calculated initial tau in dust-pulse '

opt(:,:,1)= 0.0
do k= 1, kd
    opt(:,:,k+1)= opt(:,:,k) + tau(:,:,k)*delp(:,:,k)
enddo

!    since this is an elevated dust layer, dont normalize with surface pressure
rnorm(:,:) = 0.
where (opt(:,:,kd+1) .gt. 0.)  rnorm(:,:)= 1.0  / ( opt(:,:,kd+1) )
do k= 1, kd+1
    opt(:,:,k)= opt(:,:,k)*rnorm(:,:)
enddo

!              optical depth of individual layers
do k= 1, kd
    tau(:,:,k)= opt(:,:,k+1)-opt(:,:,k)
enddo


return
end subroutine dust_distrib_pulse2

!=======================================================================
!=======================================================================


subroutine dust_opacity_assim( is, js, areolat, lon, lat,  &
                                  aerosol,  pref, pp, pph, delp, tau, taucc )
!
!   calculate and scale total dust opacity
!

integer, intent(in)                   ::  is, js
real, intent(in)                      ::  areolat
real, intent(in), dimension(:,:)      ::  lon    !  in radians
real, intent(in), dimension(:,:)      ::  lat    !  in radians
real, intent(in), dimension(:,:,:,:)  ::  aerosol
real, intent(in)                      ::  pref
real, intent(in), dimension(:,:,:)    ::  pp   !
real, intent(in), dimension(:,:,:)    ::  pph  !
real, intent(in), dimension(:,:,:)    ::  delp   ! k
real, intent(out), dimension(:,:,:)   ::  tau    ! opacity
real, intent(out), dimension(:,:,:,:) ::  taucc

!   LOCAL:

real, dimension(size(pph,1),size(pph,2),size(pph,3))  :: opt

real, dimension(size(aerosol,1),size(aerosol,2),size(aerosol,3),size(aerosol,4)) :: aerorp

real   ::   a, zmaxx, fac1, delareo, scale, scale2
integer  ::  ie, je, id, jd, i, j, k, kd, kzz, nt, ntrace

real :: reffec = 1.0E-6       !  Effective particle radius
real :: rhodust = 2.50E3     !   kg/m^3

logical :: mcpu0

kd= size(pp,3);   jd= size(pp,2);   id= size(pp,1); ntrace= size(aerosol,4)

!            Form predicted column-integrated dust
!        Assume basic scaling;

scale = 0.75*qext_cnst/(grav*rhodust*reffec )

mcpu0 = (mpp_pe() == mpp_root_pe())

aerorp(:,:,:,:)= max( aerosol(:,:,:,:), 0.0 )

!         sum over ntrace aerosol contributors
tau(:,:,:)= 0.0

do nt= 1, ntrace
    do k= 1, kd
        taucc(:,:,k,nt)= aerorp(:,:,k,nt)*delp(:,:,k)*scale
    enddo
    tau(:,:,:)= tau(:,:,:) + taucc(:,:,:,nt)
enddo


return
end subroutine dust_opacity_assim

!   =================================================================
!   =================================================================

subroutine season_table( areolat, is, ie, js, je, optical_depth2d, zmax  )
!
!   fetch dust distribution from dust map
!

real, intent(in)                    ::   areolat
integer, intent(in)                 ::   is, ie, js, je
real, intent(out), dimension(:,:)   ::  optical_depth2d, zmax    !  in radians

real   ::   a, zmaxx, fac1, delareo
integer  ::  iaa, iaam


do iaa= 1, size(areo_cycle)
    if( areo_cycle(iaa) > areolat ) then
        exit
    else
        cycle
    endif
enddo

iaam= iaa - 1
if( iaam < 1 )  iaam= iaam + size(areo_cycle)

delareo= areo_cycle(iaa) - areo_cycle(iaam)
if( delareo  < 0.0 )       delareo = delareo + 360.0

fac1= (areo_cycle(iaa)- areolat)/delareo

zmax(:,:)= fac1*zmax_cycle(is:ie,js:je,iaam) + (1.0-fac1)*zmax_cycle(is:ie,js:je,iaa)

optical_depth2d(:,:)= fac1*taufill_cycle(is:ie,js:je,iaam) +  &
                         (1.0-fac1)*taufill_cycle(is:ie,js:je,iaa)

return
end subroutine season_table

!   =================================================================
!   =================================================================



 subroutine aerosol_end ( days )
!
!  release memory at end of simulation
!
integer,  intent(in):: days

deallocate ( dust_indx )
deallocate ( ice_bin_indx ) !TB18c
deallocate ( radiative_active_inpt )
deallocate ( reff_dust )

end  subroutine aerosol_end




end module aerosol_mod










