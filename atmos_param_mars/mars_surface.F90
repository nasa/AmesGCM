
module mars_surface_mod
!
!  module for calculating soil temperatures and surface CO2 condensation

use     constants_mod, only: PI, KAPPA, CP_AIR, GRAV, STEFAN, co2_lheat, SECONDS_PER_DAY

use           fms_mod, only: error_mesg, FATAL,       &
                             open_namelist_file, check_nml_error, &
                             mpp_pe, mpp_root_pe, close_file,     &
                             write_version_number, stdlog,        &
                             uppercase, read_data, write_data,    &
                             field_size, field_exist

use fms2_io_mod,            only:  file_exists

use   mpp_domains_mod, only: domain2d

use        fms_io_mod, only: register_restart_field, restart_file_type, &
                             save_restart, restore_state, get_mosaic_tile_file

use  aerosol_util_mod, only: do_moment_water, do_bulk_water

use  time_manager_mod, only: time_type, get_time

use  diag_manager_mod, only: diag_axis_init, register_diag_field,  &
                             register_static_field, send_data
use initracer_mod

implicit none


!-----------------------------------------------------------------------
!---------- interfaces ------------

public :: mars_surface_init, mars_surface_end,  progts, sfc_snow, sfc_frost, &
         sfc_frost_mom, sfc_frost_blk, sfc_h2o2_chem, id_sfc_h2o2_chem, &
         cumulative_prec_blk

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

integer :: nlayers = 16                 !  number of soil layers
integer :: nlay_rst = 0                 !  number of soil layers in restart
integer :: d2is = 7   ! layer depth to ice in south
real :: gk1 = 1000.   ! soil ice ti in north
real :: gk2 = 1000.   ! soil ice ti in south

real ::  zoland=      0.01              !  surface roughness [cm]
real ::  drag_cnst=   0.003             !  surface drag constant
real ::  soil_ti=     350.0             !  Nominal background value; mks units
real ::  soil_alb=    0.25              !  Nominal background value; mks units
real ::  soil_temp=   170.0             !  Isothermal soil temperature if no file is used

real :: albedo_ice_np = 0.65             !  Northern CO2 ice albedo
real :: albedo_ice_sp = 0.43            !  Southern CO2 ice albedo
real :: emiss_ice_np = 0.81             !  Northern CO2 ice emissivity
real :: emiss_ice_sp = 0.88             !  Southern CO2 ice emissivity

real :: alpha = 0.5                   !  parameter for surface temp calculation (0.5 = semi implicit, 0.0 = fully explicit, 1.0 = fully implicit)
real :: albedo_h2o = 0.4                ! h2o ice albedo if do_waterice_albedo
real :: albedo_h2o_liquid = 0.07        ! liquid water albedo
real :: emiss_h2o = 1.0                 ! h2o ice emissivity if do_waterice_albedo 

logical :: do_co2_condensation= .true.  !  Maintain surface temperature at or above Tcrit

logical :: do_subsfc_ice=   .true.      !  Include influence of subsurface water ice on thermal diffusivity
real :: init_sfc_frost= 500.         !  Initial surface frost value (units?)
logical :: edit_subsfc_temp= .false.    !  change subsurface temperatures after initialization
integer :: subsfc_ice_case = 7          !  subsurface ice distribution scenario
logical :: use_legacy_soil = .false.    !  use legacy soil rho cp values
logical :: use_equilibrated_ts = .true. !  use equilibrated surface temperatures

real ::  rescale_sh_ti = -100.0         !  used to rescale TI in the SH when positive
real ::  np_cap_ti = -900.0             !  when positive, set ti when lat > np_cap_lat
real ::  np_cap_lat = 65.0              !  latitude boundary for np_cap_ti
real ::  np_cap_ti_max = 800.           !  Maximum np_cap_ti
real ::  frost_threshold = 1.e6         !  The threshhold above which water frost influences surface albedo
logical :: do_waterice_albedo= .false.  !  Activate change of albedo by water ice

real  :: rho_ground= 1.50*1.0E3         !  Soil density:    kg / m**3
real  :: cp_ground=  627.9              !  Soil heat capcity:     joules / kg / K

real  :: rho_ice = 916.7                !  water ice density: kg/m**3
real  :: cp_ice  = 2108.0               !  water ice heat capacity:  joules/kg/K
real  :: rho_legacy_soilice = 1781.99   !  legacy soil-ice mixed density
real  :: cp_legacy_soilice  = 1404.09   !  legacy soil-ice mixed heat capacity
real  :: rho_legacy_soil = 1481.39      !  legacy soil density
real  :: cp_legacy_soil  =  840.00      !  legacy soil heat capacity
real  :: rho_co2    = 910.0             !  CO2 ice density [kg/m^3]
real  :: cp_co2     = 1075.0            !  CO2 ice specific heat

real, dimension(:), allocatable    ::   delzg, &                        !  soil layer thicknesses [m]
                                        zgrid, zgrid_r                  !  soil boundary depths [m]

real, dimension(:,:,:),   allocatable, save  ::  tsoil,tsoil_r          !  soil temperatures
real, dimension(:,:),     allocatable, save  ::  thermal_inertia        !  map of surface thermal inertia
real, dimension(:,:),     allocatable, save  ::  t_surf                 !  surface temperature
real, dimension(:,:),     allocatable, save  ::  sfc_albedo             !  map of surface albedo
real, dimension(:,:),     allocatable, save  ::  sfc_emiss              !  map of surface emissivity
real, dimension(:,:),     allocatable, save  ::  sfc_snow               !  (CO2)
real, dimension(:,:,:),     allocatable, save  ::  sfc_frost            !  (water)
real, dimension(:,:,:),     allocatable, save  ::  sfc_frost_mom      !  (water)
real, dimension(:,:,:),     allocatable, save  ::  sfc_frost_blk      !  (water)
real, dimension(:,:,:),     allocatable, save  ::  cumulative_prec_blk  !  (cumulative precipitation bulk scheme)
real, dimension(:,:),     allocatable, save  ::  sfc_h2o2_chem          !  (h2o2)
real, dimension(:,:),     allocatable, save  ::  sfc_roughness          !  map of surface roughness
real, dimension(:,:),     allocatable, save  ::  sfc_topo               !  map of surface topography
real, dimension(:,:),     allocatable, save  ::  grs_ice                !  map of GRS ice latitudinal extent
real, dimension(:,:),     allocatable, save  ::  npcflag                !  north polar residual cap flag

real, dimension(:,:,:),   allocatable, save  ::  soil_icex, soil_icex_r

#ifdef REGOLITH
real, dimension(:,:,:),   allocatable, save  ::  soil_water
real, dimension(:,:,:),   allocatable, save  ::  soil_absorb
#endif

!--- for restart file
type(restart_file_type), pointer, save :: Surf_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()   !needed for tile restarts
type(restart_file_type), pointer, save :: Surf_restart2 => NULL()
type(restart_file_type), pointer, save :: Til_restart2 => NULL()   !needed for tile restarts
logical                                :: in_different_file = .false.
logical                                :: rst2 = .false. ! input restart has different nlayers

character(len=128) :: version='$Id: mars_surface.F90,v 1.1.2.1.2.1 2011/11/22 20:56:59 rjw Exp $'
character(len=128) :: tagname='$Name: mars_feb2012_rjw $'


integer :: id_insol
real    :: missing_value = -1.e10
character(len=12) :: mod_name = 'mars_surface'

logical :: module_is_initialized = .false.

integer, dimension(:),   allocatable, save  ::  id_frost_mom
integer  :: id_sfc_h2o2_chem
integer  ::  id_zgrid, id_subsfc, id_ts, id_thin, id_frost, id_frost_blk,id_snow, id_sflux
integer  ::  id_alb, id_cprecip_blk

logical ::   mcpu0

namelist /surface_data_nml/  zoland, drag_cnst, soil_ti, soil_temp, &
                            do_co2_condensation,                    &
                            do_subsfc_ice, subsfc_ice_case,         &
                            albedo_ice_np, albedo_ice_sp,           &
                            emiss_ice_np, emiss_ice_sp,             &
                            np_cap_ti, np_cap_lat, frost_threshold, &
                            np_cap_ti_max, soil_alb,  &
                            edit_subsfc_temp, rescale_sh_ti,        &
                            nlayers,do_waterice_albedo,use_legacy_soil, &
                            alpha, albedo_h2o, emiss_h2o,nlay_rst,d2is, &
                            gk1,gk2,use_equilibrated_ts,            &
                            init_sfc_frost,albedo_h2o_liquid

contains


subroutine mars_surface_init ( nlon, mlat, lonb, latb, lon, lat, axes, Time, phys_domain )
!
!  Initialize the surface model by reading various maps and creating the subsurface grid
!  
use tracer_manager_mod, only: query_method, get_tracer_index,  &
                                 get_number_tracers, get_tracer_names

integer, intent(in)                   :: nlon, mlat
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon , lat
integer, intent(in) :: axes(4)

type(time_type), intent(in) :: Time
    
type(domain2d),      intent(inout) :: phys_domain

!           Default soil layer thickness:  meters
real, dimension(13)  :: delzg_12 =                     &
                                        (/ 0.0030, 0.0050, 0.0075,  &
                                           0.0110, 0.0150, 0.0200,  &
                                           0.0400, 0.0800, 0.1600,  &
                                           0.3200, 0.640,  1.500, 3.50 /)

character (len=128) :: filename, fieldname, tracer_name, tname, tracer_name2, f_rst

integer  ::  unit, io, ierr, id, jd, i, j, k, is, js, ie, je
integer  ::  im, jm, km, days, fld_dims(4), nt, ndx
integer  ::  id_alb0,  id_thin0, id_emiss0, id_ruff0, id_gice0, id_npc0
logical  ::  used

mcpu0=  (mpp_pe() == mpp_root_pe())

!     ----- read namelist -----

if (file_exists('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1; do while (ierr /= 0)
        read  (unit, nml=surface_data_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'surface_data_nml')
    enddo
10     call close_file (unit)
endif

!     ----- write version info and namelist to log file -----

call write_version_number (version,tagname)
if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=surface_data_nml)


id= size(lonb,1)-1; jd= size(latb,2)-1


if(mcpu0) print *, 'Allocating storage for soil surface properties:  ', id, jd

allocate (  thermal_inertia(id,jd)  )
allocate (  sfc_albedo     (id,jd)  )
allocate (  sfc_emiss      (id,jd)  )
allocate (  t_surf         (id,jd)  )
allocate (  sfc_snow       (id,jd)  )
allocate (  sfc_frost      (id,jd,1)  )
allocate (  sfc_frost_mom      (id,jd,nice_mass)  )
allocate (  sfc_frost_blk  (id,jd,1)  )
allocate (  cumulative_prec_blk  (id,jd,3)  )
allocate (  sfc_h2o2_chem  (id,jd)  )
allocate (  sfc_roughness  (id,jd)  )
allocate (  sfc_topo       (id,jd)  )
allocate (  grs_ice        (id,jd)  )
allocate (  npcflag        (id,jd)  )

allocate ( delzg(nlayers+1) )
allocate ( zgrid(nlayers+1) )

allocate (  tsoil      (id,jd,nlayers+1)  )
allocate (  soil_icex  (id,jd,nlayers  )  )

#ifdef REGOLITH
allocate (  soil_water      (id,jd,nlayers)  )
allocate (  soil_absorb     (id,jd,nlayers)  )
#endif

!         ------------Read soil restart file --------------
!         Use this to obtain nlayers and subsurface temperatures

filename= 'INPUT/soil_temp.res.nc'
f_rst= 'INPUT/soil_temp.res.nc'
!new way of reading restarts
rst2 = (nlay_rst /=0 .and. nlay_rst /= nlayers)
if (rst2) then
    allocate(tsoil_r(id,jd,nlay_rst+1))
    allocate(soil_icex_r(id,jd,nlay_rst))
    allocate(zgrid_r(nlay_rst+1))
endif

call mars_surface_register_restart('soil_temp.res.nc',phys_domain)

if( file_exists( trim( filename ) ) ) then
    
    if (rst2) then
        call restore_state(Surf_restart2)
        if (in_different_file) call restore_state(Til_restart2)  
        km = min(nlayers,nlay_rst)
        tsoil(:,:,:km) = tsoil_r(:,:,:km)
        soil_icex(:,:,:km) = soil_icex_r(:,:,:km)
        do k= 1, km
            delzg(k)= zgrid_r(k+1)-zgrid_r(k)
        enddo

        !!! check nb layers
        if( nlayers > km ) then
            if(mcpu0) print *, 'initializing addition soil layers  '
            do k= km, nlayers
                delzg(k)= 1.8*delzg(k-1)
                tsoil    (:,:,k)= soil_temp
                soil_icex(:,:,k)= 0.0
            enddo
        endif
        tsoil(:,:,nlayers+1)= soil_temp
        delzg(nlayers+1)= 1.8 * delzg(nlayers)

        zgrid(1)= 0.0
        do k= 2, nlayers+1
            zgrid(k)= zgrid(k-1) + delzg(k-1)
        enddo
    else
        call restore_state(Surf_restart)
        if (in_different_file) call restore_state(Til_restart)  
        DO k= 1, nlayers
            delzg(k)= zgrid(k+1)-zgrid(k)
        ENDDO
        delzg(nlayers+1)= 1.8 * delzg(nlayers)
    endif
    
      
    if(mcpu0) print *, 'Have read soil restart file: nlayers= ', nlayers

    t_surf(:,:)= tsoil(:,:,1)

    do nt= 1, nice_mass
        if(mcpu0)  print*, 'surface frost mom',  nt, sfc_frost_mom(1,:,nt)
    enddo

else

    if( nlayers <= 12 ) then
        do k= 1, nlayers+1
            delzg(k)= delzg_12(k)
        enddo
    else if (nlayers <= 20) then
        delzg(1:12)= delzg_12(1:12)
        do k= 13, nlayers+1
            delzg(k)= 2.0*delzg(k-1)
        enddo
    else
        delzg(1)= 0.015
        do k = 2, nlayers+1
            delzg(k)= 1.2*delzg(k-1)
        enddo
    endif

    zgrid(1)= 0.0
    do k= 2, nlayers+1
        zgrid(k)= zgrid(k-1) + delzg(k-1)
    enddo

    if (use_equilibrated_ts) then
        where( lat < -15.0*pi/180. )
            t_surf(:,:) = 165.+55.5*0.5*(1.+cos((lat+0.261799)*180./75.))
        elsewhere
            t_surf(:,:) = 167.+52.*cos((lat+0.261799)*90./105.)
        end where
        do k=1,nlayers+1
            tsoil(:,:,k) = t_surf(:,:)
        enddo
    else
        tsoil(:,:,:)= soil_temp
        t_surf(:,:)= tsoil(:,:,1)
    endif

    sfc_snow(:,:)= 0.0
    cumulative_prec_blk(:,:,:)= 0.0
    do nt=1, nice_mass
        where( lat > 80.0*pi/180.0 )
            sfc_frost(:,:,1)= init_sfc_frost
            sfc_frost_mom(:,:,nt)= init_sfc_frost
            sfc_frost_blk(:,:,1)= init_sfc_frost
        elsewhere
            sfc_frost(:,:,1)= 0.0
            sfc_frost_mom(:,:,nt)= 0.0
            sfc_frost_blk(:,:,1)= 0.0
        end where
    enddo
    soil_icex(:,:,:)= 0.0
    sfc_h2o2_chem(:,:) = 0.0

endif

if( rst2 ) then
    deallocate(tsoil_r)
    deallocate(soil_icex_r)
    deallocate(zgrid_r)
endif

if( edit_subsfc_temp  ) then

    where ( abs(lat) >  70.0 * pi/180.0   )
        tsoil(:,:,nlayers+1)=  150.0
    end where

endif



!     ----- register diagnostic fields -----

id_zgrid = diag_axis_init('zgrid', zgrid, 'm', 'z', 'soil level depths',  &
                    direction=-1,set_name=mod_name)

id_subsfc = register_diag_field ( mod_name, 'tsoil',            &
                                 (/axes(1:2),id_zgrid/), Time, &
                                'Soil Temperature', 'K',       &
                                 missing_value=missing_value )

id_sflux = register_diag_field ( mod_name, 'sflux',            &
                                 (/axes(1:2),id_zgrid/), Time, &
                                'Soil Heat Diffusion Flux', 'W/m/m',       &
                                 missing_value=missing_value )

id_ts = register_diag_field ( mod_name, 'ts',                      &
                                 (/axes(1:2)/), Time,             &
                                'Surface Temperature', 'K',       &
                                 missing_value=missing_value )

id_thin = register_diag_field ( mod_name, 'thin',                  &
                                 (/axes(1:2)/), Time,             &
                                'Surface Thermal Inertia', 'mks', &
                                 missing_value=missing_value )

id_frost = register_diag_field ( mod_name, 'frost',                &
                                 (/axes(1:2)/), Time,             &
                                'Surface water ice for bulk microphysics', 'kg/m/m',    &
                                 missing_value=missing_value )

id_frost_blk = register_diag_field ( mod_name, 'frost_blk',                &
                                 (/axes(1:2)/), Time,             &
                                'Surface water ice for new bulk microphysics', 'kg/m/m',    &
                                 missing_value=missing_value )

id_cprecip_blk = register_diag_field ( mod_name, 'cprecip_blk',                &
                                 (/axes(1:2)/), Time,             &
                                'Cumulative precipitation bulk microphysics', 'kg/m/m',    &
                                 missing_value=missing_value )

id_sfc_h2o2_chem = register_diag_field ( mod_name, 'sfc_h2o2_chem', &
                                 (/axes(1:2)/), Time,             &
                                 'Surface H2O2', 'kg/m/m',    &
                                 missing_value=missing_value )

allocate ( id_frost_mom(nice_mass) )
ndx= ice_mass_indx(1)
call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)
tname= "frost_mom"
id_frost_mom(1) = register_diag_field ( mod_name, trim(tname),  &
         axes(1:2), Time, 'Surface ice mass for moment microphysics', 'kg/m/m', &
         missing_value=missing_value)

if (nice_mass.gt.1) then
    do nt= 2, nice_mass
        ndx= ice_mass_indx(nt)
        call get_tracer_names(model_atmos, ndx, tracer_name2)
        print*, 'len 2',len_trim(adjustl(tracer_name2))
        print*, 'len 1',len_trim(adjustl(tracer_name))

        tname= "frost_mom" // trim(tracer_name2(len_trim(adjustl(tracer_name))+1:len_trim(adjustl(tracer_name2))))
        id_frost_mom(nt) = register_diag_field ( mod_name, trim(tname),  &
             axes(1:2), time, 'surface ice mass for moment microphysics', 'kg/m/m', &
             missing_value=missing_value)

    enddo
endif

id_snow = register_diag_field ( mod_name, 'snow',                 &
                                 (/axes(1:2)/), Time,            &
                                'Surface CO2 ice', 'kg/m/m',     &
                                 missing_value=missing_value )

id_alb = register_diag_field ( mod_name, 'alb',                   &
                                 (/axes(1:2)/), Time,            &
                                'Current Surface Albedo', ' ',         &
                                 missing_value=missing_value )

id_alb0 = register_static_field ( mod_name, 'alb0',               &
                                 (/axes(1:2)/),                  &
                                'Initial Surface Albedo', ' ',         &
                                 missing_value=missing_value )

id_thin0 = register_static_field ( mod_name, 'thin0',                 &
                                 (/axes(1:2)/),                      &
                                'Initial Surface Thermal Inertia', 'mks',    &
                                 missing_value=missing_value )

id_emiss0 = register_static_field ( mod_name, 'emiss0',               &
                                 (/axes(1:2)/),                      &
                                'Initial Surface Emissivity', ' ',           &
                                 missing_value=missing_value )

id_ruff0 = register_static_field ( mod_name, 'zruff0',                &
                                 (/axes(1:2)/),                      &
                                'Initial Surface Roughness', ' ',            &
                                 missing_value=missing_value )

id_gice0 = register_static_field ( mod_name, 'gice0',                &
                                 (/axes(1:2)/),                      &
                                'Initial GRS Ice', ' ',            &
                                 missing_value=missing_value )

id_npc0  = register_static_field ( mod_name, 'npc0',                &
                                 (/axes(1:2)/),                      &
                                'Initial npcflag', ' ',            &
                                 missing_value=missing_value )



!   ------------Read input surface thermal inertia and albedo field data--------------
!
!          2 possibilities:  use resolution-specific input files
!                  or calculate new fields from high-resolution input

filename= 'INPUT/thermal_inertia.nc'
fieldname= 'thin'

if( file_exists( trim( filename ) ) ) then
#ifdef CUBE_CORE
    if(mcpu0) print *, 'Calling read_cube_sfc_field:, ', nlon, mlat
    call read_cube_sfc_field( nlon, mlat, filename, fieldname, thermal_inertia )
#else
    call read_sfc_field( nlon, mlat, lonb, latb, filename, fieldname, thermal_inertia )
#endif
    if(mcpu0) print *, 'Have read thermal inertia data file: '
else     !       else use default value from namelist
    thermal_inertia(:,:)= soil_ti
endif


if( np_cap_ti > 0.0 )  then
    if(mcpu0)  print *, 'Resetting TI north of latitude ', np_cap_lat, np_cap_ti
    where( lat > np_cap_lat*pi/180.0 )
        thermal_inertia= np_cap_ti
    end where
endif

where( lat > np_cap_lat*pi/180.0  .and. thermal_inertia > np_cap_ti_max )
    thermal_inertia= np_cap_ti_max
end where

!              Fudge to test Variations in TI in the Southern Hemisphere
if( rescale_sh_ti > 0.0 )  then
    if(mcpu0)  print *, 'Resetting TI in the SH ', rescale_sh_ti
    where( lat > -50.0*pi/180.0 .and. lat < 30.0*pi/180.0 )
        thermal_inertia= rescale_sh_ti * thermal_inertia
    end where
endif




!   ------------Read input surface albedo data--------------

filename= 'INPUT/albedo.nc'
fieldname= 'albedo'

if( file_exists( trim( filename ) ) ) then
#ifdef CUBE_CORE
    call read_cube_sfc_field( nlon, mlat, filename, fieldname, sfc_albedo )
#else
    call read_sfc_field( nlon, mlat, lonb, latb, filename, fieldname, sfc_albedo )
#endif
    if(mcpu0) print *, 'Have read albedo data file: '

else  !       else use default value  (= 0.25)
    sfc_albedo(:,:)= soil_alb
endif



!   ------------Read input surface emissivity data--------------

filename= 'INPUT/sfc_emiss.nc'
fieldname= 'emis'

if( file_exists( trim( filename ) ) ) then
#ifdef CUBE_CORE
    call read_cube_sfc_field( nlon, mlat, filename, fieldname, sfc_emiss )
#else
    call read_sfc_field( nlon, mlat, lonb, latb, filename, fieldname, sfc_emiss )
#endif
    if(mcpu0) print *, 'Have read surface emissivity data file: '

else  !       else use default value  (= 1.0 )
    sfc_emiss(:,:)= 1.0
endif


!   ------------Read input surface roughness data--------------

filename= 'INPUT/sfc_roughness.nc'
fieldname= 'zruf'

if( file_exists( trim( filename ) ) ) then
#ifdef CUBE_CORE
    call read_cube_sfc_field( nlon, mlat, filename, fieldname, sfc_roughness )
#else
    call read_sfc_field( nlon, mlat, lonb, latb, filename, fieldname, sfc_roughness )
#endif
    if(mcpu0) print *, 'Have read surface roughness data file: '

else  !       else use default value  ( from namelist )
    sfc_roughness(:,:)= zoland
    if(mcpu0) print *, 'Setting roughness to constant '
endif

!   ------------Read input surface topography data--------------
!!!!    Note that this is not necessarily the same as the topography used by the dynamics
!!   For example, this is not being filtered :
!   This field is currently used by the dust_source_sink module

filename= 'INPUT/mars_topo.nc'
fieldname= 'topo'

if( file_exists( trim( filename ) ) ) then
    if(mcpu0) print *, 'Topo is there'
#ifdef CUBE_CORE
    call read_cube_sfc_field( nlon, mlat, filename, fieldname, sfc_topo )
#else
    call read_sfc_field( nlon, mlat, lonb, latb, filename, fieldname, sfc_topo )
#endif
    if(mcpu0) print *, 'Have read surface topography data file: '

else  !       else use default value
    sfc_topo(:,:)= 0.0
    if(mcpu0) print *, 'Setting topo to 0'
endif

!   ------------Read input subsurface ice data--------------

filename= 'INPUT/GRS_interp16.nc'
fieldname= 'grs_indices'  !indices contains the latitude of ice boundary

if( file_exists( trim( filename ) ) ) then
#ifdef CUBE_CORE
    call read_cube_sfc_field( nlon, mlat, filename, fieldname, grs_ice )
#else
    call read_sfc_field( nlon, mlat, lonb, latb, filename, fieldname, grs_ice )
#endif
    if(mcpu0) print *, 'Have read GRS ice data file: '

else  !       else use default value  (= 0.0)
    grs_ice(:,:)= 0.0
endif


if (mcpu0) print*,'lon=0 ice values are: ',grs_ice(1,:)

!   ------------Read input subsurface ice data--------------

filename= 'INPUT/npcflag8.nc'
fieldname= 'npcflag'

if( (file_exists( trim( filename ) ))  ) then
#ifdef CUBE_CORE
    call read_cube_sfc_field( nlon, mlat, filename, fieldname, npcflag )
#else
    call read_sfc_field( nlon, mlat, lonb, latb, filename, fieldname, npcflag )
#endif
    if(mcpu0) print *, 'Have read npcflag data file: '
    if (.not. (file_exists( trim( f_rst ))) ) then
!!----reset the surface water ice to follow the npc flag file----!!
        do nt=1,nice_mass
            where (npcflag .gt. 0.5)
                sfc_frost_mom(:,:,nt) = init_sfc_frost
                sfc_frost(:,:,1) = init_sfc_frost
                sfc_frost_blk(:,:,1) = init_sfc_frost
            elsewhere
                sfc_frost_mom(:,:,nt) = 0.
                sfc_frost(:,:,1) = 0.
                sfc_frost_blk(:,:,1) = 0.
            end where
        enddo
    endif
else  !       else use default value  (= 0.0)
    npcflag(:,:)= 0.0
endif


if (mcpu0) then
    print*,'lon=0 npc values are: ',npcflag(1,:)
    if (all(sfc_frost_mom(:,:,1).eq.sfc_frost(:,:,1))) print*,'sfc_frost_mom = sfc_frost'
endif


!---            Write out time independent fields  (albedo and TI)

if (id_alb0 > 0) used = send_data(id_alb0, sfc_albedo, Time)

if (id_thin0 > 0) used = send_data(id_thin0, thermal_inertia, Time)

if (id_emiss0 > 0) used = send_data(id_emiss0, sfc_emiss, Time)

if (id_ruff0 > 0) used = send_data(id_ruff0, sfc_roughness, Time)

if (id_gice0 > 0) used = send_data(id_gice0, grs_ice, Time)

if (id_npc0  > 0) used = send_data(id_npc0, npcflag, Time)

if(mcpu0) print *, 'Have completed soil module initialization '


module_is_initialized  = .true.

!-----------------------------------------------------------------------

end subroutine mars_surface_init


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


subroutine progts( is, js, dt, Time, lon, lat, ps,  phalf, dnflx, tgrnd, snowin, subday, tg_dt, &
                     shflx, dsens_datmos, dsens_dsurf, evap, devap_datmos, devap_dsurf)
!
!  This is the main soil model prediction scheme
!
!  The soil heat diffusion equation is currently solved using the absolutely stable
!          T(n+1)-T(n)=  delt * D*D { T(n+1) }
!
!  Alternatively, one can formulate a Crank-Nicholson scheme:
!
!          T(n+1)-T(n)=  delt * D*D { ( T(n+1)+T(n) )/2 }
!
!  D*D{} is the difference operator for diffusion ( a tridiagonal matrix )
!
!
!  shflx is the sensible heat flux across the atmosphere/soil boundary
!  dnflx and shflx are positive when directed downward
!
!  hence sensible heat flux
!      shflx(:,:)= -cpair*dragh(:)*( tsfc(:)-theta(:) )
!
integer, intent(in)  :: is, js

!   is, ie and js, je are starting/stopping indices on the computational domain
!    they range over the lat/lon points assigned to a particular processor;
!     This allows indexing into arrays like thermal_inertia and sfc_albedo
!
!     (this is set by the physics window,  which typically is a subset of the
!        compututational domain indices)


real, intent(in)                        :: dt                   ! time step
type(time_type), intent(in)             :: Time                 ! model time
real, intent(in),    dimension(:,:)     :: lon                  ! longitude array [rad]
real, intent(in),    dimension(:,:)     :: lat                  ! latitude array [rad]

real, intent(in), dimension(:,:)        ::  dnflx               ! downwelling radiation at surface [W/m^2]
real, intent(in), dimension(:,:)        ::  ps                  ! surface pressure [Pa]
real, intent(in), dimension(:,:,:)      ::  phalf               ! layer boundary pressures [Pa]

real, intent(in), dimension(:,:,:)      ::  tgrnd               ! soil temperatures [K]
real, intent(in), dimension(:,:)        ::  snowin              ! initial CO2 ice on ground [kg/m^2]
real, intent(out), dimension(:,:)       ::  subday              ! CO2 sublimation/condensation [kg/m^2]
real, intent(out), dimension(:,:,:)     ::  tg_dt               ! soil temperature tendency [K/s]
real, intent(in), dimension(:,:)        ::  shflx               ! sensible heat flux [W/m^2]
real, intent(in), dimension(:,:)        ::  dsens_datmos, &     ! derivative of heat flux wrt atmosphere temperature
                                            dsens_dsurf         ! derivative of heat flux wrt surface temperature

real, intent(in), dimension(:,:)        ::  evap                ! evaporative latent heat flux [W/m^2]
real, intent(in), dimension(:,:)        ::  devap_datmos, &     ! derivative of heat flux wrt atmosphere temperature
                                            devap_dsurf         ! derivative of heat flux wrt surface temperature

!  --------Local -----------------
logical  :: used
integer  ::  ie, je, id, jd, i, j, k, nlay, kk

real, dimension(0:size(tgrnd,3)-1) :: adelz, bdelz, cdelz
real, dimension(0:size(tgrnd,3)  ) :: gkz

real, dimension(size(ps,1),size(ps,2),0:size(tgrnd,3)-1) :: gk3d, rhs, adiag, &
                                                bdiag, cdiag, asav, bsav, csav, &
                                          rho3d, cp3d, crd3d

real, dimension(size(ps,1),size(ps,2)) :: irflx, tcrit, soilp, zo, &
                                          msink, tsfc, gk, grs, snow_orig, &
                                          snowtmp, subday_orig

real, dimension(size(ps,1),size(ps,2)) :: coszen, frost, frost_blk, frost_mom, albedo, emiss, delp

real, dimension(size(ps,1),size(ps,2),size(tgrnd,3)) :: tg,torig,sflux,tgtmp

real  ::  delt, cpair, sdelt

real  :: gk_ice, gk_ice2, tice_thresh, gkco2

real  :: d2ice, icefac

real, dimension(size(ps,1),size(ps,2)) :: nonlatent_flux, maxlatent, snowthick, rlatent
integer, dimension(size(ps,1),size(ps,2)) :: istart
logical, dimension(size(ps,1),size(ps,2)) :: redo_calc, check1, check2
logical :: iter


mcpu0=  (mpp_pe() == mpp_root_pe())
nlay = size(tgrnd,3)-1

k = size(phalf,3)
delp(:,:)=phalf(:,:,k)-phalf(:,:,k-1)

id= size(ps,1);    jd= size(ps,2)

ie= is + id - 1;   je= js + jd - 1


delt= dt
torig=tgrnd(:,:,:)
snow_orig=snowin(:,:)
subday=0.d0
subday_orig=subday(:,:)
iter=.true.
adiag=0.d0
bdiag=0.d0
cdiag=0.d0
rhs=0.d0
msink=0.d0
tgtmp=tgrnd(:,:,:)
tg_dt=0.d0
sflux=0.d0


!          copy globally-defined thermal inertia into local array
soilp(:,:)= thermal_inertia(is:ie,js:je)
grs(:,:)  = grs_ice(is:ie,js:je)

if (use_legacy_soil) then
    rho3d = rho_legacy_soil
    cp3d  = cp_legacy_soil
!           Impose soil ice TI of 1000 mks units
    gk_ice=   (  1100.0   /(rho_legacy_soilice * cp_legacy_soilice) )**2   !use this for legacy north ice
    gk_ice2=  (  2236.995   /(rho_legacy_soilice * cp_legacy_soilice) )**2   !use this for legacy south ice
else
    rho3d = rho_ground
    cp3d = cp_ground
    gk_ice=   (  gk1   /(rho_ground * cp_ground) )**2  !original is 1000
    gk_ice2=  (  gk2   /(rho_ground * cp_ground) )**2 
endif
!         Form soil heat diffusivity from soil thermal inertia
gk(:,:)= (soilp(:,:)/(rho3d(:,:,1) * cp3d(:,:,1) ))**2       !Urata 4/23/19. legacy soil values


!           Impose soil ice TI of 1000 mks units
gkco2 = (5000.0/(rho_co2*cp_co2))**2    !set CO2 thermal conductivity to a large number

cpair =  CP_AIR

tg(:,:,:)= tgrnd(:,:,:)
tsfc(:,:)= tgrnd(:,:,1)

zo(:,:)= sfc_roughness(is:ie,js:je)       !   Surface roughness


do k= 1, nlay
    gkz(k)= 2.0/( delzg(k)+delzg(k+1) )
    gk3d(:,:,k)= gk(:,:)
enddo
gk3d(:,:,0)=gk(:,:)

if( do_subsfc_ice  ) then
!           impose soil ice ti

    select case ( subsfc_ice_case )

    case ( 7 )
        do j= 1, jd
            do i= 1, id
                if(  lat(i,j) .ge. 0.) then
                    d2ice=icedepth(lat(i,j),sum(delzg(1:4)),grs(i,j))
                    do k=2,nlay
                        if (d2ice .lt. sum(delzg(1:k))) then
                            !for smooth transition use exponential
                            !                     do kk=2,k
                            !                      icefac=min(1.,exp(5.*((float(kk)-2.)/float(k)-1.))  )
                            !                      gk3d(i,j,kk)=icefac*gk_ice+(1.-icefac)*gk(i,j)
                            !                     enddo
                            !for legacy transition
                            gk3d(i,j,k:nlay)=gk_ice    !gk_ice2*5.10
                            if (use_legacy_soil) then
                                rho3d(i,j,k:nlay)=rho_legacy_soilice
                                cp3d(i,j,k:nlay) =cp_legacy_soilice
                            endif
                            exit
                        endif
                    enddo
                else
                        d2ice=icedepth(lat(i,j),sum(delzg(1:d2is)),grs(i,j))
                        do k=2,nlay
                            if (d2ice .lt. sum(delzg(1:k))) then
                            !for smooth transition use exponential
                            !                     do kk=3,k
                            !                      icefac=min(1.,exp(5.*((float(kk)-2.)/float(k)-1.))  )
                            !                      gk3d(i,j,kk)=icefac*gk_ice2+(1.-icefac)*gk(i,j)
                            !                     enddo
                            !for legacy transition
                            gk3d(i,j,k:nlay)=gk_ice2   !gk_ice2*2.80
                            if (use_legacy_soil) then
                                rho3d(i,j,k:nlay)=rho_legacy_soilice
                                cp3d(i,j,k:nlay) =cp_legacy_soilice
                            endif
                            exit
                        endif
                    enddo
                endif
            enddo
        enddo

    end select



    adiag(:,:,1)= 0.0
    bdiag(:,:,1)= gk3d(:,:,1)*gkz(1)
    cdiag(:,:,1)= gk3d(:,:,1)*gkz(1)
    do k= 2, nlay
        adiag(:,:,k)=    gk3d(:,:,k-1)*gkz(k-1)
        bdiag(:,:,k)=    gk3d(:,:,k-1)*gkz(k-1) + gk3d(:,:,k)*gkz(k)
        cdiag(:,:,k)=                             gk3d(:,:,k)*gkz(k)
    enddo

    do k= 1, nlay
        adiag(:,:,k)=     - alpha*delt*adiag(:,:,k)/delzg(k)
        bdiag(:,:,k)= 1.0 + alpha*delt*bdiag(:,:,k)/delzg(k)
        cdiag(:,:,k)=     - alpha*delt*cdiag(:,:,k)/delzg(k)
    enddo
    k= nlay
    cdiag(:,:,k)=             - delt*gk3d(:,:,k)*gkz(k)/delzg(k)


else         !               uniform subsurface thermal inertia
    adelz(1)= 0.0
    bdelz(1)= gkz(1)/delzg(1)
    cdelz(1)= gkz(1)/delzg(1)
    do k= 2, nlay
        adelz(k)=     gkz(k-1)         /delzg(k)
        bdelz(k)=   ( gkz(k-1)+gkz(k) )/delzg(k)
        cdelz(k)=              gkz(k)  /delzg(k)
    enddo

    do k= 1, nlay
        adiag(:,:,k)=     - gk(:,:)*alpha*delt*adelz(k)
        bdiag(:,:,k)= 1.0 + gk(:,:)*alpha*delt*bdelz(k)
        cdiag(:,:,k)=     - gk(:,:)*alpha*delt*cdelz(k)
    enddo
    k= nlay
    cdiag(:,:,k)=        - gk(:,:)*delt*cdelz(k)

endif       ! -----------------  end of tridiagonal matrix formulation

!                surface ir flux
!               =====================

! TB18c : high threshold value set in namelist if no changes of albedo wanted for frost
if (do_waterice_albedo) then
    frost(:,:)= sfc_frost(is:ie,js:je,1)
    frost_mom(:,:)= sfc_frost_mom(is:ie,js:je,1)
    frost_blk(:,:)= sfc_frost_blk(is:ie,js:je,1)
else
    frost(:,:)= 0.
    frost_mom(:,:)= 0.
    frost_blk(:,:)= 0.
endif
coszen(:,:)= 0.0
if (do_moment_water) then
    call albedo_calc( is, js, lon, lat, &
                coszen, tsfc, ps, snowin, frost_mom, albedo, emiss )
elseif (do_bulk_water) then
    call albedo_calc( is, js, lon, lat, &
                coszen, tsfc, ps, snowin, frost_blk, albedo, emiss )
else
    call albedo_calc( is, js, lon, lat, &
                coszen, tsfc, ps, snowin, frost, albedo, emiss )
endif

asav=adiag
bsav=bdiag
csav=cdiag

tcrit= 3182.48 / ( 23.3494 - log( ps*1.d-2 ) )
! modify maxlatent by current surface temperature difference from critical
! if surface is very warm with snow because of precip, maxlatent will be very negative
maxlatent = snowin*co2_lheat/delt-rho3d(:,:,1)*cp3d(:,:,1)*delzg(1)*(tsfc-tcrit)/delt

!--------------------------------------------------------
!  loop back to here if surface temperature dips below CO2 frost point
100 continue

istart=1
redo_calc=.false.
snowtmp=snowin+subday
where (snowtmp .gt. 0.)
    istart= 0
end where
do k=1,nlay
    crd3d(:,:,k)=1.0/( cp3d(:,:,k) * rho3d(:,:,k) * delzg(k) )
    sflux(:,:,k)=-gk3d(:,:,k)*rho3d(:,:,k)*cp3d(:,:,k)*gkz(k)*(tg(:,:,k)-tg(:,:,k+1))
enddo


irflx(:,:)= STEFAN*emiss(:,:)*tsfc(:,:)**4
nonlatent_flux = dnflx  + shflx  - irflx


!  Add latent heat term to implicit equation
!  5 cases:
!  1 - no snow on ground, nonlatent_flux is positive
!      --> normal situation, no latent heat flux
!  2 - no snow on ground, nonlatent_flux is negative
!      --> surface temperature cools until tcrit
!      --> snow may form on ground during adjustment phase, moves to case 5
!  3 - snow on ground, nonlatent_flux is positive, but less than maxlatent
!      --> CO2 sublimes at surface
!      --> nonlatent balances latent
!      --> new tsfc = old tsfc = tcrit
!  4 - snow on ground, nonlatent_flux is positive, but more than maxlatent
!      --> all CO2 sublimes at surface
!      --> nonlatent is greater than latent
!      --> new tsfc > tcrit
!  5 - snow on ground, nonlatent_flux is negative
!      --> CO2 condenses
!      --> new tsfc = old tsfc = tcrit


if( .not. do_co2_condensation ) then
    rhs(:,:,1)= delt*crd3d(:,:,1)*( nonlatent_flux + (1.-alpha)*sflux(:,:,1) )   &
        +     tg(:,:,1) &
        +     delt*crd3d(:,:,1)*(4.0*irflx+dsens_dsurf*tsfc)
    bdiag(:,:,1)= bdiag(:,:,1) + 4.0*delt*crd3d(:,:,1)*stefan*emiss*tsfc**3
    bdiag(:,:,1)= bdiag(:,:,1) + delt*crd3d(:,:,1)*dsens_dsurf
else
    where (snowtmp .eq. 0.)
        ! case 1, 2 -  original fv3 formulation
        rhs(:,:,1)= delt*crd3d(:,:,1)*( nonlatent_flux )   +     tg(:,:,1)
        rhs(:,:,1)= rhs(:,:,1) + 4.0*delt*crd3d(:,:,1)*irflx
        !       eliminate   sensible heating contribution from tsfc on rhs  and move to bdiag  (fully implicit)
        rhs(:,:,1)= rhs(:,:,1) + delt*crd3d(:,:,1)*dsens_dsurf *tsfc
        !    add explicit contribution
        rhs(:,:,1)= rhs(:,:,1) + (1.-alpha)*delt*crd3d(:,:,1)*sflux(:,:,1)
        !    rjw     add implicit contribution to bdiag from surface radiation and sensible heat flux
        bdiag(:,:,1)= bdiag(:,:,1) + 4.0*delt*crd3d(:,:,1)*stefan*emiss*tsfc**3
        bdiag(:,:,1)= bdiag(:,:,1) + delt*crd3d(:,:,1)*dsens_dsurf
    elsewhere ((nonlatent_flux+sflux(:,:,1)) .lt. maxlatent)
        ! cases 3, 5
        ! limit the amount of condensation based on layer thickness of bottom layer
        msink = msink-(nonlatent_flux+sflux(:,:,1))*delt/co2_lheat
        subday= msink
        snowtmp= snowtmp+subday
        istart=0
        snowthick=snowtmp/rho_co2
        gk3d(:,:,0)  = 2./(snowthick+delzg(1))   !this is actually using gk3d as gkz, not gk to track thickness
        rhs(:,:,0)  = tcrit
        bdiag(:,:,0)= 1.
        rhs(:,:,1)  = tcrit
        adiag(:,:,1)= 0. !-alpha*delt*gkco2 *gk3d(:,:,0)/delzg(1)
        bdiag(:,:,1)= 1. !bdiag(:,:,1) + alpha*delt*(gkco2*gk3d(:,:,0))/delzg(1)
        cdiag(:,:,1)= 0. !
    elsewhere
    ! Case 4
        msink = msink-(nonlatent_flux+sflux(:,:,1))*delt/co2_lheat
        subday= -snowtmp
        tsfc=tg(:,:,1)-snowtmp*co2_lheat*crd3d(:,:,1)
        irflx = stefan*emiss *tsfc **4
        nonlatent_flux = dnflx  + shflx  - irflx + (1.-alpha)*sflux(:,:,1)
        rhs(:,:,1)=tsfc + (delt*crd3d(:,:,1))*( nonlatent_flux + &
                4.0*irflx+ dsens_dsurf*tsfc)
        bdiag(:,:,1)=bdiag(:,:,1)+(delt*crd3d(:,:,1))*(4.0*stefan*emiss*tsfc**3 + dsens_dsurf )
        istart=1
        snowtmp=0.
    end where  !snow
endif  !do_co2_condensation


do k= 2, nlay
    rhs(:,:,k)=  tg(:,:,k) + (1.-alpha)*delt*(sflux(:,:,k)-sflux(:,:,k-1))*crd3d(:,:,k)
enddo
!rhs(:,:,nlay)= rhs(:,:,nlay) - cdiag(:,:,nlay)*tg(:,:,nlay+1)  &
!         - (1.-alpha)*delt*sflux(:,:,nlay)*crd3d(:,:,nlay)
rhs(:,:,nlay) = tg(:,:,k) + ( adiag(:,:,nlay) + bdiag(:,:,nlay) )*crd3d(:,:,nlay)  &
              + cdiag(:,:,nlay)*crd3d(:,:,nlay)*tg(:,:,nlay+1) 

call trdslv( nlay, id*jd, adiag, bdiag, cdiag, rhs, 0, istart )

!          update soil temperatures:   note that tg(:,:,nlay+1) remains fixed
do k= 1, nlay
    tgtmp(:,:,k)= rhs(:,:,k)
enddo
where (snowtmp.eq.0.)
    tsfc= tgtmp(:,:,1)
elsewhere
    tsfc= tcrit
end where

if( do_co2_condensation ) then
    if (iter) then
        tcrit= 3182.48 / ( 23.3494 - log( ps*1.d-2 ) )
        rlatent= cp3d(:,:,1)*rho3d(:,:,1)*delzg(1) / co2_lheat
        msink(:,:)= ( tcrit(:,:) - tsfc(:,:) )*rlatent
        !       consider condensation at the surface; recalculate temperatures with condensation
        where (( msink >  0.0 ) .and. (snowtmp .eq. 0.))
            snowtmp=1.d-10
            redo_calc=.true.
        end where

        if (any(redo_calc)) then
            tg=torig
            tgtmp=tg
            tsfc =tg(:,:,1)
            adiag=asav
            bdiag=bsav
            cdiag=csav
            msink=0.d0
            subday(:,:)=subday_orig(:,:)
            iter=.false.
            where (.not.(snowtmp .eq. 1.d-10)) snowtmp=snow_orig
            where (snowtmp.eq.1.d-10) msink=snowtmp
            where (snowtmp.eq.1.d-10) subday(:,:)=snowtmp
            goto 100
        endif
    endif !iter
else
    subday= 0.0
    snowtmp= 0.0

endif

!      Finally, update final surface temperature
tgtmp(:,:,1)= tsfc(:,:)
tg_dt=tgtmp-tgrnd


if (id_ts > 0)     used = send_data ( id_ts,     tsfc,   Time, is, js )
if (id_subsfc > 0) used = send_data ( id_subsfc, tgtmp,  Time, is, js )
if (id_sflux > 0)  used = send_data ( id_sflux,  sflux,  Time, is, js )
if (id_thin > 0)   used = send_data ( id_thin,   soilp,  Time, is, js )
if (id_snow > 0)   used = send_data ( id_snow,   snowtmp,   Time, is, js )
if (id_alb > 0)    used = send_data ( id_alb,    albedo, Time, is, js )

end subroutine progts

!--------------------------------------------------------
!--------------------------------------------------------

real function icedepth(curlat,min_d,lat1)
!
!  calculate the ice transition depth as a tanh distribution
!  calculates depth as negative height, so convert to positive depth before return
implicit none
real, intent(in) :: curlat, &   !  current latitude [rad]
                    min_d, &    !  minimum depth [m]
                    lat1        !  transition latitude [rad]
real, parameter :: botdepth=40. !  maximum depth to ice
real :: beta, alpha, degrad
real, parameter :: width = 25.  !  width of transition zone
real :: pm                      !  saves the sign of the current lat

degrad=pi/180.
pm=sign(1.,curlat)
! determines width of transition zone
alpha = pm*width*degrad/(pi)                
! determines latitude of ice boundary. puts inflection point at lat-x
beta = pm*(abs(lat1)-width/1.25)*degrad  
! the following is negative   
icedepth=tanh((curlat-beta)/alpha)-1.0      
! scale to the maximum depth where ice starts, and subtract min_d, convert to positive, flatten out at 10m
icedepth=min(-(botdepth*icedepth/2.-min_d),10.)   

end function icedepth

!--------------------------------------------------------
!--------------------------------------------------------

subroutine albedo_calc( is, js, lon, lat,  cosz, &
                         ts, ps, snow, frost, albedo, emiss_sfc )
!
!   calculate the effect of ice onsurface albedo 
!       In principal, albedo can be wavelength and zenith-angle dependent
!       Would need to distinguish between direct and diffuse solar radiation
!
!
integer, intent(in)                  :: is, js
real, intent(in), dimension(:,:)     :: lon
real, intent(in), dimension(:,:)     :: lat
real, intent(in), dimension(:,:)     :: cosz       !
real, intent(in), dimension(:,:)     :: ts         !
real, intent(in), dimension(:,:)     :: ps         !
real, intent(in), dimension(:,:)     :: snow
real, intent(in), dimension(:,:)     :: frost

real, intent(out), dimension(:,:)    :: albedo     !
real, intent(out), dimension(:,:)    :: emiss_sfc  !


! Local variables

logical  :: used
integer  ::  ie, je, id, jd
real     ::  snow_threshold

real,   dimension(size(lat,1),size(lat,2)) :: albedo_ice, emiss_ice

id= size(ps,1);     jd= size(ps,2)

ie= is + id - 1;    je= js + jd - 1

!          copy globally-defined albedo into local array
albedo(:,:)= sfc_albedo(is:ie,js:je)

emiss_sfc(:,:)= sfc_emiss(is:ie,js:je)

!             np and sp albedo and emiss values obtained from namelist input
where(  lat > 0.0 )
    albedo_ice = albedo_ice_np
    emiss_ice=   emiss_ice_np
elsewhere
    albedo_ice = albedo_ice_sp
    emiss_ice=   emiss_ice_sp
end where


!           modify albedo and emissivity where co2 ice is present
snow_threshold= 2.0

where ( snow > snow_threshold )
    albedo = albedo_ice
    emiss_sfc = emiss_ice
end where

!      include an albedo dependence on frost;
!            frost_threshhold (via namelist) in m
!              frost (kg/m**2) ;  hence the 1000.0 conversion
where ( frost > frost_threshold * 1.0e3 .and. snow .le. 0. .and. ts .lt. 273.)
    albedo = albedo_h2o
    emiss_sfc = emiss_h2o
else where ( frost > frost_threshold * 1.0e3 .and. snow .le. 0. .and. ts .ge. 273.)
    albedo = albedo_h2o_liquid
    emiss_sfc = emiss_h2o
end where


end  subroutine albedo_calc


!--------------------------------------------------------
!--------------------------------------------------------


#ifdef CUBE_CORE
subroutine read_cube_sfc_field(  npx, npy, filename, field_name, sfc_field )
!
!      Interpolate a 2-D data field to the model grid
!          blon and blat are the bounding lon and latitudes of the model grid
!      This is a simplified version of the routine in fv_surf_map

use fv_arrays_mod,          only: fv_atmos_type, fv_nest_type, fv_grid_bounds_type

use fv_mars_interface_mod,  only: get_cubed_sphere_mars, grid, agrid, bd, npx_global, &
                                  is, js, ie, je, isd, jsd, ied, jed

use fv_surf_map_mod,        only: map_to_cubed_simple


integer,            intent(in)    :: npx, npy
character(len=128), intent(in)    :: filename, field_name
real,               intent(out)   :: sfc_field(is:ie, js:je)

type(fv_grid_bounds_type) :: bd2

!        local
real                :: phis(isd:ied, jsd:jed)

real*4, allocatable :: htopo(:,:)
real,   allocatable :: rtopo(:,:)
real,   allocatable :: lon1(:),  lat1(:)
integer             :: nlon, nlat, i, j, n, nstr,  fld_dims(4)
real                :: dx1, dy1

character(len=128)  :: tilefile

if(mcpu0) write(*,*) '  Cube face dims ', npx, npy

bd2= bd

if(mcpu0) write(*,*) '  Cube indices: bd struct ', bd2%is, bd2%ie, bd2%js, bd2%je

if(mcpu0) write(*,*) '  Cube indices: ie, is... ', is, ie, js, je

if(mcpu0) write(*,*) '  Cube:  npx_global  ', npx_global


if(mcpu0) write(*,*) '  Cube:  agrid lon   ', agrid(1:5,js,1)*180/pi

if(mcpu0) write(*,*) '  Cube:  agrid lat  ', agrid(1:5,1:5,2)*180/pi

bd2%is= is
bd2%ie= ie
bd2%js= js
bd2%je= je


if ( npx > 0 ) then



    nstr= len_trim( filename )
    tilefile= filename(1:nstr-2) //  'tile1.nc'

    !        If tiled surface input fields are present, check these.
    if( file_exists( trim(tilefile) ) ) then

        call field_size( trim(tilefile), trim(field_name), fld_dims )
        nlon= fld_dims(1);  nlat= fld_dims(2);

        if(mcpu0) print *, 'Reading tileX.nc:, ', nlon, nlat

        if( nlon==npx-1 .and.  nlat==npy-1 )  then
            call read_data(trim(filename), trim(field_name),  sfc_field )
            if(mcpu0) print *, 'Have read surface data ', trim(field_name), ' from tiled files: '
            return
        else
            call error_mesg ('mars_surface','tiled input surface field is the wrong dimension', FATAL)
        endif
    endif

    !         Otherwise, will need to interpolate hi-res surface field data to cube domains
    !
    call field_size( trim(filename), trim(field_name), fld_dims )

    nlon= fld_dims(1);  nlat= fld_dims(2);

    if(mcpu0) write(*,*) trim(filename), '  Mars Hi-RES dataset dims=', nlon, nlat

    allocate ( htopo(nlon,nlat) )
    allocate ( rtopo(nlon,nlat) )

    if(mcpu0) write(*,*) trim(filename),' ', trim(field_name),'  Read data topo'
    call read_data( trim(filename), trim(field_name), rtopo, no_domain=.true. )
    if(mcpu0) write(*,*) trim(filename),' ', trim(field_name),'  data topo OK'

    !          This is needed because htopo is declared as real*4
    htopo= rtopo

    allocate ( lat1(nlat+1) )
    allocate ( lon1(nlon+1) )

    if( field_exist(  trim(filename), 'lonb') ) then
        call read_data( trim(filename), 'lonb', lon1, no_domain=.true. )
        lon1(:)= lon1(:)*PI/180.0
    else
        dx1 = 2.*pi/real(nlon)
        do i=1,nlon+1
            lon1(i) = dx1 * (i-1)
        enddo
    endif

    if( field_exist(  trim(filename), 'latb') ) then
        call read_data( trim(filename), 'latb', lat1, no_domain=.true. )
        lat1(:)= lat1(:)*PI/180.0
    else
        dy1 = pi/real(nlat)
        lat1(1) = - 0.5*pi
        lat1(nlat+1) =  0.5*pi
        do j=2,nlat
            lat1(j) = -0.5*pi + dy1*(j-1)
        enddo
    endif
    if(mcpu0) write(*,*) 'latb lonb ok'


    call map_to_cubed_simple( nlon, nlat, lat1, lon1, htopo, grid, agrid,  &
                       phis, npx, npy, npx_global, bd2 )

    if(mcpu0) write(*,*) 'cubed simple ok'

    sfc_field(is:ie,js:je) = phis(is:ie,js:je)

    deallocate ( htopo )
    deallocate ( lat1, lon1 )

else

    sfc_field(:,:) = 1.0

end if

end subroutine read_cube_sfc_field
#endif CUBE_CORE

!--------------------------------------------------------
!--------------------------------------------------------


subroutine read_sfc_field(  im, jm, blon, blat, filename, field_name, sfc_field )
!
!      Interpolate a 2-D data field to the model grid if input field dimensions dont
!           match the global dimensions of the desired surface field (sfc_fld)
!           blon and blat are the bounding lon and latitudes of the model grid needed for interpolation

use  horiz_interp_mod,  only: horiz_interp

integer,            intent(IN)    :: im, jm
real,               intent(in)    :: blon(:,:), blat(:,:)
character(len=128), intent(in)    :: filename, field_name
real,               intent(out)   :: sfc_field(:,:)


real,   allocatable ::  rtopo(:,:)
real,   allocatable ::  lon1(:),  lat1(:)
real                ::  dx1, dy1
integer             ::  nlon, nlat, i, j, n, id, jd, fld_dims(4)


id= size(blon,1)-1; jd= size(blat,2)-1

call field_size( trim(filename), trim(field_name), fld_dims )

nlon= fld_dims(1);  nlat= fld_dims(2);

if(mcpu0) write(*,*) trim(filename),  ':  Input dataset dims=', nlon, nlat

if( nlon==im   .and. nlat==jm ) then

    if(mcpu0) write(*,*) '   --->  ', trim(field_name),  '  Field sizes match:  no interpolation necessary '
    call read_data( trim(filename), trim(field_name), sfc_field )

else   !     will need to interpolate input surface data field to the required lat/lon grid
!          The input file contains lat and lon arrays; however we
!      need the bounding latitude and longitudes: So inquire
!      about their presence; else use equi-spaced 'default' values

    allocate ( rtopo(nlon,nlat) )
    allocate ( lat1(nlat+1) )
    allocate ( lon1(nlon+1) )

    call read_data( trim(filename), trim(field_name), rtopo, no_domain=.true. )


    if( field_exist(  trim(filename), 'lonb') ) then
        call read_data( trim(filename), 'lonb', lon1, no_domain=.true. )
        lon1(:)= lon1(:)*PI/180.0
    elseif ( nlon == id ) then
        if(mcpu0) write(*,*) trim(field_name),  ':  Input lon dims= requested lon dims', nlon, id
        lon1(:)= blon(:,1)
    else
        dx1 = 2.*pi/real(nlon)
        do i=1,nlon+1
            lon1(i) = dx1 * (i-1)
        enddo
    endif

    if( field_exist(  trim(filename), 'latb') ) then
        call read_data( trim(filename), 'latb', lat1, no_domain=.true. )
        lat1(:)= lat1(:)*PI/180.0
    elseif ( nlat == jd ) then
        if(mcpu0) write(*,*) trim(field_name),  ':  Input lat dims= requested lon dims', nlon, id
        lat1(:)= blat(:,1)
    else
        dy1 = pi/real(nlat)
        lat1(1) = - 0.5*pi
        lat1(nlat+1) =  0.5*pi
        do j=2,nlat
            lat1(j) = -0.5*pi + dy1*(j-1)
        enddo
    endif

    call horiz_interp ( rtopo, lon1, lat1, blon, blat, sfc_field, interp_method= 'conservative' )


    deallocate ( rtopo )
    deallocate ( lat1, lon1 )

endif


end subroutine read_sfc_field

!--------------------------------------------------------
!--------------------------------------------------------

subroutine mars_surface_register_restart(fname,phys_domain)
! register restart field to be written to restart file.
character(len=*),                 intent(in) :: fname
character(len=64)                            :: fname2
type(domain2d),                intent(inout) :: phys_domain
integer :: id_restart, n, ndx
character (len=128) :: tracer_name


call get_mosaic_tile_file(fname, fname2, is_no_domain=.false., domain=phys_domain )

!default restart file: read/write definition
allocate(Surf_restart)
if(trim(fname2) == trim(fname)) then
    Til_restart => Surf_restart
    in_different_file = .false.
else
    in_different_file = .true.
    allocate(Til_restart)
endif

id_restart = register_restart_field(Til_restart, fname, 'snow', sfc_snow, domain=phys_domain,mandatory=.false.)
id_restart = register_restart_field(Til_restart, fname, 'frost', sfc_frost, domain=phys_domain,mandatory=.false.)
id_restart = register_restart_field(Til_restart, fname, 'frost_blk', sfc_frost_blk, domain=phys_domain,mandatory=.false.)
id_restart = register_restart_field(Til_restart, fname, 'cprecip_blk', cumulative_prec_blk, domain=phys_domain,mandatory=.false.)
!id_restart = register_restart_field(Til_restart, fname, 'cprecip_mconv', cumulative_prec_mconv, domain=phys_domain,mandatory=.false.)

do n=1,nice_mass
    ndx= ice_mass_indx(n)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)
    id_restart = register_restart_field(Til_restart, fname, trim(tracer_name), sfc_frost_mom(:,:,n), domain=phys_domain,mandatory=.false.)
end do

id_restart = register_restart_field(Til_restart, fname, 'tg', tsoil, domain=phys_domain,mandatory=.false.)
id_restart = register_restart_field(Til_restart, fname, 'soil_ice', soil_icex, domain=phys_domain,mandatory=.false.)
id_restart = register_restart_field(Surf_restart, fname, 'zgrid', zgrid, no_domain = .true.,mandatory=.false.)

! if the restart file used a different soil layer structure. read only
if (rst2) then
    allocate(Surf_restart2)
    if(trim(fname2) == trim(fname)) then
        Til_restart2 => Surf_restart2
    else
        allocate(Til_restart2)
    endif

    id_restart = register_restart_field(Til_restart2, fname, 'snow', sfc_snow, domain=phys_domain,mandatory=.false.)
    id_restart = register_restart_field(Til_restart2, fname, 'frost', sfc_frost, domain=phys_domain,mandatory=.false.)
    id_restart = register_restart_field(Til_restart2, fname, 'frost_blk', sfc_frost_blk, domain=phys_domain,mandatory=.false.)
    id_restart = register_restart_field(Til_restart2, fname, 'cprecip_blk', cumulative_prec_blk, domain=phys_domain,mandatory=.false.)

    do n=1,nice_mass
        ndx= ice_mass_indx(n)
        call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)
        id_restart = register_restart_field(Til_restart2, fname, trim(tracer_name), sfc_frost_mom(:,:,n), domain=phys_domain,mandatory=.false.)
    end do

    id_restart = register_restart_field(Til_restart2, fname, 'tg', tsoil_r, domain=phys_domain,mandatory=.false.)
    id_restart = register_restart_field(Til_restart2, fname, 'soil_ice', soil_icex_r, domain=phys_domain,mandatory=.false.)
    id_restart = register_restart_field(Surf_restart2, fname, 'zgrid', zgrid_r, no_domain = .true.,mandatory=.false.)

endif


end subroutine mars_surface_register_restart

!--------------------------------------------------------
!--------------------------------------------------------


subroutine mars_surface_end( days  )
! end mars surface

integer,             intent(in) :: days
!new way to save restart files
call save_restart(Surf_restart)
if(in_different_file) call save_restart(Til_restart)

deallocate (  tsoil            )
deallocate (  t_surf           )
deallocate (  sfc_snow         )
deallocate (  thermal_inertia  )
deallocate (  sfc_albedo       )
deallocate (  sfc_emiss        )
deallocate (  grs_ice          )
#ifndef SINGLE_TILE
deallocate ( zgrid, delzg )
#endif

end subroutine mars_surface_end

!--------------------------------------------------------
!--------------------------------------------------------


subroutine trdslv( n, mx2, a, b, c, x, iflg, istart )
!
!  Tridiagonal solver
!
!     if  iflg = 0  this subroutine computes the lu-decomposition
!     of the input matrix and then solves for the solution vector  x.
!
!     if  iflg = 1  the calculation of the lu-decomposition is skipped
!     and the solution vector  x  is calculated.

integer, intent(in) :: n, mx2, iflg
integer, dimension(mx2), intent(in)   :: istart
real, dimension(mx2,0:n), intent(in)    :: a
real, dimension(mx2,0:n), intent(inout) :: b, c
real, dimension(mx2,0:n), intent(out)   :: x

integer :: nm1,np1, i, j, is

nm1 = n-1

do i=1,mx2
    is=istart(i)
    np1=is+1
    !   obtain the lu-decomposition
    if( iflg .eq. 0 ) then
        b(i,is) = 1.0/b(i,is)
        c(i,is) = c(i,is)*b(i,is)

        do j=np1,n
            b(i,j) = 1.0/( b(i,j)-a(i,j)*c(i,j-1) )
            c(i,j)= c(i,j)*b(i,j)
        enddo

    end if

    !        ----come here for back-solving----

    x(i,is) = x(i,is)*b(i,is)
    do j=np1,n
        x(i,j) = ( x(i,j)-a(i,j)*x(i,j-1) )*b(i,j)
    enddo

    do j=nm1,is,-1
        x(i,j) = x(i,j)-c(i,j)*x(i,j+1)
    enddo
enddo !mx2 loop

return
end subroutine trdslv

!--------------------------------------------------------
!--------------------------------------------------------



end module mars_surface_mod


