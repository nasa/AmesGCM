module null_physics_mod

! *****************************************************************     
! Module to contain placeholders for Ames Mars physics routines     
! *****************************************************************

use fms_mod, only: error_mesg, FATAL
use time_manager_mod, only: time_type
use   mpp_domains_mod, only: domain2d

implicit none
private

! *****************************************************************     
! ************************ INTERFACES******************************     
! *****************************************************************     

public :: micro_driver, micro_driver_init, dust_update_init,   &
          dust_update, dust_update_end, sfc_dust_mass,         &
          tcol_mass, tcol_core, dust_surf_ini, &
          dust_source_sink, dust_source_init, dust_source_end, &
          sedim_driver, sedim_driver_init, coagul_main, coagul_init
public :: palmer_drag, palmer_drag_init, palmer_drag_end, &
          cg_drag_init, cg_drag_calc, cg_drag_end, &
          topo_drag, topo_drag_init, topo_drag_end
public :: cloud_physics, cloud_physics_init, cloud_physics_end, &
		  id_wcol, id_cldcol, id_cld_r, id_rsat, cldcol, wcol
public :: photochem_driver
public :: init_null_phys
public :: nltecool

integer ::  id_wcol, id_cldcol, id_cld_r, id_rsat
real, dimension(:,:,:),   allocatable, save  ::  sfc_dust_mass
real, dimension(:,:),   allocatable, save  ::  wcol  
real, dimension(:,:),   allocatable, save  ::  cldcol

real, dimension(:,:,:), allocatable :: tcol_mass,tcol_core


real    ::  dust_surf_ini=30.                   !  Initial reservoir of dust for cold cases

contains

! *****************************************************************     
! *****************************************************************

subroutine init_null_phys(is, js, ie, je)
                       
use initracer_mod, only: ndust_mass

integer, intent(in)  :: is, js, ie, je
integer  id, jd


allocate (  sfc_dust_mass (is:ie,js:je,ndust_mass)  )
allocate (  tcol_mass (is:ie,js:je,ndust_mass)  )
allocate (  tcol_core (is:ie,js:je,ndust_mass)  )
allocate (  wcol (is:ie,js:je)  )
allocate (  cldcol (is:ie,js:je)  )
    
end subroutine init_null_phys

! *****************************************************************     
! *****************************************************************

subroutine micro_driver(is,js,Time,p_half,p_full,t,tdt, &
         r,rdt,strss,tg,dtime,drg,rdt_micro,tdt_micro,checkcons)


!===============================================================
!     Input Arguments 
!===============================================================
integer, intent(in)  :: is, js
type(time_type), intent(in)             :: Time !     Model time
real*8, intent(in), dimension(:,:,:) :: p_half  !     layer boundary pressures [Pa]
real*8, intent(in), dimension(:,:,:) :: p_full  !     layer midpoint pressures [Pa]
real*8, intent(in), dimension(:,:,:) :: t       !     Temperature [K]
real*8, intent(in), dimension(:,:,:,:) :: r     !     Tracer field
real*8, intent(in), dimension(:,:) :: strss     !     Surface Stress
real*8, intent(in), dimension(:,:) :: drg       !     Drag

real*8, intent(in), dimension(:,:) :: tg        !     Ground Temperature (K)
real*8, intent(in) :: dtime                     !     Time step
real*8, intent(in), dimension(:,:,:,:) :: rdt   !     Tracer tendencies (kg/kg/s)
real*8, intent(in), dimension(:,:,:) :: tdt     !     Temperature tendencies (K/s)
logical, intent(in) :: checkcons                !     check conservation
!===============================================================
!     Output Arguments
!===============================================================
real*8, intent(out), dimension(:,:,:,:) :: rdt_micro
real*8, intent(out), dimension(:,:,:) :: tdt_micro

call error_mesg('micro_driver','The null version of micro_driver should not be called.',FATAL)


	
end subroutine micro_driver

! *****************************************************************     
! *****************************************************************

subroutine micro_driver_init(nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )
!! Arguments
integer, intent(in)                   :: nlon, mlat, nlevels
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon, lat
integer, intent(in)                   :: axes(4)
type(time_type), intent(in) :: Time
	
	
  call error_mesg('micro_driver_init','The null version of micro_driver_init should not be called.',FATAL)

end subroutine micro_driver_init

! *****************************************************************     
! *****************************************************************


subroutine ames_pbl(is,js,dt,ps,p_half,p_full,z_half,akm,akh,  &
                  u,v,t,tsurf,q_wat,snow,frost,tdtlw,tau_x,tau_y,&
                  sens,evap,&
                  udt,vdt,tdt,qdt, &
                  udt_pbl,vdt_pbl,tdt_pbl,qdt_pbl,&
                  ustar,thstar,cdm,cdh,dsens, &
                  wind,z_pbl,p_pbl,pbl_lev,rkh_out,Time_next)

     integer, intent(in) :: is, js
     type(time_type), intent(in) :: Time_next
     real, intent(in) :: dt
     real, intent(in),    dimension(:,:,:)   :: p_half, p_full, tdtlw, z_half !tdtlw is K/s
     real, intent(in),    dimension(:,:,:)   :: u, v, t
     real, intent(in),    dimension(:,:,:)   :: q_wat
     real, intent(out),   dimension(:,:,:)   :: akm, akh
     real, intent(inout), dimension(size(t,1),size(t,2)) :: snow, frost  !snow is co2, frost is h2o
     real, intent(in),    dimension(size(t,1),size(t,2)) :: tsurf, ps
     real,   dimension(size(t,1),size(t,2)) :: ustar,thstar,cdm,cdh
     real,   dimension(size(t,1),size(t,2)) :: sens, evap,  &
                                             tau_x, tau_y, dsens, &
                                             wind, z_pbl,p_pbl
     integer, intent(out), dimension(size(t,1),size(t,2)) :: pbl_lev
     real, intent(in), dimension(size(t,1),size(t,2),size(t,3))  :: tdt, udt, vdt, qdt
     real, intent(out), dimension(size(t,1),size(t,2),size(t,3)) :: tdt_pbl, udt_pbl, vdt_pbl, qdt_pbl
!    Output eddy mixing coefficient (m2/s)
     real, intent(out), dimension(size(t,1),size(t,2),2*size(t,3)+1) :: rkh_out


  call error_mesg('ames_pbl','The null version of ames_pbl should not be called.',FATAL)
	
end subroutine ames_pbl

! *****************************************************************     
! *****************************************************************

subroutine palmer_drag(delt, uwnd, vwnd, atmp, &
                       pfull, phalf, zfull, zhalf, &
                       dtaux, dtauy, Time    )
real,    intent(in) :: delt
type(time_type), intent(in) :: Time
real, intent(in), dimension(:,:,:) :: uwnd, vwnd, atmp
real, intent(in), dimension(:,:,:) :: pfull, phalf, zfull, zhalf
real, intent(out), dimension(:,:,:) :: dtaux, dtauy


call error_mesg('palmer_drag','The null version of palmer_drag should not be called.',FATAL)
	
end subroutine palmer_drag

! *****************************************************************     
! *****************************************************************

subroutine palmer_drag_init(nlon, mlat, lonb, latb, axes_phys, Time)
integer, intent(in), dimension(4) :: axes_phys
type(time_type), intent(in) :: Time
  integer, intent(in)                   :: nlon, mlat
  real,    intent(in),  dimension(:,:)  :: lonb, latb 
	
call error_mesg('palmer_drag_init','The null version of palmer_drag_init should not be called.',FATAL)


end subroutine palmer_drag_init

! *****************************************************************     
! *****************************************************************

subroutine palmer_drag_end

write(*,*)"Nothing done here"
	
end subroutine palmer_drag_end

! *****************************************************************     
! *****************************************************************

subroutine cloud_physics(is, js, lon, lat, dt, Time, &
                         p_half, p_full, tsfc, frost, t, r, rdt, drag_q)

integer, intent(in)  :: is, js
real,    intent(in)  :: dt
type(time_type), intent(in)             :: Time
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:)     :: lat
real, intent(in),    dimension(:,:,:)   :: p_half, p_full
real, intent(in),    dimension(:,:)     :: tsfc
real, intent(inout), dimension(:,:)     :: frost
real, intent(in),    dimension(:,:,:)   :: t
real, intent(in),    dimension(:,:,:,:) :: r
real, intent(inout), dimension(:,:,:,:) :: rdt
real, intent(in),    dimension(:,:)     :: drag_q


call error_mesg('cloud_physics','The null version of cloud_physics should not be called.',FATAL)

	
end subroutine cloud_physics


! *****************************************************************     
! *****************************************************************

subroutine cloud_physics_init( nlon, mlat, lonb, latb, lon, lat, axes, Time )

   integer, intent(in)                   :: nlon, mlat
   real,    intent(in),  dimension(:,:)  :: lonb, latb 
   real,    intent(in),  dimension(:,:)  :: lon, lat 

   integer, intent(in) :: axes(4)
   type(time_type), intent(in) :: Time
	

call error_mesg('cloud_physics_init','The null version of cloud_physics_init should not be called.',FATAL)

end subroutine cloud_physics_init


! *****************************************************************     
! *****************************************************************

subroutine cloud_physics_end


call error_mesg('cloud_physics_end','The null version of cloud_physics_end should not be called.',FATAL)
	
end subroutine cloud_physics_end


! *****************************************************************     
! *****************************************************************
subroutine cg_drag_init (lonb, latb, pref, Time, axes)
!
!   cg_drag_init is the constructor for cg_drag_mod.
!

real,    dimension(:,:), intent(in)      :: lonb, &     ! 2d array of model longitudes on cell corners [radians]
                                            latb        ! 2d array of model latitudes at cell corners [radians]
real,    dimension(:),   intent(in)      :: pref        ! array of reference pressures at full levels (plus
                                                        ! surface value at nlev+1), based on 1013.25hPa pstar
                                                        ! [ Pa ]
integer, dimension(4),   intent(in)      :: axes        ! data axes for diagnostics
type(time_type),         intent(in)      :: Time        ! current time (time_type)
    

call error_mesg('cg_drag_init','The null version of cg_drag_init should not be called.',FATAL)
end subroutine cg_drag_init

! *****************************************************************     
! *****************************************************************

subroutine cg_drag_calc (is, js, lat, pfull, zfull, temp, uuu, vvv,  &
                         Time, delt, gwfcng_x, gwfcng_y, dtemp)
!  
!    cg_drag_calc defines the arrays needed to calculate the convective
!    gravity wave forcing, calls gwfc to calculate the forcing, returns 
!    the desired output fields, and saves the values for later retrieval
!    if they are not calculated on every timestep.
!
!

!
integer,                intent(in)      :: is, js
real, dimension(:,:),   intent(in)      :: lat          ! array of model latitudes at cell boundaries [radians]
real, dimension(:,:,:), intent(in)      :: pfull, &     ! pressure at model full levels [ Pa ]
                                            zfull, &    ! height at model full levels [ m ]
                                            temp, &     ! temperature at model levels [ deg K ]
                                            uuu, &      ! zonal wind  [ m/s ]
                                            vvv         ! meridional wind  [ m/s ]
type(time_type),        intent(in)      :: Time         ! current time, needed for diagnostics [ time_type ]
real           ,        intent(in)      :: delt         ! physics time step [ s ]
real, dimension(:,:,:), intent(out)     :: gwfcng_x, &  ! time tendency for u eqn due to gravity-wave forcing [ m/s^2 ]
                                            gwfcng_y    ! time tendency for v eqn due to gravity-wave forcing [ m/s^2 ]
real, dimension(:,:,:), intent(out)     :: dtemp        ! time tendency of the temperature in K/s due to gravity-wave forcing [ K/s ]



call error_mesg('cg_drag_calc','The null version of cg_drag_calc should not be called.',FATAL)
end subroutine cg_drag_calc
! *****************************************************************     
! *****************************************************************

subroutine cg_drag_end
!
!    cg_drag_end is the destructor for cg_drag_mod.
!


call error_mesg('cg_drag_end','The null version of cg_drag_end should not be called.',FATAL)


end subroutine cg_drag_end

! *****************************************************************     
! *****************************************************************

subroutine topo_drag                                                   &
                                             ( delt, uwnd, vwnd, atmp, &
                                           pfull, phalf, zfull, zhalf, &
                                            dtaux, dtauy, dtemp, Time )

real,    intent(in) :: delt
type(time_type), intent(in) :: Time

! INPUT
! -----
real, intent(in), dimension(:,:,:) :: uwnd, &   ! UWND     Zonal wind (dimensioned IDIM x JDIM x KDIM)
                                    vwnd, &     ! VWND     Meridional wind (dimensioned IDIM x JDIM x KDIM)
                                    atmp        ! ATMP     Temperature at full levels (IDIM x JDIM x KDIM)
real, intent(in), dimension(:,:,:) :: pfull, &  ! PFULL    Pressure at full levels (IDIM x JDIM x KDIM)
                                    phalf, &    ! PHALF    Pressure at half levels (IDIM x JDIM x KDIM+1)
                                    zfull, &    ! ZFULL    Height at full levels (IDIM x JDIM x KDIM)
                                    zhalf       ! ZHALF    Height at half levels (IDIM x JDIM x KDIM+1)

! OUTPUT
! ------
real, intent(out), dimension(:,:,:) :: dtaux, & ! DTAUX Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
                                    dtauy, &    ! DTAUY Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
                                    dtemp       ! DTEMP Tendency of the temperature in K/s (IDIM x JDIM x KDIM)


call error_mesg('topo_drag','The null version of topo_drag should not be called.',FATAL)
end subroutine topo_drag
! *****************************************************************     
! *****************************************************************

subroutine topo_drag_init (lonb, latb, axes_phys, Time)
! initialize topographic drag module
real, intent(in), dimension(:,:) :: lonb, latb
integer, intent(in), dimension(4) :: axes_phys
type(time_type), intent(in) :: Time
    

call error_mesg('topo_drag_init','The null version of topo_drag_init should not be called.',FATAL)

end subroutine topo_drag_init

! *****************************************************************     
! *****************************************************************

subroutine topo_drag_end


call error_mesg('topo_drag_end','The null version of topo_drag_end should not be called.',FATAL)
end subroutine topo_drag_end

! *****************************************************************     
! *****************************************************************

subroutine photochem_driver(is,ie,js,je,kd,lon,lat,p_half,p_full,delp,  &
                            t,tdt,time,dt,r,rdt,rdt_pchem,do_qmass) 

use constants_mod, only: PI,rdgas,avogno,GRAV


integer, intent(in)  :: is, js, ie, je, kd
logical, intent(in)  :: do_qmass
real,    intent(in)  :: dt
type(time_type), intent(in)             :: Time
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:)     :: lat
real, intent(in),    dimension(:,:,:)   :: p_half, p_full  !(hPa)
real, intent(in),    dimension(:,:,:)   :: delp  !(kg/m/s2)
real, intent(in),    dimension(:,:,:)   :: t, tdt
real, intent(in),    dimension(:,:,:,:) :: r, rdt
real, intent(out),   dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_pchem


call error_mesg('photochem_driver','The null version of photochem_driver should not be called.',FATAL)
end subroutine photochem_driver
! *****************************************************************     
! *****************************************************************

subroutine dust_update (   is, js, lon, lat, dt, Time, &
                              p_half, p_full, p_pbl, k_pbl, tsfc,  snow, &
                              frost, stress, taudust, t, tdt, r, rdt, shflx, &
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
real, intent(in),    dimension(:,:,:,:) :: r,rdt
real, intent(in),    dimension(:,:)     :: shflx  
real, intent(out), dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_dst
real, intent(out),    dimension(:,:)     :: source_mom  
logical, intent(out), dimension(size(r,1),size(r,2),size(r,4)) :: lifting_dust

call error_mesg('dust_update','The null version of dust_update should not be called.',FATAL)
end subroutine dust_update
! *****************************************************************     
! *****************************************************************


subroutine dust_update_init( nlon, mlat, lonb, latb, lon, lat, axes, Time, phys_domain )
! 
!   initialize dust update
!

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

!! Arguments
integer, intent(in)                   :: nlon, mlat
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon, lat
integer, intent(in)                   :: axes(4)
type(time_type), intent(in) :: Time
type(domain2d),      intent(inout) :: phys_domain



call error_mesg('dust_update_init','The null version of dust_update_init should not be called.',FATAL)
end subroutine dust_update_init
! *****************************************************************     
! *****************************************************************


subroutine dust_update_end
! 
! Write soil dust accumulations to netCDF file
implicit none


call error_mesg('dust_update_end','The null version of dust_update_end should not be called.',FATAL)
end subroutine dust_update_end

! *****************************************************************     
! *****************************************************************

subroutine dust_source_sink (   is, js, lon, lat, dt, Time, &
                              p_half, p_full, tsfc,  snow, stress,k_pbl, source_mom, t, tdt, r, rdt, rdt_dst   ) 


integer, intent(in)  :: is, js
real,    intent(in)  :: dt
type(time_type), intent(in)             :: Time
real, intent(in),    dimension(:,:)     :: lon,lat
real, intent(in),    dimension(:,:,:)   :: p_half, p_full
real, intent(in),    dimension(:,:)     :: tsfc,snow,stress
real, intent(in),    dimension(:,:), optional   :: source_mom
integer, intent(in),    dimension(:,:)  :: k_pbl
real, intent(in),    dimension(:,:,:)   :: t,tdt
real, intent(in),    dimension(:,:,:,:) :: r,rdt
real, intent(out),   dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_dst

call error_mesg('dust_source_sink','The null version of dust_source_sink should not be called.',FATAL)
end subroutine dust_source_sink

! *****************************************************************     
! *****************************************************************

subroutine dust_source_init ( nlon, mlat, lonb, latb, lon, lat, axes, Time, phys_domain )
!
! This subroutine is called by hs_forcing_init 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 
use tracer_manager_mod, only: query_method, get_tracer_index,  &
                                 get_number_tracers, get_tracer_names
use horiz_interp_mod,  only: horiz_interp, horiz_interp_init, horiz_interp_new,  &
                              horiz_interp_end,  horiz_interp_del, horiz_interp_type

integer, intent(in)                   :: nlon, mlat
real,    intent(in),  dimension(:,:)  :: lonb, latb 
real,    intent(in),  dimension(:,:)  :: lon, lat 
integer, intent(in)                   :: axes(4)
type(time_type), intent(in) :: Time
type(domain2d),      intent(inout) :: phys_domain

call error_mesg('dust_source_init','The null version of dust_source_init should not be called.',FATAL)
end subroutine dust_source_init

! *****************************************************************     
! *****************************************************************

 subroutine dust_source_end(  days  )
! 
! Write soil dust accumulations to netCDF file
!

implicit none

integer,             intent(in) :: days


call error_mesg('dust_source_end','The null version of dust_source_end should not be called.',FATAL)
end subroutine dust_source_end   


! *****************************************************************     
! *****************************************************************


subroutine sedim_driver(is,js,Time,p_half,p_full,t,tdt, &
     r,rdt,tg,dtime,kd,rdt_sedim,lifting_dust,checkcons)
     !
! Calls init_sedim and sedim to perform aerosol sedimentation and output 
! tendencies for each aerosol tracers
 
integer, intent(in)  :: is, js
type(time_type), intent(in)             :: Time
real*8, intent(in), dimension(:,:,:) :: p_half
real*8, intent(in), dimension(:,:,:) :: p_full
real*8, intent(in), dimension(:,:,:) :: t
real*8, intent(in), dimension(:,:,:,:) :: r
real*8, intent(in), dimension(:,:) :: tg      !     Ground Temperature (K)
real*8, intent(in) :: dtime                   !     NDT time step
real*8, intent(in), dimension(:,:,:) :: kd    !     PBL Eddy coefficient for sedimentation (set to 0 for now)
real*8, intent(in), dimension(:,:,:,:) :: rdt !     tracer tendencies (kg/kg/s)
real*8, intent(in), dimension(:,:,:) :: tdt   !     Temperature tendencies (K/s)
logical, intent(in), dimension(size(r,1),size(r,2),size(r,4)) :: lifting_dust         ! flag is true where lifting occurs
logical, intent(in) :: checkcons
!===============================================================
!     Output Arguments
!===============================================================
real*8, intent(out), dimension(:,:,:,:) :: rdt_sedim


call error_mesg('sedim_driver','The null version of sedim_driver should not be called.',FATAL)
end subroutine sedim_driver   


! *****************************************************************     
! *****************************************************************


subroutine sedim_driver_init( nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )
! initialize module

!! Arguments
integer, intent(in)                   :: nlon, mlat, nlevels
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon, lat
integer, intent(in)                   :: axes(4)
type(time_type), intent(in) :: Time



call error_mesg('sedim_driver_init','The null version of sedim_driver_init should not be called.',FATAL)
end subroutine sedim_driver_init   


! *****************************************************************     
! *****************************************************************

subroutine coagul_main(is, js, lon, lat, dt, temp, p_half, p_full, rini, rdt_coag)

integer, intent(in)  :: is, js
real,    intent(in)  :: dt
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:)     :: lat
real, intent(in),    dimension(:,:,:)   :: temp
real, intent(in),    dimension(:,:,:)   :: p_half, p_full
real, intent(in),    dimension(:,:,:,:)   :: rini
real, intent(inout), dimension(size(rini,1),size(rini,2),size(rini,3),size(rini,4)) :: rdt_coag

call error_mesg('coagul_main','The null version of coagul_main should not be called.',FATAL)
end subroutine coagul_main   


! *****************************************************************     
! *****************************************************************


subroutine coagul_init()


call error_mesg('coagul_init','The null version of coagul_init should not be called.',FATAL)

end subroutine coagul_init

! *****************************************************************     
! *****************************************************************


subroutine nltecool(nlayer,player,tlayer,dt)

implicit none

!
! Input and output variables
!
integer     nlayer                    ! no. of atmospheric layers
real*8       player(nlayer)            ! input pressure grid [Pa]
real*8       tlayer(nlayer)            ! input temperatures [K]
real*8       dt(nlayer)                ! output temp rate [K s-1]


call error_mesg('nltecool','The null version of nltecool should not be called.',FATAL)

end subroutine nltecool
! *****************************************************************     
! *****************************************************************

subroutine blkh2o_driver(is, js, Time, p_half, p_full, t, tdt, r, rdt, stress, tg, dt,  &
                      rhouch, rdt, rkh, tdt, checkcons)

integer, intent(in)  :: is, js
type(time_type), intent(in)             :: Time
real*8, intent(in), dimension(:,:,:) :: p_half
real*8, intent(in), dimension(:,:,:) :: p_full
real*8, intent(in), dimension(:,:,:) :: t
real*8, intent(in), dimension(:,:,:,:) :: r
real, intent(in),    dimension(:,:)     :: rhouch,stress
real*8, intent(in), dimension(:,:) :: tg      !     Ground Temperature (K)
real*8, intent(in) :: dt                      !     NDT time step
real*8, intent(in), dimension(:,:,:) :: rkh    !     PBL Eddy coefficient for sedimentation (set to 0 for now)
real*8, intent(in), dimension(:,:,:,:) :: rdt !     tracer tendencies (kg/kg/s)
real*8, intent(in), dimension(:,:,:) :: tdt   !     Temperature tendencies (K/s)
logical, intent(in) :: checkcons

call error_mesg('blkh2o_driver','The null version of blkh2o_driver should not be called.',FATAL)

! *****************************************************************     
! *****************************************************************

subroutine blkh2o_driver_init(nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )


!! Arguments
integer, intent(in)                   :: nlon, mlat, nlevels
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon, lat
integer, intent(in)                   :: axes(4)
type(time_type), intent(in) :: Time

call error_mesg('blkh2o_driver_init','The null version of blkh2o_driver_init should not be called.',FATAL)

! *****************************************************************     
! *****************************************************************

subroutine co2micro_driver(is,js,Time,p_half,p_full,t,tdt,r,rdt, &
                    stress,tsurf,dt,rhouch,rdt,precip,dmass,rkh, &
                    tdt,checkcons)

integer, intent(in)  :: is, js
type(time_type), intent(in)             :: Time
real*8, intent(in), dimension(:,:,:) :: p_half
real*8, intent(in), dimension(:,:,:) :: p_full
real*8, intent(in), dimension(:,:,:) :: tsurf
real*8, intent(in), dimension(:,:,:,:) :: r
real, intent(in),    dimension(:,:)     :: rhouch, stress
real*8, intent(in), dimension(:,:) :: tg      !     Ground Temperature (K)
real*8, intent(in) :: dt                      !     NDT time step
real*8, intent(in), dimension(:,:,:) :: rkh, dmass
real*8, intent(in), dimension(:,:,:,:) :: rdt
real*8, intent(in), dimension(:,:,:) :: tdt 
logical, intent(in) :: checkcons

call error_mesg('co2micro_driver','The null version of co2micro_driver should not be called.',FATAL)

! *****************************************************************     
! *****************************************************************

subroutine co2micro_driver_init(nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )


!! Arguments
integer, intent(in)                   :: nlon, mlat, nlevels
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon, lat
integer, intent(in)                   :: axes(4)
type(time_type), intent(in) :: Time

call error_mesg('co2micro_driver_init','The null version of co2micro_driver_init should not be called.',FATAL)

! *****************************************************************     
! *****************************************************************
end module null_physics_mod
