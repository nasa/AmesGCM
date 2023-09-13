module radiation_driver_mod
!=======================================================================
!  main radiation driver
!=======================================================================
use constants_mod, only: KAPPA, CP_AIR, GRAV, STEFAN, PI, RADIAN,      &
                        diffac, seconds_per_day,  RAD_TO_DEG, &
                        RDGAS

use fms_mod, only: error_mesg, FATAL,                      &
                   open_namelist_file, check_nml_error,                &
                   mpp_pe, mpp_root_pe, close_file,                    &
                   write_version_number, stdlog,                       &
                   uppercase, read_data, write_data, field_size, field_exist

use fms2_io_mod,            only:  file_exists
use time_manager_mod, only: time_type, get_time
use field_manager_mod, only: MODEL_ATMOS, parse, find_field_index
use tracer_manager_mod, only: query_method, get_tracer_index,  &
    get_number_tracers, get_tracer_names
use diag_manager_mod, only: register_diag_field, send_data, diag_axis_init

!rjw    use radiation_util_mod, only : lw_scattering

use astronomy_mod,  only  :   astronomy_init,  mars_calender, sol88, zenith, &
                              solar_constant


use aerosol_mod,   only:  aerosol_init, aerosol_optics_init,              &
                          dust_distrib_fix, dust_distrib_bin,             &
                          cld_distrib_fix,                                &
                          nice_bins, reff_ice, ice_bin_indx

use ames_rt_interface, only: ames_rt
use rtmod_mgcm, only: ames_radsetup
use initracer_mod, only: ndust_mass, dust_mass_indx, tau_frac

implicit none
private

!------------------- interfaces ---------------------------------------

public :: radiation_driver,                                  &
          radiation_driver_init, radiation_driver_end,       &
          use_ames_lw_rad,  use_ames_sw_rad

!-------------------- namelist -----------------------------------------

logical :: use_forget_swheat= .false.           ! do forget simple visible heating
logical :: use_dust_lwheat= .false.              ! do simple dust IR heating
integer :: rad_calc_intv = 1848                 ! radiative time step

logical :: do_diurnal_avg_rad = .false.         ! do average insolation
logical :: do_rad_diagnostic_calc = .false.     ! do diagnostic calculation with specified temperature and dust

logical :: use_ames_lw_rad = .true.            ! use ames IR radiation
logical :: use_ames_sw_rad = .true.            ! use ames visible radiation

logical :: use_newton_damping = .false.         ! use newton damping for heating

integer :: fixed_calender_date = -1             ! use fixed calender date

integer :: nchan_swx = 1                        ! channel for radiation

real    :: rampup_rad= -1.0                     ! radiation ramping to avoid shocks at initialization [sols]


namelist /radiation_driver_nml/ use_forget_swheat, use_dust_lwheat,        &
                                rad_calc_intv,                             &
                                do_diurnal_avg_rad,                        &
                                use_ames_lw_rad,  use_ames_sw_rad,         &
                                use_newton_damping, fixed_calender_date,   &
                                do_rad_diagnostic_calc,                    &
                                nchan_swx,  rampup_rad

!---------- the following are specific to the GFDL SW and LW radiation code------

integer, parameter :: nchan_lw = 3
integer, parameter :: nchan_sw = 1
integer, parameter :: nradbands= nchan_sw + nchan_lw + 1

! ------------------------------------



real, dimension(:,:,:),   allocatable  ::  heatrate
real, dimension(:,:  ),   allocatable  ::  sfc_sw_flux, sfc_lw_flux
real, dimension(:,:),     allocatable  ::  tstrat
real, dimension(:,:,:,:), allocatable  ::  taudust_reff

real, parameter :: gcp= GRAV / CP_AIR
real, parameter :: qex_ref= 2.5
real, parameter :: co2heat0= 1.3/88775.0
real, parameter :: p0nonlte= 7.5e-3

real, parameter   :: missing_value = -1.e10
character(len=16), parameter :: model='radiation'



logical :: module_is_initialized = .false.
integer :: id_areo
integer :: id_insol, id_swheat, id_lwheat, id_lwdust,                  &
           id_ir_flx, id_opac, id_tdt_rad, id_solar_flx, id_vis_od,    &
           id_vis_od_dust, id_alb_d, id_opacd

integer :: id_lwheat1,id_lwheat2,id_lwheat3, &
           id_lwheat4,id_lwheat5,id_lwheat6, &
           id_lwheat7,id_lwheat8,id_lw15HR

integer :: id_ir_flx_d, id_solar_flx_d, id_swheat_d, id_lwheat_d, id_opac_d
integer :: id_irupflx_d, id_irdnflx_d, id_swupflx_d, id_swdnflx_d
integer :: id_swnetflx_d, id_irnetflx_d
integer :: id_irupflx, id_irdnflx, id_swupflx, id_swdnflx
integer :: id_irupflx_top, id_irdnflx_top, id_swupflx_top, id_swdnflx_top
integer :: id_irupflx_sfc, id_irdnflx_sfc, id_swupflx_sfc, id_swdnflx_sfc
integer :: id_swnetflx, id_irnetflx
integer :: id_taudust_VIS,id_taudust_IR
integer :: id_taucloud_VIS,id_taucloud_IR
integer :: id_trad7, id_trad23, id_trad32
integer :: id_dustref, id_cldref, id_dso
integer, dimension(:), allocatable :: id_taudust_reff_VIS, id_taudust_reff_IR



character(len=128) :: version='$Id: radiation_driver.F90,v 1.1.2.1.2.1 2011/11/22 21:51:30 rjw Exp $'
character(len=128) :: tagname='$Name: mars_feb2012_rjw $'

logical ::  mcpu0

logical :: first_rad


contains

!=====================================================================
subroutine radiation_driver ( is, js, lon, lat, dt, Time,                 &
                               p_half, p_full, z_half, tsfc, albedo,      &
                               sfc_emiss, t, r, tdt, rdt,                 &
                               swfsfc, lwfsfc, cosz, tdtlw,tdt_rad,       &
                               taudust, taucloud, taudust_mom, taudust_fix, pref)
!=======================================================
!  main radiation driver
!=======================================================
integer, intent(in)  :: is, js
real,    intent(in)  :: dt
type(time_type), intent(in)             :: Time
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:)     :: lat
real, intent(in),    dimension(:,:,:)   :: p_half, p_full
real, intent(in),    dimension(:,:,:)   :: z_half
real, intent(in),    dimension(:,:)     :: tsfc
real, intent(in),    dimension(:,:)     :: albedo
real, intent(in),    dimension(:,:)     :: sfc_emiss
real, intent(in),    dimension(:,:,:)   :: t
real, intent(in),    dimension(:,:,:,:) :: r
real, intent(in),    dimension(:,:,:,:) :: rdt
real, intent(in),    dimension(:,:,:)   :: tdt
real, intent(inout), dimension(:,:,:)   :: tdtlw
real, intent(out),   dimension(:,:)     :: swfsfc, lwfsfc
real, intent(out),   dimension(:,:)     :: cosz
real, intent(out),    dimension(:,:,:)   :: tdt_rad
real, intent(out),    dimension(:,:,:)   :: taudust, taucloud, taudust_mom, taudust_fix
real, intent(in)                        :: pref

! Local Variables

real, dimension(size(t,1),size(t,2)) :: hang, frac, coszro
real, dimension(size(t,1),size(t,2)) :: radin, trans, outflx
real, dimension(size(t,1),size(t,2)) :: flx_sfc, flx_sfc_dust
real, dimension(size(t,1),size(t,2)) :: sfc_ir_flx

real, dimension(size(t,1),size(t,2),size(t,3)) :: hsw, heatra, fluxout
real, dimension(size(t,1),size(t,2),size(t,3),8) :: lw_heating_band
real, dimension(size(t,1),size(t,2),size(t,3)) :: lw_15umHR
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: irupflx_d,irdnflx_d,swupflx_d,swdnflx_d
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: swnetflx_d,irnetflx_d
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: irupflx,irdnflx,swupflx,swdnflx
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: swnetflx,irnetflx

real, dimension(size(t,1),size(t,2),size(t,3)) :: hr_dust
real, dimension(size(t,1),size(t,2),size(t,3)) :: delp, tau, opac, delz
real, dimension(size(t,1),size(t,2),size(t,3)) :: dso
real, dimension(size(t,1),size(t,2))           :: vis_od, vis_od_dust

real, dimension(size(t,1),size(t,2),size(p_half,3)) :: tcol
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rnew
!real, dimension(size(r,1),size(r,2),size(r,4),2) :: taudust_reff


!         Here,  nradbands ought to be specified at runtime and
!    the optical properties allocated dynamically
real, dimension(size(t,1),size(t,2),size(t,3),nradbands) :: tau_spec
real, dimension(size(t,1),size(t,2),size(t,3),nradbands) :: sscat_spec
real, dimension(size(t,1),size(t,2),size(t,3),nradbands) :: gfac_spec



real, dimension(size(t,1),size(t,2),size(t,3))   :: tau_ames, darray
real, dimension(size(t,1),size(t,2),size(t,3))   :: pf_mb
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: ph_mb

real, dimension(size(t,1),size(t,2),size(t,3))   :: dustref, cldref
real, dimension(size(t,1),size(t,2),size(t,3))   :: dustref_bin, dustref_fix, cldice
real, dimension(size(t,1),size(t,2),size(t,3))   :: dustref_mom
real, dimension(size(t,1),size(t,2),size(t,3))   :: tdiag,tnew
real, dimension(size(t,1),size(t,2),size(p_half,3)) :: tcdiag
real, dimension(size(t,1),size(t,2))  :: tsdiag
real, dimension(size(t,1),size(t,2))  :: tstrat_dt, tsat

!   currently carrying 3 diagnostic brightness temperatures
!        see dimension statement in sfc_rad_flx_diag
real, dimension(size(t,1),size(t,2),3)   :: tbands


logical :: used
real    :: flux, sink, value
character(len=128) :: scheme, params

real :: fjd, r_orbit, declin, areolat, slag, dhr, areoout
real :: rsolar, fact, ramp
real :: secs, sols

integer :: days, seconds
integer :: nchan, lwchan, n
integer :: ie, je, id, jd, kd, ntrace, nice
integer :: i, j, k, kb, idjd, nt, ndx

mcpu0 = (mpp_pe() == mpp_root_pe())


id= size(t,1); jd= size(t,2); kd= size(t,3)

ntrace= size(r,4)

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
    call mars_calender( fixed_calender_date, fjd, r_orbit, declin, areoout )
else
    call mars_calender( days, fjd, r_orbit, declin, areoout )
endif

rsolar= solar_constant / r_orbit**2

areoout = areoout * RADIAN
if( areoout < 0.0 )  areoout = areoout + 360.0
areolat = modulo(areoout,360.)

!           Next get zenith angle
call sol88( declin, lat, hang, cosz, frac )

coszro(:,:)= cosz(:,:)

if( .not.do_diurnal_avg_rad  ) then
    slag= 0.0
    dhr= 0.5 * (rad_calc_intv/seconds_per_day)*pi
    call zenith( fjd, declin, slag, lon, lat, hang, dhr, cosz, frac )
endif


!           Now calculate solar heating

! Update tracers:
rnew=r+rdt*dt
! Update temperature:
tnew=t+tdt*dt

if( first_rad .or.  mod( seconds, rad_calc_intv) == 0 ) then   !--------------------------
    first_rad = .false.

    radin(:,:)= rsolar*frac(:,:)*cosz(:,:)

    do k=1,kd
        delp(:,:,k)= p_half(:,:,k+1)-p_half(:,:,k)
        delz(:,:,k)= z_half(:,:,k+1)-z_half(:,:,k)
    enddo



    tcol(:,:,1:kd)= tnew(:,:,:)
    tcol(:,:,kd+1)= tsfc(:,:)


    call dust_distrib_bin( is, js, areolat, lon, lat, pref,                 &
                           p_full, p_half, delp, rnew(:,:,:,1:ntrace),  &
                           tau_spec,  sscat_spec, gfac_spec,         &
                           dustref_bin)

    call dust_distrib_fix( is, js, areolat, lon, lat, pref,                 &
                           p_full, p_half, delp, rnew(:,:,:,1:ntrace),  &
                           tau_spec,  sscat_spec, gfac_spec,         &
                           dustref_fix)
    dustref(:,:,:) = 0.
    cldref(:,:,:) = 0.
    nice= ice_bin_indx(1)
    cldice(:,:,:)= r(:,:,:,nice)

    call cld_distrib_fix( is, js, areolat, lon, lat, pref, p_full, p_half, delp, &
                               rnew(:,:,:,1:ntrace), cldice )

    pf_mb= p_full * 0.01
    ph_mb= p_half * 0.01

    taudust_reff(is:ie,js:je,:,:) = 0.d0

    call  ames_rt(is,js,id,jd,kd,ntrace,ph_mb,pf_mb,    &
                   tnew,tsfc,rnew,trans,flx_sfc,        &
                   albedo, sfc_emiss, cosz,             &
                   dustref, dustref_bin, dustref_fix,   &
                   cldref, cldice,                      &
                   use_ames_sw_rad, use_ames_lw_rad,    &
                   r_orbit,                             &
                   heatra,hsw,outflx,                   &
                   rsolar,                              &
                   irupflx,irdnflx,swupflx,swdnflx,     &
                   swnetflx,irnetflx,                   &
                   taudust,taucloud,taudust_mom,        &
                   lw_heating_band,lw_15umHR,.false.,   &
                   tstrat(is:ie,js:je),                 &
                   tstrat_dt,                           &
                   taudust_reff(is:ie,js:je,:,:),       &
                   taudust_fix,                         &
                   tbands                               )

!taudust_reff(is:ie,js:je,:,:) = 0.01

    tstrat(is:ie,js:je)=tstrat(is:ie,js:je)+tstrat_dt(:,:)*real(rad_calc_intv)
    tsat(:,:) = 3182.48 / ( 23.3494 - log( ph_mb(:,:,1) ) )
    WHERE (tstrat(is:ie,js:je) .lt. tsat(:,:)) tstrat(is:ie,js:je) = tsat(:,:)

    do k= 1, kd 
        hsw(:,:,k)=  frac(:,:)*hsw(:,:,k)
    enddo 
    trans(:,:)=   frac(:,:)*trans(:,:)
    sfc_ir_flx= flx_sfc


!     Gradually increase the strength of solar forcing
    if ( rampup_rad > 0.0 ) then
        if( sols .lt. (rampup_rad) ) then
            ramp= 0.5*( 1.0 - cos( (sols)*pi/rampup_rad )  )
            hsw=   ramp*hsw
            trans= ramp*trans
        end if
    endif

!      Save reference dust and ice cloud opacities
!     These are in units of  opacity/Pa
!
!           write out the normalized dust field used by the radiation code

    if (id_opac > 0) then
        opac(:,:,:)= dustref(:,:,:) / delp(:,:,:)
        used = send_data ( id_opac, opac,  time, is, js )
    endif

    if (id_opacd > 0) then
        opac(:,:,:)= dustref_fix(:,:,:) / delp(:,:,:)
        used = send_data ( id_opacd, opac,  time, is, js )
    endif


    vis_od= 0.0
    vis_od_dust= 0.0
    do k= 1, kd
        vis_od(:,:)= vis_od(:,:) + dustref_bin(:,:,k)
        vis_od_dust(:,:)= vis_od_dust(:,:) + dustref_fix(:,:,k)
    enddo
    if (id_vis_od > 0)      used = send_data ( id_vis_od,      vis_od,       time, is, js )
    if (id_vis_od_dust > 0) used = send_data ( id_vis_od_dust, vis_od_dust,  time, is, js )


    !           write out the non-normalized dust field used by the radiation code
    opac(:,:,:)= dustref(:,:,:) ! / delp(:,:,:)
    if (id_dustref > 0) used = send_data ( id_dustref, opac,  time, is, js )

    if (id_dso > 0)then
        dso = (dustref*rdgas*t)/delz/p_full
        used = send_data ( id_dso, dso, time, is, js)
    endif

    !           write out the normalized cloud field
    opac(:,:,:)= cldref(:,:,:) ! / delp(:,:,:)
    if (id_cldref > 0) used = send_data ( id_cldref, opac,  time, is, js )


    if( id_areo > 0 ) used= send_data( id_areo, areoout , Time )

    !            Combine lwave and shortwave fluxes
    if (id_swheat > 0)  used = send_data ( id_swheat, hsw,    Time, is, js)
    if (id_lwheat > 0)  used = send_data ( id_lwheat, heatra, Time,  is, js)
    if (id_taudust_VIS > 0)  used = send_data ( id_taudust_VIS,taudust(:,:,1),   Time, is, js)
    if (id_taudust_IR > 0)  used = send_data ( id_taudust_IR,taudust(:,:,2),   Time, is, js)
    if (id_taucloud_VIS > 0)  used = send_data ( id_taucloud_VIS,taucloud(:,:,1),   Time, is, js)
    if (id_taucloud_IR > 0)  used = send_data ( id_taucloud_IR,taucloud(:,:,2),   Time, is, js)
    !           write out moment column dust field by effective radius
    do nt=1,ndust_mass
        ndx= dust_mass_indx(nt)
        if (id_taudust_reff_VIS(nt) > 0) used = send_data ( id_taudust_reff_VIS(nt), taudust_reff(is:ie,js:je,ndx,1),  time, is, js )
        if (id_taudust_reff_IR(nt) > 0) used = send_data ( id_taudust_reff_IR(nt), taudust_reff(is:ie,js:je,ndx,2),  time, is, js )
    enddo

    if (id_lwheat1 > 0)  used = &
          send_data ( id_lwheat1, lw_heating_band(:,:,:,1),    Time, is, js)
    if (id_lwheat2 > 0)  used = &
          send_data ( id_lwheat2, lw_heating_band(:,:,:,2),    Time, is, js)
    if (id_lwheat3 > 0)  used = &
          send_data ( id_lwheat3, lw_heating_band(:,:,:,3),    Time, is, js)
    if (id_lwheat4 > 0)  used = &
          send_data ( id_lwheat4, lw_heating_band(:,:,:,4),    Time, is, js)
    if (id_lwheat5 > 0)  used = &
          send_data ( id_lwheat5, lw_heating_band(:,:,:,5),    Time, is, js)
    if (id_lwheat6 > 0)  used = &
          send_data ( id_lwheat6, lw_heating_band(:,:,:,6),    Time, is, js)
    if (id_lwheat7 > 0)  used = &
          send_data ( id_lwheat7, lw_heating_band(:,:,:,7),    Time, is, js)
    if (id_lwheat8 > 0)  used = &
          send_data ( id_lwheat8, lw_heating_band(:,:,:,8),    Time, is, js)
    if (id_lw15HR > 0)  used = &
          send_data ( id_lw15HR, lw_15umHR(:,:,:),    Time, is, js)

!!! write out diagnostic fluxes
! ,,,
    if (id_irupflx > 0)  then
        fluxout(:,:,:) = irupflx(:,:,2:kd+1)
        used = send_data ( id_irupflx, fluxout,    Time, is, js)
    endif
    if (id_irdnflx > 0)  then
        fluxout(:,:,:) = irdnflx(:,:,2:kd+1)
        used = send_data ( id_irdnflx, fluxout,    Time, is, js)
    endif
    if (id_swupflx > 0)  then
        fluxout(:,:,:) = swupflx(:,:,2:kd+1)
        used = send_data ( id_swupflx, fluxout,    Time, is, js)
    endif
    if (id_swdnflx > 0)  then
        fluxout(:,:,:) = swdnflx(:,:,2:kd+1)
        used = send_data ( id_swdnflx, fluxout,    Time, is, js)
    endif
    if (id_swnetflx> 0)  then
        fluxout(:,:,:) = swnetflx(:,:,2:kd+1)
        used = send_data ( id_swnetflx,fluxout,   Time, is, js)
    endif
    if (id_irnetflx> 0)  then
        fluxout(:,:,:) = irnetflx(:,:,2:kd+1)
        used = send_data ( id_irnetflx,fluxout,   Time, is, js)
    endif
    if (id_irupflx_top > 0)  then
        used = send_data ( id_irupflx_top, irupflx(:,:,1),    Time, is, js)
    endif
    if (id_irdnflx_top > 0)  then
        used = send_data ( id_irdnflx_top, irdnflx(:,:,1),    Time, is, js)
    endif
    if (id_swupflx_top > 0)  then
        used = send_data ( id_swupflx_top, swupflx(:,:,1),    Time, is, js)
    endif
    if (id_swdnflx_top > 0)  then
        used = send_data ( id_swdnflx_top, swdnflx(:,:,1),    Time, is, js)
    endif
    if (id_irupflx_sfc > 0)  then
        used = send_data ( id_irupflx_top, irupflx(:,:,kd+1),    Time, is, js)
    endif
    if (id_irdnflx_sfc > 0)  then
        used = send_data ( id_irdnflx_top, irdnflx(:,:,kd+1),    Time, is, js)
    endif
    if (id_swupflx_sfc > 0)  then
        used = send_data ( id_swupflx_top, swupflx(:,:,kd+1),    Time, is, js)
    endif
    if (id_swdnflx_sfc > 0)  then
        used = send_data ( id_swdnflx_top, swdnflx(:,:,kd+1),    Time, is, js)
    endif

    tdtlw = heatra

    heatra = heatra + hsw

    swfsfc = trans
    lwfsfc = sfc_ir_flx
    tdt_rad = heatra

    !         Update global arrays which can save radiation over multiple time steps

    sfc_sw_flux(is:ie,js:je)=   swfsfc
    sfc_lw_flux(is:ie,js:je)=   lwfsfc
    heatrate   (is:ie,js:je,:)= heatra



    !   ---------------- Calculate diagnostic brightness temperatures -----------------
    if ( id_trad7  > 0 .or.  id_trad23  > 0 .or. id_trad32 > 0 )  then

    !rjw    call sfc_rad_flx_diag( p_full, p_half, dustref, cldref, tcol,  &
    !rjw                              tcol(:,:,kd+1), sfc_emiss, tbands )

        if ( id_trad7   > 0 ) used = send_data ( id_trad7,  tbands(:,:,1),  Time, is, js)
        if ( id_trad23  > 0 ) used = send_data ( id_trad23, tbands(:,:,2),  Time, is, js)
        if ( id_trad32  > 0 ) used = send_data ( id_trad32, tbands(:,:,3),  Time, is, js)

    endif


!                ---------------- Diagnostic calculation -----------------
    if( do_rad_diagnostic_calc ) then

        !            This is an example of calling a new radiation calculation
        !        In this case, the opacity is replaced by the dust opacity only,
        !       thus allowing an assessment of radiatively active water ice clouds

        !        tau_spec(:,:,:,1)= dustref(:,:,:)
        if (kd.eq.28) then
            dustref(:,:,1:6) = 0.
            dustref(:,:,7) = .00002
            dustref(:,:,8) = .00008
            dustref(:,:,9) = .00011
            dustref(:,:,10) = .00023
            dustref(:,:,11) = .0004
            dustref(:,:,12) = .00055
            dustref(:,:,13) = .00066
            dustref(:,:,14) = .00079
            dustref(:,:,15) = .00106
            dustref(:,:,16) = .00116
            dustref(:,:,17) = .0011
            dustref(:,:,18) = .001105
            dustref(:,:,19) = .00111
            dustref(:,:,20) = .001115
            dustref(:,:,21) = .001118
            dustref(:,:,22) = .001105
            dustref(:,:,23) = .001088
            dustref(:,:,24) = .001082
            dustref(:,:,25:28) = .001081
        endif
        ! do no dust
        !         dustref = 0.
        !         dustref(:,:,1:10) = 0.


        !   Set the diagnostic test temperature
        if (kd.eq.28) then
            tsdiag(:,:) = 236.623
            tdiag(:,:,1) = 150.
            tdiag(:,:,2) = 147.
            tdiag(:,:,3) = 142.
            tdiag(:,:,4) = 143.
            tdiag(:,:,5) = 141.
            tdiag(:,:,6) = 140.
            tdiag(:,:,7) = 140.7
            tdiag(:,:,8) = 146.
            tdiag(:,:,9) = 152.
            tdiag(:,:,10) = 160.5
            tdiag(:,:,11) = 169.5
            tdiag(:,:,12) = 175.
            tdiag(:,:,13) = 181.
            tdiag(:,:,14) = 190.
            tdiag(:,:,15) = 196.
            tdiag(:,:,16) = 202.
            tdiag(:,:,17) = 207.
            tdiag(:,:,18) = 213.5
            tdiag(:,:,19) = 217.5
            tdiag(:,:,20) = 220.5
            tdiag(:,:,21) = 222.5
            tdiag(:,:,22) = 223.5
            tdiag(:,:,23) = 223.
            tdiag(:,:,24) = 223.3
            tdiag(:,:,25) = 224.2
            tdiag(:,:,26) = 224.9
            tdiag(:,:,27) = 225.6
            tdiag(:,:,28) = 226.1
            tcdiag(:,:,1:28) = tdiag(:,:,:)
            tcdiag(:,:,29) = tsdiag(:,:)
        endif
    !       tdiag(:,:,:) = 180.
    !       tcdiag(:,:,:)= 180.
    !       tsdiag(:,:)= 180.

        cldref(:,:,:)= 0.0

        pf_mb= p_full * 0.01
        ph_mb= p_half * 0.01

        call  ames_rt(is,js,id,jd,kd,ntrace,ph_mb,pf_mb,  &
                    tdiag,tsdiag,rnew,trans,flx_sfc,&
                    albedo, sfc_emiss, cosz, dustref, &
                       dustref_bin, dustref_fix, cldref, cldice, &
                       use_ames_sw_rad, use_ames_lw_rad, r_orbit,&
                       heatra,hsw,outflx,rsolar, &
                       irupflx_d,irdnflx_d, &
                       swupflx_d,swdnflx_d,swnetflx_d,&
                       irnetflx_d,taudust,taucloud,taudust_mom, &
                       lw_heating_band,lw_15umHR,.true., &
                       tstrat(is:ie,js:je), tstrat_dt,taudust_reff, &
                   taudust_fix,                         &
                       tbands     ) 


    !           Save the diagnostic heating rates and surface fluxes
        if (id_swheat_d > 0)  used = send_data ( id_swheat_d, hsw,    Time, is, js)
        if (id_lwheat_d > 0)  used = send_data ( id_lwheat_d, heatra, Time,  is, js)

        if (id_ir_flx_d     > 0)  used = send_data ( id_ir_flx_d,    flx_sfc, Time, is, js)
        if (id_solar_flx_d  > 0)  used = send_data ( id_solar_flx_d, trans, Time, is, js)

    !!! write out diagnostic fluxes
    ! ,,,
        if (id_irupflx_d > 0)  then
            fluxout(:,:,:) = irupflx_d(:,:,2:kd+1)
            used = send_data ( id_irupflx_d, fluxout,    Time, is, js)
        endif
        if (id_irdnflx_d > 0)  then
            fluxout(:,:,:) = irdnflx_d(:,:,2:kd+1)
            used = send_data ( id_irdnflx_d, fluxout,    Time, is, js)
        endif
        if (id_swupflx_d > 0)  then
            fluxout(:,:,:) = swupflx_d(:,:,2:kd+1)
            used = send_data ( id_swupflx_d, fluxout,    Time, is, js)
        endif
        if (id_swdnflx_d > 0)  then
            fluxout(:,:,:) = swdnflx_d(:,:,2:kd+1)
            used = send_data ( id_swdnflx_d, fluxout,    Time, is, js)
        endif
        if (id_swnetflx_d> 0)  then
            fluxout(:,:,:) = swnetflx_d(:,:,2:kd+1)
            used = send_data ( id_swnetflx_d,fluxout,   Time, is, js)
        endif
        if (id_irnetflx_d> 0)  then
            fluxout(:,:,:) = irnetflx_d(:,:,2:kd+1)
            used = send_data ( id_irnetflx_d,fluxout,   Time, is, js)
        endif

        if (id_taudust_VIS > 0)  used = send_data ( id_taudust_VIS,taudust(:,:,1),   Time, is, js)
        if (id_taudust_IR > 0)  used = send_data ( id_taudust_IR,taudust(:,:,2),   Time, is, js)
        if (id_taucloud_VIS > 0)  used = send_data ( id_taucloud_VIS,taucloud(:,:,1),   Time, is, js)
        if (id_taucloud_IR > 0)  used = send_data ( id_taucloud_IR,taucloud(:,:,2),   Time, is, js)


        !           Write out the visible opacity field used for diagnostic purposes
        opac(:,:,:)= dustref(:,:,:) / delp(:,:,:)
        if (id_opac_d > 0) used = send_data ( id_opac_d, darray,  Time, is, js )

        if (id_alb_d > 0) used = send_data ( id_alb_d, outflx,  Time, is, js )


    endif       ! ------------------- End of diagnostic calcs ---------



else        !    use previously computed heating rates

    tdt_rad(:,:,:)= heatrate(is:ie,js:je,:)

    swfsfc(:,:) = sfc_sw_flux(is:ie,js:je)
    lwfsfc(:,:) = sfc_lw_flux(is:ie,js:je)

endif


if (id_ir_flx     > 0)  used = send_data ( id_ir_flx,    lwfsfc, Time, is, js)
if (id_solar_flx  > 0)  used = send_data ( id_solar_flx, swfsfc, Time, is, js)
if (id_tdt_rad    > 0)  used = send_data ( id_tdt_rad,   tdt_rad,    Time, is, js)


end subroutine radiation_driver


!=======================================================================
!=======================================================================


subroutine radiation_driver_init ( nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )
!=======================================================================
! initialize radiation driver
!=======================================================================


integer, intent(in)                    :: nlon, mlat, nlevels
real,    intent(in),  dimension(:,:)   :: lonb, latb
real,    intent(in),  dimension(:,:)   :: lon , lat
integer, intent(in) :: axes(4)
type(time_type), intent(in) :: Time

integer :: unit, io, ierr
integer :: i, j, ie, je, id, jd, is, js, nt, ndx
integer :: null_axis_id
integer ::  ntrace, ntprog, ntfam, ntdiag

character (len=128) :: filename, fieldname, tracer_name, tname


id= size(lonb,1)-1; jd= size(latb,2)-1

is= 1;   js= 1;
ie= is + id - 1
je= js + jd - 1

mcpu0 = (mpp_pe() == mpp_root_pe())

first_rad=.true.

!     ----- read namelist -----

if (file_exists('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
        read  (unit, nml=radiation_driver_nml, iostat=io, end=20)
        ierr = check_nml_error (io, 'radiation_driver_nml')
    enddo
20    call close_file (unit)
endif

!     ----- write version info and namelist to log file -----

call write_version_number (version,tagname)

if (mpp_pe() == mpp_root_pe())    write (stdlog(),nml=radiation_driver_nml)


allocate (    heatrate(is:ie,js:je,nlevels)  )
allocate ( sfc_sw_flux(is:ie,js:je)  )
allocate ( sfc_lw_flux(is:ie,js:je)  )
allocate ( tstrat(is:ie,js:je)  )
! a) First get the total number of tracer : ntrace
call get_number_tracers (MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfam )
allocate ( taudust_reff(is:ie,js:je,ntrace,2) )

!       Get heating rates and surface fluxes from radiation restart file
filename= 'INPUT/radiation.res.nc'
if (file_exists( trim( filename ) ) ) then
    if(mcpu0) print *,'Reading  restart file:   ',  trim( filename )
    call read_data( trim( filename ), 'heatrate'   , heatrate)
    call read_data( trim( filename ), 'sfc_sw_flux', sfc_sw_flux)
    call read_data( trim( filename ), 'sfc_lw_flux', sfc_lw_flux)
    if (field_exist( trim( filename ), 'tstrat')) then
        call read_data( trim( filename ), 'tstrat', tstrat)
    else
        tstrat = 170.
    endif
else
    heatrate= 0.0
    sfc_sw_flux= 0.0
    sfc_lw_flux= 0.0
    tstrat= 170.
endif


if( use_ames_lw_rad  .OR. use_ames_sw_rad ) then
    filename='input.nml'
    call ames_radsetup(nlay=nlevels,mcpu0ARG=mcpu0,nml_fileARG=filename)
    if(mcpu0) print *, '  Have returned from new Ames radsetup '
else
    if(mcpu0) print *, 'Initializing heating  ',  is, ie, js, je
endif




!     ----- register diagnostic fields -----
null_axis_id = diag_axis_init('scalar_axis', (/0./), 'none', 'X', 'none')

id_areo = register_diag_field ( model, 'areo', (/null_axis_id/),          &
                                       Time, 'Areocentric Longitude', 'deg', &
                                      missing_value=missing_value )

!id_insol = register_diag_field ( model, 'insol', axes(1:2),          &
!                                       Time, 'Solar insolation', 'W/m2', &
!                                      missing_value=missing_value )

id_ir_flx = register_diag_field ( model, 'sfcirflx', axes(1:2),      &
                                 Time, 'surface net downward IR flux', 'W/m2', &
                                       missing_value=missing_value )

id_solar_flx = register_diag_field ( model, 'sfcswflx', axes(1:2),   &
                                 Time, 'surface net downward VIS flux', 'W/m2', &
                                       missing_value=missing_value )

id_swheat = register_diag_field ( model, 'swheat', axes(1:3),        &
                                       Time, 'visible radiation heating rate', 'K/s', &
                                      missing_value=missing_value )

id_lwheat = register_diag_field ( model, 'lwheat', axes(1:3),        &
                                       Time, 'IR radiation heating rate', 'K/s', &
                                      missing_value=missing_value )

id_lwheat1 = register_diag_field ( model, 'lwheat1', axes(1:3),        &
                                       Time, 'IR heating band 1', 'K/s', &
                                      missing_value=missing_value )

id_lwheat2 = register_diag_field ( model, 'lwheat2', axes(1:3),        &
                                       Time, 'IR heating band 2', 'K/s', &
                                      missing_value=missing_value )

id_lwheat3 = register_diag_field ( model, 'lwheat3', axes(1:3),        &
                                       Time, 'IR heating band 3', 'K/s', &
                                      missing_value=missing_value )

id_lwheat4 = register_diag_field ( model, 'lwheat4', axes(1:3),        &
                                       Time, 'IR heating band 4', 'K/s', &
                                      missing_value=missing_value )

id_lwheat5 = register_diag_field ( model, 'lwheat5', axes(1:3),        &
                                       Time, 'IR heating band 5', 'K/s', &
                                      missing_value=missing_value )

id_lwheat6 = register_diag_field ( model, 'lwheat6', axes(1:3),        &
                                       Time, 'IR heating band 6', 'K/s', &
                                      missing_value=missing_value )

id_lwheat7 = register_diag_field ( model, 'lwheat7', axes(1:3),        &
                                       Time, 'IR heating band 7', 'K/s', &
                                      missing_value=missing_value )

id_lwheat8 = register_diag_field ( model, 'lwheat8', axes(1:3),        &
                                       Time, 'IR heating band 8', 'K/s', &
                                      missing_value=missing_value )

id_lw15HR = register_diag_field ( model, 'lw15HR', axes(1:3),        &
                                       Time, 'NLTE 15um heating', 'K/s', &
                                      missing_value=missing_value )

!id_lwdust = register_diag_field ( model, 'lwdust', axes(1:3),       &
!                                  Time, 'IR dust heating', 'K/s', &
!                                      missing_value=missing_value )

id_opac = register_diag_field ( model, 'opac', axes(1:3),            &
                               Time, 'visible dust opacity', 'op/Pa', &
                                      missing_value=missing_value )

id_opacd = register_diag_field ( model, 'opac_fix', axes(1:3),            &
                               Time, 'visible fixed dust opacity', 'op/Pa', &
                                      missing_value=missing_value )

id_opac_d = register_diag_field ( model, 'opac_d', axes(1:3),        &
                    Time, 'visible diagnostic dust opacity', 'op/Pa', &
                                      missing_value=missing_value )

id_vis_od = register_diag_field ( model, 'vis_od', axes(1:2),      &
                               Time, 'visible bin dust opacity', 'total', &
                                      missing_value=missing_value )

id_dustref = register_diag_field ( model, 'dustref', axes(1:3),      &
                               Time, 'visible dust opacity', 'op/level', &
                                      missing_value=missing_value )

id_cldref = register_diag_field ( model, 'cldref', axes(1:3),      &
                               Time, 'visible water ice cloud opacity', 'op/level', &
                                      missing_value=missing_value )

id_dso = register_diag_field ( model, 'dso', axes(1:3),      &
                               Time, 'visible density scaled dust opacity', 'm2/kg', &
                                      missing_value=missing_value )

id_vis_od_dust = register_diag_field ( model, 'vis_od_dust', axes(1:2), &
                               Time, 'visible fixed dust column opacity', 'total', &
                                      missing_value=missing_value )

id_tdt_rad = register_diag_field ( model, 'tdt_rad', axes(1:3),      &
                    Time, 'radiative temperature tendency', 'K/s', &
                                      missing_value=missing_value )

id_ir_flx_d = register_diag_field ( model, 'sfcirflx_d', axes(1:2),  &
                                 Time, 'diagnostic net surface downward IR flux', 'W/m2', &
                                       missing_value=missing_value )

id_solar_flx_d = register_diag_field ( model, 'sfcswflx_d', axes(1:2), &
                                 Time, 'diagnostic net surface downward VIS flux', 'W/m2',  &
                                       missing_value=missing_value )

id_alb_d = register_diag_field ( model, 'alb_d', axes(1:2), &
                                 Time, 'diagnostic albedo', 'None',  &
                                       missing_value=missing_value )

id_swheat_d = register_diag_field ( model, 'swheat_d', axes(1:3),    &
                                       Time, 'diagnostic VIS heating', 'K/s', &
                                      missing_value=missing_value )

id_lwheat_d = register_diag_field ( model, 'lwheat_d', axes(1:3),    &
                                       Time, 'diagnostic IR heating', 'K/s', &
                                      missing_value=missing_value )

id_irupflx_d = register_diag_field ( model, 'irupflx_d', axes(1:3),&
                                       Time, 'diagnostic upwards IR flux diag', 'W/m2', &
                                      missing_value=missing_value )

id_irdnflx_d = register_diag_field ( model, 'irdnflx_d', axes(1:3),&
                                       Time, 'diagnostic downwards IRR flux diag', 'W/m2', &
                                      missing_value=missing_value )

id_swupflx_d = register_diag_field ( model, 'swupflx_d',axes(1:3),&
                                       Time, 'diagnostic upwards VIS flux diag', 'W/m2', &
                                      missing_value=missing_value )

id_swdnflx_d = register_diag_field ( model, 'swdnflx_d', axes(1:3),&
                                       Time, 'diagnostic downwards VIS flux diag', 'W/m2', &
                                      missing_value=missing_value )

id_swnetflx_d = register_diag_field ( model, 'swnetflx_d', axes(1:3),&
                              Time, 'diagnostic net VIS flux diag', 'W/m2', &
                              missing_value=missing_value )

id_irnetflx_d = register_diag_field ( model, 'irnetflx_d', axes(1:3),&
                              Time, 'diagnostic net IR flux diag', 'W/m2', &
                              missing_value=missing_value )

id_irupflx = register_diag_field ( model, 'irupflx', axes(1:3),&
                                       Time, 'upwards IR flux', 'W/m2', &
                                      missing_value=missing_value )

id_irdnflx = register_diag_field ( model, 'irdnflx', axes(1:3),&
                                       Time, 'downwards IR flux', 'W/m2', &
                                      missing_value=missing_value )

id_swupflx = register_diag_field ( model, 'swupflx',axes(1:3),&
                                       Time, 'upwards VIS flux', 'W/m2', &
                                      missing_value=missing_value )

id_swdnflx = register_diag_field ( model, 'swdnflx', axes(1:3),&
                                       Time, 'downwards VIS flux', 'W/m2', &
                                      missing_value=missing_value )

id_swnetflx = register_diag_field ( model, 'swnetflx', axes(1:3),&
                              Time, 'net VIS flux', 'W/m2', &
                              missing_value=missing_value )

id_irnetflx = register_diag_field ( model, 'irnetflx', axes(1:3),&
                              Time, 'net IR flux', 'W/m2', &
                              missing_value=missing_value )

id_irupflx_top = register_diag_field ( model, 'irupflx_top', axes(1:2),&
                                       Time, 'ToA upwards IR flux', 'W/m2', &
                                      missing_value=missing_value )

id_irdnflx_top = register_diag_field ( model, 'irdnflx_top', axes(1:2),&
                                       Time, 'ToA downwards IR flux', 'W/m2', &
                                      missing_value=missing_value )

id_swupflx_top = register_diag_field ( model, 'swupflx_top',axes(1:2),&
                                       Time, 'ToA upwards VIS flux', 'W/m2', &
                                      missing_value=missing_value )

id_swdnflx_top = register_diag_field ( model, 'swdnflx_top', axes(1:2),&
                                       Time, 'ToA downwards VIS flux', 'W/m2', &
                                      missing_value=missing_value )

id_irupflx_sfc = register_diag_field ( model, 'irupflx_sfc', axes(1:2),&
                                       Time, 'Surface upwards IR flux', 'W/m2', &
                                      missing_value=missing_value )

id_irdnflx_sfc = register_diag_field ( model, 'irdnflx_sfc', axes(1:2),&
                                       Time, 'Surface downwards IR flux', 'W/m2', &
                                      missing_value=missing_value )

id_swupflx_sfc = register_diag_field ( model, 'swupflx_sfc',axes(1:2),&
                                       Time, 'Surface upwards VIS flux', 'W/m2', &
                                      missing_value=missing_value )

id_swdnflx_sfc = register_diag_field ( model, 'swdnflx_sfc', axes(1:2),&
                                       Time, 'Surface downwards VIS flux', 'W/m2', &
                                      missing_value=missing_value )

id_taudust_VIS = register_diag_field ( model, 'taudust_VIS', axes(1:2),&
                              Time, 'Column dust opacity VIS', 'op', &
                              missing_value=missing_value )

id_taudust_IR = register_diag_field ( model, 'taudust_IR', axes(1:2),&
                              Time, 'Column dust opacity IR', 'op', &
                              missing_value=missing_value )

id_taucloud_VIS = register_diag_field ( model, 'taucloud_VIS', axes(1:2),&
                              Time, 'Column cloud opacity VIS', 'op', &
                              missing_value=missing_value )

id_taucloud_IR = register_diag_field ( model, 'taucloud_IR', axes(1:2),&
                              Time, 'Column cloud opacity IR', 'op', &
                              missing_value=missing_value )

id_trad7 = register_diag_field ( model, 'trad7', axes(1:2),          &
                               Time, '7um brightness temp', 'K',   &
                                      missing_value=missing_value )

id_trad23 = register_diag_field ( model, 'trad23', axes(1:2),      &
                               Time, '23um brightness temp', 'K', &
                                      missing_value=missing_value )

id_trad32 = register_diag_field ( model, 'trad32', axes(1:2),      &
                               Time, '32um brightness temp', 'K', &
                                      missing_value=missing_value )

allocate( id_taudust_reff_VIS(ndust_mass) )
allocate( id_taudust_reff_IR(ndust_mass) )


do nt= 1, ndust_mass

    !Get name tracer
    ndx= dust_mass_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)
    !!! dust effective radius for background lifting
    tname= trim(tracer_name) // '_reff_vis'
    id_taudust_reff_VIS(nt) = register_diag_field ( model, trim(tname), axes(1:2), Time, &
                      'Radiation taudust VIS by R_eff', '',  &
                       missing_value=missing_value     )
    tname= trim(tracer_name) // '_reff_ir'
    id_taudust_reff_IR(nt) = register_diag_field ( model, trim(tname), axes(1:2), Time, &
                      'Radiation taudust IR by R_eff', '',  &
                       missing_value=missing_value     )

enddo


module_is_initialized  = .true.

end subroutine radiation_driver_init


!=======================================================================
!=======================================================================


subroutine sfc_rad_flx_diag( pp, pph, dust, cloud, tcol, tsfc, sfc_emiss, &
                                  tbands )
!=======================================================================
!  diagnostic surface flux
!=======================================================================

 use radiation_util_mod, only : planck_func, planck_inv_func, lw_scattering

 real, intent(in), dimension(:,:,:)   :: pp    ! kd
 real, intent(in), dimension(:,:,:)   :: pph   ! kp
 real, intent(in), dimension(:,:,:)   :: tcol  ! kd+1 values
 real, intent(in), dimension(:,:)     :: tsfc      ! surface temperature
 real, intent(in), dimension(:,:)     :: sfc_emiss
 real, intent(in), dimension(:,:,:)   :: dust, cloud

!       OUTPUT:
 real, intent(out), dimension(:,:,:) :: tbands

!       LOCAL:
integer, parameter :: NBANDS = 3
integer     :: id, jd, kd, kp, nn
real        :: cose, xwl

integer ii, jj


real, dimension(NBANDS)  :: xcenter, qex_dust, qex_ice , sscat_dust,  &
                                  sscat_ice, gfac_dust, gfac_ice

real, dimension(size(pp,1),size(pp,2),size(pp,3))    :: tau, sscat, gfac

real, dimension(size(pph,1),size(pph,2),size(pph,3)) :: bh, fuh, fdh
real, dimension(size(pp,1),size(pp,2))               :: bsfc, radiance

real, dimension(size(pp,1),size(pp,2))               :: sfc_emiss_0

real, dimension(size(pp,3) )                         :: tcol1d
real, dimension(size(pph,3) )                        :: bh1d, fuh1d, fdh1d

real, dimension(1)  ::   tsfc1d, bsfc1d, radiance1d

real, dimension(size(pp,3))    :: tau1d, sscat1d, gfac1d


id= size(pp,1); jd= size(pp,2); kd= size(pp,3)
kp= kd+1

!          Form the cose of the emission angle
!   cose= cos( (ema*2.0*pi/360.0) )
cose= 1.0

xcenter(1)= 7.8
xcenter(2)= 25.0
xcenter(3)= 32.0

!   Note that the dust and ice extinctions have been normalized
!          with respect to the visible.

!      Dust optical properties:  2.0 micron radius  @  7.8 25 mm and 32 mm
qex_dust(1)= 0.1446
qex_dust(2)= 0.3098
qex_dust(3)= 0.1500            !

sscat_dust(1)= 0.6698
sscat_dust(2)= 0.1649
sscat_dust(3)= 0.100       !

gfac_dust(1) = 0.7254
gfac_dust(2) = 0.2347
gfac_dust(3) = 0.1500       !

!   Ice optical properties:   4.0 mm    Need to calculate values for 32 microns

qex_ice(1)= 0.8372
qex_ice(2)= 0.3063      !   this seems low for 20 microns
qex_ice(3)= 0.1500      !  32 microns    check this

sscat_ice(1)= 0.7291
sscat_ice(2)= 0.5796
sscat_ice(3)= 0.5000

gfac_ice(1) = 0.8204
gfac_ice(2) = 0.3277
gfac_ice(3) = 0.3000


!        Assume vertical emission: appropriate for a nadir instrument
!!!!  cose= 1.0

do nn= 1, nbands   !  ---- begin loop over rad bands ---------


xwl= xcenter(nn)   !   microns

!      Assume vertical emission: appropriate for a nadir instrument  (ie TES)
!            For 32 micron channel (MCS), assume a slant path
    if( xwl > 30.0 ) then
        cose= 0.3746
    else
        cose= 1.000
    endif


#ifdef SKIP

!          Apply emissivity function for wavelengths >  18 microns
    if( xwl < 18.0 ) then
        sfc_emiss_0(:,:)= 1.0
    else
        sfc_emiss_0(:,:)= sfc_emiss(:,:)
    endif


    tau  (:,:,:)= dust(:,:,:) * qex_dust(nn)
    sscat(:,:,:)= sscat_dust(nn)
    gfac (:,:,:)= gfac_dust(nn)

    where( cloud(:,:,:) .gt. dust(:,:,:) )
        tau  (:,:,:) = cloud(:,:,:) * qex_ice(nn)
        sscat(:,:,:) = sscat_ice(nn)
        gfac (:,:,:) = gfac_ice(nn)
    end where

call planck_func( id*jd*kd, tcol, xwl, bh   )
call planck_func( id*jd,    tsfc, xwl, bsfc )

call lw_scattering ( id*jd, kd, tau, sscat, gfac, sfc_emiss_0, &
                           bh, bsfc, fuh, fdh, cose  )

!   set a minimum value to avoid division by very tiny number
radiance(:,:)= MAX(1.e-15,fuh(:,:,1) / (2.0*pi ))

call planck_inv_func( id*jd, radiance, xwl, tbands(:,:,nn)  )


#else


DO ii= 1, id
DO jj= 1, jd

    tcol1d(:)= tcol(ii,jj,:)
    tsfc1d= tsfc(ii,jj)

    tau1d(:)= dust(ii,jj,:) * qex_dust(nn)
    sscat1d(:)= sscat_dust(nn)
    gfac1d(:)= gfac_dust(nn)

!          Apply emissivity function for wavelengths >  18 microns
    if( xwl < 18.0 ) then
        sfc_emiss_0(1,1)= 1.0
    else
        sfc_emiss_0(1,1)= sfc_emiss(ii,jj)
    endif

    call planck_func( kd,   tcol1d, xwl, bh1d   )
    call planck_func( 1,    tsfc1d, xwl, bsfc1d )

    call lw_scattering ( 1, kd, tau1d, sscat1d, gfac1d, sfc_emiss_0, &
                               bh1d, bsfc1d, fuh1d, fdh1d, cose  )

!   set a minimum value to avoid division by very tiny number
    radiance1d= MAX(1.e-15,fuh1d(1) / (2.0*pi ))

    call planck_inv_func( 1, radiance1d, xwl, tsfc1d  )
    tbands(ii,jj,nn)= tsfc1d(1)

enddo
enddo

#endif SKIP


enddo   !  ------------- End Loop over rad bands ---------



return
end subroutine sfc_rad_flx_diag


!=======================================================================
!=======================================================================


subroutine radiation_driver_end

character (len=128) :: filename

mcpu0=  (mpp_pe() == mpp_root_pe())

filename= 'RESTART/radiation.res.nc'
if(mcpu0) print *,'Writing  restart file:   ',  trim( filename )
#ifndef ONED_MODEL
call write_data( trim( filename ), 'heatrate'   , heatrate)
call write_data( trim( filename ), 'sfc_sw_flux', sfc_sw_flux)
call write_data( trim( filename ), 'sfc_lw_flux', sfc_lw_flux)
call write_data( trim( filename ), 'tstrat', tstrat)
#endif
deallocate ( heatrate  )
deallocate ( sfc_sw_flux )
deallocate ( sfc_lw_flux )
deallocate ( tstrat )
deallocate ( taudust_reff )


return
end subroutine radiation_driver_end

!#######################################################################

end module radiation_driver_mod
