!*******************************************************************************************
!* Notices:
!*
!* Copyright © 2023 United States Government as represented by the Administrator of the
!* National Aeronautics and Space Administration.  All Rights Reserved.
!*
!* Disclaimers
!*
!* No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND,
!* EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT
!* THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY,
!* FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT
!* SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
!* THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY
!* GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
!* PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT
!* AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE
!* ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
!*
!* Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES
!* GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S
!* USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES
!* ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S
!* USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT,
!* ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.
!* RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
!*
!*******************************************************************************************
module mars_physics_mod


! *****************************************************************
!  Module that contains the main physics driver and initialization
!  for the NASA Ames Mars GCM
! *****************************************************************


use constants_mod,          only: KAPPA, CP_AIR, RDGAS, GRAV, PI, RADIAN,               &
                                  STEFAN, seconds_per_day
use astronomy_mod,          only: astronomy_init,  mars_calender
use time_manager_mod,       only: time_type, get_time

use mpp_domains_mod,        only: domain2d

use fms_mod,                only: error_mesg, FATAL,                                    &
                                  open_namelist_file, check_nml_error,                  &
                                  mpp_pe, mpp_root_pe, close_file,                      &
                                  write_version_number, stdlog,                         &
                                  uppercase, read_data, write_data, field_size

use fms2_io_mod,            only: file_exists
use time_manager_mod,       only: time_type, get_time
use diag_manager_mod,       only: register_diag_field, send_data
use field_manager_mod,      only: MODEL_ATMOS, parse, find_field_index
use tracer_manager_mod,     only: query_method, get_tracer_index,                       &
                                  get_number_tracers,get_tracer_names

use simple_rad_driver_mod,  only: simple_rad_driver,  simple_rad_driver_init,           &
                                  simple_rad_driver_end
use radiation_driver_mod,   only: radiation_driver,  radiation_driver_init,             &
                                  radiation_driver_end

use aerosol_mod,            only: aerosol_init
use initracer_mod
#ifdef fv3_turb
use vert_turb_driver_mod,   only: vert_turb_driver, vert_turb_driver_init,              &
                                  vert_turb_driver_end
use vert_diff_driver_mod,   only: vert_diff_driver_init, vert_diff_driver_end,          &
                                  surf_diff_type,                                       &
                                  vert_diff_driver_up, vert_diff_driver_down,           &
                                  vert_diff_driver
#endif
use surface_flux_mod,       only: surface_flux_2d

use mars_surface_mod,       only: mars_surface_init, mars_surface_end,                  &
                                  progts, albedo_calc,                                  &
                                  tsoil, t_surf, sfc_snow, sfc_frost,                   &
                                  id_frost, sfc_roughness, sfc_frost_mom,             &
                                  id_frost_mom, sfc_h2o2_chem, id_sfc_h2o2_chem

use update_mars_atmos_mod,  only: pass2, co2_condense, legacy_convect
use testconserv_mod

use ames_pbl_interface,     only: ames_pbl
use pblmod_mgcm,            only: amespbl_nml_read
use aerosol_util_mod,       only: do_moment_dust, do_moment_water, init_aerosol_flags, do_moment_sedim
use  sedim_mod,             only: sedim_driver,sedim_driver_init
use  dust_update_mod
use palmer_topo_drag_mod,   only: palmer_drag, palmer_drag_init, palmer_drag_end

#ifdef RELEASE
use null_physics_mod, only: init_null_phys
use null_physics_mod, only: micro_driver, micro_driver_init, &
          coagul_main, coagul_init, &
          dust_source_sink, dust_source_init, dust_source_end
use null_physics_mod, only: cg_drag_init, cg_drag_calc, cg_drag_end
use null_physics_mod, only: cloud_physics, cloud_physics_init, cloud_physics_end, &
		  id_wcol, id_cldcol, id_cld_r, id_rsat, cldcol, wcol, &
          photochem_driver
#else
use  micromod_mgcm,         only: micro_driver,micro_driver_init
use  coagulation_mod,       only: coagul_main,coagul_init 
use cloud_physics_mod,      only: cloud_physics_init, cloud_physics_end, &
                                 cloud_physics
use cg_drag_mod,            only: cg_drag_init, cg_drag_calc, cg_drag_end
use photochem_mod,          only: photochem_driver 
use dust_source_mod,        only: dust_source_sink, dust_source_init, &
                                 dust_source_end
use tagging_method_mod,     only: tagging_main

#endif


! *****************************************************************


implicit none
private


! *****************************************************************
! ************************ INTERFACES******************************
! *****************************************************************


public  :: mars_physics, mars_physics_init, mars_physics_end, do_mars_surface,          &
          do_bin_water_cycle, do_dust_source_sink
public  :: gascol
public  :: do_qmass


! *****************************************************************
! ************************ NAMELIST *******************************
! *****************************************************************

logical :: no_forcing = .false.                ! flag: skip mars_physics
logical :: do_mars_surface = .true.            ! flag: do mars surface calculations

real    :: t_zero=315., t_strat=200., delh=60., delv=10., eps=0., sigma_b=0.7

real    :: ka = -4.                            ! Newton Damping timescale. Negative = damping time in days
real    :: ks = -4.                            ! Special case of Newton timescale. Negative = damping time in days
real    :: kf = -4000.e6                       ! Rayleigh Damping timescale. Negative = damping time in days

logical :: do_conserve_energy = .false.        ! flag: adjust temperature according to wind drag
real    :: trflux = 1.e-5                      ! surface flux for optional tracer
real    :: trsink = -4.                        ! damping time for tracer

logical :: sponge_flag = .false.               ! flag: Rayleigh damping at the top of model
real    :: sponge_tau_days  = 1.0              ! damping time scale for the sponge (days)
real    :: rflevel = 0.1                       ! Raleigh damping pressure level (Pa)
real    :: sponge_pbottom = 1.0                ! bottom of sponge layer. Below this Rayleigh damping is zero (Pa)

logical :: sponge_flag2 = .false.              ! flag: Rayleigh damping at the top of model
real    :: sponge_tau_days2  = 1.0             ! damping time scale for the sponge (days)
real    :: rflevel2 = 0.1                      ! Raleigh damping pressure level (Pa)
real    :: sponge_pbottom2= 1.0                ! bottom of sponge layer. Below this Rayleight damping is zero (Pa)
real    :: rfwidth2= 30.0                      ! bottom of sponge layer. Below this Rayleight damping is zero (Pa)

logical :: do_vert_diff = .true.               ! flag: do vertical diffusion
logical :: do_ames_pbl = .true.                ! flag: do ames MY2.0 PBL

logical :: do_mars_radiation = .true.          ! flag: do mars radiation
logical :: do_simple_radiation = .false.       ! flag: do simple mars radiation
logical :: do_teq_z = .false.                  ! flag: do newton damping forcing
logical :: do_convective_adjust = .true.       ! flag: do convective adjustment
logical :: do_fv3_convect = .true.             ! flag: do original full column adjustment

logical :: do_co2_condensation_cycle = .true.  ! flag: do mass feedback from condensation
logical :: do_co2_condensation = .true.        ! flag: do atmospheric CO2 condensation
logical :: do_atmos_co2_loss = .true.          ! flag: do atmospheric CO2 loss from condensation
integer :: klevs_co2 = 5                       ! number of layers to spread CO2 loss to avoid single layer depletion

real    :: tau_diffusion = 1800.0              ! time scale for diffusion
logical :: diffusion_smooth = .true.           ! flag: smooth diffusion over tau_diffusion time scale

logical  :: freeze_tracer_fields = .false.      ! do not apply physics tendencies to tracers
logical  :: do_dust_source_sink = .false.        ! do binned dust source sink
logical  :: do_bin_water_cycle      = .false.       ! do binned water cycle
logical  :: checkcons = .false.                 ! check tracer conservation
real     :: gust_fact = 0.0                     ! gust_fact**3 is wind gust addition for dust stress lifting
logical  :: do_coagul_dst = .false.             ! Dust moment coagulation

logical :: inj_vap_pbl = .true.                ! flag: inject water vapor throughout PBL
real    :: facsubl = 1.0                       ! water vapor sublimation factor

logical :: do_pchem = .false.                  ! flag: for photochemistry, set to .true.
logical :: do_qmass = .false.                  ! flag: to use multiple gases in total atm mass, set to .true.

integer  :: GW_drag_TOG = 0                     ! toogle switch gravity wave drag
real     :: topo_drag_fac = 1.0                 ! scaling factor for topographic gravity wave


namelist /mars_physics_nml/  no_forcing, t_zero, t_strat, delh, delv, eps,              &
                             sigma_b, ka, ks, kf, do_conserve_energy,                   &
                             trflux, trsink,                                            &
                             sponge_flag,sponge_pbottom,sponge_tau_days,                &
                             tau_diffusion,  rflevel,                                   &
                             do_vert_diff, do_mars_radiation,                           &
                             do_simple_radiation, do_teq_z,                             &
                             do_mars_surface, do_convective_adjust,                     &
                             GW_drag_TOG, diffusion_smooth,                             &
                             do_co2_condensation_cycle, do_co2_condensation,            &
                             do_dust_source_sink, do_bin_water_cycle,                   &
                             gust_fact, topo_drag_fac, freeze_tracer_fields,            &
                             klevs_co2, do_atmos_co2_loss, checkcons,                   &
                             inj_vap_pbl, facsubl, do_ames_pbl,                         &
                             do_fv3_convect, do_pchem, do_qmass, do_coagul_dst,         &
                             sponge_flag2,sponge_pbottom2,sponge_tau_days2,             &
                             rflevel2, rfwidth2


!

character(len=12)  :: mod_name = 'mars_physics'
character(len=128) :: version  = '$Id: mars_physics.F90,v 1.1.2.1.2.1.4.1 2016/03/30 15:21:05 rjw Exp $'
character(len=128) :: tagname  = '$Name: mars_feb2012_rjw $'

real    :: tka, tks, vkf
real    :: trdamp

real    :: missing_value = -1.e10
logical :: module_is_initialized = .false.

! For registering output
integer :: id_teq, id_tdt, id_tdt2, id_udt, id_vdt,id_tdt_diss, id_diss_heat,           &
           id_tdt_hrad, id_tdt_pbl, id_tdt_adj, id_tdt_micro
integer :: id_difft_ames, id_diffm_ames
integer :: id_vmass, id_lheat,id_rkh
integer :: id_tsfc, id_stress, id_stress_gw, id_udt_cgwd, id_vdt_cgwd,id_tdt_cgwd,      &
           id_zpbl, id_difft,id_ppbl,id_kpbl,id_sens,id_wflux_vap, id_wflux_vap_bin,    &
           id_precip
integer :: id_udt_topo, id_vdt_topo
integer :: id_uv, id_vt, id_uw, id_dnflux
integer :: id_rdt_h2o2, id_cond_mass,id_rdt_subl_bin
integer :: id_delz

integer, dimension(:), allocatable  :: id_rdt_pbl, id_rdt_adj, id_rdt_dst,              &
                                       id_rdt_micro, id_rdt_dif, id_rdt_hrad,           &
                                       id_rdt_subl, id_rdt_sedim, id_rdt_coag
integer, dimension(:), allocatable  :: id_gascol, id_pchem
real, dimension(:,:,:), allocatable :: gascol

!! TO BE SAVED
real, dimension(:,:,:), allocatable, save   ::  diff_m_new, diff_t_new
real, dimension(:,:,:), allocatable, save   ::  teq, zteq, teqin, pteq
real, dimension(:,:), allocatable, save     ::  rhouch_save
integer, dimension(:,:),allocatable, save   ::  k_pbl_save

logical ::  mcpu0


!-----------------------------------------------------------------------


contains


!#######################################################################
!#######################################################################
!#######################################################################


subroutine mars_physics (is, js, dt, Time, lon, lat,                                    &
                         p_half, p_full,                                                &
                         z_half, z_full,                                                &
                         u, v, t, r, um, vm, tm, rm, udt, vdt, tdt, rdt,                &
                         dmass, p_ref, mask, kbot)
!
!  Main physics driver for Mars GCM.
!  Calls:
!
!    rayleigh_damping
!
!    albedo_calc (for radiation)
!
!    radiation_driver
!
!    update_water (water sublimation)
!
!    surface_flux2d (for diffusion)
!
!    ames_pbl / vert_turb and vert_diff
!
!    progts (surface temperature and CO2 condensation)
!
!    co2_condense (atmospheric CO2 condensation)
!
!    gravity wave drag
!
!    pass2 / legacy_convect (convective adjustment)
!
!    dust_update or dust_source_sink
!
!    binned cloud physics
!
!    sedim_driver and micro_driver
!

integer, intent(in)                     :: is, js               ! horizontal starting indices
real, intent(in)                        :: dt                   ! physics time step [s]
type(time_type), intent(in)             :: Time                 ! model time [s]

real, intent(in), dimension(:,:)        :: lon                  ! grid longitudes [rad]
real, intent(in), dimension(:,:)        :: lat                  ! grid latitudes [rad]
real, intent(in), dimension(:,:,:)      :: p_half,      &       ! interface pressures [Pa]
                                           p_full               ! midpoint pressures [Pa]
real, intent(in), dimension(:,:,:)      :: z_half,      &       ! interface heights above areoid [z]
                                           z_full               ! midpoint heights above areoid [z]

real, intent(in), dimension(:,:,:)      :: u,           &       ! u wind [m/s]
                                           v,           &       ! v wind [m/s]
                                           t,           &       ! temperature [K]
                                           um,          &       ! u wind [m/s]. same input as u
                                           vm,          &       ! v wind [m/s]. same input as v
                                           tm                   ! temperature [K]. same input as t
real, intent(in), dimension(:,:,:,:)    :: r,           &       ! tracer field
                                           rm                   ! tracer field. same input as r

real, intent(in), dimension(:,:,:), optional  :: mask           ! mask some fields
integer, intent(in), dimension(:,:), optional :: kbot           ! bottom layer of model

real, intent(inout), dimension(:,:,:)   :: udt,         &       ! u wind tendency [m/s/s]
                                           vdt,         &       ! v wind tendency [m/s/s]
                                           tdt                  ! temperature tendency [K/s]
real, intent(inout), dimension(:,:,:,:) :: rdt                  ! tracer tendency [*/s]

real, intent(out), dimension(:,:,:)     :: dmass                ! mass tendency over time step [Pa]
real, intent(in)                        :: p_ref


!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------


! tendencies
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_rad, rdt_pbl, rdt_adj,  &
                                                            rdt_dst, rdt_micro, rdt_dif,&
                                                            rdt_cld, rdt0, rdt_subl,    &
                                                            rdt_subl_bin, rdt_sedim,    &
                                                            rdt_dss, rdt_coag
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_pchem
real, dimension(size(r,1),size(r,2),size(r,3))  :: tdt_pbl, udt_pbl, vdt_pbl, tdt_rad, tdt_micro
real, dimension(size(r,1),size(r,2),size(r,3))  :: tdt_adj, udt_adj, vdt_adj
real, dimension(size(r,1),size(r,2),size(r,3))  :: tdt_dif, udt_dif, vdt_dif
real, dimension(size(r,1),size(r,2),size(r,3))  :: tdt_top, udt_top, vdt_top
real, dimension(size(r,1),size(r,2),size(r,3))  :: tdt0, udt0, vdt0
real, dimension(size(r,1),size(r,2),size(r,3))  :: tdt_co2, tdt_dst
real, dimension(size(r,1),size(r,2),size(r,4))  :: mass_wat_dt
real, dimension(size(r,1),size(r,2),size(r,3))  :: rdt_h2o2
real, dimension(size(r,1),size(r,2),size(r,3))  :: cond_mass

real, dimension(size(t,1),size(t,2))            :: ps,          &       ! surface pressure
                                                   diss_heat            ! rayleigh damping disssipative heat [W/m2]
real, dimension(size(t,1),size(t,2))            :: wflux_vap            ! sublimation flux [kg/m2/s]
real, dimension(size(t,1),size(t,2),size(t,3))  :: ttnd,        &       ! temperature tendency [K/s]
                                                   utnd,        &       ! u wind tendency [m/s/s]
                                                   vtnd,        &       ! v wind tendency [m/s/s]
                                                   pmass                ! layer thicknesses [Pa]
real, dimension(size(t,1),size(t,2),size(t,3))  :: diff_m,      &       ! momentum diffusivity
                                                   diff_t,      &       ! heat diffusivity
                                                   delp,        &       ! layer thicknesses [Pa]
                                                   delz                 ! layer thicknesses [m]
real, dimension(size(r,1),size(r,2),size(r,3))  :: uwnd,        &       ! u wind at next time step[m/s]
                                                   vwnd                 ! v wind at next time step[m/s]

real, dimension(size(t,1),size(t,2),2)          :: taucloud,    &       ! column cloud opacity
                                                   taudust,     &       ! column total dust opacity
                                                   taudust_mom, &       ! column moment dust opacity
                                                   taudust_fix          ! column fixed dust opacity
real, dimension(size(t,1),size(t,2),size(t,3))  :: qarray,      &       ! water vapor array [kg/kg]
                                                   qdt                  ! water vapor tendency [kg/kg/s]
real, dimension(size(t,1),size(t,2),size(t,3),size(r,4)) :: rini        ! initial tracer field

real, dimension(size(t,1),size(t,2))            :: swfsfc,      &       ! surface visible flux [W/m2]
                                                   lwfsfc               ! surface infrared flux [W/m2]
logical, dimension(size(t,1),size(t,2),size(t,3)+1) :: lmask            ! local mask field
real, dimension(size(t,1),size(t,2),size(t,3)+1)    :: outtmp           ! temporary output array

real, dimension(size(t,1),size(t,2),size(tsoil,3))  :: tgrnd,   &       ! soil temperature [K]
                                                       tg_tnd           ! soil temperature tendency [K/s]
real, dimension(size(t,1),size(t,2))    :: tsurf,               &       ! surface temperature [K]
                                           dnflux,              &       ! downward radiation at surface [W/m2]
                                           subday,              &       ! CO2 sublimation rate at surface [Pa]
                                           mratio,              &       ! ratio of condensation mass to column mass
                                           warray                       ! column to remove atmospheric mass
real, dimension(size(t,1),size(t,2))    :: snow,                &       ! surface CO2 [kg/m2]
                                           frost,               &       ! bin micro surface water ice [kg/m2]
                                           albedo,              &       ! albedo [-]
                                           sfc_emiss,           &       ! emissivity [-]
                                           cosz                         ! cosine of zenith angle
real, dimension(size(t,1),size(t,2))    :: frost_mom,         &       ! moment micro surface water ice [kg/m2]
                                           frost_eff                    ! bin or moment depending on compiler flag [kg/m2]
real, dimension(size(t,1),size(t,2), nice_mass)  :: frost_mom0, &     ! initial water ice field for conservation [kg/m2]
                                                    frost_field         ! temporary water ice field for updates [kg/m2]
real, dimension(size(t,1),size(t,2), ndust_mass) :: sfc_dst_mom0      ! initial dust field for conservation [kg/m2]
real, dimension(size(t,1),size(t,2))    :: sfc_h2o2_chem0               ! initial h2o2 concentration [kg/m2]

real, dimension(size(t,1),size(t,2))    :: shflx,               &       ! sensible heat flux [W/m2]
                                           precip,              &       ! condensed atmospheric CO2 [kg/m2]
                                           rhokd,               &       ! bottom layer air density [kg/m3]
                                           wind                         ! bottom layer wind speed [m/s]
real, dimension(size(t,1),size(t,2))    :: zo,                  &       ! z0 surface roughness [m]
                                           theta,               &       ! bottom layer potential temperature [K]
                                           zkd                          ! bottom layer thickness [m]
real, dimension(size(t,1),size(t,2))    :: dragm,               &       ! momentum drag
                                           dragh,               &       ! heat drag
                                           drag_q,              &       ! water vapor drag
                                           u_star,              &       ! friction velocity [m/s]
                                           b_star,              &       ! friction heat
                                           q_star,              &       ! friction water vapor
                                           cd_m,                &       ! momentum drag coefficient
                                           cd_t,                &       ! heat drag coefficient
                                           cd_q,                &       ! water vapor drag coefficient
                                           rhouch                       ! Rho*Cp*Cdh*Ustar
real, dimension(size(t,1),size(t,2))    :: stress,              &       ! surface stress
                                           source_mom                   ! moment dust source
real, dimension(size(t,1),size(t,2))    :: gust,                &       ! wind gust adjustment
                                           z_pbl,               &       ! height of pbl top [m]
                                           frac_land,           &       ! fraction of grid occupied by land
                                           xx2d                         ! qstar used by edt
real, dimension(size(t,1),size(t,2))    :: p_pbl                        ! pressure of pbl top [Pa]
integer, dimension(size(t,1),size(t,2)) :: k_pbl                        ! layer of pbl top

real, dimension(size(t,1),size(t,2))    :: dtau_du,             &       ! stress slope wrt u wind
                                           dtau_dv,             &       ! stress slope wrt v wind
                                           tau_x,               &       ! stress along u
                                           tau_y                        ! stress along v

real, dimension(size(t,1),size(t,2))    :: sens,                &       ! sensible heat flux [W/m2]
                                           evap,                &       ! evaporative heat flux [W/m2]
                                           dsens_datmos,        &       ! derivative of sensible heat wrt atmosphere temperature
                                           devap_datmos                 ! derivative of latent heat wrt atmosphere temperature
real, dimension(size(t,1),size(t,2))    :: dsens_dsurf,         &       ! derivative of sensible heat wrt surface temperature
                                           dedt_surf,           &       ! derivative of entrainment and diagnostic turbulence wrt surface temperature
                                           dedq_surf,           &       ! derivative of latent heat flux wrt surface temperature
                                           flux_r,              &       ! black body flux from surface [W/m2]
                                           drdt_surf                    ! derivative of black body flux from surface [W/m2]
real, dimension(size(t,1),size(t,2),size(t,3)) :: xx3d,         &       ! tdtlw used by edt code
                                                  tdtlw                 ! infrared heating [K/s]

logical, dimension(size(t,1),size(t,2)) :: convect                      ! used by entrainment in boundary layer parameterization
type(time_type)                         :: Time_next                    ! time of next time step
real                :: secs,                                    &       ! model time in seconds
                       sols                                             ! model sol number
integer             :: days,                                    &       ! model day number
                       seconds                                          ! seconds of current day
real, dimension(size(t,1),size(t,2),size(t,3))      :: tcol             ! current adjusted atmosphere temperature [K]
real, dimension(size(t,1),size(t,2),2*size(t,3)+1)  :: rkh              ! eddy mixing coefficient (m2/s)

logical, dimension(size(r,1),size(r,2),size(r,4))   :: lifting_dust     ! flag is true where lifting occurs

integer             :: ie, je, id, jd, kd, klev, i, j, k, kb, n, ntrace, ntp, nt, ndx
integer             :: nh2o_mom, nh2o_bin, nh2o2, nice_mom
logical             :: used
real                :: flux, sink, value
real                :: alpha                                            ! diffusion smoothing factor
character(len=128)  :: scheme, params
#ifdef fv3_turb
type(surf_diff_type) :: Surf_diff
#endif


!-----------------------------------------------------------------------
!    Initializations
!-----------------------------------------------------------------------


if (.not.module_is_initialized) call error_mesg ('mars_physics','mars_physics_init has not been called', FATAL)


tdt0      = 0.d0
udt0      = 0.d0
vdt0      = 0.d0
rdt0      = 0.d0
rdt       = 0.d0
rdt_pbl   = 0.d0
rdt_adj   = 0.d0
rdt_dst   = 0.d0
rdt_dss   = 0.d0
rdt_dif   = 0.d0
rdt_rad   = 0.0
rdt_micro = 0.d0
rdt_h2o2  = 0.d0
cond_mass = 0.d0
rdt_pchem = 0.d0
rdt_sedim = 0.d0
rdt_coag  = 0.d0
rdt_subl  = 0.d0
rdt_subl_bin = 0.d0
tdt       = 0.d0
tdt_pbl   = 0.d0
tdt_adj   = 0.d0
tdt_rad   = 0.d0
tdt_dif   = 0.d0
tdt_micro = 0.d0
tdt_co2   = 0.d0
udt       = 0.d0
udt_pbl   = 0.d0
udt_dif   = 0.d0
udt_adj   = 0.d0
vdt       = 0.d0
vdt_pbl   = 0.d0
vdt_dif   = 0.d0
vdt_adj   = 0.d0
udt_dif   = 0.d0
vdt_dif   = 0.d0
tg_tnd    = 0.d0

lifting_dust = .false.

id = size(lat,1)
jd = size(lat,2)
kd = size(p_full,3)
ie = is + id - 1
je = js + jd - 1

dmass(:,:,:) = 0.d0
mratio(:,:)  = 0.d0
precip(:,:)  = 0.d0


call get_number_tracers (MODEL_ATMOS, num_tracers=ntrace, num_prog=ntp)

!! Need first time to get first Ls and first tau
call get_time(Time, seconds, days)
secs = days * seconds_per_day + seconds
sols = secs / seconds_per_day

nh2o_mom = find_field_index(MODEL_ATMOS, 'vap_mass_mom')
nh2o_bin = find_field_index(MODEL_ATMOS, 'h2o_vapor')
nice_mom = find_field_index(MODEL_ATMOS, 'ice_mass_mom')

if (do_mars_surface) then
if (do_moment_water) then
    nh2o2          = find_field_index(MODEL_ATMOS, 'h2o2_mmr_gas')
    frost_mom0   = sfc_frost_mom(is:ie,js:je,:)
    sfc_dst_mom0 = sfc_dust_mass(is:ie,js:je,:)
    frost_eff(:,:) = sfc_frost_mom(is:ie,js:je,1)
    sfc_h2o2_chem0 = sfc_h2o2_chem(is:ie,js:je)
else
    frost_eff(:,:) = sfc_frost(is:ie,js:je,1)
endif
endif

if (checkcons) then
    call testneg(is, ie, js, je, kd, r+rdt*dt, nice_mom, -1.d-15, 'initest_ima')
    call testneg(is, ie, js, je, kd, r+rdt*dt, nh2o_mom, -1.d-15, 'initest_vap')
    call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0,      &
                      p_half, r+rdt*dt, sfc_frost_mom, sfc_dust_mass, 'initest')
endif

if (no_forcing) then
    if (do_vert_diff) then
        dnflux = 0.0
        tsurf  = 0.0
        call diffuse(is, js, dt, Time, lon, lat,        &
                p_half, p_full, z_half, z_full, tsurf,  &
                u, v, t, r(:,:,:,1:ntp),                &
                u, v, t, r(:,:,:,1:ntp),                &
                udt, vdt, tdt, rdt(:,:,:,1:ntp),        &
                dnflux, mask, kbot)
    endif
    return
endif


!-----------------------------------------------------------------------
!    Atmospheric and surface pressure
!-----------------------------------------------------------------------


if (present(kbot)) then
    do j=1, size(p_half,2)
        do i=1, size(p_half,1)
            kb      = kbot(i,j)
            ps(i,j) = p_half(i,j,kb+1)
        enddo
    enddo
else
    ps(:,:) = p_half(:,:,size(p_half,3))
endif

do k=1, kd
    delp(:,:,k) = p_half(:,:,k+1) - p_half(:,:,k)
    delz(:,:,k) = z_half(:,:,k) - z_half(:,:,k+1)
enddo

if (id_delz > 0) used = send_data (id_delz, delz, Time, is, js)


!-----------------------------------------------------------------------
!    Rayleigh damping of wind components near the surface and at the model top
!-----------------------------------------------------------------------


call rayleigh_damping (ps, p_full, lat, u+udt*dt, v+vdt*dt, utnd, vtnd, mask=mask)

if (do_conserve_energy) then
    ttnd = -((um+.5*utnd*dt)*utnd + (vm+.5*vtnd*dt)*vtnd)/CP_AIR
    tdt  = tdt + ttnd
    if (id_tdt_diss > 0) used = send_data (id_tdt_diss, ttnd, Time, is, js)
! vertical integral of ke dissipation
    if (id_diss_heat > 0) then
        do k = 1, size(t,3)
            pmass(:,:,k) = p_half(:,:,k+1)-p_half(:,:,k)
        enddo
        diss_heat = CP_AIR/GRAV * sum(ttnd*pmass, 3)
        if (id_diss_heat > 0) used = send_data (id_diss_heat, diss_heat, Time, is, js)
    endif
endif

udt = udt + utnd
vdt = vdt + vtnd

if (id_udt > 0) used = send_data (id_udt, utnd, Time, is, js)
if (id_vdt > 0) used = send_data (id_vdt, vtnd, Time, is, js)


!-----------------------------------------------------------------------
!    Surface conditions
!-----------------------------------------------------------------------


if (do_mars_surface) then
    snow(:,:)    = sfc_snow(is:ie,js:je)
    tsurf(:,:)   = tsoil(is:ie,js:je,1)
    tgrnd(:,:,:) = tsoil(is:ie,js:je,:)
    zo(:,:)      = sfc_roughness(is:ie,js:je)
    cosz(:,:)    = 1.0                                      ! Currently, albedo does not depend on solar zenith angle
    call albedo_calc(is, js, lon, lat, cosz, tsurf, ps,     &
                     snow, frost_eff, albedo, sfc_emiss)
else
    tsurf(:,:)     = 170.0
    albedo(:,:)    = 0.25
    sfc_emiss(:,:) = 1.0
endif


!-----------------------------------------------------------------------
!    Radiation
!-----------------------------------------------------------------------


!  Old version: This tracer tendency array allows the radiation code to
!  modify the tracer arrays. This is currently used to allow opacity "assimilation."

tdtlw(:,:,:) = 0.0
ttnd(:,:,:)  = 0.0


!  --------------call radiation or newtonian damping for temperature tendencies

if( do_mars_radiation  ) then

    call radiation_driver ( is, js, lon, lat, dt, Time, p_half, p_full, z_half,   &
                           tsurf, albedo, sfc_emiss, t, r, tdt, rdt,   &
                           swfsfc, lwfsfc, cosz, tdtlw, tdt_rad, taudust, &
                           taucloud, taudust_mom, taudust_fix, p_ref )
    dnflux = swfsfc + lwfsfc

else if (do_simple_radiation) then
    call simple_rad_driver  (is, js, lon, lat, dt, Time, p_half, p_full,                    &
                            tsurf, albedo, sfc_emiss, t+tdt*dt, r+rdt*dt, tdt_rad, rdt_rad, &
                            swfsfc, lwfsfc, cosz)
    dnflux = swfsfc + lwfsfc

else if (do_teq_z)  then
    ! thermal forcing for held & suarez (1994) benchmark calculation
    call newtonian_damping_z (lat, z_full, t+tdt*dt, zteq, teqin, teq, tdt_rad, mask)
    if (id_teq > 0) used = send_data (id_teq, teq(is:ie,js:je,:), Time, is, js)
    dnflux(:,:) = 0.0

else
    call newtonian_damping (lat, ps, p_full, t+tdt*dt, teq(is:ie,js:je,:), tdt_rad, mask)
    if (id_teq > 0) used = send_data (id_teq, teq(is:ie,js:je,:), Time, is, js)
    dnflux(:,:) = 0.0

endif

tdt = tdt + tdt_rad
rdt = rdt + rdt_rad

if (id_dnflux > 0) used = send_data (id_dnflux, dnflux, Time, is, js)


!-----------------------------------------------------------------------
!    PBL, Diffusion, CO2 cycle
!-----------------------------------------------------------------------


if (do_mars_surface) then
    rhokd(:,:) = ps(:,:) / (rdgas * t(:,:,kd))
    zkd(:,:)   = (ps(:,:) - p_full(:,:,kd)) / (grav*rhokd(:,:))

    !************************************************
    ! UPDATE WATER VAPOR : Sublimation / Condensation
    !************************************************
    frost_field(:,:,:) = sfc_frost_mom(is:ie,js:je,:)
    call update_water(is, ie, js, je, kd, lat, lon, dt, p_half, tsurf, rhouch_save,         &
                      k_pbl_save, r(:, :, :, :), rdt(:, :, :, :), nh2o_mom, frost_field,    &
                      wflux_vap, rdt_subl, sols)
    sfc_frost_mom(is:ie,js:je,:) = frost_field

    if (id_wflux_vap > 0) used = send_data (id_wflux_vap, wflux_vap, Time, is, js)

    frost_field(:,:,:) = sfc_frost(is:ie,js:je,:)
    call update_water_fv3(is, ie, js, je, kd, lat, lon, dt, p_half, tsurf, rhouch_save,     &
                          k_pbl_save, r(:, :, :, :), rdt(:, :, :, :), nh2o_bin,             &
                          frost_field(:, :, 1), wflux_vap, rdt_subl_bin)
    sfc_frost(is:ie,js:je,1) = frost_field(:,:,1)

    if (id_wflux_vap_bin > 0) used = send_data (id_wflux_vap_bin, wflux_vap, Time, is, js)

    rdt(:,:,:,:) = rdt(:,:,:,:) + rdt_subl(:,:,:,:) + rdt_subl_bin(:,:,:,:)

    if (checkcons) then
        call testneg(is, ie, js, je, kd, r+rdt*dt, nice_mom, -1.d-15, 'updatew_ima')
        call testneg(is, ie, js, je, kd, r+rdt*dt, nh2o_mom, -1.d-15, 'updatew_vap')
        call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0,      &
                          p_half, r+rdt*dt, sfc_frost_mom, sfc_dust_mass, 'updatew')
    endif


    !************************************************
    ! SURFACE FLUX: Description TB19
    !************************************************
    qarray(:,:,:) = r(:,:,:,nh2o_bin) + rdt(:,:,:,nh2o_bin) * dt
    qdt(:,:,:)    = rdt(:,:,:,nh2o_bin)
    gust(:,:)     = 1.0

    call surface_flux_2d (lon, lat,                             &
        t(:,:,kd), qarray(:,:,kd), u(:,:,kd), v(:,:,kd),        &
        p_full(:,:,kd), zkd, ps, tsurf, qarray(:,:,kd),         &
        zo, zo, zo, zo, gust,                                   &
        sens, evap, flux_r, tau_x, tau_y,                       &
        cd_m, cd_t, cd_q,                                       &
        wind, u_star, b_star, q_star,                           &
        dsens_dsurf, dedt_surf, dedq_surf, drdt_surf,           &
        dsens_datmos, devap_datmos, dtau_du,  dtau_dv,          &
        is, js, dt, Time)

    ! dragm, dragh and drag_q are potentially used by other mars physics
    dragm  = rhokd * cd_m * wind
    dragh  = rhokd * cd_t * wind
    drag_q = rhokd * cd_q * wind
    rhouch = rhokd * cd_t * u_star / sqrt(cd_m(:,:))      ! for legacy rhouch/cp

    ! For now, zero out moisture fluxes
    ! These would be used by the soil model also
    evap(:,:) = 0.0
    devap_datmos(:,:) = 0.0


    !************************************************
    ! VERTICAL DIFFUSION
    !************************************************
    ! set up vertical diffusion----this includes coupling with the surface----

    Time_next       = Time          ! This is used by the Mellor-Yamada code to compute the
                                    ! time step for the tke field
    frac_land(:,:)  = 1.0
    convect(:,:)    = .false.
    xx2d(:,:)       = 0.0
    xx3d(:,:,:)     = 0.0

    if (do_vert_diff) then          !-------vertical diffusion; including BL------

        !-----------------------------------
        !!! VERTICAL DIFFUSION WITH AMES PBL
        !-----------------------------------
        if (do_ames_pbl) then
            rkh = 0.0

            call ames_pbl(is, js, dt, ps,                               &
                   p_half, p_full, z_half,                              &
                   diff_m, diff_t, &
                   u, v, t, tsurf, r(:,:,:,:),                          & !qarray,     &
                   snow, frost_eff, tdt_rad, tau_x, tau_y,              &
                   sens, evap,                                          &
                   udt, vdt, tdt, rdt(:,:,:,:),                         &
                   udt_pbl, vdt_pbl, tdt_pbl, rdt_pbl(:,:,:,:),         &
                   u_star,b_star,cd_m,cd_t,dsens_datmos,dsens_dsurf,    &
                   wind,z_pbl,p_pbl,k_pbl,rkh,Time)

            dragm(:,:)  = rhokd(:,:) * cd_m(:,:) * wind(:,:)
            dragh(:,:)  = rhokd(:,:) * cd_t(:,:) * wind(:,:)
            drag_q(:,:) = rhokd(:,:) * cd_t(:,:) * wind(:,:)
            rhouch(:,:) = rhokd(:,:) * cd_t(:,:) * u_star(:,:) / sqrt(cd_m(:,:))    ! this is legacy rhouch/cp
            dtau_du     = -dragm
            dtau_dv     = -dragm
            diff_m_new(is:ie,js:je,:) = diff_m(:,:,:)
            diff_t_new(is:ie,js:je,:) = diff_t(:,:,:)
            lmask(:,:,1:kd) = .true.
            lmask(:,:,kd+1) = .false.

            outtmp(:,:,1:kd) = diff_t(:,:,1:kd)
            outtmp(:,:,kd+1) = 0.0
            if (id_difft_ames > 0) used = send_data (id_difft_ames, outtmp(:,:,:), Time, is, js, 1, mask=lmask)
            outtmp(is:ie,js:je,1:kd) = diff_m(:,:,1:kd)
            if (id_diffm_ames > 0) used = send_data (id_diffm_ames, outtmp(:,:,:), Time, is, js, 1, mask=lmask)

            rdt = rdt + rdt_pbl
            tdt = tdt + tdt_pbl
            udt = udt + udt_pbl
            vdt = vdt + vdt_pbl

            if (checkcons) then
                call testneg(is, ie, js, je, kd, r+rdt*dt, nice_mom, -1.d-15, 'amespbl_ima')
                call testneg(is, ie, js, je, kd, r+rdt*dt, nh2o_mom, -1.d-15, 'amespbl_vap')
                call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0,      &
                                  p_half, r+rdt*dt, sfc_frost_mom, sfc_dust_mass(:,:,:), 'amespbl')
            endif
        else   ! ames pbl
#ifdef fv3_turb

!--------------------------------------
!!! VERTICAL DIFFUSION WITH FV3 SCHEME
!--------------------------------------
            qarray  = r(:,:,:,nh2o_bin) + rdt(:,:,:,nh2o_bin) * dt
            qdt     = rdt(:,:,:,nh2o_bin)
            call vert_turb_driver(is, js, Time, Time, dt, xx3d,                     &
                                  frac_land,                                        &
                                  p_half, p_full, z_half, z_full,                   &
                                  u_star, b_star, xx2d, zo, lat, convect,           &
                                  u+udt*dt, v+vdt*dt, t+tdt*dt, qarray,             &
                                  r (:,:,:,1:ntp)+rdt(:,:,:,1:ntp)*dt,              &
                                  um, vm, tm, qarray, rm(:,:,:,1:ntp),              &
                                  udt_pbl, vdt_pbl, tdt_pbl, qdt, rdt_pbl,          &
                                  diff_t, diff_m,                                   &
                                  gust, z_pbl, p_pbl, k_pbl, mask=mask, kbot=kbot)

            rdt = rdt + rdt_pbl
            tdt = tdt + tdt_pbl
            udt = udt + udt_pbl
            vdt = vdt + vdt_pbl

! PBL eddy coefficient
            rkh(:,:,:) = 0.
!rkh(:,:,:)=diff_t(:,:,:)

            if (diffusion_smooth)  then
                alpha = dt / tau_diffusion
                diff_m_new(is:ie,js:je,:) = (diff_m_new(is:ie,js:je,:) + alpha*diff_m(:,:,:)) &
                                          / (1.0 + alpha)
                diff_t_new(is:ie,js:je,:) = (diff_t_new(is:ie,js:je,:) + alpha*diff_t(:,:,:)) &
                                          / (1.0 + alpha)
            else
                diff_m_new(is:ie,js:je,:) = diff_m(:,:,:)
                diff_t_new(is:ie,js:je,:) = diff_t(:,:,:)
            endif

            ttnd    = tdt
            rini    = r(:,:,:,:) + rdt(:,:,:,:) * dt
            qarray  = rini(:,:,:,nh2o_bin)
            qdt     = rdt(:,:,:,nh2o_bin)
            shflx   = sens

            call vert_diff_driver(is, js, Time_next, dt,                            &
                      p_half, p_full, z_full,                                       &
                      diff_m_new(is:ie,js:je,:),                                    &
                      diff_t_new(is:ie,js:je,:),                                    &
                      u+udt*dt, v+vdt*dt, t+tdt*dt, qarray, rini(:,:,:,1:ntp),      &
                      dtau_du, dtau_dv, tau_x, tau_y,                               &
                      dsens_datmos, devap_datmos,                                   &
                      shflx, evap,                                                  &
                      udt_dif, vdt_dif, tdt_dif, qdt, rdt_dif,                      &
                      mask=mask, kbot=kbot)

            rdt = rdt + rdt_dif
            tdt = tdt + tdt_dif
            udt = udt + udt_dif
            vdt = vdt + vdt_dif

            ttnd = tdt - ttnd
            if (id_difft > 0) used = send_data (id_difft, ttnd, Time, is, js)
#endif
        endif   ! ames pbl

    else    !-----------------no vertical diffusion
        sens(:,:) = 0.0 ! no sensible heat flux calculated
    endif   ! endif vertical diffusion

    !! Saving parameters
    rhouch_save(:,:) = rhouch(:,:)
    k_pbl_save(:,:)  = k_pbl(:,:)

    if (id_zpbl > 0) used = send_data (id_zpbl, z_pbl, Time, is, js)
    if (id_ppbl > 0) used = send_data (id_ppbl, p_pbl, Time, is, js)
    if (id_kpbl > 0) used = send_data (id_kpbl, float(k_pbl), Time, is, js)


!************************************************
!                 SURFACE MODEL
!************************************************
    shflx(:,:)  = -sens(:,:)
    subday(:,:) = 0.d0

    ! subday is the atmospheric mass gain or loss due to sublimation at the
    ! surface. It is to be added to atmospheric mass change within the model layers.
    ! dnflx and shflx are positive when directed downward
    ! Progts assumes downward flux has a positive sign

    if (checkcons) then
        call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0,      &
                          p_half, r+rdt*dt, sfc_frost_mom, sfc_dust_mass(:,:,:), 'preprog')
    endif

    ! Save Tsfc field at the beginning of the timestep, before it is updated by progts
    ! This is the temperature seen by the radiation code

        if (id_tsfc > 0) used = send_data (id_tsfc, tsurf, Time, is, js)


        call progts(is, js, dt, Time, lon, lat,  &
                ps, p_half, dnflux, tgrnd, snow, subday, tg_tnd, shflx, &
                dsens_datmos, dsens_dsurf, evap, devap_datmos, dedq_surf)

        snow=snow+subday
        tgrnd=tgrnd+tg_tnd
        if (checkcons) then
            call checkconserv(is,ie,js,je,kd,p_half,r,frost_mom0,sfc_dst_mom0,p_half, &
                            r+rdt*dt,sfc_frost_mom,sfc_dust_mass(:,:,:),'progtsf')
        endif

        tsoil   (is:ie,js:je,:)= tgrnd(:,:,:)
        sfc_snow(is:ie,js:je)=   snow(:,:)
        t_surf  (is:ie,js:je)=   tsoil(is:ie,js:je,1)

        tsurf(:,:)  = tsoil(is:ie,js:je,1)

        if (id_sens > 0) used = send_data (id_sens, sens, Time, is, js)

    !************************************************
    !                 CO2 CONDENSATION
    !************************************************
    ! Adjust atmospheric temperatures to account for CO2 condensation
    ! Mass changes are stored in dmass
    ! Augment surface CO2 snow (calculated in progts) with atmospheric
    ! CO2 precip calculated in co2_condense

    if (do_co2_condensation) then
        call co2_condense(is, js, dt, Time, t, tdt,  &
                          p_half, p_full, precip, dmass, tdt_co2)
        tdt = tdt + tdt_co2

        if (do_atmos_co2_loss) then
            sfc_snow(is:ie,js:je) = sfc_snow(is:ie,js:je) + precip(:,:)
        else
            dmass(:,:,:) = 0.d0
            precip(:,:)  = 0.d0
        endif

    else
        dmass(:,:,:) = 0.d0
        precip(:,:)  = 0.d0
    endif

    if (id_precip > 0) used = send_data (id_precip, precip, Time, is, js)

    if (do_co2_condensation_cycle) then
        ! Need to combine dmass and subday for boundary layer mass tendencies
        ! Distribute  mass loss due to CO2 condensation at the surface (subday)
        ! over the bottom N layers.  (Currently 5 layers; this should be based on
        ! "physical" considerations (say a vertical mixing length)

        ! For computational purposes, care is needed to prevent excessive
        ! mass loss within a single model layer:   (typically the bottom layer)
        ! grav*subday(:,:) / delp(:,:,k)  << 1.0

        ! In any case, it makes physical sense to "mix" mass loss at the surface through
        ! the "boundary layer" :   base the index on the zfull array

        ! zagl(:,:,k)=  z_full(:,:,k) - z_half(:,:,kd+1)
        ! find klev such that zagl(:,:,kmax)= , say, 5 km.

        klev = kd - klevs_co2 - 1
        warray(:,:) = 0.d0

        do k = klev, kd
            warray(:,:) = warray(:,:) + delp(:,:,k)
        enddo

        mratio(:,:) = subday(:,:) / warray(:,:)

        do k = klev, kd
            dmass(is:ie,js:je,k) = dmass(is:ie,js:je,k) + mratio(:,:)*delp(:,:,k)
        enddo

    else    ! No condensation cycle: Eliminate layer mass tendencies
        dmass(:,:,:) = 0.0
     endif

!***************** IF DO_SURFACE_MARS = FALSE
else    ! ---------- No Mars surface physics, hence just atmospheric diffusion  ------------------
     dmass(:,:,:) = 0.0
     if (do_vert_diff) then
        call diffuse(is, js, dt, Time, lon, lat,                                            &
                     p_half, p_full, z_half, z_full, tsurf,                                 &
                     u+udt*dt, v+vdt*dt, t+tdt*dt, r(:,:,:,1:ntp)+rdt(:,:,:,1:ntp)*dt,      &
                     u+udt*dt, v+vdt*dt, t+tdt*dt, r(:,:,:,1:ntp)+rdt(:,:,:,1:ntp)*dt,      &
                     udt_dif, vdt_dif, tdt_dif, rdt_dif (:,:,:,1:ntp) ,                     &
                     dnflux, mask, kbot)

        tdt = tdt + tdt_dif
        udt = udt + udt_dif
        vdt = vdt + vdt_dif

        rdt(:,:,:,1:ntp) = rdt(:,:,:,1:ntp) + rdt_dif (:,:,:,1:ntp)

        if (checkcons) then
            call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0,  &
                              p_half, r+rdt*dt, sfc_frost_mom, sfc_dust_mass(:,:,:), 'diffuse')
        endif
    endif
endif   ! do mars surface


!-----------------------------------------------------------------------
!   Photochemistry
!-----------------------------------------------------------------------


 if (do_pchem)   then   ! photochemistry is on
    mcpu0 = (mpp_pe() == mpp_root_pe())

    call photochem_driver(is, ie, js, je, kd, lon, lat, p_half, p_full, delp,               &
                           t, tdt, Time, dt, r, rdt, rdt_pchem, do_qmass)
    ! Tendencies:
    rdt = rdt + rdt_pchem

    !! test conservation tracers
    if (checkcons) then
      call ckconchem(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_h2o2_chem0, p_half,   &
                     r+rdt*dt, sfc_frost_mom, sfc_h2o2_chem, 'pkem')
    endif

  !************************************************
  ! UPDATE H2O2 : Condensation
  !************************************************
    call  h2o2cond(is, js, kd, lat, lon, dt, p_half, p_full, t, nh2o2, cond_mass, r, rdt, rdt_h2o2)

    ! Tendencies:
    rdt(:,:,:,nh2o2) = rdt(:,:,:,nh2o2) + rdt_h2o2(:,:,:)


    !! test conservation tracers
    if (checkcons) then
      call ckconchem(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_h2o2_chem0, p_half,   &
                     r+rdt*dt, sfc_frost_mom, sfc_h2o2_chem, 'H2O2')
    endif

    ! Output
    if (id_sfc_h2o2_chem > 0) used = send_data (id_sfc_h2o2_chem, sfc_h2o2_chem, Time, is, js)
    if (id_cond_mass > 0) used = send_data (id_cond_mass, cond_mass, Time, is, js)

endif


!-----------------------------------------------------------------------
!    Topo drag
!-----------------------------------------------------------------------


if (checkcons) then
    call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0, p_half,  &
                      r+rdt*dt, sfc_frost_mom, sfc_dust_mass(:,:,:), 'pretopo')
endif

if (nml_switch(GW_drag_TOG,1)) then     ! Palmer 1986
    tdt_top = 0.
    udt_top = 0.
    vdt_top = 0.
    mcpu0   = (mpp_pe() == mpp_root_pe())
    call palmer_drag(dt, u+udt*dt, v+vdt*dt, t+tdt*dt,      &
                     p_full, p_half, z_full, z_half,        &
                     udt_top, vdt_top, Time)

    udt = udt + udt_top * topo_drag_fac
    vdt = vdt + vdt_top * topo_drag_fac
endif

if (nml_switch(GW_drag_TOG,3)) then     ! Alexander 1999
    tdt_top = 0.
    udt_top = 0.
    vdt_top = 0.

    call cg_drag_calc(is, js, lat, p_full, z_full, t+tdt*dt, u+udt*dt, v+vdt*dt, Time, dt,  &
                      udt_top, vdt_top, tdt_top)
    udt =  udt + udt_top
    vdt =  vdt + vdt_top
    tdt=   tdt + tdt_top

    if (id_udt_cgwd > 0) used = send_data (id_udt_cgwd, udt_top, Time, is, js)
    if (id_vdt_cgwd > 0) used = send_data (id_vdt_cgwd, vdt_top, Time, is, js)
endif


!-----------------------------------------------------------------------
!    Convective adjustment
!-----------------------------------------------------------------------

if (checkcons) then
    call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0, p_half,  &
                      r+rdt*dt, sfc_frost_mom, sfc_dust_mass(:,:,:), 'prevadj')
endif

if (do_convective_adjust) then
    if (do_fv3_convect) then
        call pass2(is, js, dt, Time, u, v, t, r,        &
                   udt, vdt, tdt, rdt,                  &
                   udt_adj, vdt_adj, tdt_adj, rdt_adj,  &
                   delp, p_half, p_full)
    else
        call legacy_convect(is, js, dt, kd, u, v, t, r,         &
                            udt, vdt, tdt, rdt,                 &
                            udt_adj, vdt_adj, tdt_adj, rdt_adj, &
                            delp, p_half, p_full)
    endif

    udt = udt + udt_adj
    vdt = vdt + vdt_adj
    tdt = tdt + tdt_adj
    rdt = rdt + rdt_adj
endif

if (checkcons) then
    call testneg(is, ie, js, je, kd, r+rdt*dt, nice_mom, -1.d-15, 'convadj_ima')
    call testneg(is, ie, js, je, kd, r+rdt*dt, nh2o_mom, -1.d-15, 'convadj_vap')
    call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0,          &
                      p_half, r+rdt*dt, sfc_frost_mom, sfc_dust_mass(:,:,:), 'convadj')
endif


!-----------------------------------------------------------------------
!    Dust sources
!-----------------------------------------------------------------------

if (do_mars_surface) then
    snow(:,:) = sfc_snow(is:ie,js:je)
endif

!  Recalculate surface wind, including gustiness
!  Note that dragm has not be modified.
!  This is to capture the influence of gusts on dust lifting

uwnd = u + udt * dt
vwnd = v + vdt * dt

wind(:,:)   = sqrt(uwnd(:,:,kd) * uwnd(:,:,kd) + vwnd(:,:,kd) * vwnd(:,:,kd)    &
                               + gust_fact * gust(:,:) * gust(:,:))

stress(:,:) = dragm(:,:) * wind(:,:)

if (id_stress > 0) used = send_data (id_stress, stress, Time, is, js)

if (do_moment_dust) then
    call dust_update (is, js, lon, lat, dt, Time,                               &
                      p_half, p_full, p_pbl, k_pbl, tsurf, snow,                &
                      sfc_frost_mom(is:ie,js:je,1),stress, taudust_mom, taudust_fix, t,    &
                      tdt, r, rdt, sens, rdt_dst, source_mom, lifting_dust)
    rdt = rdt + rdt_dst

    if (checkcons) then
        call testneg(is, ie, js, je, kd, r+rdt*dt, nice_mom, -1.d-15, 'dustupd_ima')
        call testneg(is, ie, js, je, kd, r+rdt*dt, nh2o_mom, -1.d-15, 'dustupd_vap')
        call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0, p_half,  &
                        r+rdt*dt, sfc_frost_mom, sfc_dust_mass(:,:,:),'dustupd')
    endif

    if (do_coagul_dst) then
        call coagul_main(is, js, lon, lat, dt, t+tdt*dt, p_half, p_full, r+rdt*dt, rdt_coag)
        !! Update tracers
        rdt = rdt + rdt_coag
    endif

endif

if (do_dust_source_sink .and. ntrace > 1) then
    call dust_source_sink (is, js, lon, lat, dt, Time,                  &
                        p_half, p_full, tsurf,  snow,  stress,          &
                        k_pbl, source_mom, t, tdt, r, rdt, rdt_dss)
    rdt = rdt + rdt_dss
endif


!-----------------------------------------------------------------------
!    Clouds scheme
!-----------------------------------------------------------------------

if (do_bin_water_cycle) then
    do k = 1, kd
        tcol(:,:,k) = t(:,:,k) + tdt(:,:,k) * dt
    enddo
        rdt_cld = 0.
        frost(:,:) = sfc_frost(is:ie,js:je,1)
    call cloud_physics (is, js, lon, lat, dt, time, p_half, p_full,     &
                        tsurf, frost, tcol, r+rdt*dt, rdt_cld, drag_q)

    rdt = rdt + rdt_cld
    sfc_frost(is:ie,js:je,1) = frost(:,:)

    if (id_frost > 0) used = send_data (id_frost, frost, Time, is, js)
endif

if (do_moment_water) then
#ifndef RELEASE
    call micro_driver(is, js, Time, p_half, p_full, t, tdt, r, rdt, stress, tsurf, dt,  &
                      rhouch, rdt_micro, tdt_micro, checkcons)
    rdt = rdt + rdt_micro
    tdt = tdt + tdt_micro
#endif
    do nt = 1, nice_mass
        if (id_frost_mom(nt) > 0)  used = send_data (id_frost_mom(nt),              &
                                                       sfc_frost_mom(is:ie,js:je,nt), &
                                                       time, is, js)
    enddo

endif


if (do_moment_dust .and. do_moment_sedim) then 
    call sedim_driver(is, js, Time, p_half, p_full, t, tdt, r,     &
                                       rdt, tsurf, dt, rkh, rdt_sedim,              &
                                       lifting_dust, checkcons)
    rdt = rdt + rdt_sedim
endif


if (checkcons) then
    call testneg(is, ie, js, je, kd, r+rdt*dt, nice_mom, -1.d-15, 'microph_ima')
    call testneg(is, ie, js, je, kd, r+rdt*dt, nh2o_mom, -1.d-15, 'microph_vap')
    call checkconserv(is, ie, js, je, kd, p_half, r, frost_mom0, sfc_dst_mom0, p_half,  &
                    r+rdt*dt, sfc_frost_mom, sfc_dust_mass(:,:,:), 'microph')
endif

if (freeze_tracer_fields) then
    do  n = 1, ntrace
        rdt(:,:,:,n) =  0.0
    enddo
endif


!-----------------------------------------------------------------------
!  Get column of gas tracers
!-----------------------------------------------------------------------

do nt = 1, ntrace_gas
    ndx            = gas_indx(nt)
    gascol(:,:,nt) = 0.
    do k = 1, kd
        gascol(:,:,nt) = gascol(:,:,nt) + (r(:,:,k,ndx) + rdt(:,:,k,ndx) * dt) &
                       * ((p_half(:,:,k+1) - p_half(:,:,k)) / grav)
    enddo
enddo


!-----------------------------------------------------------------------
!    Outputs
!-----------------------------------------------------------------------
!  The surface_temperature, dust_source and cloud_physics modules
!  all introduce changes to the bottom boundary conditions for the
!  atmosphere which should be reflected in the upward sweep of the
!  diffusion

!  In particular, the fluxes of heat and water vapor are strongly coupled by
!  surface and near-surface atmospheric temperatures


!---- mass flux outputs ----
if (id_vmass > 0) then
    pmass(:,:,:) = v(:,:,:) * delp(:,:,:)
    used = send_data (id_vmass, pmass, Time, is, js)
endif

if (id_uv > 0) then
    pmass(:,:,:) = u(:,:,:) *v(:,:,:)
    used = send_data (id_uv, pmass, Time, is, js)
endif

if (id_vt > 0) then
    pmass(:,:,:) = v(:,:,:) * t(:,:,:)
    used = send_data (id_vt, pmass, Time, is, js)
endif

!if (id_uw > 0) then
!    pmass(:,:,:) = u(:,:,:) * v(:,:,:)     ! currently no w field
!    used = send_data (id_uw, pmass, Time, is, js)
!endif


!--- Tendencies for each physical process

!*** 1/ MIXING RATIOS
do n= 1, ntrace
    if (id_rdt_pbl(n) > 0)      used = send_data(id_rdt_pbl(n),                     &
                                                 rdt_pbl(is:ie,js:je,:,n),          &
                                                 time, is, js)
    if (id_rdt_adj(n) > 0)      used = send_data(id_rdt_adj(n),                     &
                                                 rdt_adj(is:ie,js:je,:,n),          &
                                                 time, is, js)
    if (id_rdt_dst(n) > 0)      used = send_data(id_rdt_dst(n),                     &
                                                 rdt_dst(is:ie,js:je,:,n),          &
                                                 time, is, js)
    if (id_rdt_micro(n) > 0)    used = send_data(id_rdt_micro(n),                   &
                                                 rdt_micro(is:ie,js:je,:,n),        &
                                                 time, is, js)
    if (id_rdt_sedim(n) > 0)    used = send_data(id_rdt_sedim(n),                   &
                                                 rdt_sedim(is:ie,js:je,:,n),        &
                                                 time, is, js)
    if (id_rdt_coag(n) > 0)     used = send_data(id_rdt_coag(n),                    &
                                                 rdt_coag(is:ie,js:je,:,n),         &
                                                 Time, is, js)
    if (id_rdt_hrad(n) > 0)     used = send_data (id_rdt_hrad(n),                   &
                                                 rdt_rad(is:ie,js:je,:,n),          &
                                                 time, is, js)
enddo
do nt = 1, nice_mass
    ndx=vapor_indx(nt)
    if (id_rdt_subl(nt) > 0) used = send_data (id_rdt_subl(nt), rdt_subl(:,:,:,ndx), time, is, js)
enddo

if (id_rdt_subl_bin > 0) used = send_data (id_rdt_subl_bin, rdt_subl_bin(:,:,:,nh2o_bin), time, is, js)

if (id_rdt_h2o2 > 0) used = send_data (id_rdt_h2o2, rdt_h2o2, Time, is, js)

do nt = 1, ntrace_gas
    ndx = gas_indx(nt)
    if (id_pchem(nt) > 0) used = send_data(id_pchem(nt),(r(:,:,:,ndx) + (dt*rdt(:,:,:,ndx))), Time, is, js)
!   if (id_pchem(nt) > 0) used = send_data(id_pchem(nt),rdt_pchem(:,:,:,ndx), Time, is, js)
enddo


!*** 2/ TEMPERATURES
if (id_tdt_pbl > 0)     used = send_data (id_tdt_pbl, tdt_pbl, Time, is, js)
if (id_tdt_adj > 0)     used = send_data (id_tdt_adj, tdt_adj, Time, is, js)
if (id_tdt_micro > 0)   used = send_data (id_tdt_micro, tdt_micro, Time, is, js)
if (id_tdt_hrad > 0)    used = send_data (id_tdt_hrad, tdt_rad, Time, is, js)
if (id_lheat > 0)       used = send_data (id_lheat, tdt_co2, Time, is, js)

!--- Other ---
do nt = 1, ntrace_gas
    if (id_gascol(nt) > 0)  used = send_data(id_gascol(nt), gascol(:,:,nt), time, is, js)
enddo

!--- END tendencies for each physical process


end subroutine mars_physics


!#######################################################################
!#######################################################################
!#######################################################################


subroutine mars_physics_init (nlon, mlat, nlevels, lonb, latb, lon, lat,        &
                          pstd, axes, Time, phys_domain)

! routine for initializing the model physics

integer, intent(in)              ::  nlon, mlat, nlevels
real, intent(in), dimension(:,:) ::  lonb, latb
real, intent(in), dimension(:,:) ::  lon, lat
real, intent(in), dimension(:)   ::  pstd
integer, intent(in)              :: axes(4)

integer, dimension(3)            :: half = (/1,2,4/)

type(time_type), intent(in)      :: Time

type(domain2d), intent(inout)    :: phys_domain

!-----------------------------------------------------------------------
integer unit, io, ierr, j, k, id, jd, is, js, ie, je, nt, ndx
integer km, fld_dims(4),ntrace,ntp

character (len=256)  :: filename, fieldname

logical              ::  lga, lgb
#ifdef fv3_turb
type(surf_diff_type) :: Surf_diff
#endif
character (len=256)  :: tname, tracer_name

mcpu0 = (mpp_pe() == mpp_root_pe())
call get_number_tracers (MODEL_ATMOS, num_tracers=ntrace, num_prog=ntp)

is = 1
js = 1

ie = is + size(lon,1) - 1
je = js + size(lon,2) - 1

id = size(lonb,1) - 1
jd = size(lonb,2) - 1

!     ----- read namelist -----

if (file_exists('input.nml')) then
    unit = open_namelist_file ()
    ierr = 1; do while (ierr /= 0)
        read  (unit, nml=mars_physics_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'mars_physics_nml')
    enddo
10     call close_file (unit)
endif

!     ----- write version info and namelist to log file -----

call write_version_number(version,tagname)
if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=mars_physics_nml)

if (mcpu0)  print *, 'Within mars_physics_init'


! *********************************************************
!                 Initialize microphysics
! *********************************************************

call init_aerosol_flags
call initmicro

#ifdef RELEASE
call init_null_phys(is, ie, js, je)
#endif


! *********************************************************
!          ----- register diagnostic fields -----
! *********************************************************

id_delz = register_diag_field (mod_name, 'delz',                                            &
            axes(1:3), Time, 'Model layer thickness in Z', 'm',                             &
            missing_value=missing_value )

!!! TRACERS MIXING RATIO TENDENCIES
allocate (id_rdt_dst(ntrace))
allocate (id_rdt_pbl(ntrace))
allocate (id_rdt_dif (ntrace))
allocate (id_rdt_adj(ntrace))
allocate (id_rdt_micro(ntrace))
allocate (id_rdt_sedim(ntrace))
allocate (id_rdt_coag(ntrace))
allocate (id_rdt_hrad(ntrace))
allocate (id_gascol(ntrace_gas))
allocate (id_pchem(ntrace_gas))
allocate (id_rdt_subl(nice_mass))
!allocate (id_rdt_subl_bin(ntrace))

!-------- Common to all tracers
do nt = 1, ntrace
    call get_tracer_names(MODEL_ATMOS, nt, tracer_name)

    tname = trim(tracer_name) // '_Tpbl'
    id_rdt_pbl(nt)  = register_diag_field (mod_name, trim(tname),                           &
                        axes(1:3), Time, 'Tracer Tendency PBL', 'kg/kg/s',                  &
                        missing_value=missing_value)

    tname = trim(tracer_name) // '_Tdif'
    id_rdt_dif(nt)  = register_diag_field (mod_name, trim(tname),                           &
                        axes(1:3), Time, 'Tracer Tendency Diffusion', 'kg/kg/s',            &
                        missing_value=missing_value)

    tname = trim(tracer_name) // '_Tadj'
    id_rdt_adj(nt)  = register_diag_field (mod_name, trim(tname),                           &
                        axes(1:3), Time, 'Tracer Tendency Convective Adjustment', 'kg/kg/s',&
                        missing_value=missing_value)

    tname = trim(tracer_name) // '_Tdst'
    id_rdt_dst(nt)  = register_diag_field (mod_name, trim(tname),                           &
                        axes(1:3), Time, 'Tracer Tendency Dust Update', 'kg/kg/s',          &
                        missing_value=missing_value)

    tname = trim(tracer_name) // '_Tmicro'
    id_rdt_micro(nt) = register_diag_field (mod_name, trim(tname),                          &
                        axes(1:3), Time, 'Tracer Tendency Moment Microphysics', 'kg/kg/s',  &
                        missing_value=missing_value)

    tname = trim(tracer_name) // '_Tsedi'
    id_rdt_sedim(nt) = register_diag_field (mod_name, trim(tname),                          &
                        axes(1:3), Time, 'Tracer Tendency Sedimentation', 'kg/kg/s',        &
                        missing_value=missing_value)

    tname = trim(tracer_name) // '_Tcoag'
    id_rdt_coag(nt) = register_diag_field (mod_name, trim(tname),                           &
                        axes(1:3), Time, 'Tracer Tendency Coagulation', 'kg/kg/s',          &
                        missing_value=missing_value)

    tname = trim(tracer_name) // '_Trad'
    id_rdt_hrad(nt) = register_diag_field (mod_name, trim(tname),                           &
                        axes(1:3), Time, 'Tracer Tendency Radiation', 'kg/kg/s',            &
                        missing_value=missing_value)
enddo


!-------- Common to all ice mass micro tracers
do nt = 1, nice_mass
    ndx = vapor_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)

    tname= trim(tracer_name) // '_Tssubl'
    id_rdt_subl(nt) = register_diag_field (mod_name, trim(tname),                           &
                        axes(1:3), Time, 'Water Moment Sublimation Tendency', 'kg/kg/s',    &
                        missing_value=missing_value)
enddo

!!! TEMPERATURE TENDENCIES
id_tdt      = register_diag_field (mod_name, 'tdt_ndamp',                                   &
                axes(1:3), Time, 'Temperature Tendency Newtonian Damping', 'K/sec',         &
                missing_value=missing_value)

id_difft    = register_diag_field (mod_name, 'difft',                                       &
                axes(1:3), Time, 'Temperature Tendency Vertical Diffusion', 'K/s',                   &
                missing_value=missing_value)

id_tdt_hrad = register_diag_field (mod_name, 'tdt_hrad',                                    &
                axes(1:3), Time, 'Temperature Tendency Radiation', 'K',                     &
                missing_value=missing_value)

id_tdt_pbl  = register_diag_field (mod_name, 'tdt_pbl',                                     &
                axes(1:3), Time, 'Temperature Tendency  Ames PBL', 'K/s',                        &
                missing_value=missing_value)

id_tdt_adj  = register_diag_field (mod_name, 'tdt_adj',                                     &
                axes(1:3), Time, 'Temperature Tendency Convective Adjustment', 'K/s',       &
                missing_value=missing_value)

id_tdt_micro = register_diag_field (mod_name, 'tdt_micro',                                  &
                axes(1:3), Time, 'Temperature Tendency Moment Microphysics', 'K/s',         &
                missing_value=missing_value)

id_lheat    = register_diag_field (mod_name, 'tdt_co2',                                     &
                axes(1:3), Time, 'Temperature Tendency CO2 Latent Heating ', 'K/s',         &
                missing_value=missing_value)

!!! DIFFUSION
id_difft_ames = register_diag_field(mod_name, 'difft_ames',                                 &
                axes(half),Time, 'Eddy Heat Mixing Coefficient', 'K/s',                     &
                missing_value=missing_value)

id_diffm_ames = register_diag_field(mod_name, 'diffm_ames',                                 &
                axes(half), Time, 'Eddy Momentum Mixing Coefficient', 'm2/s',               &
                missing_value=missing_value)

!!! OTHERS
id_teq      = register_diag_field (mod_name, 'teq',                                         &
                axes(1:3), Time, 'Equilibrium Temperature', 'K',                            &
                missing_value=missing_value)

id_udt      = register_diag_field (mod_name, 'udt_rdamp',                                   &
                axes(1:3), Time, 'Rayleigh Damping Zonal Wind', 'm/s/s',                    &
                missing_value=missing_value )

id_vdt      = register_diag_field (mod_name, 'vdt_rdamp',                                   &
                axes(1:3), Time, 'Rayleigh Damping Meridional Wind', 'm/s/s',               &
                missing_value=missing_value )

id_vmass    = register_diag_field (mod_name, 'vmass',                                       &
                axes(1:3), Time, 'Meridional Mass Flux', 'kg/s/s/s',                        &
                missing_value=missing_value)

id_uv       = register_diag_field (mod_name, 'uv',                                          &
                axes(1:3), Time, 'Horizontal Momentum Flux', 'm2/s/s',                              &
                missing_value=missing_value)

id_vt       = register_diag_field (mod_name, 'vt',                                          &
                axes(1:3), Time, 'Meridional Heat Flux', 'Km/s',                                    &
                missing_value=missing_value)

!id_uw       = register_diag_field (mod_name, 'uw',                                          &
!                axes(1:3), Time, 'Vertical Momentum Flux', 'kg/s/s/s',                      &
!                missing_value=missing_value)

!!! PBL
id_zpbl     = register_diag_field (mod_name, 'zpbl',                                        &
                axes(1:2), Time, 'Boundary Layer Depth in Z', 'm',                          &
                missing_value=missing_value )

id_kpbl     = register_diag_field (mod_name, 'kpbl',                                        &
                axes(1:2), Time, 'Boundary Layer Depth by Level', 'lev',                    &
                missing_value=missing_value )

id_ppbl     = register_diag_field (mod_name, 'ppbl',                                        &
                axes(1:2), Time, 'Boundary Layer Depth in P', 'mbar',                       &
                missing_value=missing_value )

!!! DISSIPATION
if (do_conserve_energy) then
    id_tdt_diss = register_diag_field (mod_name, 'tdt_diss_rdamp',                          &
                axes(1:3), Time, 'Dissipative Heating from Rayleigh Damping', 'K/sec',      &
                missing_value=missing_value)

    id_diss_heat = register_diag_field (mod_name, 'diss_heat_rdamp',                        &
                axes(1:2), Time,                                                            &
                'Column-Integrated Dissipative Heating from Rayleigh Damping', 'W/m/m')
endif

!!! BASIC SURFACE FIELDS
id_tsfc     = register_diag_field (mod_name, 'tsfc',                                        &
                axes(1:2), Time, 'Surface Temperature', 'K',                                &
                missing_value=missing_value)

id_dnflux   = register_diag_field (mod_name, 'dnflux',                                      &
                axes(1:2), Time, 'Net Downward Radiation Flux (IR+VIS)', 'W/m/m',           &
                missing_value=missing_value)

id_sens     = register_diag_field (mod_name, 'sens',                                        &
                axes(1:2), Time, 'Sensible Heat Flux', 'W/m/m',                             &
                missing_value=missing_value)

id_stress   = register_diag_field (mod_name, 'stress',                                      &
                axes(1:2), Time, 'Surface Stress', 'N/m/m',                                 &
                missing_value=missing_value)

id_precip   = register_diag_field (mod_name, 'precip',                                      &
                axes(1:2), Time, 'Atmospheric CO2 Snowfall per Timestep', 'kg/m/m',                      &
                missing_value=missing_value)

!!! WATER CYCLE
id_wflux_vap = register_diag_field (mod_name, 'wflux_vap',                                  &
                axes(1:2), Time, 'Water Moment Sublimation Flux', 'kg/m/m/s',               &
                missing_value=missing_value)

id_wflux_vap_bin = register_diag_field (mod_name, 'wflux_vap_bin',                          &
                axes(1:2), Time, 'Water Bulk Sublimation Flux', 'kg/m/m/s',                 &
                missing_value=missing_value)

!!! GAS COLUMNS
do nt = 1, ntrace_gas
    ndx = gas_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)

    tname = trim(tracer_name) // '_col'
    id_gascol(nt) = register_diag_field (mod_name, trim(tname),                             &
                        axes(1:2), Time, 'Gas Column Mass', 'kg/m/m',                       &
                        missing_value=missing_value)
enddo

!!! CHEMISTRY
id_rdt_h2o2 =  register_diag_field (mod_name, 'h2o2_Tcond',                                 &
                axes(1:3), Time, 'Tendency H2O2 Condensation', 'kg/kg/s',                   &
                missing_value=missing_value )

id_cond_mass = register_diag_field (mod_name, 'cond_mass',                                  &
                axes(1:3), Time, 'Condensed H2O2 Mass', 'kg/m2',                            &
                missing_value=missing_value )

!id_sfc_h2o2_chem = register_diag_field ('mars_surface', 'sfc_h202_chem',                    &
!                     (/axes(1:2)/), Time, 'Surface H2O2', 'kg/m/m',                         &
!                     missing_value=missing_value)

do nt = 1, ntrace_gas
    ndx = gas_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)
    tname = trim(tracer_name)
    id_pchem(nt) = register_diag_field (mod_name, trim(tname),                              &
                     axes(1:3), Time, trim(tname), 'mmr',                                   &
                     missing_value=missing_value)
enddo


!!! GRAVITY WAVES
if (nml_switch(GW_drag_TOG,1)) then
    if (mcpu0) print *, 'Calling palmer_drag_init ...'
    call palmer_drag_init(nlon, mlat, lonb, latb, axes, Time)

    id_stress_gw = register_diag_field (mod_name, 'base_flux',                              &
                    axes(1:2), Time, 'Surface Stress', 'N/m/m',                             &
                    missing_value=missing_value)
endif

if (nml_switch(GW_drag_TOG,3)) then
    if (mcpu0) print *,  'Calling  cg_drag_init ...'

    call cg_drag_init (lon, lat, pstd, Time=Time, axes=axes)

    id_udt_cgwd = register_diag_field (mod_name, 'udt_cgwd',                                &
                    axes(1:3), Time, 'U Tendency for c Gravity Wave Drag', 'm/s/s',         &
                    missing_value=missing_value)
    id_vdt_cgwd = register_diag_field (mod_name, 'vdt_cgwd',                                &
                    axes(1:3), Time, 'V Tendency for c Gravity Wave Drag', 'm/s/s',         &
                    missing_value=missing_value)
    id_tdt_cgwd = register_diag_field (mod_name, 'tdt_cgwd',                                &
                    axes(1:3), Time, 'Temperature Tendency for c Gravity Wave Drag',        &
                    'm/s/s', missing_value=missing_value)

    if (mcpu0) print *,'  Have returned from cg_drag_init '
endif

!-----------------------------------------------------------------------

if (do_vert_diff) then
#ifdef fv3_turb
    call vert_turb_driver_init (id, jd, nlevels, axes, Time, lga, lgb)

    call vert_diff_driver_init(Surf_diff, id, jd, nlevels, axes, Time)
#endif
    ! Allocate arrays for time-smoothing the vertical diffusivities
    allocate (diff_m_new(is:ie,js:je,nlevels))
    allocate (diff_t_new(is:ie,js:je,nlevels))
    allocate (rhouch_save(is:ie,js:je))
    allocate (k_pbl_save(is:ie,js:je))

    rhouch_save = 0.0
    k_pbl_save  = nlevels

    ! Get time-smoothed diffusivities from physics restart file
    filename = 'INPUT/physics.res.nc'
    if (file_exists(trim(filename))) then
        if (mcpu0) print *,'Reading  restart file:   ',  trim(filename)
        call read_data(trim(filename), 'diffm', diff_m_new)
        call read_data(trim(filename), 'difft', diff_t_new)
    else
        diff_m_new = 1.e-3
        diff_t_new = 1.e-3
    endif
endif


!
call amespbl_nml_read


! Initialize astronomy
call astronomy_init()

allocate (gascol (is:ie,js:je,ntrace_gas))

if (do_mars_surface) then
    call mars_surface_init(nlon, mlat, lonb, latb, lon, lat, axes, Time, phys_domain)
    if (mcpu0) print *,'Returned from surface init: '
end if

if (do_mars_radiation)  then
    if (mcpu0) print *,'Initializing aerosol and radiation'
    call aerosol_init(nlon, mlat, lonb, latb, lon, lat, axes, Time)
    call radiation_driver_init(nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time)

else if (do_simple_radiation) then
    call simple_rad_driver_init(nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time)

else if (do_teq_z) then
    allocate (teq(is:ie,js:je,nlevels))

    filename = 'INPUT/teq.nc'
    if (file_exists(trim(filename))) then
        ! first get z dimension and then allocate  zteq and teqin
        call field_size(trim(filename), 'teq', fld_dims)
        km = fld_dims(3)
        if (mcpu0) print *, 'Input Teq vertical dimension:  ', km

        allocate (teqin(is:ie,js:je,km))
        allocate (zteq(is:ie,js:je,km))

        call read_teq_z(filename, lon, lat, zteq(is:ie,js:je,:), teqin(is:ie,js:je,:))

    endif
    if (mcpu0) print *, 'Have initialized teq_z'
    if (mcpu0) print *, 'zteq_z ', is, js, zteq(is,js,:)

else
    allocate (teq(is:ie,js:je,nlevels))

    filename = 'INPUT/teq.nc'
    if (file_exists(trim(filename))) then
        call read_teq(filename, lon, lat, pstd, teq(is:ie,js:je,:))
    endif

    if (mcpu0) print *, 'Have initialized teq'
    if (mcpu0) print *, 'teq ', is,js,teq(is,js,:)

endif


if (do_bin_water_cycle) then
    call cloud_physics_init(nlon, mlat, lonb, latb, lon, lat, axes, Time)
    if (mcpu0) print *,'Returned from cloud_physics init: '
endif

if (do_dust_source_sink) then
    call dust_source_init(nlon, mlat, lonb, latb, lon, lat, axes, Time, phys_domain)
    if (mcpu0) print *,'Returned from dust_source init: '
endif

if (do_moment_dust) then
    call dust_update_init(nlon, mlat, lonb, latb, lon, lat, axes, Time, phys_domain)
    call sedim_driver_init(nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time)
#ifndef RELEASE
    if (do_moment_water) call micro_driver_init(nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time)
    if (do_coagul_dst) call coagul_init()
#endif
endif


!     ----- compute coefficients used in the newtonian damping subroutines-----
if (ka < 0.) ka = -86400. * ka
if (ks < 0.) ks = -86400. * ks
if (kf < 0.) kf = -86400. * kf

tka = 0.; if (ka > 0.) tka = 1. / ka
tks = 0.; if (ks > 0.) tks = 1. / ks
vkf = 0.; if (kf > 0.) vkf = 1. / kf

!               ----- for tracers -----
if (trsink < 0.) trsink = -86400. * trsink
trdamp = 0.0
if (trsink > 0.) trdamp = 1. / trsink

if (mcpu0) print *, 'mars_physics is Initialized ! '

module_is_initialized  = .true.

!-----------------------------------------------------------------------
end subroutine mars_physics_init


!#######################################################################
!#######################################################################
!#######################################################################


subroutine mars_physics_end(days)

! write restart files and release memory

integer, intent(in) :: days
character (len=256) :: filename

mcpu0 = (mpp_pe() == mpp_root_pe())

if (do_vert_diff) then
    filename = 'RESTART/physics.res.nc'
    if (mcpu0) print *,'Writing physics restart file: ',  trim(filename)
    call write_data(trim(filename), 'diffm', diff_m_new)
    call write_data(trim(filename), 'difft', diff_t_new)

    deallocate (diff_m_new)
    deallocate (diff_t_new)
    deallocate (rhouch_save)
    deallocate (k_pbl_save)
#ifdef fv3_turb
    call vert_turb_driver_end
    call vert_diff_driver_end
#endif

endif

if (do_mars_surface)     call mars_surface_end(days)

if (do_dust_source_sink) call dust_source_end(days)
if (do_moment_dust)     call dust_update_end

if (do_bin_water_cycle)  call cloud_physics_end

if (do_mars_radiation)   call radiation_driver_end

if (do_simple_radiation) call simple_rad_driver_end

if (mcpu0) print *,'Completed mars_physics_end: '

module_is_initialized = .false.


end subroutine mars_physics_end


!#######################################################################
!#######################################################################
!#######################################################################


subroutine newtonian_damping (lat, ps, p_full, t, teqin, tdt,  mask)

!   routine to compute thermal forcing for held & suarez (1994)
!   benchmark calculation.

real, intent(in), dimension(:,:)             :: lat, ps
real, intent(in), dimension(:,:,:)           :: p_full, t, teqin
real, intent(out), dimension(:,:,:)          :: tdt
real, intent(in), dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------
real, dimension(size(t,1),size(t,2),size(t,3)) :: tdamp

integer :: k
real    :: tcoeff, pref

do k = 1, size(t,3)
    tdamp(:,:,k) = tka
enddo

do k = 1, size(t,3)
    tdt(:,:,k) = -tdamp(:,:,k) * (t(:,:,k) - teqin(:,:,k))
enddo

if (present(mask)) then
    tdt = tdt * mask
endif

!-----------------------------------------------------------------------

end subroutine newtonian_damping


!#######################################################################
!#######################################################################
!#######################################################################


subroutine newtonian_damping_z (lat, z_full, t, zteq, teqin, teqout, tdt,  mask)

!   routine to compute thermal forcing from a Teq field specified as a function of Z
!
!               Need to re-interpolate in height as the Z field evolves

use axis_utils_mod, only: interp_1d

real, intent(in), dimension(:,:)             :: lat
real, intent(in), dimension(:,:,:)           :: z_full, t, zteq, teqin
real, intent(out), dimension(:,:,:)          :: teqout, tdt
real, intent(in), dimension(:,:,:), optional :: mask

real, dimension(size(t,1),size(t,2),size(t,3)) :: tdamp, zflip

integer :: k, kmax
real    :: tcoeff, pref

mcpu0 =  (mpp_pe() == mpp_root_pe())

kmax = size(t,3)

do k = 1, size(t,3)
    tdamp(:,:,k) = tka
enddo

! interp_1d requires a monotonically increasing axis; hence use -z_full

call interp_1d(-zteq, -z_full, teqin, teqout, "linear")

do k = 1, size(t,3)
    tdt(:,:,k) = -tdamp(:,:,k) * (t(:,:,k) - teqout(:,:,k))
enddo

if (present(mask)) then
    tdt = tdt * mask
endif

!-----------------------------------------------------------------------

end subroutine newtonian_damping_z


!#######################################################################
!#######################################################################
!#######################################################################


subroutine newtonian_damping_anal (lat, ps, p_full, t, tdt, p_ref, mask)

!   routine to compute thermal forcing for held & suarez (1994)
!   benchmark calculation.

real, intent(in), dimension(:,:)    :: lat, ps
real, intent(in), dimension(:,:,:)  :: p_full, t
real, intent(out), dimension(:,:,:) :: tdt
real, intent(in)                    :: p_ref
real, intent(in), dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------

real, dimension(size(t,1), size(t,2)) :: sin_lat, sin_lat_2, cos_lat_2, t_star,         &
                                         cos_lat_4, tstr, sigma, the, tfactr, rps, p_norm

real, dimension(size(t,1),size(t,2),size(t,3)) :: tdamp

integer :: k
real    :: tcoeff, pref

!-----------------------------------------------------------------------
!------------latitudinal constants--------------------------------------

sin_lat  (:,:) = sin(lat(:,:))
sin_lat_2(:,:) = sin_lat(:,:) * sin_lat(:,:)
cos_lat_2(:,:) = 1.0 - sin_lat_2(:,:)
cos_lat_4(:,:) = cos_lat_2(:,:) * cos_lat_2(:,:)

t_star(:,:)    = t_zero - delh*sin_lat_2(:,:) - eps*sin_lat(:,:)
tstr(:,:)      = t_strat - eps*sin_lat(:,:)

!-----------------------------------------------------------------------

tcoeff = (tks - tka) / (1.0 - sigma_b)
pref   = p_ref
rps    = 1./ps

do k = 1, size(t,3)
!  ----- compute equilibrium temperature (teq) -----
    p_norm(:,:) = p_full(:,:,k) / pref
    the(:,:)    = t_star(:,:) - delv*cos_lat_2(:,:) * log(p_norm(:,:))
    teq(:,:,k)  = the(:,:) * (p_norm(:,:))**KAPPA
    teq(:,:,k)  = max(teq(:,:,k), tstr(:,:))

!  ----- compute damping -----
    sigma(:,:) = p_full(:,:,k) * rps(:,:)
    where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
        tfactr(:,:) = tcoeff * (sigma(:,:) - sigma_b)
        tdamp(:,:,k) = tka + cos_lat_4(:,:) * tfactr(:,:)
    elsewhere
        tdamp(:,:,k) = tka
    endwhere
enddo

!*** note: if the following loop uses vector notation for all indices
!          then the code will not run ??????

do k = 1, size(t,3)
    tdt(:,:,k) = -tdamp(:,:,k) * (t(:,:,k) - teq(:,:,k))
enddo

if (present(mask)) then
    tdt = tdt * mask
    teq = teq * mask
endif

!-----------------------------------------------------------------------

end subroutine newtonian_damping_anal


!#######################################################################
!#######################################################################
!#######################################################################


subroutine rayleigh_damping (ps, p_full, lat, u, v, udt, vdt, mask)

!          rayleigh damping of wind components near surface
!               sigma_b < sigma < 1.0
!
!          Optional additional damping for p < sponge_bottom for the upper layers

real, intent(in), dimension(:,:)             :: ps
real, intent(in), dimension(:,:)             :: lat
real, intent(in), dimension(:,:,:)           :: p_full, u, v
real, intent(out), dimension(:,:,:)          :: udt, vdt
real, intent(in), dimension(:,:,:), optional :: mask

! local vars
real, dimension(size(u,1),size(u,2)) :: sigma, vfactr, rps
real, dimension(size(u,1),size(u,2)) :: ystruc, umx, vmx

integer :: i, j, k, id, jd
real    :: vcoeff, sponge_coeff

!----------------compute damping----------------------------------------

rps          = 1.0 / ps
vcoeff       = -vkf / (1.0 - sigma_b)
sponge_coeff = 1.0 / (sponge_tau_days * 86400.0)
udt(:,:,:)   = 0.0
vdt(:,:,:)   = 0.0

do k = 1, size(u,3)
    sigma(:,:) = p_full(:,:,k) * rps(:,:)
    where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
        vfactr(:,:) = vcoeff * (sigma(:,:) - sigma_b)
        udt(:,:,k)  = vfactr(:,:) * u(:,:,k)
        vdt(:,:,k)  = vfactr(:,:) * v(:,:,k)
    elsewhere
        udt(:,:,k) = 0.0
        vdt(:,:,k) = 0.0
    end where
enddo

if (sponge_flag) then
sponge_coeff = 1.0 / (sponge_tau_days * 86400.0)

do k = 1, size(u,3)
    sigma(:,:) = p_full(:,:,k) * rps(:,:)
    where (sigma(:,:) <  sponge_pbottom)
        vfactr(:,:) = 0.5 * (1.0 + tanh(1.5 * log(rflevel / p_full(:,:,k))))
        udt(:,:,k)  = udt(:,:,k) - u(:,:,k) * sponge_coeff * vfactr(:,:)
        vdt(:,:,k)  = vdt(:,:,k) - v(:,:,k) * sponge_coeff * vfactr(:,:)
    end where
enddo
endif


if (sponge_flag2) then
    sponge_coeff = 1.0 / (sponge_tau_days2 * 86400.0)
    ystruc(:,:)  = exp(-1.0 * (lat(:,:) * RADIAN / rfwidth2)**4)
    do k = 1, size(u,3)
        sigma(:,:) = p_full(:,:,k) * rps(:,:)
        !!  umx(:,:)= min(u(:,:,k), 100.0)
        !!      umx(:,:)= max(umx(:,:), -100.0)
        umx(:,:) = u(:,:,k)
        vmx(:,:) = v(:,:,k)
        where (sigma(:,:) < sponge_pbottom2)
            !!!  vfactr(:,:) = ystruc(:,:)*exp(-((log(p_full(:,:,k))-log(rflevel2))/2)**4);
            vfactr(:,:) =  0.5*ystruc(:,:)*(1.0 + tanh(0.8*log(rflevel2/p_full(:,:,k))))
            udt(:,:,k)  = udt(:,:,k) - umx(:,:) * sponge_coeff * vfactr(:,:)
            vdt(:,:,k)  = vdt(:,:,k) - vmx(:,:) * sponge_coeff * vfactr(:,:)
      end where
  enddo
endif

if (present(mask)) then
    udt = udt * mask
    vdt = vdt * mask
endif

!-----------------------------------------------------------------------

end subroutine rayleigh_damping


!#######################################################################
!#######################################################################
!#######################################################################


subroutine tracer_source_sink (flux, damp, p_half, r, rdt, kbot)
!
real, intent(in)                :: flux, damp, p_half(:,:,:), r(:,:,:)
real, intent(out)               :: rdt(:,:,:)
integer, intent(in), optional   :: kbot(:,:)
! local vars
real, dimension(size(r,1),size(r,2),size(r,3)) :: source, sink
real, dimension(size(r,1),size(r,2))           :: pmass

integer :: i, j, kb
real    :: rdamp
!-----------------------------------------------------------------------
rdamp = damp
if (rdamp < 0.) rdamp = -86400. * rdamp     ! convert days to seconds
if (rdamp > 0.) rdamp = 1. / rdamp

!------------ simple surface source and global sink --------------------
source(:,:,:) = 0.0

if (present(kbot)) then
    do j = 1, size(r,2)
        do i = 1, size(r,1)
            kb              = kbot(i,j)
            pmass(i,j)      = p_half(i,j,kb+1) - p_half(i,j,kb)
            source(i,j,kb)  = flux / pmass(i,j)
        enddo
    enddo
else
    kb              = size(r,3)
    pmass (:,:)     = p_half(:,:,kb+1) - p_half(:,:,kb)
    source(:,:,kb)  = flux / pmass(:,:)
endif

sink(:,:,:) = rdamp*r(:,:,:)
rdt(:,:,:)  = source(:,:,:) - sink(:,:,:)

!-----------------------------------------------------------------------

end subroutine tracer_source_sink


!#######################################################################
!#######################################################################
!#######################################################################


subroutine diffuse (is, js, dt, Time, lon, lat,                             &
                     p_half, p_full, z_half, z_full, tsurf,                 &
                     u, v, t, r, um, vm, tm, rm, udt, vdt, tdt, rdt,        &
                     dnflux, mask, kbot)

!   calculate the vertical diffusion and turbulence when no forcing or
!   no surface

integer, intent(in)                             :: is, js
real, intent(in)                                :: dt
type(time_type), intent(in)                     :: Time
real, intent(in), dimension(:,:)                :: lon
real, intent(in), dimension(:,:)                :: lat
real, intent(in), dimension(:,:,:)              :: p_half, p_full
real, intent(in), dimension(:,:,:)              :: z_half, z_full
real, intent(in), dimension(:,:)                :: tsurf
real, intent(in), dimension(:,:,:)              :: u, v, t, um, vm, tm
real, intent(in), dimension(:,:,:,:)            :: r, rm
real, intent(inout), dimension(:,:,:)           :: udt, vdt, tdt
real, intent(inout), dimension(:,:,:,:)         :: rdt
real, intent(in), dimension(:,:)                :: dnflux

real, intent(in), dimension(:,:,:), optional    :: mask
integer, intent(in), dimension(:,:) , optional  :: kbot

!  Local
real, dimension(size(t,1),size(t,2))            :: ps
real, dimension(size(t,1),size(t,2),size(t,3))  :: diff_m, diff_t
real, dimension(size(t,1),size(t,2),size(t,3))  :: ttnd

integer                                         ::   days, seconds
integer                                         :: i, j, k, kb, n, id, jd, kd, ie, je, ntp
logical                                         :: used
real                                            :: flux, sink, value, alpha
#ifdef fv3_turb
type(surf_diff_type)                            :: Surf_diff
#endif
real, dimension(size(t,1),size(t,2),size(t,3))  :: qarray, qdt

real, dimension(size(t,1),size(t,2))            ::  dragm, dragh, drag_q, u_star, b_star
real, dimension(size(t,1),size(t,2))            ::  stress
real, dimension(size(t,1),size(t,2))            ::  gust, z_pbl, frac_land, q_star
real, dimension(size(t,1),size(t,2))            ::  p_pbl
integer, dimension(size(t,1),size(t,2))         ::  k_pbl
real, dimension(size(t,1),size(t,2))            ::  dtau_du, dtau_dv, tau_x, tau_y
real, dimension(size(t,1),size(t,2))            ::  zo
real, dimension(size(t,1),size(t,2))            ::  sens, evap, dsens_datmos, devap_datmos

real, dimension(size(t,1),size(t,2),size(t,3))  :: tdtlw
logical, dimension(size(t,1),size(t,2))         :: convect
type(time_type)                                 :: Time_next

id  = size(lat,1)
jd  = size(lat,2)
kd  = size(p_full,3)
ntp = size(r,4)

ie  = is + size(lat,1) - 1
je  = js + size(lat,2) - 1

if (present(kbot)) then
    do j = 1, size(p_half,2)
        do i = 1, size(p_half,1)
            kb      = kbot(i,j)
            ps(i,j) = p_half(i,j,kb+1)
        enddo
    enddo
else
    ps(:,:) = p_half(:,:,size(p_half,3))
endif

Time_next = Time

qarray(:,:,:)  = 0.0
qdt(:,:,:)     = 0.0
! The following are not actually used; they are dummy arrays
frac_land(:,:) = 1.0
tdtlw(:,:,:)   = 0.0
convect(:,:)   = .false.
q_star(:,:)    = 0.0
u_star(:,:)    = 0.0
b_star(:,:)    = 0.0

zo(:,:)        = 0.01

#ifdef fv3_turb
call vert_turb_driver(is, js, Time, Time_next, dt, tdtlw,                       &
                      frac_land,                                                &
                      p_half, p_full, z_half, z_full,                           &
                      u_star, b_star, q_star, zo, lat, convect,                 &
                      u, v, t, qarray, r (:,:,:,1:ntp),                         &
                      um, vm, tm, qarray, rm(:,:,:,1:ntp),                      &
                      udt, vdt, tdt, qdt, rdt,                                  &
                      diff_t, diff_m,                                           &
                      gust, z_pbl, p_pbl, k_pbl, mask=mask,kbot=kbot)


if (diffusion_smooth) then
    alpha = dt / tau_diffusion
    diff_m_new(is:ie,js:je,:) = (diff_m_new(is:ie,js:je,:) + alpha*diff_m(:,:,:)) &
                              / (1.0 + alpha)
    diff_t_new(is:ie,js:je,:) = (diff_t_new(is:ie,js:je,:) + alpha*diff_t(:,:,:)) &
                              / (1.0 + alpha)
else
    diff_m_new(is:ie,js:je,:) = diff_m(:,:,:)
    diff_t_new(is:ie,js:je,:) = diff_t(:,:,:)
endif

dtau_du(:,:) = 0.0
dtau_dv(:,:) = 0.0
tau_x(:,:)   = 0.0
tau_y(:,:)   = 0.0

! Sensible heat flux (W/m**2)
sens(:,:)    =  0.0
dsens_datmos =  0.0

! For now, zero out moisture fluxes
evap(:,:)         = 0.0
devap_datmos(:,:) = 0.0

ttnd = tdt
call vert_diff_driver(is, js, Time_next, dt,                                    &
                      p_half, p_full,  z_full,                                  &
                      diff_m_new(is:ie,js:je,:),                                &
                      diff_t_new(is:ie,js:je,:),                                &
                      um, vm, tm, qarray, rm(:,:,:,1:ntp),                      &
                      dtau_du, dtau_dv, tau_x, tau_y,                           &
                      dsens_datmos, devap_datmos,                               &
                      sens, evap,                                               &
                      udt, vdt, tdt, qdt, rdt,                                  &
                      mask=mask, kbot=kbot)

ttnd = tdt - ttnd
#endif

end subroutine diffuse


!#######################################################################
!#######################################################################
!#######################################################################


subroutine read_teq (filename, lon, lat, phalf, tout)

! routine for initializing the radiative-convective temperature cross-section

use  horiz_interp_mod, only: horiz_interp

character(len=*), intent(in)            :: filename
real, intent(in), dimension(:,:)        :: lon, lat
real, intent(in), dimension(:)          :: phalf

real, intent(out), dimension(:,:,:)     :: tout

!-----------------------------------------------------------------------
integer                                 :: unit, io, ierr, i, j, k, nlevels, klev, kk
integer                                 :: im, jm, km, fld_dims(4)
real                                    :: frac

real, dimension(:,:,:), allocatable     :: teq_inpt
real, dimension(:,:),   allocatable     :: txy

real, dimension(:), allocatable         :: lat_inpt, pres_inpt,  presh_inpt
real, dimension(:), allocatable         :: lonb_inpt, latb_inpt

real, dimension(size(phalf)-1)          ::  pstd

!-----------------------------------------------------------------------
nlevels = size(phalf) - 1

call field_size(trim(filename), 'lat', fld_dims)

allocate(lat_inpt (fld_dims(1)))
allocate(latb_inpt(fld_dims(1)+1))

call read_data(trim(filename), 'lat',  lat_inpt,  no_domain=.true.)
call read_data(trim(filename), 'latb', latb_inpt, no_domain=.true.)

call field_size(trim(filename), 'lonb', fld_dims)

allocate(lonb_inpt (fld_dims(1)))

call read_data(trim(filename), 'lonb', lonb_inpt, no_domain=.true.)

call field_size(trim(filename), 'pfull', fld_dims)
allocate(pres_inpt(fld_dims(1)))
call read_data(trim(filename), 'pfull', pres_inpt, no_domain=.true.)
print *, 'pfull dims:  ', fld_dims
pres_inpt = 100.0 * pres_inpt
!-----------------------------------------------------------------------

call field_size(trim(filename), 'teq', fld_dims)

im = fld_dims(1); jm = fld_dims(2); km = fld_dims(3)
print *, 'Input Teq dims:  ', fld_dims

allocate(teq_inpt(im,jm,km))
allocate(txy(im,jm))

call read_data(trim(filename), 'teq', teq_inpt, no_domain=.true.)

latb_inpt(:) = latb_inpt(:) / RADIAN
lonb_inpt(:) = lonb_inpt(:) / RADIAN
!-----------------------------------------------------------------------

! formulate full model levels from input half-levels
do k = 1, nlevels
#ifdef SPEC_CORE
    pstd(k) = 0.5 * (phalf(k+1) + phalf(k))
#else
    pstd(k) = (phalf(k+1) - phalf(k)) / log(phalf(k+1) / phalf(k))
#endif SPEC_CORE
enddo
!-----------------------------------------------------------------------

! If km != nlevels  then require vertical interpolation
if (nlevels > km .or. nlevels < km)  then
    do k = 1, nlevels
        if (pres_inpt(1) > pstd(k)) then
            klev = 2; frac = 1.0

        else
            do kk = 2, km
                if (pres_inpt(kk) > pstd(k)) then
                    frac = (pres_inpt(kk) - pstd(k)) / (pres_inpt(kk) - pres_inpt(kk-1))
                    klev = kk
                    exit
                endif
            enddo

            if (kk > km)  then
                klev = km; frac = 0.0
            endif
        endif
        ! complete pressure interpolation
        do i = 1, im
            do j = 1, jm
                txy(i,j) = (1.0-frac) * teq_inpt(i,j,klev) + frac * teq_inpt(i,j,klev-1)
            enddo
        enddo

        call horiz_interp(txy(:,:), lonb_inpt, latb_inpt, lon, lat,                         &
                          tout(:,:,k), interp_method='bilinear')
    enddo      ! -----------  end loop over k

else
    do k = 1, nlevels
        call horiz_interp(teq_inpt(:,:,k), lonb_inpt, latb_inpt, lon, lat,                  &
                          tout(:,:,k), interp_method='bilinear')
    enddo
endif

deallocate (teq_inpt)
deallocate (txy)
deallocate (lat_inpt)
deallocate (latb_inpt)
deallocate (lonb_inpt)
deallocate (pres_inpt)

end subroutine read_teq


!#######################################################################
!#######################################################################
!#######################################################################


subroutine read_teq_z (filename, lon, lat, zgrd, tout)

! routine for initializing the radiative-convective temperature cross-section

use  horiz_interp_mod,  only: horiz_interp

character(len=*), intent(in)            :: filename
real, intent(in), dimension(:,:)        :: lon, lat
real, intent(out), dimension(:,:,:)     ::  zgrd, tout

!-----------------------------------------------------------------------
integer                                 :: k, im, jm, km, fld_dims(4)

real, dimension(:,:,:), allocatable     :: teq_inpt, zzf_inpt

real, dimension(:), allocatable         :: lat_inpt, pres_inpt, presh_inpt
real, dimension(:), allocatable         :: lonb_inpt, latb_inpt

call field_size(trim(filename), 'lat', fld_dims)

allocate(lat_inpt (fld_dims(1)))
allocate(latb_inpt(fld_dims(1)+1))

call read_data(trim(filename), 'lat',  lat_inpt,  no_domain=.true.)
call read_data(trim(filename), 'latb', latb_inpt, no_domain=.true.)

call field_size(trim(filename), 'lonb', fld_dims)

allocate(lonb_inpt (fld_dims(1)))

call read_data(trim(filename), 'lonb', lonb_inpt, no_domain=.true.)

latb_inpt(:) = latb_inpt(:) / RADIAN
lonb_inpt(:) = lonb_inpt(:) / RADIAN

call field_size(trim(filename), 'teq', fld_dims)

im = fld_dims(1); jm = fld_dims(2); km = fld_dims(3)

allocate(teq_inpt(im,jm,km))
call read_data(trim(filename), 'teq', teq_inpt, no_domain=.true.)

allocate(zzf_inpt(im,jm,km))
call read_data(trim(filename), 'zzf', zzf_inpt, no_domain=.true.)

! Need to flip indices so that zgrd is monotonically increasing
do k = 1, km
    call horiz_interp(teq_inpt(:,:,k), lonb_inpt, latb_inpt, lon, lat,          &
                     tout(:,:,k), interp_method='bilinear')

    call horiz_interp(zzf_inpt(:,:,k), lonb_inpt, latb_inpt, lon, lat,          &
                     zgrd(:,:,k), interp_method='bilinear')
enddo

deallocate (teq_inpt)
deallocate (zzf_inpt)
deallocate (lat_inpt)
deallocate (latb_inpt)
deallocate (lonb_inpt)

end subroutine read_teq_z


!#######################################################################
!#######################################################################
!#######################################################################


function nml_switch(int_nml, int_opt)
! toogle switch for namelist, used to replace the compiler flags
! -dtopo_drag -dcg_drag by a namelist switch: gw_drag_tog= 12,
!                                       !0: no gravity wave
!
!                                       !1: palmer 1986 scheme
!
!                                       !2: garner 2005 scheme    (not tested yet)
!
!                                       !3: alexander 1999 scheme (not tested yet)
!
! returns nml_switch = .true. if a digit of int_nml matches int_opt
!
! example: if int_nml =  0        --> no topo drag
!
!          if int_nml =  1        --> do palmer 1986 scheme
!
!          if int_nml = 12 or 21  --> do palmer 1986 scheme and garner 2005 scheme
!

implicit none

integer, intent(in)     :: int_nml          ! integer from namelist, 5 digit max
integer, intent(in)     :: int_opt          ! integer option ranging from 0->9
logical                 :: nml_switch
!---work variables---
character(5)            :: txt_nml
character(1)            :: txt_opt
integer i
nml_switch = .false.

write(txt_nml,'(i5)') int_nml
write(txt_opt,'(i1)') int_opt

do i = 1, 5
    if (txt_nml(i:i) .eq. txt_opt) then
        nml_switch = .true.
    endif
enddo

end function nml_switch


!#######################################################################
!#######################################################################
!#######################################################################


subroutine update_water(is, ie, js, je, nz, lat, lon, dt, pl, tg, drg, kpbl, r,     &
                        rdt, nh2o, qpig, wflux, rdt_subl, sols)

! update the water vapor field for the moment microphysics scheme

use constants_mod, only: grav

implicit none

real, intent(in)                        :: dt                                           ! physical time step [s]
integer, intent(in)                     :: nz,is,js,ie,je
real, intent(in), dimension(:,:)        :: lat                                          ! latitude [rad]
real, intent(in), dimension(:,:)        :: lon                                          ! longitude [rad]
real, intent(in), dimension(:,:,:)      :: pl                                           ! pressure at each half level [mbar]
real, intent(in), dimension(:,:)        :: drg                                          ! drag
integer, intent(in), dimension(:,:)     :: kpbl                                         ! level of pbl top
integer, intent(in)                     :: nh2o                                         ! index for water vapor
real, intent(in), dimension(:,:)        :: tg                                           ! ground temperature (equivalent to GT) [K]
real, intent(in)                        :: sols

real*8, intent(in),dimension(:,:,:,:)   :: r                                            ! tracer [kg/kg]
real*8, intent(in),dimension(:,:,:,:)   :: rdt                                          ! tracer tend [kg/kg/s]
real*8, intent(inout),dimension(:,:,:)  :: qpig                                         ! tracer on surface [kg/m2]
real*8, intent(out), dimension(size(pl,1),size(pl,2))                   :: wflux        ! sublimation flux [kg/m2/s]
real*8, intent(out), dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_subl     ! tendency [kg/kg/s]

!  Local variables
!  ---------------
real*8, dimension(size(r,1),size(r,2),size(r,3),size(r,4))  :: rini                     ! initial tracer [kg/kg]
real*8, dimension(size(pl,1),size(pl,2))                    :: qgnd,masspbl
real*8, dimension(size(pl,1),size(pl,2),size(pl,3)-1)       :: qpi
real*8, dimension(size(pl,1),size(pl,2))                    :: qpig_ini,wflux_subl
logical                                 :: mcpu0
integer                                 :: i,j,k,ndx,nt
character(len=15), dimension(4)         :: tagging_list

!  Treatment
!  ---------
mcpu0 = (mpp_pe() == mpp_root_pe())

!**********************************************************
!   sublimation/direct deposition of water from surface
!**********************************************************
! TB18q : no mixing in the atm ?
!  Note:  wflux > 0 is sublimation
!         wflux < 0 is condensation

masspbl(:,:) = 0.0
if (inj_vap_pbl) then
    do i = 1, size(pl,1)
        do j = 1, size(pl,2)
            do k = kpbl(i,j), nz
                masspbl(i,j) = masspbl(i,j) + (pl(i,j,k+1) - pl(i,j,k)) / grav
            enddo
        enddo
    enddo
else
    do i = 1, size(pl,1)
        do j = 1, size(pl,2)
            masspbl(i,j) = (pl(i,j,nz+1) - pl(i,j,nz)) / grav
        enddo
    enddo
endif

!! Updated tracer field and initialisations
rini(:,:,:,:)       = r(:,:,:,:) + rdt(:,:,:,:) * dt
rdt_subl(:,:,:,:)   = 0.
qgnd(:,:)           = (18.0/44.0) * 611.0 * exp(22.5 * (1.0 - (273.16 / tg(:,:))))      &
                    / pl(:,:,nz+1)

do nt = 1, nice_mass
    ndx        = vapor_indx(nt)
    qpi(:,:,:) = rini(:,:,:,ndx)
    wflux(:,:) = drg(:,:) * (qgnd(:,:) - qpi(:,:,nz))

    where (qpig(:,:,nt) .lt. 0.)
        qpig(:,:,nt) = 0.
    endwhere

    where (wflux(:,:) .gt. 0.)
        wflux(:,:) = wflux(:,:)*facsubl
    endwhere

    !! tag for sublimation
    wflux_subl(:,:) = max(wflux(:,:),0.)
#ifndef RELEASE
    tagging_list(1) = "geosource"
    tagging_list(2) = "loctime"
    tagging_list(3) = "antigeosource"
    tagging_list(4) = "cutsource"
    call tagging_main( tagging_list &
                     ,lat, lon, ndx, field2d=wflux_subl(:,:), solrun=sols)
#endif
    wflux(:,:)=min(wflux(:,:),wflux_subl(:,:))

    do i = is, ie
        do j = js, je
            ! sublimation : check we do not sublime the entire reservoir of ice
            if (wflux(i,j) .gt. 0.) then
                if (wflux(i,j)*dt .ge. qpig(i,j,nt)) then
                    wflux(i,j) = qpig(i,j,nt) / dt
                endif
                if (inj_vap_pbl) then
                    do k = kpbl(i,j), nz
                        qpi(i,j,k) = qpi(i,j,k) + dt * wflux(i,j) / masspbl(i,j)
                    enddo
                else
                    qpi(i,j,nz) = qpi(i,j,nz) + dt * wflux(i,j) * grav / (pl(i,j,nz+1) - pl(i,j,nz))
                endif
            endif

! condensation : check we do not reach negative values for vapor
            if (wflux(i,j) .lt. 0.) then
                qpi(i,j,nz) = qpi(i,j,nz) + dt * wflux(i,j) * grav / (pl(i,j,nz+1) - pl(i,j,nz))
                if (qpi(i,j,nz).lt.0.) then
                    wflux(i,j) = -rini(i,j,nz,ndx) * ((pl(i,j,nz+1) - pl(i,j,nz)) / grav) / dt
                    qpi(i,j,nz) = 0.
                endif
            endif
        enddo
    enddo

    qpig(:,:,nt)        = qpig(:,:,nt) - wflux(:,:) * dt
    rdt_subl(:,:,:,ndx) = (qpi(:,:,:) - rini(:,:,:,ndx)) / dt
enddo
return
end subroutine update_water


!#######################################################################
!#######################################################################
!#######################################################################


subroutine update_water_fv3(is, ie, js, je, nz, lat, lon, dt, pl, tg, drg, kpbl,    &
                            r, rdt, nh2o, qpig, wflux, rdt_subl)

! update the water vapor field for the bin microphysics scheme

use constants_mod, only: grav

implicit none

!  Arguments
!  ---------
real, intent(in)                            :: dt                                   ! physical time step [s]
integer, intent(in)                         :: nz,is,js,ie,je
real, intent(in), dimension(:,:)            :: lat                                  ! latitude [rad]
real, intent(in), dimension(:,:)            :: lon                                  ! longitude [rad]
real, intent(in), dimension(:,:,:)          :: pl                                   ! pressure at each half level [mbar]
real, intent(in), dimension(:,:)            :: drg                                  ! drag
integer, intent(in), dimension(:,:)         :: kpbl                                 ! level of pbl top
integer, intent(in)                         :: nh2o                                 ! index for water vapor
real, intent(in), dimension(:,:)            :: tg                                   ! ground temperature (equivalent to GT) [K]

!    Tracers :
real, intent(in),dimension(:,:,:,:)         :: r                                    ! tracer [kg/kg] bottom level
real, intent(in),dimension(:,:,:,:)         :: rdt                                  ! tracer tend [kg/kg/s] bottom level
real, intent(inout),dimension(:,:)          :: qpig                                 ! tracer on surface [kg/m2]
real, intent(out), dimension(size(pl,1),size(pl,2))                   :: wflux      ! sublimation flux [kg/m2/s]
real, intent(out), dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_subl   ! tendency [kg/kg/s]

!  Local variables
!  ---------------
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4))              :: rini       ! initial tracer [kg/kg]
real, dimension(size(pl,1),size(pl,2))                                :: qgnd, masspbl
real, dimension(size(pl,1),size(pl,2),size(pl,3)-1)                   :: qpi
real, dimension(size(pl,1),size(pl,2))                                :: qpig_ini
logical                                     :: mcpu0
integer                                     :: i, j, k, ndx, nt

!  Treatment
!  ---------
mcpu0 = (mpp_pe() == mpp_root_pe())

!**********************************************************
!   sublimation/direct deposition of water from surface
!**********************************************************
! TB18q : no mixing in the atm ?
!  Note:  wflux > 0 is sublimation
!         wflux < 0 is condensation

masspbl(:,:) = 0.0
if (inj_vap_pbl) then
    do i = 1, size(pl,1)
        do j = 1, size(pl,2)
            do k = kpbl(i,j), nz
                masspbl(i,j) = masspbl(i,j) + (pl(i,j,k+1) - pl(i,j,k)) / grav
            enddo
        enddo
    enddo

else
    do i = 1, size(pl,1)
        do j = 1, size(pl,2)
            masspbl(i,j) = (pl(i,j,nz+1) - pl(i,j,nz)) / grav
        enddo
    enddo
endif

rini(:,:,:,:)       = r(:,:,:,:) + rdt(:,:,:,:) * dt
qpi(:,:,:)          = rini(:,:,:,nh2o)
rdt_subl(:,:,:,:)   = 0.
wflux(:,:)          = 0.
qpig_ini(:,:)       = qpig(:,:)

do i = is, ie
    do j = js, je
        if (qpig(i,j) .ge. 0.) then
            qgnd(i,j)  = (18.0 / 44.0) * 611.0 * exp(22.5 * (1.0 - (273.16 / tg(i,j)))) / pl(i,j,nz+1)
            wflux(i,j) = facsubl * drg(i,j) * (qgnd(i,j) - qpi(i,j,nz))
            ! Sublimation : check we do not sublime the entire reservoir of ice
            if (wflux(i,j) .gt. 0.) then
                if (wflux(i,j) * dt .ge. qpig(i,j)) then
                    wflux(i,j) = qpig(i,j) / dt
                endif
                if (inj_vap_pbl) then
                    do k = kpbl(i,j), nz
                        qpi(i,j,k) = qpi(i,j,k) + dt * wflux(i,j) / masspbl(i,j)
                    enddo
                else
                    qpi(i,j,nz) = qpi(i,j,nz) + dt * wflux(i,j) * grav / (pl(i,j,nz+1) - pl(i,j,nz))
                endif
            endif
            ! condensation : check we do not reach negative values for vapor
            if (wflux(i,j) .lt. 0.) then
                qpi(i,j,nz) = qpi(i,j,nz) + dt * wflux(i,j) * grav / (pl(i,j,nz+1) - pl(i,j,nz))
                if (qpi(i,j,nz).lt.0.) then
                    wflux(i,j)  = -rini(i,j,nz,nh2o) * ((pl(i,j,nz+1) - pl(i,j,nz)) / grav) / dt
                    qpi(i,j,nz) = 0.
                endif
            endif
        endif
    enddo
enddo

qpig(:,:)            = qpig(:,:) - wflux(:,:) * dt
rdt_subl(:,:,:,nh2o) = (qpi(:,:,:) - rini(:,:,:,nh2o)) / dt

return
end subroutine update_water_fv3


!#######################################################################
!#######################################################################
!#######################################################################


subroutine h2o2cond(is,js,nz,lat,lon,dt,phalf,pfull,temp,nh2o2,cond_mass,r,rdt,rdt_h2o2)

use constants_mod, only: GRAV
use mars_surface_mod, only: sfc_h2o2_chem

implicit none

integer, intent(in)                     :: is,js,nz
integer, intent(in)                     :: nh2o2        ! index for h2o2
real, intent(in)                        :: dt           ! physical time step (s)
real, intent(in),    dimension(:,:)     :: lat
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:,:)   :: phalf        ! pressure[Pa]
real, intent(in),    dimension(:,:,:)   :: pfull        ! pressure[Pa]
real, intent(in),    dimension(:,:,:)   :: temp         ! temp [K]
real, intent(in),    dimension(:,:,:,:) :: r            ! [mmr]
real, intent(in),    dimension(:,:,:,:) :: rdt          ! tendancy updated [kg/kg/s]
real, intent(out),   dimension(:,:,:)   :: rdt_h2o2     ! H2O2 tendancy updated [kg/kg/s]
real, intent(out),   dimension(:,:,:)   :: cond_mass    ! condensing mass (kg m-2)
!real, intent(inout), dimension(:,:)     :: sfc_h2o2     ! H2O2 only tracer on surface (kg/m2)

! Local variables
logical                                 :: mcpu0 !debug only, identify master processor
integer                                 :: L,j,i,k,kd,iter,id,jd,ie,je
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rini            ! tendency kg/kg
real, dimension(size(temp,1),size(temp,2),size(temp,3))  :: q_h2o2          ! H2O2 mass mixing ratio [kg/kg]
real, dimension(size(temp,1),size(temp,2))               :: sfc_h2o2_ini    ! H2O2 only tracer on surface (kg/m2)
real, dimension(size(temp,1),size(temp,2),size(temp,3))  :: laymass         ! mass in a layer [kg m-2]
real, dimension(size(temp,1),size(temp,2),size(temp,3))  :: pvp_h2o2        ! H2O2 Partial vapor pressure (mmHg -> Pa)
real, dimension(size(temp,1),size(temp,2),size(temp,3))  :: qsat_h2o2       ! H2O2 mass mixing ratio at saturation (kg/kg of air)
! real, dimension(size(temp,3))                           :: cond_mass       ! condensing mass (kg m-2)
real, dimension(size(temp,1),size(temp,2),size(temp,3))  :: qcond_h2o2      ! mmr to be condensed: difference between mass mixing ratio and saturated mmr (kg/kg)
!---------------

id = size(temp,1); jd = size(temp,2)
ie = is + id - 1
je = js + jd - 1

!======================================================================!
! Purpose:
!     This routine is to calculate H2O2 staturation ratio and determine when
!     H2O2 condenses out to the surface.
!
! ** Note: pvp_h2o2 from Lindner, 1988
! ** Note: L=1 => top of atm and L=nz => surface.
!
!======================================================================!

mcpu0 = (mpp_pe() == mpp_root_pe())
!! Updated tracer field and initialisations
rini(:,:,:,:)       = r(:,:,:,:) + rdt(:,:,:,:) * dt
q_h2o2(:,:,:)       = rini(:,:,:,nh2o2)
rdt_h2o2(:,:,:)     = 0.
sfc_h2o2_ini(:,:)   = sfc_h2o2_chem(:,:)

pvp_h2o2(:,:,:)     = 0.
qcond_h2o2(:,:,:)   = 0.
qsat_h2o2(:,:,:)    = 1.

do k = 1, nz
    WHERE (temp(:,:,k) < 220)
        pvp_h2o2(:,:,k)  = (10**(11.98-3422 / temp(:,:,k))) * 133.322       ! [mmHg] -> [Pa]
        qsat_h2o2(:,:,k) = (pvp_h2o2(:,:,k) / pfull(:,:,k)) * (34.01 / 44.01)
    END WHERE
enddo ! k-loop

do k = 1, nz-1
    WHERE (q_h2o2(:,:,k) >= qsat_h2o2(:,:,k))
        laymass(:,:,k)    = (phalf(:,:,k+1) - phalf(:,:,k)) / grav
        qcond_h2o2(:,:,k) = q_h2o2(:,:,k) - qsat_h2o2(:,:,k)
        cond_mass(:,:,k)  = qcond_h2o2(:,:,k) * laymass(:,:,k)
     ! Update tracer
        q_h2o2(:,:,k)     = q_h2o2(:,:,k) - qcond_h2o2(:,:,k)
        q_h2o2(:,:,k+1)   =q_h2o2(:,:,k+1) + qcond_h2o2(:,:,k)
    END WHERE
end do ! k-loop

WHERE (q_h2o2(:,:,nz) >= qsat_h2o2(:,:,nz))
    laymass(:,:,nz)     = (phalf(:,:,nz+1) - phalf(:,:,nz)) / grav
    qcond_h2o2(:,:,nz)  = q_h2o2(:,:,nz) - qsat_h2o2(:,:,nz)
    cond_mass(:,:,nz)   = qcond_h2o2(:,:,nz) * laymass(:,:,nz)
    q_h2o2(:,:,nz)      = q_h2o2(:,:,nz) - qcond_h2o2(:,:,nz)
  ! Update surface reservoir
    sfc_h2o2_chem(:,:)  = sfc_h2o2_chem(:,:) + cond_mass(:,:,nz)
END WHERE

! Update surface reservoir tendency
! sfc_h2o2_chem(:,:) = sfc_h2o2_ini(:,:) + cond_mass(nz)

! Update tracer tendency
rdt_h2o2(:,:,:) = (q_h2o2(:,:,:) - rini(:,:,:,nh2o2)) / dt

return
end


!#######################################################################
!#######################################################################
!#######################################################################


end module mars_physics_mod
