module co2micromod_mgcm
!
!  module for bulk co2 cloud microphysics
!
use constants_mod, only: grav,cp_air,rdgas,pi,kbz=>kboltz
use initracer_mod

use   mpp_mod, only: input_nml_file
use           fms_mod, only: error_mesg, FATAL,       &
                             check_nml_error, &
                             mpp_pe, mpp_root_pe, &
                             write_version_number, stdlog,        &
                             uppercase

use       fms2_io_mod, only:  file_exists, FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                                   register_restart_field, register_axis, unlimited, &
                                   open_file, read_restart, write_restart, close_file, &
                                   register_field, read_data, write_data, register_variable_attribute, &
                                   get_global_io_domain_indices, get_variable_size, variable_exists

use time_manager_mod, only: time_type, get_time
use diag_manager_mod, only: register_diag_field, send_data,  diag_axis_init
use rtmod_mgcm, only: dtridgl   ! for call to dtridgl TB18q
use field_manager_mod, only  : MODEL_ATMOS, find_field_index 

implicit none

public :: co2micro_driver,co2micro_driver_init

real*8  ::  ccn_co2=1.e5                ! Number of co2 cloud condensation nuclei #/kg of gaseous CO2 for sedimentation


namelist /co2microphys_nml/ ccn_co2

integer ::  id_co2cld, id_co2cld_col, id_co2cld_r, id_co2_sed_dt, id_co2cld_gen, id_co2_sed_v

real, dimension(:,:),   allocatable, save  ::  co2col


!-------------------- Other -----------------------------------------
logical ::  mcpu0
logical,save :: firstcall_micro=.true.

!! effective T and Q from previous timestep for cloud scheme
real*8, dimension(:,:,:),     allocatable,save :: tlprev
real*8, dimension(:,:,:),     allocatable,save :: qpiprev,qpi_dmprev
real*8, parameter   :: missing_value = -1.e10
character(len=12) :: mod_name = 'co2cloudphys'

contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine co2micro_driver(is,ie,js,je,kd,ntp,Time,p_half,p_full,t,tdt, &
         r,rdt,strss,tg,dtime,drg,rdt_co2micro,precip,dmass,rkh,tdt_co2,checkcons)

!
!  Driver for co2 microphysics
!
use initracer_mod
use constants_mod, only: grav,rdgas
implicit none
!     Vertical structure is assumed as follows:
!     Midpoints are assumed for arrays with dimension NZ (top at index = 1)
!     Midpoints and boundaries are assumed for arrays with 2*NZ+X
!     - Midpoints are at even indices, and boundaries at odd indices
!     - The top boundary of the model is assumed to be at index = 3
!
!     ---------------------    2*L + 1
!     |                   |
!     |                   |    2*L + 2   or L
!     |                   | 
!     ---------------------    2*L + 3
!===============================================================
!     Input Arguments 
!===============================================================
integer, intent(in)  :: is, js, ie, je, kd, ntp
type(time_type), intent(in)             :: Time !     Model time
real*8, intent(in), dimension(is:ie,js:je,kd+1) :: p_half  !     layer boundary pressures [Pa]
real*8, intent(in), dimension(is:ie,js:je,kd) :: p_full  !     layer midpoint pressures [Pa]
real*8, intent(in), dimension(is:ie,js:je,kd) :: t       !     Temperature [K]
real*8, intent(in), dimension(is:ie,js:je,kd,ntp) :: r     !     Tracer field
real*8, intent(in), dimension(is:ie,js:je) :: strss     !     Surface Stress
real*8, intent(in), dimension(is:ie,js:je) :: drg       !     Drag
real*8, intent(in), dimension(is:ie,js:je,2*kd+1) :: rkh     !     PBL Eddy coefficient for sedimentation

real*8, intent(in), dimension(is:ie,js:je) :: tg        !     Ground Temperature (K)
real*8, intent(in) :: dtime                     !     Time step
real*8, intent(in), dimension(is:ie,js:je,kd,ntp) :: rdt   !     Tracer tendencies (kg/kg/s)
real*8, intent(in), dimension(is:ie,js:je,kd) :: tdt     !     Temperature tendencies (K/s)
logical, intent(in) :: checkcons                !     check conservation
!===============================================================
!     Output Arguments
!===============================================================
real*8, intent(out), dimension(is:ie,js:je,kd,ntp) :: rdt_co2micro
real*8, intent(out), dimension(is:ie,js:je,kd) :: tdt_co2
real*8, intent(out), dimension(is:ie,js:je)         ::  precip      ! CO2 snow accumulation this time step [kg/m2]
real*8, intent(out), dimension(is:ie,js:je,kd)       ::  dmass       ! The array dmass may be optionally used to modify the atmospheric mass [kg/m2]
!===============================================================
!===============================================================
!     Local Variables To Microphys
real*8, dimension(size(t,1),size(t,2)) :: taucld, ccldcol
integer n, nco2, nt, ndx
integer  :: id, jd, i, j, k, l,ilay
real*8 Rn,Rs,cst,cst2,fact
logical :: used
!     Number of layer midpoints
integer :: nz, lbot

!     **********************
!     Atmospheric variable on regular vertical grid
!     **********************
real*8, dimension(size(t,1),size(t,2)) :: deposit  ! tendency Tracer field on surface, co2 only
real*8, dimension(size(r,1),size(r,2),ntrace_mom) :: qpig  ! tracer on surface

real*8, dimension(size(t,1),size(t,2),size(t,3)) :: rho ! Atmospheric density 
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: rco2 ! Particle radius (m) 
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: tini ! Updated Atmospheric temperature
real*8, dimension(size(t,3),ntrace_mom) :: ratio_mass,ratio_nb,ratio_core, ratio_vap  ! Ratio for tagging
!     Tracers MMR (kg/kg) 
real*8, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rini   ! Updated tracer field

real*8, dimension(size(r,1),size(r,2),size(r,3)) :: rfin 


!     **********************
!     Atmospheric variable on new vertical grid
!     ********************** 
!     Temperature (K) at boundaries and midpoints
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3) :: tl,tlsub,tl0,tleff,tl_ini
!     Tendencies for Temperature (K) at boundaries and midpoints
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3) :: tlsubdt,tldteff
!     Pressure (Pa) at boundaries and midpoints
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3) :: pl
! Tracer field on specific new vertical grid
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3) :: qpi,qpisub,qpi0,qpieff,qpi_ini
!     Tendancies temperature (clouds scheme)
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: delp

! constants
real*8 aerdensco2 
real*8 dev
real*8 No_co2,Mo
real*8 g,xlhtc,cp2


! *******************************************************************************     
! ************************ Initializations **************************************     
! *******************************************************************************     
mcpu0 = (mpp_pe() == mpp_root_pe())
id= size(tl,1); jd= size(tl,2)

! Get the index of the CO2 cloud tracer
nco2= find_field_index( MODEL_ATMOS, 'co2_cloud' )

! Constants 
aerdensco2=1620.
No_co2 = ccn_co2   ! Number of cloud condensation nuclei (# / kg of gaseous CO2)
xlhtc= 5.902e+5
cp2 = 7.3594D+2
dev = 0.3087     ! gives an effective variance of 0.1  for co2 ice

! Original Vertical Levels
nz = size(t,3)

! Fields initialized to 0
tdt_co2(:,:,:) = 0.0
rdt_co2micro(:,:,:,:) = 0.0

qpi0(:,:,:) = 0.0
tl0(:,:,:) = 0.0
qpi(:,:,:) = 0.0
qpi_ini(:,:,:) = 0.0
tl(:,:,:) = 0.0
pl(:,:,:) = 0.0

qpieff(:,:,:) = 0.0
tleff(:,:,:) = 0.0
pl(:,:,:) = 0.0

rfin = 0.d0
precip(:,:) = 0.d0
deposit(:,:) = 0.d0
dmass(:,:,:) = 0.d0


! *******************************************************************************     
! ************************ Atm P, T, Q  *****************************************
! *******************************************************************************     
!! Pressure on new vertical grid : pl
!! Used for clouds microphysics :
!!   Updated temperature on original vertical grid : tini  
!!   Updated temperature on new vertical grid : tl 
!!   Initial temperature on original vertical grid : t  
!!   Initial temperature on new vertical grid : tl0
!!   Temperature tendency of all previous processes on the original vertical grid : tdt

!! Updated tracers on original vertical grid : rini
!! Updated tracers on new vertical grid : qpi
!! Initial tracers on origial vertical grid : r
!! Initial tracers on new vertical grid : qpi0


! k is levels (midpoints and boundaries)
! n is layers (midpoints)

!! Pressure at the model top
pl(:,:,3) = p_half(:,:,1)
pl(:,:,2) = pl(:,:,3)*0.5
pl(:,:,1) = pl(:,:,2)*1.e-6

!! Update field t => tini
tini(:,:,:)=t(:,:,:)+tdt(:,:,:)*dtime

!! get tl0, tl and pl
do n= 1, nz
    k= 2*n + 2
    tl(:,:,k)= tini(:,:,n)
    tl0(:,:,k)= t(:,:,n)
    pl(:,:,k)= p_full(:,:,n)
    pl(:,:,k+1)=p_half(:,:,n+1)
    if (n.eq.nz) then
        tl(:,:,k+1) = 0.5*( tg(:,:)+tini(:,:,n) )
        tl0(:,:,k+1) = 0.5*( tg(:,:)+t(:,:,n) )
    else
        tl(:,:,k+1)= 0.5*( tini(:,:,n+1)+tini(:,:,n) )
        tl0(:,:,k+1)= 0.5*( t(:,:,n+1)+t(:,:,n) )
    end if
enddo
tl(:,:,3) = tl(:,:,4)
tl0(:,:,3) = tl0(:,:,4)

!! density
do n = 1, nz
    k   = n * 2 + 2
    rho(:,:,n) = pl(:,:,k) / ( rdgas*tl(:,:,k) )
enddo

! Update tracer field
rini(:,:,:,nco2)=r(:,:,:,nco2)+(rdt(:,:,:,nco2)*dtime)


! Get tracer fields on new vertical grid : qpi, qpi0
do n=1,nz ! levels 
    k = 2*n
    ndx = nco2
    qpi(:,:,k+2) = max(rini(:,:,n,ndx),0.)
    qpi0(:,:,k+2) = max(r(:,:,n,ndx),0.)
end do

!Delp
do n=1,nz
    delp(:,:,n)= p_half(:,:,n+1)-p_half(:,:,n)
enddo

!! Saving current T and Q field before cloud scheme
qpi_ini(:,:,:)=qpi(:,:,:)
tl_ini(:,:,:)=tl(:,:,:)

!! Tracers for cloud scheme
qpisub(:,:,:)=qpi(:,:,:)
tlsub(:,:,:)=tl(:,:,:)

! Loop through I and J
do i=is,ie
    do j=js,je

        ! do condensation and sublimation
        call co2condsub(dtime,pl(i,j,:),tl(i,j,:),qpisub(i,j,:),xlhtc,cp2,nz)
     
        do n = 1, nz
          k = n * 2 + 2

          ! dmass atmospheric mass change from condensation/sublimation ! kg/m2
          dmass(i,j,n)= (qpisub(i,j,k)-qpi_ini(i,j,k)) * delp(i,j,n) / grav

          rho(i,j,n) = pl(i,j,k) / ( rdgas*tl(i,j,k) )
          Mo = qpisub(i,j,k)
          rco2(i,j,n) = ( Mo / No_co2 * 0.75 / pi / aerdensco2 )**(athird) * dexp( -0.5*dev**2. )

          rco2(i,j,n)=max(rco2(i,j,n),1.e-7)

          rco2(i,j,n)=rco2(i,j,n)*dexp ( 1.5 * dev**2. )
        enddo

        ! sedimentation of CO2 
        call sedimco2(dtime,pl(i,j,:),tl(i,j,:),rho(i,j,:),rkh(i,j,:),qpisub(i,j,:),rco2(i,j,:),aerdensco2,deposit(i,j),nz)
 
    enddo    ! j
enddo     ! i


!! Get final tendencies 
do n=1,nz
    k = 2*n
    rdt_co2micro(:,:,n,nco2)=(qpisub(:,:,k+2)-qpi_ini(:,:,k+2))/dtime
    tdt_co2(:,:,n)=(tl(:,:,k+2)-tl_ini(:,:,k+2))/dtime
end do

!Fill precip array with co2 snow that falls to surface via sedimentation 
precip(:,:) = precip(:,:) + deposit(:,:) 

fact= 1.0/grav

co2col(is:ie,js:je)= 0.0

! check tendency
do n= 1, nz
    co2col(is:ie,js:je)= co2col(is:ie,js:je) + (r(:,:,n,nco2)+dtime*(rdt(:,:,n,nco2)+rdt_co2micro(:,:,n,nco2)))*delp(:,:,n)*fact
enddo

rfin = r(:,:,:,nco2)+dtime*(rdt(:,:,:,nco2)+rdt_co2micro(:,:,:,nco2))

if (id_co2cld > 0) used = send_data ( id_co2cld, rfin, time, is, js )
if (id_co2cld_col > 0) used = send_data ( id_co2cld_col, co2col(is:ie,js:je), time, is, js )
if (id_co2cld_r > 0)  used = send_data ( id_co2cld_r, rco2, time, is, js )


return
end subroutine co2micro_driver

!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine co2condsub(dt,plc,tlc,qpic,xlhtc,cp2,nz)
!                                                              !
!     This routine updates species concentrations due          !
!     to both nucleation and condensation-induced variations.  !
!     Gain and loss rates associated to each one of these      !
!     processes are computed separately in other routines.     !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Arguments
!  ---------

      integer nz

      real*8 dt                    ! dtime time step
      real*8 plc(2*nz+3),tlc(2*nz+3) ! pressure and temperature in the column
      real*8 qpic(2*nz+3)           !
      real*8 cp2,xlhtc

!  Local
!  -----

      integer i,l

      real*8   Cste
      real*8 p,temp

      real*8 qcond,qpisav_co2,tsat,psat,tlsav,qpinew_co2

      integer ilay

      real*8   derf

      real*8 dqpi

! Treatment
! Start loop over heights


      DO 100 l = 1, nz

        ilay = 2 * l + 2

        qpisav_co2 = qpic(ilay)

        tlsav=tlc(ilay)

        psat=plc(ilay) / 100.  ! convert to mbar
        tsat=3182.48/(23.3494-log(psat))
        qcond=0.

        if (qpisav_co2 .gt. 0.) then

          tlc(ilay)=tlc(ilay)-qpisav_co2*xlhtc/cp2

          if (tlc(ilay) .lt. tsat) then
            qcond = cp2*(tsat-tlc(ilay)) /xlhtc
            tlc(ilay) = tsat
            qpinew_co2 = qcond
          else
            qcond = (-1.)*qpisav_co2
            qpinew_co2 = 0.
          endif

        else

          if (tlc(ilay) .lt. tsat) then

            qcond = cp2*(tsat-tlc(ilay))/xlhtc
            tlc(ilay) = tsat
            qpinew_co2 = qcond

          else

            qcond = 0.
            qpinew_co2 = 0.

          endif

        endif

        qpic(ilay) = qpinew_co2

100   CONTINUE


      return

end subroutine co2condsub


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sedimco2(dt,pls,tls,rhos,kd,qpis,rco2s,aerdensco2,deposit_co2,nz)
!                                                              !
!                Computing aerosol sedimentation               !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!      use constants_h, only: GRAV, RGAS
!                                                              !
!      use comp3cmn_h, only: jcmn,icmn
!      use standard_h, only: srfupflx,srfdnflx,dsig,p

      implicit none

! Arguments
! ---------

      integer nz           ! dimension of the c array

      real*8 dt              ! time step
      real*8 pls(2*nz+3),tls(2*nz+3)
      real*8 rhos(nz)         ! atmospheric density
      real*8 kd(2*nz+1)        ! Eddy mixing coefficient (m2/s)
      real*8 qpis(2*nz+3)! Integrated concentrations of particle/each bin (#/m3)
      real*8 qpi2s(2*nz+3)
      real*8 rco2s(nz)

! Local variables
! ---------------

      real*8 rhob(nz)       ! atmospheric density at layer boundaries
      real*8 dz(nz)         ! layer thickness (m)
      real*8 vf             ! fall velocity of particle (m/s)
      real*8 cour,w1        ! Courant number & corrected velocity (vf+kd/H with h scale height)
      real*8 dzbx
      real*8 sigma,theta,hc,lg,rap,cmp,w,wp
      real*8 fs(nz+1),ft(nz+1)
      real*8 as(nz),bs(nz),cs(nz),ds(nz)
      real*8 asi(nz),bsi(nz),csi(nz),dsi(nz),xsol(nz)
      real*8 cold(2),deposit_co2
      real*8 c(nz),aerdensco2,q1,q2

      integer i,l,n
      integer ilay,ilev


      theta = 0.0 


!     Layer thickness: dz
      do l = 1, nz
        dz(l) = ( pls(2*l+3)-pls(2*l+1) ) / grav / rhos(l)
      enddo

!     Compute density at the layer boundaries
      do l = 1, nz
        ilev = l * 2 + 3
        if (l.lt.nz) then
          rhob(l) = pls(ilev)  / (rdgas*tls(ilev))
        else
          rhob(l) = pls(ilev-1) / (rdgas*tls(ilev-1))
        endif
      enddo

! Loop over layers

      do 20 l = 1, nz

      ilay = 2 * l + 2

      c(l) = qpis(ilay) * rhos(l)

      if (l.eq.1) goto 20

!     Compute fall velocity
      ilev = 2 * l + 3
      call fallvel(scale_dt,aerdensco2,rco2s(l),tls(ilev),rhob(l),vf)

      dzbX = ( dz(l)+dz(l-1) ) / 2.

      w  = -1. * vf !* exp(-stdv(3)**2.)

!     Get the corrected fall velocity (virtual speed accounting for mixing)

      if (kd(2*l-1) .ne. 0.) then
        theta = 0.5 * ( w*dzbX/kd(2*l-1) + log(rhos(l-1)/rhos(l)) )
        if (theta.ne.0) then
          sigma = 1./dtanh(theta) - 1./theta
        else
          sigma = 1.
        endif
      else
        sigma = 1.
      endif

      if (c(l).eq.0.) then
        rap=10.
        if (c(l-1).eq.0.) then
          rap=1.
        endif
      else
        rap = min( max(c(l-1)/c(l),0.1), 10.)
      endif

      cour=abs(w*dt)

      if (rap.gt.0.9 .and. rap.lt.1.1 .or. cour.gt.dz(l)) then
        w1 = w
      else
        if (w.lt.0) then
          hc = dzbX / dlog(rap)
          lg = dzbX / (w*dt) * (dexp(-w*dt/hc)-1.) / (1.-rap)
          wp = w * 1.d0
          cmp= dlog(-wp) + abs(sigma) * dlog(lg)
          w1 = -dexp(cmp)
        else
          w1 = 0.
        endif
      endif

!  Fluxes at layer boundaries

      if (kd(2*l-1).ne.0.) then
        if (theta.ne.0.) then
          ft(l)=( w1 + log(rhos(l-1)/rhos(l))*kd(2*l-1)/dzbX ) / ( dexp(2.*theta) - 1. )
          fs(l) = ft(l) * dexp(2.*theta)
        else
          ft(l) = kd(2*l-1) / dzbX
          fs(l) = kd(2*l-1) / dzbX
        endif
      else
        if (w1.lt.0.)then
          ft(l) = -w1
          fs(l) = 0.
        else
          ft(l) = 0.
          fs(l) = w1
        endif
      endif

20    continue

! Boundary conditions for the fluxes

      fs(1)    =  0.
      ft(1)    =  0.
      fs(nz+1) =  0.
      ft(nz+1) = -w1

! Compute the coefficient of the continuity equation

      do l=1,nz
        cs(l) =  ft(l+1) + fs(l) - dz(l) / dt
        if ( cs(l) .gt. 0. ) goto 1010
        as(l) = -dz(l) / dt
        bs(l) = -ft(l)
        ds(l) = -fs(l+1)
      enddo

! Depending on the cs value, switch to an explicit or an implicit scheme

! Explicit case 

      cold(1)  = c(1)
      c(1) = ( cs(1)*c(1) + ds(1)*c(2) ) / as(1)

      do l = 2, nz-1
        cold(2)  = c(l)
        c(l) = ( bs(l)*cold(1) + cs(l)*c(l) + ds(l)*c(l+1) ) / as(l)
        cold(1)  = cold(2)
      enddo

! Compute the mass of co2 ice falling on the ground
      deposit_co2 = deposit_co2 + c(nz) * ft(nz+1) * dt


      c(nz) = ( bs(nz)*cold(1) + cs(nz)*c(nz) ) / as(nz)

      do l = 1, nz
        qpis(2*l+2) = c(l) / rhos(l)
      enddo

      GOTO 111

1010  continue

! Implicit case 

      do l = 1, nz
        asi(l) =  ft(l)
        bsi(l) = -( ft(l+1) + fs(l) + dz(l)/dt )
        csi(l) =  fs(l+1)
        dsi(l) = -dz(l) / dt * c(l)
      enddo

! Matrix inversion

      call dtridgl(nz,asi,bsi,csi,dsi,xsol)

      do l = 1, nz
        c(l) = xsol(l)
        qpis(2*l+2) = c(l) / rhos(l)
      enddo

! Compute the mass of water ice falling on the ground
      deposit_co2 = deposit_co2 + c(nz) * ft(nz+1) * dt


111   CONTINUE


      RETURN


end subroutine sedimco2


!**************************************************************************
!**************************************************************************

subroutine fallvel(scale,dpden,r,t,rho,vf)
!  Calculate the fall velocity of dust particles in the Martian
!  atmosphere at the model levels.  Part of the dust tracer scheme for 
!  the c-grid model.
!
!  VARIABLES:
!
!  WT        - mean thermal velocity
!  MFP       - mean free path
!  DV        - dynamic viscosity (kg m^-1 s^-1)
!  KN        - Knudsen number
!  ALPHA     - ALPHA: from (1 + ALPHA*Kn) which is the Cunningham
!              slip-flow correction
!  CONST     - the level-independent part of the gravitational settling
!              velocity:  (2*dpden*GRAV*r^2)/9
!
!
use constants_mod, only: GRAV

implicit none

real*8, intent(in) ::   t, &        ! temperature [K]
                        rho, &      ! density [kg/m^3]
                        r, &        ! particle radius [m]
                        scale, &    ! 3.0*KBAMU/44.0 (44.0 is mean molecular weight of atm)
                        dpden       ! particle density
real*8, intent(out) :: vf           ! fall velocity [m/s]
! Local variables
real*8  wt, mfp, dv, kn, alpha
real*8  const

!C======================================================================C

const = 2.0*dpden*r*r*grav/9.0

wt    = sqrt(scale*t)
dv    = (1.59E-6*(t**1.5))/(t+244.4)
mfp   = 2.0*dv/(rho*wt)
kn    = mfp/r
alpha = 1.246 + 0.42*exp(-0.87/kn)
vf = const*(1.0+alpha*kn)/dv

return
end subroutine fallvel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine co2micro_driver_init( nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )
!
!  Initialize microphysics driver
!

!! Arguments
integer, intent(in)                   :: nlon, mlat, nlevels
real,    intent(in),  dimension(:,:)  :: lonb, latb
real,    intent(in),  dimension(:,:)  :: lon, lat
integer, intent(in)                   :: axes(4)
type(time_type), intent(in) :: Time

!! Local
integer  unit, io, ierr
integer  id, jd, km, i, j, k, is, js, ie, je, nt, ndx,n
character (len=128) :: filename, fieldname, tracer_name, tname

is= 1
js= 1
id= size(lon,1)
jd= size(lat,2)
ie= is + id - 1
je= js + jd - 1



!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
read (input_nml_file, nml=co2microphys_nml, iostat=io)
ierr = check_nml_error(io,'co2microphys_nml')

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=co2microphys_nml)

! *********************************************************
!   --- Allocate T + Q fields needed for cloud scheme ---
! *********************************************************
allocate (  tlprev(is:ie,js:je,2*nlevels+3)  )
allocate (  qpiprev(is:ie,js:je,2*nlevels+3)  )
allocate (  qpi_dmprev(is:ie,js:je,2*nlevels+3)  )


allocate (  co2col(is:ie,js:je)  )

! *********************************************************
!     ----- register diagnostic fields -----
! *********************************************************


id_co2cld = register_diag_field ( mod_name, 'co2cld',  &
                                 (/axes(1:3)/), Time,           &
                                'co2 cld ', '',        &
                                 missing_value=missing_value )

id_co2cld_col = register_diag_field ( mod_name, 'co2cld_col',  &
                                 (/axes(1:2)/), Time,           &
                                'co2 cld column ', '',        &
                                 missing_value=missing_value )

id_co2cld_r = register_diag_field ( mod_name, 'co2cld_rad',  &
                                 (/axes(1:3)/), Time,           &
                                'co2 ice particle radius ', '',        &
                                 missing_value=missing_value )

id_co2_sed_dt = register_diag_field ( mod_name, 'co2_sed_dt',  &
                                 (/axes(1:3)/), Time,           &
                                'co2 cloud sedimentation tendency', '',        &
                                 missing_value=missing_value )

id_co2cld_gen = register_diag_field ( mod_name, 'dco2cld_dt',  &
                                 (/axes(1:3)/), Time,           &
                                'co2 cloud microphysics tendency', '',        &
                                 missing_value=missing_value )

id_co2_sed_v = register_diag_field ( mod_name, 'co2sedvel',  &
                                 (/axes(1:3)/), Time,           &
                                'co2 cloud sedimentation velocity', '',        &
                                 missing_value=missing_value )



end subroutine co2micro_driver_init





end module co2micromod_mgcm
