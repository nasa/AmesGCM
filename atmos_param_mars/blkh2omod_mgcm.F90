module blkh2omod_mgcm
!
! module for bulk water clouds
!
use constants_mod, only: grav,cp_air,rdgas,pi,kboltz

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

use   mpp_domains_mod, only: domain2d
use mpp_mod, only: input_nml_file

use time_manager_mod, only: time_type
use diag_manager_mod, only: register_diag_field, send_data
use rtmod_mgcm, only: dtridgl   ! for call to dtridgl TB18q
use mars_surface_mod,  only: sfc_frost_blk, cumulative_prec_blk
use field_manager_mod, only  : MODEL_ATMOS, find_field_index 

implicit none

public :: blkh2o_driver, blkh2o_driver_init

real*8  ::  ccn_blkh2o=1.e5                ! Number of cloud condensation nuclei [#/kg of gaseous CO2]
real*8  ::  prec_threshold=0.001           ! Cloud mass mixing ratio (kg/kg) for precipitation threshold

real*8, parameter :: aerdensh2o=917.0      ! density of ice [kg/m^3]
real*8, parameter :: lw  = 2.8e+6          ! latent heat of vaporization [J/kg]

namelist /blkh2oclouds_nml/ ccn_blkh2o,prec_threshold

integer ::  id_blkh2ocld, id_blkh2ocld_col, id_blkh2ocld_r, id_blkh2o_sed_dt, id_blkh2ocld_gen, id_blkh2o_sed_v
integer ::  id_cprecip, id_cprecip_rain, id_cprecip_snow

!-------------------- Other -----------------------------------------
logical ::  mcpu0
logical,save :: firstcall_micro=.true.

!! effective T and Q from previous timestep for cloud scheme
real*8, dimension(:,:,:),     allocatable,save :: tlprev
real*8, dimension(:,:,:),     allocatable,save :: qpiprev,qpi_dmprev
real*8, parameter   :: missing_value = -1.e10
character(len=15) :: mod_name = 'h2oblkcloudphys'

contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine blkh2o_driver(is,js,ie,je,kd,ntp,Time,p_half,p_full,t,tdt, &
         r,rdt,tg,dtime,rdt_blkh2o,rkh,tdt_blkh2o,checkcons)
!
! Main driver for bulk h2o clouds
!

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
real*8, intent(in), dimension(is:ie,js:je,2*kd+1) :: rkh     !     PBL Eddy coefficient for sedimentation

real*8, intent(in), dimension(is:ie,js:je) :: tg        !     Ground Temperature (K)
real*8, intent(in) :: dtime                     !     Time step
real*8, intent(in), dimension(is:ie,js:je,kd,ntp) :: rdt   !     Tracer tendencies (kg/kg/s)
real*8, intent(in), dimension(is:ie,js:je,kd) :: tdt     !     Temperature tendencies (K/s)
logical, intent(in) :: checkcons                !     check conservation
!===============================================================
!     Output Arguments
!===============================================================
real*8, intent(out), dimension(is:ie,js:je,kd,ntp) :: rdt_blkh2o
real*8, intent(out), dimension(is:ie,js:je,kd) :: tdt_blkh2o
!===============================================================
!===============================================================
!     Local Variables To Microphys
!     Water/cloud column
real*8, dimension(size(t,1),size(t,2)) :: taucld 
integer n, nt, ndx, microstep !ndx_ma,ndx_nb,ndx_cor,ndx_vap !,nma_dst,nnb_dst
integer nh2o_blk,nice_blk

integer, dimension(1) :: locma,locnb
integer  :: id, jd, i, j, k, l,ilay
real*8 Rn,Rs,cst,cst2,fact
real*8 qsat2,qevap,qcond,tnewe,tnewc,qvap
logical :: used
! for microphysics time sampling :
integer imicro ! number of microphysics timesteps
real*8 microdt ! time step used for microphysics
!     Number of layer midpoints
integer :: nz, lbot

!     **********************
!     Atmospheric variable on regular vertical grid
!     **********************
real*8, dimension(size(t,1),size(t,2)) :: h2ocol  ! h2o column abundance [kg/m2]

real*8, dimension(size(t,1),size(t,2),size(t,3)) :: rho ! Atmospheric density 
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: rh2o ! Particle radius (m) 
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: sat_ratio ! Saturation
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: tini ! Updated Atmospheric temperature

!     Tracers MMR (kg/kg) 
real*8, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rini   ! Updated tracer field

real*8, dimension(size(r,1),size(r,2),size(r,3)) :: rfin 

real*8, dimension(size(t,1),size(t,2),size(t,3)) :: delp



!     **********************
!     Atmospheric variable on new vertical grid
!     ********************** 
!     Temperature (K) at boundaries and midpoints
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3) :: tl,tl0,tl_ini,tl_fin
!     Pressure (Pa) at boundaries and midpoints
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3) :: pl
! Tracer field on specific new vertical grid
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3) :: qpi,qpisub,qpicld0,qpivap0, &
                                                        qpicld_ini,qpivap_ini,qpicld,qpivap

! constants
real*8 deposit_h2o
real*8 No_h2o,Mo
real*8 ym,prec, scale_dt


! *******************************************************************************     
! ************************ Initializations **************************************     
! *******************************************************************************     
mcpu0 = (mpp_pe() == mpp_root_pe())
id= size(tl,1); jd= size(tl,2)

! Get the indices of bulk scheme tracers
nh2o_blk = find_field_index(MODEL_ATMOS, 'vap_mass_blk')
nice_blk = find_field_index(MODEL_ATMOS, 'ice_mass_blk')

! Constants !
No_h2o = ccn_blkh2o   ! Number of cloud condensation nuclei (# / kg of gaseous CO2)
scale_dt = 3.0*(kboltz/1.66054E-27)/44.0
!dev = 0.3087     ! gives an effective variance of 0.1

! Original Vertical Levels
nz = size(t,3)

! Fields initialized to 0
tdt_blkh2o(:,:,:) = 0.0
rdt_blkh2o(:,:,:,:) = 0.0

! Arrays
qpicld0(:,:,:) = 0.0
qpivap0(:,:,:) = 0.0
tl0(:,:,:) = 0.0
qpicld(:,:,:) = 0.0
qpivap(:,:,:) = 0.0
qpicld_ini(:,:,:) = 0.0
qpivap_ini(:,:,:) = 0.0
tl(:,:,:) = 0.0
pl(:,:,:) = 0.0
rfin = 0.d0

! ******************************************************************************* 
!!                          Define P, T, Q on new grids
! ******************************************************************************* 
!! Used for clouds microphysics :
!!   Updated temperature on original vertical grid : tini  
!!   Updated temperature on new vertical grid : tl_fin 
!!   Initial temperature on original vertical grid : t  
!!   Initial temperature on new vertical grid : tl_ini

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

! Update tracer field
rini(:,:,:,nh2o_blk)=r(:,:,:,nh2o_blk)+(rdt(:,:,:,nh2o_blk)*dtime)
rini(:,:,:,nice_blk)=r(:,:,:,nice_blk)+(rdt(:,:,:,nice_blk)*dtime)

! Get tracer fields on new vertical grid : qpi, qpi0
do n=1,nz ! levels 
    k = 2*n
    qpicld(:,:,k+2) = max(rini(:,:,n,nice_blk),0.)
    qpicld0(:,:,k+2) = max(r(:,:,n,nice_blk),0.)
    qpivap(:,:,k+2) = max(rini(:,:,n,nh2o_blk),0.)
    qpivap0(:,:,k+2) = max(r(:,:,n,nh2o_blk),0.)
end do

!Delp
do n=1,nz
    delp(:,:,n)= p_half(:,:,n+1)-p_half(:,:,n)
enddo

!! Saving current T and Q field before cloud scheme
qpicld_ini(:,:,:)=qpicld(:,:,:)
qpivap_ini(:,:,:)=qpivap(:,:,:)
tl_ini(:,:,:)=tl(:,:,:)
tl_fin(:,:,:)=tl(:,:,:)


! Begin H2O Cloud Condensation/Evaporation

! Loop through I and J
do i=is,ie
    do j=js,je
        do n = 1, nz
          k = n * 2 + 2

          call watsat(tl_fin(i,j,k),pl(i,j,k),qsat2)

          if (qpicld(i,j,k) .gt. 0. .and. qpivap(i,j,k) .lt. (qsat2*0.999999)) then   ! cloud evaporation
            qevap=-2.0*qsat2           ! initial guess
            if (qpicld(i,j,k) .gt. qsat2) then
              qevap=-2.0*qpicld(i,j,k)
            end if
            tnewe=tl_fin(i,j,k)
            call findq(qpivap(i,j,k),pl(i,j,k),tnewe,qevap,lw,cp_air)

            if ((-1.*qevap) .le. qpicld(i,j,k)) then
              tl_fin(i,j,k)=tl_fin(i,j,k)+qevap*lw/cp_air
              qpicld(i,j,k)=qpicld(i,j,k)+qevap
              qpivap(i,j,k)=qpivap(i,j,k)-qevap
            else
              tl_fin(i,j,k)=tl_fin(i,j,k)-qpicld(i,j,k)*lw/cp_air 
              qpivap(i,j,k)=qpivap(i,j,k)+qpicld(i,j,k)
              qpicld(i,j,k)=0.0
            end if
          end if

          if (qpivap(i,j,k) .gt. (1.000001*qsat2)) then   ! condensation
            qcond=qpivap(i,j,k)         ! initial guess
            tnewc=tl_fin(i,j,k)
            call findq(qpivap(i,j,k),pl(i,j,k),tnewc,qcond,lw,cp_air)
            tl_fin(i,j,k)=tl_fin(i,j,k)+qcond*lw/cp_air
            qpivap(i,j,k)=qpivap(i,j,k)-qcond 
            qpicld(i,j,k)=qpicld(i,j,k)+qcond
          end if

          if (qpicld(i,j,k) .gt. prec_threshold) then
            ym = (pl(i,j,k+1)-pl(i,j,k-1))/grav  !layer mass
            prec=ym*(qpicld(i,j,k)-prec_threshold)
            sfc_frost_blk(i,j,1)=sfc_frost_blk(i,j,1)+prec
            ! Total cumulative precip bucket
            cumulative_prec_blk(i,j,1)=cumulative_prec_blk(i,j,1)+prec
            ! Check if rain or snow
            if (tl(i,j,2*nz+2) .ge. 273.) then
              cumulative_prec_blk(i,j,2)=cumulative_prec_blk(i,j,2)+prec
            else
              cumulative_prec_blk(i,j,3)=cumulative_prec_blk(i,j,3)+prec
            end if

            qpicld(i,j,k)=prec_threshold

          end if

        end do  ! condensation loop through levels


! Begin H2O Cloud sedimentation calculations

        deposit_h2o=0.

        do n = 1, nz
          k   = n * 2 + 2
!         Compute atmospheric density with updated temperature
          rho(i,j,n) = pl(i,j,k) / ( rdgas*tl_fin(i,j,k) )

          Mo = qpicld(i,j,k)

          rh2o(i,j,n) = ( Mo / No_h2o * 0.75 / pi / aerdensh2o )**(1./3.)
          rh2o(i,j,n)=max(rh2o(i,j,n),1.e-7)
        enddo

        call sedimh2o(dtime,pl(i,j,:),tl(i,j,:),rho(i,j,:),rkh(i,j,:),  & 
                      qpicld(i,j,:),rh2o(i,j,:),aerdensh2o,deposit_h2o,scale_dt,nz)

        sfc_frost_blk(i,j,1)=sfc_frost_blk(i,j,1) + deposit_h2o

        cumulative_prec_blk(i,j,1)=cumulative_prec_blk(i,j,1)+deposit_h2o
        ! check if rainfall or snowfall
        if (tl(i,j,2*nz+2) .ge. 273.) then
          cumulative_prec_blk(i,j,2)=cumulative_prec_blk(i,j,2)+deposit_h2o
        else
          cumulative_prec_blk(i,j,3)=cumulative_prec_blk(i,j,3)+deposit_h2o
        end if
 
    enddo    ! j
enddo     ! i

!! Get final tendencies 
do n=1,nz
    k = 2*n
    rdt_blkh2o(:,:,n,nice_blk)=(qpicld(:,:,k+2)-qpicld_ini(:,:,k+2))/dtime
    rdt_blkh2o(:,:,n,nh2o_blk)=(qpivap(:,:,k+2)-qpivap_ini(:,:,k+2))/dtime
    tdt_blkh2o(:,:,n)=(tl_fin(:,:,k+2)-tl_ini(:,:,k+2))/dtime
end do


fact= 1.0/grav

h2ocol(is:ie,js:je)= 0.0
do n= 1, nz
    h2ocol(is:ie,js:je)= h2ocol(is:ie,js:je) + (r(:,:,n,nice_blk)+dtime*(rdt(:,:,n,nice_blk) &
                         +rdt_blkh2o(:,:,n,nice_blk)))*delp(:,:,n)*fact
enddo

rfin = r(:,:,:,nice_blk)+dtime*(rdt(:,:,:,nice_blk)+rdt_blkh2o(:,:,:,nice_blk))

if (id_blkh2ocld > 0) used = send_data ( id_blkh2ocld, rfin, time, is, js )
if (id_blkh2ocld_col > 0) used = send_data ( id_blkh2ocld_col, h2ocol(is:ie,js:je), time, is, js )
if (id_blkh2ocld_r > 0)  used = send_data ( id_blkh2ocld_r, rh2o, time, is, js )
if (id_cprecip > 0)  used =send_data ( id_cprecip, cumulative_prec_blk(:,:,1), time, is, js)
if (id_cprecip_rain > 0)  used =send_data ( id_cprecip_rain, cumulative_prec_blk(:,:,2), time, is, js)
if (id_cprecip_snow > 0)  used =send_data ( id_cprecip_snow, cumulative_prec_blk(:,:,3), time, is, js)


return
end subroutine blkh2o_driver

!!!

!****************************************************************
      subroutine findq(qvap,p,temp,qcond,lw,cp)
!*     Start moddfied Newton-Raphson method to find             *
!*     the amount of condensed/sublimed ice required to         *
!*     reach saturation. Iterations are made to                 *
!*     solve the F(Qc)=S-1=0 equation (Qc=condensed mass,       *
!*     Saturation ratio) in the case of latent heat             *
!*     release associated to Qc..                               *
!****************************************************************
      implicit none

      real*8 qvap,p,temp,qcond

      real*8 x1,x2,xl,xh,f,df,fl,fh
      real*8 dx,dxold,tempo
      real*8 rtsafe,lw,cp

      integer i

      x1 = 0.
      x2 = qcond
      call newton(x1,qvap,p,temp,fl,df,lw,cp)
      call newton(x2,qvap,p,temp,fh,df,lw,cp)

      if (fl*fh.ge.0.) then
        print*,'root not bracketed'
        print*,'x2 ',x2,' qvap ',qvap,' T ',temp
        print*,'p ',p,' fl ', fl,' df ',df
        stop
      endif
      if (fl.lt.0.) then
        xl = x1
        xh = x2
      else
        xh = x1
        xl = x2
      endif
      rtsafe = 0.5 * (x1+x2)
      dxold  = abs(x2-x1)
      dx     = dxold
      call newton(rtsafe,qvap,p,temp,f,df,lw,cp)
      do i = 1, 500
        if ( ((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0. &
            .or. abs(2.*f).gt.abs(dxold*df) ) then
          dxold = dx
          dx    = 0.5 * (xh-xl)
          rtsafe= xl+dx
          if (xl.eq.rtsafe) goto 667
        else
          dxold = dx
          dx    = f / df
          tempo  = rtsafe
          rtsafe= rtsafe - dx
          if (tempo.eq.rtsafe) goto 667
        endif
        if (abs(dx).lt.1.e-8) goto 667
          call newton(rtsafe,qvap,p,temp,f,df,lw,cp)
          if (f.lt.0) then
            xl = rtsafe
          else
            xh = rtsafe
          endif
        enddo
        print*,'500 reached'
667     continue

        qcond = rtsafe

        return
        end subroutine findq

!*********************************************************
      subroutine newton(dq,q,press,temp,f,df,lw,cpp)
!*     subroutine called during the Newton-Raphson loop  *
!*     to get the function and its derivative for the    *
!*              condensed mass determination             *
!*********************************************************

      implicit none

      real*8 dq,q,press,temp
      real*8 newT
      real*8 qsat,lw,cpp

      real*8 f,df

      newT = temp + dq * lw / cpp
      newT = max(newT,60.)

      call watsat(newT,press,qsat)
      qsat = max(qsat,1.d-50)

      f    = (q-dq) / qsat - 1.
      df   = - 2. / qsat - (f+1.) * (6146.1*lw/cpp) / newT**2.

      return
      end subroutine newton


!****************************************************************
!*                                                              *
      subroutine sedimh2o(dt,pls,tls,rhos,kd,qpic,rh2os,aerdensh2o, &
                          deposit_h2o,scale_dt,nz)
!*                                                              *
!*                Computing aerosol sedimentation               *
!*                                                              *
!****************************************************************

      implicit none

! Arguments
! ---------

      integer nz 

      real*8 dt              ! comp3 time step
      real*8 pls(2*nz+3),tls(2*nz+3)
      real*8 rhos(nz)         ! atmospheric density
      real*8 kd(2*nz+1)        ! Eddy mixing coefficient (m2/s)
      real*8 qpic(2*nz+3)! Integrated concentrations of particle/each bin (#/m3)
      real*8 rh2os(nz)

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
      real*8 cold(2),deposit_h2o
      real*8 c(nz),aerdensh2o,q1,q2
      real*8 scale_dt

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

      c(l) = qpic(ilay) * rhos(l)

      if (l.eq.1) goto 20

!     Compute fall velocity
      ilev = 2 * l + 3
      call fallvel(scale_dt,aerdensh2o, &
                   rh2os(l),tls(ilev),rhob(l),vf)

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
          ft(l)=( w1 + log(rhos(l-1)/rhos(l))*kd(2*l-1)/dzbX ) &
                / ( dexp(2.*theta) - 1. )
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
        if ( cs(l) .gt. 0. ) goto 1000
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
        c(l) = ( bs(l)*cold(1) + cs(l)*c(l)  &
                   + ds(l)*c(l+1) ) / as(l)
        cold(1)  = cold(2)
      enddo

! Compute the mass of water ice falling on the ground
      deposit_h2o = deposit_h2o &
                         + c(nz) * ft(nz+1) * dt

      c(nz) = ( bs(nz)*cold(1) + cs(nz)*c(nz) ) / as(nz)

      do l = 1, nz
        qpic(2*l+2) = c(l) / rhos(l)
      enddo

      GOTO 1

1000  continue

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
        qpic(2*l+2) = c(l) / rhos(l)
      enddo

! Compute the mass of water ice falling on the ground
      deposit_h2o = deposit_h2o  &
                        + c(nz) * ft(nz+1) * dt

1     CONTINUE


      RETURN

end subroutine sedimh2o

!!!!



!**************************************************************************
!**************************************************************************

subroutine fallvel(scale,dpden,r,t,rho,vf)
!  Calculate the fall velocity of dust particles in the Martian
!  atmosphere at the model levels. 
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


subroutine watsat(temp,press,qsat)

implicit none

real*8 pvs,temp,press,qsat

! Vapor pressure for Pa
pvs  = 611.0*exp(22.5*(1.0-(273.16/temp)))

qsat = pvs * 18.0 / (44.0*press)

return
end subroutine watsat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine blkh2o_driver_init( nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )
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
read (input_nml_file, nml=blkh2oclouds_nml, iostat=io)
ierr = check_nml_error(io,'blkh2oclouds_nml')

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=blkh2oclouds_nml)

! *********************************************************
!   --- Allocate T + Q fields needed for cloud scheme ---
! *********************************************************
allocate (  tlprev(is:ie,js:je,2*nlevels+3)  )
allocate (  qpiprev(is:ie,js:je,2*nlevels+3)  )
allocate (  qpi_dmprev(is:ie,js:je,2*nlevels+3)  )


! *********************************************************
!     ----- register diagnostic fields -----
! *********************************************************


id_blkh2ocld = register_diag_field ( mod_name, 'blkh2ocld',  &
                                 (/axes(1:3)/), Time,           &
                                'bulk h2o cld ', '',        &
                                 missing_value=missing_value )

id_blkh2ocld_col = register_diag_field ( mod_name, 'blkh2ocld_col',  &
                                 (/axes(1:2)/), Time,           &
                                'bulk h2o cld column ', '',        &
                                 missing_value=missing_value )

id_cprecip = register_diag_field ( mod_name, 'cprecip',  &
                                 (/axes(1:2)/), Time,           &
                                'cumulative total precipitation bulk scheme', '',     &
                                 missing_value=missing_value )

id_cprecip_rain = register_diag_field ( mod_name, 'cprecip_rain',  &
                                 axes(1:2), Time,           &
                                'cumulative total rainfall bulk scheme', '',     &
                                 missing_value=missing_value )

id_cprecip_snow = register_diag_field ( mod_name, 'cprecip_snow',  &
                                 (/axes(1:2)/), Time,           &
                                'cumulative total snowfall bulk scheme', '',     &
                                 missing_value=missing_value )

id_blkh2ocld_r = register_diag_field ( mod_name, 'blkh2ocld_rad',  &
                                 (/axes(1:3)/), Time,           &
                                'bulk h2o ice particle radius ', '',        &
                                 missing_value=missing_value )

id_blkh2o_sed_dt = register_diag_field ( mod_name, 'blkh2o_sed_dt',  &
                                 (/axes(1:3)/), Time,           &
                                'bulk h2o cloud sedimentation tendency', '',        &
                                 missing_value=missing_value )

id_blkh2ocld_gen = register_diag_field ( mod_name, 'dblkh2ocld_dt',  &
                                 (/axes(1:3)/), Time,           &
                                'bulk h2o cloud microphysics tendency', '',        &
                                 missing_value=missing_value )

id_blkh2o_sed_v = register_diag_field ( mod_name, 'blkh2osedvel',  &
                                 (/axes(1:3)/), Time,           &
                                'bulk h2o cloud sedimentation velocity', '',        &
                                 missing_value=missing_value )



end subroutine blkh2o_driver_init



end module blkh2omod_mgcm
