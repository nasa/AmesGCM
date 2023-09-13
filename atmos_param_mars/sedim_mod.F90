module sedim_mod
!  module to calculate moment tracer sedimentation	
use constants_mod, only: grav,cp_air,rdgas,pi
use initracer_mod
use fms_mod, only: error_mesg, FATAL,                      &
       open_namelist_file, check_nml_error,                &
       mpp_pe, mpp_root_pe, close_file,                    &
       write_version_number, stdlog,                       &
       uppercase, read_data, write_data, field_size
use fms2_io_mod,            only:  file_exists
use time_manager_mod, only: time_type, get_time
use diag_manager_mod, only: register_diag_field, send_data,  diag_axis_init
use dust_update_mod, only: sfc_dust_mass
use testconserv_mod
use rtmod_mgcm   ! for call to dtridgl TB18q
use mars_surface_mod,  only: sfc_frost_mom
use field_manager_mod, only: find_tagging_index

implicit none

public :: sedim_driver,sedim_driver_init

logical ::  scaveng=.true.    !  Sedimentation for core dust

namelist /sedim_nml/ scaveng

!-------------------- Other -----------------------------------------
logical ::  mcpu0
integer, dimension(:),  allocatable  ::   id_dens,id_rcor,id_tends_sedi  

real*8, parameter   :: missing_value = -1.e10
character(len=11) :: mod_name = 'sedim'

contains

!=====================================================================
!=====================================================================

subroutine sedim_driver(is,js,Time,p_half,p_full,t,tdt, &
     r,rdt,tg,dtime,kd,rdt_sedim,lifting_dust,checkcons)
     !
! Calls init_sedim and sedim to perform aerosol sedimentation and output 
! tendencies for each aerosol tracers
 
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
!===============================================================
!===============================================================
!     Local Variables To Microphys
integer n, nt, ndx 
integer  :: ie, je, id, jd, i, j, k, l,ilay 
integer, dimension(ntrace_mom) :: sedimflag ! active or desactive sedimentation
logical :: used
!     Number of layer midpoints
integer :: nz
!     Surface fields
real*8, dimension(size(t,1),size(t,2),ntrace_mom) :: deposit  ! tendency Tracer field on surface
real*8, dimension(size(r,1),size(r,2),ntrace_mom) :: qpig  ! tracer on surface
real*8, dimension(size(t,1),size(t,2), nice_mass) :: frost_mom0
real*8, dimension(size(t,1),size(t,2), ndust_mass) :: sfc_dst_mom0
!     **********************
!     Atmospheric variable on regular vertical grid
!     **********************
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: rho ! Atmospheric density 
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: tini ! Updated Atmospheric temperature
!     Tracers MMR (kg/kg) 
real*8, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rini   ! Updated tracer field
!     **********************
!     Atmospheric variable on new vertical grid
!     ********************** 
!     Temperature (K), pressure (Pa) at boundaries and midpoints
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3) :: tl, pl
! Tracer field on specific new vertical grid
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3,ntrace_mom) :: qpi,qpi_ini
! Density and Radius of dust+ice particles
real*8, dimension(size(t,1),size(t,2),size(t,3),ntrace_mom) :: dens,Rcor
     
! *******************************************************************************     
! ************************ Initializations **************************************     
! *******************************************************************************     
mcpu0 = (mpp_pe() == mpp_root_pe())
id= size(tl,1); jd= size(tl,2) 
ie= is + id - 1
je= js + jd - 1
! Original Vertical Levels
nz = size(t,3)

qpi(:,:,:,:) = 0.0 
qpi_ini(:,:,:,:) = 0.0 
tl(:,:,:) = 0.0 !
pl(:,:,:) = 0.0 !
qpig(:,:,:) = 0.   ! reservoirs on ground set to 0...

! *******************************************************************************     
! ************************ Atm P, T, Q  *****************************************
! *******************************************************************************     
!!   Updated temperature / pressure / traer field on original vertical grid : tini, pfull, rini
!!   Updated temperature / pressure / tracer field on new vertical grid : tl, pl, qpi
!!   Temperature / tracer tendency of all previous processes on the original vertical grid : tdt, rdt

!! Update field t => tini = t + tdt * dt
tini(:,:,:)=t(:,:,:)+tdt(:,:,:)*dtime
!! Update Tracers fields rini = r + rdt * dt
do nt=1,ntrace_mom
    ndx=mom_indx(nt)
    rini(:,:,:,ndx)=r(:,:,:,ndx)+rdt(:,:,:,ndx)*dtime
end do
! Used to check mass conservation:
frost_mom0(:,:,:)=sfc_frost_mom(:,:,:)
sfc_dst_mom0(:,:,:)=sfc_dust_mass(:,:,:)

!! Get tl and pl

! Pressure at the model top
pl(:,:,3) = p_half(:,:,1)
pl(:,:,2) = pl(:,:,3)*0.5
pl(:,:,1) = pl(:,:,2)*1.e-6

do k= 1, nz
    n= 2*k + 2     
    tl(:,:,n)= tini(:,:,k)
    pl(:,:,n)= p_full(:,:,k)
    pl(:,:,n+1)=p_half(:,:,k+1)
    if (k.eq.nz) then
        tl(:,:,n+1) = 0.5*( tg(:,:)+tini(:,:,k) )
    else
        tl(:,:,n+1)= 0.5*( tini(:,:,k+1)+tini(:,:,k) )
    end if
enddo
tl(:,:,3) = tl(:,:,4)

!! atm density
do l = 1, nz
    ilay   = l * 2 + 2
    rho(:,:,l) = pl(:,:,ilay) / ( rdgas*tl(:,:,ilay) )
enddo

!! get tracer fields on new vertical grid : qpi
do n=1,nz ! levels 
    k = 2*n
    do nt=1,ntrace_mom
        ndx=mom_indx(nt)
        qpi(:,:,k+2,nt) = max(rini(:,:,n,ndx),0.)  
    end do
end do

! *****************************************************************     
! ************************ Sedimentation **************************     
! *****************************************************************     
dens(:,:,:,:)=0. ! aerosol density
Rcor(:,:,:,:)=0. ! sedimentation radius

!! Saving current tracer field before sedimentation scheme
qpi_ini(:,:,:,:)=qpi(:,:,:,:)  

do i=is,ie
    do j=js,je
        !! Prepare sedimentation : Get dens, Rcor and sedimflag
        call init_sedim(nz,dtime,qpi(i,j,:,:),dens(i,j,:,:),Rcor(i,j,:,:),sedimflag)
        !! sedimentation over the new (2*k) vertical grid
        call sedim(nz,dtime,pl(i,j,:),tl(i,j,:),rho(i,j,:),kd(i,j,:),qpi(i,j,:,:),Rcor(i,j,:,:),dens(i,j,:,:),deposit(i,j,:),sedimflag)
    enddo
enddo

!! Get final tendencies over the original vertical grid after sedimentation
do n=1,nz  
    k = 2*n
    do nt=1,ntrace_mom
        ndx=mom_indx(nt)
        rdt_sedim(:,:,n,ndx)=(qpi(:,:,k+2,nt)-qpi_ini(:,:,k+2,nt))/dtime
        if (n.eq.nz) then
            where(lifting_dust(:,:,ndx))
                rdt_sedim(:,:,n,ndx)=0.d0
                deposit(:,:,nt)=0.d0
            end where
        end if
    end do
end do     

! *****************************************************************     
! ****************** Update surface reservoirs ********************
! *****************************************************************     
do nt=1,ntrace_mom
    ndx=mom_indx(nt)
    do n=1,ndust_mass
        if (ndx.eq.dust_mass_indx(n)) then
            sfc_dust_mass(:,:,n)=sfc_dust_mass(:,:,n)+deposit(:,:,nt)
        endif
    enddo
    do n=1,nice_mass
        if (ndx.eq.vapor_indx(n)) then
            sfc_frost_mom(:,:,n)=sfc_frost_mom(:,:,n)+deposit(:,:,nt)
        endif
    enddo
enddo

if (checkcons) then     
    call testneg(is,ie,js,je,2*nz+2,qpi,3,-1.d-15,'qpsedim_ima')
    call testneg(is,ie,js,je,2*nz+2,qpi,6,-1.d-15,'qpsedim_vap')
    call testneg(is,ie,js,je,2*nz+2,qpi_ini,3,-1.d-15,'qpinsed_ima')
    call testneg(is,ie,js,je,2*nz+2,qpi_ini,6,-1.d-15,'qpinsed_vap')
    call testneg(is,ie,js,je,nz,rini,8,-1.d-15,'microin_ima')
    call testneg(is,ie,js,je,nz,rini,11,-1.d-15,'microin_vap')
    call testneg(is,ie,js,je,nz,rini+rdt_sedim*dtime,8,-1.d-15,'sedimen_ima')
    call testneg(is,ie,js,je,nz,rini+rdt_sedim*dtime,11,-1.d-15,'sedimen_vap')
    call checkconserv(is,ie,js,je,nz,p_half,rini,frost_mom0,sfc_dst_mom0,p_half,rini+rdt_sedim*dtime,sfc_frost_mom(:,:,:),sfc_dust_mass(:,:,:),'sedimen')
endif

! *****************************************************************     
! ************************ OUTPUT *********************************     
! *****************************************************************     

do nt=1,ntrace_mom
    ndx=mom_indx(nt)
    if( id_tends_sedi(nt) > 0 )  used = send_data( id_tends_sedi(nt), deposit(is:ie,js:je,nt)/dtime, Time, is, js )
    if( id_dens(nt) > 0 )  used = send_data( id_dens(nt), dens(:,:,:,nt), Time, is, js )
    if( id_rcor(nt) > 0 )  used = send_data( id_rcor(nt), Rcor(:,:,:,nt), Time, is, js )
enddo

return
end subroutine sedim_driver

!****************************************************************
!****************************************************************

subroutine init_sedim(nz,dt,qpi,dens,Rcor,sedimflag)
! initialize module
use constants_mod, only: pi
use initracer_mod

implicit none

!  Arguments
!  ---------

real*8, intent(in) :: dt              ! physical time step (s)
integer, intent(in) :: nz              

!    Tracers :
real*8, intent(in),dimension(:,:) :: qpi  ! tracer (kg/kg)
! Output 
real*8, intent(out),dimension(nz,size(qpi,2)) :: dens
real*8, intent(out),dimension(nz,size(qpi,2)) :: Rcor
integer, intent(out),dimension(ntrace_mom) :: sedimflag ! active or desactive sedimentation

!  Local variables
!  ---------------
integer i,j,l,n,nt   ! Loop integers
integer ilay,iq
integer nt_sed  !index for tag
logical mcpu0
real*8, dimension(size(qpi,1),size(qpi,2)) :: qpitrac
real*8 Mo,No

!  Treatment 
!  ---------
mcpu0 = (mpp_pe() == mpp_root_pe())
sedimflag(:)=0

do l = 1, nz
    ilay   = l * 2 + 2

    do n = 1, ntrace_mom 
        qpitrac(ilay,n) = max(qpi(ilay,n),1.e-31)
    enddo 

    do n = 1, ntrace_mom   
        if (n==iMa_dt) then
                dens(l,n) = aerdens(iMa_dt)
                Rcor(l,n) = ( qpitrac(ilay,iMa_dt) / (qpitrac(ilay,iNb_dt)) * 0.75 / pi / dens(l,n) )**(athird) &
                * exp( 3.* stdv(n)**2. )
                Rcor(l,n) = min( max(Rcor(l,n),aerad(1)) , aerad(nbin) )
                sedimflag(n)=1
        else if (n==iNb_dt) then
                dens(l,n) = aerdens(iMa_dt)
                Rcor(l,n) = Rcor(l,iMa_dt) * exp( -3.* stdv(n)**2. )
                sedimflag(n)=1
        else if (n==iMa_cld) then
                Mo = qpitrac(ilay,iMa_cld) + qpitrac(ilay,iMa_cor) 
                No = qpi(ilay,iNb_cld) 
                dens(l,n) = qpitrac(ilay,iMa_cld) / Mo * aerdens(iMa_cld) &
                +qpitrac(ilay,iMa_cor) / Mo * aerdens(iMa_dt)
                dens(l,n) = min(max(dens(l,n),dpden_ice),dpden_dt)   ! can not be heavier than ice, or lighter than dust
                Rcor(l,n) = ( (qpitrac(ilay,iMa_cld) + qpitrac(ilay,iMa_cor) ) / (qpitrac(ilay,iNb_cld))&
                * 0.75 / pi / dens(l,n) )**(athird) &
                * exp( 3.* stdv(n)**2. )
                Rcor(l,n) = min( max(Rcor(l,n),aerad(1)) , 100.*aerad(nbin) )     ! aerad is set to 4 : descibe population
                sedimflag(n)=1
        else if (n==iNb_cld) then
                dens(l,n) = dens(l,iMa_cld)
                Rcor(l,n) = Rcor(l,iMa_cld) * exp( -3.* stdv(n)**2. )
                sedimflag(n)=1
        else if (n==iMa_cor) then
                dens(l,n) = dens(l,iMa_cld)
                Rcor(l,n) = Rcor(l,iMa_cld)
                if (scaveng) then
                    sedimflag(n)=1
                else
                    sedimflag(n)=0
                endif
        end if

        !! For tags or additional Reff
        if (ndust_mass.gt.1) then
            do nt=2,ndust_mass
                if (mom_indx(n).eq.dust_mass_indx(nt)) then
                    if (dust_mass_indx(nt) .gt. nt_nontag) then
                        nt_sed = sedim_indx(n)
                        dens(l,n) = aerdens(iMa_dt)
                        Rcor(l,n) = Rcor(l,nt_sed)
                        sedimflag(n)=1
                    else
                        dens(l,n) = aerdens(iMa_dt)
                        Rcor(l,n) = ( qpitrac(ilay,n) / (qpitrac(ilay,n+1)) * 0.75 / pi / dens(l,n) )**(athird) &
                        * exp( 3.* stdv(n)**2. )
                        Rcor(l,n) = min( max(Rcor(l,n),aerad(1)) , aerad(nbin) )
                        sedimflag(n)=1
                    endif
                elseif (mom_indx(n).eq.dust_nb_indx(nt)) then
                    if (dust_mass_indx(nt) .gt. nt_nontag) then
                        nt_sed = sedim_indx(n)
                        dens(l,n) = aerdens(iMa_dt)
                        Rcor(l,n) = Rcor(l,nt_sed) 
                        sedimflag(n)=1
                    else
                        dens(l,n) = aerdens(iMa_dt)
                        Rcor(l,n) = Rcor(l,n-1) * exp( -3.* stdv(n)**2. )
                        sedimflag(n)=1
                    endif
                elseif (mom_indx(n).eq.dust_cor_indx(nt)) then
                    dens(l,n) = dens(l,iMa_cld)
                    Rcor(l,n) = Rcor(l,iMa_cld) 
                    sedimflag(n)=sedimflag(iMa_cor)
                endif
            enddo ! ndust_mass
        endif ! (ndust_mass.gt.1)

        if (nice_mass.gt.1) then
            do nt=2,nice_mass
                if (mom_indx(n).eq.ice_mass_indx(nt)) then
                    dens(l,n) = dens(l,ima_cld)
                    rcor(l,n) = rcor(l,ima_cld)
                    sedimflag(n)=1
                elseif (mom_indx(n).eq.ice_nb_indx(nt)) then
                    dens(l,n) = dens(l,inb_cld)
                    rcor(l,n) = rcor(l,inb_cld) 
                    sedimflag(n)=1
                endif
            enddo ! nice_mass
        endif ! (nice_mass.gt.1)

    enddo   ! ntrace micro
enddo  ! ilay

return
end

!**************************************************************************
!**************************************************************************

!****************************************************************
subroutine sedim(nz,dt,pl,tl,rho,kd,&
               qpi,Rcor,dens,deposit,sedimflag)
!                Computing aerosol sedimentation               
!

use constants_mod, only: grav,cp=>cp_air,rgas=>rdgas,pi

implicit none

!c Arguments
!c ---------
! IN 
integer, intent(in) :: nz           ! dimension of the c array
real*8, intent(in) :: dt              ! time step
real*8, intent(in), dimension(:) :: pl,tl
real*8, intent(in), dimension(:) :: rho
real*8, intent(in), dimension(:) :: kd  ! Eddy mixing coefficient (m2/s)
real*8, intent(in), dimension(:,:) :: Rcor
real*8, intent(in), dimension(:,:) :: dens  ! density of aerosol (kg/m3)
integer, intent(in) :: sedimflag(:) 

! INOUT
real*8, intent(inout), dimension(:,:) :: qpi 
! OUT
real*8, intent(out), dimension(:) :: deposit   ! amount of tracer felt on the ground during dt

!c Local variables
!c ---------------

real*8 rhob(nz)       ! atmospheric density at layer boundaries
real*8 dz(nz)         ! layer thickness (m)
real*8 vf             ! fall velocity of particle (m/s)
real*8 cour,w1        ! Courant number & corrected velocity (vf+kd/H with h scale height)
real*8 dzbX
real*8 wall(nz),vfall(nz),sgmall(nz)
real*8 sigma,theta,hc,lg,rap,cmp,w,wp 
real*8 fs(nz+1),ft(nz+1)
real*8 as(nz),bs(nz),cs(nz),ds(nz)
real*8 asi(nz),bsi(nz),csi(nz),dsi(nz),xsol(nz)
real*8 cold(2)
real*8 c(nz)

integer i,l,n,nt
integer ilay,ilev
logical mcpu0

!c Treatment
!c ---------------

mcpu0 = (mpp_pe() == mpp_root_pe())
theta = 0.0    ! added 6/26/08  Tim's bug list

!c    Water ice deposit reset to zero
do i=1,ntrace_mom
    deposit(i) = 0.
enddo

!c    Layer thickness: dz
do l = 1, nz
    dz(l) = ( pl(2*l+3)-pl(2*l+1) ) / grav / rho(l)
enddo

!c    Compute density at the layer boundaries
do l = 1, nz
    ilev = l * 2 + 3
    if (l.lt.nz) then
        rhob(l) = pl(ilev)  / (rgas*tl(ilev))
    else
        rhob(l) = pl(ilev-1) / (rgas*tl(ilev-1))
    endif
enddo

!c Loop over ntrace_mom
do n = 1, ntrace_mom

    ! Loop over layers
    if (sedimflag(n).eq.1) then
        do l = 1, nz

            ilay = 2 * l + 2
            c(l) = qpi(ilay,n) * rho(l)

            if (l.eq.1) cycle

            ! Compute fall velocity
            ilev = 2 * l + 3
            call fallvel(scale_dt,dens(l,n), &
                  Rcor(l,n),tl(ilev),rhob(l),vf)

            dzbX = ( dz(l)+dz(l-1) ) / 2.
            w  = -1. * vf * exp(-stdv(n)**2.)
            vfall(l)=vf

            ! Get the corrected fall velocity (virtual speed accounting for mixing)
            if (kd(2*l-1) .ne. 0.) then
                theta = 0.5 * ( w*dzbX/kd(2*l-1) + log(rho(l-1)/rho(l)) )
                if (theta.ne.0) then
                    sigma = 1./dtanh(theta) - 1./theta
                else
                    sigma = 1.
                endif
            else
                sigma = 1.
            endif

            sgmall(l)=sigma

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

            wall(l)=w1

            !Fluxes at layer boundaries
            if (kd(2*l-1).ne.0.) then
                if (theta.ne.0.) then
                    ft(l)=( w1 + log(rho(l-1)/rho(l))*kd(2*l-1)/dzbX ) &
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

        enddo

        !c Boundary conditions for the fluxes
        fs(1)    =  0.
        ft(1)    =  0.
        fs(nz+1) =  0.
        ft(nz+1) = -w1


        !c Compute the coefficient of the continuity equation
        !c Depending on the cs value, switch to an explicit or an implicit scheme
        do l=1,nz
            cs(l) =  ft(l+1) + fs(l) - dz(l) / dt
            if ( cs(l) .gt. 0. ) goto 1000 
            as(l) = -dz(l) / dt
            bs(l) = -ft(l)
            ds(l) = -fs(l+1)
        enddo

        !c Explicit case 
        cold(1)  = c(1)
        c(1) = ( cs(1)*c(1) + ds(1)*c(2) ) / as(1)

        do l = 2, nz-1
            cold(2)  = c(l)
            c(l) = ( bs(l)*cold(1) + cs(l)*c(l) &
                  + ds(l)*c(l+1) ) / as(l)
            cold(1)  = cold(2)
        enddo 

        !c Compute the mass of water ice falling on the ground
        !! The ice sedimented on surface is stored as surface vap deposit (iMa_vap)
        if (n.eq.iMa_cld) deposit(iMa_vap) = deposit(iMa_vap) &
                                 + c(nz) * ft(nz+1) * dt
        if (n.eq.iMa_dt)  deposit(iMa_dt)   = deposit(iMa_dt) &
                                 + c(nz) * ft(nz+1) * dt
        !! Dust core put back with dust mass
        if (n.eq.iMa_cor) deposit(iMa_dt)  = deposit(iMa_dt) &
                                 + c(nz) * ft(nz+1) * dt

        if (ndust_mass.gt.1) then
            do nt=2,ndust_mass
                if (mom_indx(n).eq.dust_mass_indx(nt)) then
                    deposit(n)  = deposit(n) + c(nz) * ft(nz+1) * dt
                endif              
                if (mom_indx(n).eq.dust_cor_indx(nt)) then  ! TB18tmp : make it better here
                    deposit(n-2)  = deposit(n-2) + c(nz) * ft(nz+1) * dt
                endif              
            enddo
        endif

        if (nice_mass.gt.1) then
            do nt=2,nice_mass
                if (mom_indx(n).eq.ice_mass_indx(nt)) then
                    deposit(n+2)  = deposit(n+2) + c(nz) * ft(nz+1) * dt
                endif              
            enddo
        endif

        c(nz) = ( bs(nz)*cold(1) + cs(nz)*c(nz) ) / as(nz)

        do l = 1, nz
            qpi(2*l+2,n) = c(l) / rho(l)
        enddo

        cycle
 
1000    continue

        !c Implicit case 

        do l = 1, nz
            asi(l) =  ft(l)
            bsi(l) = -( ft(l+1) + fs(l) + dz(l)/dt )
            csi(l) =  fs(l+1)
            dsi(l) = -dz(l) / dt * c(l)
        enddo

        !c Matrix inversion

        call dtridgl(nz,asi,bsi,csi,dsi,xsol)

        do l = 1, nz
            c(l) = xsol(l)
            qpi(2*l+2,n) = c(l) / rho(l)
        enddo

        !c Compute the mass of water ice falling on the ground

        !! Sedim is perform according to sedimflag
        if (n.eq.iMa_cld) deposit(iMa_vap) = deposit(iMa_vap) &
                                         + c(nz) * ft(nz+1) * dt
        if (n.eq.iMa_dt)  deposit(iMa_dt)   = deposit(iMa_dt)  &
                                         + c(nz) * ft(nz+1) * dt
        !! Dust core put back with dust mass
        if (n.eq.iMa_cor) deposit(iMa_dt)  = deposit(iMa_dt) &
                                         + c(nz) * ft(nz+1) * dt
        if (ndust_mass.gt.1) then
            do nt=2,ndust_mass
                if (mom_indx(n).eq.dust_mass_indx(nt)) then
                    deposit(n)  = deposit(n) + c(nz) * ft(nz+1) * dt
                endif              
                if (mom_indx(n).eq.dust_cor_indx(nt)) then
                    !! Dust core put back with dust mass TB18q : dangerous way to do it
                    deposit(n-2)  = deposit(n-2) + c(nz) * ft(nz+1) * dt
                endif              
            enddo
        endif

        if (nice_mass.gt.1) then
            do nt=2,nice_mass
                if (mom_indx(n).eq.ice_mass_indx(nt)) then
                    deposit(n+2)  = deposit(n+2) + c(nz) * ft(nz+1) * dt
                endif              
            enddo
        endif

    endif  ! sedimflag
enddo

return
end 

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
end

!****************************************************************
!****************************************************************

subroutine sedim_driver_init( nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )
! initialize module

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

! *********************************************************
!     ----- read namelist /_nml/   -----
! *********************************************************

if (file_exists('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1; do while (ierr /= 0)
        read  (unit, nml=sedim_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'sedim_nml')
    enddo
10     call close_file (unit)
endif

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=sedim_nml)

! *********************************************************
!     ----- register diagnostic fields -----
! *********************************************************
allocate ( id_tends_sedi(ntrace_mom) )
allocate ( id_dens(ntrace_mom) )
allocate ( id_rcor(ntrace_mom) )

!! Preparing output
do nt= 1, ntrace_mom
    !Get name tracer
    ndx= mom_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)
    if (mpp_pe() == mpp_root_pe()) print*,'Micro Tracer ID:',nt,ndx,tracer_name

    tname= trim(tracer_name) // '_dep' 
    id_tends_sedi(nt) = register_diag_field ( mod_name, trim(tname),  &
         axes(1:2), Time, 'Sedimentation deposition rate', 'kg/m2/s', &
         missing_value=missing_value)

    tname= trim(tracer_name) // '_dens' 
    id_dens(nt) = register_diag_field ( mod_name, trim(tname),  &
         axes(1:3), Time, 'Sedimentation density', 'kg/m3', &
         missing_value=missing_value)

    tname= trim(tracer_name) // '_rcor' 
    id_rcor(nt) = register_diag_field ( mod_name, trim(tname),  &
         axes(1:3), Time, 'Sedimentation radius', 'm', &
         missing_value=missing_value)
enddo

end subroutine sedim_driver_init

end module sedim_mod

