module micromod_mgcm
!
!  module to calculate moment microphysics condensation and nucleation
!
use constants_mod, only: grav,cp_air,rdgas,pi,kbz=>kboltz
use initracer_mod
use fms_mod, only: error_mesg, FATAL,                      &
           open_namelist_file, check_nml_error,                &
           mpp_pe, mpp_root_pe, close_file,                    &
           write_version_number, stdlog,                       &
           uppercase, read_data, write_data, field_size
use fms2_io_mod,            only:  file_exists
use time_manager_mod, only: time_type, get_time
use diag_manager_mod, only: register_diag_field, send_data,  diag_axis_init
use dust_update_mod
use testconserv_mod
use rtmod_mgcm, only: dtridgl   ! for call to dtridgl TB18q
use mars_surface_mod,  only: sfc_frost_mom

implicit none

public :: micro_driver,micro_driver_init,wcldcol,wcol

real*8    ::  microtimestep=-1.         !  timestep for microphysiq (s) 
real*8    ::  facsubl=1.                !  coefficient for sublimaiton of water ice
logical ::  latent_heat=.false.         !  changes of T during cloud formation
logical ::  makeclouds=.false.          !  clouds formation
integer ::  subeffmode=0                !  0: effective T + Q for cloud scheme are computed from initial T + Q in physics
                                        !  1: effective T + Q for cloud scheme are computed from T + Q saved after previous call to cloud scheme     
real*8    ::  mteta = 0.965                            ! Contact parameter ( m=cos(theta) orig .975 )
integer ::  mteta_case = 0              !  flag to use temperature dependent contact parameter
                                        !  0: no temperature dependence
                                        !  1: Iraci temperature dependence
                                        !  2: Trainer (2009) exponential temperature dependence
                                        !  3: Maattanen (2014) Trainer tanh temperature dependence


namelist /microphys_nml/ microtimestep,latent_heat,makeclouds,subeffmode,facsubl,mteta_case,mteta

!-------------------- Other -----------------------------------------
logical ::  mcpu0
logical,save :: firstcall_micro=.true.
integer, dimension(:),  allocatable  ::   id_tend_nucl,id_taucld,id_wcol,id_wcldcol  
integer ::   id_sat_ratio
integer ::   id_mt_3d

real*8, dimension(:,:,:), allocatable :: wcldcol, wcol

!! effective T and Q from previous timestep for cloud scheme
real*8, dimension(:,:,:),     allocatable,save :: tlprev
real*8, dimension(:,:,:,:),     allocatable,save :: qpiprev,qpi_dmprev,qpi_dnprev
real*8, parameter   :: missing_value = -1.e10
character(len=11) :: mod_name = 'microphys'

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine micro_driver(is,js,Time,p_half,p_full,t,tdt, &
         r,rdt,strss,tg,dtime,drg,rdt_micro,tdt_micro,checkcons)
!
! Main driver for moment microphysics 
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
!===============================================================
!===============================================================
!     Local Variables To Microphys
!     Water/cloud column
real*8, dimension(size(t,1),size(t,2),nice_mass) :: taucld
integer n, nt, ndx, microstep, ndx_ma,ndx_nb,ndx_cor,ndx_vap !,nma_dst,nnb_dst
integer, dimension(1) :: locma,locnb
integer  :: ie, je, id, jd, i, j, k, l,ilay 
real*8 Rn,Rs,cst,cst2
logical :: used
! for microphysics time sampling :
integer imicro ! number of microphysics timesteps
real*8 microdt ! time step used for microphysics
!     Number of layer midpoints
integer :: nz
!     Surface fields
real*8, dimension(size(r,1),size(r,2)) :: mass0,mass1,mass2,mass3  ! total mass in column 
real*8, dimension(size(r,1),size(r,2)) :: mass4,mass5,mass6,mass7  ! total mass in column 
real*8, dimension(size(t,1),size(t,2), nice_mass) :: frost_micro0
real*8, dimension(size(t,1),size(t,2), ndust_mass) :: sfc_dst_micro0
!     **********************
!     Atmospheric variable on regular vertical grid
!     **********************
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: rho ! Atmospheric density 
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: sat_ratio ! Water saturation
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: mt_3d ! runtime contact parameter
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: tini ! Updated Atmospheric temperature
real*8, dimension(size(t,3),ntrace_mom) :: ratio_mass,ratio_nb,ratio_core, ratio_vap  ! Ratio for tagging
!     Tracers MMR (kg/kg) 
real*8, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rini   ! Updated tracer field
!     Tendancies nucleation, sedimentation, update of water vapor
real*8, dimension(size(r,1),size(r,2),size(r,3),ntrace_mom) :: qdt_nucl
real*8, dimension(size(r,1),size(r,2),size(r,3),size(r,4)) :: rdt_nucl
!     Tendancies temperature (clouds scheme)
real*8, dimension(size(t,1),size(t,2),size(t,3)) :: tdt_nucl
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
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3,ntrace_mom) :: qpi,qpisub,qpi0,qpieff,qpi_ini
! Tendencies for Tracer field on specific new vertical grid
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3,ntrace_mom) :: qpisubdt,qpidteff
! Tracer field for extra dust modes
real*8, dimension(size(t,1),size(t,2),2*size(t,3)+3,ndust_mass) :: qpi_dmass0,qpi_dnum0,qpi_dmeff, &
                                                                   qpi_dneff,qpi_dmass,qpi_dnum, &
                                                                   qpidteff_m,qpidteff_n,qpi_dmini, &
                                                                   qpi_dnini,qpi_dmsub,qpi_dnsub, &
                                                                   qpisubdt_m,qpisubdt_n,qdt_nucl_m,qdt_nucl_n
logical :: multi_modes = .false.
     
! *******************************************************************************     
! ************************ Initializations **************************************     
! *******************************************************************************     
mcpu0 = (mpp_pe() == mpp_root_pe())
id= size(tl,1); jd= size(tl,2) 
ie= is + id - 1
je= js + jd - 1
! Original Vertical Levels
nz = size(t,3)

! Fields initialized to 0
tdt_micro(:,:,:) = 0.0 
rdt_micro(:,:,:,:) = 0.0 
rdt_nucl(:,:,:,:) = 0.0 

qpi0(:,:,:,:) = 0.0 
tl0(:,:,:) = 0.0 
qpi(:,:,:,:) = 0.0 
qpi_ini(:,:,:,:) = 0.0 
tl(:,:,:) = 0.0 
qpieff(:,:,:,:) = 0.0 
tleff(:,:,:) = 0.0 
pl(:,:,:) = 0.0 
qdt_nucl(:,:,:,:) = 0.d0
tdt_nucl(:,:,:) = 0.0
tldteff(:,:,:) = 0.0 
qpidteff(:,:,:,:) = 0.0  
sat_ratio(:,:,:) = 0.0 
mt_3d(:,:,:) = 0.0
qpi_dmass0(:,:,:,:) = 0.d0
qpi_dnum0(:,:,:,:) = 0.d0
qpi_dmass(:,:,:,:) = 0.d0
qpi_dnum(:,:,:,:) = 0.d0
qpi_dmeff(:,:,:,:) = 0.d0
qpi_dneff(:,:,:,:) = 0.d0
qpi_dmini(:,:,:,:) = 0.d0
qpi_dnini(:,:,:,:) = 0.d0
qdt_nucl_m = 0.d0
qdt_nucl_n = 0.d0
qpidteff_m = 0.d0
qpidteff_n = 0.d0
multi_modes = (ndust_mass .gt. 1)

! *******************************************************************************     
! ************* Timestep sampling for cloud scheme nucleacond *******************     
! *******************************************************************************     
! May be adapted for each dtime
imicro=int(dtime/microtimestep)
imicro=max(imicro,1)
microdt=dtime/dble(imicro)

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

!! Pressure at the model top
pl(:,:,3) = p_half(:,:,1)
pl(:,:,2) = pl(:,:,3)*0.5
pl(:,:,1) = pl(:,:,2)*1.e-6

!! Update field t => tini
tini(:,:,:)=t(:,:,:)+tdt(:,:,:)*dtime

!! get tl0, tl and pl
do k= 1, nz
    n= 2*k + 2     
    tl(:,:,n)= tini(:,:,k)
    tl0(:,:,n)= t(:,:,k)
    pl(:,:,n)= p_full(:,:,k)
    pl(:,:,n+1)=p_half(:,:,k+1)
    if (k.eq.nz) then
        tl(:,:,n+1) = 0.5*( tg(:,:)+tini(:,:,k) )
        tl0(:,:,n+1) = 0.5*( tg(:,:)+t(:,:,k) )
    else
        tl(:,:,n+1)= 0.5*( tini(:,:,k+1)+tini(:,:,k) )
        tl0(:,:,n+1)= 0.5*( t(:,:,k+1)+t(:,:,k) )
    end if
enddo
tl(:,:,3) = tl(:,:,4)
tl0(:,:,3) = tl0(:,:,4)

!! density
do l = 1, nz
    ilay   = l * 2 + 2
    rho(:,:,l) = pl(:,:,ilay) / ( rdgas*tl(:,:,ilay) )
enddo

!! Update Tracers fields
do nt=1,ntrace_mom
    ndx=mom_indx(nt)
    rini(:,:,:,ndx)=r(:,:,:,ndx)+rdt(:,:,:,ndx)*dtime
end do
frost_micro0(:,:,:)=sfc_frost_mom(:,:,:)
sfc_dst_micro0(:,:,:)=sfc_dust_mass(:,:,:)

!! Get tracer fields on new vertical grid : qpi, qpi0
do n=1,nz ! levels 
    k = 2*n
    do nt=1,ntrace_mom
        ndx=mom_indx(nt)
        qpi(:,:,k+2,nt) = max(rini(:,:,n,ndx),0.)  
        qpi0(:,:,k+2,nt) = max(r(:,:,n,ndx),0.)  
    end do
    do nt=1,ndust_mass
        ndx=dust_mass_indx(nt)
        if (nt .le. nt_nontag) then
            qpi_dmass0(:,:,k+2,nt) = max(r(:,:,n,ndx),0.)  
            qpi_dnum0(:,:,k+2,nt) = max(r(:,:,n,ndx+1),0.)  
            qpi_dmass(:,:,k+2,nt) = max(rini(:,:,n,ndx),0.)  
            qpi_dnum(:,:,k+2,nt) = max(rini(:,:,n,ndx+1),0.)  
        end if
    end do
end do

!! Get water vapor column  : calculated after cloud scheme? 
wcol(:,:,:) = 0.
do nt=1,nice_mass
    locma=FINDLOC(mom_indx, VALUE=vapor_indx(nt))
    ndx_vap=locma(1)
    do l = 1, nz
        ilay = 2 * l + 2
        wcol(:,:,nt) = wcol(:,:,nt) + qpi(:,:,ilay,ndx_vap) * &
                 (pl(:,:,2*l+3) - pl(:,:,2*l+1)) / grav
    enddo
enddo

! *****************************************************************     
! ************************ Clouds *********************************     
! *****************************************************************     

!! The temperature and tracer tendencies from other processes are used as stepped entry for the cloud scheme 

if (makeclouds) then

    !! Effective temperature and tracer before the cloud scheme: tleff, qpieff
    ! If first time step, we use tl0 and qpi0
    if ((subeffmode.eq.0).or.(firstcall_micro)) then  ! Using field computed after the dynamics, when physics starts
        tleff(:,:,:)=tl0(:,:,:)
        qpieff(:,:,:,:)=qpi0(:,:,:,:)
        qpi_dmeff(:,:,:,:) = qpi_dmass0(:,:,:,:)
        qpi_dneff(:,:,:,:) = qpi_dnum0(:,:,:,:)
        firstcall_micro=.false.
    else                       ! Using field computed during previous time step, after cloud sheme call
        tleff(:,:,:)=tlprev(:,:,:)
        qpieff(:,:,:,:)=qpiprev(:,:,:,:)
        qpi_dmeff(:,:,:,:) = qpi_dmprev(:,:,:,:)
        qpi_dneff(:,:,:,:) = qpi_dnprev(:,:,:,:)
    endif

    !! Get effective tendency of temperature + tracer fields corresponding to previous processes 
    do nt=1,ntrace_mom
        where (qpi(:,:,:,nt) .gt. qpieff(:,:,:,nt))
            qpidteff(:,:,:,nt) = (qpi(:,:,:,nt)-qpieff(:,:,:,nt))/dtime  
        elsewhere  
            qpidteff(:,:,:,nt)  = 0.0
            qpieff(:,:,:,nt) = qpi(:,:,:,nt)
        end where
    end do
    do nt=1,ndust_mass
        where (qpi_dmass(:,:,:,nt) .gt. qpi_dmeff(:,:,:,nt))
            qpidteff_m(:,:,:,nt) = (qpi_dmass(:,:,:,nt)-qpi_dmeff(:,:,:,nt))/dtime  
        elsewhere
            qpidteff_m(:,:,:,nt) = 0.0
            qpi_dmeff(:,:,:,nt) = qpi_dmass(:,:,:,nt)
        end where
        where (qpi_dnum(:,:,:,nt) .gt. qpi_dneff(:,:,:,nt))
            qpidteff_n(:,:,:,nt) = (qpi_dnum(:,:,:,nt)-qpi_dneff(:,:,:,nt))/dtime  
        elsewhere
            qpidteff_n(:,:,:,nt) = 0.
            qpi_dneff(:,:,:,nt) = qpi_dnum(:,:,:,nt)
        end where
    end do
    do l = 1, nz
        ilay   = l * 2 + 2
        tldteff(:,:,ilay)=(tl(:,:,ilay)-tleff(:,:,ilay))/dtime
    enddo

    !! Saving current T and Q field before cloud scheme
    qpi_ini(:,:,:,:)=qpi(:,:,:,:)  
    tl_ini(:,:,:)=tl(:,:,:)  
    qpi_dmini(:,:,:,:)=qpi_dmass(:,:,:,:)  
    qpi_dnini(:,:,:,:)=qpi_dnum(:,:,:,:)  

    !! Initialization of effective T and Q field before clouds scheme
    qpisub(:,:,:,:)=qpieff(:,:,:,:)   
    tlsub(:,:,:)=tleff(:,:,:)  
    qpi_dmsub(:,:,:,:)=qpi_dmeff(:,:,:,:)   
    qpi_dnsub(:,:,:,:)=qpi_dneff(:,:,:,:)   

! ******** Loop subtimesteps ***********
    do microstep=1,imicro
        !! initializations output tendencies from cloud scheme 
        tlsubdt(:,:,:)=0.
        qpisubdt(:,:,:,:)=0.
        qpisubdt_m(:,:,:,:)=0.
        qpisubdt_n(:,:,:,:)=0.

        !! update t and q fields with part of the tendency from previous physical processes 
        tlsub(:,:,:)=tlsub(:,:,:)+tldteff(:,:,:)*microdt  
        qpisub(:,:,:,:)=qpisub(:,:,:,:)+qpidteff(:,:,:,:)*microdt  
        qpisub=max(qpisub,0.d0)
        qpi_dmsub(:,:,:,:)=qpi_dmsub(:,:,:,:)+qpidteff_m(:,:,:,:)*microdt  
        qpi_dmsub=max(qpi_dmsub,0.d0)
        qpi_dnsub(:,:,:,:)=qpi_dnsub(:,:,:,:)+qpidteff_n(:,:,:,:)*microdt  
        qpi_dnsub=max(qpi_dnsub,0.d0)

        !! compute new atmospheric density
        do l = 1, nz
            ilay   = l * 2 + 2
            rho(:,:,l) = pl(:,:,ilay) / ( rdgas*tlsub(:,:,ilay) )
        enddo

        !! call nucleacond
        do i=is,ie
            do j=js,je
                call nucleacond(nz,microdt,pl(i,j,:),tlsub(i,j,:),rho(i,j,:),qpisub(i,j,:,:), &
                                qpi_dmsub(i,j,:,:),qpi_dnsub(i,j,:,:),multi_modes, &
                                tlsubdt(i,j,:),qpisubdt(i,j,:,:),qpisubdt_m(i,j,:,:),qpisubdt_n(i,j,:,:), &
                                sat_ratio(i,j,:),mt_3d(i,j,:))
            enddo
        enddo

        !! update t and q with subtendency of the cloud scheme 
        tlsub(:,:,:)=tlsub(:,:,:)+tlsubdt(:,:,:)*microdt
        qpisub(:,:,:,:)=qpisub(:,:,:,:)+qpisubdt(:,:,:,:)*microdt
        qpi_dmsub(:,:,:,:)=qpi_dmsub(:,:,:,:)+qpisubdt_m(:,:,:,:)*microdt
        qpi_dnsub(:,:,:,:)=qpi_dnsub(:,:,:,:)+qpisubdt_n(:,:,:,:)*microdt

        qpisub=max(qpisub,0.d0)
        qpi_dmsub=max(qpi_dmsub,0.d0)
        qpi_dnsub=max(qpi_dnsub,0.d0)


    enddo  ! microstep

    !********************************************************************************* 
    !! Get T and Q final tendencies over the original vertical grid 
    !********************************************************************************* 
    do n=1,nz  
        k = 2*n
        do nt=1,6 !ntrace_mom
            ndx=mom_indx(nt)
            qdt_nucl(:,:,n,nt) = (qpisub(:,:,k+2,nt)-qpi_ini(:,:,k+2,nt))/dtime      
        end do
        if (multi_modes) then
            do nt=2,ndust_mass
                ndx=dust_mass_indx(nt)
                if (ndx .le. nt_nontag) then
                    qdt_nucl_m(:,:,n,nt) = (qpi_dmsub(:,:,k+2,nt)-qpi_dmini(:,:,k+2,nt))/dtime
                    qdt_nucl_n(:,:,n,nt) = (qpi_dnsub(:,:,k+2,nt)-qpi_dnini(:,:,k+2,nt))/dtime
                end if
            end do
        end if
        tdt_nucl(:,:,n)=(tlsub(:,:,k+2)-tl_ini(:,:,k+2))/dtime   
    end do

!cloud tags not yet working with multi modes so skip   
#ifdef SKIP
    !********************************************************************************* 
    !**** Tags following clouds scheme / nucleation **********************************
    !********************************************************************************* 
    if (ndust_mass.gt.1) then

        do nt=2,ndust_mass

            do n=1,ntrace_mom
                if (mom_indx(n).eq.dust_mass_indx(nt)) then
                    ndx_ma=n
                endif
                if (mom_indx(n).eq.dust_nb_indx(nt)) then
                    ndx_nb=n
                endif
                if (mom_indx(n).eq.dust_cor_indx(nt)) then
                    ndx_cor=n
                endif
            enddo

            do i=is,ie
                do j=js,je

                    ratio_mass(:,:)=0.d0
                    ratio_nb(:,:)=0.d0
                    ratio_core(:,:)=0.d0
                    do l = 1, nz
                        k = 2*l

                        if ( qpi_ini(i,j,k+2,ndx_cor).gt.1.e-12 .and. qpi_ini(i,j,k+2,iMa_cor)  &
                            .gt. 1.e-12 .and. qdt_nucl(i,j,l,iMa_dt) .gt. 1.e-17 ) then

                            !! Get ratio core dust
                            ratio_core(l,nt-1)=qpi_ini(i,j,k+2,ndx_cor)/qpi_ini(i,j,k+2,iMa_cor)
                            ratio_core(l,nt-1)=min(max(ratio_core(l,nt-1),0.d0),1.d0)
                            !! Remove fraction of core and add it in dust mass ( same ratio applied on Number TB18td : ice tag)
                            qdt_nucl(i,j,l,ndx_cor)=qdt_nucl(i,j,l,iMa_cor)*ratio_core(l,nt-1)
                            qdt_nucl(i,j,l,ndx_ma)=-qdt_nucl(i,j,l,ndx_cor)
                            qdt_nucl(i,j,l,ndx_nb)=qdt_nucl(i,j,l,iNb_dt)*ratio_core(l,nt-1)

                        else if ( qpi_ini(i,j,k+2,ndx_ma).gt.1.e-12 .and. qpi_ini(i,j,k+2,iMa_dt) .gt. 1.e-12 .and. & 
                            qdt_nucl(i,j,l,iMa_dt) .lt. -1.e-17 ) then ! if equal to 0 then tendency remains to 0
                            !! If dust goes into dust core:
                            ratio_mass(l,nt-1)=qpi_ini(i,j,k+2,ndx_ma)/qpi_ini(i,j,k+2,iMa_dt)
                            ratio_mass(l,nt-1)=min(max(ratio_mass(l,nt-1),0.d0),1.d0)
                            ratio_nb(l,nt-1)=qpi_ini(i,j,k+2,ndx_nb)/qpi_ini(i,j,k+2,iNb_dt)
                            ratio_nb(l,nt-1)=min(max(ratio_nb(l,nt-1),0.d0),1.d0)
                            !! Remove fraction of dust and add it in dust core
                            qdt_nucl(i,j,l,ndx_ma)=qdt_nucl(i,j,l,iMa_dt)*ratio_mass(l,nt-1)
                            qdt_nucl(i,j,l,ndx_nb)=qdt_nucl(i,j,l,iNb_dt)*ratio_nb(l,nt-1)
                            qdt_nucl(i,j,l,ndx_cor)=-qdt_nucl(i,j,l,ndx_ma)
                        endif


                    enddo

                enddo    ! j
            enddo     ! i

        enddo ! ndustmass
    endif ! if ndustmass gt 1
    if (nice_mass.gt.1) then

        do nt=2,nice_mass

            do n=1,ntrace_mom
                if (mom_indx(n).eq.ice_mass_indx(nt)) then
                    ndx_ma=n
                endif
                if (mom_indx(n).eq.ice_nb_indx(nt)) then
                    ndx_nb=n
                endif
                if (mom_indx(n).eq.vapor_indx(nt)) then
                    ndx_vap=n
                endif
            enddo

            do i=is,ie
                do j=js,je

                    ratio_mass(:,:)=0.d0
                    ratio_nb(:,:)=0.d0
                    ratio_vap(:,:)=0.d0
                    do l = 1, nz
                        k = 2*l
                        !! If vapor goes into ice :
                        if ( qdt_nucl(i,j,l,iMa_cld) .gt. 0.d0 ) then
                            !! Get ratio vapor
                            ratio_vap(l,nt-1)=qpi_ini(i,j,k+2,ndx_vap)/qpi_ini(i,j,k+2,iMa_vap)
                            ratio_vap(l,nt-1)=min(max(ratio_vap(l,nt-1),0.d0),1.d0)

                            qdt_nucl(i,j,l,ndx_ma)=qdt_nucl(i,j,l,iMa_cld)*ratio_vap(l,nt-1)
                            qdt_nucl(i,j,l,ndx_vap)=-qdt_nucl(i,j,l,ndx_ma)
                            qdt_nucl(i,j,l,ndx_nb)=qdt_nucl(i,j,l,iNb_cld)*ratio_vap(l,nt-1)

                        else if ( qdt_nucl(i,j,l,iMa_cld) .lt. 0.d0 ) then
                            !! If ice goes into vapor:
                            ratio_mass(l,nt-1)=qpi_ini(i,j,k+2,ndx_ma)/qpi_ini(i,j,k+2,iMa_cld)
                            ratio_mass(l,nt-1)=min(max(ratio_mass(l,nt-1),0.d0),1.d0)
                            ratio_nb(l,nt-1)=qpi_ini(i,j,k+2,ndx_nb)/qpi_ini(i,j,k+2,iNb_cld)
                            ratio_nb(l,nt-1)=min(max(ratio_nb(l,nt-1),0.d0),1.d0)

                            qdt_nucl(i,j,l,ndx_ma)=qdt_nucl(i,j,l,iMa_cld)*ratio_mass(l,nt-1)
                            qdt_nucl(i,j,l,ndx_vap)=-qdt_nucl(i,j,l,ndx_ma)
                            qdt_nucl(i,j,l,ndx_nb)=qdt_nucl(i,j,l,iNb_cld)*ratio_nb(l,nt-1)
                        endif

                    enddo

                enddo    ! i
            enddo     ! j

        enddo ! nicemass
    endif ! if nicemass gt 1
#endif
    !********************************************************************************* 
    !! Update physical T and Q tendencies
    !********************************************************************************* 
    do nt=1,ntrace_mom
        ndx=mom_indx(nt)
        rdt_micro(:,:,:,ndx)=rdt_micro(:,:,:,ndx)+qdt_nucl(:,:,:,nt)
    end do
    if (multi_modes) then
        do nt=2,ndust_mass
            ndx=dust_mass_indx(nt)
            if (ndx .le. nt_nontag) then
                rdt_micro(:,:,:,ndx)=rdt_micro(:,:,:,ndx)+qdt_nucl_m(:,:,:,nt)
                rdt_micro(:,:,:,ndx+1)=rdt_micro(:,:,:,ndx+1)+qdt_nucl_n(:,:,:,nt)
            end if
        end do
    endif
    tdt_micro(:,:,:)=tdt_micro(:,:,:)+tdt_nucl(:,:,:)

    !! Update current values of T + Q field on new vertical grid
    tl=tlsub   
    !! Compute new atmospheric density
    do l = 1, nz
        ilay   = l * 2 + 2
        rho(:,:,l) = pl(:,:,ilay) / ( rdgas*tl(:,:,ilay) )
    enddo

    !! Save current values of T + Q field for next timestep
    if (subeffmode.eq.1) then 
        tlprev(:,:,:)=tl(:,:,:)
        qpiprev(:,:,:,:)=qpi(:,:,:,:)
        qpi_dmeff(:,:,:,:) = qpi_dmass(:,:,:,:)
        qpi_dneff(:,:,:,:) = qpi_dnum(:,:,:,:)
    endif

    if (checkcons) then     
        call testneg(is,ie,js,je,2*nz+2,qpi,3,-1.d-15,'qclouds_ima')
        call testneg(is,ie,js,je,2*nz+2,qpi,6,-1.d-15,'qclouds_vap')
        call testneg(is,ie,js,je,nz,rini+rdt_micro*dtime,8,-1.d-15,'wclouds_ima')
        call testneg(is,ie,js,je,nz,rini+rdt_micro*dtime,11,-1.d-15,'wclouds_vap')
        call checkconserv(is,ie,js,je,nz,p_half,rini,frost_micro0,sfc_dst_micro0,p_half,rini+ &
                rdt_micro*dtime,sfc_frost_mom(:,:,:),sfc_dust_mass(:,:,:),'wclouds')
    endif

endif  !! makeclouds

! *****************************************************************     
! ****************** Update column opacities **********************
! *****************************************************************     

!! Get clouds column and opacity
cst     = .75 / (pi*aerdens(iMa_cld)) * exp( -4.5*dev_ice**2. )
wcldcol(:,:,:) = 0.
taucld(:,:,:) = 0.

do nt=1,nice_mass
    locma=FINDLOC(mom_indx, VALUE=ice_mass_indx(nt))
    locnb=FINDLOC(mom_indx, VALUE=ice_nb_indx(nt))
    ndx_ma=locma(1)
    ndx_nb=locnb(1)

    do i=is,ie
        do j=js,je
            do l = 1, nz
                ilay = 2 * l + 2
                cst2 = (pl(i,j,2*l+3) - pl(i,j,2*l+1)) / grav * 2.0
                Rn   = ( qpi(i,j,ilay,ndx_ma)/( qpi(i,j,ilay,ndx_nb)+1.e-32) &
                *cst )**(athird)
                Rn = max(Rn, 1.e-7)
                Rs   = Rn * dexp( dev_ice**2. )
                taucld(i,j,nt) = taucld(i,j,nt) + pi * Rs**2. * qpi(i,j,ilay,ndx_nb) &
                        * cst2
                wcldcol(i,j,nt) = wcldcol(i,j,nt) + qpi(i,j,ilay,ndx_ma) * &
                     (pl(i,j,2*l+3) - pl(i,j,2*l+1)) / grav
            enddo
        enddo
    enddo
enddo

! *****************************************************************     
! ************************ Outputs ********************************     
! *****************************************************************     
do nt=1,nice_mass
    if (id_taucld(nt) > 0)   used = send_data ( id_taucld(nt), taucld(:,:,nt), Time, is, js )
    if (id_wcol(nt) > 0)   used = send_data ( id_wcol(nt), wcol(:,:,nt), Time, is, js )
    if (id_wcldcol(nt) > 0)   used = send_data ( id_wcldcol(nt), wcldcol(:,:,nt), Time, is, js )
enddo

! some tendencies
do nt=1,ntrace_mom
    if( id_tend_nucl(nt) > 0 )  used = send_data( id_tend_nucl(nt), qdt_nucl(:,:,:,nt), Time, is, js )
enddo
if( id_sat_ratio > 0 )   used = send_data( id_sat_ratio , sat_ratio(:,:,:), Time, is, js )
if( id_mt_3d > 0 )   used = send_data( id_mt_3d , mt_3d(:,:,:), Time, is, js )

return
end subroutine micro_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nucleacond(nz,dt,pl,tl,rho,qpi,qpi_dm,qpi_dn,multi_modes, &
                    tldt,qpidt,qpidt_m,qpidt_n,sat_ratio,mt_out)
!
!                                                              
!     This routine updates species concentrations due          
!     to both nucleation and condensation-induced variations.  
!     Gain and loss rates associated to each one of these      
!     processes are computed separately in other routines.     
!                                                              
!

use constants_mod, only: grav,cp=>cp_air,rgas=>rdgas,pi
use initracer_mod

implicit none

!c  Arguments
!c  ---------

integer, intent(in) :: nz
real*8, intent(in) :: dt                        ! time step
real*8, intent(in), dimension(:) :: pl          ! pressures [Pa]
real*8, intent(in), dimension(:) :: tl          ! temperatures [K]
real*8, intent(in),dimension(:) :: rho          ! atmospheric density (kg/m3)
real*8, intent(in), dimension(:,:) :: qpi       !qpi(2*nz+3,ntrace_mom)
real*8, intent(in), dimension(:,:) :: qpi_dm       !qpi(2*nz+3,ndust_mass)
real*8, intent(in), dimension(:,:) :: qpi_dn       !qpi(2*nz+3,ndust_mass)
logical, intent(in) :: multi_modes              !do multi mode calc
real*8, intent(out), dimension(:) :: tldt       !qpi(2*nz+3,ntrace_mom)
real*8, intent(out), dimension(:,:) :: qpidt    !qpi(2*nz+3,ntrace_mom)
real*8, intent(out), dimension(:,:) :: qpidt_m    !qpi(2*nz+3,ndust_mass)
real*8, intent(out), dimension(:,:) :: qpidt_n    !qpi(2*nz+3,ndust_mass)
real*8, intent(out), dimension(:) :: sat_ratio  ! Water vapor saturation ratio over ice
real*8, intent(out), dimension(:) :: mt_out !runtime contact parameter

!c  Local
!c  -----

integer i,l,nt

real*8 n_aer(nbin)   ! number concentrations of particle/each size bin
real*8 m_aer(nbin)   ! number concentrations of particle/each size bin
real*8 ph2o          ! Water vapor partial pressure (Pa) 
real*8 qsat          ! Water vapor mass mixing ratio at saturation (kg/kg of air)
real*8 rate(nbin)    ! Nucleation rate (s-1)

real*8 qpisav(ntrace_mom)
real*8 qpi_tmp(size(tl,1),ntrace_mom)
real*8   Cste
real*8 p,temp
real*8 Mo,No,dens
real*8 up,dwn,Ctot,gr,rad,seq
real*8 newvap,dN,dM
real*8 Rn,Rm,dev2
real*8   sig

real*8 newT,newS
real*8 Qcond
real*8 :: lw  = 2.8343e+6

integer ilay

real*8   derf

real*8 dqpi

real*8 :: tjs
!multiple modes
real, dimension(nbin,ndust_mass) :: m_aer_tab,rat_tab_m
real, dimension(nbin,ndust_mass) :: n_aer_tab,rat_tab_n
real, dimension(nbin) :: m_cor_tab, n_cld_tab
real :: mtot,mtotnew,ntot,ntotnew
real*8 qpisav_dm(ndust_mass), qpisav_dn(ndust_mass)
real*8 qpi_tmp_dm(size(tl,1),ndust_mass), qpi_tmp_dn(size(tl,1),ndust_mass)
integer :: ndx

!********************************************************
!                     Treatment
!********************************************************
Cste = dt * 4. * pi * dpden_ice
sat_ratio(:) = 0.
mt_out(:) = 0.
qpidt(:,:)=0.d0
qpidt_m(:,:)=0.d0
qpidt_n(:,:)=0.d0
rat_tab_m = 0.d0
rat_tab_n = 0.d0
m_aer_tab = 0.d0
n_aer_tab = 0.d0
m_cor_tab = 0.d0
n_cld_tab = 0.d0
qpisav_dm = 0.d0
qpi_tmp_dm = 0.d0
qpisav_dn = 0.d0
qpi_tmp_dn = 0.d0


!***********************
!Start loop over heights
!***********************
do l = 1, nz
    mtot = 0.d0
    mtotnew = 0.d0
    ntot = 0.d0
    ntotnew = 0.d0

    !! Initializations
    ilay = 2 * l + 2

    p      = pl(ilay)
    temp   = tl(ilay)
    
    !! Save the values of the tracer arrays before condensation 
    do i = 1, ntrace_mom
        qpisav(i) = max((0.d0),qpi(ilay,i))
        qpi_tmp(ilay,i) = qpisav(i)
    enddo
    do i = 1, ndust_mass
        ndx= dust_mass_indx(i)
        if (ndx .le. nt_nontag) then ! do not take into account tags
            qpisav_dm(i) = max((0.d0),qpi_dm(ilay,i))
            qpi_tmp_dm(ilay,i) = qpisav_dm(i)
            qpisav_dn(i) = max((0.d0),qpi_dn(ilay,i))
            qpi_tmp_dn(ilay,i) = qpisav_dn(i)
        endif
    end do
    mtot = sum(qpisav_dm(:)) + qpisav(iMa_cor)
    ntot = sum(qpisav_dn(:)) + qpisav(iNb_cld)

    !! get qsat
    call watsat(temp,p,qsat)

    !! Get the partial presure of water vapor and its saturation ratio
    ph2o      = qpi_tmp(ilay,iMa_vap) * (44./18.) * p
    sat_ratio(l) = qpi_tmp(ilay,iMa_vap) / qsat

    !================ acocunt for multiple dust modes if more than one ==========
    if (multi_modes) then
        !! Expand each dust moment into single binned distribution
        !!! Radius and Volumes defined in initracer/initmicro
        ! rads_coag, vols_coag, vrat_coag

        !!! Initial tracers
        n_aer_tab(:,:)=0.d0

        do nt=1,ndust_mass
            ndx= dust_mass_indx(nt)
            if (ndx .le. nt_nontag) then ! do not take into account tags
                Mo = qpi_tmp_dm(ilay,nt)
                No = qpi_tmp_dn(ilay,nt)
                if (No.gt.0.d0) then
                    Rn = ( Mo / No * 0.75 / pi / aerdens(iMa_dt) )**(athird) &
                        * exp( -1.5 * stdv(iMa_dt)**2.)
                else
                    Rn = 0.
                endif
                Rn = min( max(Rn,aerad(1)) , aerad(nbin) )
                Rm = Rn * exp( 3 * stdv(iMa_dt)**2. )
                Rn = 1. / Rn
                Rm = 1. / Rm
                dev2 = 1 / ( sqrt(2.) * stdv(iMa_dt) )

                do i=1,nbin ! Sum of the ndis for all dust modes
                    n_aer_tab(i,nt) = 0.5 * No * ( derf( dlog(rb(i+1)*Rn) * dev2 ) &
                            -derf( dlog(rb(i) * Rn) * dev2 ) )
                    m_aer_tab(i,nt) = 0.5 * Mo * ( derf( dlog(rb(i+1)*Rm) * dev2 ) &
                            -derf( dlog(rb(i) * Rm) * dev2 ) )
                enddo
            else
                exit
            endif
        enddo

        do nt=1,ndust_mass
            ndx= dust_mass_indx(nt)
            if (ndx .le. nt_nontag) then ! do not take into account tags
                ! Initial Ratio for each mode
                rat_tab_n(:,nt)=n_aer_tab(:,nt)/sum(n_aer_tab(:,:),2)
                rat_tab_m(:,nt)=m_aer_tab(:,nt)/sum(m_aer_tab(:,:),2)
            else
                exit
            endif
        enddo

        ! Normalization of the ratio (make sure sum=1)
        do nt=1,ndust_mass
             rat_tab_n(:,nt)=rat_tab_n(:,nt)/sum(rat_tab_n(:,:),2)
             rat_tab_m(:,nt)=rat_tab_m(:,nt)/sum(rat_tab_m(:,:),2)
        enddo
        ! Total initial distribution
        n_aer(:)=sum(n_aer_tab(:,:),2)
        m_aer(:)=sum(m_aer_tab(:,:),2)
    else
        !! Expand the dust moments into a binned distribution
        Mo = qpi_tmp(ilay,iMa_dt)
        No = qpi_tmp(ilay,iNb_dt)
        if (No.gt.0.d0) then
            Rn = ( Mo / No * 0.75 / pi / aerdens(iMa_dt) )**(athird) &
                * exp( -1.5 * stdv(iMa_dt)**2.)
        else
            Rn = 0.
        endif
        Rn = min( max(Rn,aerad(1)) , aerad(nbin) )
        Rm = Rn * exp( 3 * stdv(iMa_dt)**2. )
        Rn = 1. / Rn
        Rm = 1. / Rm
        dev2 = 1 / ( sqrt(2.) * stdv(iMa_dt) )

        do i = 1, nbin   
            n_aer(i) = 0.5 * No * ( derf( dlog(rb(i+1)*Rn) * dev2 ) &
                            -derf( dlog(rb(i) * Rn) * dev2 ) )
            m_aer(i) = 0.5 * Mo * ( derf( dlog(rb(i+1)*Rm) * dev2 ) &
                            -derf( dlog(rb(i) * Rm) * dev2 ) )
        enddo
    end if

    !! Get the rates of nucleation. =0 if sat_ratio <= 1
    call nuclea(ph2o,temp,sat_ratio(l),n_aer,rate,mt_out(l))

    if (multi_modes) then
        do nt = 1, ndust_mass
            dN = 0.d0
            dM = 0.d0
            n_aer_tab(:,nt) = rat_tab_n(:,nt)*n_aer_tab(:,nt) / ( 1. + rate(:)*dt)
            m_aer_tab(:,nt) = rat_tab_m(:,nt)*m_aer_tab(:,nt) / ( 1. + rate(:)*dt)
            dN = sum(n_aer_tab(:,nt)*rate(:)*dt)
            dM = sum(m_aer_tab(:,nt)*rate(:)*dt)
            if (nt.eq.1) then
                !! Update original qpi tracers
                qpi_tmp(ilay,iMa_dt ) = qpi_tmp(ilay,iMa_dt ) - dM
                qpi_tmp(ilay,iNb_dt ) = qpi_tmp(ilay,iNb_dt ) - dN
            end if
            qpi_tmp_dm(ilay,nt) = qpi_tmp_dm(ilay,nt) - dM
            qpi_tmp_dn(ilay,nt) = qpi_tmp_dn(ilay,nt) - dN
            qpi_tmp(ilay,iMa_cor) = qpi_tmp(ilay,iMa_cor) + dM
            qpi_tmp(ilay,iNb_cld) = qpi_tmp(ilay,iNb_cld) + dN
        end do
    else
        dN = 0.d0
        dM = 0.d0
        do i = 1, nbin
            n_aer(i) = n_aer(i) / ( 1. + rate(i)*dt )
            m_aer(i) = m_aer(i) / ( 1. + rate(i)*dt )
            dN       = dN + n_aer(i) * rate(i) * dt
            dM       = dM + m_aer(i) * rate(i) * dt
        enddo

        !! Update tracers
        qpi_tmp(ilay,iMa_dt ) = qpi_tmp(ilay,iMa_dt ) - dM
        qpi_tmp(ilay,iNb_dt ) = qpi_tmp(ilay,iNb_dt ) - dN
        qpi_tmp(ilay,iMa_cor) = qpi_tmp(ilay,iMa_cor) + dM
        qpi_tmp(ilay,iNb_cld) = qpi_tmp(ilay,iNb_cld) + dN
        qpi_tmp_dm(ilay,1) = qpi_tmp_dm(ilay,1) - dM
        qpi_tmp_dn(ilay,1) = qpi_tmp_dn(ilay,1) - dN
    endif

    Qcond = 0.d0
    if (qpi_tmp(ilay,inb_cld).ge.1.e-32.and.sat_ratio(l).ne.1.) then

        !! get radius cloud particle
        Mo   = qpisav(iMa_cld) + qpi_tmp(ilay,iMa_cor)
        if (Mo.lt.1.e-15) then
            rad = 1.e-8
        else
            No   = qpi_tmp(ilay,iNb_cld)
            dens = qpisav(iMa_cld) / Mo * aerdens(iMa_cld) &
                +qpi_tmp(ilay,iMa_cor) / Mo * aerdens(iMa_dt)
            dens = max(dens,dpden_ice)
            rad  = ( Mo / No * 0.75 / pi / dens ) **(athird)
        endif

        sig = (141. - 0.15 * temp) * 1.e-3
        seq  = exp( 2.*sig*mh2o / (dpden_ice*rgp*temp*rad) )

        tjs = ph2o/sat_ratio(l)
        call growthrate2(dt,temp,p,ph2o,tjs,seq,rad,gr)

        up  = Cste * gr * rad * No * seq + qpi_tmp(ilay,iMa_vap)
        dwn = Cste * gr * rad * No / qsat+ 1.

        Ctot = qpisav(iMa_cld) + qpi(ilay,iMa_vap)

        newvap = min(up/dwn,Ctot)

        gr = gr * ( newvap/qsat - seq )
        Qcond = max( Cste * No * rad * gr , -qpisav(iMa_cld) )

        !!**************************************
        if (latent_heat) then
            newT = temp + Qcond * lw / cp_air
            call watsat(newT,p,qsat)
            newS = ( qpi_tmp(ilay,iMa_vap) - Qcond ) / qsat
            if (sat_ratio(l).lt.1..and.newS.gt.1.01) then
                call findQ(qpi_tmp(ilay,iMa_vap),p,temp,Qcond)
            elseif(sat_ratio(l).gt.1..and.newS.lt.0.99) then
                call findQ(qpi_tmp(ilay,iMa_vap),p,temp,Qcond)
            endif
        endif
        !!**************************************

        qpi_tmp(ilay,iMa_cld) = qpisav(iMa_cld) + Qcond

        !! Check if cloud entirely sublimed 
        if (qpi_tmp(ilay,iMa_cld).le.0.d0) then
            if (multi_modes) then
                Mo = qpi_tmp(ilay,iMa_cor)
                No = qpi_tmp(ilay,iNb_cld)
                if (No.gt.0.d0) then
                    Rn = ( Mo / No * 0.75 / pi / aerdens(iMa_dt) )**(athird) &
                        * exp( -1.5 * stdv(iMa_dt)**2.)
                else
                    Rn = 0.
                endif
                Rn = min( max(Rn,aerad(1)) , aerad(nbin) )
                Rm = Rn * exp( 3 * stdv(iMa_dt)**2. )
                Rn = 1. / Rn
                Rm = 1. / Rm
                dev2 = 1 / ( sqrt(2.) * stdv(iMa_dt) )
                do i = 1, nbin   
                    n_cld_tab(i) = 0.5 * No * ( derf( dlog(rb(i+1)*Rn) * dev2 ) &
                                    -derf( dlog(rb(i) * Rn) * dev2 ) )
                    m_cor_tab(i) = 0.5 * Mo * ( derf( dlog(rb(i+1)*Rm) * dev2 ) &
                                    -derf( dlog(rb(i) * Rm) * dev2 ) )
                enddo
                do nt = 1,ndust_mass
                    qpi_tmp_dm(ilay,nt) = qpi_tmp_dm(ilay,nt) + sum(m_cor_tab(:)*rat_tab_m(:,nt))
                    qpi_tmp_dn(ilay,nt) = qpi_tmp_dn(ilay,nt) + sum(n_cld_tab(:)*rat_tab_n(:,nt))
                end do
                qpi_tmp(ilay,iMa_dt ) = qpi_tmp_dm(ilay,1)
                qpi_tmp(ilay,iNb_dt ) = qpi_tmp_dn(ilay,1)
            else
                qpi_tmp(ilay,iMa_dt ) = qpi_tmp(ilay,iMa_dt) + qpi_tmp(ilay,iMa_cor)
                qpi_tmp(ilay,iNb_dt ) = qpi_tmp(ilay,iNb_dt) + qpi_tmp(ilay,iNb_cld)
            endif
                qpi_tmp(ilay,iMa_cld) = 0.d0
                qpi_tmp(ilay,iMa_cor) = 0.d0
                qpi_tmp(ilay,iNb_cld) = 0.d0
        endif

    endif

    qpi_tmp(ilay,iMa_vap) = qpisav(iMa_vap) - (qpi_tmp(ilay,iMa_cld)-qpisav(iMa_cld))

    !! Update tendencies tracers       
    do i = 1, ntrace_mom
        qpidt(ilay,i) = max(-qpisav(i),(qpi_tmp(ilay,i)-qpisav(i)))/dt
    enddo 
    qpidt_m(ilay,:) = max(-qpisav_dm(:),(qpi_tmp_dm(ilay,:)-qpisav_dm(:)))/dt
    qpidt_n(ilay,:) = max(-qpisav_dn(:),(qpi_tmp_dn(ilay,:)-qpisav_dn(:)))/dt

    mtotnew = sum(qpi_tmp_dm(ilay,:)) + qpi_tmp(ilay,iMa_cor)
    ntotnew = sum(qpi_tmp_dn(ilay,:)) + qpi_tmp(ilay,iNb_cld)
    !! Update Temperature       
    if (latent_heat) then
        tldt(ilay)= Qcond * lw / cp_air / dt !! lw depending on T ? TB18q
    else
        tldt(ilay)=0.d0
    endif

enddo

return
end subroutine nucleacond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nuclea(ph2o,temp,sat,n,nucrate,mnew)
!
!                                                     
!   This subroutine computes the nucleation rate      
!   as given in Pruppacher & Klett (1978) in the      
!   case of water ice forming on a solid substrate.   
!     Definition refined by Keese (jgr,1989)          
!                                                     
!

implicit none

integer i

real*8 n(nbin)

real*8 nucrate(nbin)
real*8 ph2o,temp,sat

real*8 nh2o
real*8 sig          ! Water-ice/air surface tension  (N.m)
real*8 rstar        ! Radius of the critical germ (m)
real*8 gstar        ! # of molecules forming a critical embryo
real*8 x            ! Ratio rstar/radius of the nucleating dust particle
real*8 fistar       ! Activation energy required to form a critical embryo (J)
real*8 zeldov       ! Zeldovitch factor (no dim)
real*8 fshape       ! function defined at the end of the file
real*8 sinteta      ! sine of the contact angle
real*8 deltaf

real*8 scrit        ! temperature dependent critical saturation ratio for nucleation
real*8 mnew         ! temperature dependent contact parameter


call scrit_tdep(temp,scrit,mnew)

if (sat .gt. scrit) then    ! minimum condition to activate nucleation

    nh2o    = ph2o / kbz / temp
    sig = (141. - 0.15 * temp) * 1.e-3
    rstar  = 2. * sig * vo1 / (rgp*temp*dlog(sat))
    gstar  = 4. * nav * pi * (rstar**3) / (3.*vo1)

    !c       Loop over size bins
    do i=1,nbin

        if ( n(i) .eq. 0. ) then  ! no dust, no need to compute nucleation
            nucrate(i)=0.
            cycle
        endif

        x      = aerad(i) / rstar
        sig = (141. - 0.15 * temp) * 1.e-3
        call fshape2(mnew,x,fshape)
        fistar = (4./3.*pi) * sig * (rstar**2.) * fshape
        deltaf = min( max((2.*desorp-surfdif-fistar)/(kbz*temp) &
            , -100.), 100.)

        if (deltaf.eq.-100.) then
            nucrate(i) = 0.
        else
            zeldov = sqrt ( fistar / (3.*pi*kbz*temp*(gstar**2.)) )
            call fshape2(mnew,x,fshape)
            nucrate(i)= zeldov * kbz* temp * rstar &
                 * rstar * 4. * pi * ( nh2o*aerad(i) )**2. &
                 / ( fshape * nus * m0 ) &
                 * dexp (deltaf)
        endif

    end do

else ! sat <= 1

    do i=1,nbin
        nucrate(i) = 0.
    enddo

endif

return
end subroutine nuclea

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine growthrate2(timestep,t,p,ph2o,psat,seq,r,Cste)
!
!
!     Determination of the water ice crystal growth rate
!
!


IMPLICIT NONE
!
!   arguments:
!   

real*8 timestep
real*8 t    ! temperature in the middle of the layer (K)
real*8 p    ! pressure in the middle of the layer (Pa)
real*8 ph2o ! water vapor partial pressure (Pa)
real*8 psat ! water vapor saturation pressure (Pa)
real*8 r    ! crystal radius before condensation (m)
real*8 seq  ! Equilibrium saturation ratio
real*8 dr   ! crystal radius variation (m)

!   local:
!   ------

real*8 :: molco2 = 2.2e-10    !     Effective CO2 molecular radius (m)

real*8 :: molh2o = 1.2e-10    !     Effective H2O molecular radius (m)

real*8 :: Mco2 = 44.e-3       !     Molecular weight of CO2 kg.mol-1

real*8 :: Mh2o = 18.e-3       !     Molecular weight of H2O kg.mol-1

real*8 :: sigh2o = 0.12       !     surface tension of ice/vapor N.m

real*8 :: rho_i = 917.        !     Ice density kg.m-3

real*8 :: nav = 6.023e23      !     Avogadro number

real*8 :: rgp = 8.3143        !     Perfect gas constant

real*8 :: To = 273.15         !     Reference temperature, T=273,15 K

real*8 k,Lv
real*8 knudsen           ! Knudsen number (gas mean free path/particle radius)
real*8 a,Dv,lambda,Rk,Rd ! Intermediate computations for growth rate
real*8 Cste, rf

!-----------------------------------------------------------------------
!      Ice particle growth rate by diffusion/impegement of water molecules
!                r.dr/dt = (S-Seq) / (Seq*Rk+Rd)
!        with r the crystal radius, Rk and Rd the resistances due to
!        latent heat release and to vapor diffusion respectively
!-----------------------------------------------------------------------

!     - Equilibrium saturation accounting for KeLvin Effect
!       seq=exp(2*sigh2o*Mh2o/(rho_i*rgp*t*r))

!     - Thermal conductibility of CO2
k  = (0.17913 * t - 13.9789) * 4.184e-4
!     - Latent heat of h2o (J.kg-1)
Lv = (2834.3 - 0.28 * (t-To) - 0.004 * (t-To)**2 ) * 1.e+3

!     - Constant to compute gas mean free path
!     l= (T/P)*a, with a = (  0.707*8.31/(4*pi*molrad**2 * avogadro))
a = 0.707*rgp/(4 * pi* molco2**2  * nav)

!     - Compute Dv, water vapor diffusion coefficient
!       accounting for both kinetic and continuum regime of diffusion,
!       the nature of which depending on the Knudsen number.

Dv = 1./3. * sqrt( 8*kbz*t/(pi*Mh2o/nav) )* kbz * t / &
( pi * p * (molco2+molh2o)**2 * sqrt(1.+Mh2o/Mco2) )

knudsen = t / p * a / r
lambda  = (1.333+0.71/knudsen) / (1.+1./knudsen)
Dv      = Dv / (1. + lambda * knudsen)

!     - Compute Rk
Rk = Lv**2 * rho_i * Mh2o / (k*rgp*t**2.)
!     - Compute Rd
Rd = rgp * t *rho_i / (Dv*psat*Mh2o)

!     - Compute Cste=rdr/dt, then r(t+1)= sqrt(r(t)**2.+2.*Cste*dt)
Cste = 1. / (seq*Rk+Rd)

return
end subroutine growthrate2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fshape2(cost,rap,fshape)
!
!        function computing the f(m,x) factor           
! related to energy required to form a critical embryo  
!

implicit none

real*8 cost,rap,fshape
real*8 phi,a,b,c

phi = sqrt( 1. - 2.*cost*rap + rap**2 )
a = 1. + ( (1.-cost*rap)/phi )**3
b = (rap**3) * (2.-3.*(rap-cost)/phi+((rap-cost)/phi)**3)
c = 3. * cost * (rap**2) * ((rap-cost)/phi-1.)

fshape = 0.5*(a+b+c)

if (rap.gt.3000.) fshape = ((2.+cost)*(1.-cost)**2)/4.

return
end subroutine fshape2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine findQ(Qvap,p,temp,Qcond)
!
!     Start moddfied Newton-Raphson method to find             
!     the amount of condensed/sublimed ice required to         
!     reach saturation. Iterations are made to                 
!     solve the F(Qc)=S-1=0 equation (Qc=condensed mass,       
!     Saturation ratio) in the case of latent heat             
!     release associated to Qc..                               
!
implicit none

real*8 Qvap,p,temp,Qcond

real*8 x1,x2,xl,xh,f,df,fl,fh
real*8 dx,dxold,tempo
real*8 rtsafe

integer i

x1 = 0.
x2 = Qcond
call newton(x1,Qvap,p,temp,fl,df)
call newton(x2,Qvap,p,temp,fh,df)

if (fl*fh.ge.0.) then
    print*,'root not bracketed'
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
call newton(rtsafe,Qvap,p,temp,f,df)

do i = 1, 500
    if ( ((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0. &
        .or. abs(2.*f).gt.abs(dxold*df) ) then
        dxold = dx
        dx    = 0.5 * (xh-xl)
        rtsafe= xl+dx
        if (xl.eq.rtsafe) exit
    else
        dxold = dx
        dx    = f / df
        tempo  = rtsafe
        rtsafe= rtsafe - dx
        if (tempo.eq.rtsafe) exit
    endif
    if (abs(dx).lt.1.e-8) exit
    call newton(rtsafe,Qvap,p,temp,f,df)
    if (f.lt.0) then
        xl = rtsafe
    else
        xh = rtsafe
    endif
    if (i.eq.500) print*,'500 reached'
enddo

Qcond = rtsafe

return
end subroutine findQ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine newton(dq,q,press,temp,f,df)
!
!     subroutine called during the Newton-Raphson loop  
!     to get the function and its derivative for the    
!              condensed mass determination             
!

implicit none

real*8 dq,q,press,temp
real*8 newT
real*8 qsat

real*8 :: lw  = 2.8e+6
save lw

real*8 f,df

newT = temp + dq * lw / cp_air
newT = max(newT,60.)

call watsat(newT,press,qsat)
qsat = max(qsat,1.d-50)

f    = (q-dq) / qsat - 1.
df   = - 2. / qsat - (f+1.) * (6146.1*lw/cp_air) / newT**2.

return
end subroutine newton

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine watsat(temp,press,qsat)
!
!  subroutine to calculate saturation water vapor pressure
!
real*8 pvs,temp,press,qsat

pvs  = 611.0*exp(22.5*(1.0-(273.16/temp)))

qsat = pvs * 18.0 / (44.0*press)

return
end subroutine watsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine micro_driver_init( nlon, mlat, nlevels, lonb, latb, lon, lat, axes, Time )
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

! *********************************************************
!     ----- read namelist /_nml/   -----
! *********************************************************

if (file_exists('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1; do while (ierr /= 0)
        read  (unit, nml=microphys_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'microphys_nml')
    enddo
10     call close_file (unit)
endif

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=microphys_nml)

! *********************************************************
!   --- Allocate T + Q fields needed for cloud scheme ---
! *********************************************************
allocate (  tlprev(is:ie,js:je,2*nlevels+3)  )
allocate (  qpiprev(is:ie,js:je,2*nlevels+3,ntrace_mom)  )
allocate (  qpi_dmprev(is:ie,js:je,2*nlevels+3,ndust_mass)  )
allocate (  qpi_dnprev(is:ie,js:je,2*nlevels+3,ndust_mass)  )

allocate (  wcldcol(is:ie,js:je,nice_mass) )
allocate (  wcol(is:ie,js:je,nice_mass) )

! *********************************************************
!     ----- register diagnostic fields -----
! *********************************************************
allocate ( id_tend_nucl(ntrace_mom) )

allocate ( id_taucld(nice_mass) )
allocate ( id_wcol(nice_mass) )
allocate ( id_wcldcol(nice_mass) )

!! Preparing output
DO nt= 1, ntrace_mom
    !Get name tracer
    ndx= mom_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)
    if (mpp_pe() == mpp_root_pe()) print*,'Micro Tracer ID:',nt,ndx,tracer_name

    tname= trim(tracer_name) // '_Tnucl' 
    id_tend_nucl(nt) = register_diag_field ( mod_name, trim(tname),  &
         axes(1:3), Time, 'Tendency nucleation', 'kg/kg/s', &
         missing_value=missing_value)
ENDDO

DO nt= 1, nice_mass
    !Get name tracer
    ndx= ice_mass_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)

    tname= trim(tracer_name) // '_tau' 
    id_taucld(nt) = register_diag_field ( mod_name, trim(tname),  &
                                 (/axes(1:2)/), Time,           &
                                'cloud tau micro ', '',        &
                                 missing_value=missing_value )

    tname= trim(tracer_name) // '_col'
    id_wcldcol(nt) = register_diag_field ( mod_name, trim(tname),  &
                                 (/axes(1:2)/), Time,           &
                                'water ice column micro ', '',        &
                                 missing_value=missing_value )

    ndx= vapor_indx(nt)
    call get_tracer_names(MODEL_ATMOS, ndx, tracer_name)

    tname= trim(tracer_name) // '_col'
    id_wcol(nt) = register_diag_field ( mod_name, trim(tname),  &
                                 (/axes(1:2)/), Time,           &
                                'water vapor column micro ', '',        &
                             missing_value=missing_value )
ENDDO

id_sat_ratio = register_diag_field ( mod_name, 'sat_ratio',  &
           axes(1:3), Time, 'Vap saturation ratio', '', &
           missing_value=missing_value)

id_mt_3d = register_diag_field ( mod_name, 'mt_3d',  &
           axes(1:3), Time, '3D Contact Parameter', '', &
           missing_value=missing_value)

end subroutine micro_driver_init

!********************************************************
subroutine scrit_tdep(temp,scrit,m_tdep)
!     subroutine called during nucleation to find temp  *
!     dependent critical saturation ratio, and contact  *
!     parameter                                         *
!********************************************************

implicit none

real*8 temp,scrit,m_tdep


select case(mteta_case)
    case (1)
        !-----------------------------------------------------------------------
        !  Values for Arizona Test Dust Case (Iraci)
        !
        scrit = -0.0718*temp + 14.4
        m_tdep = 0.0046*temp + 0.1085
    case (2)
        !-----------------------------------------------------------------------
        !  Values for Trainer et al. (2009)
        !
        m_tdep = 0.94-6005.*exp(-0.065*temp)
        scrit = 1.0
    case (3)
        !-----------------------------------------------------------------------
        !  Values for Maattanen et al. (2014)
        !
        scrit = 1.0
        m_tdep = -0.698 + (1.635*tanh((temp/144.21)**3.239))
    case default
        scrit = 1.0
        m_tdep = mteta
end select

return
end subroutine scrit_tdep


end module micromod_mgcm

