module ames_rt_interface
! module to interface FV3 with legacy radiative transfer
use rtmod_mgcm, only: ames_rt_driver, qextv, l_nrefv, l_nspecti, scale_by_scon
use field_manager_mod,  only: MODEL_ATMOS, parse, find_field_index
use aerosol_util_mod, only: do_moment_dust, do_bulk_water
use astronomy_mod, only: semi_major_axis

implicit none
private
public :: ames_rt


contains


subroutine ames_rt(is,js,id,jd,kd,ntrace,p_half,p_full,             &
               t,tsurf,r,trans,flx_sfc,                             &
               albedo,sfc_emiss,coszro,                             &
               dustref,dustref_bin,dustref_fix,                     &
               cldref,cldco2ref,cldicebin,                          &
               dosw, dolw,                                          &
               rorbit,                                              &
               heatra,hsw,out_solar_flx,                            &
               rsolar,                                              &
               irupflx,irdnflx,swupflx,swdnflx,swnetflx,irnetflx,   &
               taudust,taucloud,tauco2cloud,taudust_mom,            &
               lw_heating_band,lw_15umHR,diag,tstrat_in,            &
               tstrat_dt,                                           &
               taudust_reff,taudust_fix,                            &
               tbands          )

! main interface routine
implicit none
integer :: is, js, nz,id,jd,kd,ntrace

real,intent(in),    dimension(id,jd,kd)   :: p_full        ! layer midpoint pressures [Pa]
real,intent(in),    dimension(id,jd,kd+1) :: p_half        ! layer interface pressures [Pa]
real,intent(in),    dimension(id,jd,kd)   :: t             ! layer midpoint temperatures [K]
real,intent(in),    dimension(id,jd,kd,ntrace) :: r        ! tracers [*/kg]
real,intent(inout), dimension(id,jd,kd)   :: dustref       ! output dust opacity from RT
real,intent(in),    dimension(id,jd,kd)   :: dustref_fix   ! input fixed dust opacity
real,intent(in),    dimension(id,jd,kd)   :: dustref_bin   ! input bin dust opacity
real,intent(inout), dimension(id,jd,kd,2)   :: cldref        ! output cloud opacity from RT
real,intent(inout), dimension(id,jd,kd)   :: cldco2ref     ! output co2 cloud opacity from RT
real, intent(in),   dimension(id,jd,kd)   :: cldicebin     ! input bin cloud opacity
real,intent(in),    dimension(id,jd)      :: albedo, &     ! surface albedo
                                            sfc_emiss, &   ! surface emissivity
                                            coszro         ! cosine of zenith angle
real,intent(in),    dimension(id,jd)      :: tstrat_in     ! input tstrat [K]
real,intent(out),    dimension(id,jd)     :: tstrat_dt     ! output tstrat tendency [K/s]
logical,intent(in) :: dosw, &                              ! do visible calc
                      dolw, &                              ! do IR calc
                      diag                                 ! print diagnostic comments
real,intent(in), dimension(size(t,1),size(t,2)) :: tsurf      ! surface temperature [K]
real,intent(out), dimension(size(t,1),size(t,2)) :: trans, &  ! transmitted radiation [W/m^2]
                                                flx_sfc, &    ! downward IR at surface [W/m^2]
                                                out_solar_flx ! outgoing visible radiation [W/m^2]
real,intent(out), dimension(size(t,1),size(t,2),size(t,3)) :: heatra, & ! IR heating rate [K/s]
                                                            hsw         ! visible heating rate [K/s]
real,intent(out), dimension(size(t,1),size(t,2),size(t,3)+1) :: irupflx, &  ! IR upwards flux [W/m^2]
                                                            irdnflx, &      ! IR downwards flux [W/m^2]
                                                            swupflx, &      ! visible upwards flux [W/m^2]
                                                            swdnflx, &      ! visible downwards flux [W/m^2]
                                                            swnetflx, &     ! visible net flux [W/m^2]
                                                            irnetflx        ! IR net flux [W/m^2]

real,intent(out), dimension(size(t,1),size(t,2),3) :: tbands            !     brightness temperatures:  t7 t23 t32


!   Opacities TB18c
real,intent(inout), dimension(size(t,1),size(t,2),2) :: taudust, &      ! total column dust opacity
                                                    tauco2cloud, &     ! total column co2 cloud opacity
                                                    taudust_mom, &     ! moment dust column dust opacity
                                                    taudust_fix     ! fixed dust column dust opacity
real,intent(inout), dimension(size(t,1),size(t,2),4) :: taucloud     ! total column h2o cloud opacity
real,intent(out), dimension(size(r,1),size(r,2),size(r,4),2) :: taudust_reff  ! moment column dust opacity by effective radius
real, intent(in) :: rsolar, &        ! total solar constant from FV3
                    rorbit          ! Mars-Sun distance

!===============================================================
!     Local Variables
!===============================================================
real,    dimension(kd)         :: tauscale_bin   !rescaled opd to qextref_ames
real,    dimension(kd)         :: tauscale_fix
real,    dimension(kd)         :: cldicecol
real, parameter :: qext_fv3 = 2.75   !reference extinction in FV3. Should eventually be changed?
real :: rsdist

integer :: nma_vap,nma_vap2,nma_cor,nma_dst,nnb_dst,nnb_cld,nma_cld

real, dimension(2*size(t,3)+3):: QH2O
real, dimension(size(t,1),size(t,2),size(t,3),L_NSPECTI) :: lw_heating_band
real, dimension(size(t,1),size(t,2),size(t,3)) :: lw_15umHR
real, dimension(L_NSPECTI,size(t,3)) :: lw_heating_spec
real, dimension(size(t,3)) :: htrt2

real, dimension(2*size(t,3)+3) :: tl        !   Absolute Temperature (K) at boundaries and midpoints
real, dimension(2*size(t,3)+3) :: pl        !   Pressure (Pa) at boundaries and midpoints

real, dimension(size(t,3)+1) :: fmnetv,fluxdni,fmneti
real :: albi
real :: swfactor,sfmars

integer :: ie, je
integer :: i,j,k,n,n2
integer :: nh2o
integer :: nice_blk
integer :: nco2
real :: dnvsol
logical :: do_tstrat = .false.
real :: tstrat_pass
real, dimension(size(r,4),2) :: taudust_tmp

real :: mwratio
real, parameter :: mwco2 = 4.41D-2
real, parameter :: mwh2o = 1.80153D-2
    
ie = size(t,1)
je = size(t,2)
nz = size(t,3)

hsw = 0.
heatra = 0.
lw_heating_band = 0.
lw_15umHR = 0.
out_solar_flx = 0.
swupflx = 0.
swdnflx = 0.
irupflx = 0.
irdnflx = 0.
swnetflx = 0.
irnetflx = 0.
flx_sfc = 0.
taudust = 0.
taucloud = 0.
tauco2cloud = 0.
taudust_reff = 0.
mwratio = mwco2/mwh2o
taudust_fix = 0.

! --Set no dust case --
!       dustref = 0.
if (scale_by_scon) then
    rsdist = 1356.1/rsolar      !This is to scale the solar flux to the GFDL value so that SOL/rsdist = GFDL band-split solar
else
    rsdist = (semi_major_axis*rorbit)**2  !This is so that rsdist is exactly the mars-sun distance squared, as expected by the Ames RT
endif

if (do_moment_dust) then
    nh2o= find_field_index( MODEL_ATMOS, 'vap_mass_mom' )
    nma_vap= find_field_index( MODEL_ATMOS, 'vap_mass_mom' )
    nma_cld= find_field_index( MODEL_ATMOS, 'ice_mass_mom' )
    nma_dst= find_field_index( MODEL_ATMOS, 'dst_mass_mom' )
    nma_cor= find_field_index( MODEL_ATMOS, 'cor_mass_mom' )
    nnb_dst= find_field_index( MODEL_ATMOS, 'dst_num_mom' )
    nnb_cld= find_field_index( MODEL_ATMOS, 'ice_num_mom' )
elseif (do_bulk_water) then
    nh2o= find_field_index( MODEL_ATMOS, 'vap_mass_blk' )
    nice_blk = find_field_index( MODEL_ATMOS, 'ice_mass_blk' )
else
    nh2o= find_field_index( MODEL_ATMOS, 'h2o_vapor' )
endif

! Get the index of the CO2 cloud tracer
nco2= find_field_index( MODEL_ATMOS, 'co2_cloud' )

do i=1,id
    do j=1,jd

        fmnetv = 0.
        fluxdni = 0.
        sfmars = 0.
        fmneti = 0.
        pl = 0.
        tl = 0.
        albi = 0.
        tauscale_bin = 0.
        tauscale_fix = 0.
        taudust_tmp = 0.

        pl(3) = p_half(i,j,1)
        pl(2) = pl(3)*0.5
        pl(1) = pl(2)*1.e-5

        do k= 1, nz
            n= 2*k + 2
            tl(n)= t(i,j,k)
            pl(n)= p_full(i,j,k)
            pl(n+1)=p_half(i,j,k+1)
            if (nh2o .ge. 1) then
                qh2o(n)=mwratio*r(i,j,k,nh2o)
            else
                qh2o(n)=0.d0
            endif
            if (k.eq.nz) then
                tl(n+1) = tsurf(i,j)
            else
                tl(n+1)= 0.5*( t(i,j,k+1)+t(i,j,k) )
            end if
        enddo


        tl(1:3)=t(i,j,1)
        albi = 1.-sfc_emiss(i,j)

        tauscale_bin(:) = max(0.,dustref_bin(i,j,:))  !no rescaling needed
        tauscale_fix(:) = max(0.,dustref_fix(i,j,:))  !no rescaling needed
	    cldicecol(:)= max( 0., cldicebin(i,j,:) )    !no rescaling needed
        if (do_tstrat) then
            tstrat_pass = tstrat_in(i,j)
        else
            tstrat_pass = tl(3)
        endif

        call ames_rt_driver(pl=pl,tl=tl,ptop=pl(3),                 &
                   tstrat=tstrat_pass,                              &
                   tsurf=tsurf(i,j),QH2Oin=QH2O(:),                 &
                   acosz=coszro(i,j), rsdist=rsdist,                &
                   albv=albedo(i,j), albi=albi,                     &
                   taud_bin_in=tauscale_bin(:),                     &
                   taud_fix_in=tauscale_fix(:),                     &
                   dosw=dosw, dolw=dolw,                            &
                   sw_heating=hsw(i,j,:),                           &
                   lw_heating=heatra(i,j,:),                        &
                   qtrace=r(i,j,:,:),cldice_bin= cldicecol,             &
                   taudust_diagARG=taudust(i,j,:),                  &
                   taucloud_diagARG=taucloud(i,j,:),                &
                   tauco2cloud_diagARG=tauco2cloud(i,j,:),                &
                   fmnetvARG=fmnetv,fluxdniARG=fluxdni,             &
                   outsolARG=out_solar_flx(i,j),                    &
                   swtotARG=sfmars,fluxupvARG=swupflx(i,j,:),       &
                   fluxdnvARG=swdnflx(i,j,:),                       &
                   fluxupiARG=irupflx(i,j,:),fmnetiARG=fmneti,      &
                   lw_heating_spec=lw_heating_spec,htrt2=htrt2,     &
                   diag=diag,taurefd_out=dustref(i,j,:),            &
                   taurefc_out=cldref(i,j,:,:),                       &
                   taurefco2c_out=cldco2ref(i,j,:),                 &
                   taudust_momARG=taudust_mom(i,j,:),               &
                   taudust_fixARG=taudust_fix(i,j,:),               &
                   tstrat_dt=tstrat_dt(i,j),                        &
                   nco2=nco2,                                       &
                   nice_blk=nice_blk,                               &
                   taudust_reffARG=taudust_tmp(:,:),                &
                   tbands= tbands(i,j,:)                )



        do k= 1, nz
            do n=1,l_nspecti !# ir spectral bands
                lw_heating_band(i,j,k,n) = lw_heating_spec(n,k)
            enddo
            lw_15umhr(i,j,k) = htrt2(k)
        enddo


        trans(i,j)   = -fmnetv(nz+1)
        hsw(i,j,:)   = hsw(i,j,:)
        flx_sfc(i,j) = fluxdni(nz+1)
        irdnflx(i,j,:) = fluxdni(:)
        swdnflx(i,j,:) = swdnflx(i,j,:)
        swupflx(i,j,:) = swupflx(i,j,:)
        swnetflx(i,j,:)= fmnetv(:)
        irnetflx(i,j,:)= fmneti(:)
        taudust_reff(i,j,:,:) = taudust_tmp(:,:)

    end do
end do


end subroutine ames_rt



end module
