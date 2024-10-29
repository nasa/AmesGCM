module rtmod_mgcm
!  Ames legacy Radiation Transfer

use constants_mod, only: grav, cp_air, pi
! The following is required to parallel read the namelist in FMS. Either requires
! full copy of source code to use module, or namelist read should be adapted on a
! model-by-model basis. (Urata 11/2017)
use fms_mod, only: error_mesg, FATAL,             &
    open_namelist_file, check_nml_error, close_file, &
    mpp_pe, mpp_root_pe,stdlog
use fms2_io_mod,            only:  file_exists
use field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only:  get_number_tracers
use initracer_mod, only: ntrace_mom, ndust_mass, dust_mass_indx, &
    dust_nb_indx, dust_cor_indx,                                     &
    ice_mass_indx, ice_nb_indx, vapor_indx, nt_nontag
use aerosol_util_mod, only: do_15band, Reff_fixed
#ifdef RELEASE
use null_physics_mod, only: nltecool
#endif

implicit none
private
public :: ames_rt_driver, ames_radsetup, settozero4, &
       settozero, dtridgl, fillpt, qextv, l_nrefv, &
       nltecool,interp1,interp2,interp3, l_nspecti, &
       radactive_dust_TOG, nbin, scale_by_scon

real*8, parameter :: scalep = 100.D+0
real*8, parameter :: cmk = 3.51D+22
integer, parameter :: nbin  = 4
real*8 vrat
real*8 aerad(nbin),rb(nbin+1),dr(nbin)
real*8 vol(nbin)


real*8, parameter ::     dev_dt = 0.63676    !  standard deviation of the dust distribution
! gives an effective variance of 0.5  for dust

real*8, parameter ::     dev_ice= 0.3087     !  standard deviation of the water ice distribution
! gives an effective variance of 0.1  for water ice

real*8, parameter ::     dev_co2ice= 0.3087     ! gives an effective variance of 0.1 for co2 ice, same as water ice

real*8, parameter ::     athird = (1. / 3.)

real*8, parameter ::     dpden_ice = 917.     !  water ice particle density
real*8, parameter ::     dpden_dt = 2.5e+3    !  dust particle density

real*8, parameter ::     dpden_co2ice = 1620.     !  co2 ice particle density

integer :: ntrace

integer :: l_nlayrad
integer :: l_nlevrad
integer :: l_levels, l_layers
integer ::    ima_dt,inb_dt, &
 ima_cld,inb_cld,ima_cor,ima_vap

integer :: l_nspecti, l_nrefi, l_npref, l_ntref, l_ngauss

integer, parameter :: L_NSPECTV =  7
integer, parameter :: L_TAUMAX  = 35
real*8, parameter  :: MAXEXP    = 35.0D0
integer, parameter :: L_PINT    = 51
integer, parameter :: L_REFH2O  = 10
integer, parameter :: L_NREFV   = 6

real*8, parameter :: UBARI = 0.5D0

integer, parameter :: nbin_rt = 20
integer, parameter :: nbin_rtblk = 29
integer, parameter :: nbin_rtco2 = 29

integer, parameter :: nratio  = 15
integer, parameter :: nratioblk = 1
integer, parameter :: nlonv   = L_NSPECTV
integer :: nloni

integer, parameter :: L_NSPECTD =  3
integer, parameter :: nlond   = L_NSPECTD


real*8, allocatable :: GWEIGHT(:)     !     These are for the Gauss-split 0.95 case

!     These are for the CO2+H2O k-coefficients

real*8, parameter :: WREFCO2(L_REFH2O) = [                       &
             9.999999D-1, 9.99999D-1, 9.9999D-1, 9.999D-1,     &
             9.99D-1, 9.9D-1, 9.0D-1, 8.0D-1, 7.0D-1, 6.0D-1 ]

real*8, parameter :: WREFH2O(L_REFH2O) = [                       &
             1.0D-7, 1.0D-6, 1.0D-5, 1.0D-4, 1.0D-3, 1.0D-2,   &
             1.0D-1, 2.0D-1, 3.0D-1, 4.0D-1                  ]

!     These are for CO2-CO2 and CO2-H2 CIA
  real*8, allocatable :: kgbar_tab(:,:)
  real*8, allocatable :: kbbar_tab(:,:)
  real*8, allocatable :: khbar_tab(:,:)

!     If the CO2 optical depth (top to the surface) is less than
!     this value, we place that Gauss-point into the "zeros"
!     channel.  NRC parameter.

real*8, parameter :: TLIMITS = 1.0E-3     !  TLIMITS - TLIMIT for solar part of the spectrum
real*8, parameter :: TLIMITI = 5.0D-3     !  TLIMITI - TLIMIT for the IR
!John's ver. has TLIMITI=1.e-3
!  real*8, parameter :: TLIMITI = 1.0E-3

real*8 PFGASREF(L_PINT)

real*8, allocatable :: WNOI(:), DWNI(:), WAVEI(:)

REAL*8 WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV)
REAL*8 SOLARF(L_NSPECTV), TAURAY(L_NSPECTV)

real*8, allocatable :: CO2I(:,:,:,:,:)
real*8, allocatable :: CO2V(:,:,:,:,:)

real*8, allocatable :: FZEROI(:), PGASREF(:), TGASREF(:), planckir(:,:)
real*8 FZEROV(L_NSPECTV)


!rjw    Arrays qexti, qscati, wi, gi, ggi   are based on an assumed
!      size distribution
!      These are filled in subroutine setrad
!   Same for qextv, qscatv, wv, gv, ggv

real*8, allocatable :: qexti(:),  qscati(:), wi(:)
real*8, allocatable :: GI(:), ggi(:)
!  real*8, allocatable :: ggi(:)

real*8 qextv(L_NSPECTV), qscatv(L_NSPECTV), wv(L_NSPECTV)
real*8 GV(L_NSPECTV), ggv(L_nspectv)
!  real*8 ggv(L_nspectv)


!rjw       declare new diagnostic optical property arrays
real*8, allocatable :: qextd(:), qscatd(:), gd(:)


real*4, parameter :: cor_ratio(nratio) = [ 0.1,0.2,0.25,0.3,     &
                      0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,       &
                      0.8,0.9,0.99 ]

real*8 rad_rt(nbin_rt),radb_rt(nbin_rt+1)
real*8 rad_rtco2(nbin_rtco2),radb_rtco2(nbin_rtco2+1)
real*8 rad_rtblk(nbin_rtblk),radb_rtblk(nbin_rtblk+1)


!    Optical property array ingested from tabulated scattering calculations

real*4, dimension(nratio,nbin_rt,nlonv)   ::  qextv_cld, qscatv_cld, gv_cld
real*4, dimension(nbin_rt,nratio)   ::  fcorrect
real*4, dimension(nratioblk,nbin_rtco2,nlonv)   ::  qextv_co2cld, qscatv_co2cld, gv_co2cld
real*4, dimension(nratioblk,nbin_rtblk,nlonv)   ::  qextv_bcld, qscatv_bcld, gv_bcld
real*4, dimension(nratioblk,nbin_rtblk,nlonv)   ::  qextv_blcld, qscatv_blcld, gv_blcld

real*4, allocatable :: qexti_cld(:,:,:), qscati_cld(:,:,:), gi_cld(:,:,:)
real*4, allocatable :: qexti_dst(:,:),   qscati_dst(:,:),   gi_dst(:,:)

real*4, allocatable :: qexti_bcld(:,:,:), qscati_bcld(:,:,:), gi_bcld(:,:,:)
real*4, allocatable :: qexti_blcld(:,:,:), qscati_blcld(:,:,:), gi_blcld(:,:,:)
real*4, allocatable :: qexti_co2cld(:,:,:), qscati_co2cld(:,:,:), gi_co2cld(:,:,:)

real*4, dimension(nbin_rt,nlonv)  ::   qextv_dst, qscatv_dst, gv_dst


! rjw      read in additional optical arrays for individual brightness temperatures
!      rjw:   should rename this to indicate its for diagnostic channels

  integer :: nwave

!rjw    Add these diagnostic optical properties, to be input from file

real*4, allocatable :: qextd_cld(:,:,:), qscatd_cld(:,:,:), gd_cld(:,:,:)
real*4, allocatable :: qextd_dst(:,:),   qscatd_dst(:,:),   gd_dst(:,:)

!rjw     The following real*4 variables are read in, and copied to the real*8 arrays above
!rjw     The following real*4 variables are read in, and copied to the real*4 arrays above
!rjw    The routine calc_opts currently requires the optical tables values as real*4

  real*4, allocatable ::  qextTb_cld(:,:,:), qscatTb_cld(:,:,:), gTb_cld(:,:,:)
  real*4, allocatable ::  qextTb_dst(:,:),   qscatTb_dst(:,:),   gTb_dst(:,:)


! from dtcommon

integer :: nlay_dt, nlev_dt, JM_dt, IM_dt
integer :: ndp_dt, not_dt

real*8  :: rgas_dt
real*8  :: area_dt, ssarea_dt, dsd_dt, threshold_dt
real*8  :: gtau_dt

logical :: specified_dt

!  Cloud Radiative properties

real*8, allocatable  :: TAUREFCLD(:), TAUREFCLD_UV(:)
real*8, allocatable  :: QEXTREFCLD(:)

real*8, allocatable  :: QXVCLD(:,:), QSVCLD(:,:), GVCLD(:,:)
real*8, allocatable  :: QXICLD(:,:), QSICLD(:,:), GICLD(:,:)

! CO2 cloud radiative properties
real*8, allocatable  :: TAUREFCO2CLD(:)
real*8, allocatable  :: QEXTREFCO2CLD(:)

real*8, allocatable  :: QXVCO2CLD(:,:), QSVCO2CLD(:,:), GVCO2CLD(:,:)
real*8, allocatable  :: QXICO2CLD(:,:), QSICO2CLD(:,:), GICO2CLD(:,:)

!  Dust Radiative properties

real*8, allocatable  :: QEXTREFDST(:)

real*8, allocatable  :: QXVDST(:,:), QSVDST(:,:),  GVDST(:,:)
real*8, allocatable  :: QXIDST(:,:), QSIDST(:,:),  GIDST(:,:)


!     diagnostic radiative properties for brightness TEMPERATURES:
!      note "d" instead of "v" or "i"
!   Note that these will be dimensioned with nlayers,  not levels

real*8, allocatable  :: QXDDST(:,:), QSDDST(:,:), GDDST(:,:)
real*8, allocatable  :: QXDCLD(:,:), QSDCLD(:,:), GDCLD(:,:)



!  NLTE Ref. Atmosphere Table
integer, parameter :: npnlte = 68
real*8 :: pnb(npnlte), ef1(npnlte), ef2(npnlte)
real*8 :: co2vmr(npnlte), o3pvmr(npnlte), n2covmr(npnlte)

logical ::  mcpu0

!  Namelist---------------------------
logical :: tmom = .false.
logical :: radactive_water= .false.
logical :: radactive_cloud = .false.
logical :: radactive_cloud_bin = .false.
logical :: radactive_cloud_bulk = .false. ! Include radiative effects of clouds from bulk scheme
logical :: radactive_co2cloud = .false.   ! Radiatively active CO2 clouds
real    :: scale_cloud_bin= 1.0
real    :: simple_cloud_radius = 1.5                !  cloud radius in microns    
real    :: bulk_ccn = 1.e5   ! number of seed nuclei for bulk h2o cloud RT ( # / kg gaseous CO2)
real    :: bulk_co2_ccn = 1.e5   ! number of seed nuclei for co2 cloud RT ( # / kg gaseous CO2)
logical :: do_nlte = .false.             ! do non-LTE correction in IR
logical :: ames_15band = .false.         ! Use 15 band version. If false, use 12 band
logical :: use_extended_cor_ks = .false. ! Use RT tables with extended temperture range and CO2 line widths appropriate for higher pressures
logical :: do_cia = .false. ! Account for collision induced absorption (CIA) opacity, needed for massive CO2 atmospheres
real*8  :: cia_co2 = 0.0     ! Atmospheric molar concentration of CO2 for collision induced absorption (CIA) calculation, do_cia must be true
real*8  :: cia_h2 = 0.0     ! Atmospheric molar concentration of H2 for collision induced absorption (CIA) calculation, do cia must be true
integer :: radactive_dust_TOG = 0   ! 0 for clear
                           ! 1 for moment tracers
                           ! 2 for bin tracers
                           ! 3 for prescribed, fixed
      
character (len=128) :: rtdata_path = 'INPUT/'
logical :: use_boxinterp12 = .true. ! use boxinterp for k-coef interpolation (ames_15band flag must be false)

logical :: do_fv3_dust_opt= .false.
logical :: do_dust_ir_scale = .false.
logical :: do_irflux_scale = .false. ! optional scaling factor to adjust upward ir flux equal to sigmaT4
logical :: scale_by_scon = .false. ! flag to scale the solar flux by solar constant namelist variable in astronomy
real :: fv3_gfac= 0.65
real :: fv3_sscat = 0.90
real :: fv3_ir_scale= 1.0
logical :: user_fixed_dust_opts = .false.
real, dimension(7) :: qxv_read = (/1.834D0, 2.296D0, 2.672D0,        &
                    2.829D0, 2.698D0, 2.452D0, 2.261D0/)
real, dimension(7) :: qsv_read = (/1.695D0, 2.031D0, 2.583D0,        &
                    2.744D0, 2.626D0, 2.225D0, 1.525D0/)
real, dimension(7) :: gv_read = (/0.551D0, 0.640D0, 0.661D0,        &
                    0.678D0, 0.690D0, 0.743D0, 0.868D0/)
real, dimension(8) :: qxi_read = (/0.008D0, 0.262D0, 0.491D0,        &
                                          1.017D0, 0.444D0, 0., 0., 0. /)
real, dimension(8) :: qsi_read = (/0.001D0, 0.037D0, 0.122D0,        &
                                          0.351D0, 0.336D0, 0., 0., 0. /)
real, dimension(8) :: gi_read = (/0.004D0, 0.030D0, 0.095D0,        &
                                          0.214D0, 0.316D0, 0., 0., 0. /)

namelist /ames_rtmod_nml/ tmom, radactive_water,                &
       radactive_cloud, radactive_cloud_bulk,                   &
       radactive_co2cloud,                                      &
       do_nlte, rtdata_path,                                    &
       ames_15band, use_extended_cor_ks,                        &
       do_cia, cia_co2, cia_h2,                                 &
       radactive_dust_TOG, use_boxinterp12,                     &
       do_fv3_dust_opt, fv3_gfac, fv3_sscat, do_dust_ir_scale,  &
       fv3_ir_scale, radactive_cloud_bin, scale_cloud_bin,      &
       simple_cloud_radius, bulk_ccn, bulk_co2_ccn,             &
       do_irflux_scale, scale_by_scon,                          &
       qxv_read, qsv_read, gv_read, qxi_read, qsi_read, gi_read, &
       user_fixed_dust_opts

!==================================


contains



!=====================================================================
!=====================================================================

subroutine ames_radsetup(nlay,imd,ind,imc,inc,imcor,imvap,mcpu0ARG,nml_fileARG)
!  initialize ames RT module
implicit none
integer, intent(in) :: nlay
integer, intent(in), optional :: imd,ind,imc,inc,imcor,imvap
logical, intent(in), optional :: mcpu0ARG
character (len=128), intent(in), optional :: nml_fileARG
logical :: mcpu0,mppio
real*8 :: PTOP
integer :: unit,io,ierr,i
character (len=128) :: nml_file
integer :: ntprog, ntfam, ntdiag
!===============================================================

! a) First get the total number of tracer : ntrace
call get_number_tracers (MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfam )

l_layers = nlay
l_levels = 2*l_layers+3
l_nlayrad= l_layers+1
l_nlevrad= l_layers+2

if(present(imd)) then
    ima_dt = imd
else
    ima_dt = dust_mass_indx(1)
endif
if(present(ind)) then
    inb_dt = ind
else
    inb_dt = dust_nb_indx(1)
endif
if(present(imc)) then
    ima_cld = imc
else
    ima_cld = ice_mass_indx(1)
endif
if(present(inc)) then
    inb_cld = inc
else
    inb_cld = ice_nb_indx(1)
endif
if(present(imcor)) then
    ima_cor = imcor
else
    ima_cor = dust_cor_indx(1)
endif
if(present(imvap)) then
    ima_vap = imvap
else
    ima_vap = vapor_indx(1)
endif
if(present(mcpu0ARG)) then
    mcpu0 = mcpu0ARG
    mppio = .true.
else
    mcpu0 = .true.
    mppio = .false.
endif
if(present(nml_fileARG)) then
    nml_file = trim(nml_fileARG)
else
    nml_file = 'mars'
endif

! read namelist
if (mppio) then
    ! This is the FV3 method for reading the namelist in parallel
    if (file_exists(trim(nml_file))) then
        unit = open_namelist_file ( )
        ierr=1
        do while (ierr /= 0)
            read  (unit, nml=ames_rtmod_nml, iostat=io, end=20)
            ierr = check_nml_error (io, 'ames_rtmod_nml')
        enddo
    20       call close_file (unit)
    endif
else
    unit = 5
    read  (unit, nml=ames_rtmod_nml)
endif

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=ames_rtmod_nml)

if (mcpu0) then
    print*,'ames_rtmod_nml is read: rtdata_path = ',trim(rtdata_path)
endif

if (ames_15band .and. use_boxinterp12) call error_mesg ('ames_radsetup','ames_15band is not compatible with boxinterp12', FATAL)
if (ames_15band .neqv. do_15band) call error_mesg ('ames_radsetup','ames_15band flag does not match the do_15band setting in aerosol_util_mod', FATAL)
if (radactive_co2cloud .and. (.not. ames_15band)) call error_mesg ('ames_radsetup','radactive_co2cloud must be run with ames_15band', FATAL)
if (radactive_cloud_bulk .and. (.not. ames_15band)) call error_mesg ('ames_radsetup','radactive_cloud_bulk must be run with ames_15band', FATAL)

call set_bands(mcpu0)
call setspv(mcpu0)
call setspi()

!    The following call Interpolate CO2 k coefficients to the finer pressure grid
!           and sets up default optical constants for dust (based on assumed size distribution)
call setrad(mcpu0)

!     The following call initializes optical arrays for dust and water ice clouds
call initcld(mcpu0)

!     The following call allocates storage for the optical property arrays
call init_radtrans()


PTOP = 10.0**PFGASREF(1)
! -----
!  Reading in NLTE Ref. Atmosphere and coefficients (ASB)
!  Pressure is log(nb) in table; converted here to Pa

!     print*,'****OPEN NLTE FILE****'
#ifndef RELEASE
open(96, file=trim(rtdata_path)//'input_file_nlte',status='old')
read (96,*)

pnb = 0.0
ef1 = 0.0
ef2 = 0.0
co2vmr = 0.0
o3pvmr = 0.0
n2covmr = 0.0

!     print*,'****READ NLTE FILE****'
do i=1,npnlte
    read(96,*) pnb(i), ef1(i), ef2(i), co2vmr(i), o3pvmr(i), n2covmr(i)

    pnb(i)=1.0e-4*exp(pnb(i))
end do
close (96)
#endif
! -----

end subroutine ames_radsetup

!=====================================================================
!=====================================================================


subroutine set_bands(mcpu0)
!  set up 12 band or 15 bands
implicit none
integer :: err
logical :: mcpu0

if (use_extended_cor_ks) then !do 15 band radiation & 16 guass points
    L_NSPECTI =  8
    L_NREFI   = 7
    L_NPREF   = 25
    L_NTREF   = 16
    L_NGAUSS  = 17
    allocate(CO2I(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS), stat=err)
    allocate(CO2V(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS), stat=err)
    allocate(GWEIGHT(L_NGAUSS), stat=err)
    allocate(kgbar_tab(L_NSPECTI,15), stat=err)
    allocate(kbbar_tab(L_NSPECTI,15), stat=err)
    allocate(khbar_tab(L_NSPECTI,15), stat=err)
!  For 16 gauss points
    GWEIGHT(:) = [                       &
             4.8083554740D-02, 1.0563099137D-01,              &
             1.4901065679D-01, 1.7227479710D-01,              &
             1.7227479710D-01, 1.4901065679D-01,              &
             1.0563099137D-01, 4.8083554740D-02,              &
             2.5307134073D-03, 5.5595258613D-03,              &
             7.8426661469D-03, 9.0670945845D-03,              &
             9.0670945845D-03, 7.8426661469D-03,              &
             5.5595258613D-03, 2.5307134073D-03,  0.0D0 ]


else if (ames_15band) then !do 15 band radiation & 32 guass points
    L_NSPECTI =  8
    L_NREFI   = 7
    L_NPREF   = 25
    L_NTREF   = 13
    L_NGAUSS  = 33
    allocate(CO2I(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS), stat=err)
    allocate(CO2V(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS), stat=err)
    allocate(GWEIGHT(L_NGAUSS), stat=err)
    ! For 32 gauss points
    GWEIGHT(:) = [                       &
         1.2897418221D-02, 2.9570423871D-02, 4.5200293049D-02, &
         5.9198761346D-02, 7.1058094688D-02, 8.0349346713D-02, &
         8.6736622146D-02, 8.9989039966D-02, 8.9989039966D-02, &
         8.6736622146D-02, 8.0349346713D-02, 7.1058094688D-02, &
         5.9198761346D-02, 4.5200293049D-02, 2.9570423871D-02, &
         1.2897418221D-02, 6.7881148529D-04, 1.5563380985D-03, &
         2.3789627921D-03, 3.1157242814D-03, 3.7398997204D-03, &
         4.2289129849D-03, 4.5650853761D-03, 4.7362652614D-03, &
         4.7362652614D-03, 4.5650853761D-03, 4.2289129849D-03, &
         3.7398997204D-03, 3.1157242814D-03, 2.3789627921D-03, &
         1.5563380985D-03, 6.7881148529D-04, 0.0D0  ]
else    !do 12 band radiation & 16 gauss points
    L_NSPECTI =  5
    L_NREFI   = 4
    L_NPREF   = 11
    L_NTREF   =  7
    L_NGAUSS  = 17
    allocate(CO2I(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS), stat=err)
    allocate(CO2V(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS), stat=err)
    allocate(GWEIGHT(L_NGAUSS), stat=err)
    !  For 16 gauss points
    GWEIGHT(:) = [                       &
          4.8083554740D-02, 1.0563099137D-01,              &
          1.4901065679D-01, 1.7227479710D-01,              &
          1.7227479710D-01, 1.4901065679D-01,              &
          1.0563099137D-01, 4.8083554740D-02,              &
          2.5307134073D-03, 5.5595258613D-03,              &
          7.8426661469D-03, 9.0670945845D-03,              &
          9.0670945845D-03, 7.8426661469D-03,              &
          5.5595258613D-03, 2.5307134073D-03,  0.0D0 ]
endif
nloni   = L_NSPECTI
allocate(WNOI(L_NSPECTI), stat=err)
allocate(DWNI(L_NSPECTI), stat=err)
allocate(WAVEI(L_NSPECTI), stat=err)
allocate(FZEROI(L_NSPECTI), stat=err)
allocate(PGASREF(L_NPREF), stat=err)
allocate(TGASREF(L_NTREF), stat=err)
allocate(qexti(L_NSPECTI), stat=err)
allocate(qscati(L_NSPECTI), stat=err)
allocate(wi(L_NSPECTI), stat=err)
allocate(planckir(L_NSPECTI,8501), stat=err)
allocate(qexti_cld(nratio,nbin_rt,nloni), stat=err)
allocate(qscati_cld(nratio,nbin_rt,nloni), stat=err)
allocate(gi_cld(nratio,nbin_rt,nloni), stat=err)
allocate(qexti_bcld(nratioblk,nbin_rtblk,nloni), stat=err)
allocate(qscati_bcld(nratioblk,nbin_rtblk,nloni), stat=err)
allocate(gi_bcld(nratioblk,nbin_rtblk,nloni), stat=err)
allocate(qexti_blcld(nratioblk,nbin_rtblk,nloni), stat=err)
allocate(qscati_blcld(nratioblk,nbin_rtblk,nloni), stat=err)
allocate(gi_blcld(nratioblk,nbin_rtblk,nloni), stat=err)
allocate(qexti_co2cld(nratioblk,nbin_rtco2,nloni), stat=err)
allocate(qscati_co2cld(nratioblk,nbin_rtco2,nloni), stat=err)
allocate(gi_co2cld(nratioblk,nbin_rtco2,nloni), stat=err)
allocate(qexti_dst(nbin_rt,nloni), stat=err)
allocate(qscati_dst(nbin_rt,nloni), stat=err)
allocate(gi_dst(nbin_rt,nloni), stat=err)
allocate(GI(L_NSPECTI), stat=err)
allocate(ggi(L_NSPECTI), stat=err)


!   rjw   #ifdef SKIP

!     New variables for optical properties at selected diagnostic wavelengths

nwave= 5    !    nlond should be <= nwave

allocate( qextd(nlond), stat=err)
allocate(qscatd(nlond), stat=err)
allocate(    gd(nlond), stat=err)

allocate( qextd_dst(nbin_rt,nlond), stat=err)
allocate(qscatd_dst(nbin_rt,nlond), stat=err)
allocate(    gd_dst(nbin_rt,nlond), stat=err)

allocate( qextd_cld(nratio,nbin_rt,nlond), stat=err)
allocate(qscatd_cld(nratio,nbin_rt,nlond), stat=err)
allocate(    gd_cld(nratio,nbin_rt,nlond), stat=err)

!   These are for input only (real*4)  nwave = 5;  nlond = 3

  allocate(  qextTb_cld(nratio,nbin_rt,nwave), stat=err)
  allocate( qscatTb_cld(nratio,nbin_rt,nwave), stat=err)
  allocate(     gTb_cld(nratio,nbin_rt,nwave), stat=err)

  allocate(  qextTb_dst(nbin_rt,nwave), stat=err)
  allocate( qscatTb_dst(nbin_rt,nwave), stat=err)
  allocate(     gTb_dst(nbin_rt,nwave), stat=err)

!   rjw   #endif SKIP

if (mcpu0) then
    print*,'15 band:',ames_15band,'gweight after initialization:',gweight
    print*,'nspecti = ',l_nspecti
endif

end subroutine set_bands

!=====================================================================
subroutine init_radtrans()
!=====================================================================

!   allocate and load cloud and dust optical properties
implicit none
integer :: err

allocate(taurefcld(l_levels+1), stat=err)
allocate(taurefcld_uv(l_levels+1), stat=err)
allocate(qextrefcld(l_levels+1), stat=err)

allocate(qxvcld(l_levels+1,l_nspectv), stat=err)
allocate(qsvcld(l_levels+1,l_nspectv), stat=err)
allocate(gvcld(l_levels+1,l_nspectv), stat=err)

allocate(qxicld(l_levels+1,l_nspecti), stat=err)
allocate(qsicld(l_levels+1,l_nspecti), stat=err)
allocate(gicld(l_levels+1,l_nspecti), stat=err)


! CO2 cloud properties

allocate(taurefco2cld(l_levels+1), stat=err)
allocate(qextrefco2cld(l_levels+1), stat=err)

allocate(qxvco2cld(l_levels+1,l_nspectv), stat=err)
allocate(qsvco2cld(l_levels+1,l_nspectv), stat=err)
allocate(gvco2cld(l_levels+1,l_nspectv), stat=err)

allocate(qxico2cld(l_levels+1,l_nspecti), stat=err)
allocate(qsico2cld(l_levels+1,l_nspecti), stat=err)
allocate(gico2cld(l_levels+1,l_nspecti), stat=err)

!  dust radiative properties

allocate(qextrefdst(l_levels+1), stat=err)

allocate(qxvdst(l_levels+1,l_nspectv), stat=err)
allocate(qsvdst(l_levels+1,l_nspectv), stat=err)
allocate(gvdst(l_levels+1,l_nspectv), stat=err)

allocate(qxidst(l_levels+1,l_nspecti), stat=err)
allocate(qsidst(l_levels+1,l_nspecti), stat=err)
allocate(gidst(l_levels+1,l_nspecti), stat=err)

!    rjw     zeros out the cloud optical arrays
call ini_optcld(qxvcld,qxicld,qsvcld,qsicld,gvcld,gicld, &
             qextrefcld,taurefcld)

!  zero out CO2 cloud optical arrays
call ini_optcld(qxvco2cld,qxico2cld,qsvco2cld,qsico2cld, &
             gvco2cld,gico2cld, &
             qextrefco2cld,taurefco2cld)

!    rjw   sets default values for dust optical arrays
call ini_optdst(gv,gi,                                   &
       qxvdst,qxidst,qsvdst,qsidst,gvdst,gidst,qextrefdst  )



!   rjw ----optical arrays for diagnostic brightness temperatures----
!        Note that these are dimensioned by layers; not levels
allocate(qxddst(l_layers,l_nspectd), stat=err)
allocate(qsddst(l_layers,l_nspectd), stat=err)
allocate( gddst(l_layers,l_nspectd), stat=err)

allocate(qxdcld(l_layers,l_nspectd), stat=err)
allocate(qsdcld(l_layers,l_nspectd), stat=err)
allocate( gdcld(l_layers,l_nspectd), stat=err)

!  rjw   Initialize the diagnostic arrays  (set vals to 0)

call ini_opt_diag( qxddst, qsddst, gddst )


end subroutine init_radtrans



!=====================================================================
subroutine ini_opt_diag(Qxd,Qsd,gdd)
!=====================================================================


integer dimv

real*8  Qxd(l_layers,L_NSPECTd)
real*8  Qsd(l_layers,L_NSPECTd)
real*8  gdd(l_layers,L_NSPECTd)

dimv= l_nspectd

call settozero(dimv*l_layers,qxd)
call settozero(dimv*l_layers,qsd)
call settozero(dimv*l_layers,gdd)


return

end subroutine ini_opt_diag


#ifdef SKIP
        qextd_dst(:,ii)= qextTb_dst(:,ii)
       qscatd_dst(:,ii)= qscatTb_dst(:,ii)
           gd_dst(:,ii)= gTb_dst(:,ii)

        qextd_cld(:,:,ii)=   qextTb_cld(:,:,ii)
       qscatd_cld(:,:,ii)=  qscatTb_cld(:,:,ii)
           gd_cld(:,:,ii)=   gTb_cld(:,:,ii)

#endif SKIP

!=====================================================================
!=====================================================================

subroutine end_radtrans()
!  deallocate arrays

implicit none


deallocate(taurefcld)
deallocate(taurefcld_uv)
deallocate(qextrefcld)
deallocate(qxvcld)
deallocate(qsvcld)
deallocate(gvcld)
deallocate(qxicld)
deallocate(qsicld)
deallocate(gicld)
deallocate(qextrefdst)
deallocate(qxvdst)
deallocate(qsvdst)
deallocate(gvdst)
deallocate(qxidst)
deallocate(qsidst)
deallocate(gidst)

deallocate(taurefco2cld)
deallocate(qextrefco2cld)
deallocate(qxvco2cld)
deallocate(qsvco2cld)
deallocate(gvco2cld)
deallocate(qxico2cld)
deallocate(qsico2cld)
deallocate(gico2cld)

deallocate(co2i)
deallocate(co2v)
deallocate(gweight)
deallocate(wnoi)
deallocate(dwni)
deallocate(wavei)
deallocate(fzeroi)
deallocate(pgasref)
deallocate(TGASREF)
deallocate(qexti)
deallocate(qscati)
deallocate(wi)
deallocate(planckir)
deallocate(qexti_cld)
deallocate(qscati_cld)
deallocate(gi_cld)
deallocate(qexti_dst)
deallocate(qscati_dst)
deallocate(gi_dst)
deallocate(GI)
deallocate(ggi)


!rjw    deallocate all the arrays associated with brightness temperatures

deallocate(qextd)
deallocate(qscatd)
deallocate(gd)

deallocate(qextd_cld)
deallocate(qscatd_cld)
deallocate(gd_cld)

deallocate(qextd_dst)
deallocate(qscatd_dst)
deallocate(gd_dst)

deallocate(qxddst)
deallocate(qsddst)
deallocate(gddst)

deallocate(qxdcld)
deallocate(qsdcld)
deallocate(gdcld)

deallocate(qextTb_dst)
deallocate(qscatTb_dst)
deallocate(gTb_dst)

deallocate(qextTb_cld)
deallocate(qscatTb_cld)
deallocate(gTb_cld)


end subroutine end_radtrans

!=====================================================================
!=====================================================================

subroutine ames_rt_driver(pl, tl, ptop,           &
                 tstrat, tsurf, QH2Oin,                 &
                 acosz, rsdist,                         &
		 albv, albi,                            &
                 taud_bin_in,                           &
		 taud_fix_in,                           &
		 dosw, dolw,                            &
		 sw_heating, lw_heating,                &
                 qtrace, cldice_bin,                    &
		 outsolARG, directonlyARG,              &
                 taudust_diagARG,                       &   !taudust (:,:,2)
		 taucloud_diagARG,                      &   !taucloud(:,:,2)
                 tauco2cloud_diagARG,                   &   !taucloud(:,:,2)
                 diffvtARG, directsolARG, detauARG,     &
                 fluxupvARG, fluxdnvARG, fmnetvARG,     &
                 fluxupiARG, fluxdniARG, fmnetiARG,     &
                 nfluxtopvARG, nfluxtopiARG,            &
		 swtotARG,                              &
                 lw_heating_spec,                       &
		 htrt2,                                 &
		 diag,                                  &   ! scalar
		 taurefd_out,                           &   ! dustref(:,:,kd)
                 taurefc_out,                           &   ! cldref (:,:,kd)
                 taurefco2c_out,                        &   ! cldco2ref (:,:,kd)
	         taudust_momARG,taudust_fixARG,         &   ! taudust_mom(:,:,2)
		 tstrat_dt, nco2,nice_blk,              &
                 taudust_reffARG,                       &
    tbands  )



implicit none


real*8, intent(in) :: pl(l_levels)
real*8, intent(in) :: tl(l_levels)
real*8, intent(in) :: ptop
real*8, intent(in) :: tsurf
real*8, intent(in) :: tstrat
real*8, intent(in) :: qh2oin(l_layers)
real*8, intent(in) :: acosz
real*8, intent(in) :: rsdist
real*8, intent(in) :: albv
real*8, intent(in) :: albi
real*8, intent(in) :: taud_bin_in(l_layers)
real*8, intent(in) :: taud_fix_in(l_layers)
real*8, intent(in) :: cldice_bin (l_layers)
logical, intent(in) :: dosw
logical, intent(in) :: dolw
real*8, intent(inout) :: sw_heating(l_layers)
real*8, intent(inout) :: lw_heating(l_layers)
real*8, intent(in) :: qtrace(l_layers,ntrace)
logical, intent(in), optional :: directonlyarg
logical, intent(in) :: diag
integer, intent(in) :: nco2
integer, intent(in) :: nice_blk
real*8, intent(out), optional :: taudust_diagarg(2)
real*8, intent(out), optional :: taudust_momarg(2)
real*8, intent(out), optional :: taudust_fixarg(2)
real*8, intent(out), optional :: taucloud_diagarg(4)
real*8, intent(out), optional :: tauco2cloud_diagarg(2)
real*8, intent(out), optional :: outsolarg
real*8, intent(out), optional :: diffvtarg
real*8, intent(out), optional :: directsolarg
real*8, intent(inout), optional :: detauarg(l_nspectv,l_ngauss)
real*8, intent(out), optional :: fluxupvarg(l_nlayrad)
real*8, intent(out), optional :: fluxdnvarg(l_nlayrad)
real*8, intent(out), optional :: fmnetvarg(l_nlayrad)
real*8, intent(out), optional :: fmnetiarg(l_nlayrad)
real*8, intent(out), optional :: fluxupiarg(l_nlayrad)
real*8, intent(out), optional :: fluxdniarg(l_nlayrad)
real*8, intent(out), optional :: nfluxtopvarg
real*8, intent(out), optional :: nfluxtopiarg
real*8, intent(out), optional :: swtotarg
real*8, intent(out), optional :: lw_heating_spec(l_nspecti,l_layers)
real*8, intent(out), optional :: htrt2(l_layers) ![k/s]
real*8, intent(out), optional :: taurefd_out(l_layers)
real*8, intent(out), optional :: taurefc_out(l_layers,2)
real*8, intent(out), optional :: taurefco2c_out(l_layers)
real*8, intent(out) :: tstrat_dt
real*8, intent(out), optional :: taudust_reffarg(ntrace,2)
real*8, intent(out), optional :: tbands(3)

! local versions of optional arguments
logical :: directonly
real*8 :: taudust_diag(2)
real*8 :: taureff_split(ntrace,2)
real*8 :: taucloud_diag(4)
real*8 :: tauco2cloud_diag(2)
real*8 :: outsol
real*8 :: diffvt
real*8 :: directsol
real*8 :: detau(l_nspectv,l_ngauss)
real*8 :: fluxupv(l_nlayrad)
real*8 :: fluxdnv(l_nlayrad)
real*8 :: fmnetv(l_nlayrad)
real*8 :: fluxupi(l_nlayrad)
real*8 :: fluxdni(l_nlayrad)
real*8 :: fluxupis(l_nspecti,l_nlayrad)
real*8 :: fluxdnis(l_nspecti,l_nlayrad)
real*8 :: fmneti(l_nlayrad), fmnetis(l_nspecti,l_nlayrad)
real*8 :: nfluxtopv
real*8 :: nfluxtopi, nfluxtopis(l_nspecti)
real*8 :: swtot

real*8 :: acoszt

real*8, parameter :: pnlteref = 0.001 !mbar

real*8  :: lw_heatingin(l_layers)
real*8  :: alpha(l_layers), q4(l_layers)
real*8  :: heatingir15(l_layers)
real*8  :: plevnlte(l_levels)
real*8  :: plevnlte2(l_levels)
real*8  :: tl2(l_levels)
real*8  :: htrt(l_levels) ![k/s]
real*8  :: htrt3(l_levels) ![k/s]
real*8  :: plev(l_levels), pmid(l_levels)
real*8  :: plev2(l_levels)
real*8  :: delp(l_layers)
real*8  :: tlev(l_levels), tmid(l_levels)
real*8  :: qh2o(l_levels)
real*8  :: tauref(l_levels+1)
real*8  :: tauref_mom(l_levels+1)
real*8  :: tauref_mom_temp(l_levels+1)
real*8  :: taureftemp(l_levels+1)
real*8  :: tauref_bin(l_levels+1)
real*8  :: tauref_fix(l_levels+1)
real*8  :: taucum(l_levels)
real*8  :: cldice_bin_scale(L_LAYERS)

real*8  :: taugsurf(l_nspectv,l_ngauss-1)
real*8  :: dtauv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: tauv(l_nlevrad,l_nspectv,l_ngauss)
real*8  :: taucumv(l_levels,l_nspectv,l_ngauss)
real*8  :: dtaui(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: taucumi(l_levels,l_nspecti,l_ngauss)
real*8  :: sol(l_nspectv)

real*8  :: cosbv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: wbarv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: wbari(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: cosbi(l_nlayrad,l_nspecti,l_ngauss)

real*8  :: taugsurfi(l_nspecti,l_ngauss-1)
real*8  :: cumtauv(l_nspectv,l_ngauss)

real*8  :: gcp, deltap
integer :: nw, m, l, k, n, nn, ng, ndx, nt
real*8  :: tausurf

real*8,  dimension(l_levels+1) :: tauref_save, taurefcld_save

real*8 :: swfac
real*8 :: swtot0, irtot0

logical :: diag2


!    Some new arrays

integer, parameter  ::  nbrite = 3

real*8 :: taurefcol(l_layers), taurefcldcol(l_layers)

real*8 :: tau1d(l_layers),  sscat1d(l_layers), gfac1d(l_layers), tauxx(l_layers)
real*8 :: stemp1(l_layers),  stemp2(l_layers)


!rjw      New array for tracking the expanded particle size distribution
real*8 ::  surf2d(nbin_rt,l_layers)
real*8 :: surfcld(nbin_rt,l_layers)
real*8 :: surf1d(nbin_rt)

real*8 :: tcol(l_layers+1)

real*8  :: sfc_emiss, xwl

real*8 ::  xcenter(nbrite)

integer ::  nxi, nxiarr(3)

!==================================================


!    added brightness temperature diagnostics
xcenter(1)= 7.8
xcenter(2)= 25.0
xcenter(3)= 32.0



diag2 = acosz.ge.0.98

gcp= grav/cp_air

if (present(directonlyARG)) then
    directonly = directonlyARG
else
    directonly = .false.
endif
if (present(detauARG)) then
    detau = detauARG
else
    detau = 0.d0
endif

do NW=1,L_NSPECTV
    SOL(nw) = SOLARF(NW)/RSDIST
end do

swtot = sum(solarf)

call fillpt(pl, ptop, tsurf, tstrat, tl, &
         plev, tlev, pmid, tmid, delp)

if (directonly) then
    if (acosz .ge. 1.e-4) then
        call dsolflux(SOL,ACOSZ,DETAU,DIRECTSOL)
    else
        directsol = 0.0
    endif
    if (present(directsolARG)) directsolARG=directsol
    return
endif


!  Fill QPI with water information
qh2o = 0.d0
if(radactive_water .eqv. .true.) then
    m = ima_vap
    do  l = 1, l_layers
        k = 2*l+2
        qh2o(k)   = qh2oin(l)
        qh2o(k+1) = qh2o(k)
    end do
else
    qh2o(4:)   = 1.0e-7
end if

!  Fill the TAUREF array, the dust column density for each GCM sub-layer
!  Calculates
tauref = 0.0
taureff_split = 0.0
taureftemp = 0.0
tauref_mom = 0.0
tauref_mom_temp = 0.0
tauref_bin = 0.0
tauref_fix = 0.0
taucum = 0.0
taudust_diag = 0.0

!rjw     add new argument :   surf2d

call opt_dst(qtrace,pl, &
           qxvdst,qxidst,qsvdst,qsidst,gvdst,gidst,    &
           qextrefdst,taureftemp,  surf2d,             &
           taudust_diag,taureff_split)

if (present(taudust_reffarg)) taudust_reffarg(:,:) = taureff_split(:,:)
if (present(taudust_momarg)) taudust_momarg = taudust_diag
if (present(taudust_fixarg)) then
    call ini_optdst(GV,GI,Qxvdst,Qxidst,Qsvdst,Qsidst,gvdst,gidst,Qextrefdst, surf2d, .true.)
    taudust_fixarg(1) = sum(taud_fix_in)
    taudust_fixarg(2) = taudust_fixarg(1)/Qextrefdst(1) * ( Qxidst(1,l_nrefi)-Qsidst(1,l_nrefi) )
endif

if(.not.(nml_switch(radactive_dust_tog,0))) then
    if(radactive_dust_tog.eq.1) then
        tauref_mom=taureftemp
        tauref = tauref_mom
    endif
    if(radactive_dust_tog.eq.2) then
        call ini_optdst(GV,GI,Qxvdst,Qxidst,Qsvdst,Qsidst,gvdst,gidst,Qextrefdst, surf2d, .true.)
        do k= 1, l_layers
            n= 2*k + 2
            tauref_bin(n  )= taud_bin_in(k)*(plev(n)-plev(n-1))/delp(k)
            tauref_bin(n+1)= taud_bin_in(k)*(plev(n+1)-plev(n))/delp(k)
        enddo
        tauref = tauref_bin
    endif
    if(radactive_dust_tog.eq.3) then
        call ini_optdst(GV,GI,Qxvdst,Qxidst,Qsvdst,Qsidst,gvdst,gidst,Qextrefdst, surf2d, .true.)
        do k= 1, l_layers
            n= 2*k + 2
            tauref_fix(n  )= taud_fix_in(k)*(plev(n)-plev(n-1))/delp(k)
            tauref_fix(n+1)= taud_fix_in(k)*(plev(n+1)-plev(n))/delp(k)
        enddo
        tauref = tauref_fix
    endif

!    tauref = tauref_mom + tauref_bin + tauref_fix

    do k = 1,3
        tauref(k) = 0.0
        taucum(k) = 0.0
    end do

    do l=1,l_layers
        do nn=1,2
            k=2*l+1+nn
            taucum(k) = taucum(k-1) + tauref(k)
        end do
        n= 2*l + 2
        if (present(taurefd_out)) taurefd_out(l) = tauref(n)/((plev(n)-plev(n-1))/delp(l))
    end do
    taudust_diag(1)=taucum(l_levels)
    taudust_diag(2)= taudust_diag(1)*( qxidst(l_levels,l_nrefi) - qsidst(l_levels,l_nrefi) )/qextrefdst(l_levels)

else
    taudust_diag(:)=0.
    tauref(:)=0.
endif  ! end check on dust cases


!  Fill special bottom radiation level to zero.
!
tauref(l_levels+1) = 0.0
tausurf       = taucum(l_levels)

mcpu0 = (mpp_pe() == mpp_root_pe())

if (radactive_cloud) then
    call opt_cld(qtrace,pl, &
        qxvcld,qxicld,qsvcld,qsicld,gvcld,gicld, &
        qextrefcld,taurefcld,  &
        taucloud_diag,taurefcld_uv)
    if (present(taurefc_out)) then
        do l=1,l_layers
            n= 2*l + 2
            taurefc_out(l,1) = taurefcld(n)/((plev(n)-plev(n-1))/delp(l))
            taurefc_out(l,2) = taurefcld_uv(n)/((plev(n)-plev(n-1))/delp(l))
        end do
    end if

else if( radactive_cloud_bin ) then

    cldice_bin_scale(:)= scale_cloud_bin*cldice_bin(:)

    call opt_cld_simple(cldice_bin_scale, pl,               &
        qxvcld,qxicld,qsvcld,qsicld,gvcld,gicld, &
        qextrefcld,taurefcld, surfcld,                 &
        taucloud_diag, taurefcld_uv    )

    if (present(taurefc_out)) then
        do l=1,l_layers
            n= 2*l + 2
            taurefc_out(l,1) = taurefcld(n)/((plev(n)-plev(n-1))/delp(l))
            taurefc_out(l,2) = taurefcld_uv(n)/((plev(n)-plev(n-1))/delp(l))
        end do
    end if

else if( radactive_cloud_bulk ) then

    call opt_cld_blk(qtrace,pl,tlev, &
        qxvcld,qxicld,qsvcld,qsicld,gvcld,gicld, &
        qextrefcld,taurefcld,  &
        taucloud_diag,nice_blk,mcpu0, taurefcld_uv)

    if (present(taurefc_out)) then
        do l=1,l_layers
            n= 2*l + 2
            taurefc_out(l,1) = taurefcld(n)/((plev(n)-plev(n-1))/delp(l))
            taurefc_out(l,2) = taurefcld_uv(n)/((plev(n)-plev(n-1))/delp(l))
        end do
    end if

else
    do k=1,l_levels
        taurefcld(k) = 0.0d0
    end do
    if (present(taurefc_out)) taurefc_out = 0.0d0
endif


!CO2 cloud radiative effects
if( radactive_co2cloud) then

    call opt_co2cld(qtrace,pl,               &
        qxvco2cld,qxico2cld,qsvco2cld,qsico2cld, &
        gvco2cld,gico2cld, &
        qextrefco2cld,taurefco2cld,                 &
        tauco2cloud_diag, nco2)


    if (present(taurefco2c_out)) then
        do l=1,l_layers
            n= 2*l + 2
            taurefco2c_out(l) = taurefco2cld(n)/((plev(n)-plev(n-1))/delp(l))
        end do
    end if
endif

if (dosw) then
    !  set up, and solve for, the solar (visible) fluxes, if the sun
    !  is up
    if(acosz.ge.1.0e-4) then
        if (ames_15band .or. use_boxinterp12) then
            call optcv15(dtauv,tauv,taucumv,plev,     &
                     qxvdst,qsvdst,gvdst,wbarv,cosbv,  &
                      tauref,tmid,pmid,taugsurf,qh2o,   &
                      qextrefcld,taurefcld,qxvcld,qsvcld,gvcld, &
                      qextrefco2cld,taurefco2cld,  &
                      qxvco2cld,qsvco2cld,gvco2cld)
        else
            call optcv(dtauv,tauv,taucumv,plev,     &
                     qxvdst,qsvdst,gvdst,wbarv,cosbv,  &
                      tauref,tmid,pmid,taugsurf,qh2o,   &
                      qextrefcld,taurefcld,qxvcld,qsvcld,gvcld)
        endif

        do nw=1,l_nspectv
            do ng=1,l_ngauss
                cumtauv(nw,ng) = taucumv(l_levels,nw,ng)
            enddo
        enddo

        call sfluxv(dtauv,tauv,taucumv,albv,wbarv,cosbv,      &
                   acosz,sol,nfluxtopv,fmnetv,     &
                   fluxupv,fluxdnv,diffvt,taugsurf, &
                   detau,diag)

        call dsolflux(sol,acosz,detau,directsol)
        !          combined surface/atmosphere albedo at the top of the atmosphere
        outsol =   fluxupv(1)/fluxdnv(1)

    else
    !  If the sun is down, no solar flux, nor downward flux. . .
                DIRECTSOL = 0.0
                NFLUXTOPV = 0.0
                fmnetv(:) = 0.0
                outsol    = 0.0
                fluxupv = 0.0
                fluxdnv = 0.0
                diffvt = 0.0
    end if    ! condition on sun elevation
endif  !end sw check

if (dolw) then !do ir calculation
    if (ames_15band .or. use_boxinterp12) then
        call optci15(dtaui,taucumi,plev,tlev,    &
                     qextrefdst,qxidst,qsidst,gidst,cosbi,wbari,  &
                     tauref,tmid,pmid,taugsurfi,qh2o,     &
                     qextrefcld,taurefcld,qxicld,qsicld,gicld,     &
                     qextrefco2cld,taurefco2cld,    &
                     qxico2cld,qsico2cld,gico2cld)
    else
        call optci(dtaui,taucumi,plev,    &
                     qextrefdst,qxidst,qsidst,gidst,cosbi,wbari,  &
                     tauref,tmid,pmid,taugsurfi,qh2o,     &
                     qextrefcld,taurefcld,qxicld,qsicld,gicld)
    endif

    call sfluxi(plev,tlev,dtaui,taucumi,albi,  &
                  cosbi,wbari,nfluxtopi,fmneti,     &
                  fluxupi,fluxdni,taugsurfi,nfluxtopis, &
                  fmnetis,fluxupis,fluxdnis,diag)

endif !end lw check


!rjw -----------------------------------------
!  diagnostic brightness temperatures

  nxiarr(1)= 5
  nxiarr(2)= 2
  nxiarr(3)= 2

  sfc_emiss= 1.0 - albi

DO k= 1, l_layers
   nn= 2*k + 2
   tcol(k)= tl(nn)
enddo
!  include a bottom layer
tcol(l_layers+1)= tsurf    !    tcol(l_layers)

!      Calculate optical properties at diagnostic wavelengths :   nlond of these
!       Use surf2d(nbin_rt,l_layers)  === expansion of the particle size distribution
!       qextd_dst(nbin_rt,nlond)      -->   qxddst(l_layers,nlond)   etc.

    DO k= 1, l_layers
        nn= 2*k + 2

        surf1d(:)= surf2d(:,k)
        call calc_opts(nlon=nlond,nbin=nbin_rt, &
                        qext_in=qextd_dst,qscat_in=qscatd_dst,g_in=gd_dst, &
                        surf=surf1d,   &
                        qx_out=qxddst(k,:),qs_out=qsddst(k,:),g_out=gddst(k,:)  )

        if( radactive_cloud_bin ) then
           surf1d(:)= surfcld(:,k)
           call calc_opts(nlon=nlond,nbin=nbin_rt, &
                        qext_in=qextd_cld(1,:,:),qscat_in=qscatd_cld(1,:,:), &
                        g_in=gd_cld(1,:,:),   &
                        surf=surf1d,          &
                        qx_out=qxdcld(k,:),qs_out=qsdcld(k,:),g_out=gdcld(k,:)  )
         endif

    enddo

  DO nxi= 1, nlond
    nw= 3               !   ******   3 versions of the 32 micron brightness temp
    xwl= xcenter(nw)

!           extract layer optical thicknesses:

    DO k= 1, l_layers
         nn= 2*k + 2
         taurefcol(k)=    qxddst(k,nw) * ( tauref(nn)    + tauref(nn+1) )    / ( qextrefdst(nn)+1.e-50)
         taurefcldcol(k)= qxdcld(k,nw) * ( taurefcld(nn) + taurefcld(nn+1) ) / ( qextrefcld(nn)+1.e-50)
    enddo

    DO k= 1, l_layers
        nn= 2*k + 2
        tau1d(k)=   taurefcol(k)
        sscat1d(k)= qsddst(k,nw) / ( qxddst(k,nw)+1.e-50)
        gfac1d(k)=  gddst(k,nw)

        stemp1(k)= sscat1d(k)
        stemp2(k)= qsdcld(k,nw) / ( qxdcld(k,nw)+1.e-50)

        !!!!! tauxx(k)= gfac1d(k)

    !!!    sscat1d(k)= qsidst(nn,nxi) / qxidst(nn,nxi)
    !!!     gfac1d(k)=  gidst(nn,nxi)

!   Will eventually need to properly combine the optical properties of dust and clouds

    if( nxi==1 ) then                   ! dust only
        tau1d(k)=   taurefcol(k)
        sscat1d(k)= qsddst(k,nw) / ( qxddst(k,nw)+1.e-50)
        gfac1d(k)=  gddst(k,nw)
     else if( nxi==2 ) then              ! cloud only
        tau1d(k)=   taurefcldcol(k)
        sscat1d(k)= qsdcld(k,nw) / ( qxdcld(k,nw)+1.e-50)
        gfac1d(k)=  gdcld(k,nw)
     else                              ! dust and cloud
        tau1d(k)=   taurefcol(k)  + taurefcldcol(k)
        sscat1d(k)= ( taurefcol(k)*stemp1(k) + taurefcldcol(k)*stemp2(k) )  &
                             / ( tau1d(k) + 1.e-50 )
        gfac1d(k)= ( taurefcol(k)*gddst(k,nw) + taurefcldcol(k)*gdcld(k,nw) )  &
                            / ( tau1d(k) + 1.e-50 )

    endif

  enddo

    call tbright_diag(  l_layers, tau1d, sscat1d, gfac1d, tcol, tsurf, sfc_emiss, xcenter(nw), tbands(nxi) )

enddo

!rjw -------------



!             Heating rates:  based on flux divergence
k=1
deltap= plev(5)-plev(3)
if (dolw) lw_heating(k)= ( fmneti(k+1) - fmneti(k) )/deltap
if (dosw) sw_heating(k)= ( fmnetv(k+1) - fmnetv(k) )/deltap

do nw = 1,L_NSPECTI
    if (dolw) lw_heating_spec(nw,k)=  &
                 ( fmnetis(nw,k+1) - fmnetis(nw,k) )/deltap
end do

if (dolw .or. dosw) then
    do k= 2, l_layers
        if (dolw) lw_heating(k)= ( fmneti(k+1) - fmneti(k) )/delp(k)
        if (dosw) sw_heating(k)= ( fmnetv(k+1) - fmnetv(k) )/delp(k)
    enddo

    ! spectral dependence calculation for NLTE (ASB)
    do nw = 1,L_NSPECTI
        do k= 2, l_layers
            if (dolw) lw_heating_spec(nw,k)=  &
                       ( fmnetis(nw,k+1) - fmnetis(nw,k) )/delp(k)
        end do
    end do
endif

if (dolw) lw_heating= gcp * lw_heating / scalep
if (dolw) lw_heating_spec= gcp * lw_heating_spec / scalep
if (dosw) sw_heating= gcp * sw_heating / scalep

#ifndef RELEASE
if (do_nlte) then
! longwave (IR 15 micron cooling) NLTE correction

    do n=1,l_levels
        plevnlte(n)=plev(n)*100.0  ! converting from mbar to Pa
    end do
    plevnlte2 = plevnlte(l_levels:1:-1)  ! Reverse array for nltecool
    tl2 = tl(l_levels:1:-1)

    call nltecool(l_levels,plevnlte2,tl2,htrt3)

    htrt = htrt3(l_levels:1:-1)


    if (ames_15band) then
        do k=1,l_layers
            n=2*k+2
            heatingir15(k) = lw_heating_spec(3,k)+lw_heating_spec(4,k)+ &
                         lw_heating_spec(5,k)
            alpha(k) = 1.0/(1.0+(plev(n)/pnlteRef)**4.0)
            q4(k) = (alpha(k)*htrt(n))+((1-alpha(k))*heatingir15(k))
            lw_heating(k) = lw_heating(k)-heatingir15(k)+q4(k)
            htrt2(k) = htrt(n)
        end do
    else

        lw_heatingIN = lw_heating
        do k=1,l_layers
            n=2*k+2
            heatingir15(k) = lw_heating_spec(3,k)
            alpha(k) = 1.0/(1.0+(plev(n)/pnlteRef)**4.0)
            q4(k) = (alpha(k)*htrt(n))+((1-alpha(k))*heatingir15(k))
            lw_heating(k) = lw_heating(k)-heatingir15(k)+q4(k)
            htrt2(k) = htrt(n)
        end do
    endif !15 band
endif !IR NLTE correction
#endif

! shortwave (solar heating) NLTE correction
irtot0 = gcp*(FMNETI(1) - NFLUXTOPI)/scalep
swtot0 = gcp*(FMNETV(1) - NFLUXTOPV)/scalep
tstrat_dt = irtot0+swtot0*2.2e4*plev(2)/ &
         (1.0+2.2e4*plev(2))
do k=1,l_layers
    n=2*k+2
    sw_heating(k) = sw_heating(k)*2.2e4*plev(n)/ &
         (1.0+2.2e4*plev(n))
enddo

if (present(taudust_diagARG)) taudust_diagARG = taudust_diag
if (present(taucloud_diagARG)) taucloud_diagARG = taucloud_diag
if (present(tauco2cloud_diagARG)) tauco2cloud_diagARG = tauco2cloud_diag
if (present(detauARG)) detauARG = detau
if (present(outsolARG)) outsolARG = outsol
if (present(diffvtARG)) diffvtARG = diffvt
if (present(directsolARG)) directsolARG = directsol
if (present(fluxupvARG)) fluxupvARG = fluxupv
if (present(fluxdnvARG)) fluxdnvARG = fluxdnv
if (present(fmnetvARG)) fmnetvARG = fmnetv
if (present(fmnetiARG)) fmnetiARG = fmneti
if (present(fluxupiARG)) fluxupiARG = fluxupi
if (present(fluxdniARG)) fluxdniARG = fluxdni
if (present(nfluxtopvARG)) nfluxtopvARG = nfluxtopv
if (present(nfluxtopiARG)) nfluxtopiARG = nfluxtopi
if (present(swtotARG)) swtotARG = swtot

return
end subroutine ames_rt_driver



subroutine tbright_diag(  kd, tau1d, sscat1d, gfac1d, tcol, tsfc, sfc_emiss, xwl, tbright )

!=======================================================================
!  diagnostic brightness temperature
!=======================================================================

 use radiation_util_mod, only : planck_func, planck_inv_func, lw_scattering

 integer ::  kd

 real, intent(in), dimension(kd+1)   :: tcol      ! kd+1 values
 real, intent(in)                    :: tsfc      ! surface temperature
 real, intent(in)                    :: sfc_emiss, xwl
 real, intent(in), dimension(1,kd)   :: tau1d

 real*8   :: sscat1d(1,kd)
 real*8   :: gfac1d(1,kd)



!       OUTPUT:
 real, intent(out)     :: tbright

!       LOCAL:
integer, parameter :: NBANDS = 3
integer     :: id, jd,  kp, nn
real        :: cose

real, dimension(NBANDS)  :: xcenter


real, dimension( kd+1 )                       :: bh1d, fuh1d, fdh1d

real, dimension(1)  ::   sfc_emiss_0, tsfc1d, bsfc1d, radiance1d


kp= kd+1

!          Form the cose of the emission angle
!   cose= cos( (ema*2.0*pi/360.0) )
cose= 1.0

xcenter(1)= 7.8
xcenter(2)= 25.0
xcenter(3)= 32.0


#ifdef SKIP
!   Note that the dust and ice extinctions have been normalized
!          with respect to the visible.

!      Dust optical properties:  2.0 micron radius  @  7.8 25 mm and 32 mm

!  qex_dust(1)= 0.1446
qex_dust(1)= 0.0
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

#endif SKIP


!   xwl= xcenter(nn)   !   microns

!      Assume vertical emission: appropriate for a nadir instrument  (ie TES)
!            For 32 micron channel (MCS), assume a slant path
    if( xwl > 30.0 ) then
        cose= 0.3746
    else
        cose= 1.000
    endif

    tsfc1d= tsfc

!          Apply emissivity function for wavelengths >  18 microns
    if( xwl < 18.0 ) then
        sfc_emiss_0= 1.0
    else
        sfc_emiss_0= sfc_emiss
    endif

    call planck_func( kd+1, tcol,   xwl, bh1d   )
    call planck_func( 1,    tsfc1d, xwl, bsfc1d )

    call lw_scattering ( 1, kd, tau1d, sscat1d, gfac1d, sfc_emiss_0, &
                               bh1d, bsfc1d, fuh1d, fdh1d, cose  )

!   set a minimum value to avoid division by very tiny number
    radiance1d= MAX(1.e-15,fuh1d(1) / (2.0*pi ))

!    Convert to a brightness temperature

    call planck_inv_func( 1, radiance1d, xwl, tsfc1d  )

    tbright= tsfc1d(1)



return
end subroutine tbright_diag





#ifdef SKIP





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




#endif SKIP



!=====================================================================
!=====================================================================

subroutine ini_optcld(Qxv,Qxi,Qsv,Qsi,gvini,giini, &
             Qextref,TAUREFCLDini)
!  initialize cloud optical properties

implicit none

!  Arguments
!  ---------

real*8  Qxv(L_LEVELS+1,L_NSPECTV)
real*8  Qsv(L_LEVELS+1,L_NSPECTV)
real*8  gvini(L_LEVELS+1,L_NSPECTV)

real*8  Qxi(L_LEVELS+1,L_NSPECTI)
real*8  Qsi(L_LEVELS+1,L_NSPECTI)
real*8  giini(L_LEVELS+1,L_NSPECTI)

real*8  Qextref(L_LEVELS+1)
real*8  TAUREFCLDini(L_LEVELS+1)

integer nlonv,nloni

!  Local variables
!  ---------------

integer i,k

! Initialyze various variables
! ----------------------------

nlonv = l_nspectv
nloni = l_nspecti

do k = 1, l_levels+1
    qextref(k)    = 1.
enddo

call settozero(l_levels+1,taurefcldini)

call settozero(nlonv*(l_levels+1),qxv)
call settozero(nlonv*(l_levels+1),qsv)
call settozero(nlonv*(l_levels+1),gvini)

call settozero(nloni*(l_levels+1),qxi)
call settozero(nloni*(l_levels+1),qsi)
call settozero(nloni*(l_levels+1),giini)

return

end subroutine ini_optcld

!=====================================================================
!=====================================================================

subroutine ini_optdst(GrefV,GrefI, &
               Qxv,Qxi,Qsv,Qsi,gvini,giini,Qextref,surf2d, fixed )
!  initialize dust optical properties

implicit none

!  Arguments
!  ---------

real*8, intent(in) ::  GrefV(L_NSPECTV)
real*8, intent(in) ::  GrefI(L_NSPECTI)
logical, intent(in), optional :: fixed

real*8, intent(out) ::    Qxv(L_LEVELS+1,L_NSPECTV)
real*8, intent(out) ::    Qsv(L_LEVELS+1,L_NSPECTV)
real*8, intent(out) ::  gvini(L_LEVELS+1,L_NSPECTV)

real*8, intent(out) ::    Qxi(L_LEVELS+1,L_NSPECTI)
real*8, intent(out) ::    Qsi(L_LEVELS+1,L_NSPECTI)
real*8, intent(out) ::  giini(L_LEVELS+1,L_NSPECTI)

real*8, intent(out) ::  Qextref(L_LEVELS+1)

real*8,  intent(out), optional  ::  surf2d(nbin_rt,l_layers)


!  Local variables
!  ---------------

integer i,k

real*8 dev2
real*8 rs
real*8 surf(nbin_rt)

real*8 derf

! Initialyze various variables
! ----------------------------
if(.not.(present(fixed))) then

    do k = 1, l_levels+1
        qextref(k) = qextv(l_nrefv)
    enddo

    do i = 1, nlonv
        do k = 1, l_levels+1
            qxv(k,i) = qextv(i)
            qsv(k,i) = qscatv(i)
            gvini(k,i)  = ggv(i)

        enddo
    enddo
    do i = 1, nloni
        do k = 1, l_levels+1
            qxi(k,i) = qexti(i)
            qsi(k,i) = qscati(i)
            giini(k,i)  = ggi(i)

        enddo
    enddo

    if( present(surf2d) )  surf2d(:,:)= 0.0


else   !     --------------
    if (user_fixed_dust_opts) then
        qextref(:) = qxv_read(l_nrefv)

        do i = 1, nlonv
            do k = 1, l_levels+1
                qxv(k,i) = qxv_read(i)
                qsv(k,i) = qsv_read(i)
                gvini(k,i)  = gv_read(i)
            enddo
        enddo
        do i = 1, nloni
            do k = 1, l_levels+1
                qxi(k,i) = qxi_read(i)
                qsi(k,i) = qsi_read(i)
                giini(k,i)  = gi_read(i)
            enddo
        enddo

        if( present(surf2d) )  surf2d(:,:)= 0.0
    else
        
        ! Calculate Qext, Qscat, g for fixed dust
        dev2 = 1. / ( sqrt(2.)*dev_dt )
        Rs = 0.
        surf = 0.

        Rs = min( max(reff_fixed,1.e-7) , 50.e-6 )
        Rs = 1. / Rs

        do i = 1, nbin_rt
            surf(i) = 0.5 * ( derf( dlog(radb_rt(i+1)*Rs) * dev2 )     &
             -derf( dlog(radb_rt(i)  *Rs) * dev2 ) )
        enddo

    !      Collect size-bin expansion for each layer :   surf2d(nbin_rt,l_layers)

        if( present(surf2d) ) then
            do k= 1, l_layers
              surf2d(:,k)= surf(:)
            enddo
        endif

        call calc_opts(nlon=nlonv,nbin=nbin_rt, &
                        qext_in=qextv_dst,qscat_in=qscatv_dst,g_in=gv_dst, &
                        surf=surf,qx_out=qxv(1,:),qs_out=qsv(1,:),g_out=gvini(1,:))
        do k=1,l_levels+1
            qxv(k,:)=qxv(1,:)
            qsv(k,:)=qsv(1,:)
            gvini(k,:)=gvini(1,:)
        enddo

        call calc_opts(nlon=nloni,nbin=nbin_rt, &
                        qext_in=qexti_dst,qscat_in=qscati_dst,g_in=gi_dst, &
                        surf=surf,qx_out=qxi(1,:),qs_out=qsi(1,:),g_out=giini(1,:))

        do k=1,l_levels+1
            qxi(k,:)=qxi(1,:)
            qsi(k,:)=qsi(1,:)
            giini(k,:)=giini(1,:)
        enddo

        qextref(:) = qxv(:,l_nrefv)
        
    endif ! user_fixed_dust_opts

endif


return

end subroutine ini_optdst

!=====================================================================
!=====================================================================

subroutine fillpt(pl,ptrop,tg,tstrat,tl,plev,tlev,pmid,   &
                tmid, delp)
!  Put the T & P GCM arrays onto the NRC grid:  PLEV, PMID, TLEV, TMID
!  PMID and TMID are the pressure and temperature at the GCM layer
!  mid-points.  PLEV and TLEV are the pressures and temperatures at
!  the GCM layer boundaries, i.e. at GCM levels.


implicit none

integer :: k, l
real*8  :: pl(l_levels), plev(l_levels), pmid(l_levels)
real*8  :: tstrat, tl(l_levels), tlev(l_levels), tmid(l_levels)
real*8  :: ptrop, tg
real*8  :: delp(l_layers)

!======================================================================C

!  Fill the new radiation code variables.
!  PLEV and TLEV are the pressure and tempertures on a vertical grid
!  that the new radiation code uses.

plev(2) = ptrop/2.0d0
plev(1) = 1.e-5*plev(2)  !ptrop/2.0d0

do k=3,l_levels
    plev(k) = pl(k)
end do

do k=1,3
    tlev(k) = tstrat
end do

do k=4,l_levels-1,2
    tlev(k) = tl(k)
end do

do k=5,l_levels-2,2
    tlev(k) = tlev(k+1) + (tlev(k-1)-tlev(k+1))*                   &
              dlog(plev(k)/plev(k+1))/                             &
              dlog(plev(k-1)/plev(k+1))
end do

!  Temperature of the bottom level is the ground temperature.

tlev(l_levels) = tg

!  fill the pmid & tmid arrays used by optci and optcv subroutines.
!  tmid and pmid used to get the index for co2 k-coefficient
!  interpolation.

tmid(:) = 0.d0
pmid(:) = 0.d0

tmid(1) = tlev(2)
tmid(2) = tlev(2)
pmid(1) = ptrop/2.0  !plev(1)
pmid(2) = ptrop/2.0  !plev(2)

do l=1,l_layers
    tmid(2*l+1) = tlev(2*l+1)
    tmid(2*l+2) = tlev(2*l+1)
    pmid(2*l+1) = plev(2*l+1)
    pmid(2*l+2) = plev(2*l+1)
    delp(l)=    plev(2*l+3)-plev(2*l+1)
end do

tmid(l_levels) = tlev(l_levels)
pmid(l_levels) = plev(l_levels)



return
end subroutine fillpt

!=====================================================================
!=====================================================================

subroutine opt_dst(qtrace,pl,   &
                Qxv,Qxi,Qsv,Qsi,gv,gi,   &
                Qextref,TAUREF, surf2d,  &
                taudust,taudust_split)
!  calculate dust opacities

implicit none

!  Arguments
!  ---------

real*8  qtrace(l_layers,ntrace)
real*8  pl(l_levels)

real*8  qxv(l_levels+1,l_nspectv)
real*8  qsv(l_levels+1,l_nspectv)
real*8  gv(l_levels+1,l_nspectv)

real*8  qxi(l_levels+1,l_nspecti)
real*8  qsi(l_levels+1,l_nspecti)
real*8  gi(l_levels+1,l_nspecti)

real*8  qextref(l_levels+1)
real*8  tauref(l_levels+1)

real*8  surf2d(nbin_rt,l_layers)

real*8  taudust(2)
real*8  taudust_split(ntrace,2)


!  local variables
!  ---------------

integer i,iwav,j,k,l
integer nn

logical do_it(ndust_mass)

real*8 dens,cst,dev,dev2
real*8 rn,rs,ao,mo,no
real*8 surf(nbin_rt)
real*8 mass

real*8 derf

integer ndx,nd
real*8 Ntot(l_layers)
real*8 qxv_split(ndust_mass,l_levels+1,l_nspectv)
real*8 qsv_split(ndust_mass,l_levels+1,l_nspectv)
real*8 qxi_split(ndust_mass,l_levels+1,l_nspecti)
real*8 qsi_split(ndust_mass,l_levels+1,l_nspecti)
real*8 tautmp,sratio
real*8 Ao_split(ndust_mass)
real*8 Rs_split(ndust_mass)
real*8 surf_split(ndust_mass,nbin_rt)
! Initialize various variables
! ----------------------------
do_it = .false.
do k = 1, l_levels+1
    qextref(k) = 1.
    tauref(k)  = 0.
enddo

do i = 1, nlonv
    do k = 1, l_levels+1
        qxv(k,i) = qextref(k)
        qsv(k,i) = qextref(k) * 0.99
        gv(k,i)  = 0.

        qxv_split(:,k,i) = qextref(k)
        qsv_split(:,k,i) = qextref(k) * 0.99
    enddo
enddo
do i = 1, nloni
    do k = 1, l_levels+1
        qxi(k,i) = qextref(k)
        qsi(k,i) = qextref(k) * 0.99
        gi(k,i)  = 0.

        qxi_split(:,k,i) = qextref(k)
        qsi_split(:,k,i) = qextref(k) * 0.99
    enddo
enddo
Ntot=0.
do nd = 1,ndust_mass
    ndx = dust_mass_indx(nd)
    if (ndx .gt. nt_nontag) exit
    Ntot(:)=Ntot(:)+qtrace(:,ndx+1)
enddo

dev  = dev_dt
dens = dpden_dt
cst  = 0.75 / (pi*dens)
dev2 = 1. / ( sqrt(2.)*dev )

!  Treatment
!  ---------

taudust(:) = 0.
taudust_split(:,:) = 0.

DO L = 1, L_LAYERS
    Ao = 0.
    Ao_split = 0.
    surf = 0.
    surf_split=0.
    Rs_split = 0.
    do nd = 1, ndust_mass
        ndx = dust_mass_indx(nd)
        if (ndx .gt. nt_nontag) exit      ! exit dust_mass loop if entering tag tracers

        do_it(nd) = qtrace(l,ndx) .gt. 1.e-12  .and.  qtrace(l,ndx+1) .gt. 1.

        if (do_it(nd)) then

            !     Get the cross-section mean radius (Rs) of the log-normal distribution
            Mo = qtrace(l,ndx)          ! Mass mixing ratio
            No = qtrace(l,ndx+1)        ! Number mixing ratio

            Rs = ( Mo/No*cst )**(athird) * dexp( -0.5*dev**2. )

            !     Get the total cross sectional area Ao of particles
            Ao_split(nd) = No * pi * Rs**2.
            Ao = Ao + Ao_split(nd)

            ! *********************************************************************
            !     Define the cross-section weighted distribution, i.e. surface/size
            !    bin.  Change Rs to Reff.  MAK 6 May 2008.

            Rs = Rs * dexp ( 1.5 * dev**2. )

            ! *********************************************************************

            Rs = min( max(Rs,1.e-7) , 50.e-6 )
            Rs = 1. / Rs
            Rs_split(nd) = Rs

        endif
    end do
    do nd = 1, ndust_mass
        ndx = dust_mass_indx(nd)
        if (ndx .gt. nt_nontag) exit      ! exit dust_mass loop if entering tag tracers

        do_it(nd) = qtrace(l,ndx) .gt. 1.e-12 .and.  qtrace(l,ndx+1) .gt. 1.

        if (do_it(nd)) then
            sratio = Ao_split(nd)/Ao
            do i = 1, nbin_rt
                surf_split(nd,i) = 0.5 * ( derf( dlog(radb_rt(i+1)*Rs_split(nd)) * dev2 )     &
                             -derf( dlog(radb_rt(i)  *Rs_split(nd)) * dev2 ) )
                surf(i) = surf(i) + sratio*surf_split(nd,i)
            enddo
        endif    ! end condition on do_it
    enddo    ! end loop over dust_mass

!      Collect size-bin expansion for each layer :   surf2d(nbin_rt,l_layers)
    surf2d(:,l)= surf(:)

    if (ANY(do_it)) then
    !     get the average values of <qext>, <qscat>, and <g> for the whole distribution.
    do nn = 1, 2
        k = 2*l + 1 + nn
        call calc_opts(nlon=nlonv,nbin=nbin_rt, &
                    qext_in=qextv_dst,qscat_in=qscatv_dst,g_in=gv_dst, &
                    surf=surf,qx_out=qxv(k,:),qs_out=qsv(k,:),g_out=gv(k,:), &
                    surf_split=surf_split,qx_split=qxv_split(:,k,:), &
                    qs_split=qsv_split(:,k,:))

        call calc_opts(nlon=nloni,nbin=nbin_rt, &
                    qext_in=qexti_dst,qscat_in=qscati_dst,g_in=gi_dst, &
                    surf=surf,qx_out=qxi(k,:),qs_out=qsi(k,:),g_out=gi(k,:), &
                    surf_split=surf_split,qx_split=qxi_split(:,k,:), &
                    qs_split=qsi_split(:,k,:))


        qextref(k) = qxv(k,l_nrefv)
        mass       = 100. * (pl(k) - pl(k-1)) / grav
        tauref(k)  = ao * qextref(k) * mass

        !       for diagnostics: opacity at ref wavelengths (in the vis and ir)
        taudust(1) = taudust(1) + tauref(k)
        taudust(2) = taudust(2) + tauref(k)/qextref(k) *      &
                      ( qxi(k,l_nrefi)-qsi(k,l_nrefi) )

        do nd = 1, ndust_mass
            tautmp = 0.d0
            if (do_it(nd)) then
                ndx = dust_mass_indx(nd)
                if (ndx .gt. nt_nontag) exit      ! exit dust_mass loop if entering tag tracers
                tautmp = Ao_split(nd) * qxv_split(nd,k,l_nrefv) * mass

                taudust_split(ndx,1) = taudust_split(ndx,1) + tautmp
                taudust_split(ndx,2) = taudust_split(ndx,2) +  &
                                    tautmp / qxv_split(nd,k,l_nrefv) * &
                                    ( qxi_split(nd,k,l_nrefi)-qsi_split(nd,k,l_nrefi) )
            endif
        enddo
    enddo
    endif
enddo    ! end loop over layers

return

end subroutine opt_dst


!=====================================================================

subroutine opt_cld(qtrace,pl,   &
                Qxv,Qxi,Qsv,Qsi,gv,gi, &
                Qextrefcld,TAUREFCLD, &
                taucloud,taurefcld_uv)
!=====================================================================

!  calculate cloud opacities


implicit none

!  Arguments
!  ---------


real*8  qtrace(:,:)
real*8  pl(l_levels)

real*8  qxv(l_levels+1,l_nspectv)
real*8  qsv(l_levels+1,l_nspectv)
real*8  gv(l_levels+1,l_nspectv)

real*8  qxi(l_levels+1,l_nspecti)
real*8  qsi(l_levels+1,l_nspecti)
real*8  qbi(l_levels+1,l_nspecti)
real*8  gi(l_levels+1,l_nspecti)

real*8  qextrefcld(l_levels+1)
real*8  taurefcld(l_levels+1)
real*8  taurefcld_uv(l_levels+1)

real*8  taucloud(4)


!  local variables
!  ---------------

integer i,iwav,j,k,l,irap
integer nn

real*8 dens,cst,dev,dev2
real*8 rn,rs,ao,mo,no
real*8 surf(nbin_rt)
real*8 mass

real*4 mantletocore
logical do_it

real*8 derf

! Initialize various variables
! ----------------------------
do_it = .false.
do k = 1, l_levels+1
    qextrefcld(k) = 1.
    taurefcld(k) = 0.
    taurefcld_uv(k) = 0.
enddo

do i = 1, nlonv
    do k = 1, l_levels+1
        qxv(k,i) = qextrefcld(k)
        qsv(k,i) = qextrefcld(k) * 0.99
        gv(k,i)  = 0.
    enddo
enddo

do i = 1, nloni
    do k = 1, l_levels+1
        qxi(k,i) = qextrefcld(k)
        qsi(k,i) = qextrefcld(k) * 0.99
        gi(k,i)  = 0.
    enddo
enddo

!  Treatment
!  ---------

dev= dev_ice
dev2 = 1. / ( sqrt(2.)*dev )

taucloud(:) = 0.

do l = 1, l_layers

    !     do not do the computations if the amount of cloud is too low
    do_it = qtrace(l,ima_cld)+qtrace(l,ima_cor).gt. 1.e-7   &
    .and.  qtrace(l,inb_cld) .gt. 1.

    if (do_it) then

        !     Determine the ratio of the dust core radius over that of the ice mantle
        mantletocore =   qtrace(l,iMa_cor)/dpden_dt      &
                /( qtrace(l,iMa_cld)/dpden_ice     &
                  +qtrace(l,iMa_cor)/dpden_dt )
        mantletocore = mantletocore**(athird)

        !     Find the index to which corresponds the optical properties of the
        !     core to mantle radius ratio. Those properties were determined off-line
        !     using the Toon and Ackerman coated spheres code.
        irap = nratio
        do i = 1, nratio
            if (mantletocore.lt.cor_ratio(i) .and. i.ne.1) then
                irap = i - 1
                exit
            elseif (mantletocore.eq.cor_ratio(i)) then
                irap = i
                exit
            elseif (mantletocore.lt.cor_ratio(i) .and. i.eq.1) then
                irap = 1
                exit
            endif
        enddo

        !     Get the cross-section mean radius (Rs) of the log-normal distribution
        mo = qtrace(l,ima_cld)       &
          +qtrace(l,ima_cor)       ! mass mixing ratio
        no = qtrace(l,inb_cld)       ! number mixing ratio

        dens =  qtrace(l,ima_cld) / mo * dpden_ice    &
             +qtrace(l,ima_cor) / mo * dpden_dt

        cst  = 0.75 / (pi*dens)
        rs = ( mo/no*cst )**(athird) * dexp( -0.5*dev**2. )

        !     get the total cross sectional area ao of the water ice particles
        ao = no * pi * rs**2.

        ! *********************************************************************

        !     define the cross-section weighted distribution, i.e. surface/size
        !   bin.   change rs to reff.  mak 6 may 2008.

        rs = rs * dexp ( 1.5 * dev**2. )

        ! *********************************************************************

        rs = min( max(rs,1.e-7) , 100.e-6 )
        rs = 1. / rs

        do i = 1, nbin_rt
            surf(i) = 0.5 * ( derf( dlog(radb_rt(i+1)*rs) * dev2 )     &
                         -derf( dlog(radb_rt(i)  *Rs) * dev2 ) )
        enddo

        !     get the average values of <qext>, <qscat>, and <g> for the whole distribution.
        do nn = 1, 2
            k = 2*l + 1 + nn

            do iwav = 1, nlonv
                qxv(k,iwav) = 0.
                qsv(k,iwav) = 0.
                do i = 1, nbin_rt
                    qxv(k,iwav) = qxv(k,iwav)+ surf(i) * qextv_cld(irap,i,iwav)
                    qsv(k,iwav) = qsv(k,iwav)+ surf(i) * qscatv_cld(irap,i,iwav)
                    gv(k,iwav)  = gv(k,iwav) + surf(i) * gv_cld(irap,i,iwav)
                enddo
                qsv(k,iwav) = min( qsv(k,iwav) , 0.99999*qxv(k,iwav) )
            enddo

            do iwav = 1, nloni
                qxi(k,iwav) = 0.
                qsi(k,iwav) = 0.
                qbi(k,iwav) = 0.
                do i = 1, nbin_rt
                    qxi(k,iwav) = qxi(k,iwav)+ surf(i) * qexti_cld(irap,i,iwav)
                    qsi(k,iwav) = qsi(k,iwav)+ surf(i) * qscati_cld(irap,i,iwav)
                    gi(k,iwav)  = gi(k,iwav) + surf(i) * gi_cld(irap,i,iwav)
                    qbi(k,iwav) = qbi(i,iwav)+ surf(i) * fcorrect(i,irap) &
                                * (qexti_cld(irap,i,iwav)-qscati_cld(irap,i,iwav))
                enddo
                qsi(k,iwav) = min( qsi(k,iwav) , 0.99999*qxi(k,iwav) )
            enddo

            qextrefcld(k) = qxv(k,l_nrefv)
            mass          = 100. * (pl(k) - pl(k-1)) / grav
            taurefcld(k)  = ao * qextrefcld(k) * mass
            taurefcld_uv(k) = ao * qxv(k,l_nspectv) * mass

            !       for diagnostics: cloud opacity at ref wavelengths (in the vis and ir)
            taucloud(1) = taucloud(1) + taurefcld(k)
            taucloud(2) = taucloud(2)              &
                                   +taurefcld(k) / qextrefcld(k) *     &
                                ( qxi(k,l_nrefi) - qsi(k,l_nrefi) )
            taucloud(3) = taucloud(3) + taurefcld(k) &
                                / qextrefcld(k) * qbi(k,l_nrefi)
            taucloud(4) = taucloud(4) + taurefcld_uv(k)
        enddo

    endif    ! condition on cloud amount

enddo    ! end loop over layers

return

end subroutine opt_cld

!=====================================================================
!=====================================================================
subroutine opt_cld_simple(cldice_bin,pl,       &
                    Qxv,Qxi,Qsv,Qsi,gv,gi,   &
                    Qextrefcld,TAUREFCLD,    &
                    surfcld,                 &
                    taucloud, taurefcld_uv    )
!
!    Simple cloud opacity calculation
!

implicit none

!  Arguments
!  ---------

real*8,  intent(in), dimension(l_layers)  ::  cldice_bin
real*8,  intent(in), dimension(l_levels)  ::  pl

real*8  Qxv(L_LEVELS+1,L_NSPECTV)
real*8  Qsv(L_LEVELS+1,L_NSPECTV)
real*8  gv (L_LEVELS+1,L_NSPECTV)

real*8  Qxi(L_LEVELS+1,L_NSPECTI)
real*8  Qsi(L_LEVELS+1,L_NSPECTI)
real*8  gi (L_LEVELS+1,L_NSPECTI)
real*8  qbi(l_levels+1,l_nspecti)

real*8  Qextrefcld(L_LEVELS+1)
real*8  TAUREFCLD (L_LEVELS+1)
real*8  TAUREFCLD_uv (L_LEVELS+1)
real*8  surfcld( nbin_rt,l_layers )

real*8  taucloud(4)


!  Local variables
!  ---------------

integer i,iwav,j,k,l,nn, irap, ibin

real*8 Rn,Rs,ao,mo,no, Mass, ssize
real*8 cst,dev

! Initialyze various variables
! ----------------------------

do k = 1, l_levels+1
    qextrefcld(k) = 1.
    taurefcld(k) = 0.
    taurefcld_uv(k) = 0.
enddo

do i = 1, nlonv
    do k = 1, l_levels+1
        qxv(k,i) = qextrefcld(k)
        qsv(k,i) = qextrefcld(k) * 0.99
        gv(k,i)  = 0.
    enddo
enddo

do i = 1, nloni
    do k = 1, l_levels+1
        qxi(k,i) = qextrefcld(k)
        qsi(k,i) = qextrefcld(k) * 0.99
        gi(k,i)  = 0.
    enddo
enddo

!  Treatment
!  ---------

dev= dev_ice

taucloud(1) = 0.
taucloud(2) = 0.

irap = 1

ssize= simple_cloud_radius * 1.0e-6

do nn= 1, nbin_rt
    if( rad_rt(nn) > ssize ) then
        ibin= nn
        exit
    endif
enddo

ao= 0.75 / ( dpden_ice * ssize )

surfcld(:,:)= 0.0
surfcld(ibin,:)= 1.0


do l = 1, l_layers

!     get the average values of <qext>, <qscat>, and <g> for the whole distribution

    do nn = 1, 2
        k = 2*l + 1 + nn

        do iwav = 1, nlonv
            qxv(k,iwav) = qextv_cld (irap,ibin,iwav)
            qsv(k,iwav) = qscatv_cld(irap,ibin,iwav)
            gv (k,iwav)  = gv_cld   (irap,ibin,iwav)
            qsv(k,iwav) = min( qsv(k,iwav) , 0.99999*qxv(k,iwav) )
        enddo

        do iwav = 1, nloni
            qxi(k,iwav) = qexti_cld (irap,ibin,iwav)
            qsi(k,iwav) = qscati_cld(irap,ibin,iwav)
            gi (k,iwav) = gi_cld    (irap,ibin,iwav)
            qsi(k,iwav) = min( qsi(k,iwav) , 0.99999*qxi(k,iwav) )
        enddo

        qextrefcld(k) = qxv(k,l_nrefv)
        mass          = 100. * (pl(k) - pl(k-1)) / grav
        taurefcld(k)  = ao * cldice_bin(l) * qextrefcld(k) * mass
        taurefcld_uv(k) = ao * cldice_bin(l) * qxv(k,l_nspectv) * mass

        !       for diagnostics: cloud opacity at ref wavelengths (in the vis and ir)
        taucloud(1) = taucloud(1) + taurefcld(k)
        taucloud(2) = taucloud(2) + taurefcld(k) / qextrefcld(k) *     &
                                               ( qxi(k,l_nrefi) - qsi(k,l_nrefi) )
    enddo

enddo    ! end loop over layers

return

end subroutine opt_cld_simple


!=====================================================================

subroutine opt_cld_blk(qtrace,pl,tl,   &
                Qxv,Qxi,Qsv,Qsi,gv,gi, &
                Qextrefcld,TAUREFCLD, &
                taucloud,nice,mcpu0,taurefcld_uv)
!=====================================================================

!  calculate cloud opacities for bulk scheme


implicit none

!  Arguments
!  ---------


real*8  qtrace(:,:)
real*8  pl(l_levels)
real*8  tl(l_levels)

real*8  qxv(l_levels+1,l_nspectv)
real*8  qsv(l_levels+1,l_nspectv)
real*8  gv(l_levels+1,l_nspectv)

real*8  qxi(l_levels+1,l_nspecti)
real*8  qsi(l_levels+1,l_nspecti)
real*8  gi(l_levels+1,l_nspecti)
real*8  qbi(l_levels+1,l_nspecti)

real*8  qextrefcld(l_levels+1)
real*8  taurefcld(l_levels+1)
real*8  taurefcld_uv(l_levels+1)

real*8  taucloud(4)

integer nice

logical mcpu0

!  local variables
!  ---------------

integer i,iwav,j,k,l,irap,irad
integer nn

real*8 dens,cst,dev,dev2
real*8 rn,rs,ao,mo,no
real*8 surf(nbin_rtblk)
real*8 mass

real*4 mantletocore
logical do_it

real*8 derf

! Initialize various variables
! ----------------------------
do_it = .false.
do k = 1, l_levels+1
    qextrefcld(k) = 1.
    taurefcld(k) = 0.
    taurefcld_uv(k) = 0.
enddo

do i = 1, nlonv
    do k = 1, l_levels+1
        qxv(k,i) = qextrefcld(k)
        qsv(k,i) = qextrefcld(k) * 0.99
        gv(k,i)  = 0.
    enddo
enddo

do i = 1, nloni
    do k = 1, l_levels+1
        qxi(k,i) = qextrefcld(k)
        qsi(k,i) = qextrefcld(k) * 0.99
        gi(k,i)  = 0.
    enddo
enddo

!  Treatment
!  ---------

dev= dev_ice
dev2 = 1. / ( sqrt(2.)*dev )

taucloud(1) = 0.
taucloud(2) = 0.

!  Assume smallest core to ice ratio for bulk scheme
irap= nratioblk

do l = 1, l_layers

    !     Get the cross-section mean radius (Rs) of the log-normal distribution
    mo = qtrace(l,nice)  
    no = bulk_ccn       ! number mixing ratio

    cst  = 0.75 / (pi * dpden_ice)
    rs = ( mo/no*cst )**(athird) !* dexp( -0.5*dev**2. )

    do_it = Rs .gt. 1.e-7

    if (do_it) then

        !     get the total cross sectional area ao of the water ice particles
        ao = no * pi * rs**2.

        ! assume single particle size for bulk scheme 

        !! *********************************************************************
        !
        !!     define the cross-section weighted distribution, i.e. surface/size
        !!   bin.   change rs to reff.  mak 6 may 2008.
        !
        !rs = rs * dexp ( 1.5 * dev**2. )
        !
        !! *********************************************************************

        rs = min( max(rs,1.e-7) , 100.e-6 )
!        rs = 1. / rs
!
!        do i = 1, nbin_rtblk
!            surf(i) = 0.5 * ( derf( dlog(radb_rtblk(i+1)*rs) * dev2 )     &
!                         -derf( dlog(radb_rtblk(i)  *Rs) * dev2 ) )
!        enddo


        irad=nbin_rtblk-1
        do i = 1, nbin_rtblk
          if (Rs.le.rad_rtblk(i) .and. i.ne.1) then
            irad = i - 1
            exit
          elseif (Rs.lt.rad_rtblk(i) .and. i.eq.1) then
            irad = 1
            exit
          endif
        enddo

        do nn = 1, 2
            k = 2*l + 1 + nn

            if(tl(k).lt.273.) then

            ! Use optical constants for ice
            do iwav = 1, nlonv
              qxv(k,iwav) = 0.
              qsv(k,iwav) = 0.
              qxv(k,iwav) = qextv_bcld(irap,irad+1,iwav)+ &
               (qextv_bcld(irap,irad,iwav)-  &
               qextv_bcld(irap,irad+1,iwav))*  &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              qsv(k,iwav) = qscatv_bcld(irap,irad+1,iwav)+  &
               (qscatv_bcld(irap,irad,iwav)-   & 
               qscatv_bcld(irap,irad+1,iwav))*   &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              gv(k,iwav)  = gv_bcld(irap,irad+1,iwav)+   &
               (gv_bcld(irap,irad,iwav)-gv_bcld(irap,irad+1,iwav))*   & 
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1)) 
              qsv(k,iwav) = min( qsv(k,iwav) , 0.99999*qxv(k,iwav) )
            enddo

            do iwav = 1, nloni
              qxi(k,iwav) = 0.
              qsi(k,iwav) = 0.
              qxi(k,iwav) = qexti_bcld(irap,irad+1,iwav)+  &
               (qexti_bcld(irap,irad,iwav)-  & 
               qexti_bcld(irap,irad+1,iwav))*  &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              qsi(k,iwav) = qscati_bcld(irap,irad+1,iwav)+  &
               (qscati_bcld(irap,irad,iwav)-  &
               qscati_bcld(irap,irad+1,iwav))*  &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              gi(k,iwav)  = gi_bcld(irap,irad+1,iwav)+  &
               (gi_bcld(irap,irad,iwav)-gi_bcld(irap,irad+1,iwav))*  &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              qsi(k,iwav) = min( qsi(k,iwav) , 0.99999*qxi(k,iwav) )
            enddo


            else
           
            !Use optical constants for liquid
            do iwav = 1, nlonv
              qxv(k,iwav) = 0.
              qsv(k,iwav) = 0.
              qxv(k,iwav) = qextv_blcld(irap,irad+1,iwav)+ &
               (qextv_blcld(irap,irad,iwav)-  &
               qextv_blcld(irap,irad+1,iwav))*  &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              qsv(k,iwav) = qscatv_blcld(irap,irad+1,iwav)+  &
               (qscatv_blcld(irap,irad,iwav)-   &
               qscatv_blcld(irap,irad+1,iwav))*   &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              gv(k,iwav)  = gv_blcld(irap,irad+1,iwav)+   &
               (gv_blcld(irap,irad,iwav)-gv_blcld(irap,irad+1,iwav))*   &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              qsv(k,iwav) = min( qsv(k,iwav) , 0.99999*qxv(k,iwav) )
            enddo

            do iwav = 1, nloni
              qxi(k,iwav) = 0.
              qsi(k,iwav) = 0.
              qxi(k,iwav) = qexti_blcld(irap,irad+1,iwav)+  &
               (qexti_blcld(irap,irad,iwav)-  &
               qexti_blcld(irap,irad+1,iwav))*  &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              qsi(k,iwav) = qscati_blcld(irap,irad+1,iwav)+  &
               (qscati_blcld(irap,irad,iwav)-  &
               qscati_blcld(irap,irad+1,iwav))*  &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              gi(k,iwav)  = gi_blcld(irap,irad+1,iwav)+  &
               (gi_blcld(irap,irad,iwav)-gi_blcld(irap,irad+1,iwav))*  &
               log(rs/rad_rtblk(irad+1))/log(rad_rtblk(irad)/rad_rtblk(irad+1))
              qsi(k,iwav) = min( qsi(k,iwav) , 0.99999*qxi(k,iwav) )
            enddo

            end if

            qextrefcld(k) = qxv(k,l_nrefv)
            mass          = 100. * (pl(k) - pl(k-1)) / grav
            taurefcld(k)  = ao * qextrefcld(k) * mass
            taurefcld_uv(k)  = ao * qxv(k,l_nspectv) * mass


            !       for diagnostics: cloud opacity at ref wavelengths (in the vis and ir)
            taucloud(1) = taucloud(1) + taurefcld(k)
            taucloud(2) = taucloud(2)              &
                                   +taurefcld(k) / qextrefcld(k) *     &
                                ( qxi(k,l_nrefi) - qsi(k,l_nrefi) )
!            taucloud(3) = taucloud(3) + taurefcld(k) &
!                                / qextrefcld(k) * qbi(k,l_nrefi)
!            taucloud(4) = taucloud(4) + taurefcld_uv(k)
 
        enddo        

    endif    ! condition on cloud amount

enddo    ! end loop over layers


return

end subroutine opt_cld_blk


!=====================================================================
!=====================================================================
subroutine opt_co2cld(qtrace,pl,         &
                    Qxv,Qxi,Qsv,Qsi,gv,gi,   &
                    Qextrefco2cld,TAUREFCO2CLD,    &
                    tauco2cloud, nco2 )
!
!    CO2 cloud opacity calculation
!

implicit none

!  Arguments
!  ---------

real*8,  intent(in), dimension(l_levels)  ::  pl

real*8  qtrace(:,:)

real*8  Qxv(L_LEVELS+1,L_NSPECTV)
real*8  Qsv(L_LEVELS+1,L_NSPECTV)
real*8  gv (L_LEVELS+1,L_NSPECTV)

real*8  Qxi(L_LEVELS+1,L_NSPECTI)
real*8  Qsi(L_LEVELS+1,L_NSPECTI)
real*8  gi (L_LEVELS+1,L_NSPECTI)

real*8  Qextrefco2cld(L_LEVELS+1)
real*8  TAUREFCO2CLD (L_LEVELS+1)

real*8  tauco2cloud(2)

integer nco2

!  Local variables
!  ---------------

integer i,iwav,j,k,l,nn,irap

real*8 rn,rs,ao,mo,no,mass
real*8 dens,cst,dev,dev2

real*8  surf(nbin_rtco2)

real*4 mantletocore
logical do_it

real*8 derf


! Initialyze various variables
! ----------------------------

do k = 1, l_levels+1
    qextrefco2cld(k) = 1.
    taurefco2cld(k) = 0.
enddo


do i = 1, nlonv
    do k = 1, l_levels+1
        qxv(k,i) = qextrefco2cld(k)
        qsv(k,i) = qextrefco2cld(k) * 0.99
        gv(k,i)  = 0.
    enddo
enddo

do i = 1, nloni
    do k = 1, l_levels+1
        qxi(k,i) = qextrefco2cld(k)
        qsi(k,i) = qextrefco2cld(k) * 0.99
        gi(k,i)  = 0.
    enddo
enddo



!  Treatment
!  ---------

dev= dev_co2ice
dev2 = 1. / ( sqrt(2.)*dev )

tauco2cloud(1) = 0.
tauco2cloud(2) = 0.

!  Assume smallest core to ice ratio for bulk co2 cloud scheme
irap= nratioblk


do l = 1, l_layers

        !     Get the cross-section mean radius (Rs) of the log-normal distribution
    mo = qtrace(l,nco2)      
    no = bulk_co2_ccn

    cst  = 0.75 / (pi*dpden_co2ice)
    rs = ( mo/no*cst )**(athird) * dexp( -0.5*dev**2. )
    do_it = Rs .gt. 1.e-7

    if (do_it) then
        !     get the total cross sectional area ao of the water ice particles
        ao = no * pi * rs**2.

        ! *********************************************************************

        !     define the cross-section weighted distribution, i.e. surface/size
        !   bin.   change rs to reff.  mak 6 may 2008.

        rs = rs * dexp ( 1.5 * dev**2. )

        ! *********************************************************************

        rs = min( max(rs,1.e-7) , 100.e-6 )
        rs = 1. / rs

        do i = 1, nbin_rtco2
            surf(i) = 0.5 * ( derf( dlog(radb_rtco2(i+1)*rs) * dev2 )     &
                         -derf( dlog(radb_rtco2(i)  *Rs) * dev2 ) )
        enddo

        !     get the average values of <qext>, <qscat>, and <g> for the whole distribution.
        do nn = 1, 2
            k = 2*l + 1 + nn

            do iwav = 1, nlonv
                qxv(k,iwav) = 0.
                qsv(k,iwav) = 0.
                do i = 1, nbin_rtco2
                    qxv(k,iwav) = qxv(k,iwav)+ surf(i) * qextv_co2cld(irap,i,iwav)
                    qsv(k,iwav) = qsv(k,iwav)+ surf(i) * qscatv_co2cld(irap,i,iwav)
                    gv(k,iwav)  = gv(k,iwav) + surf(i) * gv_co2cld(irap,i,iwav)
                enddo
                qsv(k,iwav) = min( qsv(k,iwav) , 0.99999*qxv(k,iwav) )
            enddo

            do iwav = 1, nloni
                qxi(k,iwav) = 0.
                qsi(k,iwav) = 0.
                do i = 1, nbin_rtco2
                    qxi(k,iwav) = qxi(k,iwav)+ surf(i) * qexti_co2cld(irap,i,iwav)
                    qsi(k,iwav) = qsi(k,iwav)+ surf(i) * qscati_co2cld(irap,i,iwav)
                    gi(k,iwav)  = gi(k,iwav) + surf(i) * gi_co2cld(irap,i,iwav)
                enddo
                qsi(k,iwav) = min( qsi(k,iwav) , 0.99999*qxi(k,iwav) )
            enddo

            qextrefco2cld(k) = qxv(k,l_nrefv)
            mass          = 100. * (pl(k) - pl(k-1)) / grav
            taurefco2cld(k)  = ao * qextrefco2cld(k) * mass

            !       for diagnostics: cloud opacity at ref wavelengths (in the vis and ir)
            tauco2cloud(1) = tauco2cloud(1) + taurefco2cld(k)
            tauco2cloud(2) = tauco2cloud(2)              &
                                   +taurefco2cld(k) / qextrefco2cld(k) *     &
                                ( qxi(k,l_nrefi) - qsi(k,l_nrefi) )

        enddo


    endif    ! condition on cloud amount

enddo    ! end loop over layers

return

end subroutine opt_co2cld

!=====================================================================
!=====================================================================

subroutine optcv(dtauv,tauv,taucumv,plev,          &
               qxvdst,qsvdst,gvdst,wbarv,cosbv,        &
               tauref,tmid,pmid,taugsurf,qh2o,  &
               qextrefcld,taurefcld,qxvcld,qsvcld,gvcld)
!
! this subroutine sets the optical constants in the visible
! it calcualtes for each layer, for each specral interval in the visible
! layer: wbar, dtau, cosbar
! level: tau
!
! tauv(l,nw,ng) is the cumulative optical depth at the top of radiation code
! layer l. nw is spectral wavelength interval, ng the gauss point index.
!
!     tlev(l) - temperature at the layer boundary
!     plev(l) - pressure at the layer boundary (i.e. level)
!     co2v(nt,nps,nw,ng) - visible co2 k-coefficients
!
!----------------------------------------------------------------------c


implicit none


real*8  :: dtauv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: dtaukv(l_levels+1,l_nspectv,l_ngauss)
real*8  :: tauv(l_nlevrad,l_nspectv,l_ngauss)
real*8  :: taucumv(l_levels,l_nspectv,l_ngauss)
real*8  :: plev(l_levels)
real*8  :: tmid(l_levels), pmid(l_levels)



real*8  :: cosbv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: wbarv(l_nlayrad,l_nspectv,l_ngauss)


!     for dust
real*8  :: qxvdst(l_levels+1,l_nspectv)
real*8  :: qsvdst(l_levels+1,l_nspectv)
real*8  :: gvdst(l_levels+1,l_nspectv)
real*8  :: qextref(l_levels+1)
real*8  :: tauref(l_levels+1)
real*8  :: tauref_save(l_levels+1)

!     for clouds
real*8  :: qxvcld(l_levels+1,l_nspectv)
real*8  :: qsvcld(l_levels+1,l_nspectv)
real*8  :: gvcld(l_levels+1,l_nspectv)
real*8  :: qextrefcld(l_levels+1)
real*8  :: taurefcld(l_levels+1)
real*8  :: taurefcld_save(l_levels+1)

real*8  :: tcloud(l_levels,l_nspectv)

real*8  :: taureflk(l_levels+1,l_nspectv)
real*8  :: taucldk(l_levels+1,l_nspectv)

integer :: l, nw, ng, k, ng1(l_nspectv), lk
integer :: mt(l_levels), mp(l_levels), np(l_levels)
real*8  :: ans, taugas
real*8  :: tray(l_levels,l_nspectv)
real*8  :: taeros(l_levels,l_nspectv)
real*8  :: dpr(l_levels), u(l_levels)
real*8  :: lcoef(4), lkcoef(l_levels,4)

real*8  :: taugsurf(l_nspectv,l_ngauss-1), trayaer

!  reference wavelength is (now) bin #2 - put into qextref
!      real*8 qextref

!  water mixing ratio stuff

real*8  :: qh2o(l_levels),  wratio(l_levels)
real*8  :: kcoef(4)
integer :: nh2o(l_levels)

!======================================================================C

!  Save old tauref values

do k=1,l_levels+1
    tauref_save(k) = tauref(k)
    taurefcld_save(k) = taurefcld(k)
end do

!  determine the total gas opacity throughout the column, for each
!  spectral interval, nw, and each gauss point, ng.
!  calculate the continuum opacities, i.e., those that do not depend on
!  ng, the gauss index.

do ng=1,l_ngauss-1
    do nw=1,l_nspectv
        taugsurf(nw,ng) = 0.0d0
    end do
end do
tray = 0.0
taeros = 0.0
tcloud = 0.0
tauv = 0.0
dtauv = 0.0
taucumv = 0.0
dtaukv = 0.0

do k=2,l_levels
    dpr(k) = plev(k)-plev(k-1)
    u(k)   = cmk*dpr(k)

    call tpindex(pmid(k),tmid(k),qh2o(k),pfgasref,tgasref, &
             lcoef,mt(k),mp(k),nh2o(k),wratio(k))

    do lk=1,4
        lkcoef(k,lk) = lcoef(lk)
    end do

    qextref(k) = qxvdst(k,l_nrefv)

    tauref(k)    = tauref(k) / qextref(k)
    taurefcld(k) = taurefcld(k) / qextrefcld(k)

    do nw=1,l_nspectv
        tray(k,nw)   = tauray(nw)*dpr(k)
        taeros(k,nw) = tauref(k)    * qxvdst(k,nw)
        tcloud(k,nw) = taurefcld(k) * qxvcld(k,nw)
    end do
end do

!  TRAYAER is Tau RAYleigh scattering, plus AERosol opacity

do k=2,l_levels

    do nw=1,l_nspectv

        trayaer = tray(k,nw) + taeros(k,nw) + tcloud(k,nw)

        do ng=1,l_ngauss-1

            !           now compute taugas

            !  interpolate between water mixing ratios
            !  wratio = 0.0 if the requested water amount is equal to, or outside the
            !  the range of water amount data.

            kcoef(1) = co2v(mt(k),mp(k),nh2o(k),nw,ng) + wratio(k)*    &
                      (co2v(mt(k),mp(k),nh2o(k)+1,nw,ng) -             &
                       co2v(mt(k),mp(k),nh2o(k),nw,ng))

            kcoef(2) = co2v(mt(k),mp(k)+1,nh2o(k),nw,ng) + wratio(k)*  &
                      (co2v(mt(k),mp(k)+1,nh2o(k)+1,nw,ng) -           &
                       co2v(mt(k),mp(k)+1,nh2o(k),nw,ng))

            kcoef(3) = co2v(mt(k)+1,mp(k)+1,nh2o(k),nw,ng) + wratio(k)*&
                      (co2v(mt(k)+1,mp(k)+1,nh2o(k)+1,nw,ng) -         &
                       co2v(mt(k)+1,mp(k)+1,nh2o(k),nw,ng))

            kcoef(4) = co2v(mt(k)+1,mp(k),nh2o(k),nw,ng) + wratio(k)*  &
                      (co2v(mt(k)+1,mp(k),nh2o(k)+1,nw,ng) -           &
                       co2v(mt(k)+1,mp(k),nh2o(k),nw,ng))

            !  interpolate the co2 k-coefficients to the requested t,p


            ans = lkcoef(k,1)*kcoef(1) + lkcoef(k,2)*kcoef(2) +        &
                  lkcoef(k,3)*kcoef(3) + lkcoef(k,4)*kcoef(4)

            taugas          = u(k)*ans
            taugsurf(nw,ng) = taugsurf(nw,ng) + taugas
            dtaukv(k,nw,ng) = taugas + trayaer
        end do

        !  now fill in the "clear" part of the spectrum (ng = l_ngauss)
        !  which holds continuum opacity only

        ng = l_ngauss
        dtaukv(k,nw,ng) = taeros(k,nw)+tray(k,nw)+tcloud(k,nw)
    end do
end do

!  Now the full treatment for the layers, where besides the opacity
!  we need to calculate the scattering albedo and asymmetry factors
!  for each layer

do nw=1,l_nspectv
    do k=2,l_levels
        taureflk(k,nw) = tauref(k)    * qsvdst(k,nw)
        taucldk(k,nw)  = taurefcld(k) * qsvcld(k,nw)
    enddo
enddo

do nw=1,l_nspectv

    !  First, the special "clear" channel

    ng = l_ngauss
    do l=1,l_layers
        k              = 2*l+1

        dtauv(l,nw,ng) = dtaukv(k,nw,ng)+dtaukv(k+1,nw,ng)
        cosbv(l,nw,ng) = ( gvdst(k,nw)  * taureflk(k,nw) +           &
                         gvdst(k+1,nw)* taureflk(k+1,nw) +         &
                         gvcld(k,nw)  * taucldk(k,nw)  +           &
                         gvcld(k+1,nw)* taucldk(k+1,nw) ) /        &
                       ( tray(k,nw)     + tray(k+1,nw) +           &
                         taureflk(k,nw) + taureflk(k+1,nw) +       &
                         taucldk(k,nw)  + taucldk(k+1,nw) )

        wbarv(l,nw,ng) = ( taureflk(k,nw) + taureflk(k+1,nw) +       &
                         taucldk(k,nw)  + taucldk(k+1,nw)  +       &
                        (tray(k,nw)+tray(k+1,nw))*0.9999)/         &
                         dtauv(l,nw,ng)
    end do

    !  special bottom layer

    l              = l_nlayrad
    k              = 2*l+1
    dtauv(l,nw,ng) = dtaukv(k,nw,ng)
    cosbv(l,nw,ng) = ( gvdst(k,nw) * taureflk(k,nw) +              &
                   gvcld(k,nw) * taucldk(k,nw) ) /             &
                 ( tray(k,nw)  + taureflk(k,nw) +              &
                   taucldk(k,nw) )

    wbarv(l,nw,ng) = (taureflk(k,nw) + taucldk(k,nw) +             &
                  tray(k,nw)*0.9999)/dtauv(l,nw,ng)

    !  . . .now the other gauss points, if needed.

    do ng=1,l_ngauss-1
        if(taugsurf(nw,ng).gt.tlimits) then
            do l=1,l_layers
                k              = 2*l+1
                dtauv(l,nw,ng) = dtaukv(k,nw,ng)+dtaukv(k+1,nw,ng)
                cosbv(l,nw,ng) = cosbv(l,nw,l_ngauss)
                wbarv(l,nw,ng) = ( taureflk(k,nw) + taureflk(k+1,nw) +     &
                                   taucldk(k,nw)  + taucldk(k+1,nw) +      &
                                  (tray(k,nw)+tray(k+1,nw))*0.9999)/       &
                                   dtauv(l,nw,ng)
            end do


            !  special bottom layer

            l              = l_nlayrad
            k              = 2*l+1
            dtauv(l,nw,ng) = dtaukv(k,nw,ng)
            cosbv(l,nw,ng) = cosbv(l,nw,l_ngauss)
            wbarv(l,nw,ng) = ( taureflk(k,nw) + taucldk(k,nw) +          &
                             tray(k,nw)*0.9999 ) / dtauv(l,nw,ng)
        endif
    end do
end do     ! nw spectral loop

!     TOTAL EXTINCTION OPTICAL DEPTHS

do nw=1,l_nspectv
    ng = l_ngauss
    tauv(1,nw,ng) = 0.0d0
    do l=1,l_nlayrad
        tauv(l+1,nw,ng) = tauv(l,nw,ng)+dtauv(l,nw,ng)
    end do

    taucumv(1,nw,ng)=0.0d0
    do k=2,l_levels
        taucumv(k,nw,ng)=taucumv(k-1,nw,ng)+dtaukv(k,nw,ng)
    end do

    do ng=1,l_ngauss-1
        tauv(1,nw,ng)=0.0d0
        do l=1,l_nlayrad
            tauv(l+1,nw,ng)=tauv(l,nw,ng)+dtauv(l,nw,ng)
        end do

        taucumv(1,nw,ng)=0.0d0
        do k=2,l_levels
            taucumv(k,nw,ng)=taucumv(k-1,nw,ng)+dtaukv(k,nw,ng)
        end do
    end do
end do

!  Restore old tauref values

do k=1,l_levels+1
    tauref(k) = tauref_save(k)
    taurefcld(k) = taurefcld_save(k)
end do

return
end subroutine optcv

!=====================================================================
!=====================================================================

subroutine optcv15(dtauv,tauv,taucumv,plev,          &
               qxvdst,qsvdst,gvdst,wbarv,cosbv,        &
               tauref,tmid,pmid,taugsurf,qh2o,  &
               qextrefcld,taurefcld,qxvcld,qsvcld,gvcld, &
               qextrefco2cld,taurefco2cld,  &
               qxvco2cld,qsvco2cld,gvco2cld)
! this subroutine sets the optical constants in the visible for 15 bands
! it calcualtes for each layer, for each specral interval in the visible
! layer: wbar, dtau, cosbar
! level: tau
!
! tauv(l,nw,ng) is the cumulative optical depth at the top of radiation code
! layer l. nw is spectral wavelength interval, ng the gauss point index.
!
!     tlev(l) - temperature at the layer boundary
!     plev(l) - pressure at the layer boundary (i.e. level)
!     co2v(nt,nps,nw,ng) - visible co2 k-coefficients
!
!----------------------------------------------------------------------C

implicit none

real*8  :: dtauv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: dtaukv(l_levels+1,l_nspectv,l_ngauss)
real*8  :: tauv(l_nlevrad,l_nspectv,l_ngauss)
real*8  :: taucumv(l_levels,l_nspectv,l_ngauss)
real*8  :: plev(l_levels)
real*8  :: tmid(l_levels), pmid(l_levels), lpmid(l_levels)
real*8  :: cosbv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: wbarv(l_nlayrad,l_nspectv,l_ngauss)

!     for dust
real*8  :: qxvdst(l_levels+1,l_nspectv)
real*8  :: qsvdst(l_levels+1,l_nspectv)
real*8  :: gvdst(l_levels+1,l_nspectv)
real*8  :: qextref(l_levels+1)
real*8  :: tauref(l_levels+1)
real*8  :: tauref_save(l_levels+1)

!     for clouds
real*8  :: qxvcld(l_levels+1,l_nspectv)
real*8  :: qsvcld(l_levels+1,l_nspectv)
real*8  :: gvcld(l_levels+1,l_nspectv)
real*8  :: qextrefcld(l_levels+1)
real*8  :: taurefcld(l_levels+1)
real*8  :: taurefcld_save(l_levels+1)

real*8  :: tcloud(l_levels,l_nspectv)

real*8  :: taureflk(l_levels+1,l_nspectv)
real*8  :: taucldk(l_levels+1,l_nspectv)

!   CO2 clouds
real*8  :: qxvco2cld(l_levels+1,l_nspectv)
real*8  :: qsvco2cld(l_levels+1,l_nspectv)
real*8  :: gvco2cld(l_levels+1,l_nspectv)
real*8  :: qextrefco2cld(l_levels+1)
real*8  :: taurefco2cld(l_levels+1)
real*8  :: taurefco2cld_save(l_levels+1)

real*8  :: tco2cloud(l_levels,l_nspectv)

real*8  :: tauco2cldk(l_levels+1,l_nspectv)

integer :: l, nw, ng, k, ng1(l_nspectv), lk
integer :: mt(l_levels,2), mp(l_levels,2), mw(l_levels)
integer :: mtt(2), mpt(2), mn

real*8  :: ans, taugas
real*8  :: tray(l_levels,l_nspectv)
real*8  :: taeros(l_levels,l_nspectv)
real*8  :: dpr(l_levels), u(l_levels)

real*8  :: taugsurf(l_nspectv,l_ngauss-1), trayaer

!  Reference wavelength is (now) bin #2 - put into qextref
!      real*8 QextREF

!  Water mixing ratio stuff

real*8  :: qh2o(l_levels),  wratio(l_levels)
integer :: nh2o(l_levels)

logical :: pinter(l_levels), tinter(l_levels)

!======================================================================C

!  Save old tauref values

do k=1,l_levels+1
    tauref_save(k) = tauref(k)
    taurefcld_save(k) = taurefcld(k)
    taurefco2cld_save(k) = taurefco2cld(k)
end do

!  determine the total gas opacity throughout the column, for each
!  spectral interval, nw, and each gauss point, ng.
!  calculate the continuum opacities, i.e., those that do not depend on
!  ng, the gauss index.

do ng=1,l_ngauss-1
    do nw=1,l_nspectv
        taugsurf(nw,ng) = 0.0d0
    end do
end do

tray = 0.0
taeros = 0.0
tcloud = 0.0
dtauv = 0.0
taucumv = 0.0
dtaukv = 0.0
tco2cloud = 0.0

do k=2,l_levels
    dpr(k)   = plev(k)-plev(k-1)
    u(k)     = cmk*dpr(k)
    lpmid(k) = log10(pmid(k))

    if (use_boxinterp12) then
        call tpindex15(lpmid(k),tmid(k),qh2o(k),pfgasref,tgasref, &
                 mtt,mpt,mw(k),pinter(k),tinter(k))
    else
        call tpindex15(lpmid(k),tmid(k),qh2o(k),pgasref,tgasref, &
                 mtt,mpt,mw(k),pinter(k),tinter(k))
    endif

    do mn=1,2
        mt(k,mn) = mtt(mn)
        mp(k,mn) = mpt(mn)
    end do

    qextref(k) = qxvdst(k,l_nrefv)

    tauref(k)    = tauref(k) / qextref(k)
    taurefcld(k) = taurefcld(k) / qextrefcld(k)
    taurefco2cld(k) = taurefco2cld(k) / qextrefco2cld(k)

    do nw=1,l_nspectv
        tray(k,nw)   = tauray(nw)*dpr(k)
        taeros(k,nw) = tauref(k)    * qxvdst(k,nw)
        tcloud(k,nw) = taurefcld(k) * qxvcld(k,nw)
        tco2cloud(k,nw) = taurefco2cld(k) * qxvco2cld(k,nw)
    end do
end do

!  TRAYAER is Tau RAYleigh scattering, plus AERosol opacity

do k=2,l_levels

    do nw=1,l_nspectv

        trayaer = tray(k,nw) + taeros(k,nw) + tcloud(k,nw) + tco2cloud(k,nw)

        do ng=1,l_ngauss-1

            !           now compute taugas
            if (use_boxinterp12) then
                call boxinterp(co2v(mt(k,1),mp(k,1),mw(k),nw,ng),          &
                     co2v(mt(k,2),mp(k,1),mw(k),nw,ng),                    &
                     co2v(mt(k,1),mp(k,2),mw(k),nw,ng),                    &
                     co2v(mt(k,2),mp(k,2),mw(k),nw,ng),tgasref(mt(k,1)),   &
                     tgasref(mt(k,2)),pfgasref(mp(k,1)),pfgasref(mp(k,2)),   &
                     tmid(k),lpmid(k),tinter(k),pinter(k),ans)
            else
                call boxinterp(co2v(mt(k,1),mp(k,1),mw(k),nw,ng),          &
                     co2v(mt(k,2),mp(k,1),mw(k),nw,ng),                    &
                     co2v(mt(k,1),mp(k,2),mw(k),nw,ng),                    &
                     co2v(mt(k,2),mp(k,2),mw(k),nw,ng),tgasref(mt(k,1)),   &
                     tgasref(mt(k,2)),pgasref(mp(k,1)),pgasref(mp(k,2)),   &
                     tmid(k),lpmid(k),tinter(k),pinter(k),ans)
            endif

            taugas          = u(k)*10.0d0**ans
            taugsurf(nw,ng) = taugsurf(nw,ng) + taugas
            dtaukv(k,nw,ng) = taugas + trayaer
        end do

        !  Now fill in the "clear" part of the spectrum (NG = L_NGAUSS)
        !  Which holds continuum opacity only

        ng = l_ngauss
        dtaukv(k,nw,ng) = taeros(k,nw)+tray(k,nw)+tcloud(k,nw)+tco2cloud(k,nw)

    end do
end do

!  Now the full treatment for the layers, where besides the opacity
!  we need to calculate the scattering albedo and asymmetry factors
!  for each layer

do nw=1,l_nspectv
    do k=2,l_levels
        taureflk(k,nw) = tauref(k)    * qsvdst(k,nw)
        taucldk(k,nw)  = taurefcld(k) * qsvcld(k,nw)
        tauco2cldk(k,nw) = taurefco2cld(k) * qsvco2cld(k,nw)
    enddo
enddo

do nw=1,l_nspectv

    !  First, the special "clear" channel
    ng = l_ngauss
    do l=1,l_layers
        k              = 2*l+1
        dtauv(l,nw,ng) = dtaukv(k,nw,ng)+dtaukv(k+1,nw,ng)
        cosbv(l,nw,ng) = ( gvdst(k,nw)  * taureflk(k,nw) +          &
                         gvdst(k+1,nw)* taureflk(k+1,nw) +         &
                         gvcld(k,nw)  * taucldk(k,nw)  +           &
                         gvcld(k+1,nw)* taucldk(k+1,nw) +          &
                         gvco2cld(k,nw)  * tauco2cldk(k,nw)  +     &
                         gvco2cld(k+1,nw)* tauco2cldk(k+1,nw) ) /  &
                       ( tray(k,nw)     + tray(k+1,nw) +           &
                         taureflk(k,nw) + taureflk(k+1,nw) +       &
                         taucldk(k,nw)  + taucldk(k+1,nw)  +       &
                         tauco2cldk(k,nw)  + tauco2cldk(k+1,nw) )

        wbarv(l,nw,ng) = ( taureflk(k,nw) + taureflk(k+1,nw) +     &
                         taucldk(k,nw)  + taucldk(k+1,nw)  +       &
                         tauco2cldk(k,nw)  + tauco2cldk(k+1,nw)  +       &
                        (tray(k,nw)+tray(k+1,nw))*0.9999)/         &
                         dtauv(l,nw,ng)
    end do

    !  special bottom layer

    l              = l_nlayrad
    k              = 2*l+1
    dtauv(l,nw,ng) = dtaukv(k,nw,ng)
    cosbv(l,nw,ng) = ( gvdst(k,nw) * taureflk(k,nw) +              &
                   gvcld(k,nw) * taucldk(k,nw) +               &
                   gvco2cld(k,nw) * tauco2cldk(k,nw)  ) /             &
                 ( tray(k,nw)  + taureflk(k,nw) +              &
                   taucldk(k,nw) + tauco2cldk(k,nw) )

    wbarv(l,nw,ng) = (taureflk(k,nw) + taucldk(k,nw) + tauco2cldk(k,nw) +             &
                  tray(k,nw)*0.9999)/dtauv(l,nw,ng)

    !  . . .now the other gauss points, if needed.

    do ng=1,l_ngauss-1
        do l=1,l_layers
            k              = 2*l+1
            dtauv(l,nw,ng) = dtaukv(k,nw,ng)+dtaukv(k+1,nw,ng)
            cosbv(l,nw,ng) = cosbv(l,nw,l_ngauss)
            wbarv(l,nw,ng) = ( taureflk(k,nw) + taureflk(k+1,nw) +     &
                               taucldk(k,nw)  + taucldk(k+1,nw) +      &
                               tauco2cldk(k,nw)  + tauco2cldk(k+1,nw) +      &
                              (tray(k,nw)+tray(k+1,nw))*0.9999)/       &
                               dtauv(l,nw,ng)
        end do

        !  special bottom layer

        l              = l_nlayrad
        k              = 2*l+1
        dtauv(l,nw,ng) = dtaukv(k,nw,ng)
        cosbv(l,nw,ng) = cosbv(l,nw,l_ngauss)
        wbarv(l,nw,ng) = ( taureflk(k,nw) + taucldk(k,nw) +         &
                           tauco2cldk(k,nw)  +                      &
                         tray(k,nw)*0.9999 ) / dtauv(l,nw,ng)
    end do

end do     ! nw spectral loop

!     TOTAL EXTINCTION OPTICAL DEPTHS

do nw=1,l_nspectv
    ng = l_ngauss
    tauv(1,nw,ng) = 0.0d0
    do l=1,l_nlayrad
        tauv(l+1,nw,ng) = tauv(l,nw,ng)+dtauv(l,nw,ng)
    end do

    taucumv(1,nw,ng)=0.0d0
    do k=2,l_levels
        taucumv(k,nw,ng)=taucumv(k-1,nw,ng)+dtaukv(k,nw,ng)
    end do

    do ng=1,l_ngauss-1
        tauv(1,nw,ng)=0.0d0
        do l=1,l_nlayrad
            tauv(l+1,nw,ng)=tauv(l,nw,ng)+dtauv(l,nw,ng)
        end do

        taucumv(1,nw,ng)=0.0d0
        do k=2,l_levels
            taucumv(k,nw,ng)=taucumv(k-1,nw,ng)+dtaukv(k,nw,ng)
        end do
    end do
end do


!  Restore old tauref values

do k=1,l_levels+1
    tauref(k) = tauref_save(k)
    taurefcld(k) = taurefcld_save(k)
    taurefco2cld(k) = taurefco2cld_save(k)
end do

return
end subroutine optcv15


!=====================================================================
!=====================================================================

subroutine optci(dtaui,taucumi,plev,       &
               qrefv,qxidst,qsidst,gidst,cosbi,wbari,tauref,   &
               tmid,pmid,taugsurf,qh2o,                &
               qextrefcld,taurefcld,qxicld,qsicld,gicld)
! THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE INFRARED
! IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE IR
! LAYER: WBAR, DTAU, COSBAR
! LEVEL: TAU
!
! Qrefv is the extinction coefficient at the reference (visible)
! wavelength - 0.67 microns.
!
! TAUI(L,LW) is the cumulative optical depth at level L (or alternatively
! at the *bottom* of layer L), LW is the spectral wavelength interval.
!
!     TLEV(L) - Temperature at the layer boundary (i.e. level)
!     PLEV(L) - Pressure at the layer boundary (i.e. level)
!     CO2_KI(NT,NP,NW,NG) - IR CO2 k-coefficients
!                           CO2_K(temp,Pres,Waveln,gauss)
!                           currently: CO2_K(7,11,5,17)
!
!----------------------------------------------------------------------C


implicit none


real*8  :: dtaui(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: dtauki(l_levels+1,l_nspecti,l_ngauss)
real*8  :: taui(l_nlevrad,l_nspecti,l_ngauss)
real*8  :: taucumi(l_levels,l_nspecti,l_ngauss)
real*8  :: taugas
real*8  :: plev(l_levels)
real*8  :: tmid(l_levels), pmid(l_levels)



real*8  :: cosbi(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: wbari(l_nlayrad,l_nspecti,l_ngauss)

integer :: l, nw, ng, k, lk
integer :: mt(l_levels), mp(l_levels), np(l_levels)
real*8  :: ans, taurefl
real*8  :: taeros(l_levels,l_nspecti)
real*8  :: dpr(l_levels), u(l_levels), tauac
real*8  :: lcoef(4), lkcoef(l_levels,4)

!     for dust
real*8  :: qxidst(l_levels+1,l_nspecti)
real*8  :: qsidst(l_levels+1,l_nspecti)
real*8  :: gidst(l_levels+1,l_nspecti)
real*8  :: qrefv(l_levels+1)
real*8  :: tauref(l_levels+1)
real*8  :: tauref_save(l_levels+1)

!     for clouds
real*8  :: qxicld(l_levels+1,l_nspecti)
real*8  :: qsicld(l_levels+1,l_nspecti)
real*8  :: gicld(l_levels+1,l_nspecti)
real*8  :: qextrefcld(l_levels+1)
real*8  :: taurefcld(l_levels+1)
real*8  :: taurefcld_save(l_levels+1)

!!!   rjw
!   for diagnostic brightness temperature calculations: dust
!   Note that these are dimensioned by l_layers,  not levels.



real*8  :: tcloud(l_levels,l_nspecti),taureflcld

real*8  :: taureflk(l_levels+1,l_nspecti)
real*8  :: taucldk( l_levels+1,l_nspecti)

! fraction of zeros in each spectral interval, as a function of t, p

real*8  :: dt, tt
real*8  :: taugsurf(l_nspecti,l_ngauss-1)

!  water mixing ratio variables

real*8  :: qh2o(l_levels), wratio(l_levels)
real*8  :: kcoef(4)
integer :: nh2o(l_levels)

logical  :: pinter(L_LEVELS), tinter(L_LEVELS)
real*8  :: kcoe
!======================================================================C

dtauki(l_levels+1,:,:)   = 0.0d0
taureflk(l_levels+1,:) = 0.0d0
taucldk(l_levels+1,:)  = 0.0d0

dpr = 0.0
lkcoef = 0.0
taeros = 0.0
tcloud = 0.0
u = 0.0
dtauki = 0.0
taureflk = 0.0
taucldk = 0.0
taucumi = 0.0
dtaui = 0.0

!  Save old tauref values

do k=1,L_LEVELS+1
    tauref_save(k) = tauref(k)
    taurefcld_save(k) = taurefcld(k)
end do

!  Determine the total gas opacity throughout the column, for each
!  spectral interval, NW, and each Gauss point, NG.

taugsurf(:,:l_ngauss-1) = 0.0d0

do k=2,l_levels
    dpr(k) = plev(k)-plev(k-1)
    u(k)   = cmk*dpr(k)

    call tpindex(pmid(k),tmid(k),qh2o(k),pfgasref,tgasref ,&
             lcoef,mt(k),mp(k),nh2o(k),wratio(k))

    do lk=1,4
        lkcoef(k,lk) = lcoef(lk)
    end do

    tauref(k)    = tauref(k)    / qrefv(k)
    taurefcld(k) = taurefcld(k) / qextrefcld(k)

    do nw=1,l_nspecti
        taeros(k,nw) = tauref(k)    * qxidst(k,nw)
        tcloud(k,nw) = taurefcld(k) * qxicld(k,nw)
    end do
end do

do k=2,l_levels
    do nw=1,l_nspecti
        do ng=1,l_ngauss-1

            !           now compute taugas

            !  interpolate between water mixing ratios
            !  wratio = 0.0 if the requested water amount is equal to, or outside the
            !  the range of water amount data.

            kcoef(1) = co2i(mt(k),mp(k),nh2o(k),nw,ng) + wratio(k)*    &
                  (co2i(mt(k),mp(k),nh2o(k)+1,nw,ng) -             &
                   co2i(mt(k),mp(k),nh2o(k),nw,ng))

            kcoef(2) = co2i(mt(k),mp(k)+1,nh2o(k),nw,ng) + wratio(k)*  &
                  (co2i(mt(k),mp(k)+1,nh2o(k)+1,nw,ng) -           &
                   co2i(mt(k),mp(k)+1,nh2o(k),nw,ng))

            kcoef(3) = co2i(mt(k)+1,mp(k)+1,nh2o(k),nw,ng) + wratio(k)*&
                  (co2i(mt(k)+1,mp(k)+1,nh2o(k)+1,nw,ng) -         &
                   co2i(mt(k)+1,mp(k)+1,nh2o(k),nw,ng))

            kcoef(4) = co2i(mt(k)+1,mp(k),nh2o(k),nw,ng) + wratio(k)*  &
                  (co2i(mt(k)+1,mp(k),nh2o(k)+1,nw,ng) -           &
                   co2i(mt(k)+1,mp(k),nh2o(k),nw,ng))


            !  interpolate the co2 k-coefficients to the requested t,p


            ans = lkcoef(k,1)*kcoef(1) + lkcoef(k,2)*kcoef(2) +        &
              lkcoef(k,3)*kcoef(3) + lkcoef(k,4)*kcoef(4)

            taugas          = u(k)*ans
            taugsurf(nw,ng) = taugsurf(nw,ng) + taugas
            dtauki(k,nw,ng) = taugas+taeros(k,nw)+tcloud(k,nw)
        end do

        !  now fill in the "clear" part of the spectrum (ng = l_ngauss)
        !  which holds continuum opacity only

        ng              = l_ngauss
        dtauki(k,nw,ng) = taeros(k,nw)+tcloud(k,nw)
    end do
end do

!  Now the full treatment for the layers, where besides the opacity
!  we need to calculate the scattering albedo and asymmetry factors
!  for each layer

do nw=1,l_nspecti
    do k=2,l_levels+1
        taureflk(k,nw) = tauref(k)    * qsidst(k,nw)
        taucldk(k,nw)  = taurefcld(k) * qsicld(k,nw)
    enddo
enddo

do nw=1,l_nspecti

    !  first, the special "clear" channel

    ng = l_ngauss

    do l=1,l_nlayrad
        k              = 2*l+1
        dtaui(l,nw,ng) = dtauki(k,nw,ng) + dtauki(k+1,nw,ng) + 1.d-50
        if(dtaui(l,nw,ng) .gt. 1.0e-9) then
            wbari(l,nw,ng) = (taureflk(k,nw)+ taureflk(k+1,nw) +       &
                              taucldk(k,nw) + taucldk(k+1,nw)   ) /    &
                              dtaui(l,nw,ng)
        else
            wbari(l,nw,ng) = 0.0d0
            dtaui(l,nw,ng) = 1.0e-9
        endif

        tauac = taureflk(k,nw)+ taureflk(k+1,nw) + taucldk(k,nw) +   &
              taucldk(k+1,nw)

        if(tauac .gt. 0.0) then
            cosbi(l,nw,ng) = ( gidst(k,nw)   * taureflk(k,nw) +        &
                               gidst(k+1,nw) * taureflk(k+1,nw) +      &
                               gicld(k,nw)   * taucldk(k,nw) +         &
                               gicld(k+1,nw) * taucldk(k+1,nw) ) /     &
                              (taureflk(k,nw)+ taureflk(k+1,nw) +      &
                               taucldk(k,nw) + taucldk(k+1,nw)   )
        else
            cosbi(l,nw,ng) = 0.0d0
        end if
    end do

    !  . . .now the other gauss points, if needed.

    do ng=1,l_ngauss-1
        if(taugsurf(nw,ng) .gt. tlimiti) then  !like john's version
            do l=1,l_nlayrad
                k              = 2*l+1
                dtaui(l,nw,ng) = dtauki(k,nw,ng)+dtauki(k+1,nw,ng)+1.d-50
                if(dtaui(l,nw,ng) .gt. 1.0e-9) then
                    wbari(l,nw,ng) = (taureflk(k,nw)+ taureflk(k+1,nw) +     &
                                    taucldk(k,nw) + taucldk(k+1,nw)   ) /  &
                                    dtaui(l,nw,ng)
                else
                    wbari(l,nw,ng) = 0.0d0
                    dtaui(l,nw,ng) = 1.0e-9
                endif

                cosbi(l,nw,ng) = cosbi(l,nw,l_ngauss)
            end do
        endif
    end do

end do     ! nw spectral loop

!     TOTAL EXTINCTION OPTICAL DEPTHS

do nw=1,l_nspecti
    ng = l_ngauss
    taui(1,nw,ng) = 0.0d0
    do l=1,l_nlayrad
        taui(l+1,nw,ng) = taui(l,nw,ng)+dtaui(l,nw,ng)
    end do

    taucumi(1,nw,ng)=0.0d0
    do k=2,l_levels
        taucumi(k,nw,ng)=taucumi(k-1,nw,ng)+dtauki(k,nw,ng)
    end do

    do ng=1,l_ngauss-1
        !  with no check, the low values have non-zero numbers order 1.e-3
        !  with the check, the new code sets to 0. this matches john's version
        !  which is better? urata 03/24/2018
        taui(1,nw,ng)=0.0d0
        do l=1,l_nlayrad
            taui(l+1,nw,ng)=taui(l,nw,ng)+dtaui(l,nw,ng)
        end do

        taucumi(1,nw,ng)=0.0d0
        do k=2,l_levels
            taucumi(k,nw,ng)=taucumi(k-1,nw,ng)+dtauki(k,nw,ng)
        end do
    end do
end do

!  Restore old tauref values

do k=1,l_levels+1
    tauref(k) = tauref_save(k)
    taurefcld(k) = taurefcld_save(k)
end do

return
end subroutine optci


!=====================================================================
!=====================================================================

subroutine optci15(dtaui,taucumi,plev,tlev,       &
               qrefv,qxidst,qsidst,gidst,cosbi,wbari,tauref,   &
               tmid,pmid,taugsurf,qh2o,                       &
               qextrefcld,taurefcld,qxicld,qsicld,gicld,      &
               qextrefco2cld,taurefco2cld,                    &
               qxico2cld,qsico2cld,gico2cld)
!
! THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE INFRARED
! IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE IR
! LAYER: WBAR, DTAU, COSBAR
! LEVEL: TAU
!
! Qrefv is the extinction coefficient at the reference (visible)
! wavelength - 0.67 microns.
!
! TAUI(L,LW) is the cumulative optical depth at level L (or alternatively
! at the *bottom* of layer L), LW is the spectral wavelength interval.
!
!     TLEV(L) - Temperature at the layer boundary (i.e. level)
!     PLEV(L) - Pressure at the layer boundary (i.e. level)
!     CO2_KI(NT,NP,NW,NG) - IR CO2 k-coefficients
!                           CO2_K(temp,Pres,Waveln,gauss)
!                           currently: CO2_K(7,11,5,17)
!
!----------------------------------------------------------------------C


implicit none

real*8  :: dtaui(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: dtauki(l_levels+1,l_nspecti,l_ngauss)
real*8  :: taui(l_nlevrad,l_nspecti,l_ngauss)
real*8  :: taucumi(l_levels,l_nspecti,l_ngauss)
real*8  :: taugas
real*8  :: plev(l_levels)
real*8  :: tlev(l_levels)
real*8  :: tmid(l_levels), pmid(l_levels), lpmid(l_levels)

real*8  :: cosbi(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: wbari(l_nlayrad,l_nspecti,l_ngauss)

integer :: l, nw, ng, k, lk
integer :: mt(l_levels,2), mp(l_levels,2), mw(l_levels)
integer :: mtt(2), mpt(2), mn

real*8  :: ans, taurefl
real*8  :: taeros(l_levels,l_nspecti)
real*8  :: dpr(l_levels), u(l_levels), tauac

!     for dust
real*8  :: qxidst(l_levels+1,l_nspecti)
real*8  :: qsidst(l_levels+1,l_nspecti)
real*8  :: gidst(l_levels+1,l_nspecti)
real*8  :: qrefv(l_levels+1)
real*8  :: tauref(l_levels+1)
real*8  :: tauref_save(l_levels+1)

!     for clouds
real*8  :: qxicld(l_levels+1,l_nspecti)
real*8  :: qsicld(l_levels+1,l_nspecti)
real*8  :: gicld(l_levels+1,l_nspecti)
real*8  :: qextrefcld(l_levels+1)
real*8  :: taurefcld(l_levels+1)
real*8  :: taurefcld_save(l_levels+1)

real*8  :: tcloud(l_levels,l_nspecti),taureflcld

real*8  :: taureflk(l_levels+1,l_nspecti)
real*8  :: taucldk(l_levels+1,l_nspecti)

!     for CO2 clouds
real*8  :: qxico2cld(l_levels+1,l_nspecti)
real*8  :: qsico2cld(l_levels+1,l_nspecti)
real*8  :: gico2cld(l_levels+1,l_nspecti)
real*8  :: qextrefco2cld(l_levels+1)
real*8  :: taurefco2cld(l_levels+1)
real*8  :: taurefco2cld_save(l_levels+1)

real*8  :: tco2cloud(l_levels,l_nspecti)

real*8  :: tauco2cldk(l_levels+1,l_nspecti)

! fraction of zeros in each spectral interval, as a function of T, P

real*8  :: dt, tt
real*8  :: taugsurf(l_nspecti,l_ngauss-1)

!  water mixing ratio variables

real*8  :: qh2o(l_levels), wratio(l_levels)

logical  :: pinter(l_levels), tinter(l_levels)
real*8  :: kcoe

! For CO2 and H2 CIA calculations (CO2-CO2 CIA based on Wordsworth et al., 2010; CO2-H2 CIA based on Turbet et al., 2020)

real*8  :: wnoi(L_NSPECTI)
real*8  :: dwni(L_NSPECTI)
real*8  :: dtaucia(L_LEVELS+1,L_NSPECTI)
real*8  :: dtau_co2cia(L_LEVELS+1,L_NSPECTI)
real*8  :: dtau_h2n2cia(L_LEVELS+1,L_NSPECTI)
real*8  :: taucia(L_NSPECTI)
real*8  :: kgbar,kbbar,khbar
real*8  :: tmp,pmp
real*8  :: plength,namgCO2,namgH2
real*8  :: kcoeff(l_levels,l_nspecti,l_ngauss)
integer :: j

!======================================================================C

dtauki(l_levels+1,:,:)   = 0.0d0
taureflk(l_levels+1,:) = 0.0d0
taucldk(l_levels+1,:)  = 0.0d0
tauco2cldk(l_levels+1,:)  = 0.0d0

dpr = 0.0

taeros = 0.0
tcloud = 0.0
u = 0.0
dtauki = 0.0
taureflk = 0.0
taucldk = 0.0
tauco2cldk = 0.0
taucumi = 0.0
dtaui = 0.0
dtaucia = 0.0d0


!  save old tauref values

do k=1,l_levels+1
    tauref_save(k) = tauref(k)
    taurefcld_save(k) = taurefcld(k)
    taurefco2cld_save(k) = taurefco2cld(k)
end do

!  Determine the total gas opacity throughout the column, for each
!  spectral interval, NW, and each Gauss point, NG.

taugsurf(:,:l_ngauss-1) = 0.0d0

do k=2,l_levels
    dpr(k)   = plev(k)-plev(k-1)
    u(k)     = cmk*dpr(k)
    lpmid(k) = log10(pmid(k)) !log10(pmid(k))

    if (use_boxinterp12) then
        call tpindex15(lpmid(k),tmid(k),qh2o(k),pfgasref,tgasref, &
                 mtt,mpt,mw(k),pinter(k),tinter(k))
    else
        call tpindex15(lpmid(k),tmid(k),qh2o(k),pgasref,tgasref, &
                 mtt,mpt,mw(k),pinter(k),tinter(k))
    endif

    do mn=1,2
        mt(k,mn) = mtt(mn)
        mp(k,mn) = mpt(mn)
    end do

    tauref(k)    = tauref(k)    / qrefv(k)
    taurefcld(k) = taurefcld(k) / qextrefcld(k)
    taurefco2cld(k) = taurefco2cld(k) / qextrefco2cld(k)

    do nw=1,l_nspecti
        taeros(k,nw) = tauref(k)    * qxidst(k,nw)
        tcloud(k,nw) = taurefcld(k) * qxicld(k,nw)
        tco2cloud(k,nw) = taurefco2cld(k) * qxico2cld(k,nw)

        !Begin CIA calculation
        if (do_cia) then

            tmp = .5*(tlev(k)+tlev(k-1))
            pmp = .5*(plev(k)+plev(k-1))

            call kinter(nw,tmp,kgbar_tab,kgbar)
            call kinter(nw,tmp,kbbar_tab,kbbar)
            call kinter(nw,tmp,khbar_tab,khbar)

            dtau_co2cia(k,nw) = 0.
            dtau_h2n2cia(k,nw) = 0.
            namgCO2=0.
            namgH2=0.

            if(k.gt.3) then
                plength=50.8*tmp*(dpr(k)/pmp) !path length in meters

                ! define amounts of CO2 and H2 in amagats where cia_co2 and cia_h2 are the molar concentration of CO2 and H2 (e.g. mols H2 / mols of both CO2+H2) 
                namgCO2=cia_co2*(pmp/tmp)*(273.15/1013.25)
                namgH2=cia_h2*(pmp/tmp)*(273.15/1013.25)

                dtau_co2cia(k,nw)=namgCO2*namgCO2*100.*(kgbar+kbbar)*plength
                dtau_h2n2cia(k,nw)=namgH2*namgCO2*100.*khbar*plength

            end if

            dtaucia(k,nw) = dtau_co2cia(k,nw)+dtau_h2n2cia(k,nw)
        else
            dtaucia(k,nw) = 0.0
        endif  !End CIA tau calculation
        

    end do
end do

do k=2,l_levels
    do nw=1,l_nspecti
        do ng=1,l_ngauss-1

            !           now compute taugas
            if (use_boxinterp12) then
                call boxinterp(co2i(mt(k,1),mp(k,1),mw(k),nw,ng),          &
                     co2i(mt(k,2),mp(k,1),mw(k),nw,ng),                    &
                     co2i(mt(k,1),mp(k,2),mw(k),nw,ng),                    &
                     co2i(mt(k,2),mp(k,2),mw(k),nw,ng),tgasref(mt(k,1)),   &
                     tgasref(mt(k,2)),pfgasref(mp(k,1)),pfgasref(mp(k,2)),   &
                     tmid(k),lpmid(k),tinter(k),pinter(k),ans)
            else
                call boxinterp(co2i(mt(k,1),mp(k,1),mw(k),nw,ng),          &
                     co2i(mt(k,2),mp(k,1),mw(k),nw,ng),                    &
                     co2i(mt(k,1),mp(k,2),mw(k),nw,ng),                    &
                     co2i(mt(k,2),mp(k,2),mw(k),nw,ng),tgasref(mt(k,1)),   &
                     tgasref(mt(k,2)),pgasref(mp(k,1)),pgasref(mp(k,2)),   &
                     tmid(k),lpmid(k),tinter(k),pinter(k),ans)
            endif

            taugas          = u(k)*10.0d0**ans

            taugsurf(nw,ng) = taugsurf(nw,ng) + taugas
            dtauki(k,nw,ng) = taugas+taeros(k,nw)+tcloud(k,nw)+  &  
                              tco2cloud(k,nw)+dtaucia(k,nw)
        end do

        !  now fill in the "clear" part of the spectrum (ng = l_ngauss)
        !  which holds continuum opacity only

        ng              = l_ngauss
        dtauki(k,nw,ng) = taeros(k,nw)+tcloud(k,nw)+tco2cloud(k,nw)+dtaucia(k,nw)
    end do

end do

!  Now the full treatment for the layers, where besides the opacity
!  we need to calculate the scattering albedo and asymmetry factors
!  for each layer

do nw=1,l_nspecti
    do k=2,l_levels+1
        taureflk(k,nw) = tauref(k)    * qsidst(k,nw)
        taucldk(k,nw)  = taurefcld(k) * qsicld(k,nw)
        tauco2cldk(k,nw)  = taurefco2cld(k) * qsico2cld(k,nw)
    enddo
enddo

do nw=1,l_nspecti

!  First, the special "clear" channel

    ng = l_ngauss

    do l=1,l_nlayrad
        k              = 2*l+1
        dtaui(l,nw,ng) = dtauki(k,nw,ng) + dtauki(k+1,nw,ng) + 1.d-50
        if(dtaui(l,nw,ng) .gt. 1.0e-9) then
            wbari(l,nw,ng) = (taureflk(k,nw)+ taureflk(k+1,nw) +       &
                              taucldk(k,nw) + taucldk(k+1,nw)  +       &
                              tauco2cldk(k,nw) + tauco2cldk(k+1,nw)     ) /    &
                              dtaui(l,nw,ng)
        else
            wbari(l,nw,ng) = 0.0d0
            dtaui(l,nw,ng) = 1.0d-9
        endif

        tauac = taureflk(k,nw)+ taureflk(k+1,nw) + taucldk(k,nw) +   &
              taucldk(k+1,nw) + tauco2cldk(k,nw) + tauco2cldk(k+1,nw) 

        if(tauac .gt. 0.0) then
            cosbi(l,nw,ng) = ( gidst(k,nw)   * taureflk(k,nw) +        &
                               gidst(k+1,nw) * taureflk(k+1,nw) +      &
                               gicld(k,nw)   * taucldk(k,nw) +         &
                               gicld(k+1,nw) * taucldk(k+1,nw)  +      &
                               gico2cld(k,nw)   * tauco2cldk(k,nw) +         &  
                               gico2cld(k+1,nw) * tauco2cldk(k+1,nw) ) /     &
                              (taureflk(k,nw)+ taureflk(k+1,nw) +      &
                               taucldk(k,nw) + taucldk(k+1,nw)  +      &
                               tauco2cldk(k,nw) + tauco2cldk(k+1,nw)   )
        else
            cosbi(l,nw,ng) = 0.0d0
        end if
    end do

    !  . . .now the other gauss points, if needed.

    do ng=1,l_ngauss-1

        do l=1,l_nlayrad
            k              = 2*l+1
            dtaui(l,nw,ng) = dtauki(k,nw,ng)+dtauki(k+1,nw,ng)+1.d-50
            if(dtaui(l,nw,ng) .gt. 1.0e-9) then
                wbari(l,nw,ng) = (taureflk(k,nw)+ taureflk(k+1,nw) +     &
                                taucldk(k,nw) + taucldk(k+1,nw)    +     &
                                tauco2cldk(k,nw) + tauco2cldk(k+1,nw) ) /  &
                                dtaui(l,nw,ng)
            else
                wbari(l,nw,ng) = 0.0d0
                dtaui(l,nw,ng) = 1.0d-9
            endif

            cosbi(l,nw,ng) = cosbi(l,nw,l_ngauss)
        end do
    end do

end do     ! nw spectral loop

!     total extinction optical depths

do nw=1,l_nspecti
    ng = l_ngauss
    taui(1,nw,ng) = 0.0d0
    do l=1,l_nlayrad
        taui(l+1,nw,ng) = taui(l,nw,ng)+dtaui(l,nw,ng)
    end do

    taucumi(1,nw,ng)=0.0d0
    do k=2,l_levels
        taucumi(k,nw,ng)=taucumi(k-1,nw,ng)+dtauki(k,nw,ng)
    end do

    do ng=1,l_ngauss-1
        taui(1,nw,ng)=0.0d0
        do l=1,l_nlayrad
            taui(l+1,nw,ng)=taui(l,nw,ng)+dtaui(l,nw,ng)
        end do

        taucumi(1,nw,ng)=0.0d0
        do k=2,l_levels
            taucumi(k,nw,ng)=taucumi(k-1,nw,ng)+dtauki(k,nw,ng)
        end do
    end do
end do

!  Restore old tauref values

do k=1,l_levels+1
    tauref(k) = tauref_save(k)
    taurefcld(k) = taurefcld_save(k)
    taurefco2cld(k) = taurefco2cld_save(k)
end do

return
end subroutine optci15


                                     
!=====================================================================
!=====================================================================

subroutine kinter(nw,temp,kbar_tab,kbar)

! This routine interpolates the CIA coefficients to the current
! atmospheric temperature
! -------------------------------------------------------------
! kbar_tab is a 2D array whose rows correspond to wavenumber
! and whose columns correspond to temperature
! -------------------------------------------------------------

    implicit none

    integer, parameter :: nwave=8
    integer, parameter :: ntemp=15
    real*8  kbar_tab(nwave,ntemp)
    real*8 kbar
    integer nw,nt

    real*8 t_tab(ntemp)
    real*8 temp
    integer kt
    real*8 logk1,logk2,logk

    data t_tab/100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,&
                       600.,650.,700.,750.,800./

! limit the bounds of temperature if we go beyond the table
    if(temp.lt.100.) temp=100.
    if(temp.gt.800.) temp=800.

    do nt=1,ntemp-1
        if(temp.ge.t_tab(nt).and.temp.le.t_tab(nt+1)) go to 100
    end do
100 kt=nt

    logk1=log(kbar_tab(nw,kt))
    logk2=log(kbar_tab(nw,kt+1))
    logk=logk1+(temp-t_tab(kt))*(logk2-logk1)/(t_tab(kt+1)-t_tab(kt))
    kbar = exp(logk)
    if(kbar_tab(nw,kt).eq.0.and.kbar_tab(nw,kt+1).eq.0.) kbar=0.

    return
end subroutine kinter

!=====================================================================
!=====================================================================

subroutine sfluxv(dtauv,tauv,taucumv,rsfv,wbarv,cosbv,      &
                ubar0,sol,nfluxtopv,fmnetv,            &
                fluxupv,fluxdnv,diffvt,taugsurf,        &
                detau1,diag)
!  calculate visible fluxes

implicit none

real*8  :: fmnetv(l_nlayrad)
real*8  :: taucumv(l_levels,l_nspectv,l_ngauss)
real*8  :: tauv(l_nlevrad,l_nspectv,l_ngauss)
real*8  :: dtauv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: fmupv(l_nlayrad), fmdv(l_nlayrad)
real*8  :: cosbv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: wbarv(l_nlayrad,l_nspectv,l_ngauss)
real*8  :: sol(l_nspectv)
real*8  :: fluxupv(l_nlayrad), fluxdnv(l_nlayrad)
real*8  :: nfluxtopv, fluxup, fluxdn

!  delta-eddington surface tau

real*8  :: detau1(l_nspectv,l_ngauss)

integer l, ng, nw, ng1
real*8  rsfv, ubar0, f0pi, btop, bsurf, taumax, eterm

real*8 diffv, diffvt

real*8 taugsurf(l_nspectv,l_ngauss-1), fzero

logical :: diag, diag2

!======================================================================C

taumax = l_taumax

diag2 = ubar0.ge.0.98

nfluxtopv = 0.0

do l=1,l_nlayrad
    fmnetv(l)  = 0.0
    fluxupv(l) = 0.0
    fluxdnv(l) = 0.0
end do

diffvt = 0.0

!     we now enter a major loop over spectral intervals in the visible
!     to calculate the net flux in each spectral interval

do nw=1,l_nspectv

    f0pi = sol(nw)

    fzero = fzerov(nw)
    if (.not.(fzero.ge.0.99)) then
        do ng=1,l_ngauss-1

            call getdetau(dtauv(1,nw,ng),tauv(1,nw,ng),                &
                        taucumv(1,nw,ng),wbarv(1,nw,ng),             &
                        cosbv(1,nw,ng),ubar0,detau1(nw,ng))


            if(taugsurf(nw,ng) .lt. tlimits) then
                fzero = fzero + (1.0-fzerov(nw))*gweight(ng)
                cycle
            end if

            !         set up the upper and lower boundary conditions on the visible

            btop = 0.0

            !         loop over the nterms beginning here

            !  detau(j1d,i1d,nw,ng) is the scaled optical depth at the surface

            eterm = min(detau1(nw,ng)/ubar0,maxexp)
            bsurf = rsfv*ubar0*sol(nw)*exp(-eterm)


            !         we can now solve for the coefficients of the two stream
            !         call a subroutine that solves  for the flux terms
            !         within each interval at the midpoint wavenumber
            !
            !         fuw and fdw are working flux arrays that will be used to
            !         return fluxes for a given nt

            call gfluxv(dtauv(1,nw,ng),tauv(1,nw,ng),taucumv(1,nw,ng),   &
                      wbarv(1,nw,ng),cosbv(1,nw,ng),ubar0,f0pi,rsfv,   &
                      btop,bsurf,fmupv,fmdv,diffv,fluxup,fluxdn,       &
                      detau1(nw,ng))

            !         now calculate the cumulative visible net flux

            nfluxtopv = nfluxtopv+(fluxup-fluxdn)*gweight(ng)*           &
                                (1.0-fzerov(nw))
            do l=1,l_nlayrad
                fmnetv(l)=fmnetv(l)+( fmupv(l)-fmdv(l) )*                  &
                                     gweight(ng)*(1.0-fzerov(nw))
                fluxupv(l) = fluxupv(l) + fmupv(l)*gweight(ng)*            &
                             (1.0-fzerov(nw))
                fluxdnv(l) = fluxdnv(l) + fmdv(l)*gweight(ng)*             &
                             (1.0-fzerov(nw))
            end do

            !         the diffuse component of the downward solar flux

            diffvt = diffvt + diffv*gweight(ng)*(1.0-fzerov(nw))

        end do   ! the gauss loop

    endif
    !       special 17th gauss point

    ng = l_ngauss

    !       set up the upper and lower boundary conditions on the visible

    btop = 0.0

    call getdetau(dtauv(1,nw,ng),tauv(1,nw,ng),                &
                  taucumv(1,nw,ng),wbarv(1,nw,ng),             &
                  cosbv(1,nw,ng),ubar0,detau1(nw,ng))

    !       loop over the nterms beginning here

    eterm = min(detau1(nw,ng)/ubar0,maxexp)
    bsurf = rsfv*ubar0*sol(nw)*exp(-eterm)

    !       we can now solve for the coefficients of the two stream
    !       call a subroutine that solves  for the flux terms
    !       within each interval at the midpoint wavenumber
    !
    !       fuw and fdw are working flux arrays that will be used to
    !       return fluxes for a given nt

    call gfluxv(dtauv(1,nw,ng),tauv(1,nw,ng),taucumv(1,nw,ng),     &
                wbarv(1,nw,ng),cosbv(1,nw,ng),ubar0,f0pi,rsfv,     &
                btop,bsurf,fmupv,fmdv,diffv,fluxup,fluxdn,         &
                detau1(nw,ng))

    !       now calculate the cumulative visible net flux

    nfluxtopv = nfluxtopv+(fluxup-fluxdn)*fzero
    do l=1,l_nlayrad
        fmnetv(l)=fmnetv(l)+( fmupv(l)-fmdv(l) )*fzero
        fluxupv(l) = fluxupv(l) + fmupv(l)*fzero
        fluxdnv(l) = fluxdnv(l) + fmdv(l)*fzero
    end do


    !       the diffuse component of the downward solar flux

    diffvt = diffvt + diffv*fzero

enddo
!     *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE VISIBLE*****

return
end subroutine sfluxv



!=====================================================================
!=====================================================================

subroutine sfluxi(plev,tlev,dtaui,taucumi,rsfi,  &
                cosbi,wbari,nfluxtopi,fmneti,          &
                fluxupi,fluxdni,taugsurf,nfluxtopis, &
                  fmnetis,fluxupis,fluxdnis,diag)
!  calculate IR fluxes

implicit none

integer :: nlevrad, l, nw, ng, nts, ntt

real*8  :: tlev(l_levels), plev(l_levels)
real*8  :: taucumi(l_levels,l_nspecti,l_ngauss)
real*8  :: fmneti(l_nlayrad),fmnetis(l_nspecti,l_nlayrad)
real*8  :: dtaui(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: fmupi(l_nlayrad), fmdi(l_nlayrad)
real*8  :: cosbi(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: wbari(l_nlayrad,l_nspecti,l_ngauss)
real*8  :: nfluxtopi, nfluxtopis(l_nspecti)
real*8  :: ftopup

real*8  :: rsfi, tsurf, bsurf, ttop, btop, tautop
real*8  :: planck, pltop
real*8  :: fluxupi(l_nlayrad), fluxdni(l_nlayrad)
real*8  :: taugsurf(l_nspecti,l_ngauss-1), fzero
real*8  :: fluxupis(l_nspecti,l_nlayrad)
real*8  :: fluxdnis(l_nspecti,l_nlayrad)

real*8  :: bsurftot, stfour

logical :: diag

!======================================================================C

nlevrad = l_nlevrad

nfluxtopi = 0.0

do l=1,l_nlayrad
    fmneti(l)  = 0.0
    fluxupi(l) = 0.0
    fluxdni(l) = 0.0
end do
nfluxtopis(:) = 0.0
fmnetis(:,:)  = 0.0
fluxupis(:,:) = 0.0
fluxdnis(:,:) = 0.0

!     we now enter a major loop over spectral intervals in the infrared
!     to calculate the net flux in each spectral interval

ttop  = tlev(2)
tsurf = tlev(l_levels)

nts   = tsurf*10.0d0-499.0d0
ntt   = ttop *10.0d0-499.0d0

!calculate optional scale factor to make spectrally integrated upward IR flux at the surface equal to sigma*T^4
bsurftot=0.
do nw=1,l_nspecti
    bsurftot=bsurftot+dwni(nw)*(1.-rsfi)*planckir(nw,nts)
end do
bsurftot=3.14159*bsurftot
stfour=5.67e-8*tsurf**4.

do nw=1,l_nspecti

    !       surface emissions - independent of gauss points

    ! check whether to apply the scaling factor to bsurf
    if( do_irflux_scale) then
        bsurf = (stfour/bsurftot)*(1.-rsfi)*planckir(nw,nts)
    else
        bsurf = (1.-rsfi)*planckir(nw,nts)
    endif
    pltop = planckir(nw,ntt)

    !  if fzeroi(nw) = 1, then the k-coefficients are zero - skip to the
    !  special gauss point at the end.

    fzero = fzeroi(nw)
    if (.not.(fzero.ge.0.99)) then

        do ng=1,l_ngauss-1

            if(taugsurf(nw,ng).lt. tlimiti) then
                fzero = fzero + (1.0-fzeroi(nw))*gweight(ng)
                cycle
            end if

            !         set up the upper and lower boundary conditions on the ir
            !         calculate the downwelling radiation at the top of the model
            !         or the top layer will cool to space unphysically

            tautop = dtaui(1,nw,ng)*plev(2)/(plev(4)-plev(2))
            btop   = (1.0-exp(-tautop/ubari))*pltop

            !         we can now solve for the coefficients of the two stream
            !         call a subroutine that solves  for the flux terms
            !         within each interval at the midpoint wavenumber

            call gfluxi(nlevrad,tlev,nw,dwni(nw),dtaui(1,nw,ng),         &
                      taucumi(1,nw,ng),                                &
                      wbari(1,nw,ng),cosbi(1,nw,ng),rsfi,btop,   &
                      bsurf,ftopup,fmupi,fmdi)

            !         now calculate the cumulative ir net flux

            nfluxtopi = nfluxtopi+ftopup*dwni(nw)*gweight(ng)*           &
                                 (1.0-fzeroi(nw))
            nfluxtopis(nw) = nfluxtopis(nw) +                            &
                           ftopup*dwni(nw)*gweight(ng)*(1.0-fzeroi(nw))

            do l=1,l_nlevrad-1

                !           correct for the wavenumber intervals

                fmneti(l)  = fmneti(l)+(fmupi(l)-fmdi(l))*dwni(nw)*        &
                                        gweight(ng)*(1.0-fzeroi(nw))
                fmnetis(nw,l) = fmnetis(nw,l) +                            &
                                          (fmupi(l)-fmdi(l))*dwni(nw)*     &
                                           gweight(ng)*(1.0-fzeroi(nw))
                fluxupi(l) = fluxupi(l) + fmupi(l)*dwni(nw)*gweight(ng)*   &
                                          (1.0-fzeroi(nw))
                fluxupis(nw,l) = fluxupis(nw,l) +                          &
                                 fmupi(l)*dwni(nw)*gweight(ng)*            &
                                 (1.0-fzeroi(nw))
                fluxdni(l) = fluxdni(l) + fmdi(l)*dwni(nw)*gweight(ng)*    &
                                          (1.0-fzeroi(nw))
                fluxdnis(nw,l) = fluxdnis(nw,l) +                          &
                                 fmdi(l)*dwni(nw)*gweight(ng)*             &
                                 (1.0-fzeroi(nw))
            end do

        end do       !end ngauss loop

    endif

    !      special 17th gauss point

    ng     = l_ngauss

    tautop = dtaui(1,nw,ng)*plev(2)/(plev(4)-plev(2))
    btop   = (1.0-exp(-tautop/ubari))*pltop

    !      we can now solve for the coefficients of the two stream
    !      call a subroutine that solves  for the flux terms
    !      within each interval at the midpoint wavenumber

    call gfluxi(nlevrad,tlev,nw,dwni(nw),dtaui(1,nw,ng),            &
              taucumi(1,nw,ng),                                &
              wbari(1,nw,ng),cosbi(1,nw,ng),rsfi,btop,   &
              bsurf,ftopup,fmupi,fmdi)

    !      now calculate the cumulative ir net flux

    nfluxtopi = nfluxtopi+ftopup*dwni(nw)*fzero
    nfluxtopis(nw) = nfluxtopis(nw) + ftopup*dwni(nw)*fzero

    do l=1,l_nlevrad-1
        !        correct for the wavenumber intervals

        fmneti(l)  = fmneti(l)+(fmupi(l)-fmdi(l))*dwni(nw)*fzero
        fmnetis(nw,l) = fmnetis(nw,l)+(fmupi(l)-fmdi(l))*dwni(nw)*fzero
        fluxupi(l) = fluxupi(l) + fmupi(l)*dwni(nw)*fzero
        fluxupis(nw,l) = fluxupis(nw,l) + fmupi(l)*dwni(nw)*fzero
        fluxdni(l) = fluxdni(l) + fmdi(l)*dwni(nw)*fzero
        fluxdnis(nw,l) = fluxdnis(nw,l) + fmdi(l)*dwni(nw)*fzero
    end do

enddo      !End Spectral Interval LOOP

!! *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE INFRARED****

return
end subroutine sfluxi



!=====================================================================
!=====================================================================

subroutine gfluxv(dtdel,tdel,taucumin,wdel,cdel,ubar0,f0pi,rsf,  &
                btop,bsurf,fmidp,fmidm,diffv,fluxup,fluxdn,    &
                detau)
!  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
!  FOR THE VISIBLE  FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
!  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
!  MEASURED FROM THE TOP OF EACH LAYER.  (DTAU) TOP OF EACH LAYER HAS
!  OPTICAL DEPTH TAU(N).IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
!  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
!  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
! THIS SUBROUTINE DIFFERS FROM ITS IR COUNTERPART IN THAT HERE WE SOLVE FOR
! THE FLUXES DIRECTLY USING THE GENERALIZED NOTATION OF MEADOR AND WEAVOR
! J.A.S., 37, 630-642, 1980.
! THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY
! VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
!
! NLL           = NUMBER OF LEVELS (NAYER + 1) THAT WILL BE SOLVED
! NAYER         = NUMBER OF LAYERS (NOTE DIFFERENT SPELLING HERE)
! WAVEN         = WAVELENGTH FOR THE COMPUTATION
! DTDEL(NLAYER) = ARRAY OPTICAL DEPTH OF THE LAYERS
! TDEL(NLL)     = ARRAY COLUMN OPTICAL DEPTH AT THE LEVELS
! WDEL(NLEVEL)  = SINGLE SCATTERING ALBEDO
! CDEL(NLL)     = ASYMMETRY FACTORS, 0=ISOTROPIC
! UBARV         = AVERAGE ANGLE,
! UBAR0         = SOLAR ZENITH ANGLE
! F0PI          = INCIDENT SOLAR DIRECT BEAM FLUX
! RSF           = SURFACE REFLECTANCE
! BTOP          = UPPER BOUNDARY CONDITION ON DIFFUSE FLUX
! BSURF         = REFLECTED DIRECT BEAM = (1-RSFI)*F0PI*EDP-TAU/U
! FP(NLEVEL)    = UPWARD FLUX AT LEVELS
! FM(NLEVEL)    = DOWNWARD FLUX AT LEVELS
! FMIDP(NLAYER) = UPWARD FLUX AT LAYER MIDPOINTS
! FMIDM(NLAYER) = DOWNWARD FLUX AT LAYER MIDPOINTS
! added Dec 2002
! DIFFV         = downward diffuse solar flux at the surface
!
!======================================================================!

implicit none

integer, parameter :: nlp = 500     ! must be larger than nlevel



real*8 em, ep
real*8 w0(l_nlayrad), cosbar(l_nlayrad), dtau(l_nlayrad)
real*8 tau(l_nlevrad), wdel(l_nlayrad), cdel(l_nlayrad)
real*8 dtdel(l_nlayrad), tdel(l_nlevrad)
real*8 fmidp(l_nlayrad), fmidm(l_nlayrad)
real*8 lamda(nlp), alpha(nlp), xk1(nlp), xk2(nlp)
real*8 g1(nlp), g2(nlp), g3(nlp), gama(nlp), cp(nlp), cm(nlp)
real*8 cpm1(nlp)
real*8 cmm1(nlp), e1(nlp), e2(nlp), e3(nlp), e4(nlp), exptrm(nlp)
real*8 fluxup, fluxdn
real*8 factor, taucumin(l_levels), taucum(l_levels+2)
real*8  :: detau

integer nayer, l, k
real*8  ubar0, f0pi, rsf, btop, bsurf, g4, denom, am, ap
real*8  taumax, taumid, cpmid, cmmid
real*8  diffv

!======================================================================C

nayer  = l_nlayrad
taumax = l_taumax    !default is 35.0

!  delta-eddington scaling
factor    = 1.0 - wdel(1)*cdel(1)**2

tau(1)    = tdel(1)*factor
taucum(1) = 0.0
taucum(2) = taucumin(2)*factor
taucum(3) = taucum(2) +(taucumin(3)-taucumin(2))*factor

do l=1,l_nlayrad-1
    factor      = 1.0 - wdel(l)*cdel(l)**2
    w0(l)       = wdel(l)*(1.0-cdel(l)**2)/factor
    cosbar(l)   = cdel(l)/(1.0+cdel(l))
    dtau(l)     = dtdel(l)*factor
    tau(l+1)    = tau(l)+dtau(l)
    k           = 2*(l+1)
    taucum(k)   = tau(l+1)
    taucum(k+1) = taucum(k) + (taucumin(k+1)-taucumin(k))*factor
end do

!  Bottom layer

l             = l_nlayrad
factor        = 1.0 - wdel(l)*cdel(l)**2
w0(l)         = wdel(l)*(1.0-cdel(l)**2)/factor
cosbar(l)     = cdel(l)/(1.0+cdel(l))
dtau(l)       = dtdel(l)*factor
tau(l+1)      = tau(l)+dtau(l)
taucum(2*l+1) = tau(l+1)
detau         = taucum(2*l+1)



!     we go with the quadrature approach here.  the "sqrt(3)" factors
!     are the ubarv term.

do l=1,l_nlayrad
    !       set of constants determined by dom

    g1(l)    = (sqrt(3.0)*0.5)*(2.0- w0(l)*(1.0+cosbar(l)))
    g2(l)    = (sqrt(3.0)*w0(l)*0.5)*(1.0-cosbar(l))
    g3(l)    = 0.5*(1.0-sqrt(3.0)*cosbar(l)*ubar0)
    lamda(l) = sqrt(g1(l)**2 - g2(l)**2)
    gama(l)  = (g1(l)-lamda(l))/g2(l)
end do

do l=1,l_nlayrad
    g4    = 1.0-g3(l)
    denom = lamda(l)**2 - 1./ubar0**2

    !       there is a potential problem here if w0=0 and ubarv=ubar0
    !       then denom will vanish. this only happens physically when
    !       the scattering goes to zero
    !       prevent this with an if statement

    if ( denom .eq. 0.) then
        denom=1.e-10
    end if

    am = f0pi*w0(l)*(g4   *(g1(l)+1./ubar0) +g2(l)*g3(l) )/denom
    ap = f0pi*w0(l)*(g3(l)*(g1(l)-1./ubar0) +g2(l)*g4    )/denom

    !       cpm1 and cmm1 are the cplus and cminus terms evaluated
    !       at the top of the layer, that is lower   optical depth tau(l)

    cpm1(l) = ap*exp(-min(tau(l)/ubar0,maxexp))
    cmm1(l) = am*exp(-min(tau(l)/ubar0,maxexp))


!       cp and cm are the cplus and cminus terms evaluated at the
!       bottom of the layer.  that is at higher optical depth tau(l+1)

    cp(l) = ap*exp(-min(tau(l+1)/ubar0,maxexp))
    cm(l) = am*exp(-min(tau(l+1)/ubar0,maxexp))

end do

!     now calculate the exponential terms needed
!     for the tridiagonal rotated layered method

do l=1,l_nlayrad
    exptrm(l) = min(taumax,lamda(l)*dtau(l))  ! clipped exponential
    ep = exp(exptrm(l))

    em        = 1.0/ep
    e1(l)     = ep+gama(l)*em
    e2(l)     = ep-gama(l)*em
    e3(l)     = gama(l)*ep+em
    e4(l)     = gama(l)*ep-em
end do

call dsolver(nayer,gama,cp,cm,cpm1,cmm1,e1,e2,e3,e4,btop,        &
           bsurf,rsf,xk1,xk2)

!     now we calculate the fluxes at the midpoints of the layers.

do l=1,l_nlayrad-1
    exptrm(l) = min(taumax,lamda(l)*(taucum(2*l+1)-taucum(2*l)))

    ep = exp(exptrm(l))

    em    = 1.0/ep
    g4    = 1.0-g3(l)
    denom = lamda(l)**2 - 1./ubar0**2

    if ( denom .eq. 0.) then
        denom=1.e-10
    end if

    am = f0pi*w0(l)*(g4   *(g1(l)+1./ubar0) +g2(l)*g3(l) )/denom
    ap = f0pi*w0(l)*(g3(l)*(g1(l)-1./ubar0) +g2(l)*g4    )/denom

    taumid   = taucum(2*l+1)

    cpmid = ap*exp(-min(taumid/ubar0,maxexp))
    cmmid = am*exp(-min(taumid/ubar0,maxexp))


    fmidp(l) = xk1(l)*ep + gama(l)*xk2(l)*em + cpmid
    fmidm(l) = xk1(l)*ep*gama(l) + xk2(l)*em + cmmid

    !       add the direct flux to the downwelling term

    fmidm(l)= fmidm(l)+ubar0*f0pi*exp(-min(taumid/ubar0,maxexp))

end do

!     flux at the top layer

ep    = 1.0
em    = 1.0
g4    = 1.0-g3(1)
denom = lamda(1)**2 - 1./ubar0**2


if ( denom .eq. 0.) then
    denom=1.e-10
end if

am = f0pi*w0(1)*(g4   *(g1(1)+1./ubar0) +g2(1)*g3(1) )/denom
ap = f0pi*w0(1)*(g3(1)*(g1(1)-1./ubar0) +g2(1)*g4    )/denom

cpmid  = ap
cmmid  = am

fluxup = xk1(1)*ep + gama(1)*xk2(1)*em + cpmid
fluxdn = xk1(1)*ep*gama(1) + xk2(1)*em + cmmid

!     add the direct flux to the downwelling term

fluxdn = fluxdn+ubar0*f0pi*exp(-min(taucum(1)/ubar0,maxexp))

!     this is for the "special" bottom layer, where we take
!     dtau instead of dtau/2.

l     = l_nlayrad
exptrm(l) = min(taumax,lamda(l)*(taucum(l_levels)-               &
                               taucum(l_levels-1)))

ep    = exp(exptrm(l))
em    = 1.0/ep
g4    = 1.0-g3(l)
denom = lamda(l)**2 - 1./ubar0**2

if ( denom .eq. 0.) then
    denom=1.e-10
end if

am = f0pi*w0(l)*(g4   *(g1(l)+1./ubar0) +g2(l)*g3(l) )/denom
ap = f0pi*w0(l)*(g3(l)*(g1(l)-1./ubar0) +g2(l)*g4    )/denom

taumid   = min(taucum(l_levels),taumax)
cpmid    = ap*exp(-min(taumid/ubar0,maxexp))
cmmid    = am*exp(-min(taumid/ubar0,maxexp))

fmidp(l) = xk1(l)*ep + gama(l)*xk2(l)*em + cpmid
fmidm(l) = xk1(l)*ep*gama(l) + xk2(l)*em + cmmid

!  save the diffuse downward flux for tempgr calculations

diffv = fmidm(l)

!     add the direct flux to the downwelling term

fmidm(l)= fmidm(l)+ubar0*f0pi*exp(-min(taumid/ubar0,maxexp))

return
end subroutine gfluxv


!=====================================================================
!=====================================================================

subroutine gfluxi(nll,tlev,nw,dw,dtau,taucum,w0,cosbar,    &
                rsf,btop,bsurf,ftopup,fmidp,fmidm)
                !
!  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
!  FOR THE INFRARED FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
!  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
!  MEASURED FROM THE TOP OF EACH LAYER.  THE TOP OF EACH LAYER HAS
!  OPTICAL DEPTH ZERO.  IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
!  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
!  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
! THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY
! VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
!
! NLL            = NUMBER OF LEVELS (NLAYERS + 1) MUST BE LESS THAT NL (101)
! TLEV(L_LEVELS) = ARRAY OF TEMPERATURES AT GCM LEVELS
! WAVEN          = WAVELENGTH FOR THE COMPUTATION
! DW             = WAVENUMBER INTERVAL
! DTAU(NLAYER)   = ARRAY OPTICAL DEPTH OF THE LAYERS
! W0(NLEVEL)     = SINGLE SCATTERING ALBEDO
! COSBAR(NLEVEL) = ASYMMETRY FACTORS, 0=ISOTROPIC
! UBARI          = AVERAGE ANGLE, MUST BE EQUAL TO 0.5 IN IR
! RSF            = SURFACE REFLECTANCE
! BTOP           = UPPER BOUNDARY CONDITION ON IR INTENSITY (NOT FLUX)
! BSURF          = SURFACE EMISSION = (1-RSFI)*PLANCK, INTENSITY (NOT FLUX)
! FP(NLEVEL)     = UPWARD FLUX AT LEVELS
! FM(NLEVEL)     = DOWNWARD FLUX AT LEVELS
! FMIDP(NLAYER)  = UPWARD FLUX AT LAYER MIDPOINTS
! FMIDM(NLAYER)  = DOWNWARD FLUX AT LAYER MIDPOINTS
!
!----------------------------------------------------------------------C


implicit none

integer, parameter :: nlp = 500    ! must be larger than nlevel

integer :: nll, nlayer, l, nw, nt, nt2
real*8  :: term, cpmid, cmmid
real*8  :: planck
real*8  :: em,ep
real*8  :: cosbar(l_nlayrad), w0(l_nlayrad), dtau(l_nlayrad)
real*8  :: taucum(l_levels), dtauk
real*8  :: tlev(l_levels)
real*8  :: waven, dw, rsf
real*8  :: btop, bsurf, fmidp(l_nlayrad), fmidm(l_nlayrad)
real*8  :: b0(nlp),b1(nlp),alpha(nlp),lamda(nlp),xk1(nlp),xk2(nlp)
real*8  :: gama(nlp),cp(nlp),cm(nlp),cpm1(nlp),cmm1(nlp),e1(nlp)
real*8  :: e2(nlp),e3(nlp),e4(nlp)

real*8  :: ftopup, fluxup, fluxdn
real*8  :: taumax = l_taumax

!======================================================================C

!     we go with the hemispheric constant approach in the infrared

if (nll .gt. nlp) stop 'parameter nl too small in glufv'

nlayer = l_nlayrad

do l=1,l_nlayrad-1
    alpha(l) = sqrt( (1.0-w0(l))/(1.0-w0(l)*cosbar(l)) )
    lamda(l) = alpha(l)*(1.0-w0(l)*cosbar(l))/ubari

    nt2   = tlev(2*l+2)*10.0d0-499
    nt    = tlev(2*l)*10.0d0-499

    b1(l) = (planckir(nw,nt2)-planckir(nw,nt))/dtau(l)
    b0(l) = planckir(nw,nt)
end do

!     take care of special lower layer

l        = l_nlayrad
alpha(l) = sqrt( (1.0-w0(l))/(1.0-w0(l)*cosbar(l)) )
lamda(l) = alpha(l)*(1.0-w0(l)*cosbar(l))/ubari

nt    = tlev(2*l+1)*10.0d0-499
nt2   = tlev(2*l)*10.0d0-499
b1(l) = (planckir(nw,nt)-planckir(nw,nt2))/dtau(l)
b0(l) = planckir(nw,nt2)

do l=1,l_nlayrad
    gama(l) = (1.0-alpha(l))/(1.0+alpha(l))
    term    = ubari/(1.0-w0(l)*cosbar(l))

    !       cp and cm are the cplus and cminus terms evaluated at the
    !       bottom of the layer.  that is at dtau optical depth

    cp(l) = b0(l)+b1(l)*dtau(l) +b1(l)*term
    cm(l) = b0(l)+b1(l)*dtau(l) -b1(l)*term

    !       cpm1 and cmm1 are the cplus and cminus terms evaluated
    !       at the top of the layer, that is zero optical depth

    cpm1(l) = b0(l)+b1(l)*term
    cmm1(l) = b0(l)-b1(l)*term
end do

!     now calculate the exponential terms needed
!     for the tridiagonal rotated layered method
!     warning if dtau(j) is greater than about 35 (vax)
!     we clip it to avoid overflow.

do l=1,l_nlayrad

    !       clip the exponential here.

    ep    = exp( min((lamda(l)*dtau(l)),taumax))
    em    = 1.0/ep
    e1(l) = ep+gama(l)*em
    e2(l) = ep-gama(l)*em
    e3(l) = gama(l)*ep+em
    e4(l) = gama(l)*ep-em
end do

!     double precision tridiagonal solver

call dsolver(nlayer,gama,cp,cm,cpm1,cmm1,e1,e2,e3,e4,btop,       &
           bsurf,rsf,xk1,xk2)

!     now we calculate the fluxes at the midpoints of the layers.

do l=1,l_nlayrad-1
    dtauk = taucum(2*l+1)-taucum(2*l)
    ep    = exp(min(lamda(l)*dtauk,taumax)) ! clipped exponential
    em    = 1.0/ep
    term  = ubari/(1.-w0(l)*cosbar(l))

    !       cp and cm are the cplus and cminus terms evaluated at the
    !       bottom of the layer.  that is at dtau  optical depth

    cpmid    = b0(l)+b1(l)*dtauk +b1(l)*term
    cmmid    = b0(l)+b1(l)*dtauk -b1(l)*term
    fmidp(l) = xk1(l)*ep + gama(l)*xk2(l)*em + cpmid
    fmidm(l) = xk1(l)*ep*gama(l) + xk2(l)*em + cmmid

    !       for flux we integrate over the hemisphere treating intensity constant

    fmidp(l) = fmidp(l)*pi
    fmidm(l) = fmidm(l)*pi
end do

!     and now, for the special bottom layer

l    = l_nlayrad

ep   = exp(min((lamda(l)*dtau(l)),taumax)) ! clipped exponential
em   = 1.0/ep
term = ubari/(1.-w0(l)*cosbar(l))

!     cp and cm are the cplus and cminus terms evaluated at the
!     bottom of the layer.  that is at dtau  optical depth

cpmid    = b0(l)+b1(l)*dtau(l) +b1(l)*term
cmmid    = b0(l)+b1(l)*dtau(l) -b1(l)*term
fmidp(l) = xk1(l)*ep + gama(l)*xk2(l)*em + cpmid
fmidm(l) = xk1(l)*ep*gama(l) + xk2(l)*em + cmmid

!     for flux we integrate over the hemisphere treating intensity constant

fmidp(l) = fmidp(l)*pi
fmidm(l) = fmidm(l)*pi

!     flux at the top level

ep   = 1.0
em   = 1.0
term = ubari/(1.0-w0(1)*cosbar(1))

!     cp and cm are the cplus and cminus terms evaluated at the
!     bottom of the layer.  that is at dtau  optical depth

cpmid  = b0(1)+b1(1)*term
cmmid  = b0(1)-b1(1)*term

fluxup = xk1(1)*ep + gama(1)*xk2(1)*em + cpmid
fluxdn = xk1(1)*ep*gama(1) + xk2(1)*em + cmmid

!     for flux we integrate over the hemisphere treating intensity constant

ftopup = (fluxup-fluxdn)*pi

return
end subroutine gfluxi

!=====================================================================
!=====================================================================

subroutine getdetau(dtdel,tdel,taucumin,wdel,cdel,ubar0,detau)
!
!   calculate delta eddington tau
!
!======================================================================!

implicit none

real*8  :: w0(l_nlayrad), cosbar(l_nlayrad), dtau(l_nlayrad)
real*8  :: tau(l_nlevrad), wdel(l_nlayrad), cdel(l_nlayrad)
real*8  :: dtdel(l_nlayrad), tdel(l_nlevrad)
real*8  :: factor, taucumin(l_levels), taucum(l_levels)
real*8  :: detau

integer :: nayer, l, k
real*8  :: ubar0

!======================================================================c

nayer  = l_nlayrad

!  delta-eddington scaling

factor    = 1.0d0 - wdel(1)*cdel(1)**2

tau(1)    = tdel(1)*factor
taucum(1) = 0.0d0
taucum(2) = taucumin(2)*factor
taucum(3) = taucum(2) +(taucumin(3)-taucumin(2))*factor

do l=1,l_nlayrad-1
    factor      = 1.0d0 - wdel(l)*cdel(l)**2
    w0(l)       = wdel(l)*(1.0d0-cdel(l)**2)/factor
    cosbar(l)   = cdel(l)/(1.0d0+cdel(l))
    dtau(l)     = dtdel(l)*factor
    tau(l+1)    = tau(l)+dtau(l)
    k           = 2*(l+1)
    taucum(k)   = tau(l+1)
    taucum(k+1) = taucum(k) + (taucumin(k+1)-taucumin(k))*factor
end do

!  bottom layer

l             = l_nlayrad
factor        = 1.0d0 - wdel(l)*cdel(l)**2
w0(l)         = wdel(l)*(1.0d0-cdel(l)**2)/factor
cosbar(l)     = cdel(l)/(1.0d0+cdel(l))
dtau(l)       = dtdel(l)*factor
tau(l+1)      = tau(l)+dtau(l)
taucum(2*l+1) = tau(l+1)
detau         = taucum(2*l+1)

return
end subroutine getdetau



!=====================================================================
!=====================================================================

subroutine tpindex(pw,tw,qh2o,pref,tref,lcoef,mt,mp,     &
                 nh2o,wratio)
!
!    Get the TI, UI values for a 2-dimensional interpolation
!    based on the following (The interpolation is done in interpco2):
!    Interpolate the CO2 K-coefficients to the current P,T values.
!    The CO2 coefficients are given on a P,T grid:
!    P = {1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 1E+1, 1E+2, 1E+3, 1E+4},
!    T = {50, 100, 150, 200, 250, 300, 350}.
!
!    The interpolation is the usual interpolation in 2-dimensions given
!    in "Numerical Recipes", where the "X" are P, the "Y" are
!    T, and the F(X,Y) are the CO2 K-coefficients.
!
!     The interpolating box is designated as follows:
!
!           (PL,TU)                        (PRR,TU)
!
!                          (TW,PW)
!
!
!           (PL,TL)                        (PRR,TL)
!
!     PL  - Pressure left
!     PRR - Pressure right
!     TL  - Temperature lower
!     TU  - Temperature upper
!     PW  - Pressure wanted
!     TW  - Temperature wanted
!
!
!  OUTPUT PARAMETERS
!    TI                 - Interpolation term (pressure)
!    UI                 - Interpolation term (temperature)
!    MT                 - Temperature index (bottom left Temperature)
!                         of bounding box
!    MP                 - Pressure index (bottom left pressure)
!                         of bounding box
!
!----------------------------------------------------------------------C

implicit none

real*8  :: Tref(L_NTREF)        ! The temperature grid array.
real*8  :: pref(L_PINT)         ! The pressure grid array.

integer :: MT, MP, N, M, NP, NH2O
real*8  :: PW, &                ! The pressure to interpolate to
            TW, &               ! The temperature to interpolate to
            Qh2o, &
            wratio
real*8  :: PWL, LCOEF(4), T, U

!======================================================================C

!     Get the upper and lower Temperature-grid indicies that bound the
!     requested temperature.  If the requested temperature is outside
!     the T-grid, set up to extrapolate from the appropriate end.

if(tw.le.tref(1)) then
    mt = 1
else
    do n=1,l_ntref-1
        if((tw.gt.tref(n) .and. tw.le.tref(n+1)).or.(n.eq.l_ntref-1)) then
            mt = n
            exit
        else
            cycle
        end if
    end do
end if

u = (tw-tref(mt))/(tref(mt+1)-tref(mt))

!     get the upper and lower pressure-grid indicies that bound the
!     requested pressure.  if the requested pressure is outside
!     the p-grid, set up to extrapolate from the appropiate end.

pwl = log10(pw)

do n=2,l_pint-1
    if((pwl.le.pref(n)).or.(n.eq.l_pint-1)) then
        mp = n-1
        exit
    else
        cycle
    end if
end do

t = (pwl-pref(mp))/(pref(mp+1)-pref(mp))

!  fill the interpolation coeficients:

lcoef(1) = (1.0-t)*(1.0-u)
lcoef(2) = t*(1.0-u)
lcoef(3) = t*u
lcoef(4) = (1.0-t)*u

!  get the indicies for water abundance.  there are 10 sets of
!  k-coefficients with differing amounts of water vs. co2.

if(qh2o.le.wrefh2o(1)) then
    nh2o   = 1
    wratio = 0.0d0
elseif(qh2o.ge.wrefh2o(l_refh2o)) then
    nh2o   = l_refh2o
    wratio = 0.0d0
else
    do n=2,l_refh2o
        if((qh2o.ge.wrefh2o(n-1) .and. qh2o.lt.wrefh2o(n)).or.(n.eq.l_refh2o)) then
            nh2o   = n-1
            wratio = (qh2o - wrefh2o(n-1))/(wrefh2o(n) - wrefh2o(n-1))
            exit
        else
            cycle
        end if
    end do
end if

return
end subroutine tpindex

!=====================================================================
!=====================================================================

subroutine dsolver(nl,gama,cp,cm,cpm1,cmm1,e1,e2,e3,e4,btop,     &
                 bsurf,rsf,xk1,xk2)
! DOUBLE PRECISION VERSION OF SOLVER
!* THIS SUBROUTINE SOLVES FOR THE COEFFICIENTS OF THE    *
!* TWO STREAM SOLUTION FOR GENERAL BOUNDARY CONDITIONS   *
!* NO ASSUMPTION OF THE DEPENDENCE ON OPTICAL DEPTH OF   *
!* C-PLUS OR C-MINUS HAS BEEN MADE.                      *
implicit none

integer, parameter :: nmax=201

integer :: l, nl, i, lm2, lm1, n
real*8  :: gama(nl),cp(nl),cm(nl),cpm1(nl),cmm1(nl),xk1(nl)
real*8  :: xk2(nl),e1(nl),e2(nl),e3(nl),e4(nl)
real*8  :: af(nmax),bf(nmax),cf(nmax),df(nmax),xk(nmax)
real*8  :: btop, bsurf, rsf

!*********************************************************
!* NL     = NUMBER OF LAYERS IN THE MODEL                *
!* CP     = C-PLUS EVALUATED AT TAO=0 (TOP)              *
!* CM     = C-MINUS EVALUATED AT TAO=0 (TOP)             *
!* CPM1   = C-PLUS  EVALUATED AT TAOSTAR (BOTTOM)        *
!* CMM1   = C-MINUS EVALUATED AT TAOSTAR (BOTTOM)        *
!* EP     = EXP(LAMDA*DTAU)                              *
!* EM     = 1/EP                                         *
!* E1     = EP + GAMA *EM                                *
!* E2     = EP - GAMA *EM                                *
!* E3     = GAMA*EP + EM                                 *
!* E4     = GAMA*EP - EM                                 *
!* BTOP   = THE DIFFUSE RADIATION INTO THE MODEL AT TOP  *
!* BSURF  = THE DIFFUSE RADIATION INTO THE MODEL AT      *
!*          THE BOTTOM: INCLUDES EMMISION AND REFLECTION *
!*          OF THE UNATTENUATED PORTION OF THE DIRECT    *
!*          BEAM. BSTAR+RSF*FO*EXP(-TAOSTAR/U0)          *
!* RSF    = REFLECTIVITY OF THE SURFACE                  *
!* XK1    = COEFFICIENT OF THE POSITIVE EXP TERM         *
!* XK2    = COEFFICIENT OF THE NEGATIVE EXP TERM         *
!*********************************************************

!======================================================================C

l=2*nl

!     ************mixed coefficents**********
!     this version avoids singularities assoc.
!     with w0=0 by solving for xk1+xk2, and xk1-xk2.

af(1) = 0.0
bf(1) = gama(1)+1.
cf(1) = gama(1)-1.
df(1) = btop-cmm1(1)
n     = 0
lm2   = l-2

!     even terms

do i=2,lm2,2
    n     = n+1
    af(i) = (e1(n)+e3(n))*(gama(n+1)-1.)
    bf(i) = (e2(n)+e4(n))*(gama(n+1)-1.)
    cf(i) = 2.0*(1.-gama(n+1)**2)
    df(i) = (gama(n+1)-1.) * (cpm1(n+1) - cp(n)) +                 &
              (1.-gama(n+1))* (cm(n)-cmm1(n+1))
end do

n   = 0
lm1 = l-1
do i=3,lm1,2
    n     = n+1
    af(i) = 2.0*(1.-gama(n)**2)
    bf(i) = (e1(n)-e3(n))*(1.+gama(n+1))
    cf(i) = (e1(n)+e3(n))*(gama(n+1)-1.)
    df(i) = e3(n)*(cpm1(n+1) - cp(n)) + e1(n)*(cm(n) - cmm1(n+1))
end do

af(l) = e1(nl)-rsf*e3(nl)
bf(l) = e2(nl)-rsf*e4(nl)
cf(l) = 0.0
df(l) = bsurf-cp(nl)+rsf*cm(nl)

call dtridgl(l,af,bf,cf,df,xk)

!     ***unmix the coefficients****

do n=1,nl
    xk1(n) = xk(2*n-1)+xk(2*n)
    xk2(n) = xk(2*n-1)-xk(2*n)

    !       now test to see if xk2 is really zero to the limit of the
    !       machine accuracy  = 1 .e -30
    !       xk2 is the coefficeint of the growing exponential and must
    !       be treated carefully

    if(xk2(n) .eq. 0.0) cycle
    if (abs (xk2(n)/(xk(2*n-1)+1.e-20)) .lt. 1.e-30) xk2(n)=0.0

enddo

return
end subroutine dsolver



!=====================================================================
!=====================================================================

subroutine dtridgl(l,af,bf,cf,df,xk)
!*    this subroutine solves a system of tridiagional matrix
!*    equations. the form of the equations are:
!*    a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i)
!*    where i=1,l  less than 103.

implicit none

!     double prescision version of tridgl

integer, parameter :: nmax = 201
integer :: l, i

real*8  :: af(l),bf(l),cf(l),df(l),xk(l)
real*8  :: as(nmax),ds(nmax), x, xkb

!======================================================================C

as(l) = af(l)/bf(l)
ds(l) = df(l)/bf(l)

do i=2,l
    x         = 1./(bf(l+1-i) - cf(l+1-i)*as(l+2-i))
    as(l+1-i) = af(l+1-i)*x
    ds(l+1-i) = (df(l+1-i)-cf(l+1-i)*ds(l+2-i))*x
end do

xk(1)=ds(1)
do i=2,l
    xkb   = xk(i-1)
    xk(i) = ds(i)-as(i)*xkb
end do

return
end subroutine dtridgl

!=====================================================================
!=====================================================================

subroutine settozero(nsize,array)
! set r8 input array to zero
implicit none

integer nsize
real*8 array(nsize)

integer i

do i=1,nsize
    array(i) = 0.
enddo

return

end subroutine settozero

!=====================================================================
!=====================================================================

subroutine settozero4(nsize,array)
! set r4 input array to 0
implicit none

integer nsize
real*4 array(nsize)

integer i

do i=1,nsize
    array(i) = 0.0
enddo

return

end subroutine settozero4
!=====================================================================
!=====================================================================

subroutine initcld(mcpu0)
!
!  Initialize the variables used in the cloud scheme

implicit none


integer i,j,l,nw,ii
logical mcpu0

real*8 :: rmin  = 0.1e-6
real*8 :: rmax  = 10.e-6
real*8 :: rbmin = 0.0001e-6
real*8 :: rbmax = 1.e-2

real*8 vrat_rt
real*8 vrat_rtblk
real*8 vrat_rtco2

real*8 factor

real*8 :: vistoir = 2.75

real*8 a0

!======================================================================C

!     Definition of the size grid
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     aerad is the primary radius grid used for microphysics computation.
!     the grid spacing is based on a volume ratio between two
!     consecutive bins; i.e. vrat.
!     rb defines the boundary values for each aerad bin.

!     Volume ratio between two adjacent bins
vrat = log(rmax/rmin) / float(nbin-1) *3.
vrat = exp(vrat)

rb(1)      = rbmin
aerad(1)   = rmin
vol(1)     = 4./3. * pi * rmin**3.

do i=1,nbin-1
    aerad(i+1)  = aerad(i) * vrat**(athird)
    vol(i+1)    = vol(i) * vrat
enddo

do i=1,nbin
    rb(i+1)= ( (2.*vrat) / (vrat+1.) )**(athird) * aerad(i)
    dr(i)  = rb(i+1) - rb(i)
enddo
rb(nbin+1) = rbmax
dr(nbin)   = rb(nbin+1) - rb(nbin)

if (mcpu0) then
    print*, ' '
    print*, ' Here are the size bins (lower bnd, center values & dr):'
    print*, ' -------------------------------------------------------'
    do i=1,nbin
        write(*,'(i2,3x,3(e12.6,4x))') i,rb(i), aerad(i),dr(i)
    enddo
    write(*,'(i2,3x,e12.6)') nbin+1,rb(nbin+1)
    print*, ' ----------------------------------------------------'
endif

!     Now, initialize the water ice cloud radiative properties.
!     Here, we use a different size grid (more refined) called rad_rt.

rad_rt(1)      = 1.e-7
rad_rt(nbin_rt)= 50.e-6
radb_rt(1)     = rbmin
radb_rt(nbin_rt+1)= rbmax

vrat_rt = log(rad_rt(nbin_rt)/rad_rt(1)) / float(nbin_rt-1) *3.
vrat_rt = exp(vrat_rt)

do i = 1, nbin_rt-1
    rad_rt(i+1)  = rad_rt(i) * vrat_rt**(athird)
    radb_rt(i+1)=((2.*vrat_rt) / (vrat_rt+1.))**(athird) * rad_rt(i)
enddo

! Bulk scheme properties, 1 nratio, 29 particle sizes
rad_rtblk(1)      = 1.e-7
rad_rtblk(20)=  50.e-6          !same first 20 bins as normal, extends 9 sizes larger
radb_rtblk(1)     = rbmin
radb_rtblk(nbin_rtblk+1)= rbmax

vrat_rtblk=log(rad_rtblk(20)/rad_rtblk(1)) / float(19) *3.
vrat_rtblk=exp(vrat_rtblk)

do i = 1, nbin_rtblk-1
    rad_rtblk(i+1) = rad_rtblk(i) * vrat_rtblk**(athird)
    radb_rtblk(i+1)=((2.*vrat_rtblk) / (vrat_rtblk+1.))**(athird) * rad_rtblk(i)
enddo


rad_rtco2(1)      = 1.e-7
rad_rtco2(nbin_rtco2)= 260.e-6
radb_rtco2(1)     = rbmin
radb_rtco2(nbin_rtco2+1)= rbmax

vrat_rtco2=log(rad_rtco2(nbin_rtco2)/rad_rtco2(1)) / float(nbin_rtco2-1) *3.
vrat_rtco2=exp(vrat_rtco2)

do i = 1, nbin_rtco2-1
    rad_rtco2(i+1) = rad_rtco2(i) * vrat_rtco2**(athird)
    radb_rtco2(i+1)=((2.*vrat_rtco2) / (vrat_rtco2+1.))**(athird) * rad_rtco2(i)
enddo

!     Read in the data files the Qext,Qscat and g values for each
!     size bins and for each spectral interval.

if (ames_15band .OR. use_extended_cor_ks) then
    open(60,file=trim(rtdata_path)//'waterCoated_vis_JD_15bands.dat')
    open(61,file=trim(rtdata_path)//'waterCoated_ir_JD_15bands.dat')
    open(62,file=trim(rtdata_path)//'Dust_vis_wolff2010_JD_15bands.dat')
    open(63,file=trim(rtdata_path)//'Dust_ir_wolff2010_JD_15bands.dat')

    do j = 1, nratio
        do i = 1, nbin_rt
            read(60,'(7(e12.7,x))') (qextv_cld(j,i,l), l=1,nlonv)
            read(60,'(7(e12.7,x))') (qscatv_cld(j,i,l), l=1,nlonv)
            read(60,'(7(e12.7,x))') (gv_cld(j,i,l), l=1,nlonv)
            read(61,'(8(e12.7,x))') (qexti_cld(j,i,l), l=1,nloni)
            read(61,'(8(e12.7,x))') (qscati_cld(j,i,l), l=1,nloni)
            read(61,'(8(e12.7,x))') (gi_cld(j,i,l), l=1,nloni)
        enddo
    enddo

    do i = 1, nbin_rt
        read(62,'(7(e10.3,2x))') (qextv_dst(i,l), l=1,nlonv)
        read(62,'(7(e10.3,2x))') (qscatv_dst(i,l), l=1,nlonv)
        read(62,'(7(e10.3,2x))') (gv_dst(i,l), l=1,nlonv)
        read(63,'(8(e10.3,2x))') (qexti_dst(i,l), l=1,nloni)
        read(63,'(8(e10.3,2x))') (qscati_dst(i,l), l=1,nloni)
        read(63,'(8(e10.3,2x))') (gi_dst(i,l), l=1,nloni)
        !        factor  = qextv_dst(i,L_NREFV) / qexti_dst(i,L_NREFI)
        factor = 1.0
        do l = 1, nloni
            qexti_dst(i,l)  = qexti_dst(i,l)  * factor
            qscati_dst(i,l) = qscati_dst(i,l) * factor
        enddo
    enddo
    close(60)
    close(61)
    close(62)
    close(63)

else
    open(60,file=trim(rtdata_path)//'waterCoated_vis_JD_12bands.dat')
    open(61,file=trim(rtdata_path)//'waterCoated_ir_JD_12bands.dat')
    open(62,file=trim(rtdata_path)//'Dust_vis_wolff2010_JD_12bands.dat')
    open(63,file=trim(rtdata_path)//'Dust_ir_wolff2010_JD_12bands.dat')

    do j = 1, nratio
        do i = 1, nbin_rt
            read(60,'(7(e12.7,x))') (qextv_cld(j,i,l), l=1,nlonv)
            read(60,'(7(e12.7,x))') (qscatv_cld(j,i,l), l=1,nlonv)
            read(60,'(7(e12.7,x))') (gv_cld(j,i,l), l=1,nlonv)
            read(61,'(5(e12.7,x))') (qexti_cld(j,i,l), l=1,nloni)
            read(61,'(5(e12.7,x))') (qscati_cld(j,i,l), l=1,nloni)
            read(61,'(5(e12.7,x))') (gi_cld(j,i,l), l=1,nloni)
        enddo
    enddo

    do i = 1, nbin_rt
        read(62,'(7(e11.5,x))') (qextv_dst(i,l), l=1,nlonv)
        read(62,'(7(e11.5,x))') (qscatv_dst(i,l), l=1,nlonv)
        read(62,'(7(e11.5,x))') (gv_dst(i,l), l=1,nlonv)
        read(63,'(5(e11.5,x))') (qexti_dst(i,l), l=1,nloni)
        read(63,'(5(e11.5,x))') (qscati_dst(i,l), l=1,nloni)
        read(63,'(5(e11.5,x))') (gi_dst(i,l), l=1,nloni)
        !        factor  = qextv_dst(i,6) / (vistoir*qexti_dst(i,4))
        factor  = 1.0
        do l = 1, nloni
            qexti_dst(i,l)  = qexti_dst(i,l)  * factor
            qscati_dst(i,l) = qscati_dst(i,l) * factor
        enddo
    enddo

  close(60)
  close(61)
  close(62)
  close(63)
endif


if (radactive_cloud_bulk) then
    open(60,file=trim(rtdata_path)//'h2oIceCoated_vis_KS_15bands29Radii_1nratio.dat')
    open(61,file=trim(rtdata_path)//'h2oIceCoated_ir_KS_15bands29Radii_1nratio.dat')
    open(62,file=trim(rtdata_path)//'h2oLiqCoated_vis_KS_15bands29Radii_1nratio.dat')
    open(63,file=trim(rtdata_path)//'h2oLiqCoated_ir_KS_15bands29Radii_1nratio.dat')

!   read in optical properties for bulk clouds both ice and liquid 
    do j = 1, nratioblk
        do i = 1, nbin_rtblk
            read(60,'(7(e12.7,x))') (qextv_bcld(j,i,l), l=1,nlonv)
            read(60,'(7(e12.7,x))') (qscatv_bcld(j,i,l), l=1,nlonv)
            read(60,'(7(e12.7,x))') (gv_bcld(j,i,l), l=1,nlonv)
            read(61,'(8(e12.7,x))') (qexti_bcld(j,i,l), l=1,nloni)
            read(61,'(8(e12.7,x))') (qscati_bcld(j,i,l), l=1,nloni)
            read(61,'(8(e12.7,x))') (gi_bcld(j,i,l), l=1,nloni)
            read(62,'(7(e12.7,x))') (qextv_blcld(j,i,l), l=1,nlonv)
            read(62,'(7(e12.7,x))') (qscatv_blcld(j,i,l), l=1,nlonv)
            read(62,'(7(e12.7,x))') (gv_blcld(j,i,l), l=1,nlonv)
            read(63,'(8(e12.7,x))') (qexti_blcld(j,i,l), l=1,nloni)
            read(63,'(8(e12.7,x))') (qscati_blcld(j,i,l), l=1,nloni)
            read(63,'(8(e12.7,x))') (gi_blcld(j,i,l), l=1,nloni)
        enddo
    enddo
    close(60)
    close(61)
    close(62)
    close(63)
endif

if (radactive_co2cloud) then
    open(64,file=trim(rtdata_path)//'co2Coated_vis_MK_15bands29radii_1nratio.dat')
    open(65,file=trim(rtdata_path)//'co2Coated_ir_MK_15bands29radii_1nratio.dat')
    
    !read in co2 cloud files
    do j = 1, nratioblk
        do i = 1, nbin_rtco2
            read(64,'(7(e12.7,x))') (qextv_co2cld(j,i,l), l=1,nlonv)
            read(64,'(7(e12.7,x))') (qscatv_co2cld(j,i,l), l=1,nlonv)
            read(64,'(7(e12.7,x))') (gv_co2cld(j,i,l), l=1,nlonv)
            read(65,'(8(e12.7,x))') (qexti_co2cld(j,i,l), l=1,nloni)
            read(65,'(8(e12.7,x))') (qscati_co2cld(j,i,l), l=1,nloni)
            read(65,'(8(e12.7,x))') (gi_co2cld(j,i,l), l=1,nloni)
        enddo
    enddo
    close(64)
    close(65)
endif

if (do_cia) then

    open(70,file=trim(rtdata_path)//'kgbar_8band_100_800.dat',form='formatted')
    open(71,file=trim(rtdata_path)//'kbbar_8band_100_800.dat',form='formatted')
    open(72,file=trim(rtdata_path)//'khbar_8band_T20_100_800.dat',form='formatted') !Turbet et al. (2020)
!    open(72,file=trim(rtdata_path)//'khbar_8band_W17_100_800.dat',form='formatted') ! Wordsworth et al. (2017)
    do nw=1,L_NSPECTI
        read(70,656) (kgbar_tab(nw,j),j=1,15)
        read(71,656) (kbbar_tab(nw,j),j=1,15)
        read(72,656) (khbar_tab(nw,j),j=1,15)
  656   format(1x,15(1pe13.3))
    end do
    close(70)
    close(71)
    close(72)

endif

open(60,file=trim(rtdata_path)//'waterCoatedQ_Tbright_MK.dat')
open(61,file=trim(rtdata_path)//'DustQ_Tbright_MK.dat')

do j= 1, nratio
    do i = 1, nbin_rt
      read(60,'(5(e12.7,x))') (qextTb_cld(j,i,l),l=1,5)
      read(60,'(5(e12.7,x))') (qscatTb_cld(j,i,l),l=1,5)
      read(60,'(5(e12.7,x))') (gTb_cld(j,i,l),l=1,5)
  !!!   print*,qextTb_cld(j,i,1),qextTb_cld(j,i,2),qextTb_cld(j,i,3),&
  !!!          qextTb_cld(j,i,4),qextTb_cld(j,i,5)
    enddo
enddo

do i = 1, nbin_rt
      read(61,'(5(e12.7,x))') (qextTb_dst(i,l),l=1,5)
      read(61,'(5(e12.7,x))') (qscatTb_dst(i,l),l=1,5)
      read(61,'(5(e12.7,x))') (gTb_dst(i,l),l=1,5)
!!!     print*,qextTb_dst(i,1),qextTb_dst(i,2),qextTb_dst(i,3),&
!!!            qextTb_dst(i,4),qextTb_dst(i,5)
enddo

close(60)
close(61)

open(60,file=trim(rtdata_path)//'QabsFCorrect.dat', access='sequential',form='unformatted')
read(60) fcorrect
close(60)

!  rjw    Copy input arrays into working arrays ( real*4 --> real*4)

 do ii= 1, nlond
        qextd_dst(:,ii)= qextTb_dst(:,ii)
       qscatd_dst(:,ii)= qscatTb_dst(:,ii)
           gd_dst(:,ii)= gTb_dst(:,ii)

        qextd_cld(:,:,ii)=   qextTb_cld(:,:,ii)
       qscatd_cld(:,:,ii)=  qscatTb_cld(:,:,ii)
           gd_cld(:,:,ii)=   gTb_cld(:,:,ii)
 enddo


if (mcpu0) then
    print*, ' '
    print*, ' Diagnostic optical properties:  ext scat gfac  '
    do i= 1, nbin_rt
        factor= qscatd_dst(i,3) / qextd_dst(i,3)
        print *, qextd_dst(i,3), factor, gd_dst(i,1), gi_dst(i,5)
    enddo


endif


return
end subroutine initcld


!=====================================================================
!=====================================================================

subroutine setrad(mcpu0)
!     
!     Set up values used by the radiation code, such as the CO2 gas
!     absorption coefficients.  True constants are defined, and the
!     time-independent quantities used by the radiation code are
!     calculated.
!
!
!     AEROSOL RADIATIVE OPTICAL CONSTANTS
!     Values are at the wavelenght interval center
!
!     MIE SCATTERING - Size distribution weighted
!     Qextv    - Extinction efficiency - in the visible.
!     Qscatv   - Scattering efficiency - in the visible.
!     WV       - Single scattering albedo - in the visible.
!     GV       - Asymmetry parameter - in the visible.
!     ggv       - Asymmetry parameter - in the visible.

!     Qexti    - Extinction efficiency - in the infrared.
!     Qscati   - Scattering efficiency - in the infrared.
!     WI       - Single scattering albedo - in the infrared.
!     GI       - Asymmetry parameter - in the infrared.
!     ggi       - Asymmetry parameter - in the infrared.

!----------------------------------------------------------------------C

implicit none

integer :: n, ns
logical :: mcpu0
integer :: nt, nw, ng, err
!----------------------------------------------------------------------

!  Visible dust properties:  M. Wolff Planck-weighted values (T=6000K)
!  Log-normal size distribution:  Reff = 1.5 microns, Veff = 0.5

!     VISIBLE wavelengths
!     Qext - M. Wolff values (order is increasing waveNUMBER)
real*8 :: qev1(L_NSPECTV) =  [ 1.834D0, 2.296D0, 2.672D0,        &
                    2.829D0, 2.698D0, 2.452D0, 2.261D0   ]
!     Qscat - M. Wolff values
real*8 :: qsv1(L_NSPECTV) =  [ 1.695D0, 2.031D0, 2.583D0,        &
                    2.744D0, 2.626D0, 2.225D0, 1.525D0   ]
!     G - M. Wolff values
real*8 :: gv1(L_NSPECTV)  =  [ 0.551D0, 0.640D0, 0.661D0,        &
                    0.678D0, 0.690D0, 0.743D0, 0.868D0   ]

!     And now the INFRARED
!     M. Wolff Planck-weighted values (T=215K)
!     INFRARED wavelengths.  (The order is increasing waveNUMBER.)

real*8, allocatable :: qei1(:), qsi1(:), gi1(:)


! These are John's values for dust optical properties
!      real*8, dimension(L_NSPECTV) :: qev1,qsv1,gv1
!      real*8, dimension(L_NSPECTI) :: qei1,qsi1,gi1
!      real*8, parameter :: sscat=0.9
!      real*8, parameter :: gfac=0.65
!
!      data qev1 / 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5 /
!      data qsv1 / 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2 /
!      data gv1  / .65D0, .65D0, .65D0, .65D0, .65D0, .65D0, .65D0 /
!      data qei1 / .193D0, .867D0, 1.209D0, 2.173D0, 0.638D0 /
!      data qsi1 / .027D0, .319D0, .558D0, 1.136D0, .237D0 /
!      data gi1  / .024D0, .127D0, .288D0, .423D0, .548D0  /

!=======================================================================

!     Set the reference pressure and temperature arrays.  These are
!     the pressures and temperatures at which we have k-coefficients.

allocate(qei1(L_NSPECTI), stat=err)
allocate(qsi1(L_NSPECTI), stat=err)
allocate(gi1(L_NSPECTI), stat=err)

if (ames_15band .OR. use_extended_cor_ks) then
    !     Qext for the IR
    qei1(:) =  [ 0.158D0, 0.467D0, 0.447D0,        &
               0.418D0, 0.414D0, 0.698D0, 1.026D0, 0.765D0   ]
    !     Qsca for M. Wolff      INFRARED wavelengths
    qsi1(:) =  [ 0.021D0, 0.064D0, 0.074D0,        &
               0.114D0, 0.174D0, 0.339D0, 0.271D0, 0.660D0   ]
    !     g for M. Wolff values  INFRARED wavelengths
    gi1(:)  =  [ 0.020D0, 0.052D0, 0.088D0,        &
               0.112D0, 0.134D0, 0.181D0, 0.241D0, 0.375D0   ]
else
    !     Qext for the IR
    qei1(:) =  [ 0.008D0, 0.262D0, 0.491D0,        &
                                          1.017D0, 0.444D0   ]
    !     Qsca for M. Wolff      INFRARED wavelengths
    qsi1(:) =  [ 0.001D0, 0.037D0, 0.122D0,        &
                                          0.351D0, 0.336D0   ]
    !     g for M. Wolff values  INFRARED wavelengths
    gi1(:)  =  [ 0.004D0, 0.030D0, 0.095D0,        &
                                          0.214D0, 0.316D0   ]
endif

if (use_extended_cor_ks) then
    pgasref( 1) = 1.000D-8
    pgasref( 2) = 3.162D-8
    pgasref( 3) = 1.000D-7
    pgasref( 4) = 3.162D-7
    pgasref( 5) = 1.000D-6
    pgasref( 6) = 3.162D-6
    pgasref( 7) = 1.000D-5
    pgasref( 8) = 3.162D-5
    pgasref( 9) = 1.000D-4
    pgasref(10) = 3.162D-4
    pgasref(11) = 1.000D-3
    pgasref(12) = 3.162D-3
    pgasref(13) = 1.000D-2
    pgasref(14) = 3.162D-2
    pgasref(15) = 1.000D-1
    pgasref(16) = 3.162D-1
    pgasref(17) = 1.000D0
    pgasref(18) = 3.162D0
    pgasref(19) = 1.000D+1
    pgasref(20) = 3.162D+1
    pgasref(21) = 1.000D+2
    pgasref(22) = 3.162D+2
    pgasref(23) = 1.000D+3
    pgasref(24) = 3.162D+3
    pgasref(25) = 1.000D+4

    tgasref(1) = 50.000D0
    tgasref(2) = 100.000D0
    tgasref(3) = 150.000D0
    tgasref(4) = 200.000D0
    tgasref(5) = 250.000D0
    tgasref(6) = 300.000D0
    tgasref(7) = 350.000D0
    tgasref(8)  = 400.000D0
    tgasref(9)  = 450.000D0
    tgasref(10) = 500.000D0
    tgasref(11) = 550.000D0
    tgasref(12) = 600.000D0
    tgasref(13) = 650.000D0
    tgasref(14) = 700.000D0
    tgasref(15) = 750.000D0
    tgasref(16) = 800.000D0
 
else if (ames_15band) then
    pgasref( 1) = 1.000D-8
    pgasref( 2) = 3.162D-8
    pgasref( 3) = 1.000D-7
    pgasref( 4) = 3.162D-7
    pgasref( 5) = 1.000D-6
    pgasref( 6) = 3.162D-6
    pgasref( 7) = 1.000D-5
    pgasref( 8) = 3.162D-5
    pgasref( 9) = 1.000D-4
    pgasref(10) = 3.162D-4
    pgasref(11) = 1.000D-3
    pgasref(12) = 3.162D-3
    pgasref(13) = 1.000D-2
    pgasref(14) = 3.162D-2
    pgasref(15) = 1.000D-1
    pgasref(16) = 3.162D-1
    pgasref(17) = 1.000D0
    pgasref(18) = 3.162D0
    pgasref(19) = 1.000D+1
    pgasref(20) = 3.162D+1
    pgasref(21) = 1.000D+2
    pgasref(22) = 3.162D+2
    pgasref(23) = 1.000D+3
    pgasref(24) = 3.162D+3
    pgasref(25) = 1.000D+4

    tgasref( 1) =  50.0D0
    tgasref( 2) =  75.0D0
    tgasref( 3) = 100.0D0
    tgasref( 4) = 125.0D0
    tgasref( 5) = 150.0D0
    tgasref( 6) = 175.0D0
    tgasref( 7) = 200.0D0
    tgasref( 8) = 225.0D0
    tgasref( 9) = 250.0D0
    tgasref(10) = 275.0D0
    tgasref(11) = 300.0D0
    tgasref(12) = 325.0D0
    tgasref(13) = 350.0D0
else
    pgasref( 1) = 1.0E-6
    pgasref( 2) = 1.0E-5
    pgasref( 3) = 1.0E-4
    pgasref( 4) = 1.0E-3
    pgasref( 5) = 1.0E-2
    pgasref( 6) = 1.0E-1
    pgasref( 7) = 1.0
    pgasref( 8) = 1.0E+1
    pgasref( 9) = 1.0E+2
    pgasref(10) = 1.0E+3
    pgasref(11) = 1.0E+4

    tgasref(1)  =  50.0
    tgasref(2)  = 100.0
    tgasref(3)  = 150.0
    tgasref(4)  = 200.0
    tgasref(5)  = 250.0
    tgasref(6)  = 300.0
    tgasref(7)  = 350.0
endif

!     Fill the (VISIBLE) arrays Qextv, Qscatv, WV, GV
!   Original Ames
do n=1,l_nspectv
    qextv(n)  = qev1(n)
    qscatv(n) = qsv1(n)
    if(qscatv(n).ge.qextv(n)) then
        qscatv(n) = 0.99999*qextv(n)
    end if
    wv(n)     = qscatv(n)/qextv(n)
    gv(n)     = gv1(n)
    ggv(n)    = gv1(n)

end do
!  FV3-like-----------------------------------
!      Do N=1,L_NSPECTV
!        Qextv(n) = qev1(n)
!        Qscatv(n) = Qextv(n)*sscat
!      If(Qscatv(n).GE.Qextv(n)) then
!        Qscatv(n) = 0.99999*Qextv(n)
!      END IF
!        WV(n) = sscat
!        GV(n) = gfac
!      END DO
!----------------------------------------------

!     Fill the (INFRARED) arrays Qexti, Qscati, WI, GI
!   Original Ames
do n=1,l_nspecti
    qexti(n)  = qei1(n)
    qscati(n) = qsi1(n)
    if(qscati(n).ge.qexti(n)) then
        qscati(n) = 0.99999*qexti(n)
    end if
    wi(n)     = qscati(n)/qexti(n)
    gi(n)     = gi1(n)
    ggi(n)    = gi1(n)
    if (.not. user_fixed_dust_opts) then
        qxi_read(n)=qei1(n)
        qsi_read(n)=qsi1(n)
        if(qscati(n).ge.qexti(n)) then
            qsi_read(n) = 0.99999*qexti(n)
        end if
        gi_read(n) = gi1(n)
    endif

end do
!  FV3-like-------------------------------------
!      Do N=1,L_NSPECTI
!        Qexti(n) = qei1(n)
!        Qscati(n)= qsi1(n)*Qexti(n)/qev1(1)
!        If(Qscati(n).ge.Qexti(n)) then
!          Qscati(n) = 0.99999*Qexti(n)
!        END IF
!        WI(n) = Qscati(n)/Qexti(n)
!        GI(n) = gi1(n)
!      END DO

!     Interpolate CO2 k coefficients to the finer pressure grid.

if (ames_15band .OR. use_extended_cor_ks) then
    call initinterp(PGASREF,mcpu0)
else
    call laginterp(PGASREF,PFGASREF)
endif

deallocate(qei1)
deallocate(qsi1)
deallocate(gi1)

return
end subroutine setrad


!=====================================================================
!=====================================================================

subroutine setspi()
!
!     Set up the spectral intervals in the infrared.
!
!
!**********************************************************************C


implicit none

!     BWNI - Bin wavenumber of the edges of the IR spectral bins
!     units are inverse centimeters.  Dimension needs to be changed
!     if the number of IR bins changes.

REAL*8, allocatable  :: BWNI(:)

real*8  :: a, b, ans, y, bpa, bma, T
real*8  :: wn1, wn2
integer :: n, nw, nt, m, err

!  C1 and C2 values from Goody and Yung (2nd edition)  MKS units
!  These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4

real*8 :: c1 = 3.741832D-16      ! W m^-2
real*8 :: c2 = 1.438786D-2       ! m K

real*8 :: x(12) = [ -0.981560634246719D0,  -0.904117256370475D0, &
              -0.769902674194305D0,  -0.587317954286617D0, &
              -0.367831498998180D0,  -0.125233408511469D0, &
               0.125233408511469D0,   0.367831498998180D0, &
               0.587317954286617D0,   0.769902674194305D0, &
               0.904117256370475D0,   0.981560634246719D0 ]

real*8 :: w(12) = [  0.047175336386512D0,   0.106939325995318D0, &
               0.160078328543346D0,   0.203167426723066D0, &
               0.233492536538355D0,   0.249147045813403D0, &
               0.249147045813403D0,   0.233492536538355D0, &
               0.203167426723066D0,   0.160078328543346D0, &
               0.106939325995318D0,   0.047175336386512D0 ]

!======================================================================C

allocate(bwni(l_nspecti+1),stat=err)
!     bin wavenumber - wavenumber [cm^(-1)] at the edges of the ir
!     spectral bins.

if (ames_15band .OR. use_extended_cor_ks) then
    bwni( 1) =   10.000d0       ! 1000.0 microns
    bwni( 2) =  300.300d0       !   33.3 microns
    bwni( 3) =  549.451d0       !   18.2 microns
    bwni( 4) =  636.943d0       !   15.7 microns
    bwni( 5) =  704.225d0       !   14.2 microns
    bwni( 6) =  775.194d0       !   12.9 microns
    bwni( 7) =  970.874d0       !   10.3 microns
    bwni( 8) = 1282.051d0       !    7.8 microns
    bwni( 9) = 2222.222d0       !    4.5 microns
else
    bwni( 1) =   10.000d0       ! 1000.0 microns
    bwni( 2) =  166.667d0       !   60.0 microns
    bwni( 3) =  416.667d0       !   24.0 microns
    bwni( 4) =  833.333d0       !   12.0 microns
    bwni( 5) = 1250.000d0       !    8.0 microns
    bwni( 6) = 2222.222d0       !    4.5 microns
endif

!     Set up mean wavenumbers and wavenumber deltas.  Units of
!     wavenumbers is cm^(-1); units of wavelengths is microns.

do m=1,l_nspecti
    wnoi(m)  = 0.5*(bwni(m+1)+bwni(m))
    dwni(m)  = bwni(m+1)-bwni(m)
    wavei(m) = 1.0e+4/wnoi(m)
end do

!  for each ir wavelength interval, compute the integral of b(t), the
!  planck function, divided by the wavelength interval, in cm-1.  the
!  integration is in mks units, the final answer is the same as the
!  original planck.f; w m^-2 wavenumber^-1, where wavenumber is in cm^-1.

do nw=1,l_nspecti
    a = 1.0d-2/bwni(nw+1)
    b = 1.0d-2/bwni(nw)
    bpa = (b+a)/2.0
    bma = (b-a)/2.0
    do nt=500,9000   !3500 is in john's version
        t   = dble(nt)/1.0d+1
        ans = 0.0d0
        do m=1,12
            y    = bma*x(m)+bpa
            ans  = ans + w(m)*c1/(y**5*(exp(c2/(y*t))-1.0d0))
        end do
        planckir(nw,nt-499) = ans*bma/(pi*dwni(nw))
    end do
end do

deallocate(bwni)

return
end subroutine setspi


!=====================================================================
!=====================================================================

subroutine setspv(mcpu0)
!
!     Set up the spectral intervals in the visible (solar).
!
!**********************************************************************C


implicit none

!     bwnv - bin wavenumber of the edges of the visible spectral bins
!     units are inverse centimeters.  dimension needs to be changed
!     if the number of visible bins changes.

real*8  :: sum, wl
integer :: n, m
logical :: mcpu0

!     p0      - rayleigh scattering reference pressure in pascals.
!     grav    - acceleration due to gravity (g) - mks

real*8  :: p0 = 9.423d+6

!     bin wavenumber - wavenumber [cm^(-1)] at the edges of the visible
!     spectral bins.  go from smaller to larger wavenumbers, the same as
!     in the ir.

!     2222.22d0    ->   4.50 microns
!     3087.37d0    ->   3.24 microns
!     4030.63d0    ->   2.48 microns
!     5370.57d0    ->   1.86 microns
!     7651.11d0    ->   1.31 microns
!     12500.00d0   ->   0.80 microns
!     25000.00d0   ->   0.40 microns
!     41666.67d0   ->   0.24 microns

real*8 :: bwnv(l_nspectv+1) = [ 2222.22d0, 3087.37d0, 4030.63d0, &
                             5370.57d0, 7651.11d0, 12500.00d0, &
                             25000.00d0, 41666.67d0 ]

!     solar flux within each spectral interval, at 1au (w/m^2)
!     sum equals 1356 w/m^2 (values from wehrli, 1985)

real*8 :: solar(l_nspectv) = [ 12.7, 24.2, 54.6, 145.9, 354.9,   &
                             657.5, 106.3 ]


!======================================================================C

!     Set up mean wavenumbers and wavenumber deltas.  Units of
!     wavenumbers is cm^(-1); units of wavelengths is microns.

do m=1,l_nspectv
    wnov(m)  = 0.5*(bwnv(m+1)+bwnv(m))
    dwnv(m)  = bwnv(m+1)-bwnv(m)
    wavev(m) = 1.0e+4/wnov(m)
end do

!     sum the solar flux, and write out the result.

sum = 0.0
do n=1,l_nspectv
    solarf(n) = solar(n)
    sum       = sum+solarf(n)
end do

if (mcpu0) write(6,'("solar flux at 1au = ",f7.2," w/m^2")') sum

!     set up the wavelength independent part of rayleigh scattering.
!     the pressure dependent part will be computed elsewhere (optcv).
!     wavev is in microns.  there is no rayleigh scattering in the ir.

do n=1,l_nspectv
    wl        = wavev(n)
    tauray(n) = (8.7/grav)*(1.527*(1.0+0.013/wl**2)/wl**4)*        &
             scalep/p0
end do

return
end subroutine setspv

!=====================================================================
!=====================================================================

subroutine laginterp(pgref,pint)
!  Lagrange interpolation (linear in log pressure) of the CO2
!  k-coefficients in the pressure domain.  Subsequent use of these
!  values will use a simple linear interpolation in pressure.

implicit none

integer :: n, nt, np, nh, ng, nw, m, i
real*8  :: co2i8(l_ntref,l_npref,l_refh2o,l_nspecti,l_ngauss)
real*8  :: co2v8(l_ntref,l_npref,l_refh2o,l_nspectv,l_ngauss)
real*8  :: pgref(l_npref)


real*8  :: x, xi(4), yi(4), ans
real*8  :: pint(l_pint), pref(l_npref), p

real*8  :: pin(l_pint) =  [                                      &
                   -6.0d0, -5.8d0, -5.6d0, -5.4d0, -5.2d0,     &
                   -5.0d0, -4.8d0, -4.6d0, -4.4d0, -4.2d0,     &
                   -4.0d0, -3.8d0, -3.6d0, -3.4d0, -3.2d0,     &
                   -3.0d0, -2.8d0, -2.6d0, -2.4d0, -2.2d0,     &
                   -2.0d0, -1.8d0, -1.6d0, -1.4d0, -1.2d0,     &
                   -1.0d0, -0.8d0, -0.6d0, -0.4d0, -0.2d0,     &
                    0.0d0,  0.2d0,  0.4d0,  0.6d0,  0.8d0,     &
                    1.0d0,  1.2d0,  1.4d0,  1.6d0,  1.8d0,     &
                    2.0d0,  2.2d0,  2.4d0,  2.6d0,  2.8d0,     &
                    3.0d0,  3.2d0,  3.4d0,  3.6d0,  3.8d0,     &
                    4.0d0                                  ]

!======================================================================!

!  Fill pint for output from this subroutine

do n=1,l_pint
    pint(n) = pin(n)
end do

!  take log of the reference pressures

do n=1,l_npref
    pref(n) = log10(pgref(n))
end do

!     get co2 k coefficients

open(20,file=trim(rtdata_path)//'CO2H2O_V_12_95_INTEL',   &
      form='unformatted')
read(20) co2v8
read(20) fzerov
close(20)

open(20,file=trim(rtdata_path)//'CO2H2O_IR_12_95_INTEL',  &
      form='unformatted')
read(20) co2i8
read(20) fzeroi
close(20)

!  Take Log10 of the values - we interpolate the log10 of the values,
!  not the values themselves.   Smallest value is 1.0E-200.

do nt=1,l_ntref
    do np=1,l_npref
        do nh=1,l_refh2o
            do ng = 1,l_ngauss

                do nw=1,l_nspectv
                    if(co2v8(nt,np,nh,nw,ng).gt.1.0d-200) then
                        co2v8(nt,np,nh,nw,ng) = log10(co2v8(nt,np,nh,nw,ng))
                    else
                        co2v8(nt,np,nh,nw,ng) = -200.0
                    end if
                end do

                do nw=1,l_nspecti
                    if(co2i8(nt,np,nh,nw,ng).gt.1.0d-200) then
                        co2i8(nt,np,nh,nw,ng) = log10(co2i8(nt,np,nh,nw,ng))
                    else
                        co2i8(nt,np,nh,nw,ng) = -200.0
                    end if
                end do

            end do
        end do
    end do
end do

!  Interpolate the values:  first the IR

do nt=1,l_ntref
    do nh=1,l_refh2o
        do nw=1,l_nspecti
            do ng=1,l_ngauss

                !  first, the initial interval (p=1e-6 to 1e-5)

                n = 1
                do m=1,5
                    x     = pint(m)
                    xi(1) = pref(n)
                    xi(2) = pref(n+1)
                    xi(3) = pref(n+2)
                    xi(4) = pref(n+3)
                    yi(1) = co2i8(nt,n,nh,nw,ng)
                    yi(2) = co2i8(nt,n+1,nh,nw,ng)
                    yi(3) = co2i8(nt,n+2,nh,nw,ng)
                    yi(4) = co2i8(nt,n+3,nh,nw,ng)
                    call lagrange(x,xi,yi,ans)
                    co2i(nt,m,nh,nw,ng) = 10.0**ans
                end do

                do n=2,l_npref-2
                    do m=1,5
                        i     = (n-1)*5+m
                        x     = pint(i)
                        xi(1) = pref(n-1)
                        xi(2) = pref(n)
                        xi(3) = pref(n+1)
                        xi(4) = pref(n+2)
                        yi(1) = co2i8(nt,n-1,nh,nw,ng)
                        yi(2) = co2i8(nt,n,nh,nw,ng)
                        yi(3) = co2i8(nt,n+1,nh,nw,ng)
                        yi(4) = co2i8(nt,n+2,nh,nw,ng)
                        call lagrange(x,xi,yi,ans)
                        co2i(nt,i,nh,nw,ng) = 10.0**ans
                    end do
                end do

                !  Now, get the last interval (P=1e+3 to 1e+4)

                n = l_npref-1

                do m=1,5
                    i     = (n-1)*5+m
                    x     = pint(i)
                    xi(1) = pref(n-2)
                    xi(2) = pref(n-1)
                    xi(3) = pref(n)
                    xi(4) = pref(n+1)
                    yi(1) = co2i8(nt,n-2,nh,nw,ng)
                    yi(2) = co2i8(nt,n-1,nh,nw,ng)
                    yi(3) = co2i8(nt,n,nh,nw,ng)
                    yi(4) = co2i8(nt,n+1,nh,nw,ng)
                    call lagrange(x,xi,yi,ans)
                    co2i(nt,i,nh,nw,ng) = 10.0**ans
                end do

                !  fill the last pressure point

                co2i(nt,l_pint,nh,nw,ng) = 10.0**co2i8(nt,l_npref,nh,nw,ng)

            end do
        end do
    end do
end do

!  Interpolate the values:  now the visible

do nt=1,l_ntref
    do nh=1,l_refh2o
        do nw=1,l_nspectv
            do ng=1,l_ngauss

                !  first, the initial interval (p=1e-6 to 1e-5)

                n = 1
                do m=1,5
                    x     = pint(m)
                    xi(1) = pref(n)
                    xi(2) = pref(n+1)
                    xi(3) = pref(n+2)
                    xi(4) = pref(n+3)
                    yi(1) = co2v8(nt,n,nh,nw,ng)
                    yi(2) = co2v8(nt,n+1,nh,nw,ng)
                    yi(3) = co2v8(nt,n+2,nh,nw,ng)
                    yi(4) = co2v8(nt,n+3,nh,nw,ng)
                    call lagrange(x,xi,yi,ans)
                    co2v(nt,m,nh,nw,ng) = 10.0**ans
                end do

                do n=2,l_npref-2
                    do m=1,5
                        i     = (n-1)*5+m
                        x     = pint(i)
                        xi(1) = pref(n-1)
                        xi(2) = pref(n)
                        xi(3) = pref(n+1)
                        xi(4) = pref(n+2)
                        yi(1) = co2v8(nt,n-1,nh,nw,ng)
                        yi(2) = co2v8(nt,n,nh,nw,ng)
                        yi(3) = co2v8(nt,n+1,nh,nw,ng)
                        yi(4) = co2v8(nt,n+2,nh,nw,ng)
                        call lagrange(x,xi,yi,ans)
                        co2v(nt,i,nh,nw,ng) = 10.0**ans
                    end do
                end do

                !  Now, get the last interval (P=1e+3 to 1e+4)

                n = l_npref-1

                do m=1,5
                    i     = (n-1)*5+m
                    x     = pint(i)
                    xi(1) = pref(n-2)
                    xi(2) = pref(n-1)
                    xi(3) = pref(n)
                    xi(4) = pref(n+1)
                    yi(1) = co2v8(nt,n-2,nh,nw,ng)
                    yi(2) = co2v8(nt,n-1,nh,nw,ng)
                    yi(3) = co2v8(nt,n,nh,nw,ng)
                    yi(4) = co2v8(nt,n+1,nh,nw,ng)
                    call lagrange(x,xi,yi,ans)
                    co2v(nt,i,nh,nw,ng) = 10.0**ans
                end do

                !  fill the last pressure point

                co2v(nt,l_pint,nh,nw,ng) = 10.0**co2v8(nt,l_npref,nh,nw,ng)

            end do
        end do
    end do
end do

!for boxintep use log10 values
if (use_boxinterp12) then
    co2i=log10(co2i)
    co2v=log10(co2v)
endif

return
end subroutine laginterp


!=====================================================================
!=====================================================================

subroutine lagrange(x, xi, yi, ans)
!  Lagrange interpolation - Polynomial interpolation at point x
!  xi(1) <= x <= xi(4).  Yi(n) is the functional value at XI(n).

implicit none

real*8 :: x, xi(4), yi(4), ans
real*8 :: fm1, fm2, fm3, fm4

!======================================================================!

fm1 = x - xi(1)
fm2 = x - xi(2)
fm3 = x - xi(3)
fm4 = x - xi(4)

!  get the "answer" at the requested x

ans = fm2*fm3*fm4*yi(1)/                                         &
          ((xi(1)-xi(2))*(xi(1)-xi(3))*(xi(1)-xi(4)))  +   &
    fm1*fm3*fm4*yi(2)/                                         &
              ((xi(2)-xi(1))*(xi(2)-xi(3))*(xi(2)-xi(4)))  +   &
    fm1*fm2*fm4*yi(3)/                                         &
              ((xi(3)-xi(1))*(xi(3)-xi(2))*(xi(3)-xi(4)))  +   &
    fm1*fm2*fm3*yi(4)/                                         &
              ((xi(4)-xi(1))*(xi(4)-xi(2))*(xi(4)-xi(3)))

return
end subroutine lagrange

!=====================================================================
!=====================================================================

subroutine dsolflux(sol,ubar0,detau,directsol)
!
!  calculate the direct surface solar flux (solar flux at the surface due
!  to the direct solar beam.
!

implicit none

real*8 ::  sol(l_nspectv), ubar0, directsol
real*8 :: factor
real*8 :: detau(l_nspectv,l_ngauss)

integer nw, ng

!======================================================================c

directsol = 0.0d0

do nw=1,l_nspectv
    factor = ubar0*sol(nw)

    do ng=1,l_ngauss-1
        if(detau(nw,ng) .le. 5.0) then
            directsol = directsol + factor*   &
                 exp(-detau(nw,ng)/ubar0)*  &
                 gweight(ng)*(1.0d0-fzerov(nw))
        end if
    end do
    ng        = l_ngauss
    if(detau(nw,ng) .le. 5.0) then
        directsol = directsol + factor*    &
           exp(-detau(nw,ng)/ubar0)*fzerov(nw)
    endif
end do

return
end subroutine dsolflux

!======================================================================C
!======================================================================C

subroutine boxinterp(ft1p1,ft2p1,ft1p2,ft2p2,tref1,tref2,pref1,  &
         pref2,tmid,pmid,tinter,pinter,ans)
!
!   Calculate 2d interpolation from pressure and temperature
!----------------------------------------------------------------------!
!                   T2 .FT2P1                    .FT2P2
!
!
!                   T1 .FT1P1                    .FT1P2
!                      P1                        P2
!----------------------------------------------------------------------!

implicit none

real*8 :: ft1p1, ft2p1, ft1p2, ft2p2, tref1, tref2
real*8 :: pref1, pref2, tmid, pmid, ans1, ans2, ans
logical :: tinter, pinter

!======================================================================C

if(.not.tinter .and. .not.pinter) then
    ans = ft1p1
elseif(.not.tinter .and. pinter) then
    ans = ft1p1 + (ft1p2 - ft1p1)*(pmid - pref1)/(pref2 - pref1)
elseif(tinter .and. .not.pinter) then
    ans = ft1p1 + (ft2p1 - ft1p1)*(tmid - tref1)/(tref2 - tref1)
elseif(tinter .and. pinter) then
    ans1 = ft1p1 + (ft2p1 - ft1p1)*(tmid - tref1)/(tref2 - tref1)
    ans2 = ft1p2 + (ft2p2 - ft1p2)*(tmid - tref1)/(tref2 - tref1)
    ans  = ans1 + (ans2 - ans1)*(pmid - pref1)/(pref2 - pref1)
endif

return
end subroutine boxinterp

!======================================================================C
!======================================================================C

subroutine initinterp(pgref,mcpu0)
!  Set up for interpolation (linear in log pressure) of the CO2 
!  k-coefficients in the pressure domain.  Subsequent use of these
!  values will use a simple linear interpolation in pressure.


implicit none

integer :: n, nt, np, nh, ng, nw, m, i, k, mn
real*8  :: co2i8(l_ntref,l_npref,l_refh2o,l_nspecti,l_ngauss)
real*8  :: co2v8(l_ntref,l_npref,l_refh2o,l_nspectv,l_ngauss)
real*8  :: pgref(l_npref)

logical :: mcpu0

real*8  :: x, xi(3), yi(3), ans
real*8  :: p

!======================================================================!

!  Take log of the reference pressures

do n=1,L_NPREF
    pgref(n) = LOG10(PGREF(n))
!    if(mcpu0) print*,'n, pgref(n)',n, pgref(n)
end do

!     Get CO2 k coefficients

!!   For 16 gauss points
!      open(20,file=trim(rtdata_path)//'CO2H2O_V_2013_16', status='old',      &
!              form='unformatted')
!      read(20) co2v8
!      read(20) fzerov
!      close(20)
!
!      open(20,file=trim(rtdata_path)//'CO2H2O_IR_2013_16', status='old',    &
!              form='unformatted')
!      read(20) co2i8
!      read(20) fzeroi
!      close(20)
!

if (use_extended_cor_ks) then
    open(20,file=trim(rtdata_path)//'CO2H2O_V_15B_800K_v4', &
          status='old', form='unformatted')
    read(20) co2v8
    read(20) fzerov
    close(20)

    open(20,file=trim(rtdata_path)//'CO2H2O_IR_15B_800K_v4', &
          status='old', form='unformatted')
    read(20) co2i8
    read(20) fzeroi
    close(20)

     
else

!!   For 32 gauss points
    open(20,file=trim(rtdata_path)//'CO2H2O_V_2013_32', status='old',      &
      form='unformatted')
    read(20) co2v8
    read(20) fzerov
    close(20)

    open(20,file=trim(rtdata_path)//'CO2H2O_IR_2013_32', status='old',    &
      form='unformatted')
    read(20) co2i8
    read(20) fzeroi
    close(20)
endif

!  Take Log10 of the values - we interpolate the log10 of the values,
!  not the values themselves.   Smallest value is 1.0E-200.

do nt=1,l_ntref
    do np=1,l_npref
        do nh=1,l_refh2o
            do ng = 1,l_ngauss

                do nw=1,l_nspectv
                    if(co2v8(nt,np,nh,nw,ng).gt.1.0d-200) then
                        co2v(nt,np,nh,nw,ng) = log10(co2v8(nt,np,nh,nw,ng))
                    else
                        co2v(nt,np,nh,nw,ng) = -200.0
                    end if
                end do

                do nw=1,l_nspecti
                    if(co2i8(nt,np,nh,nw,ng).gt.1.0d-200) then
                        co2i(nt,np,nh,nw,ng) = log10(co2i8(nt,np,nh,nw,ng))
                    else
                        co2i(nt,np,nh,nw,ng) = -200.0
                    end if
                end do

            end do
        end do
    end do
end do


! "zero" out 17th (or 33rd) Gauss point

      do nt=1,L_NTREF
        do nh=1,L_REFH2O
          do nw=1,L_NSPECTV
            do np=1,L_NPREF
              co2v(nt,np,nh,nw,L_NGAUSS) = -200.0D0
            end do
          end do
          do nw=1,L_NSPECTI
            do np=1,L_NPREF
              co2i(nt,np,nh,nw,L_NGAUSS) = -200.0D0
            end do
          end do
        end do
      end do


      return
      end subroutine initinterp

!======================================================================C
!======================================================================C

subroutine tpindex15(pw,tw,qh2o,pref,tref,MT,MP,MW,pinter, &
                         tinter)
!    For the 15 band model:
!    Get the TI, UI values for a 2-dimensional interpolation
!    based on the following (The interpolation is done in interpco2):
!    Interpolate the CO2 K-coefficients to the current P,T values.
!    The CO2 coefficients are given on a P,T grid:
!    P = {1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 1E+1, 1E+2, 1E+3, 1E+4},
!    T = {50, 100, 150, 200, 250, 300, 350}.
!
!----------------------------------------------------------------------C

implicit none

integer :: mt(2), mp(2), mw, n
real*8  :: pw, tw, qh2o
real*8, intent(in)  :: pref(:), tref(:)
logical :: pinter, tinter

!======================================================================C

!     Get the upper and lower Temperature-grid indicies that bound the
!     requested temperature.  If the requested temperature is outside
!     the T-grid, set up to extrapolate from the appropriate end.

if(tw.le.tref(1)) then
    mt(1) = 1
    mt(2) = mt(1)
    tinter = .false.
elseif(tw.ge.tref(l_ntref)) then
    mt(1) = l_ntref
    mt(2) = mt(1)
    tinter = .false.
else
    do n=1,l_ntref-1
        if(tw.ge.tref(n) .and. tw.lt.tref(n+1)) then
            mt(1) = n
            mt(2) = n + 1
            tinter = .true.
            exit
        endif
    end do
end if

!     Get the upper and lower Pressure-grid indicies that bound the
!     requested pressure.  If the requested pressure is outside
!     the P-grid, set up to extrapolate from the appropiate end.
if (use_boxinterp12) then
    if(pw.le.pref(1)) then
        mp(1) = 1
        mp(2) = mp(1)
        pinter = .false.
    elseif(pw.ge.pref(l_pint)) then
        mp(1) = l_pint
        mp(2) = mp(1)
        pinter = .false.
    else
        do n=1,l_pint-1
            if(pw.ge.pref(n) .and. pw.lt.pref(n+1)) then
                mp(1) = n
                mp(2) = n + 1
                pinter = .true.
                exit
            endif
        end do
    end if
else
    if(pw.le.pref(1)) then
        mp(1) = 1
        mp(2) = mp(1)
        pinter = .false.
    elseif(pw.ge.pref(l_npref)) then
        mp(1) = l_npref
        mp(2) = mp(1)
        pinter = .false.
    else
        do n=1,l_npref-1
            if(pw.ge.pref(n) .and. pw.lt.pref(n+1)) then
                mp(1) = n
                mp(2) = n + 1
                pinter = .true.
                exit
            endif
        end do
    end if
endif

!  Get the indicies for water abundance.  There are 10 sets of
!  k-coefficients with differing amounts of water vs. CO2.

if(qh2o.ge.wrefh2o(l_refh2o)) then
    mw = l_refh2o
elseif(qh2o.le.wrefh2o(1)) then
    mw = 1
else
    do n=1,l_refh2o-1
        if(qh2o.ge.wrefh2o(n) .and. qh2o.lt.wrefh2o(n+1)) then
            mw = n
            exit
        end if
    end do
end if

return
end subroutine tpindex15

!======================================================================C
!======================================================================C

!subroutine nltecool(nlayer,player,tlayer,dt)

!***********************************************************************
!***********************************************************************

subroutine interp1(escout1,p1,nlayer1,escin1,pin1,np)
!
! subroutine to perform linear interpolation in pressure from 1D profile
! escin(nl) sampled on pressure grid pin(nl) to profile
! escout(nlayer) on pressure grid p(nlayer).
!
real*8 escout1(nlayer1),p1(nlayer1)
real*8 escin1(np),pin1(np),wm1,wp1
integer nl1,nlayer1,n11,n,nm1,np1,np

nl1 = np
nm1 = 1
np1 = 1
wm1 = 1.
wp1 = 0.


do n11=1,nlayer1
    if(p1(n11) .gt. 800.0 .or. p1(n11) .lt. 5.0e-9) then
        escout1(n11) = 0.0
    else
        do n = 1,nl1-1
            if (p1(n11).le.pin1(n).and.p1(n11).ge.pin1(n+1)) then
                nm1=n
                np1=n+1
                wm1=abs(pin1(np1)-p1(n11))/(pin1(nm1)-pin1(np1))
                wp1=1.0 - wm1
            endif
        enddo
        escout1(n11) = escin1(nm1)*wm1 + escin1(np1)*wp1
    endif
enddo
return
end

!***********************************************************************
!***********************************************************************

subroutine interp2(escout2,p2,nlayer2,escin2,pin2,np)
!
! subroutine to perform linear interpolation in pressure from 1D profile
! escin(nl) sampled on pressure grid pin(nl) to profile
! escout(nlayer) on pressure grid p(nlayer).
!
real*8 escout2(nlayer2),p2(nlayer2)
real*8 escin2(np),pin2(np),wm2,wp2
integer nl2,nlayer2,n12,n,nm2,np2,np


nl2 = np
nm2 = 1
np2 = 1
wm2 = 1.
wp2 = 0.

do n12=1,nlayer2
    if(p2(n12) .gt. 800.0 .or. p2(n12) .lt. 5.0e-9) then
        escout2(n12) = 0.0
    else
        do n = 1,nl2-1
            if (p2(n12).le.pin2(n).and.p2(n12).ge.pin2(n+1)) then
                nm2=n
                np2=n+1
                wm2=abs(pin2(np2)-p2(n12))/(pin2(nm2)-pin2(np2))
                wp2=1.0 - wm2
            endif
        enddo

        escout2(n12) = escin2(nm2)*wm2 + escin2(np2)*wp2
    endif

enddo

return

end

!***********************************************************************
!***********************************************************************

subroutine interp3(esco1,esco2,esco3,p3,nlayer3,  &
   esci1,esci2,esci3,pin3,np)
!
! subroutine to perform 3 simultaneous linear interpolations in pressure frm
! 1D profiles esci1-3(nl) sampled on pressure grid pin(nl)to 1D profiles
! esco1-3(nlayer) on pressure grid p(ngrid,nlayer).
!
real*8 esco1(nlayer3),esco2(nlayer3),esco3(nlayer3),p3(nlayer3)
real*8 esci1(np),    esci2(np),    esci3(np), pin3(np),wm3,wp3
integer nl3,nlayer3,n13,n,nm3,np3,np

nl3 = np
nm3 = 1
np3 = 1
wm3 = 1.
wp3 = 0.


do n13=1,nlayer3
    if(p3(n13) .gt. 800.0 .or. p3(n13) .lt. 5.0e-9) then
        esco1(n13)=0.0
        esco2(n13)=0.0
        esco3(n13)=0.0
    else
        do n = 1,nl3-1
            if (p3(n13).le.pin3(n).and.p3(n13).ge.pin3(n+1)) then
                nm3=n
                np3=n+1
                wm3=abs(pin3(np3)-p3(n13))/(pin3(nm3)-pin3(np3))
                wp3=1.0 - wm3
            endif
        enddo
        esco1(n13) = esci1(nm3)*wm3 + esci1(np3)*wp3
        esco2(n13) = esci2(nm3)*wm3 + esci2(np3)*wp3
        esco3(n13) = esci3(nm3)*wm3 + esci3(np3)*wp3
    endif

enddo

return
end

!======================================================================C
!======================================================================C
function nml_switch(int_nml,int_opt)
!returns nml_switch = .true. if a digit of int_nml matches int_opt
implicit none
!---input----
integer, intent(in) :: int_nml !integer from namelist, 5 digit max
integer, intent(in) :: int_opt !integer option ranging from 0->9
!---output---
logical :: nml_switch
!---work var-ables---
character(5) :: txt_nml
character(1) :: txt_opt
integer i
!initialization:
nml_switch = .false.

write(txt_nml,'(i5)') int_nml
write(txt_opt,'(i1)') int_opt
do i =1,5
    if (txt_nml(i:i) .eq. txt_opt) then
        nml_switch = .true.
    endif
enddo

end function nml_switch

!======================================================================C
!======================================================================C

subroutine calc_opts(nlon,nbin,qext_in,qscat_in,g_in, &
                    surf,qx_out,qs_out,g_out,  &
                    surf_split,qx_split,qs_split,g_split)

integer, intent(in) :: nlon, nbin
real*8, intent(in), dimension(nbin) :: surf
real*4, intent(in), dimension(nbin,nlon) :: qext_in, qscat_in, g_in
real*8, intent(out), dimension(nlon) :: qx_out, qs_out, g_out
real*8, intent(in), optional, dimension(ndust_mass,nbin) :: surf_split
real*8, intent(out), optional, dimension(ndust_mass,nlon) :: qx_split,qs_split,g_split

integer :: iwav, i, nd, ndx
logical :: do_split

do_split = present(qx_split) .or. present(qs_split) .or. present(g_split)

do iwav = 1, nlon
    qx_out(iwav) = 0.
    qs_out(iwav) = 0.
     g_out(iwav)= 0.0
    do i = 1, nbin
        qx_out(iwav) = qx_out(iwav)+ surf(i) * qext_in(i,iwav)
        qs_out(iwav) = qs_out(iwav)+ surf(i) * qscat_in(i,iwav)
        g_out(iwav)  = g_out(iwav) + surf(i) * g_in(i,iwav)
        if (do_split) then
            do nd = 1, ndust_mass
                ndx = dust_mass_indx(nd)
                if (ndx .gt. nt_nontag) exit
                if (present(qx_split)) qx_split(nd,iwav) = qx_split(nd,iwav)+ surf_split(nd,i) * qext_in(i,iwav)
                if (present(qs_split)) qs_split(nd,iwav) = qs_split(nd,iwav)+ surf_split(nd,i) * qscat_in(i,iwav)
                if (present(g_split)) g_split(nd,iwav) = g_split(nd,iwav)+ surf_split(nd,i) * g_in(i,iwav)
            enddo
        endif
    enddo
    qs_out(iwav) = min( qs_out(iwav) , 0.99999*qx_out(iwav) )
    if (present(qs_split)) then
       do nd = 1, ndust_mass
           qs_split(nd,iwav) = min( qs_split(nd,iwav), 0.99999*qx_split(nd,iwav) )
       enddo
    endif
enddo


end subroutine calc_opts

end module rtmod_mgcm
