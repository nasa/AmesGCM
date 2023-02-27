module initracer_mod
!
!  module to initialize tracers
!
	
use constants_mod, only: grav,cp=>cp_air,rgas=>rdgas,pi,kboltz
use field_manager_mod, only: MODEL_ATMOS, parse, find_field_index, num_tags, find_tagging_index
use tracer_manager_mod, only: query_method, get_tracer_index,  &
                       get_number_tracers, get_tracer_names
use fms_mod, only: error_mesg, FATAL, file_exist,                      &
           open_namelist_file, check_nml_error,                &
           mpp_pe, mpp_root_pe, close_file,                    &
           write_version_number, stdlog,                       &
           uppercase, read_data, write_data, field_size
use aerosol_util_mod, only: Reff_backgd, do_15band

implicit none

integer,public :: iMa_dt,iNb_dt,iMa_cor
integer,public :: iMa_vap,iMa_cld,iNb_cld
integer, public, parameter :: nbin  = 4                  ! number of bin tracers
integer, public, parameter :: naer  = 5                  ! number of moment tracers
integer, public :: ndust_mass                            ! number of dust mass aerosols
integer, public :: nice_mass                             ! number of ice mass aerosols
integer, public :: ntrace_micro                          ! number of dust mass aerosols
integer, public :: ntrace_gas                            ! number of gas tracer ouside co2 (argon...)
integer, public :: nt_nontag                             ! number of non tag tracers
real*8, public, parameter :: dpden_dt = 2.5E+3, &        ! density of dust [kg/m^3]
                dpden_ice = 917.                         ! density of ice [kg/m^3]
real*8, public, parameter :: athird = 1./3.              ! 1/3
real*8, public :: aerdens(naer), &                       ! density of all aerosols [kg/m^3]
                scale_dt                                 ! scale dust
real*8, public :: dev_dt, &                              ! standard deviation of dust distribution
                dev_ice, &                               ! standard deviation of ice distribution
                aerad(nbin), &                           ! radii of size distribution bins
                rb(nbin+1)                               ! radius boundaries
real*8, public, parameter :: mh2o = 18.01e-3, &          ! molar weight of water [kg]
                rgp = 8.3143, &                          ! perfect gas constant
                vo1 = mh2o / dpden_ice , &               ! volume of water molecule [m^3]
                nav = 6.02e23, &                         ! Avogadro number
!                kbz = rgp / nav, &                       ! Boltzmann constant
                mteta = 0.965                            ! Contact parameter ( m=cos(theta) orig .975 )
real*8, public, parameter :: desorp = 0.288e-19, &       ! Activation energy for desorption of water on a dust-like substrate (J/molecule)
                surfdif = desorp / 10., &                ! Estimated activation energy for surface diffusion of water molecules (J/molecule)
                nus = 1.e+13, &                          ! Jump frequency of a water molecule (s-1)
                m0 = mh2o / nav, &                       ! Weight of a water molecule (kg)
                a0 = 4.52E-10, &
                d_nuc = a0 / sqrt(2.)
real*8, public, parameter :: mco2 = 44.01e-3, &          ! Molecular weight of co2  (kg mol-1) 
                mco20 = mco2 / nav, &                    ! Weight of a CO2 molecule (kg)
                dynvis0 = 14.8e-6                        ! dynamic viscosity CO2 at 293.15 K (Pa.s)
                
real*8, public, parameter :: small_mass = 1.e-15, &      ! Small dust mass mixing ratio
                small_num = 1.e-3                        ! Small dust number mixing ratio


! Specific coagulation
real*8, public :: r1_coag = 0.01e-6
real*8, public :: rn_coag = 40.e-6
integer, public, parameter :: nres_coag = 25 
real*8, public :: vrat_coag,rads_coag(nres_coag),vols_coag(nres_coag),deltar_coag(nres_coag)

public :: initmicro
public :: dust_mass_indx,dust_nb_indx,dust_cor_indx,micro_indx,gas_indx,sedim_indx
public :: ice_mass_indx,ice_nb_indx,vapor_indx
public :: reff_mom, tau_frac, qext_dt


integer,  dimension(:),  allocatable  ::   dust_mass_indx,dust_nb_indx,dust_cor_indx,micro_indx,gas_indx,sedim_indx
integer,  dimension(:),  allocatable  ::   ice_mass_indx,ice_nb_indx,vapor_indx
real*8,   dimension(:),  allocatable  ::   stdv
real*8,   dimension(:),  allocatable  ::   reff_mom, tau_frac
real*8,   dimension(:,:),  allocatable  ::   qext_dt

contains

!=====================================================================
!=====================================================================

subroutine initmicro
!
!initialize microphysics
!

implicit none
integer  err, n 
real*8 :: rmin  = 0.1e-6
real*8 :: rmax  = 10.e-6
real*8 :: rbmin = 0.0001e-6
real*8 :: rbmax = 1.e-2

real*8 :: vistoir = 2.75
real*8 vrat_rt

integer i,j,k,l,nt

real*8 vrat
real*8 dr(nbin)
real*8 vol(nbin)

logical :: fullcomp3

integer, parameter :: nbin_rt = 20
integer, parameter :: nratio  = 15

real*8 rad_rt(nbin_rt),radb_rt(nbin_rt+1)

real*4, parameter :: cor_ratio(nratio) = [ 0.1,0.2,0.25,0.3,     &
                      0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,     &
                      0.8,0.9,0.99 ]

character (len=128) :: rtdata_path = 'INPUT/'

real*8 factor

integer ::  ntrace, ntprog, ntfam, ntdiag, ntdustmass, &
            ntdustnb, ntmicro, ntdustcor, ntgas, ntreff
integer ::  nticemass, nticenb, ntvapor, flag, ndx
character(len=128) :: scheme, params
character (len=128) :: filename, fieldname, tracer_name
logical mcpu0

real :: reff_rd
integer, parameter :: nlonv = 7
real*4 qextv_cld(nratio,nbin_rt,nlonv)
real*4 qscatv_cld(nratio,nbin_rt,nlonv)
real*4 gv_cld(nratio,nbin_rt,nlonv)

real*4, allocatable :: qexti_cld(:,:,:), qscati_cld(:,:,:), gi_cld(:,:,:), &
   qexti_dst(:,:), qscati_dst(:,:), gi_dst(:,:)

real*4 qextv_dst(nbin_rt,nlonv),qscatv_dst(nbin_rt,nlonv)
real*4 gv_dst(nbin_rt,nlonv) 
real :: dev2
real :: surf(nbin_rt), Rs, qtmp
integer :: nloni, l_nrefi
integer, parameter :: l_nrefv = 6
real, parameter :: qex_tmp = 2.938

!C======================================================================C
mcpu0 = (mpp_pe() == mpp_root_pe())

!***********************************************************************
!C    1) Set the various tracer index for the moment scheme
!***********************************************************************

!c  For dust: Mass and Number 
iMa_dt  = 1
iNb_dt  = 2 
!c  For water ice: Mass and Number 
iMa_cld = 3 
iNb_cld = 4 
!c  For dust core: only Mass
iMa_cor = 5
!c  For water vapor: only Mass
iMa_vap = 6 

!***********************************************************************
!C    2) Some general constants 
!***********************************************************************
!     fill variables that will be passed to the cloud scheme via cldcommon.h.

!C    Densities dust and ice

aerdens(iMa_dt) = dpden_dt
aerdens(iNb_dt) = dpden_dt
aerdens(iMa_cld)= dpden_ice
aerdens(iNb_cld)= dpden_ice
aerdens(iMa_cor)= dpden_ice


!c    nu_eff = exp( dev^2) - 1
!c    Standard deviation of the dust distribution
dev_dt = 0.63676    ! gives an effective variance of 0.5  for dust
!c    Standard deviation of the water ice distribution
dev_ice= 0.3087     ! gives an effective variance of 0.1  for water ice
!!    dev_ice= 0.05     ! gives an effective variance of 0.02  for water ice

scale_dt = 3.0*(kboltz/1.66054E-27)/44.0


!***********************************************************************
!     3) Definition of the size grid and dust moment distribution
!***********************************************************************
!c     aerad is the primary radius grid used for microphysics computation.
!c     the grid spacing is based on a volume ratio between two
!c     consecutive bins; i.e. vrat.
!c     rb defines the boundary values for each aerad bin.

!c     Volume ratio between two adjacent bins
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


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

if (do_15band) then
    nloni = 8
    l_nrefi = 7
allocate(qexti_cld(nratio,nbin_rt,nloni), stat=err)
allocate(qscati_cld(nratio,nbin_rt,nloni), stat=err)
allocate(gi_cld(nratio,nbin_rt,nloni), stat=err)
allocate(qexti_dst(nbin_rt,nloni), stat=err)
allocate(qscati_dst(nbin_rt,nloni), stat=err)
allocate(gi_dst(nbin_rt,nloni), stat=err)
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
    nloni = 5
    l_nrefi = 4
allocate(qexti_cld(nratio,nbin_rt,nloni), stat=err)
allocate(qscati_cld(nratio,nbin_rt,nloni), stat=err)
allocate(gi_cld(nratio,nbin_rt,nloni), stat=err)
allocate(qexti_dst(nbin_rt,nloni), stat=err)
allocate(qscati_dst(nbin_rt,nloni), stat=err)
allocate(gi_dst(nbin_rt,nloni), stat=err)
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

!***********************************************************************
!     3bis) Specific coagulation: size grid and dust moment distribution
!***********************************************************************
vrat_coag=(rn_coag/r1_coag)**(3./(nres_coag-1))

do i = 1, nres_coag
     rads_coag(i)  = r1_coag*vrat_coag**((i-1)/3.)
     vols_coag(i)  = 4./3.*pi*r1_coag**3*vrat_coag**(i-1)
enddo
! diameter width
deltar_coag(:)=2.*rads_coag(:)*2.**(1/3.)*(vrat_coag**(1./3.)-1)/(1+vrat_coag)**(1./3.)

!***********************************************************************
!     4) Get indexes of all tracers
!***********************************************************************

! a) First get the total number of tracer : ntrace
call get_number_tracers (MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfam )
nt_nontag = ntrace - num_tags
if (mcpu0) print*,'number tracers = ',ntrace, ' number tags = ',num_tags,' number nontag = ',nt_nontag

! b) First get the number of dust_mass to check if there are more than one distribution : ndust_mass
!! Doing the same for gas : ntrace_gas
ndust_mass= 0
nice_mass= 0
ntrace_gas= 0
do n = 1, ntrace
    if (query_method('longname', model_atmos, n, scheme, params)) then
        if( trim(scheme)== 'dust_mass' ) then
            ndust_mass= ndust_mass + 1
        endif
        if( trim(scheme)== 'ice_mass' ) then
            nice_mass= nice_mass + 1
        endif
    endif
        
    if (query_method('type2', model_atmos, n, scheme, params)) then
        if( trim(scheme)== 'gas' ) then
            ntrace_gas= ntrace_gas + 1
        endif
    endif
enddo ! ntrace

! c) Get total number of microphys tracers 
ntrace_micro=3*ndust_mass+3*nice_mass

! d) Allocate table of microphysical tracers
if ( ndust_mass > 0 ) then
    if(mcpu0)  print *, 'Allocating dust mass:  ', ndust_mass
    allocate ( dust_mass_indx(ndust_mass) )
    allocate ( dust_nb_indx(ndust_mass) )
    allocate ( dust_cor_indx(ndust_mass) )
    allocate ( micro_indx(ntrace_micro) )
    allocate ( sedim_indx(ntrace_micro) )
    allocate ( stdv(ntrace_micro) )
    allocate ( reff_mom(ntrace) )
    allocate ( tau_frac(ntrace) )
    allocate ( qext_dt(ntrace,2) )
endif
if ( nice_mass > 0 ) then
    if(mcpu0)  print *, 'Allocating ice mass:  ', nice_mass
    allocate ( ice_mass_indx(nice_mass) )
    allocate ( ice_nb_indx(nice_mass) )
    allocate ( vapor_indx(nice_mass) )
endif
if ( ntrace_gas > 0 ) then
    allocate ( gas_indx(ntrace_gas) )
endif

! e) Fill some important tables
stdv(iMa_dt) = dev_dt
stdv(iNb_dt) = dev_dt
stdv(iMa_cld)= dev_ice
stdv(iNb_cld)= dev_ice
stdv(iMa_cor)= dev_ice

! f) Get all indexes
ntdustmass= 0
ntdustnb= 0
ntdustcor= 0
nticemass= 0
nticenb= 0
ntvapor= 0
ntmicro= 0
ntgas= 0
ntreff= 0

reff_rd=0.0
tau_frac=0.d0
reff_mom=Reff_backgd  !initialize reff_mom to background

do n = 1, ntrace
    !! check longname
    if (query_method('longname', MODEL_ATMOS, n, scheme, params)) then
        if(mcpu0)  print *, 'Initracer: Tracer', n, trim(scheme), trim(params) 

        call get_tracer_names(MODEL_ATMOS, n, tracer_name)
        if(mcpu0) print *, 'Field', n, trim(tracer_name), '  ', trim(scheme), &
                                       '  ', trim(params) 

        if( trim(scheme)== 'dust_mass' ) then
            ntdustmass= ntdustmass + 1
            dust_mass_indx(ntdustmass)= n 
            if(mcpu0)  print *, 'Setting dustmass index ',  trim(scheme) 
            if(mcpu0)  then
                print*, 'the sedimentation index is', find_tagging_index(n)
            endif
        elseif( trim(scheme)== 'dust_number' ) then
            ntdustnb= ntdustnb + 1
            dust_nb_indx(ntdustnb)= n 
            if(mcpu0)  print *, 'Setting dustnb index ',  trim(scheme) 
        elseif( trim(scheme)== 'dust_core_mass' ) then
            ntdustcor= ntdustcor + 1
            dust_cor_indx(ntdustcor)= n 
            if(mcpu0)  print *, 'Setting dustcor index ',  trim(scheme) 
        endif 

        if( trim(scheme)== 'ice_mass' ) then
            nticemass= nticemass + 1
            ice_mass_indx(nticemass)= n 
            if(mcpu0)  print *, 'Setting icemass index ',  trim(scheme) 
        elseif( trim(scheme)== 'ice_number' ) then
            nticenb= nticenb + 1
            ice_nb_indx(nticenb)= n 
            if(mcpu0)  print *, 'Setting icenb index ',  trim(scheme) 
        elseif( trim(scheme)== 'vapor_mass' ) then
            ntvapor= ntvapor + 1
            vapor_indx(ntvapor)= n 
            if(mcpu0)  print *, 'Setting vapor index ',  trim(scheme) 
        endif 

    endif

    !! check type
    if (query_method('type', MODEL_ATMOS, n, scheme, params)) then
        call get_tracer_names(MODEL_ATMOS, n, tracer_name)

        if( trim(scheme)== 'microphys' ) then
            ntmicro= ntmicro + 1
            micro_indx(ntmicro)= n 
            if(mcpu0)  print *, 'Setting micro index ',  trim(scheme) 
        endif

    endif

    if (query_method('type2', MODEL_ATMOS, n, scheme, params)) then
        call get_tracer_names(MODEL_ATMOS, n, tracer_name)

        if( trim(scheme)== 'gas' ) then
            ntgas= ntgas + 1
            gas_indx(ntgas)= n
            if(mcpu0)  print *, 'Setting gas index ',  trim(scheme)
        endif

    endif
    
    !! check effective radius
    if (query_method('radius_eff', MODEL_ATMOS, n, scheme, params)) then
        ntreff= ntreff + 1
        !! check effective radius. This overrides the namelist value for background Reff
        if ( parse(params,'radius',reff_rd) == 1 ) then
            reff_mom(n) = reff_rd
        else
            reff_mom(n) = Reff_backgd
        endif
        !! check if this dust is injected
        if ( parse(params,'tau_frac',reff_rd) == 1 ) then
            tau_frac(n) = reff_rd
        else
            tau_frac(n) = 0.d0
        endif
    endif

enddo ! ntrace
 
! g) Prints and security checks

if(mcpu0)  print *, 'TB18 numbers',ntmicro,ntrace_micro,ntrace_gas
if(mcpu0)  print *, 'TB18 micro_ind',micro_indx
if(mcpu0)  print *, 'TB18 gas_ind',gas_indx
if(mcpu0)  print *, 'TB18 dust_mass_ind',dust_mass_indx
if(mcpu0)  print *, 'TB18 dust_nb_ind',dust_nb_indx
if(mcpu0)  print *, 'TB18 dust_cor_ind',dust_cor_indx
if(mcpu0)  print *, 'TB18 ice_mass_ind',ice_mass_indx
if(mcpu0)  print *, 'TB18 ice_nb_ind',ice_nb_indx
if(mcpu0)  print *, 'TB18 vapor_ind',vapor_indx
if(mcpu0)  print *, 'Reff = ',reff_mom,'tau_frac = ',tau_frac

if (mcpu0) then
    if (ntdustmass.ne.ntdustnb.or.nticemass.ne.nticenb) then
        print *, 'stop number of tracer mass and number different'
        stop 'in initracer'
    endif 
    if (ntdustmass.ne.ntdustcor.or.nticemass.ne.ntvapor) then
        print *, 'stop number of tracer mass and cor/vapor different',ntdustmass,ntdustcor,nticemass,ntvapor
        stop 'in initracer'
    endif 
    if (ntmicro.ne.ntrace_micro) then
        print *, 'stop number of micro tracer different',ntmicro,ntrace_micro
        stop 'in initracer'
    endif 
endif

! h) TAGS : filling some important tables
sedim_indx=0
do n=1,ntrace_micro
    if (ndust_mass.gt.1) then
        do nt=2,ndust_mass
            if (micro_indx(n).eq.dust_mass_indx(nt)) then
                stdv(n) = stdv(iMa_dt)
                if (micro_indx(n).gt.nt_nontag) then
                   sedim_indx(n) = find_tagging_index(dust_mass_indx(nt))
                endif
            endif
            if (micro_indx(n).eq.dust_nb_indx(nt)) then
                stdv(n) = stdv(iNb_dt)
                if (micro_indx(n).gt.nt_nontag) then
                    sedim_indx(n) = find_tagging_index(dust_nb_indx(nt))
                endif
            endif
            if (micro_indx(n).eq.dust_cor_indx(nt)) then
                stdv(n) = stdv(iMa_cor)
                if (micro_indx(n).gt.nt_nontag) then
                    sedim_indx(n) = find_tagging_index(dust_cor_indx(nt))
                endif
            endif
        enddo
    endif
    if (nice_mass.gt.1) then
        do nt=2,nice_mass
            if (micro_indx(n).eq.ice_mass_indx(nt)) then
                stdv(n) = stdv(iMa_cld)
            endif
            if (micro_indx(n).eq.ice_nb_indx(nt)) then
                stdv(n) = stdv(iNb_cld)
            endif
        enddo
    endif
enddo

do nt = 1,ntrace_micro
    if (sedim_indx(nt).gt.0) then
        sedim_indx(nt) = findloc(micro_indx,sedim_indx(nt),1) 
    endif
enddo
if(mcpu0)  print *, 'sedimentation indexes are: ',sedim_indx

! i) Make sure tau_frac adds to 1 (fraction of dust bin that contributes to background)
if (sum(tau_frac(:nt_nontag)).ne.1.0) then
    if(mcpu0)  print *, 'tau_frac does not sum to 1. forcing first dust_mass bin to 1'
    n=dust_mass_indx(1)
    tau_frac(:) = 0.d0
    tau_frac(n) = 1.d0
endif


! j) Calculate Qext for lifting based on qext from table
dev2 = 1. / ( sqrt(2.)*dev_dt )
qext_dt = 0.
Rs = 0.
surf = 0.
qtmp = 0.

do nt =1,ndust_mass
    n=dust_mass_indx(nt)
    Rs = min( max(reff_mom(n),1.e-7) , 50.e-6 )
    Rs = 1. / Rs

    do i = 1, nbin_rt
        surf(i) = 0.5 * ( derf( dlog(radb_rt(i+1)*Rs) * dev2 )     &
         -derf( dlog(radb_rt(i)  *Rs) * dev2 ) )
    enddo

    qtmp = 0.
    do i = 1, nbin_rt
        qext_dt(n,1) = qext_dt(n,1) + surf(i) * qextv_dst(i,l_nrefv)
        qext_dt(n,2) = qext_dt(n,2) + surf(i) * qexti_dst(i,l_nrefi)
        qtmp = qtmp + surf(i) * qscati_dst(i,l_nrefi)     ! calculate IR Qscat
    enddo
    qtmp = min( qtmp , 0.99999*qext_dt(n,2) )
    qext_dt(n,2) = qext_dt(n,2) - qtmp                    ! convert IR Qext to Qabs
enddo

if(mcpu0)  print *, 'qext_vis = ',qext_dt(:,1)
if(mcpu0)  print *, 'qext_ir = ',qext_dt(:,2)

deallocate(qexti_cld)
deallocate(qscati_cld)
deallocate(gi_cld)
deallocate(qexti_dst)
deallocate(qscati_dst)
deallocate(gi_dst)

return
end subroutine initmicro

!=====================================================================
!=====================================================================

end module initracer_mod

