module palmer_topo_drag_mod
!
! OROGRAPHIC DRAG -- Palmer (1986)
!
!
!  Calculates horizontal velocity tendency due to topographic drag
!

use          mpp_mod, only: input_nml_file
use          fms_mod, only: open_namelist_file,            &
                            close_file, error_mesg, FATAL, NOTE,       &
                            mpp_pe, mpp_root_pe, stdout, stdlog,       &
                            check_nml_error, write_version_number
use       mpp_io_mod, only: mpp_open, mpp_close,                       &
                            mpp_read, mpp_write, mpp_write_meta,       &
                            MPP_NETCDF, MPP_MULTI,                     &
                            MPP_SINGLE, MPP_RDONLY, MPP_OVERWR,        &
                            mpp_get_info, mpp_get_fields,              &
                            mpp_get_atts, mpp_get_axis_data,           &
                            mpp_get_axes, axistype, fieldtype
use       fms_io_mod, only: read_data, field_size
use       fms_io_mod, only: register_restart_field, restart_file_type
use       fms_io_mod, only: save_restart, restore_state
use fms2_io_mod,      only: file_exists
use    constants_mod, only: Grav, Cp_Air, Rdgas, Pi
use horiz_interp_mod, only: horiz_interp_type, horiz_interp_init, &
                            horiz_interp_new, horiz_interp, horiz_interp_del
use diag_manager_mod, only: register_static_field, register_diag_field, send_data
use time_manager_mod, only: time_type
use mars_surface_mod, only: read_sfc_field !read orography file

implicit none

private

character(len=128) :: version = '$Id: palmer_topo_drag.F90,v  2018/06/15 20:56:59 AK$'
character(len=128) :: tagname = '$Name: mars_june2018_ak $'

character(len=*), parameter :: module='palmer_drag'

type(axistype),  save :: Axes(3)
type(fieldtype), save :: Fields(10)
type(restart_file_type), save :: Topo_restart

integer :: id_udt_topo, id_vdt_topo,  id_tdt_topo,                     &
           id_sfc_stressGW, id_topo_SD

logical :: module_is_initialized = .false.

real :: missing_value = -999.

real, dimension(:,:),  allocatable  ::  topo_SD !standard deviation[m] 

logical :: mcpu0 !debug only, identify master processor

! parameters:

character(len=64) :: restart_read='INPUT/topo_drag.res.nc'
character(len=64) :: restart_write = 'RESTART/topo_drag.res.nc'

! axis arrays for restart file


! parameters in namelist (palmer_drag_nml):

real ::   KAP=5e-6               !Legitimately Tunable Parameter (KAP) COllins et al (1997)  
logical ::do_IR_damping=.true.   !Do infrared  damping 
integer ::grid_res=24            !grid resolution. e.g 24 for c24

NAMELIST /palmer_drag_nml/                                             &
  KAP, do_IR_damping,grid_res

public palmer_drag, palmer_drag_init, palmer_drag_end, IR_damping

contains

!#######################################################################

subroutine palmer_drag                                                   &
                                             ( delt, uwnd, vwnd, atmp, &
                                           pfull, phalf, zfull, zhalf, &
                                            dtaux, dtauy, Time         ) 
!
!Orographic gravity wave drag scheme from palmer et al. 1986
! "Alleviation of a systematic westerly bias in general circulation and 
!  numerical weather prediction models through an orographic gravity wave drag parametrization"
!This code includes Eckermann et al 2011 damping rates for vertical wavelenght of 10.km
!Adapted for FV3 by Alex Kling from James Murphy's implementation 
!Oct 2018
!                                                                     
real,    intent(in) :: delt         ! DELT     Time step
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
real, intent(out), dimension(:,:,:) :: dtaux, & ! DTAUX  Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
                                    dtauy       ! DTAUY  Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
real,    dimension(size(zhalf,1),size(zhalf,2)) :: sfc_stressGW ! sfc_stressGW  surface stress for diagnostics (IDIM x JDIM)

! Work arrays and variables
integer :: i, idim
integer :: j, jdim
integer :: k, kdim

real, dimension(size(zfull,3)+1)::stress
real, dimension(size(zfull,3)) :: u,v,t,pz
real, dimension(size(zfull,3)) :: geot,du,dv,uav,tavg,rho,dtheta,theta
real, dimension(size(zfull,3)) :: dz, thtavg, dstress,dzp
real, dimension(size(zfull,3)) :: dudt,dvdt,uacc
real, dimension(size(zfull,3)) :: Upar !wind component parallel to the surface stress

real psfc
integer NLAY
integer nl
real wgw               !vertical phase/group velocity
real stressnl1 
real pln               !logarithm of pressure (p in mbar)
real damper            !IR damping



real Nsfc
!  Derived Model Variables
!  Brunt-Vaisala Frequency (N), Scale Height (H), Vertical Wind Shear (U
real, dimension(size(zfull,3)) :: N, H, UZ

!  Isentrope Displacement (dh)
real, dimension(size(zfull,3)) :: dh
!  Isentrope Displacement diagnosed from Subcritical Rimin(nl) value
real, dimension(size(zfull,3)) :: dhsat
!
!  Surface Drag (stress) (TS)
!  Surface Topographic Variance (VAR) 
real VAR, TS 

real xRi,rhosfc
real, dimension(size(zfull,3)) :: Ri,Rimin,eps

real, dimension(size(zfull,3)) :: DP !thickness of a layer 

real thetaWind !wind direction from xaxis (purely western winds  --> is 0 degree)
    
idim = size(pfull,1) !12 for c24
jdim = size(pfull,2) !12 for c24
kdim = size(pfull,3) !36 for L36
! NOTE: ,size(phalf,3) !37 for L36
NLAY= kdim

mcpu0 = (mpp_pe() == mpp_root_pe()) 

!=====extract column vectors =====
do j=1,jdim
    do i=1,idim
        u(:)=uwnd(i,j,:)
        v(:)=vwnd(i,j,:)
        t(:)=atmp(i,j,:)
        psfc=pfull(i,j,kdim) 
        pz(:)=pfull(i,j,:) !in Pa [Alex]
        VAR = topo_SD(i,j) * topo_SD(i,j)
        if(VAR.gt.1.0E6) VAR = 1.0E6     !Limit the standard deviation to 1000 meters (1000**2= 1.0E6 m2)

        !    Determine the direction of the surface stress (theta= angle between X-axis and  sfc wind)

        if    ((v(nlay)==0.) .AND.  (u(nlay)==0.)) then  !special case where u=v=0 for atan2()
            theta=0.
        else
            thetaWind=atan2(v(nlay),u(nlay))
        endif

        do nl=nlay,2,-1 !prepare  column variables
          
            dzp(nl) =zhalf(i,j,nl-1) - zhalf(i,j,nl) !compute thickness of layer
            dz(nl) = zfull(i,j,nl-1) - zfull(i,j,nl) !distance between layer midpoint
            dp(nl) = phalf(i,j,nl+1) - phalf(i,j,nl) 

! component of the wind parallel to the surface stress. From Palmer 1986, Upar is derived from the dot product of the vectors: Upar= U.Ts/|Ts|      
!      From the linerarity of the surface stress with respect to the sfc wind,  we can also use :    Upar= U.Usfc/|Usfc| 

            Upar(nl)=(u(nl)*u(nlay)+v(nl)*v(nlay))/sqrt(u(nlay)**2+v(nlay)**2)


!     Calculate potential temperature values at layer midpoints, 
            THETA(nl) = t(nl) * (psfc/pz(nl))**(Rdgas/Cp_Air)
        enddo

        ! first layer  
        dzp(1) = zhalf(i,j,2) - zhalf(i,j,1)
        dz(1)=0. 
        dp(1) = phalf(i,j,2) - phalf(i,j,1)
        THETA(1) = t(1) * (psfc/pz(1))**(Rdgas/Cp_Air)
        Upar(1)=sqrt(u(nlay)**2+v(nlay)**2)


!======================================================
     
    
        do nl=nlay,2,-1

!    Calculate the potential temperature at layer boundaries, via a simp
!   averaging of the potential temperatures at the neighboring layer mid
!    This can be imrpoved by using the layer boundary pressure (known
!   via  psfc  and  sigma) to directly calculate the pressure term
!   in the potential temperature calculation.. need a temperature at
!   the layer boundary
!
            thtavg(nl) = 0.5 * (theta(nl)+theta(nl-1))
!
!    Calculate the average temperature at layer boundaries, via a simple
!   averaging of the temperatures at the neighboring layer midpoints
!
            tavg(nl) = 0.5 * (t(nl)+t(nl-1))
!

!   Calculate the density values at layer boundaries using the average o
!  the pressures at the two neighboring layer midpoints and the average
!  temperature at that same layer boundary;  this is a poor way of calcu
!  this average density and will be improved upon using independently
!  determined values of pressure and temperature at layer boundaries
!
            rho(nl) = 0.5*(pz(nl)+pz(nl-1))/(Rdgas*Tavg(nl))
!
!  Difference in potential temperature values across layer boundaries,
! calcultes by calculating the difference between the two neighboring
! layer midpoint potential temperatures;  this should be  (nl-1)-(nl)
!
            dtheta(nl) = theta(nl-1) - theta(nl)
        enddo  !end  l=nlay,2,-1

        uav(1) = 0.0

!
!      Calculate the surface drag magnitude, using lowest layer midpoint
!     density and wind speed values, and N (Brunt-Vaisala frequency) val
!     at the interface between the bottom two atmosphere layers

!      TS = KAP * rhosfc * Nsfc * abs(U(NLAY)) * VAR !original, zonal only

        rhosfc = psfc / (Rdgas * T(nlay))
        rhosfc = rho(nlay)
        Nsfc = (grav/theta(nlay))*(dtheta(nlay)/dz(nlay))
        if(dtheta(nlay).gt.0.0) then
            Nsfc = SQRT(Nsfc)
        else
            Nsfc = 0.00000001
        endif
        Nsfc = 0.01

        TS = KAP * rhosfc * Nsfc * sqrt(U(NLAY)**2+V(NLAY)**2) * VAR 


!     Set STRESS(NLAY+1) value equal to the surface drag value
        stress(nlay+1) = TS

!
!  THIS IS THE PRIMARY  LOOP  OVER VERTICAL LAYERS, WITHIN
!  WHICH WAVE VERTICAL PROPAGATION AND BREAKING IS DIAGNOSED


        do nl=nlay,1,-1
            eps(nl) = 0.0
            dhsat(nl) = 0.0
            dh(nl) = 0.0
            Ri(nl)=0.0
            Rimin(nl)=0.0
            stress(nl) = 0.0
            dstress(nl) = 0.0
            uacc(nl) = 0.0
        enddo
!

        do nl=nlay,2,-1
             

!    Calculate the average wind speed at the top of each model layer
!   (but not the top layer), using a simple averaging of the zonal wind
!   speeds at the two neighboring layer midpoints
!
!         uav(nl) = 0.5 * (u(nl)+u(nl-1)) !use only the zonal component
            uav(nl) = 0.5 * (Upar(nl)+Upar(nl-1)) !component of the wind parallel to the surface stress


!   Calculate the difference in zonal wind speed across a layer boundary
!
            du(nl) =   u(nl-1) - u(nl)
            dv(nl) =   v(nl-1) - v(nl)

!    DO NOT IMPLEMENT WAVE PROPAGATION (more properly, 'BREAKING')
!    WITHIN THE BOTOOM THREE MODEL LAYERS
!     try excluding the bottom four layers  2/12/2010
!
            if(nl.ge.nlay-3) then
                stress(nl) = stress(nl+1)
                dstress(nl) = 0.0


!   The  ELSE below accounts for model layers above the bottom three
!  layers... in these upper layers wave breaking is permitted
            else
!
!  Determine if a CRITICAL LEVEL situation is occurring, and if it is
! deposit all of the stress within the offending layer (the layer within
! which the sign of the zonal wind speed is oppositie the sign of the
! stress magnitude at the bottom of that layer
!
!  DIAGNOSE CRITICAL LEVEL OCCURRENCE BASED UPON EITHER LAYER MIDPOINT
!  ZONAL WIND (Upar(nl)) OR AVERAGE ZONAL WIND AT THE TOP OF THE LAYER
!  (uav(nl)) BEING OPPOSITE IN SIGN TO THE STRESS VALUE AT THE BOTTOM
!  OF THE LAYER
!

                if(Upar(nl)*stress(nl+1) .lt. 0.0 .or.&
                &     uav(nl)*stress(nl+1) .lt. 0.0) then

!--------Zonal component only-----------------
!      if(u(nl)*stress(nl+1) .lt. 0.0 .or.&
!     &     uav(nl)*stress(nl+1) .lt. 0.0) then
!----------------------------------------------
!     stress at the layer top is set equal to zero
                    stress(nl) = 0.0

!     dstress at the layer top is set equal to the stress value
!   at the bottom of the layer
!
                    dstress(nl) = stress(nl+1)

!
!   for this layer within which a CRITICAL LEVEL has arisen and
!   wave momentum has been deposited, calculate the wind acceleration
!   (units of m/s/s) by dividing the  dstress value by the mass (kg) 
!   in the layer; a NEGATIVE sign indicates a westward acceleration,
!   and a positive value indicates an eastward accaleration
!

                    uacc(nl)=-dstress(nl)/(dp(nl)/grav)

!       set dh(nl) and dhsat(nl) values to -9.99 in this CRITICAL LAYER
!     occurrence to easily help identify such layers
                    dh(nl) = -9.99
                    dhsat(nl) = -9.99
!     the  goto  statement below sends the code to the next layer,
!   but in actuality there is no further need to continue upwards 
!   in this column since there is no wave energy to consider...
!   THIS CAN BE IMPROVED AND MADE MORE EFFICIENT    

                    goto 1000 
                endif


!    the above  endif  ends the treatment of CRITICAL LEVELS
!
!  Calculate the Brunt-Vaisala frequency at the top boundary of
! layer  nl
!
                N(nl) = (grav/thtavg(nl)) * (theta(nl-1)-theta(nl))/dz(nl)
!  THERE IS CURRENLY NO 'CHECK' FOR NEGATIVE N(nl) VALUES, WHICH
!  would arise in the presence of a superadiabatic lapse rate.. which
!  will not occur in the Ames model output but could, I believe, arise
!  the no-hydrostatic GITM..
!
                if(N(nl).gt.0.0) then
                    N(nl) = sqrt(N(nl))
                else
                    N(nl) = 0.00000001
                endif
!
!  Estimate the isentropic vertical displacement at the top of layer  nl
! using the calculated stress value at the value at the lnect lower laye
! boundary and the density, Brunt-Vaisala frequency and average zonal wi
! speed at the layer boundary of interest
!
!    account for radiative damping
                if (do_IR_damping) then
                     call IR_damping (pz(nl),damper,wgw)
                else
                    wgw=1.
                    damper=0.
                endif


                stressnl1 = stress(nl+1)*exp(-(dzp(nl)/wgw)*damper/88775.5)

                dh(nl) = sqrt(abs(stressnl1/(rho(nl)*kap*N(nl)*uav(nl))))

!  Calculate the Richardson Number at the top boundary of Layer  nl
!
!   MAKE Ri  and Rimin  vectors with NLAY elements
!
!    try using N(nl) value calculated above in Ri(nl) calculation
!        Ri(nl) = N(nl) / ((du(nl)/dz(nl))**2) !zonal wind only
                Ri(nl) = N(nl) / ((du(nl)/dz(nl))**2+(dv(nl)/dz(nl))**2)
!
!  Calculate the 'minimum' Richardson Number at the top of Layer  nl
! using the Brunt-Vaisala frequency, estimated isentropic vertical displ
! average zonal wind speed, and Richardson Number already calculated at 
! layer boundary of interest
!

                if(ri(nl).lt.0.0) ri(nl) = 0.0
                rimin(nl) = ri(nl) * (1.0-(n(nl)*dh(nl)/abs(uav(nl)))) /&
                & (1.0 + ((ri(nl)**0.5)*n(nl)*dh(nl)/abs(uav(nl))))**2.0
!
!  If Rimin is greater than 0.25, than there is no wave breaking within
! layer  nl  and thus the stress value at the top of Layer  nl  is the s
! as the stress value at the bottom of Layre  nl  and the dstress value 
! Layer  nl  is zero, and there is no wind acceleration within that laye
!
                if(Rimin(nl).gt.0.25) then

! Do IR damping--------
                    if (do_IR_damping) then
                     call IR_damping (pz(nl),damper,wgw)
                    else
                        wgw=1.
                        damper=0.
                    endif

                    stress(nl) = stress(nl+1)*exp(-(dzp(nl)/wgw)*damper/88775.5)
                    dstress(nl) = stress(nl+1)*(1-exp(-30000.0*damper/88775.5))
                else
!
!   if Rimin is less than 0.25, than wave breaking is occurring and some
!  (acceleration.. really deceleration for westerly flow) is being appli
!  the zonal wind in Layer  nl;
!
!  with the below inclusion of xRi I am NOT permitting
!  the value of  Ri  to be less than 0.25 when I calculate
!  the value of  epsILON... I am not sure if this is the
!  proper way to be proceeding here
!
                    xRi = max (Ri(nl),0.25000)
!
                    if(xri.eq.0.25000) then
                        eps(nl) = 0.0
                    else
                        eps(nl)=(xri**(-0.5))*(1.0+2.0*xri**0.5)*&
                        &  ((2*(xri**0.250)*(1.0+2.0*(xri**0.50))**(-0.50))-1.0)
!
!  eps  will not be less than zero, since if it did have a 
!  negative value it would result in a negative  dh  value
!
                        eps(nl) = max(eps(nl),0.0)
                    
                    endif
!  Calculate a more representative value of the isentropic displacement
! at the top of Layer  nl, and then use that value to calculate a new va
! of the stress at the top of Layer  nl, and then a value (dstress) that
! is the change in the stress across the layer.  it is the magnitude of 
! dstress value that determines the magnitude of the wind acceleration t
! results
!
!
                    dhsat(nl)= eps(nl) / (N(nl)/ABS(uav(nl)))
!
!   Determine the stress value at a layer's top boundary if wave breakin
!  is diagnosed within that layer
!
                    stress(nl) = eps(nl) * eps(nl) * kap * rho(nl)*(uav(nl)**3)/N(nl)


!  If the stress at the top of layer  nl  is greater than the stress val
! at the bottom of Layer  nl, then the stress value at the top of Layer 
! is set equal to zero and the  dstress value for layer  nl  is set equa
! to the stress value at the bottom of Layer  nl
                    if(stress(nl).ge.0.0.and.stress(nl+1).gt.0.0) then
                        if(stress(nl).lt.stress(nl+1)) then
                            dstress(nl) = stress(nl+1)-stress(nl)
                        endif
                    endif
!  add the following  IF(stress(nl).le.0.0.and.stress(nl+1).le.0.0) then
                    if(stress(nl).le.0.0.and.stress(nl+1).lt.0.0) then
                        if(stress(nl).gt.stress(nl+1)) then
                            dstress(nl) = stress(nl+1) - stress(nl)
                        endif
                    endif
!
                    if(stress(nl).eq.0.0) goto 1000

!
!  If the stress value at the top of Layer  nl  is greater than
! zero but less than the stress value at the bottom of layer  nl,
! the stress value at the top of Layer  nl  is as calculated before 
! entering this  if  statement and the dstress value is set equal to
! stress(nl) - stress(nl+1)
!
                endif
            endif !     Endif NL.LE.NLAY-3
!
!    Calculate the layer midoint wind accaleration value based upon
!   the layer's value of  dstress  and the mass within the layer

            uacc(nl) = -dstress(nl)/(dp(nl)/grav) !
      
        enddo ! end primary nl=nlay,2,-1 loop

!
!  Now, deal with the top model layer; if the stress value at the
! bottom of this top layer is zero, the layer's  dstress  value is
! set equal to the stress value at the bottom of that layer, and an
! acceleration is then calculated for the top layer;  no wave energy
! or momentum flux leaks out of the top of the model
!
        stress(1) = 0.0
        dstress(1) = 0.0
        uacc(1) = 0.0

        if(abs(stress(2)).gt.1.0E-9) then
            dstress(1) = stress(2)
            uacc(1)= - dstress(1) / (dp(1)/grav)

        endif

1000     continue


!==========save outputs===================

!-----Reproject the decceleration on the x and y axis 
        dtaux(i,j,:)=uacc*cos(thetaWind)
        dtauy(i,j,:)=uacc*sin(thetaWind) 
        sfc_stressGW(i,j)=TS      !sfc stress
    enddo !end i=1,idim
enddo  !end j=1,jdim

call diag ( dtaux, dtauy, sfc_stressGW, Time )

end subroutine palmer_drag

!===============================
subroutine IR_damping (press,damper,wgw)
!
! IR damping rates from Eckermann et al. 2011, as implemented by Jim Murphy 
! Scale-dependant infrared radiative damping rates on Mars and their role in the 
! deposition of gravity-wave momentum flux, Icarus 211 (2011) 429-442 
!
![Alex Kling]
! NOTE: Testing showed that vertical wavelengths <7.5-10km provide a strong damping that quickly inhibits GW's growth.
!       While that damping is physical in conditions where such short wavelengths are generated, there is no interest
!       in selecting such wavelenghts  for  this  monocromatic scheme since it makes it non-effective. 
!       Additionally, the pressure-log fits for the shorter wavelenghts initially considered were suitable for the 
!       lower atmosphere but provided non-physical (negative) values for low pressures, which was an issue for high-top 
!       simulations. Therefore, only the 10km vertical wavelength was retained and the magnitude of the damping is primarily 
!       controlled by the KAP parameter.      
       

real, intent(in) :: press            !pressure [Pa]
real, intent(out):: damper           !damping rate for the stress within the atmospheric layer
real, intent(out):: wgw              !vertical phase/group velocity
real :: pln                          !logarithm of the pressure (p in mbar)

mcpu0 = (mpp_pe() == mpp_root_pe())
damper = 0.4888441 
wgw = 500.0/800.0

pln=log(press/100.) !conversion from Pascal to mbar

!Use 10 km vertical wavelength from Eckermann et al. 2011
wgw = 500.0/800.0
damper = 0.74464 + 0.088451*pln +                 &
         0.0087765*pln**2 - 0.059745*pln**3 -   &
         0.014485*pln**4 + 0.0032386*pln**5 +   &
         0.0020174*pln**6 + 0.00014703*pln**7 - &
         0.000020184*pln**8 - 0.0000019778*pln**9
end subroutine IR_damping


!=======================================================================
!=======================================================================
subroutine diag ( dtaux, dtauy, sfc_stressGW, Time )
! send data to diagnostics

real, intent(in), dimension(:,:,:) :: dtaux
real, intent(in), dimension(:,:,:) :: dtauy
real, intent(in), dimension(:,:)   :: sfc_stressGW
type(time_type), intent(in) :: Time

logical :: used
if ( id_udt_topo > 0 ) then
    used = send_data ( id_udt_topo, dtaux, Time )
endif

if ( id_vdt_topo > 0 ) then
    used = send_data ( id_vdt_topo, dtauy, Time )
endif

if ( id_sfc_stressGW > 0 ) then
    used = send_data ( id_sfc_stressGW, sfc_stressGW, Time )
endif


end subroutine diag

!=======================================================================
!=======================================================================

subroutine palmer_drag_init(nlon, mlat, lonb, latb, axes_phys, Time) 
! initialize module

integer, intent(in), dimension(4) :: axes_phys
type(time_type), intent(in) :: Time
  integer, intent(in)                   :: nlon, mlat
  real,    intent(in),  dimension(:,:)  :: lonb, latb 
logical :: used
character (len=128) :: filename, fieldname
integer, dimension(3)::resolution_avail
integer id,jd
!---------------
integer :: io, ierr,unit_nml
! read namelist
#ifdef INTERNAL_FILE_NML
read (input_nml_file, nml=palmer_drag_nml, iostat=io)
ierr = check_nml_error(io,'palmer_drag_nml')
#else   
unit_nml = open_namelist_file ( )
ierr = 1
do while ( ierr /= 0 )
    read( unit_nml, nml = palmer_drag_nml, iostat = io, end = 10 )
    ierr = check_nml_error (io, 'palmer_drag_nml')
end do
10 call close_file ( unit_nml )
#endif

! write version number and namelist to logfile

call write_version_number (version, tagname)
if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=palmer_drag_nml)

!-----------------




id= size(lonb,1)-1; jd= size(latb,2)-1

allocate ( topo_SD(id,jd)  )

mcpu0 = (mpp_pe() == mpp_root_pe()) 
  
 
filename= 'INPUT/palmer_drag_input_c24_c48_c96.nc' 
!--- find nearest resolution and write to fieldname. Available resolutions are c24, c48, c96
! For example, requesting a c50 grid would return  fieldname='SD48'
! grid_res is defined in the namelist


if(mcpu0) then 
    print *, 'KAP=, ',KAP 
    print *, 'grid_res=, c',grid_res 
endif

resolution_avail=(/24,48,96/)
write(fieldname,'(a2,i2)')'SD',resolution_avail(minloc(abs(resolution_avail-grid_res))) 
!----
if( file_exists( trim( filename ) ) ) then 
    call read_sfc_field( nlon, mlat, lonb, latb, filename, fieldname, topo_SD ) 
    if(mcpu0) print *, 'I successfully read ',trim(fieldname),' in palmer_drag_input_c24_c48_c96.nc' 
else         
    if(mcpu0) print *, ' In palmer_drag_init, ',filename, 'not found'  
endif

! register diagnostics 


!static field containing the standard deviation
id_topo_SD = register_static_field ('mars_physics', 'topo_SD',            &
                                 (/axes_phys(1:2)/),                  &
                                'topography standard deviation', 'm',   &
                                 missing_value=missing_value )

if (id_topo_SD > 0) used = send_data(id_topo_SD, topo_SD, Time)

!sfc stress and decceleration tendencies
id_sfc_stressGW = &
register_diag_field ('mars_physics', 'sfc_stressGW', axes_phys(1:2), Time,         &
                 'sfc base stress for gravity wave drag', 'N/m/m',      &
                    missing_value=missing_value )

id_udt_topo =  register_diag_field ('mars_physics', 'udt_topo', axes_phys(1:3),  &
                      Time, 'utend for topo wave drag',   'm/s2',   &
                       missing_value=missing_value               )

id_vdt_topo =  register_diag_field ('mars_physics', 'vdt_topo', axes_phys(1:3),  &
                      Time, 'vtend for topo wave drag',   'm/s2',   &
                       missing_value=missing_value               )


module_is_initialized = .true.

end subroutine palmer_drag_init



!=====================================================================
subroutine palmer_drag_end

write(*,*)"Nothing done here"

end subroutine palmer_drag_end



endmodule palmer_topo_drag_mod
