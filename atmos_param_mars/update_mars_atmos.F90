module update_mars_atmos_mod
!
!
!  Module to calculate atmospheric CO2 condensation as well as
!  convective adjustment. There are two convective adjustment schemes:
!  Bottom to top of model adjustment (pass2) and bottom to top of unstable
!  atmosphere (legacy_convect)

use constants_mod,    only: grav, cp_air, kappa, co2_lheat
use time_manager_mod, only: time_type
use fms_mod,          only: error_mesg, FATAL,       &
                         open_namelist_file, check_nml_error, &
                         mpp_pe, mpp_root_pe, close_file,     &
                         write_version_number, stdlog,        &
                         uppercase, read_data, write_data, field_size
implicit none

!---------- interfaces ------------

public :: pass2, co2_condense, legacy_convect

!---------- shared variables ------
logical :: mix_tracers = .true.   !hard-coded logical to apply convective adjustment to tracers
logical :: mix_momentum = .true.  !hard-coded logical to apply convective adjustment to momentum


contains

subroutine co2_condense( is, js, dt, Time, temp, dt_t, &
                        p_half, p_full, precip, dmass, dt_co2 )
!
!
!  Compute CO2 condensation when  T < Tcrit
!
!  Assume that CO2 ice is immediatiately deposited
!  at the surface, augmenting the snow already present 
!
!  Tcrit is currently formulated as:
!   Tcrit= 3182.48 / ( 23.3494 - log(pres) )  ;   pres in mb 

integer, intent(in)                       ::  is, js      ! horizontal dimensions
real, intent(in)                          ::  dt          ! time step
type(time_type), intent(in)               ::  Time        ! model tim
real, intent(in), dimension(:,:,:)        ::  temp        ! temperature at midpoints [K] 
real, intent(in), dimension(:,:,:)        ::  dt_t        ! dt_t is required temperature tendency to return the atmopshere back to Tcrit 
real, intent(in), dimension(:,:,:)        ::  p_full      ! layer midpoint pressure [Pa]
real, intent(in), dimension(:,:,:)        ::  p_half      ! layer interface pressure [Pa]

real, intent(out), dimension(:,:)         ::  precip      ! CO2 snow accumulation this time step [kg/m2]
real, intent(out), dimension(:,:,:)       ::  dmass       ! The array dmass may be optionally used to modify the atmospheric mass [kg/m2]
real, intent(out), dimension(:,:,:)       ::  dt_co2      ! temperature tendency from CO2 condensation [K/s]

! Local arrays
real, dimension( size(temp,1),size(temp,2) )  :: tcrit, tnew, tdiff, delp 
real     ::  fxx
integer  ::  k

!=====================================================================

fxx= (cp_air/co2_lheat) / grav


precip(:,:)= 0.0

do k= 1, size(temp,3)

    tcrit(:,:)= 3182.48 / ( 23.3494 - log( p_full(:,:,k)*1.e-2 ) )
    tnew(:,:)= temp(:,:,k) + dt_t(:,:,k)*dt
    delp(:,:)= p_half(:,:,k+1)-p_half(:,:,k)
    tdiff(:,:)= MAX( 0.0, tcrit(:,:) - tnew(:,:)  )

    dmass(:,:,k)= fxx * tdiff(:,:) * delp(:,:)

!       Limit the extent of the mass loss 
!       This is for computational reasons. 
    dmass(:,:,k)= MIN( dmass(:,:,k), 0.1*delp(:,:) )
    precip(:,:)= precip(:,:) + dmass(:,:,k)
    
!         The required temperature tendency to return the atmopshere back to Tcrit
    dt_co2(:,:,k)= tdiff(:,:)/dt 

enddo 

end subroutine co2_condense

!============================================================================
!============================================================================

subroutine pass2( is, js, dt, Time, u, v, temp, r, dt_u, dt_v, dt_t, dt_r, &
                                                dt_u_adj, dt_v_adj, dt_t_adj, &
                                          dt_r_adj, delp, p_half, p_full )
!
!
!  Compute dry convective adjustment
!
!  Dry Lapse rate:  dT/dp = <  (1/g)(RT/p)(g/cp) =  (T/p)*kappa
!    
!


integer, intent(in)                   ::  is, js      ! horizontal dimensions
real, intent(in)                      ::  dt          ! time step
type(time_type), intent(in)           ::  Time        ! model time
real, intent(in), dimension(:,:,:)    ::  u           ! u wind at midpoints [m/s]
real, intent(in), dimension(:,:,:)    ::  v           ! v wind at midpoints [m/s]
real, intent(in), dimension(:,:,:)    ::  temp        ! temperature at midpoints [K]
real, intent(in), dimension(:,:,:,:)  ::  r           ! tracers
real, intent(in), dimension(:,:,:)    ::  dt_u        ! u wind tendency [m/s/s]
real, intent(in), dimension(:,:,:)    ::  dt_v        ! v wind tendency [m/s/s]
real, intent(in), dimension(:,:,:)    ::  dt_t        ! temperature tendency [K/s]
real, intent(in), dimension(:,:,:,:)  ::  dt_r        ! tracer tendency [/s]
real, intent(in), dimension(:,:,:)    ::  delp        ! layer pressure thickness [Pa]
real, intent(in), dimension(:,:,:)    ::  p_full      ! layer midpoint pressure [Pa]
real, intent(in), dimension(:,:,:)    ::  p_half      ! layer interface pressure [Pa]

real, intent(out), dimension(size(temp,1),size(temp,2),size(temp,3))    ::  dt_u_adj        ! output u wind tendency [m/s/s]
real, intent(out), dimension(size(temp,1),size(temp,2),size(temp,3))    ::  dt_v_adj        ! output v wind tendency [m/s/s]
real, intent(out), dimension(size(temp,1),size(temp,2),size(temp,3))    ::  dt_t_adj        ! output temperature tendency [K/s]
real, intent(out), dimension(size(r,1),size(r,2),size(r,3),size(r,4))   ::  dt_r_adj        ! output tracer tendency [/s]

! Local variables
real, dimension( size(temp,1),size(temp,2) )  :: dtdp,  &              ! d(theta)/d(P)
                                                 alf,   &              ! (1+dtdp)/(1-dtdp)
                                                 rfrac, &              ! inverse of dp(k)+dp(k+1)*alf
                                                 rfrac2                ! inverse of dp(k)+dp(k+1)
real, dimension( size(temp,1),size(temp,2),size(temp,3) )  :: tnew, &  ! calculated temperature
                                                              unew, &  ! calculated u wind
                                                              vnew     ! calculated v wind
real, dimension( size(r,1),size(r,2),size(r,3),size(r,4) ) :: rnew     ! calculated tracers
integer  ::  i, j, id, jd, k, kd, kmax, ntrace, nt   
real :: test         ! stability test

!=====================================================================

id= size(temp,1)
jd= size(temp,2)
kd= size(temp,3)
ntrace= size(r,4)

dt_r_adj=0.
dt_u_adj=0.
dt_v_adj=0.
dt_t_adj=0.
tnew(:,:,:)= temp(:,:,:) + dt_t(:,:,:)*dt
unew(:,:,:)= u(:,:,:) + dt_u(:,:,:)*dt
vnew(:,:,:)= v(:,:,:) + dt_v(:,:,:)*dt
rnew(:,:,:,:)= r(:,:,:,:) + dt_r(:,:,:,:)*dt 
kmax= kd-1
#ifdef SUPERADIABAT
kmax= kd-2   !  this allows the bottom temperature layer to remain superadiabatic
#endif SUPERADIABAT


do k= kmax, 1, -1 
    dtdp(:,:)= 0.5*kappa*(p_full(:,:,k+1)-p_full(:,:,k))/p_half(:,:,k+1)
    alf(:,:)= ( 1.0 + dtdp(:,:) )/( 1.0 - dtdp(:,:) )
    rfrac (:,:)= 1.0 / ( delp(:,:,k) + delp(:,:,k+1)*alf(:,:) )
    rfrac2(:,:)= 1.0 / ( delp(:,:,k) + delp(:,:,k+1)          )

    do i= 1, id
        do j= 1, jd 

            test = (dtdp(i,j)-1.0)*tnew(i,j,k+1) + (dtdp(i,j)+1.0)*tnew(i,j,k)

            if( test .lt. 0.0 ) then
                tnew(i,j,k) = (delp(i,j,k  )*tnew(i,j,k  ) +  &
                     delp(i,j,k+1)*tnew(i,j,k+1)   )* rfrac(i,j)
                tnew(i,j,k+1)= tnew(i,j,k) * alf(i,j) 

                if( mix_tracers ) then
                    rnew(i,j,k,:) = (delp(i,j,k  )*rnew(i,j,k,  :) + &
                         delp(i,j,k+1)*rnew(i,j,k+1,:) )* rfrac2(i,j)
                    rnew(i,j,k+1,:) = rnew(i,j,k,:) 
                endif

                if( mix_momentum ) then
                    unew(i,j,k) = (delp(i,j,k  )*unew(i,j,k  ) + &
                         delp(i,j,k+1)*unew(i,j,k+1) )* rfrac2(i,j)
                    unew(i,j,k+1) = unew(i,j,k) 

                    vnew(i,j,k) = (delp(i,j,k  )*vnew(i,j,k  ) +  &
                         delp(i,j,k+1)*vnew(i,j,k+1) )* rfrac2(i,j)
                    vnew(i,j,k+1) = vnew(i,j,k) 
                endif
            endif
        enddo    ! j loop
    enddo    ! i loop
enddo    ! k loop

dt_t_adj(:,:,:)=  ( tnew(:,:,:) - (temp(:,:,:)+dt_t(:,:,:)*dt) ) / dt 

if( mix_tracers ) dt_r_adj(:,:,:,:)= ( rnew(:,:,:,:)-(r(:,:,:,:)+ dt_r(:,:,:,:)*dt) ) / dt

if( mix_momentum ) then
    dt_u_adj(:,:,:)=  ( unew(:,:,:) - (u(:,:,:)+ dt_u(:,:,:)*dt) ) / dt 
    dt_v_adj(:,:,:)=  ( vnew(:,:,:) - (v(:,:,:)+ dt_v(:,:,:)*dt) ) / dt 
endif

end subroutine pass2

!============================================================================
!============================================================================

subroutine legacy_convect( is, js, dt, nz, u, v, temp, &
	                       r, dt_u, dt_v, dt_t, dt_r, &
                           dt_u_adj, dt_v_adj, dt_t_adj, &
                           dt_r_adj, delp, p_half, p_full )
!  purpose
!      convect checks for (superadiabatic) potential temperature
!      instabilities between adjacent layers for a given 'pi' grid
!      point. any instabilities are resolved by the process of
!      atmospheric convection. the extent and the potential temperature
!      of the convective zone are also determined. (the convective
!      zone starts at the surface and includes all layers that exchange
!      heat with the surface through convection.)
!
!  in convect, the temperature at the midpoint of a layer is used
!  as an approximation for the temperature of the whole layer.
!  for the atmosphere to be stable at a point, the potential
!  temperature must be non-decreasing as a function of altitude
!  (non-increasing as a function of the level index k for even k.)
!  if the atmosphere is unstable, convect adjusts the potential
!  temperature values until the atmosphere is stable. when convect
!  finds a set of adjacent layers that lacks stability, it changes
!  the potential temperature for all layers in the set to an
!  energy-weighted average value in a way that conserves the total
!  heat energy of the layers in the set. any adjustments made that
!  involve levels 4 or 6 (the top two layers of the troposphere)
!  include an adjustment to the potential temperature of the
!  stratosphere ( teta(3) ) to keep teta(3) a linear extrapolation
!  in log of pressure from teta(4) and teta(6). as specified in the
!  note of 7-7-82, there are three types of adjustments (types a, b,
!  and c) which are used depending on the involvement of the
!  stratosphere.  the algorithm used in convect resembles
!  mathematical induction.  start with set of the top two layers
!  of the troposphere.  make any adjustments required to produce
!  stablility for the set.  add next layer of the troposphere to
!  the set.  make any adjustments required to produce stablility
!  for the set.  continue adding layers and stablizing resultant
!  set until the set  contains all of the troposphere.
!
!  if the bottom layer gets adjusted in this process, then there
!  will be a (vertical) interval from the surface of the planet to
!  what is called 'the top of the convective zone' over which the
!  potential temperatures are constant. in this case, this interval
!  is the convective zone. (see below in code for the other case.)
!  (if two adjacent layers have the same potential temperature, the
!  adjustment made is academic and only serves to advance the
!  algorithm in determining the extent of the convective zone. thus
!  the borderline stable case is handled as unstable.)
!

integer, intent(in)                   :: is, js, nz    ! dimension lengths
real, intent(in)                      :: dt            ! time step [s]
real, intent(in), dimension(:,:,:)    ::  u            ! u wind at midpoints [m/s]
real, intent(in), dimension(:,:,:)    ::  v            ! v wind at midpoints [m/s]
real, intent(in), dimension(:,:,:)    ::  temp         ! temperature at midpoints [K]
real, intent(in), dimension(:,:,:,:)  ::  r            ! tracers
real, intent(in), dimension(:,:,:)    ::  dt_u         ! u wind tendency [m/s/s]
real, intent(in), dimension(:,:,:)    ::  dt_v         ! v wind tendency [m/s/s]
real, intent(in), dimension(:,:,:)    ::  dt_t         ! temperature tendency [K/s]
real, intent(in), dimension(:,:,:,:)  ::  dt_r         ! tracer tendency [/s]
real, intent(in), dimension(:,:,:)    ::  delp         ! layer pressure thickness [Pa]
real, intent(in), dimension(:,:,:)    ::  p_full       ! layer midpoint pressure [Pa]
real, intent(in), dimension(:,:,:)    ::  p_half       ! layer interface pressure [Pa]

real, intent(out), dimension(size(temp,1),size(temp,2),size(temp,3))    ::  dt_u_adj        ! output u wind tendency [m/s/s]
real, intent(out), dimension(size(temp,1),size(temp,2),size(temp,3))    ::  dt_v_adj        ! output v wind tendency [m/s/s]
real, intent(out), dimension(size(temp,1),size(temp,2),size(temp,3))    ::  dt_t_adj        ! output temperature tendency [K/s]
real, intent(out), dimension(size(r,1),size(r,2),size(r,3),size(r,4))   ::  dt_r_adj        ! output tracer tendency [/s]

! local variables
real, dimension(size(temp,1),size(temp,2),size(temp,3))    ::  unew        ! calculated u wind
real, dimension(size(temp,1),size(temp,2),size(temp,3))    ::  vnew        ! calculated v wind
real, dimension(size(temp,1),size(temp,2),size(temp,3))    ::  tnew        ! calculated temperature
real, dimension(size(r,1),size(r,2),size(r,3),size(r,4))   ::  rnew        ! calculated tracers
real  :: wtc(size(temp,1),size(temp,2),2*nz+3)
real  :: om(size(temp,1),size(temp,2),2*nz+3), ym(size(temp,1),size(temp,2),2*nz+3)
real  :: teta(size(temp,1),size(temp,2),2*nz+3)
real  :: pl(size(temp,1),size(temp,2),2*nz+3)
real  :: tl(size(temp,1),size(temp,2),2*nz+3)
real  :: upi(size(temp,1),size(temp,2),2*nz+3), vpi(size(temp,1),size(temp,2),2*nz+3)
real  :: plogadj(size(temp,1),size(temp,2),2*nz+3)
real  :: qpi(size(temp,1),size(temp,2),2*nz+3,size(r,4)), qcon(size(temp,1),size(temp,2),size(r,4))

integer :: i, j, k, m, kk, kj, mkj, km2, kp3, ntrace, id, jd
integer :: l_levels, l_levelm1, l_levelm2, l_levelm3, l_levelm4
real    :: sumy, sumuv, pcexp
real, dimension(size(temp,1),size(temp,2)) :: pcon
real, dimension(size(temp,1),size(temp,2)) :: tecon

!=====================================================================

l_levels = 2*nz+3
l_levelm1= 2*nz+2
l_levelm2= 2*nz+1
l_levelm3= 2*nz
l_levelm4= 2*nz-1
ntrace=size(r,4)
upi = 0.d0
vpi = 0.d0
qpi(:,:,:,:) = 0.d0
teta = 0.d0
ym = 0.d0
om = 0.d0
plogadj=0.d0
id=size(temp,1)
jd=size(temp,2)

tnew(:,:,:)= temp(:,:,:) + dt_t(:,:,:)*dt
unew(:,:,:)= u(:,:,:) + dt_u(:,:,:)*dt
vnew(:,:,:)= v(:,:,:) + dt_v(:,:,:)*dt
rnew(:,:,:,:)= r(:,:,:,:) + dt_r(:,:,:,:)*dt 

pl(:,:,3) = p_half(:,:,1)
pl(:,:,2) = pl(:,:,3)*0.5
pl(:,:,1) = pl(:,:,2)*1.e-6
ym(:,:,2) = p_half(:,:,1)/grav
do m=1,nz  
    k = 2*m
    upi(:,:,k+2) = unew(:,:,m)
    vpi(:,:,k+2) = vnew(:,:,m)
    qpi(:,:,k+2,:) = rnew(:,:,m,:)
    pl(:,:,k+2)= p_full(:,:,m)
    pl(:,:,k+3)= p_half(:,:,m+1)
    ym(:,:,k+2)=(pl(:,:,k+3)-pl(:,:,k+1))/grav
end do

pl=pl/100.   !convert pa to mbar
call potemp1(nz,tnew,plogadj,pl,om,teta)

!   wtc is a weighting factor for convection.  wtc * teta is
!   proportional to the heat energy of the layer.
wtc(:,:,3) = om(:,:,3) * ym(:,:,2)
do k = 4, l_levelm1, 2
    wtc(:,:,k) = om(:,:,k)*ym(:,:,k)
end do

!   each time through the following loop adds another layer to the
!   set and then makes the atmosphere stable for that set.
!   if levels k and k+2 are stable, the whole set of layers is
!   stable. this is because if k = 4, there are no other layers, and
!   if k > 4, the previous time through this loop left teta(k)
!   non-decreasing in altitude for levels 4 to k.

do i=1,id
    do j=1,jd
        do k=4,l_levelm3,2

!       if stable, no adjustments are needed so by-pass the next
!       section
            if (.not.(teta(i,j,k).gt.teta(i,j,k+2)))   then

!       instability found.  make adjustments required to make
!       teta(k) non-decreasing in altitude from level k+2 to the
!       tropopause.
                sumy  = wtc(i,j,k+2)
                sumuv = ym(i,j,k+2)
                do m = 1,ntrace
                    qcon(i,j,m) = qpi(i,j,k+2,m)
                end do
                tecon(i,j) = teta(i,j,k+2)
                kp3   = k+3

!       work upwards from level k+2 to the tropopause including one more
!       layer each time through this loop. the first time through this
!       loop, kj = k and tecon = teta(k+2); so we are adjusting the
!       instability we found above. thus the first time through this
!       loop we  goto 500 . after this teta(k) and teta(k-2) may no
!       longer be stable, so we keep working until we reach the
!       tropopause or until we perform an adjustment that does not cause
!       any new instability.
!
!       since each time through this loop we include one more layer, the
!       loop must handle two cases, the case of no more instabilities
!       found, and the case of new instabilities found.

                km2   = k-2
                do  mkj = 2, km2, 2
                    kj = k-(mkj-2)
!       case of no more instabilities found. only need to check for
!       pcon adjustment. then exit the mkj loop.
                    if (.not.(teta(i,j,kj).le.tecon(i,j)))   then
!       check teta at layer boundary (instead of layer midpoint) at
!       top of convective zone. adjust pcon if teta not
!       nondecreasing here.
                        if (teta(i,j,kj+1).ge.tecon(i,j))   cycle     ! k loop
                        if (k+2.ne.l_levelm1)      cycle     ! k loop
                        pcexp = 2.0*(tecon(i,j)-teta(i,j,kj+1))/(teta(i,j,kj)-teta(i,j,kj+1))
                        pcon(i,j)  = pl(i,j,kj+1)*(pl(i,j,kj)/pl(i,j,kj+1))**pcexp

                        if (pcon(i,j).lt.pl(i,j,kj-1)) then
                            pcon(i,j) = pl(i,j,kj-1)
                        endif

                        cycle     ! k loop
                    endif

!       case of new instabilities found. three subcases
!       according to kj.
                    if (teta(i,j,kj).ne.tecon(i,j))  then
                        tecon(i,j) = (teta(i,j,kj)*wtc(i,j,kj)+ &
                                    tecon(i,j)*sumy)/(wtc(i,j,kj)+sumy)
                    endif
                    do m = 1, ntrace
                        if(qpi(i,j,kj,m).ne.qcon(i,j,m)) then
                            qcon(i,j,m) = (qpi(i,j,kj,m)*ym(i,j,kj)+ &
                                        qcon(i,j,m)*sumuv)/(ym(i,j,kj)+sumuv)
                        endif
                    end do

                    do kk = kj, kp3
                        teta(i,j,kk) = tecon(i,j)
                        do m = 1,ntrace
                            qpi(i,j,kk,m)  = qcon(i,j,m)
                        end do
                    end do

                    pcon(i,j)  = pl(i,j,kj-1)
                    sumy   = sumy+wtc(i,j,kj)
                    sumuv = sumuv+ym(i,j,kj)
                end do
                cycle     ! k loop
            end if    ! stable check

!       case of no adjustment required when level k+2 added to
!       the set. if level k+2 is the midpoint of the bottom layer,
!       check for a 'shallow' convective layer. otherwise repeat
!       the k loop for the next layerr.
            if (k+2.ne.l_levelm1)   cycle     ! k loop

!       check teta at layer boundary (instead of layer midpoint)
!       at top of convective zone. adjust pcon if teta not
!       nondecreasing here.
            if (.not.(teta(i,j,l_levelm2).gt.teta(i,j,l_levelm1)))   then
                tecon(i,j) = teta(i,j,l_levelm1)
                pcexp = 2.0*(tecon(i,j)-teta(i,j,l_levelm2))/  &
                    (teta(i,j,l_levelm3)-teta(i,j,l_levelm2))
                    pcon(i,j)  = pl(i,j,l_levelm2)*(pl(i,j,l_levelm3)/pl(i,j,l_levelm2))**pcexp

                if (pcon(i,j).lt.pl(i,j,l_levelm4)) then
                    pcon(i,j) = pl(i,j,l_levelm4)
                endif
                cycle    ! k loop
            end if

!       shallow convective layer not present. whole atmosphere
!       stable. no convection found.
!       assume a little convection near the ground.
            tecon(i,j) = 0.5*(teta(i,j,2*nz+3)+teta(i,j,l_levelm1))
            pcon(i,j)  = pl(i,j,l_levelm1)
        end do    ! k loop   
    end do    ! j loop
end do    ! i loop


do k=1,nz
    m=2*k+2
    tnew(:,:,k)=TETA(:,:,m)*OM(:,:,m)
    rnew(:,:,k,:)=QPI(:,:,m,:)
end do

dt_t_adj(:,:,:)=  ( tnew(:,:,:) - (temp(:,:,:)+dt_t(:,:,:)*dt) ) / dt 
if( mix_tracers ) dt_r_adj(:,:,:,:)= ( rnew(:,:,:,:)-(r(:,:,:,:)+ dt_r(:,:,:,:)*dt) ) / dt
if( mix_momentum ) then
    dt_u_adj(:,:,:)=  ( unew(:,:,:) - (u(:,:,:)+ dt_u(:,:,:)*dt) ) / dt 
    dt_v_adj(:,:,:)=  ( vnew(:,:,:) - (v(:,:,:)+ dt_v(:,:,:)*dt) ) / dt
endif

return
end subroutine legacy_convect

!======================================================================    
!======================================================================    

subroutine potemp1(nz,t,plog_diff,p_levels,om,teta)
!     
!
!      potemp calculates the potential temperature and several other
!      quantities at each level for a given grid point.
!
!
implicit none

integer, intent(in)  :: nz                     ! number of vertical layers
real, intent(in)     :: t(:,:,:)               ! temperature field [K]
real, intent(in)     :: p_levels(:,:,:)        ! pressure field at interfaces and midpoints [Pa]

real, intent(out)    :: om(:,:,:)              ! potential temperature conversion factor
real, intent(out)    :: teta(:,:,:)            ! potential temperature [K]
real, intent(out)    :: plog_diff(:,:,:)       ! difference in log(p) between levels

!   local variables
integer :: k, l, l_levels, l_levelm1, l_levelm3, l_levelm2, l_levelm4
real  :: tl(size(t,1),size(t,2),2*nz+3)        ! temperature on levels
real  :: slope(size(t,1),size(t,2))            ! slope of potential temperature over log(p)

!=====================================================================

l_levels = 2*nz+3	
l_levelm1= 2*nz+2
l_levelm2= 2*nz+1
l_levelm3= 2*nz
l_levelm4= 2*nz-1

!     calculate pressure at all levels
om(:,:,l_levels) = 1.0
do  k = 3, l_levelm1
    om(:,:,k) = (p_levels(:,:,k)/p_levels(:,:,l_levels))**kappa
end do
om(:,:,2) = (p_levels(:,:,2)/p_levels(:,:,l_levels))**kappa

!     temperature and potential temperature at the layer midpoints.
tl(:,:,2) = 1.05*t(:,:,1)
do  l = 1, nz
    tl(:,:,2*l+2) = t(:,:,l)
end do

do  k=2,l_levelm1,2
    teta(:,:,k) = tl(:,:,k)/om(:,:,k)
end do

!     interpolate to find potential temperature at layer boundaries.
plog_diff(:,:,2)=log(2.0)
do  k = 3, l_levelm1
    plog_diff(:,:,k) = log(p_levels(:,:,k+1))-log(p_levels(:,:,k))
end do

do  k=3,l_levelm4,2
    slope(:,:)   = (teta(:,:,k+1)-teta(:,:,k-1))/(plog_diff(:,:,k)+plog_diff(:,:,k-1))
    teta(:,:,k) = teta(:,:,k-1)+slope*plog_diff(:,:,k-1)
end do

!     extrapolate for the top and bottom levels
k       = l_levelm2
slope(:,:)   = (teta(:,:,k-1)-teta(:,:,k-3))/(plog_diff(:,:,k-2)+plog_diff(:,:,k-3))
teta(:,:,k) = teta(:,:,k-1)+slope(:,:)*plog_diff(:,:,k-1)
k       = l_levels
slope(:,:)   = (teta(:,:,k-1)-teta(:,:,k-2))/plog_diff(:,:,k-2)
teta(:,:,k) = teta(:,:,k-1)+slope(:,:)*plog_diff(:,:,k-1)

return
end subroutine potemp1



end module update_mars_atmos_mod 



