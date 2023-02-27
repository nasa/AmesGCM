module pblmod_mgcm

use constants_mod, only: grav,cp=>cp_air,rgas=>rdgas
use fms_mod, only: error_mesg, FATAL, file_exist,                      &
                   open_namelist_file, check_nml_error,                &
                   mpp_pe, mpp_root_pe, close_file,                    &
                   write_version_number, stdlog

implicit none

private
public :: initpbl, pbl_driver, amespbl_nml_read


!      Altered for tracer scheme
!      parameter (nvar=3+ntrace)


real*8, parameter :: KAPA = 0.25684D0    !     A thermodynamic constant.

logical, parameter :: latent_heat = .false.

!     Set up flags for water sublimation and injecting
!     Should be set once for model and not changed
logical, parameter :: do_sublimation = .false.   !flag to sublimate water within pblmod
logical, parameter :: do_injection   = .false.   !flag to inject previously sublimated water within pblmod

real*8  :: z0_pbl, epsl0_pbl, vk_pbl, alphl0_pbl, ric_pbl
real*8  :: dtmu_pbl, rmu1mu_pbl
real*8  :: rmu

logical :: mcpu0

!  Namelist---------------------------
logical :: htflux_recalc = .false.     ! If true, calculate heat flux at the end of newpbl with updated temperature

namelist /ames_pblmod_nml/ htflux_recalc

!==================================

contains

!=====================================================================
!=====================================================================

subroutine initpbl(nz,dtime,du_pbl1)
! Initialize module
implicit none
integer nz, err, n
real*8 dtime  !nc3*dt(s)
real*8, dimension(2*nz+1) :: du_pbl1

!     Numerical Constants that need to be defined for newpbl

rmu        = 1.0
z0_pbl     = 0.01
epsl0_pbl  = 0.1 
vk_pbl     = 0.4
alphl0_pbl = 0.1

!     Calculated Constants

ric_pbl = 0.195

!     speed constants

dtmu_pbl   = rmu * dtime
rmu1mu_pbl = (1. - rmu) / rmu
do n=2,2*nz
    du_pbl1(n+1) = 0.
enddo

mcpu0 = (mpp_pe() == mpp_root_pe())

return
end subroutine initpbl

subroutine amespbl_nml_read()
integer :: unit, io, ierr

mcpu0 = (mpp_pe() == mpp_root_pe())

!     ----- read namelist -----

if (file_exist('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
        read  (unit, nml=ames_pblmod_nml, iostat=io, end=20)
        ierr = check_nml_error (io, 'ames_pblmod_nml')
    enddo
20    call close_file (unit)
endif

return 
end subroutine amespbl_nml_read

!=====================================================================
!=====================================================================

subroutine pbl_driver(nvar,nz,tl,pl,uwind,vwind,qvap,qg,gt,psurf, &
                dtime,h2oflux,qrad1, &
                polarcap,ttend,utend,vtend,qtend,qgtend, &
                rhouch,htflux,strx,stry,rkh_out,rkm_out,sup,sdn, &
                du_out,latheat,co2g, &
                ustar,thstar,cdm,cdh,dsensa,dsenss )
! main pbl driver
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
integer, intent(in) :: nvar !     Number of tracer to be mixed in the PBL
integer, intent(in) :: nz   !     Number of layer midpoints
real*8, dimension(2*nz+3), intent(in) :: tl !     Absolute Temperature (K) at boundaries and midpoints
real*8, dimension(2*nz+3), intent(in) :: pl !     Pressure (Pa) at boundaries and midpoints
real*8, dimension(nz), intent(in) :: uwind  !     U wind (m/s) at midpoints
real*8, dimension(nz), intent(in) :: vwind  !     V wind (m/s) at midpoints
real*8, dimension(:,:), intent(in) :: qvap  !     tracer Vapor MMR (kg/kg) at midpoints
real*8, dimension(nz), intent(in) :: qrad1  !     Radiation Heating Rates (K/s) at midpoints
real*8, intent(in) :: qg    !     Surface water Ice (kg/m2)
real*8, intent(in) :: gt    !     Ground Temperature (K)
real*8, intent(in) :: psurf !     Surface pressure minus model top pressure (Pa)
real*8, intent(in) :: dtime !     PBL time step
real*8, intent(inout), optional :: co2g !     Input surface CO2 ice for roughness (kg/m2)
real*8, intent(inout) :: h2oflux    !     Surface water vapor flux (kg/m2)
logical, intent(inout) :: polarcap  !     If over residual cap
!===============================================================
!     Output Arguments
!===============================================================
real*8, dimension(nz), intent(out) :: ttend !     Temperature tendency (K/s)
real*8, dimension(nz), intent(out) :: utend !     U wind tendency (m/s/s)
real*8, dimension(nz), intent(out) :: vtend !     V wind tendency (m/s/s)
real*8, dimension(nz,size(qvap,2)), intent(out) :: qtend    !     vapor MMR tendency (kg/kg/s)
real*8, intent(out) :: qgtend   !     Surface water ice tendency (kg/m2/s)
!===============================================================
!     Optional Arguments
!===============================================================
real*8, intent(inout), optional :: rhouch   !     Output Rho*Cp*Cdh*Ustar
real*8, intent(inout), optional :: strx     !     Output horizontal wind stress in X direction (zonal)
real*8, intent(inout), optional :: stry     !     Output horizontal wind stress in Y direction (meridional)
real*8, intent(inout), optional :: sup      !     Output upwards water vapor flux (kg/kg/s)
real*8, intent(inout), optional :: sdn      !     Output downwards water vapor flux (kg/kg/s)
real*8, intent(inout), optional :: htflux   !     Output sensible heat flux (W/m2)
real*4, intent(inout), optional :: latheat  !     Output latent heat flux (W/m2) for diagnostics only
!     The following two variables are calculated at boundaries
!     down to the top boundary of the second layer from the surface
real*8, intent(inout), optional :: du_out(2*nz+1)   !     Output shear (/s)
real*8, intent(inout), optional :: rkh_out(2*nz+1)  !     Output heat eddy mixing coefficient (m2/s)
real*8, intent(inout), optional :: rkm_out(2*nz+1)  !     Output momentum eddy mixing coefficient (m2/s)
!test variables...
real*8, intent(inout) :: ustar,thstar,cdm,cdh,dsensa,dsenss
!===============================================================
!===============================================================
!     Local Variables To NewPBL
!     First dimension variables passed from GCM
real*8 om(2*nz+3),dxi_pbl(2*nz+1)
real*8 upi(2*nz+3),vpi(2*nz+3),teta(2*nz+3),     &
       qrad(2*nz+3)
       
real*8 qpi(2*nz+3,size(qvap,2))
real*8 qpig
integer k, l, n, l_levels
real*8 :: dumz0
real*8 :: tgtmp

!===============================================================
!===============================================================

if (.not.present(latheat) ) latheat=0.
if (.not.present(rhouch) ) rhouch=0.
if (.not.present(strx) ) strx=0.
if (.not.present(stry) ) stry=0.
if (.not.present(rkh_out) ) rkh_out=0.
if (.not.present(sup) ) sup=0.
if (.not.present(sdn) ) sdn=0.
if (.not.present(co2g) ) co2g=0.
if (.not.present(htflux) ) htflux=0.
if (.not.present(du_out) ) du_out=0.


!  Modify surface roughness when CO2 ice is on the ground.
!  Arya, pg 164.  Update added 7/21/99

dumz0 = z0_pbl
IF(co2g .gt. 0.0) dumz0 = 1.0E-4

!  Do the same for water ice if thick enough. 
!  Slab ice roughness  has been measured to be 0.01 mm.
if(qg .GT. 100. .or. polarcap )  dumz0 = 1.0E-4

if(.not.do_injection) h2oflux=0.d0

qrad(:) = 0.d0
ttend = 0.d0
utend = 0.d0
vtend = 0.d0
qtend(:,:) = 0.d0
qgtend= 0.d0
upi = 0.d0
vpi = 0.d0
qpi(:,:) = 0.d0


teta = 0.d0
om = 0.d0

do n=1,nz  
    k = 2*n
    upi(k+2) = uwind(n)
    vpi(k+2) = vwind(n)
    qpi(k+2,:) = qvap(n,:)
    qrad(k+2)= qrad1(n)
end do	

qpig = qg


dxi_pbl(:) = 0.d0
do n = 1,2*nz
    dxi_pbl(n) = (pl(n+3)-pl(n+1))/(psurf)
enddo


l_levels = 2*nz+3
!     om is the conversion factor from teta to t.
om(l_levels) = 1.0
do  k = 3, l_levels-1
    om(k) = (pl(k)/pl(l_levels))**kapa
enddo
om(2) = (pl(2)/pl(l_levels))**kapa

tgtmp=gt/om(l_levels)

!     potential temperature
do  k=2,l_levels
    teta(k) = tl(k)/om(k)
enddo



!does the RT heating need to be converted to potemp?

call newpbl(nvar,dumz0, epsl0_pbl, vk_pbl, alphl0_pbl,           &
                  ric_pbl, dtmu_pbl,                      &
                  rmu1mu_pbl,dxi_pbl,du_out,   &
                  tgtmp, &
                  psurf,dtime,pl,om,upi,vpi,teta, &
                  htflux,strx,stry,rhouch,rkh_out,rkm_out,qrad,qpi,qpig, &
                  latheat,polarcap,sup,sdn,h2oflux,&
                  nz &
                  ,ustar,thstar,cdm,cdh,dsensa,dsenss)

do n=1,nz  
    k = 2*n

    utend(n) = (upi(k+2) - uwind(n))/dtime 
    vtend(n) = (vpi(k+2) - vwind(n))/dtime
    qtend(n,:) = (qpi(k+2,:) - qvap(n,:))/dtime
    ttend(n) = (teta(k+2)*om(k+2) - tl(k+2))/dtime

end do
qgtend = (qpig - qg)/dtime


return
end subroutine pbl_driver

!=====================================================================
!=====================================================================

subroutine newpbl(nvar,z0, epsl0, vk, alphl0, ric, dtmu,    &
               rmu1mu,dxi,du,tg, &
               pi,dt,pl,om,upi,vpi,teta, &
               htflux,strx,stry,rhouch,rkh,rkm,qrad,qpi,qpig, &
               latent,polarcap,sup,sdn,h2oflux,n, &
               ustar,thstar,cdm,cdh,dsensa,dsenss)
!
!     LegacyMellor-Yamada 2.0 turbulence 
!
implicit none

integer, intent(in) :: nvar   !! Number of tracer to be mixed in PBL (gas type)
integer :: n

!     First dimension variables passed from GCM

real*8 om(2*n+3), rhouch
real*8 upi(2*n+3),vpi(2*n+3),teta(2*n+3),     &
       pl(2*n+3),z(2*n+1),qrad(2*n+3)

!  May 2002  -  dust tracer variables (QPI)

real*8, intent(inout),dimension(:,:) :: qpi
real*8 qpig

!     Now dimension variables local to newpbl

real*8 r1roxi(n),r2roxi(n),r1roro(n),r2roro(n),vect(n)
real*8 mat(n,3)
real*8 dxi(2*n+1),q0(2*n+1),  &
       ro(2*n+1),rnum(2*n+1), &
       vg(2*n+1),xinum(2*n+1),rnumnum(2*n+1),ronum(2*n+1), &
       rkm(2*n+1),rkh(2*n+1), &
       dudz
real*8 du(2*n+1)
real*8 psi(2*n+1,nvar),   &
       pc(2*n+1,nvar),qc(2*n+1,nvar),rc(2*n+1,nvar)


!     Local variables
real*8 qsat    ! mass mixing ratio of water vapor near the surface at temp=tg
real*8 h2oflux, qtmp, qold
real*8 checkm(n,nvar),checkv(n)

real*4 latent
real*8   coef

logical polarcap
integer polarflag

integer :: i, j, k, l
real*8  :: cdh, w, dphi, thstar, ustar, tg, dt,  rmu1mu
real*8  :: htflux, stry, strx, cdm, alphl0, epsl0, pi
real*8  :: dtmu, ric, vk, alpha, tbar, hscale, rlnzz, z0
real*8  :: sup, sdn, dsensa, dsenss




!======================================================================C

psi = 0.
latent = 0.
xinum = 0.

do k = 2,2*n,2
    psi(k,1) = upi(k+2)
    psi(k,2) = vpi(k+2)
    psi(k,3) = teta(k+2)
    psi(k,4:nvar) = qpi(k+2,:)
end do


!     Calculate heights

z(2*n+1) = 0.0

do k = 2*n, 2, -2
    tbar = 0.0

    do l = k,2*n,2
        tbar = tbar + teta(l+2)*om(l+2)* ((pl(l+3) - pl(l+1))/(pi))
    end do
    tbar   = tbar / ((pl(2*n+3) - pl(k+1))/(pi))

    hscale = rgas*tbar/grav
    z(k)   = hscale * log( pl(2*n+3) / pl(k+2) )
    z(k-1) = hscale * log( pl(2*n+3) / pl(k+1) )
end do

rlnzz = log(z(2*n) / z0)

do k = 1,2*n+1
    ro(k)   = + pl(k+2) / rgas / (om(k+2)*teta(k+2))
    rnum(k) = - pi / ( ro(k)*grav )
end do

do k = 2, 2*n
    ronum(k) = ro(k) * rnum(k)
end do

do j = 1, n-1
    r1roxi(j)   = -dtmu * ronum(2*j+1) / ronum(2*j) /     &
                     dxi(2*j+2)   / dxi(2*j+1)
    r2roxi(j)   = -dtmu * ronum(2*j+1) / ronum(2*j+2) /   &
                     dxi(2*j)     / dxi(2*j+1)
    r1roro(j+1) = ronum(2*j)   / ronum(2*j+2)
    r2roro(j)   = ronum(2*j+2) / ronum(2*j)
end do

do k = 3, 2*n-1, 2
    xinum(k)   = 1. / dxi(k) / rnum(k)
    rnumnum(k) = rnum(k) * rnum(k)
end do


call eddycoef(nvar,ric,vk,alphl0,epsl0,dxi,q0,du,  &
           xinum,z,rkm,rkh,psi,dt,n)

call bndcond(nvar,tg,vk,rlnzz,z0,ustar,   &
          thstar,ric,dphi,w,z,psi,cdh,cdm,n)


alpha = atan2( psi(2*n,2), psi(2*n,1) )
strx = + ro(2*n+1)*ustar*ustar*cos(alpha)
stry = + ro(2*n+1)*ustar*ustar*sin(alpha)

htflux = - ro(2*n+1)*cp*ustar*thstar
dsenss = ro(2*n+1)*cp*cdh*ustar
dsensa = - ro(2*n+1)*cp*cdh*ustar/om(2*n+3)

!  Compute RHOUCH for the new method of determing Tg in TEMPGR:

RHOUCH = RO(2*N+1)*cp*CDH*USTAR

call scldef(nvar,ro,rnum,rkm,rkh,rnumnum,pc,qc,rc,psi,qrad,n)
call watsat_pbl(tg,pl(2*n+3)/100.,qsat)

!     Set sublimation switch
if (do_sublimation) then
    coef = 1.0D0
    qold = qpi(2*n+2,1)
else
    coef = 0.0D0
endif

!     loop over variables

do i = 1, nvar
    call matrix(nvar,i,dt,dtmu,rmu1mu,ustar,vk,rlnzz,tg,    &
           r1roro,r2roro,r1roxi,r2roxi,rnum,dxi,ro,vect,   &
           pc,qc,rc,psi,mat,cdh,cdm,qsat,coef,n)

    if (do_sublimation) then
    !       Check the amount of subliming water ice.  - Feb 2003
    !       If it exceeds the available surface water, change it.
        if (i .eq. 4) then

            !       Make a preliminary computation of vect(n) (the true vect(n) is computed in solve.f) 
            !       After scaling, vect(n) becomes qpi of water in the first layer.
            checkm(1,2) = mat(1,2) 
            checkv(1)   = vect(1)
            do j = 1, n-1
                checkm(j+1,2) = mat(j+1,2) - mat(j+1,1)*mat(j,3)/checkm(j,2)
                checkv(j+1)   = vect(j+1) - mat(j+1,1)*checkv(j)/checkm(j,2)
            end do
            qtmp = checkv(n) / checkm(n,2) / ro(2*n) 
            !         Use a temporary value of qpi (qtmp) to compute the near surface water flux*dt    
            !         (the flux accounts through dtmu for the type of scheme used; i.e. expl. or impl.) 
            !         Uses "coef" switch so sublimation won't happen if switch is off. Urata 3/11/19
            h2oflux = coef * dt * ro(2*n) / rnum(2*n) * cdh * ustar    &
                  * (qsat - dtmu/dt*qtmp - (1.0D0-dtmu/dt)*qold)
            !         If subliming more than available (qpig), change vect(n).
            !         This time, no boundary condition on the ground but water vapor flux is passed like 
            !         a mass mixing ratio tendency in the first layer. 

            !For some reason if polarcap is T, it will always evaluate to T
            !in gfortran when inside the f90 module. Still an issue as of 5/6/2016.
            !Work around this by setting polarcap to an integer flag

            if (polarcap) then
                polarflag = 1
            else
                polarflag = 0
            endif

            if((polarflag.eq.0) .and. (h2oflux.gt.qpig)) then
                h2oflux  = qpig
            end if

            mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i)
            vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) +       & 
                      (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) +  &
                      rc(2*n,i) * dt + ro(2*n) * h2oflux /   &
                      ( (pl(2*n+3)-pl(2*n+1))/grav)
        endif  !nvar=4 (water vapor)
    endif  !do_sublimation
    call solve(nvar,i,vect,mat,psi,n)

end do

call descale (nvar,ro,rnum,rkm,rkh,rnumnum,psi,n)

do k = 2,2*n,2
    upi(k+2)  = psi(k,1)
    vpi(k+2)  = psi(k,2)
    teta(k+2) = psi(k,3)
    qpi(k+2,:) = psi(k,4:nvar)
end do

! Improve energy balance by calculating heat flux with new teta from solved matrix
if (htflux_recalc) then
    htflux = - ro(2*n+1)*cp*ustar*cdh*((psi(2*n,3)*om(2*n+2)) - &
                                   (tg*om(2*n+3)))
    dsenss = ro(2*n+1)*cp*cdh*ustar
    dsensa = - ro(2*n+1)*cp*cdh*ustar/om(2*n+3)
endif

if (do_sublimation) then
!  updating water ice budget on the surface - Feb 2003

    qpig = qpig - h2oflux 
    if(.not.polarcap .and. qpig.lt.0.0) then
        qpig = 0.0D0
    end if

    if (h2oflux .gt. 0.) then
        sup = h2oflux/dt
        sdn = 0.
    else
        sdn = (-1.) * h2oflux/dt
        sup = 0.
    endif

    if(latent_heat)  latent = 2.8e+6*(-h2oflux)/dt
endif !do_sublimation

return
end subroutine newpbl


!=====================================================================
!=====================================================================

subroutine watsat_pbl(temp,press,qsat)

real*8 pvs,temp,press,qsat

! Bob's vapor pressure
pvs  = 6.11*exp(22.5*(1.0-(273.16/temp)))

qsat = pvs * 18.0 / (44.0*press)

return
end subroutine watsat_pbl


!=====================================================================
!=====================================================================

subroutine eddycoef(nvar,ric,vk,alphl0,epsl0,dxi,q0,du,   &
                 xinum,z,rkm,rkh,psi,dt,n)
! A new scheme for determining the eddy mixing coefficients
! 2 (count 'em) fiddles...errr...parameterizations are present:
! dudvdz smoothed in time to remove oscillations in lowest model level
! Kh and Km have minimum values: From Arya, Into to Micromet, p164-166
!
! Apart from all this the scheme is completely kosher
!

implicit none
integer :: n
integer, intent(in) :: nvar   !! Number of tracer to be mixed in PBL (gas type)

!     global variables:

real*8 dxi(2*n+1),q0(2*n+1),du(2*n+1),   &
       xinum(2*n+1),z(2*n+1),rkm(2*n+1),rkh(2*n+1)
real*8 psi(2*n+1,nvar)

!     implicit none

integer :: k
real*8  :: ric, vk, alphl0, epsl0, rl0, dudvdz, ri, dt, rim, rkh0
real*8  :: rkm0, rkmin, rl, beta, dudz, dvdz, dthdz

!======================================================================C

beta = 0.
dudvdz = 0.
ri = 0.
rim = 0.
rkh0 = 0.
rkm0 = 0.
rkmin = 0.
rl = 0.
dudz = 0.
dvdz = 0.
dthdz = 0.

! Specify constants

rl0    = 150.0

do k = 3, 2*n-1, 2


    ! Calculate mixing length and beta

    rl        = rl0 * vk * z(k) / (rl0 + vk * z(k))
    beta      = (dxi(k-1) + dxi(k+1)) / (dxi(k-1) *   &
                psi(k-1,3) + dxi(k+1) * psi(k+1,3))

    ! Calculate shears and Richardson number

    dudz      = (psi(k+1,1) - psi(k-1,1)) * xinum(k)
    dvdz      = (psi(k+1,2) - psi(k-1,2)) * xinum(k)
    dthdz     = (psi(k+1,3) - psi(k-1,3)) * xinum(k)
    dudvdz    = dudz*dudz + dvdz*dvdz
    ri        = beta*grav*dthdz / (dudvdz+1.0E-9)

    du(k) = du(k)-(du(k)-dudvdz)*dt/1.0E4
    rim       = beta*grav*dthdz / (du(k)+1.0E-9)

    q0(k)     = rim

    ! Neutrally stable mixing coefficients = q(Ri=0)*l*l*S(Ri=0)

    rkh0      = sqrt(du(k)/0.153) * rl * rl * 0.493
    rkm0      = sqrt(du(k)/0.153) * rl * rl * 0.393

    ! Variation of Km and Kh with Ri (from Arya p.164)

    if(rim .le. 0.0)then

        rkh(k)  = rkh0*((1.0-15.0*rim)**0.50)
        rkm(k)  = rkm0*((1.0-15.0*rim)**0.25)

    else

        rkh(k)  = rkh0*(1.0-rim/ric)
        rkm(k)  = rkm0*(1.0-rim/ric)

    endif

    ! set limiting value to RKM
    rkmin = 0.001
    !        if(z(k) .lt. 300.0)rkmin=0.1

    rkh(k)  = max(rkh(k),rkmin)
    rkm(k)  = max(rkm(k),rkmin)

end do

return
end subroutine eddycoef


!=====================================================================
!=====================================================================

subroutine bndcond(nvar,tg,vk,rlnzz,z0,ustar,    &
                thstar,ric,dphi,w,z,psi,cdh,cdm,n)
! New routine calculates ustar and thstar from rib, 
! based on Savijarvi, Icarus Vol 117 p121, section 2.1 (rib < 0)
! and Hourdin et al, JGR Vol 10 0 p5505, equations 5 and 6 (rib > 0)

implicit none
integer n, iters
integer, intent(in) :: nvar   !! Number of tracer to be mixed in PBL (gas type) + U,V,Heat

real*8 z(2*n+1), psi(2*n+1,nvar)

real*8  :: tg, vk, rlnzz, z0, ustar, thstar, ric, dphi, w, cdh 
real*8  :: cdm, rib, frih, frim

real*8  :: ricrit, beta, zeta_trans, ritrans, const
real*8, parameter :: drag_min=1.e-5
real*8 :: rich_crit
logical, parameter :: legacy_stability = .true.

!======================================================================C

w  = sqrt( psi(2*n,1)*psi(2*n,1) + psi(2*n,2)*psi(2*n,2) )

!     For now use psi at midpoint of bottom layer.
!     Later interpolate to surface

dphi = (psi(2*n,3) - tg)
rib  = ((grav * z(2*n)) / ( psi(2*n,3)*w*w + 1.0e-09))*dphi

frim=0.
frih=0.


! FV3 parameters----------------
ricrit=10.
rich_crit=0.95*ricrit
beta=1./ricrit
zeta_trans=0.5
ritrans=zeta_trans/(1.+5.*zeta_trans)

if ( rib .ge. 0.0 ) then
    if (legacy_stability) then
        ! NEWPBL version
        frih   = 1.0/(1.0+(15.0*rib/sqrt(1.0+5.0*rib)))
        frim   = 1.0/(1.0+(10.0*rib/sqrt(1.0+5.0*rib)))
    else
        ! FV3 version
        !   weakly stable
        if (rib .lt. ritrans) then
            frih = (1.-5.*rib)**2
        else
        !strongly stable
            frih = ((1.-beta*rib)/(1.+(5.-beta)*zeta_trans))**2
        endif
        frim=frih
    endif
else if ( rib .lt. 0.0 ) then

! LMD formulation
 
!           frih   = 1.0 - 15.0*rib/(1+75.0*((vk/rlnzz)**2)*
!     &              sqrt(-rib*exp(rlnzz)))
!           frim   = 1.0 - 15.0*rib/(1+75.0*((vk/rlnzz)**2)*
!     &              sqrt(-rib*(1.0+exp(rlnzz))))

! Original, incorrect version

!       frih   = sqrt(1.0-16.0*rib)
!       frim   = sqrt(1.0-64.0*rib)

    ! Corrected 5-30-06
    if (legacy_stability) then
        ! NEWPBL version
        frih = sqrt(1.0-64.0*rib)
        frim = sqrt(1.0-16.0*rib)
    else
        ! FV3 version
        const=-256.*rib**2
        call solve_unstable(fh,fhprime,frih,frih,iters,.false.,const)
        call solve_unstable(fm,fmprime,frim,frim,iters,.false.,const)
    endif
end if

if (legacy_stability) then
    cdh    = sqrt(frih)*vk/rlnzz
    cdm    = frim*((vk/rlnzz)**2)
else
    cdh    = sqrt(frih)*vk/rlnzz
    cdm    = frim*((vk/rlnzz)**2)
    if (rib .ge. rich_crit) then
        cdh = sqrt(drag_min)
        cdm = drag_min
    else
        cdh = max(cdh,sqrt(drag_min))
        cdm = max(cdm,drag_min)
    endif
endif

ustar  = sqrt(cdm)*w 
thstar = cdh*dphi

return
end subroutine bndcond


!=====================================================================
!=====================================================================

subroutine scldef(nvar,ro,rnum,rkm,rkh,rnumnum,pc,qc,rc,psi,qrad,n)
! scale the fields by number density
implicit none
integer n
integer, intent(in) :: nvar   !! Number of tracer to be mixed in PBL (gas type) + U,V,Heat

real*8 ro(2*n+1), rnum(2*n+1), rkm(2*n+1), rkh(2*n+1)
real*8 rnumnum(2*n+1)

real*8 psi(2*n+1,nvar), pc(2*n+1,nvar), qc(2*n+1,nvar)
real*8 rc(2*n+1,nvar), qrad(2*n+3)

integer :: k, m
 
!======================================================================C
 
do k = 2, 2*n
    ro(k) = ro(k) * rnum(k)
end do

!  dust tracer modification - May 2002

do m=1,nvar
    do k=2,2*n,2
        psi(k,m) = ro(k)*psi(k,m)
    end do
end do

do k = 3, 2*n-1, 2
    rkm(k) = rkm(k) / rnumnum(k) 
    rkh(k) = rkh(k) / rnumnum(k) 
end do

!  dust tracer modification - May 2002

do m=1,nvar
    do k=2,2*n,2
        pc(k,m) = 0.0
        qc(k,m) = 0.0
        rc(k,m) = 0.0
    end do
end do

do k=2,2*n,2
    rc(k,3) = ro(k)*qrad(k+2)
    rc(k,4:nvar) = 0.
end do

do k = 3, 2*n-1, 2
    pc(k, 1) = rkm(k)
    pc(k, 2) = rkm(k)
    pc(k, 3) = rkh(k)
    pc(k,4:nvar) = rkh(k)
end do

return
end subroutine scldef




!=====================================================================
!=====================================================================

subroutine matrix(nvar,i,dt,dtmu,rmu1mu,ustar,vk,rlnzz,tg,          &
           r1roro,r2roro,r1roxi,r2roxi,rnum,dxi,ro,vect, &
           pc,qc,rc,psi,mat,cdh,cdm,qsat,coef,n)

! fill the tridiagonal matrix

implicit none

!     local areas: none
!     local variables:
integer, intent(in) :: nvar   !! Number of tracer to be mixed in PBL (gas type) + U,V,Heat
real*8 qsat   ! mass mixing ratio of water vapor near the surface at temp=tg
real*8 coef
integer n

!     global arrays:

real*8 mat(n,3)
real*8 r1roro(n),r2roro(n),r1roxi(n),r2roxi(n),vect(n),    &
    rnum(2*n+1),dxi(2*n+1),ro(2*n+1),        &
    pc(2*n+1,nvar),qc(2*n+1,nvar),rc(2*n+1,nvar),       &
    psi(2*n+1,nvar)

integer :: i, j
real*8  :: dtmu, ustar, cdm, cdh, tg, vk, rlnzz, rmu1mu, dt

!======================================================================C

do j = 1, n-1
    mat(j+1, 1) = r1roxi(j)  *  pc(2*j+1,i)      
    mat(j, 3)   = r2roxi(j)  *  pc(2*j+1,i)
end do

mat(1, 1) = 0.0
mat(n, 3) = 0.0

do j = 2, n-1
    mat(j, 2) = 1. - r1roro(j) * mat(j,1) -         &
              r2roro(j) * mat(j,3)      -    &
              dtmu * qc(2*j,i) 
end do

mat(1, 2) = 1. -  r2roro(1) * mat(1,3) - dtmu  *  qc(2,i) 

do j = 2, n-1
    vect(j) = -rmu1mu * mat(j,1) * psi(2*j-2,i) +    &
            (1. + rmu1mu * (1. - mat(j,2)))   *    &
            psi(2*j,i) - rmu1mu * mat(j,3)    *    &
    psi(2*j+2,i) + rc(2*j,i) * dt 
end do

vect(1) = (1. + rmu1mu * (1. - mat(1,2)))   *     &
        psi(2,i) - rmu1mu   * mat(1,3)    *    &
        psi(4,i) + rc(2,i)   * dt        


if(i.eq.1 .or. i.eq.2) then

    !     Surface stress boundary condition in u and v

    mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i) -   &
            ( ustar / rnum(2*n) ) *    &
            ( sqrt(cdm) ) * ( dtmu / dxi(2*n) )

    vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) +   &
            (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) +    &
            rc(2*n,i) * dt


elseif (i .eq. 3) then

    !     Surface heat flux boundary condition in theta
    
    mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i) -  &
            ( ustar / rnum(2*n) ) * ( cdh ) *  &
            ( dtmu / dxi(2*n) )

    vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) + &
            (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) + &
            rc(2*n,i) * dt - ( ustar / rnum(2*n) ) * &
            ( cdh ) *     &
            ( dt * ro(2*n)* tg/dxi(2*n)) 


elseif (coef .gt. 0.) then

    !     Water vapor flux condition (sublimation from surface)
    
    mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i) - &
            ( ustar / rnum(2*n) ) * ( cdh ) * coef *    &
            ( dtmu / dxi(2*n) )

    vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) + &
            (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) + &
            rc(2*n,i) * dt - ( ustar / rnum(2*n) ) * &
            ( cdh ) * coef *   &
            ( dt * ro(2*n)* qsat/dxi(2*n))

else

    !     Surface conditions for tracer quantities 

    mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i)

    vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) +              &
    (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) +  &
    rc(2*n,i) * dt

end if

return
end subroutine matrix


!=====================================================================
!=====================================================================

subroutine solve(nvar,i,vect,mat,psi,n)

!    solve the tridiagonal matrix

implicit none
integer n

!     local arrays:

!     global arrays:
integer, intent(in) :: nvar   !! Number of tracer to be mixed in PBL (gas type) + U,V,Heat

real*8  :: mat(n,3)
real*8  :: vect(n),psi(2*n+1,nvar)
integer :: i, j

!======================================================================C

do j = 1, n-1
    mat(j+1, 2) = mat(j+1,2) - mat(j+1,1) * mat(j,3) / mat(j,2)
    vect(j+1)   = vect(j+1)  - mat(j+1,1) * vect(j)  / mat(j,2)
end do

vect(n) = vect(n) / mat(n,2)

do j = n-1, 1, -1
    vect(j) = (vect(j) - mat(j,3) * vect(j+1)) / mat(j,2)
end do

do j = 1, n
    psi(2*j, i) = vect (j)
end do

return
end subroutine solve




!=====================================================================
!=====================================================================

subroutine descale(nvar,ro,rnum,rkm,rkh,rnumnum,psi,n)  

!     descale the fields

implicit none
integer n
integer, intent(in) :: nvar   !! Number of tracer to be mixed in PBL (gas type) + U,V,Heat

!     local arrays: none

!     global arrays:

real*8 ro(2*n+1),rnum(2*n+1),rkm(2*n+1),rkh(2*n+1),   &
       rnumnum(2*n+1) 
real*8 psi(2*n+1,nvar)

integer :: k, m

!======================================================================C

do m=1,nvar
    do k = 2, 2*n, 2
        psi(k,m) = psi(k,m) / ro(k)
    end do
end do

do k = 2, 2*n
    ro(k) = ro(k) / rnum(k)
end do

do k = 3, 2*n-1, 2
    rkm(k) = rkm(k) * rnumnum(k)
    rkh(k) = rkh(k) * rnumnum(k)
end do

return
end subroutine descale


!=====================================================================


subroutine solve_unstable(f, fp, x0, x, iters, debug, c)
! Estimate the zero of f(x) using Newton's method. 
 
implicit none
real*8, intent(in) :: x0,&          ! the initial guess
                     c
real*8, external :: f, &          ! the function to find a root of
                    fp              ! function returning the derivative f'
logical, intent(in) :: debug        ! logical, prints iterations if debug=.true.
real*8, intent(out) :: x            ! the estimate x satisfying f(x)=0 (assumes Newton converged!) 
integer, intent(out) :: iters       ! the number of iterations iters
integer, parameter :: maxiter = 20
real*8, parameter :: tol=1.d-14

! Declare any local variables:
real(kind=8) :: deltax, fx, fxprime
integer :: k


! initial guess
x = x0

if (debug) then
    print 11, x
    11     format('Initial guess: x = ', e22.15)
endif

! Newton iteration to find a zero of f(x) 

do k=1,maxiter

    ! evaluate function and its derivative:
    fx = f(x,c)
    fxprime = fp(x)

    if (abs(fx) < tol) then
        exit  ! jump out of do loop
    endif

    ! compute Newton increment x:
    deltax = fx/fxprime

    ! update x:
    x = x - deltax

    if (debug) then
        print 12, k,x
        12         format('After', i3, ' iterations, x = ', e22.15)
    endif

enddo


if (k > maxiter) then
! might not have converged

    fx = f(x,c)
endif 

! number of iterations taken:
iters = k-1


end subroutine solve_unstable

!=====================================================================


real*8 function fh(x,c)
implicit none
real*8, intent(in) :: x,c

fh=x**3-2.*x**2+x+c

end function fh

!=====================================================================

real*8 function fhprime(x)
implicit none
real*8, intent(in) :: x

fhprime=3.*x**2-4.*x+1.

end function fhprime

!=====================================================================

real*8 function fm(x,c)
implicit none
real*8, intent(in) :: x,c

fm=x**5-2.*x**3+x+c

end function fm

!=====================================================================

real*8 function fmprime(x)
implicit none
real*8, intent(in) :: x

fmprime=5.*x**4-6.*x**2+1.

end function fmprime


end module pblmod_mgcm
