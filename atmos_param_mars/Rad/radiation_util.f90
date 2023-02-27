module radiation_util_mod
!-----------------------------------------------------------------------
! module to define radiation types and other utilities
!-----------------------------------------------------------------------
use constants_mod, only : pi, grav, diffac

implicit none
private

! interfaces

public surface_type

type surface_type
    real, dimension(:,:),   pointer ::  asfc=>NULL(),   &
                                        land=>NULL(),  &
                                        asfc_vis_dir=>NULL(), &
                                        asfc_nir_dir=>NULL(), &
                                        asfc_vis_dif=>NULL(), &
                                        asfc_nir_dif=>NULL()
end type surface_type

public atmos_input_type

type atmos_input_type
     real, dimension(:,:,:), pointer :: press=>NULL(),   &
                                        temp=>NULL(), &
                                        rh2o=>NULL(),  &
                                        zfull=>NULL(),  &
                                        pflux=>NULL(), &
                                        tflux=>NULL(),  &
                                        deltaz=>NULL(),  &
                                        phalf=>NULL(),   &
                                        rel_hum=>NULL(), &
                                        cloudtemp=>NULL(),   &
                                        clouddeltaz=>NULL(), &
                                        cloudvapor=>NULL(), &
                                        aerosoltemp=>NULL(), &
                                        aerosolvapor=>NULL(), &
                                        aerosolpress=>NULL(), &
                                        aerosolrelhum=>NULL()
     real, dimension(:,:),   pointer :: tsfc=>NULL(),   &
                                        psfc=>NULL()              
end type atmos_input_type

public :: scatter_calc, planck_function, calc_plank_values,            &
          co2_trans_hourdin, co2_dust_trans, calc_15mm_fluxes,         &
          lw_scattering, planck_table_eval, planck_func, planck_inv_func 

contains

!-----------------------------------------------------------------------

subroutine scatter_calc( id, kd, tau, sscat, gfac, fscat, &
                                  coszro, albedo, heat, grnflx, outflx )
!
!  calculate flux and heating
!
!  INPUT:
integer, intent(in) :: id, kd
real, dimension(id),    intent(in) :: albedo, coszro
real, dimension(id,kd), intent(in) :: tau, sscat, gfac, fscat 

!  OUTPUT:
real, dimension(id),    intent(out) :: grnflx, outflx
real, dimension(id,kd), intent(out) :: heat 

!  LOCAL ARRAYS:

real, dimension(id,kd+1) :: dnflx, upflx, netflx, optdp, direct
real, dimension(id,kd)   :: dirct, sscatp, gfacp, taup
real, dimension(id)      :: zenith, seczen, denom, diffuse, zen2
real, dimension(id,kd)   :: tr, rf, tr1, rf1
real, dimension(id,kd+1) :: tdn, rdn1, rup, rup1 

integer :: k, kp, i

kp= kd+1

!             Similarity scaling
!           Rescale optical depths within each layer:

sscatp(:,:)= sscat(:,:)*(1.0-fscat(:,:)) / (1.0-sscat(:,:)*fscat(:,:))
gfacp(:,:)=  (gfac(:,:)-fscat(:,:))/(1.0-fscat(:,:))

taup(:,:)= (1.0-sscat(:,:)*fscat(:,:)) * tau(:,:)


!             Calculate integrated column optical depth
optdp(:,1)= 0.0
do k=1,kd
    optdp(:,k+1)= optdp(:,k) + taup(:,k)
enddo

zenith(:)= coszro(:)
seczen(:)= 35.0/sqrt( 1.224e3*zenith(:)*zenith(:) + 1.0 )
zen2(:)= 1.0/seczen(:)

do i=1,id
    if( zenith(i) .gt. 0.0 )  then
        direct(i,1)= 1.0
        do k=1,kd
            dirct (i,k  )= exp( -taup (i,k  )*seczen(i) )
            direct(i,k+1)= exp( -optdp(i,k+1)*seczen(i) )
        enddo 
    else
        dirct (i,:)= 0.0
        direct(i,:)= 0.0
    end if 
enddo  

!     Form transmission and reflectivities for the individual layers

call ecoeffs( id, kd, taup, sscatp, gfacp, zen2, dirct, tr1, rf1, tr, rf )

!          Now march upward and downward, combining levels:

!           -------------First, do the upward marching--------

!    Combine levels by marching upwards from the surface

rup (:,kp)= albedo(:)   !   surface reflectivity  :  direct beam
rup1(:,kp)= albedo(:)   !   surface reflectivity  :  diffuse radiation

do k=kd,1,-1
    diffuse(:)= tr(:,k)-dirct(:,k)                 ! total transmission - direct
    denom(:)= 1.0/( 1.0 - rf1(:,k)*rup1(:,k+1) )   !  factor for multiple refections
    rup1(:,k)= rf1(:,k) + ( tr1(:,k)*rup1(:,k+1)*tr1(:,k) )*denom(:)

    rup(:,k)= rf(:,k) + tr1(:,k)*( diffuse(:)*rup1(:,k+1) +       &
                                dirct(:,k)*rup(:,k+1) )*denom(:)
enddo

!             -----------Now, do the downward marching----------

tdn (:,1)= 1.0   ! total transmission through the top layer:  only direct beam is involved
rdn1(:,1)= 0.0   ! reflection of top layer to radiation from below:  only diffuse radiation

tdn (:,2)=  tr(:,1)
rdn1(:,2)= rf1(:,1)

do k=2,kd
    diffuse(:)= tdn(:,k)-direct(:,k)           !   % total transmission - direct through top k-1 levels
    denom(:)= 1.0/( 1.0 - rdn1(:,k)*rf1(:,k) ) !   % factor for multiple refections

    rdn1(:,k+1)= rf1(:,k) +  ( tr1(:,k)*rdn1(:,k)*tr1(:,k) )*denom(:)

    tdn(:,k+1)= direct(:,k)*tr(:,k) + tr1(:,k)*( diffuse(:) +     &
                        direct(:,k)*rdn1(:,k)*rf(:,k) )*denom(:)
enddo


do i=1,id

    if( zenith(i) .gt. 0.0 )  then

        do k=1,kp
            diffuse(i)= tdn(i,k) - direct(i,k)
            denom(i)= 1.0 / ( 1.0 - rdn1(i,k)*rup1(i,k) )
            upflx(i,k)= ( direct(i,k)*rup(i,k) +  diffuse(i)*rup1(i,k) ) * denom(i)

            dnflx(i,k)= direct(i,k) + ( diffuse(i) +                &
                       direct(i,k)*rup(i,k)*rdn1(i,k) ) * denom(i)

            netflx(i,k)= upflx(i,k)-dnflx(i,k)
        enddo

        outflx(i)= upflx(i,1)
        grnflx(i)= -netflx(i,kp)
        do k=1,kd
            heat(i,k)=  netflx(i,k+1) - netflx(i,k)
        enddo


    else

        outflx(i)= 0.0
        grnflx(i)= 0.0
        heat(i,:)= 0.0

    endif 

enddo


return
end subroutine scatter_calc

!-----------------------------------------------------------------------

subroutine ecoeffs(  id, kd, tau, w, g, zen, dirct, tr1, rf1, tr, rf  )
!
! transmission and reflection coefficients
!
!  INPUT:
integer, intent(in) :: id, kd
real, dimension(id,kd), intent(in) :: tau, w, g, dirct
real, dimension(id),    intent(in) :: zen 

!  OUTPUT:
real, dimension(id,kd), intent(out) :: tr1, rf1, tr, rf 

!  LOCAL:
real, dimension(id,kd) :: lam, u, denom, alpha, gam 

integer :: i, k
real :: expp, expm, en

!              dirct is the direct beam transmission through the layer

lam(:,:)= sqrt( 3.0*(1.0-w(:,:))*(1.0-w(:,:)*g(:,:)) )
u(:,:)= 1.5 * ( 1.0-w(:,:)*g(:,:) )/lam(:,:)

do k=1,kd      
    denom(:,k)= 1.0/( 1.0 - (lam(:,k)*zen(:))**2 )
    alpha(:,k)= 0.75*zen(:)*w(:,k)*                               &
               ( 1.0 + g(:,k)*(1.0-w(:,k)) )*denom(:,k)
    gam(:,k)= 0.5*w(:,k)*( 1.0 + 3.0*g(:,k)*                      &
                (1.0-w(:,k))*zen(:)*zen(:) )*denom(:,k)
enddo

!       tr1 and rf1 are diffuse layer transmission and reflection: independent of zenith angle

!       tr and rf are transmission and reflection for direct beam radiation
do i=1,id

    if( zen(i) .gt. 0.0 ) then
        do k=1,kd
            expp= exp( lam(i,k)*tau(i,k) )
            expm= exp(-lam(i,k)*tau(i,k) )
            en= ((u(i,k) + 1.0)**2)*expp - ((u(i,k) - 1.0)**2)*expm 

            tr1(i,k)= 4.0*u(i,k)/en
            rf1(i,k)= (u(i,k) + 1.0)*(u(i,k) - 1.0)*(expp - expm)/en

            tr(i,k)=  (alpha(i,k)+gam(i,k))*tr1(i,k) +              &
                (  (alpha(i,k)-gam(i,k))*rf1(i,k)                &
                 - (alpha(i,k)+gam(i,k)-1.0) )*dirct(i,k)

            rf(i,k)= (alpha(i,k)-gam(i,k))*tr1(i,k)*dirct(i,k) +    &
                  (alpha(i,k)+gam(i,k))*rf1(i,k) -               &
                  (alpha(i,k)-gam(i,k))

        enddo  
    else
        tr1(i,:)= 1.0
        rf1(i,:)= 0.0
        tr (i,:)= 1.0
        rf (i,:)= 0.0
    endif
enddo

return
end subroutine ecoeffs

!-----------------------------------------------------------------------

subroutine pade_eval( pc, idim, ueq, td )
!
! pade approximation
!
!  INPUT:
integer, intent(in) :: idim
real, dimension(5),    intent(in) :: pc
real, dimension(idim), intent(in) :: ueq

!  OUTPUT:
real, dimension(idim), intent(out) :: td

!  LOCAL:
real, dimension(idim) :: top, denom

!             Pade approximants:  
top(:)  =  pc(1) + ueq(:)*( pc(2) + ueq(:)*( pc(3)      )   )
denom(:)=  pc(1) + ueq(:)*( pc(4) + ueq(:)*( pc(5) + ueq(:)))

td(:)= top(:) / denom(:)

return
end subroutine pade_eval

!-----------------------------------------------------------------------

subroutine planck_function( temp, idim, planck )
!
! planck function
!
integer, intent(in) :: idim
real, dimension(idim), intent(in)    :: temp
real, dimension(idim), intent(inout) :: planck

real, parameter :: h=6.6262e-27  !  Planck Constant:     erg-sec
real, parameter :: rk=1.3806e-16 !  Boltzman constant:   erg/K
real, parameter :: c=3.0e10      !  Speed of light:      cm/sec
!      gnu= 660.0 
!      gnu= 645.0    
real, parameter :: gnu=620.0     !  this gives the closest correspondence
!                                         to the actual value over T = [150,250]     
integer :: i
real :: con1, con2, freq3, arg

!     c1= 2*h*c^2  Sometimes this has a factor of pi included
!                  The pi factor comes from integrating the normal component
!                  of intensity(radiance) over the entire solid angle.               

con1=  2*h*c
con2=  h*c/rk    !  cm-K

!      arg= h*c*gnu / ( rk*temp)
arg= h*c*gnu / rk

freq3= gnu*gnu*gnu

!             Units:   ergs/sec/cm^2 / cm
do i=1,idim
    planck(i)= con1*freq3  / (  exp( arg/temp(i) ) - 1.0 )
enddo

planck= planck * pi * c * ( 865.0-500.0 )

planck= planck * 1.e-3     !  conversion from cgs to mks units

return
end subroutine planck_function

!-----------------------------------------------------------------------

subroutine planck_func( idim, tt, ww, radiance )  
!
!
!   c1 = constant 1 = 1.19106e8 = 2hc^2 
!       where h = planck constant and c = speed of light
!   c2 = constant 2 = 1.43883e4 = hc/k microns-K  
!
!               ww should be given in microns
!

integer idim
real tt(idim),   radiance(idim)
real  ::   c1, c2, c2prime, ww  

data c1, c2  /1.19106e8,1.43883e4/


c2prime=  c2/ww
radiance(:) = ( c1 / ww**5) /  ( exp( c2prime/tt(:) ) - 1.0 ) 

return

END subroutine planck_func

!-----------------------------------------------------------------------

subroutine planck_inv_func( idim, rad, ww, tbright  )  
!
!
!   c1 = constant 1 = 1.19106e8 = 2hc^2 
!      where h = plancks first radiation const in W/(m2-sr-(cm-1)**4)
!      constant and c = speed of light
!   c2 = constant 2 = 1.43883e4 = hc/k 
!      where h, c, and k are defined as above
!

!         Note that ww is in microns : which affects the scaling of c1 

integer idim
real  ::   c1, c2, c1prime, ww  
real tbright(idim),   rad(idim)

data c1,  c2   /1.19106e8,1.43883e4/

c1prime=  c1/(ww**5)

tbright(:) =  ( c2/ww ) /   log( 1.0 + c1prime/rad(:) ) 

return

END subroutine planck_inv_func

!-----------------------------------------------------------------------

subroutine co2_trans_hourdin( id, kd, pp, pph, pstd, tcol, tu1, td1, tnu, tnd, td )
!
!   hourdin co2 transmission
!     Warning:  be sure to properly distinguish between 
!          td(:,kd,kp)       Contribution to downward flux at surface
!               and   
!          td(:,kp,kk) kk=1:kd   Contribution to upward flux from near-surface
!

!   INPUT:
integer, intent(in) :: id, kd 
real, dimension(id,kd),        intent(in)  :: pp, tcol
real, dimension(id,kd+1),      intent(in)  :: pph
real, intent(in) :: pstd

!   OUTPUT:
real, dimension(id,kd+1),      intent(out) :: tu1, tnu, td1, tnd 
real, dimension(id,kd+1,kd+1), intent(out) :: td


!   LOCAL:
real, dimension(id,kd+1) :: ubar, upbar, ubarw, upbarw 
real, dimension(id,kd)   ::  delp, delp2 
real, dimension(kd+1)    :: bwt1, bwt2, bwt3,  bdwt1, bdwt2, bdwt3
real, dimension(id)      :: ueq, ueqw, updiff, udiff, tdc, tdw
real, dimension(id,kd+1) ::  ubara, upbara, ubarwa, upbarwa
real, dimension(id,kd+1) ::  ubard, upbard, ubarwd, upbarwd
real, dimension(id)      ::  ubars, upbars, ubarws, upbarws
real, dimension(5) :: pcoef, pcoefw, pdcoef, pdcoefw
real, dimension(4) :: acent, awing

!            Pade Coefficients:   ao, a1, a2, b1, b2    200 K
!             From table 3   Hourdin (1992) 

!             Band center
data  pcoef / 2.882E-5, 1.708E-2, -3.397E-2, 1.454E-2, 5.438E-1 /
data pdcoef / 2.924E-5, 1.732E-2, -3.442E-2, 1.475E-2, 5.511E-1 /

!             Band wings 
data  pcoefw  / 2.893E-2,  1.906,  3.841,  1.895,   6.004 /
data pdcoefw  / 2.832E-2,  1.893,  3.792,  1.881,   5.977 /

!             temperature dependency of absorber amounts:  Table 2
data acent / 0.694E-1,  0.328E-3, 0.544E-2,  0.596E-5  /
data awing / 0.275E-1, -0.705E-3, 0.356E-1, -0.326E-4  /

!            Weighting of u and up for equivalent absorber amounts: Table 2

real, parameter :: c1c=0.005, c2c=0.01
real, parameter :: c1w=0.015, c2w=0.1
real, parameter :: grav_cgs= grav*100.
real, parameter :: tzero= 200.0
real, parameter :: co2_frac= 0.95

integer :: i, k, kp, kk
real :: factor_1, factor_2, deltat, exf1, exf2, dpa,             &
      dpa2, dpd, dpd2

kp= kd+1

factor_1= diffac/grav_cgs
factor_2= factor_1 / pstd

do k=1,kp
    bwt1(k)= (135./365.)
    bwt2(k)= (70./365.)
    bwt3(k)= (160./365.)

    bdwt1(k)= (135./365.)
    bdwt2(k)= (70./365.)
    bdwt3(k)= (160./365.)
enddo

do k=1,kd
    delp (:,k)=   pph(:,k+1)-pph(:,k)
    delp2(:,k)= ( pph(:,k+1)*pph(:,k+1)-pph(:,k)*pph(:,k) )*0.5
enddo
 
!   Integrate downwards to obtain path lengths at half and full levels

ubar  (:,1)= 0.0
upbar (:,1)= 0.0
ubarw (:,1)= 0.0
upbarw(:,1)= 0.0

do k=1,kd
    do i=1,id

        !             calculate path lengths in the band center:

        deltat = tcol(i,k) - tzero
        exf1= exp( deltat*( acent(1) + acent(2)*deltat ) )*co2_frac
        exf2= exp( deltat*( acent(3) + acent(4)*deltat ) )*co2_frac

        ubar(i,k+1) = ubar(i,k)  + factor_1*exf1 * delp (i,k)
        upbar(i,k+1)= upbar(i,k)  + factor_2*exf2 * delp2(i,k)


        !              Nearest neighbor path; ascending
        dpa= pp(i,k)-pph(i,k) 
        dpa2= dpa*( pp(i,k)+pph(i,k) )*0.5

        ubara(i,k)= factor_1*exf1 * dpa
        upbara(i,k)= factor_2*exf2 * dpa2 

        !              Nearest neighbor path; descending
        dpd= pph(i,k+1)-pp(i,k) 
        dpd2= dpd * ( pph(i,k+1)+pp(i,k) )*0.5

        ubard(i,k+1)= factor_1*exf1 * dpd
        upbard(i,k+1)= factor_2*exf2 * dpd2


        !                Repeat for the wings:                 
        deltat = tcol(i,k) - tzero
        exf1= exp( deltat*( awing(1) + awing(2)*deltat ) )*co2_frac
        exf2= exp( deltat*( awing(3) + awing(4)*deltat ) )*co2_frac

        ubarw(i,k+1) = ubarw(i,k)  + factor_1*exf1 * delp (i,k)
        upbarw(i,k+1)= upbarw(i,k)  + factor_2*exf2 * delp2(i,k)


        !              Nearest neighbor path; ascending
        ubarwa(i,k)= factor_1*exf1 * dpa
        upbarwa(i,k)= factor_2*exf2 * dpa2

        !              Nearest neighbor path; descending
        ubarwd(i,k+1)= factor_1*exf1 * dpd
        upbarwd(i,k+1)= factor_2*exf2 * dpd2

    enddo
enddo 

!  Final absorber amounts:   to between pp(kd) and psfc 
do i=1,id
    deltat = tcol(i,kd) - tzero
    dpd= 0.5*( pph(i,kp)-pp(i,kd) )
    dpd2= 0.5*( (pp(i,kd)+dpd)**2 - pp(i,kd)**2 )

    exf1= exp( deltat*( acent(1) + acent(2)*deltat ) )*co2_frac
    exf2= exp( deltat*( acent(3) + acent(4)*deltat ) )*co2_frac

    ubars(i) = ubar(i,kp)  - factor_1*exf1 * dpd
    upbars(i)= upbar(i,kp)  - factor_2*exf2 * dpd2


    exf1= exp( deltat*( awing(1) + awing(2)*deltat ) )*co2_frac
    exf2= exp( deltat*( awing(3) + awing(4)*deltat ) )*co2_frac

    ubarws(i) = ubarw(i,kp)  - factor_1*exf1 * dpd
    upbarws(i)= upbarw(i,kp)  - factor_2*exf2 * dpd2
enddo 


!                    Downward paths:
!           use top temperature for pade coeffs

td1(:,1)= 1.0
tu1(:,kp)= 1.0

do k=2,kp
    !                equivalent absorber amount:
    ueq(:) =  sqrt(  upbar(:,k) ) + c1c*( ubar (:,k)**c2c )
    !                       wings: 
    ueqw(:) = sqrt( upbarw(:,k) ) + c1w*( ubarw(:,k)**c2w )

    !                pade approximants:  
    call pade_eval( pcoef,  id, ueq,  tdc )
    call pade_eval( pcoefw, id, ueqw, tdw )

    !             finally, combine with planck function weighting:
    td1(:,k)= tdc(:)*bwt2(1) + tdw(:)*( bwt1(1)+bwt3(1) )  
enddo


!              Upward paths:
!           Use surface temperature for pade coeffs
do k=1,kd
    updiff(:)= upbar(:,kp)-upbar(:,k)
    udiff(:)=   ubar(:,kp)- ubar(:,k)
    ueq(:) = sqrt( updiff(:) ) + c1c*( udiff(:)**c2c )

    !                   wings: 
    !                equivalent absorber amount:
    updiff(:)= upbarw(:,kp)-upbarw(:,k)
    udiff(:)=   ubarw(:,kp)- ubarw(:,k)
    ueqw(:) = sqrt( updiff(:) ) + c1w*( udiff(:)**c2w )

    !                pade approximants:  
    call pade_eval( pcoef , id, ueq , tdc )
    call pade_eval( pcoefw, id, ueqw, tdw )

    !          finally, combine with planck function weighting:
    tu1(:,k)= tdc(:)*bwt2(kp) + tdw(:)*( bwt1(kp)+bwt3(kp) )
enddo

!          Compute nearest-neighbor transmissions

tnu(:,kp)= 0.0
do k=1,kd
    ueq (:)= sqrt( upbara (:,k) ) + c1c*( ubara (:,k)**c2c )
    ueqw(:)= sqrt( upbarwa(:,k) ) + c1w*( ubarwa(:,k)**c2w )

    call pade_eval( pdcoef,  id, ueq,  tdc )
    call pade_eval( pdcoefw, id, ueqw, tdw )

    tnu(:,k)= tdc(:)*bdwt2(k) + tdw(:)*( bdwt1(k)+bdwt3(k) )
enddo

tnd(:,1)= 1.0
do k=2,kp
    ueq (:)= sqrt( upbard (:,k) ) + c1c*( ubard (:,k)**c2c )
    ueqw(:)= sqrt( upbarwd(:,k) ) + c1w*( ubarwd(:,k)**c2w )

    call pade_eval( pdcoef,  id, ueq,  tdc )
    call pade_eval( pdcoefw, id, ueqw, tdw )

    tnd(:,k)= tdc(:)*bdwt2(k) + tdw(:)*( bdwt1(k)+bdwt3(k) )
enddo

! ----------------  Differential Paths  ----------------:                      

do k=1,kp
    td(:,k,k)= 1.0
enddo

do k=1,kp 
    do kk=1,k-1
        !                  center:
        udiff(:)  =   ubar(:,k) -  ubar(:,kk)
        updiff(:) =   upbar(:,k) - upbar(:,kk)
        !                equivalent absorber amount:
        ueq(:) = sqrt( updiff(:) ) + c1c*udiff(:)**c2c

        !                contribution from the wings:  
        udiff(:)  =   ubarw(:,k) -  ubarw(:,kk)
        updiff(:) =  upbarw(:,k) - upbarw(:,kk) 
        ueqw(:) = sqrt( updiff(:) ) + c1w*udiff(:)**c2w
        !                pade approximants:  
        call pade_eval( pdcoef,  id, ueq,  tdc )
        call pade_eval( pdcoefw, id, ueqw, tdw )

        !              planck function weighting:  
        td(:,k,kk)=  tdc(:)*bdwt2(k) +                             &
                                     tdw(:)*( bdwt1(k) + bdwt3(k) )

        td(:,kk,k)= tdc(:)*bdwt2(kk) +                             &
                                   tdw(:)*( bdwt1(kk) + bdwt3(kk) )
    enddo 
enddo

!          Final special differential path length: near-surface to up
!                    ie   td(:,kp,kk);  kk= kd, 1, -1
k= kp
do kk=1,kd
    !                  center:
    udiff(:)  =    ubars(:) -  ubar(:,kk)
    updiff(:) =   upbars(:) - upbar(:,kk)
    !                equivalent absorber amount:
    ueq(:) = sqrt( updiff(:) ) + c1c*udiff(:)**c2c

    !                wings:  
    udiff(:)  =   ubarws(:) -  ubarw(:,kk)
    updiff(:) =  upbarws(:) - upbarw(:,kk) 
    ueqw(:) = sqrt( updiff(:) ) + c1w*udiff(:)**c2w

    call pade_eval( pdcoef,  id, ueq,  tdc )
    call pade_eval( pdcoefw, id, ueqw, tdw )

    !              planck function weighting:  
    td(:,k,kk)=  tdc(:)*bdwt2(k) + tdw(:)*( bdwt1(k) + bdwt3(k) )

enddo 

!             warning:  be sure to proper distinguish between 
!               td(:,kd,kp)  and   td(:,kp,kk) kk=1:kd

return
end subroutine co2_trans_hourdin

!-----------------------------------------------------------------------

subroutine co2_dust_trans( id, kd, tau_od, pp, pph, tcol, &
                                               td, tu1, td1, tnu, tnd )
!
!  co2 dust transmission
!
!     Warning:  be sure to properly distinguish between 
!          td(:,kd,kp)       Contribution to downward flux at surface
!               and   
!          td(:,kp,kk) kk=1:kd   Contribution to upward flux from near-surface
!

!   INPUT:
integer, intent(in)  :: id, kd 
real, dimension(id,kd),   intent(in) :: pp, tcol, tau_od
real, dimension(id,kd+1), intent(in) :: pph

!   INPUT/OUTPUT:
real, dimension(id,kd+1),      intent(inout)  ::  tu1, tnu, td1
real, dimension(id,kd+1),      intent(inout)  ::  tnd
real, dimension(id,kd+1,kd+1), intent(inout)  ::  td

!   LOCAL:
real, dimension(id)           ::  dfacu, dfacd
real, dimension(id,kd+1)      ::  taut
real, dimension(id,kd+1,kd+1) ::  dep
real, dimension(id,kd)        ::  delpi

integer :: k, kp, kk

kp=kd+1

do k=1,kd
    delpi(:,k)= 1.0 / ( pph(:,k+1)-pph(:,k) )
enddo

taut(:,1)= 0.0
do k=1,kd
    taut(:,k+1)= taut(:,k) + tau_od(:,k)
enddo

!            now form differential paths 
do k =1,kp
    do kk=1,kp
        dep(:,k,kk)=  abs( taut(:,k)-taut(:,kk) )
    enddo
enddo


td(:,1:kd,:) = td(:,1:kd,:) * exp( -dep(:,1:kd,: ) )

!          td(:,kp,1:kd)  is the tranmission from the bottom half-level to level k 
td(:,kp,1:kd-1)= td(:,kp,1:kd-1) * exp( -dep(:,kp,1:kd-1) )
td(:,kp,kd) =    td(:,kp,kd)     * exp( -dep(:,kp,kd )*0.5 )


td1(:,:  ) = td1(:,:  ) * exp( -dep(:,:,1 ) )
tu1(:,:  ) = tu1(:,:  ) * exp( -dep(:,:,kp) )

!           contributions from nearest neighbors
do k=1,kd
    dfacu(:)=  ( pp (:,k)-pph(:,k  ) )*tau_od(:,k)*delpi(:,k)
    dfacd(:)=  ( pph(:,k+1)-pp (:,k) )*tau_od(:,k)*delpi(:,k)
    tnu(:,k)   = tnu(:,k  ) * exp( -diffac * dfacu(:) )
    tnd(:,k+1) = tnd(:,k+1) * exp( -diffac * dfacd(:) )
enddo

return
end subroutine co2_dust_trans

!-----------------------------------------------------------------------

subroutine calc_plank_values( id, kd, sigma, emiss, tcol, tsfc,  &
                                    bsource, bsource_sfc, bfl )
!
!  calculate planck values
!
!   INPUT:
integer, intent(in) :: id, kd 
real, intent(in)    :: sigma
real, dimension(id),      intent(in)  :: tsfc, emiss
real, dimension(id,kd+1), intent(in)  :: tcol

!   OUTPUT:
real, dimension(id,kd+1), intent(out) :: bsource, bfl 
real, dimension(id),      intent(out) :: bsource_sfc

!   LOCAL:
integer, dimension(id,kd+1) :: ind_temp
integer, dimension(id)      :: ind_temp1d
real, dimension(id,kd+1)    :: tfl

integer :: k, kp

kp=kd+1

!           Planck function on full model levels 

call planck_function( tcol, id*kp, bsource )
call planck_function( tsfc, id,    bsource_sfc )
bsource_sfc(:)=  bsource_sfc(:)*emiss(:)


!         Calculate temperatures at flux levels

tfl(:,1)= tcol(:,1)
do k=2,kd
    tfl(:,k)= 0.5*( tcol(:,k) + tcol(:,k-1) )
enddo
tfl(:,kp)= tcol(:,kp)

call planck_function( tfl, id*kp, bfl )

return
end subroutine calc_plank_values

!-----------------------------------------------------------------------

subroutine calc_15mm_fluxes( id, kd, emiss, td, tu1, td1, tnu, tnd, &
                                    bsource, bsource_sfc, bfl, fup, fdn )
!
!  calculate 15 micron fluxes
!     Warning:  be sure to properly distinguish between 
!          td(:,kd,kp)       Contribution to downward flux at surface
!               and   
!          td(:,kp,kk) kk=1:kd   Contribution to upward flux from near-surface
!
!   INPUT:
integer, intent(in) :: id, kd 

real, dimension(id),           intent(in) :: emiss, bsource_sfc 
real, dimension(id,kd+1),      intent(in) :: tu1, tnu, td1, tnd 
real, dimension(id,kd+1),      intent(in) :: bsource, bfl 
real, dimension(id,kd+1,kd+1), intent(in) :: td

!   OUTPUT:
real, dimension(id,kd+1),     intent(out) :: fup, fdn

!   LOCAL:
real, dimension(id,kd) :: b_diff 
real, dimension(id)    :: bsurf_diff

integer :: k, kp, kk

kp= kd+1

do k=1,kd
    b_diff(:,k)= bsource(:,k) - bsource(:,k+1)
enddo
bsurf_diff(:)= bsource_sfc(:) - bsource(:,kp)


fdn(:,1)= 0.0
do k=2,kp
    fdn(:,k)=   bsource(:,1)*td1(:,k) - bfl(:,k)

    do kk=1,k-2
        fdn(:,k)= fdn(:,k) - b_diff(:,kk)*td(:,kk+1,k)
    enddo
    !     now include nearest-neighbor contribution (from above)
    fdn(:,k)= fdn(:,k) -                                         &
               (bsource(:,k-1)-bfl(:,k))*0.5*(tnd(:,k)+1.0)
enddo


!     --------Now, formulate upward flux
!            Be sure to include contribution from "reflection" of 
!           downward flux by surface with non-unit emissivity-----------

!                 Warning:  be sure to proper distinguish between 
!          td(:,kd,kp)       Contribution to downward flux at surface
!               and   
!    td(:,kp,kk) kk=1:kd   Contribution to upward flux from near-surface
!                                                    layer 

fup(:,kp)= bsource_sfc(:) - (1.0-emiss(:))*fdn(:,kp)

bsurf_diff(:)= bsurf_diff(:) - (1.0-emiss(:))*fdn(:,kp)

do k=1,kd
    fup(:,k)=  bsurf_diff(:)*tu1(:,k) + bfl(:,k)  

    do kk=k,kd
        fup(:,k)= fup(:,k) - b_diff(:,kk)*td(:,kk+1,k)
    enddo
    !        now include nearest-neighbor contribution (from below)
    fup(:,k)= fup(:,k) -                                        &
                 (bfl(:,k)-bsource(:,k))*0.5*(tnu(:,k)+1.0)
enddo

return
end subroutine calc_15mm_fluxes

!-----------------------------------------------------------------------

subroutine lw_scattering( id, kd, tau, sscat, gfac, emis, &
                                        bh, bsol, fuh, fdh, coszen )
!
! Calculation of upward (fuh) and downward (fdh) fluxes
! 2-stream approximation; hemispheric average;
! with the source function technique
! Fluxes are calculated at the interfaces of kd layers
! Layers are ordered from the top of the atmosphere downward
!

!    INPUT:
integer, intent(in)                    :: id, kd
real,dimension(id,kd),   intent(in)    :: tau, &    ! optical depth: k'th layer
                                        sscat       ! single scattering albedo: k'th layer
real,dimension(id,kd),   intent(inout) :: gfac      ! asymmetry parameter: k'th layer
real,dimension(id,kd+1), intent(in)    :: bh        ! Planck function: top of k'th lyer
real,dimension(id),      intent(in)    :: bsol, &   ! Planck function at the surface
                                        emis        ! surface emissivity
real,                    intent(in)    :: coszen

!    OUTPUT:
real,dimension(id,kd+1), intent(out) :: fuh, &     ! upward flux at the top of the k'th layer
                                     fdh        ! downward flux at the top of the k'th layer,

!   LOCAL:

real, dimension(id,kd)   :: b0, b1, e1, e2, e3, e4
real, dimension(id,kd)   :: cut, cdt, cub, cdb
real, dimension(id,kd*2) :: adiag, bdiag, ddiag, rhs, y
real, dimension(id,kd)   :: grg, grh, grj, grk, alpha1, alpha2,    &
                         sigma1, sigma2
real, dimension(id,kd)   :: beta, gama1, gama2, grgama,            &
                         alambda, denom, exfact

!   Gaussian points and weights

integer, parameter :: ngauss = 8

real :: xg(ngauss), wt(ngauss)
real :: exfact2(id,ngauss), exfactg(id,kd), accum(id,ngauss) 
real :: gri(id,kd,ngauss)

DATA xg /1.9855071751231860E-2 , 0.1016667612931866E+0,            &
  	  0.2372337950418355E+0 , 0.4082826787521751E+0,            &
  	  0.5917173212478250E+0 , 0.7627662049581645E+0,            &
  	  0.8983332387068134E+0 , 0.9801449282487682E+0 /

DATA wt /5.0614268145185310E-2 , 0.1111905172266872E+0,            &
  	  0.1568533229389437E+0 , 0.1813418916891810E+0,            &
  	  0.1813418916891810E+0 , 0.1568533229389437E+0,            &
  	  0.1111905172266872E+0 , 5.0614268145185310E-2 /

real, parameter :: eps=1.0e-10  ! small number to avoid infinities
real, parameter :: amu1= 0.5    ! cosine of emission angle (for hemispheric mean)

integer :: i, k, kp, ng, j1, j2, iflag

kp= kd+1

do k=1,kd
    do i=1,id
        if( abs(gfac(i,k)) > 1.0 ) then
            gfac(i,k)= 0.0
        endif
    enddo
enddo

beta(:,:)= 1.0 - gfac(:,:)

!            From table 1:  Hemispheric Mean Scheme
gama1(:,:)= 2.0 - sscat(:,:)*(1.0 + gfac(:,:))
gama2(:,:)= sscat(:,:)*beta(:,:)

!             Equations (21) and (22)
alambda(:,:)= sqrt( gama1(:,:)*gama1(:,:) -                    &
             gama2(:,:)*gama2(:,:)  )
grgama(:,:)= (gama1(:,:)-alambda(:,:)) / ( gama2(:,:) + eps )


!        The Planck function is assumed to vary linearly with tau
!                  within each layer

do i=1,id
    do k=1,kd
        if (tau(i,k) .gt. 0.00001) then
            b0(i,k)= bh(i,k)
            b1(i,k)= (bh(i,k+1)-bh(i,k))/tau(i,k)
        else
            b0(i,k)= 0.5*( bh(i,k)+bh(i,k+1) )
            b1(i,k)= 0.0
        endif
    enddo
enddo

!      Particular solutions to the 2-stream equations: Equation (27)
!      These are evaluated at the top and bottom of the Kth layer
!      hence c(u/d)(t/b) ==> cut, cub, cdt, cdb

denom(:,:)= 1.0 / ( gama1(:,:)+gama2(:,:) )

cut(:,:)= 2.0*pi*amu1*( b0(:,:) + b1(:,:) * denom(:,:) )
cdt(:,:)= 2.0*pi*amu1*( b0(:,:) - b1(:,:) * denom(:,:) )

cub(:,:)= 2.0*pi*amu1*( b0(:,:) + b1(:,:)* ( tau(:,:) +        &
                                            denom(:,:) )  )

cdb(:,:)= 2.0*pi*amu1*( b0(:,:) + b1(:,:)* ( tau(:,:) -        &
                                            denom(:,:) )  )

exfact(:,:) = exp( -alambda(:,:) * tau(:,:) )

!               Equation (44)
e1(:,:)= 1.0 + grgama(:,:) * exfact(:,:)
e2(:,:)= 1.0 - grgama(:,:) * exfact(:,:)
e3(:,:)= grgama(:,:) + exfact(:,:)
e4(:,:)= grgama(:,:) - exfact(:,:)

!      Formulate diagonal elements and rhs terms for 
!        the tridiagonal system of equations;
!      Odd terms (equation 41) and even terms equation (42)
!           with appropriate boundary conditions

adiag(:,1)  =   0.0
bdiag(:,1)  =  e1(:,1)
ddiag(:,1)  = -e2(:,1)
rhs(:,1)= -cdt(:,1)

do k=1,kd-1
    j1= 2*k+1
    j2= 2*k
    adiag(:,j1)= e2(:,k)*e3(:,k  ) - e4(:,k)*e1(:,k  )
    bdiag(:,j1)= e1(:,k)*e1(:,k+1) - e3(:,k)*e3(:,k+1)
    ddiag(:,j1)= e3(:,k)*e4(:,k+1) - e1(:,k)*e2(:,k+1)

    rhs(:,j1)= e3(:,k)*( cut(:,k+1)-cub(:,k  ) ) +               &
                             e1(:,k)*( cdb(:,k  )-cdt(:,k+1)   )

    adiag(:,j2)= e2(:,k+1)*e1(:,k  ) - e3(:,k  )*e4(:,k+1)
    bdiag(:,j2)= e2(:,k  )*e2(:,k+1) - e4(:,k  )*e4(:,k+1)
    ddiag(:,j2)= e1(:,k+1)*e4(:,k+1) - e2(:,k+1)*e3(:,k+1)

    rhs(:,j2)= e2(:,k+1)*( cut(:,k+1)-cub(:,k  ) ) +             &
                            e4(:,k+1)*( cdb(:,k  )-cdt(:,k+1)  )      

enddo

adiag(:,kd*2)= e1(:,kd) - (1.0-emis(:))*e3(:,kd)
bdiag(:,kd*2)= e2(:,kd) - (1.0-emis(:))*e4(:,kd)
ddiag(:,kd*2)= 0.0
rhs(:,kd*2)= emis(:)*pi*bsol(:) - cub(:,kd) +                   &
                     (1.0-emis(:))*cdb(:,kd)

!                 call tridiagonal solver
iflag = 0
call trdslv( kd*2, id, adiag, bdiag, ddiag, rhs, iflag )

y(:,:)= rhs(:,:)

!
! One can formulate the upward and downward fluxes
! from these terms:

! fuh(k)= e1(k)*y(2*k-1) + e2(k)*y(2*k) + ...

! fdh(k)= e3(k)*y(2*k-1) + e4(k)*y(2*k) + ...

! Instead, it is preferable to invoke the source
! function technique as described by Toon

! Here, we use the 2-stream solutions to evaluate the
! source term of the radiative transfer equation and then
! integrate this equation explicitly for a series of zenith
! angles:
! Finally, the upward and downward fluxes are evaluated
! by Gaussian quadrature


!           Formulate the following terms: Table 3 in Toon et al.

grg(:,:)= 1.0/amu1 - alambda(:,:)
grh(:,:)= grgama(:,:) * ( alambda(:,:) + 1.0/amu1 )
grj(:,:)= grh(:,:)
grk(:,:)= grg(:,:)

do k=1,kd
    grg(:,k)= grg(:,k)*( y(:,2*k-1) + y(:,2*k) )
    grh(:,k)= grh(:,k)*( y(:,2*k-1) - y(:,2*k) )
    grj(:,k)= grj(:,k)*( y(:,2*k-1) + y(:,2*k) )
    grk(:,k)= grk(:,k)*( y(:,2*k-1) - y(:,2*k) )
enddo

alpha1(:,:)= 2.0*pi*( b0(:,:) + b1(:,:)*( denom(:,:) - amu1) ) 

sigma1(:,:)= 2.0*pi*( b0(:,:) - b1(:,:)*( denom(:,:) - amu1) )

alpha2(:,:)= 2.0*pi*b1(:,:)
sigma2(:,:)= 2.0*pi*b1(:,:)

!      Note that alpha2 and  sigma2 are identical


!            First, integrate downwards:
!           Do these for gaussian zenith angles, xg(1:ngauss)

!         exfact (:,:) = exp( -tau(:,k)*alambda(:,k) );
!         exfactg(:,:) = exp( -tau(:,:)/xg(ng) );


do ng= 1, ngauss
    exfactg(:,:)= exp( -tau(:,:)/xg(ng) )
    gri(:,:,ng)=                                                 &
          grj(:,:)*( 1.0-exfact(:,:)*exfactg(:,:) ) /           &
                         ( alambda(:,:)*xg(ng)+1.0 )            &
       +  grk(:,:)*( exfactg(:,:)-exfact(:,:) )  /              &
                         ( alambda(:,:)*xg(ng)-1.0)             &
       +  sigma1(:,:)*( 1.0-exfactg(:,:) )                      &
       +  sigma2(:,:)*( xg(ng)*exfactg(:,:) + tau(:,:)-xg(ng) )
enddo

!    No incident flux at the top boundary
fdh(:,:)= 0.0
accum(:,1:ngauss)= 0.0

do k=1,kd
    do ng= 1, ngauss
        exfact2(:,ng)= exp( -tau(:,k)/xg(ng) )
        accum(:,ng)= accum(:,ng) * exfact2(:,ng) + gri(:,k,ng)
    enddo  
    !                 gaussian quadrature
    do ng= 1, ngauss
        fdh(:,k+1)= fdh(:,k+1) + xg(ng)*wt(ng)*accum(:,ng)
    enddo  

enddo

if( coszen > 0.0 ) then  !  solve for upward intensity at a specified zenith angle

    do ng= 1, 1
        exfactg(:,:)= exp( -tau(:,:)/coszen )
        gri(:,:,ng)=                                         &
            grh(:,:)*( 1.0-exfact(:,:)*exfactg(:,:) ) / &
                           ( alambda(:,:)*coszen+1.0 )      &
         +  grg(:,:)*( exfactg(:,:)-exfact(:,:) )  /    &
                          ( alambda(:,:)*coszen-1.0)        &
         +  alpha1(:,:)*( 1.0-exfactg(:,:) )            &
        +  alpha2(:,:)*( coszen - (tau(:,:)+coszen)*exfactg(:,:) )
    enddo

    !        apply a surface reflection boundary condition
    !       assume diffuse reflection of ir radiation at bottom boundary

    fuh(:,1:kd)= 0.0
    fuh(:,kp)= 2.0*pi*emis(:)*bsol + 2.0*(1.0-emis(:))*fdh(:,kp)

    ng= 1
    do k= kd, 1, -1
        exfact2(:,ng)= exp( -tau(:,k)/coszen )
        fuh(:,k) = fuh(:,k+1) * exfact2(:,ng) + gri(:,k,ng)
    enddo


else   ! ----------- Integrate upwards from the surface to obtain the upward flux component


!    Integrate upwards from the surface to obtain the upward flux component

    do ng=1,ngauss
        exfactg(:,:)= exp( -tau(:,:)/xg(ng) )
        gri(:,:,ng)=                                                &
              grh(:,:)*( 1.0-exfact(:,:)*exfactg(:,:) ) /           &
                            ( alambda(:,:)*xg(ng)+1.0 )             &
           +  grg(:,:)*( exfactg(:,:)-exfact(:,:) )  /              &
                            ( alambda(:,:)*xg(ng)-1.0)              &
           +  alpha1(:,:)*( 1.0-exfactg(:,:) )                      &
           +  alpha2(:,:)*( xg(ng) - (tau(:,:)+xg(ng))*exfactg(:,:) )
    enddo

    !        apply a surface reflection boundary condition
    !       assume diffuse reflection of ir radiation at bottom boundary

    fuh(:,1:kd) = 0.0
    fuh(:,kp)= pi*emis(:)*bsol(:) + (1.0-emis(:))*fdh(:,kp)

    do ng=1,ngauss
        accum(:,ng)= fuh(:,kp)/amu1
    enddo

    do k=kd,1,-1
        do ng=1,ngauss
            exfact2(:,ng)= exp( -tau(:,k)/xg(ng) )
            accum(:,ng)= accum(:,ng) * exfact2(:,ng) + gri(:,k,ng)
        enddo  
        !                 gaussian quadrature
        do ng= 1, ngauss
            fuh(:,k)= fuh(:,k) + xg(ng)*wt(ng)*accum(:,ng)
        enddo  
    enddo

end if 

return
end subroutine lw_scattering


!-----------------------------------------------------------------------

subroutine newton_damping( id, kd, tcol, qcol, heatra )
!
!  calculate simple newtonian damping
!
!    INPUT:
integer, intent(in) :: id, kd
real, dimension(id,kd), intent(in) :: tcol
real, dimension(kd),    intent(in) :: qcol

!    input/output:
real, dimension(id,kd), intent(inout) :: heatra

real, parameter :: tnot= 140.0  
real, parameter :: qnot= 0.9447012000e-05

integer :: k
real :: tbar, vfac

do k=1,3
    tbar= sum( tcol(:,k), 1 )/float(id)
    vfac= ( ( log(qcol(k))/log(qnot) )**6  )/5.0e4
    heatra(:,k)= heatra(:,k) - vfac*(tbar-tnot)
enddo

return
end subroutine newton_damping

!-----------------------------------------------------------------------

subroutine trdslv( n, mx2, a, b, c, x, iflg )
!
!  tridiagonal solver
!
!     if  iflg = 0  this subroutine computes the lu-decomposition
!     of the input matrix and then solves for the solution vector  x.

!     if  iflg = 1  the calculation of the lu-decomposition is skipped
!     and the solution vector  x  is calculated.
!

integer, intent(in) :: n, mx2, iflg
real, dimension(mx2,n), intent(in)    :: a
real, dimension(mx2,n), intent(inout) :: b, c   
real, dimension(mx2,n), intent(out)   :: x   

integer :: nm1, j     

nm1 = n-1

!   obtain the lu-decomposition
if( iflg .eq. 0 ) then

    b(:,1) = 1.0/b(:,1)
    c(:,1) = c(:,1)*b(:,1)

    do j=2,n
        b(:,j) = 1.0/( b(:,j)-a(:,j)*c(:,j-1) )
        c(:,j)= c(:,j)*b(:,j)
    enddo

end if

!        ----come here for back-solving----

x(:,1) = x(:,1)*b(:,1)
do j=2,n
    x(:,j) = ( x(:,j)-a(:,j)*x(:,j-1) )*b(:,j)
enddo

do j=nm1,1,-1
    x(:,j) = x(:,j)-c(:,j)*x(:,j+1)
enddo

return
end subroutine trdslv

!-----------------------------------------------------------------------

subroutine planck_table_eval( id, ntab, planck_frac, index, temp, planck )
!
!  evaluate planck table
!

integer, intent(in) :: id, ntab
integer, dimension(id),   intent(in)  :: index
real,    dimension(ntab), intent(in)  :: planck_frac
real,    dimension(id),   intent(in)  :: temp
real,    dimension(id),   intent(out) :: planck

integer :: i

do i=1,id 
    planck(i)= planck_frac(index(i)) * temp(i)**4
enddo 

return
end subroutine planck_table_eval

!#######################################################################

end module radiation_util_mod
