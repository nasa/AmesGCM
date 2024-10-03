module coagulation_mod

use constants_mod, only: pi,grav,avogno,wtmco2,rdgas,kbz=>kboltz
use initracer_mod
use fms_mod, only: error_mesg, FATAL, file_exist,                      &
                   open_namelist_file, check_nml_error,                &
                   mpp_pe, mpp_root_pe, close_file,                    &
                   write_version_number, stdlog,                       &
                   uppercase, read_data, write_data, field_size

implicit none
!private

!!! Constants
! kbz = 1.38064852e-23
! dynvis0=14.8e-6 ! dyn viscosity CO2 at 293.15 K (Pa.s)
! avogno : 6.02214076e23 # mol-1
! wtmco2 mco2=44.01e-3 ! molar mass co2 kg/mol
! mco20=wtmco2/avogno ! mass of one air molecule
! grav 3.711 # gravity

!------------------- interfaces ---------------------------------------
public ::  coagul_init

! Namelist --------------------
integer :: coal_kg=1  ! coalescence option for gravitational coagulation
                    ! 0: default, E=1, 1= sticking coefficent as in Jacobson 2005
logical :: kernel_b=.true. ! brownian
logical :: kernel_g=.true. ! gravitational
logical :: kernel_de=.true. ! diffusion enhancement
logical :: kernel_ti=.false. ! turbulent inertia motion
real    ::  effb=1.  !  Sticking efficiency for B coagulation
real    ::  coal_fac=1.  ! Overall coalescence factor


namelist /coagulation_nml/ coal_kg,kernel_b,kernel_g,kernel_de,kernel_ti,effb,coal_fac

! Other --------------------
logical ::  mcpu0
logical :: module_is_initialized = .false.
integer,parameter :: table_numt = 15
integer,parameter :: table_nump = 25
integer,parameter :: table_numm = 25
real, save :: table_temp(table_numt)
real, save :: table_pres(table_nump)
real, save :: table_mfp(table_numm)
real, save :: table_b(table_numt,table_nump,nres_coag,nres_coag)
real, save :: table_g(table_numt,table_nump,nres_coag,nres_coag)
real, save :: table_de(table_numt,table_nump,nres_coag,nres_coag)
real, save :: table_ti(table_numt,table_nump,nres_coag,nres_coag)

contains

!#######################################################################
!#######################################################################
!#######################################################################

!======================================================================
subroutine coagul_main(is, js, lon, lat, dt, temp, p_half, p_full, rini, rdt_coag)

integer, intent(in)  :: is, js
real,    intent(in)  :: dt
real, intent(in),    dimension(:,:)     :: lon
real, intent(in),    dimension(:,:)     :: lat
real, intent(in),    dimension(:,:,:)   :: temp
real, intent(in),    dimension(:,:,:)   :: p_half, p_full
real, intent(in),    dimension(:,:,:,:)   :: rini
real, intent(inout), dimension(size(rini,1),size(rini,2),size(rini,3),size(rini,4)) :: rdt_coag

!======================================================================
!   Main routine : operations for adding tags 

!-----------------------------------------------------------------------
!   Local variables 
!-----------------------------------------------------------------------
integer  :: ie, je, id, jd, kd, i, j, k, l, ii, jj, ndx,nt
real dens,dev,cst,sig0
real mfp0,temp0,rad0,rad10,rad20,rho0,pres0
real term1,term2,kernel,norm
real, dimension(size(temp,1),size(temp,2),size(temp,3)) :: r0,rho, rn, rn_tot, mtot, mtotnew
real, dimension(size(temp,1),size(temp,2),size(temp,3),ndust_mass) :: rn_new,rn_ini,rm_ini,rm_new
real, dimension(size(temp,1),size(temp,2),size(temp,3),nres_coag) :: ndis,ndis_new 
real, dimension(size(temp,1),size(temp,2),size(temp,3),nres_coag,ndust_mass) :: ndis_tab,rat_tab,ndis_tab_new
!for boxinterp
integer :: t1,t2,p1,p2,m1,m2

mcpu0 = (mpp_pe() == mpp_root_pe()) 

id= size(temp,1); jd= size(temp,2); kd= size(temp,3)
ie= is + id - 1
je= js + jd - 1

!-----------------------------------------------------------------------
!*** 1) Initialisations
!-----------------------------------------------------------------------

!!! Particles properties 
dens    = dpden_dt  ! density dust 
sig0     = dev_dt    ! Standard deviation of the dust distribution
cst     = .75 / (pi*dens) * exp( -4.5*sig0**2. )   ! fixed distribution

!!! Atmospheric state
do k=1,kd
     rho(:,:,k) = p_full(:,:,k) / ( rdgas*temp(:,:,k) )
enddo 

!!! Radius and Volumes defined in initracer/initmicro
! rads_coag, vols_coag, vrat_coag

!!! Initial tracers
rdt_coag(:,:,:,:)=0.d0
ndis_tab(:,:,:,:,:)=0.d0
ndis_new(:,:,:,:)=0.d0
mtot(:,:,:)=0.d0

do nt=1,ndust_mass
     ndx= dust_mass_indx(nt)
     if (ndx .le. nt_nontag) then ! do not take into account tags
       rm_ini(:,:,:,nt)=max(rini(:,:,:,ndx),1e-15)
       rn(:,:,:)=max(rini(:,:,:,ndx+1),1e-15)

       !!! Get total mass to ensure mass conservation
       mtot(:,:,:)=mtot(:,:,:)+rm_ini(:,:,:,nt)

       !!! Initial Lognormal distribution
       ! sig0 remains constant
       r0(:,:,:)=(3/4.*rm_ini(:,:,:,nt)/(pi*rn(:,:,:)*dens))**(1/3.)*exp(-1.5*sig0**2)

       do i=1,nres_coag ! Sum of the ndis for all dust modes
          ndis_tab(:,:,:,i,nt)=rn(:,:,:)*rho(:,:,:)*deltar_coag(i)/(2.*rads_coag(i)*(2*pi)**0.5*sig0)*exp(-0.5*(log(rads_coag(i)/r0(:,:,:)))**2/sig0**2)
       enddo
       rn_ini(:,:,:,nt)=sum(ndis_tab(:,:,:,:,nt),4)/rho(:,:,:)
       rn_new(:,:,:,nt)=rn_ini(:,:,:,nt)
     endif
enddo

do nt=1,ndust_mass
     ndx= dust_mass_indx(nt)
     if (ndx .le. nt_nontag) then ! do not take into account tags
       ! Initial Ratio for each mode
       rat_tab(:,:,:,:,nt)=ndis_tab(:,:,:,:,nt)/sum(ndis_tab(:,:,:,:,:),5)
     endif
enddo

! Normalization of the ratio (make sure sum=1)
do nt=1,ndust_mass
     rat_tab(:,:,:,:,nt)=rat_tab(:,:,:,:,nt)/sum(rat_tab(:,:,:,:,:),5)
enddo
! Total initial distribution
ndis(:,:,:,:)=sum(ndis_tab(:,:,:,:,:),5)
rn_tot(:,:,:)=sum(ndis(:,:,:,:),4)/rho(:,:,:)


!-----------------------------------------------------------------------
!*** 2) Coagulation
!-----------------------------------------------------------------------
! Spatial loop => 1D
#ifdef FULLCOAG
do i=is,ie
  do j=js,je
    do l= 1, kd  
      if (rn_tot(i,j,l).gt.1000.) then
       temp0=temp(i,j,l)
       mfp0=mfp(temp0,rho(i,j,l))

       ! Bin loop
       do k=1,nres_coag

         ! Term 1
         term1=0.
         do jj=1,k
          do ii=1,k-1
             kernel=0.
             if (kernel_b) kernel=kernel+betab(temp0,mfp0,rads_coag(ii),rads_coag(jj))
             if (kernel_g) kernel=kernel+betag(temp0,mfp0,rads_coag(ii),rads_coag(jj))
             if (kernel_de) kernel=kernel+betade(temp0,mfp0,rho(i,j,l),rads_coag(ii),rads_coag(jj))
             if (kernel_ti) kernel=kernel+betati(temp0,mfp0,rho(i,j,l),rads_coag(ii),rads_coag(jj))
             term1=term1+frac(ii,jj,k)*vols_coag(ii)*kernel*ndis_new(i,j,l,ii)*ndis(i,j,l,jj)
          enddo
         enddo
         term1=term1*dt/vols_coag(k)

         ! Term 2
         term2=0.
         do jj=1,nres_coag
            kernel=0.
            if (kernel_b) kernel=kernel+betab(temp0,mfp0,rads_coag(k),rads_coag(jj))
            if (kernel_g) kernel=kernel+betag(temp0,mfp0,rads_coag(k),rads_coag(jj))
            if (kernel_de) kernel=kernel+betade(temp0,mfp0,rho(i,j,l),rads_coag(k),rads_coag(jj))
            if (kernel_ti) kernel=kernel+betati(temp0,mfp0,rho(i,j,l),rads_coag(k),rads_coag(jj))
            term2=term2+(1-frac(k,jj,k))*kernel*ndis(i,j,l,jj)
         enddo
         term2=term2*dt

         ndis_new(i,j,l,k)=(ndis(i,j,l,k)+coal_fac*term1)/(1.+coal_fac*term2)

       enddo

       norm=sum(ndis(i,j,l,:)*vols_coag(:))/sum(ndis_new(i,j,l,:)*vols_coag(:))
       ndis_new(i,j,l,:)=ndis_new(i,j,l,:)*norm

      endif
    enddo
  enddo
enddo
#else
! Spatial loop => 1D
do i=is,ie
  do j=js,je
    do l= 1, kd  
      if (rn_tot(i,j,l).gt.1000.) then
       temp0=temp(i,j,l)
       mfp0=mfp(temp0,rho(i,j,l))
       pres0=p_full(i,j,l)
!find the corners for interpolation       
       t2=findval(table_temp,temp0)
       if (t2==0) t2=table_numt
       if (t2==1) t2=2
       t1=t2-1
       m2=findval(table_mfp,mfp0)
       if (m2==0) m2=table_numm
       if (m2==1) m2=2
       m1=m2-1
       p2=findval(table_pres,pres0)
       if (p2==0) p2=table_numt
       if (p2==1) p2=2
       p1=p2-1
       
       ! Bin loop
       do k=1,nres_coag

         ! Term 1
         term1=0.
         do jj=1,k
          do ii=1,k-1
             kernel=0.
             if (kernel_b) kernel=kernel+boxinterp(table_b(t1,m1,ii,jj),       &
                                                   table_b(t1,m2,ii,jj),       &
                                                   table_b(t2,m1,ii,jj),       &
                                                   table_b(t2,m2,ii,jj),       &
                                                   table_temp(t1),table_temp(t2),log10(table_mfp(m1)),log10(table_mfp(m2)),temp0,log10(mfp0))
             if (kernel_g) kernel=kernel+boxinterp(table_g(t1,m1,ii,jj),       &
                                                   table_g(t1,m2,ii,jj),       &
                                                   table_g(t2,m1,ii,jj),       &
                                                   table_g(t2,m2,ii,jj),       &
                                                   table_temp(t1),table_temp(t2),log10(table_mfp(m1)),log10(table_mfp(m2)),temp0,log10(mfp0))
             if (kernel_de) kernel=kernel+boxinterp(table_de(t1,p1,ii,jj),    &
                                                   table_de(t1,p2,ii,jj),     &
                                                   table_de(t2,p1,ii,jj),     &
                                                   table_de(t2,p2,ii,jj),     &
                                                   table_temp(t1),table_temp(t2),log10(table_pres(p1)),log10(table_pres(p2)),temp0,log10(pres0))
!             if (kernel_ti) kernel=kernel+betati(temp0,mfp0,rho(i,j,l),ii,jj)  !not implemented yet
             term1=term1+frac(ii,jj,k)*vols_coag(ii)*kernel*ndis_new(i,j,l,ii)*ndis(i,j,l,jj)
          enddo
         enddo
         term1=term1*dt/vols_coag(k)

         ! Term 2
         term2=0.
         do jj=1,nres_coag
            kernel=0.
            if (kernel_b) kernel=kernel+boxinterp(table_b(t1,m1,k,jj),       &
                                                  table_b(t1,m2,k,jj),       &
                                                  table_b(t2,m1,k,jj),       &
                                                  table_b(t2,m2,k,jj),       &
                                                  table_temp(t1),table_temp(t2),log10(table_mfp(m1)),log10(table_mfp(m2)),temp0,log10(mfp0))
             if (kernel_g) kernel=kernel+boxinterp(table_g(t1,m1,k,jj),       &
                                                   table_g(t1,m2,k,jj),       &
                                                   table_g(t2,m1,k,jj),       &
                                                   table_g(t2,m2,k,jj),       &
                                                   table_temp(t1),table_temp(t2),log10(table_mfp(m1)),log10(table_mfp(m2)),temp0,log10(mfp0))
             if (kernel_de) kernel=kernel+boxinterp(table_de(t1,p1,k,jj),       &
                                                   table_de(t1,p2,k,jj),       &
                                                   table_de(t2,p1,k,jj),       &
                                                   table_de(t2,p2,k,jj),       &
                                                   table_temp(t1),table_temp(t2),log10(table_pres(p1)),log10(table_pres(p2)),temp0,log10(pres0))
!            if (kernel_ti) kernel=kernel+betati(temp0,mfp0,rho(i,j,l),k,jj)  !not implemented yet
            term2=term2+(1-frac(k,jj,k))*kernel*ndis(i,j,l,jj)
         enddo
         term2=term2*dt

         ndis_new(i,j,l,k)=(ndis(i,j,l,k)+coal_fac*term1)/(1.+coal_fac*term2)
       enddo

      !! Volume conservation
      norm=sum(ndis(i,j,l,:)*vols_coag(:))/sum(ndis_new(i,j,l,:)*vols_coag(:))
      ndis_new(i,j,l,:)=ndis_new(i,j,l,:)*norm

      endif ! rn_tot

    enddo
  enddo
enddo
#endif FULLCOAG
!-----------------------------------------------------------------------
!*** 3) new moments mass and number
!-----------------------------------------------------------------------
mtotnew(:,:,:)=1.d-15
do nt=1,ndust_mass
    ndx= dust_mass_indx(nt)
    if (ndx .le. nt_nontag) then ! do not take into account tags
       !! Estimate new distribution in each mode
       ndis_tab_new(:,:,:,:,nt)=ndis_new(:,:,:,:)*rat_tab(:,:,:,:,nt)     

       do i=is,ie
        do j=js,je
         do l= 1, kd  
          !! New number and mass mixing ratio in each mode
          rn_new(i,j,l,nt)=sum(ndis_tab_new(i,j,l,:,nt))/rho(i,j,l)
          rm_new(i,j,l,nt)=sum(ndis_tab_new(i,j,l,:,nt)*vols_coag(:))*dens/rho(i,j,l)
         enddo
        enddo
       enddo

       !! New total mass
       mtotnew(:,:,:)=mtotnew+rm_new(:,:,:,nt)
    endif
enddo

do nt=1,ndust_mass
    ndx= dust_mass_indx(nt)
    if (ndx .le. nt_nontag) then ! do not take into account tags
       !! mass conservation
       rm_new(:,:,:,nt)=rm_new(:,:,:,nt)/mtotnew(:,:,:)*mtot(:,:,:)
       !! New tendancy
       rdt_coag(:,:,:,ndx)=(rm_new(:,:,:,nt)-rm_ini(:,:,:,nt))/dt
       rdt_coag(:,:,:,ndx+1)=(rn_new(:,:,:,nt)-rn_ini(:,:,:,nt))/dt
       !! check
       where (rini(:,:,:,ndx)+rdt_coag(:,:,:,ndx)*dt.le.0.)
         rdt_coag(:,:,:,ndx)=-rini(:,:,:,ndx)/dt+1.e-15
       end where
       where (rini(:,:,:,ndx+1)+rdt_coag(:,:,:,ndx+1)*dt.le.0.)
         rdt_coag(:,:,:,ndx+1)=-rini(:,:,:,ndx+1)/dt+1.e-14
       end where
    endif
enddo


end subroutine coagul_main


!======================================================================
!======================================================================
real function dynvis(temp)
  implicit none
  ! dynamic viscosity following sutherland's formula for CO2
  ! Pa.s = kg m-1 s-1
  real, intent(in) :: temp
  dynvis=1.37e-5*(273.15 + 222.)/(temp + 222.)*(temp/273.15)**(3./2.)
end function dynvis

real function kinvis(temp,rho)
  implicit none
  real, intent(in) :: temp,rho
  kinvis=dynvis(temp)/rho
end function kinvis

real function thervel(temp,massm)
  implicit none
  ! thermal velocity of an air molecule or a dust particle
  ! m s-1
  real, intent(in) :: temp, massm
  thervel=(8*kbz*temp/(pi*massm))**0.5
end function thervel

real function mfp(temp,rho)
  implicit none
  ! thermal velocity of an air molecule 
  ! m s-1
  real, intent(in) :: temp, rho
  mfp=2.*kinvis(temp,rho)/thervel(temp,mco20) ! mean free path (m)
end function mfp

real function cun(temp,mfp,rad)
  implicit none
  !! Cunningham slip flow correction
  real, intent(in) :: temp, mfp, rad
  cun=1.+mfp/rad*(1.246+0.42*exp(-0.87/(mfp/rad)))
end function cun

real function diff(temp,mfp,rad)
  implicit none
  !! diffusion coefficient
  real, intent(in) :: temp, mfp, rad
  diff=kbz*temp/(6.*pi*rad*dynvis(temp))*cun(temp,mfp,rad)
end function diff

real function lambda(temp, mfp, rad)
  implicit none
  !! particle mean free path 
  real, intent(in) :: temp,mfp,rad
  real massm,dens
  dens    = dpden_dt  ! density dust 
  massm=dens*4./3.*pi*rad**3
  lambda=2.*diff(temp, mfp, rad)/(pi*thervel(temp,massm))
end function lambda

real function delta(temp,mfp,rad)
  implicit none
  !! mean distance 
  real, intent(in) :: temp,mfp,rad
  real :: ltmp
  
  ltmp=lambda(temp,mfp,rad)
  delta=( (2.*rad+ltmp)**3 - (4.*rad**2+ltmp**2)**1.5 ) /(6*rad*ltmp) - 2*rad
end function delta

real function fallv(temp,mfp,rad)
  implicit none
  !! fall velocity
  real, intent(in) :: temp,mfp,rad
  real dens
  dens    = dpden_dt  ! density dust 
  fallv=2./9.*rad**2*dens*grav/dynvis(temp)*cun(temp,mfp,rad)
end function fallv

real function reyn(temp,mfp,rad,rho)
  implicit none
  !! Reynold number
  real, intent(in) :: temp,rad,mfp,rho
  reyn=2.*rad*fallv(temp,mfp,rad)/kinvis(temp,rho)
end function reyn

real function stokes(temp,mfp,rad1,rad2)
  implicit none
  !! Stokes number
  real, intent(in) :: temp,mfp,rad1,rad2
  real :: f1
  
  f1=fallv(temp,mfp,rad1)
  stokes=f1*abs(fallv(temp,mfp,rad2)-f1)/(rad2*grav)
end function stokes

real function schmidt(temp,mfp,rho,rad)
  implicit none
  !! Schmidt number
  real, intent(in) :: temp,mfp,rho,rad
  schmidt=kinvis(temp,rho)/diff(temp,mfp,rad)
end function schmidt

real function coal(temp,mfp,rad1,rad2,rho)
  implicit none
  !! Coalescence efficiency
  real, intent(in) :: temp,mfp,rad1,rad2,rho
  real :: coalea,coalev
  real :: stmp,rtmp
  
  stmp=stokes(temp,mfp,rad1,rad2)
  rtmp=reyn(temp,mfp,rad2,rho)
  coalea=stmp**2/(stmp+0.5)**2
  if (stmp>1.214) then
    coalev=(1+0.75*log(2.*stmp)/(stmp-1.214))**(-2)
  else
    coalev=0.
  endif
  coal=(60.*coalev+coalea*rtmp)/(60.+rtmp)
end function coal

real function betab(temp,mfp,rad1,rad2)
  implicit none
  real, intent(in) :: temp,mfp,rad1,rad2
  real :: m1,m2,num,den1,den2
  real :: dens
  real :: d1,d2
  
  dens    = dpden_dt  ! density dust 
  d1=diff(temp,mfp,rad1)
  d2=diff(temp,mfp,rad2)
  !! kernel transition regime
  ! Mass of 1 particle 1 and 2
  m1=dens*4./3.*pi*rad1**3
  m2=dens*4./3.*pi*rad2**3
  num=4.*pi*(rad1+rad2)*(d1+d2)
  den1=(rad1+rad2)/(rad1+rad2+(delta(temp,mfp,rad1)**2+delta(temp,mfp,rad2)**2)**0.5)
  den2=4./effb*(d1+d2)/((thervel(temp,m1)**2+thervel(temp,m2)**2)**0.5*(rad1+rad2))
  betab=num/(den1+den2)
end function betab

real function betag(temp,mfp,rad1,rad2)
  implicit none
  real, intent(in) :: temp,mfp,rad1,rad2
  real w1,w2,ee
  real :: radt
  !! kernel gravitation
  ! fall vel
  w1=fallv(temp,mfp,rad1)
  w2=fallv(temp,mfp,rad2)
  radt=(rad1+rad2)**2
  if (coal_kg.eq.0) then ! E=1
    ee=1.
  elseif (coal_kg.eq.1) then
    ee=1.5*(min(rad1,rad2))**2/radt
  elseif (coal_kg.eq.2) then
    ee=0.25*(min(rad1,rad2))**2/radt
  endif
  !b=coal(temp,r1,r2)*pi*(r1+r2)**2*abs(w1-w2)
  betag=ee*pi*radt*abs(w1-w2)
end function betag

real function betade(temp,mfp,rho,rad1,rad2)
  implicit none
  real, intent(in) :: temp,mfp,rad1,rad2,rho
  real rd,sc
  !! kernel diffusion enhancement
  rd=reyn(temp,mfp,max(rad1,rad2),rho)
  sc=schmidt(temp,mfp,rho,min(rad1,rad2))
  if (rd.le.1.) then
    betade=betab(temp,mfp,rad1,rad2)*0.45*rd**(1/3.)*sc**(1/3.)
  else
    betade=betab(temp,mfp,rad1,rad2)*0.45*rd**(1/2.)*sc**(1/3.)
  endif
end function betade

real function betade_tabfunc(temp,pres,rad1,rad2)
  implicit none
  real, intent(in) :: temp,pres,rad1,rad2
  real rd,sc,rho,mf
  !! kernel diffusion enhancement
  rho=pres/(temp*rdgas)
  mf=mfp(temp,rho)
  rd=reyn(temp,mf,max(rad1,rad2),rho)
  sc=schmidt(temp,mf,rho,min(rad1,rad2))
  if (rd.le.1.) then
    betade_tabfunc=betab(temp,mf,rad1,rad2)*0.45*rd**(1/3.)*sc**(1/3.)
  else
    betade_tabfunc=betab(temp,mf,rad1,rad2)*0.45*rd**(1/2.)*sc**(1/3.)
  endif
end function betade_tabfunc

real function betati(temp,mfp,rho,rad1,rad2)
  implicit none
  real, intent(in) :: temp,mfp,rad1,rad2,rho
  real w1,w2,eps
  ! fall vel
  w1=fallv(temp,mfp,rad1)
  w2=fallv(temp,mfp,rad2)
  betati=eps**(0.75)*pi/(grav*kinvis(temp,rho)**(0.25))*(rad1+rad2)**2*abs(w1-w2)
end function betati

real function betati_tabfunc(temp,pres,rad1,rad2)
  implicit none
  real, intent(in) :: temp,rad1,rad2
  real w1,w2,eps,mf,rho,pres
  rho=pres/(temp*rdgas)
  mf=mfp(temp,rho)
  ! fall vel
  w1=fallv(temp,mf,rad1)
  w2=fallv(temp,mf,rad2)
  betati_tabfunc=eps**(0.75)*pi/(grav*kinvis(temp,rho)**(0.25))*(rad1+rad2)**2*abs(w1-w2)
end function betati_tabfunc

real function frac(i,j,k)
  implicit none
  integer, intent(in) :: i,j,k
  real vint
  ! intermediate volume
  vint=vols_coag(i)+vols_coag(j)
  if (k.lt.nres_coag) then
     if ( (vint.ge.vols_coag(k)) .and. (vint.lt.vols_coag(k+1))) then
         frac=(vols_coag(k+1)-vint)/(vols_coag(k+1)-vols_coag(k))*vols_coag(k)/vint
         return
     endif
  endif
  if ((k.eq.nres_coag) .and. (vint.ge.vols_coag(k))) then
     frac=1.
     return
  endif
  if (k.gt.1) then
     if ( (vint.lt.vols_coag(k)) .and. (vint.gt.vols_coag(k-1))) then
        frac=(-vols_coag(k-1)+vint)/(-vols_coag(k-1)+vols_coag(k))*vols_coag(k)/vint
        return
     endif
  endif

  frac=0.
  return
end function frac

! *********************************************************
! *********************************************************

subroutine coagul_init()

!! Local
integer  unit, io, ierr 

! *********************************************************
!     ----- read namelist /dust_update_nml/   -----
! *********************************************************

if (file_exist('input.nml')) then
      unit = open_namelist_file ( )
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=coagulation_nml, iostat=io, end=10)
         ierr = check_nml_error (io, 'coagulation_nml')
      enddo
10     call close_file (unit)
endif

if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=coagulation_nml)

#ifndef FULLCOAG
call make_tables()
#endif

module_is_initialized  = .true.

end subroutine coagul_init

! *********************************************************
! *********************************************************

subroutine make_tables()
implicit none
integer :: i,j,k,l
real :: t1,t2,p1,p2

mcpu0 = (mpp_pe() == mpp_root_pe())

do i = 1,table_numt
    table_temp(i)=50.+25.*real(i)
end do
do j = 1,table_nump
    table_pres(j)=(1.e-8)*10.**(real(j)/2.)
end do
do j = 1,table_numm
    table_mfp(j)=(1.e-8)*10.**(real(j)/2.)
end do
do i = 1,table_numt
    do j = 1,table_numm
        do k = 1,nres_coag
            do l = 1,nres_coag
                table_b(i,j,k,l)=betab(table_temp(i),table_mfp(j),rads_coag(k),rads_coag(l))
                table_g(i,j,k,l)=betag(table_temp(i),table_mfp(j),rads_coag(k),rads_coag(l))
            end do
        end do
    end do
end do
    

do i = 1,table_numt
    do j = 1,table_nump
        do k = 1,nres_coag
            do l = 1,nres_coag
                table_de(i,j,k,l)=betade_tabfunc(table_temp(i),table_pres(j),rads_coag(k),rads_coag(l))
!                table_ti(i,j,k,l)=betati_tabfunc()
            end do
        end do
    end do
end do

if (mcpu0) print*,'Finished making coagulation tables: table_b(10,10,4,:) =',table_b(10,10,4,:)
end subroutine make_tables

!======================================================================C
!======================================================================C

real function boxinterp(ft1p1,ft2p1,ft1p2,ft2p2,tref1,tref2,pref1,  &
         pref2,tmid,pmid)
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

real :: ft1p1, ft2p1, ft1p2, ft2p2, tref1, tref2
real :: pref1, pref2, tmid, pmid, ans1, ans2, ans

!======================================================================C

ans1 = ft1p1 + (ft2p1 - ft1p1)*(tmid - tref1)/(tref2 - tref1)
ans2 = ft1p2 + (ft2p2 - ft1p2)*(tmid - tref1)/(tref2 - tref1)
boxinterp  = ans1 + (ans2 - ans1)*(pmid - pref1)/(pref2 - pref1)

end function boxinterp

!======================================================================C
!======================================================================C

integer function findval(array,value)
implicit none
real :: array(:)
real :: value

findval=minloc(array,dim=1,mask = (array > value))
    
end function findval

end module coagulation_mod






