module ames_pbl_interface
! module to interface FV3 to ames PBL module        
use pblmod_mgcm, only : pbl_driver,initpbl
use field_manager_mod,  only: MODEL_ATMOS, parse, find_field_index
use diag_manager_mod, only: send_data
use constants_mod, only: seconds_per_day,grav
!use vert_turb_driver_mod, only: id_z_pbl
use time_manager_mod, only: time_type
use initracer_mod
implicit none
private
public :: ames_pbl

contains

subroutine ames_pbl(is,js,dt,ps,p_half,p_full,z_half,akm,akh,  &
                  u,v,t,tsurf,q,snow,frost,tdtlw,tau_x,tau_y,&
                  sens,evap,&
                  udt,vdt,tdt,qdt, &
                  udt_pbl,vdt_pbl,tdt_pbl,qdt_pbl,&
                  ustar,thstar,cdm,cdh,dsens_datm,dsens_dsrf, &
                  wind,z_pbl,p_pbl,pbl_lev,rkh_out,Time_next)

implicit none
integer, intent(in) :: is, js
type(time_type), intent(in) :: Time_next
real, intent(in) :: dt
real, intent(in),    dimension(:,:,:)   :: p_half, &                                ! layer interface pressure [Pa]
                                        p_full, &                                   ! layer midpoint pressure [Pa]
                                        z_half, &                                   ! layer interface heights [m]
                                        tdtlw                                       ! radiative heating [K/s]
real, intent(in),    dimension(:,:,:)   :: u, &                                     ! u wind [m/s]
                                        v, &                                        ! v wind [m/s]
                                        t                                           ! temperature [K]
real, intent(in),    dimension(:,:,:,:)   :: q                                      ! tracers [*/kg]
real, intent(in),    dimension(size(t,1),size(t,2)) :: tsurf, &                     ! surface temperature
                                                    ps                              ! surface pressure
real, intent(in), dimension(size(t,1),size(t,2),size(t,3))  :: tdt, &               ! temperature tendency [K/s]
                                                            udt, &                  ! u wind tendency [m/s/s]
                                                            vdt                     ! v wind tendency [m/s/s]
real, intent(in), dimension(:,:,:,:)  :: qdt                                        ! tracer tendency [*/kg/s]
real, intent(inout), dimension(size(t,1),size(t,2)) :: snow, &                      ! co2 surface ice
                                                    frost                           ! h2o surface ice
integer, intent(out), dimension(size(t,1),size(t,2)) :: pbl_lev                     ! pbl top index
real, intent(out),   dimension(:,:,:)   :: akm, &                                   ! momentum mixing coefficient
                                        akh                                         ! heat mixing coefficient
real, intent(out), dimension(size(t,1),size(t,2),size(t,3)) :: tdt_pbl, &           ! pbl temperature tendency [K/s]
                                                            udt_pbl, &              ! pbl u wind tendency [m/s/s]
                                                            vdt_pbl                 ! pbl v wind tendency [m/s/s]
real, intent(out), dimension(size(q,1),size(q,2),size(q,3),size(q,4)) :: qdt_pbl    ! pbl tracer tendency [*/kg/s]
real, intent(out), dimension(size(t,1),size(t,2),2*size(t,3)+1) :: rkh_out          ! Output eddy mixing coefficient (m2/s)
                                                    
! local variables
real,   dimension(size(t,1),size(t,2)) :: ustar,ustar2,thstar,thstar2,cdm,cdm2,cdh,cdh2
real,   dimension(size(t,1),size(t,2)) :: sens, evap, sens_tnd, &
                                     evap_tnd, tau_x, tau_y, dsens_datm, dsens_dsrf, &
                                     wind, z_pbl,p_pbl

logical :: polarcap     !    If over residual cap
!===============================================================
!     Output Arguments from pbl_driver
!===============================================================
real :: qgtend  !    Surface water ice tendency (kg/m2/s)
!===============================================================
!     Optional Arguments for pbl_driver
!===============================================================
real :: rhouch      !    Output Rho*Cp*Cdh*Ustar
real :: strx        !    Output horizontal wind stress in X direction (zonal)
real :: stry        !    Output horizontal wind stress in Y direction (meridional)
real :: sup         !    Output upwards water vapor flux (kg/kg/s)
real :: sdn         !    Output downwards water vapor flux (kg/kg/s)
real*4 :: latheat   !    Output latent heat flux (W/m2) for diagnostics only
!    The following two variables are calculated at boundaries
!     down to the top boundary of the second layer from the surface
real :: du_out(2*size(t,3)+1)   !    Output shear (/s)
real :: rkm_out(2*size(t,3)+1)  !    Output eddy mixing coefficient (m2/s)

!===============================================================
!     Local Variables
!===============================================================	
!    local, updated fields
real, dimension(size(t,1),size(t,2),size(t,3)) :: u_bl, v_bl, t_bl, qrad
real, dimension(size(q,1),size(q,2),size(q,3),ntrace_gas) :: q_bl
integer :: nz           !    Number of layer midpoints
real, dimension(2*size(t,3)+3) :: tl    !    Absolute Temperature (K) at boundaries and midpoints
real, dimension(2*size(t,3)+3) :: pl    !    Pressure (Pa) at boundaries and midpoints
integer :: ie, je
integer :: i,j,k,n,n2,nt,ndx
integer :: nh2o,nvar
real :: sens_old, evap_old, dsens_olda, dsens_olds
real, dimension(size(t,3)) :: tttmp,uttmp,vttmp
real, dimension(size(q,3),ntrace_gas) :: qttmp
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: zhalf_ag
real ::tautmpx,tautmpy
real, dimension(size(t,3)+1) :: rkh_avg
real, parameter :: kcrit = 0.01
logical :: used

ie = size(t,1)
je = size(t,2)
nz = size(t,3)

qdt_pbl=0.
tdt_pbl=0.
udt_pbl=0.
vdt_pbl=0.
qrad=0.  !pass in 0 radiative heating rate since it is included in tdt
do k=1,nz+1
    zhalf_ag(:,:,k)=z_half(:,:,k)-z_half(:,:,nz+1)
enddo

cdm2 = cdm
thstar2 = thstar
cdh2 = cdh
ustar2 = ustar

rkh_avg(:)=0.

!convert cdh to cdh in newpbl
cdh2(:,:) = cdh(:,:)/sqrt(cdm(:,:))
!convert thstar to thstar in newpbl
thstar2(:,:) = thstar(:,:)/(-grav/tsurf(:,:))

call initpbl(size(t,3),dt,du_out)

! Grid vertical level of the PBL height
pbl_lev(:,:)=nz
p_pbl(:,:)=p_full(:,:,nz)

t_bl=t+tdt*dt
u_bl=u+udt*dt
v_bl=v+vdt*dt

!! Filling tracer field for gas type only
do nt= 1, ntrace_gas
    ndx= gas_indx(nt)
    q_bl(:,:,:,nt)=q(:,:,:,ndx)+qdt(:,:,:,ndx)*dt
enddo
! Number of gas variables to be mixed in PBL
nvar=ntrace_gas+3 ! u,v,heat

rkh_out(:,:,:) = 0.
do i=1,ie
    do j=1,je
        tttmp(:) = 0.
        uttmp(:) = 0.
        vttmp(:) = 0.
        qttmp(:,:) = 0.
        rkm_out(:) = 0.
        du_out(:) = 0.
        pl(3) = p_half(i,j,1)
        pl(2) = pl(3)*0.5
        pl(1) = pl(2)*1.e-6
        do k= 1, nz
            n= 2*k + 2
            tl(n)= t_bl(i,j,k)         
            pl(n)= p_full(i,j,k)
            pl(n+1)=p_half(i,j,k+1)
            if (k.eq.nz) then
                tl(n+1) = 0.5*( tsurf(i,j)+t_bl(i,j,k) )
            else
                tl(n+1)= 0.5*( t_bl(i,j,k+1)+t_bl(i,j,k) )
            end if
        enddo 
        tl(3) = tl(4)

        polarcap = frost(i,j) .gt. 1.e-5
        sens_old = sens(i,j)
        evap_old = 0.d0      
        dsens_olda = dsens_datm(i,j)
        dsens_olds = dsens_dsrf(i,j)

        call pbl_driver(nvar,nz,tl,pl,u_bl(i,j,:),v_bl(i,j,:),q_bl(i,j,:,:), &
                   frost(i,j),  &
                   tsurf(i,j),  &
                   ps(i,j), &
                   dt,evap_old,  &
                   qrad(i,j,:),  &
                   polarcap,tttmp(:),uttmp(:), &
                   vttmp(:),qttmp(:,:),qgtend, &
                   rhouch,sens_old,tautmpx,tautmpy,rkh_out(i,j,:),rkm_out,sup,sdn, &
                   du_out,latheat,snow(i,j),ustar2(i,j),thstar2(i,j),cdm2(i,j), &
                   cdh2(i,j),dsens_olda,dsens_olds )
        do k=2,nz
            n2=2*k-1   
            akh(i,j,k) = rkh_out(i,j,n2)
            akm(i,j,k) = rkm_out(n2)
            rkh_avg(k) = 0.5*(akh(i,j,k)+akm(i,j,k))
        end do

        tau_x(i,j)=tautmpx
        tau_y(i,j)=tautmpy

        sens(i,j) = sens_old
        dsens_datm(i,j) = dsens_olda
        dsens_dsrf(i,j) = dsens_olds

        if (rkh_avg(nz).gt.kcrit) then
            k=nz
            do while (k.gt.2 .and. rkh_avg(k-1).gt.kcrit)
                k=k-1
            enddo
            pbl_lev(i,j)=k-1
            z_pbl(i,j) = zhalf_ag(i,j,k)+ &
            (zhalf_ag(i,j,k-1)-zhalf_ag(i,j,k))* &
            (rkh_avg(k)-kcrit) / &
            (rkh_avg(k)-rkh_avg(k-1))
        else
            z_pbl(i,j) = 0.
        endif

        p_pbl(i,j)=p_full(i,j,pbl_lev(i,j))

        cdm(i,j) = cdm2(i,j)
        cdh(i,j) = cdh2(i,j)*sqrt(cdm2(i,j))
        ustar(i,j)=ustar2(i,j)
        thstar(i,j)=thstar2(i,j)*(-grav/tsurf(i,j))
        wind(i,j) = ustar(i,j)/sqrt(cdm(i,j))

        udt_pbl(i,j,:) = uttmp(:)
        vdt_pbl(i,j,:) = vttmp(:)
        tdt_pbl(i,j,:) = tttmp(:)

        do nt= 1, ntrace_gas
            ndx= gas_indx(nt)
            qdt_pbl(i,j,:,ndx) = qttmp(:,nt)
        enddo

    end do
end do
!if (id_z_pbl > 0) then
!    used = send_data(id_z_pbl , z_pbl , Time_next, is, js )
!endif

end subroutine ames_pbl

end module ames_pbl_interface
