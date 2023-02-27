module mars_physics_update_mod

!!!!      use mpp_domains_mod, only: mpp_update_domains
!!!!      use mp_mod,          only: domain

!!!!      use fv_diagnostics_mod, only: id_prec
    use diag_manager_mod,   only: send_data
    use time_manager_mod,   only: time_type

    use constants_mod,    only: grav, kappa, rdgas, pi, radian

    use mars_physics_mod,   only: mars_physics, do_qmass

    use fms_mod,          only:  mpp_pe, mpp_root_pe, error_mesg, FATAL,       &
                             open_namelist_file, check_nml_error, &
                             mpp_pe, mpp_root_pe, close_file,     &
                             write_version_number, stdlog,        &
                             uppercase, read_data, write_data, field_size

    use fms2_io_mod, only: file_exists

    use field_manager_mod,  only: MODEL_ATMOS, parse, find_field_index

    use tracer_manager_mod, only: query_method, get_tracer_index,  &
                              get_number_tracers,get_tracer_names
    use initracer_mod



implicit none
!-----------------------------------------------------------------------

private

public :: mars_physics_update

   logical :: mcpu0 !debug only, identify master processor


contains

!-----------------------------------------------------------------------

 subroutine mars_physics_update(npx, npy, npz, is, ie, js, je, ng, nq,      &
                     u_dt, v_dt, t_dt, q_dt, ua, va, pt, q, phis, &
                     pe, delp, peln,                              &
                     dt, grid, ak, bk, qratio, rayf, master, &
                     p_ref, Time, time_total)


      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq
      real   , INTENT(IN   ) :: dt
      logical, INTENT(IN   ) :: rayf, master
      real   , INTENT(IN   ) :: grid(is-ng:ie+ng,js-ng:je+ng, 1:2)
      real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
      real   , INTENT(IN   ) :: phis(is-ng:ie+ng,js-ng:je+ng)

      type(time_type), intent(in) :: Time
      real, INTENT(IN), optional:: time_total

      real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)
      real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )

! Tendencies:
      real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)

      real, INTENT(INOUT):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: va(is-ng:ie+ng,js-ng:je+ng,npz)

      real, INTENT(OUT):: qratio(is:ie,js:je,npz)
      real, intent(in) :: p_ref

! Local

    integer  i, j, k, isw, iew, jsw, jew, nx_win, ny_win
    integer  iloc, jloc, beglat, nlev, nt_prog, ncnst 

    real  ::  rrg, ginv, tvm

    integer nco2,nar,nh2o,unit,io,ierr
    real, dimension(is-ng:ie+ng,js-ng:je+ng,npz, nq)::  q1     !inidvidual tracer mmr

    logical used

!       Need to keep track of layer mass tendencies; These need to be
!       polar filtered prior to use in update_fv_physics

  real, dimension(is:ie,js:je,npz) :: dmass

  real, dimension(is:ie,js:je,npz) :: dt_delp

  real p_full(is:ie, js:je, npz)
  real p_half(is:ie, js:je, npz+1)

  real z_full(is:ie, js:je, npz)
  real z_half(is:ie, js:je, npz+1)

  real rlon(is:ie, js:je)
  real rlat(is:ie, js:je)


  mcpu0 = (mpp_pe() == mpp_root_pe())

    ginv= 1.0/grav
    rrg = rdgas / grav        

    nlev= npz
    nt_prog= nq
    ncnst= nq

!-----------------------------------

!**      Label tracers - get tracer names

   if (do_qmass) then
    nco2   = find_field_index( MODEL_ATMOS, 'co2_mmr_gas' )
    nar    = find_field_index( MODEL_ATMOS, 'ar_mmr_gas' )
   endif

!-----------------------------------


#ifdef SKIP
    if( mpp_pe()==0 ) then 
     print *, 'Mars Phys: ',  is, ie, js, je
     print *, 'Latitudes: i= 1',  grid(1,:,2)*radian 
     print *, 'Latitudes: i= 5',  grid(5,:,2)*radian 
     print *, 'Latitudes: i=18',  grid(8,:,2)*radian 
    endif
#endif SKIP 

!        These are the "A" grid longitudes and latitudes
   rlon(is:ie,js:je)= grid(is:ie,js:je,1)
   rlat(is:ie,js:je)= grid(is:ie,js:je,2)


!          Calculate p_half, p_full, z_half, z_full
!  Note that pe and peln are defined with reversed j and k indices !

  DO j= js, je 
       do k=1,nlev+1
          do i=is,ie
             p_half(i,j,k) = pe(i,k,j)
          enddo
       enddo

       do k=1,nlev
          do i=is,ie
             p_full(i,j,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo

      do i=is,ie
          z_half(i,j,nlev+1) = phis(i,j) * ginv
      enddo

      do k=nlev,1,-1
          do i= is, ie
             tvm = rrg*pt(i,j,k)
             z_full(i,j,k) = z_half(i,j,k+1) + &
                            tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
             z_half(i,j,k) = z_half(i,j,k+1) + &
                            tvm*(peln(i,k+1,j)-peln(i,k,j))
          enddo
       enddo
  ENDDO



!    Distinguish between computational and Physics windows
!    Original code had the physics window ny_win = 1   (1 latitude at a time) 
!    I will make these the same. 


  nx_win=  ie + 1 - is 
  ny_win=  je + 1 - js
!!!!  ny_win=  1

  DO jsw= js, je, ny_win    !  ----------  Begin big loop over latitudes  -----------

!            print *, 'Entering j loop:', mpp_pe(), js, j, je 

     jew= jsw + (ny_win-1)
     jloc= jsw + 1 - js

     isw= is 
     iloc= 1

     call mars_physics( isw-is+1, jsw-js+1, dt, Time,          &
            rlon(is:ie,jsw:jew),      rlat(is:ie,jsw:jew),               &
          p_half(is:ie,jsw:jew,:),  p_full(is:ie,jsw:jew,:),             &
          z_half(is:ie,jsw:jew,:),  z_full(is:ie,jsw:jew,:),             &
          ua(is:ie,jsw:jew,1:nlev), va(is:ie,jsw:jew,1:nlev),            &
          pt(is:ie,jsw:jew,1:nlev),  q(is:ie,jsw:jew,1:nlev,1:ncnst),    &
          ua(is:ie,jsw:jew,1:nlev), va(is:ie,jsw:jew,1:nlev),            &
          pt(is:ie,jsw:jew,1:nlev),  q(is:ie,jsw:jew,1:nlev,1:nt_prog),  &
          u_dt(is:ie,jsw:jew,:), v_dt(is:ie,jsw:jew,:),                  &
          t_dt(is:ie,jsw:jew,:), q_dt(is:ie,jsw:jew,:,1:nt_prog),        &
          dmass(is:ie,jsw:jew,:), p_ref                   )


!-------------------------------------------
!-------------------------------------------
  if (do_qmass)   then ! mass of multiple tracers is on

!-------------------------------------------
!  -- Update/change tracers [mmr] --
!** NOTE ** fv_phys.F90 scales/normalizes ALL tracers by "qratio"
!-------------------------------------------
!   Accumulate adjustments to layer depths resulting from CO2 condensation/sublimation
! Calculate Mean Mass of Atmosphere
!          delp(i,j,k) = delp(i,j,k) - grav*dmass(i,beglat,k)
!-------------------------------------------


    q1(:,jsw:jew,:,nco2) = q(:,jsw:jew,:,nco2)
    q1(:,jsw:jew,:,nar) = q(:,jsw:jew,:,nar)


    do k = 1,nlev
       do i = is,ie
          q1(i,jsw:jew,k,nco2)    = q1(i,jsw:jew,k,nco2)-(dmass(i,jsw:jew,k)*(grav/delp(i,jsw:jew,k)))
          q1(i,jsw:jew,k,nar)     = q1(i,jsw:jew,k,nar)

          q_dt(i,jsw:jew,k,nco2)  = q_dt(i,jsw:jew,k,nco2)-((q(i,jsw:jew,k,nco2) - q1(i,jsw:jew,k,nco2))/dt)
          q_dt(i,jsw:jew,k,nar)   = q_dt(i,jsw:jew,k,nar)-((q(i,jsw:jew,k,nar) - q1(i,jsw:jew,k,nar))/dt)


          dt_delp(i,jsw:jew,k)    = delp(i,jsw:jew,k)*((q1(i,jsw:jew,k,nco2)-q(i,jsw:jew,k,nco2))+  &
                                                      (q1(i,jsw:jew,k,nar)-q(i,jsw:jew,k,nar)))

        enddo
     enddo

!-------------------------------------------
  else  ! total mass = surface pressure
!-------------------------------------------
!   Accumulate adjustments to layer depths resulting from CO2 condensation/sublimation
!          delp(i,j,k) = delp(i,j,k) - grav*dmass(i,beglat,k)
 
!     In the case where do_co2_condensation_cycle= F, dmass== 0.0
     do k= 1,nlev
        do i=is,ie
           dt_delp(i,jsw:jew,k) =  -grav*dmass(i,jsw:jew,k)
        enddo
     enddo

  endif ! end of qmass loop
!-------------------------------------------
!-------------------------------------------
  ENDDO    !   loop over latitudes




   do k=1,nlev
         do j= js, je
            do i= is, ie
               qratio(i,j,k)= 1.0 + dt_delp(i,j,k) / delp(i,j,k)
               delp(i,j,k)= delp(i,j,k) * qratio(i,j,k)  
            enddo
         enddo
    enddo

!-----------------------------------------------------------------------------


  end subroutine mars_physics_update
  
end module mars_physics_update_mod
