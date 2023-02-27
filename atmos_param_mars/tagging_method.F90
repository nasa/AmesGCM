module tagging_method_mod
! Module to manage tagged tracers
use constants_mod, only: PI,RADIUS  

use fms_mod, only: error_mesg, FATAL,                      &
                   open_namelist_file, check_nml_error,                &
                   mpp_pe, mpp_root_pe, close_file,                    &
                   write_version_number, stdlog,                       &
                   uppercase, read_data, write_data, field_size

use mars_surface_mod,  only:  sfc_frost_micro, sfc_roughness, sfc_topo, &
                              sfc_albedo,sfc_emiss,thermal_inertia 

implicit none

!------------------- interfaces ---------------------------------------
public ::  tagging_main


!----------------------------------------------------------------------
logical ::  mcpu0

contains

!#######################################################################
!#######################################################################
!#######################################################################

!======================================================================
subroutine tagging_main(methods_array, lat, lon, ndx, field2D, field3D, solrun, &
                    sources, inidust, taudust_curr, dustfield, pres, flag, stress, source)
!
!   Main routine : operations for adding tags 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index,get_field_info 
use tracer_manager_mod, only: query_method, get_tracer_index,  &
                                get_number_tracers, get_tracer_names

CHARACTER(len = *), DIMENSION(:)        :: methods_array
integer, intent(in)                      :: ndx
real,    intent(inout),  dimension(:,:), optional  :: field2D 
real,    intent(inout),  dimension(:,:,:), optional  :: field3D 
real,    intent(in),  dimension(:,:)     :: lat,lon 
real, intent(in), optional               :: solrun   ! current sol
real, intent(in), dimension(:,:), optional         :: sources 
real, intent(in), dimension(:,:), optional         :: inidust 
real, intent(in), dimension(:,:), optional         :: taudust_curr  ! tau_current (vis or ir)
real, intent(in),  dimension(:,:,:), optional  :: dustfield 
real, intent(in),  dimension(:,:,:), optional  :: pres 
CHARACTER(len = *), intent(in), optional ::flag
real,    intent(in),  dimension(:,:), optional  :: stress
real,    intent(in),  dimension(:,:), optional  :: source
!-----------------------------------------------------------------------
!   Local variables 
!-----------------------------------------------------------------------
integer  :: id, jd
character(len=1024) :: scheme, params
character (len=1024) :: fld_type, fld_name
integer :: model, num_methods, nbmax, nbmethods
integer :: i
real, dimension(size(lat,1),size(lat,2))   :: sum_field2D,ini_field2D,new_field2D
character(len=1) :: num
character(len=20) :: txt

mcpu0 = (mpp_pe() == mpp_root_pe()) 

id= size(lat,1); jd= size(lat,2)

!!! Read possible methods in methods_array
nbmethods=size(methods_array,1)


do i=1,nbmethods
    txt=trim(methods_array(i))
    if ((txt(1:3).ne."add").and.(txt(1:3).ne."cpl")) then
        if (query_method(txt, MODEL_ATMOS, ndx, scheme, params)) then
             if(present(field2D)) call listmethods2D(ndx, field2D, id, jd, lat, lon, txt, trim(scheme), trim(params), solrun, sources, inidust, flag, stress, source, taudust_curr)
             if(present(field3D)) call listmethods3D(ndx, field3D, id, jd, lat, lon, txt, trim(scheme), trim(params), taudust_curr, dustfield, pres, solrun, flag, stress, source)
        endif
    endif  ! neither add nor cpl methods
enddo

!********************************************
!!! Other possible sets of methods : ADDITION
!********************************************
if (present(field2D)) then
    sum_field2D(:,:)=0.
    ini_field2D(:,:)=field2D(:,:)
endif
nbmax=5

do i =1,nbmax
    write(num,'(i1)') i
    if (query_method('add'//trim(num), MODEL_ATMOS, ndx, scheme, params)) then
        if (present(field2D)) then
            new_field2D(:,:)=ini_field2D(:,:)
            call listmethods2D(ndx, new_field2D, id, jd, lat, lon, trim(scheme), trim(params), trim(params), solrun, sources, inidust, flag, stress, source, taudust_curr)
            sum_field2D=sum_field2D+new_field2D
            field2D(:,:)=min(sum_field2D,ini_field2D)
        endif
    endif 
enddo

!********************************************
!!! Other possible sets of methods : COUPLING
!********************************************
if (present(field2D)) then
    sum_field2D(:,:)=0.
    ini_field2D(:,:)=field2D(:,:)
endif
nbmax=5

do i =1,nbmax
    write(num,'(i1)') i
    if (query_method('cpl'//trim(num), model_atmos, ndx, scheme, params)) then
        if (present(field2d)) then
            new_field2d(:,:)=ini_field2d(:,:)
            call listmethods2d(ndx, new_field2d, id, jd, lat, lon, trim(scheme), trim(params), trim(params), solrun, sources, inidust, flag, stress, source)
            sum_field2d=sum_field2d+new_field2d
            where (sum_field2d(:,:).lt.i*field2d(:,:))
                field2d(:,:)=0.
            endwhere
        endif
    endif ! method coupling
enddo ! loop on cplx


end subroutine tagging_main

!======================================================================
!======================================================================
!======================================================================
subroutine listmethods2D(ndx, field, id, jd, lat, lon, method, scheme, &
                    params, solrun, sources, inidust, flag, stress, source, taudst)
!
!   parse 2d tag methods

integer, intent(in)                   :: ndx
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:)  :: field 
real,    intent(in),  dimension(:,:)  :: lat,lon 
character(len=*),  intent(in)  :: params, scheme
character(len=*),  intent(in)  :: method 
real, intent(in), optional               :: solrun   ! current sol
real, intent(in), dimension(:,:), optional         :: sources 
real, intent(in), dimension(:,:), optional         :: inidust 
character(len=*),  intent(in), optional  :: flag
real,    intent(in),  dimension(:,:), optional  :: stress
real,    intent(in),  dimension(:,:), optional  :: source
real,    intent(in),  dimension(:,:), optional  :: taudst

select case (trim(method))
    CASE ("geosource")
        call geo_source(ndx, field, id, jd, lat, lon, trim(params), stress,source, taudst,flag)
    CASE ("loctime")
        call loctime_source(ndx, field, solrun, id, jd, lat, lon, trim(params))
    CASE ("cutsource")
        call cutsource(ndx, field, trim(scheme), flag)
    CASE ("solsource")
        call solsource(ndx, field, solrun, id, jd, lat, lon, trim(params))
    CASE ("antigeosource")
        call anti_geo_source(ndx, field, id, jd, lat, lon, trim(params), stress,source)
    CASE ("custom_sources")
        call custom_sources(ndx, field, sources, id, jd, lat, lon, trim(scheme), trim(params))
    CASE ("inidust")
        call inidust_only(ndx, field, inidust, id, jd, lat, lon, trim(params))
    CASE DEFAULT
        write(*,*) 'BUG Method 2D ',trim(method),' is not implemented in tagging_method.F90'
        stop
end select

end subroutine listmethods2D

!======================================================================
!======================================================================
!======================================================================
subroutine listmethods3D(ndx, field, id, jd, lat, lon, method, scheme, params, taudust_curr, dustfield, pres &
                     , solrun, flag, stress, source)
!
!   parse 3d tag methods
!
integer, intent(in)                   :: ndx
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:,:)  :: field 
real,    intent(in),  dimension(:,:)  :: lat,lon 
character(len=*),  intent(in)  :: params, scheme 
character(len=*), intent(in) :: method
real, intent(in), dimension(:,:), optional         :: taudust_curr
real, intent(in),  dimension(:,:,:), optional  :: dustfield 
real, intent(in),  dimension(:,:,:), optional  :: pres 
real, intent(in), optional               :: solrun   ! current sol
character(len=*),  intent(in), optional  :: flag
real,    intent(in),  dimension(:,:), optional  :: stress
real,    intent(in),  dimension(:,:), optional  :: source

select case (trim(method))
    CASE ("envir")
        call environment(ndx, field, id, jd, lat, lon, trim(params), taudust_curr, dustfield, pres)
    CASE ("activatestorm") 
        call activate_storm(ndx, field, id, jd, lat, lon, solrun, pres, dustfield, trim(params),stress,source,taudust_curr,flag)
    CASE DEFAULT
        write(*,*) 'BUG Method 3D ',trim(method),' is not implemented in tagging_method.F90'
        stop
end select

end subroutine listmethods3D


!======================================================================
!======================================================================
!======================================================================
subroutine geo_source(ndx, field, id, jd, lat, lon, params,stress,source,taudst,flag)
!
!   Set source to 0 if outside required box (tag options)

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

integer, intent(in)                   :: ndx
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:)  :: field 
real,    intent(in),  dimension(:,:)  :: lat,lon 
character(len=*),  intent(in)  :: params 
real,    intent(in),  dimension(:,:), optional  :: stress
real,    intent(in),  dimension(:,:), optional  :: source
real,    intent(in),  dimension(:,:), optional  :: taudst
character(len=*),  intent(in), optional  :: flag 

!   Local variables 
!-----------------------------------------------------------------------
integer,parameter :: nbopt=12
character*5   :: geo_opt(nbopt) = (/ 'longi','latit','topog','emiss','albed','therm','rough','frost','stres','sourc','opaci','dista' /)
REAL :: geo_arr(nbopt*2)
real, dimension(size(lat,1),size(lat,2)) :: tab,dist,a,c,dlon,dlat
integer  l,i,j  
real    ::   value,activ,plat,plon

! Initialisations

mcpu0 = (mpp_pe() == mpp_root_pe()) 

geo_arr(:)=-10000    ! Table of min and max for each criterion

! If Flag is set to 0, we deactivate the tag
activ=1.
if (parse(params,trim(flag),value) == 1) then
   activ=value
endif

if (activ.eq.1) then
 !! Get point if present
 if (parse(params,trim('point'),value) == 1) then
   plat=value
   plon=value
   if (parse(params,trim('point'//'2'),value) == 1) then
          plon=value
   endif
   dlon(:,:) = plon*pi/180. - lon(:,:)
   dlat(:,:) = plat*pi/180. - lat(:,:)
   a(:,:) = (sin(dlat(:,:)/2.))**2 + cos(lat(:,:)) * cos(plat*pi/180.) * (sin(dlon(:,:)/2.))**2
   c(:,:) = 2. * atan2( sqrt(a(:,:)), sqrt(1.-a(:,:)) )
   dist(:,:) = RADIUS/1000. * c(:,:)
 endif
 !! Loop for all flags
 DO l=1,nbopt
     if (parse(params,trim(geo_opt(l)),value) == 1) then
      !if(mcpu0)  print *, 'TB18 params=',params
      !if(mcpu0)  print *, 'TB18 VAL=',value
      geo_arr(2*l-1)=value
      geo_arr(2*l)=value
      if (parse(params,trim(geo_opt(l)//'2'),value) == 1) then
          geo_arr(2*l)=value
      endif
      !if(mcpu0)  print *, 'TB18 VAL2=',value

      select case (l)
             CASE (1)
                tab=lon(:,:)*180./pi
             CASE (2)
                tab=lat(:,:)*180./pi
             CASE (3)
                tab=sfc_topo(:,:)
             CASE (4)
                tab=sfc_emiss(:,:)
             CASE (5)
                tab=sfc_albedo(:,:)
             CASE (6)
                tab=thermal_inertia(:,:)
             CASE (7)
                tab=sfc_roughness(:,:)
             CASE (8)
                tab=sfc_frost_micro(:,:,1)
             CASE (9)
                if (present(stress)) then
                   tab=stress(:,:)
                else
                   cycle
                endif
             CASE (10)
                if (present(source)) then
                   tab=source(:,:)
                else
                   cycle
                endif 
             CASE (11)
                if (present(taudst)) then
                   tab=taudst(:,:)
                else
                   cycle
                endif 
             CASE (12)
                tab=dist(:,:)
             CASE DEFAULT
                write(*,*) 'BUG in Tagging method Geosource: No surface properties selected'
                stop
      end select

      where((tab(:,:).gt.geo_arr(2*l)).or.(tab(:,:).lt.geo_arr(2*l-1))) 
                field(:,:)= 0.       
      endwhere

     endif
 ENDDO
endif

end subroutine geo_source

!======================================================================
!======================================================================
!======================================================================
subroutine loctime_source(ndx, field, solrun, id, jd, lat, lon, params)
!
!   Set source to 0 if outside required local time 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

integer, intent(in)                   :: ndx
real, intent(in)                   :: solrun   ! current sol
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:)  :: field 
real,    intent(in),  dimension(:,:)  :: lat,lon 
character(len=*),  intent(in) :: params 

!   Local variables 
!-----------------------------------------------------------------------
integer,parameter :: nbopt=1
character*4   :: loctime_opt(nbopt) = (/ 'time' /)
REAL :: loctime_arr(nbopt*2)
real, dimension(size(lat,1),size(lat,2)) :: tab
integer  l  
real    ::   value 

mcpu0 = (mpp_pe() == mpp_root_pe()) 
  
loctime_arr(:)=-10000.    ! Table of min and max for each criterion

do l=1,nbopt
    if (parse(params,trim(loctime_opt(l)),value) == 1) then
        loctime_arr(2*l-1)=value
        loctime_arr(2*l)=value
        if (parse(params,trim(loctime_opt(l)//'2'),value) == 1) then
            loctime_arr(2*l)=value
        endif

        select case (l)
            CASE (1)
                tab(:,:)=modulo( modulo(solrun+0.5,1.)*24. - (180.-modulo(lon(:,:)*180./pi,360.))*12./180.  , 24.)
            CASE DEFAULT
                write(*,*) 'BUG in Tagging method Loctime: wrong local time'
            stop
        end select
        where ((tab(:,:).gt.loctime_arr(2*l)).or.(tab(:,:).lt.loctime_arr(2*l-1))) 
            field(:,:)= 0.      
        endwhere
    endif
enddo

end subroutine loctime_source

!======================================================================
!======================================================================
!======================================================================
subroutine anti_geo_source(ndx, field, id, jd, lat, lon, params, stress, source)
!
!   Set source to 0 if IN the selected area 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

integer, intent(in)                   :: ndx
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:)  :: field 
real,    intent(in),  dimension(:,:)  :: lat,lon 
character(len=*),  intent(in) :: params 
real,    intent(in),  dimension(:,:), optional  :: stress,source

!   Local variables 
!
integer,parameter :: nbopt=10
character*5   :: geo_opt(nbopt) = (/ 'longi','latit','topog','emiss','albed','therm','rough','frost','stres','sourc' /)
REAL :: geo_arr(nbopt*2)
real, dimension(size(lat,1),size(lat,2)) :: tab
integer  l  
real    ::   value 

mcpu0 = (mpp_pe() == mpp_root_pe()) 

geo_arr(:)=-10000    ! Table of min and max for each criterion

! For each criteria, we get the range of values accepted
DO l=1,nbopt
    if (parse(params,trim(geo_opt(l)),value) == 1) then
        geo_arr(2*l-1)=value
        geo_arr(2*l)=value
        if (parse(params,trim(geo_opt(l)//'2'),value) == 1) then
            geo_arr(2*l)=value
        endif

        select case (l)
            CASE (1)
                tab=lon(:,:)*180./pi
            CASE (2)
                tab=lat(:,:)*180./pi
            CASE (3)
                tab=sfc_topo(:,:)
            CASE (4)
                tab=sfc_emiss(:,:)
            CASE (5)
                tab=sfc_albedo(:,:)
            CASE (6)
                tab=thermal_inertia(:,:)
            CASE (7)
                tab=sfc_roughness(:,:)
            CASE (8)
                tab=sfc_frost_micro(:,:,1)
            CASE (9)
                if (present(stress)) then
                    tab=stress(:,:)
                else
                    cycle
                endif
            CASE (10)
                if (present(source)) then
                    tab=source(:,:)
                else
                    cycle
                endif
            CASE DEFAULT
                write(*,*) 'TB18 BUG in antigeosource'
                stop
        end select

        where((tab(:,:).le.geo_arr(2*l)).and.(tab(:,:).gt.geo_arr(2*l-1)))
            field(:,:)= 0.         
        endwhere

    endif
ENDDO

end subroutine anti_geo_source

!======================================================================
!======================================================================
!======================================================================
subroutine solsource(ndx, field, solrun, id, jd, lat, lon, params)
!
!   Set source to 0 if outside required local time 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

integer, intent(in)                   :: ndx
real, intent(in)                   :: solrun   ! current sol
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:)  :: field 
real,    intent(in),  dimension(:,:)  :: lat,lon 
character(len=*),  intent(in) :: params 

!   Local variables 
!
integer,parameter :: nbopt=1
character*3   :: loctime_opt(nbopt) = (/ 'sol' /)
REAL :: loctime_arr(nbopt*2)
real, dimension(size(lat,1),size(lat,2)) :: tab
integer  l  
real    ::   value 

mcpu0 = (mpp_pe() == mpp_root_pe()) 
  
loctime_arr(:)=-10000.    ! Table of min and max for each criterion

do l=1,nbopt
    if (parse(params,trim(loctime_opt(l)),value) == 1) then
        loctime_arr(2*l-1)=value
        loctime_arr(2*l)=value
        if (parse(params,trim(loctime_opt(l)//'2'),value) == 1) then
            loctime_arr(2*l)=value
        endif

        select case (l)
            CASE (1)
                tab(:,:)=solrun
            CASE DEFAULT
                write(*,*) 'BUG in Tagging method Loctime: wrong local time'
                stop
        end select
        where ((tab(:,:).gt.loctime_arr(2*l)).or.(tab(:,:).lt.loctime_arr(2*l-1))) 
            field(:,:)= 0.      
        endwhere
    endif
enddo

end subroutine solsource

!======================================================================
!======================================================================
!======================================================================
subroutine cutsource(ndx, field, scheme, flag)
!
!   Set source to 0 if outside required local time 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

integer, intent(in)                   :: ndx
real,    intent(inout),  dimension(:,:)  :: field 
character(len=*),  intent(in)  :: scheme 
character(len=*),  intent(in), optional  :: flag 

character(len=1),  parameter :: semicolon         = ";"
integer :: control_array(10,2)
integer icount,l,ltrec
character(len=3)  :: valname
character(len=20)  :: control_str


mcpu0 = (mpp_pe() == mpp_root_pe())

ltrec= len_trim(scheme)
control_str=trim(adjustl(scheme))
control_array(:,1)=1
control_array(:,2)=ltrec
icount=1

do l= 1, ltrec ! Check where we have semicolon in the string
    if (control_str(l:l) == semicolon ) then
        control_array(icount,2) = l-1 ! End string
        control_array(icount+1,1) = l+1 ! Start new string
        icount = icount + 1
    endif
enddo


do l = 1,icount
    valname = trim(adjustl(control_str(control_array(l,1):control_array(l,2))))
    if (valname=='all' .or. valname==flag) then
        field(:,:)=0.
    endif
enddo

end subroutine cutsource


!======================================================================
!======================================================================
!======================================================================
subroutine custom_sources(ndx, field, sources, id, jd, lat, lon, scheme, params)
!
!   Set source to 0 if outside required local time 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

integer, intent(in)                   :: ndx
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:)  :: field 
real,    intent(in),  dimension(:,:)  :: lat,lon 
real,    intent(in),  dimension(:,:)  :: sources
character(len=*),  intent(in)  :: params,scheme 

!   Local variables 

real    ::   value,val1,val2 

mcpu0 = (mpp_pe() == mpp_root_pe())

val1=0.5
val2=1.

! Simple case

if (trim(scheme)=='pos') then
    where(sources(:,:).le.0.)
        field(:,:)=0.
    endwhere
elseif (trim(scheme)=='neg') then
    where(sources(:,:).gt.0.)
        field(:,:)=0.
    endwhere
elseif (trim(scheme)=='mul') then
    field(:,:)=field(:,:)*sources(:,:)
elseif (trim(scheme)=='add') then
    field(:,:)=field(:,:)+sources(:,:)
elseif (trim(scheme)=='dib') then
    field(:,:)=field(:,:)/sources(:,:)
elseif (trim(scheme)=='sub') then
    field(:,:)=field(:,:)-sources(:,:)
endif

! Case with values
if (parse(params,'thresh',value) == 1) then
    val1=value
    if (parse(params,'thresh'//'2',value) == 1) then
        val2=value
        where((sources(:,:).lt.val1).or.(sources(:,:).gt.val2))
            field(:,:)=0.
        endwhere
    else
        where(sources(:,:).lt.val1)
            field(:,:)=0.
        endwhere
    endif
endif

if (parse(params,'nothresh',value) == 1) then
    val1=value
    if (parse(params,'nothresh'//'2',value) == 1) then
        val2=value
        where((sources(:,:).gt.val1).and.(sources(:,:).lt.val2))
            field(:,:)=0.
        endwhere
    else
        where(sources(:,:).gt.val1)
            field(:,:)=0.
        endwhere
    endif
endif

end subroutine custom_sources


!======================================================================
!======================================================================
!======================================================================
subroutine inidust_only(ndx, field, inidust, id, jd, lat, lon, params)
!
!   Set source to 0 if outside required local time 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 
use tracer_manager_mod, only: query_method, get_tracer_index,  &
                                get_number_tracers, get_tracer_names

integer, intent(in)                   :: ndx
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:)  :: field 
real,    intent(in),  dimension(:,:)  :: lat,lon 
real,    intent(in),  dimension(:,:)  :: inidust
character(len=*),  intent(in)  :: params 

!   Local variables 

real    ::   value 

mcpu0 = (mpp_pe() == mpp_root_pe()) 

if (parse(params,"pos",value) == 1) then
    where(inidust(:,:).le.value)
        field(:,:)=0.
    endwhere
endif

if (parse(params,"neg",value) == 1) then
    where(inidust(:,:).gt.value)
        field(:,:)=0.
    endwhere
endif

end subroutine inidust_only


!======================================================================
!======================================================================
!======================================================================
subroutine environment(ndx, field, id, jd, lat, lon, params,taudst,dustfield,pres)
!
!   Set source to 0 if IN the selected area 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

integer, intent(in)                   :: ndx
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:,:)  :: field 
character(len=*),  intent(in) :: params 
real,    intent(in),  dimension(:,:)  :: lat,lon 
real,    intent(in),  dimension(:,:), optional  :: taudst
real,    intent(in),  dimension(:,:,:), optional  :: dustfield
real,    intent(in),  dimension(:,:,:), optional  :: pres

!   Local variables 

integer,parameter :: nbopt=2
character*5   :: env_opt(nbopt) = (/ 'opaci' , 'press' /)
LOGICAL :: case3D
REAL :: env_arr(nbopt*2)
real, dimension(size(field,1),size(field,2),size(field,3))   :: newfield, tab3D
real, dimension(size(field,1),size(field,2))   :: tab

integer  l,i,j,k,kd

real    ::   value 

kd=size(field,3)
mcpu0 = (mpp_pe() == mpp_root_pe()) 

env_arr(:)=-10000    ! Table of min and max for each criterion

! For each criteria, we get the range of values accepted
DO l=1,nbopt
    if (parse(params,trim(env_opt(l)),value) == 1) then
        env_arr(2*l-1)=value
        env_arr(2*l)=value
        if (parse(params,trim(env_opt(l)//'2'),value) == 1) then
            env_arr(2*l)=value
        endif

        select case (l)
            CASE (1)
                tab=taudst(:,:)
                case3D=.false.
            CASE (2)
                tab3D=pres(:,:,:)
                case3D=.true.
            CASE DEFAULT
                write(*,*) 'TB18 BUG in environment'
                stop
        end select


        if (case3D) then
            do i=1,id
                do j=1,jd
                    do k=1,kd
                        where((tab3D(i,j,:).ge.env_arr(2*l-1)).and.(tab3D(i,j,:).le.env_arr(2*l))) 
                            field(i,j,:)=dustfield(i,j,:)
                        endwhere
                    enddo
                enddo
            enddo
        else
            do i=1,id
                do j=1,jd
                    if((tab(i,j).ge.env_arr(2*l-1)).and.(tab(i,j).le.env_arr(2*l))) then
                        field(i,j,:)=dustfield(i,j,:)
                    endif
                enddo
            enddo
        endif

   endif
ENDDO

end subroutine environment

!======================================================================
!======================================================================
!======================================================================
subroutine activate_storm(ndx, field, id, jd, lat, lon, solrun, pres, dustfield, params,stress,source,taudst,flag)
!
!   Set source to 0 if IN the selected area 

use  field_manager_mod, only: MODEL_ATMOS, parse, find_field_index 

integer, intent(in)                   :: ndx
integer, intent(in)                   :: id, jd
real,    intent(inout),  dimension(:,:,:)  :: field 
real, intent(in)                   :: solrun   ! current sol
character(len=*),  intent(in) :: params 
real,    intent(in),  dimension(:,:)  :: lat,lon 
real,    intent(in),  dimension(:,:,:), optional  :: dustfield
real,    intent(in),  dimension(:,:,:), optional  :: pres
real,    intent(in),  dimension(:,:), optional  :: stress
real,    intent(in),  dimension(:,:), optional  :: source
real,    intent(in),  dimension(:,:), optional  :: taudst
character(len=*),  intent(in), optional  :: flag 

!   Local variables 
integer,parameter :: nbopt=13
character*5   :: geo_opt(nbopt) = (/ 'longi','latit','topog','emiss','albed','therm','rough','frost','stres','sourc','opaci','dista','press' /)
integer,parameter :: nbops=2
character*3   :: ops_opt(nbops) = (/ 'mul','add' /)
LOGICAL :: case3D
REAL :: geo_arr(nbopt*2)
REAL :: ops_arr(nbops)
REAL :: sol_arr(2)
real, dimension(size(field,1),size(field,2),size(field,3))   :: tmpfield3D, tab3D
real, dimension(size(field,1),size(field,2))   :: tmpfield2D
real, dimension(size(lat,1),size(lat,2)) :: tab,dist,a,c,dlon,dlat
real :: nb2d,nb3d,activ,plat,plon,mulval,addval
integer  l,i,j,k,kd
real    ::   value 

kd=size(field,3)
mcpu0 = (mpp_pe() == mpp_root_pe()) 

geo_arr(:)=-10000.    ! Table of min and max for each criterion
sol_arr(:)=-10000.    ! Table of min and max for each criterion
dist(:,:)=100.

! If Flag is set to 0, we deactivate the tag
activ=1.
if(mcpu0)  print *, 'TB21 par=',params    
if (parse(params,trim(flag),value) == 1) then
   activ=value
endif
! If sol is outside required sols for storm activation, we deactivate the tag
if (parse(params,trim('sol'),value) == 1) then
   sol_arr(1)=value
   sol_arr(2)=value
   if (parse(params,trim('sol'//'2'),value) == 1) then
      sol_arr(2)=value
   endif
   if(mcpu0)  print *, 'TB21 SOLVAL=',sol_arr    
   if ((solrun.gt.sol_arr(2)).or.(solrun.lt.sol_arr(1))) then
      activ=0.      
   endif
endif

if(mcpu0)  print *, 'TB21 ACTIV=',activ    
if (activ.eq.1.) then
 if(mcpu0)  print *, 'TB21 ACTIVATION='   
 tmpfield2D=0.
 tmpfield3D=0.
 nb2d=0.
 nb3d=0.
 !! Loop for all flags
 DO l=1,nbopt
     if(mcpu0)  print *, 'TB21 Loop=',l   
     if (parse(params,trim(geo_opt(l)),value) == 1) then
      geo_arr(2*l-1)=value
      geo_arr(2*l)=value
      if (parse(params,trim(geo_opt(l)//'2'),value) == 1) then
          geo_arr(2*l)=value
      endif

      case3D=.false.
      select case (l)
             CASE (1)
                tab=lon(:,:)*180./pi
             CASE (2)
                tab=lat(:,:)*180./pi
             CASE (3)
                tab=sfc_topo(:,:)
             CASE (4)
                tab=sfc_emiss(:,:)
             CASE (5)
                tab=sfc_albedo(:,:)
             CASE (6)
                tab=thermal_inertia(:,:)
             CASE (7)
                tab=sfc_roughness(:,:)
             CASE (8)
                tab=sfc_frost_micro(:,:,1)
             CASE (9)
                if (present(stress)) then
                   tab=stress(:,:)
                else
                   cycle
                endif
             CASE (10)
                if (present(source)) then
                   tab=source(:,:)
                else
                   cycle
                endif 
             CASE (11)
                if (present(taudst)) then
                   tab=taudst(:,:)
                else
                   cycle
                endif 
             CASE (12)
                !! Get point if present
                if (parse(params,trim('point'),value) == 1) then
                  plat=value
                  plon=value
                  if (parse(params,trim('point'//'2'),value) == 1) then
                    plon=value
                  endif
                  dlon(:,:) = plon*pi/180. - lon(:,:)
                  dlat(:,:) = plat*pi/180. - lat(:,:)
                  a(:,:) = (sin(dlat(:,:)/2.))**2 + cos(lat(:,:)) * cos(plat*pi/180.) * (sin(dlon(:,:)/2.))**2
                  c(:,:) = 2. * atan2( sqrt(a(:,:)), sqrt(1.-a(:,:)) )
                  dist(:,:) = RADIUS/1000. * c(:,:)
                  if(mcpu0)  print *, 'TB21 dist=',dist    
                endif
                tab=dist(:,:)
             CASE (13)
                tab3D=pres(:,:,:)
                case3D=.true.
             CASE DEFAULT
                write(*,*) 'BUG in Tagging method active storm: No properties selected'
                stop
      end select

      if (case3D) then
            nb3d=nb3d+1.
            do i=1,id
                do j=1,jd
                        where((tab3D(i,j,:).ge.geo_arr(2*l-1)).and.(tab3D(i,j,:).le.geo_arr(2*l))) 
                            tmpfield3D(i,j,:)=tmpfield3D(i,j,:)+1.
                        endwhere
                enddo
            enddo
      else
            nb2d=nb2d+1.
            where((tab(:,:).le.geo_arr(2*l)).and.(tab(:,:).ge.geo_arr(2*l-1))) 
                tmpfield2D(:,:)= tmpfield2D(:,:)+1.       
            endwhere
      endif
      !if(mcpu0)  print *, 'TB21 tmp2D=',tmpfield2D    
      !if(mcpu0)  print *, 'TB21 tmp3D=',tmpfield3D    
      if(mcpu0)  print *, 'TB21 nbs=',nb2d,nb3d    

     endif ! parse value
 ENDDO ! loop geo opt

 ! Defining the type of operation
 mulval=1.
 addval=0.
 DO l=1,nbops
     if (parse(params,trim(ops_opt(l)),value) == 1) then
      ops_arr(l)=value
      select case (l)
             CASE (1)
                mulval=value
             CASE (2)
                addval=value
             CASE DEFAULT
                write(*,*) 'BUG in Tagging method active storm: No properties selected'
                stop
      end select
     endif ! parse value
 ENDDO ! loop ops opt

 do k=1,kd
    where(tmpfield3D(:,:,k).eq.nb3d.and.tmpfield2D(:,:).eq.nb2d)
      field(:,:,k)=dustfield(:,:,k)*mulval+addval
    endwhere
 enddo
endif ! activ

end subroutine activate_storm


end module tagging_method_mod












