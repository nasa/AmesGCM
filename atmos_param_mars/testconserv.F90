module testconserv_mod
! Module to test tracer conservation

use constants_mod, only: grav
use field_manager_mod,  only: MODEL_ATMOS, parse, find_field_index
implicit none


public ::  comparetrac,testneg,checkconserv,checkco2

contains

!****************************************************************************************
!****************************************************************************************
!****************************************************************************************
!****************************************************************************************

subroutine comparetrac(is,ie,js,je,nz,field,ind1,ind2,incert,flag)
!
!     Purpose
!     
!     compare values in a  4D field for different tracers
!  
!
!

!     Arguments
INTEGER, intent(in) :: nz,is,ie,js,je,ind1,ind2
REAL, intent(in), dimension(:,:,:,:) :: field
REAL, intent(in) :: incert
character (len=11), intent(in), optional :: flag
! LOCAL VARIABLES
INTEGER l,i,j
!-----------------------------------------------------------------------

do i=is,ie
    do j=js,je
        do l = 1, nz
            IF ((field(i,j,l,ind1)-field(i,j,l,ind2)).lt.0.d0) THEN
                print*, 'TB18 compare issue in ',flag,' :',i,j,l,field(i,j,l,ind1),field(i,j,l,ind2)
                stop
            ENDIF
        enddo
    enddo
enddo
end subroutine comparetrac

!****************************************************************************************
!****************************************************************************************
!****************************************************************************************
!****************************************************************************************

subroutine testneg(is,ie,js,je,nz,field,ind,incert,flag)
!
!     Purpose
!     
!     Test if tracer field is lower than incert
!
!

!     Arguments
INTEGER, intent(in) :: nz,is,ie,js,je,ind
REAL, intent(in), dimension(:,:,:,:) :: field
REAL, intent(in) :: incert
character (len=11), intent(in), optional :: flag
! LOCAL VARIABLES
INTEGER l,i,j
!-----------------------------------------------------------------------

do i=is,ie
    do j=js,je
        do l = 1, nz
            if (field(i,j,l,ind).lt.incert) then
                print*, 'TB18 Negative field in ',flag,' :',i,j,l,field(i,j,l,ind)
                stop
            endif
        enddo
    enddo
enddo

end subroutine testneg


subroutine ckconchem(is,ie,js,je,nz,pres0,r0,surfwat0,surfh2o20,pres1,r1,surfwat1,surfh2o21,flag)
!==================================================================
!     Purpose
!     -------
!     check conservation all tracers in the column (mass)
!  
!==================================================================
!-----------------------------------------------------------------------
!     Arguments
integer, intent(in) :: nz,is,ie,js,je
real, intent(in), dimension(:,:,:) :: pres0,pres1
real, intent(in), dimension(:,:,:,:) :: r0,r1
real, intent(in), dimension(:,:,:) :: surfwat0,surfwat1
real, intent(in), dimension(:,:) :: surfh2o20,surfh2o21
character (len=4), intent(in), optional :: flag
! LOCAL VARIABLES
integer l,i,j,nh2o,nice,nice_surf,id,jd
integer nco2,nco,no,no1d,no2,no3,  &
        nh,nh2,noh,nho2,nh2o2,nn2
real :: colh2o20,colh2o21,colwat0,colwat1,atmwat0,atmwat1,atmh2o20,atmh2o21
real :: incert_h2o2,incert_wat
!-----------------------------------------------------------------------

!! Choice tracer water and dust NDX
nh2o   = find_field_index( MODEL_ATMOS, 'vap_mass_micro' )
nice= find_field_index( MODEL_ATMOS, 'ice_mass_micro' )
nh2o2  = find_field_index( MODEL_ATMOS, 'h2o2_mmr_gas' )

!! Choices index for surface water and H2O2 
nice_surf=1
!ndst_surf=2
incert_wat=1.e-12
incert_h2o2=1.e-12
do i=is,ie
    id=i-is+1
    do j=js,je
        jd=j-js+1
        colh2o20=surfh2o20(id,jd)
        colh2o21=surfh2o21(i,j)
        colwat0=surfwat0(id,jd,nice_surf)
        colwat1=surfwat1(i,j,nice_surf)

        atmwat0=0.
        atmwat1=0.
        atmh2o20=0.
        atmh2o21=0.
        do l = 1, nz
            atmwat0=atmwat0+r0(i,j,l,nh2o)*(pres0(i,j,l+1)-pres0(i,j,l))/grav
            atmwat0=atmwat0+r0(i,j,l,nice)*(pres0(i,j,l+1)-pres0(i,j,l))/grav
            atmwat1=atmwat1+r1(i,j,l,nh2o)*(pres1(i,j,l+1)-pres1(i,j,l))/grav
            atmwat1=atmwat1+r1(i,j,l,nice)*(pres1(i,j,l+1)-pres1(i,j,l))/grav

            atmh2o20=atmh2o20+r0(i,j,l,nh2o2)*(pres0(i,j,l+1)-pres0(i,j,l))/grav
            atmh2o21=atmh2o21+r1(i,j,l,nh2o2)*(pres1(i,j,l+1)-pres1(i,j,l))/grav
        enddo
        colh2o20=colh2o20+atmh2o20
        colh2o21=colh2o21+atmh2o21
        colwat0=colwat0+atmwat0
        colwat1=colwat1+atmwat1

        if (abs(colh2o21-colh2o20).gt.incert_h2o2) THEN
            print*, 'AB20 conserv issue H2O2 in ',flag,' :',i,j,colh2o21,colh2o20,abs(colh2o21-colh2o20)
            print*, 'AB20 atm',atmh2o21,atmh2o20,abs(atmh2o21-atmh2o20)
            print*, 'AB20 surf',surfh2o21(i,j),surfh2o20(id,jd),abs(surfh2o21(i,j)-surfh2o20(id,jd))
            !stop
         endif
         if (abs(colwat1-colwat0).gt.incert_wat) THEN
            print*, 'AB20 conserv issue WAT in ',flag,' : ',i,j,colwat1,colwat0,abs(colwat1-colwat0)
            print*, 'AB20 atm',atmwat1,atmwat0,abs(atmwat1-atmwat0)
            print*, 'AB20 surf',surfwat1(i,j,nice_surf),surfwat0(id,jd,nice_surf),abs(surfwat1(i,j,nice_surf)-surfwat0(id,jd,nice_surf))
            !stop
         endif 
    enddo
enddo
end subroutine ckconchem


subroutine checkconserv(is,ie,js,je,nz,pres0,r0,surfwat0,surfdust0,pres1,r1,surfwat1,surfdust1,flag)
!
!     Purpose
!     
!     check conservation all tracers in the column (mass)
!  
!

!     Arguments
INTEGER, intent(in) :: nz,is,ie,js,je
REAL, intent(in), dimension(:,:,:) :: pres0,pres1
REAL, intent(in), dimension(:,:,:,:) :: r0,r1
REAL, intent(in), dimension(:,:,:) :: surfwat0,surfwat1
REAL, intent(in), dimension(:,:,:) :: surfdust1,surfdust0
character (len=7), intent(in), optional :: flag
! LOCAL VARIABLES
INTEGER l,i,j,nh2o,nice,ndst,ncor,ndst_surf,nice_surf,id,jd
REAL :: coldust0,coldust1,colwat0,colwat1,atmwat0,atmwat1,atmdust0,atmdust1
REAL :: incert_dst,incert_wat
!-----------------------------------------------------------------------

!! Choice tracer water and dust NDX
nh2o= find_field_index( MODEL_ATMOS, 'vap_mass_micro' )
nice= find_field_index( MODEL_ATMOS, 'ice_mass_micro' )
ncor= find_field_index( MODEL_ATMOS, 'cor_mass_micro' )
ndst= find_field_index( MODEL_ATMOS, 'dst_mass_micro' )


!! Choices index for surface water and dust NT
ndst_surf=1
nice_surf=1
!ndst_surf=2
incert_wat=1.e-12
incert_dst=1.e-12
do i=is,ie
    id=i-is+1
    do j=js,je
        jd=j-js+1
        coldust0=surfdust0(id,jd,ndst_surf)
        coldust1=surfdust1(i,j,ndst_surf)
        colwat0=surfwat0(id,jd,nice_surf)
        colwat1=surfwat1(i,j,nice_surf)

        atmwat0=0.
        atmwat1=0.
        atmdust0=0.
        atmdust1=0.
        do l = 1, nz
            atmwat0=atmwat0+r0(i,j,l,nh2o)*(pres0(i,j,l+1)-pres0(i,j,l))/grav
            atmwat0=atmwat0+r0(i,j,l,nice)*(pres0(i,j,l+1)-pres0(i,j,l))/grav
            atmwat1=atmwat1+r1(i,j,l,nh2o)*(pres1(i,j,l+1)-pres1(i,j,l))/grav
            atmwat1=atmwat1+r1(i,j,l,nice)*(pres1(i,j,l+1)-pres1(i,j,l))/grav

            atmdust0=atmdust0+r0(i,j,l,ndst)*(pres0(i,j,l+1)-pres0(i,j,l))/grav
            atmdust0=atmdust0+r0(i,j,l,ncor)*(pres0(i,j,l+1)-pres0(i,j,l))/grav
            atmdust1=atmdust1+r1(i,j,l,ndst)*(pres1(i,j,l+1)-pres1(i,j,l))/grav
            atmdust1=atmdust1+r1(i,j,l,ncor)*(pres1(i,j,l+1)-pres1(i,j,l))/grav
        enddo
        coldust0=coldust0+atmdust0
        coldust1=coldust1+atmdust1
        colwat0=colwat0+atmwat0
        colwat1=colwat1+atmwat1

        IF (abs(coldust1-coldust0).gt.incert_dst) THEN
            print*, 'TB18 conserv issue DUST in',flag,' :',i,j,coldust1,coldust0,abs(coldust1-coldust0)
            print*, 'TB18 atm',atmdust1,atmdust0,abs(atmdust1-atmdust0)
            print*, 'TB18 surf',surfdust1(i,j,ndst_surf),surfdust0(id,jd,ndst_surf),abs(surfdust1(i,j,ndst_surf)-surfdust0(id,jd,ndst_surf))
        ENDIF
        IF (abs(colwat1-colwat0).gt.incert_wat) THEN
            print*, 'TB18 conserv issue WAT in ',flag,' : ',i,j,colwat1,colwat0,abs(colwat1-colwat0)
            print*, 'TB18 atm',atmwat1,atmwat0,abs(atmwat1-atmwat0)
            print*, 'TB18 surf',surfwat1(i,j,nice_surf),surfwat0(id,jd,nice_surf),abs(surfwat1(i,j,nice_surf)-surfwat0(id,jd,nice_surf))
        ENDIF
    enddo
enddo
end subroutine checkconserv

subroutine checkco2(pres,atmosco2,dmass,surfco2,newsnow,flag)
!
!     Purpose
!     
!     check conservation all co2 in the column (mass)
!
!

!     Arguments
REAL, intent(in), dimension(:,:) :: pres
REAL, intent(in), dimension(:,:,:) :: atmosco2, dmass
REAL, intent(in), dimension(:,:) :: surfco2, newsnow
character (len=7), intent(in), optional :: flag
! LOCAL VARIABLES
INTEGER l,i,j,id,jd
REAL :: incert_co2=1.d-20
REAL, dimension(size(pres,1),size(pres,2)) :: masstot0,masstot1
!-----------------------------------------------------------------------


masstot0(:,:)=SUM(atmosco2(:,:,:),dim=3)+surfco2(:,:)*grav
masstot1(:,:)=masstot0+newsnow*grav-SUM(dmass(:,:,:),dim=3)*grav


do i=1,size(pres,1)
    do j=1,size(pres,2)
        IF (ABS(masstot0(i,j)-masstot1(i,j)).gt.incert_co2) THEN
            print*, ' conserv issue CO2 vs CO2 after ',flag,' :',masstot1(i,j),masstot0(i,j),abs(masstot0(i,j)-masstot1(i,j))
            print*,' newsnow, dmass, newsnow-dmass = ',newsnow(i,j)*grav,SUM(dmass(i,j,:))*grav,abs(newsnow(i,j)-sum(dmass(i,j,:)))*grav
        ENDIF
    end do
end do


end subroutine checkco2

end module testconserv_mod
