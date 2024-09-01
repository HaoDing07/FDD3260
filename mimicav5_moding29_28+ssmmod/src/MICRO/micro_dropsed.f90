#include "ctrparam.h"

! ================================================================
!
! MICRO_DROPSED.F90:

!	Subroutines processing droplet sedimentation
!
! ================================================================
!
module micro_dropsed
!
  USE gridno
  USE shared_data
  USE shared_nuclei
  USE shared_surf
  USE shared_hydro
  USE shared_diag
  USE shared_pressure
  USE shared_thermo
  USE sbpack
  USE kessler
  USE advection_lw
!	
  IMPLICIT NONE

  private

  public :: calc_sedimentation

  contains
!
!----------------------------------------------------------!
! Calculating terminal velocity in Stokes regime (Bergot,2016)!
!----------------------------------------------------------!
!
!
  function cal_vdrop (h,q,n,pxx)

  integer :: h
  real :: q, n, vdrop_mean, lambda_drop, lambda_c, pxx, cal_vdrop
  real :: n0

  cal_vdrop = 0.

    if (h == drop .AND. with_dropsed) then 

      n0 = max(1.e6,min(5.e7,n))

    !!!define alpha_drop and nu_drop in cm.nml and mimica.f90
      vdrop_mean = alpha_drop*(1+alpha_drop)/(lambda_drop**2) 

      cal_vdrop = max(vdrop_mean, 1.e-4)

      lambda_c = (pi/6*rhow*gamma(nu_drop+3/2)/gamma(nu_drop)*n0/1.29/q)**(1/3)
    
      lambda_drop = min(max(lambda_c, 1.e3),1.e4)

    end if

  return

  end function

!----------------------------------------------------------!
!     Initialise terminal velocity in the loop             !
!----------------------------------------------------------!
!
!  ==================================================
  subroutine calc_sedimentation ( dref, hydromtrl, qtm, hydromtrm, dens)
!  ==================================================
!
  integer :: i, j, k, kp, h
  real, dimension(ip_start:,jp_start:,:), intent(in) :: dref
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: dens
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: qtm
  type (hydrometeor), dimension(1:nhydro), intent(in) :: hydromtrl
  type (hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrm
  ! type (aero_3d), dimension(1:nmode) :: aero3dm
  real :: dt1
  real, dimension(1:nz) :: fup, fupn, fup2, fupn2, p, pn, pp, vp, vpn

  if (verbose > 1) write(7,*)'Strating calc_sedimentation'

!
  do h = 1,1
    do j = jt_start, jt_end
      do i = it_start, it_end
        if ( maxval(hydromtrl(h)%q) > qmin(h) .and. maxval(hydromtrl(h)%n) > xnmin(h) .and.with_dropsed ) then
          vp = 0.0 
          vpn = 0.0  
          do k = 1, nz
            if ( hydromtrl(h)%q(i,j,k) > qmin(h) .and. hydromtrl(h)%n(i,j,k) > xnmin(h) ) then
              vp(k) = -cal_vdrop (h, hydromtrl(h)%q(i,j,k), hydromtrl(h)%n(i,j,k), pxx(k,h))
              vpn(k) = -cal_vdrop (h, hydromtrl(h)%q(i,j,k), hydromtrl(h)%n(i,j,k), pxx(k,h))
            endif
          enddo
        endif
      enddo
    enddo
  enddo 

!
!----------------------------------------------------------!
!                    Store diagnostics                     !
!----------------------------------------------------------!
!  

dt1 = 1./dt0
do h = 1,1
  do j=jt_start,jt_end
    do i=it_start,it_end

      p = 0.
      pn = 0.

      fup = hydromtrl(h)%q(i,j,:)
call tridiag( nz, qmin(h), dref(i,j,:), -dt0*dz1w*vp, fup, fup2 ) 
      p = (fup - fup2)*dt1
      qtm(i,j,:) = qtm(i,j,:) - p
      hydromtrm(h)%q(i,j,:) = hydromtrm(h)%q(i,j,:) - p

      if ( moments == 2 ) then
        fupn = hydromtrl(h)%n(i,j,:)
call tridiag( nz, xnmin(h), dref(i,j,:), -dt0*dz1w*vpn, fupn, fupn2 ) 
        pn = (fupn - fupn2)*dt1
        hydromtrm(h)%n(i,j,:) = hydromtrm(h)%n(i,j,:) - pn
      endif

      ! #ifndef ANELASTIC
      ! rhom(i,j,:) = rhom(i,j,:) - dref(i,j,:)*p
      ! #endif
    
      if (out_micro) then
        if (out_diagq) diag(7)%qt(i,j,:) = diag(7)%qt(i,j,:) - dref(i,j,:)*cint*p
        if (h == drop .AND. out_diagl .AND. with_dropsed) diag(7)%qc(i,j,:) = diag(7)%qc(i,j,:) - dref(i,j,:)*cint*p
        if (h == drop .AND. with_dropsed) diag(10)%micro(drop)%q(i,j,:) = diag(10)%micro(drop)%q(i,j,:) - dref(i,j,:)*cint*p
        if ( moments == 2 ) then
          if (h == drop .AND. with_dropsed) diag(10)%micro(drop)%n(i,j,:) = diag(10)%micro(drop)%n(i,j,:) - dref(i,j,:)*cint*pn
        endif 

        if (out_diagl .AND. lmicro>0 .AND. with_dropsed) diag(7)%qc(i,j,:) = diag(7)%qc(i,j,:) - dens(i,j,:)*cint*pp
        if (out_micro .AND. lmicro>0 .AND. with_dropsed) diag(10)%micro(drop)%q(i,j,:) = diag(10)%micro(drop)%q(i,j,:) - dens(i,j,:)*cint*pp

      endif 

    enddo
  enddo
enddo


if (verbose > 1) write(7,*) 'Terminating calc_sedimentation'
end subroutine calc_sedimentation



! !  ==================================================
! subroutine calc_sedimentation_k ( dens, rhom, qtm, hydromtrm )
!   !  ==================================================
  
!   real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: dens
!   real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: rhom, qtm
!   type(hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrm
! !
!   integer :: i, j, k, h
!   real :: dt1, vp(1:nz), pp(1:nz), fup(1:nz), fup2(1:nz)
! !
!   if (verbose > 1) call write_debug('Strating calc_sedimentation_k')
! !
!   dt1 = 1./dt0
! !
! !----------------------------------------------------------!
! !                  Solve precipitation                     !
! !----------------------------------------------------------!
! !
!   do j=jt_start,jt_end
!     do i=it_start,it_end
!       ! if ( maxval(qr(i,j,:)) > qmin(drop) ) then
! !
!       do k = 1, nz
!         vp(k) = calc_sedimentation_k( dens(i,j,k), thermo%T(i,j,k) )
!       enddo
! !
! !  Solve tridiag system
! !
!       fup = qr(i,j,:) + cint*dt0*hydromtrm(drop)%q(i,j,:)
! !
!       call tridiag ( nz, qmin(drop), dens(i,j,:), -dt0*dz1w*vp, fup, fup2 ) 
!       !call bott ( nz, qmin(rain), dens(i,j,:), -dt0*dz1w*vp, dzw, fup, fup2 ) 
! !
!       pp = (fup - fup2)*dt1
! !	  
! !  Update tendencies
! !
!       hydromtrm(drop)%q(i,j,:) = hydromtrm(drop)%q(i,j,:) - pp
!       qtm(i,j,:)  = qtm(i,j,:) - pp
! !
! #ifndef ANELASTIC
!       rhom(i,j,:) = rhom(i,j,:) - dens(i,j,:)*pp
! #endif
! !
! !  diagnostics

!       if (out_diagl .AND. lmicro>0 .AND. with_dropsed) diag(7)%qc(i,j,:) = diag(7)%qc(i,j,:) - dens(i,j,:)*cint*pp
!       ! if (out_diagq ) diag(7)%qt(i,j,:) = diag(7)%qt(i,j,:) - dens(i,j,:)*cint*pp
!       if (out_diagr .AND. lmicro>0 .AND. with_dropsed) diag(6)%qr(i,j,:) = diag(6)%qr(i,j,:) - dens(i,j,:)*cint*pp
!       if (out_micro .AND. lmicro>0 .AND. with_dropsed) diag(10)%micro(drop)%q(i,j,:) = diag(10)%micro(drop)%q(i,j,:) - dens(i,j,:)*cint*pp
! !
!       endif
!     end do
!   end do    
!
!   if (verbose > 1) call write_debug('Terminating calc_sedimentation_k')
! !
!   end subroutine calc_sedimentation_k

end module micro_dropsed
