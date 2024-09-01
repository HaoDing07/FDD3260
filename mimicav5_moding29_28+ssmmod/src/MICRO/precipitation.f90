!
!------------------------------------------------------------------------!
! This file is part of MIMICA                                            !
!                                                                        !
! Copyright 2017-2021 Julien Savre                                       ! 
!                                                                        !
! This program is free software: you can redistribute it and/or modify   !
! it under the terms of the GNU General Public License as published by   !
! the Free Software Foundation, either version 3 of the License, or      !
! (at your option) any later version.                                    !
!                                                                        !
! This program is distributed in the hope that it will be useful,        !
! but WITHOUT ANY WARRANTY; without even the implied warranty of         !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          !
! GNU General Public License for more details.                           !
!                                                                        !
! You should have received a copy of the GNU General Public License      !
! along with this program.  If not, see <http://www.gnu.org/licenses/>.  !
!------------------------------------------------------------------------!
!

#include "ctrparam.h"

! ================================================================
!
!  PRECIPITATION.F                   
!
!  Purpose:
!	Subroutines for calculating precipitation of hydrometeors
!
!  Author
!	Julien Savre, MISU
!
! ================================================================
!
  module precipitation
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
  
  public :: calc_precip, calc_precip_k, cal_vdrop, calc_sedimentation

  contains
!
!  ==================================================
  subroutine calc_precip_k ( dens, qr, rhom, qtm, hydromtrm )
!  ==================================================

  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: dens, qr
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: rhom, qtm
  type(hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrm
!
  integer :: i, j, k, h
  real :: dt1, vp(1:nz), pp(1:nz), fup(1:nz), fup2(1:nz)
!
  if (verbose > 1) call write_debug('Strating calc_precip_k')
!
  dt1 = 1./dt0
!
!----------------------------------------------------------!
!                  Solve precipitation                     !
!----------------------------------------------------------!
!
  do j=jt_start,jt_end
    do i=it_start,it_end
      if ( maxval(qr(i,j,:)) > qmin(rain) ) then
!
      do k = 1, nz
        vp(k) = calc_vel_k( dens(i,j,k), thermo%T(i,j,k), qr(i,j,k) )
      enddo
!
!  Solve tridiag system
!
      fup = qr(i,j,:) + cint*dt0*hydromtrm(rain)%q(i,j,:)
!
      call tridiag ( nz, qmin(rain), dens(i,j,:), -dt0*dz1w*vp, fup, fup2 ) 
      !call bott ( nz, qmin(rain), dens(i,j,:), -dt0*dz1w*vp, dzw, fup, fup2 ) 
!
      pp = (fup - fup2)*dt1
!	  
!  Update tendencies
!
      hydromtrm(rain)%q(i,j,:) = hydromtrm(rain)%q(i,j,:) - pp
      qtm(i,j,:)  = qtm(i,j,:) - pp
!
#ifndef ANELASTIC
      rhom(i,j,:) = rhom(i,j,:) - dens(i,j,:)*pp
#endif
!
!  diagnostics
!
   !! implement droplet_sedimentation
      !if (out_diagl .AND. lmicro>0 .AND. with_dropsed) diag(7)%qc(i,j,:) = diag(7)%qc(i,j,:) - dens(i,j,:)*cint*pp
   !!
      if (out_diagq) diag(7)%qt(i,j,:) = diag(7)%qt(i,j,:) - dens(i,j,:)*cint*pp
      if (out_diagr) diag(6)%qr(i,j,:) = diag(6)%qr(i,j,:) - dens(i,j,:)*cint*pp
      if (out_micro) diag(10)%micro(rain)%q(i,j,:) = diag(10)%micro(rain)%q(i,j,:) - dens(i,j,:)*cint*pp
!
      endif
    end do
  end do    
!
  if (verbose > 1) call write_debug('Terminating calc_precip_k')
!
  end subroutine calc_precip_k
!
!  ==================================================
   subroutine calc_precip ( dref, hydromtrl, rhom, qtm, hydromtrm, aero3dm )
!  ==================================================
!
  integer  :: i, j, k, kp, h
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: dref
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: rhom, qtm
  type (hydrometeor), dimension(1:nhydro), intent(in) :: hydromtrl
  type (hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrm
  type (aero_3d), dimension(1:nmode) :: aero3dm
!  
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz,1:nhydro) :: precipn, veln
!
!  Internal
!
  real :: dt1, nmax, nmini, qmax, qmini
  real :: lambda, fl, Fold, Fold_n, Flux, Flux_n, ptot, pntot
  real, dimension(1:nz) :: fup, fupn, fup2, fupn2, vp, vpn, p, pn
!
  if (verbose > 1) call write_debug('Strating calc_precip')
!
  dt1 = 1./dt0
  precipn  = 0.0
!
!----------------------------------------------------------!
!                  Solving precipitation                   !
!----------------------------------------------------------!
!
!  Starting main hydrometeor loop
!
  !! implement droplet_sedimentation
  !do h = 1, nhydro
  !!
  do h = 2, nhydro
    do j=jt_start,jt_end
      do i=it_start,it_end
        if ( maxval(hydromtrl(h)%q(i,j,:)) > qmin(h) .and. maxval(hydromtrl(h)%n(i,j,:)) > xnmin(h) ) then
!
!  Calculate precip velocities
!
        vp = 0.
        vpn = 0.
        do k = 1, nz
          if ( hydromtrl(h)%q(i,j,k) > qmin(h) .and. hydromtrl(h)%n(i,j,k) > xnmin(h) ) then
            vp(k) = cal_vp (h, hydromtrl(h)%q(i,j,k), hydromtrl(h)%n(i,j,k), pxx(k,h))
            vpn(k) = cal_vpn (h, hydromtrl(h)%q(i,j,k), hydromtrl(h)%n(i,j,k), pxx(k,h))
          endif
        enddo
!
!  Solve precipitation implicitly: q
!
        fup = hydromtrl(h)%q(i,j,:) + cint*dt0*hydromtrm(h)%q(i,j,:)
!
	call tridiag ( nz, qmin(h), dref(i,j,:), -dt0*dz1w*vp, fup, fup2 ) 
        !call bott ( nz, qmin(h), dref(i,j,:), -dt0*dz1w*vp, dzw, fup, fup2 )  
!
	p = (fup - fup2)*dt1
!
	qtm(i,j,:) = qtm(i,j,:) - p
!
	hydromtrm(h)%q(i,j,:) = hydromtrm(h)%q(i,j,:) - p
! 
!  Solve precipitation implicitly: n
!
	if ( moments == 2 ) then
          fupn = hydromtrl(h)%n(i,j,:) + cint*dt0*hydromtrm(h)%n(i,j,:)
!
  	  call tridiag ( nz, xnmin(h), dref(i,j,:), -dt0*dz1w*vpn, fupn, fupn2 ) 
          !call bott ( nz, xnmin(h), dref(i,j,:), -dt0*dz1w*vpn, dzw, fupn, fupn2 )  
!
          pn = (fupn - fupn2)*dt1
!
  	  hydromtrm(h)%n(i,j,:) = hydromtrm(h)%n(i,j,:) - pn
!
	  precipn(i,j,:,h) = pn
	endif
!	  
!  Other precipitation tendencies
!
	if ( lmicro > 3 .and. h > ice )     &
	    hydromtrm(h)%w(i,j,:) = hydromtrm(h)%w(i,j,:) - p*hydromtrl(h)%w(i,j,:) / max(hydromtrl(h)%q(i,j,:),qmin(h))
!
#ifndef ANELASTIC
	rhom(i,j,:) = rhom(i,j,:) - dref(i,j,:)*p
#endif
!
!  micro diagnostics
!
	if ( out_diagq ) diag(7)%qt(i,j,:) = diag(7)%qt(i,j,:) - dref(i,j,:)*cint*p
        if ( h == rain .and. out_diagr ) diag(6)%qr(i,j,:) = diag(6)%qr(i,j,:) - dref(i,j,:)*cint*p
        if ( h == ice .and. out_diagi ) diag(6)%qi(i,j,:) = diag(6)%qi(i,j,:) - dref(i,j,:)*cint*p
        !! implement droplet_sedimentation
        !if ( h==drop .AND. out_diagl .AND. with_dropsed) diag(7)%qc(i,j,:) = diag(7)%qc(i,j,:) - dref(i,j,:)*cint*p
        !!
        !
        if ( out_micro ) then
        !! implement droplet_sedimentation
          !if (h == drop .AND. with_dropsed) diag(10)%micro(drop)%q(i,j,:) = diag(10)%micro(drop)%q(i,j,:) - dref(i,j,:)*cint*p
        !!
          if (h == rain) diag(10)%micro(rain)%q(i,j,:) = diag(10)%micro(rain)%q(i,j,:) - dref(i,j,:)*cint*p
          if (h == ice) diag(10)%micro(ice)%q(i,j,:) = diag(10)%micro(ice)%q(i,j,:) - dref(i,j,:)*cint*p
          if (h == grau) diag(10)%micro(grau)%q(i,j,:) = diag(10)%micro(grau)%q(i,j,:) - dref(i,j,:)*cint*p
          if (h == snow) diag(10)%micro(snow)%q(i,j,:) = diag(10)%micro(snow)%q(i,j,:) - dref(i,j,:)*cint*p
          if (h == hail) diag(10)%micro(hail)%q(i,j,:) = diag(10)%micro(hail)%q(i,j,:) - dref(i,j,:)*cint*p
	  if ( moments == 2 ) then
          !! implement droplet_sedimentation
            !if (h == drop .AND. with_dropsed) diag(10)%micro(drop)%n(i,j,:) = diag(10)%micro(drop)%n(i,j,:) - dref(i,j,:)*cint*pn
          !!       
            if (h == rain) diag(10)%micro(rain)%n(i,j,:) = diag(10)%micro(rain)%n(i,j,:) - dref(i,j,:)*cint*pn
            if (h == ice) diag(10)%micro(ice)%n(i,j,:) = diag(10)%micro(ice)%n(i,j,:) - dref(i,j,:)*cint*pn
            if (h == grau) diag(10)%micro(grau)%n(i,j,:) = diag(10)%micro(grau)%n(i,j,:) - dref(i,j,:)*cint*pn
            if (h == snow) diag(10)%micro(snow)%n(i,j,:) = diag(10)%micro(snow)%n(i,j,:) - dref(i,j,:)*cint*pn
            if (h == hail) diag(10)%micro(hail)%n(i,j,:) = diag(10)%micro(hail)%n(i,j,:) - dref(i,j,:)*cint*pn
	  endif
        endif
!
#ifdef NUC_CNT
        veln(i,j,:,h) = vpn(:)
#endif
!
!----------------------------------------------------------!
!                       Terminating                        !
!----------------------------------------------------------!
!
!  Terminate loops
!
        endif
      enddo
    enddo
  enddo
!
!  Scavenging of aerosols & chemicals   
!
  if (aero_flg%reg) then 
!
!  Liquid phase chemicals
!
#ifdef AQCHEM_ENABLE
        call aqprecip   ( pup, pdn )
#endif
!
!  Solid phase chemicals
!
#ifdef SOLIDCHEM_ENABLE
        call solprecip   ( pup, pdn )
#endif
!
!  Bulk CCN
!
#ifdef AERO_ENABLE
        call aero_precip ( aero3d, hydromtrl, precipn, aero3dm )
!#else 
!        call ccn_precip ( nuc, hydromtrl, precipn )
#endif
!
!  INs
!
#ifdef NUC_CNT
        call in_precip ( nucin2, hydromtrl, veln )
!
!       call in_scavenge ( nucin2, hydromtrl, pressure )
#endif
!
!  Impaction scavenging
!
#ifdef AERO_ENABLE
        if ( aero_flg%imp_scv ) call impaction_scav ( dref, aero3d, hydromtrl, aero3dm )
#endif
!
  endif
!
  if (verbose > 1) call write_debug('Terminating calc_precip')
!
!----------------------------------------------------------!
!      	
return
end
!
!  ==================================================
   subroutine ccn_precip ( nuc_tmp, hydromtr_tmp, precip )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Calculates ice nuclei sources/sinks due to precipitation
!
! ================================================================
!
  integer :: i, j, k, h, hl
  real :: ql
!
  type(nuclei) :: nuc_tmp  
  type(hydrometeor), dimension(1:nhydro)  :: hydromtr_tmp
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz,1:nhydro) :: precip
!
!----------------------------------------------------------!
!       Ice nuclei sources/sinks due to precipitation      !
!----------------------------------------------------------!
!
      hl = rain
      if ( lmicro > 1 ) hl = ice
      if ( lmicro > 2 ) hl = snow
!
      do j = jt_start, jt_end
	do i = it_start, it_end
	  do k = 1, nz
!
	    ql = 0.
	    do h = 1, hl
	      ql = ql + hydromtr_tmp(h)%n(i,j,k)
	    enddo
	    if ( ql > xnmin(rain) ) nuc_tmp%ccn(i,j,k) = nuc_tmp%ccn(i,j,k) - sum(precip(i,j,k,rain:hl)) * cint*dt0
!
	  enddo
	enddo	  
      enddo
!
  return
  end
!
!  ==================================================
   subroutine aero_precip ( aero_tmp, hydromtr_tmp, precip, aero3dm )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Calculates activated aerosols precipitation
!
! ================================================================
!
  integer :: i, j, k, h, hl, im
  real :: dt1, nl, pp, paero
!
  type (aero_3d), dimension(:)  :: aero_tmp
  type(hydrometeor), dimension(1:nhydro)  :: hydromtr_tmp
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz,1:nhydro) :: precip
  type (aero_3d), dimension(1:nmode) :: aero3dm
!
!----------------------------------------------------------!
!       Ice nuclei sources/sinks due to precipitation      !
!----------------------------------------------------------!
!
      dt1 = 1./dt0
!
      hl = rain
      if ( lmicro > 1 ) hl = ice
      if ( lmicro > 2 ) hl = snow
!
!  precipitation rates
!
      do j = jt_start, jt_end
	do i = it_start, it_end
          do k = 1, nz-1
            nl = 0.
            pp = 0.
	    do h = 1, hl
	      nl = nl + hydromtr_tmp(h)%n(i,j,k)
	    enddo	  
            if (nl > xnmin(rain)) pp = sum(precip(i,j,k,rain:hl))/nl
!
	    do im = 1, nmode
              paero = min( max(-pp*aero_tmp(im)%ma(i,j,k),-dt1*aero_tmp(im)%ma(i,j,k)), max(dt1*(aero_tmp(im)%ma(i,j,k+1)-aero_tmp(im)%ma(i,j,k)),0.) )
!
	      aero3dm(im)%ma(i,j,k) = aero3dm(im)%ma(i,j,k) + paero
	      if (ldiag .and. out_diaga) then
	        diag(8)%aero(im)%m(i,j,k) = diag(8)%aero(im)%m(i,j,k) + cint*paero
	        diag(9)%aero(im)%m(i,j,k) = diag(9)%aero(im)%m(i,j,k) + cint*paero
	      endif
	    enddo
!
	  enddo
	enddo	  
      enddo
!
  return
  end
!
#ifdef NUC_CNT
!
!  ==================================================
   subroutine in_precip ( nucin_tmp, hydromtr_tmp, vpn )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Calculates ice nuclei sources/sinks due to precipitation
!
! ================================================================
!
  integer :: i, j, k, m, h, hl
  real :: ql, qi, ul(1:nz), ui(1:nz) 
!
  type(hydrometeor), dimension(1:nhydro)  :: hydromtr_tmp
  type(nuclei_3d), dimension(3) :: nucin_tmp  
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz,1:nhydro) :: vpn
  real, dimension(1:nz,1:3)  :: nn, nnc, mm, mmc, nni, nnic, mmi, mmic
  real, dimension(1:nz,1:3)  :: pnucl, pnuclm, pnuclc, pnuclmc, pnuci, pnucim, pnucic, pnucimc
!
!----------------------------------------------------------!
!       Ice nuclei sources/sinks due to precipitation      !
!----------------------------------------------------------!
!
      hl = rain
      if ( lmicro > 1 ) hl = ice
      if ( lmicro > 2 ) hl = snow
!
      do j = jt_start, jt_end
	do i = it_start, it_end
	  ul = 0.0
	  ui = 0.0
	  pnucl = 0.
	  pnuci = 0.
	  pnuclm = 0.
	  pnucim = 0.
	  pnuclc = 0.
	  pnucic = 0.
	  pnuclmc = 0.
	  pnucimc = 0.
!
	  do h = 1, 3
	    nn(:,h) = nucin_tmp(h)%mode(2)%n(i,j,:)
	    nnc(:,h) = nucin_tmp(h)%mode(2)%nc(i,j,:)
	    nni(:,h) = nucin_tmp(h)%mode(3)%n(i,j,:)
	    nnic(:,h) = nucin_tmp(h)%mode(3)%nc(i,j,:)
!
#ifdef NUC_CNT1
	    mm(:,h) = nucin_tmp(h)%mode(2)%m(i,j,:)
	    mmc(:,h) = nucin_tmp(h)%mode(2)%mc(i,j,:)
	    mmi(:,h) = nucin_tmp(h)%mode(3)%m(i,j,:)
	    mmic(:,h) = nucin_tmp(h)%mode(3)%mc(i,j,:)
#endif
	  enddo
!
!  Weighted precipitation velocities
!
	  do k = 1, nz
	    ql = hydromtr_tmp(rain)%n(i,j,k) + hydromtr_tmp(drop)%n(i,j,k)
	    qi = sum( hydromtr_tmp(ice:hl)%n(i,j,k) )
!
  	    if ( ql > xnmin(rain) ) ul(k) = hydromtr_tmp(rain)%n(i,j,k) * vpn(i,j,k,rain) / ql
	    if ( qi > xnmin(ice) ) ui(k) = sum( hydromtr_tmp(ice:hl)%n(i,j,k) * vpn(i,j,k,ice:h1) ) / ql
	  enddo
!
!----------------------------------------------------------!
!          Calculate precipitation fluxes in drops         !
!----------------------------------------------------------!
!
	  do h = 1, 3
!
!  Solve precipitation implicitly: Fine mode
!
	    call tridiag ( nz, xin_min, den0, -dt0*dz1w*ul, nn(:,h), nn(:,h) ) 
            !call bott ( nz, xin_min, den0, -dt0*dz1w*ul, dzw, nn(:,h), nn(:,h) )  
!
            pnucl(:,h) = nucin_tmp(h)%mode(2)%n(i,j,:) - nn(:,h)
!
#ifdef NUC_CNT1
	    call tridiag ( nz, qin_min, den0, -dt0*dz1w*ul, mm(:,h), mm(:,h) ) 
            !call bott ( nz, qin_min, den0, -dt0*dz1w*ul, dzw, mm(:,h), mm(:,h) )  
!
            pnuclm(:,h) = nucin_tmp(h)%mode(2)%m(i,j,:) - mm(:,h)
#endif
!
!  Solve precipitation implicitly: Coarse mode
!
	    call tridiag ( nz, xin_min, den0, -dt0*dz1w*ul, nnc(:,h), nnc(:,h) ) 
            !call bott ( nz, xin_min, den0, -dt0*dz1w*ul, dzw, nnc(:,h), nnc(:,h) )  
!
            pnuclc(:,h) = nucin_tmp(h)%mode(2)%nc(i,j,:) - nnc(:,h)
!
#ifdef NUC_CNT1
	    call tridiag ( nz, qin_min, den0, -dt0*dz1w*ul, mmc(:,h), mmc(:,h) ) 
            !call bott ( nz, qin_min, den0, -dt0*dz1w*ul, dzw, mmc(:,h), mmc(:,h) )  
!
            pnuclmc(:,h) = nucin_tmp(h)%mode(2)%mc(i,j,:) - mmc(:,h)
#endif
!
	  enddo
!
!----------------------------------------------------------!
!          Calculate precipitation fluxes in ice           !
!----------------------------------------------------------!
!
	  do h = 1, 3
!
!  Solve precipitation implicitly: Fine mode
!
	    call tridiag ( nz, xin_min, den0, -dt0*dz1w*ui, nni(:,h), nni(:,h) ) 
            !call bott ( nz, xin_min, den0, -dt0*dz1w*ui, dzw, nni(:,h), nni(:,h) )  
!
            pnuci(:,h) = nucin_tmp(h)%mode(3)%n(i,j,:) - nni(:,h)
!
#ifdef NUC_CNT1
	    call tridiag ( nz, qin_min, den0, -dt0*dz1w*ui, mmi(:,h), mmi(:,h) ) 
            !call bott ( nz, qin_min, den0, -dt0*dz1w*ui, dzw, mmi(:,h), mmi(:,h) )  
!
            pnucim(:,h) = nucin_tmp(h)%mode(3)%m(i,j,:) - mmi(:,h)
#endif
!
!  Solve precipitation implicitly: Coarse mode
!
	    call tridiag ( nz, xin_min, den0, -dt0*dz1w*ui, nnic(:,h), nnic(:,h) ) 
            !call bott ( nz, xin_min, den0, -dt0*dz1w*ui, dzw, nnic(:,h), nnic(:,h) )  
!
            pnucic(:,h) = nucin_tmp(h)%mode(3)%nc(i,j,:) - nnic(:,h)
!
#ifdef NUC_CNT1
	    call tridiag ( nz, qin_min, den0, -dt0*dz1w*ui, mmic(:,h), mmic(:,h) ) 
            !call bott ( nz, qin_min, den0, -dt0*dz1w*ui, dzw, mmic(:,h), mmic(:,h) )  
!
            pnucimc(:,h) = nucin_tmp(h)%mode(3)%mc(i,j,:) - mmic(:,h)
#endif
!
	  enddo
!
!----------------------------------------------------------!
!                	 Update INs                        !
!----------------------------------------------------------!
!
          do k = 1, nz-1
	    do h = 1, 3
!
!  All nuclei
!
	      nucin_tmp(h)%mode(1)%n(i,j,:)  = nucin_tmp(h)%mode(1)%n(i,j,:)  - (pnucl(k,h) + pnuci(k,h))
	      nucin_tmp(h)%mode(1)%nc(i,j,:) = nucin_tmp(h)%mode(1)%nc(i,j,:) - (pnuclc(k,h) + pnucic(k,h))
#ifdef NUC_CNT1
	      nucin_tmp(h)%mode(1)%m(i,j,:)  = nucin_tmp(h)%mode(1)%m(i,j,:)  - (pnuclm(k,h) + pnucim(k,h))
	      nucin_tmp(h)%mode(1)%mc(i,j,:) = nucin_tmp(h)%mode(1)%mc(i,j,:) - (pnuclmc(k,h) + pnucimc(k,h))
#endif
!
!  Immersed nuclei
!
	      nucin_tmp(h)%mode(2)%n(i,j,:)  = nucin_tmp(h)%mode(2)%n(i,j,:)  - pnucl(k,h)
	      nucin_tmp(h)%mode(2)%nc(i,j,:) = nucin_tmp(h)%mode(2)%nc(i,j,:) - pnuclc(k,h)
#ifdef NUC_CNT1
	      nucin_tmp(h)%mode(2)%m(i,j,:)  = nucin_tmp(h)%mode(2)%m(i,j,:)  - pnuclm(k,h)
	      nucin_tmp(h)%mode(2)%mc(i,j,:) = nucin_tmp(h)%mode(2)%mc(i,j,:) - pnuclmc(k,h)
#endif
!
!  Frozen nuclei
!
	      nucin_tmp(h)%mode(3)%n(i,j,:)  = nucin_tmp(h)%mode(3)%n(i,j,:)  - pnuci(k,h)
	      nucin_tmp(h)%mode(3)%nc(i,j,:) = nucin_tmp(h)%mode(3)%nc(i,j,:) - pnucic(k,h)
#ifdef NUC_CNT1
	      nucin_tmp(h)%mode(3)%m(i,j,:)  = nucin_tmp(h)%mode(3)%m(i,j,:)  - pnucim(k,h)
	      nucin_tmp(h)%mode(3)%mc(i,j,:) = nucin_tmp(h)%mode(3)%mc(i,j,:) - pnucimc(k,h)
#endif
	    enddo
	  enddo
!
	enddo	  
      enddo
!
  return
  end
!
#endif
!
! ===============================================
  subroutine tridiag ( nzp, nlim, dd, mu, nn0, nn1 )  
! ===============================================

  integer :: k, nzp, ierr
  real :: nlim, nn0(1:nzp), nn1(1:nzp), nn2(1:nzp), mu(1:nzp), dd(1:nzp)
  real :: ak(1:nzp),dk(1:nzp),ck(1:nzp)
  real :: c = 1.
!
  nn2 = nn0
  nn1 = nn0
!
! solve tri-diagonal system a_k x_k-1 + d_k x_k + c_k x_k+1 = b_k.
!
  ak(1) = 0.
  ck(1) = dd(2)/dd(1) * 0.5*(mu(1)+mu(2))
  dk(1) = 1. - mu(1)
!
  do k=2,nzp-1
    ak(k) = 0.
    ck(k) = dd(k+1)/dd(k) * 0.5*(mu(k)+mu(k+1))
    dk(k) = 1. - 0.5*(mu(k-1)+mu(k))
  enddo
!
  ak(nzp) = 0.
  ck(nzp) = mu(nzp)
  dk(nzp) = 1. - 0.5*(mu(nzp-1)+mu(nzp))
!
!  Solve tridiag system
!
  call TRIDAG (ak,dk,ck,nn2,nn1,nzp,ierr)
!
!  Apply 2nd order filter to smooth out oscillations
!
!  nn1(1) = nn(1) + c/4.*(-nn(1) + nn(2))
!  do k=2,nzp-1
!    nn1(k) = nn(k) + c/4.*(nn(k-1) - 2.*nn(k) + nn(k+1))
!  enddo
!  nn1(nzp) = nn(nzp) + c/4.*(nn(nzp-1) - nn(nzp))
!
  end subroutine tridiag
!
  SUBROUTINE TRIDAG(A,B,C,R,U,N,CODE)
  
  !*****************************************************************
  ! Solves for a vector U of length N the tridiagonal linear set
  ! M U = R, where A, B and C are the three main diagonals of matrix
  ! M(N,N), the other terms are 0. R is the right side vector.
  !*****************************************************************

  INTEGER :: CODE,N,J
  REAL :: BET,GAM(N),A(N),B(N),C(N),R(N),U(N)
!
  IF(B(1) == 0.) THEN
    CODE=1
    RETURN
  END IF

  BET=B(1)
  U(1)=R(1)/BET
  DO J=2,N                    !Decomposition and forward substitution
    GAM(J)=C(J-1)/BET
    BET=B(J)-A(J)*GAM(J)
    IF(BET == 0.) THEN            !Algorithm fails
      CODE=2
      RETURN
    END IF
    U(J)=(R(J)-A(J)*U(J-1))/BET
  END DO

  DO J=N-1,1,-1                     !Back substitution
    U(J)=U(J)-GAM(J+1)*U(J+1)
  END DO
  
  CODE=0
  RETURN
  END
!
! ===============================================
  subroutine bott ( nzp, nlim, dd, mu, dzz, nn0, nn1 )  
! ===============================================

  USE advection_lw

  IMPLICIT NONE

  integer :: k, nzp 
  real :: nlim, nn0(1:nzp), nn1(1:nzp), nn2(1:nzp), nn(1:nzp), mu(1:nzp), dd(1:nzp), dzz(1:nzp)
  real :: dt1, vp, I, Im(1:nzp), Flux(1:nzp), vpn(1:nzp), dn0(1:nzp)
  real :: c=1.
!
  dt1 = 1./dt0
!
!  Corrected fluxes
!
  Im(1) = max( 0., dd(1)*nn0(1)*abs(mu(1)) )
  do k = 1, nzp-1
    vp = 0.5*(mu(k+1) + mu(k))
    Im(k+1) = max( 0., dd(k+1)*nn0(k+1)*abs(vp) - 0.5*(dd(k+1)*nn0(k+1) - dd(k)*nn0(k))*abs(vp)*(1. - abs(vp)) )
  enddo
!
!  Integration
!
  nn1 = nn0
  Flux(1) = -dzz(1)*dt1 * Im(1)/max(dd(1)*nn0(1),nlim) * dd(1)*nn0(1)
  do k = 1, nzp-1
    I = max( dd(k+1)*nn0(k+1), Im(k+1) + nlim )
    Flux(k+1) = -dzz(k+1)*dt1 * Im(k+1)/I * dd(k+1)*nn0(k+1)
    nn1(k) = nn0(k) - dt0*(Flux(k+1) - Flux(k))/dzz(k)/dd(k)
  enddo
!
!  Apply filter to smooth out oscillations
!
!  nn1(1) = nn(1) + c/16.*(-2.*nn(1) + 3.*nn(2) - nn(3))
!  nn1(2) = nn(2) + c/16.*(3.*nn(1) - 6.*nn(2) + 4.*nn(3) - nn(4))
!  do k=3,nzp-2
!    nn1(k) = nn(k) + c/16.*(-nn(k-2) + 4.*nn(k-1) - 6.*nn(k) + 4.*nn(k+1) - nn(k+2))
!  enddo
!  nn1(nzp-1) = nn(nzp-1) + c/16.*(3.*nn(nzp) - 6.*nn(nzp-1) + 4.*nn(nzp-2) - nn(nzp-3))
!  nn1(nzp) = nn(nzp) + c/16.*(-2.*nn(nzp) + 3.*nn(nzp-1) - nn(nzp-2))
!
  end subroutine bott
!

  !
!----------------------------------------------------------!
! Calculating terminal velocity in Stokes regime (Bergot,2016)!
!----------------------------------------------------------!
!
!
  function cal_vdrop (h,q,n,pxx)

    integer :: h
    real :: q, n, vdrop_mean, lambda_drop, lambda_c, pxx, cal_vdrop
    real :: n0, q0, eta_vis, d2_expect
  
    cal_vdrop = 0.
    eta_vis = 1.827e-5

      if (h == drop .AND. with_dropsed) then 
  

        ! n0 = max(1.e6,n)
        n0 = n
        ! n0 = 300.e6
        q0 = q
        ! q0 = max(q,1.e-12)
        if (q0/=0 .and. n0/=0) then
        ! lambda_c = (pi/6*1000*gamma(nu_drop+3/alpha_drop)/gamma(nu_drop)*n0/1.29/q)**(1/3)
          lambda_c = (4.0/3.0*pi*1000.0/1.29*n0/q0*gamma(nu_drop+3.0/alpha_drop)/gamma(nu_drop))**(1.0/3.0)
        
          ! lambda_drop = min(max(lambda_c, 1.0),10.0)
          lambda_drop = lambda_c

          ! d2_expect = lambda_drop**((alpha_drop-1)/alpha_drop)/alpha_drop*(alpha_drop*nu_drop+1)/alpha_drop*gamma((alpha_drop*nu_drop+1)/alpha_drop)/gamma(nu_drop)

          !!!define alpha_drop and nu_drop in cm.nml and mimica.f90
          ! vdrop_mean = 1.0/18.0*9.8*d2_expect*(1000-1.29)/eta_vis
          vdrop_mean = 2.0/9.0*9.8*(1000.0-1.29)/eta_vis/(lambda_drop**2.0)*gamma(nu_drop+5.0/alpha_drop)/gamma(nu_drop+3.0/alpha_drop)
    
          ! cal_vdrop = max(vdrop_mean, 1.e-4)
          cal_vdrop = vdrop_mean
          
        end if
        
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
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz,1:nhydro) :: precipn
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
              if ( hydromtrl(h)%q(i,j,k) > qmin(h) .and. hydromtrl(h)%n(i,j,k) > xnmin(h) .and.with_dropsed) then
                vp(k) = cal_vdrop (h, hydromtrl(h)%q(i,j,k), hydromtrl(h)%n(i,j,k), pxx(k,h))
                vpn(k) = cal_vdrop (h, hydromtrl(h)%q(i,j,k), hydromtrl(h)%n(i,j,k), pxx(k,h))
              endif
              ! if(vp(k)>0) write(7,*) 'The vp=', vp(k) 
              ! if(vp(k)>0) write(7,*) 'In cal_vdrop q=',hydromtrl(h)%q(i,j,k)
              ! if(vp(k)>0) write(7,*) 'In cal_vdrop n=',hydromtrl(h)%n(i,j,k)
            enddo
          endif
        enddo
      enddo
    enddo 
  
  !
  !----------------------------------------------------------!
  !   Solving sedimentation loop and Store diagnostics       !
  !----------------------------------------------------------!
  !  
  
  dt1 = 1./dt0
  precipn  = 0.0

  do h = 1,1
    do j=jt_start,jt_end
      do i=it_start,it_end
  
        p = 0.
        pn = 0.
  
        fup = hydromtrl(h)%q(i,j,:) + cint*dt0*hydromtrm(h)%q(i,j,:)
        call tridiag( nz, qmin(h), dref(i,j,:), -dt0*dz1w*vp, fup, fup2 ) 
        p = (fup - fup2)*dt1
        qtm(i,j,:) = qtm(i,j,:) - p
        hydromtrm(h)%q(i,j,:) = hydromtrm(h)%q(i,j,:) - p
  
        if ( moments == 2 ) then
          fupn = hydromtrl(h)%n(i,j,:) + cint*dt0*hydromtrm(h)%n(i,j,:)
        call tridiag( nz, xnmin(h), dref(i,j,:), -dt0*dz1w*vpn, fupn, fupn2 ) 
          pn = (fupn - fupn2)*dt1
          hydromtrm(h)%n(i,j,:) = hydromtrm(h)%n(i,j,:) - pn

        precipn(i,j,:,h) = pn
        endif
  
        ! #ifndef ANELASTIC
        ! rhom(i,j,:) = rhom(i,j,:) - dref(i,j,:)*p
        ! #endif

        if (out_diagq) diag(7)%qt(i,j,:) = diag(7)%qt(i,j,:) - dref(i,j,:)*cint*p
        if (h == drop .AND. out_diagl .AND. with_dropsed) diag(7)%qc(i,j,:) = diag(7)%qc(i,j,:) - dref(i,j,:)*cint*p
        if (out_micro) then
          if (h == drop .AND. with_dropsed) diag(10)%micro(drop)%q(i,j,:) = diag(10)%micro(drop)%q(i,j,:) - dref(i,j,:)*cint*p
        endif
        if ( moments == 2 ) then
          if (h == drop .AND. with_dropsed) diag(10)%micro(drop)%n(i,j,:) = diag(10)%micro(drop)%n(i,j,:) - dref(i,j,:)*cint*pn
        endif 
  
          if (out_diagl .AND. lmicro>0 .AND. with_dropsed) diag(7)%qc(i,j,:) = diag(7)%qc(i,j,:) - dens(i,j,:)*cint*pp
          if (out_micro .AND. lmicro>0 .AND. with_dropsed) diag(10)%micro(drop)%q(i,j,:) = diag(10)%micro(drop)%q(i,j,:) - dens(i,j,:)*cint*pp
  
  
      enddo
    enddo
  enddo
  
  
  if (verbose > 1) write(7,*) 'Terminating calc_sedimentation'
  end subroutine calc_sedimentation
  
end module precipitation
