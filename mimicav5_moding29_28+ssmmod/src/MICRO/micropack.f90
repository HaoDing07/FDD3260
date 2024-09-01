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
!  MICROPACK:
!	Package of subroutines common to all microphysics schemes
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================
!
  module micropack
!
  USE shared_data
  USE shared_diag
  USE shared_hydro
  USE shared_thermo
  USE shared_pressure
  USE shared_nuclei
  USE shared_rad
  USE shared_aerosol_new
  USE thermodynamics
  USE aerosols
!
  IMPLICIT NONE
!
  private
!
  public :: resolved_saturation, limit_hydro, init_micro
  public :: hydro_sizes, reflectivity, one_moment
!
  contains
!
!  ==================================================
   subroutine resolved_saturation ( k, cd, cdr, cdh, depi, depg, deps, deph, pp, Tem, ww, dTdt, cp, cv, exn,    &
   	sw, qsw, si, qsi, qv, qc, qr, qi, qg, qs, qh, dqc, dqr, dqi, dqg, dqs, dqh, des )
!  ==================================================
!
  integer :: k
  real :: l, dtc, tt, cd, cdr, cdh, depi, depg, deps, deph, dTdt
  real :: dqc, dqr, dqi, dqg, dqs, dqh, des, dqa, dql
  real :: qv, qc, qr, qi, qg, qs, qh, qsw, qsi
  real :: ww, cp, cv, sw, si, dqhl, dqhi, Tem, pp, exn
!
  real :: fbeta, ebeta, ebetai, gamma, gammai, gamma2, gammai2, alpha, beta, betai, mk_tmp(10)
!         
!----------------------------------------------------------!
!            Calculate resolved supersaturation            !
!----------------------------------------------------------!
!
	dqc = 0. 
	dqr = 0. 
	dqi = 0. 
	dqg = 0. 
	dqs = 0. 
	dqh = 0. 
	dqhl = 0. 
	dqhi = 0.
!
        alpha  = g / ((cp_a-cv_a)*Tem)
        beta   = cal_dqsdt(Tem,pp) !qsw*cal_flv(Tem) / ((cp_v-cv_v)*Tem*Tem)
        betai  = cal_dqsidt(Tem,pp) !qsi*cal_fls(Tem) / ((cp_v-cv_v)*Tem*Tem)
!
	gamma = 1. + beta*cal_flv(Tem)/cp
	gammai = 1. + betai*cal_fls(Tem)/cp
	gamma2 = 1. + beta*cal_fls(Tem)/cp
	gammai2 = 1. + betai*cal_flv(Tem)/cp
!
!  First: Semi-implicit method
!
        if ( dtcon <= 0. ) then
!
!  Warm microphysics only
!
	if ( lmicro == 1 ) then
!
	  mk_tmp(1) = gamma*(cd + cdr)
	  if (mk_tmp(1) > 1.e-12) then
	    mk_tmp(2) = 1. / mk_tmp(1)
	    mk_tmp(3) = -beta*(dTdt - g*ww/cp_a) !+ qsw*alpha*ww) 
!
	    ebeta = exp(-dt0/mk_tmp(2))
	    sw = mk_tmp(2)*mk_tmp(3) + mk_tmp(2)/dt0*(sw - mk_tmp(2)*mk_tmp(3))*(1. - ebeta)
!
            dqc = cd*sw
            dqr = cdr*sw
          endif  
!
!  Cold microphysics case
!
	else  if ( lmicro > 1 ) then
	  mk_tmp(1) = gamma*(cd + cdr) !+ cdh
	  mk_tmp(2) = gammai*(depi + depg + deps) !+ deph
	  mk_tmp(3) = mk_tmp(1) + gamma2*(depi + depg + deps)
	  mk_tmp(4) = mk_tmp(2) + gammai2*(cd + cdr)
!
	  if (mk_tmp(3) > 1.e-12 .and. mk_tmp(4) > 1.e-12) then
	    mk_tmp(5) = 1. / mk_tmp(3) 
	    mk_tmp(6) = 1. / mk_tmp(4) 
!
	    mk_tmp(7) = -beta*(dTdt - g*ww/cp_a) !+ qsw*alpha*ww) 
	    mk_tmp(8) = -betai*(dTdt - g*ww/cp_a) !+ qsi*alpha*ww) 
	    mk_tmp(9) = mk_tmp(7) - gamma2*(depi + depg + deps)*(qsw - qsi)
	    mk_tmp(10)= mk_tmp(8) - gammai2*(cd + cdr)*(qsi - qsw)
!
            ebeta = exp(-dt0/mk_tmp(5))
            ebetai = exp(-dt0/mk_tmp(6))
            sw = mk_tmp(9)*mk_tmp(5) + mk_tmp(5)/dt0*(sw - mk_tmp(9)*mk_tmp(5))*(1. - ebeta)    
            si = mk_tmp(10)*mk_tmp(6) + mk_tmp(6)/dt0*(si - mk_tmp(10)*mk_tmp(6))*(1. - ebetai)             
            !si = sw + qsw - qsi
!
            dqc = cd*sw
            dqr = cdr*sw
            dqi = depi*si
            dqg = depg*si
            dqs = deps*si
            !dqhl = cdh*sw/gamma
            !dqhi = deph*si/gammai
!
          else if (mk_tmp(1) > 1.e-12) then
	    mk_tmp(5) = 1. / mk_tmp(1)
	    mk_tmp(6) = -beta*(dTdt - g*ww/cp_a) !+ qsw*alpha*ww) 
!
	    ebeta = exp(-dt0/mk_tmp(5))
	    sw = mk_tmp(5)*mk_tmp(6) + mk_tmp(5)/dt0*(sw - mk_tmp(5)*mk_tmp(6))*(1. - ebeta)
!
            dqc = cd*sw
            dqr = cdr*sw
!
          else if (mk_tmp(2) > 1.e-12) then
	    mk_tmp(5) = 1. / mk_tmp(2)
	    mk_tmp(6) = -betai*(dTdt - g*ww/cp_a) !+ qsi*alpha*ww) 
!
	    ebetai = exp(-dt0/mk_tmp(5))
	    si = mk_tmp(5)*mk_tmp(6) + mk_tmp(5)/dt0*(si - mk_tmp(5)*mk_tmp(6))*(1. - ebetai)
!
            dqi = depi*si
            dqg = depg*si
            dqs = deps*si
	  endif
!
	endif
!
!  Else: Explicit sub-stepping
!
        else
!
        dtc = min(dt0,max(dtcon,1.e-8))
        dtc = dt0 / real(ceiling(dt0/dtc))
!
        if ( lmicro == 1 ) then
!
          mk_tmp(1) = gamma*(cd + cdr)
	  mk_tmp(2) = -beta*(dTdt - g*ww/cp_a) !+ qsw*alpha*ww) 
!
          tt = 0.
          do while (tt < dt0)
            sw = sw + (mk_tmp(2) - sw*mk_tmp(1))*dtc   
            dqc = dqc + cd*sw * dtc/dt0
            dqr = dqr + cdr*sw * dtc/dt0
            tt = tt + dtc
          enddo
!
        else if ( lmicro > 1 ) then
!
          mk_tmp(1) = gamma*(cd + cdr) !+ cdh
	  mk_tmp(2) = gammai*(depi + depg + deps) !+ deph
	  mk_tmp(3) = mk_tmp(1) + gamma2*(depi + depg + deps) 
	  mk_tmp(4) = mk_tmp(2) + gammai2*(cd + cdr)
	  mk_tmp(7) = -beta*(dTdt - g*ww/cp_a) !+ qsw*alpha*ww) 
	  mk_tmp(8) = -betai*(dTdt - g*ww/cp_a) !+ qsw*alpha*ww) 
!
          tt = 0.
          do while (tt < dt0)
	    mk_tmp(9) = mk_tmp(7) - gamma2*(depi + depg + deps)*(si - sw)
	    mk_tmp(10)= mk_tmp(8) - gammai2*(cd + cdr)*(sw - si)
!
            sw = sw + (mk_tmp(9) - sw*mk_tmp(3))*dtc
            si = si + (mk_tmp(10) - si*mk_tmp(4))*dtc
            dqc = dqc + cd*sw * dtc/dt0
            dqr = dqr + cdr*sw * dtc/dt0
            dqi = dqi + depi*si * dtc/dt0
            dqg = dqg + depg*si * dtc/dt0
            dqs = dqs + deps*si * dtc/dt0
            tt = tt + dtc
            !dqhl = dqhl + cdh*sw * dtc/dt0
            !dqhi = dqhi + deph*si * dtc/dt0
          enddo
!
        endif
!
	endif
!
!  Variations and limit
!
	if (qc > qc_min) dqc = max(dqc, -qc/dt0)
	if (qr > qr_min) dqr = max(dqr, -qr/dt0)
	if (qi > qi_min) dqi = max(dqi, -qi/dt0)
	if (qg > qg_min) dqg = max(dqg, -qg/dt0)
	if (qs > qs_min) dqs = max(dqs, -qs/dt0)
	if (qh > qh_min) dqhl = max(dqhl, -0.0066*qh/dt0)	! Hail evaporation is limited to the outter liquid layer of size dr=0.001Dh
	if (qh > qh_min) dqhi = max(dqhi, -qh/dt0)
	if (qh > qh_min) dqh = dqhl + dqhi
!
        dqa = dqc + dqr + dqi + dqg + dqs + dqh
	if ( qv - dqa*dt0 < 1.e-9 ) then
	  dql = (qv + 1.e-9)/(dqa*dt0)
	  dqc = dqc*dql
	  dqr = dqr*dql
	  dqi = dqi*dql
	  dqg = dqg*dql
	  dqs = dqs*dql
	  dqh = dqh*dql
	endif
!
!  Energy
!
#ifdef ISENTROPIC
	des = (cal_flv(Tem)*(dqc + dqr + dqhl) + cal_fls(Tem)*(dqi + dqg + dqs + dqhi)) / (cp*exn)
#endif
!
!----------------------------------------------------------!
!
return
end
!
! ==========================================
  subroutine limit_hydro ( qt, es, hydromtrl, aero3dl )
! ==========================================
!
!----------------------------------------------------------!
!                                                          !
!  Limit hydrometeors				           !
!                                                          !
!----------------------------------------------------------!
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: qt, es
  type(hydrometeor), dimension(:), intent(inout)  :: hydromtrl
  type(aero_3d), dimension(:), intent(inout) :: aero3dl
!
  integer  :: i, j, k, h, l, im, nb
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: qv, sumh
  real :: dq, qsum
  real, parameter :: small=1.e-14
!
  if (verbose > 1) call write_debug('Start limit_hydro')
!
!----------------------------------------------------------!
!           	    Limit hydrometeors            	   !
!----------------------------------------------------------!
!
#ifdef M3D
  nb = 4
#else
  nb = 2
#endif
!
!  Limit aerosols first
!
#ifdef AERO_ENABLE
  if (aero_flg%any) then
    do k=1,nz
      do j=jt_start,jt_end
        do i=it_start,it_end
          if ( sum(hydromtrl%n(i,j,k)) < xnmin(drop) .or. sum(hydromtrl%q(i,j,k)) < qmin(drop) ) then
            do h = 1, nmode
              aero3dl(h)%m(i,j,k) = aero3dl(h)%m(i,j,k) + aero3dl(h)%ma(i,j,k)
              aero3dl(h)%ma(i,j,k) = 0.
            enddo
          endif
        end do
      end do
    end do
  endif
#endif
!
!  Limit hydrometeors
!
  do h = 1, nhydro
    do k=1,nz
      if ( (minval(hydromtrl(h)%n(:,:,k)) < xnmin(h) .and. minval(hydromtrl(h)%n(:,:,k)) /= 0.) .or.     &
           (minval(hydromtrl(h)%q(:,:,k)) < qmin(h) .and. minval(hydromtrl(h)%q(:,:,k)) /= 0.) ) then
!
        do j=jt_start,jt_end
          do i=it_start,it_end
            if ( (hydromtrl(h)%n(i,j,k) < xnmin(h) .and. hydromtrl(h)%n(i,j,k) /= 0.) .or.        &
                 (hydromtrl(h)%q(i,j,k) < qmin(h) .and. hydromtrl(h)%q(i,j,k) /= 0.) ) then
                  ! write(7,*) 'TTTTEST', thermo_prop%cp(i,j,k), thermo%exn(i,j,k)
!
#ifdef ISENTROPIC
	      if ( h <= rain ) then
  	        es(i,j,k) = es(i,j,k) - cal_flv(thermo%T(i,j,k))*hydromtrl(h)%q(i,j,k) / (thermo_prop%cp(i,j,k)*thermo%exn(i,j,k))
	      else
	        es(i,j,k) = es(i,j,k) - cal_fls(thermo%T(i,j,k))*hydromtrl(h)%q(i,j,k) / (thermo_prop%cp(i,j,k)*thermo%exn(i,j,k))
	      endif
#endif
	      hydromtrl(h)%q(i,j,k) = 0.
!
	      if ( h /= drop .or. lndrop == 1 ) then
#ifdef AERO_ENABLE
	        if (aero_flg%reg) call regeneration ( i, j, k, pressure%dens(i,j,k), sum(hydromtrl%n(i,j,k)), aero3dl, abs(hydromtrl(h)%n(i,j,k))/dt0 )
#endif
                hydromtrl(h)%n(i,j,k) = 0.
	      endif
!
	      if ( lmicro > 3 .and. h > ice ) hydromtrl(h)%w(i,j,k) = 0.
            endif
          end do
        end do
!
      endif
    end do
  end do
!
!----------------------------------------------------------!
!                    Other limitations               	   !
!----------------------------------------------------------!
!
!  Limit if qv < 0
!
  call get_qv (qt, hydromtrl, qv, sumh=sumh)
!
  where ( qv < 0. ) qt = sumh + small
!
!  Limit if more rimed water than ice
!
  if ( lmicro > 3 .and. h >= grau ) then
    do h = grau, hail
      where ( hydromtrl(h)%q > qmin(h) .and. hydromtrl(h)%q < hydromtrl(h)%w .or. hydromtrl(h)%w < qmin(h) ) 
        hydromtrl(h)%w = hydromtrl(h)%q
      end where
    enddo
  endif
!
!----------------------------------------------------------!
!                     Limit aerosols     	           !
!----------------------------------------------------------!
!
  call limit_aerochem (	nuc				&
#ifdef AERO_ENABLE
			, aero3dl			&
#endif
#ifdef NUC_CNT
			, nucin2			&
#endif  
#ifdef CHEM_ENABLE
			, gas2				&
#endif
#ifdef AQCHEM_ENABLE	  
			, aqc2, aqr2			&
#endif
#ifdef SOLIDCHEM_ENABLE
			, solidi2			&
#endif
					)
! 
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate limit_hydro')
!
  end subroutine limit_hydro
!  
!  ==========================================
  subroutine limit_aerochem ( nucl		  &
#ifdef AERO_ENABLE
				  , aerol	  &
#endif
#ifdef NUC_CNT
				  , nucin	  &
#endif  
#ifdef CHEM_ENABLE
				  , gas 	  &
#endif
#ifdef AQCHEM_ENABLE	  
				  , aqc, aqr	  &
#endif
#ifdef SOLIDCHEM_ENABLE
				  , solidi	  &
#endif
							)
!  ==========================================
!
!----------------------------------------------------------!
!                                                          !
!  		Limit aerosols and chemicals		   !
!                                                          !
!----------------------------------------------------------!
!
USE shared_hydro
USE shared_aerosol_new
USE shared_nuclei
USE shared_gas
USE shared_aq
USE shared_solid
!
IMPLICIT NONE
!
!  Input/Output
!
  integer  :: i, j, k, h, m, im
!
  type(nuclei) :: nucl
!
#ifdef NUC_CNT
  type(nuclei_3d), dimension(3) :: nucin	
#endif
!
#ifdef AERO_ENABLE
  type(aero_3d), dimension(1:nmode) :: aerol	
#endif
!                                                          !
!----------------------------------------------------------!
!  
    do k=1,nz
      do j=jt_start,jt_end
	do i=it_start,it_end
!
!  Aerosols
!
#if (defined AERO_ENABLE)
	  do h = 1, nmode
	    if (aerol(h)%n(i,j,k) < naero_min .or. aerol(h)%m(i,j,k) < qaero_min) then
              aerol(h)%ma(i,j,k) = aerol(h)%ma(i,j,k) + max(aerol(h)%m(i,j,k),0.)
	      aerol(h)%n(i,j,k) = 0.0
	      aerol(h)%m(i,j,k) = 0.0
	    endif
	  enddo
#endif
!
!  Chemistry
!
#ifdef CHEM_ENABLE
	  if (gas(i,j,k)%o3    < 0.0) gas(i,j,k)%o3    = 0.0
	  if (gas(i,j,k)%co    < 0.0) gas(i,j,k)%co    = 0.0
	  if (gas(i,j,k)%h2o2  < 0.0) gas(i,j,k)%h2o2  = 0.0
	  if (gas(i,j,k)%xno   < 0.0) gas(i,j,k)%xno   = 0.0
	  if (gas(i,j,k)%xno2  < 0.0) gas(i,j,k)%xno2  = 0.0
	  if (gas(i,j,k)%xn2o5 < 0.0) gas(i,j,k)%xn2o5 = 0.0
	  if (gas(i,j,k)%hno3  < 0.0) gas(i,j,k)%hno3  = 0.0
	  if (gas(i,j,k)%ch2o  < 0.0) gas(i,j,k)%ch2o  = 0.0
	  if (gas(i,j,k)%ch3o2h< 0.0) gas(i,j,k)%ch3o2h= 0.0
	  if (gas(i,j,k)%so2   < 0.0) gas(i,j,k)%so2   = 0.0
	  if (gas(i,j,k)%h2so4 < 0.0) gas(i,j,k)%h2so4 = 0.0
	  if (gas(i,j,k)%dms   < 0.0) gas(i,j,k)%dms   = 0.0
	  if (gas(i,j,k)%nh3   < 0.0) gas(i,j,k)%nh3   = 0.0
#endif	
!
#ifdef AQCHEM_ENABLE
	  if (aqc(i,j,k)%o3	< 0.0) aqc(i,j,k)%o3	 = 0.0
	  if (aqc(i,j,k)%civ	< 0.0) aqc(i,j,k)%civ	 = 0.0
	  if (aqc(i,j,k)%h2o2	< 0.0) aqc(i,j,k)%h2o2   = 0.0
	  if (aqc(i,j,k)%xnv	< 0.0) aqc(i,j,k)%xnv	 = 0.0
	  if (aqc(i,j,k)%ch2o	< 0.0) aqc(i,j,k)%ch2o   = 0.0
	  if (aqc(i,j,k)%ch3o2h < 0.0) aqc(i,j,k)%ch3o2h = 0.0
	  if (aqc(i,j,k)%siv	< 0.0) aqc(i,j,k)%siv	 = 0.0
	  if (aqc(i,j,k)%svi	< 0.0) aqc(i,j,k)%svi	 = 0.0
	  if (aqc(i,j,k)%nh4	< 0.0) aqc(i,j,k)%nh4   = 0.0
	
	  if (aqr(i,j,k)%o3	< 0.0) aqr(i,j,k)%o3	 = 0.0
	  if (aqr(i,j,k)%civ	< 0.0) aqr(i,j,k)%civ	 = 0.0
	  if (aqr(i,j,k)%h2o2	< 0.0) aqr(i,j,k)%h2o2   = 0.0
	  if (aqr(i,j,k)%xnv	< 0.0) aqr(i,j,k)%xnv    = 0.0
	  if (aqr(i,j,k)%ch2o	< 0.0) aqr(i,j,k)%ch2o   = 0.0
	  if (aqr(i,j,k)%ch3o2h < 0.0) aqr(i,j,k)%ch3o2h = 0.0
	  if (aqr(i,j,k)%siv	< 0.0) aqr(i,j,k)%siv    = 0.0
	  if (aqr(i,j,k)%svi	< 0.0) aqr(i,j,k)%svi    = 0.0
	  if (aqr(i,j,k)%nh4	< 0.0) aqr(i,j,k)%nh4    = 0.0
#endif
!
#ifdef SOLIDCHEM_ENABLE
	  if (solidi(i,j,k)%o3	   < 0.0) solidi(i,j,k)%o3     = 0.0
	  if (solidi(i,j,k)%h2o2   < 0.0) solidi(i,j,k)%h2o2   = 0.0
	  if (solidi(i,j,k)%xnv	   < 0.0) solidi(i,j,k)%xnv    = 0.0
	  if (solidi(i,j,k)%ch2o   < 0.0) solidi(i,j,k)%ch2o   = 0.0
	  if (solidi(i,j,k)%ch3o2h < 0.0) solidi(i,j,k)%ch3o2h = 0.0
	  if (solidi(i,j,k)%siv	   < 0.0) solidi(i,j,k)%siv    = 0.0
	  if (solidi(i,j,k)%svi	   < 0.0) solidi(i,j,k)%svi    = 0.0
	  if (solidi(i,j,k)%nh4    < 0.0) solidi(i,j,k)%nh4    = 0.0
#endif
!
!----------------------------------------------------------!
!  			Limit INs			   !
!----------------------------------------------------------!
!
#ifdef NUC_CNT
	  do h = 1, 3
!
#ifdef NUC_CNT1
	    do m = 1, 3
	      if (nucin(h)%mode(m)%n(i,j,k) < xin_min .or. nucin(h)%mode(m)%m(i,j,k) < qin_min) then
	        nucin(h)%mode(m)%n(i,j,k) = 0.0
	        nucin(h)%mode(m)%m(i,j,k) = 0.0
	      endif
!
 	      if (nucin(h)%mode(m)%nc(i,j,k) < xin_min .or. nucin(h)%mode(m)%mc(i,j,k) < qin_min) then
	        nucin(h)%mode(m)%nc(i,j,k) = 0.0
 	        nucin(h)%mode(m)%mc(i,j,k) = 0.0
	      endif
	    enddo
!
	    if ( (nucin(h)%mode(2)%n(i,j,k) > xin_min .and. nucin(h)%mode(2)%m(i,j,k) > qin_min) .and. (hydromtr(drop)%q(i,j,k) < qmin(drop) .and. hydromtr(rain)%q(i,j,k) < qmin(rain)) ) then
	      nucin(h)%mode(2)%n(i,j,k) = 0.0
	      nucin(h)%mode(2)%m(i,j,k) = 0.0
	    endif
!	    
	    if ( (nucin(h)%mode(2)%nc(i,j,k) > xin_min .and. nucin(h)%mode(2)%mc(i,j,k) > qin_min) .and. (hydromtr(drop)%q(i,j,k) < qmin(drop) .and. hydromtr(rain)%q(i,j,k) < qmin(rain)) ) then
	      nucin(h)%mode(2)%nc(i,j,k) = 0.0
	      nucin(h)%mode(2)%mc(i,j,k) = 0.0
	    endif
!	    
	    if ( (nucin(h)%mode(3)%n(i,j,k) > xin_min .and. nucin(h)%mode(3)%m(i,j,k) > qin_min) .and. hydromtr(ice)%q(i,j,k) < qmin(ice) ) then
	      nucin(h)%mode(3)%n(i,j,k) = 0.0
	      nucin(h)%mode(3)%m(i,j,k) = 0.0
	    endif
!	    
	    if ( (nucin(h)%mode(3)%nc(i,j,k) > xin_min .and. nucin(h)%mode(3)%mc(i,j,k) > qin_min) .and. hydromtr(ice)%q(i,j,k) < qmin(ice) ) then
	      nucin(h)%mode(3)%nc(i,j,k) = 0.0
	      nucin(h)%mode(3)%mc(i,j,k) = 0.0
	    endif
#else
	    do m = 1, 3
	      if (nucin(h)%mode(m)%n(i,j,k) < xin_min) then
	        nucin(h)%mode(m)%n(i,j,k) = 0.0
	      endif
!
 	      if (nucin(h)%mode(m)%nc(i,j,k) < xin_min) then
	        nucin(h)%mode(m)%nc(i,j,k) = 0.0
	      endif
	    enddo
!
	    if ( nucin(h)%mode(2)%n(i,j,k) > xin_min .and. (hydromtr(drop)%q(i,j,k) < qmin(drop) .and. hydromtr(rain)%q(i,j,k) < qmin(rain)) ) then
	      nucin(h)%mode(2)%n(i,j,k) = 0.0
	    endif
!	    
	    if ( nucin(h)%mode(2)%nc(i,j,k) > xin_min .and. (hydromtr(drop)%q(i,j,k) < qmin(drop) .and. hydromtr(rain)%q(i,j,k) < qmin(rain)) ) then
	      nucin(h)%mode(2)%nc(i,j,k) = 0.0
	    endif
!	    
	    if ( nucin(h)%mode(3)%n(i,j,k) > xin_min .and. hydromtr(ice)%q(i,j,k) < qmin(ice) ) then
	      nucin(h)%mode(3)%n(i,j,k) = 0.0
	    endif
!	    
	    if ( nucin(h)%mode(3)%nc(i,j,k) > xin_min .and. hydromtr(ice)%q(i,j,k) < qmin(ice) ) then
	      nucin(h)%mode(3)%nc(i,j,k) = 0.0
	    endif
#endif
!
	  enddo	  
#endif
!
	end do
      end do
    end do
!                                                          !
!----------------------------------------------------------!
!
  return
  end
!
! ==========================================
  subroutine hydro_sizes ( hydromtrl )
! ==========================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculate mean hydrometeor sizes		           !
!                                                          !
!----------------------------------------------------------!
!
  type(hydrometeor), dimension(:), intent(inout)  :: hydromtrl
  integer  :: i, j, k, h
!  
!----------------------------------------------------------!
!
  if (out_dc .and. lmicro > 0) hydromtrl(drop)%d = 0.
  if (out_dr .and. lmicro > 0) hydromtrl(rain)%d = 0.
  if (out_di .and. lmicro > 1) hydromtrl(ice)%d = 0.
  if (out_dg .and. lmicro > 2) hydromtrl(grau)%d = 0.
  if (out_ds .and. lmicro > 2) hydromtrl(snow)%d = 0.
  if (out_dh .and. lmicro > 3) hydromtrl(hail)%d = 0.
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
	if ( hydromtrl(drop)%n(i,j,k) > xnmin(drop) .and. hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. out_dc ) then
          hydromtrl(drop)%d(i,j,k) = 1000.*cal_avd ( drop, hydromtrl(drop)%q(i,j,k), hydromtrl(drop)%n(i,j,k) )
	endif
!        
	if ( hydromtrl(rain)%n(i,j,k) > xnmin(rain) .and. hydromtrl(rain)%q(i,j,k) > qmin(rain) .and. out_dr ) then
          hydromtrl(rain)%d(i,j,k) = 1000.*cal_avd ( rain, hydromtrl(rain)%q(i,j,k), hydromtrl(rain)%n(i,j,k) )
	endif
!        
	if ( lmicro > 1 ) then
	  if ( hydromtrl(ice)%n(i,j,k) > xnmin(ice) .and. hydromtrl(ice)%q(i,j,k) > qmin(ice) .and. out_di ) then
            hydromtrl(ice)%d(i,j,k) = 1000.*cal_avd ( ice, hydromtrl(ice)%q(i,j,k), hydromtrl(ice)%n(i,j,k) )
          endif
	endif
!       
        if ( lmicro > 2 ) then
	  if ( hydromtrl(grau)%n(i,j,k) > xnmin(grau) .and. hydromtrl(grau)%q(i,j,k) > qmin(grau) .and. out_dg ) then
            hydromtrl(grau)%d(i,j,k) = 1000.*cal_avd ( grau, hydromtrl(grau)%q(i,j,k), hydromtrl(grau)%n(i,j,k) )
          endif
!        
	  if ( hydromtrl(snow)%n(i,j,k) > xnmin(snow) .and. hydromtrl(snow)%q(i,j,k) > qmin(snow) .and. out_ds ) then
            hydromtrl(snow)%d(i,j,k) = 1000.*cal_avd ( snow, hydromtrl(snow)%q(i,j,k), hydromtrl(snow)%n(i,j,k) )
          endif
        endif
!       
        if ( lmicro > 3 ) then
	  if ( hydromtrl(hail)%n(i,j,k) > xnmin(hail) .and. hydromtrl(hail)%q(i,j,k) > qmin(hail) .and. out_dh ) then
            hydromtrl(hail)%d(i,j,k) = 1000.*cal_avd ( hail, hydromtrl(hail)%q(i,j,k), hydromtrl(hail)%n(i,j,k) )
          endif
        endif
      enddo
    enddo
  enddo  
!  
!----------------------------------------------------------!
!
  return
  end subroutine
!
! ==========================================
  subroutine reflectivity ( hydromtrl, Z )
! ==========================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculate reflectivity			           !
!                                                          !
!----------------------------------------------------------!
!
  type(hydrometeor), dimension(:), intent(in)  :: hydromtrl
  real, dimension(ip_start:,jp_start:,:), intent(out) :: Z
!
  integer  :: i, j, k, h
  real :: d, czliq, aliq, cliq
  real, parameter :: Z0 = 1.e-18
!  
!----------------------------------------------------------!
!
  Z = 1.e-21
  cliq = hydrop(drop)%cm
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
	if ( hydromtrl(drop)%n(i,j,k) > xnmin(drop) .and. hydromtrl(drop)%q(i,j,k) > qmin(drop) ) then
          Z(i,j,k) = Z(i,j,k) + hydromtrl(drop)%n(i,j,k)*hydroc(drop)%czm*(cliq*cal_lambda(drop,hydromtrl(drop)%q(i,j,k),hydromtrl(drop)%n(i,j,k))**(-hydrop(drop)%am/hydrop(drop)%mu))**6.
	endif
!        
	if ( hydromtrl(rain)%n(i,j,k) > xnmin(rain) .and. hydromtrl(rain)%q(i,j,k) > qmin(rain) ) then
          Z(i,j,k) = Z(i,j,k) + hydromtrl(rain)%n(i,j,k)*hydroc(rain)%czm*(cliq*cal_lambda(rain,hydromtrl(rain)%q(i,j,k),hydromtrl(rain)%n(i,j,k))**(-hydrop(ice)%am/hydrop(rain)%mu))**6.
	endif
!        
	if ( lmicro > 1 ) then
	  if ( hydromtrl(ice)%n(i,j,k) > xnmin(ice) .and. hydromtrl(ice)%q(i,j,k) > qmin(ice) ) then
            Z(i,j,k) = Z(i,j,k) + hydromtrl(ice)%n(i,j,k)*hydroc(ice)%czm*(cliq*cal_lambda(ice,hydromtrl(ice)%q(i,j,k),hydromtrl(ice)%n(i,j,k))**(-hydrop(drop)%am/hydrop(ice)%mu))**6.
          endif
	endif
!       
        if ( lmicro > 2 ) then
	  if ( hydromtrl(grau)%n(i,j,k) > xnmin(grau) .and. hydromtrl(grau)%q(i,j,k) > qmin(grau) ) then
            Z(i,j,k) = Z(i,j,k) + hydromtrl(grau)%n(i,j,k)*hydroc(grau)%czm*(cliq*cal_lambda(grau,hydromtrl(grau)%q(i,j,k),hydromtrl(grau)%n(i,j,k))**(-hydrop(drop)%am/hydrop(grau)%mu))**6.
          endif
!        
	  if ( hydromtrl(snow)%n(i,j,k) > xnmin(snow) .and. hydromtrl(snow)%q(i,j,k) > qmin(snow) ) then
            Z(i,j,k) = Z(i,j,k) + hydromtrl(snow)%n(i,j,k)*hydroc(snow)%czm*(cliq*cal_lambda(snow,hydromtrl(snow)%q(i,j,k),hydromtrl(snow)%n(i,j,k))**(-hydrop(drop)%am/hydrop(snow)%mu))**6.
          endif
        endif
!       
        if ( lmicro > 3 ) then
	  if ( hydromtrl(hail)%n(i,j,k) > xnmin(hail) .and. hydromtrl(hail)%q(i,j,k) > qmin(hail) ) then
            Z(i,j,k) = Z(i,j,k) + hydromtrl(hail)%n(i,j,k)*hydroc(hail)%czm*(cliq*cal_lambda(hail,hydromtrl(hail)%q(i,j,k),hydromtrl(hail)%n(i,j,k))**(-hydrop(drop)%am/hydrop(hail)%mu))**6.
          endif
        endif
!
	Z(i,j,k) = 10.*log10( Z(i,j,k) / Z0 )	! Converting to dBz
      enddo
    enddo
  enddo  
!  
!----------------------------------------------------------!
!
  return
  end subroutine
!
! ==========================================
  subroutine one_moment ( hydromtrl )
! ==========================================
!
!----------------------------------------------------------!
!                                                          !
!  One-moment SB scheme: calculate number concentrations   !
!	diagnostically (except from drops and ice)	   !
!                                                          !
!----------------------------------------------------------!
!
  type(hydrometeor), dimension(:), intent(inout)  :: hydromtrl
  integer  :: i, j, k, h
!  
!----------------------------------------------------------!
!
  if (lmicro > 0) hydromtrl(rain)%n = 0.
  if (lmicro > 2) hydromtrl(grau)%n = 0.
  if (lmicro > 2) hydromtrl(snow)%n = 0.
  if (lmicro > 3) hydromtrl(hail)%n = 0.
!
  do k=1,nz
    do j=jp_start,jp_end
      do i=ip_start,ip_end
!        
	if ( hydromtrl(rain)%q(i,j,k) > qmin(rain) ) then
          hydromtrl(rain)%n(i,j,k) = cal_n1 ( rain, hydromtrl(rain)%q(i,j,k), hydromtrl(rain)%n(i,j,k) )
	endif
!       
        if ( lmicro > 2 ) then
	  if ( hydromtrl(grau)%q(i,j,k) > qmin(grau) ) then
            hydromtrl(grau)%n(i,j,k) = cal_n1 ( grau, hydromtrl(grau)%q(i,j,k), hydromtrl(grau)%n(i,j,k) )
          endif
!        
	  if ( hydromtrl(snow)%q(i,j,k) > qmin(snow) ) then
            hydromtrl(snow)%n(i,j,k) = cal_n1 ( snow, hydromtrl(snow)%q(i,j,k), hydromtrl(snow)%n(i,j,k) )
          endif
        endif
!       
        if ( lmicro > 3 ) then
	  if ( hydromtrl(hail)%q(i,j,k) > qmin(hail) ) then
            hydromtrl(hail)%n(i,j,k) = cal_n1 ( hail, hydromtrl(hail)%q(i,j,k), hydromtrl(hail)%n(i,j,k) )
          endif
        endif
!
      enddo
    enddo
  enddo  
!  
!----------------------------------------------------------!
!
  return
  end subroutine
!
! ==========================================
  subroutine init_micro
! ==========================================
	
  integer :: h, l, i, j, k
!
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Start init_micro')
!
!  Allocations
!
  allocate( pxx(1:nz,1:nhydro), qmin(1:nhydro), lmax(1:nhydro), lmin(1:nhydro), xnmin(1:nhydro) ) 
!
  if (out_dc .and. lmicro > 0) allocate ( hydromtr(drop)%d(ip_start:ip_end,jp_start:jp_end,1:nz) )
  if (out_dr .and. lmicro > 0) allocate ( hydromtr(rain)%d(ip_start:ip_end,jp_start:jp_end,1:nz) )
  if (out_di .and. lmicro > 1) allocate ( hydromtr(ice)%d(ip_start:ip_end,jp_start:jp_end,1:nz) )
  if (out_dg .and. lmicro > 2) allocate ( hydromtr(grau)%d(ip_start:ip_end,jp_start:jp_end,1:nz) )
  if (out_ds .and. lmicro > 2) allocate ( hydromtr(snow)%d(ip_start:ip_end,jp_start:jp_end,1:nz) )
  if (out_dh .and. lmicro > 3) allocate ( hydromtr(hail)%d(ip_start:ip_end,jp_start:jp_end,1:nz) )
!
!  Limit concentrations
!
  if (lmicro > 0) qmin(drop) = qc_min
  if (lmicro > 0) qmin(rain) = qr_min
  if (lmicro > 1) qmin(ice) = qi_min
  if (lmicro > 2) qmin(grau) = qg_min
  if (lmicro > 2) qmin(snow) = qs_min
  if (lmicro > 3) qmin(hail) = qh_min
  if (lmicro > 0) lmin(drop) = lc_min
  if (lmicro > 0) lmin(rain) = lr_min
  if (lmicro > 1) lmin(ice) = li_min
  if (lmicro > 2) lmin(grau) = lg_min
  if (lmicro > 2) lmin(snow) = ls_min
  if (lmicro > 3) lmin(hail) = lh_min
  if (lmicro > 0) lmax(drop) = lc_max
  if (lmicro > 0) lmax(rain) = lr_max
  if (lmicro > 1) lmax(ice) = li_max
  if (lmicro > 2) lmax(grau) = lg_max
  if (lmicro > 2) lmax(snow) = ls_max
  if (lmicro > 3) lmax(hail) = lh_max
  if (lmicro > 0) xnmin(drop) = xnc_min
  if (lmicro > 0) xnmin(rain) = xnr_min
  if (lmicro > 1) xnmin(ice) = xni_min
  if (lmicro > 2) xnmin(grau) = xng_min 
  if (lmicro > 2) xnmin(snow) = xns_min 
  if (lmicro > 3) xnmin(hail) = xnh_min 
!               
!  Set ice habit
!
  if (ice_habit == 'PLA') then
    ihab = 1
  else if (ice_habit == 'DEN') then
    ihab = 2
  else if (ice_habit == 'BEN') then
    ihab = 3
  else if (ice_habit == 'COL') then
    ihab = 4
  else if (ice_habit == 'HCO') then
    ihab = 5
  else if (ice_habit == 'BUL') then
    ihab = 6
  else if (ice_habit == 'ASS') then
    ihab = 7
  else if (ice_habit == 'PAS') then
    ihab = 8
  else if (ice_habit == 'ISD') then
    ihab = 9
  else
    ihab = 10
  endif
!
!  Initialize hydro properties
!
  allocate ( hydrop(1:nhydro), hydroc(1:nhydro) )
!
  if (lmicro == 1) then
    hydrop = (/ hydro_param(0.124, 3.75e5, 0.333, 0.667, 1., 1., 1., 0., 0., 0.5, 0., 4.2e-15, 2.e-10),	                &	! cloud
	        hydro_param(0.124, 159., 0.333, 0.266, -2./3., 1./3., 0.78, 0.308, 0., 0.5, 5.e5, 6.54e-11, 5.e-6) /)  		! rain 
!
  else if (lmicro == 2) then
    hydrop = (/ hydro_param(0.124, 3.75e5, 0.333, 0.667, 1., 1., 1., 0., 0., 0.5, 0., 4.2e-15, 2.e-10),	                &	! cloud
	        hydro_param(0.124, 159., 0.333, 0.266, -2./3., 1./3., 0.78, 0.308, 0., 0.5, 5.e5, 6.54e-11, 5.e-6),	& 	! rain
	        hydro_param(0.217, 317., 0.302, 0.363, 1., 1./3., 0.86, 0.28, 0., 0.318, 0., 4.42e-14, 5.e-7) /)  		! ice
!
  else if (lmicro == 3) then
    hydrop = (/ hydro_param(0.124, 3.75e5, 0.333, 0.667, 1., 1., 1., 0., 0., 0.5, 0., 4.2e-15, 2.e-10),	        &	! cloud
	        hydro_param(0.124, 159., 0.333, 0.266, -2./3., 1./3., 0.78, 0.308, 0., 0.5, 5.e5, 6.54e-11, 5.e-6),	& 	! rain
	        hydro_param(0.217, 317., 0.302, 0.363, 1., 1./3., 0.86, 0.28, 0., 0.318, 0., 4.42e-14, 5.e-7),   	&	! ice
	        hydro_param(0.190, 40., 0.323, 0.230, -2./3., 1./3., 0.78, 0.308, 0., 0.5, 4.5e6, 6.54e-11, 1.e-4),      &  	! graupel
	        hydro_param(8.156, 27.7, 0.526, 0.216, 1./3., 1./2., 0.78, 0.308, 0., 0.318, 1.8e8, 6.54e-11, 5.e-6) /)	        ! snow
!
  else if (lmicro == 4) then
    hydrop = (/ hydro_param(0.124, 3.75e5, 0.333, 0.667, 1., 1., 1., 0., 0., 0.5, 0., 4.2e-15, 2.e-10),	                &	! cloud
	        hydro_param(0.124, 159., 0.333, 0.266, -2./3., 1./3., 0.78, 0.308, 0., 0.5, 5.e5, 6.54e-11, 5.e-6),	& 	! rain
	        hydro_param(0.217, 317., 0.302, 0.363, 1., 1./3., 0.86, 0.28, 0., 0.318, 0., 4.42e-14, 5.e-7),   	&	! ice
	        hydro_param(0.190, 40., 0.323, 0.230, -2./3., 1./3., 0.78, 0.308, 0., 0.5, 4.5e6, 6.54e-11, 1.e-4),      &  	! graupel
	        hydro_param(8.156, 27.7, 0.526, 0.216, 1./3., 1./2., 0.78, 0.308, 0., 0.318, 1.8e8, 6.54e-11, 5.e-6),	&       ! snow
	        hydro_param(0.131, 56.97, 0.333, 0.197, -2./3., 1./3., 0.78, 0.308, 0., 0.5, 2.e3, 1.73e-9, 1.e-2) /)	        ! hail
  endif
!
  if (lmicro > 1 .and. ihab <= 9) hydrop(ice) = icep(ihab)
!
!  Initialize microphysical constants
!
  do h = 1, nhydro
!
!  One moment parameter
!
    hydroc(h)%clam1 = 0.
    if (h/=drop .and. h/=ice) hydroc(h)%clam1 = hydrop(h)%mu / gamma((hydrop(h)%nu+2.)/hydrop(h)%mu) / hydrop(h)%a0
!
!  Lambda parameter
!
    hydroc(h)%clam = gamma((hydrop(h)%nu+1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu+2.)/hydrop(h)%mu)
!
!  Parameters for mean diameters and volumes
!
    hydroc(h)%cmm = gamma((hydrop(h)%nu + 2.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu)
    hydroc(h)%cdm = gamma((hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu)
    hydroc(h)%cvm = gamma((3.*hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) 
    hydroc(h)%cem = gamma((3.*hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((2.*hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) 
    hydroc(h)%czm = gamma((6.*hydrop(drop)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu)
    hydroc(h)%cr6 = gamma((hydrop(h)%nu + 3.)/hydrop(h)%mu)*gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 2.)/hydrop(h)%mu)**2.
!
!  Parameters for terminal fall speeds
!
    hydroc(h)%cvpq = gamma((hydrop(h)%nu + hydrop(h)%btv + 2.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu)
    hydroc(h)%cvpn = gamma((hydrop(h)%nu + hydrop(h)%btv + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu)
!
!  Parameters of the ventilation factor
!
    hydroc(h)%aventn = hydrop(h)%av*gamma((hydrop(h)%nu + hydrop(h)%am)/hydrop(h)%mu)/gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(hydrop(h)%am-1.)
    hydroc(h)%bventn = hydrop(h)%bv*gamma((hydrop(h)%nu + 1.5*hydrop(h)%am + 0.5*hydrop(h)%btv)/hydrop(h)%mu)/gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(1.5*hydrop(h)%am+0.5*hydrop(h)%btv-1.)
    hydroc(h)%cventn = hydrop(h)%cv*gamma((hydrop(h)%nu + 2.*hydrop(h)%am + hydrop(h)%btv)/hydrop(h)%mu)/gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(2.*hydrop(h)%am+hydrop(h)%btv-1.)
!
    hydroc(h)%aventq = hydrop(h)%av*gamma((hydrop(h)%nu + hydrop(h)%am + 1.)/hydrop(h)%mu)/gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(hydrop(h)%am)
    hydroc(h)%bventq = hydrop(h)%bv*gamma((hydrop(h)%nu + 1.5*hydrop(h)%am + 0.5*hydrop(h)%btv + 1.)/hydrop(h)%mu)/gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(1.5*hydrop(h)%am+0.5*hydrop(h)%btv)
    hydroc(h)%cventq = hydrop(h)%cv*gamma((hydrop(h)%nu + 2.*hydrop(h)%am + hydrop(h)%btv + 1.)/hydrop(h)%mu)/gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(2.*hydrop(h)%am+hydrop(h)%btv)
!
!  Parameters for collision integrals
!
    hydroc(h)%deltaq = gamma((2.*hydrop(h)%am + hydrop(h)%nu + 2.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(2.*hydrop(h)%am+1.)
    hydroc(h)%deltan = gamma((2.*hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(2.*hydrop(h)%am)
!
    hydroc(h)%thetaq = gamma((2.*hydrop(h)%btv + 2.*hydrop(h)%am + hydrop(h)%nu + 2.)/hydrop(h)%mu) / gamma((2.*hydrop(h)%am + hydrop(h)%nu + 2.)/hydrop(h)%mu) * hydroc(h)%clam**(2.*hydrop(h)%btv)
    hydroc(h)%thetan = gamma((2.*hydrop(h)%btv + 2.*hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((2.*hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(h)%clam**(2.*hydrop(h)%btv)
!
    do l = 1, nhydro
      hydroc(h)%deltaijq(l) = gamma((hydrop(l)%am + hydrop(l)%nu + 2.)/hydrop(l)%mu) / gamma((hydrop(l)%nu + 1.)/hydrop(l)%mu) * gamma((hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(l)%clam**(hydrop(l)%am+1.) * hydroc(h)%clam**(hydrop(h)%am)
      hydroc(h)%deltaijn(l) = gamma((hydrop(l)%am + hydrop(l)%nu + 1.)/hydrop(l)%mu) / gamma((hydrop(l)%nu + 1.)/hydrop(l)%mu) * gamma((hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(l)%clam**(hydrop(l)%am) * hydroc(h)%clam**(hydrop(h)%am)
!
      hydroc(h)%thetaijq(l) = gamma((hydrop(l)%btv + hydrop(l)%am + hydrop(l)%nu + 2.)/hydrop(l)%mu) / gamma((hydrop(l)%am + hydrop(l)%nu + 2.)/hydrop(l)%mu) * gamma((hydrop(h)%btv + hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(l)%clam**(hydrop(l)%btv) * hydroc(h)%clam**(hydrop(h)%btv)
      hydroc(h)%thetaijn(l) = gamma((hydrop(l)%btv + hydrop(l)%am + hydrop(l)%nu + 1.)/hydrop(l)%mu) / gamma((hydrop(l)%am + hydrop(l)%nu + 1.)/hydrop(l)%mu) * gamma((hydrop(h)%btv + hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%am + hydrop(h)%nu + 1.)/hydrop(h)%mu) * hydroc(l)%clam**(hydrop(l)%btv) * hydroc(h)%clam**(hydrop(h)%btv)
    enddo
!
!  Parameters for impaction scavenging
!
    hydroc(h)%cimp = gamma((2.*hydrop(h)%am + hydrop(h)%btv + hydrop(h)%nu + 1.)/hydrop(h)%mu) / gamma((hydrop(h)%nu + 1.)/hydrop(h)%mu)
  enddo
!
  if (verbose > 1) call write_debug('Terminate init_micro')
!  
!----------------------------------------------------------!
!
  return
  end subroutine
!  
  end module micropack
