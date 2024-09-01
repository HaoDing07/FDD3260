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
!  kessler.F90:
!
!  Purpose:
!	Original Kessler microphysics
!
!  Author:
!	Julien Savre, MIM, Ludwig-Maximilian Universitat, Munchen
!
! ================================================================

  module kessler

!  
!  ---------------------------------------------
!
  USE gridno
  USE shared_data
  USE shared_pressure
  USE shared_hydro
  USE shared_thermo
  USE shared_surf
  USE shared_diag
  USE micropack

  IMPLICIT NONE
  
  ! General parameters
  real, parameter :: ar = 523.6
  real, parameter :: br = 3.
  real, parameter :: as = 2.5e-2
  real, parameter :: bs = 2.
  real, parameter :: kmr = 7.4866
  real, parameter :: kms = 0.3684
  real, parameter :: bmr = 0.75
  real, parameter :: bms = 0.6667
  real, parameter :: kdr = 0.2427
  real, parameter :: kds = 3.83877
 
  ! Autoconversion parameters
  real, parameter :: k1r = 1.e3
  real, parameter :: k1s = 1.e3
  ! qauto defined in cm.nml
  
  ! Accretion parameters
  real, parameter :: k2r  = 0.705
  real, parameter :: k2s  = 253.73
  real, parameter :: effr = 0.6283
  real, parameter :: effs = 0.0471
  real, parameter :: nr0 = 1.e7
  real, parameter :: ns0 = 1.e7
  
  ! Fall-speed parameters
  real, parameter :: k3r = 92.107
  real, parameter :: k3s = 6.5444
  
  ! Phase-changes parameters
  real, parameter :: k4r = 0.7795
  real, parameter :: k4s = 26.492
  real, parameter :: k5r = 532.975  
  real, parameter :: k5s = 1416.311
  real, parameter :: avr = 0.78
  real, parameter :: avs = 0.65
  real, parameter :: bvr = 0.27
  real, parameter :: bvs = 0.39
  
  private
  
  public :: micro_k, calc_avd_k, calc_n_k, calc_vel_k
  
  contains
!
!  ==========================================
  subroutine micro_k ( dref, pressurel, statel, hydromtrl, esm, hydromtrm ) 
!  ==========================================

    real, dimension(ip_start:,jp_start:,:), intent(in) :: dref
    type(atm_state), intent(in) :: statel
    type(hydrometeor), dimension(1:nhydro), intent(in) :: hydromtrl
    type(atm_pressure), intent(in) :: pressurel
!    
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: esm
    type(hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrm
!
    integer :: i, j, k, h
    real :: dt1, qctmp, nctmp, qrtmp, estmp, ttmp, dtmp, dql, sumh, fice
    real :: cnc, auq, acc, eva, qv, qsw, qsi, exn, pp
!
    if (verbose > 1) call write_debug('Strating micro_k')
!   
!----------------------------------------------------------!
!                     Kessler model                        !
!----------------------------------------------------------!
!
    cnc = 0.
    if (lndrop == 1) cnc = 1.
    dt1 = 1./dt0
!
    if (lmicro > 0) then
!
    do k=1,nz
      do j=jt_start,jt_end
        do i=it_start,it_end
!
	  auq = 0.
	  acc = 0.
	  eva = 0.
!
#ifdef ANELASTIC
	  pp = p0(k)
#else
          pp  = pressurel%p(i,j,k) + p0(k)
#endif
!
          sumh = 0.
          do h = 1, nhydro
	    sumh = sumh + hydromtrl(h)%q(i,j,k)
	  enddo
	  qctmp = hydromtrl(drop)%q(i,j,k)
	  nctmp = hydromtrl(drop)%n(i,j,k)
	  qrtmp = hydromtrl(rain)%q(i,j,k)
	  estmp = statel%es(i,j,k)
	  ttmp  = thermo%T(i,j,k)
	  dtmp = pressurel%dens(i,j,k)
	  exn = thermo%exn(i,j,k)
 
    	  fice = 1.
          if (lfreeze > 0) fice = max(min( 1., (ttmp - tii) / (t00 - tii) ), 0.)
!
!  Autoconversion from cloud to rain: c + c --> r
!
          call aucr_k ( fice, dtmp, ttmp, hydromtrl(drop)%q(i,j,k), hydromtrl(rain)%q(i,j,k), auq )
!
          if ( out_micro ) then
	    diag(2)%micro(drop)%q(i,j,k) = diag(2)%micro(drop)%q(i,j,k) - cint*dref(i,j,k)*min(auq*dt0,qctmp)
	    diag(2)%micro(rain)%q(i,j,k) = diag(2)%micro(rain)%q(i,j,k) + cint*dref(i,j,k)*min(auq*dt0,qctmp)
          endif
!
          qrtmp = qrtmp + min(auq*dt0,qctmp)
          qctmp = qctmp - min(auq*dt0,qctmp)
          if (hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. lndrop > 0) nctmp = nctmp - cnc*min(auq*hydromtrl(drop)%n(i,j,k)/hydromtrl(drop)%q(i,j,k)*dt0,nctmp)
!  
!  Collection of cloud drop by rain: c + r --> r
!
	  if ( hydromtrl(rain)%q(i,j,k) > qmin(rain) ) then
            call ccr_k ( fice, dtmp, hydromtrl(drop)%q(i,j,k), hydromtrl(rain)%q(i,j,k), acc )
!
            if ( out_micro ) then
	      diag(3)%micro(drop)%q(i,j,k) = diag(3)%micro(drop)%q(i,j,k) - cint*dref(i,j,k)*min(acc*dt0,qctmp)
	      diag(3)%micro(rain)%q(i,j,k) = diag(3)%micro(rain)%q(i,j,k) + cint*dref(i,j,k)*min(acc*dt0,qctmp)
            endif
!
            qrtmp = qrtmp + min(acc*dt0,qctmp)
            qctmp = qctmp - min(acc*dt0,qctmp)
            if (hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. lndrop > 0) nctmp = nctmp - cnc*min(acc*hydromtrl(drop)%n(i,j,k)/hydromtrl(drop)%q(i,j,k)*dt0,nctmp)
!  
!  Evaporation/condensation of precipitable water: r <--> v
!	  
	    qv  = statel%qt(i,j,k) - sumh
  	    qsw = cal_qsw( ttmp, pp )
  	    qsi = cal_qsi( ttmp, pp )
!
	    call evap_k ( fice, dtmp, ttmp, qv/qsw, qv/qsi, hydromtrl(rain)%q(i,j,k), eva )
!
            eva = min(eva,hydromtrl(rain)%q(i,j,k)*dt1)
            if ( out_micro ) then
	      diag(1)%micro(rain)%q(i,j,k) = diag(1)%micro(rain)%q(i,j,k) + cint*dref(i,j,k)*max(eva*dt0,-qrtmp)
            endif
!
	    eva = max(eva*dt0,-qrtmp)
#ifdef ISENTROPIC
            estmp = estmp + dref(i,j,k) * cal_flv(ttmp) * eva / (thermo_prop%cp(i,j,k)*exn)
#endif
            qrtmp = qrtmp + eva
	  endif
!
!----------------------------------------------------------!
!                          Updates                         !
!----------------------------------------------------------!
!
          hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) + (qctmp - hydromtrl(drop)%q(i,j,k))*dt1
          hydromtrm(rain)%q(i,j,k) = hydromtrm(rain)%q(i,j,k) + (qrtmp - hydromtrl(rain)%q(i,j,k))*dt1
          hydromtrm(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) + cnc*(nctmp - hydromtrl(drop)%n(i,j,k))*dt1
!
#ifdef ISENTROPIC
          esm(i,j,k) = esm(i,j,k) + (estmp - statel%es(i,j,k))*dt1
#endif
!
!  Store diagnostics 
!
          if (ldiag) then
#ifdef ISENTROPIC
            if (out_diagt) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + cint*dref(i,j,k)*(estmp - statel%es(i,j,k))*dt1
#endif
            if (out_diagr) diag(4)%qr(i,j,k) = diag(4)%qr(i,j,k) + cint*dref(i,j,k)*eva
          endif
!
        enddo
      enddo
    enddo
!
    endif
!
    if (verbose > 1) call write_debug('Terminating micro_k')

  end subroutine micro_k
!
!  ==================================================  
  subroutine aucr_k ( ff, rho, tt, qc, qr, auq )
!  ==================================================

    real, intent(in) :: rho, tt, qc, qr, ff
    real, intent(out) :: auq
    real :: auqr, auqi, psi, ks, dd
!
!----------------------------------------------------------!
!
!  Kessler autoconversion
!
!    auqr = max((ff*qc - qauto)/k1r, 0.)
!
!  Berry autoconversion
!
    auqr = 0.
    auqi = 0.
!
    if ( ff*qc > qmin(drop) ) then
      psi = 1000.*rho*ff*qc
      dd = 0.146 - 0.05964*log(xn_ccn0/2000.e6)
!
      auqr = 1.67e-5*psi*psi*(5. + 0.036*xn_ccn0/(1.e6*dd*psi))**(-1.) / rho
    endif
!
!  Snow autoconversion
!
    if ( tt > 273.15 - 15. ) then
      ks = k1s - 0.05333*k1s*(273.15 - tt)
    else if (tt <= 273.15 - 15. .and. tt > 273.15 - 30.) then
      ks = k1s - 0.05333*k1s*(tt - 273.15 + 30.)
    else
      ks = k1s
    endif
!
    if ( (1.-ff)*qc > qmin(drop) ) auqi = (1.-ff)*qc/ks
!
    auq = min(auqr + auqi, qc/dt)
  
  end subroutine aucr_k
!
!  ==================================================  
  subroutine ccr_k ( ff, rho, qc, qr, acc )
!  ==================================================

    real, intent(in) :: rho, ff, qc, qr
    real, intent(out) :: acc
    real :: mr, ms, clr, cli
!
!----------------------------------------------------------!
!
    clr = 0.
    cli = 0.
!
    mr = kmr*(rho*ff*qr/nr0)**bmr
    ms = kms*(rho*(1.-ff)*qr/ns0)**bms
!
    if ( ff*qr > qmin(rain) ) clr = ff*qr*k2r*effr*qc*mr**(-0.1667)
    if ( (1.-ff)*qr > qmin(rain) ) cli = (1.-ff)*qr*k2s*effs*qc*ms**(0.125)
    
    acc = min(clr + cli, qc/dt)
    
  end subroutine ccr_k
!
!  ==================================================  
  subroutine evap_k ( ff, rho, tt, ssw, ssi, qr, eva )
!  ==================================================

    real, intent(in) :: ff, rho, tt, qr, ssw, ssi
    real, intent(out) :: eva
    real :: gw, gi, mr, ms, fr, fs, evar, evai
!
!----------------------------------------------------------!
!
    evar = 0.
    evai = 0.
!
    gw = 1.e-7 / (2.2*tt/cal_esw(tt) + 220./tt)
    gi = 1.e-7 / (2.2*tt/cal_esi(tt) + 220./tt)
!
    mr = kmr*(rho*ff*qr/nr0)**bmr
    ms = kms*(rho*(1.-ff)*qr/ns0)**bms
!
    fr = avr + bvr*k5r*mr**(0.25)
    fs = avs + bvs*k5s*ms**(0.3125)
!
    if ( ff*qr > qmin(rain) ) evar = ff*qr*k4r*fr*gw*(ssw - 1.)*mr**(-0.6667)
    if ( (1.-ff)*qr > qmin(rain) ) evai = (1.-ff)*qr*k4s*fs*gi*(ssi - 1.)*ms**(-0.5)
!    
    eva = max(evar + evai, -qr/dt)
    
  end subroutine evap_k
!
!  ==================================================
  real function calc_vel_k ( rho, tt, qr )
!  ==================================================

    real, intent(in) :: rho, qr, tt
    real :: velqr, velqi, ff, mr, ms
!
!----------------------------------------------------------!
!
    velqr = 0.
    velqi = 0.
!
    ff = 1.
    if ( lfreeze > 0 ) ff = max(min( 1., (tt - tii) / (t00 - tii) ), 0.)
!
    if ( ff*qr > qmin(rain) ) velqr = max(k3r*(rho*ff*qr/nr0)**(0.125), 0.)
    if ( (1.-ff)*qr > qmin(rain) ) velqi = max(k3s*(rho*(1.-ff)*qr/ns0)**(0.0833), 0.)
    
    calc_vel_k = ff*velqr + (1.-ff)*velqi
    
    return
  
  end function calc_vel_k
!
!  ==================================================
  real function calc_avd_k ( ind, rho, tt, qr )
!  ==================================================

    real, intent(in) :: rho, qr, tt
    integer, intent(in) :: ind
    real :: kd, b, n, q
!
!----------------------------------------------------------!
!
    if (ind == 1) then
      q = qr
      n = xn_ccn0
      b = br
      kd = kdr
    else if (ind == 2) then
      q = qr
      if (lfreeze > 0) q = qr*max(min( 1., (tt - tii) / (t00 - tii) ), 0.)
      n = nr0
      b = br
      kd = kdr
    else if (ind == 3) then
      q = 0.
      if (lfreeze > 0) q = qr*(1.-max(min( 1., (tt - tii) / (t00 - tii) ), 0.))
      n = ns0
      b = bs
      kd = kds
    endif
!
    if ( q > qmin(drop) ) then
      calc_avd_k = kd*(rho*q/n)**(1./(b+1.))
    else
      calc_avd_k = 0.
    endif
    
    return
  
  end function calc_avd_k
!
!  ==================================================
  real function calc_n_k ( rho, tt, q )
!  ==================================================

    real, intent(in) :: rho, q, tt
    real :: ff, mr, ms, nr, ns
!
!----------------------------------------------------------!
!
    ff = 1.
    if (lfreeze > 0) ff = max(min( 1., (tt - tii) / (t00 - tii) ), 0.)
!
    mr = kmr*(rho*ff*q/nr0)**bmr         
    ms = kms*(rho*(1.-ff)*q/ns0)**bms      
!
    nr = 0.
    ns = 0.
    if ( ff*q > qmin(rain) ) nr = rho*ff*q/mr
    if ( (1.-ff)*q > qmin(rain) ) ns = rho*(1.-ff)*q/ms
! 
    calc_n_k = nr + ns 
!    
    return
  
  end function calc_n_k

  end module kessler
