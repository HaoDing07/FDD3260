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
!  MICRO_DIAG:
!	Micro_diag treats all the diagnostic microphysical
!       processes that need to be calculated after advection
!       and other physical phenomena. This is the case for
!       example for CCN nucleation or diffusional growth 
!       through saturation adjustment.
!
!	Subroutines called are found in micropack2.f90
!
!  Author:
!	Julien Savre, MISU
!
! ==============================================================

  module micro_diagnostics
!  
!  ---------------------------------------------
!
  USE shared_all
  USE allocation
  USE thermodynamics
  USE micropack
  USE sbpack
  USE kessler
  USE activation
  USE freezing
  USE aerosols 
  USE boundary_conditions
 
  IMPLICIT NONE
  
  private
  
  public :: diagnostic_micro, sat_adj, init_activ 
  
  contains
  
!	==========================================
	subroutine diagnostic_micro ( pressurel, windl, statel, hydromtrl, aero3dl ) 
!	==========================================
!
  type (atm_winds), intent(inout) :: windl
  type (atm_state), intent(inout) :: statel
  type (atm_pressure), intent(inout) :: pressurel
  type (hydrometeor), dimension(:), intent(inout) :: hydromtrl
  type(aero_3d), dimension(1:nmode), intent(inout)  :: aero3dl
!
  integer :: i, j, k, h
  real  :: xfkd, xfmu, exn, cpm, cvm, ss, ssi, pp, tnew, qtnew, nnew
  real  :: qnew(1:nhydro), ql, qi
  real, allocatable, dimension(:,:,:) :: tem, rh, dref
  real, allocatable, dimension(:,:,:) :: dqc, dnc, des
!
  if (verbose > 0) call write_debug('Start diagnostic_micro')
!
!----------------------------------------------------------!
!                  	Initialization                     !
!----------------------------------------------------------!
!
!  General allocations
!
  call alloc ( dref )
!
!  Reference density 
!
#if (defined CONSERVATIVE) && (!defined ANELASTIC)
  dref = pressurel%dens
#elif (defined CONSERVATIVE) && (defined ANELASTIC)
  do k = 1, nz
    dref(:,:,k) = den0(k)
  enddo
#else
  dref = 1.
#endif
!
!----------------------------------------------------------!
!                 Activation and Nucleation                !
!----------------------------------------------------------!
! 
#ifdef SEIFERT
  call activate_and_nucleate ( dref, windl, pressurel, statel, hydromtrl, aero3dl )
#endif
!
!----------------------------------------------------------!
!                  Saturation adjustment                   !
!----------------------------------------------------------!
!
#if (defined SAT_ADJ) || (!defined SEIFERT)
!
!  Other allocations
!
  call alloc ( tem )
  call alloc ( rh )
  call alloc ( dqc )
  call alloc ( dnc )
  call alloc ( des )
!
!  Recalculate temperature and RH
!
  call get_temperature ( pressurel, statel, hydromtrl, tem )
!
  call get_rh ( tem, pressurel%p, statel, hydromtrl, rh )
!
!  Starting saturation adjustment loops
!         
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
!
#ifdef ANELASTIC
	pp = p0(k)
#else
	pp = pressurel%p(i,j,k) + p0(k)
#endif
!
	tnew = tem(i,j,k)
	qtnew = statel%qt(i,j,k)
	nnew = hydromtrl(drop)%n(i,j,k)
	do h = 1, nhydro
  	  qnew(h) = hydromtrl(h)%q(i,j,k)
	enddo
!
        ss = rh(i,j,k)/100.
        ssi = ss*cal_qsw(tnew,pp)/cal_qsi(tnew,pp)
!
!  Saturation adjustment
!
#ifndef SEIFERT
        if ( (qnew(drop) > qmin(drop) .and. ss < 1.) .or. (ss >= 1. .and. tnew >= tii .and. nnew >= xnmin(drop)) .or. (ssi >= 1. .and. tnew < t00 .and. nnew >= xnmin(drop)) ) then
#else
        if ( (qnew(drop) > qmin(drop) .and. ss < 1.) .or. (ss >= 1. .and. nnew >= xnmin(drop)) ) then
#endif
!
          call sat_adj( pp, tnew, qtnew, qnew )
!
        endif
!
	dqc(i,j,k) = (qnew(drop) - hydromtrl(drop)%q(i,j,k)) / dt0
!
!  Tendencies
!
        xfkd = cal_xfkd (tnew, pp)
	xfmu = cal_xfmu (tnew)
!
  	if ( lndrop == 1 .and. dqc(i,j,k) < 0. ) then
          call sbevap ( k, drop, qnew(drop), nnew, xfkd, xfmu, dqc(i,j,k), dnc(i,j,k) )
        endif
!
#ifdef ISENTROPIC
        if ( lmicro > 1 ) qi = sum(qnew(ice:nhydro))
        call get_cp_cv ( qtnew, qnew(drop)+qnew(rain), qi, cpm, cvm )
   	exn = (pp/pref)**((cp_a-cv_a)/cp_a)
!
    	des(i,j,k) = cal_flv(tnew)*dqc(i,j,k) / (cpm*exn)
#endif
!
!  For Kessler, also diagnose number of precipitating particles
!
#ifndef SEIFERT
        if ( qnew(rain) > qmin(rain) ) then
          hydromtrl(rain)%n(i,j,k) = calc_n_k ( dref(i,j,k), tnew, qnew(rain) )
	else
	  hydromtrl(rain)%n(i,j,k) = 0.0
	endif
#endif
!
      end do
    end do
  end do
!
!----------------------------------------------------------!
!                      Updates and BC                      !
!----------------------------------------------------------!
!
!  Updates
!
#ifdef ISENTROPIC
    statel%es = statel%es + des*dt0
#endif
!
    hydromtrl(drop)%q = hydromtrl(drop)%q + dqc*dt0
    hydromtrl(drop)%n = hydromtrl(drop)%n + dnc*dt0
!
!  BCs
!
#ifdef ISENTROPIC
    call statebc ( statel%es )
#endif
!
    call hydrobc ( hydromtrl )
!
!----------------------------------------------------------!
!               Regeneration after adjustment              !
!----------------------------------------------------------!
!
    if ( aero_flg%reg ) then
      do k=1,nz
        do j=jt_start,jt_end
          do i=it_start,it_end
#ifdef NUC_CNT
	    call in_evap ( i, j, k, hydromtrl, nucin2, Tnew-273.15, dnc(i,j,k), 0., 0., 0., 0., 0., 0., 0., 0., 0. )
#endif
!
#ifdef AERO_ENABLE
            call regeneration ( i, j, k, dref(i,j,k), sum(hydromtrl%n(i,j,k)), aero3d, dnc(i,j,k) )
!
#if (defined AQCHEM_ENABLE) || (defined SOLIDCHEM_ENABLE)
            call evscaveng ( i, j, k, hydromtrl, dqc(i,j,k), 0., 0. )
#endif
#endif
          end do
        end do
      end do
    endif
!
!  BC
!
#ifdef AERO_ENABLE
    call aerobc ( aero3d )
#endif
!
!----------------------------------------------------------!
!                         Finalize                         !
!----------------------------------------------------------!
!
!  Diagnostics
!
    if (ldiag) then
      if (out_diagt.and.lmicro>0) diag(4)%pt = diag(4)%pt + dref*des
      if (out_diagl.and.lmicro>0) diag(4)%qc = diag(4)%qc + dref*dqc
      if (out_diagt.and.lmicro>0) diag(9)%pt = diag(9)%pt + dref*des
      if (out_diagl.and.lmicro>0) diag(9)%qc = diag(9)%qc + dref*dqc
!
      if (out_micro) then
        do k = 1, nz
          do j = jt_start, jt_end
	    do i = it_start, it_end
              diag(1)%micro(drop)%q(i,j,k) = diag(1)%micro(drop)%q(i,j,k) + dref(i,j,k)*dqc(i,j,k)
              diag(1)%micro(drop)%n(i,j,k) = diag(1)%micro(drop)%n(i,j,k) + dref(i,j,k)*dnc(i,j,k)
            enddo
	  enddo
        enddo
      endif
    endif
!
! Deallocate
!
    call dealloc ( rh )
    call dealloc ( dqc )
    call dealloc ( dnc )
    call dealloc ( des )
#endif
!
    call dealloc ( dref )
    call dealloc ( tem )
!
! ------------------------------------------------------------------
!
  if (verbose > 0) call write_debug('Terminate diagnostics_micro')
!	
  return
  end subroutine diagnostic_micro
!
!  =================================================
  subroutine sat_adj ( pp, tt, qt, qall )
!  ==================================================
!
! === Saturation adjustment: calculates liquid water mixing 
! === ratio diagnostically.
!
  integer  :: k
  integer, parameter :: kmax=6
  real, parameter :: small=1.e-8
!
  real :: qt, qall(nhydro), qs, ql, qi, qv, qcold, tt, pp, cpm, cvm
  real :: dqlold, dql, dqi, dte, ff
!
!----------------------------------------------------------!
!
!  Iterative method: update qc and qv only
!  We assume that diffusional growth as calculated by the saturation adjustment 
!  only affects cloud droplets, qr is held constant. Saturation adjustmen is thus
!  only done where qc > 0. During the overall first step, qc must be initialized.
!
  qcold = qall(drop)
  ql = qall(drop) + qall(rain)
  qi = sum(qall(ice:nhydro))
  qv = qt - sum(qall)
!
  call get_cp_cv (qt, ql, qi, cpm, cvm)
!
!  Saturation adjustment over water only
!
  k  = 1
  dqlold = 1.
  do while ( abs(dqlold) > small .and. k < kmax )
!
#ifdef SEIFERT
    qs = cal_qsw(tt, pp)
#else
    ff = min(max((t00 - tt) / (t00 - tii),0.),1.)
    qs = (1.-ff)*cal_qsw(tt, pp) + ff*cal_qsi(tt, pp)
#endif
!
    dql = (qv - qs) / (1. + cal_flv(tt)*cal_flv(tt)*cal_qsw(tt, pp) / (cpm*(cp_v-cv_v)*tt*tt))
!
    dql = max(dql, -qall(drop))
    dte = cal_flv(tt)*dql/cpm
!
    qall(drop) = qall(drop) + dql
    qv = qv - dql
    tt = tt + dte 
!
    k = k + 1
    dqlold = dql
    ql = qall(drop) + qall(rain)
    call get_cp_cv (qt, ql, qi, cpm, cvm)
!
  enddo
!
  dql = qall(drop) - qcold
!
!----------------------------------------------------------!
!
  return
  end subroutine sat_adj
!
! ==========================================
  subroutine init_activ ( pressurel, statel, hydromtrl, nucl, aero3dl )
! ==========================================
!
! ================================================================
!
!  Purpose:
!	Simple CCN activation scheme following Khvorostyanov & Curry (JAS, 2006)
!	Activation is calculated from state 1 (n-1 in LeapFrog), that is from the
!	state prior to integration.
!
! ================================================================
!
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
  type (atm_state) :: statel
  type (atm_pressure) :: pressurel  
  type (nuclei) :: nucl
  type (aero_3d), dimension(1:nmode)  :: aero3dl
!
  integer :: i, j, k, im, h
  real  :: dnc, dqc, dtt, m0
  real  :: sumh, p, tem, exn, qv, Told 
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: dref
!
  if (verbose > 1) call write_debug('Start init_activ')
!
!----------------------------------------------------------!
!                        Starting                          !
!----------------------------------------------------------!
!
!  Reference density
!
#ifdef CONSERVATIVE
#if (!defined ANELASTIC)
  dref = pressurel%dens
#else
  do k = 1, nz
    dref(:,:,k) = den0(k)
  enddo
#endif
#else
  dref = 1.
#endif
!
!  Update thermo
!
  call equation_of_state ( statel, hydromtrl, pressure_in=pressurel%p, thermo_out=thermo, conservative=.false. )
!
!----------------------------------------------------------!
!                    Initial activation                    ! 
!----------------------------------------------------------!
!
  m0 = 4./3.*pi*rho_l*(6.e-6)**3.
!
  do k=1,nz-1
    do j=jp_start,jp_end
      do i=ip_start,ip_end
!
#ifdef ANELASTIC
        p = p0(k)
#else
        p = pressurel%p(i,j,k) + p0(k)
#endif
        exn = (p/pref)**((cp_a-cv_a)/cp_a)
!
        sumh = 0.
        do h = 1, nhydro
          sumh = sumh + hydromtrl(h)%q(i,j,k)
        enddo
        qv = statel%qt(i,j,k) - sumh
        tem = thermo%T(i,j,k)
!
#ifdef AERO_ENABLE
	dnc = 0.
	do im = 1, nmode
	  dnc = dnc + aero3dl(im)%n(i,j,k)
	enddo
#else
	dnc = nucl%ccn(i,j,k)
#endif
!
	Told = 1e9
	do while ( abs(Told-Tem) > 0.01 )
	  Told = Tem
          dqc = min(m0*dnc,qv-cal_qsw(Tem,p))
          dtt = cal_flv(tem)/(cp_a)*dqc
	  Tem = Tem + dtt
	enddo
!
!----------------------------------------------------------!
!                         Updates                          !
!----------------------------------------------------------!
!
!  Activate
!
          if ( dqc > qmin(drop) .and. dnc > xnmin(drop) ) then
	    hydromtrl(drop)%n(i,j,k) = hydromtrl(drop)%n(i,j,k) + dnc
	    hydromtrl(drop)%q(i,j,k) = hydromtrl(drop)%q(i,j,k) + dqc
#ifdef ISENTROPIC
            statel%es(i,j,k) = statel%es(i,j,k) + dtt/exn
#endif
!
!  Diagnostics
!
#ifdef ISENTROPIC
            if (ldiag .and. out_diagt) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + dref(i,j,k)*dtt/exn/dt0
#endif
	    if (ldiag .and. out_diagl) diag(8)%qc(i,j,k) = diag(8)%qc(i,j,k) + dref(i,j,k)*dqc/dt0
            if (ldiag .and. out_micro) then
    	      diag(9)%micro(drop)%q(i,j,k) = diag(9)%micro(drop)%q(i,j,k) + dref(i,j,k)*dqc/dt0
	      diag(9)%micro(drop)%n(i,j,k) = diag(9)%micro(drop)%n(i,j,k) + dref(i,j,k)*dnc/dt0
            endif
!
	    if (ldiag .and. out_lwp) surf%cond(i,j) = surf%cond(i,j) + dref(i,j,k)*max(dqc,0.)*dzw(k)/dt0
	  endif
!
      enddo
    enddo
  enddo
!
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate init_activ')
!
  return
  end
!
! ===================================================
  subroutine activate_and_nucleate ( dref, windl, pressurel, statel, hydromtrl, aero3dl )
! ===================================================

  integer :: i, j, k, h
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,nz), intent(in) :: dref
  type(atm_winds), intent(in) :: windl
  type(atm_state), intent(inout) :: statel 
  type(atm_pressure), intent(inout) :: pressurel
  type(hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrl 
  type (aero_3d), dimension(1:nmode), intent(inout)  :: aero3dl
!
  real, dimension(:,:,:), allocatable :: ww
  type(atm_state) :: statem
  type(hydrometeor), dimension(1:nhydro) :: hydromtrm
  type (aero_3d), dimension(:), allocatable  :: aero3dm
#ifdef NUC_CNT
  type(nuclei_3d) :: nucinm
#endif
!
  if (verbose > 0) call write_debug('Enter activate_and_nucleate')
!
!-----------------------------------------------------------!
!
  call alloc ( ww )
  call alloc ( hydromtrm )
  call alloc ( statem ) 
#ifdef AERO_ENABLE
  allocate ( aero3dm(1:nmode) )
  call alloc ( aero3dm )
#endif
#ifdef NUC_CNT
  call alloc ( nucinm )	
#endif
!
!  Vertical wind
!
  do k=1,nz-1
    do j=jp_start,jp_end
      do i=ip_start,ip_end
        ww(i,j,k) = 0.5*(windl%w(i,j,k+1) + windl%w(i,j,k))
      end do
    end do
  end do
  ww(:,:,nz) = 0.5*(0. + windl%w(:,:,nz))
!
!  Update thermo
!
  call equation_of_state ( statel, hydromtrl, pressure_in=pressurel%p, thermo_out=thermo, conservative=.false. )
!
!----------------------------------------------------------!
!                Activation and nucleation                 !
!----------------------------------------------------------!
!
!  Activation
!
  if ( lndrop == 1 ) then 
#ifndef AERO_ENABLE
    call CCN_activ ( dref, thermo%T, ww, rad%dtnet, pressurel, statel, hydromtrl, nuc, statem, hydromtrm )
#else 
    call aero_activ ( dref, thermo%T, ww, rad%dtnet, pressurel, statel, hydromtrl, aero3dl, statem, hydromtrm, aero3dm )
#endif
!
#ifdef NUC_CNT
    call in_activ ( dref, thermo%T, pressurel, statel, hydromtrl, nucinl, nucinm )
#endif
  endif
!
!  Freezing
!
  if ( lmicro > 1 .and. lfreeze > 0 .and. time > ice_delay ) then
    call freez ( dref, ww, pressurel, statel, hydromtrl, statem, hydromtrm )
  endif
!
!----------------------------------------------------------!
!                          Updates                         !
!----------------------------------------------------------!
!
#ifdef ISENTROPIC
  statel%es = statel%es + dt0*statem%es
#endif
!
  do h = 1, nhydro
    hydromtrl(h)%q = hydromtrl(h)%q + dt0*hydromtrm(h)%q
  enddo
!
#ifdef SEIFERT
  if ( moments == 2 ) then
    if ( lndrop == 1 ) hydromtrl(drop)%n = hydromtrl(drop)%n + dt0*hydromtrm(drop)%n
!
    do h = 2, nhydro
      hydromtrl(h)%n = hydromtrl(h)%n + dt0*hydromtrm(h)%n
    enddo
  endif
#endif
!
#ifdef AERO_ENABLE
  if ( lndrop == 1 .and. aero_flg%any ) then
    do h = 1, nmode
      aero3dl(h)%n = aero3dl(h)%n + dt0*aero3dm(h)%n
!
      aero3dl(h)%m = aero3dl(h)%m + dt0*aero3dm(h)%m
!
      aero3dl(h)%ma = aero3dl(h)%ma + dt0*aero3dm(h)%ma
    enddo
  endif
#endif
!
!  Boundary conditions
!
#ifdef ISENTROPIC
  call statebc ( statel )
#endif
!
  call hydrobc ( hydromtrl )
!
#ifdef AERO_ENABLE
  if ( lndrop == 1 .and. aero_flg%any ) call aerobc ( aero3dl )
#endif
!
!  Limit hydrometeors
!
  call limit_hydro ( statel%qt, statel%es, hydromtrl, aero3dl )
!
!  One moment
!
#ifdef SEIFERT
  if ( moments == 1 ) call one_moment ( hydromtrl )
#endif
!
!  Deallocate
!
  call dealloc ( ww )
  call dealloc ( statem )
  call dealloc ( hydromtrm )
#ifdef AERO_ENABLE
  call dealloc ( aero3dm )
  deallocate ( aero3dm )
#endif
#ifdef NUC_CNT
  call dealloc ( nucinm )	
#endif
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate activate_and_nucleate')
! 
  return
!
  end subroutine activate_and_nucleate
!
end module micro_diagnostics
