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
!  MICROPACK2:
!	Microphysics2 treats evaporation/condensation and
!       sublimation/deposition processes using prognostic
!       supersaturation to advance qv, pt, qi and ni.
!
!  Author:
!	Julien Savre, MISU
!
! ==============================================================
!
  module seifert_beheng_ice
!
  USE shared_all
  USE sbpack
  USE micropack
  USE aerosols
  USE averages
  USE allocation
  USE thermodynamics
#ifdef SPMD
  USE mpi
#endif
!
  IMPLICIT NONE

  private
 
  logical :: with_cond = .true., no_cp = .false.

  public :: sb_ice

  contains

!	==========================================
	subroutine sb_ice ( dref, ww, dTdt, qt, pressurel, hydromtrl, aero3dl, statem, hydromtrm, aero3dm ) 
!	==========================================
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: dref, ww, dTdt, qt
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(1:nhydro), intent(in) :: hydromtrl
  type (aero_3d), dimension(1:nmode), intent(inout) :: aero3dl
!
  type (hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrm  
  type (atm_state), intent(inout) :: statem
  type (aero_3d), dimension(1:nmode) :: aero3dm
!
  integer :: i, j, k, h, im, iel
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: evapr
  real :: qc, qr, qi, qg, qs, qh, nc, nr, ni, ng, ns, nh, Tem, pp, qv,  &
	     qc0, qr0, qi0, qg0, qs0, qh0, nc0, nr0, ni0, ng0, ns0, nh0,  &
	     sw, si, es0, xfkd, xfmu, cd, cdr, frimi, frims, frimr, 	&
	     cli, clin, clii, ais, aisn, cri, crin, crii,       &
	     cis, ciss, cisn, clg, clgg, clgn, cssn, css,       &
	     cls, clss, clsn, crg, crgg, crgn, csg, csgg, csgn, &
	     crs, crsn, crss, xmir, xmirn, xmgr, xmgrn, xmsr,   &
	     xmsrn, xmggn, xmhhn, xmiin, xmssn,			&
	     depi, depg, deps, dqc, dqr, cghs, cghsn,		&
	     dqi, dqg, dqs, dnc,  dnr,  dni, dng, dns, xmrim, 	&
	     sheli, sheri, shels, shers, shelg, sherg, rho,  	&
	     exn, esi, esw, qsw, qsi, des, dp, hm, hmn,   	&
	     cc, evqc, evqr, evnc, evnr, sbqi, sbqg, 		&
	     sbqs, sbni, sbng, sbns, shegh, cdh, deph,   	&
	     sbqh, sbnh, xmhr, xmhrn, xmhh, clh, clhh, clhn,    &
	     crh, crhh, crhn, cih, cihh, cihn, csh, cshh, cshn, &
	     hf, hfs, hfsn, shelh, dnh, dqh, q, Tdg, Tdh, sumdq,&
             ql, cp, cv, gamma, evm(nz), x 
! 
#if ( defined SPMD )      
  integer :: ierr, status(MPI_STATUS_SIZE)
#endif
!
  if (verbose > 1) call write_debug('Starting sb_ice')
!
!----------------------------------------------------------!
!                      External loop                       !
!----------------------------------------------------------!
!
  x = 1.e9
  evapr = 0.0
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
        if ( (hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. hydromtrl(drop)%n(i,j,k) > xnmin(drop)) .or. &
             (hydromtrl(rain)%q(i,j,k) > qmin(rain) .and. hydromtrl(rain)%n(i,j,k) > xnmin(rain)) .or. &
             (hydromtrl(ice)%q(i,j,k) > qmin(ice) .and. hydromtrl(ice)%n(i,j,k) > xnmin(ice)) .or. &
             (lmicro > 2 .and. hydromtrl(grau)%q(i,j,k) > qmin(grau) .and. hydromtrl(grau)%n(i,j,k) > xnmin(grau)) .or. &
             (lmicro > 2 .and. hydromtrl(snow)%q(i,j,k) > qmin(snow) .and. hydromtrl(snow)%n(i,j,k) > xnmin(snow)) ) then   
!
#ifdef ANELASTIC
	pp = p0(k)
        rho = den0(k)
#else
        pp  = pressurel%p(i,j,k) + p0(k)
        rho = pressurel%dens(i,j,k)
#endif
!
	qv  = qt(i,j,k)
        ql  = 0.
        qi  = 0.
	if (lmicro > 0) ql = hydromtrl(drop)%q(i,j,k) + hydromtrl(rain)%q(i,j,k)
	if (lmicro > 1) qi = hydromtrl(ice)%q(i,j,k)
	if (lmicro > 2) qi = qi + hydromtrl(grau)%q(i,j,k) + hydromtrl(snow)%q(i,j,k)
	if (lmicro > 3) qi = qi + hydromtrl(hail)%q(i,j,k)
        qv = qv - ql - qi
!
	Tem = thermo%T(i,j,k)
	exn = thermo%exn(i,j,k)
	es0 = cal_esw (273.15)
	esw = cal_esw (Tem)
	esi = cal_esi (Tem)
        qsw = cal_qsw (Tem, pp)
        qsi = cal_qsi (Tem, pp)
	sw  = qv - qsw
	si  = qv - qsi
!
        call get_cp_cv (qv, ql, qi, cp, cv)
!
	xfkd = cal_xfkd (Tem, pp)
	xfmu = cal_xfmu (Tem)
	gamma = 1./qv + cal_flv(Tem)**2./((cp_v-cv_v)*cp*Tem*Tem) 
!
!----------------------------------------------------------!
!                Initialize microphysics                   !
!----------------------------------------------------------!
!
	cli   = 0.0
	clii  = 0.0
	clin  = 0.0
	cls   = 0.0
	clss  = 0.0
	clsn  = 0.0
	cri   = 0.0
	crii  = 0.0
	crin  = 0.0
	crs   = 0.0
	crss  = 0.0
	crsn  = 0.0
	ais   = 0.0
	aisn  = 0.0
	cis   = 0.0
	ciss  = 0.0
	cisn  = 0.0
	clg   = 0.0
	clgg  = 0.0
	clgn  = 0.0
	crg   = 0.0
	crgg  = 0.0
	crgn  = 0.0
	csg   = 0.0
	csgg  = 0.0
	csgn  = 0.0
	css   = 0.0
	cssn  = 0.0
	cghs  = 0.0
	cghsn = 0.0
	clh   = 0.0
	clhh  = 0.0
	clhn  = 0.0
	crh   = 0.0
	crhh  = 0.0
	crhn  = 0.0
	cih   = 0.0
	cihh  = 0.0
	cihn  = 0.0
	csh   = 0.0
	cshh  = 0.0
	cshn  = 0.0
	hf    = 0.0
	hfs   = 0.0
	hfsn  = 0.0
	sheli = 0.0
	sheri = 0.0
	shels = 0.0
	shers = 0.0
	shelg = 0.0
	sherg = 0.0
	shelh = 0.0
	shegh = 0.0
	hm    = 0.0
	hmn   = 0.0
	cd    = 0.0
	cdr   = 0.0
	cdh   = 0.0
	depg  = 0.0
	depi  = 0.0  
	deps  = 0.0
	deph  = 0.0
	dqc   = 0.0
	dqr   = 0.0
	dqi   = 0.0
	dqg   = 0.0
	dqs   = 0.0
	dqh   = 0.0
	des   = 0.0
	dp    = 0.0
	qc    = 0.0
	qr    = 0.0
	qi    = 0.0
	qg    = 0.0
	qs    = 0.0
        qh    = 0.0
	nc    = 0.0
	nr    = 0.0
	ni    = 0.0
	ng    = 0.0
	ns    = 0.0
        nh    = 0.0
	dnc   = 0.0
	dnr   = 0.0
	dni   = 0.0
	dng   = 0.0
	dns   = 0.0
	dnh   = 0.0
	evqc  = 0.0
	evnc  = 0.0
	evqr  = 0.0
	evnr  = 0.0
	sbqi  = 0.0
	sbni  = 0.0
	sbqg  = 0.0
	sbng  = 0.0
	sbqs  = 0.0
	sbns  = 0.0
	sbqh  = 0.0
	sbnh  = 0.0
	xmir  = 0.0
	xmirn = 0.0
	xmiin = 0.0
	xmgr  = 0.0
	xmgrn = 0.0
	xmggn = 0.0
	xmsr  = 0.0
	xmsrn = 0.0
	xmssn = 0.0
	xmhr  = 0.0
	xmhrn = 0.0
	xmhhn = 0.0
	frimi = 0.0
	frims = 0.0
	frimr = 0.0
	xmrim = 0.0
	Tdg   = 0.0
	Tdh   = 0.0
!
        qc = hydromtrl(drop)%q(i,j,k)
        qr = hydromtrl(rain)%q(i,j,k)
	qi = hydromtrl(ice)%q(i,j,k)
        nc = hydromtrl(drop)%n(i,j,k)
        nr = hydromtrl(rain)%n(i,j,k)
	ni = hydromtrl(ice)%n(i,j,k)
	if (lmicro > 2) qg = hydromtrl(grau)%q(i,j,k)
	if (lmicro > 2) qs = hydromtrl(snow)%q(i,j,k)
	if (lmicro > 3) qh = hydromtrl(hail)%q(i,j,k)
	if (lmicro > 2) ng = hydromtrl(grau)%n(i,j,k)
	if (lmicro > 2) ns = hydromtrl(snow)%n(i,j,k)
	if (lmicro > 3) nh = hydromtrl(hail)%n(i,j,k)
!
	qc0 = qc
	qr0 = qr
	qi0 = qi
	qg0 = qg
	qs0 = qs
	qh0 = qh
	nc0 = nc
	nr0 = nr
	ni0 = ni
	ng0 = ng
	ns0 = ns
	nh0 = nh
!
!----------------------------------------------------------!
!               Cold: deposition/sublimation               !
!----------------------------------------------------------!
!
!  Evaporation:
!
        if (with_cond) then
!
#ifndef SAT_ADJ
          if(qc0 > qmin(drop) .and. nc0 > xnmin(drop))then
            call sbcond (k, drop, qsw, esw, Tem, xfkd, xfmu, qc, nc, cd)
          endif
#endif
! 
          if(qr0 > qmin(rain) .and. nr0 > xnmin(rain)) then
            call sbcond (k, rain, qsw, esw, Tem, xfkd, xfmu, qr, nr, cdr)
          endif
! 
!  Deposition of ice:
! 
          if (qi0 > qmin(ice) .and. ni0 > xnmin(ice)) then
            call sbdep (k, ice, qsi, esi, Tem, xfkd, xfmu, qi, ni, depi)
          endif
! 
!  Deposition of graupel and snow:
! 
          if (lmicro > 2) then
	    if (qg0 > qmin(grau) .and. ng0 > xnmin(grau)) then
              call sbdep (k, grau, qsi, esi, Tem, xfkd, xfmu, qg, ng, depg)
            endif
! 
            if (qs0 > qmin(snow) .and. ns0 > xnmin(snow)) then
              call sbdep (k, snow, qsi, esi, Tem, xfkd, xfmu, qs, ns, deps)
            endif
          endif
!
!  Resolved saturation
!
        call resolved_saturation ( k, cd, cdr, 0., depi, depg, deps, 0., pp, Tem, ww(i,j,k), dTdt(i,j,k), cp, cv,                     &
				   exn, sw, qsw, si, qsi, qv, qc, qr, qi, qg, qs, qh, evqc, evqr, sbqi, sbqg, sbqs, sbqh, des )
!
!  Evaporation
!
#ifndef SAT_ADJ
        if(evqc < 0.0 .and. lndrop == 1 .and. moments == 2)then  
          call sbevap (k, drop, qc, nc, xfkd, xfmu, evqc, evnc)
	endif
#endif
!
	if (evqr < 0.0 .and. moments == 2) then
          evapr(i,j,k) = evqr
          call sbevap (k, rain, qr, nr, xfkd, xfmu, evqr, evnr)
        endif
!
!  Evaporation of ice
!
	if (sbqi < 0.0 .and. moments == 2) then
          call sbevap (k, ice, qi, ni, xfkd, xfmu, sbqi, sbni)
	endif
!
!  Evaporation of graupel
!
	if (sbqg < 0.0 .and. lmicro > 2 .and. moments == 2) then
          call sbevap (k, grau, qg, ng, xfkd, xfmu, sbqg, sbng)	  
  	endif
!
!  Evaporation of snow
!
	if (sbqs < 0.0 .and. lmicro > 2 .and. moments == 2) then
          call sbevap (k, snow, qs, ns, xfkd, xfmu, sbqs, sbns)
	endif
!
!  Evaporation of hail
!
	if (sbqh < 0.0 .and. lmicro > 3 .and. moments == 2) then
          call sbevap (k, hail, qh, nh, xfkd, xfmu, sbqh, sbnh)
	endif
!
	if (out_diagt) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + dref(i,j,k)*cint*des
	if (out_diagl) diag(4)%qc(i,j,k) = diag(4)%qc(i,j,k) + dref(i,j,k)*cint*evqc
	if (out_diagr) diag(4)%qr(i,j,k) = diag(4)%qr(i,j,k) + dref(i,j,k)*cint*evqr
	if (out_diagi) diag(4)%qi(i,j,k) = diag(4)%qi(i,j,k) + dref(i,j,k)*cint*sbqi
!
        if (ldtmic .and. (cd + cdr) > 1.e-12) x = max(min(x,1./(gamma*(cd+cdr))),1.e-12)
!
        endif
!
!----------------------------------------------------------!
!            Full cold Microphysical processes	           !
!----------------------------------------------------------!
!	
        if ( lmicro > 2 ) then
!
          if (lrime) frimr = 1.
!
!  Collection of cloud drop by ice (riming): c + i --> g or r
!
          if (qc > qmin(drop) .and. nc > xnmin(drop) .and. qi > qmin(ice) .and. ni > xnmin(ice)) then
            call sbcollec (ice, drop, rho, Tem, qi, ni, qc, nc, clii, cli, clin)
!
	    call sbrime (ice, Tem, qi, ni, frimi)
!
	    if (Tem < 273.15) sheli = 1.
!
	    call hallett_mossop (Tem, qc, nc, sheli*frimi*cli, hm, hmn)
          endif
!
!  Collection of cloud drop by snow (riming): c + s --> g or r
!
          if (qc > qmin(drop) .and. nc > xnmin(drop) .and. qs > qmin(snow) .and. ns > xnmin(snow)) then
            call sbcollec (snow, drop, rho, Tem, qs, ns, qc, nc, clss, cls, clsn)
!
	    call sbrime (snow, Tem, qs, ns, frims)
!
	    if (Tem < 273.15) shels = 1.
!
	    call hallett_mossop (Tem, qc, nc, shels*frims*cls, hm, hmn)
          endif
!
!  Collection of cloud drop by graupel: c + g --> g or r
!
          if (qc > qmin(drop) .and. nc > xnmin(drop) .and. qg > qmin(grau) .and. ng > xnmin(grau))then
            call sbcollec (grau, drop, rho, Tem, qg, ng, qc, nc, clgg, clg, clgn)
!
	    if (Tem < 273.15) shelg = 1.
!
	    call hallett_mossop (Tem, qc, nc, shelg*clg, hm, hmn)
          endif
!
!  Collection of rain by graupel: r + g --> g or r
!
          if (qr > qmin(rain) .and. nr > xnmin(rain) .and. qg > qmin(grau) .and. ng > xnmin(grau))then
            call sbcollec (rain, grau, rho, Tem, qr, nr, qg, ng, crg, crgg, crgn)
!
	    if (Tem < 273.15) sherg = 1.
          endif
!
!  Collection of rain by ice (riming): r + i --> g or r
!
          if (qr > qmin(rain) .and. nr > xnmin(rain) .and. qi > qmin(ice) .and. ni > xnmin(ice)) then
            call sbcollec (rain, ice, rho, Tem, qr, nr, qi, ni, cri, crii, crin)
!
	    if (Tem < 273.15) sheri = 1.
          endif
!
!  Collection of rain by snow (riming): r + s --> g or r
!
          if (qr > qmin(rain) .and. nr > xnmin(rain) .and. qs > qmin(snow) .and. ns > xnmin(snow)) then
            call sbcollec (rain, snow, rho, Tem, qr, nr, qs, ns, crs, crss, crsn)
!
	    if (Tem < 273.15) shers = 1.
          endif
!
!  autoconversion from ice to snow: i + i --> s
!
          if(qi > qmin(ice) .and. ni > xnmin(ice))then
            call sbcollec (ice, ice, rho, Tem, qi, ni, qi, ni, ais, ais, aisn)
          endif
!
! Collection of ice by snow: i + s --> s
!
          if(qi > qmin(ice) .and. ni > xnmin(ice) .and. qs > qmin(snow) .and. ns > xnmin(snow))then
	    call sbcollec (snow, ice, rho, Tem, qs, ns, qi, ni, ciss, cis, cisn)
          endif	
!
!  Self-collection of snow: s + s --> s
!
	  if (qs > qmin(snow) .and. ns > xnmin(snow)) then
	    call sbcollec (snow, snow, rho, Tem, qs, ns, qs, ns, css, css, cssn)
	  endif
!
! Collection of snow by graupel: g + s --> g
!
          if(qs > qmin(snow) .and. ns > xnmin(snow) .and. qg > qmin(grau) .and. ng > xnmin(grau))then
            call sbcollec (grau, snow, rho, Tem, qg, ng, qs, ns, csgg, csg, csgn)
          endif	
!
        endif
!
!----------------------------------------------------------!
!          	    Hail Microphysics	                   !
!----------------------------------------------------------!
!		
        if ( lmicro > 3 ) then
!
!  Conversion graupel to hail: g + l --> h
!
          if (Tem < 273.15 .and. qg > qmin(grau) .and. ng > xnmin(grau))then
	    call fcaltdi (k, grau, Tem, qc + qr, esi, xfkd, xfmu, qg, ng, shelg*clg + sherg*crg, 0., Tdg, shegh)
!
	    if ( shegh > 0. ) call cshed (grau, Tem, qc+qr, qg, ng, shelg*clg + sherg*crg, 0., cghs, cghsn) 
	  endif
!
!  Collection of cloud droples: h + c --> h
!
          if (qc > qmin(drop) .and. nc > xnmin(drop) .and. qh > qmin(hail) .and. nh > xnmin(hail))then
            call sbcollec (hail, drop, rho, Tem, qh, nh, qc, nc, clhh, clh, clhn)
	  endif
!
!  Collection of rain: h + r --> h
!
          if (qr > qmin(rain) .and. nr > xnmin(rain) .and. qh > qmin(hail) .and. nh > xnmin(hail))then
            call sbcollec (hail, rain, rho, Tem, qh, nh, qr, nr, crhh, crh, crhn)
	  endif	  
!
!  Collection of ice: h + i --> h
!
          if (qi > qmin(ice) .and. ni > xnmin(ice) .and. qh > qmin(hail) .and. nh > xnmin(hail))then
            call sbcollec (hail, ice, rho, Tem, qh, nh, qi, ni, cihh, cih, cihn)
	  endif	  
!
!  Collection of snow: h + s --> h
!
          if (qs > qmin(snow) .and. ns > xnmin(snow) .and. qh > qmin(hail) .and. nh > xnmin(hail))then
            call sbcollec (hail, snow, rho, Tem, qh, nh, qs, ns, cshh, csh, cshn)
	  endif	  
! 
!  Dry or wet growth
!
          if (qh > qmin(hail) .and. nh > xnmin(hail)) then
	    call fcaltdi (k, hail, Tem, qc + qr, esi, xfkd, xfmu, qh, nh, clh + crh, cih + csh, Tdh, shelh)
!
	    hf = clh + crh + cih + csh
!
            call sbcond (k, hail, qsw, esw, Tem, xfkd, xfmu, qh, nh, cdh)
!
            call sbdep (k, hail, qsi, esi, Tem, xfkd, xfmu, qh, nh, deph)
!
	    cdh  = shelh*cdh/qsw 
	    deph = (1. - shelh)*deph/qsi
!
	    hf = clh + crh + cih + csh
!
	    if ( shelh > 0. ) call cshed (hail, Tem, qc+qr, qh, nh, clh + crh, cih + csh, hfs, hfsn)
          endif
!
	endif
!
!----------------------------------------------------------!
!                Phase transitions: melting                !
!	   Includes melting enhancement by riming          !
!----------------------------------------------------------!
!
!  Melting of ice to rain: i --> r
!
	if (Tem >= 273.15 .and. qi > qmin(ice) .and. ni > xnmin(ice)) then
          call sbmelt (k, ice, qi, ni, Tem, esw, es0, xfkd, xfmu, xmir, xmirn, xmiin)
!
          call sbrimemelt (Tem, clii, 0., xmrim)
!
          xmir = xmir + xmrim
  	endif
!
!  Melting of graupel and snow to rain: g --> r; s --> r
!
        if (lmicro > 2) then
  	  if (Tem >= 273.15 .and. qg > qmin(grau) .and. ng > xnmin(grau)) then
	    call sbmelt (k, grau, qg, ng, Tem, esw, es0, xfkd, xfmu, xmgr, xmgrn, xmggn)
!  
            call sbrimemelt (Tem, clgg, crgg, xmrim)
!
	    xmgr = xmgr + xmrim
  	  endif
!
          if (Tem >= 273.15 .and. qs > qmin(snow) .and. ns > xnmin(snow)) then
	    call sbmelt (k, snow, qs, ns, Tem, esw, es0, xfkd, xfmu, xmsr, xmsrn, xmssn)
!
            call sbrimemelt (Tem, clss, crss, xmrim)
!
	    xmsr = xmsr + xmrim 
  	  endif
	endif
!
!  Melting of hail to rain: h --> r
!
        if (lmicro > 3) then
  	  if (Tem >= 273.15 .and. qh > qmin(hail) .and. nh > xnmin(hail)) then
	    call sbmelt (k, hail, qh, nh, Tem, esw, es0, xfkd, xfmu, xmhr, xmhrn, xmhhn)
!  
            call sbrimemelt (Tem, clhh, crhh, xmrim)
!
	    xmhr = xmhr + xmrim
  	  endif
	endif
!
!  Energy
!
#ifdef ISENTROPIC
	if (Tem >= 273.15) des = des - cal_flm(Tem)*(xmir + xmgr + xmsr + xmhr + (1.-sheri)*crii) / (cp*exn)
#endif
!
!----------------------------------------------------------!
!             Update microphysics source terms             !
!----------------------------------------------------------!
!
!  Q 
!
        dqc = evqc - cli - cls - clg - clh
	dqr = evqr + xmir + xmgr + xmsr + xmhr - sheri*cri + (1.-sheri)*crii - shers*crs + (1.-shers)*crss + (1.-sheli)*(cli + clii) + (1.-shels)*(cls + clss) + (1.-shelg)*clg - sherg*crg - crh + shegh*cghs + shelh*hfs
	dqi = sbqi - xmir + sheli*(1.-frimi)*(cli + clii) - clii + sheri*(1.-frimr)*(cri + crii) - crii - 2.*ais - cis + hm - cih
	dqg = sbqg - xmgr + sheli*frimi*(cli + clii) + shels*frims*(cls + clss) + frimr*(sheri*(cri + crii) + shers*(crs + crss)) + csg + (1.-shegh)*(shelg*clg + sherg*crg) - shegh*(shelg*clgg + sherg*crgg) - hm
	dqs = sbqs - xmsr + shels*(1.-frims)*(cls + clss) - clss + shers*(1.-frimr)*(crs + crss) - crss + 2.*ais + cis - csg - csh
	dqh = sbqh - xmhr + shegh*(shelg*(clg + clgg) + sherg*(crg + crgg) - cghs) + hf - shelh*hfs
!
        hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) + dqc
        hydromtrm(rain)%q(i,j,k) = hydromtrm(rain)%q(i,j,k) + dqr
        hydromtrm(ice)%q(i,j,k)  = hydromtrm(ice)%q(i,j,k) + dqi
        if (lmicro > 2) then
          hydromtrm(grau)%q(i,j,k) = hydromtrm(grau)%q(i,j,k) + dqg
          hydromtrm(snow)%q(i,j,k) = hydromtrm(snow)%q(i,j,k) + dqs
        endif
        if (lmicro > 3) then
          hydromtrm(hail)%q(i,j,k) = hydromtrm(hail)%q(i,j,k) + dqh
        endif
!
!  N 
!
	if ( moments == 2 ) then
	  dnc = evnc - clin - clsn - clgn - clhn
	  dnr = evnr + xmirn + xmgrn + xmsrn + xmhrn - sheri*crin - shers*crsn - sherg*crgn + (1.-sheli)*clin + (1.-shels)*clsn + (1.-shelg)*clgn - crhn + shegh*cghsn + shelh*hfs
	  dni = sbni - xmiin + (sheli*(1.-frimi) - 1.)*clin + (sheri*(1.-frimr) - 1.)*crin - 2.*aisn - cisn + hmn - cihn
	  dng = sbng - xmggn + sheli*frimi*clin + shels*frims*clsn + frimr*(sheri*crin + shers*crsn) - shegh*(shelg*clgn + sherg*crgn)
	  dns = sbns - xmssn + (shels*(1.-frims) - 1.)*clsn + (shers*(1.-frimr) - 1.)*crsn + aisn - cssn - csgn - cshn
	  dnh = sbnh - xmhhn + shegh*(shelg*clgn + sherg*crgn)
!
          if (lndrop == 1) hydromtrm(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) + dnc
          hydromtrm(rain)%n(i,j,k) = hydromtrm(rain)%n(i,j,k) + dnr
          hydromtrm(ice)%n(i,j,k)  = hydromtrm(ice)%n(i,j,k) + dni
          if (lmicro > 2) then
            hydromtrm(grau)%n(i,j,k) = hydromtrm(grau)%n(i,j,k) + dng 
            hydromtrm(snow)%n(i,j,k) = hydromtrm(snow)%n(i,j,k) + dns
          endif
          if (lmicro > 3) then
            hydromtrm(hail)%n(i,j,k) = hydromtrm(hail)%n(i,j,k) + dnh
          endif
	endif
!
!  Liquid water mass W
!
	if ( lmicro > 3 .and. moments == 2 ) then
	  dqg = (min(sbqg,0.) - xmgr - shegh*(shelg*clgg + sherg*crgg))*hydromtrl(grau)%w(i,j,k)/max(qg,qg_min) + (1.-shegh)*(shelg*clg + sherg*crg) + sheli*frimi*cli + shels*frims*cls + sheri*cri + shers*crs
	  dqs = (min(sbqs,0.) - xmsr - shels*frims*clss - shers*crss - csg - csh)*hydromtrl(snow)%w(i,j,k)/max(qs,qs_min) + sheli*(1.-frimi)*cli + shels*(1.-frims)*cls
	  dqh = (min(sbqh,0.) - xmhr)*hydromtrl(hail)%w(i,j,k)/max(qh,qh_min) + hf - shelh*hfs
!
          hydromtrm(grau)%w(i,j,k) = hydromtrm(grau)%w(i,j,k) + max(dqg,-hydromtrl(grau)%w(i,j,k)/dt0)
          hydromtrm(snow)%w(i,j,k) = hydromtrm(snow)%w(i,j,k) + max(dqs,-hydromtrl(snow)%w(i,j,k)/dt0)
          hydromtrm(hail)%w(i,j,k) = hydromtrm(hail)%w(i,j,k) + max(dqh,-hydromtrl(hail)%w(i,j,k)/dt0)
	endif
!
!  Energy 
!
#ifdef ISENTROPIC
	statem%es(i,j,k) = statem%es(i,j,k) + des
#endif
!
!----------------------------------------------------------!
!                    Store diagnostics                     !
!----------------------------------------------------------!
!
	if ( ldiag .and. out_lwp ) surf%cond(i,j) = surf%cond(i,j) + dref(i,j,k)*cint*(max(evqc,0.) + max(evqr,0.) + max(sbqi,0.) + max(sbqg,0.) + max(sbqs,0.))*dzw(k)
	if ( ldiag .and. out_lwp ) surf%evap(i,j) = surf%evap(i,j) + dref(i,j,k)*cint*(min(evqc,0.) + min(evqr,0.) + min(sbqi,0.) + min(sbqg,0.) + min(sbqs,0.))*dzw(k)
!
        if ( ldiag .and. out_micro ) then
    	  diag(1)%micro(drop)%q(i,j,k) = diag(1)%micro(drop)%q(i,j,k) + dref(i,j,k)*cint*evqc
    	  diag(4)%micro(drop)%q(i,j,k) = diag(4)%micro(drop)%q(i,j,k) - dref(i,j,k)*cint*cli
    	  diag(5)%micro(drop)%q(i,j,k) = diag(5)%micro(drop)%q(i,j,k) - dref(i,j,k)*cint*clg
    	  diag(6)%micro(drop)%q(i,j,k) = diag(6)%micro(drop)%q(i,j,k) - dref(i,j,k)*cint*cls
!
	  diag(1)%micro(rain)%q(i,j,k) = diag(1)%micro(rain)%q(i,j,k) + dref(i,j,k)*cint*evqr
	  diag(4)%micro(rain)%q(i,j,k) = diag(4)%micro(rain)%q(i,j,k) + dref(i,j,k)*cint*(xmir + xmgr + xmsr)
	  diag(5)%micro(rain)%q(i,j,k) = diag(5)%micro(rain)%q(i,j,k) + dref(i,j,k)*cint*((1.-sheli)*(cli + clii) + (1.-sheri)*crii - sheri*cri)
	  diag(6)%micro(rain)%q(i,j,k) = diag(6)%micro(rain)%q(i,j,k) + dref(i,j,k)*cint*((1.-shelg)*clg - sherg*crg)
	  diag(7)%micro(rain)%q(i,j,k) = diag(7)%micro(rain)%q(i,j,k) + dref(i,j,k)*cint*((1.-shels)*(cls + clss) + (1.-shers)*crss - shers*crs)
!
    	  diag(1)%micro(ice)%q(i,j,k) = diag(1)%micro(ice)%q(i,j,k) + dref(i,j,k)*cint*sbqi
    	  diag(2)%micro(ice)%q(i,j,k) = diag(2)%micro(ice)%q(i,j,k) - dref(i,j,k)*cint*xmir
    	  diag(3)%micro(ice)%q(i,j,k) = diag(3)%micro(ice)%q(i,j,k) - 2.*dref(i,j,k)*cint*ais
    	  diag(4)%micro(ice)%q(i,j,k) = diag(4)%micro(ice)%q(i,j,k) + dref(i,j,k)*cint*(sheli*(1.-frimi)*(cli + clii) - clii)
    	  diag(5)%micro(ice)%q(i,j,k) = diag(5)%micro(ice)%q(i,j,k) + dref(i,j,k)*cint*(sheri*(1.-frimr)*(cri + crii) - crii)
    	  diag(6)%micro(ice)%q(i,j,k) = diag(6)%micro(ice)%q(i,j,k) - dref(i,j,k)*cint*cis
    	  diag(7)%micro(ice)%q(i,j,k) = diag(7)%micro(ice)%q(i,j,k) + dref(i,j,k)*cint*hm
!	  
	  if ( lmicro > 2 ) then
	    diag(1)%micro(grau)%q(i,j,k) = diag(1)%micro(grau)%q(i,j,k) + dref(i,j,k)*cint*sbqg
    	    diag(2)%micro(grau)%q(i,j,k) = diag(2)%micro(grau)%q(i,j,k) - dref(i,j,k)*cint*xmgr
    	    diag(3)%micro(grau)%q(i,j,k) = diag(3)%micro(grau)%q(i,j,k) + dref(i,j,k)*cint*(sheli*frimi*(cli + clii) + shels*frims*(cls + clss))
    	    diag(4)%micro(grau)%q(i,j,k) = diag(4)%micro(grau)%q(i,j,k) + dref(i,j,k)*cint*(sheri*frimr*(cri + crii) + shers*frimr*(crs + crss))
    	    diag(5)%micro(grau)%q(i,j,k) = diag(5)%micro(grau)%q(i,j,k) + dref(i,j,k)*cint*((1.-shegh)*(shelg*clg + sherg*crg) - shegh*(shelg*clgg + sherg*crgg))
    	    diag(6)%micro(grau)%q(i,j,k) = diag(6)%micro(grau)%q(i,j,k) + dref(i,j,k)*cint*csg
	    diag(7)%micro(grau)%q(i,j,k) = diag(7)%micro(grau)%q(i,j,k) - dref(i,j,k)*cint*hm
!	    
	    diag(1)%micro(snow)%q(i,j,k) = diag(1)%micro(snow)%q(i,j,k) + dref(i,j,k)*cint*sbqs
    	    diag(2)%micro(snow)%q(i,j,k) = diag(2)%micro(snow)%q(i,j,k) - dref(i,j,k)*cint*xmsr
	    diag(3)%micro(snow)%q(i,j,k) = diag(3)%micro(snow)%q(i,j,k) + 2.*dref(i,j,k)*cint*ais
	    diag(4)%micro(snow)%q(i,j,k) = diag(4)%micro(snow)%q(i,j,k) + dref(i,j,k)*cint*(shels*(1.-frims)*(cls + clss) - clss)
	    diag(5)%micro(snow)%q(i,j,k) = diag(5)%micro(snow)%q(i,j,k) + dref(i,j,k)*cint*(shers*(1.-frimr)*(crs + crss) - crss)
	    diag(6)%micro(snow)%q(i,j,k) = diag(6)%micro(snow)%q(i,j,k) - dref(i,j,k)*cint*csg
	    diag(7)%micro(snow)%q(i,j,k) = diag(7)%micro(snow)%q(i,j,k) + dref(i,j,k)*cint*cis
	  endif
!
	  if ( lmicro > 3 ) then
	    diag(1)%micro(hail)%q(i,j,k) = diag(1)%micro(hail)%q(i,j,k) + dref(i,j,k)*cint*sbqh
    	    diag(2)%micro(hail)%q(i,j,k) = diag(2)%micro(hail)%q(i,j,k) - dref(i,j,k)*cint*xmhr
    	    diag(3)%micro(hail)%q(i,j,k) = diag(3)%micro(hail)%q(i,j,k) + dref(i,j,k)*cint*shegh*(shelg*(clg + clgg) + sherg*(crg + crgg))
    	    diag(4)%micro(hail)%q(i,j,k) = diag(4)%micro(hail)%q(i,j,k) - dref(i,j,k)*cint*shegh*cghs
    	    diag(5)%micro(hail)%q(i,j,k) = diag(5)%micro(hail)%q(i,j,k) + dref(i,j,k)*cint*(clh + crh)
    	    diag(6)%micro(hail)%q(i,j,k) = diag(6)%micro(hail)%q(i,j,k) + dref(i,j,k)*cint*(cih + csh)
    	    diag(7)%micro(hail)%q(i,j,k) = diag(7)%micro(hail)%q(i,j,k) - dref(i,j,k)*cint*shelh*hfs
    	    diag(8)%micro(hail)%q(i,j,k) = diag(8)%micro(hail)%q(i,j,k) + cint*Tdg
	  endif
	  
	  if ( moments == 2 ) then
	    diag(1)%micro(drop)%n(i,j,k) = diag(1)%micro(drop)%n(i,j,k) + dref(i,j,k)*cint*evnc
	    diag(4)%micro(drop)%n(i,j,k) = diag(4)%micro(drop)%n(i,j,k) - dref(i,j,k)*cint*clin
	    diag(5)%micro(drop)%n(i,j,k) = diag(5)%micro(drop)%n(i,j,k) - dref(i,j,k)*cint*clgn
	    diag(6)%micro(drop)%n(i,j,k) = diag(6)%micro(drop)%n(i,j,k) - dref(i,j,k)*cint*clsn
!
	    diag(1)%micro(rain)%n(i,j,k) = diag(1)%micro(rain)%n(i,j,k) + dref(i,j,k)*cint*evnr
	    diag(4)%micro(rain)%n(i,j,k) = diag(4)%micro(rain)%n(i,j,k) + dref(i,j,k)*cint*(xmirn + xmgrn + xmsrn)
	    diag(5)%micro(rain)%n(i,j,k) = diag(5)%micro(rain)%n(i,j,k) + dref(i,j,k)*cint*((1.-sheli)*clin - sheri*crin)
	    diag(6)%micro(rain)%n(i,j,k) = diag(6)%micro(rain)%n(i,j,k) + dref(i,j,k)*cint*((1.-shelg)*clgn - sherg*crgn)
	    diag(7)%micro(rain)%n(i,j,k) = diag(7)%micro(rain)%n(i,j,k) + dref(i,j,k)*cint*((1.-shels)*clsn - shers*crsn)
!
	    diag(1)%micro(ice)%n(i,j,k) = diag(1)%micro(ice)%n(i,j,k) + dref(i,j,k)*cint*sbni
	    diag(2)%micro(ice)%n(i,j,k) = diag(2)%micro(ice)%n(i,j,k) - dref(i,j,k)*cint*xmirn
	    diag(3)%micro(ice)%n(i,j,k) = diag(3)%micro(ice)%n(i,j,k) - 2.*dref(i,j,k)*cint*aisn
	    diag(4)%micro(ice)%n(i,j,k) = diag(4)%micro(ice)%n(i,j,k) + dref(i,j,k)*cint*(sheli*(1.-frimi) - 1.)*clin
	    diag(5)%micro(ice)%n(i,j,k) = diag(5)%micro(ice)%n(i,j,k) + dref(i,j,k)*cint*(sheri*(1.-frimr) - 1.)*crin
	    diag(6)%micro(ice)%n(i,j,k) = diag(6)%micro(ice)%n(i,j,k) - dref(i,j,k)*cint*cisn
	    diag(7)%micro(ice)%n(i,j,k) = diag(7)%micro(ice)%n(i,j,k) + dref(i,j,k)*cint*hmn
!	  
	    if ( lmicro > 2 ) then
	      diag(1)%micro(grau)%n(i,j,k) = diag(1)%micro(grau)%n(i,j,k) + dref(i,j,k)*cint*sbng
	      diag(2)%micro(grau)%n(i,j,k) = diag(2)%micro(grau)%n(i,j,k) - dref(i,j,k)*cint*xmgrn
	      diag(3)%micro(grau)%n(i,j,k) = diag(3)%micro(grau)%n(i,j,k) + dref(i,j,k)*cint*(sheli*frimi*clin + shels*frims*clsn)
	      diag(4)%micro(grau)%n(i,j,k) = diag(4)%micro(grau)%n(i,j,k) + dref(i,j,k)*cint*(sheri*frimr*crin + shers*frimr*crsn)
!	      
	      diag(1)%micro(snow)%n(i,j,k) = diag(1)%micro(snow)%n(i,j,k) + dref(i,j,k)*cint*sbns
	      diag(2)%micro(snow)%n(i,j,k) = diag(2)%micro(snow)%n(i,j,k) - dref(i,j,k)*cint*xmsrn
	      diag(3)%micro(snow)%n(i,j,k) = diag(3)%micro(snow)%n(i,j,k) + dref(i,j,k)*cint*aisn       
	      diag(4)%micro(snow)%n(i,j,k) = diag(4)%micro(snow)%n(i,j,k) + dref(i,j,k)*cint*(shels*(1.-frims) - 1.)*clsn
	      diag(5)%micro(snow)%n(i,j,k) = diag(4)%micro(snow)%n(i,j,k) + dref(i,j,k)*cint*(shers*(1.-frimr) - 1.)*crsn
	      diag(7)%micro(snow)%n(i,j,k) = diag(5)%micro(snow)%n(i,j,k) - dref(i,j,k)*cint*csgn
	      diag(8)%micro(snow)%n(i,j,k) = diag(7)%micro(snow)%n(i,j,k) - dref(i,j,k)*cint*cssn
	    endif
!
	    if ( lmicro > 3 ) then
	      diag(1)%micro(hail)%n(i,j,k) = diag(1)%micro(hail)%n(i,j,k) + dref(i,j,k)*cint*sbnh
	      diag(2)%micro(hail)%n(i,j,k) = diag(2)%micro(hail)%n(i,j,k) - dref(i,j,k)*cint*xmhrn
	      diag(3)%micro(hail)%n(i,j,k) = diag(3)%micro(hail)%n(i,j,k) + dref(i,j,k)*cint*shegh*(shelg*clgn + sherg*crgn)
	      diag(7)%micro(hail)%n(i,j,k) = diag(7)%micro(hail)%n(i,j,k) + dref(i,j,k)*cint*(shegh*cghs + shelh*hfs)
	    endif
	  endif
	endif
!
!----------------------------------------------------------!
!          Regeneration of aerosols and chemicals          !
!----------------------------------------------------------!	
!
!  particles release during evaporation
!
        if ( lndrop == 1 .and. aero_flg%reg ) then
!
!  Ice nuclei
!
#ifdef NUC_CNT
	  call in_evap ( i, j, k,  hydromtrl, nucin2, Tem-273.15, evnc, evnr, sbni, sbng, sbns, sbnh, xmirn, xmgrn, xmsrn, xmhrn )
#endif
!
!  Aerosols
!
#ifdef AERO_ENABLE
  	  call regeneration ( i, j, k, dref(i,j,k), sum(hydromtrl%n(i,j,k)), aero3dl, evnc+evnr+sbni+sbng+sbns, aero3dm )
!
!  Chemicals
!	     
#if (defined AQCHEM_ENABLE) || (defined SOLIDCHEM_ENABLE)
          call evscaveng ( i, j, k, hydromtrl, evqc, evqr, evqi )
#endif
!
!  Liquid phase chemicals
!
#ifdef AQCHEM_ENABLE
          call aqscaveng ( i, j, k, hydromtrl, auq, clr, cli, clg, cirr, crg, csrr, cls )
#endif
!
!  Solid phase chemicals
!
#ifdef SOLIDCHEM_ENABLE
          call solscaveng ( i, j, k, hydromtrl, cli, aig, ais, cir, cis )
#endif
#endif
!
	endif
!
!----------------------------------------------------------!
!                     End i j k loops                      !
!----------------------------------------------------------!
!
        endif
      end do
    end do
  end do
!
  if (ldtmic) then
#ifdef SPMD
    call MPI_ALLReduce (x, dtmic, 1, REALTYPE, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
    dtmic = x
#endif 
  endif
!
!  Homogenize evaporation
!
  if (no_cp) then
    call horav ( evapr, evm )
!
    do k = 1, nz
      do j = jt_start, jt_end
        do i = it_start, it_end
	  Tem = thermo%T(i,j,k)
	  exn = thermo%exn(i,j,k)
          statem%es(i,j,k) = statem%es(i,j,k) - cal_flv(Tem)*(evapr(i,j,k) - evm(k)) / (cp_a*exn)
          statem%qt(i,j,k) = statem%qt(i,j,k) + (evapr(i,j,k) - evm(k))
        enddo
      enddo
    enddo
  endif
!
  if (verbose > 1) call write_debug('Terminating sb_ice')
!
! ------------------------------------------------------------------
!	
return
end
!
end module seifert_beheng_ice
