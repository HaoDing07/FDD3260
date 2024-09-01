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
!  MICROPHYSICS1.F                   
!
!  Purpose:
!	Subroutines for calculating the first stage of microphysics 
!       processes
!
!  Author
!	Julien Savre, MISU
!
! ================================================================
!
  module seifert_beheng_warm
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

  public :: sb_warm

  contains

!	==========================================
	subroutine sb_warm ( dref, ww, dTdt, qt, pressurel, hydromtrl, aero3dl, statem, hydromtrm, aero3dm ) 
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
  integer :: i, j, k, h
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: evapr
  real :: xfkd, xfmu, aunc, aunsc, aunr, auq, cccn, clr, clrn,  &
	     crrn, cd, cdr, dqc, dqr, evqc, evqr, evnc,     &
	     evnr, dnc, dnr, dqi, dqg, dqs, dqh, pp, Tem,   &
	     exn, rho, esw, qsw, sw, qc, nc, qr, nr, des,   &
	     dp, q, qv, ql, qi, cp, cv, gamma, evm(nz), x
! 
#if ( defined SPMD )      
  integer :: ierr
#endif
!
  if (verbose > 1) call write_debug('Starting sb_warm')
!
!----------------------------------------------------------!
!                     Starting loops                       !
!----------------------------------------------------------!
!
  x = 1.e9
  evapr = 0.0
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
	if ( (hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. hydromtrl(drop)%n(i,j,k) > xnmin(drop))       &
        .or. (hydromtrl(rain)%q(i,j,k) > qmin(rain) .and. hydromtrl(rain)%n(i,j,k) > xnmin(rain)) ) then   
!
!----------------------------------------------------------!
!                 Define thermo properties                 !
!----------------------------------------------------------!
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
        ql = 0.
        qi = 0.
	if (lmicro > 0) ql = hydromtrl(drop)%q(i,j,k) + hydromtrl(rain)%q(i,j,k)
	if (lmicro > 1) qi = hydromtrl(ice)%q(i,j,k)
	if (lmicro > 2) qi = qi + hydromtrl(grau)%q(i,j,k) + hydromtrl(snow)%q(i,j,k)
	if (lmicro > 3) qi = qi + hydromtrl(hail)%q(i,j,k)
        qv = qv - ql - qi
!
	Tem = thermo%T(i,j,k)
	exn = thermo%exn(i,j,k)
	esw = cal_esw (Tem)
	qsw = cal_qsw (Tem, pp)
	sw  = qv - qsw
!
        call get_cp_cv (qv, ql, qi, cp, cv)
!
        xfkd = cal_xfkd (Tem, pp)
	xfmu = cal_xfmu (Tem)
	gamma = 1./qv + cal_flv(Tem)**2./((cp_v-cv_v)*cp*Tem*Tem) 
!
!----------------------------------------------------------!
!       	  Initialize Microphysics	           !
!----------------------------------------------------------!
!
	auq    = 0.0
	aunc   = 0.0
	aunsc  = 0.0
	cccn   = 0.0
	aunr   = 0.0
	crrn   = 0.0
	clr    = 0.0
	clrn   = 0.0
	cd     = 0.0
	cdr    = 0.0
	evqc   = 0.0
	evqr   = 0.0
	evnc   = 0.0
	evnr   = 0.0
	dqc    = 0.0
	dqr    = 0.0
	des    = 0.0
	dp     = 0.0
	dnc    = 0.0
	dnr    = 0.0
!
	qc = hydromtrl(drop)%q(i,j,k)
	qr = hydromtrl(rain)%q(i,j,k)
	nc = hydromtrl(drop)%n(i,j,k)
	nr = hydromtrl(rain)%n(i,j,k)
!
!----------------------------------------------------------!
!             Warm: evaporation/condensation               !
!----------------------------------------------------------!
!
      if ( lmicro == 1 .or. lfreeze == 0 ) then
!
      if (with_cond) then
!
#ifndef SAT_ADJ
	if( qc > qmin(drop) .and. nc > xnmin(drop) )then
	  call sbcond (k, drop, qsw, esw, Tem, xfkd, xfmu, qc, nc, cd)
	endif
#endif
! 
	if( qr > qmin(rain) .and. nr > xnmin(rain) ) then
	  call sbcond (k, rain, qsw, esw, Tem, xfkd, xfmu, qr, nr, cdr)
	endif
!
!  Calculate tendencies with resolved saturation
!
      call resolved_saturation ( k, cd, cdr, 0., 0., 0., 0., 0., pp, Tem, ww(i,j,k), dTdt(i,j,k), cp, cv, &
                exn, sw, qsw, 0., 0., qv, qc, qr, 0., 0., 0., 0., evqc, evqr, dqi, dqg, dqs, dqh, des )
!
!  Evaporation
!
      if ( moments == 2 ) then
#ifndef SAT_ADJ
	if(evqc < 0.0 .and. lndrop == 1)then  
	  call sbevap (k, drop, qc, nc, xfkd, xfmu, evqc, evnc)
	endif
#endif
!
	if (evqr < 0.0) then
          evapr(i,j,k) = evqr
	  call sbevap (k, rain, qr, nr, xfkd, xfmu, evqr, evnr)
	endif
      endif
!
!  Diagnostics
!
      if (out_diagt) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + dref(i,j,k)*cint*des
      if (out_diagl) diag(4)%qc(i,j,k) = diag(4)%qc(i,j,k) + dref(i,j,k)*cint*evqc
      if (out_diagr) diag(4)%qr(i,j,k) = diag(4)%qr(i,j,k) + dref(i,j,k)*cint*evqr
!
      if (ldtmic .and. (cd + cdr) > 1.e-12) x = max(min(x,1./(gamma*(cd+cdr))),1.e-12)
!
      endif
      endif
!
!----------------------------------------------------------!
!       	 Warm Microphysical processes	           !
!                   Seifert & Beheng model                 !
!----------------------------------------------------------!
!
!  Autoconversion from cloud to rain: c + c --> r
!
	if( qc > qmin(drop) .and. nc > xnmin(drop) )then
          if ( auto == 1 ) then
	    call kess ( rho, qc, nc, auq, aunc, aunr )
          else if ( auto == 2 ) then
	    call liu_daum ( k, rho, qc, nc, nr, auq, aunc, aunr )
	  else
	    if ( ldrizz ) call aucr ( k, rho, qc, qr, nc, nr, auq, aunc, aunr )
	  endif
	endif
!
!  Self-collection of cloud drops without autoconversion: c + c --> c
!  Note: aunsc includes aunc, no need to count autoconversion twice
!
	if( lndrop == 1 .and. qc > qmin(drop) .and. nc > xnmin(drop) )then
          call scc ( k, rho, qc, nc, aunsc )
        endif 
!
!  Collection of cloud drop by rain: c + r --> r
!
	if ( qc > qmin(drop) .and. nc > xnmin(drop) .and. qr > qmin(rain) .and. nr > xnmin(rain) ) then
	  call ccr ( k, rho, qc, qr, nc, nr, clr, clrn )
	endif
!
!  Self-collection, breakup, & spontaneous breakup 
!  of rain drop: r + r --> r
!
	if( moments == 2 .and. qr > qmin(rain) .and. nr > xnmin(rain) )then
	  call scr ( k, rho, qr, nr, crrn )
	endif
!
!----------------------------------------------------------!
!             Update microphysics source terms             !
!----------------------------------------------------------!
!
!  Total tendencies    
!
	dqc = evqc - min(auq + clr, qc/dt0)
	dqr = evqr + min(auq + clr, qc/dt0)
	if ( moments == 2 ) then
          if ( lndrop == 1 ) dnc = evnc - min(aunsc + clrn, nc/dt0)
	  dnr = evnr + min(aunr - crrn, nc/dt0)
	endif
!
	hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) + dqc
	hydromtrm(rain)%q(i,j,k) = hydromtrm(rain)%q(i,j,k) + dqr 
	if ( moments == 2 ) then
          if (lndrop == 1) hydromtrm(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) + dnc
	  hydromtrm(rain)%n(i,j,k) = hydromtrm(rain)%n(i,j,k) + dnr 
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
	if ( ldiag .and. out_lwp ) surf%cond(i,j) = surf%cond(i,j) + dref(i,j,k)*cint*(max(evqc,0.) + max(evqr,0.))*dzw(k)
	if ( ldiag .and. out_lwp ) surf%evap(i,j) = surf%evap(i,j) + dref(i,j,k)*cint*(min(evqc,0.) + min(evqr,0.))*dzw(k)
!
        if ( ldiag .and. out_micro ) then
	  diag(1)%micro(drop)%q(i,j,k) = diag(1)%micro(drop)%q(i,j,k) + dref(i,j,k)*cint*evqc
	  diag(2)%micro(drop)%q(i,j,k) = diag(2)%micro(drop)%q(i,j,k) - dref(i,j,k)*cint*auq
	  diag(3)%micro(drop)%q(i,j,k) = diag(3)%micro(drop)%q(i,j,k) - dref(i,j,k)*cint*clr
!
          diag(1)%micro(rain)%q(i,j,k) = diag(1)%micro(rain)%q(i,j,k) + dref(i,j,k)*cint*evqr
	  diag(2)%micro(rain)%q(i,j,k) = diag(2)%micro(rain)%q(i,j,k) + dref(i,j,k)*cint*auq
	  diag(3)%micro(rain)%q(i,j,k) = diag(3)%micro(rain)%q(i,j,k) + dref(i,j,k)*cint*clr
!
	  if (moments == 2) then
	    diag(1)%micro(drop)%n(i,j,k) = diag(1)%micro(drop)%n(i,j,k) + dref(i,j,k)*cint*evnc
	    diag(2)%micro(drop)%n(i,j,k) = diag(2)%micro(drop)%n(i,j,k) - dref(i,j,k)*cint*aunsc
	    diag(3)%micro(drop)%n(i,j,k) = diag(3)%micro(drop)%n(i,j,k) - dref(i,j,k)*cint*clrn
!
            diag(1)%micro(rain)%n(i,j,k) = diag(1)%micro(rain)%n(i,j,k) + dref(i,j,k)*cint*evnr
	    diag(2)%micro(rain)%n(i,j,k) = diag(2)%micro(rain)%n(i,j,k) + dref(i,j,k)*cint*aunr
	    diag(3)%micro(rain)%n(i,j,k) = diag(3)%micro(rain)%n(i,j,k) - dref(i,j,k)*cint*crrn	 
	  endif
	endif
!
!----------------------------------------------------------!
!          Regeneration of aerosols and chemicals          !
!----------------------------------------------------------!	
!
        if ( lmicro == 1 .and. aero_flg%reg) then
!
!  Ice nuclei
!
#ifdef NUC_CNT
	  call in_evap ( i, j, k, hydromtrl, nucin2, Tem-273.15, dnc, dnr, 0., 0., 0., 0., 0., 0., 0., 0. )
#endif
!
!  Aerosols
!
#ifdef AERO_ENABLE
	  call regeneration ( i, j, k, dref(i,j,k), sum(hydromtrl%n(i,j,k)), aero3dl, dnc+dnr, aero3dm )
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
  if (verbose > 1) call write_debug('Terminating sb_warm')
!
!----------------------------------------------------------!
!
return
end

end module seifert_beheng_warm
