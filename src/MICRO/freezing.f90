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
!  FREEZING:
!	Freezing routines
!
!  Author:
!	Julien Savre, LMU Munich
!
! ==============================================================
!
  module freezing
!
USE gridno
USE shared_all
USE thermodynamics
!
  IMPLICIT NONE

  private
  
  public :: freez, ice_nuc, foch, focdw

  contains
!
!	==========================================
	subroutine freez ( dref, ww, pressurel, statel, hydromtrl, statem, hydromtrm ) 
!	==========================================
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dref, ww
  type(atm_state) :: statel, statem
  type(atm_pressure) :: pressurel
  type(hydrometeor), dimension(1:nhydro) :: hydromtrl, hydromtrm
!
  integer :: i, j, k, h, l
  real    :: Tem, rho, pp, qv, ql, qi, cp, cv, exn, ssw, ssi,    &
	     e, Td, fci, fcni, fdi, fdni, fcd, fcnd, fch, fcnh,  &
	     frh, frnh, frg, frgn, dqcdt, dqrdt, des
!
  if (verbose > 0) call write_debug('Start freez')
!
!----------------------------------------------------------!
!                     Starting loop                        !
!----------------------------------------------------------!
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
!
	if ( ((hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. hydromtrl(drop)%n(i,j,k) > xnmin(drop))       &
         .or. (hydromtrl(rain)%q(i,j,k) > qmin(rain) .and. hydromtrl(rain)%n(i,j,k) > xnmin(rain)))      & 
         .and. thermo%T(i,j,k) < 273.15 ) then   
!
!----------------------------------------------------------!
!                  Local initializations                   !
!----------------------------------------------------------!
!
	fch  = 0.0
	fcnh = 0.0
	frh  = 0.0
	frnh = 0.0
	fdi  = 0.0
	fdni = 0.0
	fci  = 0.0
	fcni = 0.0
	fcd  = 0.0
	fcnd = 0.0
	frg  = 0.0
	frgn = 0.0
	des  = 0.0
!
#ifdef ANELASTIC
	pp = p0(k)
        rho = den0(k)
#else
        pp = pressurel%p(i,j,k)+p0(k)
        rho = pressurel%dens(i,j,k)
#endif
!
	qv  = statel%qt(i,j,k)
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
	ssw = qv/cal_qsw(Tem, pp)
	ssi = qv/cal_qsi(Tem, pp)
        e   = pp * qv/0.622 / (1. + qv/0.622)
	Td  = Tem 
!
        call get_cp_cv (qv, ql, qi, cp, cv)
!
        dqcdt = (hydromtrl(drop)%q(i,j,k) - hydromtr2(drop)%q(i,j,k)) / dt0
        dqrdt = (hydromtrl(rain)%q(i,j,k) - hydromtr2(rain)%q(i,j,k)) / dt0
!
!----------------------------------------------------------!
!             Advanced freezing parameterization           !
!----------------------------------------------------------!
!
#if (defined NUC_CNT)
!
!  Heterogeneous freezing from CNT model
!
        call freez_cnt ( i, j, k, nucin2, hydromtrl(drop)%q(i,j,k), hydromtrl(drop)%n(i,j,k), hydromtrl(ice)%n(i,j,k),    &
	                 Tem, ssw, ssi, e, fci, fcni, fcd, fcnd )    
!
!  Homogeneous freezing of cloud droplets
!
	if ( Tem < 238. .and. hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. hydromtrl(drop)%n(i,j,k) > xnmin(drop) ) then
          call foch ( i, j, k, drop, hydromtrl(drop)%q(i,j,k), hydromtrl(drop)%n(i,j,k), hydromtrl(ice)%n(i,j,k),       &
	              Tem, rad%dtnet(i,j,k), ww(i,j,k), ssw, ssi, fch, fcnh )
	endif
!
!  Homogeneous freezing of rain drops: deactivated
!
        if (lmicro > 2) then
	  if ( Tem < 238. .and. hydromtrl(rain)%q(i,j,k) > qmin(rain) .and. hydromtrl(rain)%n(i,j,k) > xnmin(rain) ) then	    
            call foch ( i, j, k, rain, hydromtrl(rain)%q(i,j,k), hydromtrl(rain)%n(i,j,k), hydromtrl(grau)%n(i,j,k)    & 
	                Tem, rad%dtnet(i,j,k), ww(i,j,k), ssw, ssi, frh, frnh )
    	  endif
	endif
!
!  Updates
!
        hydromtrm(ice)%q(i,j,k)  = hydromtrm(ice)%q(i,j,k)  + min((fci + fcd + fch),hydromtrl(drop)%q(i,j,k)/dt0)
        hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) - min((fci + fcd + fch),hydromtrl(drop)%q(i,j,k)/dt0)
        hydromtrm(ice)%n(i,j,k)  = hydromtrm(ice)%n(i,j,k)  + min((fcni + fcnd + fcnh),hydromtrl(drop)%n(i,j,k)/dt0)
        if (lndrop == 1) hydromtrm(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) - min((fcni + fcnd + fcnh),hydromtrl(drop)%n(i,j,k)/dt0)
!
        if (lmicro > 2) then
	  hydromtrm(grau)%q(i,j,k) = hydromtrm(grau)%q(i,j,k) + min(frh,hydromtrl(rain)%q(i,j,k)/dt0) 
          hydromtrm(rain)%q(i,j,k) = hydromtrm(rain)%q(i,j,k) - min(frh,hydromtrl(rain)%q(i,j,k)/dt0)
          hydromtrm(grau)%n(i,j,k) = hydromtrm(grau)%n(i,j,k) + min(frnh,hydromtrl(rain)%n(i,j,k)/dt0)
          hydromtrm(rain)%n(i,j,k) = hydromtrm(rain)%n(i,j,k) - min(frnh,hydromtrl(rain)%n(i,j,k)/dt0)
	endif
!
#ifdef ISENTROPIC
	des = cal_flm(Tem)*(min(fci + fcd + fch,hydromtrl(drop)%q(i,j,k)/dt0) + min(frh,hydromtrl(rain)%q(i,j,k)/dt0)) / (cp*exn)
        statem%es(i,j,k) = statem%es(i,j,k) + des 
#endif
!
!----------------------------------------------------------!
!                     Simple freezing                      !
!----------------------------------------------------------!
!
#else
!
!  Deposition nucleation
!
          if ( ssi > 1.02 .and. Tem < 268.15 ) then
	    if ( lfreeze > 1 ) then
              call cooper ( Tem, ssi, hydromtrl(ice)%n(i,j,k), fdi, fdni )
            endif
          endif
!
!  Droplet freezing
!
!
	  if ( hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. hydromtrl(drop)%n(i,j,k) > xnmin(drop) ) then
            if ( lfreeze == 1 ) then 
	      call ice_nuc ( Tem, ssi, hydromtrl(drop)%q(i,j,k), hydromtrl(drop)%n(i,j,k), hydromtrl(ice)%n(i,j,k), nuc%in(i,j,k), fci, fcni )
            else if ( lfreeze == 2 ) then
              call bigg ( drop, Tem, hydromtrl(drop)%q(i,j,k), hydromtrl(drop)%n(i,j,k), hydromtrl(ice)%n(i,j,k), fci, fcni)
	    else if ( lfreeze == 3 ) then
              call focdw ( drop, hydromtrl(drop)%q(i,j,k), hydromtrl(drop)%n(i,j,k), hydromtrl(ice)%q(i,j,k), hydromtrl(ice)%n(i,j,k), 	&
	  		nuc%in(i,j,k)/nuc%ccn(i,j,k), Tem, rad%dtnet(i,j,k), ww(i,j,k), dqcdt, ssw, fci, fcni )
    	    endif
    	  endif
!
!  Rain freezing
!
          if ( hydromtrl(rain)%q(i,j,k) > qmin(rain) .and. hydromtrl(rain)%n(i,j,k) > xnmin(rain) .and. lmicro > 2 ) then
            if ( lfreeze == 2 ) then
              call bigg ( rain, Tem, hydromtrl(rain)%q(i,j,k), hydromtrl(rain)%n(i,j,k), hydromtrl(grau)%n(i,j,k), frg, frgn)
            else if ( lfreeze == 3 ) then	    
              call focdw ( rain, hydromtrl(rain)%q(i,j,k), hydromtrl(rain)%n(i,j,k), hydromtrl(grau)%q(i,j,k), hydromtrl(grau)%n(i,j,k), 	&
	  	  	nuc%in(i,j,k)/nuc%ccn(i,j,k), Tem, rad%dtnet(i,j,k), ww(i,j,k), dqrdt, ssw, frg, frgn )
    	    endif
	  endif    
!	
!  Homogeneous freezing of cloud droplets
!
	  if ( Tem < 235.15 .and. hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. hydromtrl(drop)%n(i,j,k) > xnmin(drop) ) then
            call foch ( i, j, k, drop, hydromtrl(drop)%q(i,j,k), hydromtrl(drop)%n(i,j,k), hydromtrl(ice)%n(i,j,k),     &
	                Tem, rad%dtnet(i,j,k), ww(i,j,k), ssw, ssi, fch, fcnh )
	  endif
!
	  if ( Tem < 235.15 .and. hydromtrl(rain)%q(i,j,k) > qmin(rain) .and. hydromtrl(rain)%n(i,j,k) > xnmin(rain) ) then
            call foch ( i, j, k, rain, hydromtrl(rain)%q(i,j,k), hydromtrl(rain)%n(i,j,k), hydromtrl(grau)%n(i,j,k),     &
	                Tem, rad%dtnet(i,j,k), ww(i,j,k), ssw, ssi, frh, frnh )
	  endif
!
!  freezing
!
          hydromtrm(ice)%q(i,j,k)  = hydromtrm(ice)%q(i,j,k) + fdi +  min(fci + fch,hydromtrl(drop)%q(i,j,k)/dt0)
          hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) - min(fci + fch,hydromtrl(drop)%q(i,j,k)/dt0)
          hydromtrm(ice)%n(i,j,k)  = hydromtrm(ice)%n(i,j,k) + fdni + min(fcni + fcnh,hydromtrl(drop)%n(i,j,k)/dt0)
          if (lndrop == 1) hydromtrl(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) - min(fcni + fcnh,hydromtrl(drop)%n(i,j,k)/dt0)
!
          if (lmicro > 2) then
	    hydromtrm(grau)%q(i,j,k) = hydromtrm(grau)%q(i,j,k) + min(frg + frh,hydromtrl(rain)%q(i,j,k)/dt0) 
            hydromtrm(rain)%q(i,j,k) = hydromtrm(rain)%q(i,j,k) - min(frg + frh,hydromtrl(rain)%q(i,j,k)/dt0)
            hydromtrm(grau)%n(i,j,k) = hydromtrm(grau)%n(i,j,k) + min(frgn + frnh,hydromtrl(rain)%n(i,j,k)/dt0)
            hydromtrm(rain)%n(i,j,k) = hydromtrm(rain)%n(i,j,k) - min(frgn + frnh,hydromtrl(rain)%n(i,j,k)/dt0)
	  endif
!
#ifdef ISENTROPIC
	  des = (cal_fls(Tem)*fdi + cal_flm(Tem)*(min(fci + fch,hydromtrl(drop)%q(i,j,k)/dt0) + min(frg + frh,hydromtrl(rain)%q(i,j,k)/dt0))) / (cp*exn)
          statem%es(i,j,k) = statem%es(i,j,k) + des 
#endif
!
#endif
!
!----------------------------------------------------------!
!  		    Micro diagnostics			   !
!----------------------------------------------------------!
!	
        if (ldiag) then
          if (out_diagt) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + dref(i,j,k)*cint*des
          if (out_diagl) diag(8)%qc(i,j,k) = diag(8)%qc(i,j,k) - dref(i,j,k)*cint*(fci + fch)
          if (out_diagi) diag(8)%qi(i,j,k) = diag(8)%qi(i,j,k) + dref(i,j,k)*cint*(fdi + fci + fch)
          if (out_diagr) diag(8)%qr(i,j,k) = diag(8)%qr(i,j,k) - dref(i,j,k)*cint*(frg + frh)
!
          if (out_diagt) diag(9)%pt(i,j,k) = diag(9)%pt(i,j,k) + dref(i,j,k)*cint*des
          if (out_diagl) diag(9)%qc(i,j,k) = diag(9)%qc(i,j,k) - dref(i,j,k)*cint*(fci + fch)
          if (out_diagi) diag(9)%qi(i,j,k) = diag(9)%qi(i,j,k) + dref(i,j,k)*cint*(fdi + fci + fch)
          if (out_diagr) diag(9)%qr(i,j,k) = diag(9)%qr(i,j,k) - dref(i,j,k)*cint*(frg + frh)
!
          if (out_micro) then
    	    diag(8)%micro(drop)%q(i,j,k) = diag(8)%micro(drop)%q(i,j,k) - dref(i,j,k)*cint*(fci + fch)
	    if (lndrop == 1) diag(8)%micro(drop)%n(i,j,k) = diag(8)%micro(drop)%n(i,j,k) - dref(i,j,k)*cint*(fcni + fcnh)
    	    diag(8)%micro(ice)%q(i,j,k) = diag(8)%micro(ice)%q(i,j,k) + dref(i,j,k)*cint*(fdi + fci + fch)
	    diag(8)%micro(ice)%n(i,j,k) = diag(8)%micro(ice)%n(i,j,k) + dref(i,j,k)*cint*(fdni + fcni + fcnh)
!
	    if (lmicro > 2) then
       	      diag(8)%micro(rain)%q(i,j,k) = diag(8)%micro(rain)%q(i,j,k) - dref(i,j,k)*cint*(frg + frh)
	      diag(8)%micro(rain)%n(i,j,k) = diag(8)%micro(rain)%n(i,j,k) - dref(i,j,k)*cint*(frgn + frnh)
    	      diag(8)%micro(grau)%q(i,j,k) = diag(8)%micro(grau)%q(i,j,k) + dref(i,j,k)*cint*(frg + frh)
	      diag(8)%micro(grau)%n(i,j,k) = diag(8)%micro(grau)%n(i,j,k) + dref(i,j,k)*cint*(frgn + frnh)
            endif
	  endif
	endif
!
!----------------------------------------------------------!
!                     End i j k loops                      !
!----------------------------------------------------------! 
!
        endif
!
      end do
    end do
  end do
!
  if (verbose > 0) call write_debug('Terminate freez')
!
!----------------------------------------------------------!
!	
return
end
!
!  ==================================================
   subroutine ice_nuc ( t, si, qc, nc, ni, xin, fr, frn )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Calculation of ice nucleation
!
! ================================================================
!
  real :: t, si, qc, nc, ni, xin, fr, frn
  real :: dt1, xnun, xmini
!
!----------------------------------------------------------!
!                        Starting                          !
!----------------------------------------------------------!
!
  dt1  = 1./dt0
!
  xnun = 0.
  xmini = hydrop(ice)%xmin !(dmini/hydrop(ice)%cm)**(1./hydrop(ice)%am)
!
!  Immersion
!
  select case (trim(casename))
    case('MOCCHA')
      if ( qc >= 1.e-6 .and. si > 1.05 ) then
        xnun = max(xin - ni, 0.)*dt1
      endif
!
    case ('MIDLAT')
      if ( (qc >= 1.e-7 .and. T <= 267.15) .or. (si > 1.02 .and. T <= 248.15) ) then
        xnun = 100000. 
      endif
!
    case default
      if ( qc >= 1.e-7 .and. T <= 268.15 .and. si > 1.05 ) then
        xnun = max(xin - ni, 0.)*dt1
      endif
!
  end select
!
!  Update tendencies
!
  frn = min(xnun, (nc - ni)*dt1)
  fr  = min(frn*xmini, qc*dt1)
!
!----------------------------------------------------------!
!
return
end
!
!  ==================================================
   subroutine cooper ( T, si, ni, fr, frn )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Cooper simple immersion freezing
!
! ================================================================
!
  real :: T, si, ni, xin, fr, frn
  real :: dt1, xnun, xmini
!
!----------------------------------------------------------!
!                        Starting                          !
!----------------------------------------------------------!
!
  dt1  = 1./dt0
  xmini = hydrop(ice)%xmin
!
  frn = 0.
  fr = 0.
!
  if ( T <= 268.15 .and. si > 1.02 ) then
    !xin = 117.*exp(0.125*max(T00 - T,35.))
    xin = 5.*exp(-0.304*max(T - T00,-35.))
!
    xnun = max(xin - ni, 0.)*dt1
!
!  Update tendencies
!
    frn = xnun
    fr  = frn*xmini
  endif
!
!----------------------------------------------------------!
!
return
end
!
!  ==================================================
   subroutine bigg ( h, T, qc, nc, ni, fr, frn )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!       Deposition nucleation from Bigg (1953)
!
! ================================================================
!
  integer :: h
  real :: T, mc, qc, nc, ni, xin, fr, frn
  real :: dt1, xnun, xmini, jb, ff
  real :: a = 0.2, b = 0.65
!
!----------------------------------------------------------!
!                        Starting                          !
!----------------------------------------------------------!
!
  dt1 = 1./dt0
  mc = cal_avm(h, qc, nc)
!
  frn = 0.
  fr = 0.
!
  if ( T <= 268.15 .and. qc > qmin(h) .and. nc > xnmin(h) ) then
    Jb = a*exp(-b*max(T - T00,-35.) - 1.)
    xnun = dt1*nc*(1. - exp(-mc*Jb*dt0))
!
    frn = max(min(xnun, (nc - ni)/dt0), 0.)
    fr = hydroc(h)%cr6*mc*frn
  endif
!
!----------------------------------------------------------!
!
return
end
!
!  ==================================================
   subroutine foch ( i, j, k, h, qc, nc, ni, Tem, dtrad, ww, sw, si, fci, fcin )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Homogeneous freezing of cloud droplets (constants in shared_data):
!         - Follows mostly Barahona, ACP 2014
!
! ================================================================
!
  integer  :: h, i, j, k, l
  real     :: nc, qc, ni, mc, fci, fcin
  real     :: sw, si, dTdt, Tem, dtrad, ww
  real     :: dt1, a, d, vc, dgact, dfcr, J0, aw, xnun
!
!----------------------------------------------------------!
!            Calculate critical free energy                !
!            of homogeneous germ formation		   !
!----------------------------------------------------------!
!
  dt1 = 1./dt0
  mc  = cal_avm( h, qc, nc )
!
  aw    = cal_aw( 1, Tem, sw, qc, nc )
  dgact = cal_dgact( Tem )
  dTdt = dtrad - g*ww/cp_a
  dfcr  = 16./3.*pi*cal_sigmaiw(Tem)**3.*(mwm/(rhoi*kb*Tem*log(aw*si/sw)))**2.
!
!----------------------------------------------------------!
!         Calculate freezing and nucleation rates          !
!----------------------------------------------------------!
!
!  Logarithmic calculation to limit floating point errors
!
  fcin = 0.
  fci = 0.
!
  if ( Tem <= 268.15 .and. qc > qmin(h) .and. nc > xnmin(h) ) then
    a   = log(kbhp) + log(Na) + log(Tem/Mw)
    J0  = -(dgact + dfcr) / (kb*Tem)
    J0  = exp(a + J0) 
    xnun = dt1*nc*(1. - exp(-mc*J0*dt0))
!
    fcin = max(min(xnun, (nc - ni)/dt0), 0.)
    fci = hydroc(h)%cr6*mc*fcin
  endif
!
!----------------------------------------------------------!
!         Frozen INs are added to frozen particles         !
!----------------------------------------------------------!
!
#ifdef NUC_CNT
  do l = 1, 3
      nucin2(l)%mode(3)%n(i,j,k)  = nucin2(l)%mode(3)%n(i,j,k) + fcin * nucin2(l)%mode(2)%n(i,j,k)*dt0
      nucin2(l)%mode(3)%nc(i,j,k) = nucin2(l)%mode(3)%nc(i,j,k) + fcin * nucin2(l)%mode(2)%nc(i,j,k)*dt0
      nucin2(l)%mode(2)%n(i,j,k)  = nucin2(l)%mode(2)%n(i,j,k) - fcin * nucin2(l)%mode(2)%n(i,j,k)*dt0
      nucin2(l)%mode(2)%nc(i,j,k) = nucin2(l)%mode(2)%nc(i,j,k) - fcin * nucin2(l)%mode(2)%nc(i,j,k)*dt0
#ifdef NUC_CNT1
      nucin2(l)%mode(3)%m(i,j,k)  = nucin2(l)%mode(3)%m(i,j,k) + fcin * nucin2(l)%mode(2)%m(i,j,k)*dt0
      nucin2(l)%mode(3)%mc(i,j,k) = nucin2(l)%mode(3)%mc(i,j,k) + fcin * nucin2(l)%mode(2)%mc(i,j,k)*dt0
      nucin2(l)%mode(2)%m(i,j,k)  = nucin2(l)%mode(2)%m(i,j,k) - fcin * nucin2(l)%mode(2)%m(i,j,k)*dt0
      nucin2(l)%mode(2)%mc(i,j,k) = nucin2(l)%mode(2)%mc(i,j,k) - fcin * nucin2(l)%mode(2)%mc(i,j,k)*dt0
#endif
  enddo
#endif
!
!----------------------------------------------------------!
!
return
end
!
!  ==================================================
   subroutine focdw ( h, qc, nc, qi, ni, fin, Tem, dtrad, ww, dqcdt, sw, fcdw, fcndw )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Immersion freezing of cloud droplets following Diehl and Wurzler 2004:
!         - Freezing efficiencies for illite and soot based on lab data
!	  - Freezing point depression obtained using ice water activity from Koop et al., Nature 2000
!
! ================================================================
!
  integer  :: h, i
  real     :: nc, qc, qi, ni, nin, fin, Tem, dTdt, dqcdt, dtrad, ww, sw, fcdw, fcndw
  real     :: d, v, dtf, tvar, fall, ff(3)
!
#ifndef NUC_CNT
!
!----------------------------------------------------------!
!            Calculate immersion freezing
!----------------------------------------------------------!
!
  dtf = cal_fpd ( 1, Tem, sw, qc, nc )
!
  v = pi/6.*cal_avd3(h, qc, nc)
  do i = 1, 3
    ff(i) = V*nucprop(i)%Beff*exp(-nucprop(i)%a*(Tem - T00 - dtf))
  enddo
!
!----------------------------------------------------------!
!    Frozen INs are added to immersed frozen particles     !
!----------------------------------------------------------!
!
  fall = 0.
  fcndw = 0.
  fcdw = 0.
  dTdt = dtrad - g*ww/cp_a
  do i = 1, 3
    tvar = max(dqcdt/qc - nucprop(i)%a*dTdt,0.)
    fall = fall + nucprop(i)%frac * max(min(1. - exp(-ff(i)*tvar*dt0),1.),0.)
  enddo
!
  nin = nc*fin
!
  fcndw = max(min(fall*nin, nin - ni), 0.) / dt0
  fcdw  = fcndw * qc / nc
!
!----------------------------------------------------------!
!
#endif
!
return
end
!
!  ==================================================
   subroutine focphil (qc, nc, Tem, siw, sw, si, rho, fcp, fcnp)
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Immersion freezing of cloud droplets following Phillips et al. 2008
!	 	- So far, only valid for monodisperse INs
!		- no parameterization for contact freezing
!
! ================================================================
!
USE gridno
USE shared_data
USE shared_hydro
USE shared_nuclei
!
  IMPLICIT NONE
!
  real     :: fcp, fcnp, dnid1, dnid2, dnibc1, dnibc2
  real     :: nc, qc, Tc, Tem, siw, si, sw, rho
  real     :: dfoc, fc, dtc, dsi, dsw, si0, hx, xit, xsi, sig, kappa, mux
  real     :: nin1, nin11, nin12, nin1d, nin1bc, omega1d, omega1bc
  real     :: ninad1, ninad2, ninabc1, ninabc2
  real     :: wmin=4.42e-14, alphad=2./3., alphabc=1./3., fffc=0.
!
!----------------------------------------------------------!
!              Calculate nucleation rates		   !
!----------------------------------------------------------!
!
  Tc = Tem - 273.15
  fcnp = 0.
  fcp = 0.
!
  xsi = -1.0261 + 3.1656e-3*Tc + 5.3938e-4*Tc*Tc + 8.2584e-6*Tc*Tc*Tc
  kappa = cal_delta_ab (1., 0., Tc, -35., -25.)
  dsw  = cal_delta_ab (0., 1., sw, 0.97, 1.)
  dfoc = cal_delta_ab (0., 1., fffc, 0.1, 1.)
  sig  = cal_delta_ab (0., 1., 2.-fffc, 1., 2.)*cal_delta_ab (1., 0., fffc, 0.1, 1.)
!
!  Reference IN available surface for normalization of parameterization
!
  omega1d  = 5.e-7
  omega1bc = 1.e-7
!
!----------------------------------------------------------!
!           Calculate active INs concentration		   !
!----------------------------------------------------------!
!
  if ( Tc >= -25. ) then
    nin1   = 154.492*exp(12.96*(si - 1.) - 0.639)
  else if ( Tc <= -35. ) then
    nin1   = 2631.58*exp(12.96*(si - 1.1))**0.3
  else if ( Tc < -25. .and. Tc > -30. ) then
    nin11   = 154.492*exp(12.96*(siw - 1.) - 0.639)
    nin12   = 2631.58*exp(12.96*(siw - 1.1))**0.3
    nin1    = min(nin11,nin11) * (min(nin12,nin11)/min(nin11,nin11))**kappa 
    nin1    = min(nin1,nin11)
  else if ( Tc < -30. .and. Tc > -35. ) then
    nin11   = 154.492*exp(12.96*(siw - 1.) - 0.639)
    nin12   = 2631.58*exp(12.96*(siw - 1.1))**0.3
    nin1    = min(nin11,nin12) * (min(nin12,nin12)/min(nin11,nin12))**kappa 
    nin1    = min(nin1,nin12)
  endif
!
!----------------------------------------------------------!
!             Number of nucleated INs for dust		   !
!----------------------------------------------------------!
!
  if ( sw < 1. ) then
    si0 = 1. + 10.**(xsi)
    dsi = cal_delta_ab (0., 1., si, si0, si0+0.1)
    dtc = cal_delta_ab (1., 0.15, Tc, -40., -35.)
    fc  = 0.5*dtc*dsi
!
    hx = min(fc + (1. - fc)*dsw, 1.)
  else
    hx = 1.
  endif 
!
  xit = cal_delta_ab (1., 0., Tc, -30., -10.)
  mux = hx * xit * alphad*nin1/omega1d
!
!  Number of activated dust particles (consider all the dust particles, that
!  is interstitial+immersed+frozen)
!
!  ninad1 = nucin_tmp%species(1)%mode(1)%n*(1. - exp(-pi*nucp(1)%d**2.*mux)) 
!  ninad2 = nucin_tmp%species(1)%mode(1)%nc*(1. - exp(-pi*nucp(1)%dc**2.*mux))
!
!----------------------------------------------------------!
!              Number of nucleated INs for BC		   !
!----------------------------------------------------------!
!
  if ( sw < 1. ) then
    si0 = 1.3 + dfoc * (1.2*siw - 1.3)
    dsi = cal_delta_ab (0., 1., si, si0, si0+0.1)
    dtc = cal_delta_ab (1., 0., Tc, -50., -40.)
    fc  = 0.5*dtc*dsi
!  
    hx = min(fc + (1. - fc)*dsw, 1.)
  else
    hx = 1.
  endif
!
  xit = cal_delta_ab (1., 0., Tc, -25., -15.)
  mux = hx * xit * sig*alphabc*nin1/omega1bc
!
!  Number of activated BC particles (consider all the BC particles, that
!  is interstitial+immersed+frozen)
!
!  ninabc1 = nucin_tmp%species(2)%mode(1)%n*(1. - exp(-pi*nucp(2)%d**2.*mux)) 
!  ninabc2 = nucin_tmp%species(2)%mode(1)%nc*(1. - exp(-pi*nucp(2)%dc**2.*mux)) 
!
!----------------------------------------------------------!
!              Calculate nucleation rates		   !
!----------------------------------------------------------!
!
!  Ice production rates
!
!  dnid1 = max(ninad1 - nucin_tmp%species(1)%mode(5)%n,0.)/dt
!  dnid2 = max(ninad2 - nucin_tmp%species(1)%mode(5)%nc,0.)/dt
!  dnibc1 = max(ninabc1 - nucin_tmp%species(2)%mode(5)%n,0.)/dt
!  dnibc2 = max(ninabc2 - nucin_tmp%species(2)%mode(5)%nc,0.)/dt
!
  fcnp = fcnp + dnid1 + dnid2 + dnibc1 + dnibc2
  if ( sw >= 1. ) then
    fcp = fcp + fcnp * qc/nc
  else
    fcp = fcp + fcnp * wmin
  endif
!
return
end
!
end module freezing
