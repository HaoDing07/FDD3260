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
  module activation
!
  USE gridno
  USE shared_data
  USE shared_nuclei
  USE shared_aerosol_new
  USE shared_thermo
  USE shared_hydro
  USE shared_pressure
  USE shared_wind
  USE shared_rad
  USE shared_diag
  USE aeroactivate
  USE thermodynamics
  USE allocation

  ! USE time_step
!
  IMPLICIT NONE

  private
  
  public :: aero_activ, CCN_activ



  contains
!
! ==========================================
  subroutine CCN_activ ( dref, T, w, dTdt, pressurel, statel, hydromtrl, nucl, statem, hydromtrm )
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
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: T, dref, w, dTdt
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl, hydromtrm
  type (atm_state) :: statel, statem
  type (atm_pressure) :: pressurel  
  type (nuclei) :: nucl
!
  integer :: i, j, k, h, nt
  real    :: s, tem, dnc, dqc, des, sumh, qv, pp, exn
  real    :: r0, rd0, rw0, sd, xxx, xnc00, ak, s0, b
!
  if (verbose > 0) call write_debug('Start CCN_activ_simple')
!
!----------------------------------------------------------!
!                        Starting                          !
!----------------------------------------------------------!
!
  sd  = xnc0_s
  b   = xnc0_k
  rd0 = 0.5*xnc0_d
!
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
!
        sumh = 0.
        do h = 1, nhydro
          sumh = sumh + hydromtrl(h)%q(i,j,k)
        enddo
        qv = statel%qt(i,j,k) - sumh
        tem = T(i,j,k)
#ifdef ANELASTIC
        pp = p0(k)
#else
        pp = pressurel%p(i,j,k) + p0(k)
#endif
        exn = (pp/pref)**((cp_a-cv_a)/cp_a)
!
        s = min(qv*cal_qsw1(tem,pp) - 1., 0.2)
        if (s > 1.e-5) then
          dnc = 0.
          dqc = 0.
          des = 0.
!
!  Activation spectrum
!
          ak  = 2.*cal_sigmawv(tem) / ((cp_v-cv_v)*tem*rho_l)
          rw0 = sqrt(b / ak) * rd0**(1. + beta)
	  r0  = sqrt(3.*b / ak) * rd0**(1. + beta)
	  s0  = sqrt(4.*ak*ak*ak/(27.*b)) * rd0**(-(1. + beta))
!
	  xxx = log(s0/s) / ((1.+beta)*sqrt(2.)*log(sd))
          xnc00 = 0.5 * nucl%ccn(i,j,k) * (1. - calerf( xxx, 0 ))
!
          dnc = max(xnc00 - hydromtrl(drop)%n(i,j,k),0.) 
!
#ifndef SAT_ADJ
          dqc = 4./3.*pi*rho_l*(1.e-6)**3. * dnc/dt0
          des = cal_flv(tem)/(exn*cp_a)*dqc
#endif
!
!----------------------------------------------------------!
!                         Updates                          !
!----------------------------------------------------------!
!
!  Activate
!
          if ( dqc > qmin(drop) .and. dnc > xnmin(drop) ) then
	    hydromtrm(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) + dnc/dt0
	    hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) + dqc
#ifdef ISENTROPIC
            statem%es(i,j,k) = statem%es(i,j,k) + des
#endif
!
!  Diagnostics
!
            ! write(7,*) 'in activation',i,j,k,cint,dqc,dref(i,j,k),diag(9)%micro(drop)%q(i,j,k)
            if (ldiag .and. out_diagt) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + dref(i,j,k)*cint*des
	    if (ldiag .and. out_diagl) diag(8)%qc(i,j,k) = diag(8)%qc(i,j,k) + dref(i,j,k)*cint*dqc
            if (ldiag .and. out_micro) then
    	      diag(9)%micro(drop)%q(i,j,k) = diag(9)%micro(drop)%q(i,j,k) + dref(i,j,k)*cint*dqc
	      diag(9)%micro(drop)%n(i,j,k) = diag(9)%micro(drop)%n(i,j,k) + dref(i,j,k)*cint*dnc/dt0
            endif
	  endif
        endif
!
      enddo
    enddo
  enddo
!
!----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate CCN_activ_simple')
!
  return
  end
!	  
! ================================================================
  subroutine aero_activ ( dref, T, w, dTdt, pressurel, statel, hydromtrl, aerol, statem, hydromtrm, aero3dm )
! ================================================================
!
!  Purpose:
!	CCN activation with aerosol module
!

  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: T     !< temperature
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: dref  !< reference density
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: w     !< vertival wind speed
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: dTdt  !< Temperature rate of change
  type (aero_3d), dimension(1:nmode)  :: aerol
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl, hydromtrm
  type (atm_state) :: statel, statem
  type (atm_pressure) :: pressurel  
  type (aero_3d), dimension(1:nmode) :: aero3dm
!
  integer  :: i, j, k, h, im
  real  :: pp, exn
  real  :: Tem, dtem, ww
  real  :: sumh, qv, s
  real  :: nc
!
  real  :: dnc, dqc, des
  real, dimension(nmode) :: remnr                     !< number of activated aerosols per mode
  real, dimension(nmode) :: remmass                   !< mass of activated aerosols per mode
!
#ifdef AERO_ENABLE
!
  if (verbose > 0)  call write_debug('Start aero_activ')
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!
!
  if ( sorder > 1 ) then
    cdiag = 0.5
  else
    cdiag = 1.
  endif
!
  if ( split_mic ) then
    cint = 1.
  else
    cint = cdiag
  endif

  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
!
#ifdef ANELASTIC
        pp = p0(k)
#else
        pp = pressurel%p(i,j,k) + p0(k)
#endif
!
        sumh = 0.
        do h = 1, nhydro
          sumh = sumh + hydromtrl(h)%q(i,j,k)
        enddo
        qv = statel%qt(i,j,k) - sumh
        tem = T(i,j,k)
        dtem = dTdt(i,j,k)
        ww = w(i,j,k)
        exn = (pp/pref)**((cp_a-cv_a)/cp_a)
        nc = hydromtrl(drop)%n(i,j,k)
!
!----------------------------------------------------------!
!                       Activation                         !
!----------------------------------------------------------!
!
        dnc = 0.
        dqc = 0.
        remnr = 0.
        remmass = 0.
!
        call activ_local ( pp, tem, ww, dtem, qv, nc, aerol%n(i,j,k), aerol%m(i,j,k), dnc, dqc, remnr, remmass )
!
!----------------------------------------------------------!
!                        Updates                           !
!----------------------------------------------------------!
       
      !   if ( sorder > 1 ) then
      !     cdiag = 0.5
      !   else
      !     cdiag = 1.
      !   endif
      ! !
      !   if ( split_mic ) then
      !     cint = 1.
      !   else
      !     cint = cdiag
      !   endif

!  Update concentrations  
!
        hydromtrm(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) + dnc/dt0
        hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) + dqc
#ifdef ISENTROPIC
        des = cal_flv(tem)/(exn*cp_a)*dqc
        statem%es(i,j,k) = statem%es(i,j,k) + des
#endif
!
        if (aero_flg%act_scv .and. sum(remnr)>1e-5) then
          do im = 1, nmode
            ! The term dnc/sum(remnr) makes sure, that only activated aerosol get depleted
            ! in case activation is limited by supersaturation
            aero3dm(im)%n(i,j,k)  = aero3dm(im)%n(i,j,k) - dnc/sum(remnr) * remnr(im)/dt0
            aero3dm(im)%m(i,j,k)  = aero3dm(im)%m(i,j,k) - dnc/sum(remnr) * remmass(im)/dt0
            aero3dm(im)%ma(i,j,k) = aero3dm(im)%ma(i,j,k) + dnc/sum(remnr) * remmass(im)/dt0
          enddo
        endif
!
!  Diagnostics
	      ! write(7,*) 'in aero_activ',cint,dqc,dref(i,j,k),diag(9)%micro(drop)%q(i,j,k)
	if (ldiag .and. out_diagl) diag(8)%qc(i,j,k) = diag(8)%qc(i,j,k) + dref(i,j,k)*cint*dqc 
#ifdef ISENTROPIC
	if (ldiag .and. out_diagt) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + dref(i,j,k)*cint*des
#endif
!
        if (ldiag .and. out_micro) then
    	  diag(9)%micro(drop)%q(i,j,k) = diag(9)%micro(drop)%q(i,j,k) + dref(i,j,k)*cint*dqc 
	  diag(9)%micro(drop)%n(i,j,k) = diag(9)%micro(drop)%n(i,j,k) + dref(i,j,k)*cint*dnc/dt0
          if (out_diaga .and. aero_flg%act_scv) then
            do im = 1, nmode
	      diag(5)%aero(im)%n(i,j,k) = diag(5)%aero(im)%n(i,j,k) - dref(i,j,k)*cint*remnr(im)/dt0
	      diag(5)%aero(im)%m(i,j,k) = diag(5)%aero(im)%m(i,j,k) + dref(i,j,k)*cint*remmass(im)/dt0
	      diag(9)%aero(im)%n(i,j,k) = diag(9)%aero(im)%n(i,j,k) - dref(i,j,k)*cint*remnr(im)/dt0
	      diag(9)%aero(im)%m(i,j,k) = diag(9)%aero(im)%m(i,j,k) + dref(i,j,k)*cint*remmass(im)/dt0
            enddo
	  endif
	endif
!
!  Chemicals  
!
#ifdef AQCHEM_ENABLE     
      aqc1(i,j,k)%svi = aqc1(i,j,k)%svi + sum(remmass)
#endif
!
      enddo ! end i loop
    enddo ! end k loop
  enddo ! end j loop
!
!----------------------------------------------------------!
!
  if (verbose > 0)  call write_debug('Terminate aero_activ')
!
#endif
!
  return
  end
!
! ==================================================
  subroutine activ_local ( p, Tem, w, dtem, qv, n, num, mass, dnc, dqc, remnr, remmass )
! ==================================================
!
  real :: p, Tem, dtem, w, qv, n, dnc, dqc
  real, dimension(nmode) :: num, mass, remnr, remmass
!
  ! Local variables
  integer  :: im, nt                                  !< integer indizes
  real  :: gg                                         !<Diffusion coefficient
  real  :: ak                                         !<
  real  :: alpha,gamm                                 !<Saturation change rate
  real  :: s, s1                                      !<Saturation
!
#ifdef AERO_ENABLE
  real, dimension(nmode) :: rgeo                      !< geometric radius of aerosol
  real, dimension(nmode) :: fact                      !< activated fraction
  real, dimension(nmode) :: mact                      !< activated mass
  real, dimension(nmode) :: rwet                      !< wet radius
  real, dimension(nmode) :: rdry                      !< dry radius
  real :: dtl                                         !< sub time step helper
  real :: m0                                          !< mass of a 1 micron droplet
  real :: ds                                          !< change in supersaturation
  real :: dtt                                         !< Temperature change in sub-timestep
  real :: tt                                          !< sub time step helper
!
  ! Arrays for sub time step
  real  :: sc(100,nmode)                               !< equilibrium supersaturation for given wet radius
  real  :: dq(100,nmode)                               !< total change in water 
  real  :: rd(100,nmode)                               !< dry radius of activated particles
  real  :: rl(100,nmode)                               !< wet radius of activated particles
  real  :: dr(100,nmode)                               !< growth in radius 
  real  :: dn(100,nmode)                               !< change in numer density
!
!----------------------------------------------------------!
!
!  Local variables
!
  dtl = dt0
!
  ! Growth coefficient Eq. A4 Ghan et al., 2011
  gg = 1. / ((Rw*tem)/(cal_esw(tem)*cal_xfkd(tem,p)) + cal_flv(tem)/(Kt*tem)*(cal_flv(tem)/(Rw*tem) - 1.))
  ! Curvature effect
  ak  = 2.*cal_sigmawv(Tem) / ((cp_v-cv_v)*Tem*rhow) 
  ! Rate supersaturation change
  alpha = cal_flv(tem)*(dtem - g*w/cp_a)/(cp_v-cv_v)/tem/tem + g*w/(cp_a-cv_a)/tem
  gamm = 1./qv + cal_flv(Tem)*cal_flv(Tem)/(cp_a*(cp_v-cv_v)*Tem*Tem)  
!
!  Starting saturation
!
  s = min(qv*cal_qsw1(tem,p) - 1., 0.2)
  if ( .not.lkin ) s = min( max(s, s + dt0*(s+1.)*alpha), 0.2 )
!
  do im = 1, nmode
    rgeo(im) = cal_rgeo (aero(im)%size%sigma, aero(im)%mix%rho, mass(im), num(im))
  enddo
!
!----------------------------------------------------------!
!             Calculate activated fraction                 !
!----------------------------------------------------------!
!
  dq = 0.
  dn = 0.
  dr = 0.
  rl = 1.e9
  rd = 0.
  sc = 0.
  dqc = 0.
  dnc = 0.
!
  tt = 0.
  nt = 1
  ds = 1.
  s1 = 0.
!
  do while ( abs(ds) > 1.e-9 .and. tt < dt0 .and. nt < 20 )
      remnr = 0.
      remmass = 0.
      fact = 0.
      mact = 0.
      rwet = 1.e9
      rdry = 0.
!
      do im = 1, nmode
        if ( s > 1.e-5 .and. num(im) > naero_min .and. mass(im) > qaero_min ) then
           call activate (s, s1, rgeo(im), aero(im)%size%sigma, aero(im)%mix%rho, ak, aero(im)%mix%kappa,   &
			  fact(im), mact(im), rwet=rwet(im), ract=rdry(im))
!
           ! Number and mass of activated particles
           remnr(im) = fact(im)*num(im)
           remmass(im) = mact(im)*num(im)
	endif
      enddo
!
      if (.not.aero_flg%act_scv .and. sum(fact) > 1.e-9) remnr = fact/sum(fact)*max(sum(remnr) - n,0.)
      ! Number of activated droplets
      dnc = sum(remnr)
      s1 = s
!
!  Kinetic growth
!
#ifndef SAT_ADJ
      if (lkin) then
!
        dtt = dtem
        ds = (s + 1.)*alpha
        do im = 1, nmode
          ! Number of activated particles during iteration
          dn(nt,im) = remnr(im) - sum(dn(1:nt,im))
          ! Mean wet size of activated aerosols during iteration
          rl(nt,im) = rwet(im)
          ! Mean dry size of activated aerosols during iteration
          rd(nt,im) = rdry(im) 
!
          ! Critical supersaturation for all bins
          sc(1:nt,im) = ak/rl(1:nt,im) - aero(im)%mix%kappa*rd(1:nt,im)**3./rl(1:nt,im)**3. 
          ! Rate of change of wet particle size in each bin
          dr(1:nt,im) = gg*max(s - sc(1:nt,im),0.)/rl(1:nt,im)/rho_l
          ! Rate of change of water mass in each bin
          dq(1:nt,im) = 4.*pi*gg*rl(1:nt,im)*max(s - sc(1:nt,im),0.)*dn(1:nt,im)
!
          ! Temperature change during iteration
          dtt = dtt + cal_flv(Tem)/cp_a*sum(dq(1:nt,im))
          ! Supersaturation change during iteration
          ds  = ds - (s  + 1.)*gamm*sum(dq(1:nt,im))
          ! Total water change during iteration  
          dqc = dqc + sum(dq(1:nt,im)) 
        enddo
!
        ! Determine stable small time-step
        dtl = min(0.5*(s+1.)/abs(ds), 0.2, dt0-tt)
        ! Advance properties before next iteration
        tem = tem + dtt*dtl
        s   = s   + ds*dtl
        rl  = rl  + dr*dtl
        ! Update time and step
        tt = tt + dtl
        nt = nt + 1
      else
	! else, all activated droplets have the size of wet CCNs
        m0 = 0.
        do im = 1, nmode
          m0 = m0 + 4./3.*pi*rhow*sum(rwet**3.)
        enddo
        dqc = min(m0*dnc,qv-cal_qsw(Tem,p))/dt0
	dnc = dt0 * dqc / m0
	tt = dt0
      endif
#else
      ! else, all activated droplets have the same size: 6micron
      m0 = 4./3.*pi*rho_l*(6.e-6)**3.	!Minimum cloud droplet size of 6mum
      dqc = min(m0*dnc,qv-cal_qsw(Tem,p))/dt0
      dnc = dt0 * dqc / m0
      tt = dt0
#endif
!
  enddo    !< end while loop
#endif
!
!----------------------------------------------------------!
!
  return
  end
!
! ==================================================
  subroutine in_activ ( dref, T, rh, pressurel, statel, hydromtrl, nucin_tmp )
! ==================================================
!
! ================================================================
!
  integer :: i, j, k, h
  real    :: s, tem, dnc, dqc
  real    :: r0, r0c, rd0, rd0c, sd, sdc, s0, s0c
  real    :: xxx, xnc00, xnc00c, ak, b
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: rh, T, dref
  type(nuclei_3d), dimension(3) :: nucin_tmp
  type (atm_state) :: statel
  type (atm_pressure) :: pressurel  
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
!
#ifdef NUC_CNT
!
!----------------------------------------------------------!
!             Calculate ice nuclei activation              !
!----------------------------------------------------------!
!
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
!
!  Thermo properties (limit supersaturation to 5%)
!
        tem = T(i,j,k)
        s   = min(rh(i,j,k)/100. - 1.,0.05)
        ak  = 2.*cal_sigmawv(tem) / ((cp_v-cv_v)*tem*rho_l)
!
!  Aerosol properties: limit standard deviation by 1.2
!
	do h = 1, 3
	  if ( s > 0.0 .and. nucin_tmp(h)%mode(1)%n(i,j,k) + nucin_tmp(h)%mode(1)%nc(i,j,k) > xin_min ) then
!
	    if ( h == 1 ) then
	      b = 0.01
	    else if ( h == 2 ) then
	      b = 0.1
	    else if ( h == 3 ) then
	      b = 0.25
	    endif
!
!  Activation spectrum
!
#ifdef NUC_CNT1
	    rd0  = cal_rgeo ( nucp(h)%s, nucp(h)%rhop, nucin_tmp(h)%mode(1)%m(i,j,k), nucin_tmp(h)%mode(1)%n(i,j,k) )
	    rd0c = cal_rgeo ( nucp(h)%sc, nucp(h)%rhop, nucin_tmp(h)%mode(1)%mc(i,j,k), nucin_tmp(h)%mode(1)%nc(i,j,k) )
#else
	    rd0 = nucp(h)%d * exp(log(nucp(h)%s)*log(nucp(h)%s))
	    rd0c = nucp(h)%dc * exp(log(nucp(h)%sc)*log(nucp(h)%sc))
#endif
!
	    r0 = sqrt(3.*b / ak) * rd0**1.5
	    s0 = 2./3.*ak/r0
	    xxx = log(s/s0) / (1.5*sqrt(2.)*log(nucp(h)%s))
            xnc00 = 0.5 * nucin_tmp(h)%mode(1)%n(i,j,k) * (1. + min(max(calerf(xxx,0),-1.),1.))
!
	    r0c = sqrt(3.*b / ak) * rd0c**1.5
	    s0c = 2./3.*ak/r0c
	    xxx = log(s/s0c) / (1.5*sqrt(2.)*log(nucp(h)%sc))
            xnc00c = 0.5 * nucin_tmp(h)%mode(1)%nc(i,j,k) * (1. + min(max(calerf(xxx,0),-1.),1.))
!      
!  Activate
!	
	    dnc = max(xnc00 - nucin_tmp(h)%mode(2)%n(i,j,k),0.)
            nucin_tmp(h)%mode(2)%n(i,j,k) = nucin_tmp(h)%mode(2)%n(i,j,k) + dnc
!
#ifdef NUC_CNT1
	    dqc = max(4./3.*pi*nucp(h)%rhop*r0**3.*xnc00 - nucin_tmp(h)%mode(2)%m(i,j,k),0.)
            nucin_tmp(h)%mode(2)%m(i,j,k) = nucin_tmp(h)%mode(2)%m(i,j,k) + dqc
#endif
!	
	    dnc = max(xnc00c - nucin_tmp(h)%mode(2)%nc(i,j,k),0.)
            nucin_tmp(h)%mode(2)%nc(i,j,k) = nucin_tmp(h)%mode(2)%nc(i,j,k) + dnc
!	    
#ifdef NUC_CNT1
	    dqc = max(4./3.*pi*nucp(h)%rhop*r0c**3.*xnc00c - nucin_tmp(h)%mode(2)%mc(i,j,k),0.)
            nucin_tmp(h)%mode(2)%mc(i,j,k) = nucin_tmp(h)%mode(2)%mc(i,j,k) + dqc
#endif
!
	  endif
	enddo
!
      enddo
    enddo
  enddo
!
!----------------------------------------------------------!
!
#endif
!
return
end
!
end module activation
