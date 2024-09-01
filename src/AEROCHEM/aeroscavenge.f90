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
!  AEROSCAVENGE.F                   
!
!  Purpose:
!	Subroutines for calculating scavenging of aerosols
!
!  Author
!	Julien Savre
!       MISU
!
! ================================================================
!
!	==========================================
	subroutine impaction_scav ( dref, aero_tmp, hydromtr_tmp, aero3dm )
!     	==========================================
!
USE gridno
USE shared_aq
USE shared_solid
USE shared_gas
USE shared_data
USE shared_nuclei
USE shared_hydro
USE shared_thermo
USE shared_pressure
USE shared_aerosol_new
USE shared_diag
!
IMPLICIT NONE
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dref
  type (aero_3d), dimension(1:nmode)  :: aero_tmp
  type(hydrometeor), dimension(1:nhydro)  :: hydromtr_tmp
  type (aero_3d), dimension(1:nmode) :: aero3dm
!
#ifdef AERO_ENABLE
!
  integer :: im, i, j, k
  real :: p, rho, Tem, xfmu
  real :: mass, num, dens, dav, mav, c_a(nhydro), cnacc, remnr, remmass
!
!----------------------------------------------------------!
!                     External loops                       !
!----------------------------------------------------------!
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end 
!
        do im = 1, nmode
          c_a  = 0.
!
          num  = aero_tmp(im)%n(i,j,k) 
          mass = aero_tmp(im)%m(i,j,k)
          dens = aero(im)%mix%rho
          dav  = 2.*cal_rav( aero(im)%size%sigma, dens, mass, num )
          mav  = cal_mav( aero(im)%size%sigma, dens, mass, num )
!
          p = p0(k)
          rho = den0(k)
          Tem = thermo%T(i,j,k)
          xfmu = cal_xfmu( Tem )
!
!----------------------------------------------------------!
!            Aerosol collection by hydrometeors            !
!              Collection efficiency is 1.e-2              !
!----------------------------------------------------------!
!
!  Scavenging of aerosols by rain
!
          if (hydromtr_tmp(rain)%n(i,j,k) > xnmin(rain) .and. lmicro > 0) 		&
            call contact ( rain, k, hydromtr_tmp(rain)%q(i,j,k), hydromtr_tmp(rain)%n(i,j,k), p, rho, &
			Tem, xfmu, dav, dens, c_a(2) )

!
!  By ice particles
!
          if(hydromtr_tmp(ice)%n(i,j,k) > xnmin(ice) .and. lmicro > 1)		&
            call contact ( ice, k, hydromtr_tmp(ice)%q(i,j,k), hydromtr_tmp(ice)%n(i,j,k), p, rho, &
			Tem, xfmu, dav, dens, c_a(3) )
!
!  By graupel
!
          if(hydromtr_tmp(grau)%n(i,j,k) > xnmin(grau) .and. lmicro > 2)		&
            call contact ( grau, k, hydromtr_tmp(grau)%q(i,j,k), hydromtr_tmp(grau)%n(i,j,k), p, rho, &
			Tem, xfmu, dav, dens, c_a(4) )
!
!  By snow
!
          if(hydromtr_tmp(snow)%n(i,j,k) > xnmin(snow) .and. lmicro > 2)		&
            call contact ( snow, k, hydromtr_tmp(snow)%q(i,j,k), hydromtr_tmp(snow)%n(i,j,k), p, rho, &
			Tem, xfmu, dav, dens, c_a(5) )
!
!----------------------------------------------------------!
!                      Scavenging                          !
!----------------------------------------------------------!
!
          cnacc = 1. - exp(-dt0*sum(c_a))
!   
          remnr   = cnacc * num
          remmass = remnr * mav 
!
          aero3dm(im)%n(i,j,k) = aero3dm(im)%n(i,j,k) - remnr / dt0
          aero3dm(im)%m(i,j,k) = aero3dm(im)%m(i,j,k) - remmass / dt0
          aero3dm(im)%ma(i,j,k) = aero3dm(im)%ma(i,j,k) + remmass / dt0
!      
          diag(7)%aero(im)%n(i,j,k) = diag(7)%aero(im)%n(i,j,k) - cint*remnr / dt0
          diag(7)%aero(im)%m(i,j,k) = diag(7)%aero(im)%m(i,j,k) + cint*remmass / dt0
          diag(9)%aero(im)%n(i,j,k) = diag(9)%aero(im)%n(i,j,k) - cint*remnr / dt0
          diag(9)%aero(im)%m(i,j,k) = diag(9)%aero(im)%m(i,j,k) + cint*remmass / dt0
!
        enddo
!
!----------------------------------------------------------!
!                       End loops                          !
!----------------------------------------------------------!
!
      enddo
    enddo
  enddo
!
#endif
!
!----------------------------------------------------------!
!
return
end
!
!  ==================================================
   subroutine contact ( h, k, q, n, p, rho, Tem, xfmu, dd, densp, Kf )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Calculates impactions between drops and ice nuclei
!
! ================================================================
!
USE gridno
USE shared_data
USE shared_hydro
!
IMPLICIT NONE
!
  integer :: h, k, i
  real :: p, Tem, xfmu, rho, q, n, m, dd, densp
  real :: bb, csc, Vc, d, v, Dc, vp, d2v
  real :: Re, Sc, St, Sts, tau, lair, Kn, alpha, E
!
  real :: Kf
!
!----------------------------------------------------------!
!           Calculate collisions between droplets	   !
!        and all aerosols due to Brownian motions and      !
!              thermophoresis (no ventilation)	           !
!----------------------------------------------------------!
!
  m   = cal_avm(h, q, n)
  d   = hydrop(h)%cm*m**(hydrop(h)%am)
  v   = hydrop(h)%ctv*m**(hydrop(h)%btv)
  Dc  = cal_avd(h, q, n)
  Vc  = cal_vpn(h, q, n, 1.)
  d2v = hydrop(h)%cm*hydrop(h)%cm*hydrop(h)%ctv*hydroc(h)%cimp*cal_lambda(h,q,n)**(-(2.*hydrop(h)%am+hydrop(h)%btv)/hydrop(h)%mu)
!
  Kf = 0.
  lair = 21.55e-6*xfmu*sqrt(Tem)/p
!
  if (dd > 1.e-9 .and. Dc > 1.e-6) then
!
!  Brownian motions
!
    Kn = lair / dd
    alpha = 1.257 + 0.4*exp(-0.55/Kn)
    csc  = 1. + 2.*alpha*Kn
    bb  =  csc*kb*Tem / (3.*pi*xfmu*dd)
!
!  Non dimensional numbers
!
    Re  = 0.5*rho*Vc*Dc / xfmu
    Sc  = xfmu / bb / rho
!
!  Inertial impaction
!
    vp  = g*csc*(densp - rho)*dd**2. / (18.*xfmu)
    tau = vp / g
    St  = 2.*tau * (Vc - vp) / Dc
    Sts = (1.2 + 0.08333*log(1. + Re)) / (1. + log(1. + Re))
!
!  Collision efficiency
!
    E  = 4./(Re*Sc) * (1. + 0.4*sqrt(Re)*Sc**(1./3.) + 0.16*sqrt(Re*Sc))	&
       + 4.*dd/Dc * (xfmu/0.0018 + (1. + 2.*sqrt(Re))*dd/Dc)
!
    if (with_collision .and. St > Sts ) E = E + ((St - Sts) / (St - Sts + 2./3.))**(1.5) * sqrt(densp/rho)

	if (with_collision == .false.) E = 0.0
!
    Kf  = pi/4. * max(min(E,1.),0.) * d2v * n 
  endif
!
!----------------------------------------------------------!
!
return
end
!
#ifdef NUC_CNT
!
!  ==================================================
   subroutine in_scavenge ( nucin_tmp, hydromtr_tmp, pressure_tmp )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Calculates ice nuclei scavenging by impaction
!
! ================================================================
!
USE gridno
USE shared_data
USE shared_aerosol_new
USE shared_nuclei
USE shared_hydro
USE shared_thermo
USE shared_pressure
!
IMPLICIT NONE
!
  integer :: i, j, k, h
  real    :: pp, qv, rho, Tem, xfmu, xfkd, e, esw, rhov, rhovs, fn, fnc
  type(nuclei_3d), dimension(3) :: nucin_tmp
  type (hydrometeor), dimension(1:nhydro) :: hydromtr_tmp
  type (atm_pressure) :: pressure_tmp
!
  real, dimension(1:3) :: Kcont, Kcontc, Kconti, Kcontci, Kfd, Kfcd, Kfr, Kfcr, Kfi, Kfci, Kfg, Kfcg, Kfs, Kfcs, Kfh, Kfch
  real, dimension(1:3) :: dd, ddc, nn, nnc, mm, mmc
!
!----------------------------------------------------------!
!    Calculate impaction scavenging using Slinn formula    !
!----------------------------------------------------------!
!
    do k = 1, nz
      do j = jt_start, jt_end
	do i = it_start, it_end
          Kcont  = 0.
          Kcontc = 0.
          Kconti  = 0.
          Kcontci = 0.
          Kfd	 = 0.
          Kfcd   = 0.
          Kfr	 = 0.
          Kfcr   = 0.
          Kfi	 = 0.
          Kfci   = 0.
          Kfg	 = 0.
          Kfcg   = 0.
          Kfs	 = 0.
          Kfcs   = 0.
          Kfh	 = 0.
          Kfch   = 0.
!
#ifdef ANELASTIC
	  pp = p0(k)
#else
          pp  = pressure_tmp%p(i,j,k) + p0(k)
#endif
!
	  Tem = thermo%T(i,j,k) 
	  qv  = statel%qt(i,j,k)
	  if (lmicro > 0) qv = qv - hydromtrl(drop)%q(i,j,k) - hydromtrl(rain)%q(i,j,k)
	  if (lmicro > 1) qv = qv - hydromtrl(ice)%q(i,j,k)
	  if (lmicro > 2) qv = qv - hydromtrl(grau)%q(i,j,k) - hydromtrl(snow)%q(i,j,k)
	  if (lmicro > 3) qv = qv - hydromtrl(hail)%q(i,j,k)
          rho = pressure_tmp%dens(i,j,k)
!
	  e = pp * qv/0.622 / (1. + qv/0.622)
	  esw = cal_esw (Tem)
	  rhov = rho*0.622*e/pp
	  rhovs = rho*0.622*esw/pp
!
	  xfkd = cal_xfkd (Tem, pp)
	  xfmu = cal_xfmu (Tem)
!
          do h = 1, 3
#ifdef NUC_CNT1
            dd(h) = cal_logmom (1., nucp(h)%s, nucin_tmp(h)%mode(2)%m(i,j,k), nucin_tmp(h)%mode(2)%n(i,j,k), nucp(h)%rhop)
            ddc(h)= cal_logmom (1., nucp(h)%sc, nucin_tmp(h)%mode(2)%mc(i,j,k), nucin_tmp(h)%mode(2)%nc(i,j,k), nucp(h)%rhop)
#else
	    dd(h) = nucp(h)%d * exp(1.5*log(nucp(h)%s)*log(nucp(h)%s))
	    ddc(h) = nucp(h)%dc * exp(1.5*log(nucp(h)%sc)*log(nucp(h)%sc))
#endif
          enddo
!
!  Scavenging of aerosols by cloud drops
!
!          if (hydromtr_tmp(drop)%n(i,j,k) > xnmin(drop) .and. lmicro > 0) then
!              call contact (drop, k, hydromtr_tmp(drop)%q(i,j,k), hydromtr_tmp(drop)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, dd, nucp(h)%rhop, Kfd )	 
!              call contact (drop, k, hydromtr_tmp(drop)%q(i,j,k), hydromtr_tmp(drop)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, ddc, nucp(h)%rhop, Kfcd )	 
!	  endif
!
!  By rain
!
!          if (hydromtr_tmp(rain)%n(i,j,k) > xnmin(rain) .and. lmicro > 0) then
!              call contact (rain, k, hydromtr_tmp(rain)%q(i,j,k), hydromtr_tmp(rain)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, dd, nucp%rhop, Kfr )      
!              call contact (rain, k, hydromtr_tmp(rain)%q(i,j,k), hydromtr_tmp(rain)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, ddc, nucp%rhop, Kfcr )    
!	  endif
!
!  By ice particles
!
!	  if ( lmicro > 1 ) then
!            if (hydromtr_tmp(ice)%n(i,j,k) > xnmin(ice)) then
!              call contact (ice, k, hydromtr_tmp(ice)%q(i,j,k), hydromtr_tmp(ice)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, dd, nucp%rhop, Kfi )	  
!              call contact (ice, k, hydromtr_tmp(ice)%q(i,j,k), hydromtr_tmp(ice)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, ddc, nucp%rhop, Kfci )    
!	    endif
!	  endif
!
!  By graupel
!
!	  if ( lmicro > 2 ) then
!            if (hydromtr_tmp(grau)%n(i,j,k) > xnmin(grau)) then
!              call contact (grau, k, hydromtr_tmp(grau)%q(i,j,k), hydromtr_tmp(grau)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, dd, nucp%rhop, Kfg )      
!              call contact (grau, k, hydromtr_tmp(grau)%q(i,j,k), hydromtr_tmp(grau)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, ddc, nucp%rhop, Kfcg )    
!	    endif
!
!  By snow
!
!            if (hydromtr_tmp(snow)%n(i,j,k) > xnmin(snow)) then
!              call contact (snow, k, hydromtr_tmp(snow)%q(i,j,k), hydromtr_tmp(snow)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, dd, nucp%rhop, Kfs )      
!              call contact (snow, k, hydromtr_tmp(snow)%q(i,j,k), hydromtr_tmp(snow)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, ddc, nucp%rhop, Kfcs )    
!	    endif
!	  endif
!
!  By hail
!
!	  if ( lmicro > 3 ) then
!            if (hydromtr_tmp(hail)%n(i,j,k) > xnmin(hail)) then
!              call contact (hail, k, hydromtr_tmp(hail)%q(i,j,k), hydromtr_tmp(hail)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, dd, nucp%rhop, Kfh )      
!              call contact (hail, k, hydromtr_tmp(hail)%q(i,j,k), hydromtr_tmp(hail)%n(i,j,k), pp, rho, Tem, Tem, rhov, rhovs, xfmu, xfkd, ddc, nucp%rhop, Kfch )    
!	    endif
!	  endif
!
!----------------------------------------------------------!
!                      Update INs         	           !
!----------------------------------------------------------!
!
	  do h = 1, 3
	    nn(h) = nucin_tmp(h)%mode(1)%n(i,j,k) - nucin_tmp(h)%mode(2)%n(i,j,k) - nucin_tmp(h)%mode(3)%n(i,j,k)
	    nnc(h) = nucin_tmp(h)%mode(1)%nc(i,j,k) - nucin_tmp(h)%mode(2)%nc(i,j,k) - nucin_tmp(h)%mode(3)%nc(i,j,k)
#ifdef NUC_CNT1
	    mm(h) = nucin_tmp(h)%mode(1)%m(i,j,k) - nucin_tmp(h)%mode(2)%m(i,j,k) - nucin_tmp(h)%mode(3)%m(i,j,k)
	    mmc(h) = nucin_tmp(h)%mode(1)%mc(i,j,k) - nucin_tmp(h)%mode(2)%mc(i,j,k) - nucin_tmp(h)%mode(3)%mc(i,j,k)
#endif
!
!  Liquid hydrometeors
!
	    fn = Kfd(h)*hydromtr_tmp(drop)%n(i,j,k) + Kfr(h)*hydromtr_tmp(rain)%n(i,j,k)
	    fnc = Kfcd(h)*hydromtr_tmp(drop)%n(i,j,k) + Kfcr(h)*hydromtr_tmp(rain)%n(i,j,k)
	    Kcont(h) = 1. - exp(-fn*dt0)
	    Kcontc(h) = 1. - exp(-fnc*dt0)
!
            nucins(h)%mode(2)%n(i,j,k)  = nucins(h)%mode(2)%n(i,j,k) + nn(h) * Kcont(h) / dt0
            nucins(h)%mode(2)%nc(i,j,k) = nucins(h)%mode(2)%nc(i,j,k) + nnc(h) * Kcontc(h) / dt0
#ifdef NUC_CNT1
            nucins(h)%mode(2)%m(i,j,k)  = nucins(h)%mode(2)%m(i,j,k) + mm(h) * Kcont(h) / dt0
            nucins(h)%mode(2)%mc(i,j,k) = nucins(h)%mode(2)%mc(i,j,k) + mmc(h) * Kcontc(h) / dt0
#endif
!
!  Ice hydrometeors
!
	    fn = Kfi(h)*hydromtr_tmp(ice)%n(i,j,k) 
	    if ( lmicro > 2 ) fn = fn + Kfg(h)*hydromtr_tmp(snow)%n(i,j,k) + Kfs(h)*hydromtr_tmp(grau)%n(i,j,k)
	    fnc = Kfci(h)*hydromtr_tmp(ice)%n(i,j,k)
	    if ( lmicro > 2 ) fnc = fnc + Kfcg(h)*hydromtr_tmp(snow)%n(i,j,k) + Kfcs(h)*hydromtr_tmp(grau)%n(i,j,k)
	    Kconti(h) = 1. - exp(-fn*dt0)
	    Kcontci(h) = 1. - exp(-fnc*dt0)
!
            nucins(h)%mode(3)%n(i,j,k)  = nucins(h)%mode(3)%n(i,j,k) + nn(h) * Kconti(h) / dt0
            nucins(h)%mode(3)%nc(i,j,k) = nucins(h)%mode(3)%nc(i,j,k) + nnc(h) * Kcontci(h) / dt0
#ifdef NUC_CNT1
            nucins(h)%mode(3)%m(i,j,k)  = nucins(h)%mode(3)%m(i,j,k) + mm(h) * Kconti(h) / dt0
            nucins(h)%mode(3)%mc(i,j,k) = nucins(h)%mode(3)%mc(i,j,k) + mmc(h) * Kcontci(h) / dt0
#endif
	  enddo
!
	enddo
      enddo
    enddo
!
return
end
!
!  ==================================================
   subroutine in_evap ( i, j, k, hydromtr_tmp, nucin_tmp, Tc, dnc, dnr, dni, dng, dns, dnh,	&
		        xmirn, xmgrn, xmsrn, xmhrn )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Calculates release of ice nuclei due to evaporation/sublimation/melting
!       Warning: melting sink of an ice particle is xmirn > 0 while the evaporative
!       sink is dni < 0 (evaporation of cloud droplet is dnc < 0)
!
!
! ================================================================
!
USE gridno
USE shared_data
USE shared_nuclei
USE shared_hydro
!
IMPLICIT NONE
!
  integer :: i,j,k, h
  real    :: dnc, dnr, dni, dng, dns, dnh, dnin
  real    :: Tc, xmirn, xmgrn, xmsrn, xmhrn
  real, dimension(1:nhydro) :: q, n
!  
  type(hydrometeor), dimension(1:nhydro) :: hydromtr_tmp
  type(nuclei_3d), dimension(3) :: nucin_tmp
!
!----------------------------------------------------------!
!           Calculate ice nuclei release due to            !
!                 evaporation/sublimation        	   !
!----------------------------------------------------------!
!
	do h = 1, nhydro
	  q(h) = hydromtr_tmp(h)%q(i,j,k)
	  n(h) = hydromtr_tmp(h)%n(i,j,k)
	enddo
!
!  In liquid drops
!
	if ( lmicro > 0 )  then
	  if ( n(drop) > xnmin(drop) .and. dnc < 0. ) then
	    do h = 1, 3	
	      nucins(h)%mode(2)%n(i,j,k) = nucins(h)%mode(2)%n(i,j,k) + nucin_tmp(h)%mode(2)%n(i,j,k) * min(dnc,0.) / n(drop)
	      nucins(h)%mode(2)%nc(i,j,k) = nucins(h)%mode(2)%nc(i,j,k) + nucin_tmp(h)%mode(2)%nc(i,j,k) * min(dnc,0.) / n(drop)
#ifdef NUC_CNT1
	      nucins(h)%mode(2)%m(i,j,k) = nucins(h)%mode(2)%m(i,j,k) + nucin_tmp(h)%mode(2)%m(i,j,k) * min(dnc,0.) / n(drop)
	      nucins(h)%mode(2)%mc(i,j,k) = nucins(h)%mode(2)%mc(i,j,k) + nucin_tmp(h)%mode(2)%mc(i,j,k) * min(dnc,0.) / n(drop)
#endif
	    enddo
	  endif
!
	  if ( n(rain) > xnmin(rain) .and. dnr < 0. ) then
	    do h = 1, 3	
	      nucins(h)%mode(2)%n(i,j,k) = nucins(h)%mode(2)%n(i,j,k) + nucin_tmp(h)%mode(2)%n(i,j,k) * min(dnr,0.) / n(rain)
	      nucins(h)%mode(2)%nc(i,j,k) = nucins(h)%mode(2)%nc(i,j,k) + nucin_tmp(h)%mode(2)%nc(i,j,k) * min(dnr,0.) / n(rain)
#ifdef NUC_CNT1
	      nucins(h)%mode(2)%m(i,j,k) = nucins(h)%mode(2)%m(i,j,k) + nucin_tmp(h)%mode(2)%m(i,j,k) * min(dnr,0.) / n(rain)
	      nucins(h)%mode(2)%mc(i,j,k) = nucins(h)%mode(2)%mc(i,j,k) + nucin_tmp(h)%mode(2)%mc(i,j,k) * min(dnr,0.) / n(rain)
#endif
	    enddo
	  endif			       
	endif
!
!  In ice crystals
!
	if ( lmicro > 1 )  then
	  if ( n(ice) > xnmin(ice) .and. dni < 0. ) then
	    do h = 1, 3	
	      nucins(h)%mode(3)%n(i,j,k) = nucins(h)%mode(3)%n(i,j,k) + nucin_tmp(h)%mode(3)%n(i,j,k) * min(dni,0.) / n(ice)
	      nucins(h)%mode(3)%nc(i,j,k) = nucins(h)%mode(3)%nc(i,j,k) + nucin_tmp(h)%mode(3)%nc(i,j,k) * min(dni,0.) / n(ice)
#ifdef NUC_CNT1
	      nucins(h)%mode(3)%m(i,j,k) = nucins(i,j,k)%species(h)%mode(3)%m(i,j,k) + nucin_tmp(h)%mode(3)%m(i,j,k) * min(dni,0.) / n(ice)
	      nucins(h)%mode(3)%mc(i,j,k) = nucins(i,j,k)%species(h)%mode(3)%mc(i,j,k) + nucin_tmp(h)%mode(3)%mc(i,j,k) * min(dni,0.) / n(ice)
#endif
	    enddo
	  endif
	endif
!
	if ( lmicro > 2 )  then
	  if ( n(grau) > xnmin(grau) .and. dng < 0. ) then
	    do h = 1, 3	
	      nucins(h)%mode(3)%n(i,j,k) = nucins(h)%mode(3)%n(i,j,k) + nucin_tmp(h)%mode(3)%n(i,j,k) * min(dng,0.) / n(grau)
	      nucins(h)%mode(3)%nc(i,j,k) = nucins(h)%mode(3)%nc(i,j,k) + nucin_tmp(h)%mode(3)%nc(i,j,k) * min(dng,0.) / n(grau)
#ifdef NUC_CNT1
	      nucins(h)%mode(3)%m(i,j,k) = nucins(h)%mode(3)%m(i,j,k) + nucin_tmp(h)%mode(3)%m(i,j,k) * min(dng,0.) / n(grau)
	      nucins(h)%mode(3)%mc(i,j,k) = nucins(h)%mode(3)%mc(i,j,k) + nucin_tmp(h)%mode(3)%mc(i,j,k) * min(dng,0.) / n(grau)
#endif
	    enddo
	  endif
!	  
	  if ( n(snow) > xnmin(snow) .and. dns < 0. ) then
	    do h = 1, 3	
	      nucins(h)%mode(3)%n(i,j,k) = nucins(h)%mode(3)%n(i,j,k) + nucin_tmp(h)%mode(3)%n(i,j,k) * min(dns,0.) / n(snow)
	      nucins(h)%mode(3)%nc(i,j,k) = nucins(h)%mode(3)%nc(i,j,k) + nucin_tmp(h)%mode(3)%nc(i,j,k) * min(dns,0.) / n(snow)
#ifdef NUC_CNT1
	      nucins(h)%mode(3)%m(i,j,k) = nucins(h)%mode(3)%m(i,j,k) + nucin_tmp(h)%mode(3)%m(i,j,k) * min(dns,0.) / n(snow)
	      nucins(h)%mode(3)%mc(i,j,k) = nucins(h)%mode(3)%mc(i,j,k) + nucin_tmp(h)%mode(3)%mc(i,j,k) * min(dns,0.) / n(snow)
#endif
	    enddo
	  endif
	endif
!	  
	if ( lmicro > 3 )  then
	  if ( n(hail) > xnmin(hail) .and. dnh < 0. ) then
	    do h = 1, 3	
	      nucins(h)%mode(3)%n(i,j,k) = nucins(h)%mode(3)%n(i,j,k) + nucin_tmp(h)%mode(3)%n(i,j,k) * min(dnh,0.) / n(hail)
	      nucins(h)%mode(3)%nc(i,j,k) = nucins(h)%mode(3)%nc(i,j,k) + nucin_tmp(h)%mode(3)%nc(i,j,k) * min(dnh,0.) / n(hail)
#ifdef NUC_CNT1
	      nucins(h)%mode(3)%m(i,j,k) = nucins(h)%mode(3)%m(i,j,k) + nucin_tmp(h)%mode(3)%m(i,j,k) * min(dnh,0.) / n(hail)
	      nucins(h)%mode(3)%mc(i,j,k) = nucins(h)%mode(3)%mc(i,j,k) + nucin_tmp(h)%mode(3)%mc(i,j,k) * min(dnh,0.) / n(hail)
#endif
	    enddo
	  endif
	endif
!
!----------------------------------------------------------!
!           Calculate ice nuclei release due to            !
!                 melting of ice particles                 !
!----------------------------------------------------------!
!
	if ( Tc >= 0. )  then
	  if ( lmicro > 1 ) then
	    if ( n(ice) > xnmin(ice) ) then
	      do h = 1, 3
	        nucins(h)%mode(3)%n(i,j,k) = nucins(h)%mode(3)%n(i,j,k) - nucin_tmp(h)%mode(3)%n(i,j,k) * xmirn / n(ice)
	        nucins(h)%mode(3)%nc(i,j,k) = nucins(h)%mode(3)%nc(i,j,k) - nucin_tmp(h)%mode(3)%nc(i,j,k) * xmirn / n(ice)
#ifdef NUC_CNT1
	        nucins(h)%mode(3)%m(i,j,k) = nucins(h)%mode(3)%m(i,j,k) - nucin_tmp(h)%mode(3)%m(i,j,k) * xmirn / n(ice)
	        nucins(h)%mode(3)%mc(i,j,k) = nucins(h)%mode(3)%mc(i,j,k) - nucin_tmp(h)%mode(3)%mc(i,j,k) * xmirn / n(ice)
#endif
	      enddo
	    endif
	  endif
!
	  if ( lmicro > 2 ) then
	    if ( n(grau)+n(snow) > xnmin(grau)+xnmin(snow) ) then
	      do h = 1, 3
	        nucins(h)%mode(3)%n(i,j,k) = nucins(h)%mode(3)%n(i,j,k) - nucin_tmp(h)%mode(3)%n(i,j,k) * (xmgrn + xmsrn) / (n(grau) + n(snow))
	        nucins(h)%mode(3)%nc(i,j,k) = nucins(h)%mode(3)%nc(i,j,k) - nucin_tmp(h)%mode(3)%nc(i,j,k) * (xmgrn + xmsrn) / (n(grau) + n(snow))
#ifdef NUC_CNT1
	        nucins(h)%mode(3)%m(i,j,k) = nucins(h)%mode(3)%m(i,j,k) - nucin_tmp(h)%mode(3)%m(i,j,k) * (xmgrn + xmsrn) / (n(grau) + n(snow))
	        nucins(h)%mode(3)%mc(i,j,k) = nucins(h)%mode(3)%mc(i,j,k) - nucin_tmp(h)%mode(3)%mc(i,j,k) * (xmgrn + xmsrn) / (n(grau) + n(snow))
#endif
	      enddo
	    endif
	  endif
!	  
	  if ( lmicro > 3 ) then
	    if ( n(hail) > xnmin(hail) ) then
	      do h = 1, 3
	        nucins(h)%mode(3)%n(i,j,k) = nucins(h)%mode(3)%n(i,j,k) - nucin_tmp(h)%mode(3)%n(i,j,k) * xmhrn / n(hail)
	        nucins(h)%mode(3)%nc(i,j,k) = nucins(h)%mode(3)%nc(i,j,k) - nucin_tmp(h)%mode(3)%nc(i,j,k) * xmhrn / n(hail)
#ifdef NUC_CNT1
	        nucins(h)%mode(3)%m(i,j,k) = nucins(h)%mode(3)%m(i,j,k) - nucin_tmp(h)%mode(3)%m(i,j,k) * xmhrn / n(hail)
	        nucins(h)%mode(3)%mc(i,j,k) = nucins(h)%mode(3)%mc(i,j,k) - nucin_tmp(h)%mode(3)%mc(i,j,k) * xmhrn / n(hail)
#endif
	      enddo
	    endif
	  endif
	endif
!
!----------------------------------------------------------!
!
return
end
!
#endif
!
!  ==================================================
   subroutine coaero (k,h,lambda,n,can)
!  ==================================================
!
! ================================================================
!  Collection of aerosols by hydrometeors (Mass of aerosol 
!    collected per sec. per unit mass of hydro. (unit s))
!
!  Assumes Seifert and Beheng microphysics.
! ================================================================
!
  integer  :: k, h
  real     :: lambda, n
  real     :: eff=0.01
  real     :: can
!
!----------------------------------------------------------!  
!
!  can = pi/4. * eff * den0(k) * cal_vpn(h,lambda,pxx(k,h)) * cal_avd2(h,lambda) 
!  can = can * n
!
  can = max(0.0,can)
!
return
end
