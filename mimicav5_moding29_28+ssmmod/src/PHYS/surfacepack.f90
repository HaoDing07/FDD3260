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
!  SURFACEPACK:
!	Package of subroutines related to turbulent mixing of all quantities 
!       (wind, pt, qv, hydrometeors) created to simplify and clarify 
!       the code (renamed from sgspack)
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================	 

  module surfacemod
  
  USE gridno
  USE shared_data
  USE shared_state
  USE shared_wind
  USE shared_thermo
  USE shared_pressure
  USE shared_hydro
  USE shared_surf
  USE shared_turbu
  USE thermodynamics
  USE averages
  
  IMPLICIT NONE
  
  public :: calc_surface
  
  CONTAINS

  subroutine calc_surface ( pressurel, windl, statel, hydromtrl )

! ----------------------------------------------
! --- Subroutine calculating friction velocity
! ----------------------------------------------
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_winds), intent(in) :: windl
  type (atm_state), intent(in) :: statel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
!
  integer :: i, j, k
  real  :: windm, l0, lnz, zeta, x, zeff, qv1, qc1, ptv1, tstar, qvstar, ustar,		&
  	   ptvstar, phim, phih, umean, zr, pt1, l0old, delta, zh, 			&
	   windm_mean, wt_mean, wq_mean, u0m, uq0m, ut0m
  real  :: flux_t, flux_q, flux_b, t0_g, qv0_g, exn0_g, ptv0_g, l, sstl
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: pt, q, qv
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: usurf
!
!  Parameters for stability functions
!
  real, parameter :: small=1.e-12, large=1.e10, betam=5.3, betah=5.3, gamma1=16., gamma2=16.
!
  if (verbose > 1) call write_debug('Starting calc_surface')
!
!----------------------------------------------------------!
!                     Initialisations                      !
!----------------------------------------------------------!
!
  q = 0.0
  if (lmicro > 0) then
    q = q + hydromtrl(drop)%q + hydromtrl(rain)%q
  endif
  if (lmicro > 1) then
    q = q + hydromtrl(ice)%q
  endif
  if (lmicro > 2) then
    q = q + hydromtrl(grau)%q + hydromtrl(snow)%q
  endif
  if (lmicro > 3) then
    q = q + hydromtrl(hail)%q
  endif
!
  qv = statel%qt - q
  call get_pt ( pressurel, statel, hydromtrl, pt )
!
!-----------------------------------------------------------!
!            Roughness length and surface winds             !
!-----------------------------------------------------------!
!
  zh = 0.5*dz*fdz0(1)
  if (zrough <= 0.) then
    call horav (surf%ustar, umean)
    zr = max(0.0001,(0.016/g)*umean**2) 
  else
    zr = zrough
  end if
  lnz = log(zh/zr)
!
!  Surface wind speed
!
  usurf = 0.
#ifndef COLUMN
  do j=jt_start,jt_end
    do i=it_start,it_end
      usurf(i,j) = (0.5*(windl%u(i+1,j,1)+windl%u(i,j,1)))**2.
#ifdef MODEL_3D
      usurf(i,j) = usurf(i,j) + (0.5*(windl%v(i,j+1,1)+windl%v(i,j,1)))**2.
#endif
    enddo
  enddo
#endif
  usurf = max( min_w, sqrt(usurf) )
!
!  Other initializations
!
  exn0_g = (psurf/pref)**((cp_a-cv_a)/cp_a)
  ustar = 0.0
  qvstar = 0.0
  tstar = 0.0
!
  if (isurf == 4) then
    call horav (usurf, u0m)
    call horav (usurf*pt(:,:,1), ut0m)
    call horav (usurf*qv(:,:,1), uq0m)
  endif
!
!-----------------------------------------------------------!
!                   Main loop over surface                  !
!-----------------------------------------------------------!
!
  do j=jt_start,jt_end
!
!  Latitude dependent surface values
!
#ifdef CHANNEL
    l = abs(real(j-3) - 0.5*real(ny-5))*dx/111100.
    sstl = 273.15 + 0.5*(sst-273.5)*( (cos(pi*(ctr_lat+l)/60.))**2. + (cos(pi*(ctr_lat+l)/60.))**4. )  
    qv0_g = ssm*cal_qsw(sstl, psurf)
    t0_g  = sstl/exn0_g
#else
    qv0_g = ssm*cal_qsw(sst, psurf)
    t0_g  = sst/exn0_g
#endif
    ptv0_g = t0_g*(1. + epsm*qv0_g)
!
    do i=it_start,it_end
!
      qc1  = q(i,j,1)
      qv1  = qv(i,j,1)	 
      pt1  = pt(i,j,1)
      ptv1 = pt1*(1. + epsm*qv1 - qc1)
      windm = usurf(i,j)
!
!-----------------------------------------------------------!
!              First case: Fixed surface fluxes             !
!-----------------------------------------------------------!
!
      if (isurf == 0) then
!
!  Initialise fluxes
!
	flux_t = shf / (pressurel%dens(i,j,1)*thermo_prop%cp(i,j,1))
        flux_q = lhf / (pressurel%dens(i,j,1)*cal_flv(pt1*exn0_g))
	flux_b = g*(flux_t/Pt0(1) + epsm*flux_q)
!
	if ( ust > 0. ) then
	  ustar = ust
	else
          ustar = vk*windm/lnz
        endif
!
	tstar = flux_t / ustar
	qvstar = flux_q / ustar
        ptvstar = tstar*(1. + epsm*qv1) + epsm*pt1*qvstar
!
!  Subiterations: find ustar, ptstar and qvstar
!
	if ( momsf > 0 .and. ust <= 0. ) then
	  k = 1
	  l0old = 1000.
	  delta = 1.
          do while ( k <= 6 )
            l0 = -(ustar**3.) / (vk*flux_b + small)
	    l0 = sign(max(abs(l0), small), l0)
	    l0 = max(min(l0,large), -large)
!
            zeta = zh / l0
            if (flux_b == 0.) then
	      phim = 1.
   	      phih = abs(pran)
!
              x = 0.
            else if (zeta > 0.) then
	      phim = 1. + betam*zeta
	      phih = abs(pran) + betah*zeta
!
	      x = 1. - phim
            else 
              phim = (1. - gamma1*zeta)**(-0.25)
              phih = abs(pran)*(1. - gamma2*zeta)**(-0.25)
!
              x = 2.*log(0.5*(1. + 1./phim)) + log(0.5*(1. + 1./(phim**2.))) - 2.*atan(1./phim) + 0.5*pi
            endif
	    delta = abs(l0 - l0old)
	    l0old = l0
	    k = k+1
!
            ustar  = vk*windm / (lnz - x)
	    tstar  = flux_t / ustar
	    qvstar = flux_q / ustar
!
       	    ptvstar = tstar*(1. + epsm*qv1) + epsm*pt1*qvstar
	  enddo
	endif
!
!-----------------------------------------------------------!
!           Second case: Fixed surface quantities           !
!    Caution: no function implemented for unstable case     !
!-----------------------------------------------------------!
!
      else if (isurf == 1) then
!
!  Initialize ustar, ptstar, qvstar
!
        ustar = vk*windm / lnz
  	tstar = -vk*(pt1 - t0_g) / lnz
	qvstar= -vk*(qv1 - qv0_g) / lnz
        ptvstar = tstar*(1. + epsm*qv1) + epsm*pt1*qvstar
!
!  Subiterations: find ustar, ptstar and qvstar
!
	k = 1
	l0old = 1000.
	delta = 1.
        do while ( k <= 6 )
          l0 = -ustar**2.*Pt0(1)/(vk*g*ptvstar + small)
	  l0 = sign(max(abs(l0), small), l0)
	  l0 = max(min(l0,large),-large)
!
!  Limit stability parameter to avoid excessive +/- fluxes
!
	  if (zh/l0 < -6.) then
	    zeff = -6.*l0
	  else if (zh/l0 > 12.) then
	    zeff = 12.*l0
	  else
	    zeff = zh
	  endif
!
          zeta = zeff / l0
	  if (ptvstar == 0.) then
	    phim = 1.
	    phih = abs(pran)  
!
          else if (zeta > 0.) then
	    phim = 1. + betam*zeta
	    phih = abs(pran) + betah*zeta
!
	    x = 1. - phim
	    ustar = vk*windm / (lnz  - x)
!
	    tstar  = -vk*(pt1 - t0_g) / (abs(pran)*lnz + betah*zeta)
	    qvstar = -vk*(qv1 - qv0_g) / (abs(pran)*lnz + betah*zeta)
!
          else if (zeta < 0.) then
            phim = (1. - gamma1*zeta)**(-0.25)
            phih = abs(pran) *(1. - gamma2*zeta)**(-0.25)
!
            x = 2.*log(0.5*(1. + 1./phim)) + log(0.5*(1. + 1./(phim**2.))) - 2.*atan(1./phim) + 0.5*pi
	    ustar = vk*windm / (lnz - x)
!
            tstar = -vk*(pt1 - t0_g) / (abs(pran)*(lnz - 2.*log(0.5*(1.+1./(phih**2.)))))
            qvstar = -vk*(qv1 - qv0_g) / (abs(pran)*(lnz - 2.*log(0.5*(1.+1./(phih**2.)))))
          endif
	  delta = abs(l0 - l0old)
	  l0old = l0
	  k = k+1
!
          ptvstar = tstar*(1. + epsm*qv1) + epsm*pt1*qvstar
	enddo
!
!-----------------------------------------------------------!
!              Third case: Fixed SST and drag	            !
!-----------------------------------------------------------!
!
      else if (isurf == 2) then
!
!  Near-neutral approximation with fixed drag coefficient
!
	ustar  = sqrt(c_dm)*windm
        tstar  = -c_ds/sqrt(c_dm)*(pt1 - t0_g) 
        qvstar = -c_ds/sqrt(c_dm)*(qv1 - qv0_g) 
        ptvstar = tstar*(1. + epsm*qv1) + epsm*pt1*qvstar
!
!-----------------------------------------------------------!
!         Fourth case: bulk according to Grabowski	    !
!-----------------------------------------------------------!
!
      else if (isurf == 3) then
!
!  Friction velocity is obtained from surface buoyancy
!
	ustar  = sqrt( 0.7*600.*g*((t0_g - pt1)/Pt0(1) + 0.61*(qv0_g - qv1)) )
	ustar  = max(min_w,sqrt(ustar**2. + windm**2.))
        tstar  = -c_ds*(pt1 - t0_g) 
        qvstar = -c_ds*(qv1 - qv0_g) 
        ptvstar = tstar*(1. + epsm*qv1) + epsm*pt1*qvstar
!
!-----------------------------------------------------------!
!         Fifth case: Use of flux-equivalent sst	    !
!-----------------------------------------------------------!
!
      else if (isurf == 4) then
!
!  Calculate sst from imposed flux and averaged first-layer properties
!
	flux_t = shf / (pressurel%dens(i,j,1)*thermo_prop%cp(i,j,1))
        flux_q = lhf / (pressurel%dens(i,j,1)*cal_flv(pt1*exn0_g))
        t0_g  = (flux_t + c_ds*ut0m) / (c_ds*u0m)
        qv0_g = (flux_q + c_ds*uq0m) / (c_ds*u0m)
!
	ustar  = max(min_w,sqrt(windm**2.))
        tstar  = -c_ds*(pt1 - t0_g) 
        qvstar = -c_ds*(qv1 - qv0_g) 
        ptvstar = tstar*(1. + epsm*qv1) + epsm*pt1*qvstar
!
      endif
!
!-----------------------------------------------------------!
!                   Assemble surface fluxes                 !
!-----------------------------------------------------------!
!
      surf%ustar(i,j) = ustar
      surf%mflux(i,j) = pressurel%dens(i,j,1)*ustar*ustar
      surf%qvflux(i,j) = pressurel%dens(i,j,1)*ustar*qvstar
      surf%esflux(i,j) = pressurel%dens(i,j,1)*ustar*tstar
!
#ifndef ISENTROPIC
      surf%esflux(i,j) = exn0_g*thermo_prop%cp(i,j,1)*surf%esflux(i,j) + flv00*surf%qvflux(i,j)
#endif
!
    end do
  end do
!
  call localflux ( surf%esflux, surf%qvflux )
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating calc_surface')
!
return
end
!
!  ===================================================
   subroutine localflux (esflux, qvflux)
!  ===================================================
!
  real, dimension(ip_start:,ip_end:), intent(inout) :: esflux, qvflux
!
  integer :: i, j
  real :: Lx, Ly, xx, yy, esg, qvg
  real, parameter :: sigma=2000.
!
!------------------------------------------------------------!
!
  select case (trim(casename))
!
    case ('SINGLE')
      call horav ( esflux, esg )
      call horav ( qvflux, qvg )
!
      Lx = real(nx-5)*dx/2.
      Ly = real(ny-5)*dy/2.
!
      do j = jt_start, jt_end
        do i = it_start, it_end
          xx = real(i-3)*dx
          yy = real(j-3)*dy
          esflux(i,j) = esflux(i,j) + 2.*esg*exp( -((xx-Lx)**2. + (yy-Ly)**2.) / sigma**2. )
          qvflux(i,j) = qvflux(i,j) + 2.*qvg*exp( -((xx-Lx)**2. + (yy-Ly)**2.) / sigma**2. )
        enddo
      enddo
!
  end select
!
return
end
!
end module surfacemod

