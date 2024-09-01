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
!  RADIAG.F90:
!
!  Purpose:
!      A front subroutine for preparing data and calling
!            Qiang Fu's radiation model
!
!  Author:
!      Chien Wang
!      MIT Joint Program for Science and Policy of Global Change
!
! ================================================================

  module radiationmod
  
  USE shared_all      
  USE shared_rad
  USE radia
  USE fuliou
  USE thermodynamics
  USE averages
#ifdef SPMD
  USE mpi
#endif
!
  IMPLICIT NONE
!
  private
!
  public :: rad_source, rad_front, rad_forc, init_rad

  CONTAINS
!
!      =============================================
      subroutine rad_source ( dens, ess )
!      =============================================
!
      integer :: i, j, k
      real, dimension(ip_start:,jp_start:,:), intent(in) :: dens
      real, dimension(ip_start:,jp_start:,:), intent(inout) :: ess
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dref
!
!---------------------------------------------------------------------!
!                    Calculate radiative source                       !
!---------------------------------------------------------------------!
!
#if (defined CONSERVATIVE)
  dref = dens
#else
  dref = 1.
#endif
!
  do k=1,nz
    do j=jp_start,jp_end
      do i=ip_start,ip_end
#ifdef ISENTROPIC
	ess(i,j,k) = ess(i,j,k) + dref(i,j,k)*rad%dtnet(i,j,k) / thermo%exn(i,j,k)
#else
    	ess(i,j,k) = ess(i,j,k) + dref(i,j,k)*thermo_prop%cp(i,j,k) * rad%dtnet(i,j,k)
#endif
      end do
    end do
  end do
!
#ifdef ISENTROPIC
  if (out_diagt) diag(7)%pt = diag(7)%pt + cdiag*dref*rad%dtnet / thermo%exn
#else
  if (out_diagt) diag(7)%pt = diag(7)%pt + cdiag*dref*rad%dtnet * thermo_prop%cp
#endif
!
!---------------------------------------------------------------------!
!
return
end subroutine rad_source

!      =============================================
      subroutine rad_forc 
!      =============================================
!
  integer  :: i, j, k, l, ki, k1, flag
  real     :: F1, F2, krad, frad1, frad2
  real     :: x, z, zi, cool
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: lwp, frad
!
  if (verbose > 0) call write_debug('Starting rad_forc')
!
!-----------------------------------------------------------------!
!
  rad%dtnet = 0.0
  frad = 0.0
!  
  select case (casename)
!
  case ('DYCOMS','ISDAC','SMOKE')
!
    if (casename == 'DYCOMS') then
      F1 = 22.
      F2 = 70.
      krad = 85.
    else if (casename == 'ISDAC') then
      F1 = 15.
      F2 = 72.
      krad = 170.    
    else if (casename == 'SMOKE') then
      F1 = 0.
      F2 = 60.
      krad = 0.02    
    endif
!
!  Radiative forcing for DYCOMS II 
!   
    do k = 1, nz
      lwp(:,:,k) = fdz0(k)*den0(k)*dz*hydromtr2(drop)%q(:,:,k)
    enddo
!
    do j = jt_start, jt_end
      do i = it_start, it_end
!
        ki = nz 
        flag = 0
        z = 0.5*dz*fdz0(1)
        zi = 500.0
        do k = 1, nz
          k1 = min(k+1,nz)
          if (state2%qt(i,j,k) >= 0.008) then
            ki = k
            zi = z
          endif
!    
          frad1 = 0.0
          frad2 = 0.0
          do l = 1, k
            frad1 = frad1 + lwp(i,j,l)
          enddo
          do l = k1, nz
            frad2 = frad2 + lwp(i,j,l)
          enddo
!
          frad2 = F2*exp(-krad*frad2)
          frad1 = F1*exp(-krad*frad1)
          frad(i,j,k) = frad1 + frad2
          if (k < nz) z = z + dz*0.5*(fdz0(k)+fdz0(k+1))
        enddo
!
        if (casename == 'DYCOMS') then
          z = zi
          do k = ki, nz
            frad(i,j,k) = frad(i,j,k) + 1.12*cp*Ddiv*(0.25*(z-zi)**(4./3.) + zi*(z-zi)**(1./3.))
            if (k < nz) z = z + dz*0.5*(fdz0(k)+fdz0(k+1))
          enddo
        endif
!
      enddo
    enddo
!      
   do k = 2, nz-1
      rad%dtnet(:,:,k) = -(frad(:,:,k+1) - frad(:,:,k-1)) / (dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1))) / (den0(k) * cp)
   enddo
   rad%dtnet(:,:,nz) = -(frad(:,:,nz) - frad(:,:,nz-1)) / (dz*(0.5*fdz0(nz)+0.5*fdz0(nz-1))) / (den0(nz) * cp)
   rad%dtnet(:,:,1) = -(frad(:,:,2) - frad(:,:,1)) / (dz*(0.5*fdz0(2)+0.5*fdz0(1))) / (den0(1) * cp)
!
!  Radiative convective equilibrium with fixed cooling
!
  case ('RCE1.5')
!
    cool = -1.5			!cooling rate, K/day
    cool = cool/86400.		!cooling rate, K/sec
!
    do k = 1, nz
      if (p0(k) > 40000.) then
        rad%dtnet(:,:,k) = cool
      else if (p0(k) <= 40000. .and. p0(k) > 20000.) then
        rad%dtnet(:,:,k) = cool * (p0(k) - 20000.) / 20000.
      else if (p0(k) <= 20000.) then
        rad%dtnet(:,:,k) = 0.
      endif
    enddo
!
  case ('RCE2','RCE2C')
!
    cool = -2.			!cooling rate, K/day
    cool = cool/86400.		!cooling rate, K/sec
!
    do k = 1, nz
      if (p0(k) > 40000.) then
        rad%dtnet(:,:,k) = cool * (p0(k)/pref)**((cp_a-cv_a)/cp_a)
      else if (p0(k) <= 40000. .and. p0(k) > 20000.) then
        rad%dtnet(:,:,k) = cool * (p0(k)/pref)**((cp_a-cv_a)/cp_a) * (p0(k) - 20000.) / 20000.
      else if (p0(k) <= 20000.) then
        rad%dtnet(:,:,k) = 0.
      endif
    enddo
!
  case ('RCE4','RCE4C')
!
    cool = -4.			!cooling rate, K/day
    cool = cool/86400.		!cooling rate, K/sec
!
    do k = 1, nz
      if (p0(k) > 40000.) then
        rad%dtnet(:,:,k) = cool * (p0(k)/pref)**((cp_a-cv_a)/cp_a)
      else if (p0(k) <= 40000. .and. p0(k) > 20000.) then
        rad%dtnet(:,:,k) = cool * (p0(k)/pref)**((cp_a-cv_a)/cp_a) * (p0(k) - 20000.) / 20000.
      else if (p0(k) <= 20000.) then
        rad%dtnet(:,:,k) = 0.
      endif
    enddo
!
  case ('RCE8','RCE8C')
!
    cool = -8.			!cooling rate, K/day
    cool = cool/86400.		!cooling rate, K/sec
!
    do k = 1, nz
      if (p0(k) > 40000.) then
        rad%dtnet(:,:,k) = cool * (p0(k)/pref)**((cp_a-cv_a)/cp_a)
      else if (p0(k) <= 40000. .and. p0(k) > 20000.) then
        rad%dtnet(:,:,k) = cool * (p0(k)/pref)**((cp_a-cv_a)/cp_a) * (p0(k) - 20000.) / 20000.
      else if (p0(k) <= 20000.) then
        rad%dtnet(:,:,k) = 0.
      endif
    enddo
!
  case ('RCE12','RCE12C')
!
    cool = -12.			!cooling rate, K/day
    cool = cool/86400.		!cooling rate, K/sec
!
    do k = 1, nz
      if (p0(k) > 40000.) then
        rad%dtnet(:,:,k) = cool * (p0(k)/pref)**((cp_a-cv_a)/cp_a)
      else if (p0(k) <= 40000. .and. p0(k) > 20000.) then
        rad%dtnet(:,:,k) = cool * (p0(k)/pref)**((cp_a-cv_a)/cp_a) * (p0(k) - 20000.) / 20000.
      else if (p0(k) <= 20000.) then
        rad%dtnet(:,:,k) = 0.
      endif
    enddo
!
  case ('RCEMIP','RCERH')
!
    do k = 1, nz
      if (z0(k) <= 10000.) then
        rad%dtnet(:,:,k) = -0.15e-4
      else if (z0(k) > 10000. .and. z0(k) <= 17500.) then
        rad%dtnet(:,:,k) = -0.15e-4 + (z0(k) - 10000.) * 0.3e-4/7500.
      else if (z0(k) > 17500. .and. z0(k) <= 22500.) then
        rad%dtnet(:,:,k) = 0.15e-4 - (z0(k) - 17500.) * 0.2e-4/5000.
      else if (z0(k) > 22500.) then
        rad%dtnet(:,:,k) = -0.05e-4
      endif
    enddo
!
!  Very idealized test case with fixed warming
!
  case ('ADJUST_1D')
!
    rad%dtnet = 0.
    x = (it_start-3)*dx
    do i = it_start, it_end
      if (x >= 0.5*dx*real(nx)-1000. .and. x <= 0.5*dx*real(nx)+1000.) rad%dtnet(i,:,:) = 0.1
      x = x+dx
    enddo
!
  case ('ADJUST_2D')
!
    rad%dtnet = 0.
    z = 0.
    do k = 1, nz
      if (z <= 3000. .and. z >= 1000.) rad%dtnet(:,:,k) = 0.001388889 * (p0(k)/pref)**((cp_a-cv_a)/cp_a)
      z=z+dz
    enddo
!
  case default
    rad%dtnet = 0.
!
  end select
!
! Clustered RCE cases
!
  if ( ANY( casename(1:6) == (/'RCE2C ','RCE4C ','RCE8C ','RCE12C'/) ) ) then
    do j = jp_start, jp_end
      do i = ip_start, ip_end
        if ( sqrt( real(i-nx/2)**2./real(nx/4)**2. + real(j-ny/2)**2./real(ny/4)**2. ) > 1. ) then
	  rad%dtnet(i,j,:) = rad%dtnet(i,j,:) / 2.
	endif
      enddo
    enddo
  endif
!
!  Store radiative flux in W/m2 per grid cell
!
  do k = 1, nz
    frad(:,:,k) = thermo_prop%cp(:,:,k)*pressure2%dens(:,:,k)*rad%dtnet(:,:,k)/dz1(k)
  enddo
!
!---------------------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating rad_forc')
!
  return
  end subroutine rad_forc
!
!  =============================================
  subroutine rad_front
!  =============================================
!
    use typedef_aerosol_new, only : nelem 
    use shared_data, only : nmode
    !USE mo_salsa_optical_properties, ONLY : initialize_optical_properties

      logical :: McICA = .FALSE.      
      real    :: hour, lat, radabove
      real    :: l, sc, u0, lamb, decl, epsi, func, sumh
      real, dimension(nz) :: fluxsu, fluxsd, fluxiru, fluxird
      real, dimension(ip_start:ip_end,jp_start:jp_end) :: twp, tolr, tssr, toar
      real, dimension(ip_start:ip_end,jp_start:jp_end) :: tosr, tisr !outgoing and incoming shortwave radiation
      real, dimension(ip_start:ip_end,jp_start:jp_end) :: tolr_cs, tosr_cs !outgoing longwave/shortwave radiation
!
      integer, save :: count, i, j, k, h, kk, ierr, iscal, im,ispec
      real, dimension(1:nz) :: pres, den

      real, dimension (nmode,nelem,nv) :: volc                !Volume concentration for each of the aerosols species in the mode
      real, dimension (nmode,nv) :: voltot                    !Total volume for each aerosol mode [] 
      real, dimension (nmode,nv) :: aero1d_num               !Aerosol number concentration in column ij 
      integer, parameter :: iba = 1                          ! wavenumber band containing shortest wavelenghts 
      real, parameter :: minSolarZenithCosForVisible = 1.e-3
!      
      real, dimension (nv) :: tau_aerosol,tau_cloud,tau_rain
      real, dimension (nv)   :: w_cloud,w_rain,w_aer,dzz1
      real, dimension (nv,4) :: ww_cloud, ww_rain, ww_aer

      if (verbose > 0) call write_debug('Start rad_front')
!
!---------------------------------------------------------------------!
!                    Define constants and  arrays                     !
!---------------------------------------------------------------------!
!
#ifdef RAD_ENABLE
!
!  Reinitialize rad
!
      rad%dtnet  = 0.0
      rad%dtsw   = 0.0
      rad%dtlw   = 0.0
      rad%frad   = 0.0
      rad%fluxs  = 0.0  
      rad%fluxir = 0.0

      twp  = 0.0
      tolr = 0.0
      tssr = 0.0 
      toar = 0.0
      tosr = 0.0
      tisr =0.0
      tolr_cs=0.0
      tosr_cs=0.0

      tau_aerosol = 0.0
      tau_cloud = 0.0
      tau_rain =0.0
!
      sc   = SolarConstant      	! Solar constant
      hour = time/3600. + t_local	! local hour
!
! --- Special cases
      select case ( trim(casename) )
      
      case ('RCEMIP','COLUMN','RCERH')
	sc = 551.58			! Halved solar constant
	u0 = cos(pi*42.05/180.)		! Fixed solar zenith angle of 42.05deg
      
      case default    
     	hour = 2.*pi*(hour - 12.)/24.				! Hour angle, local time (hour angle is 0 at noon)
        epsi = 23.4439*pi/180.0
        lamb = 2.*pi*real(j_day - 79.479522)/365.24
        decl = asin(sin(epsi) * sin(lamb))				! Sinus of declination angle
!
#ifdef CHANNEL
        l = abs(real(j-3) - 0.5*real(ny-5))*dx/111100.
#else
        l = 0.
#endif
        lat = (ctr_lat+l)/180.0*pi						! Calculate latitude
        u0 = sin(decl)*sin(lat) + cos(decl)*cos(lat)*cos(hour)		! Calculate solar zenith angle
      
      end select 
!
! --- Calculate thermodynamics properties      
      call get_cp_cv ( state2%qt, hydromtr2 )
!
!-----Initialize aerosol optical properties 
!-----(probably this should be moved to another module )
!      call initialize_optical_properties()

!  Start Loop over entire domain
!    

      do j = jt_start,jt_end      
        do i = it_start,it_end
	  fus = 0.0
	  fds = 0.0
	  fuir = 0.0
	  fdir = 0.0
	  taus = 0.0
	  tauir = 0.0
!
!---------------------------------------------------------------------!
!        Initialize atmospheric column for radiation model            !
!---------------------------------------------------------------------!
!
          do k = 1, nz
            kk = nv - k + 1
	    sumh = 0.
	    do h = 1, nhydro
	      sumh = sumh + hydromtr2(h)%q(i,j,k)
	    enddo
!
            pt(kk) = thermo%T(i,j,k)
            ph(kk) = max(state2%qt(i,j,k) - sumh,1.e-11)
	    po(kk) = 1.e-21
!
#ifdef ANELASTIC
            pres(k) = 0.01*p0(k)
	    den(k) = 1000.*den0(k)
#else
            pres(k) = 0.01*(pressure2%p(i,j,k) + p0(k))
	    den(k) = 1000.*pressure2%dens(i,j,k)
#endif
!
!  Ozone profile
!
#ifdef CHEM_ENABLE
	    if (rad_o3 == 1) then
              po(kk) = 48./mair * gas(i,j,k)%o3*1.e-9     !ppbv to kg/kg	
	    else
              po(kk) = 48./mair * gas0(k)%o3*1.e-9     !ppbv to kg/kg	
	    endif
#else
	    select case ( trim(casename) )
	      case ('RCEMIP','COLUMN')
	    	if (rad_o3 == 1) then
	          iscal = min(1,nscal)
	    	  po(kk) = 48./mair * state2%scal(i,j,k,iscal)*1.e-6      !ppmv to kg/kg	 
	        else if (rad_o3 == 0) then
	          iscal = min(1,nscal)
	    	  po(kk) = 48./mair * scal0(k,iscal)*1.e-6           !ppmv to kg/kg	
	    	endif
	    end select
#endif
	  enddo
!
!  Calculate hydrometeor mixing ratios for rad (in g/m3)
!
	  pp(nv1) = psurf * 0.01
	  plwc = 0.0
	  prwc = 0.0
	  piwc = 0.0
	  pgwc = 0.0
	  pre = 0.0
	  pde = 0.0
!
          do k = 1, nz
            kk = nv - (k-1)
            pp(kk) = pres(k)
!
            if (lmicro > 0) plwc(kk) = den(k)*hydromtr2(drop)%q(i,j,k)
            if (lmicro > 0) prwc(kk) = den(k)*hydromtr2(rain)%q(i,j,k)
            if (lmicro > 1) piwc(kk) = den(k)*hydromtr2(ice)%q(i,j,k)
            if (lmicro > 2) pgwc(kk) = den(k)*hydromtr2(grau)%q(i,j,k)
            if (lmicro > 3) pgwc(kk) = pgwc(kk) + den(k)*hydromtr2(hail)%q(i,j,k)
!
	    if (lmicro > 0) twp(i,j) = twp(i,j) + max(plwc(kk) + prwc(kk) + piwc(kk) + pgwc(kk),0.0)*dz*fdz0(k)
!
!  Calculate hydrometeor sizes
!
!  Cloud liquid
!
	    if (lmicro > 0) then
              if (plwc(kk) > 1.e-5 .and. hydromtr2(drop)%q(i,j,k) > qmin(drop) .and. hydromtr2(drop)%n(i,j,k) > xnmin(drop)) then
                pre(kk) = 1.e6*cal_reff (drop, hydromtr2(drop)%q(i,j,k), hydromtr2(drop)%n(i,j,k))
                pre(kk) = min(max(pre(kk),4.18),31.23)
              endif
            endif
!
!  Ice
!
#ifdef SEIFERT
	    if (lmicro > 1) then
              if (lfreeze > 0 .and. piwc(kk) > 1.e-5 .and. hydromtr2(ice)%q(i,j,k) > qmin(ice) .and. hydromtr2(ice)%n(i,j,k) > xnmin(ice)) then
                pde(kk) = 1.e6 * 2.*cal_reff (ice, hydromtr2(ice)%q(i,j,k), hydromtr2(ice)%n(i,j,k))
              endif
            endif
#endif
          enddo
!
#ifdef AERO_RADIA
!
voltot=0
volc=0 
aero1d_num=0

do k = 1, nz
  kk = nv - k + 1
  do im=1,nmode
    voltot(im,kk)=(aero3d(im)%m(i,j,k)/aero(im)%mix%rho) !* (100*pres(kk)/(R_air*pt(kk)))     !units [ m^3/m^3 ]

    aero1d_num(im,kk)=aero3d(im)%n(i,j,k) !*(100*pres(kk)/(R_air*pt(kk)))                     !units [ # /m^3 ]

    do ispec = 1, nelem  !number of aerosol species(same as nspec in SALSA):
      volc(im,ispec,kk)=voltot(im,kk) * aero(im)%init%frac(ispec)                            !units [ m^3 /m^3 ]
    enddo
  enddo
enddo
!
!---------------------------------------------------------------------!
!                Call main solver: infra red + visible                !
!     All the fluxes calculated by the radiation solver, that is      !
!    		     fd, fus, fdir, fuir, are positive		      !
!---------------------------------------------------------------------!
!
      if(out_ocs) then
!       !Computing clear sky fluxes
        call rad_ir( sst, emi, pp, pt, ph, po, fdir, fuir, &
                  volc=volc, voltot=voltot, aero1d_num=aero1d_num)  
!
        if (rad_sw == 1) then 
          call rad_vis( alb, u0, sc, pp, pt, ph, po, fds, fus, &
                          volc=volc, voltot=voltot, aero1d_num=aero1d_num, tau_aerosol=tau_aerosol)
        endif
!
        tolr_cs(i,j) = fuir(1)                           
        tosr_cs(i,j) = fus(1)
!
        !Reseting all fluxes to zero
        fdir=0.0
        fuir=0.0
        fds=0.0
        fus=0.0
      endif
! 
!--------------------------------------------------------------------
!
      call rad_ir( sst, emi, pp, pt, ph, po, fdir, fuir, &
                   volc=volc, voltot=voltot, aero1d_num=aero1d_num, &
		   plwc=plwc, pre=pre, piwc=piwc, pde=pde, 	 &
	   	   prwc=prwc, pgwc=pgwc) 
!
      if (rad_sw == 1) then 
        call rad_vis( alb, u0, sc, pp, pt, ph, po, fds, fus, &
                      volc=volc, voltot=voltot, aero1d_num=aero1d_num,	&
		      plwc=plwc, pre=pre, piwc=piwc, pde=pde, 	&
		      prwc=prwc, pgwc=pgwc , 	&
                      tau_aerosol=tau_aerosol,tau_cloud=tau_cloud, &
                      tau_rain=tau_rain)

        if(u0 > minSolarZenithCosForVisible) then   
            rad%tau_aer_surf(i,j) = SUM(tau_aerosol)
            rad%tau_cloud_surf(i,j) = SUM(tau_cloud)
            rad%tau_rain_surf(i,j) = SUM(tau_rain)
        endif
      endif
!
        if(u0 <= minSolarZenithCosForVisible) then
          call thicks(pp, pt, ph, dzz1)

          call cloud_water(iba, pre, plwc, dzz1, tau_cloud, w_cloud, ww_cloud)
          call rain_water( iba, prwc, dzz1, tau_rain, w_rain, ww_rain )
          call aero_rad(iba, volc,voltot,aero1d_num, dzz1, tau_aerosol, w_aer, ww_aer)

          rad%tau_aer_surf(i,j) = SUM(tau_aerosol)
          rad%tau_cloud_surf(i,j) = SUM(tau_cloud)
          rad%tau_rain_surf(i,j) = SUM(tau_rain)
      endif
!                                       
#else                                      
!
!---------------------------------------------------------------------!
!                Call main solver: infra red + visible                !
!     All the fluxes calculated by the radiation solver, that is      !
!    		     fd, fus, fdir, fuir, are positive		      !
!---------------------------------------------------------------------!
!
  if (out_ocs) then
    !Computing clear sky fluxes
      call rad_ir( sst, emi, pp, pt, ph, po, fdir, fuir)  
    !
    if (rad_sw == 1) call rad_vis( alb, u0, sc, pp, pt, ph, po, fds, fus)
    
    tolr_cs(i,j) = fuir(1)                           
    tosr_cs(i,j) = fus(1)

    !Reseting all fluxes to zero
    fdir=0.0
    fuir=0.0
    fds=0.0
    fus=0.0
  endif
!
!--------------------------------------------------------------------
!
  call rad_ir( sst, emi, pp, pt, ph, po, fdir, fuir, &
               plwc=plwc, pre=pre, piwc=piwc, pde=pde, 	 &
               prwc=prwc, pgwc=pgwc)  
!
  if (rad_sw == 1) then
    call rad_vis( alb, u0, sc, pp, pt, ph, po, fds, fus, &
              plwc=plwc, pre=pre, piwc=piwc, pde=pde, 	&
              prwc=prwc, pgwc=pgwc,tau_cloud=tau_cloud, &
              tau_rain=tau_rain ) 
!                        
    if(u0 > minSolarZenithCosForVisible) then   
      rad%tau_cloud_surf(i,j) = SUM(tau_cloud)
      rad%tau_rain_surf(i,j) = SUM(tau_rain)
    endif
  endif

  if(u0 <= minSolarZenithCosForVisible) then
    call thicks(pp, pt, ph, dzz1)

    call cloud_water(iba, pre, plwc, dzz1, tau_cloud, w_cloud, ww_cloud)
    call rain_water( iba, prwc, dzz1, tau_rain, w_rain, ww_rain )

    rad%tau_cloud_surf(i,j) = SUM(tau_cloud)
    rad%tau_rain_surf(i,j) = SUM(tau_rain)
  endif
!
#endif
!
!---------------------------------------------------------------------!
!                      Update fluxes and dtnet                        !
!---------------------------------------------------------------------!
!                                       
!  Fluxsu, Fluxsd:
!
      do k = 1, nz
        kk = nv - (k-1)
        fluxsu (k) = fus (kk)
        fluxsd (k) = fds (kk)
        fluxiru(k) = fuir(kk)
        fluxird(k) = fdir(kk)
      end do
      tolr(i,j) = fuir(1)
      tosr(i,j) = fus(1)
      tisr(i,j) = fds(1)
      toar(i,j) = fds(1) - fus(1) + fdir(1) - fuir(1)
      tssr(i,j) = fds(nv) - fus(nv) + fdir(nv) - fuir(nv)
!
!  Dtnet & frad net, Sign convention:
!	Fluxes should be considered positive when directed downward and
!	negative when they are directed upward. A - minus is thus necessary
!	in front of the upward fluxes (which have positive values when calculated 
!	by the radiation solver)
!
      do k = 1, nz
	rad%fluxs(i,j,k)  = fluxsd(k) - fluxsu(k)
	rad%fluxir(i,j,k) = fluxird(k) - fluxiru(k)
	rad%frad(i,j,k) = rad%fluxs(i,j,k) + rad%fluxir(i,j,k)

      end do
      radabove = fds(nv-nz) - fus(nv-nz) + fdir(nv-nz) - fuir(nv-nz)
!
      rad%dtnet(i,j,1) = ( rad%frad(i,j,2) - rad%frad(i,j,1) ) / ( 0.5*dz*(fdz0(1)+fdz0(2))*thermo_prop%cp(i,j,1)*den(1)/1000. )
      rad%dtsw(i,j,1) = ( rad%fluxs(i,j,2) - rad%fluxs(i,j,1) ) / ( 0.5*dz*(fdz0(1)+fdz0(2))*thermo_prop%cp(i,j,1)*den(1)/1000. )
      rad%dtlw(i,j,1) = ( rad%fluxir(i,j,2) - rad%fluxir(i,j,1) ) / ( 0.5*dz*(fdz0(1)+fdz0(2))*thermo_prop%cp(i,j,1)*den(1)/1000. )
!
      do k = 2, nz-1
        rad%dtnet(i,j,k) = ( rad%frad(i,j,k+1) - rad%frad(i,j,k-1) ) / ( dz*(0.5*fdz0(k-1)+fdz0(k)+0.5*fdz0(k+1))*thermo_prop%cp(i,j,k)*den(k)/1000. )
        rad%dtsw(i,j,k) = ( rad%fluxs(i,j,k+1) - rad%fluxs(i,j,k-1) ) / ( dz*(0.5*fdz0(k-1)+fdz0(k)+0.5*fdz0(k+1))*thermo_prop%cp(i,j,k)*den(k)/1000. )
        rad%dtlw(i,j,k) = ( rad%fluxir(i,j,k+1) - rad%fluxir(i,j,k-1) ) / ( dz*(0.5*fdz0(k-1)+fdz0(k)+0.5*fdz0(k+1))*thermo_prop%cp(i,j,k)*den(k)/1000. )
      enddo
      rad%dtnet(i,j,nz) = ( rad%frad(i,j,nz) - rad%frad(i,j,nz-1) ) / ( 0.5*dz*(fdz0(nz-1)+fdz0(nz))*thermo_prop%cp(i,j,nz)*den(nz)/1000. )
      rad%dtsw(i,j,nz) = ( rad%fluxs(i,j,nz) - rad%fluxs(i,j,nz-1) ) / ( 0.5*dz*(fdz0(nz-1)+fdz0(nz))*thermo_prop%cp(i,j,nz)*den(nz)/1000. )  
      rad%dtlw(i,j,nz) = ( rad%fluxir(i,j,nz) - rad%fluxir(i,j,nz-1) ) / ( 0.5*dz*(fdz0(nz-1)+fdz0(nz))*thermo_prop%cp(i,j,nz)*den(nz)/1000. )  
!
!  Smoothly decrease dtnet to 0 above zdec
!
      if ( iradx < 0. ) then
        do k = 1, nz
          func = max(1. - exp(4.*(z0(k) - zdec)/zdec),0.)
          rad%dtnet(i,j,k) = func * rad%dtnet(i,j,k)

          rad%dtsw(i,j,k) = func * rad%dtsw(i,j,k)
          rad%dtlw(i,j,k) = func * rad%dtlw(i,j,k)
	enddo
      endif
!  
!  End i, j loops
!
    end do
  end do
!
!  Store radiative flux in W/m2 per grid cell
!
  if (out_srad) then
    surf%lwf = rad%fluxir(:,:,1)
    surf%swf = rad%fluxs(:,:,1)
  endif
!  
  if (out_olr) surf%olr = tolr
  if (out_olr) surf%toa = toar
  if (out_osr) surf%osr = tosr
  if (out_osr) surf%isr = tisr
  if (out_ocs) surf%olr_cs= tolr_cs
  if (out_ocs) surf%osr_cs= tosr_cs

!
!  Calculate olr & ssr
!
  call horav (tssr, ssr)
  call horav (tolr, olr)
  call horav (toar, toa)
  call horav (tolr, olr_cs, xq=twp, xc=1.e-1, filter='lower')
  call horav (toar, toa_cs, xq=twp, xc=1.e-1, filter='lower')
!
#endif
!
!---------------------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate rad_front')
!
  return
  end subroutine rad_front
!
!  ====================================
  subroutine init_rad (nn, zp)
!  ====================================
    USE mo_salsa_optical_properties, ONLY : initialize_optical_properties

      integer, intent (in) :: nn
      real, intent (in) :: zp(nn)
!
#ifdef RAD_ENABLE
!
      real, allocatable  :: sp(:), st(:), sh(:), so(:), sl(:), si(:)
!
      integer :: k, ns, index
      logical :: blend
      real    :: pa, pb, ptest, test, dp1, dp2, dp3, Tsurf, tt
     
!
!  Reading standard atmosphere for radiation from sounding_rad.dat
!
!
      if (with_radmod .and. time == dt0 ) then
        write(7,*) 'Reading 1st Background Sounding for radiation: sounding_rad_time1.dat'
        
        open ( unit = 08, file = './sounding_rad_time1.dat', status = 'old' )
!
        read (08,*) Tsurf, ns
        allocate ( sp(ns), st(ns), sh(ns), so(ns), sl(ns), si(ns))
        do k=1,ns
          read (08, *) sp(k), st(k), sh(k), so(k), sl(k), si(k)
        enddo
!
        close (08)
      endif
      
      if (with_radmod .and. time == t_radmod) then 
        write(7,*) 'Reading 2nd Background Sounding for radiation: sounding_rad_time2.dat'
        deallocate (pp,fds,fus,fdir,fuir,taus,tauir)
        deallocate (pt,ph,po,pre,pde,plwc,prwc,piwc,pgwc)
        open ( unit = 09, file = './sounding_rad_time2.dat', status = 'old' )
!
        read (09,*) Tsurf, ns
        allocate ( sp(ns), st(ns), sh(ns), so(ns), sl(ns), si(ns))
        do k=1,ns
          read (09, *) sp(k), st(k), sh(k), so(k), sl(k), si(k)
        enddo
!
        close (09)

      endif

!
! identify what part, if any, of background sounding to use
!
      ptop = 0.01*zp(nn)
!
      if (sp(2) < ptop) then 
       pa = sp(1)
       pb = sp(2)
       k = 3
       do while (sp(k) < ptop)
          pa = pb
          pb = sp(k)
          k  = k+1
       end do
       k=k-1           ! identify first level above top of input
       blend = .True. 
      else
       blend = .False.
      endif
!
! if blend is true then the free atmosphere above the sounding will be
! specified based on the specified background climatology, here the 
! pressure levels for this part of the sounding are determined
!
      if (blend) then
        dp1 = pb-pa
        dp2 = ptop - pb
        dp3 = zp(nn-1) - zp(nn)
!        if (dp1 > 2.*dp2) k = k-1 ! first level is too close, blend from prev
        npts  = k
        norig = k
        ptest = sp(k)
        test = ptop-ptest
        do while (test > 2.*dp3)
          ptest = (ptest+ptop)*0.5
          test  = ptop-ptest
          npts  = npts + 1
        end do
        nv1 = npts + nz
      else
       nv1 = nz+1
       npts = nz
      endif
      nv = nv1-1
!     
      allocate (pp(nv1),fds(nv1),fus(nv1),fdir(nv1),fuir(nv1),taus(nv1),tauir(nv1))
      allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),prwc(nv),piwc(nv),pgwc(nv))
      pp=0.; pt=0.; ph=0.; po=0.; plwc=0.; prwc=0.; piwc=0.; pgwc=0.
!
!  fill in with reference atmosphere above domain top
!
      if (blend) then
       pp(1:norig) = sp(1:norig)
       pt(1:norig) = st(1:norig)
       ph(1:norig) = sh(1:norig)
       po(1:norig) = so(1:norig)
       plwc(1:norig) = sl(1:norig)
       piwc(1:norig) = si(1:norig)
       prwc(1:norig) = 0.
       pgwc(1:norig) = 0.
!
!  Interpolate in model domain
!       
       do k=norig+1,npts
          pp(k) = (ptop + pp(k-1))*0.5
          index = getindex(sp,ns,pp(k))
          pt(k) =  intrpl(sp(index),st(index),sp(index+1),st(index+1),pp(k))
          ph(k) =  intrpl(sp(index),sh(index),sp(index+1),sh(index+1),pp(k))
          po(k) =  intrpl(sp(index),so(index),sp(index+1),so(index+1),pp(k))
          plwc(k) = intrpl(sp(index),sl(index),sp(index+1),sl(index+1),pp(k))
          piwc(k) = intrpl(sp(index),si(index),sp(index+1),si(index+1),pp(k))
   	  prwc(k) = 0.
	  pgwc(k) = 0.
       end do
!       
!  Simple parameterization to determine ice crystal size in cirrus clouds (Lohmann and Karcher, 2002)
!
       pre(1:nv) = 0.0
       pde(1:nv) = 0.0
       do k = 1, npts
	 tt = pt(k) - 273.
	 if (piwc(k) > 0.) pde(k) = max(0.5*(326.3 + 12.4*tt + 0.2*tt*tt + 0.001*tt*tt*tt),20.)
       enddo
      end if
#ifdef ANELASTIC 
      piwc=0.
#else
      if (lfreeze == 0) piwc=0.
#endif
!
       do k = npts+1, nv
         po(k) = po(npts)
       enddo
!
       deallocate(sp, st, sh, so, sl, si)
!
#endif

#ifdef AERO_RADIA
!-----Initialize aerosol optical properties       
call initialize_optical_properties()

#endif

!
!----------------------------------------------------------!
!               Calculate radiative fluxes                 !
!----------------------------------------------------------!
!
       if ( ntime == 0 ) then
#ifdef RAD_ENABLE
         call rad_front
         call rad_front  
#else
         call rad_forc
#endif
       endif
!
      return
      
!  ===============================

      contains

      integer function getindex(x,n,xval)

      integer, intent (in)  :: n
      real,    intent (in)  :: x(n),xval

      integer :: ia, ib

      ia=1
      ib=n
      if (xval < x(1)) then
       getindex = 1
      elseif (xval > x(n)) then
       getindex = n-1
      else
       getindex = (ia+ib)/2
       do while (getindex /= ia .or. ib /= getindex+1)
          getindex = (ia+ib)/2
          if ((xval-x(getindex)) >= 0.0) then
             ia = getindex
          else
             ib = getindex
          end if
       end do
      endif

      end function getindex
!
!  ===============================
    real function intrpl(x1,y1,x2,y2,x)
 
    real, intent (in)  :: x1,y1,x2,y2,x

    real :: slope

    slope  = (y2-y1)/(x2 - x1 + epsilon(1.))
    intrpl = y1+slope*(x-x1)

    end function intrpl
!
    end subroutine init_rad
    
    end module radiationmod
