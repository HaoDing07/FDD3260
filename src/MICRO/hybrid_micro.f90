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

  module hybrid
  
  USE shared_data
  USE shared_hydro
  USE shared_pressure
  USE shared_thermo
  USE shared_surf
  USE shared_diag
  USE micropack
  
  IMPLICIT NONE
  
  ! General gamma distribution parameters
  real, parameter :: am(4)=(/523.5988,523.5988,52.36,171./)
  real, parameter :: bm(4)=(/3.,3.,3.,3.1/)
  real, parameter :: cap(4)=(/2.,2.,2.,2./)
  real, parameter :: gamp1(4)=(/24.,1.,1.,1./)			!gamma function value G(nu+1)
  real, parameter :: gampb(4)=(/5040.,6.,6.,6.81262/)	 	!gamma function value G(nu+bm+1)
  real :: nu(4)=(/4.,0.,0.,0./)
  logical :: diag_nu=.false.
  
  ! General parameters for the SCE solution
  integer :: nbin(4)
  real, parameter :: alpha=1.414, st=50.e-6
  real, parameter :: smin(4)=(/2.e-6,25.e-6,2.e-6,10.e-6/) 
  real, parameter :: smax(4)=(/100.e-6,5.e-3,1.e-3,1.e-2/) 
  real, parameter :: gmin(4)=(/1.e-3,1.e-3,1.e-3,1.e-3/)
  
  ! General options
  character(len=20) :: fallvel='best'
  character(len=20) :: kernel='hall'
  
  ! Other parameters
  real, parameter :: nmin=1.e-12

  ! g-distribution type
  type :: g_distribution
  	  real :: m0, m1
  	  real :: dm0, dm1
	  real :: vm0, vm1
  	  real :: dnauto, dmauto
          real, allocatable :: n(:)
  end type g_distribution
  
  private
  
  public :: micro_hybrid, calc_precip_hybrid
  
  contains
!
!----------------------------------------------------------!
  
!      =================================================                            
  subroutine micro_hybrid ( flag, cint, dref, qt, pressurel, hydromtrl, esm, hydromtrm, ww, dTdt )
!      =================================================                            

  integer, intent(in) :: flag
  real, intent(in) :: cint
  real, dimension(ip_start:,jp_start:,:), intent(in) :: dref, ww, dTdt, qt
  type (hydrometeor), dimension(1:nhydro), intent(in) :: hydromtrl
  type (atm_pressure), intent(in) :: pressurel
!
  type (hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrm
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: esm
!
  integer :: i, j, k, h
!
  real :: lambda, n0, dmax, qv
  real :: maxh, rho, qc, nc, qr, nr, qi, ni, qg, ng, sw, si, q, Tem, pp
  real :: auq, clr, aunc, aunr, clrn, crrn, auni, cig, cign
  real :: rcic, rcii, rrir, rrii, rcg, rrg, rcicn, rciin, rrirn, rriin, rcgn, rrgn
  real :: cdc, cdr, cdnc, cdnr, dei, deg, deni, deng, dqc, dqr, dqi, dqg, dqs, dqh, des, dnc, dnr, dni, dng
!
  type(g_distribution) :: gdistc, gdistr, gdisti, gdistg
!
!----------------------------------------------------------!
!                     Starting loops                       !
!----------------------------------------------------------!
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
!
        maxh = 0.
	do h = 1, rain
	  maxh = max(maxh,hydromtrl(h)%n(i,j,k))
	enddo
!
	if (maxh > 1.e-3) then
!
!----------------------------------------------------------!
!              	   Various initializations                 !
!----------------------------------------------------------!
!
	auq=0.; clr=0.; aunc=0.; aunr=0.; clrn=0.; crrn=0.; cig=0.; cign=0.
	rcic=0.; rcii=0.; rcicn=0.; rciin=0.; rrir=0.; rrii=0.; rrirn=0.; rriin=0.; rcg=0.; rcgn=0.; rrg=0.; rrgn=0.
	dqc=0.; dqr=0.; dqi=0.; dqg=0.; dqs=0.; dqh=0.; des=0.; dnc=0.; dnr=0.; dni=0.; dng=0.
	qc=0.; qr=0.; qi=0.; qg=0.
!
	rho = pressurel%dens(i,j,k)
!
	q = qt(i,j,k)
	if (lmicro > 0) qc = hydromtrl(drop)%q(i,j,k)
	if (lmicro > 0) qr = hydromtrl(rain)%q(i,j,k)
	if (lmicro > 1) qi = hydromtrl(ice)%q(i,j,k)
	if (lmicro > 2) qg = hydromtrl(grau)%q(i,j,k)
	if (lmicro > 0) nc = hydromtrl(drop)%n(i,j,k)
	if (lmicro > 0) nr = hydromtrl(rain)%n(i,j,k)
	if (lmicro > 1) ni = hydromtrl(ice)%n(i,j,k)
	if (lmicro > 2) ng = hydromtrl(grau)%n(i,j,k)
!
!  Define g-distributions
!
!  Cloud drops
!
        if( qc > qmin(drop) .and. nc > xnmin(drop) ) then
	  if ( diag_nu ) nu(drop) = 3. + min(exp(8.*nc/xn_ccn0),30.)
 	  lambda = (am(drop)*gampb(drop)*nc / (gamp1(drop)*qc))**(1./bm(drop))
  	  n0 = nc*lambda**(nu(drop)+1.) / gamp1(drop)
	  dmax = min(smax(drop), -log(gmin(drop)/n0) / lambda)
          nbin(drop) = ceiling(1. + bm(drop)*log(dmax/smin(drop)) / log(alpha))
!
          allocate ( gdistc%n(1:nbin(drop)) ) 
!	
  	  call get_gdistribution_froms ( drop, lambda, n0, gdistc )
	endif
!
!  Rain
!
        if( qr > qmin(rain) .and. nr > xnmin(rain) ) then
 	  lambda = (am(rain)*gampb(rain)*nr / (gamp1(rain)*qr))**(1./bm(rain))
  	  n0 = nr*lambda**(nu(rain)+1.) / gamp1(rain)
	  dmax = min(smax(rain), -log(gmin(rain)/n0) / lambda)
	  nbin(rain) = ceiling(1. + bm(rain)*log(dmax/smin(rain)) / log(alpha))
!
          allocate ( gdistr%n(1:nbin(rain)) ) 
!
	  call get_gdistribution_froms ( rain, lambda, n0, gdistr )
	endif
!
!  Ice
!
        if(lmicro > 1 .and. qi > qmin(ice) .and. ni > xnmin(ice) ) then
 	  lambda = (am(ice)*gampb(ice)*ni / (gamp1(ice)*qi))**(1./bm(ice))
  	  n0 = ni*lambda**(nu(ice)+1.) / gamp1(ice)
	  dmax = min(smax(ice), -log(gmin(ice)/n0) / lambda)
	  nbin(ice) = ceiling(1. + bm(ice)*log(dmax/smin(ice)) / log(alpha))
!
          allocate ( gdisti%n(1:nbin(ice)) ) 
!
	  call get_gdistribution_froms ( ice, lambda, n0, gdisti )
	endif
!
!  Graupel
!
        if(lmicro > 2 .and. qg > qmin(grau) .and. ng > xnmin(grau) ) then
 	  lambda = (am(grau)*gampb(grau)*ng / (gamp1(grau)*qg))**(1./bm(grau))
  	  n0 = ng*lambda**(nu(grau)+1.) / gamp1(grau)
	  dmax = min(smax(grau), -log(gmin(grau)/n0) / lambda)
	  nbin(grau) = ceiling(1. + bm(grau)*log(dmax/smin(grau)) / log(alpha))
!
          allocate ( gdistg%n(1:nbin(grau)) ) 
!
	  call get_gdistribution_froms ( grau, lambda, n0, gdistg )
	endif
!
!----------------------------------------------------------!
!       	      Warm Microphysics	                   !
!----------------------------------------------------------!
!
!  Rain auto-conversion and droplet self-collection
!
        if ( qc > qmin(drop) .and. nc > xnmin(drop) ) then
	  call self_collection ( drop, gdistc )
!
	  auq  = gdistc%dmauto
	  aunc = gdistc%dnauto + gdistc%dm0
	  aunr = gdistc%dnauto
	endif
!
!  Rain self-collection
!
        if ( qr > qmin(rain) .and. nr > xnmin(rain) ) then
	  call self_collection ( rain, gdistr )
!
	  crrn = gdistr%dm0
	endif
!
!  Droplet collection by rain
!
        if ( qc > qmin(drop) .and. nc > xnmin(drop) .and. qr > qmin(rain) .and. nr > xnmin(rain) ) then
	  call collection ( drop, rain, gdistc, gdistr )
!
	  clr  = gdistc%dm1
	  clrn = gdistc%dm0
	endif
!
!----------------------------------------------------------!
!       	      Ice Microphysics	                   !
!----------------------------------------------------------!
!
!  Ice only: self-collection
!
	if ( lmicro > 1 ) then
          if ( qi > qmin(ice) .and. ni > xnmin(ice) ) then
	    call self_collection ( ice, gdisti )
!
	    auni = gdisti%dm0
	  endif
	endif
!
!  If graupel is included
!
	if ( lmicro > 2 ) then
!
!  Droplet riming
!
          if ( qi > qmin(ice) .and. ni > xnmin(ice) .and. qc > qmin(drop) .and. nc > xnmin(drop) ) then
	    call collection_third ( drop, ice, grau, gdistc, gdisti, gdistg )
!
	    rcii  = gdisti%dm1
	    rcic  = gdistc%dm1
	    rciin = gdisti%dm0
	    rcicn = gdistc%dm0
	  endif
!
          if ( qg > qmin(grau) .and. ng > xnmin(grau) .and. qc > qmin(drop) .and. nc > xnmin(drop) ) then
	    call collection ( drop, grau, gdistc, gdistg )
!
	    rcg  = gdistc%dm1
	    rcgn = gdistc%dm0
	  endif
!
!  Rain riming
!
          if ( qi > qmin(ice) .and. ni > xnmin(ice) .and. qr > qmin(rain) .and. nr > xnmin(rain) ) then
	    call collection_third ( rain, ice, grau, gdistr, gdisti, gdistg )
!
	    rrii  = gdisti%dm1
	    rrir  = gdistr%dm1
	    rriin = gdisti%dm0
	    rrirn = gdistr%dm0
	  endif
!
          if ( qg > qmin(grau) .and. ng > xnmin(grau) .and. qr > qmin(rain) .and. nr > xnmin(rain) ) then
	    call collection ( rain, grau, gdistr, gdistg )
!
	    rrg  = gdistr%dm1
	    rrgn = gdistr%dm0
	  endif
!
!  Graupel-ice collisions
!
          if ( qg > qmin(grau) .and. ng > xnmin(grau) .and. qi > qmin(ice) .and. ni > xnmin(ice) ) then
	    call collection ( ice, grau, gdisti, gdistg )
!
	    cig  = gdisti%dm1
	    cign = gdisti%dm0
	  endif
!
	endif
!
!----------------------------------------------------------!
!               Condensation/evaporation                   !
!----------------------------------------------------------!
!
	Tem = thermo%T(i,j,k)
	qv  = qt(i,j,k)
	if (lmicro > 0) qv = qv - hydromtrl(drop)%q(i,j,k) - hydromtrl(rain)%q(i,j,k)
	if (lmicro > 1) qv = qv - hydromtrl(ice)%q(i,j,k)
	if (lmicro > 2) qv = qv - hydromtrl(grau)%q(i,j,k) - hydromtrl(snow)%q(i,j,k)
	if (lmicro > 3) qv = qv - hydromtrl(hail)%q(i,j,k)
	pp = pressure%p(i,j,k) + p0(k)
	sw  = qv - cal_qsw (Tem, pp)
	si  = qv - cal_qsi (Tem, pp)
!
!  Calculate evaporation/condensation rates
!
#ifndef SAT_ADJ
	if ( qc > qmin(drop) .and. nc > xnmin(drop) )then
          call cal_cond_gauss ( drop, den0(k), Tem, pp, gdistc, cdc, cdnc )
        endif
#endif
! 
        if ( qr > qmin(rain) .and. nr > xnmin(rain) ) then
          call cal_cond_gauss ( rain, den0(k), Tem, pp, gdistr, cdr, cdnr )
        endif
! 
        if ( lmicro > 1 .and. qi > qmin(ice) .and. ni > xnmin(ice) ) then
          call cal_cond_gauss ( ice, den0(k), Tem, pp, gdisti, dei, deni )
        endif
! 
        if ( lmicro > 2 .and. qg > qmin(grau) .and. ng > xnmin(grau) ) then
          call cal_cond_gauss ( grau, den0(k), Tem, pp, gdistg, deg, deng )
        endif
!
!  Complete tendencies (with resolved saturation)
!
        call resolved_saturation ( cdc, cdr, 0., dei, deg, 0., 0., pp, Tem, ww(i,j,k), dTdt(i,j,k), thermo_prop%cp(i,j,k), 		& 
		thermo_prop%cv(i,j,k), thermo%exn(i,j,k), sw, si, q, qc, qr, qi, qg, 0., 0., dqc, dqr, dqi, dqg, dqs, dqh, des )
!
	if (qc > qmin(drop) .and. lndrop == 1) dnc = min(max(dqc*nc/qc, -nc/dt0),0.)
	if (qr > qmin(rain)) dnr = min(max(dqr*nr/qr, -nr/dt0),0.)
	if (lmicro > 1 .and. qi > qmin(ice)) dni = min(max(dei*ni/qi, -ni/dt0),0.)
	if (lmicro > 2 .and. qg > qmin(grau)) dng = min(max(deg*ng/qg, -ng/dt0),0.)
!
!----------------------------------------------------------!
!             Update microphysics source terms             !
!----------------------------------------------------------!
!
        hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) - auq - clr + dqc
        hydromtrm(rain)%q(i,j,k) = hydromtrm(rain)%q(i,j,k) + auq + clr + dqr
!
        hydromtrm(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) - aunc - clrn + dnc
        hydromtrm(rain)%n(i,j,k) = hydromtrm(rain)%n(i,j,k) + aunr - crrn + dnr
!
	if (lmicro > 1) then
          hydromtrm(ice)%q(i,j,k) = hydromtrm(ice)%q(i,j,k) + dqi
          hydromtrm(ice)%n(i,j,k) = hydromtrm(ice)%n(i,j,k) - auni + dni
	endif
!
	if (lmicro > 2) then
          hydromtrm(drop)%q(i,j,k) = hydromtrm(drop)%q(i,j,k) - rcic - rcg
          hydromtrm(rain)%q(i,j,k) = hydromtrm(rain)%q(i,j,k) - rrir - rrg
          hydromtrm(ice)%q(i,j,k)  = hydromtrm(ice)%q(i,j,k)  - rcii - rrii - cig
          hydromtrm(grau)%q(i,j,k) = hydromtrm(grau)%q(i,j,k) + rcii + rcic + rrii + rrir + rcg + rrg + cig + dqg
!
          hydromtrm(drop)%n(i,j,k) = hydromtrm(drop)%n(i,j,k) - rcicn - rcgn
          hydromtrm(rain)%n(i,j,k) = hydromtrm(rain)%n(i,j,k) - rrirn - rrgn
          hydromtrm(ice)%n(i,j,k)  = hydromtrm(ice)%n(i,j,k)  - rciin - rriin - cign
          hydromtrm(grau)%n(i,j,k) = hydromtrm(grau)%n(i,j,k) + rcicn + rrirn + dng
	endif
!
#ifdef ISENTROPIC
	esm(i,j,k) = esm(i,j,k) + des
#endif
!
!  Diagnostics
!
	if (flag == 1 .and. ldiag) then
	  if (lmicro > 0) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + dref(i,j,k)*des
	  if (lmicro > 0) diag(4)%qc(i,j,k) = diag(4)%qc(i,j,k) + dref(i,j,k)*dqc
	  if (lmicro > 0) diag(4)%qr(i,j,k) = diag(4)%qr(i,j,k) + dref(i,j,k)*dqr
	  if (lmicro > 1) diag(4)%qi(i,j,k) = diag(4)%qi(i,j,k) + dref(i,j,k)*dqi
	else if (ldiag) then
	  if (lmicro > 0) diag(4)%pt(i,j,k) = diag(4)%pt(i,j,k) + (1.-cint)*dref(i,j,k)*des
	  if (lmicro > 0) diag(4)%qc(i,j,k) = diag(4)%qc(i,j,k) + (1.-cint)*dref(i,j,k)*dqc
	  if (lmicro > 0) diag(4)%qr(i,j,k) = diag(4)%qr(i,j,k) + (1.-cint)*dref(i,j,k)*dqr
	  if (lmicro > 1) diag(4)%qi(i,j,k) = diag(4)%qi(i,j,k) + (1.-cint)*dref(i,j,k)*dqi
	endif
!
!----------------------------------------------------------!
!                     End i j k loops                      !
!----------------------------------------------------------!
!
	if ( allocated(gdistc%n) ) deallocate (gdistc%n)
	if ( allocated(gdistr%n) ) deallocate (gdistr%n)
	if ( allocated(gdisti%n) ) deallocate (gdisti%n)
	if ( allocated(gdistg%n) ) deallocate (gdistg%n)
!
        endif
!
      end do
    end do
  end do
!
  if (lndrop == 0) hydromtrm(drop)%n = 0.
!
!----------------------------------------------------------!
!
  return
  end subroutine micro_hybrid
!
!
!  ==================================================
  subroutine calc_precip_hybrid ( dens, hydromtrl, qtm, hydromtrm )
!  ==================================================

    type(hydrometeor), dimension(1:nhydro), intent(in) :: hydromtrl
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: dens
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: qtm
    type(hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrm
!
    integer :: i, j, k
    real :: lambda, n0, dmax
    real :: dt1, qr, qrold, nr, nrold, pq, pn, qmax, qmini, nmax, nmini
    type(g_distribution) :: gdist_old, gdist
!
!----------------------------------------------------------!
!                   Initialize parameters                  !
!----------------------------------------------------------!
!
  dt1 = 1./dt0
  surf%precip = 0.0
!
!----------------------------------------------------------!
!                Estimating precipitation                  !
!----------------------------------------------------------!
!
  do j=jt_start,jt_end
    do i=it_start,it_end
!
      qrold = hydromtrl(rain)%q(i,j,nz)
      nrold = hydromtrl(rain)%n(i,j,nz)
      gdist_old%vm0 = 0.
      gdist_old%vm1 = 0.
!
      do k = nz-1, 1, -1
!	
        gdist%vm0 = 0.
	gdist%vm1 = 0.
	if (lmicro > 0) qr = hydromtrl(rain)%q(i,j,k)
	if (lmicro > 0) nr = hydromtrl(rain)%n(i,j,k)
!
        if ( qr > qmin(rain) .and. nr > xnmin(rain) ) then
!
!  Initialise size distribution
!
 	  lambda = (am(rain)*gampb(rain)*nr / (gamp1(rain)*qr))**(1./bm(rain))
  	  n0 = nr*lambda**(nu(rain)+1.) / gamp1(rain)
	  dmax = min(smax(rain), -log(gmin(rain)/n0)/lambda)
	  nbin(rain) = ceiling(1. + bm(rain)*log(dmax/smin(rain)) / log(alpha))
!
          allocate ( gdist%n(1:nbin(rain)) ) 
!
  	  call get_gdistribution_froms ( rain, lambda, n0, gdist )
!
!  Calculate precipitation velocities and fluxes
!
          call cal_vel_gauss( rain, den0(1)/den0(k), qr, nr, gdist )
!
	  pn = (dens(i,j,k+1)*nrold*gdist_old%vm0 - dens(i,j,k)*nr*gdist%vm0)*dz1(k)	 
	  pq = (dens(i,j,k+1)*qrold*gdist_old%vm1 - dens(i,j,k)*qr*gdist%vm1)*dz1(k)	 
  	  pn = pn / dens(i,j,k)
  	  pq = pq / dens(i,j,k)
!
!  Limit fluxes to avoid numerical instabilities
!
      	  qmax = max(qrold, qr)
          qmini = min(qrold, qr)
          pq = max(min(pq, (qr - qmini)*dt1), (qr - qmax)*dt1)
!
      	  nmax = max(nrold, nr)
          nmini = min(nrold, nr)
          pn = max(min(pn, (nr - nmini)*dt1), (nr - nmax)*dt1)
!
!  Update and diagnostics
!
	  hydromtrm(rain)%q(i,j,k) = hydromtrm(rain)%q(i,j,k) - pq
	  hydromtrm(rain)%n(i,j,k) = hydromtrm(rain)%n(i,j,k) - pn
  	  qtm(i,j,k)               = qtm(i,j,k) - pq
!
          if (k == 1) surf%precip(i,j) = -86400.*den0(1)*qr*gdist%vm1	! Precipitation fluxes mm/day
!
          qrold = qr
	  nrold = nr
	  gdist_old%vm0 = gdist%vm0	
	  gdist_old%vm1 = gdist%vm1
!
	  if ( allocated(gdist%n) ) deallocate (gdist%n)
!
	endif
!
      enddo
    end do
  end do    
!
!----------------------------------------------------------!
!  
  return
  end subroutine calc_precip_hybrid
!  
!      =================================================                            
  subroutine self_collection ( lhydro, gdist )
!      =================================================                            

  integer, intent(in) :: lhydro
  type(g_distribution), intent(inout) :: gdist
!
  integer :: i, j, k
  real :: dy, y0, yi, yj, mi, mj
  real :: dni, dnj, dna, dmi, dmj, dma, dnm0, dnm1, dnauto, dmauto
!
  y0 = log(smin(lhydro))
  dy = log(alpha) / bm(lhydro)
!
!  Considering collected drops i
!
  gdist%dm0 = 0.
  gdist%dm1 = 0.
  gdist%dnauto = 0.
  gdist%dmauto = 0.
!
  yi = y0
  do i = 1, nbin(lhydro)
!
    if ( gdist%n(i) > nmin ) then 
    mi = am(lhydro)*exp(bm(lhydro)*yi)
!
!  Collected by collector drops j 
!
    yj = y0 + real(i-1)*dy
    do j = i, nbin(lhydro)
      if ( gdist%n(j) > nmin .and. gdist%n(i) > nmin ) then 
        mj = am(lhydro)*exp(bm(lhydro)*yj)
!
        call gauss ( i, j, dy, dy, lhydro, lhydro, gdist, gdist, mi, mj, dni, dnj, dna, dmi, dmj, dma )
!
	gdist%dm0 = gdist%dm0 + dni
	gdist%dm1 = gdist%dm1 + dmi
	gdist%dnauto = gdist%dnauto + dna
	gdist%dmauto = gdist%dmauto + dma
      endif
      yj = yj + dy
    enddo
!
    endif
    yi = yi + dy
  enddo
!
!----------------------------------------------------------!
!  
  end subroutine self_collection
  
!      =================================================                            
  subroutine collection ( lhydroi, lhydroj, gdisti, gdistj )
!      =================================================                            

  integer, intent(in) :: lhydroi, lhydroj
  type(g_distribution), intent(inout) :: gdisti, gdistj
!
  integer :: i, j, k
  real :: dyi, dyj, y0i, y0j, yi, yj, mi, mj
  real :: dni, dnj, dna, dmi, dmj, dma
!
  y0i = log( smin(lhydroi) )
  y0j = log( smin(lhydroj) )
  dyi = log(alpha) / bm(lhydroi)
  dyj = log(alpha) / bm(lhydroj)
!
!  Considering collected drops i
!
  gdisti%dm0 = 0.
  gdisti%dm1 = 0.
  gdistj%dm0 = 0.
  gdistj%dm1 = 0.
!
  yi = y0i
  do i = 1, nbin(lhydroi)
!
    if ( gdisti%n(i) > nmin ) then 
    mi = am(lhydroi)*exp(bm(lhydroi)*yi)
!
!  Collected by collector drops j 
!
    yj = y0j + real(i-1)*dyj
    do j = i, nbin(lhydroj)
      mj = am(lhydroj)*exp(bm(lhydroj)*yj)
      if ( gdistj%n(j) > nmin .and. gdisti%n(i) > nmin .and. mj >= mi ) then 
!
        call gauss ( i, j, dyi, dyj, lhydroi, lhydroj, gdisti, gdistj, mi, mj, dni, dnj, dna, dmi, dmj, dma )
!
	gdisti%dm0 = gdisti%dm0 + dni
	gdisti%dm1 = gdisti%dm1 + dmi
	gdistj%dm1 = gdistj%dm1 + dmi
      endif
      yj = yj + dyj
    enddo
!
    endif
    yi = yi + dyi
  enddo
!
!----------------------------------------------------------!
!  
  end subroutine collection
  
!      =================================================                            
  subroutine collection_third ( lhydroi, lhydroj, lhydrok, gdisti, gdistj, gdistk )
!      =================================================                            

  integer, intent(in) :: lhydroi, lhydroj, lhydrok
  type(g_distribution), intent(inout) :: gdisti, gdistj, gdistk
!
  integer :: i, j, k
  real :: dyi, dyj, y0i, y0j, yi, yj, mi, mj
  real :: dni, dnj, dna, dmi, dmj, dma
!
  y0i = log( smin(lhydroi) )
  y0j = log( smin(lhydroj) )
  dyi = log(alpha) / bm(lhydroi)
  dyj = log(alpha) / bm(lhydroj)
!
!  Considering collected drops i
!
  gdisti%dm0 = 0.
  gdisti%dm1 = 0.
  gdistj%dm0 = 0.
  gdistj%dm1 = 0.
  gdistk%dm0 = 0.
  gdistk%dm1 = 0.
!
  yi = y0i
  do i = 1, nbin(lhydroi)
!
    if ( gdisti%n(i) > nmin ) then 
    mi = am(lhydroi)*exp(bm(lhydroi)*yi)
!
!  Collected by collector drops j 
!
    yj = y0j + real(i-1)*dyj
    do j = i, nbin(lhydroj)
      mj = am(lhydroj)*exp(bm(lhydroj)*yj)
      if ( gdistj%n(j) > nmin .and. gdisti%n(i) > nmin .and. mj >= mi ) then 
!
        call gauss ( i, j, dyi, dyj, lhydroi, lhydroj, gdisti, gdistj, mi, mj, dni, dnj, dna, dmi, dmj, dma )
!
	gdisti%dm0 = gdisti%dm0 + dni
	gdistj%dm0 = gdistj%dm0 + dnj
	gdistk%dm0 = gdistk%dm0 + dnj
!
	gdisti%dm1 = gdisti%dm1 + dmi
	gdistj%dm1 = gdistj%dm1 + dmj
	gdistk%dm1 = gdistj%dm1 + dmi + dmj
      endif
      yj = yj + dyj
    enddo
!
    endif
    yi = yi + dyi
  enddo
!
!----------------------------------------------------------!
!  
  end subroutine collection_third

!      =================================================                            
  subroutine gauss ( i, j, dyi, dyj, lhydroi, lhydroj, gdisti, gdistj, mi, mj, dni, dnj, dna, dmi, dmj, dma )
!      =================================================                            

  integer, intent(in) :: i, j, lhydroi, lhydroj
  real, intent(in) :: mi, mj, dyi, dyj
  type(g_distribution), intent(in) :: gdisti, gdistj
  real, intent(out) :: dni, dnj, dna, dmi, dmj, dma
!
  integer :: k
  real :: cmj, cpj, cj, cmi, cpi, ci
  real :: mt, mj1, mj3, mi1, mi3, dni1, dni2, dni3, dnj1, dnj2, dnj3
  real :: F11, F12, F13, F21, F22, F23, F31, F32, F33
  real, parameter :: y1=-sqrt(3./5.), y2=0., y3=sqrt(3./5.)
  real, parameter :: w1=5./18., w2=8./18., w3=5./18.
!
!  Estimate coefficients for piecewise linear approximations
!
  if ( j == 1 ) then
    cmj = 0.5*(gdistj%n(j+1) - gdistj%n(j))
    cpj = 0.5*(gdistj%n(j+1) - gdistj%n(j))
    cj  = gdistj%n(j)
  else if ( j == size(gdistj%n) ) then
    cmj = 0.5*(gdistj%n(j) - gdistj%n(j-1))
    cpj = 0.5*(gdistj%n(j) - gdistj%n(j-1))
    cj  = gdistj%n(j)
  else
    cmj = 0.5*(gdistj%n(j) - gdistj%n(j-1))
    cpj = 0.5*(gdistj%n(j+1) - gdistj%n(j))
    cj  = gdistj%n(j)
  endif
!
  if ( i == 1 ) then
    cmi = 0.5*(gdisti%n(i+1) - gdisti%n(i))
    cpi = 0.5*(gdisti%n(i+1) - gdisti%n(i))
    ci  = gdisti%n(i)
  else if ( i == size(gdisti%n) ) then
    cmi = 0.5*(gdisti%n(i) - gdisti%n(i-1))
    cpi = 0.5*(gdisti%n(i) - gdisti%n(i-1))
    ci  = gdisti%n(i)
  else
    cmi = 0.5*(gdisti%n(i) - gdisti%n(i-1))
    cpi = 0.5*(gdisti%n(i+1) - gdisti%n(i))
    ci  = gdisti%n(i)
  endif
!
!  Point locations in m coordinate
!
  mj1 = mj*alpha**(0.5*y1)
  mj3 = mj*alpha**(0.5*y3) 
!
  mi1 = mi*alpha**(0.5*y1)
  mi3 = mi*alpha**(0.5*y3) 
!
!  Estimate function value at 3 locations
!
  F11 = ker(lhydroi, mi1, lhydroj, mj1) * (cmj*y1 + cj) * (cmi*y1 + ci)
  F12 = ker(lhydroi, mi1, lhydroj, mj ) * gdistj%n(j)   * (cmi*y1 + ci)
  F13 = ker(lhydroi, mi1, lhydroj, mj3) * (cpj*y3 + cj) * (cmi*y1 + ci)
!
  F21 = ker(lhydroi, mi, lhydroj, mj1) * (cmj*y1 + cj) * gdisti%n(i)
  F22 = ker(lhydroi, mi, lhydroj, mj ) * gdistj%n(j)   * gdisti%n(i)
  F23 = ker(lhydroi, mi, lhydroj, mj3) * (cpj*y3 + cj) * gdisti%n(i)
!
  F31 = ker(lhydroi, mi3, lhydroj, mj1) * (cmj*y1 + cj) * (cpi*y3 + ci)
  F32 = ker(lhydroi, mi3, lhydroj, mj ) * gdistj%n(j)   * (cpi*y3 + ci)
  F33 = ker(lhydroi, mi3, lhydroj, mj3) * (cpj*y3 + cj) * (cpi*y3 + ci)
!
!  Gauss quadrature for number
!
  dni1 = (w1*F11/mj1 + w2*F12/mj + w3*F13/mj3)*dyj
  dni2 = (w1*F21/mj1 + w2*F22/mj + w3*F23/mj3)*dyj
  dni3 = (w1*F31/mj1 + w2*F32/mj + w3*F33/mj3)*dyj
!
  dni = (w1*dni1/mi1 + w2*dni2/mi + w3*dni3/mi3)*dyi
!
  dnj1 = (w1*F11/mi1 + w2*F21/mi + w3*F31/mi3)*dyi
  dnj2 = (w1*F12/mi1 + w2*F22/mi + w3*F32/mi3)*dyi
  dnj3 = (w1*F13/mi1 + w2*F23/mi + w3*F33/mi3)*dyi
!
  dnj = (w1*dnj1/mj1 + w2*dnj2/mj + w3*dnj3/mj3)*dyj
!
!  Gauss quadrature for mass
!
  dmi = (w1*dni1 + w2*dni2 + w3*dni3)*dyi
!
  dmj = (w1*dnj1 + w2*dnj2 + w3*dnj3)*dyj
!
!  Auto-conversion
!
  mt = am(lhydroj)*exp(bm(lhydroj)*log(st))
!
  dna = 0.
  dma = 0.
  if ( mi+mj >= mt ) then
    dna = dnj
    dma = dmi + dmj
  endif
!
!----------------------------------------------------------!
!  
  end subroutine gauss

!      =================================================                            
  subroutine cal_vel_gauss ( lhydro, den_ratio, q, n, gdist )
!      =================================================                            

  integer, intent(in) :: lhydro
  real, intent(in) :: den_ratio, q, n
  type(g_distribution), intent(inout) :: gdist
!
  integer :: i
  logical :: yes
  real :: cmi, cpi, ci
  real :: yi, dy, F1, F2, F3, m, m1, m3, d, d1, d3, vn, vm
  real, parameter :: y1=-sqrt(3./5.), y2=0., y3=sqrt(3./5.)
  real, parameter :: w1=5./18., w2=8./18., w3=5./18.
!
  gdist%vm0 = 0.
  gdist%vm1 = 0.
!
  yes = .false.
  yi = log(smin(lhydro))
  dy = log(alpha) / bm(lhydro)
  do i = 1, nbin(lhydro)
    if ( gdist%n(i) > nmin ) then 
!
!  Estimate coefficients for piecewise linear approximations
!
    if ( i == 1 ) then
      cmi = 0.5*(gdist%n(i+1) - gdist%n(i))
      cpi = 0.5*(gdist%n(i+1) - gdist%n(i))
      ci  = gdist%n(i)
    else if ( i == size(gdist%n) ) then
      cmi = 0.5*(gdist%n(i) - gdist%n(i-1))
      cpi = 0.5*(gdist%n(i) - gdist%n(i-1))
      ci  = gdist%n(i)
    else
      cmi = 0.5*(gdist%n(i) - gdist%n(i-1))
      cpi = 0.5*(gdist%n(i+1) - gdist%n(i))
      ci  = gdist%n(i)
    endif
!
!  Point locations in m and d coordinate
!
    m = am(lhydro)*exp(bm(lhydro)*yi)
    m1 = m*alpha**(0.5*y1)
    m3 = m*alpha**(0.5*y3) 
!
    d  = (m /am(lhydro))**(1./bm(lhydro))
    d1 = (m1/am(lhydro))**(1./bm(lhydro))
    d3 = (m3/am(lhydro))**(1./bm(lhydro))
!
!  Estimate function value at 3 locations
!
    F1 = fall_speed ( lhydro, d1 ) * (cmi*y1 + ci)
    F2 = fall_speed ( lhydro, d  ) * gdist%n(i)
    F3 = fall_speed ( lhydro, d3 ) * (cpi*y3 + ci)
!
!  Gauss quadrature for number
!
    vn = (w1*F1/m1 + w2*F2/m + w3*F3/m3)*dy
    vm = (w1*F1    + w2*F2   + w3*F3   )*dy
!
    gdist%vm0 = gdist%vm0 + vn
    gdist%vm1 = gdist%vm1 + vm
!
    yes = .true.
    endif
!
    yi = yi + dy
  enddo
!
!  Rescale velocities
!
  if ( yes) then
    gdist%vm0 = -gdist%vm0 / n * sqrt(den_ratio)
    gdist%vm1 = -gdist%vm1 / q * sqrt(den_ratio)
  endif
!
!----------------------------------------------------------!
!  
  end subroutine cal_vel_gauss
  
!      =================================================                            
  real function ker ( lhydro1, m1, lhydro2, m2 )
!      =================================================                            

  integer, intent(in) :: lhydro1, lhydro2
  real, intent(in) :: m1, m2
  real :: d1, d2
!  
  select case (trim(kernel))
!
    case ('hall')
      d1 = (m1/am(lhydro1))**(1./bm(lhydro1))
      d2 = (m2/am(lhydro2))**(1./bm(lhydro2))
!
      ker = hall_ker ( lhydro1, d1, m1, lhydro2, d2, m2 )
!
    case('golovin')
      ker = golovin_ker ( m1, m2 )
!
    case('long')
      d1 = (m1/am(lhydro1))**(1./bm(lhydro1))
      d2 = (m2/am(lhydro2))**(1./bm(lhydro2))
!
      ker = long_ker ( d1, d2, m1, m2 )
!
  end select
!
  contains
!
!----------------------------------------------------------!
!  
  real function hall_ker ( lhydro1, d1, m1, lhydro2, d2, m2 )

  integer, intent(in) :: lhydro1, lhydro2
  real, intent(in) :: d1, d2, m1, m2
!
  integer :: k
  real :: v1, v2, kc, e
  real :: e1(20)=(/0.0001,0.0001,0.0001,0.014,0.017,0.019,0.022,0.027,0.03,0.033,0.035,0.037,0.038,0.038,0.037,0.036,0.035,0.032,0.029,0.027/)
  real :: e2(20)=(/0.0001,0.0001,0.005,0.016,0.022,0.03,0.043,0.052,0.064,0.072,0.079,0.082,0.08,0.076,0.067,0.057,0.048,0.04,0.033,0.027/)
  real :: e3(20)=(/0.0001,0.002,0.02,0.04,0.085,0.17,0.27,0.4,0.5,0.55,0.58,0.59,0.58,0.54,0.51,0.49,0.47,0.45,0.47,0.52/)
  real :: e4(20)=(/0.001,0.07,0.28,0.5,0.62,0.68,0.74,0.78,0.8,0.8,0.8,0.78,0.77,0.76,0.77,0.77,0.78,0.79,0.95,1.4/)
  real :: e5(20)=(/0.005,0.4,0.6,0.7,0.78,0.83,0.86,0.88,0.9,0.9,0.9,0.9,0.89,0.88,0.88,0.89,0.92,1.01,1.3,2.3/)
  real :: e6(20)=(/0.05,0.43,0.64,0.77,0.84,0.87,0.89,0.9,0.91,0.91,0.91,0.91,0.91,0.92,0.93,0.95,1.,1.03,1.7,3./)
  real :: e7(20)=(/0.2,0.58,0.75,0.84,0.88,0.9,0.92,0.94,0.95,0.95,0.95,0.95,0.96,0.96,0.97,1.,1.02,1.04,2.3,4./)
  real :: e8(20)=(/0.5,0.79,0.91,0.95,0.97,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
  real :: e9(20)=(/0.77,0.93,0.97,0.98,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
  real :: e10(20)=(/0.87,0.96,0.98,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
  real :: e11(20)=(/0.97,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
!
  v1 = fall_speed ( lhydro1, d1 )
  v2 = fall_speed ( lhydro2, d2 )
!
  if ( lhydro1 == ice .and. (lhydro2 == ice .or. lhydro1 == grau) ) then
    kc = 0.5*(m1*m2)/(m1+m2)*(v2 - v1)**2.
    e = exp(-20654.*kc/pi/d1**2.)
  else 
    k = int( 20.*max(min(d1/d2,1.),0.05) )
    if ( d2 < 40.e-6 ) then
      e = e1(k)
    else if ( d2 >= 40.e-6 .and. d2 < 60.e-6 ) then
      e = e2(k)  
    else if ( d2 >= 60.e-6 .and. d2 < 80.e-6 ) then
      e = e3(k)  
    else if ( d2 >= 80.e-6 .and. d2 < 100.e-6 ) then
      e = e4(k)  
    else if ( d2 >= 100.e-6 .and. d2 < 120.e-6 ) then
      e = e5(k)  
    else if ( d2 >= 120.e-6 .and. d2 < 140.e-6 ) then
      e = e6(k)  
    else if ( d2 >= 140.e-6 .and. d2 < 200.e-6 ) then
      e = e7(k)  
    else if ( d2 >= 200.e-6 .and. d2 < 300.e-6 ) then
      e = e8(k)  
    else if ( d2 >= 300.e-6 .and. d2 < 400.e-6 ) then
      e = e9(k)  
    else if ( d2 >= 400.e-6 .and. d2 < 600.e-6 ) then
      e = e10(k)  
    else
      e = e11(k)  
    endif
  endif
!
  hall_ker = 0.25*pi*e*(d1 + d2)**2.*abs(v2 - v1)
!
  end function hall_ker
!
!----------------------------------------------------------!
!  
  real function golovin_ker ( m1, m2 )

  real, intent(in) :: m1, m2
  real, parameter :: b=1.5
!
  golovin_ker = b*(m1 + m2)
!
  end function golovin_ker
!
!----------------------------------------------------------!
!  
  real function long_ker ( d1, d2, m1, m2 )

  real, intent(in) :: d1, d2, m1, m2
  real, parameter :: k1=9.44e9, k2=5.78
!
  if ( max(d1,d2) <= 100.e-6 ) then
    long_ker = k1*(m1**2. + m2**2.)
  else
    long_ker = k2*(m1 + m2)
  endif
!
  end function long_ker
!
!----------------------------------------------------------!
!  
  end function ker
 
!      =================================================                            
  real function fall_speed ( lhydro, d )
!      =================================================                            

  integer, intent(in) :: lhydro
  real, intent(in) :: d
!
  select case (lhydro)
!
  case (drop,rain)
    select case (trim(fallvel))
!
!  Atlas et al. 1973
!
    case ('atlas73')
      fall_speed = max(9.65 - 10.3*exp(-600.*d),0.)
!
!  Atlas and Ulbrich 1977
!
    case ('au77')
      fall_speed = 386.58*d**(2./3.)
!
!  Best 1950
!
    case default
      fall_speed = 9.58*(1. - exp(-(100.*d/0.171)**1.147))
    end select
!
  case (ice)
    fall_speed = 11.72*d**(0.41)
!
  case (grau)
    fall_speed = 199.05*d**(0.8)
!
  end select
!
!----------------------------------------------------------!
!  
  end function fall_speed
!
!  ==================================================
   subroutine cal_cond_gauss ( lhydro, den, Tem, pp, gdist, cd, cdn )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Condensation of cloud drops, no ventilation factor
!
! ================================================================
!
  integer, intent(in)  :: lhydro
  type(g_distribution), intent(in) :: gdist
  real, intent(in)     :: Tem, pp, den
  real, intent(out)    :: cd, cdn
!
  integer :: i
  real :: xfkd, xfmu, Fv, Fv1, Fv3, Re, Sc, ff, es, qs, dq
  real :: yi, dy, cmi, cpi, ci
  real :: m, m1, m3, d, d1, d3, F1, F2, F3
  real, parameter :: y1=-sqrt(3./5.), y2=0., y3=sqrt(3./5.)
  real, parameter :: w1=5./18., w2=8./18., w3=5./18.
!
!----------------------------------------------------------!
!
  cd  = 0.
  cdn = 0.
!
!  Global parameters
!
  xfkd = cal_xfkd (Tem, pp)
  xfmu = cal_xfmu (Tem)
! 
  if (lhydro <= 2) then  
    es = cal_esw (Tem)
    qs = cal_qsw (Tem, pp)
    ff = 2.*pi / (qs * ((Rw*Tem)/(es*xfkd) + cal_flv(Tem)/(Kt*Tem)*(cal_flv(Tem)/(Rw*Tem) - 1.)))
  else
    es = cal_esi (Tem)
    qs = cal_qsi (Tem, pp)
    ff = 2.*pi / cap(lhydro) / (qs * ((Rw*Tem)/(es*xfkd) + cal_fls(Tem)/(Kt*Tem)*(cal_fls(Tem)/(Rw*Tem) - 1.)))
  endif
!
  Sc = xfmu / diff
!
  yi = log(smin(lhydro))
  dy = log(alpha) / bm(lhydro)
!
  do i = 1, nbin(lhydro)
    if ( gdist%n(i) > nmin ) then 
!
!  Estimate coefficients for piecewise linear approximations
!
    if ( i == 1 ) then
      cmi = 0.5*(gdist%n(i+1) - gdist%n(i))
      cpi = 0.5*(gdist%n(i+1) - gdist%n(i))
      ci  = gdist%n(i)
    else if ( i == size(gdist%n) ) then
      cmi = 0.5*(gdist%n(i) - gdist%n(i-1))
      cpi = 0.5*(gdist%n(i) - gdist%n(i-1))
      ci  = gdist%n(i)
    else
      cmi = 0.5*(gdist%n(i) - gdist%n(i-1))
      cpi = 0.5*(gdist%n(i+1) - gdist%n(i))
      ci  = gdist%n(i)
    endif  
!
!  Point locations in m and d coordinate
!
    m = am(lhydro)*exp(bm(lhydro)*yi)
    m1 = m*alpha**(0.5*y1)
    m3 = m*alpha**(0.5*y3) 
!
    d1 = (m1/am(lhydro))**(1./bm(lhydro))
    d  = (m /am(lhydro))**(1./bm(lhydro))
    d3 = (m3/am(lhydro))**(1./bm(lhydro))
!
!  Ventilation factor
!
    Re = fall_speed( lhydro, d ) * sqrt(den/den0(1)) * d / xfmu 
    if ( Sc**(1./3.)*sqrt(Re) < 1.4 ) then
      Fv1 = 1. + 0.108 * Sc**(2./3.) * ( fall_speed(lhydro,d1)*d1/xfmu )
      Fv  = 1. + 0.108 * Sc**(2./3.) * Re
      Fv3 = 1. + 0.108 * Sc**(2./3.) * ( fall_speed(lhydro,d3)*d3/xfmu )
    else
      Fv1 = 0.78 + 0.308 * Sc**(1./3.) * sqrt( fall_speed(lhydro,d1)*d1/xfmu )
      Fv  = 0.78 + 0.308 * Sc**(1./3.) * sqrt(Re)
      Fv3 = 0.78 + 0.308 * Sc**(1./3.) * sqrt( fall_speed(lhydro,d3)*d3/xfmu )
    endif
!
!  Local condensation rates including ventilation
!
    F1 = ff * d1 * Fv1 * (cmi*y1 + ci)
    F2 = ff * d *  Fv  * gdist%n(i)
    F3 = ff * d3 * Fv3 * (cpi*y3 + ci)
!
!  Gauss quadrature for condensation
!
    dq = (w1*F1/m1 + w2*F2/m + w3*F3/m3)*dy
    cd  = cd  + dq
    cdn = cdn + dq/m
!
    endif
!
    yi = yi + dy
  enddo
!
return
end
!
!      =================================================                            
  subroutine get_gdistribution_froms ( lhydro, lambda, n0, gdist )
!      =================================================                            

  integer, intent(in) :: lhydro
  real, intent(in) :: lambda, n0
  type(g_distribution), intent(inout) ::gdist
!
  integer :: k
  real :: m
!
  m = am(lhydro)*smin(lhydro)**bm(lhydro)
  do k = 1, nbin(lhydro)
    gdist%n(k) = bm(lhydro)*m*m*dmass( nu(lhydro), am(lhydro), bm(lhydro), n0, lambda, m )
    m = alpha*m
  enddo
!
!----------------------------------------------------------!
!  
  end subroutine get_gdistribution_froms
!
!      =================================================                            
  real function dmass ( nu, am, bm, n0, lambda, m )
!      =================================================                            

  real, intent(in) :: nu, am, bm, lambda, n0, m
  real :: lambdam, m0m, num, mum
!
!  Convert size distribution parameters to mass distribution parameters
!
  lambdam = lambda * am**(-1./bm)
  m0m = n0 * am**(-(nu + 1.)/bm) / bm
  num = (nu + 1.)/bm - 1.
  mum = 1. / bm
!
!  Estimate mass distribution at mass m
!
  dmass = m0m * m**num * exp(-lambdam*m**mum)
  
  end function dmass
  
  end module hybrid
