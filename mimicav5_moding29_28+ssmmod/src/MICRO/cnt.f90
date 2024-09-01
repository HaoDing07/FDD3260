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
!  LOOKUP
!	New routine for Classical Nucleation Theory
!
!  Author:
!	Julien Savre, MISU
!
! ==============================================================
!
#ifdef NUC_CNT
!
!      =====================================================
  subroutine cnt_init
!      =====================================================

      USE gridno
      USE shared_data
      USE shared_nuclei

      IMPLICIT NONE 
      
      integer :: h, l
!
! --- Initialize INs for detailed ice nucleation
!
  lmode(:,1) = ldep
  lmode(:,2) = limm
  lmode(:,3) = lcon
!
!  Contact angle
!
  do h = 1, 3
    do l = 1, 3
      if (nucprop(h)%theta(l) > 180.) lmode(h,l) = .FALSE.
    enddo
    if (nucprop(h)%frac == 0. .and. nucprop(h)%fracc == 0.) lmode(h,:) = .FALSE.
  enddo
!
return
end
!
!  ==================================================
   subroutine freez_cnt ( i, j, k, nucin_tmp, qc, nc, ni, Tem, Sw, Si, e, fci, fcni, fcd, fcnd )
!
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Freezing of liquid drops --> immersion and contact
!
! ================================================================
!
USE shared_data
USE shared_aerosol_new
USE shared_thermo
USE shared_nuclei
!
IMPLICIT NONE
!
  integer  :: i, j, k, h, j, i, l
  real     :: qc, nc, ni, nin, dni
  real     :: fci, fcni, fcd, fcnd
  real     :: frac, fracc, fracm, fracmc
  real     :: Tem, Tc, Sw, Si, e
  real, dimension(3)  :: nn, nnc, mm, mmc
  real, parameter     :: wmin = 4.42e-14
  type(nuclei_3d), dimension(3) :: nucin_tmp
!
  real, dimension(3,2) :: ncont, phif, phifm
  real, dimension(3)   :: dd, ddc, theta00
!
!----------------------------------------------------------!
!
!  Detailed ice nucleation process accounting for homogeneous and heterogeneous nucleation.
!  foch, foci and focc (homogeneous, immersion and contact freezing).
!
  Tc   = Tem - T00
  dd   = 0.
  ddc  = 0.
  fcd  = 0.
  fcnd = 0.
  fci  = 0.
  fcni = 0.
!
!----------------------------------------------------------!
!               Calculate nucleation rates                 !
!----------------------------------------------------------!
!
  do l = 1, 2
    if ( lmode(1,l).or.lmode(2,l).or.lmode(3,l) ) then 	
    phif = 0.0
    phifm = 0.0
!
!----------------------------------------------------------!
!           Total nucleation rates: deposition             !
!----------------------------------------------------------!
!
  if (l == 1 .and. Sw < 1. .and. Si > 1.) then
!
    do h = 1, 3
      nn(h)  = max(nucin_tmp(h)%mode(1)%n(i,j,k) - nucin_tmp(h)%mode(2)%n(i,j,k) - nucin_tmp(h)%mode(3)%n(i,j,k), 0.)
      nnc(h) = max(nucin_tmp(h)%mode(1)%nc(i,j,k) - nucin_tmp(h)%mode(2)%nc(i,j,k) - nucin_tmp(h)%mode(3)%nc(i,j,k), 0.)
#ifdef NUC_CNT1
      mm(h)  = max(nucin_tmp(h)%mode(1)%m(i,j,k)- nucin_tmp(h)%mode(2)%m(i,j,k) - nucin_tmp(h)%mode(3)%m(i,j,k), 0.)
      mmc(h) = max(nucin_tmp(h)%mode(1)%mc(i,j,k) - nucin_tmp(h)%mode(2)%mc(i,j,k) - nucin_tmp(h)%mode(3)%mc(i,j,k), 0.)
#endif
    enddo
!
!  Initialize sizes
!
    dd = 0.
    ddc = 0.
    do h = 1, 3
#ifdef NUC_CNT1
      if ( nn(h) > xin_min .and. mm(h) > qin_min ) dd(h) = sqrt( cal_logmom( 2., nucp(h)%s, mm(h), nn(h), nucp(h)%rhop ) )
      if ( nnc(h) > xin_min .and. mmc(h) > qin_min ) ddc(h) = sqrt( cal_logmom( 2., nucp(h)%sc, mmc(h), nnc(h), nucp(h)%rhop ) )
#else
      if ( nn(h) > xin_min ) dd(h) = nucp(h)%d * exp(2.*log(nucp(h)%s)*log(nucp(h)%s))
      if ( nnc(h) > xin_min ) ddc(h) = nucp(h)%dc * exp(2.*log(nucp(h)%sc)*log(nucp(h)%sc))
#endif
    enddo
!
!  Rates
!
    call calc_rates ( i, j, k, l, nucin_tmp, nucprop(:)%theta(l), nucprop(:)%sigma(l),   &
    		      dd, ddc, Tc, Sw, Si, e, theta00, phif, phifm )
!
!  Calculate rates
!
    do h = 1, 3
      frac  = phif(h,1)
      fracm = phifm(h,1) 
      fracc  = phif(h,2) 
      fracmc = phifm(h,2)
      dni = (nn(h)*frac + nnc(h)*fracc) / dt0
!
      nucin_tmp(h)%mode(3)%n(i,j,k)  = nucin_tmp(h)%mode(3)%n(i,j,k)  + nn(h) * frac
      nucin_tmp(h)%mode(3)%nc(i,j,k) = nucin_tmp(h)%mode(3)%nc(i,j,k) + nnc(h) * fracc
#ifdef NUC_CNT1
      nucin_tmp(h)%mode(3)%m(i,j,k)  = nucin_tmp(h)%mode(3)%m(i,j,k)  + mm(h) * fracm
      nucin_tmp(h)%mode(3)%mc(i,j,k) = nucin_tmp(h)%mode(3)%mc(i,j,k) + mmc(h) * fracmc 
#endif
!
      fcnd = fcnd + dni
    enddo
    fcnd = max(min(fcnd,(sum(nn+nnc) - ni)/dt0),0.)
    fcd  = fcnd * wmin
!
  endif
!
!----------------------------------------------------------!
!           Total nucleation rates: immersion              !
!----------------------------------------------------------!
!
  if (l == 2 .and. Sw > 1. .and. qc > qmin(drop) .and. nc > xnmin(drop)) then	
!
    do h = 1, 3
      nn(h)  = nucin_tmp(h)%mode(2)%n(i,j,k) 
      nnc(h) = nucin_tmp(h)%mode(2)%nc(i,j,k)
#ifdef NUC_CNT1
      mm(h)  = nucin_tmp(h)%mode(2)%m(i,j,k) 
      mmc(h) = nucin_tmp(h)%mode(2)%mc(i,j,k)
#endif
    enddo
!
!  Initialize sizes
!
    dd = 0.
    ddc = 0.
    do h = 1, 3
#ifdef NUC_CNT1
      if ( nn(h) > xin_min .and. mm(h) > qin_min ) dd(h) = sqrt( cal_logmom( 2., nucp(h)%s, mm(h), nn(h), nucp(h)%rhop ) )
      if ( nnc(h) > xin_min .and. mmc(h) > qin_min ) ddc(h) = sqrt( cal_logmom( 2., nucp(h)%sc, mmc(h), nnc(h), nucp(h)%rhop ) )
#else
      if ( nn(h) > xin_min ) dd(h) = nucp(h)%d * exp(2.*log(nucp(h)%s)*log(nucp(h)%s))
      if ( nnc(h) > xin_min ) ddc(h) = nucp(h)%dc * exp(2.*log(nucp(h)%sc)*log(nucp(h)%sc))
#endif
    enddo
!
!  Rates
!
    call calc_rates ( i, j, k, l, nucin_tmp, nucprop(:)%theta(l), nucprop(:)%sigma(l),   &
       	  	      dd, ddc, Tc, Sw, Si, e, theta00, phif, phifm )
!
!  Calculate rates
!
    do h = 1, 3
      if ( nn(h) + nnc(h) > xin_min ) then
        frac  = phif(h,1)
        fracm = phifm(h,1) 
        fracc  = phif(h,2) 
        fracmc = phifm(h,2)
        dni = (nn(h)*frac + nnc(h)*fracc) / dt0
      
!		
        nucin_tmp(h)%mode(3)%n(i,j,k)  = nucin_tmp(h)%mode(3)%n(i,j,k)  + nucin_tmp(h)%mode(2)%n(i,j,k) * frac
        nucin_tmp(h)%mode(3)%nc(i,j,k) = nucin_tmp(h)%mode(3)%nc(i,j,k) + nucin_tmp(h)%mode(2)%nc(i,j,k) * fracc
        nucin_tmp(h)%mode(2)%n(i,j,k)  = nucin_tmp(h)%mode(2)%n(i,j,k)  - nucin_tmp(h)%mode(2)%n(i,j,k) * frac
        nucin_tmp(h)%mode(2)%nc(i,j,k) = nucin_tmp(h)%mode(2)%nc(i,j,k) - nucin_tmp(h)%mode(2)%nc(i,j,k) * fracc
!
#ifdef NUC_CNT1
        nucin_tmp(h)%mode(3)%m(i,j,k)  = nucin_tmp(h)%mode(3)%m(i,j,k)  + nucin_tmp(h)%mode(2)%m(i,j,k) * fracm
        nucin_tmp(h)%mode(3)%mc(i,j,k) = nucin_tmp(h)%mode(3)%mc(i,j,k) + nucin_tmp(h)%mode(2)%mc(i,j,k) * fracmc 
        nucin_tmp(h)%mode(2)%m(i,j,k)  = nucin_tmp(h)%mode(2)%m(i,j,k)  - nucin_tmp(h)%mode(2)%m(i,j,k) * fracm
        nucin_tmp(h)%mode(2)%mc(i,j,k) = nucin_tmp(h)%mode(2)%mc(i,j,k) - nucin_tmp(h)%mode(2)%mc(i,j,k) * fracmc 
#endif
!
        fcni = fcni + dni
      endif
    enddo
    fcni = max(min(fcni,(sum(nn+nnc) - ni)/dt0),0.)
    fci  = fcni * qc / nc
!
  endif
!
!----------------------------------------------------------!
!        		  Ending       	                   !
!----------------------------------------------------------!
!
    endif
  enddo
!
return
end
!
!  ==================================================
   subroutine calc_rates ( i, j, k, flag, nucin_tmp, thetam, thetas, dd, ddc, Tc, Sw, Si, e, theta00, eff, effm )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Freezing of cloud drops
!
! ================================================================
!
USE shared_data
USE shared_nuclei
USE shared_thermo
!
IMPLICIT NONE
!
  integer  :: i, j, k, h, flag
!
  real     :: Tc, Sw, Si, e
  real     :: thetalim, dm, sm, dplus, dminus
  real     :: ff, ffc, ffm, ffmc
  real, dimension(2) 	:: err, prob, freez
  real, dimension(3)    :: dd, ddc, nn, nnc
  real, dimension(3)    :: thetam, thetas, theta00, theta00c
  real, dimension(3,2) 	:: eff, effm
  type(nuclei_3d), dimension(3) :: nucin_tmp
  external :: focd, foci
!
!----------------------------------------------------------!
!
!  Detailed ice nucleation process accounting for homogeneous and heterogeneous nucleation.
!  foch, foci and focc (homogeneous, immersion and contact freezing).
! 
  theta00 = 0.
  theta00c = 0.
  if ( theta_low ) then		
    do h = 1, 3
      if ( flag == 1 ) then
        nintot  = nucin_tmp(h)%mode(1)%n(i,j,k)  - nucin_tmp(h)%mode(2)%n(i,j,k)
        nintotc = nucin_tmp(h)%mode(1)%nc(i,j,k) - nucin_tmp(h)%mode(2)%nc(i,j,k)
      else 
        nintot  = nucin_tmp(h)%mode(2)%n(i,j,k) + nucin_tmp(h)%mode(3)%n(i,j,k)
        nintotc = nucin_tmp(h)%mode(2)%nc(i,j,k) + nucin_tmp(h)%mode(3)%nc(i,j,k)	  
      endif
!
      if ( nintot > xnmin(ice) .and. nucin_tmp(h)%mode(3)%n(i,j,k) > xnmin(ice) ) then
        frac = min(max(nucin_tmp(h)%mode(3)%n(i,j,k)/nintot,0.000001),0.999999)
        call normal_cdf_inv ( frac,  thetam(h), thetas(h), theta00(h) )
        theta00(h) = min(max(theta00(h),0.),180.)
      endif
!
      if  ( nintotc > xnmin(ice) .and. nucin_tmp(h)%mode(3)%nc(i,j,k) > xnmin(ice) ) then
       frac = min(max(nucin_tmp(h)%mode(3)%nc(i,j,k)/nintotc,0.000001),0.999999)
        call normal_cdf_inv ( frac, thetam(h), thetas(h), theta00c(h) )
        theta00c(h) = min(max(theta00c(h),0.),180.)
      endif
    enddo
  endif
!
  eff   = 0.
  effm  = 0.
!
!  Integrate using gauss-legendre quadrature
!
  do h = 1, 3
    if ( lmode(h,flag) ) then
!
!  Fine IN mode
!
      thetalim = theta00(h)
!      
      if ( flag == 1 .and. dd(h) > 0.0 ) call gauss_legendre ( focd, thetam(h), thetas(h), theta00(h), Tc, e, Si, dd(h), nucp(h)%s, ff, ffm )
!      
      if ( flag >= 2 .and. dd(h) > 0.0 ) call gauss_legendre ( foci, thetam(h), thetas(h), theta00(h), Tc, e, Si, dd(h), nucp(h)%s, ff, ffm )
!
      eff(h,1) = ff
      effm(h,1) = ffm
!
!  Coarse IN mode
!
      thetalim = theta00c(h)
!      
      if ( flag == 1 .and. ddc(h) > 0.0 ) call gauss_legendre ( focd, thetam(h), thetas(h), theta00c(h), Tc, e, Si, ddc(h), nucp(h)%sc, ffc, ffmc )
!      
      if ( flag >= 2 .and. ddc(h) > 0.0 ) call gauss_legendre ( foci, thetam(h), thetas(h), theta00c(h), Tc, e, Si, ddc(h), nucp(h)%sc, ffc, ffmc )
!      
      eff(h,2) = ffc
      effm(h,2) = ffmc
    endif
  enddo	
!
return
end
!
!  ==================================================
   subroutine foci ( theta, Tc, e, si, d, sigma, J, Jm )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Heterogeneous freezing of cloud droplets in the immersion mode:
!         -  Theoretically based parameterization: CNT extended to heterogeneous nucleation
!            & all the temperature dependent parameters are taken from PK97
!
! ================================================================
!
USE shared_data
USE shared_nuclei
USE shared_thermo
!
IMPLICIT NONE
!
  real     :: d, theta, sigma
  real     :: si, e, Tc, Tem, sigmaiw, dfcr, rcr, zeldo, dgact
  real     :: J0, J, Jm, a, fac, f, m
!
!----------------------------------------------------------!
!
!  Nucleation in the immersion mode should take the same form as homogeneous nucleation
!  except that the nucleation rate is calculated as a surface rate multiplied by the immersed
!  aerosol surface area. 
!  Assumptions for nucleation rates: immersion freezing occurs upon unsoluble fraction of activated aerosols
!
!----------------------------------------------------------!
!     Calculate critical free energy of germ formation     !
!         in immersion mode & Zeldovitch factor		   !
!----------------------------------------------------------!
!
  Tem   = Tc + 273.15
  sigmaiw = cal_sigmaiw (Tem)
  dgact = cal_dgact (Tem)
!
!  Critical radius and Gibbs energy
!
  rcr  = 2.*sigmaiw*(mwm/(rhoi*kb*Tem*log(Si)))
  dfcr  = 4./3.*pi*rcr*rcr*cal_sigmaiw(Tem)
!
!  Shape parameter
!
  m 	= theta * pi/180.
  f 	= 0.25 * (2. + cos(m)) * (1. - cos(m))**2.
!
  a = log(2.*nmol) + log(kbhp*Tem) + log(mwm/rhoi) - 2.*log(rcr) + 0.5*log(sigmaiw*f/(kb*Tem))
!
!----------------------------------------------------------!
!              Freezing probability of dust                !
!----------------------------------------------------------!
!
  J0    = - (dgact + dfcr*f) / (kb*Tem)
  J0    = exp(J0 + a)
!
!----------------------------------------------------------!
!            Integrate over size distributions             !
!----------------------------------------------------------!
!
  J  = pi * d**2. * J0
!
  Jm = pi * d**2. * exp(6.*log(sigma)*log(sigma)) * J0
!
return
end
!
!  ==================================================
   subroutine focd ( theta, Tc, e, si, d, sigma, J, Jm )
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Heterogeneous freezing of cloud droplets in the deposition mode:
!         -  Theoretically based parameterization: CNT extended to heterogeneous nucleation
!            & all the temperature dependent parameters are taken from PK97
!
! ================================================================
!
USE shared_data
USE shared_nuclei
USE shared_thermo
!
IMPLICIT NONE
!
  real     :: d, sigma, rhop, theta
  real     :: si, e, Tc, Tem, sigmaiv, dfcr, dfads, rcr
  real     :: J00, J0, J, Jm, frac, f, m
!
!----------------------------------------------------------!
!
!  Nucleation in the deposition mode is defined as in Khvorostyanov & Curry 2000.
!  Newly nucleated crystals are given the minimum allowed mass wmin.
!
!----------------------------------------------------------!
!            Calculate critical free energy                !
!         of germ formation in deposition mode		   !
!----------------------------------------------------------!
!
  Tem   = Tc + 273.15
  sigmaiv = cal_sigmaiv (Tem)
!
  rcr   = 4.*sigmaiv / (rhoi*Rw*Tem*log(si))
  dfcr  = 1./3.*pi*sigmaiv * rcr**2. 
  J00   = e*e / (4.*pi*rhow*kb*Tem*nuv) * sqrt(sigmaiv/(kb*Tem))	! Includes Zeldovitch factor
!
!----------------------------------------------------------!
!               Freezing probability of dust               !
!----------------------------------------------------------!
!
  m 	= theta * pi/180.							! Shape factor
  f 	= 0.25 * (2. + cos(m)) * (1. - cos(m))**2.
  dfads = cal_dfads (3, Tem)
  J0    = J00 / sqrt(f) * exp(- (dfads + dfcr*f) / (kb*Tem))
!
!----------------------------------------------------------!
!            Integrate over size distributions             !
!	          Fine and coarse modes			   !
!----------------------------------------------------------!
!
  J  = pi * d**2. * J0
!
  Jm = pi * d**2. * exp(6.*log(sigma)*log(sigma)) * J0
!
  return
  end
!
!  ==================================================
  subroutine gauss_legendre ( func, thetam, thetas, thetalim, Tc, e, S, dd, sigma, ff, ffm )
!  ==================================================
!
  USE shared_data
  USE shared_nuclei
!
  IMPLICIT NONE 
  
  real :: Tc, e, S, dd, thetam, thetas, thetalim, sigma, ff, ffm
  external :: func
  
  integer :: i
  real  :: J, Jm, ftheta, thet, thetaplus, thetaminus

  ! Points and weights for 10th order method
  real, dimension(10), parameter  :: x = (/ -0.973908239224106,-0.865060566454542,-0.679410688379483,-0.433395397408358,-0.148874339007507, &
        				    0.148874339007507,0.433395397408358,0.679410688379483,0.865060566454542,0.973908239224106 /)
  real, dimension(10), parameter  :: w = (/ 0.066667031481273176,0.14945422671662059,0.21908574324601154,0.26926671836765848,0.29552422471242434, &
					    0.29552422471242434,0.26926671836765848,0.21908574324601154,0.14945422671662059,0.066667031481273176 /)
!
!----------------------------------------------------------!
!
    thetaplus = 180.
    thetaminus = thetalim
!    
    ff = 0.
    ffm = 0.
    do i = 1, 10
!    
      thet = 0.5*(thetaplus - thetaminus)*x(i) + 0.5*(thetaplus + thetaminus) 
!    
      call func ( thet, Tc, e, S, dd, sigma, J, Jm )
!    
      if ( theta_low ) then
        if (thet > thetalim) then
          ftheta = cal_npdf (thetam, thetas, thet)
        else
          ftheta = 0.
        endif
      else 
        ftheta = cal_npdf (thetam, thetas, thet)
      endif
!      
      ff = ff + w(i)*J*ftheta
!      
      ffm = ffm + w(i)*Jm*ftheta
!      
    enddo
!
    ff = 0.5*(thetaplus - thetaminus)*ff
!
    ffm = 0.5*(thetaplus - thetaminus)*ffm
!
    ff = 1. - exp(- ff * dt0)
!
    ffm = 1. - exp(- ffm * dt0)
!
    if ( ff < 1.e-12 ) ff = 0.
    if ( ffm < 1.e-12 ) ffm = 0.
!
!----------------------------------------------------------!
!  
  return
  end
!
#endif
