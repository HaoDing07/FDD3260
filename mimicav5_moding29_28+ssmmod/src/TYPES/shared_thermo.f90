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
!  SHARED_THERMO.F:                   
!
!  Purpose:
!	Module containing diagnostic thermodynamics variables		  
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================

	module shared_thermo

	use gridno
	use shared_data
	use typedef_thermo
	use shared_aerosol_new
	
	IMPLICIT NONE
		
	SAVE
	
	type (thermodyn) :: thermo, thermo_mean
	type (thermoprop) :: thermo_prop

!
!  Definition of thermo parameters
!
	real, parameter :: cp_a  = 1004.		! Heat capacity dry air
	real, parameter :: cv_a  = 716. 		! Heat capacity dry air
	real, parameter :: cp_v  = 1827.		! Heat capacity water vapour
	real, parameter :: cv_v  = 1388.5		! Heat capacity water vapour
	real, parameter :: cp_l  = 4185.		! Heat capacity liquid water
	real, parameter :: cv_l  = 4185.		! Heat capacity liquid water
	real, parameter :: cp_i  = 2102.		! Heat capacity ice water
	real, parameter :: cv_i  = 2102.		! Heat capacity ice water
	real, parameter :: epsm  = 0.60243		! Coefficient for virtual temperature
	real, parameter :: eps = 0.624052		! (cp_a-cv_a)/(cp_v-cv_v)
	real, parameter :: delta_aw = 0.305		! Theoretical water activity delta (koop et al.)
	
	real, parameter :: rho_l = 1000.		! Liquid water density
	real, parameter :: rho_i = 920.			! Bulk ice density

	real, parameter :: fls00 = 2.838889e6		! Latent heat of sublimation at 0K
	real, parameter :: flv00 = 2.495800e6		! Latent heat of vaporization at 0K
	real, parameter :: fiv = 2.839e6		! Latent heat of sublimation
	real, parameter :: flv = 2.501e6		! Latent heat of vaporization
	real, parameter :: flm = 338000.		! Latent heat of melting

        real, parameter :: t00 = 273.15, tii = 263.15	! limit temperatures
	real, parameter :: thom0 = 235.03		! Theoretical homogeneous freezing temp.

!
!  Defining public elementary thermodynamic functions
!
	public :: cal_td
	public :: cal_tdew
	public :: cal_esw
	public :: cal_esi
	public :: cal_qsw
	public :: cal_qsw1
	public :: cal_qsi
	public :: cal_qsi1
	public :: cal_dqsdt
	public :: cal_dqsidt
	public :: cal_dqsdp
	public :: cal_dqsidp
	public :: cal_flv
	public :: cal_fls
	public :: cal_flm
	public :: cal_sigmawv
	public :: cal_sigmaiw
	public :: cal_sigmaiv
	public :: cal_dfact
	public :: cal_dgact
	public :: cal_dfads
	public :: cal_xfkd
	public :: cal_xfmu
	public :: cal_aw
	public :: cal_fpd

	CONTAINS

!	==========================================

	real function cal_td( T, p, qt, qh )

	real	:: T, p, qt, qh
	real    :: a

	if (T >= t00) then
	  a = log(max(qt-qh,1.e-15)/cal_qsw( T, p )) + 17.62*(T - t00)/(T - 30.03)
	  cal_td = t00 + a*(t00 - 30.03)/(17.62 - a)
	else
	  a = log(max(qt-qh,1.e-15)/cal_qsi( T, p )) + 21.875*(T - t00)/(T - 7.65)
	  cal_td = t00 + a*(t00 - 7.65)/(21.875 - a)
	endif
	
	return                                                                   
	end function

!	==========================================

	real function cal_tdew( T, rh )

	real	:: T, rh
	real    :: g, tc

	tc = T - t00
	g = log( 0.01*min(max(rh,1.e-6),100.) ) + 17.62*tc/(tc + 243.12)

	cal_tdew = t00 + 243.12*g / (17.62 - g)
	
	return                                                                   
	end function

!	==========================================

	real function cal_esw( T )

	real	:: T, a

	a = 17.62*(T - t00)/(T - 30.03)
	
	cal_esw = 611.2 * exp(a)
	
!	cal_esw = exp(54.842763 - 6763.22/T - 4.21*log(T) + 0.000367*T + tanh(0.0415*(T-218.8))*(53.878 - 1331.22/T - 9.44523*log(T) + 0.014025*T))

	return                                                                   
	end function

!	==========================================

	real function cal_esi( T )

	real :: T, a

	a = 21.875*(T - t00)/(T - 7.65)
	
	cal_esi = 611.2 * exp(a)
	
	return 
	end function  
	
!	==========================================

	real function cal_qsw( T, p )

	real	:: T, p, es

	es = cal_esw (T)
!
	cal_qsw = eps * es / (p - es) 

	return                                                                   
	end function

!	==========================================

	real function cal_qsi( T, p )

	real :: T, p, ei

	ei = cal_esi (T)
!
	cal_qsi = eps * ei / (p - ei) 
	
	return 
	end function  
	
!	==========================================

	real function cal_qsw1( T, p )

	real	:: T, p, es

	es = cal_esw (T)
!
	cal_qsw1 = (p - es) / (eps * es)

	return                                                                   
	end function

!	==========================================

	real function cal_qsi1( T, p )

	real :: T, p, ei

	ei = cal_esi (T)
!
	cal_qsi1 = (p - ei) / (eps * ei)
	
	return 
	end function  

!	==========================================

	real function cal_dqsdp( T, p )

	real :: T, p, es
!
	es = cal_esw (T)

	cal_dqsdp = - eps * es/(p - es)**2.

	return                                                                   
	end function

!	==========================================

	real function cal_dqsidp( T, p )

	real :: T, p, es
!
	es = cal_esi (T)
	
	cal_dqsidp = - eps * es/(p - es)**2.

	return                                                                   
	end function

!	==========================================

	real function cal_dqsdt( T, p )

	real :: T, p, qs, desdt, es
!
!  Clausius-Clapeyron
!	
	es = cal_esw (T)

	desdt = 4283.7744 / (T - 30.03)**2. * es

	cal_dqsdt = desdt * eps*p/(p - es)**2.

	return                                                                   
	end function

!	==========================================

	real function cal_dqsidt( T, p )

	real :: T, p, qi, desdt, es
!
!  Clausius-Clapeyron
!	
	es = cal_esi (T)

	desdt = 5807.8125 / (T - 7.65)**2. * es
	
	cal_dqsidt = desdt * eps*p/(p - es)**2.

	return                                                                   
	end function

!	==========================================
	!> Calculate Latent heat of vaporization 
	real function cal_flv( T )

	real :: T
!	
	if (cst_cp == 1) then
	  cal_flv = flv00
	else
	  cal_flv = flv00 - (cp_l - cp_v)*(T - 273.15)
	endif

	return                                                                   
	end function

!	==========================================

	real function cal_fls( T )

	real :: T
!	
	if (cst_cp == 1) then
	  cal_fls = fls00
	else
	  cal_fls = fls00 - (cp_i - cp_v)*(T - 273.15)
	endif

	return                                                                   
	end function

!	==========================================

	real function cal_flm( T )

	real :: T
!	
	if (cst_cp == 1) then
	  cal_flm = fls00 - flv00
	else
	  cal_flm = cal_fls( T ) - cal_flv( T )
	endif

	return                                                                   
	end function

!	==========================================

	real function cal_sigmaiw( T )

	real :: T
!
!  Surface tension of ice/water interface (from Barahona, ACP 2014)
!	
	cal_sigmaiw = 3.04e-4*T - 0.04919

	return                                                                   
	end function

	!==========================================
	!> Surface tension of water/vapor interface
	!> Following Khain, A. P., & Pinsky, M. (2018). Physical Processes in Clouds and Cloud Modeling. Cambridge: Cambridge University Press. https://doi.org/10.1017/9781139049481
    !> 
	!> Note: Use -5.5°C as threshold which calculation is used
	!>       (there the smallest difference between the two formulations occur)
	!> 
	real function cal_sigmawv( T )

	real :: T, Tc
!
!  Surface tension of water/vapor interface
!	
	Tc = T - t00
	if (Tc <=-5.5) then
	  cal_sigmawv = (75.93 + Tc*(0.115 + Tc*(6.818e-2 + Tc*(6.511e-3 + Tc*(2.933e-4 + Tc*(6.283e-6 + Tc*5.285e-8)))))) * 1.e-3
	else if (Tc > -5.5) then
	  cal_sigmawv = 7.566165e-2 - 1.55e-4*Tc
	endif

	return                                                                   
	end function

!	==========================================

	real function cal_sigmaiv( T )

	real :: T
!
!  Surface tension of vapor/ice interface
!	
	cal_sigmaiv = cal_sigmaiw (T) + cal_sigmawv (T)

	return                                                                   
	end function

!	==========================================

	real function cal_dfact( T )

	real :: T, Tc
!
!  Activation energy for water molecule diffusion across water/ice interface
!	
	Tc = T - t00
	cal_dfact = 0.694e-19 * (1. + 0.0027*(Tc + 30.))	
!	cal_dfact = 3.85606e-20 * exp(-8.423e-3*Tc + 6.384e-4*Tc*Tc + 7.891e-6*Tc*Tc*Tc) 	! Fornea et al 2009
	
	return                                                                   
	end function

!	==========================================

	real function cal_dgact( T )

	real :: T, Tc, E
!
!  Activation energy for homogeneous freezing, from Zobrist et al. 2007
!	
	Tc = T - 118.
	E = 892.
!
	cal_dgact = kb * E * (T/Tc)**2.
	
	return                                                                   
	end function

!	==========================================

	real function cal_dfads( k, T )

	integer :: k
	real	:: T, Tc
	real, dimension(3)    :: A = (/0., 2.1e-19, 2.e-19/)		! Illite: 2.3e-19, soot: 2.1e-19, Natural dust: 2e-19
!
!  Activation energy for water molecule adsorption for deposition freezing		! Welti et al. 2014
!	
	Tc = T - t00
!	cal_dfads = 4.4e-20
	cal_dfads = A(k) * exp(0.03678*Tc)
	
	return                                                                   
	end function
	
!	================================================

	real function cal_xfkd (T, p)

	USE shared_data
	
	IMPLICIT NONE
	
	real :: T, p

!	cal_xfkd = 10./p * (1.53e-3*max(150.0,T) - 0.192)
	cal_xfkd = 3.e-5

	return
	end function

!	================================================

	real function cal_xfmu (T)

	USE shared_data
	
	IMPLICIT NONE
	
	real     :: T

!	cal_xfmu = 1.72e-5*(3.93e2/(T + 1.2e2)) * (T/2.73e2)**(3./2.)
	cal_xfmu = 1.409e-5

	return
	end function
!
!  ==================================================
  real function cal_aw ( h, T, sw, qc, nc )
!  ==================================================

  integer :: h
  real  :: sw, T, qc, nc
  real  :: Tm, b, aw, m
!
!  Water activity is computed using bulk approach from from Khvorostyanov and Curry, JGR 99, 06. 
!  Requires specification of bulk aerosol properties in cm.nml
!
  if (lndrop == 1) then
    m  = qc/nc
    aw = 1. / (1. + pi/6.*rhow*aerob(h)%b*xnc0_d**(2.*(1. + aerob(h)%beta))/m)
    cal_aw = max(aw,0.)
  else
    cal_aw = 1.
  endif

  return								   
  end function
!
!  ==================================================
  real function cal_fpd ( h, T, sw, qc, nc )
!  ==================================================

  integer :: h
  real  :: sw, T, qc, nc
  real  :: Tm, b, aw, m
  real  :: cal_f, cal_fprime, epsilon
!
!  Calculates the freezing point depression considering ice water activity from
!  Koop et al., Nature 2000. Pressure and curvature corrections are neglected.
!  Water activity is computed using bulk approach from from Khvorostyanov and Curry, JGR 99, 06. 
!
  if (lndrop == 1) then
    Tm = Thom0
    aw = cal_aw( h, T, sw, qc, nc )
!  
    epsilon = 1.
    do while (abs(epsilon) > 1.e-3)
      cal_f = 210368. + (131.438-8.314*log(aw-delta_aw))*Tm - 3323730./Tm - 41729.1*log(Tm)
      cal_fprime = (131.438-8.314*log(aw-delta_aw)) + 3323730./Tm**2. - 41729.1/Tm
!
      epsilon = - cal_f / cal_fprime
      Tm = Tm + epsilon
    enddo
  
    cal_fpd = Tm - Thom0
  else
    cal_fpd = 0.
  endif

  return								   
  end function

  end module shared_thermo
