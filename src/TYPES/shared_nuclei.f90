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
!  SHARED_NUCLEI.F:                   
!
!  Purpose:
!	Module containing condensation nuclei	  
!
!  Author
!	Julien Savre, MISU
!
! ================================================================

	module shared_nuclei

	use gridno
	use shared_data
	use typedef_nuclei
		
	SAVE
!
!  General parameters
!	
	type (nuclei) :: nuc2, nuc, nuc_mean
!
	logical :: ldep = .FALSE., limm = .TRUE., lcon = .FALSE., lcont = .FALSE., theta_low = .FALSE., lmode(3,4)
	real :: dmini=10.e-6	! Minimum size of an ice crystal
!
!  Nuclei structure
!
	type (nuclei_3d), dimension(3) :: nucin2, nucin, nucins, nucin_mean
	type (nuclei0), dimension(3) :: nucin0
!
!  Nuclei size mode parameters
!
	type(nuc_param), dimension(3)  								&
		 :: nucp = (/ 									& 
		 		nuc_param(28.e6, 0.094e-6, 1.5,		& ! dust (fine mode)
		 		          2.e6, 0.4e-6, 2.3,		& ! dust (coarse mode)	
					  2100.),			&			
		 		nuc_param(28.e6, 0.094e-6, 1.5,		& ! bc (fine mode)
		 		          2.e6, 0.4e-6, 2.3,		& ! bc (coarse mode)	
					  1000.),			&	
		 		nuc_param(28.e6, 0.094e-6, 1.5,		& ! bio (fine mode)
		 		          2.e6, 0.4e-6, 2.3,		& ! bio (coarse mode)	
					  1000.)	/)
!
#ifdef NUC_CNT
!
!  Nuclei properties for classical nucleation theory
!
	type(nuc_prop), dimension(3) 							& !
		 :: nucprop = (/ 							& ! 3 aerosols types:
!					(frac, fracc, thetad, thetai, thetac)
				 nuc_prop(0.02, 0.1,  					& ! fractions, dust
				 	  (/14.2, 180., 132./),				& ! Contact angles
					  (/1.2, 17.9, 20./) ),				& ! Standard deviations
				 nuc_prop(0.04, 0.04,	 				& ! BC
				 	  (/15.2, 270., 93./),				&
					  (/1., 30.4, 11.5/) ),		        	&
		 		 nuc_prop(0., 0., 					& ! bio-aerosols
				 	  (/360., 360., 360./),				&
					  (/0.1, 0.1, 0.1/) )    	  /)	                
#else
!
!  Nuclei properties for Diehl & Wurzler model (immersion)
!
	type(nuc_prop), dimension(3) 							& !
		 :: nucprop = (/ 							& ! 3 aerosols types:
!					(frac, fracc, Beff, a, ac, bc, theta)
				 nuc_prop(0.33, 0., 61.9, 1., 0.1014, 0.3277, 14.),		& ! dust (20.11: from Broadley et al., 2012)
				 nuc_prop(0.67, 0., 0.00291, 1., 0.0614, 0.5730, 40.),		& ! BC (from Diehl & wurzler, 2004)
		 		 nuc_prop(0., 0., 0., 1., 0., 0., 180.)	/)		 	  ! bio-particles (inactive)
#endif
!
!  List of elementary functions
!
	public :: cal_mu			! Calculate mu parameter for log-normal
	public :: cal_sigma			! Calculate sigma parameter for log-normal
	public :: cal_mean			! Calculate mean of a log-normal
	public :: cal_npdf			! Calculate normal distribution
	public :: cal_cumul			! Cumulative function for normal distribution
	public :: cal_delta_ab			! delta function from Phillips et al. (cubic interpolation)
	public :: erfx				! Error function

	CONTAINS

! ============================
	
	function cal_mu (d,s)
	
	real     :: d, s, cal_mu
	
	! d is the modal diameter (d at max of distribution)
	
	if (s > 0.) then
	  cal_mu = log(d) + log(s)*log(s)
	else
	  cal_mu = log(d)
	endif
	
	return
	end function cal_mu

! ============================
	
	function cal_sigma (s)

	real     :: s, cal_sigma
	
	! s is the geometric standard deviation
	
	if (s > 0.) then
  	  cal_sigma = log(s)
	else
	  cal_sigma = 0.
	endif
	
	return
	end function cal_sigma

! ============================
	
	function cal_mean (d,s)
	
	real     :: s, d, cal_mean
	real     :: sig, mu
	
	mu  = cal_mu (d,s)
	sig = cal_sigma (s)
	
	cal_mean = exp(mu + 0.5*sig*sig)
	
	return
	end function cal_mean
	
! ============================
	
	function cal_npdf (d,s,x)
	
	real     :: d, s, x, cal_npdf
	real     :: sigma2, mu
	
	sigma2 = s*s
	mu     = d
	
	cal_npdf = 1. / sqrt(2.*pi*sigma2) * exp (- (x - mu)**2. / (2.*sigma2))
	
	return
	end function cal_npdf	

! ============================
	
	function cal_cumul (d,s,x)
	
	real     :: s, d, x, cal_cumul
	real     :: sig, mu
	
	mu  = cal_mu (d,s)
	sig = cal_sigma (s)
	
	cal_cumul = 0.5*(1. + erfx((log(x) - mu) / sqrt(2.*sig*sig)))
	
	return
	end function cal_cumul
	
! ============================
		
	function cal_delta_ab (a,b,x,x1,x2)

	IMPLICIT NONE
	
	real     :: aa, bb, a, b, x, x1, x2
	real     :: a0, a1, a2, a3
	real     :: cal_delta_ab
		
	aa = 6.*(a-b)/(x2 - x1)**3.
	bb = a + aa*x1**3./6. - aa*x1*x1*x2/2.
	
	a0 = bb
	a1 = aa*x1*x2
	a2 = -aa*(x1+x2)/2.
	a3 = aa/3.

	if ( x > x1 .and. x < x2 ) then
	  cal_delta_ab = a0 + a1*x + a2*x*x + a3*x*x*x
	else if ( x >= x2 ) then
	  cal_delta_ab = b
	else if ( x <= x1 ) then
	  cal_delta_ab = a
	endif
	
	return
	end function
		
! ============================
	
	function erfx (x)
	
	real     :: x, t, xx, erfx
	real, parameter :: p0=0.3275911,          &
	                   p1=0.254829592,        &
			   p2=-0.284496736,       &
			   p3=1.421413741,        &
			   p4=-1.453152027,       &
			   p5=1.061405429
		
	t = 1. / (1. + p0*abs(x))	   
	xx = 1. - (p1*t + p2*t*t + p3*t*t*t + p4*t*t*t*t + p5*t*t*t*t*t) * exp(-x*x)
	erfx = sign(1.,x)*xx
	
	return
	end function	

! ============================
	
      real function cal_inverror (p)
      
      IMPLICIT NONE
      
      real :: p,p_low,p_high
      real :: a1,a2,a3,a4,a5,a6
      real :: b1,b2,b3,b4,b5
      real :: c1,c2,c3,c4,c5,c6
      real :: d1,d2,d3,d4
      real :: z,q,r
      
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1.-p_low
      
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
      
201   q=sqrt(-2.*log(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/			&
     	((((d1*q+d2)*q+d3)*q+d4)*q+1.)
      goto 204
      
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
      
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/		&
     	(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
      
302   if((p.gt.p_high).and.(p.lt.1.)) goto 203
203   q=sqrt(-2.*log(1.-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/			&
     	 ((((d1*q+d2)*q+d3)*q+d4)*q+1.)

204   cal_inverror=z
      
      return
      end function				
	
	end module shared_nuclei
