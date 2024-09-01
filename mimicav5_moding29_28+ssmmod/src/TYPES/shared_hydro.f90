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
!  SHARED_HYDRO.F:                   
!
!  Purpose:
!	Module containing hydrometeor variables (replaces hydro12.h) 
!
!  Author
!	Julien Savre, MISU
!
! ================================================================

	module shared_hydro

	use gridno
	use shared_data
	use typedef_hydrometeor
		
	SAVE
!
	type (hydrometeor), dimension(:), allocatable :: hydromtr, hydromtrs, hydromtrp, hydromtr2, hydromtr_mean
!
	type(hydro_const), dimension(:), allocatable :: hydroc
	type(hydro_param), dimension(:), allocatable :: hydrop

	type(hydro_param), dimension(9)  							&
		 :: icep = (/ 									&
!					   (cm, ctv, am, btv, nu, mu, cap, a0)
		 	        hydro_param(0.88, 16.53, 0.3584, 0.2222, -1./3., 1./3., 1., -0.0604, 0.028, 0.287, 0., 4.42e-14, 1.e-7), 		&	! ice (plates, PLA)
		 	        hydro_param(5.1635, 11.055, 0.4367, 0.21, -1./3., 1./3., 1., 0.0355, 0.0355, 0.287, 0., 4.42e-14, 1.e-7), 		&	! ice (dendrites, DEN)
		 	        hydro_param(7.39, 61.857, 0.3584, 0.2222, -1./3., 1./3., 1., 0.0355, 0.0355, 0.287, 0., 4.42e-14, 1.e-7), 		&	! ice (broad dendrites, BDE)
		 	        hydro_param(13.447, 107.074, 0.5555, 0.271, -1./3., 1./3., 1., 0.0309, 0.1447, 0.3183, 0., 4.42e-14, 1.e-7), 		&	! ice (column, COL)
		 	        hydro_param(2.032, 30.376, 0.4955, 0.271, -1./3., 1./3., 1., 0.0309, 0.1447, 0.3183, 0., 4.42e-14, 1.e-7), 		&	! ice (hollow column, HCO)
		 	        hydro_param(3.6384, 610.22, 0.44, 0.3524, -1./3., 1./3., 1., 0.3005, -0.0022, 0.3183, 0., 4.42e-14, 1.e-7), 		&	! ice (bullet rosette C2a, BUL)
		 	        hydro_param(4.4066, 21.512, 0.4762, 0.195, -1./3., 1./3., 0.86, 0.28, 0., 0.3183, 0., 4.42e-14, 1.e-7), 		&	! ice (polycrystal assemblage S3, ASS)
		 	        hydro_param(3.1, 2.153, 0.4545, 0.0545, -1./3., 1./3., 0.86, 0.28, 0., 0.3183, 0., 4.42e-14, 1.e-7), 		&	! ice (plate assemblage P7a, PAS)
		 	        hydro_param(0.2828, 6.38, 0.3333, 0.16667, -1./3., 1./3., 0.86, 0.28, 0., 0.3183, 0., 4.42e-14, 1.e-7) 	/)	 	! ice (default, ISDAC, DEF)
!	
!  List of functions included
!			
	interface cal_lambda
	    MODULE PROCEDURE cal_lambda_1, cal_lambda_3
	end interface cal_lambda
!			
	interface cal_avm 
	    MODULE PROCEDURE cal_avm_1, cal_avm_3
	end interface cal_avm
!			
	interface cal_vp
	    MODULE PROCEDURE cal_vp_1, cal_vp_3
	end interface cal_vp
!			
	interface cal_vpn
	    MODULE PROCEDURE cal_vpn_1, cal_vpn_3
	end interface cal_vpn

	public :: cal_lambda
	public :: cal_avm, cal_avd, cal_avd3
        public :: cal_vp, cal_vpn
	public :: cal_reff
	public :: cal_n1
	public :: cal_gamma
	public :: cal_ingamma


	CONTAINS

! ============================
	
	function cal_lambda_1 (h, q, n)

	IMPLICIT NONE
	
	integer  :: h
	real     :: q, n, lambda, cal_lambda_1

        lambda = 0.
        
	if (moments == 2 .or. h == drop .or. h == ice) then
          if ( q > qmin(h) .and. n > xnmin(h) ) lambda = (hydroc(h)%clam*cal_avm(h,q,n))**(-hydrop(h)%mu)
	else
          if ( q > qmin(h) ) lambda = (hydroc(h)%clam1*q)**(-hydrop(h)%mu/(hydrop(h)%nu+2.))
	endif
    
        if (h == drop .and. lambda > 0.) lambda=min(max(lambda,lmin(drop)),lmax(drop))
        if (h == ice .and. lambda > 0.) lambda=min(max(lambda,lmin(ice)),lmax(ice))
        if (h == rain .and. lambda > 0.) lambda=min(max(lambda,lmin(rain)),lmax(rain))
        if (h == grau .and. lambda > 0.) lambda=min(max(lambda,lmin(grau)),lmax(grau))
        if (h == snow .and. lambda > 0.) lambda=min(max(lambda,lmin(snow)),lmax(snow))
        if (h == hail .and. lambda > 0.) lambda=min(max(lambda,lmin(hail)),lmax(hail))

	cal_lambda_1 = lambda
	
	end function cal_lambda_1

! ============================
		
	function cal_lambda_3 (h, q, n) 

	IMPLICIT NONE
	
	integer  :: h
	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: q, n 
	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: lambda, cal_lambda_3

        lambda = 0.

	if (moments == 2 .or. h == drop .or. h == ice) then
	  where ( q > qmin(h) .and. n > xnmin(h) ) lambda = (hydroc(h)%clam*cal_avm(h,q,n))**(-hydrop(h)%mu)
	else
	  where ( q > qmin(h) ) lambda = (hydroc(h)%clam1*q)**(-hydrop(h)%mu/(hydrop(h)%nu+2.))
	endif
    
        if (h == drop) lambda=min(max(lambda,lmin(drop)),lmax(drop))
        if (h == ice) lambda=min(max(lambda,lmin(ice)),lmax(ice))
        if (h == rain) lambda=min(max(lambda,lmin(rain)),lmax(rain))
        if (h == grau) lambda=min(max(lambda,lmin(grau)),lmax(grau))
        if (h == snow) lambda=min(max(lambda,lmin(snow)),lmax(snow))
        if (h == hail) lambda=min(max(lambda,lmin(hail)),lmax(hail))

	cal_lambda_3 = lambda
	
	end function cal_lambda_3

! ============================
	
	function cal_avm_1 (h, q, n)

	IMPLICIT NONE
	
	integer  :: h
	real     :: q, n, avm, cal_avm_1

        avm = 0.

        if (moments == 2 .or. h == drop .or. h == ice) then
          if ( q > qmin(h) .and. n > xnmin(h) ) avm = min(max(q/n, hydrop(h)%xmin), hydrop(h)%xmax)
        else
          if ( q > qmin(h) ) avm = min(max(hydroc(h)%cmm*(hydroc(h)%clam1*q)**(1./(hydrop(h)%nu+2.)), hydrop(h)%xmin), hydrop(h)%xmax)
        endif

        cal_avm_1 = avm
	
	return
	end function

! ============================
	
	function cal_avm_3 (h, q, n)

	IMPLICIT NONE
	
	integer  :: h
	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: q, n, avm, cal_avm_3

        avm = 0.

        if (moments == 2 .or. h == drop .or. h == ice) then
          where ( q > qmin(h) .and. n > xnmin(h) ) avm = min(max(q/n, hydrop(h)%xmin), hydrop(h)%xmax)
        else
          where ( q > qmin(h) ) avm = min(max(hydroc(h)%cmm*(hydroc(h)%clam1*q)**(1./(hydrop(h)%nu+2.)), hydrop(h)%xmin), hydrop(h)%xmax)
        endif

        cal_avm_3 = avm
	
	return
	end function

! ============================
	
	function cal_avd (h, q, n)

	IMPLICIT NONE
	
	integer  :: h
	real     :: q, n, avd, cal_avd
        
        avd = 0.

        if (moments == 2 .or. h == drop .or. h == ice) then
          if ( q > qmin(h) .and. n > xnmin(h) ) avd = hydroc(h)%cdm*hydroc(h)%clam**(hydrop(h)%am) * hydrop(h)%cm*cal_avm(h, q, n)**(hydrop(h)%am)
        else
          if ( q > qmin(h) ) avd = hydroc(h)%cdm*hydrop(h)%cm*(hydroc(h)%clam1*q)**(hydrop(h)%am/(hydrop(h)%nu+2.))      
        endif

        cal_avd = avd

	return
	end function

! ============================
	
	function cal_avd3 (h, q, n)

	IMPLICIT NONE
	
	integer  :: h
	real     :: q, n, avd3, cal_avd3

        avd3 = 0.

        if (moments == 2 .or. h == drop .or. h == ice) then
          if ( q > qmin(h) .and. n > xnmin(h) ) avd3 = hydroc(h)%cvm*hydroc(h)%clam**(3.*hydrop(h)%am) * (hydrop(h)%cm*cal_avm(h, q, n)**(hydrop(h)%am))**3.
        else
          if ( q > qmin(h) ) avd3 = hydroc(h)%cvm*(hydrop(h)%cm*(hydroc(h)%clam1*q)**(hydrop(h)%am/(hydrop(h)%nu+2.)))**3.      
        endif

        cal_avd3 = avd3
	
	return
	end function

! ============================
	
	function cal_reff (h, q, n)

	IMPLICIT NONE
	
	integer  :: h
	real     :: q, n, reff, cal_reff
        
        reff = 0.

        if (moments == 2 .or. h == drop .or. h == ice) then
          if ( q > qmin(h) .and. n > xnmin(h) ) reff = 0.5*hydroc(h)%cem*hydroc(h)%clam**(hydrop(h)%am) * hydrop(h)%cm*cal_avm(h, q, n)**(hydrop(h)%am)
        else
          if ( q > qmin(h) ) reff = 0.5*hydroc(h)%cem*hydrop(h)%cm*(hydroc(h)%clam1*q)**(hydrop(h)%am/(hydrop(h)%nu+2.))
        endif

        cal_reff = reff

	return
	end function

! ============================	

  function cal_n1 (h, q, n)

  IMPLICIT NONE

  integer :: h
  real    :: q, n, n1, cal_n1

  n1 = 0.

  if (moments == 2 .or. h == drop .or. h == ice) then
    n1 = n
  else
    if ( q > qmin(h) ) n1 = hydroc(h)%clam/hydroc(h)%clam1 * (hydroc(h)%clam1*q)**((hydrop(h)%nu+1.)/(hydrop(h)%nu+2.))
  endif

  cal_n1 = n1

  return
  end function

! ============================

  function cal_vp_3 (h, q, n, pxx)

  IMPLICIT NONE

  integer :: i, j, h, k
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: cal_vp_3
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: q, n, n0
  real :: m, pxx

  cal_vp_3 = 0.

  if (moments == 2 .or. h == drop .or. h == ice) then
    do j = jt_start, jt_end
      do i = it_start, it_end
!        if (h == rain) then
!          n0 = max(250000.,min(2.e7,n*(3141./min(5.e-6,max(2.6e-10,q/n)))**(1./3.)))
!          cal_vp_3 = max((9.65 - 10.3*(1.+600./min(max((3141.*n0/q)**(1./4.),1000.),10000.))**(-4.)) * pxx, 1.e-3)
        if ( q(i,j) > qmin(h) .and. n(i,j) > xnmin(h) ) then
          cal_vp_3(i,j) = hydroc(h)%cvpq*hydroc(h)%clam**(hydrop(h)%btv+1) * hydrop(h)%ctv*cal_avm(h,q(i,j),n(i,j))**hydrop(h)%btv * pxx
        endif
      enddo
    enddo
  else
    do j = jt_start, jt_end
      do i = it_start, it_end
        if ( q(i,j) > qmin(h) ) then
          cal_vp_3(i,j) = hydroc(h)%cvpq*hydrop(h)%ctv*(hydroc(h)%clam1*q(i,j))**(hydrop(h)%btv/(hydrop(h)%nu+2.)) * pxx
        endif
      enddo
    enddo
  endif

  return
  end function

! ============================

  function cal_vpn_3 (h, q, n, pxx)
  
  IMPLICIT NONE

  integer :: i, j, h, k
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: cal_vpn_3
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: q, n, n0
  real :: m, pxx

  cal_vpn_3 = 0.

  if (moments == 2 .or. h == drop .or. h == ice) then
    do j = jt_start, jt_end
      do i = it_start, it_end
!        if (h == rain) then
!          n0 = max(250000.,min(2.e7,n*(3141./min(5.e-6,max(2.6e-10,q/n)))**(1./3.)))
!          cal_vpn_3 = max((9.65 - 10.3*(1.+600./min(max((3141.*n0/q)**(1./4.),1000.),10000.))**(-1.)) * pxx, 1.e-3)
        if ( q(i,j) > qmin(h) .and. n(i,j) > xnmin(h) ) then
          cal_vpn_3(i,j) = hydroc(h)%cvpn*hydroc(h)%clam**(hydrop(h)%btv) * hydrop(h)%ctv*cal_avm(h,q(i,j),n(i,j))**hydrop(h)%btv * pxx
        endif
      enddo
    enddo
  else
    do j = jt_start, jt_end
      do i = it_start, it_end
        if ( q(i,j) > qmin(h) ) then
          cal_vpn_3(i,j) = hydroc(h)%cvpn*hydrop(h)%ctv*(hydroc(h)%clam1*q(i,j))**(hydrop(h)%btv/(hydrop(h)%nu+2.)) * pxx
        endif
      enddo
    enddo
  endif

  return
  end function

! ============================

  function cal_vp_1 (h, q, n, pxx)

  IMPLICIT NONE

  integer :: h
  real :: q, n, m, pxx, cal_vp_1
  real :: n0

  cal_vp_1 = 0.

  if (moments == 2 .or. h == drop .or. h == ice) then
!    if (h == rain) then
!      n0 = max(250000.,min(2.e7,n*(3141./min(5.e-6,max(2.6e-10,q/n)))**(1./3.)))
!      cal_vp_1 = max((9.65 - 10.3*(1.+600./min(max((3141.*n0/q)**(1./4.),1000.),10000.))**(-4.)) * pxx, 1.e-3)
    if ( q > qmin(h) .and. n > xnmin(h) ) then
      cal_vp_1 = hydroc(h)%cvpq*hydroc(h)%clam**(hydrop(h)%btv+1) * hydrop(h)%ctv*cal_avm(h,q,n)**hydrop(h)%btv * pxx
    endif
  else
    if ( q > qmin(h) ) cal_vp_1 = hydroc(h)%cvpq*hydrop(h)%ctv*(hydroc(h)%clam1*q)**(hydrop(h)%btv/(hydrop(h)%nu+2.)) * pxx
  endif

  return
  end function

! ============================

  function cal_vpn_1 (h, q, n, pxx)
  
  IMPLICIT NONE

  integer :: h
  real :: q, n, m, pxx, cal_vpn_1
  real :: n0

  cal_vpn_1 = 0.

  if (moments == 2 .or. h == drop .or. h == ice) then
!    if (h == rain) then
!      n0 = max(250000.,min(2.e7,n*(3141./min(5.e-6,max(2.6e-10,q/n)))**(1./3.)))
!      cal_vpn_1 = max((9.65 - 10.3*(1.+600./min(max((3141.*n0/q)**(1./4.),1000.),10000.))**(-1.)) * pxx, 1.e-3)
    if ( q > qmin(h) .and. n > xnmin(h) ) then 
      cal_vpn_1 = hydroc(h)%cvpn*hydroc(h)%clam**(hydrop(h)%btv) * hydrop(h)%ctv*cal_avm(h,q,n)**hydrop(h)%btv * pxx
    endif
  else
    if ( q > qmin(h) ) cal_vpn_1 = hydroc(h)%cvpn*hydrop(h)%ctv*(hydroc(h)%clam1*q)**(hydrop(h)%btv/(hydrop(h)%nu+2.)) * pxx
  endif

  return
  end function

! ============================
	
	function cal_gamma (x)
	
	! Gamma function calculator

	IMPLICIT NONE
	
	real     :: x, xx, cal_gamma
	real, parameter :: p0=1.000000000190015,        &
	                   p1=76.18009172947146,        &
			   p2=-86.50532032941677,       &
			   p3=24.01409824083091,        &
			   p4=-1.231739572450155,       &
			   p5=1.208650973866179e-3,     &
			   p6=-5.395239384953e-6

	xx = p0 + p1/(x + 1.) + p2/(x + 2.) + p3/(x + 3.) + p4/(x + 4.) + p5/(x + 5.) + p6/(x + 6.)
	cal_gamma = sqrt(2.*pi) * xx/x * (x + 5.5)**(x+0.5) * exp(-(x + 5.5))
	
	return
	end function	

! ============================
	
	function cal_ingamma (h,z)
	
	! Lower Incomplete Gamma function

	IMPLICIT NONE
	
	integer  :: h
	real     :: x, z, cal_ingamma

	x = real (hydrop(h)%nu + 1)
	cal_ingamma = exp(-z) * z**x/x * (1. + z/(x + 1.) * (1. + z/(x + 2.) * (1. + z/(x + 3.) * (1. + z/(x + 4.) * (1. + z/(x + 5.) * (1. + z/(x + 6.)))))))
	
	return
	end function	

			
	end module shared_hydro
