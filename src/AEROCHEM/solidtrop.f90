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
!  SOLIDTROP.F:
!
!  Purpose:
!	A package for calculating solid phase reactions in troposphere.			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
! ================================================================

!	==================================================
	subroutine solidtrop ( ifss, gas, solidi       )
!	==================================================

	USE gridno
	USE typedef_gas
	USE typedef_solid
	USE shared_data
	USE shared_state
	USE shared_thermo
	USE shared_hydro

	IMPLICIT NONE
	
! --------------------------------------------------
!  Brief Description of Dummy Variables:
!
!	dtr:	time step of integration in second
!	ifss:	1 for using steady-state assumption
!		0 for not using
!       ktrop:  layer label of tropopause
!
! --------------------------------------------------
!  Units of Variables:
!
!	Solid phase:
!		inside:  ppb(m)
!		outside: ppb(m)
!	Gaseous phase:
!		inside:	 ppb(m)
!		outside: ppb(m)
!
! --------------------------------------------------
	
	integer :: ifss, i, j, k, ind
	real    :: pres, den, r_den, Temp, surf
	logical :: l_qi
	
	type (solid_chemical) :: gamma	! reaction possibility

#ifdef CHEM_ENABLE
	type (gas_chemical),                              &
       dimension(ip_start:ip_end,jp_start:jp_end,1:nz)    &
     			       :: gas
#else
        type (gas_chemical) :: gas
#endif

#ifdef SOLIDCHEM_ENABLE
	type (solid_chemical),                             &
       dimension(ip_start:ip_end,jp_start:jp_end,1:nz)     &
     			       :: solidi
#else
        type (solid_chemical) :: solidi
#endif

#ifdef SOLIDCHEM_ENABLE

	! --- note: currently 8 diagnostic variables
	integer				:: mnv
	real, dimension (nz,8)		:: tdiag, t2diag
	real, dimension (nz,8)		:: gdiag, g2diag
	real, dimension (nz)		:: wdiag

! ===============================================================

!#include "patch_size.h"

	if (lmicro > 1) then
!
	tdiag(:,:) = 0.0
	do k=1,ktrop
	do j=jt_start,jt_end
	do i=it_start,it_end
	  tdiag(k, 1) = tdiag(k, 1) + solidi(i,j,k)%o3
	  tdiag(k, 2) = tdiag(k, 2) + solidi(i,j,k)%h2o2
	  tdiag(k, 3) = tdiag(k, 3) + solidi(i,j,k)%xnv
	  tdiag(k, 4) = tdiag(k, 4) + solidi(i,j,k)%ch2o
	  tdiag(k, 5) = tdiag(k, 5) + solidi(i,j,k)%ch3o2h
	  tdiag(k, 6) = tdiag(k, 6) + solidi(i,j,k)%siv
	  tdiag(k, 7) = tdiag(k, 7) + solidi(i,j,k)%svi
	  tdiag(k, 8) = tdiag(k, 8) + solidi(i,j,k)%nh4
	end do
	end do
	end do

	main_loop: do k=1,ktrop

	do j=jt_start,jt_end
	do i=it_start,it_end

	  l_qi = hydromtr1(i,j,k,ice)%q.gt.qmin(ice).and.   &
     		hydromtr1(i,j,k,ice)%n.gt.xnmin(ice)
	  
	if ( l_qi ) then
	
	  pres = 100000.*state2(i,j,k)%p**(1.005e3/2.87e2)
	  den  = pres/(2.87e2*thermo%T(i,j,k))
	  r_den= 1./den
	  Temp = thermo%T(i,j,k)

	  !
	  ! --- Surface area (m^2/m^3)
	  !	derived based on bullet crystal with d = 0.25L
	  !	Am = 0.0088
	  !
	  surf = 100.406*den*hydromtr1(i,j,k,ice)%q

	  call calc_cc0 (awO3,   Temp, surf, gamma_o3,                &
     			gas(i,j,k)%o3, solidi(i,j,k)%o3)
         call calc_cc0 (awH2O2, Temp, surf, gamma_h2o2, 	      &
     			gas(i,j,k)%h2o2, solidi(i,j,k)%h2o2)
         call calc_cc0 (awHNO3, Temp, surf, gamma_xnv,  	      &
     			gas(i,j,k)%hno3, solidi(i,j,k)%xnv)
         ! --- Note: N(v) includes N2O5
         call calc_cc0 (awN2O5, Temp, surf, gamma_xnv,  	      &
     			gas(i,j,k)%xn2o5, solidi(i,j,k)%xnv)
         call calc_cc0 (awCH2O, Temp, surf, gamma_ch2o, 	      &
     			gas(i,j,k)%ch2o, solidi(i,j,k)%ch2o)
         call calc_cc0 (awCH3O2H, Temp, surf, gamma_ch3o2h, 	      &
     			gas(i,j,k)%ch3o2h, solidi(i,j,k)%ch3o2h)
         call calc_cc0 (awSO2,  Temp, surf, gamma_siv,  	      &
     			gas(i,j,k)%so2, solidi(i,j,k)%siv)
         call calc_cc0 (awH2SO4,Temp, surf, gamma_svi,  	      &
     			gas(i,j,k)%h2so4, solidi(i,j,k)%svi)
         call calc_cc0 (awNH4, Temp, surf, gamma_nh4, 		      &
     			gas(i,j,k)%nh3, solidi(i,j,k)%nh4)
	end if

	end do
	end do

	end do main_loop
!	----------------

	t2diag(:,:) = 0.0
	do k=1,ktrop
	do j=jt_start,jt_end
	do i=it_start,it_end
	  t2diag(k, 1) = t2diag(k, 1) + solidi(i,j,k)%o3
	  t2diag(k, 2) = t2diag(k, 2) + solidi(i,j,k)%h2o2
	  t2diag(k, 3) = t2diag(k, 3) + solidi(i,j,k)%xnv
	  t2diag(k, 4) = t2diag(k, 4) + solidi(i,j,k)%ch2o
	  t2diag(k, 5) = t2diag(k, 5) + solidi(i,j,k)%ch3o2h
	  t2diag(k, 6) = t2diag(k, 6) + solidi(i,j,k)%siv
	  t2diag(k, 7) = t2diag(k, 7) + solidi(i,j,k)%svi
	  t2diag(k, 8) = t2diag(k, 8) + solidi(i,j,k)%nh4
	end do
	end do
	end do

	tdiag(:,:) = t2diag(:,:) - tdiag(:,:)

	endif
#endif

	return
	 end

!	===========================================================
	subroutine calc_cc0 ( aw, Temp, surf, gamma, gasc, solidc )
!	===========================================================

	USE shared_data

!
! === Program for calculating dgas due to gaseous adsorption 
!	on solid surface
!	dc/dt = -vt*gamma*surf*c/4
!
!	Note: both gasc and solidc is in ppb(m)
!
	
	real :: aw, Temp, surf, gamma, gasc, solidc, vt, dgas
	
	vt = 145.534*sqrt(Temp/aw)	! m/s
	
	dgas   = gasc*( exp(-vt*gamma*surf*0.25*dt) - 1.0 )
	gasc   = gasc + dgas
	solidc = solidc - dgas
	
	return
	  end subroutine calc_cc0
	  
	
