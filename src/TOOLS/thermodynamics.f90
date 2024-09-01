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
!  THERMODYNAMICS.F:                   
!
!  Purpose:
!	Module to calculate thermodynamic properties: the main equation_of_state
!   	routine takes either the prognostic density or pressure as an input
!		(plus the other prognostic variables pt, qt and hydrometeors), and
!		returns diagnostic properties used in the code (pressure, density, 
!		temperature, saturation mixing ratio...). 
!		Outputs are optional so that only the desired diagnostic properties are
!		computed and output when equation_of_state is called. 	  
!
!  Author
!	Julien Savre, MIM, Ludwig-Maximilian Univeristät, Munchen
!		
! ================================================================

	module thermodynamics

	use gridno
	use shared_data
	USE shared_state
	use shared_thermo
	USE shared_hydro
	USE shared_wind
	USE shared_diag
	USE shared_pressure
	USE averages
!	
	private
	
	interface equation_of_state
	    MODULE PROCEDURE equation_of_state
	end interface equation_of_state
	
	interface get_cp_cv
            MODULE PROCEDURE get_cp_cv_3d, get_cp_cv_0d, get_cp_cv_1d
	end interface get_cp_cv
	
	interface get_temperature
            MODULE PROCEDURE get_temperature_3d
	end interface get_temperature
	
	public :: equation_of_state, get_cp_cv, buoyancy
	public :: get_temperature, get_pt, get_ptv, get_mse, get_qv, get_rh, get_rhi, get_tdew, get_supersaturation

	CONTAINS
	
!	==========================================	
	  subroutine equation_of_state ( state_in, hydromtr_in, pressure_in, density_in,        &
	  		pressure_out, density_out, mse_out, pt_out, thermo_out, 		&
			conservative, with_pt )
!	==========================================
!
	  IMPLICIT NONE
!
	  logical, intent(in), optional :: conservative, with_pt
	  type(atm_state), intent(inout) :: state_in
	  type(hydrometeor), dimension(:), intent(inout) :: hydromtr_in
!
	  real, dimension(ip_start:,jp_start:,:), optional, intent(inout) :: pressure_in, density_in
!
	  real, dimension(ip_start:,jp_start:,:), optional, intent(out) :: pressure_out, density_out, mse_out, pt_out
			
	  type (thermodyn), optional, intent(inout) :: thermo_out
!
!  Local variables
!
	  logical :: cons, mse_in, pt_in
	  integer :: i, j, k, h
	  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: pp, exn, qv, tt
!
	  if (verbose > 1) call write_debug('Starting thermodynamics.f90')
!
!-----------------------------------------------------------!
!	 First thing: calculate thermo properties	    !
!-----------------------------------------------------------!
!
	  mse_in = .false.
	  pt_in = .false.
#ifndef ISENTROPIC
  	  mse_in = .true.
#else
	  pt_in = .true.
#endif
!
	  if (present(with_pt)) then
	    if (with_pt) then
	      pt_in=.true.
  	      mse_in=.false.
	    endif
	  endif
!
  	  if (present(conservative)) then
	    cons = conservative
	  else
	    cons = .false.
	  endif
!
	  if (cons .and. .not.present(density_in)) print*, 'ERROR IN THERMODYNAMICS: Conserved variables are used but no density provided'
!
  	  if (cons) then
	    state_in%es = state_in%es / density_in
	    state_in%qt = state_in%qt / density_in
	    do h = 1, nhydro
              hydromtr_in(h)%q = hydromtr_in(h)%q / density_in
	    enddo
	  endif
!
 	  call get_cp_cv ( state_in%qt, hydromtr_in )
!
!-----------------------------------------------------------!
!			Initializations		  	    !
!-----------------------------------------------------------!
!
!  Full pressure
!
#ifdef ANELASTIC
  	    do k = 1, nz
	      pp(:,:,k) = p0(k)
            enddo
#else
  	    if (present(pressure_in)) then
  	      do k = 1, nz
    	        pp(:,:,k) = pressure_in(:,:,k) + p0(k)
  	      enddo
  	    endif	  
#endif
!
!  Water vapour mixing ratio
!
	    qv = state_in%qt
!
  	    if (lmicro > 0) qv = qv - (hydromtr_in(drop)%q+hydromtr_in(rain)%q)
!
  	    if (lmicro > 1) qv = qv - hydromtr_in(ice)%q
!
  	    if (lmicro > 2) qv = qv - (hydromtr_in(grau)%q + hydromtr_in(snow)%q)
!
  	    if (lmicro > 3) qv = qv - hydromtr_in(hail)%q
!
!  Temperature from mse
!
 	    if (mse_in) then
	      do k = 1, nz
	        tt(:,:,k) = state_in%es(:,:,k) - flv00*qv(:,:,k) - g*z0(k)
	      enddo
	      if (lmicro > 1) tt = tt - (flv00-fls00)*hydromtr_in(ice)%q
	      if (lmicro > 2) tt = tt - (flv00-fls00)*(hydromtr_in(grau)%q + hydromtr_in(snow)%q)
	      if (lmicro > 3) tt = tt - (flv00-fls00)*hydromtr_in(hail)%q
	      tt = 273.15 + tt / thermo_prop%cp
	    endif
!
!-----------------------------------------------------------!
!		   Get Exner pressure			    !
!-----------------------------------------------------------!
!
!  Compressible equations: need to consider full pressure
!  Anelastic approximation (Lipps and Hemler): exner pressure 
!     remains equal to its base state values for consistency
!
#ifdef ANELASTIC
  	    exn = (pp / pref)**( (cp_a-cv_a)/cp_a )
#else
  	    if (present(pressure_in)) then
  	      exn = (pp / pref)**( (cp_a-cv_a)/cp_a )
  	    else if (present(density_in).and.mse_in) then
	      exn = (density_in * (thermo_prop%cp-thermo_prop%cv) * tt / pref)**((cp_a-cv_a)/cp_a)
  	    else if (present(density_in).and.pt_in) then
	      exn = (density_in * (thermo_prop%cp-thermo_prop%cv) * state_in%es / pref)**((cp_a-cv_a)/cv_a)
	    endif
#endif
!
!-----------------------------------------------------------!
!		  	Get temperature	  	  	    !
!-----------------------------------------------------------!
!	    
	    if (pt_in) then
	      tt = state_in%es * exn
	    endif
!
!  Get pressure from density
!
#ifndef ANELASTIC
	    if (present(density_in)) then
  	      pp = density_in*(thermo_prop%cp - thermo_prop%cv)*tt
	    endif
#endif
!
!-----------------------------------------------------------!
!		    Other energy variable		    !
!-----------------------------------------------------------!
!
	    if (pt_in.and.present(mse_out)) then
	      do k = 1, nz
	        mse_out(:,:,k) = thermo_prop%cp(:,:,k) * (tt(:,:,k) - 273.15) + flv00*qv(:,:,k) + g*z0(k)
		if (lmicro > 1) mse_out(:,:,k) = mse_out(:,:,k) + (flv00-fls00)*hydromtr_in(ice)%q(:,:,k)
		if (lmicro > 2) mse_out(:,:,k) = mse_out(:,:,k) + (flv00-fls00)*(hydromtr_in(grau)%q(:,:,k) + hydromtr_in(snow)%q(:,:,k))
		if (lmicro > 3) mse_out(:,:,k) = mse_out(:,:,k) + (flv00-fls00)*hydromtr_in(hail)%q(:,:,k)
	      enddo
	    else if (mse_in.and.present(pt_out)) then
	      pt_out = tt / exn
	    endif
!
!-----------------------------------------------------------!
!	 	 Calculate thermo structure		    !
!-----------------------------------------------------------!
!
	    if ( present(thermo_out) ) then
!
!  Temperature
!
 	      thermo_out%T = tt
!
!  Potential temperature
!
	      thermo_out%pt = tt / exn
!
!  Store exner if required
!
	      thermo_out%exn = exn
!
!  Store relative humidity
!
              if ( out_sat ) then
                call get_rh ( tt, pp, state_in, hydromtr_in, thermo%rh )
                call get_rhi ( tt, pp, state_in, hydromtr_in, thermo%rhi )
              endif
!
!  Sound speed
!
#if (!defined ANELASTIC) || (!defined ISENTROPIC)
	      if (present(pressure_in)) then
  	    	thermo_prop%cs2 = cp_a/cv_a * tt * (thermo_prop%cp - thermo_prop%cv)
  	      else if (present(density_in)) then
  	    	thermo_prop%cs2 = cp_a/cv_a * pp / density_in
  	      endif
#endif
!
	    endif
!
!-----------------------------------------------------------!
!		Output pressure or density		    !
!-----------------------------------------------------------!
!
	    if (present(pressure_out).and.present(density_in)) then
	      do k = 1, nz
	    	pressure_out(:,:,k) = pp(:,:,k) - p0(k)
	      enddo
	    endif
!
	    if (present(density_out).and.present(pressure_in)) then
	      density_out = pp / (tt * (thermo_prop%cp-thermo_prop%cv))
	    endif
!
!-----------------------------------------------------------!
!
!  Back to conserved variables
!
  	  if (cons) then
	    state_in%es = state_in%es * density_in
	    state_in%qt = state_in%qt * density_in
	    do h = 1, nhydro
              hydromtr_in(h)%q = hydromtr_in(h)%q * density_in
	    enddo
	  endif
!
	  if (verbose > 1) call write_debug('Terminating thermodynamics.f90')
!	  	  
	  end subroutine equation_of_state
!
!	==========================================	
	  subroutine get_cp_cv_3d ( qt, hydromtrl )
!	==========================================	
!
	  real, dimension(ip_start:,jp_start:,:), intent(in) :: qt
	  type(hydrometeor), dimension(:), intent(in) :: hydromtrl
!
	  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: qv, qc
!
!-----------------------------------------------------------!
!			Initializations			    !
!-----------------------------------------------------------!
!
	    qv = qt
	    qc = 0.0
!
	    if (lmicro > 0) then
	      qv = qv - hydromtrl(drop)%q - hydromtrl(rain)%q
	      qc = qc + hydromtrl(drop)%q + hydromtrl(rain)%q
	    endif   
!
	    if (lmicro > 1) then
	      qv = qv - hydromtrl(ice)%q
	      qc = qc + hydromtrl(ice)%q
	    endif   
!
	    if (lmicro > 2) then
	      qv = qv - hydromtrl(grau)%q - hydromtrl(snow)%q
	      qc = qc + hydromtrl(grau)%q + hydromtrl(snow)%q
	    endif   
!
	    if (lmicro > 3) then
	      qv = qv - hydromtrl(hail)%q
	      qc = qc + hydromtrl(hail)%q
	    endif   
!
!-----------------------------------------------------------!
!			Calculate cp/cv		   	    !
!-----------------------------------------------------------!
!
	if (cst_cp == 1) then
	  thermo_prop%cp = cp_a
	  thermo_prop%cv = cv_a
	else
	  thermo_prop%cp = (1.-qt)*cp_a + cp_v*qv
	  thermo_prop%cv = (1.-qt)*cv_a + cv_v*qv
!	
	  if (lmicro == 1) then
	    thermo_prop%cp = thermo_prop%cp + cp_l*(hydromtrl(drop)%q + hydromtrl(rain)%q)
	    thermo_prop%cv = thermo_prop%cv + cv_l*(hydromtrl(drop)%q + hydromtrl(rain)%q)

	  endif 	  
!	
	  if (lmicro == 2) then
	    thermo_prop%cp = thermo_prop%cp + cp_i*hydromtrl(ice)%q 
	    thermo_prop%cv = thermo_prop%cv + cv_i*hydromtrl(ice)%q 
	  endif
!	
	  if (lmicro == 3) then
	    thermo_prop%cp = thermo_prop%cp + cp_i*(hydromtrl(grau)%q + hydromtrl(snow)%q)
	    thermo_prop%cv = thermo_prop%cv + cv_i*(hydromtrl(grau)%q + hydromtrl(snow)%q)
	  endif
!	
	  if (lmicro == 4) then
	    thermo_prop%cp = thermo_prop%cp + cp_i*hydromtrl(hail)%q
	    thermo_prop%cv = thermo_prop%cv + cv_i*hydromtrl(hail)%q
	  endif
	endif
!
!-----------------------------------------------------------!
!
	  end subroutine get_cp_cv_3d
!
!	==========================================	
	  subroutine get_cp_cv_0d ( qt, ql, qi, cp, cv )
!	==========================================	
!
	  real, intent(in) :: qt, ql, qi
	  real, intent(out) :: cp, cv
	  real  :: qv, qc
!
!-----------------------------------------------------------!
!			Initializations			    !
!-----------------------------------------------------------!
!
        qv = qt
        qc = 0.
!
        if (lmicro > 0) then
          qv = qv - ql
          qc = qc + ql
        endif	
!
        if (lmicro > 1) then
          qv = qv - qi
          qc = qc + qi
        endif	
!
!-----------------------------------------------------------!
!	  	       Calculate cp/cv			    !
!-----------------------------------------------------------!
!
	if (cst_cp == 1) then
	  cp = cp_a
	  cv = cv_a
	else
	  cp = (1.-qt)*cp_a + cp_v*qv
	  cv = (1.-qt)*cv_a + cv_v*qv
!	
	  if (lmicro > 0) then
	    cp = cp + cp_l*ql
	    cv = cv + cv_l*ql

	  endif 	  
!	
	  if (lmicro > 1) then
	    cp = cp + cp_i*qi
	    cv = cv + cv_i*qi
	  endif   
	endif
!
!-----------------------------------------------------------!
!
	  end subroutine get_cp_cv_0d
!
!
	  subroutine get_cp_cv_1d ( qt, ql, qi, cp, cv )
!
	  real, dimension(nz), intent(in) :: qt, ql, qi
	  real, dimension(nz), intent(out) :: cp, cv
	  real, dimension(nz)  :: qv, qc
!
!-----------------------------------------------------------!
!	  	       Initializations	  	  	    !
!-----------------------------------------------------------!
!
        qv = qt
        qc = 0.
! 
        if (lmicro > 0) then
          qv = qv - ql
          qc = qc + ql
        endif   
!
        if (lmicro > 1) then
          qv = qv - qi
          qc = qc + qi
        endif   
!
!-----------------------------------------------------------!
!		      Calculate cp/cv			    !
!-----------------------------------------------------------!
!
	if (cst_cp == 1) then
	  cp = cp_a
	  cv = cv_a
	else
	  cp = (1.-qv-qc)*cp_a + cp_v*qv
	  cv = (1.-qv-qc)*cv_a + cv_v*qv
!	
	  if (lmicro > 0) then
	    cp = cp + cp_l*ql
	    cv = cv + cv_l*ql
	  endif 	  
!	
	  if (lmicro > 1) then
	    cp = cp + cp_i*qi
	    cv = cv + cv_i*qi
	  endif   
	endif
!
!-----------------------------------------------------------!
!
	end subroutine get_cp_cv_1d
!
!	==========================================	
  subroutine buoyancy ( dens, statel, hydromtrl, buoy )
!	==========================================	
!
  logical :: cons
  integer :: k, h
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: bb, pt, den1, sumh
!
  type (atm_state), intent(inout) :: statel
  type (hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrl
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: dens
  real, dimension(ip_start:,jp_start:,:), intent(out) :: buoy
!
!----------------------------------------------------------!
!
!  Initializations
!
  buoy = 0.
!
!  Calculate buoyancy from density
!
#ifndef ANELASTIC
!
  do k = 2, nz
    buoy(:,:,k) = -0.5*g*(dens(:,:,k) + dens(:,:,k-1) - den0(k) - den0(k-1))
  enddo
  buoy(:,:,1) = -g*(dens(:,:,1) - den0(1))
!
#else
!
!  Calculate buoyancy from virtual potential temperature
!
#ifdef CONSERVATIVE
  cons = .true.
  den1 = 1./dens
#else
  cons = .false.
  den1 = 1.
#endif
!
#ifndef ISENTROPIC
  call equation_of_state ( statel, hydromtrl, density_in=dens, pt_out=pt, conservative=cons )
#else
  pt = den1*statel%es
#endif
!
  bb=0.0
  sumh = 0.
  do h = 1, nhydro
    sumh = sumh + hydromtrl(h)%q
  enddo
!
  do k = 1, nz
    bb(:,:,k) = dens(:,:,k)*g*( pt(:,:,k)*( 1. + epsm*den1(:,:,k)*(statel%qt(:,:,k) - sumh(:,:,k))    &
    	      - den1(:,:,k)*sumh(:,:,k) ) / (pt0(k)*(1. + epsm*qv0(k))) - 1. )
  enddo  
!
  do k = 2, nz
    buoy(:,:,k) = 0.5*(bb(:,:,k-1) + bb(:,:,k))
  enddo
  buoy(:,:,1) = bb(:,:,1)
!
#endif
!
return
end
!
!	==========================================	
  subroutine get_temperature_3d ( pressurel, statel, hydromtrl, tem )
!	==========================================	

  integer :: k
  type (atm_state), intent(in) :: statel
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
!
  real, dimension(ip_start:,jp_start:,:), intent(out) :: tem
!
!----------------------------------------------------------!
!
!  Isentropic case
!
#ifdef ISENTROPIC
!
  do k = 1, nz
#ifdef ANELASTIC
    tem(:,:,k) = statel%es(:,:,k)*(p0(k) / pref)**( (cp_a-cv_a)/cp_a )
#else
    tem(:,:,k) = statel%es(:,:,k)*((pressurel%p(:,:,k) + p0(k)) / pref)**( (cp_a-cv_a)/cp_a )
#endif  
  enddo
!
#else
!
!  Prognostic MSE
!
  call get_cp_cv ( statel%qt, hydromtrl )
!
  do k = 1, nz
    tem(:,:,k) = statel%es(:,:,k) - g*z0(k) - flv00*statel%qt(:,:,k)  
  enddo
!
  if (lmicro > 0) tem = tem + flv00*(hydromtrl(drop)%q + hydromtrl(rain)%q)
  if (lmicro > 1) tem = tem + fls00*hydromtrl(ice)%q
  if (lmicro > 2) tem = tem + fls00*(hydromtrl(grau)%q + hydromtrl(snow)%q)
  if (lmicro > 3) tem = tem + fls00*hydromtrl(hail)%q
!
  tem = 273.15 + tem / thermo_prop%cp
!
#endif
!
return 
end subroutine get_temperature_3d
!
!	==========================================	
  subroutine get_pt ( pressurel, statel, hydromtrl, pt )
!	==========================================	

  integer :: k
  type (atm_state), intent(in) :: statel
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:,jp_start:,:), intent(out) :: pt
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: exn
!
!----------------------------------------------------------!
!
!  Isentropic case
!
#ifdef ISENTROPIC
!
  pt = statel%es
!
#else
!
!  Prognostic MSE
!
#ifdef ANELASTIC
  do k = 1, nz
    exn(:,:,k) = (p0(k) / pref)**( (cp_a-cv_a)/cp_a )
  enddo
#else
  exn = (pressurel%p / pref)**( (cp_a-cv_a)/cp_a )
#endif
!
  call get_cp_cv ( statel%qt, hydromtrl )
!
  do k = 1, nz
    pt(:,:,k) = statel%es(:,:,k) - flv00*statel%qt(:,:,k) - g*z0(k)  
  enddo
!
  if (lmicro > 0) pt = pt + flv00*(hydromtrl(drop)%q + hydromtrl(rain)%q)
!  
  if (lmicro > 1) pt = pt + fls00*hydromtrl(ice)%q
!  
  if (lmicro > 2) pt = pt + fls00*(hydromtrl(grau)%q + hydromtrl(snow)%q)
!
  if (lmicro > 3) pt = pt + fls00*hydromtrl(hail)%q
!
  pt = pt / thermo_prop%cp + 273.15
!
  pt = pt / exn
!
#endif
!
return 
end subroutine get_pt
!
!	==========================================	
  subroutine get_ptv ( pressurel, statel, hydromtrl, ptv )
!	==========================================	

  integer :: k
  type (atm_state), intent(in) :: statel
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:,jp_start:,:), intent(out) :: ptv
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: exn, pt
!
!----------------------------------------------------------!
!
!  Isentropic case
!
#ifdef ISENTROPIC
!
  pt = statel%es
!
#else
!
!  Prognostic MSE
!
#ifdef ANELASTIC
  do k = 1, nz
    exn(:,:,k) = (p0(k) / pref)**( (cp_a-cv_a)/cp_a )
  enddo
#else
  exn = (pressurel%p / pref)**( (cp_a-cv_a)/cp_a )
#endif
!
  call get_cp_cv ( statel%qt, hydromtrl )
!
  do k = 1, nz
    pt(:,:,k) = statel%es(:,:,k) - flv00*statel%qt(:,:,k) - g*z0(k)  
  enddo
!
  if (lmicro > 0) pt = pt + flv00*(hydromtrl(drop)%q + hydromtrl(rain)%q)
!  
  if (lmicro > 1) pt = pt + fls00*hydromtrl(ice)%q
!  
  if (lmicro > 2) pt = pt + fls00*(hydromtrl(grau)%q + hydromtrl(snow)%q)
!
  if (lmicro > 3) pt = pt + fls00*hydromtrl(hail)%q
!
  pt = pt / thermo_prop%cp + 273.15
!
  pt = pt / exn
!
#endif
!
!  Virtual potential temperature
!
  ptv = pt*( 1. + epsm*statel%qt )
  if (lmicro > 0) ptv = ptv - pt*(1. + epsm)*( hydromtrl(drop)%q + hydromtrl(rain)%q )
  if (lmicro > 1) ptv = ptv - pt*(1. + epsm)*hydromtrl(ice)%q
  if (lmicro > 2) ptv = ptv - pt*(1. + epsm)*( hydromtrl(grau)%q + hydromtrl(snow)%q )
  if (lmicro > 3) ptv = ptv - pt*(1. + epsm)*hydromtrl(hail)%q
!
return 
end subroutine get_ptv
!
!	==========================================	
  subroutine get_mse ( pressurel, statel, hydromtrl, mse )
!	==========================================	

  integer :: k
  type (atm_state), intent(in) :: statel
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:,jp_start:,:), intent(out) :: mse
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: exn
!
!----------------------------------------------------------!
!
!  Isentropic case
!
#ifndef ISENTROPIC
!
  mse = statel%es
!
#else
!
!  Prognostic MSE
!
#ifdef ANELASTIC
  do k = 1, nz
    exn(:,:,k) = (p0(k) / pref)**( (cp_a-cv_a)/cp_a )
  enddo
#else
  exn = (pressurel%p / pref)**( (cp_a-cv_a)/cp_a )
#endif
!
  call get_cp_cv ( statel%qt, hydromtrl )
!
  do k = 1, nz
    mse(:,:,k) = thermo_prop%cp(:,:,k)*(exn(:,:,k)*statel%es(:,:,k) - 273.15) + flv00*statel%qt(:,:,k) + g*z0(k)
  enddo
!
  if (lmicro > 0) mse = mse - flv00*(hydromtrl(drop)%q + hydromtrl(rain)%q)
!
  if (lmicro > 1) mse = mse - fls00*hydromtrl(ice)%q
!
  if (lmicro > 2) mse = mse - fls00*(hydromtrl(grau)%q + hydromtrl(snow)%q)
!
  if (lmicro > 3) mse = mse - fls00*hydromtrl(hail)%q
!
#endif
!
return 
end subroutine get_mse
!
!	==========================================	
  subroutine get_qv ( qt, hydromtrl, qv, sumh )
!	==========================================	
!  
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: qt
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: qv
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout), optional :: sumh
!
!----------------------------------------------------------!
!
  if (present(sumh)) then
    sumh = 0.
    if (lmicro > 0) sumh = sumh + (hydromtrl(drop)%q+hydromtrl(rain)%q)
    if (lmicro > 1) sumh = sumh + hydromtrl(ice)%q
    if (lmicro > 2) sumh = sumh + (hydromtrl(grau)%q + hydromtrl(snow)%q)
    if (lmicro > 3) sumh = sumh + hydromtrl(hail)%q
    qv = qt - sumh
  else
    qv = qt
    if (lmicro > 0) qv = qv - (hydromtrl(drop)%q+hydromtrl(rain)%q)
    if (lmicro > 1) qv = qv - hydromtrl(ice)%q
    if (lmicro > 2) qv = qv - (hydromtrl(grau)%q + hydromtrl(snow)%q)
    if (lmicro > 3) qv = qv - hydromtrl(hail)%q
  endif
!
return 
end subroutine get_qv
!
!	==========================================	
  subroutine get_tdew ( tem, pressurel, statel, hydromtrl, tdew )
!	==========================================	

  integer :: i,j,k
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: rh
!  
  type (atm_state), intent(in) :: statel
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: tem
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: tdew
!
!----------------------------------------------------------!
!
      call get_rh ( tem, pressurel%p, statel, hydromtrl, rh )
!
      do k = 1, nz
	do j = jt_start, jt_end
	  do i = it_start, it_end
	    tdew(i,j,k) = cal_tdew ( tem(i,j,k), rh(i,j,k) )
	  enddo
	enddo
      enddo
!
return 
end subroutine get_tdew
!
!	==========================================	
  subroutine get_supersaturation ( pressurel, statel, hydromtrl, sat )
!	==========================================	

  integer :: i,j,k
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: qv, tem
!  
  type (atm_state), intent(in) :: statel
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: sat
!
!----------------------------------------------------------!
!
  call get_qv ( statel%qt, hydromtrl, qv )
  call get_temperature ( pressurel, statel, hydromtrl, tem )
!
#ifdef ANELASTIC
      do k = 1, nz
	sat(:,:,k) = qv(:,:,k)*cal_qsw1_2d(tem(:,:,k), p0(k)) - 1.
      enddo
#else
      do k = 1, nz
	sat(:,:,k) = qv(:,:,k)*cal_qsw1_2d(tem(:,:,k), p(:,:,k)) - 1.
      enddo
#endif
!
return 
end subroutine get_supersaturation
!
!	==========================================	
  subroutine get_rh ( tem, p, statel, hydromtrl, rh )
!	==========================================	

  integer :: i,j,k
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: qv
!  
  type (atm_state), intent(in) :: statel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: tem, p
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: rh
!
!----------------------------------------------------------!
!
      call get_qv ( statel%qt, hydromtrl, qv )
!
#ifdef ANELASTIC
      do k = 1, nz
	rh(:,:,k) = 100. * qv(:,:,k)*cal_qsw1_2d (tem(:,:,k), p0(k))
      enddo
#else
      do k = 1, nz
	rh(:,:,k) = 100. * qv(:,:,k)*cal_qsw1_2d (tem(:,:,k), p(:,:,k))
      enddo
#endif
!
return 
end subroutine get_rh
!
!	==========================================	
  subroutine get_rhi ( tem, p, statel, hydromtrl, rhi )
!	==========================================	

  integer :: i,j,k
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: qv
!  
  type (atm_state), intent(in) :: statel
  type (hydrometeor), dimension(:), intent(in) :: hydromtrl
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: tem, p
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: rhi
!
!----------------------------------------------------------!
!
      call get_qv ( statel%qt, hydromtrl, qv )
!
#ifdef ANELASTIC
      do k = 1, nz
	rhi(:,:,k) = 100. * qv(:,:,k)*cal_qsi1_2d (tem(:,:,k), p0(k))
      enddo
#else
      do k = 1, nz
	rhi(:,:,k) = 100. * qv(:,:,k)*cal_qsi1_2d (tem(:,:,k), p(:,:,k))
      enddo
#endif
!
return 
end subroutine get_rhi
!
	function cal_qsi1_2d( T, p ) result(qsi1)

	real :: T(ip_start:ip_end,jp_start:jp_end), qsi1(ip_start:ip_end,jp_start:jp_end)
	real :: p, ei
        integer :: i, j
	
	do j = jp_start, jp_end
	  do i = ip_start, ip_end
	    ei = cal_esi ( T(i,j) )
	    qsi1(i,j) = (p - ei) / (eps * ei) 
	  enddo
	enddo
	end function  
	
	function cal_qsw1_2d( T, p ) result(qsw1)

	real :: T(ip_start:ip_end,jp_start:jp_end), qsw1(ip_start:ip_end,jp_start:jp_end)
	real	:: p, es
        integer :: i, j

	do j = jp_start, jp_end
	  do i = ip_start, ip_end
	    es = cal_esw ( T(i,j) )
	    qsw1(i,j) = (p - es) / (eps * es)
	  enddo
	enddo
	end function
	  
  end module thermodynamics
