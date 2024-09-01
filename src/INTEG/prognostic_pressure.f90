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
!  prognostic_pressure.f90:
!
!  Purpose:
!	Calculates tendencies for a prognostic pressure equation.
!   This includes an advection term and velocity divergence term.
!
!  Author:
!	Julien Savre, Meteorology Institute, LMU Munich
!
! ================================================================
  
  module pressure_prognos
        
!
!-----------------------------------------------------------!
!
USE shared_all
USE allocation
USE gradients
USE advection
USE integpack
USE funcpack
USE averages
USE thermodynamics
USE boundary_conditions
!
IMPLICIT NONE
  
  logical :: limit_p=.false.
  real, parameter :: adiv = 0.2

  private

  public :: prognostic_pressure

  contains

!	===================================================
  subroutine prognostic_pressure ( im, windold )						   		  
!	===================================================
!
  integer, intent(in) :: im
  type(atm_winds), intent(in) :: windold
!  
  integer :: i, j, k, it
  real :: tau1, damp(nz)
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: gpx, gpy, gpz
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: zero, div, gdivx, gdivy, gdivz
  
  type(atm_winds) :: wind_tmp, winda, windv
  type(atm_pressure) :: pressurea
  type(atm_state) :: statea
  type(hydrometeor), dimension(1:nhydro) :: hydromtra
  real, dimension(:,:,:), allocatable :: buoy
!
  if (verbose > 0) call write_debug('Enter prognostic_pressure')
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!
!
  call alloc ( wind_tmp )
  call alloc ( winda )
  call alloc ( windv )
  call alloc ( pressurea )
  call alloc ( statea )
  call alloc ( hydromtra )
  call alloc ( buoy )
!
  tau1 = 1. / real(ntau)
  dtau = dt0*tau1
  zero = 0.
!
!----------------------------------------------------------!
!               Start small time-step loop                 !
!----------------------------------------------------------!
!
  do it = 1, ntau		!<----------------------------------------------------------
!
    wind_tmp = wind
    pressurea = pressures
    statea = states
    winda = 0.
    buoy = 0.
!
!  Get velocity difference
!
    call get_velocities( pressure, wind_tmp )
!
    windv%u = wind_tmp%u - wind2%u
#ifdef MODEL_3D
    windv%v = wind_tmp%v - wind2%v
#endif
    windv%w = wind_tmp%w - wind2%w
!
!----------------------------------------------------------!
!            Time-splitting: Solve for scalars             !
!----------------------------------------------------------!
!
!  Pressure gradient term for pressure and energy 
!
#ifndef ISENTROPIC
    if (with_adv) call div_pressure ( pressure, wind_tmp, pressurea%p, statea%es )
!
    if (out_diagt) diag(4)%pt = diag(4)%pt + cdiag*statea%es*tau1
#endif
!
!----------------------------------------------------------!
!         Advect and integrate energy and density          !
!     Advection is here performed without flux limiters    !
!----------------------------------------------------------!
!
    if (with_adv) then
!
      call advect ( dtau, dens=pressure%dens, statel=state, densa=pressurea%dens, statea=statea, &
      		    hydromtra=hydromtra, limiter=.false., full=.false. )
!    
      if (out_diagp) diag(8)%p = diag(8)%p + cdiag*pressurea%dens*tau1
      if (out_diagt) diag(9)%pt = diag(9)%pt + cdiag*statea%es*tau1
!
    endif
!
!----------------------------------------------------------!
!       	   Update full variables	           !
!----------------------------------------------------------!
!
!  Integrate all scalars
!
    call integ_scal ( im, dtau, pressureold=pressure, stateold=state, hydromtrold=hydromtr,     &
   	 	      pressurel=pressure, statel=state, hydromtrl=hydromtr,             	&
   	 	      pressurea=pressurea, statea=statea, hydromtra=hydromtrs )
!
!  Boundary conditions
!
    call pressurebc ( pressure )
!
    call statebc ( state )
!
    call hydrobc ( hydromtr )
!
    call equation_of_state ( state, hydromtr, density_in=pressure%dens, pressure_out=pressure%p, thermo_out=thermo, conservative=.true. )
!
!----------------------------------------------------------!
!            Time-splitting: Solve for momentum            !
!----------------------------------------------------------!
!
!  Pressure gradient
!
    call gradp ( pressure%p, wind%w, gpx, gpy, gpz )
!
!  Buoyancy
!    
    if (with_buoy) call buoyancy ( pressure%dens, state, hydromtr, buoy )  	
!
!  Divergence gradient for divergence damping
!
    call caldiv ( wind, div )
!
    call gradp ( div, zero, gdivx, gdivy, gdivz )
!
    call horav (thermo_prop%cs2, damp)
!
!  Assemble momentum tendencies
!
    do k = 1,nz
      do j=jt_start,jt_end
        do i=it_start,it_end
   	  winda%u(i,j,k) = winds%u(i,j,k) - gpx(i,j,k) + adiv*dtau*damp(k)*gdivx(i,j,k)
!
#ifdef MODEL_3D
   	  winda%v(i,j,k) = winds%v(i,j,k) - gpy(i,j,k) + adiv*dtau*damp(k)*gdivy(i,j,k) 
#endif
!     
   	  winda%w(i,j,k) = winds%w(i,j,k) - gpz(i,j,k) + adiv*dtau*damp(k)*gdivz(i,j,k) + buoy(i,j,k)
        enddo
      enddo
    enddo 
!
    winda%w(:,:,1) = 0. 
!
!  Diagnostics
!
    if (out_diagu) diag(6)%u = diag(6)%u - cdiag*gpx
    if (out_diagu) diag(8)%u = diag(8)%u + cdiag*adiv*dtau*damp(k)*gdivx(i,j,k)
#ifdef MODEL_3D
    if (out_diagv) diag(6)%v = diag(6)%v - cdiag*gpy
    if (out_diagv) diag(8)%v = diag(8)%v + cdiag*adiv*dtau*damp(k)*gdivy(i,j,k)
#endif
    if (out_diagw) diag(6)%w = diag(6)%w - cdiag*gpz
    if (out_diagw) diag(7)%w = diag(7)%w + cdiag*buoy*tau1
    if (out_diagw) diag(8)%w = diag(8)%w + cdiag*adiv*dtau*damp(k)*gdivz(i,j,k)
!
!  Integrate momentum
!
    if ( with_mom ) call integ_mom ( im, dtau, wind, wind, winda )
!
!  BC
!
    call windbc ( wind )
!  
enddo		!<------------------------------------------------------------------------
!
!  Diagnostics
!
  thermo%buoy = buoy
!
  if (out_diagu) diag(9)%u = diag(9)%u + cdiag*(wind%u - windold%u)/dt0
#ifdef MODEL_3D
  if (out_diagv) diag(9)%v = diag(9)%v + cdiag*(wind%v - windold%v)/dt0
#endif
  if (out_diagw) diag(9)%w = diag(9)%w + cdiag*(wind%w - windold%w)/dt0
!
  call dealloc ( wind_tmp )
  call dealloc ( winda )
  call dealloc ( windv )
  call dealloc ( pressurea )
  call dealloc ( statea )
  call dealloc ( hydromtra )
  call dealloc ( buoy )
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate prognostic_pressure')
!
return                  
end subroutine prognostic_pressure                                                   
!
!	===================================================
  subroutine div_pressure ( pressurel, windl, ps, ess )				     
!	===================================================
!
  integer :: i, j, k
!
  type(atm_pressure), intent(in) :: pressurel
  type(atm_winds), intent(inout) :: windl
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: ps, ess
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ugp, ugrcv, pp, rcv, div 
!
!----------------------------------------------------------!
!                  Calculate divergence                    !
!----------------------------------------------------------!
!
  ugp = 0.0
  ugrcv = 0.0
!
  do k = 1, nz
    pp(ip_start:ip_end,jp_start:jp_end,k) = pressurel%p(ip_start:ip_end,jp_start:jp_end,k) + p0(k)
!
#ifdef ISENTROPIC
    rcv(ip_start:ip_end,jp_start:jp_end,k) = (cp_a-cv_a) / cv_a
#else
    rcv(ip_start:ip_end,jp_start:jp_end,k) = (thermo_prop%cp(ip_start:ip_end,jp_start:jp_end,k)-thermo_prop%cv(ip_start:ip_end,jp_start:jp_end,k))  &
    	/ thermo_prop%cv(ip_start:ip_end,jp_start:jp_end,k)
#endif
  enddo
!
  call caldiv ( windl, div )
!
  call u_gradp ( windl, pp, ugp )
!
#ifndef ISENTROPIC
  if (cst_cp == 0) call u_gradp ( windl, rcv, ugrcv )
#endif
!
!----------------------------------------------------------!
!                   Pressure tendency                      !
!----------------------------------------------------------!
!
#ifndef CONSERVATIVE
  do j=jt_start,jt_end
    do i=it_start,it_end
      do k=1,nz
        ps(i,j,k) = ps(i,j,k) + rcv(i,j,k)*ugp(i,j,k) + pp(i,j,k)*ugrcv(i,j,k)
!
	if (out_diagp) diag(3)%p(i,j,k) = diag(3)%p(i,j,k) + cdiag*(rcv(i,j,k)*ugp(i,j,k) + pp(i,j,k)*ugrcv(i,j,k))*dtau/dt0
      enddo
    enddo
  enddo  
#endif
!
!----------------------------------------------------------!
!                     Energy tendency                      !
!----------------------------------------------------------!
!
#ifndef ISENTROPIC
  do j=jt_start,jt_end
    do i=it_start,it_end
      do k=1,nz
#ifdef CONSERVATIVE
	ess(i,j,k) = ess(i,j,k) + ugp(i,j,k)
	if (out_diagt) diag(8)%pt(i,j,k) = diag(8)%pt(i,j,k) + cdiag*ugp(i,j,k)*dtau/dt0
#else
        ess(i,j,k) = ess(i,j,k) - div(i,j,k) * pp(i,j,k) / pressurel%dens(i,j,k)
	if (out_diagt) diag(8)%pt(i,j,k) = diag(8)%pt(i,j,k) - cdiag*div(i,j,k)*pp(i,j,k) / pressurel%dens(i,j,k)
#endif
      enddo
    enddo
  enddo  
#endif
!
!  BC
!
  call statebc ( ps )
!
  call statebc ( ess )
!
!-----------------------------------------------------------!
!
return                  
end subroutine div_pressure

end module pressure_prognos                                                 
