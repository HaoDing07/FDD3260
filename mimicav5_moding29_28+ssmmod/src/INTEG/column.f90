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
!  time_step.f90:                   
!
!  Purpose:
!	Used by modelctl.f90: contains the main routines called to 
!       integrate the dynamics
!
!  Author
!	Julien Savre
!	MIM, Ludwig Maximilian Universiat, Munich
!
! ================================================================

  module column_step
!
!----------------------------------------------------------!
!
  USE shared_all
  USE shared_tend
  USE columnmod
!  
  USE micropack
  USE micro_diagnostics
  USE radiationmod
  USE thermodynamics
  USE funcpack
  USE integpack
  USE allocation
  USE diagnostics
!
  IMPLICIT NONE
!
  private

  public :: columnstep

  CONTAINS

! ===========================================
  subroutine columnstep (per_radstep)
! ===========================================
!
  integer :: i, k, h, im
  real    :: xx, zz
  logical :: per_radstep
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: wold
!
#include "chemdef.h"
!
  if (verbose > 0) call write_debug('Start stepping')
!
!----------------------------------------------------------!
!                 Various initialisations                  !
!----------------------------------------------------------!
!
  dt = dt0
!
!  Reset diagnostics
!
  if (ldiag) call reset ( diag )
!
!  Reset stored tendencies
!
  if ( sorder > 1 ) then
    nstep = 3
    tend = 0.
    htend = 0.
#ifdef AERO_ENABLE
    atend = 0.
#else
    ntend = 0.
#endif
#ifdef NUC_CNT
    itend = 0.
#endif
  else
    nstep = 1
  endif
!
!  Initialise local variables
!
  wind%w(1,1,:) = w0
!
  pressure%p = pressure2%p
!
  pressure%dens = pressure2%dens
!
  state%es = state2%es
!
  state%qt = state2%qt
!
  state%scal = state2%scal
!
  hydromtr = hydromtr2
!
!----------------------------------------------------------!
!                     Integrate column                     !
!----------------------------------------------------------!
!
!
!  Radiative flux
!
  if (with_rad .and. per_radstep) then
#ifdef RAD_ENABLE
    call rad_front
#else
    call rad_forc
#endif
  endif
!
!  Conserved variables
!  
#if (defined CONSERVATIVE)
  call get_conserved ( pressure, state, hydromtr )
#endif
!
  do im = 1, nstep
    states = 0.
    hydromtrs = 0.
!
    call advance_column
!
    call integrate_all ( im )
  enddo
!
#if (defined CONSERVATIVE)
  call get_scalars( pressure, state, hydromtr )
#endif
!
!----------------------------------------------------------!
!                 Diagnostic microphysics                  !
!----------------------------------------------------------!
!
  if ( with_mic .and. lmicro > 0 ) then
    call diagnostic_micro ( pressure, wind, state, hydromtr, aero3d )
  endif
!
!-----------------------------------------------------------!
!		    Additional diagnostics	            !
!-----------------------------------------------------------!
!
  call thermo_diagnostics
!
  call special_diagnostics
!
!----------------------------------------------------------!
!   	      		 Updates		     	   !
!----------------------------------------------------------!
!
  call update_prognostic
!
!----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate stepping')
!
return
end subroutine columnstep
!
! ===============================================
  subroutine integrate_all ( im )
! ===============================================
!
! --- Update velotities using new pressure
!
  integer  :: k, im
!
  type(atm_winds) :: windold
  type(atm_pressure) :: pressureold
  type(atm_state) :: stateold
  type(hydrometeor), dimension(1:nhydro) :: hydromtrold
!
  if (verbose > 1) call write_debug('Start integrate_all')
!
!----------------------------------------------------------!
!           Density weighted variables at time t           !
!----------------------------------------------------------!
!
!  Allocate
!
  call alloc ( stateold )
  call alloc ( hydromtrold )
!
!  Natural variables at time t
!
  stateold%es = state2%es
!
  stateold%qt = state2%qt
!
  hydromtrold = hydromtr2
!
  stateold%scal = state2%scal
!
#if (defined CONSERVATIVE)
  call get_conserved ( pressure, stateold, hydromtrold )
#endif
!
!-----------------------------------------------------------!
!                       Time integration                    !
!-----------------------------------------------------------!
! 
  call integ_scal ( im, dt0, stateold=stateold, hydromtrold=hydromtrold,       &
   	 	    statel=state, hydromtrl=hydromtr, statea=states, hydromtra=hydromtrs )
!
!  Deallocate
!
  call dealloc ( stateold )
  call dealloc ( hydromtrold )
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate integrate_all')
!
return
end
!
! ===============================================
  subroutine update_prognostic 
! ===============================================
!
! --- Update velotities using new pressure
!
  integer  :: h
!
  if (verbose > 1) call write_debug('Start update_prognostic')
!
!-----------------------------------------------------------!
!                   Update data structures                  !
!-----------------------------------------------------------!
!
  wind2%w = wind%w
!
  pressure2%p = pressure%p
!
  pressure2%dens = pressure%dens
!
  state2%es = state%es
!
  state2%qt = state%qt
!
  hydromtr2 = hydromtr
!
  do h = 1, nscal
    state2%scal(:,:,:,h) = state%scal(:,:,:,h)
  enddo
!
#ifdef ISENTROPIC
  call equation_of_state ( state2, hydromtr2, pressure_in=pressure2%p, thermo_out=thermo, conservative=.false. )
#else
  call equation_of_state ( state2, hydromtr2, pressure_in=pressure2%p, thermo_out=thermo, conservative=.false. )
#endif 
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate update_prognostic')
!
return
end
!
end module column_step
