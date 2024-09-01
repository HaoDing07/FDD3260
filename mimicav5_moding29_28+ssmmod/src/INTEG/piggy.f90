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
!  piggy.f90:                   
!
!  Purpose:
!	Interface for piggy-backing: a different microphysics is solved
!	with main dynamics
!
!  Author
!	Julien Savre
!	MIM, Ludwig Maximilian Universiat, Munich
!
! ================================================================

  module piggyback

!
!----------------------------------------------------------!
!
  USE shared_all
  USE shared_tend
  USE momentummod
  USE scalarsmod
  USE tracersmod
  USE subgrid
  USE micro
  USE micropack
  USE micro_diagnostics
  USE pressure_solver
  USE pressure_prognos
  USE boundary_conditions
  USE thermodynamics
  USE radiationmod
  USE sources
  USE funcpack
  USE integpack
  USE allocation
  USE diagnostics
  USE netcdfmod
  USE netcdfslice
!
  IMPLICIT NONE
!
!  Space to define default micro parameters
!
  integer :: lfreeze0
  real :: xauto0 
!
  private

  public :: piggy

  CONTAINS

! ===========================================
  subroutine piggy (per_radstep)
! ===========================================
!
  integer :: im
  logical :: per_radstep
!
  if (verbose > 0) call write_debug('Start piggy')
!
!----------------------------------------------------------!
!                 Various initialisations                  !
!----------------------------------------------------------!
!
  dt = dt0
!
!
  if (split_mic) then
    cint = 1.
    ! cint = 0.5
  else
    cint = cdiag
    ! cint = 1.
  endif
!
!  Reset stored tendencies
!
  if ( sorder > 1 ) then
    nstep = 2
    tend = 0.
    htend = 0.
  else
    nstep = 1
  endif
!
!  Reset diagnostics
!
  if (ldiag) call reset ( diag )
!
!  Piggy-backing: set alternative micro parameter
!
  call set_piggy
!
! Data at time level n
!
  pressure = pressure2
!
  state = statep
!
  hydromtr = hydromtrp
!
#ifdef AERO_ENABLE
  aero3d = aero3d2
#endif
!
!----------------------------------------------------------!
!            Integrate microphysics and Dynamics           !
!----------------------------------------------------------!
!
!  Radiative fluxes
!
  if (with_rad .and. per_radstep) then
#ifdef RAD_ENABLE
    call rad_front
#else
    call rad_forc
#endif
  endif
!
!  Calculate SGS tendencies
!
  if (with_dif) call subgrid_interface ( pressure, wind, statep, hydromtrp )
!
!  Density weighting
!
#if (defined CONSERVATIVE)
  call get_conserved ( pressure, state, hydromtr, with_aero=.true. )
#endif
!
  do im = 1, nstep
!
!  Initialize tendencies
!
    states = 0.
    hydromtrs = 0.
!
    wind_adv = wind2
!
!  Scalar tendencies
!
    if (with_scal) call scalars	
!
!  Integrate
!  
    call integrate_all ( im )
!
  enddo
!
#if (defined CONSERVATIVE)
  call get_scalars( pressure, state, hydromtr, with_aero=.true. )
#endif
!
!  Limit hydrometeors
!
  call limit_hydro ( state%qt, state%es, hydromtr, aero3d )
!
!  Diagnostic microphysics step 
!
  if ( with_mic .and. lmicro > 0 ) then
    if ( split_mic ) call microphysics ( pressure, wind2, statep, hydromtrp, aero3d )
!
    call diagnostic_micro ( pressure, wind2, statep, hydromtrp, aero3d )
  endif
!
!----------------------------------------------------------!
!                         Updates                          !
!----------------------------------------------------------!
!
!  BCs 
!
  call statebc ( state )
!
  call hydrobc ( hydromtr )
!
! Data at time level n
!
  statep = state
!
  hydromtrp = hydromtr
!
!  Piggy-backing: reset alternative micro parameter
!
  call reset_piggy
!
!----------------------------------------------------------!
!                         Outputs                          !
!----------------------------------------------------------!
!
  call output_pig
!
!----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate piggy')
!
return
end
!
! ===============================================
  subroutine integrate_all ( im )
! ===============================================
!
! --- Update velotities using new pressure
!
  integer  :: im
  type(atm_pressure) :: pressureold
  type(atm_state) :: stateold
  type(hydrometeor), dimension(1:nhydro) :: hydromtrold
!
  if (verbose > 1) call write_debug('Start integrate_all')
!
!  Allocate
!
  call alloc ( pressureold )
  call alloc ( stateold )
  call alloc ( hydromtrold )
!
!----------------------------------------------------------!
!           Density weighted variables at time t           !
!----------------------------------------------------------!
!
  pressureold = pressure2
!
  stateold = statep
!
  hydromtrold = hydromtrp
!
#if (defined CONSERVATIVE)
  call get_conserved ( pressureold, stateold, hydromtrold )
#endif
!
!-----------------------------------------------------------!
!               Diagnostic pressure integration             !
!-----------------------------------------------------------!
!
  call integ_scal ( im, dt0, stateold=stateold, hydromtrold=hydromtrold,       &
   	 	    statel=state, hydromtrl=hydromtr, statea=states, hydromtra=hydromtrs )
!
!  BCs 
!
  if ( sorder > 1 ) then
    call statebc ( state )
    call hydrobc ( hydromtr )
  endif
!
!  Deallocate
!
  call dealloc ( pressureold )
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
  subroutine set_piggy
! ===============================================
!
!  Store reference micro parameters
!
  lfreeze0 = lfreeze
  !xauto0 = xauto
!
!  Define new micro parameters
!
  lfreeze = 1
  !xauto = 2.6e-10
!
return
end
!
! ===============================================
  subroutine reset_piggy 
! ===============================================
!
!  Reset to reference micro parameters
!
  lfreeze = lfreeze0
  !xauto = xauto0
!
return
end

end module piggyback
