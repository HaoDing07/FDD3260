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

  module time_step

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
  USE radiationmod
  USE thermodynamics
  USE sources
  USE funcpack
  USE integpack
  USE allocation
  USE averages
!
#ifdef SPMD
  USE mpi
#endif
!
  IMPLICIT NONE
!
  private

  public :: stepping

  CONTAINS

! ===========================================
  subroutine stepping (per_radstep)
! ===========================================
!
  integer :: i, j, k, h, im
  real    :: x, xx, zz
  real, dimension(1:nz) :: ptvav
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: buoy, ptv, windold
  logical :: per_radstep
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
  if ( sorder > 1 ) then
    cdiag = 0.5
  else
    cdiag = 1.
  endif
!
  if ( split_mic ) then
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
!  Fixed momentum field for spiral test case
!
  if (trim(casename) == 'SPIRAL') then
    zz = 0.
    do k = 3, nz-2
      xx = 0.
      do i = it_start, it_end
        wind2%u(i,:,k) = sin(pi*(xx-dx/2.))*sin(pi*(xx-dx/2.))*sin(2.*pi*zz)*cos(pi*time/5.)
        wind2%w(i,:,k) = -sin(pi*(zz-dz/2.))*sin(pi*(zz-dz/2.))*sin(2.*pi*xx)*cos(pi*time/5.)
        xx = xx + dx
      enddo
      zz = zz + dz*fdz0(k)
    enddo
    call windbc ( wind2 )
  endif
!
!  CP intensity
!  
  if (out_cp) then
    column%cpint0 = 0.0
!
    call get_ptv ( pressure, state, hydromtr, ptv )
    call horav ( ptv, ptvav )
!
    buoy = 0.
    do k = 1, nz
      buoy(:,:,k) = g*(ptv(:,:,k) - ptvav(k)) / ptvav(k)
    enddo
!
    do j = jp_start, jp_end
      do i = ip_start, ip_end
	zz = 0.
        k = 1
	do while ( zz <= zcp )
          column%cpint(i,j) = column%cpint(i,j) - 2.*buoy(i,j,k)*dz*fdz0(k)
	  zz = zz + dz*fdz0(k)
          k = k + 1
	enddo
        column%cpint0(i,j) = sqrt( max(column%cpint0(i,j),0.) )
      enddo
    enddo
  endif
!
! Data at time level n
!
  pressure%p = pressure2%p
!
#ifdef ANELASTIC
  pressure%dens = pressure2%dens
#else
  do k = 1, nz
    pressure%dens(:,:,k) = den0(k)
  enddo
#endif
!
  wind = wind2
!
  state = state2
!
  hydromtr = hydromtr2
!
  nuc = nuc2
#ifdef AERO_ENABLE
  aero3d = aero3d2
#endif
#ifdef CNT_NUC
  nucin = nucin2
#endif
!
!  Reset diagnostics
!
  if (ldiag) call reset ( diag )
!
  if (ldiag) call reset ( column )
!
!----------------------------------------------------------!
!            Integrate microphysics and Dynamics           !
!----------------------------------------------------------!
!
!----------------------------------------------------------!
!                   First step: physics                    !
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
!  SGS tendencies
!
  if ( with_dif ) call subgrid_interface ( pressure, wind, state, hydromtr )
!
!----------------------------------------------------------!
!                  Second step: dynamics                   !
!----------------------------------------------------------!
!
!  Density weighting
!
  call get_momentum ( pressure, wind )
!
#if (defined CONSERVATIVE)
  call get_conserved ( pressure, state, hydromtr, with_aero=.true. )
#endif
!
!  Dynamics loop
!
  do im = 1, nstep
!
!  Initialize tendencies
!
    windold = wind%w
    winds = 0.
    states = 0.
    pressures = 0.
    hydromtrs = 0.
#ifdef AERO_ENABLE
    call reset ( aero3ds )
#endif
#ifdef NUC_CNT
    call reset ( nucins )
#endif
!
!  Get advective velocity
!
    if ( with_adv ) call get_advective ( pressure, wind )
!
!  Scalars
!
    if ( with_scal ) call scalars
!
    if ( with_scal .and. nscal > 0 ) call tracers
!
!  Aerosols
!
#ifdef SEIFERT
    if ( aero_flg%tran ) call advc_m ( im )
#endif
!
!  Momentum
!
    if ( with_mom ) call momentum
!
!  Integrate
!  
    call integrate_all ( im )
!
!  Solve for anelastic pressure 
!
#ifdef ANELASTIC
    if ( with_mom .and. nsubp > 0 ) call anelastic_pressure ( windold )
#endif
!
  enddo
!
!  Invert density weighting
!
  call get_velocities( pressure, wind )
!
#if (defined CONSERVATIVE)
  call get_scalars( pressure, state, hydromtr, with_aero=.true. )
#endif
!
!----------------------------------------------------------!
!                        Microphysics                      !
!----------------------------------------------------------!
!
  if ( with_mic .and. lmicro > 0 ) then
!
!  Split microphysics
!
    if ( split_mic ) call microphysics ( pressure, wind, state, hydromtr, aero3d )
!
!  Diagnostic microphysics
!  
    call diagnostic_micro ( pressure, wind, state, hydromtr, aero3d )
!
!  One moment
!
#ifdef SEIFERT
    if ( moments == 1 ) call one_moment ( hydromtr )
#endif
  endif
!
!----------------------------------------------------------!
!   	                 Updates 			   !
!----------------------------------------------------------!
!
!  Update
!
  call update_prognostic
!
!  Reset scalars 
!
  do h = 1, nscal
    call scalar_reset ( h, state2%scal(:,:,:,h), zbl, (sca_set==1.or.sca_set==3) )
  enddo
!
!----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate stepping')
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
  integer  :: k, im
  character(len=1) :: c(5)
!
  type(atm_winds) :: windold
  type(atm_pressure) :: pressureold
  type(atm_state) :: stateold
  type(hydrometeor), dimension(1:nhydro) :: hydromtrold
!
  if (verbose > 1) call write_debug('Start integrate_all')
!
!  Allocate
!
  call alloc ( windold )
  call alloc ( pressureold )
  call alloc ( stateold )
  call alloc ( hydromtrold )
!
!----------------------------------------------------------!
!           Density weighted variables at time t           !
!----------------------------------------------------------!
!
!  Natural variables at time t
!
  pressureold = pressure2
!
  stateold = state2
!
  hydromtrold = hydromtr2
!
  windold = wind2
!
!  Density weighting
!
  call get_momentum ( pressureold, windold )
!
#if (defined CONSERVATIVE)
  call get_conserved ( pressureold, stateold, hydromtrold )
#endif
!
!-----------------------------------------------------------!
!               Prognostic pressure integration             !
!-----------------------------------------------------------!
!
#ifndef ANELASTIC
!
  call prognostic_pressure ( im, windold )
!
#else
!
!-----------------------------------------------------------!
!               Diagnostic pressure integration             !
!-----------------------------------------------------------!
!
  if (with_mom) call integ_mom ( im, dt0, windold, wind, winds )
! 
  call integ_scal ( im, dt0, pressureold=pressureold, stateold=stateold, hydromtrold=hydromtrold,       &
   	 	    pressurel=pressure, statel=state, hydromtrl=hydromtr,                               &
   	 	    pressurea=pressures, statea=states, hydromtra=hydromtrs )
!
#endif
!
!  BCs 
!
  if ( sorder > 1 ) then
    call windbc ( wind )
!
    call pressurebc ( pressure )
!
    call statebc ( state )
!
    call hydrobc ( hydromtr )
  endif
!
!  Limit hydrometeors
!
  if ( lmicro > 0 ) call limit_hydro ( state%qt, state%es, hydromtr, aero3d )
!
!  For debugging purposes
!
!  c = (/'C','R','I','G','S'/)
!  call debug_minmax( wind%u, 'U' )
!  call debug_minmax( wind%v, 'V' )
!  call debug_minmax( wind%w, 'W' )
!  call debug_minmax( state%qt, 'QT' )
!  call debug_minmax( state%qt, 'QT' )
!  do k = 1, nhydro
!    call debug_minmax( hydromtr(k)%q, 'Q'//c(k) )
!    call debug_minmax( hydromtr(k)%n, 'N'//c(k) )
!  enddo
!
!  Deallocate
!
  call dealloc ( windold )
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
  subroutine update_prognostic 
! ===============================================
!
! --- Update velotities using new pressure
!
  integer  :: k, h
!
  if (verbose > 1) call write_debug('Start update_prognostic')
!
!-----------------------------------------------------------!
!                    Boundary conditions                    !
!-----------------------------------------------------------!
!
  call windbc ( wind )
!
  call pressurebc ( pressure )
!
  call statebc ( state )
!
  call hydrobc ( hydromtr )
!
#ifdef AERO_ENABLE
  call aerobc ( aero3d )
#endif
!
!-----------------------------------------------------------!
!                  Data at time level n                     !
!-----------------------------------------------------------!
!
  pressure2 = pressure
!
  wind2 = wind
!
  state2 = state
!
  hydromtr2 = hydromtr
!
  nuc2 = nuc 
#ifdef AERO_ENABLE
  aero3d2 = aero3d
#endif
#ifdef CNT_NUC
  nucin2 = nucin
#endif
!
!-----------------------------------------------------------!
!                      Thermodynamics                       !
!-----------------------------------------------------------!
!
  call equation_of_state ( state2, hydromtr2, pressure_in=pressure2%p, thermo_out=thermo, conservative=.false. )
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate update_prognostic')
!
return
end

end module time_step
