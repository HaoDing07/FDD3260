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
!  UD.F:
!
!  Purpose:
!	Scalars pt, qv, q, n: All scalar tendencies are calculated first
!	including advection, sgs mixing, damping and microphysics,
!	then scalars are integrated using either forward euler or RK.
!	As an input, all tendency terms already contain SGS mixing tendencies.
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================

  module scalarsmod
  
  USE shared_all
  USE allocation
  USE advection
  USE micro
  USE radiationmod
  USE sources
  USE subgrid
  USE integpack
  USE funcpack
  USE boundary_conditions
!
  IMPLICIT NONE
!
  private
!
  public :: scalars

  CONTAINS

!	===================================================
  subroutine scalars     
!	===================================================
!    
  integer :: k, h
  character(len=1) :: car1
  real :: dref0(nz), es0(nz), zero(nz)
!
  type(atm_state) :: statea, statetmp
  type(atm_pressure) :: pressurea, pressuretmp
  type(hydrometeor), dimension(1:nhydro) :: hydromtra, hydromtrtmp
!
  if (verbose > 0) call write_debug('Entering scalars')
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!
!
!  Allocate and initialize
!
  call alloc ( statea )
  call alloc ( pressurea )
  call alloc ( hydromtra )
!
!  Reference soundings
!
#if (defined CONSERVATIVE)
  dref0 = den0
#else
  dref0 = 1.
#endif
!
#ifdef ISENTROPIC
  es0 = pt0
#else
  es0 = mse0
#endif
!
  zero = 0.
!
!----------------------------------------------------------!
!                        Sources                           !
!----------------------------------------------------------!
!
!  Turbulent diffusion    
!
  if ( with_dif ) call eff0 ( pressure2%dens, state2, hydromtr2 )
!
!  Radiation 
!
  if ( with_rad ) call rad_source ( pressure2%dens, states%es )
!
!  Nudging
!
  if ( with_nudg ) then
    call nudge_x ( 'pt ', dref0*es0, state%es, states%es )
!
    call nudge_x ( 'qt ', dref0*qt0, state%qt, states%qt )
  endif
!
!  LS advection
!
  if ( with_lsadv ) then
    call ls_advec ( 1, state%es, states%es )
!
    call ls_advec ( 2, state%qt, states%qt )
  endif
!
!  Other sources
!
  if ( with_lssrc ) then
    call other_source ( 1, pressure%dens, states%es )
!
    call other_source ( 2, pressure%dens, states%qt )
  endif
!
!  Large-scale subsidence and upwelling
!
  call advection_1d ( w0 )
!
!----------------------------------------------------------!
!                   Explicit microphysics                  !
!----------------------------------------------------------!
!
  if ( .not.split_mic .and. with_mic .and. lmicro > 0 ) then
    call microphysics ( pressure, wind, state, hydromtr, aero3d )
  endif
!
!----------------------------------------------------------!
!    		      Scalar advection	  		   !
!----------------------------------------------------------!
!
  if (with_adv) then
    call advect ( dt0, pressure%dens, statel=state, hydromtrl=hydromtr,   	&
      	          densa=pressurea%dens, statea=statea, hydromtra=hydromtra )
  endif
!
!----------------------------------------------------------!
!                    Update tendencies                     !
!----------------------------------------------------------!
!
  pressures%dens = pressurea%dens + pressures%dens
!
  states%es = statea%es + states%es
!
  states%qt = statea%qt + states%qt
!
  do h = 1, nhydro
    hydromtrs(h)%q = hydromtra(h)%q + hydromtrs(h)%q
!
    if (moments == 2) hydromtrs(h)%n = hydromtra(h)%n + hydromtrs(h)%n
!
    if ( lmicro > 3 .and. h > ice ) hydromtrs(h)%w = hydromtra(h)%w + hydromtrs(h)%w
  enddo
!
!----------------------------------------------------------!
!	     Additional sources and relaxations            !
!----------------------------------------------------------!
!
!  Apply sponge layers
!
  if (sponge == 1) then
    call damps ( 'pt', state%es, dref0*es0, states%es )
  endif
!
!----------------------------------------------------------!
!          Pressure contribution (for MSE only)            !
!----------------------------------------------------------!
!
#ifndef ISENTROPIC
  do k = 1, nz-1
    state%es(:,:,k) = state%es(:,:,k) + 0.5*dt0*stab0(k)*thermo_prop%cs2(:,:,k)*(wind%w(:,:,k+1) + wind%w(:,:,k))
  enddo
!
  if (out_diagt) then
    do k = 1, nz-1
      diag(8)%pt(:,:,k) = diag(8)%pt(:,:,k) + 0.5*cdiag*stab0(k)*thermo_prop%cs2(:,:,k)*(wind%w(:,:,k+1) + wind%w(:,:,k))
    enddo
  endif
#endif
!
!----------------------------------------------------------!
!               	  Terminate   	                   !
!----------------------------------------------------------!
!
!  1 moment microphysics
!
  if ( lmicro > 0 .and. lndrop == 0 ) hydromtrs(drop)%n = 0.0 
!
!  Diagnostics
!  
  if (ldiag) then
    if (out_diagt) diag(9)%pt = diag(9)%pt + cdiag*states%es
    if (out_diagq) diag(9)%qt = diag(9)%qt + cdiag*states%qt
    if (out_diagl.and.lmicro>0) diag(9)%qc = diag(9)%qc + cdiag*hydromtrs(drop)%q
    if (out_diagr.and.lmicro>0) diag(9)%qr = diag(9)%qr + cdiag*hydromtrs(rain)%q
    if (out_diagi.and.lmicro>1) diag(9)%qi = diag(9)%qi + cdiag*hydromtrs(ice)%q
  endif
!
!  Deallocate
!
  call dealloc ( statea ) 
  call dealloc ( pressurea ) 
  call dealloc ( hydromtra )
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating scalars')
!
return                  
end
!
end module scalarsmod
