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
!  microphysics.F:
!
!  Purpose:
!	Interface to calculate split microphysics
!
!  Author:
!	Julien Savre, MIM, Ludwig-Maximilian Universität, Munchen
!
! ================================================================

  module micro

!  
!  ---------------------------------------------
!
  USE shared_all
  USE kessler
  USE seifert_beheng_warm
  USE seifert_beheng_ice
  USE micropack
  USE micro_diagnostics
  USE thermodynamics
  USE precipitation
  USE aerosols
  USE allocation
  USE averages
  USE funcpack
  USE boundary_conditions

  ! USE micro_dropsed

#ifdef SPMD
  USE mpi
#endif

  IMPLICIT NONE
  
  private
  
  public :: microphysics
  
  contains
  
! ===================================================
  subroutine microphysics ( pressurel, windl, statel, hydromtrl, aero3dl )
! ===================================================
!
  type(atm_winds), intent(in) :: windl
  type(atm_state), intent(inout) :: statel 
  type(atm_pressure), intent(inout) :: pressurel
  type(hydrometeor), dimension(1:nhydro), intent(inout) :: hydromtrl 
  type(aero_3d), dimension(1:nmode), intent(inout)  :: aero3dl
!
  real, allocatable, dimension(:,:,:) :: dTdt, ww, dref
  type(atm_state) :: statem
  type(atm_pressure) :: pressurem
  type(hydrometeor), dimension(1:nhydro) :: hydromtrm
  type(aero_3d), dimension(:), allocatable  :: aero3dm
! 
  integer :: i, j, k, h
  real :: dtl, x, zz
!
#if ( defined SPMD )      
  integer :: ierr, status(MPI_STATUS_SIZE)
#endif
!
  if (verbose > 0) call write_debug('Enter microphysics')
!
!----------------------------------------------------------!
!                     initializations         	           !
!----------------------------------------------------------!
!
!  Allocate Tendencies
!
  call alloc ( dref )
  call alloc ( dTdt )
  call alloc ( ww )
!
!  Reference density 
!
#ifdef CONSERVATIVE
#if (!defined ANELASTIC)
  dref = pressurel%dens
#else
  do k = 1, nz
    dref(:,:,k) = den0(k)
  enddo
#endif
#else
  dref = 1.
#endif
!
!  Update thermo
!
#if (defined CONSERVATIVE)
  if ( .not.split_mic ) call get_scalars( pressurel, statel, hydromtrl, with_aero=.true. )
#endif
!
  call equation_of_state ( statel, hydromtrl, pressure_in=pressurel%p, thermo_out=thermo, conservative=.false. )
!
!  Temperature tendencies for condensation
!
  if ( .not.split_mic ) then
    do k=1,nz-1
      do j=jp_start,jp_end
        do i=ip_start,ip_end
          ww(i,j,k) = 0.5*(windl%w(i,j,k+1) + windl%w(i,j,k)) / dref(i,j,k)
          dTdt(i,j,k) = states%es(i,j,k) / (dref(i,j,k)*thermo%exn(i,j,k))
        end do
      end do
    end do
  else
    ww = 0.
    dTdt = 0.
  endif 
!
!----------------------------------------------------------!
!                 Microphysics tendencies                  !
!----------------------------------------------------------!
!
!  More allocations
!
  call alloc ( statem )
  call alloc ( pressurem )
  call alloc ( hydromtrm )
#ifdef AERO_ENABLE
  allocate ( aero3dm(1:nmode) )
  call alloc ( aero3dm )
#endif
! 
!  Microphysics sources
! 
#ifdef SEIFERT
    if ( lmicro > 0 ) &
        call sb_warm ( dref, ww, dTdt, statel%qt, pressurel, hydromtrl, aero3dl, statem, hydromtrm, aero3dm ) 
!
    if ( lmicro > 1 .and. lfreeze > 0 ) &
        call sb_ice ( dref, ww, dTdt, statel%qt, pressurel, hydromtrl, aero3dl, statem, hydromtrm, aero3dm )
#else
    call micro_k( dref, pressurel, statel, hydromtrl, statem%es, hydromtrm )
#endif
! 
!  Precipitation
!
#ifdef SEIFERT
    call calc_precip ( dref, hydromtrl, pressurem%dens, statem%qt, hydromtrm, aero3dm ) 
#else
    call calc_precip_k ( dref, hydromtrl(rain)%q, pressurem%dens, statem%qt, hydromtrm ) 
#endif

!  Droplet_sedimentation
!
#ifdef SEIFERT
    call calc_sedimentation ( dref, hydromtrl, statem%qt, hydromtrm, dref)
#endif
!
!----------------------------------------------------------!
!       Update tendencies for explicit microphysics        ! 
!----------------------------------------------------------!
!
  if ( .not.split_mic ) then
!
#ifndef ANELASTIC
    pressures%dens = pressurem%dens + pressures%dens
#endif
!
    states%es = dref*statem%es + states%es
!
    states%qt = dref*statem%qt + states%qt
!
    do h = 1, nhydro
      hydromtrs(h)%q = dref*hydromtrm(h)%q + hydromtrs(h)%q
!
      if (moments == 2) hydromtrs(h)%n = dref*hydromtrm(h)%n + hydromtrs(h)%n
!
      if ( lmicro > 3 .and. h > ice ) hydromtrs(h)%w = dref*hydromtrm(h)%w + hydromtrs(h)%w
    enddo
!
#if (defined CONSERVATIVE)
    call get_conserved ( pressurel, statel, hydromtrl, with_aero=.true. )
#endif
!
!----------------------------------------------------------!
!          Else, update microphysical variables            ! 
!----------------------------------------------------------!
!
  else
!
#ifndef ANELASTIC
    pressurel%dens = pressurel%dens + dt0*pressurem%dens
#endif
!
#ifdef ISENTROPIC
    statel%es = statel%es + dt0*statem%es
#endif
!
    statel%qt = statel%qt + dt0*statem%qt
!
    do h = 1, nhydro
      hydromtrl(h)%q = hydromtrl(h)%q + dt0*hydromtrm(h)%q
!
#ifdef SEIFERT
      if ( moments == 2 .and. (lndrop == 1 .or. h > drop) ) hydromtrl(h)%n = hydromtrl(h)%n + dt0*hydromtrm(h)%n
!
      if ( lmicro > 3 .and. h > ice ) hydromtrl(h)%w = hydromtrl(h)%w + dt0*hydromtrm(h)%w
#endif
    enddo
!
#ifdef AERO_ENABLE
    if ( aero_flg%act_scv .or.  aero_flg%reg .or. aero_flg%imp_scv ) then
      do h = 1, nmode
        aero3dl(h)%n = aero3dl(h)%n + dt0*aero3dm(h)%n
!
        aero3dl(h)%m = aero3dl(h)%m + dt0*aero3dm(h)%m
!
        aero3dl(h)%ma = aero3dl(h)%ma + dt0*aero3dm(h)%ma
      enddo
    endif
#endif
!
!  Boundary conditions
!
    call statebc ( statel )
!
    call hydrobc ( hydromtrl )
!
#ifdef AERO_ENABLE
    call aerobc ( aero3dl )
#endif
!
!  Limit hydrometeors 
!
    call limit_hydro ( statel%qt, statel%es, hydromtrl, aero3dl )
!
  endif
!
!----------------------------------------------------------!
!                 	   Finalize   	                   !
!----------------------------------------------------------!
!
!  Store microphysics tendencies 
!
  if (out_diagl) diag(7)%qc = diag(7)%qc + dref*cint*hydromtrm(drop)%q
  if (out_diagr) diag(7)%qr = diag(7)%qr + dref*cint*hydromtrm(rain)%q
  if (out_diagi .and. lmicro > 1) diag(7)%qi = diag(7)%qi + dref*cint*hydromtrm(ice)%q
!
  if (out_diagt) diag(9)%pt = diag(9)%pt + cint*dref*statem%es
  if (out_diagq) diag(9)%qt = diag(9)%qt + cint*dref*statem%qt
  if (out_diagl) diag(9)%qc = diag(9)%qc + cint*dref*hydromtrm(drop)%q
  if (out_diagr) diag(9)%qr = diag(9)%qr + cint*dref*hydromtrm(rain)%q
  if (out_diagi .and. lmicro > 1) diag(9)%qi = diag(9)%qi + cint*dref*hydromtrm(ice)%q
!
!  Deallocate local arrays
!
  call dealloc ( dref )
  call dealloc ( dTdt )
  call dealloc ( ww )
!
  call dealloc ( statem )
  call dealloc ( pressurem )
  call dealloc ( hydromtrm )
#ifdef AERO_ENABLE
  call dealloc ( aero3dm )
  deallocate ( aero3dm )
#endif
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate microphysics')
! 
  return
!
  end subroutine microphysics
!
end module micro
