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
!	Integrate u,v,w: first calculate wind tendencies including
!			advection, sgs mixing, damping, buoyancy... (but no
!			pressure gradient) then integrate using either forward Euler
!			or RK. If pressure is prognostic, pressure advection is
!			calculated here, and a time-split method is used to 
!			integrate the coupled pressure-momentum system.
!			As an input, all tendency terms already contain SGS 
!			mixing tendencies.
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================

  module momentummod

  USE shared_all
  USE allocation
  USE gradients
  USE advection
  USE windadv
  USE subgrid
  USE funcpack
  USE thermodynamics
  USE boundary_conditions
!
IMPLICIT NONE
!
  private

  public :: momentum

  CONTAINS

!	===================================================
  subroutine momentum    
!	===================================================
!
  integer :: i, j, k, h
  real :: windm, cb, zero(nz)
!
!  IF parallel run
!
  real, dimension(:,:,:), allocatable :: b1, buoy, gpx, gpy, gpz, pp 
  type(atm_winds) :: winda
  type (atm_state) :: state_tmp
  type (atm_winds) :: wind_tmp
  type (hydrometeor), dimension(1:nhydro) :: hydromtr_tmp
!
  if (verbose > 0) call write_debug('Entering momentum')
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!
!
  if ( imp_buoy ) then
    cb = 0.
  else
    cb = 0.5
  endif
!
!  Allocate and initialize
!
  call alloc ( pp )
  call alloc ( b1 )
  call alloc ( buoy )
  call alloc ( gpx )
  call alloc ( gpy )
  call alloc ( gpz )
  call alloc ( winda )
!
!----------------------------------------------------------!
!                   Turbulent diffusion                    !
!----------------------------------------------------------!
!
  if ( with_dif ) call efud ( pressure2, wind2, state2, hydromtr2 ) 
!
!----------------------------------------------------------!
!                        Nudging                           !
!----------------------------------------------------------!
!
  if ( with_nudg ) call nudge_uv
!
!----------------------------------------------------------!
!                        Coriolis                          !
!----------------------------------------------------------!
!     
  if ( with_cor ) call coriolis
!
!----------------------------------------------------------!
!                Explicit pressure gradient                !
!----------------------------------------------------------!
!
#ifdef ANELASTIC
  do k = 1, nz
    pp(:,:,k) = pressure%p(:,:,k) / den0(k)
  enddo
!
  call gradp ( pp, wind%w, gpx, gpy, gpz )
!
  do k = 1, nz
    winds%u(:,:,k) = winds%u(:,:,k) - den0(k)*gpx(:,:,k)
#ifdef MODEL_3D
    winds%v(:,:,k) = winds%v(:,:,k) - den0(k)*gpy(:,:,k)
#endif
    winds%w(:,:,k) = winds%w(:,:,k) - avden0(k)*gpz(:,:,k)
!
    if (ldiag) then
      if (out_diagu) diag(7)%u(:,:,k) = diag(7)%u(:,:,k) - cdiag*den0(k)*gpx(:,:,k)
#ifdef MODEL_3D
      if (out_diagv) diag(7)%v(:,:,k) = diag(7)%v(:,:,k) - cdiag*den0(k)*gpy(:,:,k)
#endif
      if (out_diagw) diag(7)%w(:,:,k) = diag(7)%w(:,:,k) - cdiag*avden0(k)*gpz(:,:,k)
    endif
  enddo
#endif
!
!----------------------------------------------------------!
!                         Buoyancy                         !
!----------------------------------------------------------!
!
#ifdef ANELASTIC
  if ( with_buoy ) then
    call buoyancy ( pressure%dens, state, hydromtr, b1 )
!
    buoy = cb*b1
!
    if (out_diagw) diag(6)%w = diag(6)%w + cdiag*cb*b1 
  endif
!
!  Implicit buoyancy
!
  if ( with_buoy .and. imp_buoy ) then
    call alloc ( state_tmp )
    call alloc ( hydromtr_tmp )
!
    state_tmp = state + cdiag*dt0*states
    do h = 1, nhydro
      hydromtr_tmp(h)%q = hydromtr(h)%q + cdiag*dt0*hydromtrs(h)%q
    enddo
!
    call buoyancy ( pressure%dens, state_tmp, hydromtr_tmp, b1 )
!
    buoy = buoy + (1. - cb)*b1
!
    if (out_diagw) diag(6)%w = diag(6)%w + cdiag*(1. - cb)*b1 
!
    call dealloc ( state_tmp )
    call dealloc ( hydromtr_tmp )
  endif
!
  thermo%buoy = buoy
#endif
!
!----------------------------------------------------------!
!                         Advection                        !
!----------------------------------------------------------!
!
  if (with_adv) call qnadv_uvw ( winda, wind )
!
!----------------------------------------------------------!
!     Assembling tendencies with advection and buoyancy    !
!----------------------------------------------------------!
!
  do k=1,nz
    do j=jt_start,jt_end
      do i=it_start,it_end
	winds%u(i,j,k) = - winda%u(i,j,k) + winds%u(i,j,k)
!
#ifdef MODEL_3D 
	winds%v(i,j,k) = - winda%v(i,j,k) + winds%v(i,j,k)
#endif
!
        winds%w(i,j,k) = - winda%w(i,j,k) + winds%w(i,j,k) + buoy(i,j,k)
      enddo
    enddo
  enddo
!
!----------------------------------------------------------!
!                	Sponge layer                       !
!----------------------------------------------------------!
!
  if (sponge == 1) then
    call damps ( 'u ', wind%u, den0*u00, winds%u )
!
#ifdef MODEL_3D
    call damps ( 'v ', wind%v, den0*v00, winds%v )
#endif
!
    zero = 0.0
    call damps ( 'w ', wind%w, zero, winds%w )
  endif
!
!  Divergence damping
!  
  if ( lddamp ) call damp_div
!
  if (ldiag) then
    if (out_diagu) diag(9)%u = diag(9)%u + cdiag*winds%u
#ifdef MODEL_3D
    if (out_diagv) diag(9)%v = diag(9)%v + cdiag*winds%v
#endif
    if (out_diagw) diag(9)%w = diag(9)%w + cdiag*winds%w
  endif
!
!----------------------------------------------------------!
!               	  Terminate   	                   !
!----------------------------------------------------------!
!
!  Deallocate
!
  call dealloc ( pp )
  call dealloc ( b1 )
  call dealloc ( buoy )
  call dealloc ( gpx )
  call dealloc ( gpy )
  call dealloc ( gpz )
  call dealloc ( winda )
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating momentum')
!
return                  
end                                                    

end module momentummod
