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

  module windadv

  USE gridno
  USE allocation
  USE shared_data	
  USE shared_wind
  USE shared_pressure
  USE shared_diag
  USE funcpack
  USE boundary_conditions
!
  IMPLICIT NONE
!
  private

  real, parameter :: avisc = 1./16.

  public :: qnadv_uvw

  CONTAINS

!  ===========================================
   subroutine qnadv_uvw  ( winda, windl )
!  ===========================================
!
! ---------------------------------
! --- Calculate advection term of u
! ---------------------------------
!
  type(atm_winds), intent(in) :: windl
  type(atm_winds), intent(inout) :: winda
!
! ---------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Enter qnadv_uvw')
!
!  u advection        
!
  call uadv_fd ( winda%u, windl%u, wind_adv%u, wind_adv%v, wind_adv%w )
!
!  v advection    
!
#ifdef MODEL_3D 
  call vadv_fd ( winda%v, windl%v, wind_adv%u, wind_adv%v, wind_adv%w )
#endif
!
!  w advection    
!
  call wadv_fd ( winda%w, windl%w, wind_adv%u, wind_adv%v, wind_adv%w )
!
! ---------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate qnadv_uvw')
!
return
end
!
!  ===========================================
   subroutine uadv_fd  ( ua, u, ut, vt, wt )
!  ===========================================
!
! ---------------------------------
! --- Calculate advection term of u
! ---------------------------------
!
  integer :: i, j, k
!
  real, dimension (ip_start:,jp_start:,:), intent(inout) :: ua
  real, dimension (ip_start:,jp_start:,:), intent(in) :: u, ut, vt, wt
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: uref, uo
!
! ---------------------------------------------------------!
!
  ua = 0.0
  uo = 0.0
  uref = u
!
!  Even time steps
!
  if (mod(ntime,2) == 0) then
!
!  Calculate u advection along x
!
  do k = 1, nz
    do j = jt_start, jt_end
      call advux_fd (uref(:,j,k),ut(:,j,k),ua(:,j,k))
    enddo
  enddo
!
  if (out_diagu) diag(1)%u(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(1)%u(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*ua(ip_start:ip_end,jp_start:jp_end,1:nz)
  uo = ua
!
#ifdef ADV_SPLIT
  uref(ip_start:ip_end,jp_start:jp_end,1:nz) = uref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*ua(ip_start:ip_end,jp_start:jp_end,1:nz)
  ua(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Calculate u advection along y
!
#ifdef MODEL_3D
  do k = 1, nz
    do i = it_start, it_end
      call advuy_fd (uref(i,:,k),vt(i-1,:,k),vt(i,:,k),vt(i-2,:,k),vt(i+1,:,k),ua(i,:,k))
    enddo
  enddo
!
  if (out_diagu) diag(2)%u(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(2)%u(ip_start:ip_end,jp_start:jp_end,1:nz) -  cdiag*(ua(ip_start:ip_end,jp_start:jp_end,1:nz) &
          - uo(ip_start:ip_end,jp_start:jp_end,1:nz))
  uo = ua
!
#ifdef ADV_SPLIT
  uref(ip_start:ip_end,jp_start:jp_end,1:nz) = uref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*ua(ip_start:ip_end,jp_start:jp_end,1:nz)
  ua(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
#endif
!
!  Calculate u advection along z
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      call advuz_fd (uref(i,j,:),wt(i-1,j,:),wt(i,j,:),wt(i-2,j,:),wt(i+1,j,:),ua(i,j,:))
    enddo
  enddo
!
  if (out_diagu) diag(3)%u(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(3)%u(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(ua(ip_start:ip_end,jp_start:jp_end,1:nz) &
          - uo(ip_start:ip_end,jp_start:jp_end,1:nz))
!
#ifdef ADV_SPLIT
  uref(ip_start:ip_end,jp_start:jp_end,1:nz) = uref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*ua(ip_start:ip_end,jp_start:jp_end,1:nz)
  ua(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Odd time-steps
!
  else
!
!  Calculate u advection along z
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      call advuz_fd (uref(i,j,:),wt(i-1,j,:),wt(i,j,:),wt(i-2,j,:),wt(i+1,j,:),ua(i,j,:))
    enddo
  enddo
!
  if (out_diagu) diag(3)%u(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(3)%u(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*ua(ip_start:ip_end,jp_start:jp_end,1:nz)
  uo = ua
!
#ifdef ADV_SPLIT
  uref(ip_start:ip_end,jp_start:jp_end,1:nz) = uref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*ua(ip_start:ip_end,jp_start:jp_end,1:nz)
  ua(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Calculate u advection along y
!
#ifdef MODEL_3D
  do k = 1, nz
    do i = it_start, it_end
      call advuy_fd (uref(i,:,k),vt(i-1,:,k),vt(i,:,k),vt(i-2,:,k),vt(i+1,:,k),ua(i,:,k))
    enddo
  enddo
!
  if (out_diagu) diag(2)%u(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(2)%u(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(ua(ip_start:ip_end,jp_start:jp_end,1:nz) &
          - uo(ip_start:ip_end,jp_start:jp_end,1:nz))
  uo = ua
!
#ifdef ADV_SPLIT
  uref(ip_start:ip_end,jp_start:jp_end,1:nz) = uref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*ua(ip_start:ip_end,jp_start:jp_end,1:nz)
  ua(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
#endif
!
!  Calculate u advection along x
!
  do k = 1, nz
    do j = jt_start, jt_end
      call advux_fd (uref(:,j,k),ut(:,j,k),ua(:,j,k))
    enddo
  enddo
!
  if (out_diagu) diag(1)%u(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(1)%u(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(ua(ip_start:ip_end,jp_start:jp_end,1:nz) &
          - uo(ip_start:ip_end,jp_start:jp_end,1:nz))
!
#ifdef ADV_SPLIT
  uref(ip_start:ip_end,jp_start:jp_end,1:nz) = uref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*ua(ip_start:ip_end,jp_start:jp_end,1:nz)
  ua(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
  endif
!
#ifdef ADV_SPLIT
  ua(ip_start:ip_end,jp_start:jp_end,1:nz) = -(uref(ip_start:ip_end,jp_start:jp_end,1:nz) - u(ip_start:ip_end,jp_start:jp_end,1:nz)) / dt0
#endif
!
  if (verbose > 2) call write_debug('Terminate uadv_fd')
!
return
end
!
!   ===========================================
    subroutine vadv_fd ( va, v,  ut, vt, wt )
!   ===========================================

! ---------------------------------
! --- Calculate advection term of v
! ---------------------------------
!
  integer :: i, j, k, j1
!
  real, dimension (ip_start:,jp_start:,:), intent(inout) :: va
  real, dimension (ip_start:,jp_start:,:), intent(in) :: v, ut, vt, wt
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: vref, vo
!
! -------------------------------------------------------
!
  va = 0.0
  vo = 0.0
  vref = v
!
!  Even time steps
!
  if (mod(ntime,2) == 0) then
!
!  Calculate v advection along x
!
  do k = 1, nz
    do j = jt_start, jt_end
      call advvx_fd (vref(:,j,k), ut(:,j-1,k), ut(:,j,k), ut(:,j-2,k), ut(:,j+1,k), va(:,j,k))
    enddo
  enddo
!
  if (out_diagv) diag(1)%v(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(1)%v(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*va(ip_start:ip_end,jp_start:jp_end,1:nz)
  vo = va
!
#ifdef ADV_SPLIT
  vref(ip_start:ip_end,jp_start:jp_end,1:nz) = vref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*va(ip_start:ip_end,jp_start:jp_end,1:nz)
  va(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Calculate v advection along y
!
  do k = 1, nz
    do i = it_start, it_end
      call advvy_fd (vref(i,:,k), vt(i,:,k), va(i,:,k))
    enddo
  enddo
!
  if (out_diagv) diag(2)%v(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(2)%v(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(va(ip_start:ip_end,jp_start:jp_end,1:nz) &
          - vo(ip_start:ip_end,jp_start:jp_end,1:nz))
  vo = va
!
#ifdef ADV_SPLIT
  vref(ip_start:ip_end,jp_start:jp_end,1:nz) = vref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*va(ip_start:ip_end,jp_start:jp_end,1:nz)
  va(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Calculate v advection along z
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      call advvz_fd (vref(i,j,:), wt(i,j-1,:), wt(i,j,:), wt(i,j-2,:), wt(i,j+1,:), va(i,j,:))
    enddo
  enddo
!
  if (out_diagv) diag(3)%v(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(3)%v(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(va(ip_start:ip_end,jp_start:jp_end,1:nz) & 
          - vo(ip_start:ip_end,jp_start:jp_end,1:nz))
!
#ifdef ADV_SPLIT
  vref(ip_start:ip_end,jp_start:jp_end,1:nz) = vref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*va(ip_start:ip_end,jp_start:jp_end,1:nz)
  va(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Odd time-steps
!
  else
!
!  Calculate v advection along z
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      call advvz_fd (vref(i,j,:), wt(i,j-1,:), wt(i,j,:), wt(i,j-2,:), wt(i,j+1,:), va(i,j,:))
    enddo
  enddo
!
  if (out_diagv) diag(3)%v(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(3)%v(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*va(ip_start:ip_end,jp_start:jp_end,1:nz)
  vo = va
!
#ifdef ADV_SPLIT
  vref(ip_start:ip_end,jp_start:jp_end,1:nz) = vref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*va(ip_start:ip_end,jp_start:jp_end,1:nz)
  va(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Calculate v advection along y
!
  do k = 1, nz
    do i = it_start, it_end
      call advvy_fd (vref(i,:,k), vt(i,:,k), va(i,:,k))
    enddo
  enddo
!
  if (out_diagv) diag(2)%v(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(2)%v(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(va(ip_start:ip_end,jp_start:jp_end,1:nz) & 
          - vo(ip_start:ip_end,jp_start:jp_end,1:nz))
  vo = va
!
#ifdef ADV_SPLIT
  vref(ip_start:ip_end,jp_start:jp_end,1:nz) = vref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*va(ip_start:ip_end,jp_start:jp_end,1:nz)
  va(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Calculate v advection along x
!
  do k = 1, nz
    do j = jt_start, jt_end
      call advvx_fd (vref(:,j,k), ut(:,j-1,k), ut(:,j,k), ut(:,j-2,k), ut(:,j+1,k), va(:,j,k))
    enddo
  enddo
!
  if (out_diagv) diag(1)%v(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(1)%v(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(va(ip_start:ip_end,jp_start:jp_end,1:nz) & 
          - vo(ip_start:ip_end,jp_start:jp_end,1:nz))
!
#ifdef ADV_SPLIT
  vref(ip_start:ip_end,jp_start:jp_end,1:nz) = vref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*va(ip_start:ip_end,jp_start:jp_end,1:nz)
  va(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
  endif
!
#ifdef ADV_SPLIT
  va(ip_start:ip_end,jp_start:jp_end,1:nz) = -(vref(ip_start:ip_end,jp_start:jp_end,1:nz) - v(ip_start:ip_end,jp_start:jp_end,1:nz)) / dt0
#endif
!
  if (verbose > 2) call write_debug('Terminate vadv_fd')
!
return
end
!
!  ===========================================
   subroutine wadv_fd ( wa, w, ut, vt, wt )
!  ===========================================

! ---------------------------------
! --- Calculate advection term of w
! ---------------------------------
!
IMPLICIT NONE
!
  integer :: i, j, k
!
  real, dimension (ip_start:,jp_start:,:), intent(inout) :: wa
  real, dimension (ip_start:,jp_start:,:), intent(in) :: w, ut, vt, wt
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: wref, wo
!
! --------------------------------------------------------
!
  wa = 0.0
  wo = 0.0
  wref = w
!
!  Even time steps
!
  if (mod(ntime,2) == 0) then
!
!  Calculate w advection along x
!
  do j = jt_start, jt_end
    call advwx_fd (wref(:,j,:), ut(:,j,:), wa(:,j,:))
  enddo
!
  if (out_diagw) diag(1)%w(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(1)%w(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*wa(ip_start:ip_end,jp_start:jp_end,1:nz)
  wo = wa
!
#ifdef ADV_SPLIT
  wref(ip_start:ip_end,jp_start:jp_end,1:nz) = wref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*wa(ip_start:ip_end,jp_start:jp_end,1:nz)
  wa(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Calculate w advection along y
!
#ifdef MODEL_3D
  do i = it_start, it_end
    call advwy_fd (wref(i,:,:), vt(i,:,:), wa(i,:,:))
  enddo
!
  if (out_diagw) diag(2)%w(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(2)%w(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(wa(ip_start:ip_end,jp_start:jp_end,1:nz) &
          - wo(ip_start:ip_end,jp_start:jp_end,1:nz))
  wo = wa
!
#ifdef ADV_SPLIT
  wref(ip_start:ip_end,jp_start:jp_end,1:nz) = wref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*wa(ip_start:ip_end,jp_start:jp_end,1:nz)
  wa(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
#endif
!
!  Calculate w advection along z
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      call advwz_fd (wref(i,j,:), wt(i,j,:), wa(i,j,:))
    enddo
  enddo
!
  if (out_diagw) diag(3)%w(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(3)%w(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(wa(ip_start:ip_end,jp_start:jp_end,1:nz) &
          - wo(ip_start:ip_end,jp_start:jp_end,1:nz))
!
#ifdef ADV_SPLIT
  wref(ip_start:ip_end,jp_start:jp_end,1:nz) = wref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*wa(ip_start:ip_end,jp_start:jp_end,1:nz)
  wa(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Odd time-steps
!
  else
!
!  Calculate w advection along z
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      call advwz_fd (wref(i,j,:), wt(i,j,:), wa(i,j,:))
    enddo
  enddo
!
  if (out_diagw) diag(3)%w(ip_start:ip_end,jp_start:jp_end,1:nz) =diag(3)%w(ip_start:ip_end,jp_start:jp_end,1:nz)  - cdiag*wa(ip_start:ip_end,jp_start:jp_end,1:nz)
  wo = wa
!
#ifdef ADV_SPLIT
  wref(ip_start:ip_end,jp_start:jp_end,1:nz) = wref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*wa(ip_start:ip_end,jp_start:jp_end,1:nz)
  wa(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
!  Calculate w advection along y
!
#ifdef MODEL_3D
  do i = it_start, it_end
    call advwy_fd (wref(i,:,:), vt(i,:,:), wa(i,:,:))
  enddo
!
  if (out_diagw) diag(2)%w(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(2)%w(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(wa(ip_start:ip_end,jp_start:jp_end,1:nz) &
          - wo(ip_start:ip_end,jp_start:jp_end,1:nz))
  wo = wa
!
#ifdef ADV_SPLIT
  wref(ip_start:ip_end,jp_start:jp_end,1:nz) = wref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*wa(ip_start:ip_end,jp_start:jp_end,1:nz)
  wa(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
#endif
!
!  Calculate w advection along x
!
  do j = jt_start, jt_end
    call advwx_fd (wref(:,j,:), ut(:,j,:), wa(:,j,:))
  enddo
!
  if (out_diagw) diag(1)%w(ip_start:ip_end,jp_start:jp_end,1:nz) = diag(1)%w(ip_start:ip_end,jp_start:jp_end,1:nz) - cdiag*(wa(ip_start:ip_end,jp_start:jp_end,1:nz) &
         - wo(ip_start:ip_end,jp_start:jp_end,1:nz))
!
#ifdef ADV_SPLIT
  wref(ip_start:ip_end,jp_start:jp_end,1:nz) = wref(ip_start:ip_end,jp_start:jp_end,1:nz) - dt0*wa(ip_start:ip_end,jp_start:jp_end,1:nz)
  wa(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
#endif
!
  endif
!
#ifdef ADV_SPLIT
  wa(ip_start:ip_end,jp_start:jp_end,1:nz) = -(wref(ip_start:ip_end,jp_start:jp_end,1:nz) - w(ip_start:ip_end,jp_start:jp_end,1:nz)) / dt0
#endif
!
  if (verbose > 2) call write_debug('Terminate wadv_fd')
!
return
end
!
  subroutine advux_fd (u,ut,ua)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (ip_start:ip_end) :: ua, u, ut, Flux, Adv
!
  integer  :: i
  real     :: fr, fr2, fr4, fr5, dx1, dx2, uu, uu2
  real     :: c1, c2, mu
!                                                          !
!----------------------------------------------------------!
!
  dx1 = 1./dx
  dx2 = 1./(2.*dx)
  if (mom_ord == 4) then
    c1 = 4./3.; c2 = 1./3.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
  Flux = 0.0
  Adv = 0.0
!
!  Calculate fluxes: 4th order
!
  do i = it_start-2, it_end
    uu  = (ut(i+1) + ut(i))*0.5 
    uu2 = (ut(i+2) + ut(i))*0.5
    fr  = (u(i+1) + u(i))*0.5
    fr2 = (u(i+2) + u(i))*0.5
    Flux(i) = uu*fr
    Adv(i)  = uu2*fr2 
  enddo
!  
!  Calculate advection
!
  do i = it_start, it_end
    ua(i) = ua(i) + c1*(Flux(i) - Flux(i-1))*dx1 - c2*(Adv(i) - Adv(i-2))*dx2 
  enddo     
!
!  Artificial viscosity
!
  if (lavisc) then
    do i = it_start, it_end
      mu = abs(ut(i))*dx1
      ua(i) = ua(i) - avisc*mu*(-u(i+2) + 4.*u(i+1) - 6.*u(i) + 4.*u(i-1) - u(i-2))
    enddo
  endif
!
return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advux_fd
!
  subroutine advuy_fd (u,vt1,vt2,vt3,vt4,ua)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (jp_start:jp_end) :: u, ua, vt1, vt2, vt3, vt4, Flux
!
  integer  :: j
  real     :: fr, fr2, fr4, fr5, vv
  real     :: dy1, dy3, c1, c2, mu
!                                                          !
!----------------------------------------------------------!
!
  dy1 = 1./dy
  dy3 = 1./(3.*dy)
  if (mom_ord == 4) then
    c1 = 9./8.; c2 = 1./8.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
  Flux = 0.0
!
!  Calculate fluxes
!
  do j = jt_start-2, jt_end+1
    vv  = c1*(vt1(j+1)+vt2(j+1))*0.5 - c2*(vt3(j+1)+vt4(j+1))*0.5
    fr  = (u(j+1) + u(j))*0.5
    Flux(j) = vv*fr
  enddo
!
!  Boundary conditions
!
#ifdef CHANNEL
    if (jp_start == 1) then
      Flux(jj1) = 0.
      Flux(jj1+1) = 0.5*(vt1(jj1+2)+vt2(jj1+2))*(u(jj1+2) - u(jj1+1))
    endif
! 
    if (jp_end == ny) then
      Flux(jj2) = 0.
      Flux(jj2-1) = 0.5*(vt1(jj2)+vt2(jj2))*(u(jj2) - u(jj2-1))
    endif
#endif
!
!  Calculate advection
!
    do j = jt_start, jt_end
      ua(j) = ua(j) + c1*(Flux(j) - Flux(j-1))*dy1 - c2*(Flux(j+1) - Flux(j-2))*dy3
    enddo
!
!  Artificial viscosity
!
  if (lavisc) then
    do j = jt_start, jt_end
      mu = 0.25*abs(vt1(j)+vt2(j)+vt1(j+1)+vt2(j+1))*dy1
      ua(j) = ua(j) - avisc*mu*(-u(j+2) + 4.*u(j+1) - 6.*u(j) + 4.*u(j-1) - u(j-2))
    enddo
  endif
!
return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advuy_fd
!
  subroutine advuz_fd (u,wt1,wt2,wt3,wt4,ua)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (1:nz) :: u, ua, wt1, wt2, wt3, wt4, Flux
!
  integer  :: k
  real     :: fr, fr2, fr4, fr5, fk2, fk1, ww, dt1
  real     :: c1, c2, mu
!                                                          !
!----------------------------------------------------------!
!
  Flux = 0.0
  if (mom_ord == 4) then
    c1 = 9./8.; c2 = 1./8.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
!  Fluxes
!
  do k=1,nz-1
    ww  = c1*(wt1(k+1)+wt2(k+1))*0.5 - c2*(wt3(k+1)+wt4(k+1))*0.5
    fr  = (dzw(k+1)*u(k+1) + dzw(k)*u(k))/(dzw(k+1) + dzw(k))
    Flux(k) = ww*fr
  enddo
!
  if (bcl(6) == 'nnf') Flux(nz) = 0.
  if (bcl(6) == 'ope') Flux(nz) = Flux(nz-1)
!
!  Main domain core	
!
  do k = 3,nz-1
    ua(k) = ua(k) + c1*(Flux(k) - Flux(k-1))*dz1w(k) - c2*(Flux(k+1) - Flux(k-2))/(dzw(k-1)+dzw(k)+dzw(k+1))
  enddo 
!
!  Boundary conditions
!
  if (bcl(5) == 'nnf') ua(1) = ua(1) + (Flux(1) - 0.)*dz1w(1)
  if (bcl(5) == 'nnf') ua(2) = ua(2) + c1*(Flux(2) - Flux(1))*dz1w(1) - c2*(Flux(3) - 0.)/(dzw(1)+dzw(2)+dzw(3))
  if (bcl(5) == 'ope') ua(1) = ua(1) + (Flux(1) - Flux(1))*dz1w(1)
  if (bcl(5) == 'ope') ua(2) = ua(2) + c1*(Flux(2) - Flux(1))*dz1w(1) - c2*(Flux(3) - Flux(1))/(dzw(1)+dzw(2)+dzw(3))
!
  if (bcl(6) == 'nnf') ua(nz) = ua(nz) + (0. - Flux(nz-1))*dz1w(nz)
  if (bcl(6) == 'ope') ua(nz) = ua(nz) + (Flux(nz-1) - Flux(nz-1))*dz1w(nz)
!
!  Artificial viscosity
!
  if (lavisc) then
    do k = 3, nz-2
      mu = 0.25*abs(wt1(k)+wt2(k)+wt1(k+1)+wt2(k+1))*dz1w(k)
      ua(k) = ua(k) - avisc*mu*(-u(k+2) + 4.*u(k+1) - 6.*u(k) + 4.*u(k-1) - u(k-2))
    enddo
  endif
!
return
!
!----------------------------------------------------------!
!
  end subroutine advuz_fd
!
  subroutine advvx_fd (v,ut1,ut2,ut3,ut4,va)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (ip_start:ip_end) :: v, va, ut1, ut2, ut3, ut4, Flux
!
  integer  :: i
  real     :: fr, fr2, fr4, fr5, uu
  real     :: dx1, dx3, c1, c2, mu
!
!----------------------------------------------------------!
!
  dx1 = 1./dx
  dx3 = 1./(3.*dx)
  if (mom_ord == 4) then
    c1 = 9./8.; c2 = 1./8.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
  Flux = 0.0
!
!  Calculate fluxes
!
  do i = it_start-2, it_end+1
    uu = c1*(ut1(i+1)+ut2(i+1))*0.5 - c2*(ut3(i+1)+ut4(i+1))*0.5
    fr = (v(i+1) + v(i))*0.5
    Flux(i) = uu*fr
  enddo
!  
!  Calculate advection
!
  do i = it_start, it_end
    va(i) = va(i) + c1*(Flux(i) - Flux(i-1))*dx1 - c2*(Flux(i+1) - Flux(i-2))*dx3
  enddo
!
!  Artificial viscosity
!
  if (lavisc) then
    do i = it_start, it_end
      mu = 0.25*abs(ut1(i)+ut2(i)+ut1(i+1)+ut2(i+1))*dx1
      va(i) = va(i) - avisc*mu*(-v(i+2) + 4.*v(i+1) - 6.*v(i) + 4.*v(i-1) - v(i-2))
    enddo
  endif
!
return
! 
!----------------------------------------------------------!
!
  end subroutine advvx_fd 
!
  subroutine advvy_fd (v,vt,va)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (jp_start:jp_end) :: v, va, vt, Flux, Adv
!
  integer  :: j
  real     :: fr, fr2, fr3, fr4, vv, vv2
  real     :: dy1, dy2, c1, c2, mu
!                                                          !
!----------------------------------------------------------!
!
  dy1 = 1./dy
  dy2 = 1./(2.*dy)
  if (mom_ord == 4) then
    c1 = 4./3.; c2 = 1./3.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
  Flux = 0.0
  Adv = 0.0
!
!  Calculate fluxes
!
  do j=jt_start-2,jt_end
    vv = (vt(j+1) + vt(j))*0.5
    vv2 = (vt(j+2) + vt(j))*0.5
    fr = (v(j+1) + v(j))*0.5
    fr2 = (v(j+2) + v(j))*0.5
    Flux(j) = vv*fr
    Adv(j) = vv2*fr2
  enddo  
!
!  Boundary conditions
!
#ifdef CHANNEL
    if (jp_start == 1) then
      Flux(jj1) = 0.25*(vt(jj1)+0.)*(vt(jj1) - 0.)
    endif
! 
    if (jp_end == ny) then
      Flux(jj2) = 0.25*(0.+vt(jj2))*(0. - vt(jj2))
    endif
#endif
!
!  Calculate advection
!
  do j = jt_start, jt_end
    va(j) = va(j) + c1*(Flux(j) - Flux(j-1))*dy1 - c2*(Adv(j) - Adv(j-2))*dy2
  enddo
!
!  Artificial viscosity
!
  if (lavisc) then
    do j = jt_start, jt_end
      mu = abs(vt(j))*dy1
      va(j) = va(j) - avisc*mu*(-v(j+2) + 4.*v(j+1) - 6.*v(j) + 4.*v(j-1) - v(j-2))
    enddo
  endif
!
return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advvy_fd
!
  subroutine advvz_fd (v,wt1,wt2,wt3,wt4,va)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (1:nz) :: v, va, wt1, wt2, wt3, wt4, Flux
!
  integer  :: k
  real     :: fr, fr2, fr4, fr5, fk1, fk2, ww
  real     :: c1, c2, mu
!                                                          !
!----------------------------------------------------------!
!
  Flux = 0.0
  if (mom_ord == 4) then
    c1 = 9./8.; c2 = 1./8.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
!  Fluxes
!
  do k=1,nz-1
    ww  = c1*(wt1(k+1)+wt2(k+1))*0.5 - c2*(wt3(k+1)+wt4(k+1))*0.5
    fr  = (dzw(k+1)*v(k+1) + dzw(k)*v(k))/(dzw(k+1) + dzw(k))
    Flux(k) = ww*fr
  enddo
!
  if (bcl(6) == 'nnf') Flux(nz) = 0.
  if (bcl(6) == 'ope') Flux(nz) = Flux(nz-1)
!
!  Main domain core
!
  do k = 3, nz-1
    va(k) = va(k) + c1*(Flux(k) - Flux(k-1))*dz1w(k) - c2*(Flux(k+1) - Flux(k-2))/(dzw(k-1)+dzw(k)+dzw(k+1))
  enddo  
!
!  Boundary conditions
!
  if (bcl(5) == 'nnf') va(1) = va(1) + (Flux(1) - 0.)*dz1w(1)
  if (bcl(5) == 'nnf') va(2) = va(2) + c1*(Flux(2) - Flux(1))*dz1w(2) - c2*(Flux(3) - 0.)/(dzw(1)+dzw(2)+dzw(3))
  if (bcl(5) == 'ope') va(1) = va(1) + (Flux(1) - Flux(1))*dz1w(1)
  if (bcl(5) == 'ope') va(2) = va(2) + c1*(Flux(2) - Flux(1))*dz1w(2) - c2*(Flux(3) - Flux(1))/(dzw(1)+dzw(2)+dzw(3))
!
  if (bcl(6) == 'nnf') va(nz) = va(nz) + (0. - Flux(nz-1))*dz1w(nz)
  if (bcl(6) == 'ope') va(nz) = va(nz) + (Flux(nz-1) - Flux(nz-1))*dz1w(nz)
!
!  Artificial viscosity
!
  if (lavisc) then
    do k = 3, nz-2
      mu = 0.25*abs(wt1(k)+wt2(k)+wt1(k+1)+wt2(k+1))*dz1w(k)
      va(k) = va(k) - avisc*mu*(-v(k+2) + 4.*v(k+1) - 6.*v(k) + 4.*v(k-1) - v(k-2))
    enddo
  endif
!
  return
! 
!----------------------------------------------------------!
!
  end subroutine advvz_fd 
!
  subroutine advwx_fd (w,ut,wa)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (ip_start:ip_end,1:nz) :: w, wa, ut, Flux
!
  integer  :: i, k 
  real     :: fr, fr2, fr4, fr5, uu
  real     :: dx1, dx3, c1, c2, mu
!                                                          !
!----------------------------------------------------------!
!
  dx1 = 1./dx
  dx3 = 1./(3.*dx)
  if (mom_ord == 4) then
    c1 = 9./8.; c2 = 1./8.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
  Flux = 0.0
!
!  Calculate fluxes
!
  do i=it_start-2,it_end+1
    uu  = (ut(i+1,2) + ut(i+1,1))*0.5
    fr  = (w(i+1,2) + w(i,2))*0.5
    Flux(i,2) = uu*fr
!
    uu  = (ut(i+1,nz) + ut(i+1,nz-1))*0.5
    fr  = (w(i+1,nz) + w(i,nz))*0.5
    Flux(i,nz) = uu*fr
  enddo
!
  do k = 3, nz-1
    do i=it_start-2,it_end+1
      uu  = c1*(ut(i+1,k)+ut(i+1,k-1))*0.5 - c2*(ut(i+1,k+1)+ut(i+1,k-2))*0.5
      fr  = (w(i+1,k) + w(i,k))*0.5
      Flux(i,k) = uu*fr
    enddo
  enddo
!  
!  Calculate advection
!
  do k = 2, nz
    do i = it_start, it_end
      wa(i,k) = wa(i,k) + c1*(Flux(i,k) - Flux(i-1,k))*dx1 - c2*(Flux(i+1,k) - Flux(i-2,k))*dx3
    enddo
  enddo
!
!  Artificial viscosity
!
  if (lavisc) then
    do k = 2, nz
      do i = it_start, it_end
        mu = 0.25*abs(ut(i,k)+ut(i,k-1)+ut(i+1,k)+ut(i+1,k-1))*dx1
        wa(i,k) = wa(i,k) - avisc*mu*(-w(i+2,k) + 4.*w(i+1,k) - 6.*w(i,k) + 4.*w(i-1,k) - w(i-2,k))
      enddo
    enddo
  endif
!
return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advwx_fd
!
  subroutine advwy_fd (w,vt,wa)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (jp_start:jp_end,1:nz) :: w, wa, vt, Flux, Adv
!
  integer  :: k, j
  real     :: fr, fr2, fr4, fr5, vv
  real     :: dy1, dy3, c1, c2, mu
!                                                          !
!----------------------------------------------------------!
!
  dy1 = 1./dy
  dy3 = 1./(3.*dy)
  if (mom_ord == 4) then
    c1 = 9./8.; c2 = 1./8.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
  Flux = 0.0
!
!  Calculate fluxes
!
  do j=jt_start-2,jt_end+1
    vv  = (vt(j+1,2) + vt(j+1,1))*0.5
    fr  = (w(j+1,2) + w(j,2))*0.5
    Flux(j,2) = vv*fr
!
    vv  = (vt(j+1,nz) + vt(j+1,nz-1))*0.5
    fr  = (w(j+1,nz) + w(j,nz))*0.5
    Flux(j,nz) = vv*fr
  enddo
!
  do k = 3, nz-1
    do j=jt_start-2, jt_end+1
      vv  = c1*(vt(j+1,k)+vt(j+1,k-1))*0.5 - c2*(vt(j+1,k+1)+vt(j+1,k-2))*0.5
      fr  = (w(j+1,k) + w(j,k))*0.5
      Flux(j,k) = vv*fr
    enddo
  enddo
!
!  Boundary conditions
!
#ifdef CHANNEL
    if (jp_start == 1) then
      Flux(jj1) = 0.
      Flux(jj1+1) = 0.5*(vt1(jj1+2)+vt2(jj1+2))*(w(jj1+2) - w(jj1+1))
    endif
! 
    if (jp_end == ny) then
      Flux(jj2) = 0.
      Flux(jj2-1) = 0.5*(vt1(jj2)+vt2(jj2))*(w(jj2) - w(jj2-1))
    endif
#endif
!
!  Calculate advection
!
    do k = 2, nz
      do j = jt_start, jt_end
        wa(j,k) = wa(j,k) + c1*(Flux(j,k) - Flux(j-1,k))*dy1 - c2*(Flux(j+1,k) - Flux(j-2,k))*dy3
      enddo
    enddo
!
!  Artificial viscosity
!
  if (lavisc) then
    do k = 2, nz
      do j = jt_start, jt_end
        mu = 0.25*abs(vt(j,k)+vt(j,k-1)+vt(j+1,k)+vt(j+1,k-1))*dy1
        wa(j,k) = wa(j,k) - avisc*mu*(-w(j+2,k) + 4.*w(j+1,k) - 6.*w(j,k) + 4.*w(j-1,k) - w(j-2,k))
      enddo
    enddo
  endif
!
return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advwy_fd 
!
  subroutine advwz_fd (w,wt,wa)
!
!----------------------------------------------------------!
!                                                          !
!  Calculates gradients with centered 4th order finite     !
!  differences : array                                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (1:nz) :: w, wa, wt, Flux, Adv
!
  integer  :: k
  real     :: fr, fr2, ww, ww2
  real     :: c1, c2, mu
!                                                          !
!----------------------------------------------------------!
!
  Flux = 0.0
  Adv = 0.0
  if (mom_ord == 4) then
    c1 = 4./3.; c2 = 1./3.
  else if (mom_ord == 2) then
    c1 = 1.; c2 = 0.
  endif
!
!  Rest of the domain
!
  do k = 1, nz-2
    ww  = (wt(k+1) + wt(k))*0.5
    ww2 = (dzw(k+1)*wt(k+2) + dzw(k)*wt(k))/(dzw(k+1) + dzw(k))
    fr  = (w(k+1) + w(k))*0.5
    fr2 = (dzw(k+1)*w(k+2) + dzw(k)*w(k))/(dzw(k+1) + dzw(k))
    Flux(k) = ww*fr
    Adv(k) = ww2*fr2
  enddo
!
  Flux(nz-1) = 0.25*(wt(nz) + wt(nz-1))*(w(nz) + w(nz-1)) 
  if (bcl(6) == 'nnf') Adv(nz-1) = (0. + dzw(nz-1)*wt(nz-1))*(0. + dzw(nz-1)*w(nz-1))/(dzw(nz) + dzw(nz-1))**2.
  if (bcl(6) == 'ope') Adv(nz-1) = (dzw(nz)*wt(nz) + dzw(nz-1)*wt(nz-1))*(dzw(nz)*w(nz) + dzw(nz-1)*w(nz-1))/(dzw(nz) + dzw(nz-1))**2.
!
!  Assemble fluxes
!
  do k = 3,nz-1
    wa(k) = wa(k) + c1*(Flux(k) - Flux(k-1))*dz1(k-1) - c2*(Adv(k) - Adv(k-2))/(dzw(k)+dzw(k-1))
  enddo
!
!  Boundary conditions
!
  !if (bcl(5) == 'nnf') wa(1) = wa(1) + c1*(Flux(1) + Flux(1))*dz1w(1)
  if (bcl(5) == 'nnf') wa(2) = wa(2) + c1*(Flux(2) - Flux(1))*dz1(1) - c2*(Adv(2) - 0.)/(dzw(2)+dzw(1))
  if (bcl(5) == 'ope') wa(2) = wa(2) + c1*(Flux(2) - Flux(1))*dz1(1) - c2*(Adv(2) - Flux(1))/(dzw(2)+dzw(1))
  if (bcl(6) == 'nnf') wa(nz) = wa(nz) + (0. - Flux(nz-1))*dz1(nz-1) 
  if (bcl(6) == 'ope') wa(nz) = wa(nz) + (Flux(nz-1) - Flux(nz-1))*dz1(nz-1) 
!
!  Artificial viscosity
!
  if (lavisc) then
    do k = 3, nz-2
      mu = abs(wt(k))*dz1(k-1)
      wa(k) = wa(k) - avisc*mu*(-w(k+2) + 4.*w(k+1) - 6.*w(k) + 4.*w(k-1) - w(k-2))
    enddo
  endif
!
return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advwz_fd
! 
end module windadv
