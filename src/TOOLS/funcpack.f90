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
!  FUNCPACK:
!	Package of functions
!
!  Purpose:
!	Put functions and utility subroutines together.
!
!  Author:
!	Julien Savre, Ludwig-Maximilian-Universitat, Munich
!
! ================================================================

  module funcpack

  USE shared_all
  USE nesting
  USE averages
  USE gradients
  USE thermodynamics
  USE boundary_conditions
!
  IMPLICIT NONE
!
  private

  interface get_momentum
    MODULE PROCEDURE get_momentum_1, get_momentum_3
  end interface get_momentum

  interface get_velocities
    MODULE PROCEDURE get_velocities_1, get_velocities_3
  end interface get_velocities

  public :: get_momentum, get_velocities, get_conserved, get_scalars, get_advective
  public :: damps, damp_div
  public :: cal_p1

  CONTAINS
!
! ===============================================
  subroutine get_advective ( pressurel, windl )
! ===============================================
!
! --- Calculate advective velocity
!
!  Input/Output
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_winds), intent(inout) :: windl
!  
  wind_adv%u = windl%u 
#ifdef MODEL_3D
  wind_adv%v = windl%v
#endif
  wind_adv%w = windl%w
!
!  Get velocities
!
  call get_velocities ( pressurel, wind_adv ) 
!
  wind_adv%u = wind_adv%u - u0shift
#ifdef MODEL_3D
  wind_adv%v = wind_adv%v - v0shift
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine get_advective
!
! ===============================================
  subroutine get_momentum_1 ( pressurel, windl )
! ===============================================
!
! --- Calculate momentum from density & velocity
!
!  Input/Output
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_winds), intent(inout) :: windl
!  
  integer :: k
  integer :: i, j
  real, dimension(3,ip_start:ip_end,jp_start:jp_end,1:nz) :: avden
!
  if (verbose > 2) call write_debug('Starting get_momentum')
!
!----------------------------------------------------------!
!                     momentum = rho*U                     !
!----------------------------------------------------------!
!
!  Compressible case
!
#ifndef ANELASTIC
!
  avden = 1.
!
!  X
!
  do i = ip_start+1, ip_end
    avden(1,i,jp_start:jp_end,:) = 0.5*(pressurel%dens(i,jp_start:jp_end,:)+pressurel%dens(i-1,jp_start:jp_end,:))
  enddo
!
!  Y
!
#ifdef MODEL_3D
  do j = jp_start+1, jp_end
    avden(2,ip_start:ip_end,j,:) = 0.5*(pressurel%dens(ip_start:ip_end,j,:)+pressurel%dens(ip_start:ip_end,j-1,:))
  enddo
#endif
!
!  Z
!
  do k = 2, nz
    avden(3,ip_start:ip_end,jp_start:jp_end,k) = (fdz0(k)*pressurel%dens(ip_start:ip_end,jp_start:jp_end,k)+fdz0(k-1)*pressurel%dens(ip_start:ip_end,jp_start:jp_end,k-1)) / (fdz0(k)+fdz0(k-1))
  enddo
  avden(3,ip_start:ip_end,jp_start:jp_end,1) = pressurel%dens(ip_start:ip_end,jp_start:jp_end,1)
!
!  Get momentum
!
  windl%u = windl%u*avden(1,:,:,:)
!
#ifdef MODEL_3D
  windl%v = windl%v*avden(2,:,:,:)
#endif
!
  windl%w = windl%w*avden(3,:,:,:)
!
!  BC
!
  call windbc ( windl )
!
#else
!
!  Anelastic case
!
  do k = 1, nz
    windl%u(:,:,k) = windl%u(:,:,k) * den0(k)
#ifdef MODEL_3D
    windl%v(:,:,k) = windl%v(:,:,k) * den0(k)
#endif
    windl%w(:,:,k) = windl%w(:,:,k) * avden0(k)
  enddo
!
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating get_momentum')
!
  end subroutine get_momentum_1
!
! ===============================================
  subroutine get_momentum_3 ( pressurel, u, v, w )
! ===============================================
!
! --- Calculate momentum from density & velocity
!
!  Input/Output
!
  type (atm_pressure), intent(in) :: pressurel
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: u, v, w
!  
  integer :: k
  integer :: i, j
  real, dimension(3,ip_start:ip_end,jp_start:jp_end,1:nz) :: avden
!
  if (verbose > 2) call write_debug('Starting get_momentum')
!
!----------------------------------------------------------!
!                     momentum = rho*U                     !
!----------------------------------------------------------!
!
!  Compressible case
!
#ifndef ANELASTIC
!
  avden = 1.
!
!  X
!
  do i = ip_start+1, ip_end
    avden(1,i,jp_start:jp_end,:) = 0.5*(pressurel%dens(i,jp_start:jp_end,:)+pressurel%dens(i-1,jp_start:jp_end,:))
  enddo
!
!  Y
!
#ifdef MODEL_3D
  do j = jp_start+1, jp_end
    avden(2,ip_start:ip_end,j,:) = 0.5*(pressurel%dens(ip_start:ip_end,j,:)+pressurel%dens(ip_start:ip_end,j-1,:))
  enddo
#endif
!
!  Z
!
  do k = 2, nz
    avden(3,ip_start:ip_end,jp_start:jp_end,k) = (fdz0(k)*pressurel%dens(ip_start:ip_end,jp_start:jp_end,k)+fdz0(k-1)*pressurel%dens(ip_start:ip_end,jp_start:jp_end,k-1)) / (fdz0(k)+fdz0(k-1))
  enddo
  avden(3,ip_start:ip_end,jp_start:jp_end,1) = pressurel%dens(ip_start:ip_end,jp_start:jp_end,1)
!
!  Get momentum
!
  u = u*avden(1,:,:,:)
!
#ifdef MODEL_3D
  v = v*avden(2,:,:,:)
#endif
!
  w = w*avden(3,:,:,:)
!
!  BC
!
  call windbc ( u, v, w )
!
#else
!
!  Anelastic case
!
  do k = 1, nz
    u(:,:,k) = u(:,:,k) * den0(k)
#ifdef MODEL_3D
    v(:,:,k) = v(:,:,k) * den0(k)
#endif
    w(:,:,k) = w(:,:,k) * avden0(k)
  enddo
!
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating get_momentum')
!
  end subroutine get_momentum_3
!
! ===============================================
  subroutine get_velocities_1 ( pressurel, windl )
! ===============================================
!
! --- Update velotities using new pressure
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_winds), intent(inout) :: windl
!  
  integer :: k
  integer :: i, j
  real, dimension(3,ip_start:ip_end,jp_start:jp_end,1:nz) :: avden1
!
  if (verbose > 2) call write_debug('Starting get_velocities')
!
!----------------------------------------------------------!
!                Back from momentum rho*U                  !
!----------------------------------------------------------!
!
!  Compressible case
!
#ifndef ANELASTIC
!
  avden1 = 1.
!
!  X
!
  do i = ip_start+1, ip_end
    avden1(1,i,jp_start:jp_end,:) = 2./(pressurel%dens(i,jp_start:jp_end,:)+pressurel%dens(i-1,jp_start:jp_end,:))
  enddo
!
!  Y
!
#ifdef MODEL_3D
  do j = jp_start+1, jp_end
    avden1(2,ip_start:ip_end,j,:) = 2./(pressurel%dens(ip_start:ip_end,j,:)+pressurel%dens(ip_start:ip_end,j-1,:))
  enddo
#endif
!
!  Z
!
  do k = 2, nz
    avden1(3,ip_start:ip_end,jp_start:jp_end,k) = (fdz0(k)+fdz0(k-1))/(fdz0(k)*pressurel%dens(ip_start:ip_end,jp_start:jp_end,k)+fdz0(k-1)*pressurel%dens(ip_start:ip_end,jp_start:jp_end,k-1))
  enddo
  avden1(3,ip_start:ip_end,jp_start:jp_end,1) = 1./pressurel%dens(ip_start:ip_end,jp_start:jp_end,1)
!
!  Get velocities
!
  windl%u = windl%u*avden1(1,:,:,:)
!
#ifdef MODEL_3D
  windl%v = windl%v*avden1(2,:,:,:)
#endif
!
  windl%w = windl%w*avden1(3,:,:,:)
!
!  BC
!
  call windbc ( windl )
!
#else
!
!  Anelastic case
!
  do k = 1, nz
    windl%u(:,:,k) = windl%u(:,:,k) / den0(k)
#ifdef MODEL_3D
    windl%v(:,:,k) = windl%v(:,:,k) / den0(k)
#endif
    windl%w(:,:,k) = windl%w(:,:,k) / avden0(k)
  enddo
!
#endif
!
return
!
!----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating get_velocities')
!
  end subroutine get_velocities_1
!
! ===============================================
  subroutine get_velocities_3 ( pressurel, u, v, w )
! ===============================================
!
! --- Update velotities using new pressure
!
  type (atm_pressure), intent(in) :: pressurel
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: u, v, w
!  
  integer :: k
  integer :: i, j
  real, dimension(3,ip_start:ip_end,jp_start:jp_end,1:nz) :: avden1
!
  if (verbose > 2) call write_debug('Starting get_velocities')
!
!----------------------------------------------------------!
!                Back from momentum rho*U                  !
!----------------------------------------------------------!
!
!  Compressible case
!
#ifndef ANELASTIC
!
  avden1 = 1.
!
!  X
!
  do i = ip_start+1, ip_end
    avden1(1,i,jp_start:jp_end,:) = 2./(pressurel%dens(i,jp_start:jp_end,:)+pressurel%dens(i-1,jp_start:jp_end,:))
  enddo
!
!  Y
!
#ifdef MODEL_3D
  do j = jp_start+1, jp_end
    avden1(2,ip_start:ip_end,j,:) = 2./(pressurel%dens(ip_start:ip_end,j,:)+pressurel%dens(ip_start:ip_end,j-1,:))
  enddo
#endif
!
!  Z
!
  do k = 2, nz
    avden1(3,ip_start:ip_end,jp_start:jp_end,k) = (fdz0(k)+fdz0(k-1))/(fdz0(k)*pressurel%dens(ip_start:ip_end,jp_start:jp_end,k)+fdz0(k-1)*pressurel%dens(ip_start:ip_end,jp_start:jp_end,k-1))
  enddo
  avden1(3,ip_start:ip_end,jp_start:jp_end,1) = 1./pressurel%dens(ip_start:ip_end,jp_start:jp_end,1)
!
!  Get velocities
!
  u = u*avden1(1,:,:,:)
!
#ifdef MODEL_3D
  v = v*avden1(2,:,:,:)
#endif
!
  w = w*avden1(3,:,:,:)
!
!  BC
!
  call windbc ( u, v, w )
!
#else
!
!  Anelastic case
!
  do k = 1, nz
    u(:,:,k) = u(:,:,k) / den0(k)
#ifdef MODEL_3D
    v(:,:,k) = v(:,:,k) / den0(k)
#endif
    w(:,:,k) = w(:,:,k) / avden0(k)
  enddo
!
#endif
!
return
!
!----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating get_velocities')
!
end subroutine get_velocities_3
!
! ===============================================
  subroutine get_conserved ( pressurel, statel, hydromtrl, with_aero )
! ===============================================
!
! --- Calculate momentum from density & velocity
!
!  Input/Output
!
  integer :: k, h
  logical, intent(in), optional :: with_aero
!
  type (atm_state), intent(inout) :: statel
  type (atm_pressure), intent(in) :: pressurel
  type (hydrometeor), dimension(:), intent(inout) :: hydromtrl
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dref
!
  if (verbose > 2) call write_debug('Starting get_conserved')
!
!----------------------------------------------------------!
!                 Density weighted scalars                 !
!----------------------------------------------------------!
!
  dref = pressurel%dens
!
  statel%es = statel%es * dref 
!
  statel%qt = statel%qt * dref 
!
  do h = 1, nscal
    statel%scal(:,:,:,h) = statel%scal(:,:,:,h) * dref 
  enddo
!
  do h = 1, nhydro
    hydromtrl(h)%q = hydromtrl(h)%q * dref
    if ( lndrop > 0 ) hydromtrl(h)%n = hydromtrl(h)%n * dref 
    if ( lmicro > 3 .and. h > ice ) hydromtrl(h)%w = hydromtrl(h)%w * dref 
  enddo
!
  if ( present(with_aero) .and. aero_flg%any ) then
#ifdef AERO_ENABLE
    do h = 1, nmode
      aero3d(h)%n  = aero3d(h)%n * dref 
      aero3d(h)%m  = aero3d(h)%m * dref 
      aero3d(h)%ma = aero3d(h)%ma * dref 
    enddo
#else
    nuc%ccn = nuc%ccn * dref 
    nuc%in  = nuc%in * dref 
#endif
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating get_conserved')
!
  end subroutine get_conserved
!
! ===============================================
  subroutine get_scalars ( pressurel, statel, hydromtrl, with_aero )
! ===============================================
!
! --- Update velotities using new pressure
!
  integer  :: k, h
  logical, intent(in), optional :: with_aero
!
  type (atm_pressure), intent(in) :: pressurel
  type (atm_state), intent(inout) :: statel
  type (hydrometeor), dimension(:), intent(inout) :: hydromtrl
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dref
!
  if (verbose > 2) call write_debug('Starting get_scalars')
!
!----------------------------------------------------------!
!           Back from density weighted scalars             !
!----------------------------------------------------------!
!
  dref = pressurel%dens
!
  statel%es = statel%es / dref
!
  statel%qt = statel%qt / dref
!
  do h = 1, nscal
    statel%scal(:,:,:,h) = statel%scal(:,:,:,h) / dref 
  enddo
!
  do h = 1, nhydro
    hydromtrl(h)%q = hydromtrl(h)%q / dref 
    if ( lndrop > 0 ) hydromtrl(h)%n = hydromtrl(h)%n / dref 
    if ( lmicro > 3 .and. h > ice ) hydromtrl(h)%w = hydromtrl(h)%w / dref 
  enddo
!
  if ( present(with_aero) .and. aero_flg%any ) then
#ifdef AERO_ENABLE
    do h = 1, nmode
      aero3d(h)%n  = aero3d(h)%n / dref 
      aero3d(h)%m  = aero3d(h)%m / dref 
      aero3d(h)%ma = aero3d(h)%ma / dref
    enddo
#else
    nuc%ccn = nuc%ccn / dref 
    nuc%in  = nuc%in / dref 
#endif
  endif
!
return
!
!----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating get_scalars')
!
end subroutine get_scalars
!
! ==================================================
  subroutine damps (c2,x1,x0,xs)
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Damp top layer of scalar x                              !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  character(len=2) :: c2
  real, dimension(nz) :: x0
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x1,xs,xd
!
!----------------------------------------------------------!
!
  xd(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.
!
!  Vertical damping layers
!
  call damps_v (x1,x0,xd)
!
!  Horizontal damping layers
!
  if (dxdamp>0. .or. dydamp>0.) then
    call damps_h (c2,x1,x0,xd)
  endif
!
!  Update
!
  xs = xs + xd
!
!  Diagnostics
!
  if (out_diagu .and. c2=='u ') diag(8)%u = diag(8)%u + cdiag*xd
  if (out_diagv .and. c2=='v ') diag(8)%v = diag(8)%v + cdiag*xd
  if (out_diagw .and. c2=='w ') diag(8)%w = diag(8)%w + cdiag*xd
  if (out_diagt .and. c2=='pt') diag(8)%pt = diag(8)%pt + cdiag*xd
  if (out_diagq .and. c2=='qt') diag(8)%qt = diag(8)%qt + cdiag*xd
  if (out_diags .and. c2=='sc') diag(8)%sca(:,:,:,1) = diag(8)%sca(:,:,:,1) + cdiag*xd
!
!----------------------------------------------------------!
!
  contains

  subroutine damps_v (x1,x0,xs)
!
IMPLICIT NONE  
!
  integer  :: k
  real     :: H, func
  real, dimension(nz) :: x0
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x1,xs
!
  H = zmax - zdamp
  do k = 2, nz
    if (z0(k) >= zdamp) then
      func = 1./tdamp * max(min(cos(0.5*pi*(zmax - z0(k))/H)**2.,1.),0.)
      xs(it_start:it_end,jt_start:jt_end,k) = xs(it_start:it_end,jt_start:jt_end,k) 	&
      	   - func*(x1(it_start:it_end,jt_start:jt_end,k) - x0(k))
    endif
  end do
!
  end subroutine damps_v
!
  subroutine damps_h (c2,x1,x0,xs)
!
IMPLICIT NONE  
!
  integer :: i1, i2, i3, i4, j1, j2, j3, j4, i, j, k
  character(len=2) :: c2
#ifdef NESTING
  integer :: h
  character(len=1) :: c1
#endif
!
  real :: xi, yi, lx, ly, lxs, lxe, lys, lye, func
  real, dimension(nz) :: x0
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x1,xs,x3d
!
!  Sponge layer position
!
#ifdef NESTING
  lx = real(nx-5)*dx
  i1 = it_start
  i2 = i1
  i4 = it_end
  i3 = i4
  do i = it_start, it_end 
    xi = real(i-3)*dx
    if (xi < dxdamp) i2 = i2+1
    if (lx - xi < dxdamp) i3 = i3-1
  enddo
!

  j1=1; j2=1; j3=1; j4=1
#ifdef MODEL_3D
  ly = real(ny-5)*dy
  j1 = jt_start
  j2 = j1
  j4 = jt_end
  j3 = j4
  do j = jt_start, jt_end 
    yi = real(j-3)*dy
    if (yi < dydamp) j2 = j2+1
    if (ly - yi < dydamp) j3 = j3-1
  enddo
#endif
!
  if (nest_run) then
    x3d(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.
!
    select case (c2)
      case ('u ')
        call load_nest (i1, i2, i3, i4, j1, j2, j3, j4, ftime, unest1, unest2, x3d)
      case ('v ')
        call load_nest (i1, i2, i3, i4, j1, j2, j3, j4, ftime, vnest1, vnest2, x3d)
      case ('w ')
        call load_nest (i1, i2, i3, i4, j1, j2, j3, j4, ftime, wnest1, wnest2, x3d)
      case ('pt')
        call load_nest (i1, i2, i3, i4, j1, j2, j3, j4, ftime, ptnest1, ptnest2, x3d)
      case ('qt')
        call load_nest (i1, i2, i3, i4, j1, j2, j3, j4, ftime, qtnest1, qtnest2, x3d)
      case ('q1','q2','q3','q4','q5')
        c1=c2(2:2)
        read(c1,'(i1)') h
        call load_nest (i1, i2, i3, i4, j1, j2, j3, j4, ftime, qnest1(:,:,:,h), qnest2(:,:,:,h), x3d)
      case ('n1','n2','n3','n4','n5')
        c1=c2(2:2)
        read(c1,'(i1)') h
        call load_nest (i1, i2, i3, i4, j1, j2, j3, j4, ftime, nnest1(:,:,:,h), nnest2(:,:,:,h), x3d)
    end select
  endif
#endif
!
!  Spong layer damping in X
!
  lx = real(nx-5)*dx
  lxs = real(it_start-3)*dx
  lxe = real(it_end-3)*dx
!
  if ( lxs <= dxdamp ) then
    do i = it_start, it_end 
      xi = real(i-3)*dx
      if ( xi <= dxdamp ) then
        func = 1./tdamp*max(min(0.25*(1. - cos(pi*(xi/dxdamp - 1.)))**2.,1.),0.)
        do k = 1, nz
          xs(i,:,k) = xs(i,:,k) - func*(x1(i,:,k) - x0(k))
        enddo
      endif
    enddo
  endif
!
  if ( lxe >= lx - dxdamp ) then
    do i = it_start, it_end 
      xi = real(i-3)*dx
      if ( xi >= lx - dxdamp ) then
        func = 1./tdamp*max(min(0.25*(1. - cos(pi*((lx-xi)/dxdamp - 1.)))**2.,1.),0.)
        do k = 1, nz
          xs(i,:,k) = xs(i,:,k) - func*(x1(i,:,k) - x0(k))   
        enddo
      endif
    enddo
  endif
!
!  Spong layer damping in Y
!
#ifdef MODEL_3D
  ly = real(ny-5)*dy
  lys = real(jt_start-3)*dy
  lye = real(jt_end-3)*dy
!
  if ( lys <= dydamp ) then
    do j = jt_start, jt_end 
      yi = real(j-3)*dy
      if ( yi <= dydamp ) then
        func = 1./tdamp*max(min(0.25*(1. - cos(pi*(yi/dydamp - 1.)))**2.,1.),0.)
        do k = 1, nz
          xs(:,j,k) = xs(:,j,k) - func*(x1(:,j,k) - x0(k))
        enddo
      endif
    enddo
  endif
!
  if ( lye >= ly - dydamp ) then
    do j = jt_start, jt_end 
      yi = real(j-3)*dy
      if ( yi >= ly - dydamp ) then
        func = 1./tdamp*max(min(0.25*(1. - cos(pi*((ly-yi)/dydamp - 1.)))**2.,1.),0.)
        do k = 1, nz
          xs(:,j,k) = xs(:,j,k) - func*(x1(:,j,k) - x0(k))
        enddo
      endif
    enddo
  endif
#endif
!
  end subroutine damps_h
!
!----------------------------------------------------------!
!
  end subroutine damps 
!
! ==========================================
  subroutine cal_p1 ( w, wold )
! ==========================================
!
! === Calculates hydrostatic pressure correction
!
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: w, wold
!
  integer :: k
  real, dimension(nz)  :: wav, pp
!
!----------------------------------------------------------!
!
!  Hydrostatic pressure correction
!
  p1 = 0.0
!
  call horav ( (w-wold)/dt0, wav )
  do k = nz-1, 1, -1
    p1(k) = p1(k+1) + wav(k+1)*0.5*dz*(fdz0(k)+fdz0(k+1))
  enddo
!
!  Correct for mass conservation
!
  if (p_mcons) call p_conservation
!
!  Diagnostics 
!  
  if (out_diagw) then
    do k = 1, nz
      diag(8)%w(:,:,k) = diag(8)%w(:,:,k) - cdiag*wav(k)
      diag(9)%w(:,:,k) = diag(9)%w(:,:,k) - cdiag*wav(k)
    enddo
  endif
!
!  Correct velocity
!
  call horav ( w, wav )
!
  do k = 1, nz
    w(:,:,k) = w(:,:,k) - wav(k)
  enddo
!
  return
  end subroutine cal_p1
!
!  ===========================================
   subroutine p_conservation
!  ===========================================
!
! ---------------------------------
! --- Correct pressure obtained from FFTW solver: using FFTW, the
!	pressure is calculated to a constant of integration. The
!       constant can be determined here by invoking conservation of
!	the domain' total internal energy.
! ---------------------------------
!
  integer :: k
  logical :: cons
!
  real :: pcorr, ptot
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: pc, tem
!
!----------------------------------------------------------!
!          Calculate mass conservation correction	   !
!----------------------------------------------------------!
!
#ifdef CONSERVATIVE
    cons = .true.
#else
    cons = .false.
#endif
!
    call get_temperature ( pressure, state, hydromtr, tem )
!
    do k = 1, nz
      pc(:,:,k) = den0(k)*(thermo_prop%cp(:,:,k) - thermo_prop%cv(:,:,k))*(tem(:,:,k) - t0(k))
    enddo
! 
    call totalav ( pressure%p, ptot )
!
    call totalav ( pc, pcorr )
!
!  Apply pressure correction
!
    p1 = p1 - ptot + pcorr
!
! ---------------------------------------------------------!
!
  return
  end subroutine p_conservation
!
! ==================================================
  subroutine damp_div
! ==================================================
!
! ---------------------------------
! --- Calculate divergence damping 
! ---------------------------------
!
  integer :: k
  real :: H, adiv
  real, dimension(nz)  :: cuv
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: gdivx, gdivy, gdivz, div, zero
!
! ---------------------------------------------------------!
!
  zero = 0.0
  cuv = 0.0
!
!  Height dependent coefficient 
!
  adiv = 0.2*dx*dy/dt0
  H = zmax - zdamp
  do k = 1, nz
    if (z0(k) >= zdamp) then
      cuv(k) = adiv*(1. + 15.*max(min(cos(0.5*pi*(zmax - z0(k))/H)**2.,1.),0.))
    else
      cuv(k) = adiv
    endif
  enddo
!
!  Divergence damping 
!
  call caldiv ( wind, div, horizontal=.true. )
!
  call gradp ( div, zero, gdivx, gdivy, gdivz )
!
  do k = 1, nz
    winds%u(:,:,k) = winds%u(:,:,k) + cuv(k)*gdivx(:,:,k)
#ifdef MODEL_3D
    winds%v(:,:,k) = winds%v(:,:,k) + cuv(k)*gdivy(:,:,k)
#endif
  enddo
!
!  Diagnostics
!
  if (ldiag) then
    do k = 1, nz
      if (out_diagu) diag(8)%u(:,:,k) = diag(8)%u(:,:,k) + cdiag*cuv(k)*gdivx(:,:,k)
#ifdef MODEL_3D
      if (out_diagv) diag(8)%v(:,:,k) = diag(8)%v(:,:,k) + cdiag*cuv(k)*gdivy(:,:,k)
#endif
    enddo
  endif
!
! ---------------------------------------------------------!
!
return
end subroutine damp_div

end module funcpack
