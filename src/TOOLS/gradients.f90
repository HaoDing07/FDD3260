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
!  Purpose:
!	Package of routines calculating scalar gradients at face centers
!   based one several (FV) schemes
!
!  Author:
!	Julien Savre, MIM-LMU
!
! ================================================================
!
  module gradients
!
  USE gridno
  USE allocation
  USE shared_data	
  USE shared_diag
  USE shared_wind
  USE shared_pressure
  USE boundary_conditions
!
  IMPLICIT NONE
!
  private
!
  real, parameter :: small=1.e-18
!
  interface grad_wind
    MODULE PROCEDURE grad_wind_2d, grad_wind_3d
  end interface grad_wind
!
  public :: grad_quick, grad_cus, grad_wind, grad_1d, gradp, cal_dij, caldiv, calconv, cal_rot, u_gradp
!
  contains
!
! ================================================================
!
! ======================================================
  subroutine grad_quick (x, gradx, grady, gradz)
! ======================================================
!
!----------------------------------------------------------!
!                                                          !
!  Gradient calculation based on QUICK method.       	   !
!                                                          !
!----------------------------------------------------------!
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x 
  real, dimension(ip_start:,jp_start:,:), optional, intent(out) :: gradx, grady, gradz
!
  integer :: i, j, k, i1, i2, j1, j2
!
  real :: dx1, dy1
  real, dimension(ip_start:ip_end) :: xx, gx
  real, dimension(jp_start:jp_end) :: xy, gy
  real, dimension(1:nz) :: xz, gz
! 
!----------------------------------------------------------!
!
!  Initialization
!
  dx1 = 1./dx
#ifdef MODEL_3D
  dy1 = 1./dy
#endif
!
  i1=it_start-1
  i2=it_end
  j1=1
  j2=1
#ifdef MODEL_3D
  j1=jt_start-1
  j2=jt_end
#endif
!
  gx = 0.
  gy = 0.
  gz = 0.
!
!----------------------------------------------------------!
!      		    Compute gradients    		   !
!----------------------------------------------------------!
!
!  Gradient in x
!
  if (present(gradx)) then
!
  do k=1,nz
    do j=jt_start,jt_end
      xx = x(ip_start:ip_end,j,k)   
!
      do i=i1,i2
        gx(i) = (xx(i+1) - xx(i))*dx1
      enddo
!
!  Move gradient by 1 cell when copying
!
      gradx(i1+1:i2+1,j,k) = gx(i1:i2)
    enddo
  enddo
!
  call gradbc ( gradx(:,:,:) )
!
  endif
!
!  Gradient in y
!
#ifdef MODEL_3D
  if (present(grady)) then
!
  do k=1,nz
    do i=it_start,it_end
      xy = x(i,jp_start:jp_end,k)   
!
      do j=j1,j2
        gy(j) = (xy(j+1) - xy(j))*dy1
      enddo
!
!  Move gradient by 1 cell when copying
!
      grady(i,j1+1:j2+1,k) = gy(j1:j2)
!
#ifdef CHANNEL
      if (jp_start == 1) then
        grady(i,j1+1,k) = 0.
      endif
      if (jp_end == ny) then
        grady(i,j2+1,k) = 0.
      endif
#endif
!
    enddo
  enddo
!
  call gradbc ( grady(:,:,:) )
!
  endif
#endif
!
!  Gradient in z
!
  if (present(gradz)) then
!
  do j=jt_start,jt_end
    do i=it_start,it_end
      xz = x(i,j,1:nz)
!
      do k=1,nz-1
        gz(k) = (xz(k+1) - xz(k))*dz1(k)
      enddo
!
!  Move gradient by 1 cell when copying
!
      gradz(i,j,2:nz) = gz(1:nz-1)
    enddo
  enddo
!
  call gradbc ( gradz(:,:,:) )
!
  endif
!  
  return
!
  end subroutine grad_quick
!
! ======================================================
  subroutine grad_cus (x, gradx, grady, gradz)
! ======================================================
!
!----------------------------------------------------------!
!                                                          !
!  Gradient calculation using cubic interpolations.        !
!  (cubic upwind scheme)				   !
!                                                          !
!----------------------------------------------------------!
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x 
  real, dimension(ip_start:,jp_start:,:), optional, intent(out) :: gradx, grady, gradz
!
  integer :: i, j, k, i1, i2, j1, j2
!
  real    :: dx1, dy1
  real    :: ffip1, ffip2, ffim1, deltac, deltad, deltau, deltadd
  real    :: phi0, phi1, phi2, phi3
  
  real, dimension(ip_start:ip_end) :: xx, gx
  real, dimension(jp_start:jp_end) :: xy, gy
  real, dimension(1:nz) :: xz, gz
! 
!----------------------------------------------------------!
!
!  Initialization
!
  dx1 = 1./dx
  dy1 = 1./dy
!
  i1=it_start-1
  i2=it_end
  j1=1
  j2=1
#ifdef MODEL_3D
  j1=jt_start-1
  j2=jt_end
#endif
!
  gx = 0.
  gy = 0.
  gz = 0.
!
!  Compute gradient
!
!  Gradient in x
!
  if (present(gradx)) then
!
  do k=1,nz
    do j=jt_start,jt_end
      xx = x(ip_start:ip_end,j,k)   
!
      do i=i1,i2
    	ffip1  = (xx(i+1) - xx(i))
	ffip2  = (xx(i+2) - xx(i+1))
	ffim1  = (xx(i) - xx(i-1))
!
        call coeff_cus (xx(i), ffip1, ffim1, ffip2, dx, dx, dx, dx, phi0, phi1, phi2, phi3)
!
        gx(i) = (phi1 + phi2 + 0.75*phi3)*dx1
      enddo
!
      gradx(i1+1:i2+1,j,k) = gx(i1:i2)
    enddo
  enddo
!
  call gradbc ( gradx(:,:,:) )
!
  endif
!
!  Gradient in y
!
#ifdef MODEL_3D
  if (present(grady)) then
!
  do k=1,nz
    do i=it_start,it_end
      xy = x(i,jp_start:jp_end,k)   
!
      do j=j1,j2
    	ffip1  = (xy(j+1) - xy(j))
	ffip2  = (xy(j+2) - xy(j+1))
	ffim1  = (xy(j) - xy(j-1))
!
        call coeff_cus (xy(j), ffip1, ffim1, ffip2, dy, dy, dy, dy, phi0, phi1, phi2, phi3)
!
        gy(j) = (phi1 + phi2 + 0.75*phi3)*dy1
      enddo
!
      grady(i,j1+1:j2+1,k) = gy(j1:j2)
!
#ifdef CHANNEL
      if (jp_start == 1) then
        grady(i,j1+1,k) = 0.
      endif
      if (jp_end == ny) then
        grady(i,j2+1,k) = 0.
      endif
#endif
!
    enddo
  enddo
!
  call gradbc ( grady(:,:,:) )
!
  endif
#endif
!
!  Gradient in z
!
  if (present(gradz)) then
!
  do j=j1,j2
    do i=i1,i2
      xz = x(i,j,1:nz)
!
!  Top boundary
!
      k = nz-1
      deltac  = dz * fdz0(k)
      deltad  = 0.5 * dz * (fdz0(k+1) + fdz0(k))
      deltadd = dz * fdz0(k+1)
      deltau  = 0.5 * dz * (fdz0(k) + fdz0(k-1))
!
      ffip1  = (xz(k+1) - xz(k))
      ffip2  = 0.
      ffim1  = (xz(k) - xz(k-1))
!
      call coeff_cus (xz(k), ffip1, ffim1, ffip2, deltac, deltau, deltad, deltadd, phi0, phi1, phi2, phi3)

      gz(k) = (phi1 + phi2 + 0.75*phi3)*dz1(k)
!
!  Bottom boundary (level 2)
!
      k = 1
      deltac  = dz * fdz0(k)
      deltad  = 0.5 * dz * (fdz0(k+1) + fdz0(k))
      deltadd = 0.5 * dz * (fdz0(k+2) + fdz0(k+1))
      deltau  = dz * fdz0(k)
!
      ffip1  = (xz(k+1) - xz(k))
      ffip2  = (xz(k+2) - xz(k+1))
      ffim1  = 0.
!
      call coeff_cus (xz(k), ffip1, ffim1, ffip2, deltac, deltau, deltad, deltadd, phi0, phi1, phi2, phi3)
!
      gz(k) = (phi1 + phi2 + 0.75*phi3)*dz1(k)
!
!  Domain core
!
      do k=2,nz-2
        deltac  = dz * fdz0(k)
        deltad  = 0.5 * dz * (fdz0(k+1) + fdz0(k))
        deltadd = 0.5 * dz * (fdz0(k+2) + fdz0(k+1))
        deltau  = 0.5 * dz * (fdz0(k) + fdz0(k-1))
!
    	ffip1  = (xz(k+1) - xz(k))
	ffip2  = (xz(k+2) - xz(k+1))
	ffim1  = (xz(k) - xz(k-1))
!
        call coeff_cus (xz(k), ffip1, ffim1, ffip2, deltac, deltau, deltad, deltadd, phi0, phi1, phi2, phi3)
!
        gz(k) = (phi1 + phi2 + 0.75*phi3)*dz1(k)
      enddo
!
      gradz(i,j,2:nz) = gz(1:nz-1)
    enddo
  enddo
!
  call gradbc ( gradz(:,:,:) )
!
  endif
!  
  return

  contains
!
  subroutine coeff_cus (phic, ffip1, ffim1, ffip2, deltac, deltau, deltad, deltadd, phi0, phi1, phi2, phi3)
!
!----------------------------------------------------------!
!                                                          !
!  Calculate coefficients for cubic upwind scheme.   	   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real  :: phic, ffip1, ffip2, ffim1
  real  :: deltac, deltau, deltad, deltadd
  real  :: phi0, phi1, phi2, phi3
  real  :: g1, g2, g3
!                                                          !
!----------------------------------------------------------!
!
  g1 = deltac*deltac / (deltau+deltad+deltadd)
  g2 = 1. / (deltadd*deltad*(deltadd+deltad))
  g3 = 1. / (deltau*deltad*(deltau+deltad))
!
  phi0 = phic
  phi3 = g1*((deltac*deltad*ffip2 - deltac*deltadd*ffip1)*g2 - (deltac*deltau*ffip1 - deltac*deltad*ffim1)*g3)
  phi1 = (deltac*deltau*deltau*ffip1 + deltac*deltad*deltad*ffim1)*g3 - phi3*deltad*deltau / (deltac*deltac)
  phi2 = (deltac*deltac*deltau*ffip1 - deltac*deltac*deltad*ffim1)*g3 - phi3*(deltad - deltau) / deltac
!
  return
  end
!
  end subroutine grad_cus
!
! ======================================================
  subroutine grad_1d (xz, gz)
! ======================================================
!
!----------------------------------------------------------!
!                                                          !
!  Gradient calculation based on QUICK method.       	   !
!                                                          !
!----------------------------------------------------------!
!
  real, dimension(:), intent(in) :: xz
  real, dimension(:), intent(out) :: gz
!
  integer :: k
  real :: ffim1, ffip1, ffip2, deltac, deltad, deltau
  real :: phi1, phi2
! 
!----------------------------------------------------------!
!      		    Compute gradients    		   !
!----------------------------------------------------------!
!
  gz = 0.
!
!  Top boundary
!
  k = nz-1
  deltad  = 0.5 * dz * (fdz0(k+1) + fdz0(k))
  deltau  = 0.5 * dz * (fdz0(k) + fdz0(k-1))
  deltac  = (deltad + deltau)*(deltad*deltad + deltau*deltau) - (deltad - deltau)*(deltad*deltad - deltau*deltau)
!
  ffip1 = (xz(k+1) - xz(k))
  ffip2 = 0.
  ffim1 = (xz(k) - xz(k-1))
!
  phi1 = ((ffip1 + ffim1)*(deltad*deltad + deltau*deltau) - (ffip1 - ffim1)*(deltad*deltad - deltau*deltau)) / deltac
  phi2 = ((ffip1 - ffim1)*(deltad + deltau) - (ffip1 + ffim1)*(deltad - deltau)) / deltac
!
  gz(k) = phi1
!
!  Bottom boundary (level 2)
!
  k = 1
  deltad  = 0.5 * dz * (fdz0(k+1) + fdz0(k))
  deltau  = dz * fdz0(k)
  deltac  = (deltad + deltau)*(deltad*deltad + deltau*deltau) - (deltad - deltau)*(deltad*deltad - deltau*deltau)
!
  ffip1 = (xz(k+1) - xz(k))
  ffip2 = (xz(k+2) - xz(k+1))
  ffim1 = 0.
!
  phi1 = ((ffip1 + ffim1)*(deltad*deltad + deltau*deltau) - (ffip1 - ffim1)*(deltad*deltad - deltau*deltau)) / deltac
  phi2 = ((ffip1 - ffim1)*(deltad + deltau) - (ffip1 + ffim1)*(deltad - deltau)) / deltac
!
  gz(k) = phi1
!
!  Domain core
!
  do k=2,nz-2
    deltad  = 0.5 * dz * (fdz0(k+1) + fdz0(k))
    deltau  = 0.5 * dz * (fdz0(k) + fdz0(k-1))
    deltac  = (deltad + deltau)*(deltad*deltad + deltau*deltau) - (deltad - deltau)*(deltad*deltad - deltau*deltau)
!
    ffip1  = (xz(k+1) - xz(k))
    ffip2  = (xz(k+2) - xz(k+1))
    ffim1  = (xz(k) - xz(k-1))
!
    phi1 = ((ffip1 + ffim1)*(deltad*deltad + deltau*deltau) - (ffip1 - ffim1)*(deltad*deltad - deltau*deltau)) / deltac
    phi2 = ((ffip1 - ffim1)*(deltad + deltau) - (ffip1 + ffim1)*(deltad - deltau)) / deltac
!
    gz(k) = phi1
  enddo
!  
  return
!
  end subroutine grad_1d
!
! ======================================================
  subroutine cal_dij ( windl, dij, rij )
! ======================================================		  
!
! ---------------------------------------------------------
! --- Subroutine for calculations of the deformation tensor
! --- Output: dij is cell centered 
! ---------------------------------------------------------
!
  type(atm_winds), intent(in) :: windl
  real, dimension(ip_start:,jp_start:,:,:,:), optional, intent(inout) :: dij, rij
!
  integer :: i1, i2, j1, j2, i, j, k
  real, dimension(:,:,:), allocatable :: gradux, graduy, graduz, gradvx, gradvy, gradvz, gradwx, gradwy, gradwz
!
  if (verbose > 2) call write_debug('Starting cal_dij')
!
!----------------------------------------------------------!
!
  call alloc ( gradux )
  call alloc ( graduy )
  call alloc ( graduz )
#ifdef MODEL_3D
  call alloc ( gradvx )
  call alloc ( gradvy )
  call alloc ( gradvz )
#endif
  call alloc ( gradwx )
  call alloc ( gradwy )
  call alloc ( gradwz )
!
  if (present(dij)) dij = 0.0
  if (present(rij)) rij = 0.0
!
!----------------------------------------------------------!
!                    Calculate gradient 		   !
!----------------------------------------------------------!
!
#ifdef MODEL_3D
  call grad_wind_3d (windl, gradux, graduy, graduz, gradvx, gradvy, gradvz, gradwx, gradwy, gradwz)
#else
  call grad_wind_2d (windl, gradux, graduz, gradwx, gradwz)
#endif
!
!----------------------------------------------------------!
!         Diagonal terms (already cell centered)	   !
!----------------------------------------------------------!
!
  if (present(dij)) then
    dij(:,:,:,1,1) = gradux
#ifdef MODEL_3D
    dij(:,:,:,2,2) = gradvy
#endif
    dij(:,:,:,3,3) = gradwz
  endif
!
!----------------------------------------------------------!
!                     Non-diagonal terms      	           !
!----------------------------------------------------------!
!
  i1=it_start-1
  i2=it_end+1
  j1=1
  j2=1
#ifdef MODEL_3D
  j1=jt_start-1
  j2=jt_end+1
#endif
!
!  Symmetric part
!
  if (present(dij)) then
!
  do k = 2, nz-1
    do j = j1, j2
      do i = i1, i2
        dij(i,j,k,1,3) = 0.5*( 0.25*(graduz(i,j,k)+graduz(i+1,j,k)+graduz(i,j,k-1)+graduz(i+1,j,k-1)) &
                       + 0.25*(gradwx(i-1,j,k+1)+gradwx(i,j,k+1)+gradwx(i-1,j,k)+gradwx(i,j,k)) )
        dij(i,j,k,3,1) = dij(i,j,k,1,3)
!
#ifdef MODEL_3D
        dij(i,j,k,1,2) = 0.5*( 0.25*(graduy(i,j,k)+graduy(i+1,j,k)+graduy(i,j-1,k)+graduy(i+1,j-1,k)) &
                       + 0.25*(gradvx(i-1,j+1,k)+gradvx(i,j+1,k)+gradvx(i-1,j,k)+gradvx(i,j,k)) )
        dij(i,j,k,2,1) = dij(i,j,k,1,2) 
!
        dij(i,j,k,2,3) = 0.5*( 0.25*(gradvz(i,j,k)+gradvz(i,j+1,k)+gradvz(i,j,k-1)+gradvz(i,j+1,k-1)) & 
                       + 0.25*(gradwy(i,j-1,k+1)+gradwy(i,j,k+1)+gradwy(i,j-1,k)+gradwy(i,j,k)) )
        dij(i,j,k,3,2) = dij(i,j,k,2,3) 
#endif
      end do
    end do
  end do
!
  endif
!
!  Antisymmetric part
!
  if (present(rij)) then
!
  do k = 2, nz-1
    do j = j1, j2
      do i = i1, i2
        rij(i,j,k,1,3) = 0.5*( 0.25*(graduz(i,j,k)+graduz(i+1,j,k)+graduz(i,j,k-1)+graduz(i+1,j,k-1)) &
                       - 0.25*(gradwx(i-1,j,k+1)+gradwx(i,j,k+1)+gradwx(i-1,j,k)+gradwx(i,j,k)) )
        rij(i,j,k,3,1) = -rij(i,j,k,1,3)
!
#ifdef MODEL_3D
        rij(i,j,k,1,2) = 0.5*( 0.25*(graduy(i,j,k)+graduy(i+1,j,k)+graduy(i,j-1,k)+graduy(i+1,j-1,k)) &
                       - 0.25*(gradvx(i-1,j+1,k)+gradvx(i,j+1,k)+gradvx(i-1,j,k)+gradvx(i,j,k)) )
        rij(i,j,k,2,1) = -rij(i,j,k,1,2) 
!
        rij(i,j,k,2,3) = 0.5*( 0.25*(gradvz(i,j,k)+gradvz(i,j+1,k)+gradvz(i,j,k-1)+gradvz(i,j+1,k-1)) & 
                       - 0.25*(gradwy(i,j-1,k+1)+gradwy(i,j,k+1)+gradwy(i,j-1,k)+gradwy(i,j,k)) )
        rij(i,j,k,3,2) = -rij(i,j,k,2,3) 
#endif
      end do
    end do
  end do
!
  endif
!
  call dealloc ( gradux )
  call dealloc ( graduy )
  call dealloc ( graduz )
#ifdef MODEL_3D
  call dealloc ( gradvx )
  call dealloc ( gradvy )
  call dealloc ( gradvz )
#endif
  call dealloc ( gradwx )
  call dealloc ( gradwy )
  call dealloc ( gradwz )
!
  if (verbose > 2) call write_debug('Terminating cal_dij')
!
  return
  end subroutine
!
! ======================================================
  subroutine cal_rot ( windl, rotx, roty, rotz )
! ======================================================		  
!
! ---------------------------------------------------------
! --- Subroutine for calculations of the deformation tensor
! --- Output: dij is cell centered 
! --------------------------------------------------------
!
  type(atm_winds), intent(in) :: windl
  real, dimension(ip_start:,jp_start:,:), intent(out) :: rotx, roty, rotz
  
  real, dimension(:,:,:), allocatable :: gradux, graduy, graduz, gradvx, gradvy, gradvz, gradwx, gradwy, gradwz
!
!----------------------------------------------------------!
!
  call alloc ( gradux )
  call alloc ( graduy )
  call alloc ( graduz )
#ifdef MODEL_3D
  call alloc ( gradvx )
  call alloc ( gradvy )
  call alloc ( gradvz )
#endif
  call alloc ( gradwx )
  call alloc ( gradwy )
  call alloc ( gradwz )
!
  rotx = 0.0
  roty = 0.0
  rotz = 0.0
!
!----------------------------------------------------------!
!                    Calculate gradient 		   !
!----------------------------------------------------------!
!
#ifdef MODEL_3D
  call grad_wind_3d (windl, gradux, graduy, graduz, gradvx, gradvy, gradvz, gradwx, gradwy, gradwz)
#else
  call grad_wind_2d (windl, gradux, graduz, gradwx, gradwz)
#endif
!
!----------------------------------------------------------!
!          	   Calculate rotational			   !
!----------------------------------------------------------!
!
  rotx = gradwy - gradvz
  roty = graduz - gradwx
  rotz = gradvx - graduy
!
  call dealloc ( gradux )
  call dealloc ( graduy )
  call dealloc ( graduz )
#ifdef MODEL_3D
  call dealloc ( gradvx )
  call dealloc ( gradvy )
  call dealloc ( gradvz )
#endif
  call dealloc ( gradwx )
  call dealloc ( gradwy )
  call dealloc ( gradwz )
!
  return
  end subroutine cal_rot
!
! ======================================================
  subroutine grad_wind_2d (windl, gradux, graduz, gradwx, gradwz)
! ======================================================
!
!----------------------------------------------------------!
!                                                          !
!  Gradient calculation for wind vector using 4th order    !
!  central finite difference formula on staggered grid     !
!							   !
!  Careful! Diagonal terms (grad(1,1) grad(2,2) grad(3,3)) !
!  are computed at cell centers while off diagonal terms   !
!  are defined at face centers (at the velocity node).	   !
!                                                          !
!----------------------------------------------------------!
!
  integer :: i, j, k, i1, i2, j1, j2
  real :: dx1
!
  type(atm_winds), intent(in) :: windl
  real, dimension(ip_start:,jp_start:,1:), intent(out) 	&
  	:: gradux, graduz, gradwx, gradwz
!
  real, dimension(ip_start:ip_end) :: uu, uw, gu1, gw1
  real, dimension(1:nz) :: wu, ww, gu3, gw3
! 
!----------------------------------------------------------!
!
!  Initialization
!
  dx1 = 1./dx
!
  i1=it_start
  i2=it_end
  j1=1
  j2=1
!
  gu1 = 0.; gu3 = 0.
  gw1 = 0.; gw3 = 0.
!
!  Gradients in x
!
  do k=1,nz
    do j=j1,j2
      uu = windl%u(:,j,k)
      uw = windl%w(:,j,k)
!
      do i=i1,i2
    	gu1(i) = 9./8.*(uu(i+1) - uu(i))*dx1 - 1./24.*(uu(i+2) - uu(i-1))*dx1
!
    	gw1(i) = 9./8.*(uw(i+1) - uw(i))*dx1 - 1./24.*(uw(i+2) - uw(i-1))*dx1
      enddo
!
      gradux(:,j,k) = gu1
      gradwx(:,j,k) = gw1
    enddo
  enddo
!
!  Gradient in z
!
!  Bottom boundary (level 2)
!
  do j=j1,j2
    do i=i1,i2
      wu = windl%u(i,j,:)
      ww = windl%w(i,j,:)
!
!  Main domain core
!
      do k=2,nz-2
    	gu3(k)  = 9./8.*(wu(k+1) - wu(k))*dz1(k) - 1./8.*(wu(k+2) - wu(k-1)) / ( dz*(0.5*fdz0(k-1)+fdz0(k)+fdz0(k+1)+0.5*fdz0(k+2)) )
!
    	gw3(k)  = 9./8.*(ww(k+1) - ww(k))*dz1w(k) - 1./8.*(ww(k+2) - ww(k-1)) / ( dz*(fdz0(k-1)+fdz0(k)+fdz0(k+1)) )
      enddo
!
!  boundaries
!
      gu3(1) = (wu(2) - wu(1))*dz1(1)
!
      gu3(nz-1) = (wu(nz) - wu(nz-1))*dz1(nz) 
!
      gw3(nz-1) = 9./8.*(ww(nz) - ww(nz-1))*dz1w(nz-1) - 1./8.*(0. - ww(nz-2)) / ( dz*(fdz0(nz)+fdz0(nz-1)+fdz0(nz-2)) )
!
      if (bcl(6) == 'nnf') then
        gu3(nz) = gu3(nz-1)
!
        gw3(nz) = (0. - ww(nz))*dz1w(nz)
      endif
!
      graduz(i,j,:) = gu3(:)
      gradwz(i,j,:) = gw3(:)
    enddo
  enddo
!
!  Apply boundary conditions to gradient
!
  call gradbc ( gradux )
!
  call gradbc ( gradwz )
!
  call gradbc ( graduz )
!
  call gradbc ( gradwx )
!
!  Diagnostics
!
  if (ldiag .and. out_grad) diag(1)%grad%u = gradux 
  if (ldiag .and. out_grad) diag(3)%grad%u = graduz 
!
  if (ldiag .and. out_grad) diag(1)%grad%w = gradwx 
  if (ldiag .and. out_grad) diag(3)%grad%w = gradwz
!
  return
!
  end subroutine grad_wind_2d
!
! ======================================================
  subroutine grad_wind_3d (windl, gradux, graduy, graduz, gradvx, gradvy, gradvz, gradwx, gradwy, gradwz)
! ======================================================
!
!----------------------------------------------------------!
!                                                          !
!  Gradient calculation for wind vector using 4th order    !
!  central finite difference formula on staggered grid     !
!							   !
!  Careful! Diagonal terms (grad(1,1) grad(2,2) grad(3,3)) !
!  are computed at cell centers while off diagonal terms   !
!  are defined at the corners                              !
!                                                          !
!----------------------------------------------------------!
!
  integer :: i, j, k, i1, i2, j1, j2
  real :: dx1, dy1
!
  type(atm_winds), intent(in) :: windl
  real, dimension(ip_start:,jp_start:,1:), intent(out) 	&
  	:: gradux, graduy, graduz, gradvx, gradvy, gradvz, gradwx, gradwy, gradwz
!
  real, dimension(ip_start:ip_end) :: uu, uv, uw, gu1, gv1, gw1
  real, dimension(jp_start:jp_end) :: vu, vv, vw, gu2, gv2, gw2
  real, dimension(1:nz) :: wu, wv, ww, gu3, gv3, gw3
! 
!----------------------------------------------------------!
!
!  Initialization
!
  dx1 = 1./dx
  dy1 = 1./dy
!
  i1=it_start
  i2=it_end
  j1=jt_start
  j2=jt_end
!
  gu1 = 0.; gu2 = 0.; gu3 = 0.
  gv1 = 0.; gv2 = 0.; gv3 = 0.
  gw1 = 0.; gw2 = 0.; gw3 = 0.
!
!  Gradients in x
!
  do k=1,nz
    do j=j1,j2
      uu = windl%u(:,j,k)
      uv = windl%v(:,j,k)
      uw = windl%w(:,j,k)
!
      do i=i1,i2
    	gu1(i) = 9./8.*(uu(i+1) - uu(i))*dx1 - 1./24.*(uu(i+2) - uu(i-1))*dx1
!
    	gv1(i) = 9./8.*(uv(i+1) - uv(i))*dx1 - 1./24.*(uv(i+2) - uv(i-1))*dx1
!
    	gw1(i) = 9./8.*(uw(i+1) - uw(i))*dx1 - 1./24.*(uw(i+2) - uw(i-1))*dx1
      enddo
!
      gradux(:,j,k) = gu1
      gradvx(:,j,k) = gv1
      gradwx(:,j,k) = gw1
    enddo
  enddo
!
!  Gradients in y
!
  do k=1,nz
    do i=i1,i2
      vu = windl%u(i,:,k)
      vv = windl%v(i,:,k)
      vw = windl%w(i,:,k)
!
      do j=j1,j2
    	gu2(j) = 9./8.*(vu(j+1) - vu(j))*dy1 - 1./24.*(vu(j+2) - vu(j-1))*dy1
!
    	gv2(j) = 9./8.*(vv(j+1) - vv(j))*dy1 - 1./24.*(vv(j+2) - vv(j-1))*dy1
!
    	gw2(j) = 9./8.*(vw(j+1) - vw(j))*dy1 - 1./24.*(vw(j+2) - vw(j-1))*dy1
      enddo
!
#ifdef CHANNEL
      if (jp_start == 1) then
        gu2(j1) = (vu(j1+1) - vu(j1))*dy1 
!
        gv2(j1) = (vv(j1+1) - 0.)*dy1
!
        gw2(j1) = (vw(j1+1) - vw(j1))*dy1
!
        gv2(j1+1)  = 9./8.*(vv(j1+2) - vv(j1+1))*dy1 - 1./24.*(vv(j1+3) - 0.)*dy1
      endif
!
      if (jp_end == ny) then
        gu2(j2) = gu2(j2-1) 
!
        gv2(j2) = (0. - vv(j2))*dy1
!
        gw2(j2) = gw2(j2-1) 
!
        gv2(j2-1)  = 9./8.*(wv(j2) - wv(j2-1))*dy1 - 1./24.*(0. - wv(j2-2))*dy1
      endif
#endif
!
      graduy(i,:,k) = gu2
      gradvy(i,:,k) = gv2
      gradwy(i,:,k) = gw2
    enddo
  enddo
!    
!  Gradient in z
!
  do j=j1,j2
    do i=i1,i2
      wu = windl%u(i,j,:)
      wv = windl%v(i,j,:)
      ww = windl%w(i,j,:)
!
!  Bottom boundary
!
      do k=2,nz-2
    	gu3(k)  = 9./8.*(wu(k+1) - wu(k))*dz1(k) - 1./8.*(wu(k+2) - wu(k-1))/(dz*(0.5*fdz0(k+2)+fdz0(k)+fdz0(k-1)+0.5*fdz0(k-1)))
!
    	gv3(k)  = 9./8.*(wv(k+1) - wv(k))*dz1(k) - 1./8.*(wv(k+2) - wv(k-1))/(dz*(0.5*fdz0(k+2)+fdz0(k)+fdz0(k+1)+0.5*fdz0(k-1)))
!
    	gw3(k)  = 9./8.*(ww(k+1) - ww(k))*dz1w(k) - 1./8.*(ww(k+2) - ww(k-1)) / ( dz*(fdz0(k-1)+fdz0(k)+fdz0(k+1)) )
      enddo
!
      gu3(1) = (wu(2) - wu(1))*dz1(1) 
!
      gv3(1) = (wv(2) - wv(1))*dz1(1)
!
      gw3(1) = (ww(2) - 0.)*dz1w(1)
!
!  Top boundary
!
      gu3(nz-1) = (wu(nz) - wu(nz-1))*dz1(nz-1)
!
      gv3(nz-1) = (wv(nz) - wv(nz-1))*dz1(nz-1)
!
      gw3(nz-1) = 9./8.*(ww(nz) - ww(nz-1))*dz1w(nz-1) - 1./8.*(0. - ww(nz-2)) / ( dz*(fdz0(nz)+fdz0(nz-1)+fdz0(nz-2)) )
!
      if (bcl(6) == 'nnf') then
        gu3(nz) = gu3(nz-1) 
!
        gv3(nz) = gv3(nz-1)
!
        gw3(nz) = (0. - ww(nz))*dz1w(nz)
      endif
!
      graduz(i,j,:) = gu3(:)
      gradvz(i,j,:) = gv3(:)
      gradwz(i,j,:) = gw3(:)
    enddo
  enddo
!
!  Apply boundary conditions to gradient
!
  call gradbc ( gradux )
!
  call gradbc ( gradvy )
!
  call gradbc ( gradwz )
!
  call gradbc ( graduz )
!
  call gradbc ( gradwx )
!
  call gradbc ( gradvx )
!
  call gradbc ( graduy )
!
  call gradbc ( gradvz )
!
  call gradbc ( gradwy )
!
!  Diagnostics
!
  if (ldiag .and. out_grad) diag(1)%grad%u = gradux 
  if (ldiag .and. out_grad) diag(2)%grad%u = graduy 
  if (ldiag .and. out_grad) diag(3)%grad%u = graduz 
!
  if (ldiag .and. out_grad) diag(1)%grad%v = gradvx 
  if (ldiag .and. out_grad) diag(2)%grad%v = gradvy
  if (ldiag .and. out_grad) diag(3)%grad%v = gradvz
!
  if (ldiag .and. out_grad) diag(1)%grad%w = gradwx 
  if (ldiag .and. out_grad) diag(2)%grad%w = gradwy
  if (ldiag .and. out_grad) diag(3)%grad%w = gradwz
!  
  return

  end subroutine grad_wind_3d
!
!  ===========================================
   subroutine caldiv ( windl, div, horizontal  )
!  ===========================================
!
! ---------------------------------
! --- Calculate pressure gradients
! ---------------------------------
!
  integer :: i, j, k
  real  :: dx1, dy1
!
  logical, intent(in), optional :: horizontal
  type(atm_winds), intent(inout) :: windl
  real, dimension(ip_start:,jp_start:,:), intent(out) :: div
!
  if (verbose > 2) call write_debug('Starting caldiv')
!
! ---------------------------------------------------------!
!
  dx1 = 1./dx
#ifdef MODEL_3D
  dy1 = 1./dy 
#endif
!
  div = 0.0
!
  call windbc ( windl )
!
!----------------------------------------------------------!
!             Calculate velocity divergence                !
!----------------------------------------------------------!
!
  do i = it_start,it_end
    div(i,jt_start:jt_end,1:nz) = div(i,jt_start:jt_end,1:nz) + (windl%u(i+1,jt_start:jt_end,1:nz) - windl%u(i,jt_start:jt_end,1:nz))*dx1
  enddo
!
#ifdef MODEL_3D
  do j = jt_start,jt_end
    div(it_start:it_end,j,1:nz) = div(it_start:it_end,j,1:nz) + (windl%v(it_start:it_end,j+1,1:nz) - windl%v(it_start:it_end,j,1:nz))*dy1
  enddo
#endif
!
  if ( .not.present(horizontal) ) then
!
    do k = 1, nz-1
      div(it_start:it_end,jt_start:jt_end,k) = div(it_start:it_end,jt_start:jt_end,k) 		&
    	+ (windl%w(it_start:it_end,jt_start:jt_end,k+1) - windl%w(it_start:it_end,jt_start:jt_end,k))*dz1w(k)
    enddo
!
!  BC
!
    if (bcl(6) == 'nnf') then
      div(it_start:it_end,jt_start:jt_end,nz) = div(it_start:it_end,jt_start:jt_end,nz) + 	&
      	(0. - windl%w(it_start:it_end,jt_start:jt_end,nz))*dz1w(nz)
    endif
!
  endif
!
! BC
!
  call statebc ( div ) 
!
! ---------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating caldiv')
!
return
end
!
!  ===========================================
   subroutine calconv ( windl, data, conv, horizontal  )
!  ===========================================
!
! ---------------------------------
! --- Calculate pressure gradients
! ---------------------------------
!
  integer :: i, j, k
  real  :: dx1, dy1
!
  logical, intent(in), optional :: horizontal
  type(atm_winds), intent(inout) :: windl
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: data
  real, dimension(ip_start:,jp_start:,:), intent(out) :: conv
!
  if (verbose > 2) call write_debug('Starting calconv')
!
! ---------------------------------------------------------!
!
  dx1 = 1./dx
#ifdef MODEL_3D
  dy1 = 1./dy 
#endif
!
  conv = 0.0
!
!----------------------------------------------------------!
!               Calculate data convergence                 !
!----------------------------------------------------------!
!
  do i = it_start,it_end
    conv(i,jt_start:jt_end,1:nz) = conv(i,jt_start:jt_end,1:nz)         &
            - ( 0.5*windl%u(i+1,jt_start:jt_end,1:nz)*(data(i+1,jt_start:jt_end,1:nz)+data(i,jt_start:jt_end,1:nz)) -      &
                0.5*windl%u(i,jt_start:jt_end,1:nz)*(data(i,jt_start:jt_end,1:nz)+data(i-1,jt_start:jt_end,1:nz)) )*dx1
  enddo
!
#ifdef MODEL_3D
  do j = jt_start,jt_end
    conv(it_start:it_end,j,1:nz) = conv(it_start:it_end,j,1:nz)         &
            - ( 0.5*windl%v(it_start:it_end,j+1,1:nz)*(data(it_start:it_end,j+1,1:nz)+data(it_start:it_end,j,1:nz)) -      &
                0.5*windl%v(it_start:it_end,j,1:nz)*(data(it_start:it_end,j,1:nz)+data(it_start:it_end,j-1,1:nz)) )*dy1
  enddo
#endif
!
  if ( .not.present(horizontal) ) then
    conv(it_start:it_end,jt_start:jt_end,1) = conv(it_start:it_end,jt_start:jt_end,1)         &
    	        - ( 0.5*windl%w(it_start:it_end,jt_start:jt_end,2)*(data(it_start:it_end,jt_start:jt_end,2)+data(it_start:it_end,jt_start:jt_end,1)) - 0.0 )*dz1w(1)
    do k = 2, nz-1
      conv(it_start:it_end,jt_start:jt_end,k) = conv(it_start:it_end,jt_start:jt_end,k)         &
    	        - ( 0.5*windl%w(it_start:it_end,jt_start:jt_end,k+1)*(data(it_start:it_end,jt_start:jt_end,k+1)+data(it_start:it_end,jt_start:jt_end,k)) -       &
                    0.5*windl%w(it_start:it_end,jt_start:jt_end,k)*(data(it_start:it_end,jt_start:jt_end,k)+data(it_start:it_end,jt_start:jt_end,k-1)) )*dz1w(k)
    enddo
  endif
!
! BC
!
  call statebc ( conv ) 
!
! ---------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating calconv')
!
return
end
!
!  ===========================================
   subroutine u_gradp ( windl, p, ugp )
!  ===========================================
!
! ---------------------------------
! --- Calculate u.grad(p)
! ---------------------------------
!
  integer :: i, j, k
  real  :: dx1, dy1
!
  type(atm_winds), intent(in) :: windl
  real, dimension(ip_start:,jp_start:,:), intent(in) :: p
  real, dimension(ip_start:,jp_start:,:), intent(out) :: ugp
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: gpx, gpy, gpz
!
! ---------------------------------------------------------!
!
  dx1 = 1./(2.*dx)
  dy1 = 1./(2.*dy)
!
  ugp = 0.0
  gpx = 0.0
  gpy = 0.0
  gpz = 0.0
!
!  Calculate cell centered pressure gradient
!
  do k = 1,nz
    do j = jt_start,jt_end
      do i = it_start,it_end
        gpx(i,j,k) = (p(i+1,j,k) - p(i-1,j,k))*dx1
      enddo
    enddo
  enddo
!
#ifdef MODEL_3D
  do k = 1,nz
    do j = jt_start,jt_end
      do i = it_start,it_end
        gpy(i,j,k) = (p(i,j+1,k) - p(i,j-1,k))*dy1	
      enddo
    enddo
  enddo
#endif
!
  do k = 2,nz-1
    do j = jt_start,jt_end
      do i = it_start,it_end
        gpz(i,j,k) = (p(i,j,k+1) - p(i,j,k-1)) / (dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
      enddo
    enddo
  enddo
!
  if (bcl(5) == 'ope') gpz(:,:,1) = 0.5*(p(:,:,1) + p(:,:,2)) * dz1w(1)
  if (bcl(6) == 'ope') gpz(:,:,nz) = -0.5*(p(:,:,nz) + p(:,:,nz-1)) * dz1w(nz)
  if (bcl(6) == 'nnf') gpz(:,:,nz) = -0.5*(p(:,:,nz) - p(:,:,nz-1)) * dz1(nz-1)
!
!  Calculate u.grad(p) 
!
  do j = jt_start,jt_end
    do i = it_start,it_end
      do k = 1,nz-1
        ugp(i,j,k) = 0.5*(windl%u(i+1,j,k) + windl%u(i,j,k)) * gpx(i,j,k)	      &
#ifdef MODEL_3D
		   + 0.5*(windl%v(i,j+1,k) + windl%v(i,j,k)) * gpy(i,j,k)	      &
#endif
		   + 0.5*(windl%w(i,j,k+1) + windl%w(i,j,k)) * gpz(i,j,k)
      enddo
!
      if (bcl(6) == 'ope') then
        ugp(i,j,nz) = 0.5*(windl%u(i+1,j,nz) + windl%u(i,j,nz))*gpx(i,j,nz) 	&
      		    + 0.5*(windl%v(i,j+1,k) + windl%v(i,j,k))*gpy(i,j,k) 	&
		    + 0.5*windl%w(i,j,nz)*gpz(i,j,nz)      
      endif
    enddo
  enddo
!
! ---------------------------------------------------------!
!
return
end subroutine u_gradp
!
! ======================================================
   subroutine gradp ( p, w, gpx, gpy, gpz )
! ======================================================
!
! ---------------------------------
! --- Calculate pressure gradients
! ---------------------------------
!
  integer :: i, j, k
  real  :: dx1, dy1
!
  real, dimension(ip_start:,jp_start:,:) :: p
  real, dimension(ip_start:,jp_start:,:), intent(in) :: w
  real, dimension(ip_start:,jp_start:,:), intent(out), optional :: gpx, gpy, gpz
!
  if (verbose > 2) call write_debug('Starting gradp')
!
! ---------------------------------------------------------!
!
!  Initializations
!
  dx1 = 1./dx
#ifdef MODEL_3D
  dy1 = 1./dy
#endif
!
  call statebc ( p )
!
!  Calculate gradient
!
  if ( present(gpx) ) then
!
  gpx(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
  do k = 1,nz
    do j = jt_start,jt_end
      do i = it_start,it_end
        gpx(i,j,k) = (p(i,j,k) - p(i-1,j,k))*dx1
      enddo
    enddo
  enddo
!
  endif
!
#ifdef MODEL_3D
  if ( present(gpy) ) then
!
  gpy(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
  do k = 1,nz
    do j = jt_start,jt_end
      do i = it_start,it_end
        gpy(i,j,k) = (p(i,j,k) - p(i,j-1,k))*dy1	
      enddo
    enddo
  enddo
!
#ifdef CHANNEL
  if (jp_start == 1) gpy(:,jt_start,:) = 0.
#endif
!
  endif
#endif
!
  if ( present(gpz) ) then
! 
  gpz(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
  do k = 2,nz
    do j = jt_start,jt_end
      do i = it_start,it_end
        gpz(i,j,k) = (p(i,j,k) - p(i,j,k-1))*dz1(k-1)
      enddo
    enddo
  enddo
!
  if (bcl(5) == 'nnf') then
#if (defined ANELASTIC)
    gpz(it_start:it_end,jt_start:jt_end,1) = w(it_start:it_end,jt_start:jt_end,1) / (dt0*avden0(1))
#else
    gpz(it_start:it_end,jt_start:jt_end,1) = 0. !w(it_start:it_end,jt_start:jt_end,1) / dt0
#endif
  endif
!
  endif
!
! ---------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating gradp')
!
return
end
  
  end module gradients
