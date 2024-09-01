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
!  ADVQ_lw.F90:
!
!  Purpose:
!	Calculating advection of all scalars using an upwind biased 3rd order 
!	TVD method based on finite differences. The method follows Hundsdorfer et al.
!       from  JCP 1995.
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================
!
  module advection_lw_nl
!  
!  ---------------------------------------------
!
!
USE gridno
USE shared_data	
!
IMPLICIT NONE
  
  private
  
  public :: advsx_lw_nl, advsy_lw_nl, advsz_lw_nl
  
  contains

  subroutine advsx_lw_nl (i1,i2,x,u,xa,adv)
!
!----------------------------------------------------------!
!                                                          !
!  Scalar advection along x using MUSCL scheme		   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:ip_end) :: u
  real, dimension(ip_start:ip_end) :: x, xa, Flux
!
  logical  :: adv
  integer  :: i, i1, i2
  real     :: ffip1, fr, fl, ffip2, ffim1
  real     :: dx1, uu, xx, mu, Fold
!                                                          !
!----------------------------------------------------------!
!
  dx1 = 1./dx
  Flux = 0.0
!
!  Calculate fluxes
!
  do i=i1-1,i2
    uu = u(i+1)
    mu = abs(uu)*dt0*dx1
    xx = x(i)
    ffip1  = x(i+1) - xx
!
    fl = xx + 0.5*ffip1*(1. - mu)
    fr = x(i+1) - 0.5*ffip1*(1. - mu)
!
    Flux(i) = 0.5*fl*(uu+abs(uu)) + 0.5*fr*(uu-abs(uu))
  enddo
!
!  Compute advection
!
  Fold = Flux(i1-1)
  do i=i1,i2
    xa(i) = xa(i) + (Flux(i) - Fold)*dx1
    Fold = Flux(i)
  enddo
!
  if (adv) then
    uu = u(i2+1)
    do i=i2,i1,-1
      xa(i) = xa(i) - x(i)*(uu - u(i))*dx1
      uu = u(i)
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsx_lw_nl
!
!
  subroutine advsy_lw_nl (j1,j2,x,v,xa,adv)
!
!----------------------------------------------------------!
!                                                          !
!  Scalar advection using MUSCL scheme			   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(jp_start:jp_end) :: v
  real, dimension(jp_start:jp_end) :: x, xa, Flux
!
  logical  :: adv
  integer  :: j, j1, j2
  real     :: ffip1, fr, fl, ffim1, ffip2
  real     :: dy1, vv, xx, mu, Fold
!                                                          !
!----------------------------------------------------------!
!
  dy1 = 1./dy
  Flux = 0.0
!
!  Calculate fluxes
!
  do j=j1-1,j2
    vv = v(j+1)
    mu = abs(vv)*dt0*dy1
    xx = x(j)
    ffip1 = x(j+1) - xx
!
    fr = x(j+1) - 0.5*ffip1*(1. - mu)
    fl = xx + 0.5*ffip1*(1. - mu) 
!
    Flux(j) = 0.5*fl*(vv+abs(vv)) + 0.5*fr*(vv-abs(vv))
  enddo
!
#ifdef CHANNEL
  if (jp_start == 1) then
    Flux(j1-1) = 0.
    Flux(j1) = 0.5*v(j1+1)*( (1.+sign(1.,v(j1+1)))*(0.5*(x(j1+1)+x(j1)) - 0.25*v(j1+1)*(1.+sign(1.,v(j1+1)))*(x(j1+1)-x(j1))*dt0*dy1)  & 
	     - (1.-sign(1.,v(j1+1)))*(0.5*(x(j1+1)+x(j1)) + 0.25*v(j1+1)*(1.+sign(1.,v(j1+1)))*(x(j1+1)-x(j1))*dt0*dy1) )
  endif
!
  if (jp_end == ny) then
    Flux(j2) = 0.
    Flux(j2-1) = 0.5*v(j2)*( (1.+sign(1.,v(j2)))*(0.5*(x(j2)+x(j2-1)) - 0.25*v(j2)*(1.+sign(1.,v(j2)))*(x(j2)-x(j2-1))*dt0*dy1)  & 
	       - (1.-sign(1.,v(j2)))*(0.5*(x(j2)+x(j2-1)) + 0.25*v(j2)*(1.+sign(1.,v(j2)))*(x(j2)-x(j2-1))*dt0*dy1) )
  endif
#endif
!
!  Compute advection
!
  Fold = Flux(j1-1)
  do j=j1,j2
    xa(j) = xa(j) + (Flux(j) - Fold)*dy1 
    Fold = Flux(j)
  enddo
!
  if (adv) then
    vv = v(j2+1)
    do j=j2,j1,-1
      xa(j) = xa(j) - x(j)*(vv - v(j))*dy1
      vv = v(j)
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsy_lw_nl
!
!
  subroutine advsz_lw_nl (k1,k2,x,w,xa,adv)
!
!----------------------------------------------------------!
!                                                          !
!  Scalar advection using MUSCL scheme			   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(1:nz) :: w
  real, dimension(1:nz) :: x, xa , Flux
!
  logical  :: adv
  integer  :: k, k1, k2
  real     :: deltac, deltad, ffip1, fr, fl, ffip2, ffim1
  real     :: mu, ww, xx, Fold
  real	   :: bc51, bc61
!                                                          !
!----------------------------------------------------------!
!
  Flux = 0.
!
!  Boundary flags
!
  bc51=0.; bc61=0.
!
  if (bcl(5) == 'per') then
    bc51 = 1.
  endif
!
  if (bcl(6) == 'per') then
    bc61 = 1.
  endif
!
!  Main core
!
  do k=1,nz-1
    ww = w(k+1)
    mu = abs(ww)*dt0*dz1(k)
    xx = x(k)
    ffip1  = (x(k+1) - xx) * dz1(k)
!
    fl = xx + 0.5*ffip1*(1. - mu)/dz1w(k) 
    fr = x(k+1) - 0.5*ffip1*(1. - mu)/dz1w(k)
!
    Flux(k) = 0.5*fl*(ww+abs(ww)) + 0.5*fr*(ww-abs(ww))
  enddo
!
!  Compute advection
!
  Fold = Flux(1)
  do k=2,nz-1
    xa(k) = xa(k) + (Flux(k) - Fold)*dz1w(k) 
    Fold = Flux(k)
  enddo
!
  if (adv) then
    ww = w(k2)
    do k=nz-1,2,-1
      xa(k) = xa(k) - x(k)*(ww - w(k))*dz1w(k)
      ww = w(k)
    enddo
  endif
!
!  Boundary conditions
!
  if (bcl(5) == 'nnf') xa(1) = xa(1) + (Flux(1) - 0.)*dz1w(1)
  if (bcl(5) == 'nnf' .and. adv) xa(1) = xa(1) - x(1)*(w(2) - 0.)*dz1w(1)
!
  if (bcl(6) == 'nnf') xa(nz) = xa(nz) + (0. - Flux(nz-1))*dz1w(nz)
  if (bcl(6) == 'nnf' .and. adv) xa(nz) = xa(nz) - x(nz)*(0. - w(nz))*dz1w(nz)
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsz_lw_nl
!
end module advection_lw_nl
