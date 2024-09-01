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
!  ADVQ_MUSCL.F:
!
!  Purpose:
!	Calculating advection of all scalars using a MUSCL type
!   Finite Volume scheme, using linear interpolations.
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================
!
  module advection_muscl
!  
!  ---------------------------------------------
!
!
USE gridno
USE shared_data	
!
IMPLICIT NONE
  
  private
  
  character(len=10)  :: muscllimiter='minmod'
  real :: eps = 1.e-18
  
  public :: muscllimiter, advsx_muscl, advsy_muscl, advsz_muscl
  
  contains

  subroutine advsx_muscl (i1,i2,x,u,xa,adv)
!
!----------------------------------------------------------!
!                                                          !
!  Scalar advection along x using MUSCL scheme		   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension (ip_start:ip_end) :: u
  real, dimension(ip_start:ip_end) :: x, xa, Flux
!
  logical  :: adv
  integer  :: i, i1, i2
  real     :: slope_p, slope_m, slope_p2, phil, phir, xx
  real     :: dx1, uu, uold, Fold, Flux_r, Flux_l
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
    xx = x(i)
!
    slope_p  = (x(i+1) - xx) * dx1
    slope_p2 = (x(i+2) - x(i+1)) * dx1
    slope_m  = (xx - x(i-1)) * dx1
!
    if ( dx*abs(slope_m) > lim_tol*abs(xx)+eps .or. dx*abs(slope_p2) > lim_tol*abs(xx)+eps ) then
      call limiter_muscl (xx, slope_p, slope_p2, slope_m, phir, phil)
!
      Flux_r = x(i+1) - 0.5*phir*slope_p2*dx
      Flux_l = xx + 0.5*phil*slope_m*dx
    else
      Flux_r = x(i+1) - 0.5*slope_p*dx
      Flux_l = xx + 0.5*slope_p*dx      
    endif
!
    Flux(i) = 0.5*uu*( Flux_l*(1.+sign(1.,uu)) + Flux_r*(1.-sign(1.,uu)) )
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
    uold = u(i2+1)
    do i=i2,i1,-1
      xa(i) = xa(i) - x(i)*(uold - u(i))*dx1
      uold = u(i)
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsx_muscl
!
!
  subroutine advsy_muscl (j1,j2,x,v,xa,adv)
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
  real     :: slope_p, slope_m, slope_p2, phil, phir, xx
  real     :: dy1, Flux_l, Flux_r, vv, vold, Fold
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
    xx = x(j)
!
    slope_p  = (x(j+1) - xx) * dy1
    slope_p2 = (x(j+2) - x(j+1)) * dy1
    slope_m  = (xx - x(j-1)) * dy1
!
    if ( dy*abs(slope_m) > lim_tol*abs(xx)+eps .or. dy*abs(slope_p2) > lim_tol*abs(xx)+eps ) then
      call limiter_muscl (xx, slope_p, slope_p2, slope_m, phir, phil)
!
      Flux_r = x(j+1) - 0.5*phir*slope_p2*dy
      Flux_l = xx + 0.5*phil*slope_m*dy
    else
      Flux_r = x(j+1) - 0.5*slope_p*dy
      Flux_l = xx + 0.5*slope_p*dy     
    endif
!
    Flux(j) = 0.5*vv*( Flux_l*(1.+sign(1.,vv)) + Flux_r*(1.-sign(1.,vv)) )
  enddo
!
#ifdef CHANNEL
  if (jp_start == 1) then
    Flux(j1-1) = 0.
    Flux(j1) = 0.5*v(j1+1)*( (1.+sign(1.,v(j1+1)))*(0.5*(x(j1+1)+x(j1)) - 0.25*str*v(j1+1)*(1.+sign(1.,v(j1+1)))*(x(j1+1)-x(j1))*dt0*dy1)  & 
	     - (1.-sign(1.,v(j1+1)))*(0.5*(x(j1+1)+x(j1)) + 0.25*str*v(j1+1)*(1.+sign(1.,v(j1+1)))*(x(j1+1)-x(j1))*dt0*dy1) )
  endif
!
  if (jp_end == ny) then
    Flux(j2) = 0.
    Flux(j2-1) = 0.5*v(j2)*( (1.+sign(1.,v(j2)))*(0.5*(x(j2)+x(j2-1)) - 0.25*str*v(j2)*(1.+sign(1.,v(j2)))*(x(j2)-x(j2-1))*dt0*dy1)  & 
	     - (1.-sign(1.,v(j2)))*(0.5*(x(j2)+x(j2-1)) + 0.25*v(j2)*str*(1.+sign(1.,v(j2)))*(x(j2)-x(j2-1))*dt0*dy1) )
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
    vold = v(j2+1)
    do j=j2,j1,-1
      xa(j) = xa(j) - x(j)*(vold - v(j))*dy1
      vold = v(j)
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsy_muscl
!
!
  subroutine advsz_muscl (k1,k2,x,w,xa,adv)
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
  real, dimension(1:nz) :: x, xa, Flux
!
  logical  :: adv
  integer  :: k, k1, k2
  real     :: slope_p, slope_p2, slope_m, phil, phir
  real     :: ww, wold, Flux_l, Flux_r, Fold
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
!  Calculate fluxes
!
  k = 1
  ww = w(2)
  slope_p  = (x(2) - x(1)) * dz1(1)
  slope_p2 = (x(3) - x(2)) * dz1(2)
  slope_m = bc51*(x(nz-4) - x(nz-5))*dz1(nz-5) + (1.-bc51)*slope_p
!
  if ( fdz0(2)*dz*abs(slope_p2) > lim_tol*abs(x(1))+eps ) then
    call limiter_muscl (x(1), slope_p, slope_p2, slope_m, phir, phil)
!
    Flux_l = x(1) + 0.5*phil*slope_m/dz1w(1)
    Flux_r = x(2) - 0.5*phir*slope_p2/dz1w(1)
  else
    Flux_l = x(1) + 0.5*slope_p/dz1w(1)
    Flux_r = x(2) - 0.5*slope_p/dz1w(1)
  endif
!
  Flux(1) = 0.5*ww*( Flux_l*(1.+sign(1.,ww)) + Flux_r*(1.-sign(1.,ww)) )
!
!  Main core
!
  do k=2,nz-2
    ww = w(k+1)
    slope_p  = (x(k+1) - x(k)) * dz1(k)
    slope_p2 = (x(k+2) - x(k+1)) * dz1(k+1)
    slope_m  = (x(k) - x(k-1)) * dz1(k-1)
!
    if ( fdz0(k-1)*dz*abs(slope_m) > lim_tol*abs(x(k))+eps .or. fdz0(k+1)*dz*abs(slope_p2) > lim_tol*abs(x(k))+eps ) then
      call limiter_muscl (x(k), slope_p, slope_p2, slope_m, phir, phil)
!
      Flux_r = x(k+1) - 0.5*phir*slope_p2/dz1w(k)
      Flux_l = x(k) + 0.5*phil*slope_m/dz1w(k)
    else
      Flux_r = x(k+1) - 0.5*slope_p/dz1w(k)
      Flux_l = x(k) + 0.5*slope_p/dz1w(k)
    endif
!
    Flux(k) = 0.5*ww*( Flux_l*(1.+sign(1.,ww)) + Flux_r*(1.-sign(1.,ww)) )
  enddo
!
!  Top boundary
!
  k = nz-1
  ww = w(k+1)
  slope_p  = (x(k+1) - x(k)) * dz1(k)
  slope_p2 = bc61*(x(6) - x(5))*dz1(5) + (1.-bc61)*slope_p
  slope_m  = (x(k) - x(k-1)) * dz1(k-1)
!
  if ( fdz0(k-1)*dz*abs(slope_m) > lim_tol*abs(x(nz-1))+eps ) then
    call limiter_muscl (x(k), slope_p, slope_p2, slope_m, phir, phil)
!
    Flux_r = x(k+1) - 0.5*phir*slope_p2/dz1w(k)
    Flux_l = x(k) + 0.5*phil*slope_m/dz1w(k)
  else
    Flux_r = x(k+1) - 0.5*slope_p/dz1w(k)
    Flux_l = x(k) + 0.5*slope_p/dz1w(k)
  endif
!
  Flux(k) = 0.5*ww*( Flux_l*(1.+sign(1.,ww)) + Flux_r*(1.-sign(1.,ww)) )
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
    wold = w(k2)
    do k=nz-1,2,-1
      xa(k) = xa(k) - x(k)*(wold - w(k))*dz1w(k)
      wold = w(k)
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
  end subroutine advsz_muscl
!
!
  subroutine limiter_muscl (xx, ffip1, ffip2, ffim1, phir, phil)
!
!----------------------------------------------------------!
!                                                          !
!  FLux limiter_muscls using prescribed functions.   	   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real  :: ffip1, ffim1, ffip2, xx
  real  :: phir, phil
  real  :: rp1, rp2, rp11, rp21
  real  :: kappa
!
  character(len=10) :: limiter
!                                                          !
!----------------------------------------------------------!
!
  limiter = muscllimiter
!
  rp1 = 1.
  rp2 = 1.
  rp11 = 1.
  rp21 = 1.
!
  if ( abs(ffim1) > lim_tol*abs(xx) .and. abs(ffip1) < abs(ffim1)/lim_tol ) rp1 = ffip1 / ffim1
  if ( abs(ffip2) > lim_tol*abs(xx) .and. abs(ffip1) < abs(ffip2)/lim_tol ) rp2 = ffip1 / ffip2
!
  phil = rp1
  phir = rp2
!
  select case (trim(limiter))
    case ('VL')		  
      phil = (rp1 + abs(rp1)) / (1. + abs(rp1))
      phir = (rp2 + abs(rp2)) / (1. + abs(rp2))
    case ('MC')	
      phil = max(0.,min(2.*rp1,0.5*(rp1+1.),2.))
      phir = max(0.,min(2.*rp2,0.5*(rp2+1.),2.))
    case ('minmod')
      phil = max(0.,min(1.,rp1))
      phir = max(0.,min(1.,rp2))
    case default
      phil = 0.5*(rp1+1.)
      phir = 0.5*(rp2+1.)
  end select
!
  return
  end
  
end module advection_muscl
