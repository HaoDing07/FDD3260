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
  module advection_lw
!  
!  ---------------------------------------------
!
!
USE gridno
USE shared_data	
!
IMPLICIT NONE
  
  private
  
  character(len=10), save :: lwlimiter='MC'
  real :: eps = 1.e-18
  
  public :: lwlimiter, advsx_lw, advsy_lw, advsz_lw
  
  contains

  subroutine advsx_lw (i1,i2,x,u,xa,adv)
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
  real     :: ffip1, ffim1, ffip2, phil, phir, fr, fl
  real     :: dx1, uu, xx, Fold, xp1, xp2, mu
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
    ffip1 = x(i+1) - xx
    ffim1 = xx - x(i-1)
    ffip2 = x(i+2) - x(i+1)
    mu = abs(uu)*dt0*dx1
!
    if ( abs(ffim1) > lim_tol*abs(xx)+eps .or. abs(ffip2) > lim_tol*abs(xx)+eps ) then
      call limiter_lw (xx, ffip1, ffip2, ffim1, phir, phil)
!
      fl = xx + 0.5*phil*ffim1*(1. - mu)
      fr = x(i+1) - 0.5*phir*ffip2*(1. - mu)
    else
      fl = xx + 0.5*ffip1*(1. - mu)
      fr = x(i+1) - 0.5*ffip1*(1. - mu)
    endif
!
    Flux(i) = 0.5*uu*( fl*(sign(1.,uu)+1.) - fr*(sign(1.,uu)-1.) )
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
  end subroutine advsx_lw
!
!
  subroutine advsy_lw (j1,j2,x,v,xa,adv)
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
  real     :: ffip1, ffim1, ffip2, phil, phir, fr, fl
  real     :: dy1, vv, xx, Fold, xp1, xp2, mu
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
    ffip1 = x(j+1) - xx
    ffim1 = xx - x(j-1)
    ffip2 = x(j+2) - x(j+1)
    mu = abs(vv)*dt0*dy1
!
    if ( abs(ffim1) > lim_tol*abs(xx)+eps .or. abs(ffip2) > lim_tol*abs(xx)+eps ) then
      call limiter_lw (xx, ffip1, ffip2, ffim1, phir, phil)
!
      fl = xx + 0.5*phil*ffim1*(1. - mu)
      fr = x(j+1) - 0.5*phir*ffip2*(1. - mu)
    else
      fl = xx + 0.5*ffip1*(1. - mu)
      fr = x(j+1) - 0.5*ffip1*(1. - mu)
    endif
!
    Flux(j) = 0.5*vv*( fl*(sign(1.,vv)+1.) - fr*(sign(1.,vv)-1.) )
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
  end subroutine advsy_lw
!
!
  subroutine advsz_lw (k1,k2,x,w,xa,adv)
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
  real     :: ffip1, ffim1, ffip2, phil, phir, fr, fl
  real     :: ww, xx, Fold, xp2, xp1, mu
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
!  Handle level 1
!
  k = 1
  ww = w(2)
  mu = abs(ww)*dt0*dz1(1)
  ffip1 = (x(2) - x(1)) * dz1(1)
  ffip2 = (x(3) - x(2)) * dz1(2)
  ffim1 = bc51*(x(nz-4) - x(nz-5))*dz1(nz-5) + (1.-bc51)*ffip1
!
  if ( dz*abs(ffim1) > lim_tol*abs(x(1))+eps .or. dz*abs(ffip2) > lim_tol*abs(x(1))+eps ) then
    call limiter_lw (x(1), ffip1, ffip2, ffim1, phir, phil)
!
    fl = x(1) + 0.5*phil*ffim1*(1. - mu)/dz1w(1)
    fr = x(2) - 0.5*phir*ffip2*(1. - mu)/dz1w(1)
  else
    fl = x(1) + 0.5*ffip1*(1. - mu)/dz1w(1)
    fr = x(2) - 0.5*ffip1*(1. - mu)/dz1w(1)
  endif
!
  Flux(1) = 0.5*ww*( fl*(sign(1.,ww)+1.) - fr*(sign(1.,ww)-1.) )
!
!  Main core
!
  do k=2,nz-2
    xx = x(k)
    ww = w(k+1)
    xp1 = x(k+1)
    xp2 = x(k+2)
!
    ffim1 = (xx - x(k-1)) * dz1(k-1)
    ffip1 = (xp1 - xx) * dz1(k)
    ffip2 = (xp2 - xp1) * dz1(k+1)
    mu = abs(ww)*dt0*dz1(k)
!
    if ( dz*abs(ffim1) > lim_tol*abs(xx)+eps .or. dz*abs(ffip2) > lim_tol*abs(xx)+eps ) then
      call limiter_lw (xx, ffip1, ffip2, ffim1, phir, phil)
! 
      fl = xx + 0.5*phil*ffim1*(1. - mu)/dz1w(k)
      fr = xp1 - 0.5*phir*ffip2*(1. - mu)/dz1w(k)
    else
      fl = xx + 0.5*ffip1*(1. - mu)/dz1w(k)
      fr = xp1 - 0.5*ffip1*(1. - mu)/dz1w(k)
    endif
!
    Flux(k) = 0.5*ww*( fl*(sign(1.,ww)+1.) - fr*(sign(1.,ww)-1.) )
  enddo
!
!  Handle level nz-1
!
  k = nz-1
  ww = w(nz)
  xx = x(k)
  ffim1 = (xx - x(k-1)) * dz1(k-1) 
  ffip1 = (x(k+1) - xx) * dz1(k)
  ffip2 = bc61*(x(6) - x(5))*dz1(5) + (1.-bc61)*ffip1
  mu = abs(ww)*dt0*dz1(k)
!
  if ( dz*abs(ffim1) > lim_tol*abs(xx)+eps .or. dz*abs(ffip2) > lim_tol*abs(xx)+eps ) then
    call limiter_lw (xx, ffip1, ffip2, ffim1, phir, phil)
!
    fl = xx + 0.5*phil*ffim1*(1. - mu)/dz1w(k)
    fr = x(k+1) - 0.5*phir*ffip2*(1. - mu)/dz1w(k)
  else
    fl = xx + 0.5*ffip1*(1. - mu)/dz1w(k)
    fr = x(k+1) - 0.5*ffip1*(1. - mu)/dz1w(k)
  endif
!
  Flux(k) = 0.5*ww*( fl*(sign(1.,ww)+1.) - fr*(sign(1.,ww)-1.) )
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
  end subroutine advsz_lw
!
  subroutine limiter_lw (x, ffip1, ffip2, ffim1, phir, phil)
!
!----------------------------------------------------------!
!                                                          !
!  	FLux limiter_lws using prescribed functions.      !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real  :: x, ffip1, ffim1, ffip2
  real  :: phir, phil
  real  :: Kpr, Kpl
  real  :: rp1, rp2, rp11, rp21
  real  :: kappa
!
  character(len=10) :: limiter
!                                                          !
!----------------------------------------------------------!
!
  limiter = lwlimiter
!
  rp1 = 1.
  rp2 = 1.
  rp11 = 1.
  rp21 = 1.
!
  if ( abs(ffim1) > lim_tol*abs(x)+eps ) rp1 = ffip1 / ffim1
  if ( abs(ffip2) > lim_tol*abs(x)+eps ) rp2 = ffip1 / ffip2
!
  phil = rp1
  phir = rp2
!
  select case (limiter)
    case ('koren3rd')    !Equivalent to QUICK 
      kappa = 1./3.
      Kpl = 0.5*(1.-kappa) + 0.5*(1.+kappa)*rp1
      Kpr = 0.5*(1.-kappa) + 0.5*(1.+kappa)*rp2
      phil = max(min(2.*rp1,Kpl,2.),0.)
      phir = max(min(2.*rp2,Kpr,2.),0.)
    case ('VL')		  
      phil = (rp1 + abs(rp1)) / (1. + abs(rp1))
      phir = (rp2 + abs(rp2)) / (1. + abs(rp2))
    case ('MC')	
      phil = max(min(2.*rp1,0.5*(rp1+1.),2.),0.)
      phir = max(min(2.*rp2,0.5*(rp2+1.),2.),0.)
    case ('minmod') 
      phil = max(min(1.,rp1),0.)
      phir = max(min(1.,rp2),0.)
!
!  Linear schemes (unlimited)
!
    case ('upwind')
      phil = 0.
      phir = 0.
    case ('LW')
      phil = rp1
      phir = rp2
    case ('WB')
      phil = 1.
      phir = 1.
    case ('Fromm')
      phil = 0.5*(1. + rp1)
      phir = 0.5*(1. + rp2)
    case ('3rd')
      kappa = 1./3.
      phil = 0.5*(1.-kappa) + 0.5*(1.+kappa)*rp1
      phir = 0.5*(1.-kappa) + 0.5*(1.+kappa)*rp2
    case ('4th')
      kappa = 1./6.  
      if ( abs(ffip1) > lim_tol*abs(x) ) rp11 = ffip2 / ffip1
      if ( abs(ffip1) > lim_tol*abs(x) ) rp21 = ffim1 / ffip1
      phil = rp1 + kappa - kappa*rp1*rp11
      phir = rp2 + kappa - kappa*rp2*rp21
  end select
!
  return
  end
  
end module advection_lw
