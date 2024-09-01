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
!  ADVQ_WENO.F:
!	subroutine qnadv_weno(xmin,u,v,w,x,s)
!
!  Purpose:
!	Calculating advection of all scalars using a 3rd order WENO scheme.
!	Second order ENO reconstructions of the scalar's face values are
!	computed on 3 different stencils (from i-1 to i+2) and the total
!	flux is evaluated by weighting these 3 quantities depending on
!	the smoothness of the solution within the 3 chosen stencils. In 
!	the vertical direction, Lax-Friedrichs flux splitting is employed to
!	yield an upwind biased method.
!
!  Author:
!	Julien Savre
!       MISU
!
! ================================================================
!

  module advection_weno
!  
!  ---------------------------------------------
!
!
USE gridno
USE shared_data	
!
IMPLICIT NONE
  
  private
  
  public :: advsx_weno, advsy_weno, advsz_weno
  
  contains

  subroutine advsx_weno (i1,i2,x,u,xa,adv)
!
!----------------------------------------------------------!
!                                                          !
!  Scalar advection along x using WENO method		   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:ip_end) :: u
  real, dimension(ip_start:ip_end) :: x, xa, Flux
  real, dimension(ip_start:ip_end,1:4) :: ww
!
  logical  :: adv
  integer  :: i, i1, i2, ii1, ii2
  real     :: ffip1p, ffip1m, ffim1, ffip2
  real     :: dx1, uu, Flux_l, Flux_r, uold, Fold
!                                                          !
!----------------------------------------------------------!
!
  dx1 = 1./dx
  ii1 = max(i1-1,2)
  ii2 = min(i2,nx-2)
!
!  Calculate weights
!
  call weightx (ii1, ii2, x, ww)
!
!  Calculate fluxes at the faces
!
  do i=ii1,ii2
    uu = u(i+1)
    ffim1 = (-1./2.*x(i-1) + 3./2.*x(i))
    ffip1p = (1./2.*x(i) + 1./2.*x(i+1))
    ffip1m = (1./2.*x(i) + 1./2.*x(i+1))
    ffip2 = (3./2.*x(i+1) - 1./2.*x(i+2))
!
!    ffim1 = 1./3.*x(i-2,j,k) - 7./6.*x(i-1,j,k) + 11./6.*x(i,j,k)
!    ffip1 = -1./6.*x(i-1,j,k) + 5./6.*x(i,j,k) + 1./3.*x(i+1,j,k)
!    ffip2 = 1./3.*x(i,j,k) + 5./6.*x(i+1,j,k) - 1./6.*x(i+1,j,k)
!    ffip3 = 11./6.*x(i+1,j,k) - 7./6.*x(i+2,j,k) + 1./3.*x(i+3,j,k)
!
    call limiter_weno (ffim1, ffip1m, ffip1p, ffip2, x(i))
!
    Flux_l = ww(i,1)*ffim1 + ww(i,2)*ffip1p
    Flux_r = ww(i,3)*ffip1m + ww(i,4)*ffip2
!
    Flux(i) = 0.5*Flux_l*(uu+abs(uu)) + 0.5*Flux_r*(uu-abs(uu))
  enddo
!
!  Assemble tendencies
!
  Fold = Flux(ii1)
  do i=ii1+1,ii2
    xa(i) = xa(i) + (Flux(i) - Fold)*dx1
    Fold = Flux(i)
  enddo
!
  if (adv) then
    uold = u(ii2+1)
    do i=ii2,ii1+1,-1
      xa(i) = xa(i) - x(i)*(uold - u(i))*dx1
      uold = u(i)
    enddo
  endif     
!                                                          !
!----------------------------------------------------------!
!
  return
  end subroutine advsx_weno
!
!
  subroutine advsy_weno (j1,j2,x,v,xa,adv)
!
!----------------------------------------------------------!
!                                                          !
!  Scalar advection along y using WENO method              !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(jp_start:jp_end) :: v
  real, dimension(jp_start:jp_end) :: x, xa, Flux
  real, dimension(jp_start:jp_end,4) :: ww
!
  logical  :: adv
  integer  :: j, j1, j2, jj1, jj2
  real     :: ffim1, ffip1p, ffip1m, ffip2
  real     :: dy1, Flux_l, Flux_r, vv, vold, Fold
!                                                          !
!----------------------------------------------------------!
!
  dy1 = 1./dy 
  jj1 = max(j1-1,2)
  jj2 = min(j2,ny-2)
!
!  Calculate weights
!
  call weighty (jj1, jj2, x, ww)
!
!  Calculate fluxes at the faces
!
  do j=jj1,jj2  
    vv = v(j+1)
    ffim1 = (-1./2.*x(j-1) + 3./2.*x(j))
    ffip1p = (1./2.*x(j) + 1./2.*x(j+1))
    ffip1m = (1./2.*x(j) + 1./2.*x(j+1))
    ffip2 = (3./2.*x(j+1) - 1./2.*x(j+2))
!
    call limiter_weno (ffim1, ffip1m, ffip1p, ffip2, x(j))
!
    Flux_l = ww(j,1)*ffim1 + ww(j,2)*ffip1p
    Flux_r = ww(j,3)*ffip1m + ww(j,4)*ffip2
!
    Flux(j) = 0.5*Flux_l*(vv+abs(vv)) + 0.5*Flux_r*(vv-abs(vv))
  enddo      
!
!  Assemble tendencies
!
  Fold = Flux(jj1)
  do j=jj1+1,jj2
    xa(j) = xa(j) + (Flux(j) - Fold)*dy1
    Fold = Flux(j)
  enddo      
!
  if (adv) then
    vold = v(jj2+1)
    do j=jj2,jj1+1,-1
      xa(j) = xa(j) - x(j)*(vold - v(j))*dy1
      vold = v(j)
    enddo
  endif   
!                                                          !
!----------------------------------------------------------!
!
  return
  end subroutine advsy_weno
!
!
  subroutine advsz_weno (k1,k2,x,w,xa,adv)
!
!----------------------------------------------------------!
!                                                          !
!  Scalar advection in the vertical direction using WENO   !
!  method and Lax-Friedrichs flux splitting		   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(1:nz) :: w
  real, dimension(1:nz) :: x, xa, Flux
  real, dimension(1:nz,4) :: ww
!
  logical  :: adv
  integer  :: k, k1, k2 
  real     :: ffip1p, ffip1m, ffim1, ffip2
  real     :: Flux_l, Flux_r, wold, www, Fold
!                                                          !
!----------------------------------------------------------!
!
  Flux = 0.
!
!  Calculate weights
!
  call weightz (2, nz-2, x, ww)
!
!  Bottom boundary
!
  k = 1  
  www	= w(k+1)
  ffip1p = (1./2.*x(k) + 1./2.*x(k+1))
  ffip1m = (1./2.*x(k) + 1./2.*x(k+1))
  ffip2 = (3./2.*x(k+1) - 1./2.*x(k+2))
!
  Flux_l = ffip1p
  Flux_r = ww(k,3)*ffip1m + ww(k,4)*ffip2
!
  Flux(k) = 0.5*Flux_l*(www+abs(www)) + 0.5*Flux_r*(www-abs(www))
!
!  Calculate fluxes at the faces
!
  do k=k1,k2-2
    www   = w(k+1)
    ffim1 = (-1./2.*x(k-1) + 3./2.*x(k))
    ffip1p = (1./2.*x(k) + 1./2.*x(k+1))
    ffip1m = (1./2.*x(k) + 1./2.*x(k+1))
    ffip2 = (3./2.*x(k+1) - 1./2.*x(k+2))
!
    call limiter_weno (ffim1, ffip1m, ffip1p, ffip2, x(k))
!
    Flux_l = ww(k,1)*ffim1 + ww(k,2)*ffip1p
    Flux_r = ww(k,3)*ffip1m + ww(k,4)*ffip2
!
    Flux(k) = 0.5*Flux_l*(www+abs(www)) + 0.5*Flux_r*(www-abs(www))
  enddo
!
!  Top boundary
!
  k = nz-1  
  www	= w(k+1)
  ffim1 = (-1./2.*x(k-1) + 3./2.*x(k))
  ffip1p = (1./2.*x(k) + 1./2.*x(k+1))
  ffip1m = (1./2.*x(k) + 1./2.*x(k+1))
!
  Flux_l = ww(k,1)*ffim1 + ww(k,2)*ffip1p
  Flux_r = ffip1m
!
  Flux(k) = 0.5*Flux_l*(www+abs(www)) + 0.5*Flux_r*(www-abs(www))
!
!  Compute advection
!
  Fold = Flux(1)
  do k=2,k2-1
    xa(k) = xa(k) + (Flux(k) - Fold)*dz1w(k)
    Fold = Flux(k)
  enddo
!
  if (adv) then
    wold = w(k2)
    do k=k2-1,2,-1
      xa(k) = xa(k) - x(k)*(wold - w(k))*dz1w(k)
      wold = w(k)
    enddo
  endif
!
!  Boundary conditions
!
  if (bcl(5) == 'nnf') xa(1) = xa(1) + (Flux(1) - 0.)*dz1w(1)
  if (adv .and. bcl(5) == 'nnf') xa(1) = xa(1) - x(1)*(w(2) - 0.)*dz1w(1)
!
  if (bcl(6) == 'nnf') xa(nz) = xa(nz) + (0. - Flux(nz-1))*dz1w(nz)
  if (bcl(6) == 'nnf' .and. adv) xa(nz) = xa(nz) - x(nz)*(0. - w(nz))*dz1w(nz)
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsz_weno  
!
!
  subroutine weightx (i1, i2, x, ww)
!
!----------------------------------------------------------!
!                                                          !
!  Calculate WENO weights for x advection		   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:ip_end) :: x
  real, dimension(ip_start:ip_end,4) :: ww  
!
  integer ::  i, i1, i2
  real    ::  eta(4)
  real	  ::  sumw, w1, w2, w3, w4
  real    ::  g1=1./3., g2=2./3., eps=1.e-15, exp=2., thres=20.
!                                                          !
!----------------------------------------------------------!
!
  ww(:,1) = g1
  ww(:,2) = g2
  ww(:,3) = g2
  ww(:,4) = g1 
!
  do i=i1,i2  
    eta(1) = (x(i) - x(i-1))**2.  
    eta(2) = (x(i+1) - x(i))**2.
    eta(3) = (x(i+2) - x(i+1))**2.	
!
    w1 = g1 / (eps + eta(1))**exp
    w2 = g2 / (eps + eta(2))**exp
    w3 = g2 / (eps + eta(2))**exp
    w4 = g1 / (eps + eta(3))**exp
!
    if (max(eta(1),eta(2)) / (min(eta(1),eta(2))+eps) > thres) then
      sumw = w1 + w2
      ww(i,1) = w1/sumw
      ww(i,2) = w2/sumw
    endif
!
    if (max(eta(2),eta(3)) / (min(eta(2),eta(3))+eps) > thres) then
      sumw = w3 + w4
      ww(i,3) = w3/sumw
      ww(i,4) = w4/sumw
    endif
  enddo
!
  return
  end subroutine weightx
!
!
  subroutine weighty (j1, j2, x, ww)
!
!----------------------------------------------------------!
!                                                          !
!  Calculate WENO weights for y advection		   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(jp_start:jp_end) :: x
  real, dimension(jp_start:jp_end,4) :: ww  
!
  integer  ::  j, j1, j2
  real, dimension(4) ::  eta
  real  ::  sumw, w1, w2, w3, w4 
  real  ::  g1=1./3., g2=2./3., eps=1.e-15, exp=2., thres=20.
!                                                          !
!----------------------------------------------------------!
!
  ww(:,1) = g1
  ww(:,2) = g2
  ww(:,3) = g2
  ww(:,4) = g1
!
  do j=j1,j2
    eta(1) = (x(j) - x(j-1))**2.  
    eta(2) = (x(j+1) - x(j))**2.
    eta(3) = (x(j+2) - x(j+1))**2.
!
    w1 = g1 / (eps + eta(1))**exp
    w2 = g2 / (eps + eta(2))**exp
    w3 = g2 / (eps + eta(2))**exp
    w4 = g1 / (eps + eta(3))**exp
!
    if (max(eta(1),eta(2)) / (min(eta(1),eta(2))+eps) > thres) then
      sumw = w1 + w2
      ww(j,1) = w1/sumw
      ww(j,2) = w2/sumw
    endif
! 
    if (max(eta(2),eta(3)) / (min(eta(2),eta(3))+eps) > thres) then
      sumw = w3 + w4	    
      ww(j,3) = w3/sumw
      ww(j,4) = w4/sumw
    endif
  enddo
!
  return
  end subroutine weighty
!
!
  subroutine weightz (k1, k2, x, ww)
!
!----------------------------------------------------------!
!                                                          !
!  Calculate WENO weights for z advection		   !
!  The eta exponents to calculate weights are set to 3    !
!  instead of 2 to increase solution's smoothness.         !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(1:nz) :: x
  real, dimension(1:nz,4) :: ww  
!
  integer  :: k, k1, k2
  real, dimension(4) ::  eta
  real  ::  sumw, w1, w2, w3, w4
  real  ::  g1=1./3., g2=2./3., eps=1.e-15, exp=2., thres=20.
!                                                          !
!----------------------------------------------------------!
!
  ww(:,1) = g1
  ww(:,2) = g2
  ww(:,3) = g2
  ww(:,4) = g1 
!
  do k=k1+1,k2-1
    eta(1) = (x(k) - x(k-1))**2.  
    eta(2) = (x(k+1) - x(k))**2.
    eta(3) = (x(k+2) - x(k+1))**2.
!
    w1 = g1 / (eps + eta(1))**exp
    w2 = g2 / (eps + eta(2))**exp
    w3 = g2 / (eps + eta(2))**exp
    w4 = g1 / (eps + eta(3))**exp
!
    if (max(eta(1),eta(2)) / (min(eta(1),eta(2))+eps) > thres) then
      sumw = w1 + w2
      ww(k,1) = w1/sumw
      ww(k,2) = w2/sumw
    endif
!
    if (max(eta(2),eta(3)) / (min(eta(2),eta(3))+eps) > thres) then
      sumw = w3 + w4	    
      ww(k,3) = w3/sumw
      ww(k,4) = w4/sumw
    endif
  enddo
!
  return
  end subroutine weightz
!
!
  subroutine limiter_weno (ffim1, ffip1m, ffip1p, ffip2, x)
!
!----------------------------------------------------------!
!                                                          !
!  Simple limiter (rescaling) to preserve positivity       !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real  :: ffim1, ffip1m, ffip1p, ffip2, x
  real  :: maxp, minp, phi, eps=1.e-15
!                                                          !
!----------------------------------------------------------!
!
  maxp = max(ffim1, ffip1m, ffip1p, ffip2)
  minp = min(ffim1, ffip1m, ffip1p, ffip2)
!
  phi = min(abs((0. - x)/(minp - x + eps)), 1.)
!
  ffim1  = phi*(ffim1 - x) + x
  ffip1m = phi*(ffip1m - x) + x
  ffip1p = phi*(ffip1p - x) + x
  ffip2  = phi*(ffip2 - x) + x
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine limiter_weno
  
end module advection_weno
