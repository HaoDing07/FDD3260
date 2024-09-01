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
!  ADVQ_QUICK_NOLIM.F90:
!
!  Purpose:
!	Calculating advection of all scalars using a QUICK type
!       Finite Volume scheme, using quadratic polynomial interpolations.
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================
!
  module advection_quick_nl
!  
!  ---------------------------------------------
!
!
USE gridno
USE shared_data	
!
IMPLICIT NONE
  
  private
  
  public :: advsx_quick_nl, advsy_quick_nl, advsz_quick_nl
  
  contains

  subroutine advsx_quick_nl (i1,i2,x,u,xa,adv)
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
  real     :: ffip1, ffip2, ffim1
  real     :: Flux_r, Flux_l, phil, phir, xx
  real     :: cc, dx1, uu, uold, Fold, rp1, rp2
!                                                          !
!----------------------------------------------------------!
!
  cc = 1./4.
  dx1 = 1./dx
  Flux = 0.0
!
!  Calculate fluxes
!
  do i=i1-1,i2
    uu = u(i+1)
!
    xx = x(i)
    ffip1  = (x(i+1) - xx)
    ffip2  = (x(i+2) - x(i+1))
    ffim1  = (xx - x(i-1))
!
    phil = ffip1 - cc*(ffip1 - ffim1)
    phir = ffip1 + cc*(ffip2 - ffip1)
!
    Flux_l = xx + 0.5*phil
    Flux_r = x(i+1) - 0.5*phir
!
    Flux(i) = 0.5*Flux_l*(uu+abs(uu)) + 0.5*Flux_r*(uu-abs(uu))
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
  end subroutine advsx_quick_nl
!
!
  subroutine advsy_quick_nl (j1,j2,x,v,xa,adv)
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
  real     :: ffip1, ffim1, ffip2, phil, phir, xx 
  real     :: cc, dy1, Flux_l, Flux_r, vv, vold, Fold
! 
!----------------------------------------------------------!
!
  cc = 1./4.
  dy1 = 1./dy
  Flux = 0.0
!
!  Calculate fluxes
!
  do j=j1-1,j2
    vv = v(j+1)
    xx = x(j)
    ffip1  = (x(j+1) - xx)
    ffip2  = (x(j+2) - x(j+1))
    ffim1  = (xx - x(j-1))
!
    phil = ffip1 - cc*(ffip1 - ffim1)
    phir = ffip1 + cc*(ffip2 - ffip1)
!
    Flux_l = xx + 0.5*phil
    Flux_r = x(j+1) - 0.5*phir
!
    Flux(j) = 0.5*Flux_l*(vv+abs(vv)) + 0.5*Flux_r*(vv-abs(vv))
  enddo
!
#ifdef CHANNEL
  if (jp_start == 1) then
    Flux(j1-1) = 0.
    Flux_l = x(j1) + 0.5*(x(j1+1) - x(j1))
    Flux_r = x(j1+1) - 0.5*(x(j1+1) - x(j1)) - 1./6.*(x(j1+2) - 2.*x(j1+1) + x(j1))
    Flux(j1) = 0.5*Flux_l*(v(j1+1)+abs(v(j1+1))) + 0.5*Flux_r*(v(j1+1)-abs(v(j1+1)))
  endif
!
  if (jp_end == ny) then
    Flux(j2) = 0.
    Flux_l = x(j2-1) + 0.5*(x(j2) - x(j2-1)) - 1./6.*(x(j2) - 2.*x(j2-1) + x(j2-2))
    Flux_r = x(j2) - 0.5*(x(j2) - x(j2-1))
    Flux(j2-1) = 0.5*Flux_l*(v(j2)+abs(v(j2))) + 0.5*Flux_r*(v(j2)-abs(v(j2)))
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
  end subroutine advsy_quick_nl
!
!
  subroutine advsz_quick_nl (k1,k2,x,w,xa,adv)
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
  real     :: wold, Fold, ffip1, ffip2, ffim1
  real     :: cc, ww, Flux_l, Flux_r, phil, phir, rp1, rp2, xx
  real	   :: bc51, bc61
!                                                          !
!----------------------------------------------------------!
!
  cc = 1./4.
  Flux = 0.0
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
  xx = x(1)
  ww = w(2)
  ffip1 = (x(2) - xx) * dz1(1)
  ffip2 = (x(3) - x(2)) * dz1(2)
  ffim1 = bc51*(x(nz-4) - x(nz-5))*dz1(nz-5) + (1.-bc51)*ffip1
!
  phil = ffip1 - cc*(ffip1 - ffim1)
  phir = ffip1 + cc*(ffip2 - ffip1)
!
  Flux_l = xx + 0.5*phil/dz1w(1)
  Flux_r = x(2) - 0.5*phir/dz1w(1)
!
  Flux(1) = 0.5*Flux_l*(ww+abs(ww)) + 0.5*Flux_r*(ww-abs(ww))
!
!  Main core
!
  do k=2,nz-2
    xx = x(k)
    ww = w(k+1)
    ffip1 = (x(k+1) - xx) * dz1(k)
    ffip2 = (x(k+2) - x(k+1)) * dz1(k+1)
    ffim1 = (xx - x(k-1)) * dz1(k-1)
!
    phil = ffip1 - cc*(ffip1 - ffim1)
    phir = ffip1 + cc*(ffip2 - ffip1)
!
    Flux_l = xx + 0.5*phil/dz1w(k)
    Flux_r = x(k+1) - 0.5*phir/dz1w(k)
!
    Flux(k) = 0.5*Flux_l*(ww+abs(ww)) + 0.5*Flux_r*(ww-abs(ww))
  enddo
!
!  Top boundary
!
  k = nz-1
  xx = x(nz-1)
  ww = w(nz)
  ffip1 = (x(k+1) - xx) * dz1(k)
  ffim1 = (xx - x(k-1)) * dz1(k-1)
  ffip2 = bc61*(x(6) - x(5))*dz1(5) + (1.-bc61)*ffip1
!
  phil = ffip1 - cc*(ffip1 - ffim1)
  phir = ffip1 + cc*(ffip2 - ffip1)
!
  Flux_l = xx + 0.5*phil/dz1w(k)
  Flux_r = x(k+1) - 0.5*phir/dz1w(k)
!
  Flux(k) = 0.5*Flux_l*(ww+abs(ww)) + 0.5*Flux_r*(ww-abs(ww))
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
  end subroutine advsz_quick_nl
  
end module advection_quick_nl
