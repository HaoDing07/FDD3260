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
!  ADVQ_ppm.F90:
!
!  Purpose:
!	Calculating advection of all scalars using a PPM type
!   Finite Volume scheme, using quadratic polynomial interpolations.
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================
!
  module advection_ppm
!  
!  ---------------------------------------------
!
!
USE gridno
USE shared_data	
!
IMPLICIT NONE
  
  private
  
  character(len=10), save :: ppmlimiter='PPM'
  
  public :: ppmlimiter, advsx_ppm, advsy_ppm, advsz_ppm
  
  contains

  subroutine advsx_ppm (i1,i2,x,u,xa,adv)
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
  real, dimension(ip_start:ip_end) :: x, xa, Flux_r, Flux_l
!
  logical  :: adv
  integer  :: i, i1, i2, ii1, ii2
  real     :: ffip1, ffip2, ffim1, ffim2, phi2, phi1, phir, phil, sl, slp1
  real     :: mup, muc, dx1, uold, xmin, xplus, Flux, Flux_old
!                                                          !
!----------------------------------------------------------!
!
  dx1 = 1./dx
  ii1 = max(i1-1,2)
  ii2 = min(i2,nx-2)
!
!  Left boundary
!
  mup = abs(u(ii1+1))*dt0*dx1
  muc = abs(u(ii1))*dt0*dx1
!
  ffip1  = (x(ii1+1) - x(ii1))
  ffip2  = (x(ii1+2) - x(ii1+1))
  ffim1  = (x(ii1) - x(ii1-1))
  ffim2  = (x(ii1) - x(ii1-1))
  xmin   = 0.5*(x(ii1+1)+x(ii1))
  xplus  = x(ii1) + 0.5*ffip1
!
  call limiter_ppm (x(ii1), x(ii1+1), x(ii1-1), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
  Flux_r(ii1) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
  Flux_l(ii1) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
!  Domain core
!
  xmin = xplus
  do i=ii1+1,ii2
    mup = abs(u(i+1))*dt0*dx1
    muc = abs(u(i))*dt0*dx1
!
    ffip1 = (x(i+1) - x(i))
    ffip2 = (x(i+2) - x(i+1))
    ffim1 = (x(i) - x(i-1))
    ffim2 = (x(i-1) - x(i-2))
!
    sl = 0.5*(ffip1 + ffim1)
    slp1 = 0.5*(ffip2 + ffip1)
    xplus = x(i) + 0.5*ffip1 - 1./6.*(slp1 - sl)
!
    call limiter_ppm (x(i), x(i+1), x(i-1), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
    Flux_r(i) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
    Flux_l(i) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
    xmin = xplus
  enddo
!
!  Right boundary
!
  mup = abs(u(ii2+2))*dt0*dx1
  muc = abs(u(ii2+1))*dt0*dx1
!
  ffip1  = (x(ii2+2) - x(ii2+1))
  ffip2  = (x(ii2+2) - x(ii2+1))
  ffim1  = (x(ii2+1) - x(ii2))
  ffim2  = (x(ii2) - x(ii2-1))
  xplus  = x(ii2+1) + 0.5*ffip1
!
  call limiter_ppm (x(ii2+1), x(ii2+2), x(ii2), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
  Flux_r(ii2+1) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
  Flux_l(ii2+1) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
!  BC: non periodic
!
  if (ii1 == 2 .and. bcl(1) /= 'per') then
    Flux_l(ii1+1) = Flux_r(ii1+1)
    Flux_r(ii1) = Flux_l(ii1+2)
  endif
!
  if (ii2 == nx-2 .and. bcl(2) /= 'per') then
    Flux_r(ii2) = Flux_l(ii2)
    Flux_l(ii2+1) = Flux_r(ii2-1)
  endif
!
!  Compute advection
!
  do i=ii1+1,ii2
    Flux_old = 0.5*Flux_r(i-1)*(u(i)+abs(u(i))) + 0.5*Flux_l(i)*(u(i)-abs(u(i)))
    Flux = 0.5*Flux_r(i)*(u(i+1)+abs(u(i+1))) + 0.5*Flux_l(i+1)*(u(i+1)-abs(u(i+1)))
 !
    xa(i) = xa(i) + (Flux - Flux_old)*dx1   
  enddo
!
  if (adv) then
    uold = u(ii2+1)
    do i=ii2,ii1+1,-1
      xa(i) = xa(i) - x(i)*(uold - u(i))*dx1
      uold = u(i)
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsx_ppm
!
!
  subroutine advsy_ppm (j1,j2,x,v,xa,adv)
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
  real, dimension(jp_start:jp_end) :: x, xa, Flux_r, Flux_l
!
  logical  :: adv
  integer  :: j, j1, j2, jj1, jj2
  real     :: ffip1, ffip2, ffim1, ffim2, phi1, phi2, Flux, Flux_old, vold
  real     :: mup, muc, dy1, phil, phir, xplus, xmin, sl, slp1
!                                                          !
!----------------------------------------------------------!
!
  dy1 = 1./dy
  jj1 = max(j1-1,2)
  jj2 = min(j2,ny-2)
!
!  Left boundary
!
  mup = abs(v(jj1+1))*dt0*dy1
  muc = abs(v(jj1))*dt0*dy1
!
  ffip1  = (x(jj1+1) - x(jj1))
  ffip2  = (x(jj1+2) - x(jj1+1))
  ffim1  = (x(jj1) - x(jj1-1))
  ffim2  = (x(jj1) - x(jj1-1))
  xmin   = 0.5*(x(jj1+1)+x(jj1))
  xplus  = x(jj1) + 0.5*ffip1
!
  call limiter_ppm (x(jj1), x(jj1+1), x(jj1-1), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
  Flux_r(jj1) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
  Flux_l(jj1) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
!  Domain core
!
  xmin = xplus
  do j=jj1+1,jj2
    mup = abs(v(j+1))*dt0*dy1
    muc = abs(v(j))*dt0*dy1
!
    ffip1 = (x(j+1) - x(j))
    ffip2 = (x(j+2) - x(j+1))
    ffim1 = (x(j) - x(j-1))
    ffim2 = (x(j-1) - x(j-2))
!
    sl = 0.5*(ffip1 + ffim1)
    slp1 = 0.5*(ffip2 + ffip1)
    xplus = x(j) + 0.5*ffip1 - 1./6.*(slp1 - sl)
!
    call limiter_ppm (x(j), x(j+1), x(j-1), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
    Flux_r(j) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
    Flux_l(j) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
    xmin = xplus
  enddo
!
!  Right boundary
!
  mup = abs(v(jj2+2))*dt0*dy1
  muc = abs(v(jj2+1))*dt0*dy1
!
  ffip1 = (x(jj2+2) - x(jj2+1))
  ffip2 = (x(jj2+2) - x(jj2+1))
  ffim1 = (x(jj2+1) - x(jj2))
  ffim2 = (x(jj2) - x(jj2-1))
  xplus = x(jj2+1) + 0.5*ffip1
!
  call limiter_ppm (x(jj2+1), x(jj2+2), x(jj2), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
  Flux_r(jj2+1) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
  Flux_l(jj2+1) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
!  BC: non periodic
!
  if (jj1 == 2 .and. bcl(3) /= 'per') then
    Flux_l(jj1) = Flux_r(jj1)
    Flux_r(jj1-1) = Flux_l(jj1+1)
  endif
!
  if (jj2 == ny-2 .and. bcl(4) /= 'per') then
    Flux_r(jj2) = Flux_l(jj2)
    Flux_l(jj2+1) = Flux_r(jj2-1)
  endif
!
!  Compute advection
!
  do j=jj1+1,jj2
    Flux_old = 0.5*Flux_r(j-1)*(v(j)+abs(v(j))) + 0.5*Flux_l(j)*(v(j)-abs(v(j)))
    Flux = 0.5*Flux_r(j)*(v(j+1)+abs(v(j+1))) + 0.5*Flux_l(j+1)*(v(j+1)-abs(v(j+1)))		    
!
    xa(j) = xa(j) + (Flux - Flux_old)*dy1   
  enddo
!
  if (adv) then
    vold = v(jj2+1)
    do j=jj2,jj1+1,-1
      xa(j) = xa(j) - x(j)*(vold - v(j))*dy1
      vold = v(j)
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsy_ppm
!
!
  subroutine advsz_ppm (k1,k2,x,w,xa,adv)
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
  real, dimension(1:nz) :: x, xa, Flux_r, Flux_l
!
  logical  :: adv
  integer  :: k, k1, k2
  real     :: ffip1, ffip2, ffim1, ffim2, phi1, phi2
  real     :: deltam1, delta, deltap1, deltap2, slp1, sl
  real     :: mup, muc, wold, Flux, Flux_old, phil, phir, xplus, xmin
  real	   :: bc51, bc61
!                                                          !
!----------------------------------------------------------!
!
  Flux_r = 0.0
  Flux_l = 0.0
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
!  Bottom boundary
!
  k = 1
  mup = abs(w(k+1))*dt0*dz1(k)
  muc = abs(w(k))*dt0*dz1(k)
!
  ffip1 = x(k+1) - x(k)
  ffip2 = x(k+2) - x(k+1)
  ffim1 = bc51*(x(nz-4) - x(nz-5))
  ffim2 = bc51*(x(nz-5) - x(nz-6))
!
  xmin  = x(k) - 0.5*ffip1
  xplus = x(k) + 0.5*ffip1
!
  call limiter_ppm (x(k), x(k+1), x(k+1), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
  Flux_r(k) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
  Flux_l(k) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
  k = 2
  xmin = xplus
  mup = abs(w(k+1))*dt0*dz1(k)
  muc = abs(w(k))*dt0*dz1(k)
!
  ffip1 = x(k+1) - x(k)
  ffip2 = x(k+2) - x(k+1)
  ffim1 = (x(k) - x(k-1))
  ffim2 = bc51*(x(nz-4) - x(nz-5))
!
  xplus = x(k) + 0.5*ffip1
!
  call limiter_ppm (x(k), x(k+1), x(k-1), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
  Flux_r(k) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
  Flux_l(k) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
!  Main core
!
  xmin = xplus
  do k=3,nz-2
    mup = abs(w(k+1))*dt0*dz1(k)
    muc = abs(w(k))*dt0*dz1(k)
!
    ffip1 = (x(k+1) - x(k))
    ffip2 = (x(k+2) - x(k+1))
    ffim1 = (x(k) - x(k-1))
    ffim2 = (x(k-1) - x(k-2))
!
    delta = dz1w(k)
    deltam1 = dz1w(k-1)
    deltap1 = dz1w(k+1)
    deltap2 = dz1w(k+2)
!
    sl = delta/(deltam1+delta+deltap1)*((2.*deltam1+delta)/(deltap1+delta)*ffip1 + (delta+2.*deltap1)/(deltam1+delta)*ffim1)
    slp1 = deltap1/(delta+deltap1+deltap2)*((2.*delta+deltap1)/(deltap2+deltap1)*ffip2 + (deltap1+2.*deltap2)/(delta+deltap1)*ffip1)
    xplus = x(k) + delta/(delta+deltap1)*ffip1 + 1./(deltam1+delta+deltap1+deltap2) 	 					     &
  	  * (2.*deltap1*delta/(delta+deltap1)*((deltam1+delta)/(2.*delta+deltap1) - (deltap2+deltap1)/(2.*deltap1+delta))*ffip1      &
  	  - delta*(deltam1+delta)/(2.*delta+deltap1)*slp1 + deltap1*(deltap1+deltap2)/(2.*deltap1+delta)*sl)
!
    call limiter_ppm (x(k), x(k+1), x(k-1), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
    Flux_r(k) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
    Flux_l(k) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
    xmin = xplus
  enddo
!
!  Top boundary
!
  k = nz-1
  mup = abs(w(k+1))*dt0*dz1(k)
  muc = abs(w(k))*dt0*dz1(k)
!
  ffip1 = x(k+1) - x(k)
  ffim1 = x(k) - x(k-1)
  ffim2 = x(k-1) - x(k-2)
  ffip2 = bc61*(x(6) - x(5)) + (1.-bc61)*ffip1
!
  xplus = x(k) + 0.5*ffip1
!
  call limiter_ppm (x(k), x(k+1), x(k-1), xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
  Flux_r(k) = phir - 0.5*mup*(phi1 + (1. - 2./3.*mup)*phi2)
  Flux_l(k) = phil + 0.5*muc*(phi1 - (1. - 2./3.*muc)*phi2)
!
  k = nz
  muc = abs(w(k))*dt0*dz1(k)
  ffim1 = x(k) - x(k-1)
  Flux_l(k) = x(k) + 0.5*ffim1*(1. - muc)
!
!  Compute advection
!
  Flux_old = 0.5*Flux_r(1)*(w(2)+abs(w(2))) + 0.5*Flux_l(2)*(w(2)-abs(w(2)))
  do k=2,k2-1
    Flux = 0.5*Flux_r(k)*(w(k+1)+abs(w(k+1))) + 0.5*Flux_l(k+1)*(w(k+1)-abs(w(k+1)))
!
    xa(k) = xa(k) + (Flux - Flux_old)*dz1w(k)
    Flux_old = Flux
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
  if (bcl(5) == 'nnf') xa(1) = xa(1) + (0.5*Flux_r(1)*(w(2)+abs(w(2))) + 0.5*Flux_l(2)*(w(2)-abs(w(2))) - 0.)*dz1w(1)
  if (bcl(5) == 'nnf' .and. adv) xa(1) = xa(1) - x(1)*(w(2) - 0.)*dz1w(1)
!
  if (bcl(6) == 'nnf') xa(nz) = xa(nz) + (0. - Flux_old)*dz1w(nz)
  if (bcl(6) == 'nnf' .and. adv) xa(nz) = xa(nz) - x(nz)*(0. - w(nz))*dz1w(nz)
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsz_ppm
!
!
  subroutine limiter_ppm (xx, xp, xm, xplus, xmin, ffip1, ffip2, ffim1, ffim2, phir, phil, phi1, phi2)
!
!----------------------------------------------------------!
!                                                          !
!  	PPM detection of local extrema and limitation.     !
!	Follows to some extent Zerroukat et al. 2005	   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real  :: ffip1, ffip2, ffim1, ffim2, xx, xp, xm, xplus, xmin
  real  :: phir, phi1, phi2, phil, slope
  real  :: gamma1, gamma2, gamma3, lambda, rp1
  real  :: eps=1.e-15, thres=10.
!                                                          !
!----------------------------------------------------------!
!
  rp1 = 1.
  phir = xplus
  phil = xmin
!
!  Selective limitation criterion (Blossey and Durran, 2008)
!
  gamma1 = 0.5*((ffip1-ffim1)**2. + (ffip1+ffim1)**2.)
  gamma2 = 0.5*((ffip2-ffip1)**2. + (ffip2+ffip1)**2.)
  gamma3 = 0.5*((ffim1-ffim2)**2. + (ffim1+ffim2)**2.)
  lambda = max(gamma1,gamma2,gamma3) / (min(gamma1,gamma2,gamma3) + eps)
!
!  PPM coeficients
!
  phi1 = phir - phil
  phi2 = 3.*(phil - 2.*xx + phir)
!
  if (ppmlimiter=='PPM' .and. lambda > thres) then
!
!  Constraints on face values: selective detection and correction of local extrema at faces
!							   + correct for positivity
!
  if ( phir < 0. .or. ((phir - xx)*(xp - phir) < 0. .and. 	&
       (ffim1*ffip2 >= 0. .or. ffim1*ffip1 <= 0. .or. ffip2*ffip1 <= 0.)) ) then
    if (abs(phir - xx) < abs(phir - xp)) then
      phir = xx
    else
      phir = xp
    endif
  endif
!
  if ( phil < 0. .or. ((xx - phil)*(phil - xm) < 0. .and. 	&
       (ffim2*ffip1 >= 0. .or. ffim2*ffim1 <= 0. .or. ffip1*ffim1 <= 0.)) ) then
    if (abs(phil - xx) < abs(phil - xm)) then
      phil = xx
    else
      phil = xm
    endif
  endif
!
!  Recalculate PPM coeficients after adjusting face values
!
  phi1 = phir - phil
  phi2 = 3.*(phil - 2.*xx + phir)
!
!  Local constraint: if there is a sub-grid extremum, we replace the parabola by a limited slope
!
  slope = -0.5*phi1/(sign(1.,phi2)*max(abs(phi2),eps))
  if (abs(ffim1) > sqrt(spacing(ffim1))) rp1 = ffip1/ffim1
!
  if (abs(phi2) > eps .and. (slope > -1./2. .or. slope < 1./2.)) then
    phir = xx + 0.5*ffim1*max(0.,min(2.,0.5*(rp1+1.),2.*rp1))
    phil = xx - 0.5*ffim1*max(0.,min(2.,0.5*(rp1+1.),2.*rp1))
    if ((xx - phil)*(phir - xx) < 0.) then
      phil = xx
      phir = xx
    endif      
    phi1 = phir - phil
    phi2 = 0.	 
  endif
!
  endif
!
  return
  end
  
end module advection_ppm
