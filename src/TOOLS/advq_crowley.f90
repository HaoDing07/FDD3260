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
!  ADVQ_TVD.F:
!	subroutine qnadv_tvd(xmin,u,v,w,x,s)
!
!  Purpose:
!	Calculating advection of all scalars using a MUSCL type
!   Finite Volume scheme, using linear interpolations.
!
!  Author:
!	Julien Savre
!       MISU
!
! ================================================================
!

  module advection_crow
!  
!  ---------------------------------------------
!
!
USE gridno
USE shared_data	
!
IMPLICIT NONE
  
  private
  
  public :: advsx_crow, advsy_crow, advsz_crow
  
  contains

  subroutine advsx_crow (i1,i2,x,u,xa,adv)
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
  integer  :: i, i1, i2, ii1, ii2
  real     :: fr, fr2, fr3
  real     :: dx1, dt1, uu, mu
!                                                          !
!----------------------------------------------------------!
!
  dt1 = 1./dt0
  dx1 = 1./dx
  ii1 = max(i1-1,2)
  ii2 = min(i2,nx-2)
  mu = 0.
  Flux = 0.0
!
!  Calculate fluxes
!
  do i = ii1,ii2
    uu = u(i+1)
    fr = (x(i+1) - x(i))
    fr2 = (x(i+2) - x(i+1)) - (x(i) - x(i-1))
    fr3 = (x(i+2) - x(i+1)) + (x(i) - x(i-1))
!
    Flux(i) = uu*(x(i) + 0.5*fr - 1./12.*fr2)
!
    mu = dt0*dx1*abs(uu)
    Flux(i) = Flux(i) + (1./24.*fr3 - 7./12.*fr) * dx*dt1*mu**2.
    Flux(i) = Flux(i) + (1./12.*fr2) * dx*dt1*mu**3.
    Flux(i) = Flux(i) + (1./12.*fr - 1./24.*fr3) * dx*dt1*mu**4.
  enddo
!
!  Compute advection
!
  do i=ii1+1,ii2
    xa(i) = xa(i) + (Flux(i) - Flux(i-1))*dx1
  enddo
!
  if (adv) then
    do i=ii1+1,ii2
      xa(i) = xa(i) - x(i)*(u(i+1) - u(i))*dx1
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsx_crow
!
!
  subroutine advsy_crow (j1,j2,x,v,xa,adv)
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
  integer  :: j, j1, j2, jj1, jj2
  real     :: fr, fr2, fr3
  real     :: dy1, dt1, mu, vv
!                                                          !
!----------------------------------------------------------!
!
  dt1 = 1./dt0
  dy1 = 1./dy
  jj1 = max(j1-1,2)
  jj2 = min(j2,ny-2)
  mu = 0.
  Flux = 0.0
!
!  Calculate fluxes
!
  do j = jj1, jj2
    vv  = v(j+1)
    fr  = (x(j+1) - x(j))
    fr2 = (x(j+2) - x(j+1)) - (x(j) - x(j-1))
    fr3 = (x(j+2) - x(j+1)) + (x(j) - x(j-1))
!
    Flux(j) = vv*(x(j) + 0.5*fr - 1./12.*fr2)
!
    mu = dt0*dy1*abs(vv)
    Flux(j) = Flux(j) + (1./24.*fr3 + 7./12.*fr) * dt1*dy*mu**2.
    Flux(j) = Flux(j) + (1./12.*fr2) * dt1*dy*mu**3.
    Flux(j) = Flux(j) + (1./12.*fr - 1./24.*fr3) * dt1*dy*mu**4.
  enddo
!
!  Compute advection
!
  do j=jj1+1,jj2
    xa(j) = xa(j) + (Flux(j) - Flux(j-1))*dy1
  enddo
!
  if (adv) then
    do j=jj1+1,jj2
      xa(j) = xa(j) - x(j)*(v(j+1) - v(j))*dy1
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsy_crow
!
!
  subroutine advsz_crow (k1,k2,x,w,xa,adv)
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
  real     :: fr, fr2, fr3
  real     :: ww, dt1, mu
!                                                          !
!----------------------------------------------------------!
!
  dt1 = 1./dt0
  mu = 0.
  Flux = 0.
!
!  Top boundary
!
  k = nz-1
  ww  = w(k+1)
  fr  = (x(nz) - x(nz-1))*dz1(nz-1)
  fr2 = (x(nz) - x(nz-1))*dz1(nz-1) - (x(nz-1) - x(nz-2))*dz1(nz-2)
!
  Flux(nz-1) = ww*(x(nz-1) + (0.5*fr - 1./6.*fr2)/dz1w(nz-1))
!
  mu = dt0*dz1(nz-1)*abs(ww)
  Flux(nz-1) = Flux(nz-1) - (0.5*fr)/dz1w(k2-1) * dt1/dz1w(nz-1)*mu**2.
  Flux(nz-1) = Flux(nz-1) + (1./6.*fr2)/dz1w(k2-1) * dt1/dz1w(nz-1)*mu**3.
  if (bcl(6) == 'nnf') Flux(nz) = 0.
!
!  Main core
!
  do k=k1+1,k2-2
    ww  = w(k+1)
    fr  = (x(k+1) - x(k))*dz1(k)
    fr2 = (x(k+2) - x(k+1))*dz1(k+1) - (x(k) - x(k-1))*dz1(k-1)
    fr3 = (x(k+2) - x(k+1))*dz1(k+1) + (x(k) - x(k-1))*dz1(k-1)
!
    Flux(k) = ww*(x(k) + (0.5*fr - 1./12.*fr2)/dz1w(k))
!
    mu = dt0*dz1(k)*abs(ww)
    Flux(k) = Flux(k) + (1./24.*fr3 - 7./12.*fr)/dz1w(k) * dt1/dz1w(k)*mu**2.
    Flux(k) = Flux(k) + (1./12.*fr2)/dz1w(k) * dt1/dz1w(k)*mu**3.
    Flux(k) = Flux(k) + (1./12.*fr - 1./24.*fr3)/dz1w(k) * dt1/dz1w(k)*mu**4.
  enddo
!
!  Bottom boundary
!
  ww = w(k1+1)
  fr = (x(2) - x(1))*dz1(1)
  fr2 = (x(3) - x(2))*dz1(2) - (x(2) - x(1))*dz1(1)
  if (bcl(5) == 'nnf') fr3 = (x(2) - x(1))*dz1(1) - (x(1) - x(2))*dz1(1)
  if (bcl(5) == 'ope') fr3 = (x(2) - x(1))*dz1(1)
  if (bcl(5) == 'per') fr3 = (x(2) - x(1))*dz1(1) - (x(nz-4) - x(nz-5))*dz1(nz-4)
!
  Flux(k1) = ww*(x(k1) + (0.5*fr - 1./6.*fr2)/dz1w(k1))
!
  mu = dt0*dz1(k1)*abs(ww)
  Flux(k1) = Flux(k1) - (0.5*fr)/dz1w(k1) * dt1/dz1w(k1)*mu**2.
  Flux(k1) = Flux(k1) + (1./6.*fr3)/dz1w(k1) * dt1/dz1w(k1)*mu**2.
!
!  Compute advection
!
  do k=2,k2-1
    xa(k) = xa(k) + (Flux(k) - Flux(k-1))*dz1w(k)
  enddo
!
  if (adv) then
    do k=2,k2-1
      xa(k) = xa(k) - x(k)*(w(k+1) - w(k))*dz1w(k)
    enddo
  endif
!
  if (bcl(5) == 'nnf') xa(1) = xa(1) + (Flux(1) - 0.)*dz1w(1)
  if (bcl(5) == 'nnf' .and. adv) xa(1) = xa(1) - x(1)*(w(2) - 0.)*dz1w(1)
!
  if (bcl(6) == 'nnf') xa(k2) = xa(k2) + (0. - Flux(k2-1))*dz1w(k2)
  if (bcl(6) == 'nnf' .and. adv) xa(k2) = xa(k2) - x(k2)*(0. - w(k2))*dz1w(k2)
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine advsz_crow
  
end module advection_crow
