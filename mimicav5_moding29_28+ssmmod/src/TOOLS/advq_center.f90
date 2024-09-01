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
  module advection_center
!  
!  ---------------------------------------------
!
!
USE gridno
USE shared_data	
!
IMPLICIT NONE
  
  private
  
  public :: advsx_center, advsy_center, advsz_center
  
  contains

  subroutine advsx_center (i1,i2,x,u,xa,adv)
!
!----------------------------------------------------------!
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
  real     :: dx1, uu, uold, Fold, slope_p 
!                                                          !
!----------------------------------------------------------!
!
  dx1 = 1./dx
  Flux = 0.
!
!  Calculate fluxes
!
  do i=i1-1,i2
    uu = u(i+1)
    slope_p  = 0.5*(x(i+1) + x(i)) 
    Flux(i) = uu*slope_p 
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
  end subroutine advsx_center
!
!
  subroutine advsy_center (j1,j2,x,v,xa,adv)
!
!----------------------------------------------------------!
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
  real     :: dy1, vv, vold, Fold, slope_p
!                                                          !
!----------------------------------------------------------!
!
  dy1 = 1./dy
  Flux = 0.
!
!  Calculate fluxes
!
  do j=j1-1,j2
    vv = v(j+1)
    slope_p  = 0.5*(x(j+1) + x(j))
    Flux(j) = vv*slope_p 
  enddo
!
#ifdef CHANNEL
  if (jp_start == 1) then
    Flux(j1-1) = 0.
  endif
!
  if (jp_end == ny) then
    Flux(j2) = 0.
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
  end subroutine advsy_center
!
!
  subroutine advsz_center (k1,k2,x,w,xa,adv)
!
!----------------------------------------------------------!
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
  real     :: ww, wold, Fold, slope_p, dzw(nz)
  real	   :: bc51, bc61
!                                                          !
!----------------------------------------------------------!
!
  dzw = 1./dz1
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
    slope_p  = 0.5*(x(k+1) + x(k))
    Flux(k) = ww*slope_p 
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
    wold = w(nz)
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
  end subroutine advsz_center
!
end module advection_center
