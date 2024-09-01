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
!	Interface tracking using VOF method
!
!  Author:
!	Julien Savre, MIM-LMU
!
! ================================================================
!
module vof
!
  USE gridno
  USE shared_data	
!
  IMPLICIT NONE
!
  private
!
!  LOGICAL :: lvof=.true.
  REAL, PARAMETER :: small=1.e-12, large=1.e+12
!
  interface volume
      MODULE PROCEDURE volume2d, volume3d
  end interface volume
!
  interface HEAVISIDE
      MODULE PROCEDURE HEAVISIDE0D, HEAVISIDE3D
  end interface HEAVISIDE
!
  interface LOGISTIC 
      MODULE PROCEDURE LOGISTIC0D, LOGISTIC3D
  end interface LOGISTIC 
!
  public :: volume, HEAVISIDE, LOGISTIC
!
  contains
!
! ================================================================
!
! ======================================================
  subroutine volume2d (thres, X, vol, U, W)
! ======================================================

  real, intent(in) :: thres, X(ip_start:,:)
  real, intent(in), optional :: U(ip_start:,:), W(ip_start:,:)
  real, intent(out) :: vol(ip_start:,:)
!
  integer :: i, k
!
!----------------------------------------------
!
  do k = 2, nz-1
    do i = it_start, it_end
!
      if (present(U) .and. present(W)) then
        call vof2d (k, thres, X(i-1:i+1,k-1:k+1), vol(i,k), U=U(i-1:i+2,k), W=W(i,k-1:k+2))    
      else
        call vof2d (k, thres, X(i-1:i+1,k-1:k+1), vol(i,k))
      endif
!
    enddo
  enddo
!
!----------------------------------------------
!
  return
  end subroutine volume2d
!
! ======================================================
  subroutine volume3d (thres, X, vol, U, V, W)
! ======================================================

  real, intent(in) :: thres, X(ip_start:,jp_start:,:)
  real, intent(in), optional :: U(ip_start:,jp_start:,:), V(ip_start:,jp_start:,:), W(ip_start:,jp_start:,:)
  real, intent(out) :: vol(ip_start:,jp_start:,:)
!
  integer :: i, j, k
!
!----------------------------------------------
!
  do k = 2, nz-1
    do j = jt_start, jt_end
      do i = it_start, it_end
!
        if (present(U) .and. present(V) .and. present(W)) then
          call vof3d (k, thres, X(i-1:i+1,j-1:j+1,k-1:k+1), vol(i,j,k), U=U(i-1:i+2,j,k), V=V(i,j-1:j+2,k), W=W(i,j,k-1:k+2))
        else
          call vof3d (k, thres, X(i-1:i+1,j-1:j+1,k-1:k+1), vol(i,j,k))
        endif
!
      enddo
    enddo
  enddo
!
!----------------------------------------------
!
  return
  end subroutine volume3d
!
! ======================================================
  subroutine vof3d (k, thres, X, vol, U, V, W)
! ======================================================

  integer, intent(in) :: k
  real, intent(in) :: thres, X(3,3,3)
  real, dimension(4), intent(in), optional :: U, V, W
  real, intent(out) :: vol
!
  REAL :: Xc, Xpx, Xmx, Xpy, Xmy, Xpz, Xmz, Uc, Vc, Wc
  REAL :: Xp, Yp, Zp, Xq(3), Yq(3), Zq(3)
  REAL :: dzi, H
!
!----------------------------------------------
!
  dzi = dz*fdz0(k)
  vol = 0.
!
  if ( maxval(X) > thres ) then
!
    Xp = large
    Yp = large
    Zp = large
!
    Xc = X(2,2,2)
    Xpx = X(3,2,2)
    Xmx = X(1,2,2)
    Xpy = X(2,3,2)
    Xmy = X(2,1,2)
    Xpz = X(2,2,3)
    Xmz = X(2,2,1)
!
    H = heaviside(Xc - thres)
!
!  End points position
!
    if ( (max(Xc,Xpx) >= thres .and. Xpx <= Xc) .or. (max(Xc,Xpx) >= thres .and. min(Xc,Xpx) < thres) ) then
      Xp = dx*(thres - Xc + small) / (Xpx - Xc + small)
    else if ( (max(Xc,Xmx) >= thres .and. Xmx <= Xc) .or. (max(Xc,Xmx) >= thres .and. min(Xc,Xmx) < thres) ) then
      Xp = dx*(thres - Xc + small) / (Xc - Xmx + small)
    endif
!
    if ( (max(Xc,Xpy) >= thres .and. Xpy <= Xc) .or. (max(Xc,Xpy) >= thres .and. min(Xc,Xpy) < thres) ) then
      Yp = dy*(thres - Xc + small) / (Xpy - Xc + small)
    else if ( (max(Xc,Xmy) >= thres .and. Xmy <= Xc) .or. (max(Xc,Xmy) >= thres .and. min(Xc,Xmy) < thres) ) then
      Yp = dy*(thres - Xc + small) / (Xc - Xmy + small)
    endif
!
    if ( (max(Xc,Xpz) >= thres .and. Xpz <= Xc) .or. (max(Xc,Xpz) >= thres .and. min(Xc,Xpz) < thres) ) then
      Zp = dzi*(thres - Xc + small) / (Xpz - Xc + small)
    else if ( (max(Xc,Xmz) >= thres .and. Xmz <= Xc) .or. (max(Xc,Xmz) >= thres .and. min(Xc,Xmz) < thres) ) then
      Zp = dzi*(thres - Xc + small) / (Xc - Xmz + small)
    endif
!
    xq(1) = Xp + dx/2.
    yq(1) = dy/2.
    zq(1) = dzi/2.
    xq(2) = dx/2.
    yq(2) = Yp + dy/2.
    zq(2) = dzi/2.
    xq(3) = dx/2.
    yq(3) = dy/2.
    zq(3) = Zp + dzi/2.
!
!  Advect end points if required
!
    if ( present(U) .and. present(V) .and. present(W) ) then
      Uc = 0.
      Vc = 0.
      Wc = 0.
      if ( Xp >= -dx/2. .and. Xp <= 0. ) then
        Uc = U(1) - (U(2) - U(1))*(Xp/dx)
	Vc = 0.5*(V(2) + V(3))
	Wc = 0.5*(W(2) + W(3))
      else if ( Xp > 0. .and. Xp <= dx ) then
        Uc = U(2) + (U(3) - U(2))*(Xp/dx)
	Vc = 0.5*(V(2) + V(3))
	Wc = 0.5*(W(2) + W(3))
      else if ( Xp > dx .and. Xp <= 3.*dx/2. ) then
        Uc = U(3) + (U(4) - U(3))*(Xp/dx - 1.)
	Vc = 0.5*(V(2) + V(3))
	Wc = 0.5*(W(2) + W(3))
      endif
      Xq(1) = Xq(1) + Uc*dt0
      Yq(1) = Yq(1) + Vc*dt0
      Zq(1) = Zq(1) + Wc*dt0
!
      Uc = 0.
      Vc = 0.
      Wc = 0.
      if ( Yp >= -dy/2. .and. Yp <= 0. ) then
	Uc = 0.5*(U(2) + U(3))
        Vc = V(1) - (V(2) - V(1))*(Yp/dy)
	Wc = 0.5*(W(2) + W(3))
      else if ( Yp > 0. .and. Yp <= dy ) then
	Uc = 0.5*(U(2) + U(3))
        Vc = V(2) + (V(3) - V(2))*(Yp/dy)
	Wc = 0.5*(W(2) + W(3))
      else if ( Yp > dy .and. Yp <= 3.*dy/2. ) then
	Uc = 0.5*(U(2) + U(3))
        Vc = V(3) + (V(4) - V(3))*(Yp/dy - 1.)
	Wc = 0.5*(W(2) + W(3))
      endif
      Xq(2) = Xq(2) + Uc*dt0
      Yq(2) = Yq(2) + Vc*dt0
      Zq(2) = Zq(2) + Wc*dt0
!
      Uc = 0.
      Vc = 0.
      Wc = 0.
      if ( Zp >= -dzi/2. .and. Zp <= 0. ) then
	Uc = 0.5*(U(2) + U(3))
	Vc = 0.5*(V(2) + V(3))
        Wc = W(1) - (W(2) - W(1))*(Zp/dzi)
      else if ( Zp > 0. .and. Zp <= dzi ) then
	Uc = 0.5*(U(2) + U(3))
	Vc = 0.5*(V(2) + V(3))
        Wc = W(2) + (W(3) - W(2))*(Zp/dzi)
      else if ( Zp > dzi .and. Zp <= 3.*dzi/2. ) then
	Uc = 0.5*(U(2) + U(3))
	Vc = 0.5*(V(2) + V(3))
        Wc = W(3) + (W(4) - W(3))*(Zp/dzi - 1.)
      endif
      Xq(3) = Xq(3) + Uc*dt0
      Yq(3) = Yq(3) + Vc*dt0
      Zq(3) = Zq(3) + Wc*dt0
    endif
!
!  Calculate VOF
!
    if ( abs(Xp) <= dx .and. abs(Yp) <= dy .and. abs(Zp) <= dzi ) then
      vol = VOL3D ( Xq, Yq, Zq, dx, dy, dzi )
    else if ( abs(Xp) > dx .and. abs(Yp) > dy .and. abs(Zp) > dzi .and. Xc >= thres ) then
      vol = 1.
    else if ( abs(Xp) > dx .and. abs(Yp) > dy .and. abs(Zp) > dzi .and. Xc < thres ) then
      vol = 0.
    else if ( abs(Xp) > dx .and. abs(Yp) > dy ) then
      vol = (Zp+dzi/2.)/dzi
    else if ( abs(Xp) > dx .and. abs(Zp) > dzi ) then
      vol = (Yp+dy/2.)/dy 
    else if ( abs(Zp) > dzi .and. abs(Yp) > dy ) then
      vol = (Xp+dx/2.)/dx
    else if ( abs(Xp) > dx ) then
      vol = VOL2D ( Yq(2), Yq(3), Zq(2), Zq(3), dy, dzi )
    else if ( abs(Yp) > dy ) then
      vol = VOL2D ( Xq(1), Xq(3), Zq(1), Zq(3), dx, dzi )
    else if ( abs(Zp) > dzi ) then
      vol = VOL2D ( Xq(1), Xq(2), Yq(1), Yq(2), dx, dy )
    endif
!
    vol = max(min(vol,1.),0.)
    vol = H*max(vol,(1. - vol)) + (1. - H)*min(vol,(1. - vol))
!
  endif
!
!----------------------------------------------
!
  return
!
contains
!
! ======================================================
  real function VOL3D (X, Y, Z, dxi, dyi, dzi)
! ======================================================

  Implicit none
  
  real, intent(in) :: X(3), Y(3), Z(3)
  real, intent(in) :: dxi, dyi, dzi
!
  real :: a, b, c, d
  real :: alpha, vol
  real :: x1, x2, x3, y1, y2, y3, z1, z2, z3
!
  x1 = X(1)
  y1 = Y(1)
  z1 = Z(1)
  x2 = X(2)
  y2 = Y(2)
  z2 = Z(2)
  x3 = X(3)
  y3 = Y(3)
  z3 = Z(3)
!
!  Calculate plane equation
!
  a = ( (z2-z1)/(y1*z2-y2*z1+small) - (z3-z2)/(y2*z3-y3*z2+small) ) / ( (x1*z2-x2*z1)/(y1*z2-y2*z1+small) - (x2*z3-x3*z2)/(y2*z3-y3*z2+small) ) + small
  b = ( (y2-y1)/(x1*y2-x2*y1+small) - (y3-y1)/(x1*y3-x3*y1+small) ) / ( (z1*y2-z2*y1)/(x1*y2-x2*y1+small) - (z1*y3-z3*y1)/(x1*y3-x3*y1+small) ) + small
  c = ( (x3-x2)/(z2*x3-z3*x2+small) - (x3-x1)/(z1*x3-z3*x1+small) ) / ( (y2*x3-y3*x2)/(z2*x3-z3*x2+small) - (y1*x3-y3*x1)/(z1*x3-z3*x1+small) ) + small
! 
!  Rotate axes
!
  if (a < 0.) then
    a = abs(a)
    b = b / (1. + a*dxi)
    c = c / (1. + a*dxi)
    a = a / (1. + a*dxi)
  endif
!
  if (b < 0.) then
    b = abs(b)
    a = a / (1. + b*dyi)
    c = c / (1. + b*dyi)
    b = b / (1. + b*dyi)
  endif
!
  if (c < 0.) then
    c = abs(c)
    a = a / (1. + c*dzi)
    b = b / (1. + c*dzi)
    c = c / (1. + c*dzi)
  endif
!
!  Calculate volume
!
  d = a*dxi + b*dyi + c*dzi
!
  vol = 1./(6.*a*b*c) * ( 1. - (heaviside(1. - a*dxi)*(1. - a*dxi)**3. +	  &
       				heaviside(1. - b*dyi)*(1. - b*dyi)**3. +	  &
       				heaviside(1. - c*dzi)*(1. - c*dzi)**3.)		  &
			     + (heaviside(1. - d + a*dxi)*(1. - d + a*dxi)**3. +  &
			        heaviside(1. - d + b*dyi)*(1. - d + b*dyi)**3. +  &
				heaviside(1. - d + c*dzi)*(1. - d + c*dzi)**3.) )
!
!  Volume fraction
!
  VOL3D = vol / (dxi*dyi*dzi)
!
!----------------------------------------------
!
  end function VOL3D
!
! ======================================================
  real function VOL2D (X1, X2, Z1, Z2, dxi, dzi)
! ======================================================

  Implicit none
  
  real, intent(in) :: x1, x2, z1, z2
  real, intent(in) :: dxi, dzi
!
  real :: a, b, vol
!
!  Calculate plane equation
!
  a = (z2-z1)/(x1*z2-x2*z1+small) + small
  b = (x2-x1)/(x2*z1-x1*z2+small) + small
! 
!  Rotate axes
! 
  if (a < 0.) then
    a = abs(a)
    b = b / (1. + a*dxi)
    a = a / (1. + a*dxi)
  endif
!
  if (b < 0.) then
    b = abs(b)
    a = a / (1. + b*dzi)
    b = b / (1. + b*dzi)
  endif
!
!  Calculate volume
!
  vol = 1./(2.*a*b) * ( 1. - heaviside(1. - a*dxi)*(1. - a*dxi)**2. - heaviside(1. - b*dzi)*(1. - b*dzi)**2. )
!
!  Volume fraction
!
  VOL2D = vol / (dxi*dzi)
!
!----------------------------------------------
!
  end function VOL2D
!
  end subroutine vof3d
!
! ======================================================
  subroutine vof2d (k, thres, X, vol, U, W)
! ======================================================

  integer, intent(in) :: k
  real, intent(in) :: thres, X(3,3)
  real, dimension(3), intent(in), optional :: U, W
  real, intent(out) :: vol
!
  REAL :: Xc, Xpx, Xmx, Xpz, Xmz, Uc, Wc
  REAL :: dzi, Xp, Zp, H
!
!----------------------------------------------
!
  dzi = dz*fdz0(k)
  vol = 0.
!
  if ( maxval(X) > thres ) then
    Xp = large
    Zp = large
!
    Xpx = X(3,2)
    Xmx = X(1,2)
    Xpz = X(2,3)
    Xmz = X(2,1)
    Xc = X(2,2)
    H = heaviside(Xc - thres)
!
    if ( maxval(X(:,2)) > thres .and. abs(Xpx - Xmx) > small ) &
      Xp = 2.*dx*abs(Xc - thres)/abs(Xpx - Xmx)
!
    if ( maxval(X(2,:)) > thres .and. abs(Xpz - Xmz) > small ) &
      Zp = 2.*dzi*abs(Xc - thres)/abs(Xpz - Xmz)
!
!  Advection if requested
!
    if ( present(U) .and. present(W) ) then
      if ( (thres - Xc)/(Xpx - Xc + small) >= 0. .and. Xp <= dx ) then
        Uc = U(2) + (U(3) - U(2))*(1. - abs(Xp)/dx)
        Xp = Xp + Uc*dt0
      else if ( (thres - Xc)/(Xmx - Xc + small) >= 0. .and. Xp <= dx ) then
        Uc = U(2) - (U(2) - U(1))*(1. - abs(Xp)/dx)
        Xp = Xp + Uc*dt0
      endif
!
      if ( (thres - Xc)/(Xpz - Xc + small) >= 0. .and. Zp <= dzi ) then
        Wc = W(2) + (W(3) - W(2))*(1. - abs(Zp)/dzi)
        Zp = Zp + Wc*dt0      
      else if ( (thres - Xc)/(Xmz - Xc + small) >= 0. .and. Zp <= dzi ) then
        Wc = W(2) - (W(2) - W(1))*(1. - abs(Zp)/dzi)
        Zp = Zp + Wc*dt0        
      endif
    endif
!
!  Calculate VOF
!
    if (Xp <= dx .and. Zp <= dzi) then
      vol = VOL2D ( Xp, Zp, dx, dzi )
    else
      vol = 1.
    endif
    vol = vol*H + (1. - vol)*(1. - H)    
  endif
!
!----------------------------------------------
!
  return
!
contains
!
! ======================================================
  real function VOL2D (X, Z, dxi, dzi)
! ======================================================

  Implicit none
  
  real :: X, Z
  real :: dxi, dzi
!
  real :: alpha, vol
  REAL, DIMENSION(2) :: normal
  REAL, DIMENSION(2,2) :: coor
!
!  Change coordinates
!
  coor = 0.
  coor(1,1) = X*(1. + dzi/2./Z) + dxi/2.
  coor(2,2) = Z*(1. + dxi/2./X) + dzi/2.
!
!  Calculate normal and distance
!
  alpha = 1.
  normal(1) = 1. / coor(1,1)
  normal(2) = 1. / coor(2,2)
!
!  Calculate volume
!
  vol = 1./(2.*normal(1)*normal(2)) * ( alpha**2. -	  		  &
      (heaviside(alpha - normal(1)*dxi)*(alpha - normal(1)*dxi)**2. +	  &
       heaviside(alpha - normal(2)*dzi)*(alpha - normal(2)*dzi)**2.) )
!
!  Volume fraction
!
  vol = vol / (dxi*dzi)
  VOL2D = max(min(vol,1.),0.)
!
!----------------------------------------------
!
  end function VOL2D
!
  end subroutine vof2d
!
! ======================================================
  subroutine marching_squares (k, thres, X, vol)
! ======================================================

  integer, intent(in) :: k
  real, intent(in) :: thres, X(3,3)
  real, intent(out) :: vol
!
  REAL :: Xc, X1, X2, X3, X4
  REAL :: dzi
!
!----------------------------------------------
!
  vol = 0.
  dzi = dz*fdz0(k)
!
  if ( X(2,2) > small ) then
    Xc = X(2,2)
    X1 = 0.25*(X(1,1) + X(1,2) + X(2,1) + X(2,2))
    X2 = 0.25*(X(2,1) + X(2,2) + X(3,1) + X(3,2))
    X3 = 0.25*(X(2,2) + X(2,3) + X(3,2) + X(3,3))
    X4 = 0.25*(X(1,2) + X(1,3) + X(2,2) + X(2,3))
!
    vol = vol + triangle( thres, Xc, X1, X2, dx, dzi )
!
    vol = vol + triangle( thres, Xc, X3, X4, dx, dzi )
!
    vol = vol + triangle( thres, Xc, X4, X1, dzi, dx )
!
    vol = vol + triangle( thres, Xc, X2, X3, dzi, dx )
!
    vol = max(min(vol,1.),0.)
  endif
!
!----------------------------------------------
!
  return
!
contains
!
! ======================================================
  real function triangle (thres, X1, X2, X3, dxi, dzi)
! ======================================================

  Implicit none
  
  real :: X1, X2, X3, dxi, dzi, thres
  real :: vol, volt
  real :: dx1, dx2, dx3, dz1, dz2
!
!  Find surface of sub-triangle
!
    vol = 0.
    volt = 0.25*dxi*dzi
!
    if ( X2 > thres .and. X3 > thres .and. X1 > thres ) then
      vol = volt
    else if ( .not.(X2 <= thres .and. X3 <= thres .and. X1 <= thres) ) then
    
      if ( (X2 > thres .and. X1 <= thres) .or. (X2 <= thres .and. X1 > thres) ) then
        dx1 = abs( 0.5*(thres - X2)/(X1 - X2)*dxi )
	dz1 = abs( 0.5*(thres - X2)/(X1 - X2)*dzi )
      else if ( (X2 > thres .and. X3 <= thres) .or. (X2 <= thres .and. X3 > thres) ) then
        dx1 = abs( (thres - X2)/(X3 - X2)*dxi )
	dz1 = 0.
      endif
      
      if ( (X3 > thres .and. X1 <= thres) .or. (X3 <= thres .and. X1 > thres) ) then
        dx2 = abs( 0.5*(thres - X3)/(X1 - X3)*dxi )
	dz2 = abs( 0.5*(thres - X3)/(X1 - X3)*dzi )
      else if ( (X3 > thres .and. X2 <= thres) .or. (X3 <= thres .and. X2 > thres) ) then
        dx2 = abs( (thres - X3)/(X2 - X3)*dxi )
	dz2 = 0.
      endif      
      
      dx3 = dxi
      
      vol = 0.5*dz1*dx2 + 0.5*dz2*(dx3-dx1)
    endif
!
    vol = heaviside(0.3333*(X1+X2+X3)-thres)*max(vol,volt-vol) + 	&
    	  heaviside(-(0.3333*(X1+X2+X3)-thres))*min(vol,volt-vol)
! 
    triangle = vol / (dxi*dzi)
!
  end function triangle
!
  end subroutine marching_squares
!
! ======================================================
  FUNCTION HEAVISIDE0D (X) RESULT (H)
! ======================================================
!
  IMPLICIT NONE
!
  REAL, intent(in) :: X
  REAL :: H
!
  H = max(sign(1.,X),0.)
  
  END FUNCTION HEAVISIDE0D
!
! ======================================================
  FUNCTION HEAVISIDE3D (X)  RESULT (H)
! ======================================================
!
  IMPLICIT NONE
!
  REAL, DIMENSION(ip_start:,jp_start:,:), intent(in) :: X
  REAL, DIMENSION(size(X,1),size(X,2),size(X,3)) :: H
!
  H = max(sign(1.,X),0.)
  
  END FUNCTION HEAVISIDE3D
!
! ======================================================
  FUNCTION LOGISTIC0D (X,EPS) RESULT (L)
! ======================================================
!
  IMPLICIT NONE
!
  REAL, intent(in) :: X, EPS
  REAL :: L
!
  L = 1. / (1. + exp(-X/EPS)) 
  
  END FUNCTION LOGISTIC0D
!
! ======================================================
  FUNCTION LOGISTIC3D (X,EPS)  RESULT (L)
! ======================================================
!
  IMPLICIT NONE
!
  REAL, intent(in) :: EPS
  REAL, DIMENSION(ip_start:,jp_start:,:), intent(in) :: X
  REAL, DIMENSION(size(X,1),size(X,2),size(X,3)) :: L
!
  L = 1. / (1. + exp(-X/EPS))
  
  END FUNCTION LOGISTIC3D 
  
end module vof
