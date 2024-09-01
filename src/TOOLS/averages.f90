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
!  FUNCPACK:
!	Package of functions
!
!  Purpose:
!	Put functions and utility subroutines together.
!
!  Author:
!	Chien Wang
!	Center for Clouds, Chemistry, and Climate
!
! ================================================================

  module averages

#ifdef SPMD
  USE mpi
  USE mpicomm
#endif
  use gridno
  use shared_data
  use shared_pressure, only : den0
  use boundary_conditions
!
  IMPLICIT NONE
!
  private

  interface horav
      MODULE PROCEDURE horav_all, horav_all_2df, horav_2d
  end interface horav

  interface horsum
      MODULE PROCEDURE horsum_all, horsum_2d
  end interface horsum
  
  interface vertint
      MODULE PROCEDURE vertint_all, vertint_1d
  end interface vertint

  interface max_all
      MODULE PROCEDURE max_all_3d, max_all_2d
  end interface max_all

  interface av_all
      MODULE PROCEDURE av_all_3d
  end interface

  public :: horav, horsum, totalav, totalsum, boxav, vertint, cal_fluxx, cal_flux, cal_fluxu, cal_fluxv, max_all, av_all

  CONTAINS
!
! ==================================================
  subroutine horav_all ( x, xav, xq1, xq2, xc1, xc2, filter1, filter2, xflag1, xflag2, filter3 )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates horizontal average of a 3D scalar            !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, dimension(ip_start:,jp_start:,:), intent(in), optional :: xq1, xq2
  real, dimension(nz), intent(out) :: xav
  real, intent(in), optional :: xc1, xc2
  logical, dimension(ip_start:,jp_start:), intent(in), optional :: xflag1, xflag2
  character(len=*), intent(in), optional :: filter1, filter2, filter3
!
  integer  :: i, j, k, ierr
  real :: crit, c
  real, dimension(nz) :: xav0, ii, iitot
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: f1, f2
!                                                          !
!----------------------------------------------------------!
!
!  Calculate filter 1
!
  f1(ip_start:ip_end,jp_start:jp_end,1:nz) = 1.
  if (present(xq1)) then
    if (present(xc1)) then
      crit = xc1
    else
      crit = 0.
    endif
    if (present(filter1)) then
      if (trim(filter1) == 'lower') then
        where (xq1 > crit) f1 = 0. 
      else if (trim(filter1) == 'larger') then
        where (xq1 <= crit) f1 = 0. 
      endif
    else
      where (xq1 <= crit) f1 = 0. 
    endif
    if (present(xflag1)) then
      do k = 1, nz
        where (.not.xflag1) f1(:,:,k) = 0.
      enddo
    endif
  endif
!
!  Calculate filter 2
!
  f2(ip_start:ip_end,jp_start:jp_end,1:nz) = f1(ip_start:ip_end,jp_start:jp_end,1:nz)
  if (present(xq2)) then
    if (present(xc2)) then
      crit = xc2
    else
      crit = 0.
    endif
    if (present(filter2)) then
      if (trim(filter2) == 'lower') then
        where (xq2 > crit) f2 = 0. 
      else if (trim(filter2) == 'larger') then
        where (xq2 <= crit) f2 = 0. 
      endif
    else
      where (xq2 <= crit) f2 = 0. 
    endif
    if (present(xflag2)) then
      do k = 1, nz
        where (.not.xflag2) f2(:,:,k) = 0.
      enddo
    endif
  endif
!
! Filter 3?
!
  if (present(filter3)) then
    f1(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.
    if (trim(filter3) == 'inside') then
      do k = 1, nz
        do j = jt_start, jt_end
          do i = it_start, it_end
            if (f2(i,j,k) == 1. .and. (f2(i+1,j,k)==0. .or. f2(i-1,j,k)==0. .or. f2(i,j+1,k)==0. .or. f2(i,j-1,k)==0.)) then
              f1(i,j,k) = 1.
            endif
          enddo
        enddo
      enddo
    else if (trim(filter3) == 'outside') then
      do k = 1, nz
        do j = jt_start, jt_end
          do i = it_start, it_end
            if (f2(i,j,k) == 0. .and. (f2(i+1,j,k)==1. .or. f2(i-1,j,k)==1. .or. f2(i,j+1,k)==1. .or. f2(i,j-1,k)==1.)) then
              f1(i,j,k) = 1.
            endif
          enddo
        enddo
      enddo
    endif
    f2(ip_start:ip_end,jp_start:jp_end,1:nz) = f1(ip_start:ip_end,jp_start:jp_end,1:nz)
  endif
!
!  Cumulate horizontally with filter
!
  ii = 0.0
  xav0 = 0.0
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        c = f2(i,j,k)
  	xav0(k) = xav0(k) + c*x(i,j,k)
	ii(k) = ii(k) + c
      enddo
    enddo
  enddo
!  
#ifdef SPMD
  xav = xav0
  iitot = ii
  call MPI_ALLReduce (xav0(1), xav(1), nz, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (ii, iitot, nz, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  xav = xav/max(iitot,1.)
#else
  xav = xav0/max(ii,1.)
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine horav_all 
!
! ==================================================
  subroutine horav_all_2df ( x, xav, xq1, xc1, filter1, xq2, xc2, filter2 )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates horizontal average of a 3D scalar            !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, dimension(ip_start:,jp_start:), intent(in) :: xq1
  real, dimension(ip_start:,jp_start:), intent(in), optional :: xq2
  real, intent(in) :: xc1
  real, intent(in), optional :: xc2
  character(len=*), intent(in) :: filter1
  character(len=*), intent(in), optional :: filter2
  real, dimension(nz), intent(out) :: xav
!
  integer  :: i, j, k, ierr
  real :: crit, c
  real, dimension(nz) :: xav0, ii, iitot
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: f1, f2
!                                                          !
!----------------------------------------------------------!
!
!  Calculate filter 1
!
  f1(ip_start:ip_end,jp_start:jp_end) = 1.
  if (trim(filter1) == 'lower') then
    where (xq1 > xc1) f1 = 0. 
  else if (trim(filter1) == 'larger') then
    where (xq1 <= xc1) f1 = 0. 
  endif
!
!  Calculate filter 2
!
  f2(ip_start:ip_end,jp_start:jp_end) = f1(ip_start:ip_end,jp_start:jp_end)
  if (present(xq2)) then
    if (present(xc2)) then
      crit = xc2
    else
      crit = 0.
    endif
    if (present(filter2)) then
      if (trim(filter2) == 'lower') then
        where (xq2 > crit) f2 = 0. 
      else if (trim(filter2) == 'larger') then
        where (xq2 <= crit) f2 = 0. 
      endif
    else
      where (xq2 <= crit) f2 = 0. 
    endif
  endif
!
!  Cumulate horizontally with filter
!
  ii = 0.0
  xav0 = 0.0
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        c = f2(i,j)
  	xav0(k) = xav0(k) + c*x(i,j,k)
	ii(k) = ii(k) + c
      enddo
    enddo
  enddo
!  
#ifdef SPMD
  xav = xav0
  iitot = ii
  call MPI_ALLReduce (xav0(1), xav(1), nz, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (ii, iitot, nz, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  xav = xav/max(iitot,1.)
#else
  xav = xav0/max(ii,1.)
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine horav_all_2df
!
! ==================================================
  subroutine horav_2d (x, xav, xq, xc, filter)
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates horizontal average of a 2D scalar            !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:), intent(in) :: x
  real, dimension(ip_start:,jp_start:), intent(in), optional :: xq
  real, optional :: xc
  character(len=*), intent(in), optional :: filter
  real, intent(out) :: xav
!
  integer  :: i, j, ii, itot, ierr
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: f
  real :: xav0, c
!                                                          !
!----------------------------------------------------------!
!
  f(ip_start:ip_end,jp_start:jp_end) = 1.
  if (trim(filter) == 'lower') then
    where (xq > xc) f = 0. 
  else if (trim(filter) == 'larger') then
    where (xq <= xc) f = 0. 
  endif
!
  ii = 0
  xav0 = 0.0
  do j = jt_start, jt_end
    do i = it_start, it_end
      c = f(i,j)
      xav0 = xav0 + c*x(i,j)
      ii = ii + c
    enddo
  enddo
!  
#ifdef SPMD
  xav = xav0
  itot = ii
  call MPI_ALLReduce (xav0, xav, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (ii, itot, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  xav = xav/real(max(itot,1))
#else
  xav = xav0/real(max(ii,1))
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine horav_2d
!
! ==================================================
  subroutine horsum_all (x, xav, xq1, xc1, filter1, xflag1, xq2, xc2, filter2, xflag2, scale)
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates total sum of a 3D scalar                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, dimension(ip_start:,jp_start:,:), intent(in), optional :: xq1, xq2
  real, intent(in), optional :: xc1, xc2, scale
  character(len=*), intent(in), optional :: filter1, filter2
  logical, dimension(ip_start:,jp_start:), intent(in), optional :: xflag1, xflag2
  real, dimension(nz), intent(out) :: xav
!
  real :: crit, c, fac
  real, dimension(nz) :: xav0
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: f1, f2
  integer  :: ierr, i, j, k
!                                                          !
!----------------------------------------------------------!
!
!  Scale sum
!
  if (present(scale)) then
    fac = scale
  else
    fac = 1.
  endif
!
!  Calculate filter 1
!
  f1 = 1.
!
  if (present(xq1)) then
    if (present(xc1)) then
      crit = xc1
    else
      crit = 0.
    endif
    if (present(filter1)) then
      if (trim(filter1) == 'lower') then
        where (xq1 > crit) f1 = 0. 
      else if (trim(filter1) == 'larger') then
        where (xq1 <= crit) f1 = 0. 
      endif
    else
      where (xq1 <= crit) f1 = 0. 
    endif
    if (present(xflag1)) then
      do k = 1, nz
        where (.not.xflag1) f1(:,:,k) = 0.
      enddo
    endif
  endif
!
!  Calculate filter 2
!
  f2 = f1
!
  if (present(xq2)) then
    if (present(xc2)) then
      crit = xc2
    else
      crit = 0.
    endif
    if (present(filter2)) then
      if (trim(filter2) == 'lower') then
        where (xq2 > crit) f2 = 0. 
      else if (trim(filter2) == 'larger') then
        where (xq2 <= crit) f2 = 0. 
      endif
    else
      where (xq2 <= crit) f2 = 0. 
    endif
    if (present(xflag2)) then
      do k = 1, nz
        where (.not.xflag2) f2(:,:,k) = 0.
      enddo
    endif
  endif
!
!  Cumulate horizontally with filter
!
  xav0 = 0.0
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        c = f2(i,j,k)
  	xav0(k) = xav0(k) + fac*c*x(i,j,k)
      enddo
    enddo
  enddo
!  
#ifdef SPMD
  xav = xav0
  call MPI_ALLReduce (xav0(1), xav(1), nz, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  xav = xav0
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine horsum_all
!
! ==================================================
  subroutine horsum_2d (x, xav)
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates total sum of a 2D scalar                     !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:), intent(in) :: x
  real, intent(out) :: xav
!
  integer  :: i, j, ierr
  real :: xav0
!                                                          !
!----------------------------------------------------------!
!
  xav0 = 0.0
  do j = jt_start, jt_end
    do i = it_start, it_end
      xav0 = xav0 + x(i,j)
    enddo
  enddo
!  
#ifdef SPMD
  xav = xav0
  call MPI_ALLReduce (xav0, xav, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  xav = xav0
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine horsum_2d
!
! ==================================================
  subroutine totalav ( x, xav )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates horizontal average of a scalar               !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, intent(out) :: xav
!
  integer  :: i, j, k, ii, itot, ierr
  real :: xav0
!                                                          !
!----------------------------------------------------------!
!
  ii = 0
  xav0 = 0.0
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        xav0 = xav0 + x(i,j,k)
        ii = ii + 1
      enddo
    enddo
  enddo
!  
#ifdef SPMD
  xav = xav0
  itot = ii
  call MPI_ALLReduce (xav0, xav, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (ii, itot, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  xav = xav/real(itot)
#else
  xav = xav0/real(ii)
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine totalav
!
! ==================================================
  subroutine boxav ( x, xav, m )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates horizontal average of a scalar               !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  integer, intent(in) :: m
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, dimension(ip_start:,jp_start:,:), intent(out) :: xav
!
  integer  :: i, j, k, ii, jj, ii1, jj1, ierr
  real, dimension(nx,ny) :: x2d
  real :: nn1
!                                                          !
!----------------------------------------------------------!
!
  xav = 0.
!
  nn1 = 1. / real( (2*m + 1)*(2*m + 1) )
!
  do k = 1, nz
!
#ifdef SPMD
    call collect ( x(:,:,k), x2d )
!
    call MPI_Bcast ( x2d, nx*ny, REALTYPE, 0, MPI_COMM_WORLD, ierr )
#else
    x2d = x(:,:,k)
#endif
!
    do j = jt_start, jt_end
      do i = it_start, it_end 
!
        do jj = j-m, j+m
          do ii = i-m, i+m
            ii1 = ii
            jj1 = jj
            if ( ii1 > nx ) ii1 = ii1 - nx
            if ( jj1 > ny ) jj1 = jj1 - ny
            if ( ii1 < 1 ) ii1 = ii1 + nx
            if ( jj1 < 1 ) jj1 = jj1 + ny
!
            xav(i,j,k) = xav(i,j,k) + x2d(ii1,jj1)*nn1
          enddo
        enddo
! 
      enddo
    enddo
!
  enddo
!
!  BCs
!
  call statebc ( xav )
!
  return
!   
!----------------------------------------------------------!
!
  end subroutine boxav
!
! ==================================================
  subroutine totalsum ( x, xav )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates the full 3D sum of a scalar                  !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, intent(out) :: xav
!
  integer  :: i, j, k, ierr
  real :: xav0
!                                                          !
!----------------------------------------------------------!
!
  xav0 = 0.0
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        xav0 = xav0 + x(i,j,k)
      enddo
    enddo
  enddo
!  
#ifdef SPMD
  xav = xav0
  call MPI_ALLReduce (xav0, xav, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  xav = xav0
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine totalsum
!
! ==================================================
  subroutine vertint_all ( x, xav, zlim, zbas, filter, weighted )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates vertical integral of a 3D scalar             !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, dimension(ip_start:,jp_start:), intent(out) :: xav
  real, dimension(ip_start:,jp_start:,:), intent(in), optional :: filter
  real, optional, intent(in) :: zlim, zbas
  logical, optional, intent(in) :: weighted
!
  integer  :: k
  real :: weight(nz)
!                                                          !
!----------------------------------------------------------!
!
  if (present(weighted)) then
    if (weighted) then
      weight = den0
    else
      weight = 1.
    endif
  else
    weight = 1.
  endif
!
  xav = 0.0
!
  if (present(zlim)) then
    do k = 1, nz
      if (z0(k) <= zlim) xav = xav + weight(k)*x(:,:,k)*dz*fdz0(k)
    enddo
  else if (present(zbas)) then
    do k = 1, nz
      if (z0(k) > zbas) xav = xav + weight(k)*x(:,:,k)*dz*fdz0(k)
    enddo
  else if (present(filter)) then
    do k = 1, nz
      where (filter(:,:,k) > qthres) xav = xav + weight(k)*x(:,:,k)*dz*fdz0(k)
    enddo
  else
    do k = 1, nz
      xav = xav + weight(k)*x(:,:,k)*dz*fdz0(k)
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine vertint_all 
!
! ==================================================
  subroutine vertint_1D ( x, xav, zlim, zbas, weighted )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates vertical integral of a 3D scalar             !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(:), intent(in) :: x
  real, intent(out) :: xav
  real, optional, intent(in) :: zlim, zbas
  logical, intent(in), optional :: weighted
!
  integer  :: k
  real :: zz, weight(nz)
!                                                          !
!----------------------------------------------------------!
!
  if (present(weighted)) then
    if (weighted) then
      weight = den0
    else
      weight = 1.
    endif
  else
    weight = 1.
  endif
!
  xav = 0.0
!
  if (present(zlim)) then
    do k = 1, nz
      if (z0(k) <= zlim) xav = xav + weight(k)*x(k)*dz*fdz0(k)
    enddo
  else if (present(zbas)) then
    do k = 1, nz
      if (z0(k) > zbas) xav = xav + weight(k)*x(k)*dz*fdz0(k)
    enddo
  else
    do k = 1, nz
      xav = xav + weight(k)*x(k)*dz*fdz0(k)
    enddo
  endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine vertint_1D
!
! ==================================================
  subroutine max_all_3d ( x, xmax )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Returns max value of a 3D scalar  	                   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, intent(out) :: xmax
!
  integer :: ierr
  real :: xmax0
!                                                          !
!----------------------------------------------------------!
!
  xmax0 = maxval( x(it_start:it_end,jt_start:jt_end,1:nz) )
!
#ifdef SPMD
    CALL MPI_ALLREDUCE (xmax0, xmax, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
#else
    xmax = xmax0  
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine max_all_3d 
!
! ==================================================
  subroutine max_all_2d ( x, xmax )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Returns max value of a 2D scalar  	                   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:), intent(in) :: x
  real, intent(out) :: xmax
!
  integer :: ierr
  real :: xmax0
!                                                          !
!----------------------------------------------------------!
!
  xmax0 = maxval( x(it_start:it_end,jt_start:jt_end) )
!
#ifdef SPMD
    CALL MPI_ALLREDUCE (xmax0, xmax, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
#else
    xmax = xmax0  
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine max_all_2d 
!
! ==================================================
  subroutine av_all_3d ( x, xav, xq1, xc1, filter1 )
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates horizontal average of a 3D scalar            !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x
  real, dimension(ip_start:,jp_start:,:), intent(in), optional :: xq1
  real, intent(out) :: xav
  real, intent(in), optional :: xc1
  character(len=*), intent(in), optional :: filter1
!
  integer  :: i, j, k, ierr
  real :: crit, c, ii, iitot, xav0
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: f1
!                                                          !
!----------------------------------------------------------!
!
!  Calculate filter 1
!
  f1(ip_start:ip_end,jp_start:jp_end,1:nz) = 1.
  if (present(xq1)) then
    if (present(xc1)) then
      crit = xc1
    else
      crit = 0.
    endif
    if (present(filter1)) then
      if (trim(filter1) == 'lower') then
        where (xq1 > crit) f1 = 0. 
      else if (trim(filter1) == 'larger') then
        where (xq1 <= crit) f1 = 0. 
      endif
    else
      where (xq1 <= crit) f1 = 0. 
    endif
  endif
!
!  Cumulate horizontally with filter
!
  ii = 0.0
  xav0 = 0.0
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        c = f1(i,j,k)
  	xav0 = xav0 + c*x(i,j,k)
	ii = ii + c
      enddo
    enddo
  enddo
!  
#ifdef SPMD
  xav = xav0
  iitot = ii
  call MPI_ALLReduce (xav0, xav, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (ii, iitot, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  xav = xav/max(iitot,1.)
#else
  xav = xav0/max(ii,1.)
#endif
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine av_all_3d
!
! ==================================================
  subroutine cal_fluxx (x1,x2,xx)
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates scalar corrleations and variances		   !
!  Correlations are computed at grid cell centers	   !
!  (i.e. at scalar nodes)				   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  integer  :: i, j, k
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x1, x2, xx
  real, dimension(nz) :: x1av, x2av
!                                                          !
!----------------------------------------------------------!
!
  xx = 0.
!
  call horav  (x1, x1av)
  call horav  (x2, x2av)
!
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        xx(i,j,k) = (x1(i,j,k) - x1av(k)) * (x2(i,j,k) - x2av(k))
      enddo
    enddo
  enddo
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine cal_fluxx
!
! ==================================================
  subroutine cal_flux (r,w,x,wx)
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates vertical scalar fluxes			   !
!  Fluxes are computed ad cell centers			   !
!  (i.e. at scalar nodes)				   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  integer  :: i, j, k
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: w
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: r, x, wx
  real, dimension(nz) :: wav
  real, dimension(nz) :: xav
!                                                          !
!----------------------------------------------------------!
!
  wx = 0.
!
  call horav (w, wav)
  call horav (r*x, xav)
!
  do k = 1, nz-1
    do j = jt_start, jt_end
      do i = it_start, it_end
        wx(i,j,k) =  0.5*((w(i,j,k+1)-wav(k+1))+(w(i,j,k)-wav(k))) * (r(i,j,k)*x(i,j,k)-xav(k))
      enddo
    enddo
  enddo
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine cal_flux
!
! ==================================================
  subroutine cal_fluxu (r,w,u,wu)
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates vertical scalar fluxes			   !
!  Fluxes are computed ad cell centers			   !
!  (i.e. at scalar nodes)				   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  integer  :: i, j, k
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: w
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: r, u, m, wu
  real, dimension(nz) :: wav
  real, dimension(nz) :: mav
!                                                          !
!----------------------------------------------------------!
!
  m = 0.
  wu = 0.
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      m(i,j,:) = 0.5*(u(i+1,j,:)+u(i,j,:))*r(i,j,:)
    enddo
  enddo
!
  call horav (w, wav)
  call horav (m, mav)
!
  do k = 1, nz-1
    do j = jt_start, jt_end
      do i = it_start, it_end
        wu(i,j,k) =  0.5*((w(i,j,k+1)-wav(k+1))+(w(i,j,k)-wav(k))) * (m(i,j,k)-mav(k))
      enddo
    enddo
  enddo
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine cal_fluxu
!
! ==================================================
  subroutine cal_fluxv (r,w,v,wv)
! ==================================================
!
!----------------------------------------------------------!
!                                                          !
!  Calculates vertical scalar fluxes			   !
!  Fluxes are computed ad cell centers			   !
!  (i.e. at scalar nodes)				   !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  integer  :: i, j, k
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: w
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: r, v, m, wv
  real, dimension(nz) :: wav
  real, dimension(nz) :: mav
!                                                          !
!----------------------------------------------------------!
!
  m = 0.
  wv = 0.
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      m(i,j,:) = 0.5*(v(i,j+1,:)+v(i,j,:))*r(i,j,:)
    enddo
  enddo
!
  call horav (w, wav)
  call horav (m, mav)
!
  do k = 1, nz-1
    do j = jt_start, jt_end
      do i = it_start, it_end
        wv(i,j,k) =  0.5*((w(i,j,k+1)-wav(k+1))+(w(i,j,k)-wav(k))) * (m(i,j,k)-mav(k))
      enddo
    enddo
  enddo
!
  return
!                                                          !
!----------------------------------------------------------!
!
  end subroutine cal_fluxv

  end module averages
