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
!  netcdf_interface.f90                   
!
!  Purpose:
!      Module file containing routines to output netcdf files of various dimensions                 
!
!  Author
!      Julien Savre
!      MIM, Ludwig Maximilian Universitat, Munich
!
! ================================================================

module netcdfslice

!
!-----------------------------------------------------------!
!

#ifdef SPMD
  USE mpi
  USE mpicomm
#endif

  USE shared_all
  USE averages
  USE netcdfdef
  USE netcdf
  
  implicit none
        
  private
  
  integer :: NVAR, NRECS = 1
  integer :: NZO, NYO, NXO, NTO, NZBDO
  integer :: i1, i2, j1, j2, is, js
  character (len = *), parameter :: Z_NAME = "Z"
  character (len = *), parameter :: Z_BD_NAME = "Z_BD"
  character (len = *), parameter :: Y_NAME = "Y"
  character (len = *), parameter :: X_NAME = "X"
  character (len = *), parameter :: REC_NAME = "time"
  character (len = *), parameter :: UNITS = "units", INFO = "description"
  character (len = *), parameter :: Y_UNITS = "m"
  character (len = *), parameter :: X_UNITS = "m"
  character (len = *), parameter :: Z_UNITS = "m"
  character (len = *), parameter :: ZBD_UNITS = "m"
  character (len = *), parameter :: REC_UNITS = "s"
  
  save i1, i2, j1, j2, is, js, NXO, NYO, NZO
  
  public :: write_2d_slice, write_2d_surf, write_2d_variables, write_tmp_slice
  
  contains

! ===============================================
  subroutine write_2d_slice ( cnt, lnew, pig )
! ===============================================
  
  implicit none
  
  real, intent(in) :: cnt
  logical, intent(in) :: lnew
  logical, intent(in), optional :: pig 
 
  integer :: ncidx(nslicex), ncidy(nslicey), ncidz(nslicez)
  integer :: i, ndims, ndim, include_parents=0
  integer :: startx(3,nslicex), countx(2,nslicex), dimidsx(3,nslicex), varidsx(nvarout+4,nslicex)
  integer :: starty(3,nslicey), county(2,nslicey), dimidsy(3,nslicey), varidsy(nvarout+4,nslicey)
  integer :: startz(3,nslicez), countz(2,nslicez), dimidsz(3,nslicez), varidsz(nvarout+4,nslicez)
  logical :: ex, new
!
  character(len=10) :: name
  character(len=7) :: car7
  character (len = 100) :: FILE_NAME
!
  real, allocatable, dimension(:) :: yy, xx
 
  ndims = 2
!  
  call define_indices ( i1, i2, j1, j2, is, js, NXO, NYO, NZO, NZBDO)
!
  do i = 1, nslicez
    call open_2d_z_slice ( cnt, slicez_io(i), ncidz(i), varidsz(:,i), dimidsz(:,i), nvarout, startz(:,i), countz(:,i) )
  enddo
!
  do i = 1, nslicex
    call open_2d_x_slice ( cnt, slicex_io(i), ncidx(i), varidsx(:,i), dimidsx(:,i), startx(:,i), countx(:,i) )
  enddo
!
  do i = 1, nslicey
    call open_2d_y_slice ( cnt, slicey_io(i), ncidy(i), varidsy(:,i), dimidsy(:,i), starty(:,i), county(:,i) )
  enddo

  ! Write slices      
  call write_variables (ndims,ncidx,ncidy,ncidz,varidsx,varidsy,varidsz,startx,starty,startz,countx,county,countz,cnt)
  
  ! Close z slices.
  do i = 1, nslicez
    write(car7,'(i7)')floor(slicez_io(i))
    FILE_NAME = './OUTPUT/slice_z_'//trim(adjustl(car7))//'.nc'  
    if (present(pig)) FILE_NAME = trim(FILE_NAME) // '_pig'
    inquire(FILE=trim(FILE_NAME), EXIST=ex)
    new = ( lnew .or. .not.ex )
  
#ifndef PARALLEL_OUT
    if (mypid == 0) then
#endif
      call check( nf90_close(ncidz(i)), FILE_NAME )
#ifndef PARALLEL_OUT
    endif
#endif
  enddo

  ! Close x slices.
  do i = 1, nslicex
    write(car7,'(i7)')floor(slicex_io(i))
    FILE_NAME = './OUTPUT/slice_x_'//trim(adjustl(car7))//'.nc'  
    if (present(pig)) FILE_NAME = trim(FILE_NAME) // '_pig'
    inquire(FILE=trim(FILE_NAME), EXIST=ex)
    new = ( lnew .or. .not.ex )
  
#ifndef PARALLEL_OUT
    if (mypid == 0) then
#endif
      call check( nf90_close(ncidx(i)), FILE_NAME )
#ifndef PARALLEL_OUT
    endif
#endif
  enddo

  ! Close y slices.
  do i = 1, nslicey
    write(car7,'(i7)')floor(slicey_io(i))
    FILE_NAME = './OUTPUT/slice_y_'//trim(adjustl(car7))//'.nc'  
    if (present(pig)) FILE_NAME = trim(FILE_NAME) // '_pig'
    inquire(FILE=trim(FILE_NAME), EXIST=ex)
    new = ( lnew .or. .not.ex )
  
#ifndef PARALLEL_OUT
    if (mypid == 0) then
#endif
      call check( nf90_close(ncidy(i)), FILE_NAME )
#ifndef PARALLEL_OUT
    endif
#endif
  enddo
  
  end subroutine write_2d_slice

! ===============================================
  subroutine write_2d_surf ( cnt, lnew )
! ===============================================
  
  implicit none
  
  real, intent(in) :: cnt
  logical, intent(in) :: lnew
 
  integer :: ncidz, ndims, ndim, include_parents=0
  integer :: startz(3), countz(2), dimidsz(3), varidsz(nvarsurf+4)
  logical :: ex, new
!
  character(len=10) :: name
  character(len=7) :: car7
  character (len = 100) :: FILE_NAME
!
  ndims = 2
!  
  call define_indices ( i1, i2, j1, j2, is, js, NXO, NYO, NZO, NZBDO)
  call open_2d_z_slice ( cnt, 0., ncidz, varidsz(:), dimidsz(:), nvarsurf, startz(:), countz(:), surf=.true. )
  call write_2d_variables (ndims,ncidz,varidsz,startz,countz,cnt)

! Close z slices.
  FILE_NAME = './OUTPUT/slice_surf.nc'  
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( lnew .or. .not.ex )

#ifndef PARALLEL_OUT
  if (mypid == 0) then
#endif
    call check( nf90_close(ncidz), FILE_NAME )
#ifndef PARALLEL_OUT
  endif
#endif
  
  end subroutine write_2d_surf
  
  subroutine open_2d_z_slice (cnt,height,ncid,varids,dimids,nvars,start,count,surf)
  
  implicit none
  
  integer, intent(in) :: nvars
  real, intent(in) :: cnt
  real, intent(in) :: height
  logical, intent(in), optional :: surf
  
  integer :: ndims, ndim, ncid, include_parents=0
  integer :: y_dimid, x_dimid, t_dimid, ix, iy, rec, i, NVAR
  integer :: start(3), count(2), dimids(3), varids(nvars+4)
  logical :: ex, new, lsurf
!
  character(len=10) :: name
  character(len=7) :: car7
  character (len = 100) :: FILE_NAME
!
  real, allocatable, dimension(:) :: yy, xx

  lsurf = .false.
  if (present(surf)) then
    if (surf) lsurf = .true.
  endif

  ! Create the file.
  if (lsurf) then
    FILE_NAME = './OUTPUT/slice_surf.nc'  
  else
    write(car7,'(i7)')floor(height)
    FILE_NAME = './OUTPUT/slice_z_'//trim(adjustl(car7))//'.nc'  
  endif
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( (new_run.and.ntime==0) .or. .not.ex )
  
  allocate(xx(1:NXO), yy(1:NYO))
  
  do ix = 1, NXO
    xx(ix) = (ix + is - 1)*dx
  end do
#ifdef MODEL_3D
  do iy = 1, NYO
    yy(iy) = (iy + js - 1)*dy
  end do
#endif

#ifndef PARALLEL_OUT
  if (mypid == 0) then
#endif
    !
    ! New file: define all variables
    !
    ndims = 2
    NVAR = nvars 
    if (new) then
      call check( nf90_create(FILE_NAME, cmode=or(nf90_clobber, nf90_netcdf4), ncid=ncid), FILE_NAME )
        
      ! Define the dimensions. The record dimension is defined to have
      ! unlimited length - it can grow as needed. In this example it is
      ! the time dimension.
      call check( nf90_def_dim(ncid, X_NAME, NXO, x_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, Y_NAME, NYO, y_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )

      ! Define the coordinate variables. We will only define coordinate
      ! variables for iy and ix.  Ordinarily we would need to provide
      ! an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID (x_dimid) and
      ! similarly for (y_dimid).
      call check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, varids(1)), FILE_NAME )
      call check( nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, varids(2)), FILE_NAME )
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(3)), FILE_NAME )
      call check( nf90_def_var(ncid, 'Delta_t', NF90_REAL, t_dimid, varids(4)), FILE_NAME )

      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, varids(1), UNITS, X_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), UNITS, Y_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(3), UNITS, REC_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(4), UNITS, REC_UNITS), FILE_NAME )

      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
      dimids = (/ x_dimid, y_dimid, t_dimid /)
 
      if (lsurf) then 
        call add_2d_variables (ndims,ndims+2,ncid,dimids,varids)
      else
        call add_variables (ndims,ndims+2,ncid,dimids,varids)
      endif

      ! End define mode.
      call check( nf90_enddef(ncid), FILE_NAME )
  
      ! Write the coordinate variable data.   
      call check( nf90_put_var(ncid, varids(1), xx), FILE_NAME )
      call check( nf90_put_var(ncid, varids(2), yy), FILE_NAME )
      call check( nf90_put_var(ncid, varids(3), time), FILE_NAME )
      call check( nf90_put_var(ncid, varids(4), dt0), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      count = (/ NXO, NYO /)
      start = (/ 1, 1, 1 /)
    else
  
      !
      ! Existing file: appending new dataset
      !  
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
          
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
          
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
          
      call check( nf90_inquire_dimension(ncid, dimids(ndims+1), name=name, len=NTO), FILE_NAME )

      call check( nf90_put_var(ncid, varids(ndim), time, start=(/ NTO+1 /)), FILE_NAME )
      call check( nf90_put_var(ncid, varids(ndim+1), dt0, start=(/ NTO+1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      count = (/ NXO, NYO /)
      start = (/ 1, 1, NTO+1 /)

    endif
#ifndef PARALLEL_OUT
  endif
#endif

  deallocate( xx, yy )
  
  end subroutine open_2d_z_slice

  subroutine open_2d_x_slice (cnt,height,ncid,varids,dimids,start,count)
  
  implicit none
  
  real, intent(in) :: cnt
  real, intent(in) :: height
  
  integer :: ndims, ncid, ndim, include_parents=0
  integer :: y_dimid, z_dimid, t_dimid, iz, iy, rec, i, NVAR
  integer :: start(3), count(2), dimids(3), varids(nvarout+4)
  logical :: ex, new
!
  character(len=10) :: name
  character(len=7) :: car7
  character (len = 100) :: FILE_NAME
!
  real, allocatable, dimension(:) :: yy, zz

  ! Create the file.
  write(car7,'(i7)')floor(height)
  FILE_NAME = './OUTPUT/slice_x_'//trim(adjustl(car7))//'.nc'  
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( (new_run.and.ntime==0) .or. .not.ex )
  
  allocate(zz(1:NZO), yy(1:NYO))
  
#ifdef MODEL_3D
  do iy = 1, NYO
    yy(iy) = (iy + js - 1)*dy
  end do
#endif
  zz=0.5*dz*fdz0(1)
  do iz = 2, NZO
    zz(iz) = zz(iz-1) + 0.5*(fdz0(iz-1)+fdz0(iz))*dz
  end do
  
#ifndef PARALLEL_OUT
  if (mypid == 0) then
#endif
    !
    ! New file: define all variables
    !
    ndims = 2
    NVAR = nvarout
    if (new) then
      call check( nf90_create(FILE_NAME,cmode=or(nf90_clobber, nf90_netcdf4),ncid=ncid), FILE_NAME )
        
      ! Define the dimensions. The record dimension is defined to have
      ! unlimited length - it can grow as needed. In this example it is
      ! the time dimension.
      call check( nf90_def_dim(ncid, Y_NAME, NYO, y_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, Z_NAME, NZO, z_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )

      ! Define the coordinate variables. We will only define coordinate
      ! variables for iy and ix.  Ordinarily we would need to provide
      ! an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID (x_dimid) and
      ! similarly for (y_dimid).
      call check( nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, varids(1)), FILE_NAME )
      call check( nf90_def_var(ncid, Z_NAME, NF90_REAL, z_dimid, varids(2)), FILE_NAME )
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(3)), FILE_NAME )
      call check( nf90_def_var(ncid, 'Delta_t', NF90_REAL, t_dimid, varids(4)), FILE_NAME )

      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, varids(1), UNITS, Y_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), UNITS, Z_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(3), UNITS, REC_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(4), UNITS, REC_UNITS), FILE_NAME )

      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
      dimids = (/ y_dimid, z_dimid, t_dimid /)
  
      call add_variables (ndims,ndims+2,ncid,dimids,varids)

      ! End define mode.
      call check( nf90_enddef(ncid), FILE_NAME )
  
      ! Write the coordinate variable data.   
      call check( nf90_put_var(ncid, varids(1), yy), FILE_NAME )
      call check( nf90_put_var(ncid, varids(2), zz), FILE_NAME )
      call check( nf90_put_var(ncid, varids(3), time), FILE_NAME )
      call check( nf90_put_var(ncid, varids(4), dt0), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      count = (/ NYO, NZO /)
      start = (/ 1, 1, 1 /)
    else
  
      !
      ! Existing file: appending new dataset
      !  
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
          
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
          
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
          
      call check( nf90_inquire_dimension(ncid, dimids(ndims+1), name=name, len=NTO), FILE_NAME )

      call check( nf90_put_var(ncid, varids(ndim), time, start=(/ NTO+1 /)), FILE_NAME )
      call check( nf90_put_var(ncid, varids(ndim+1), dt0, start=(/ NTO+1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      count = (/ NYO, NZO /)
      start = (/ 1, 1, NTO+1 /)

    endif
#ifndef PARALLEL_OUT
  endif
#endif

  deallocate( zz, yy )
  
  end subroutine open_2d_x_slice

  subroutine open_2d_y_slice (cnt,height,ncid,varids,dimids,start,count)
  
  implicit none
  
  real, intent(in) :: cnt
  real, intent(in) :: height
  
  integer :: ndims, ncid, ndim, include_parents=0
  integer :: x_dimid, z_dimid, t_dimid, iz, ix, rec, i, NVAR
  integer :: start(3), count(2), dimids(3), varids(nvarout+4)
  logical :: ex, new
!
  character(len=10) :: name
  character(len=7) :: car7
  character (len = 100) :: FILE_NAME
!
  real, allocatable, dimension(:) :: xx, zz

  ! Create the file.
  write(car7,'(i7)')floor(height)
  FILE_NAME = './OUTPUT/slice_y_'//trim(adjustl(car7))//'.nc'  
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( (new_run.and.ntime==0) .or. .not.ex )
  
  allocate(zz(1:NZO), xx(1:NXO))
  
#ifdef MODEL_3D
  do ix = 1, NXO
    xx(ix) = (ix + is - 1)*dx
  end do
#endif
  zz=0.5*dz*fdz0(1)
  do iz = 2, NZO
    zz(iz) = zz(iz-1) + 0.5*(fdz0(iz-1)+fdz0(iz))*dz
  end do
  
#ifndef PARALLEL_OUT
  if (mypid == 0) then
#endif
    !
    ! New file: define all variables
    !
    ndims = 2
    NVAR = nvarout
    if (new) then
      call check( nf90_create(FILE_NAME, cmode=or(nf90_clobber, nf90_netcdf4),ncid=ncid), FILE_NAME )
        
      ! Define the dimensions. The record dimension is defined to have
      ! unlimited length - it can grow as needed. In this example it is
      ! the time dimension.
      call check( nf90_def_dim(ncid, X_NAME, NXO, x_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, Z_NAME, NZO, z_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )

      ! Define the coordinate variables. We will only define coordinate
      ! variables for iy and ix.  Ordinarily we would need to provide
      ! an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID (x_dimid) and
      ! similarly for (y_dimid).
      call check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, varids(1)), FILE_NAME )
      call check( nf90_def_var(ncid, Z_NAME, NF90_REAL, z_dimid, varids(2)), FILE_NAME )
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(3)), FILE_NAME )
      call check( nf90_def_var(ncid, 'Delta_t', NF90_REAL, t_dimid, varids(4)), FILE_NAME )

      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, varids(1), UNITS, X_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), UNITS, Z_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(3), UNITS, REC_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(4), UNITS, REC_UNITS), FILE_NAME )

      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
      dimids = (/ x_dimid, z_dimid, t_dimid /)
  
      call add_variables (ndims,ndims+2,ncid,dimids,varids)

      ! End define mode.
      call check( nf90_enddef(ncid), FILE_NAME )
  
      ! Write the coordinate variable data.   
      call check( nf90_put_var(ncid, varids(1), xx), FILE_NAME )
      call check( nf90_put_var(ncid, varids(2), zz), FILE_NAME )
      call check( nf90_put_var(ncid, varids(3), time), FILE_NAME )
      call check( nf90_put_var(ncid, varids(4), dt0), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      count = (/ NXO, NZO /)
      start = (/ 1, 1, 1 /)
    else
  
      !
      ! Existing file: appending new dataset
      !  
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
          
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
          
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
          
      call check( nf90_inquire_dimension(ncid, dimids(ndims+1), name=name, len=NTO), FILE_NAME )

      call check( nf90_put_var(ncid, varids(ndim), time, start=(/ NTO+1 /)), FILE_NAME )
      call check( nf90_put_var(ncid, varids(ndim+1), dt0, start=(/ NTO+1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      count = (/ NXO, NZO /)
      start = (/ 1, 1, NTO+1 /)

    endif
#ifndef PARALLEL_OUT
  endif
#endif

  deallocate( zz, xx )
  
  end subroutine open_2d_y_slice
  
  subroutine write_tmp_slice (flag,height,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)
  
  implicit none
  
  integer, intent(in) :: flag
  real, intent(in) :: height
  real, dimension(nx,ny), intent(in) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  
  real :: yy(1:NYO), xx(1:NXO)
  integer :: i, i1, i2, j1, j2, is, js
  integer :: ncid, ndim, include_parents=0
  integer :: y_dimid, x_dimid, t_dimid, ix, iy, rec, NVAR
#ifdef MODEL_3D
  integer :: start(3), count(2), dimids(3), varids(8+3)
#else
  integer :: start(2), count(1), dimids(2), varids(8+2)
#endif 
  logical :: ex, new
  logical, save :: lfirst=.true.
!
  character(len=10) :: name
  character(len=7) :: car7
  character(len = 100) :: FILE_NAME
  
  ! Create pretend data. If this wasn't an example program, we would
  ! have some real data to write, for example, model output.
  call define_indices ( i1, i2, j1, j2, is, js, NXO, NYO, NZO, NZBDO)
  
  do ix = 1, NXO
    xx(ix) = (ix + is - 1)*dx
  end do
#ifdef MODEL_3D
  do iy = 1, NYO
    yy(iy) = (iy + js - 1)*dy
  end do
#endif

  ! Create the file.
  write(car7,'(i7)')floor(height)
  if (flag == 0) then
    FILE_NAME = './OUTPUT/slice_clt_'//trim(adjustl(car7))//'.nc'  
  endif
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( (new_run.and.lfirst) .or. .not.ex )
  
  NVAR = 2
  if (mypid == 0) then
    !
    ! New file: define all variables
    !
    if (new) then
      call check( nf90_create(FILE_NAME, cmode=or(nf90_clobber, nf90_netcdf4), ncid=ncid), FILE_NAME )
  
      ! Define the dimensions. The record dimension is defined to have
      ! unlimited length - it can grow as needed. In this example it is
      ! the time dimension.
      call check( nf90_def_dim(ncid, X_NAME, NXO, x_dimid), FILE_NAME )
#ifdef MODEL_3D
      call check( nf90_def_dim(ncid, Y_NAME, NYO, y_dimid), FILE_NAME )
#endif
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )

      ! Define the coordinate variables. We will only define coordinate
      ! variables for iy and ix.  Ordinarily we would need to provide
      ! an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID (x_dimid) and
      ! similarly for (y_dimid).
      i = 1
      call check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, varids(1)), FILE_NAME )
#ifdef MODEL_3D
      i = i + 1
      call check( nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, varids(i)), FILE_NAME )
#endif
      i = i + 1
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(i)), FILE_NAME )

      ! Assign units attributes to coordinate variables.
      i = 1
      call check( nf90_put_att(ncid, varids(1), UNITS, X_UNITS), FILE_NAME )
#ifdef MODEL_3D
      i = i + 1
      call check( nf90_put_att(ncid, varids(i), UNITS, Y_UNITS), FILE_NAME )
#endif
      i = i + 1
      call check( nf90_put_att(ncid, varids(i), UNITS, REC_UNITS), FILE_NAME )

      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
#ifdef MODEL_3D
      nvar = 3
      dimids = (/ x_dimid, y_dimid, t_dimid /)
#else
      nvar = 2
      dimids = (/ x_dimid, t_dimid /)
#endif
  
      if (flag == 0) then
        call register_var ('wmo  ', nvar, ncid, dimids, varids)
        call register_var ('buo  ', nvar, ncid, dimids, varids)
        call register_var ('ql   ', nvar, ncid, dimids, varids)
        call register_var ('lwp  ', nvar, ncid, dimids, varids)
        call register_var ('ctop ', nvar, ncid, dimids, varids)
        call register_var ('bltop', nvar, ncid, dimids, varids)
        call register_var ('cbas ', nvar, ncid, dimids, varids)
        call register_var ('cid  ', nvar, ncid, dimids, varids)
        call register_var ('uid  ', nvar, ncid, dimids, varids)
      endif

    ! End define mode.
      call check( nf90_enddef(ncid), FILE_NAME )
  
    ! Write the coordinate variable data.   
      i = 1
      call check( nf90_put_var(ncid, varids(1), xx), FILE_NAME )
#ifdef MODEL_3D
      i = i + 1
      call check( nf90_put_var(ncid, varids(i), yy), FILE_NAME )
#endif
      i = i + 1
      call check( nf90_put_var(ncid, varids(i), time), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
#ifdef MODEL_3D
      count = (/ NXO, NYO /)
      start = (/ 1, 1, 1 /)
#else
      count = (/ NXO /)
      start = (/ 1, 1 /)
#endif
    else
  
#ifdef MODEL_3D
      i = 3
#else
      i = 2
#endif

      !
      ! Existing file: appending new dataset
      !  
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
          
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
          
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
          
      call check( nf90_inquire_dimension(ncid, dimids(i), name=name, len=NTO), FILE_NAME )

      call check( nf90_put_var(ncid, varids(ndim), time, start=(/ NTO+1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
#ifdef MODEL_3D
      count = (/ NXO, NYO /)
      start = (/ 1, 1, NTO+1 /)
#else
      count = (/ NXO /)
      start = (/ 1, NTO+1 /)
#endif
    endif
  endif
        
  ! Write the pretend data. 
  if (mypid == 0) then
    call check( nf90_put_var(ncid, varids(i+1), tmp1(i1:i2,j1:j2), start=start, count=count), FILE_NAME )
    call check( nf90_put_var(ncid, varids(i+2), tmp2(i1:i2,j1:j2), start=start, count=count), FILE_NAME )
    call check( nf90_put_var(ncid, varids(i+3), tmp3(i1:i2,j1:j2), start=start, count=count), FILE_NAME )
    call check( nf90_put_var(ncid, varids(i+4), tmp4(i1:i2,j1:j2), start=start, count=count), FILE_NAME )
    call check( nf90_put_var(ncid, varids(i+5), tmp5(i1:i2,j1:j2), start=start, count=count), FILE_NAME )
    call check( nf90_put_var(ncid, varids(i+6), tmp6(i1:i2,j1:j2), start=start, count=count), FILE_NAME )
    call check( nf90_put_var(ncid, varids(i+7), tmp7(i1:i2,j1:j2), start=start, count=count), FILE_NAME )
    call check( nf90_close(ncid), FILE_NAME )
  endif
  
  lfirst = .false.
  
  end subroutine write_tmp_slice
  
  subroutine write_variables (ndims,ncidx,ncidy,ncidz,varidsx,varidsy,varidsz,startx,starty,startz,countx,county,countz,cnt)
  
    integer, intent(in) :: ncidx(:),ncidy(:),ncidz(:),ndims
    integer, dimension(:,:), intent(in) :: startx, starty, startz, countx, county, countz
    integer, dimension(:,:), intent(in) :: varidsx, varidsy, varidsz
    real, intent(in) :: cnt
    
    integer :: nvar, k, h, dim2,i,j
    real :: height2
    character(len=1) :: car1
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: tmp1, tmp2
    real, dimension(1:nz) :: xav
  
    !nvar = ndims+2
    nvar = size(varidsz,DIM=1) - nvarout
  
    if (out_u) call write_var ('u    ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, wind%u)
#ifdef MODEL_3D
    if (out_v) call write_var ('v    ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, wind%v)
#endif
    if (out_w) call write_var ('w    ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, wind%w)
    if (out_p) then
      call write_var ('p    ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, pressure%p)
      do k = 1, nz
        tmp1(:,:,k) = p0(k) + p1(k)
      enddo
      call write_var ('ptot ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
    endif
    if (out_t) call write_var ('t    ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, thermo%T)
#ifdef ISENTROPIC
    if (out_pt) call write_var ('pt   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, state%es)
    if (out_mse) call write_var ('mse  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, thermo%mse)
#else
    if (out_pt) call write_var ('pt   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, thermo%pt)
    if (out_mse) call write_var ('mse  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, state%es)
#endif
    if (out_ptv) call write_var ('pte  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, thermo%ptv)
    if (out_rho) then
      call write_var ('rho  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, pressure%dens)
      do k = 1, nz
        tmp1(:,:,k) = pressure%dens(:,:,k) - den0(k)
      enddo
      call write_var ('drho ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
    endif
    if (out_qv) call write_var ('qv   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, thermo%qv)
    if (out_qt) call write_var ('qt   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, state%qt)
    if (out_z) call write_var ('z    ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%z)
    if (lmicro>0.and.out_qc) call write_var ('qc   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(drop)%q)
    if (lmicro>0.and.out_qr) call write_var ('qr   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(rain)%q)
    if (lmicro>1.and.out_qi) call write_var ('qi   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(ice)%q)
    if (lmicro>2.and.out_qg) call write_var ('qg   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(grau)%q)
    if (lmicro>2.and.out_qs) call write_var ('qs   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(snow)%q)
    if (lmicro>3.and.out_qh) call write_var ('qh   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(hail)%q)
    if (lmicro>0.and.out_nc) call write_var ('nc   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(drop)%n)
    if (lmicro>0.and.out_nr) call write_var ('nr   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(rain)%n)
    if (lmicro>1.and.out_ni) call write_var ('ni   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(ice)%n)
    if (lmicro>2.and.out_ng) call write_var ('ng   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(grau)%n)
    if (lmicro>2.and.out_ns) call write_var ('ns   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(snow)%n)
    if (lmicro>3.and.out_nh) call write_var ('nh   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(hail)%n)
    if (lmicro>0.and.out_dc) call write_var ('dc   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(drop)%d)
    if (lmicro>0.and.out_dr) call write_var ('dr   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(rain)%d)
    if (lmicro>0.and.out_vp) call write_var ('vp   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(rain)%vp)
    if (lmicro>1.and.out_di) call write_var ('di   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(ice)%d)
    if (lmicro>2.and.out_dg) call write_var ('dg   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(grau)%d)
    if (lmicro>2.and.out_ds) call write_var ('ds   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(snow)%d)
    if (lmicro>3.and.out_dh) call write_var ('dh   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(hail)%d)
    if (lmicro>3.and.out_wg) call write_var ('wg   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(grau)%w)
    if (lmicro>3.and.out_ws) call write_var ('ws   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(snow)%w)
    if (lmicro>3.and.out_wh) call write_var ('wh   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, hydromtr(hail)%w)
    if (lmicro>0.and.out_prec) then
      tmp1 = pressure%dens*hydromtr(rain)%q*hydromtr(rain)%vp
      call write_var ('prel ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
    endif
    if (lmicro>1.and.out_prec) then
      tmp1 = 0.
      do h = ice, nhydro
        tmp1 = tmp1 + pressure%dens*hydromtr(h)%q*hydromtr(h)%vp
      enddo
      call write_var ('prei ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
    endif
    if (out_ccn) call write_var ('ccn  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nuc%ccn)
    if (out_in) call write_var ('in   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nuc%in)
    if (out_k) call write_var ('k    ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu%fkv)
    if (out_tke) call write_var ('ke   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%kres)
    if (out_tke) call write_var ('ke2  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu%ksgs)
    if (out_tke) call write_var ('wp2  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%wvar)
    if (out_tke) call write_var ('wp3  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%wske)
    if (out_tke) call write_var ('s2   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%S2)
    if (out_tke) call write_var ('n2   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%N2)
    if (out_sat) call write_var ('rh   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, thermo%rh)
    if (out_sat) call write_var ('rhi  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, thermo%rhi)
    if (out_div) call write_var ('div  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%div)
    if (out_div) call write_var ('qcon ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%qconv)
    if (out_vort) call write_var ('vox  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%vortx)
    if (out_vort) call write_var ('voy  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%vorty)
    if (out_vort) call write_var ('voz  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%vortz)
    if (out_mf) call write_var ('mf   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%mf)
    if (out_mf) call write_var ('cc   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%cc)
    if (out_mf) call write_var ('mfm  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%mfm)
    if (out_mf) call write_var ('ccm  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%ccm)
    if (out_grad) call write_var ('gu1  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(1)%grad%u)
    if (out_grad) call write_var ('gu2  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(2)%grad%u)
    if (out_grad) call write_var ('gu3  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(3)%grad%u)
    if (out_grad) call write_var ('gv1  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(1)%grad%v)
    if (out_grad) call write_var ('gv2  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(2)%grad%v)
    if (out_grad) call write_var ('gv3  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(3)%grad%v)
    if (out_grad) call write_var ('gw1  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(1)%grad%w)
    if (out_grad) call write_var ('gw2  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(2)%grad%w)
    if (out_grad) call write_var ('gw3  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(3)%grad%w)
    if (out_grad) call write_var ('gt1  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(1)%grad%pt)
    if (out_grad) call write_var ('gt2  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(2)%grad%pt)
    if (out_grad) call write_var ('gt3  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(3)%grad%pt)
    if (out_grad) call write_var ('gq1  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(1)%grad%qt)
    if (out_grad) call write_var ('gq2  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(2)%grad%qt)
    if (out_grad) call write_var ('gq3  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(3)%grad%qt)
    if (out_ent) call write_var ('activ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%activ)
    if (out_ent) call write_var ('wen  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%went)
    if (out_ent) call write_var ('wen1 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%went1)
    if (out_ent) call write_var ('wen2 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%went2)
    if (out_ent) call write_var ('wen5 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%went5)
    if (out_ent) call write_var ('wen3 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%went3)
    if (out_ent) call write_var ('wen4 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%went4)
    if (out_ent) call write_var ('den  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%dent)
    if (out_ent) call write_var ('cen  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%cent)
    if (out_ent) call write_var ('wde  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%wdet)
    if (out_ent) call write_var ('wde1 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%wdet1)
    if (out_ent) call write_var ('wde2 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%wdet2)
    if (out_ent) call write_var ('wde5 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%wdet5)
    if (out_ent) call write_var ('wde3 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%wdet3)
    if (out_ent) call write_var ('wde4 ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%wdet4)
    if (out_ent) call write_var ('dde  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%ddet)
    if (out_ent) call write_var ('cde  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, entrain%cdet)
    if (out_var) then
#ifdef ISENTROPIC
      if (out_pt) call cal_fluxx (state%es, state%es, tmp1)  
#else
      if (out_pt) call cal_fluxx (thermo%pt, thermo%pt, tmp1)  
#endif
      if (out_pt) call write_var ('ptv  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
!
      if (out_qv) call cal_fluxx (thermo%qv, thermo%qv, tmp1)  
      if (out_qv) call write_var ('qvv  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
    endif
    if (out_buoy) call write_var ('buoy ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, thermo%buoy)
    if (out_beff) call write_var ('beff ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%beff)
    if (out_dp) call write_var ('dpb  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%pb)
    if (out_dp) call write_var ('dpd  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%pd)
    if (out_dp) call write_var ('dpnh ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diagnos%pnh)
    if (out_flut) then
      if (out_u) call cal_fluxu (pressure%dens, wind%w, wind%u, tmp1)
      if (out_u) call write_var ('uw   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
#ifdef MODEL_3D
      if (out_v) call cal_fluxv (pressure%dens, wind%w, wind%v, tmp1)
      if (out_v) call write_var ('vw   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
#endif
      if (out_buoy) call write_var ('bres ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%bfres)
#ifdef ISENTROPIC
      if (out_pt) call cal_flux (pressure%dens, wind%w, state%es, tmp1)
#else
      if (out_pt) call cal_flux (pressure%dens, wind%w, thermo%pt, tmp1)
#endif
      if (out_pt) call write_var ('ptw  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
      if (out_qt) call cal_flux (pressure%dens, wind%w, state%qt, tmp1)
      if (out_qt) call write_var ('qtw  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
      if (out_ptv) call cal_flux (pressure%dens, wind%w, thermo%ptv, tmp1)
      if (out_ptv) call write_var ('ptvw ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
      if (out_sca.and.nscal>0) call cal_flux (pressure%dens, wind%w, state%scal(:,:,:,1), tmp1)
      if (out_sca.and.nscal>0) call write_var ('s1tw ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, tmp1)
    endif
    if (out_fsgs) then
      if (out_u) call write_var ('usgs ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%ufsgs)
#ifdef MODEL_3D
      if (out_v) call write_var ('vsgs ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%vfsgs)
#endif
      if (out_buoy) call write_var ('bsgs ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%bfsgs)
      if (out_pt) call write_var ('ptsgs', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%ptfsgs)
      if (out_qt) call write_var ('qtsgs', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, turbu_diag%qtfsgs)
    endif
    if (out_sca) then
      do h = 1, nscal
    	write(car1,'(i1)') h
    	call write_var ('sca'//car1//' ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, state%scal(:,:,:,h))
      enddo
    endif
    if (out_dtnet) call write_var ('dtn  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, rad%dtnet)
#ifdef RAD_ENABLE
    if (out_dtnet) call write_var ('dtsw  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, rad%dtsw)
    if (out_dtnet) call write_var ('dtlw  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, rad%dtlw)

    if (out_frad) call write_var ('fnet ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, rad%frad)
    if (out_frad) call write_var ('fsw  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, rad%fluxs)
    if (out_frad) call write_var ('flw  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, rad%fluxir)   
#endif  
  
#ifdef AERO_ENABLE
    if (out_aero) then
      do h = 1, nmode
        call write_var ('an   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, aero3d(h)%n)
        call write_var ('am   ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, aero3d(h)%m)
        call write_var ('ama  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, aero3d(h)%ma)
      enddo
    endif
#endif    

#ifdef NUC_CNT
    if (out_in) then
      call write_var ('dif  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(1)%mode(2)%n)
      call write_var ('dff  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(1)%mode(3)%n)
      call write_var ('dtf  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(1)%mode(1)%n)
      call write_var ('cif  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(2)%mode(2)%n)
      call write_var ('cff  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(2)%mode(3)%n)
      call write_var ('ctf  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(2)%mode(1)%n)
      call write_var ('dic  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(1)%mode(2)%nc)
      call write_var ('dfc  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(1)%mode(3)%nc)
      call write_var ('dtc  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(1)%mode(1)%nc)
      call write_var ('cic  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(2)%mode(2)%nc)
      call write_var ('cfc  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(2)%mode(3)%nc)
      call write_var ('ctc  ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, nucin2(2)%mode(1)%nc)
    endif
#endif    
!
    if (ldiag) then
      do h = 1, ndiag
        if (out_diagu) call write_var ('udiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%u)
        if (out_diagv) call write_var ('vdiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%v)
        if (out_diagw) call write_var ('wdiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%w)
        if (out_diagp) call write_var ('pdiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%p)
        if (out_diagt) call write_var ('tdiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%pt)
        if (out_diagtv) call write_var ('tvdiag', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%ptv)
        if (out_diagq) call write_var ('qdiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%qt)
        if (out_diagk) call write_var ('kdiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%k)
        if (lmicro>0.and.out_diagl) call write_var ('ldiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%qc)
        if (lmicro>0.and.out_diagr) call write_var ('rdiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%qr)
        if (lmicro>1.and.out_diagi) call write_var ('idiag ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%qi)
        if (lmicro>0.and.out_micro) then
	  do k = 1, nhydro
            write(car1,'(i1)')k
	    call write_var ('qmic'//car1//' ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%micro(k)%q)
	    if (moments==2) call write_var ('nmic'//car1//' ', nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%micro(k)%n)
	  enddo
	endif
	do k = 1, nscal
          write(car1,'(i1)')k
          if (out_diags) call write_var ('sdiag'//car1, nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%sca(:,:,:,k))
        enddo
#ifdef AERO_ENABLE
	do k = 1, nmode
          write(car1,'(i1)')k
          if (out_diaga) call write_var ('ndiag'//car1, nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%aero(k)%n)
          if (out_diaga) call write_var ('mdiag'//car1, nvar, ncidx, ncidy, ncidz, varidsx, varidsy, varidsz, startx, starty, startz, countx, county, countz, cnt, diag(h)%aero(k)%m)
        enddo
#endif
      enddo
    endif
!
    return
    end subroutine
!    
    subroutine write_var (var,nvar,ncidx,ncidy,ncidz,varidsx,varidsy,varidsz,startx,starty,startz,countx,county,countz,cnt,tab)
    
      integer, dimension(:), intent(in) :: ncidx, ncidy, ncidz
      character(len=5), intent(in) :: var
      integer, intent(inout) :: nvar
      integer, dimension(:,:), intent(in) :: varidsx, varidsy, varidsz
      integer, dimension(:,:), intent(in) :: startx, starty, startz, countx, county, countz
      real, intent(in) :: cnt
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: tab
!      
      integer :: i, j, k, l, ierr
      real :: data_max, data_min, maxall, minall, xx, yy
!      
      real, allocatable, dimension(:,:,:)  :: tab3d
      real, allocatable, dimension(:,:) :: slice
!
      if (verbose > 1) call write_debug('Starting writing in slice variable '//trim(var)) 
!    
      nvar=nvar+1
!
      allocate(tab3d(1:nx,1:ny,1:nz))
!
#if ( defined SPMD )  
	call collect ( tab,  tab3d )
#else
	tab3d = tab
#endif
!
!  Write slice in x-Y
!
        do l = 1, nslicez
          if ( mypid == 0 ) then 
	    allocate(slice(1:nx,1:ny))
!
	    do k = 1, nz-1
	      if (z0(k) >= slicez_io(l)) exit
	    enddo
!
	    slice(1:nx,1:ny) = tab3d(1:nx,1:ny,k)
!
            if ( minmax ) write(7,9100) trim(var)//' ', maxval( slice(i1:i2,j1:j2) ), minval( slice(i1:i2,j1:j2) )
!
            call check( nf90_put_var(ncidz(l), varidsz(nvar,l), cnt*slice(i1:i2,j1:j2), start=startz(:,l), count=countz(:,l)), var )
!
	    deallocate(slice)
	  endif
	enddo
!
!  Write slice in Y-Z
!
        do l = 1, nslicex
          if ( mypid == 0 ) then 
	    allocate(slice(1:ny,1:nz))
!
	    xx = 0.
	    do i = 3, nx-3
	      if (xx >= slicex_io(l)) exit
	      xx = xx + dx
	    enddo
!
	    slice(1:ny,1:nz) = tab3d(i,1:ny,1:nz)
!  
            if ( minmax ) write(7,9100) trim(var)//' ', maxval( slice(j1:j2,1:nz) ), minval( slice(j1:j2,1:nz) )
!                    
            call check( nf90_put_var(ncidx(l), varidsx(nvar,l), cnt*slice(j1:j2,1:nz), start=startx(:,l), count=countx(:,l)), var )
!
	    deallocate(slice)
          endif
	enddo
!
!  Write slice in X-Z
!
        do l = 1, nslicey
          if ( mypid == 0 ) then 
	    allocate(slice(1:nx,1:nz))
!
	    yy = 0.
	    do j = 3, ny-3
	      if (yy >= slicey_io(l)) exit
	      yy = yy + dy
	    enddo
!
	    slice(1:nx,1:nz) = tab3d(1:nx,j,1:nz)
!  
            if ( minmax ) write(7,9100) trim(var)//' ', maxval( slice(i1:i2,1:nz) ), minval( slice(i1:i2,1:nz) )
!                    
            call check( nf90_put_var(ncidy(l), varidsy(nvar,l), cnt*slice(i1:i2,1:nz), start=starty(:,l), count=county(:,l)), var )
!
	    deallocate(slice)
          endif
	enddo
!
        deallocate(tab3d)
!
9100 format(11x,"Max & Min ",a5,"          ",f13.4,f16.4)
    
  end subroutine
  
  subroutine write_2d_variables (ndims,ncid,varids,start,count,cnt,hov)
  
    integer, intent(in) :: ncid, ndims
    integer, dimension(:), intent(in) :: start, count
    integer, dimension(:), intent(in) :: varids
    real, intent(in) :: cnt
    logical, intent(in), optional :: hov

    character(len=1) :: car1
    integer :: nvar,h
    logical :: lhov

    lhov = .false.
    if (present(hov)) then
      if (hov) lhov = .true.
    endif
!        
    !nvar = ndims+2
    nvar = size(varids) - nvarsurf
!
!  Must be called in the correct order in which they have been registered
!
    if (out_lwp) call write_2d_var ('lwp  ', nvar, ncid, varids, start, count, cnt, column%lwp, lhov)
    if (out_lwp) call write_2d_var ('iwp  ', nvar, ncid, varids, start, count, cnt, column%iwp, lhov)
    if (out_cwp) call write_2d_var ('cwp  ', nvar, ncid, varids, start, count, cnt, column%cwp, lhov)
    if (out_cwp) call write_2d_var ('rwp  ', nvar, ncid, varids, start, count, cnt, column%rwp, lhov)  
    if (out_wvp) call write_2d_var ('wvp  ', nvar, ncid, varids, start, count, cnt, column%wvp, lhov)
    if (out_wvp) call write_2d_var ('wvpt ', nvar, ncid, varids, start, count, cnt, column%wvpt, lhov)
    if (out_wvp) call write_2d_var ('wvpb ', nvar, ncid, varids, start, count, cnt, column%wvpb, lhov)
    if (out_cmse) call write_2d_var ('cmse ', nvar, ncid, varids, start, count, cnt, column%cmse, lhov)
    if (out_cmse) call write_2d_var ('cmset', nvar, ncid, varids, start, count, cnt, column%cmset, lhov)
    if (out_cmse) call write_2d_var ('cmseb', nvar, ncid, varids, start, count, cnt, column%cmseb, lhov)
    if (out_cmfl) call write_2d_var ('cmfl ', nvar, ncid, varids, start, count, cnt, column%cmfl, lhov)
    if (out_cape) call write_2d_var ('cap  ', nvar, ncid, varids, start, count, cnt, column%cape, lhov)
    if (out_cape) call write_2d_var ('cin  ', nvar, ncid, varids, start, count, cnt, column%cin, lhov)
    if (out_cape) call write_2d_var ('lcl  ', nvar, ncid, varids, start, count, cnt, column%lcl, lhov)
    if (out_cape) call write_2d_var ('lfc  ', nvar, ncid, varids, start, count, cnt, column%lfc, lhov)
    if (out_cp) call write_2d_var ('cpi0 ', nvar, ncid, varids, start, count, cnt, column%cpint0, lhov)
    if (out_cp) call write_2d_var ('cpi  ', nvar, ncid, varids, start, count, cnt, column%cpint, lhov)
    if (out_ctop) call write_2d_var ('ctop ', nvar, ncid, varids, start, count, cnt, column%ctop, lhov)
    if (out_ctop) call write_2d_var ('cbas ', nvar, ncid, varids, start, count, cnt, column%cbas, lhov)
    if (out_zinv) call write_2d_var ('bltop', nvar, ncid, varids, start, count, cnt, column%bltop, lhov)
    if (out_zinv) call write_2d_var ('zinv ', nvar, ncid, varids, start, count, cnt, column%zinv, lhov) 
    if (out_ints) call write_2d_var ('isca ', nvar, ncid, varids, start, count, cnt, column%intsca, lhov)
    if (out_rrate) call write_2d_var ('rrate', nvar, ncid, varids, start, count, cnt, surf%precip, lhov)
    if (out_rrate) call write_2d_var ('cumul', nvar, ncid, varids, start, count, cnt, surf%cumul, lhov)
    if (out_srad) call write_2d_var ('slw  ', nvar, ncid, varids, start, count, cnt, surf%lwf, lhov)
    if (out_srad) call write_2d_var ('ssw  ', nvar, ncid, varids, start, count, cnt, surf%swf, lhov)
    if (out_olr) call write_2d_var ('olr  ', nvar, ncid, varids, start, count, cnt, surf%olr, lhov)
    if (out_olr) call write_2d_var ('toa  ', nvar, ncid, varids, start, count, cnt, surf%toa, lhov)
    if (out_sfl) call write_2d_var ('shf  ', nvar, ncid, varids, start, count, cnt, thermo_prop%cp(:,:,1)*surf%esflux, lhov)
    if (out_sfl) call write_2d_var ('lhf  ', nvar, ncid, varids, start, count, cnt, flv00*surf%qvflux, lhov)
    if (out_sst) call write_2d_var ('sst  ', nvar, ncid, varids, start, count, cnt, state%es(:,:,1), lhov)
    if (out_sst) call write_2d_var ('ssq  ', nvar, ncid, varids, start, count, cnt, thermo%qv(:,:,1), lhov)
    if (out_osr) call write_2d_var ('osr  ', nvar, ncid, varids, start, count, cnt, surf%osr, lhov)
    if (out_osr) call write_2d_var ('isr  ', nvar, ncid, varids, start, count, cnt, surf%isr, lhov)
    if (out_ocs) call write_2d_var ('olrcs', nvar, ncid, varids, start, count, cnt, surf%olr_cs, lhov)
    if (out_ocs) call write_2d_var ('osrcs', nvar, ncid, varids, start, count, cnt, surf%osr_cs, lhov)

!#ifdef AERO_RADIA    
    if (out_opthic) call write_2d_var ('opta2', nvar, ncid, varids, start, count, cnt, rad%tau_aer_surf, lhov)  
    if (out_opthic) call write_2d_var ('optc2', nvar, ncid, varids, start, count, cnt, rad%tau_cloud_surf, lhov)
    if (out_opthic) call write_2d_var ('optr2', nvar, ncid, varids, start, count, cnt, rad%tau_rain_surf, lhov) 
!#endif   
!
    return
  end subroutine
!    
  subroutine write_2d_var (var,nvar,ncid,varids,start,count,cnt,tab,hov)
    
    character(len=5), intent(in) :: var
    integer, intent(inout) :: nvar
    integer, intent(in) :: ncid
    integer, dimension(:) :: varids
    integer, dimension(:), intent(in) :: start, count
    real, intent(in) :: cnt
    real, dimension(ip_start:ip_end,jp_start:jp_end), intent(in) :: tab
    logical, intent(in) :: hov

    integer :: j, ierr
    real :: data_max, data_min, maxall, minall
    real, allocatable, dimension(:,:) :: slice
    real, dimension(nx) :: sliceav
!    
!  Find and print out min max values
!
    nvar=nvar+1
!
#ifndef PARALLEL_OUT
    allocate(slice(1:nx,1:ny))
#else
    allocate(slice(ip_start:ip_end,jp_start:jp_end))
#endif
!
    if ( minmax ) then
#if ( defined SPMD )
      data_max  = maxval( tab(it_start:it_end,jt_start:jt_end) )
      data_min  = minval( tab(it_start:it_end,jt_start:jt_end) )
      CALL MPI_ALLREDUCE (data_max, maxall, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
      CALL MPI_ALLREDUCE (data_min, minall, 1, REALTYPE, MPI_MIN, MPI_COMM_WORLD, ierr) 
#else
      maxall = maxval( tab(it_start:it_end,jt_start:jt_end) )
      minall = minval( tab(it_start:it_end,jt_start:jt_end) )
#endif
!
      if (mypid == 0) write(7,9100)trim(var)//' ',maxall,minall
    endif
!       
!  Collect slice
!
#if ( defined SPMD ) 
    call collect ( tab, slice )
#else
    slice = tab
#endif
!
!  Write slice
!
      if (.not.hov) then
!    
        if ( mypid == 0 ) then 
          call check( nf90_put_var(ncid, varids(nvar), cnt*slice(i1:i2,j1:j2), start=start, count=count), var )
        endif
!
!  Write Hovmoller
!
      else
!
        if ( mypid == 0 ) then
          sliceav = 0.0
          do j = 4, ny-2
    	    sliceav = sliceav + slice(:,j)
          enddo
          sliceav = sliceav/real(ny-5)
          call check( nf90_put_var(ncid, varids(nvar), cnt*sliceav(1:nx), start=start, count=count), var )
        endif
!
      endif
!   
    deallocate (slice)
!
9100 format(11x,"Max & Min ",a5,"          ",f13.4,f16.4)
    
  end subroutine

  subroutine check(status,name)
    integer, intent (in) :: status
    character(len=*), intent (in) :: name
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      print *, 'AN ERROR WAS DETECTED IN NETCDF INTERFACE, FILE/VARIABLE: ', trim(name)
      print *, 'FORCING MIMICA TO STOP (WITHOUT FINALIZING MPI)'
      stop
    end if
  end subroutine check  
  
end module netcdfslice
