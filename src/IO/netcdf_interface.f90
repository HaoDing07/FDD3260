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

!>
!!
!! # netcdf_interface.f90                   
!!
!! ## Purpose
!! Module file containing routines to output netcdf files of various dimensions                 
!!
!! ## Authors
!! * Julien Savre 
!!  (MIM, Ludwig Maximilian Universitat, Munich)
!! * Matthias Brakebusch
!!  (ACES, Stockholm University)
!! * Matthias Schwarz 
!!  (MISU, Stockholm University)
!!
!! ## Resources
!! Links to netCDF fortran90 documentation
!! * HTML-Version https://www.unidata.ucar.edu/software/netcdf/docs-fortran/
!! * PDF-Version available here:
!!      * /mimicav5/MANUAL/V5/Fortran-netcdf-f90.pdf
!!      * PDF-Version https://www.nag.com/market/training/fortran-workshop/netcdf-f90.pdf

!!
module netcdfmod

!
!-----------------------------------------------------------!
!

#ifdef SPMD
  USE mpi
  USE mpicomm
#endif

  USE shared_all
#ifdef LAGRANGE
  USE shared_lagrange
#endif
  USE averages
  USE netcdfslice
  USE netcdfdef
  USE netcdf
  
  implicit none
        
  private
  
  integer :: NVAR                                              !< integer pointer to variable ID
  integer :: NRECS = 1
  integer :: NZO, NYO, NXO, NTO, NVTO                          !< Integers to store length of dimensions
  integer :: i1, i2, j1, j2, is, js                            !< Indizes
  character (len = *), parameter :: Z_NAME = "Z"
  character (len = *), parameter :: VT_NAME = "VT"
  character (len = *), parameter :: Z_BD_NAME = "Z_BD"
  character (len = *), parameter :: X_BD_NAME = "X_BD"
  character (len = *), parameter :: Y_BD_NAME = "Y_BD"
  character (len = *), parameter :: Y_NAME = "Y"
  character (len = *), parameter :: X_NAME = "X"
  character (len = *), parameter :: REC_NAME = "time"
  character (len = *), parameter :: UNITS = "units", INFO = "description"
  character (len = *), parameter :: Y_UNITS = "m"
  character (len = *), parameter :: Y_BD_UNITS = "m"
  character (len = *), parameter :: X_UNITS = "m"
  character (len = *), parameter :: X_BD_UNITS = "m"
  character (len = *), parameter :: Z_UNITS = "m"
  character (len = *), parameter :: Z_BD_UNITS = "m"
  character (len = *), parameter :: REC_UNITS = "s"
  character (len = *), parameter :: AXIS = "axis"
  character (len = *), parameter :: X_AXIS = "X"
  character (len = *), parameter :: Y_AXIS = "Y"
  character (len = *), parameter :: Z_AXIS = "Z"
  
  save i1, i2, j1, j2, is, js, NXO, NYO, NZO, NVTO
  
  public :: write_dim, write_2d, write_hov, write_yav, write_1d_profile, write_lagrange
  
  contains

  ! ================================================================
  !> Write netCDF output files with maximum dimensionality (3-run --> 3d+time, 2d-run --> 2d+time)
  !! subroutine create (first call) and updates (subsequent call) netcdf
  !! output files 
  subroutine write_dim (NDIMS, cnt, lnew, fail, pig)
  !
  integer, intent(in) :: NDIMS              !< number of dimensions (depends if 2d or 3d run)
  logical, intent(in) :: fail, lnew               !< flag
  logical, intent(in), optional :: pig      !< flag
  real, intent(in) :: cnt                   !< counter
  !
  character (len = 100) :: FILE_NAME        !< name of netcdf file
  character (len = 100) :: name             !< name of netcdf variable
  character(len=3) :: npid
  integer :: idx_out_start
  integer :: ncid                           !< netCDF ID
  integer :: ndim                           !< number of netcdf dimensions
  integer :: include_parents=0              !< flag for nf90_inq_dimids
  integer :: z_dimid, y_dimid, x_dimid, t_dimid !< dimension IDs 
  integer :: iz, iy, ix, rec, i             !< integer counters
  integer :: vt_dimid                       !< dimension ID for vertices-dimension 
  integer :: start(NDIMS+1)                 !< array to tell netCDF where to put the data
  integer :: count(NDIMS)                   !< array to tell netCDF where to put the data
  integer :: dimids4vars(NDIMS+1)           !< list of dimension IDs +1 for time +1 to create variables
  integer :: dimids(NDIMS+2)                !< list of dimension IDs +1 for time +1 for z_bd
  integer :: zbd_varid, xbd_varid, ybd_varid !< varids for grid bounds
  integer :: varids(nvarout+NDIMS+1+2)       !< number of output variables + number if dimension (3or2) + 1 for time + 1 for grid cell boundaries
  integer :: start2d(NDIMS), count2d(NDIMS-1), dimids2d(NDIMS) !< reduced dimension IDs for 2d (3d run) and 1d (2d run) output
  logical :: ex                             !< flag if file exists
  logical :: new                            !< flag if new run
  !
  real, dimension(:), allocatable :: yy, xx, zz !< 1d arrays for x, z, and z values
  real, dimension(:,:), allocatable :: zz_bd !< 2d arrays for z upper and lower grid cell boundary
  !
  ! Create the file(s).
  if (fail) then
    FILE_NAME = './OUTPUT/' // trim(file_output) // '_err'
  else
    FILE_NAME = './OUTPUT/' // trim(file_output)
  endif
  !
#ifdef PARALLEL_OUT
  write(npid,'(i3)') mypid
  FILE_NAME = trim(FILE_NAME) // '_' // trim(adjustl(npid))
#endif
  !
  if (present(pig)) FILE_NAME = trim(FILE_NAME) // '_pig'
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( lnew .or. .not.ex )
  !
  ! Dimenions.
  call define_indices ( i1, i2, j1, j2, is, js, NXO, NYO, NZO, NVTO )
  !
  ! Allocate arrays
  allocate(xx(1:NXO), yy(1:NYO), zz(1:NZO))
  allocate(zz_bd(1:NVTO, 1:NZO))
  !
  ! create x-grid
  do ix = 1, NXO
    xx(ix) = (ix + is - 1)*dx
  end do
#ifdef MODEL_3D
  ! create y-grid
  do iy = 1, NYO
    yy(iy) = (iy + js - 1)*dy
  end do
#endif
  ! create z-grid
  zz=0.5*dz*fdz0(1)
  do iz = 2, NZO
    zz(iz) = zz(iz-1) + 0.5*(fdz0(iz-1)+fdz0(iz))*dz
  end do
  !
  ! z-grid cell boundaries
  zz_bd(1,1) = 0
  zz_bd(2, 1) = zz_bd(1,1) + dz * fdz0(1)
  do iz = 2, NZO
    zz_bd(1, iz) = zz_bd(2, iz-1)
    zz_bd(2, iz) = zz_bd(1, iz) + dz * fdz0(iz)
  end do
  !

#ifndef PARALLEL_OUT
  if (mypid == 0) then
#endif
    !
    ! New file: define all variables
    !
    NVAR = nvarout
    !
    if (new) then
      ! Create a net netCDF dataset
      call check( nf90_create(FILE_NAME, cmode=or(nf90_clobber, nf90_netcdf4),ncid=ncid), FILE_NAME )
      !
      ! Define the dimensions.
      call check( nf90_def_dim(ncid, X_NAME, NXO, x_dimid), FILE_NAME )
#ifdef MODEL_3D
      call check( nf90_def_dim(ncid, Y_NAME, NYO, y_dimid), FILE_NAME )
#endif
      call check( nf90_def_dim(ncid, Z_NAME, NZO, z_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, VT_NAME, NVTO, vt_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )
      !
      ! Define the coordinate variables. 
      call check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, varids(1)), FILE_NAME )
      i=1
#ifdef MODEL_3D
      call check( nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, varids(2)), FILE_NAME )
      i=2
#endif
      call check( nf90_def_var(ncid, Z_NAME, NF90_REAL, z_dimid, varids(i+1)), FILE_NAME )
      ! Define z-bounds variable   
      call check( nf90_def_var(ncid, Z_BD_NAME, NF90_REAL, (/vt_dimid, z_dimid/), varids(i+2)), FILE_NAME )
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(i+3)), FILE_NAME )
      call check( nf90_def_var(ncid, 'Delta_t', NF90_REAL, t_dimid, varids(i+4)), FILE_NAME )
!
      call check( nf90_put_att(ncid, varids(1), UNITS, X_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(1), AXIS, X_AXIS), FILE_NAME )
      i=1
#ifdef MODEL_3D
      call check( nf90_put_att(ncid, varids(2), UNITS, Y_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), AXIS, Y_AXIS), FILE_NAME )
      i=2
#endif
      call check( nf90_put_att(ncid, varids(i+1), UNITS, Z_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(i+1), AXIS, Z_AXIS), FILE_NAME )
      !
      call check( nf90_put_att(ncid, varids(i+2), "bounds", Z_BD_NAME), FILE_NAME )
      call check( nf90_put_att(ncid, varids(i+2), UNITS, Z_BD_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(i+3), UNITS, REC_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(i+4), UNITS, REC_UNITS), FILE_NAME )
      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. All of the netCDF variables we are creating
      ! share the same three/four dimensions (x/(z)/z/t). 
      ! In Fortran, the unlimited dimension must come last on the list of dimids.
      !
#ifdef MODEL_3D
      dimids4vars = (/ x_dimid, y_dimid, z_dimid, t_dimid /)
      dimids2d = (/ x_dimid, y_dimid, t_dimid /)
#else
      dimids4vars = (/ x_dimid, z_dimid, t_dimid /)
      dimids2d = (/ x_dimid, t_dimid /)
#endif
!
      call add_variables (ndims,ndims+3,ncid,dimids4vars,varids)
      !
      ! End define mode
      call check( nf90_enddef(ncid), FILE_NAME )
      !
      ! Write the coordinate variable data.  
      call check( nf90_put_var(ncid, varids(1), xx), FILE_NAME )
      i=1
#ifdef MODEL_3D
      call check( nf90_put_var(ncid, varids(2), yy), FILE_NAME )
      i=2
#endif
      call check( nf90_put_var(ncid, varids(i+1), zz), FILE_NAME )
      call check( nf90_put_var(ncid, varids(i+2), zz_bd), FILE_NAME )
      !
      call check( nf90_put_var(ncid, varids(i+3), time), FILE_NAME )
      call check( nf90_put_var(ncid, varids(i+4), dt0), FILE_NAME )
      ! These settings tell netcdf to write first timestep of data.
      if (ndims == 3) then
        count = (/ NXO, NYO, NZO /)
        start = (/ 1, 1, 1, 1 /)
        count2d = (/ NXO, NYO /)
        start2d = (/ 1, 1, 1 /)
      else
        count = (/ NXO, NZO /)
        start = (/ 1, 1, 1 /)
        count2d = (/ NXO /)
        start2d = (/ 1, 1 /)
      endif
      !
    else  ! End create new file
      !
      ! Existing file: appending new dataset
      !
      !
      ! open netcdf file
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
      !
      ! get dimension IDs from the existing netCDF file
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
      ! find all variables
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
      !
      ! get infromation about netCDF TIME dimension
      call check( nf90_inquire_dimension(ncid, dimids(ndims + 1 + 1), name=name, len=NTO), FILE_NAME )
      !
      ! update time variable
      ! index of varids --> number of dimension + z-boundary dimension + time
      call check( nf90_put_var(ncid, varids(NDIMS + 1 + 1 ), time, start=(/ NTO+1 /)), FILE_NAME )
      call check( nf90_put_var(ncid, varids(NDIMS + 1 + 2 ), dt0, start=(/ NTO+1 /)), FILE_NAME )
      ! These settings tell netcdf to write first timestep of data.
      if (ndims == 3) then
        count = (/ NXO, NYO, NZO /)
        start = (/ 1, 1, 1, NTO+1 /)
        count2d = (/ NXO, NYO /)
        start2d = (/ 1, 1, NTO+1 /)
      else
        count = (/ NXO, NZO /)
        start = (/ 1, 1, NTO+1 /)
        count2d = (/ NXO /)
        start2d = (/ 1, NTO+1 /)
      endif
    endif
!
#ifndef PARALLEL_OUT
  endif
#endif
  !
  ! Write the output data to file. 
  call write_variables (ndims,ncid,varids,start,count,cnt)

  !
  ! Close the file. This causes netCDF to flush all buffers.
#ifndef PARALLEL_OUT
  if (mypid == 0) then
#endif
    call check( nf90_close(ncid), FILE_NAME )
#ifndef PARALLEL_OUT
  endif
#endif
  !
  deallocate(xx, yy, zz, zz_bd)
  !
  end subroutine write_dim
  !
  ! ================================================================
  !> Write netCDF output files with reduced dimentionality (domain-means at all z-levels + time)
  !! subroutine create (first call) and updates (subsequent call) netcdf
  !! output files 
  !!
  subroutine write_2d ( cnt, lnew, prof, pig )
  !
  implicit none
  !
  real, intent(in) :: cnt                   !< counter
  character(len=*), intent(in) :: prof      !< param string charactarizing the profile type
  logical, intent(in) :: lnew
  logical, intent(in), optional :: pig     
  !  
  integer :: ncid                           !< netCDF ID
  integer :: ndim                           !< number of netcdf dimensions
  integer :: include_parents=0              !< flag for nf90_inq_dimids
  integer :: z_dimid, t_dimid               !< dimension IDs 
  integer :: iz, rec, i, NVAR               !< integer counters
  integer :: vt_dimid                       !< dimension ID for vertices-dimension 
  integer :: start(2)                       !<
  integer :: count(1)                       !<
  integer :: dimids4vars(2)                 !< list of dimension IDs to create variables
  integer :: dimids(3)                      !< list of dimension IDs
  integer :: varids(nvarout+4)              !< number of output variables + number if dimension (1, z-only) + 1 for time + 1 for grid cell boundaries
  logical :: ex                             !< flag if file exists
  logical :: new                            !< flag if new run
  integer :: n_dim_var
  character(len=100) :: FILE_NAME           !< name of netcdf file
  character(len=10) :: name                 !< name of netcdf variable
  real, allocatable, dimension(:) :: zz     !< 1d arrays for z values
  real, dimension(:,:), allocatable :: zz_bd!< 2d arrays for z upper and lower grid cell boundary
  ! 
  NZO  = nz
  NVTO = 2
  NVAR = nvarout
  !
  ! Allocate arrays
  allocate( zz(1:NZO) )
  allocate(zz_bd(1:NVTO, 1:NZO))
  !
  ! create z grid
  zz=0.5*dz*fdz0(1)
  do iz = 2, NZO
    zz(iz) = zz(iz-1) + 0.5*(fdz0(iz-1)+fdz0(iz))*dz
  end do
  !
  ! z-boundaries
  zz_bd(1,1) = 0
  zz_bd(2, 1) = zz_bd(1,1) + dz * fdz0(1)
  do iz = 2, NZO
    zz_bd(1, iz) = zz_bd(2, iz-1)
    zz_bd(2, iz) = zz_bd(1, iz) + dz * fdz0(iz)
  end do
  !
  ! Create the file.
  FILE_NAME = './OUTPUT/profiles_'//trim(prof)//'.nc'
  if (present(pig)) FILE_NAME = trim(FILE_NAME) // '_pig'
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( lnew .or. .not.ex )
  !
  if (mypid == 0) then
    !
    ! New file: define all variables
    !
    if (new) then
      ! Create a net netCDF dataset
      call check( nf90_create(FILE_NAME, cmode=or(nf90_clobber, nf90_netcdf4), ncid=ncid), FILE_NAME )
      !
      ! Define the dimensions.
      call check( nf90_def_dim(ncid, Z_NAME, NZO, z_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, VT_NAME, NVTO, vt_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )

      ! Define the coordinate variables. 
      call check( nf90_def_var(ncid, Z_NAME, NF90_REAL, z_dimid, varids(1)), FILE_NAME )
      call check( nf90_def_var(ncid, Z_BD_NAME, NF90_REAL, (/vt_dimid, z_dimid/), varids(2)), FILE_NAME)
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(3)), FILE_NAME )
      call check( nf90_def_var(ncid, 'Delta_t', NF90_REAL, t_dimid, varids(4)), FILE_NAME )

      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, varids(1), "axis", "Z"), FILE_NAME )
      call check( nf90_put_att(ncid, varids(1), UNITS, Z_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), "bounds", Z_BD_NAME), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), UNITS, Z_BD_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(3), UNITS, REC_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(4), UNITS, REC_UNITS), FILE_NAME )
      !
      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
      dimids4vars = (/ z_dimid, t_dimid /)
      call add_variables (1,4,ncid,dimids4vars,varids)

      ! End define mode.
      call check( nf90_enddef(ncid), FILE_NAME )
 
      ! Write the coordinate variable data.   
      call check( nf90_put_var(ncid, varids(1), zz), FILE_NAME )
      call check( nf90_put_var(ncid, varids(2), zz_bd), FILE_NAME )
      call check( nf90_put_var(ncid, varids(3), time), FILE_NAME )
      call check( nf90_put_var(ncid, varids(4), dt0), FILE_NAME )
      !
      ! These settings tell netcdf to write first timestep of data.
      start = (/ 1, 1 /)
    else ! End create new file
      !
      ! New file: define all variables
      !
      ! open netcdf file
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
      !
      ! get dimension IDs from the existing netCDF file
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
      ! find all variables
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
      !
      ! get infromation about netCDF TIME dimension
      call check( nf90_inquire_dimension(ncid, dimids(3), name=name, len=NTO), FILE_NAME )
      !
      ! update time variable
      ! index of varids --> number of dimension + z-boundary dimension + time
      call check( nf90_put_var(ncid, varids(3), time, start=(/ NTO+1 /)), FILE_NAME )
      call check( nf90_put_var(ncid, varids(4), dt0, start=(/ NTO+1 /)), FILE_NAME )
      !
      ! These settings tell netcdf to write first timestep of data.
      start = (/ 1, NTO+1 /)
    endif
  endif
  !
  ! These settings tell netcdf to write first timestep of data.
  count = (/ NZO /)
  !
  ! Write the output data to file. 
  call write_variables (1,ncid,varids,start,count,cnt,average=.true.,filter=prof)
  !
  ! Close the file. This causes netCDF to flush all buffers.
  if (mypid == 0) call check( nf90_close(ncid), FILE_NAME )
  !
  deallocate( zz )
  !
  end subroutine write_2d
  
  !====================================================================
  subroutine write_hov ( cnt )
  !====================================================================
  
  implicit none

  real, intent(in) :: cnt
!  
  integer :: ncid, ndim, include_parents=0
  integer :: x_dimid, t_dimid, ix, rec, i, NVAR
  integer :: start(2), count(1), dimids(2), varids(nvarsurf+2)
  logical :: ex, new
!  
  character(len=100) :: FILE_NAME
  character(len=10) :: name
!  
  real, allocatable, dimension(:) :: xx
! 
  NXO  = nx
  NVAR = nvarsurf

  allocate( xx(1:NXO) )
    
  ! Create pretend data. If this wasn't an example program, we would
  ! have some real data to write, for example, model output.
  xx=0.
  do ix = 2, NXO
    xx(ix) = xx(ix-1) + dx
  end do

  ! Create the file.
  FILE_NAME = './OUTPUT/hovmoller.nc'
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( (new_run.and.ntime==0) .or. .not.ex )

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
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )

      ! Define the coordinate variables. We will only define coordinate
      ! variables for iy and ix.  Ordinarily we would need to provide
      ! an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID (x_dimid) and
      ! similarly for (y_dimid).
      call check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, varids(1)), FILE_NAME )
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(2)), FILE_NAME )

      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, varids(1), UNITS, X_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), UNITS, REC_UNITS), FILE_NAME )

      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
      dimids = (/ x_dimid, t_dimid /)
  
      call add_2d_variables (1,2,ncid,dimids,varids)

    ! End define mode.
      call check( nf90_enddef(ncid), FILE_NAME )
  
    ! Write the coordinate variable data.   
      call check( nf90_put_var(ncid, varids(1), xx), FILE_NAME )
      call check( nf90_put_var(ncid, varids(2), time), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      start = (/ 1, 1 /)
    else
  
      !
      ! Existing file: appending new dataset
      !  
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
          
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
          
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
          
      call check( nf90_inquire_dimension(ncid, dimids(2), name=name, len=NTO), FILE_NAME )

      call check( nf90_put_var(ncid, varids(ndim), time, start=(/ NTO+1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      start = (/ 1, NTO+1 /)
    endif
  endif
        
  ! These settings tell netcdf to write first timestep of data.
  count = (/ NXO /)

  ! Write the pretend data. 
  call write_2d_variables (1,ncid,varids,start,count,cnt,hov=.true.)
  
  ! Close the file. This causes netCDF to flush all buffers.
  if (mypid == 0) call check( nf90_close(ncid), FILE_NAME )
  
  deallocate( xx )
  
  end subroutine write_hov

  !
  !====================================================================
  subroutine write_yav ( cnt )
  
  implicit none

  real, intent(in) :: cnt
!  
  integer :: ncid, ndim, ndims, NVAR, include_parents=0
  integer :: x_dimid, z_dimid, t_dimid, ix, iz, rec, i
  integer :: start(3), count(2), dimids(3), varids(nvarout+3)
  logical :: ex, new
!  
  character(len=100) :: FILE_NAME
  character(len=10) :: name
!  
  real, allocatable, dimension(:) :: xx, zz
! 
  NXO  = nx
  NZO  = nz
  NVAR = nvarout
  ndims = 2

  ! Create the file.
  FILE_NAME = './OUTPUT/slice_yav.nc'
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( (new_run.and.ntime==0) .or. .not.ex )

  ! Dimenions.
  call define_indices ( i1, i2, j1, j2, is, js, NXO, NYO, NZO, NVTO )

  allocate( xx(1:NXO), zz(1:NZO) )
    
  ! Create pretend data. If this wasn't an example program, we would
  ! have some real data to write, for example, model output.
  xx=0.
  do ix = 2, NXO
    xx(ix) = xx(ix-1) + dx
  end do

  ! create z-grid
  zz=0.5*dz*fdz0(1)
  do iz = 2, NZO
    zz(iz) = zz(iz-1) + 0.5*(fdz0(iz-1)+fdz0(iz))*dz
  end do

  if (mypid == 0) then
    !
    ! New file: define all variables
    !
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
      start = (/ 1, 1, 1 /)
    else
  
      !
      ! Existing file: appending new dataset
      !  
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
          
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
          
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
          
      call check( nf90_inquire_dimension(ncid, dimids(ndims+1), name=name, len=NTO), FILE_NAME )

      call check( nf90_put_var(ncid, varids(ndims+1), time, start=(/ NTO+1 /)), FILE_NAME )
      call check( nf90_put_var(ncid, varids(ndims+2), dt0, start=(/ NTO+1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      start = (/ 1, 1, NTO+1 /)
    endif
  endif
        
  ! These settings tell netcdf to write first timestep of data.
  count = (/ NXO, NZO /)

  ! Write the pretend data. 
  call write_variables (ndims,ncid,varids,start,count,cnt,yav=.true.)
  
  ! Close the file. This causes netCDF to flush all buffers.
  if (mypid == 0) call check( nf90_close(ncid), FILE_NAME )
  
  deallocate( xx, zz )
  
  end subroutine write_yav

  subroutine write_lagrange 
!    
  implicit none
!  
  integer :: ncid, ndim, include_parents=0
  integer :: p_dimid, t_dimid, ix, rec, i, NVAR
  integer :: start(2), count(1), dimids(2)
  logical :: ex, new
!  
  character(len=100) :: FILE_NAME
  character(len=10) :: name
!  
  integer, allocatable, dimension(:) :: pp, varids
! 
#ifdef LAGRANGE
  NVAR = 10
  if (lmicro > 0) NVAR = NVAR + 1
  if (lmicro > 1) NVAR = NVAR + 1
  if (nscal > 0) NVAR = NVAR + nscal
  if (out_diagw) NVAR = NVAR + 9
  if (out_diagl) NVAR = NVAR + 9  
  
  allocate( pp(1:nparc), varids(1:NVAR+2) )
  do i = 1, nparc
    pp(i) = i
  enddo

  ! Create the file.
  FILE_NAME = './OUTPUT/lagrange.nc'
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( (new_run.and.ntime==0) .or. .not.ex )

  if (mypid == 0) then
    !
    ! New file: define all variables
    !
    if (new) then
      call check( nf90_create(FILE_NAME, cmode=or(nf90_clobber, nf90_netcdf4),ncid=ncid), FILE_NAME )
  
      ! Define the dimensions. The record dimension is defined to have
      ! unlimited length - it can grow as needed. In this example it is
      ! the time dimension.
      call check( nf90_def_dim(ncid, "PARC", nparc, p_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )

      ! Define the coordinate variables. We will only define coordinate
      ! variables for iy and ix.  Ordinarily we would need to provide
      ! an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID (x_dimid) and
      ! similarly for (y_dimid).
      call check( nf90_def_var(ncid, "PARC", NF90_REAL, p_dimid, varids(1)), FILE_NAME )
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(2)), FILE_NAME )

      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, varids(1), UNITS, "#"), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), UNITS, REC_UNITS), FILE_NAME )

      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
      dimids = (/ p_dimid, t_dimid /)
  
      call add_lag_variables (ncid,dimids,varids)

    ! End define mode.
      call check( nf90_enddef(ncid), FILE_NAME )
  
    ! Write the coordinate variable data.   
      call check( nf90_put_var(ncid, varids(1), pp), FILE_NAME )
      call check( nf90_put_var(ncid, varids(2), time), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      start = (/ 1, 1 /)
    else
  
      !
      ! Existing file: appending new dataset
      !  
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
          
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
          
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
          
      call check( nf90_inquire_dimension(ncid, dimids(2), name=name, len=NTO), FILE_NAME )

      call check( nf90_put_var(ncid, varids(2), time, start=(/ NTO+1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      start = (/ 1, NTO+1 /)
    endif
  endif
        
  ! These settings tell netcdf to write first timestep of data.
  count = (/ nparc /)

  ! Write the pretend data.
  call write_lag_variables (ncid,varids,start,count)
  
  ! Close the file. This causes netCDF to flush all buffers.
  if (mypid == 0) call check( nf90_close(ncid), FILE_NAME )
  
  deallocate( pp, varids )
  
#endif
 
  end subroutine write_lagrange
  
  subroutine write_1d_profile (flag,lnew,loc,plum)
  
  implicit none
  
  logical, intent(in) :: lnew
  integer, intent(in) :: flag
  type (filtered), dimension(nz), intent(in), optional :: loc
  type (plumes), dimension(nz), intent(in), optional :: plum
  
  real :: zz(1:nz)
  integer, allocatable ::  varids(:)
  integer :: ncid, ndim, include_parents=0
  integer :: z_dimid, t_dimid, iz, rec, NVAR
  integer :: count(1), start(2), dimids(2)
  logical :: ex, new
!
  character(len=10) :: name
  character(len=100) :: FILE_NAME
  
  ! Create pretend data. If this wasn't an example program, we would
  ! have some real data to write, for example, model output.

  ! Dimenions.
  zz=0.5*dz*fdz0(1)
  do iz = 2, nz
    zz(iz) = zz(iz-1) + 0.5*(fdz0(iz-1)+fdz0(iz))*dz
  end do

  ! Create files.
  if (flag == 1) then
    FILE_NAME = './OUTPUT/profiles_wenta.nc'
  else if (flag == 2) then
    FILE_NAME = './OUTPUT/profiles_wents.nc'
  else if (flag == 3) then
    FILE_NAME = './OUTPUT/profiles_wentd.nc'
  else if (flag == 4) then
    FILE_NAME = './OUTPUT/profiles_benta.nc'
  else if (flag == 5) then
    FILE_NAME = './OUTPUT/profiles_bents.nc'
  else if (flag == 6) then
    FILE_NAME = './OUTPUT/profiles_bentd.nc'
  else if (flag == 7) then
    FILE_NAME = './OUTPUT/profiles_plumewa.nc'
  else if (flag == 8) then
    FILE_NAME = './OUTPUT/profiles_plumews.nc'
  else if (flag == 9) then
    FILE_NAME = './OUTPUT/profiles_plumewd.nc'
  else if (flag == 10) then
    FILE_NAME = './OUTPUT/profiles_plumeba.nc'
  else if (flag == 11) then
    FILE_NAME = './OUTPUT/profiles_plumebs.nc'
  else if (flag == 12) then
    FILE_NAME = './OUTPUT/profiles_plumebd.nc'
  else
    return
  endif
  
  if (flag > 6) then
    allocate(varids(1:35))
  else
    allocate(varids(1:23))
  endif
  
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = ( (new_run.and.lnew) .or. .not.ex )
  
  NVAR = 2
  if (mypid == 0) then
    !
    ! New file: define all variables
    !
    if (new) then
      call check( nf90_create(FILE_NAME, cmode=or(nf90_clobber, nf90_netcdf4),ncid=ncid), FILE_NAME )
  
      ! Define the dimensions. The record dimension is defined to have
      ! unlimited length - it can grow as needed. In this example it is
      ! the time dimension.
      call check( nf90_def_dim(ncid, Z_NAME, nz, z_dimid), FILE_NAME )
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid), FILE_NAME )

      ! Define the coordinate variables. We will only define coordinate
      ! variables for iy and ix.  Ordinarily we would need to provide
      ! an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID (x_dimid) and
      ! similarly for (y_dimid).
      call check( nf90_def_var(ncid, Z_NAME, NF90_REAL, z_dimid, varids(1)), FILE_NAME )
      call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(2)), FILE_NAME )

      ! Assign units attributes to coordinate variables.
      call check( nf90_put_att(ncid, varids(1), UNITS, Z_UNITS), FILE_NAME )
      call check( nf90_put_att(ncid, varids(2), UNITS, REC_UNITS), FILE_NAME )

      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. In Fortran, the unlimited
      ! dimension must come last on the list of dimids.
      dimids = (/ z_dimid, t_dimid /)
  
      if (flag <= 6) then
        call register_var ('een  ', nvar, ncid, dimids, varids)
        call register_var ('qen  ', nvar, ncid, dimids, varids)
        call register_var ('ten  ', nvar, ncid, dimids, varids)
        call register_var ('ven  ', nvar, ncid, dimids, varids)
        call register_var ('buen ', nvar, ncid, dimids, varids)
        call register_var ('been ', nvar, ncid, dimids, varids)      
        call register_var ('wen  ', nvar, ncid, dimids, varids)      
        call register_var ('ede  ', nvar, ncid, dimids, varids)
        call register_var ('qde  ', nvar, ncid, dimids, varids)
        call register_var ('tde  ', nvar, ncid, dimids, varids)
        call register_var ('vde  ', nvar, ncid, dimids, varids)
        call register_var ('bude ', nvar, ncid, dimids, varids)
        call register_var ('bede ', nvar, ncid, dimids, varids)      
        call register_var ('wde  ', nvar, ncid, dimids, varids)      
        call register_var ('enb  ', nvar, ncid, dimids, varids)
        call register_var ('qtb  ', nvar, ncid, dimids, varids)
        call register_var ('tb   ', nvar, ncid, dimids, varids)
        call register_var ('teb  ', nvar, ncid, dimids, varids)
        call register_var ('bub  ', nvar, ncid, dimids, varids)
        call register_var ('beb  ', nvar, ncid, dimids, varids)      
        call register_var ('wb   ', nvar, ncid, dimids, varids)      
      else 
        call register_var ('qt   ', nvar, ncid, dimids, varids)
        call register_var ('pte  ', nvar, ncid, dimids, varids)
        call register_var ('w    ', nvar, ncid, dimids, varids)
        call register_var ('qte  ', nvar, ncid, dimids, varids)
        call register_var ('tee  ', nvar, ncid, dimids, varids)
        call register_var ('we   ', nvar, ncid, dimids, varids)
        call register_var ('qtb  ', nvar, ncid, dimids, varids)
        call register_var ('teb  ', nvar, ncid, dimids, varids)
        call register_var ('wb   ', nvar, ncid, dimids, varids)
        call register_var ('mf   ', nvar, ncid, dimids, varids)
        call register_var ('dmf  ', nvar, ncid, dimids, varids)
        call register_var ('cov  ', nvar, ncid, dimids, varids)
        call register_var ('dad  ', nvar, ncid, dimids, varids)
        call register_var ('mu   ', nvar, ncid, dimids, varids)
        call register_var ('eq   ', nvar, ncid, dimids, varids)
        call register_var ('et   ', nvar, ncid, dimids, varids)
        call register_var ('ew   ', nvar, ncid, dimids, varids)
        call register_var ('eqs  ', nvar, ncid, dimids, varids)
        call register_var ('ets  ', nvar, ncid, dimids, varids)
        call register_var ('ews  ', nvar, ncid, dimids, varids)
        call register_var ('alq  ', nvar, ncid, dimids, varids)
        call register_var ('alt  ', nvar, ncid, dimids, varids)
        call register_var ('alw  ', nvar, ncid, dimids, varids)
        call register_var ('alb  ', nvar, ncid, dimids, varids)
        call register_var ('alp  ', nvar, ncid, dimids, varids)
        call register_var ('bet  ', nvar, ncid, dimids, varids)
        call register_var ('gam  ', nvar, ncid, dimids, varids)
        call register_var ('dl1  ', nvar, ncid, dimids, varids)
        call register_var ('dl2  ', nvar, ncid, dimids, varids)
        call register_var ('wp2  ', nvar, ncid, dimids, varids)
        call register_var ('wp3  ', nvar, ncid, dimids, varids)
        call register_var ('ri   ', nvar, ncid, dimids, varids)
        call register_var ('mri  ', nvar, ncid, dimids, varids)
      endif

    ! End define mode.
      call check( nf90_enddef(ncid), FILE_NAME )
  
    ! Write the coordinate variable data.   
      call check( nf90_put_var(ncid, varids(1), zz), FILE_NAME )
      call check( nf90_put_var(ncid, varids(2), time, start=(/ 1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      count = (/ nz /)
      start = (/ 1, 1 /)
    else
  
      !
      ! Existing file: appending new dataset
      !  
      call check( nf90_open(FILE_NAME, nf90_write, ncid), FILE_NAME )
          
      call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents), FILE_NAME )
          
      call check( nf90_inq_varids(ncid, nvar, varids), FILE_NAME )
          
      call check( nf90_inquire_dimension(ncid, dimids(2), name=name, len=NTO), FILE_NAME )

      call check( nf90_put_var(ncid, varids(ndim), time, start=(/ NTO+1 /)), FILE_NAME )
      
      ! These settings tell netcdf to write first timestep of data.
      count = (/ nz /)
      start = (/ 1, NTO+1 /)
    endif

  endif
        
  ! Write the pretend data. 
  if (mypid == 0) then
    if ( flag <= 6 .and. present(loc) ) then
      call check( nf90_put_var(ncid, varids(nvar+1), loc%epse, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+2), loc%qte, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+3), loc%pte, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+4), loc%ptve, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+5), loc%buoye, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+6), loc%beffe, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+7), loc%we, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+8), loc%epsd, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+9), loc%qtd, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+10), loc%ptd, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+11), loc%ptvd, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+12), loc%buoyd, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+13), loc%beffd, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+14), loc%wd, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+15), loc%epsa, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+16), loc%qta, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+17), loc%pta, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+18), loc%ptva, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+19), loc%buoya, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+20), loc%beffa, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+21), loc%wa, start=start, count=count), FILE_NAME )
    else if ( present(plum) ) then
      call check( nf90_put_var(ncid, varids(nvar+1), plum%qt, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+2), plum%ptv, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+3), plum%w, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+4), plum%qte, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+5), plum%ptve, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+6), plum%we, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+7), plum%qtb, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+8), plum%ptvb, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+9), plum%wb, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+10), plum%mfl, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+11), plum%dmfl, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+12), plum%cov, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+13), plum%dadt, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+14), plum%mu, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+15), plum%eq, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+16), plum%et, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+17), plum%ew, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+18), plum%eqs, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+19), plum%ets, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+20), plum%ews, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+21), plum%alphaq, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+22), plum%alphat, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+23), plum%alphaw, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+24), plum%alphab, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+25), plum%alphap, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+26), plum%beta, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+27), plum%gamma, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+28), plum%delta1, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+29), plum%delta2, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+30), plum%wvar, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+31), plum%wske, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+32), plum%ri, start=start, count=count), FILE_NAME )
      call check( nf90_put_var(ncid, varids(nvar+33), plum%mri, start=start, count=count), FILE_NAME )
    endif
    call check( nf90_close(ncid), FILE_NAME )
  endif

  deallocate ( varids )

  end subroutine write_1d_profile
  !
  ! ================================================================
  !> Subroutine to prepare variables for output to netCDF
  !! 
  !!  
  subroutine write_variables (ndims,ncid,varids,start,count,cnt,average,yav,filter)
    !
    integer, intent(in) :: ncid,ndims
    integer, dimension(:), intent(in) :: start, count
    integer, dimension(:), intent(in) :: varids !< list of variable IDs of netCDF file
    !
    character(len=*), optional, intent(in) :: filter
    logical, optional, intent(in) :: average, yav    !< switch for averaging on/off
    real, intent(in) :: cnt
    !
    integer :: nvar !< a counter which points at the variable ID of a variable
    integer :: k, h, i,j
    logical :: cal_mean, av
    character(len=1) :: car1
    character(len=10) :: prof
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: tmp1, tmp2
    real, dimension(1:nz) :: xav
    !
    if(present(average)) then
      cal_mean=average
    else
      cal_mean=.false.
    endif
    !
    if(present(filter)) then
      prof=filter
    else
      prof='default'
    endif
    !
    if(present(yav)) then
      av=yav
    else
      av=.false.
    endif
    !
    ! get the ID of the fist variable
    nvar = size(varids) - nvarout
    !
    ! each write_var call adds +1 to nvar to loop though all variables IDs
    if (out_u) call write_var ('u    ', nvar, ncid, varids, start, count, cnt, wind%u, cal_mean, av, prof)
#ifdef MODEL_3D
    if (out_v) call write_var ('v    ', nvar, ncid, varids, start, count, cnt, wind%v, cal_mean, av, prof)
#endif
    if (out_w) call write_var ('w    ', nvar, ncid, varids, start, count, cnt, wind%w, cal_mean, av, prof)
    if (out_p) then
      call write_var ('p    ', nvar, ncid, varids, start, count, cnt, pressure%p, cal_mean, av, prof)
      do k = 1, nz
        tmp1(:,:,k) = p0(k) + p1(k)
      enddo
      call write_var ('ptot ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
    endif
    if (out_t) call write_var ('t    ', nvar, ncid, varids, start, count, cnt, thermo%T, cal_mean, av, prof)
#ifdef ISENTROPIC
    if (out_pt) call write_var ('pt   ', nvar, ncid, varids, start, count, cnt, state%es, cal_mean, av, prof)
    if (out_mse) call write_var ('mse  ', nvar, ncid, varids, start, count, cnt, thermo%mse, cal_mean, av, prof)
#else
    if (out_pt) call write_var ('pt   ', nvar, ncid, varids, start, count, cnt, thermo%pt, cal_mean, av, prof)
    if (out_mse) call write_var ('mse  ', nvar, ncid, varids, start, count, cnt, state%es, cal_mean, av, prof)
#endif
    if (out_ptv) call write_var ('pte  ', nvar, ncid, varids, start, count, cnt, thermo%ptv, cal_mean, av, prof)
    if (out_rho) then
      call write_var ('rho  ', nvar, ncid, varids, start, count, cnt, pressure%dens, cal_mean, av, prof)
      do k = 1, nz
        tmp1(:,:,k) = pressure%dens(:,:,k) - den0(k)
      enddo
      call write_var ('drho ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
    endif
    if (out_qv) call write_var ('qv   ', nvar, ncid, varids, start, count, cnt, thermo%qv, cal_mean, av, prof)
    if (out_qt) call write_var ('qt   ', nvar, ncid, varids, start, count, cnt, state%qt, cal_mean, av, prof)
    if (out_z) call write_var ('z    ', nvar, ncid, varids, start, count, cnt, diagnos%z, cal_mean, av, prof)
    if (lmicro>0.and.out_qc) call write_var ('qc   ', nvar, ncid, varids, start, count, cnt, hydromtr(drop)%q, cal_mean, av, prof)
    if (lmicro>0.and.out_qr) call write_var ('qr   ', nvar, ncid, varids, start, count, cnt, hydromtr(rain)%q, cal_mean, av, prof)
    if (lmicro>1.and.out_qi) call write_var ('qi   ', nvar, ncid, varids, start, count, cnt, hydromtr(ice)%q, cal_mean, av, prof)
    if (lmicro>2.and.out_qg) call write_var ('qg   ', nvar, ncid, varids, start, count, cnt, hydromtr(grau)%q, cal_mean, av, prof)
    if (lmicro>2.and.out_qs) call write_var ('qs   ', nvar, ncid, varids, start, count, cnt, hydromtr(snow)%q, cal_mean, av, prof)
    if (lmicro>3.and.out_qh) call write_var ('qh   ', nvar, ncid, varids, start, count, cnt, hydromtr(hail)%q, cal_mean, av, prof)
    if (lmicro>0.and.out_nc) call write_var ('nc   ', nvar, ncid, varids, start, count, cnt, hydromtr(drop)%n, cal_mean, av, prof)
    if (lmicro>0.and.out_nr) call write_var ('nr   ', nvar, ncid, varids, start, count, cnt, hydromtr(rain)%n, cal_mean, av, prof)
    if (lmicro>1.and.out_ni) call write_var ('ni   ', nvar, ncid, varids, start, count, cnt, hydromtr(ice)%n, cal_mean, av, prof)
    if (lmicro>2.and.out_ng) call write_var ('ng   ', nvar, ncid, varids, start, count, cnt, hydromtr(grau)%n, cal_mean, av, prof)
    if (lmicro>2.and.out_ns) call write_var ('ns   ', nvar, ncid, varids, start, count, cnt, hydromtr(snow)%n, cal_mean, av, prof)
    if (lmicro>3.and.out_nh) call write_var ('nh   ', nvar, ncid, varids, start, count, cnt, hydromtr(hail)%n, cal_mean, av, prof)
    if (lmicro>0.and.out_dc) call write_var ('dc   ', nvar, ncid, varids, start, count, cnt, hydromtr(drop)%d, cal_mean, av, prof)
    if (lmicro>0.and.out_dr) call write_var ('dr   ', nvar, ncid, varids, start, count, cnt, hydromtr(rain)%d, cal_mean, av, prof)
    ! if (lmicro>0.and.out_vp) call write_var ('vp   ', nvar, ncid, varids, start, count, cnt, hydromtr(rain)%vp, cal_mean, av, prof)
    if (lmicro>0.and.out_vp) call write_var ('vdrop', nvar, ncid, varids, start, count, cnt, hydromtr(drop)%vp, cal_mean, av, prof)
    if (lmicro>1.and.out_di) call write_var ('di   ', nvar, ncid, varids, start, count, cnt, hydromtr(ice)%d, cal_mean, av, prof)
    if (lmicro>2.and.out_dg) call write_var ('dg   ', nvar, ncid, varids, start, count, cnt, hydromtr(grau)%d, cal_mean, av, prof)
    if (lmicro>2.and.out_ds) call write_var ('ds   ', nvar, ncid, varids, start, count, cnt, hydromtr(snow)%d, cal_mean, av, prof)
    if (lmicro>3.and.out_dh) call write_var ('dh   ', nvar, ncid, varids, start, count, cnt, hydromtr(hail)%d, cal_mean, av, prof)
    if (lmicro>3.and.out_wg) call write_var ('wg   ', nvar, ncid, varids, start, count, cnt, hydromtr(grau)%w, cal_mean, av, prof)
    if (lmicro>3.and.out_ws) call write_var ('ws   ', nvar, ncid, varids, start, count, cnt, hydromtr(snow)%w, cal_mean, av, prof)
    if (lmicro>3.and.out_wh) call write_var ('wh   ', nvar, ncid, varids, start, count, cnt, hydromtr(hail)%w, cal_mean, av, prof)
    if (lmicro>0.and.out_prec) then
      tmp1 = pressure%dens*hydromtr(rain)%q*hydromtr(rain)%vp
      call write_var ('prel ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
    endif
    if (lmicro>1.and.out_prec) then
      tmp1 = 0.
      do h = ice, nhydro
        tmp1 = tmp1 + pressure%dens*hydromtr(h)%q*hydromtr(h)%vp
      enddo
      call write_var ('prei ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
    endif
    if (out_ccn) call write_var ('ccn  ', nvar, ncid, varids, start, count, cnt, nuc%ccn, cal_mean, av, prof)
    if (out_in) call write_var ('in   ', nvar, ncid, varids, start, count, cnt, nuc%in, cal_mean, av, prof)
    if (out_k) call write_var ('k    ', nvar, ncid, varids, start, count, cnt, turbu%fkv, cal_mean, av, prof)
    if (out_tke) call write_var ('ke   ', nvar, ncid, varids, start, count, cnt, turbu_diag%kres, cal_mean, av, prof)
    if (out_tke) call write_var ('ke2  ', nvar, ncid, varids, start, count, cnt, turbu%ksgs, cal_mean, av, prof)
    if (out_tke) call write_var ('wp2  ', nvar, ncid, varids, start, count, cnt, turbu_diag%wvar, cal_mean, av, prof)
    if (out_tke) call write_var ('wp3  ', nvar, ncid, varids, start, count, cnt, turbu_diag%wske, cal_mean, av, prof)
    if (out_tke) call write_var ('s2   ', nvar, ncid, varids, start, count, cnt, turbu_diag%S2, cal_mean, av, prof)
    if (out_tke) call write_var ('n2   ', nvar, ncid, varids, start, count, cnt, turbu_diag%N2, cal_mean, av, prof)
    if (out_sat) call write_var ('rh   ', nvar, ncid, varids, start, count, cnt, thermo%rh, cal_mean, av, prof)
    if (out_sat) call write_var ('rhi  ', nvar, ncid, varids, start, count, cnt, thermo%rhi, cal_mean, av, prof)
    if (out_div) call write_var ('div  ', nvar, ncid, varids, start, count, cnt, diagnos%div, cal_mean, av, prof)
    if (out_div) call write_var ('qcon ', nvar, ncid, varids, start, count, cnt, diagnos%qconv, cal_mean, av, prof)
    if (out_vort) call write_var ('vox  ', nvar, ncid, varids, start, count, cnt, diagnos%vortx, cal_mean, av, prof)
    if (out_vort) call write_var ('voy  ', nvar, ncid, varids, start, count, cnt, diagnos%vorty, cal_mean, av, prof)
    if (out_vort) call write_var ('voz  ', nvar, ncid, varids, start, count, cnt, diagnos%vortz, cal_mean, av, prof)
    if (out_mf) call write_var ('mf   ', nvar, ncid, varids, start, count, cnt, diagnos%mf, cal_mean, av, prof)
    if (out_mf) call write_var ('cc   ', nvar, ncid, varids, start, count, cnt, diagnos%cc, cal_mean, av, prof)
    if (out_mf) call write_var ('mfm  ', nvar, ncid, varids, start, count, cnt, diagnos%mfm, cal_mean, av, prof)
    if (out_mf) call write_var ('ccm  ', nvar, ncid, varids, start, count, cnt, diagnos%ccm, cal_mean, av, prof)
    if (out_grad) call write_var ('gu1  ', nvar, ncid, varids, start, count, cnt, diag(1)%grad%u, cal_mean, av, prof)
    if (out_grad) call write_var ('gu2  ', nvar, ncid, varids, start, count, cnt, diag(2)%grad%u, cal_mean, av, prof)
    if (out_grad) call write_var ('gu3  ', nvar, ncid, varids, start, count, cnt, diag(3)%grad%u, cal_mean, av, prof)
    if (out_grad) call write_var ('gv1  ', nvar, ncid, varids, start, count, cnt, diag(1)%grad%v, cal_mean, av, prof)
    if (out_grad) call write_var ('gv2  ', nvar, ncid, varids, start, count, cnt, diag(2)%grad%v, cal_mean, av, prof)
    if (out_grad) call write_var ('gv3  ', nvar, ncid, varids, start, count, cnt, diag(3)%grad%v, cal_mean, av, prof)
    if (out_grad) call write_var ('gw1  ', nvar, ncid, varids, start, count, cnt, diag(1)%grad%w, cal_mean, av, prof)
    if (out_grad) call write_var ('gw2  ', nvar, ncid, varids, start, count, cnt, diag(2)%grad%w, cal_mean, av, prof)
    if (out_grad) call write_var ('gw3  ', nvar, ncid, varids, start, count, cnt, diag(3)%grad%w, cal_mean, av, prof)
    if (out_grad) call write_var ('gt1  ', nvar, ncid, varids, start, count, cnt, diag(1)%grad%pt, cal_mean, av, prof)
    if (out_grad) call write_var ('gt2  ', nvar, ncid, varids, start, count, cnt, diag(2)%grad%pt, cal_mean, av, prof)
    if (out_grad) call write_var ('gt3  ', nvar, ncid, varids, start, count, cnt, diag(3)%grad%pt, cal_mean, av, prof)
    if (out_grad) call write_var ('gq1  ', nvar, ncid, varids, start, count, cnt, diag(1)%grad%qt, cal_mean, av, prof)
    if (out_grad) call write_var ('gq2  ', nvar, ncid, varids, start, count, cnt, diag(2)%grad%qt, cal_mean, av, prof)
    if (out_grad) call write_var ('gq3  ', nvar, ncid, varids, start, count, cnt, diag(3)%grad%qt, cal_mean, av, prof)
    if (out_ent) call write_var ('activ', nvar, ncid, varids, start, count, cnt, entrain%activ, cal_mean, av, prof)
    if (out_ent) call write_var ('wen  ', nvar, ncid, varids, start, count, cnt, entrain%went, cal_mean, av, prof)
    if (out_ent) call write_var ('wen1 ', nvar, ncid, varids, start, count, cnt, entrain%went1, cal_mean, av, prof)
    if (out_ent) call write_var ('wen2 ', nvar, ncid, varids, start, count, cnt, entrain%went2, cal_mean, av, prof)
    if (out_ent) call write_var ('wen5 ', nvar, ncid, varids, start, count, cnt, entrain%went5, cal_mean, av, prof)
    if (out_ent) call write_var ('wen3 ', nvar, ncid, varids, start, count, cnt, entrain%went3, cal_mean, av, prof)
    if (out_ent) call write_var ('wen4 ', nvar, ncid, varids, start, count, cnt, entrain%went4, cal_mean, av, prof)
    if (out_ent) call write_var ('den  ', nvar, ncid, varids, start, count, cnt, entrain%dent, cal_mean, av, prof)
    if (out_ent) call write_var ('cen  ', nvar, ncid, varids, start, count, cnt, entrain%cent, cal_mean, av, prof)
    if (out_ent) call write_var ('wde  ', nvar, ncid, varids, start, count, cnt, entrain%wdet, cal_mean, av, prof)
    if (out_ent) call write_var ('wde1 ', nvar, ncid, varids, start, count, cnt, entrain%wdet1, cal_mean, av, prof)
    if (out_ent) call write_var ('wde2 ', nvar, ncid, varids, start, count, cnt, entrain%wdet2, cal_mean, av, prof)
    if (out_ent) call write_var ('wde5 ', nvar, ncid, varids, start, count, cnt, entrain%wdet5, cal_mean, av, prof)
    if (out_ent) call write_var ('wde3 ', nvar, ncid, varids, start, count, cnt, entrain%wdet3, cal_mean, av, prof)
    if (out_ent) call write_var ('wde4 ', nvar, ncid, varids, start, count, cnt, entrain%wdet4, cal_mean, av, prof)
    if (out_ent) call write_var ('dde  ', nvar, ncid, varids, start, count, cnt, entrain%ddet, cal_mean, av, prof)
    if (out_ent) call write_var ('cde  ', nvar, ncid, varids, start, count, cnt, entrain%cdet, cal_mean, av, prof)
    if (out_var) then
#ifdef ISENTROPIC
      if (out_pt) call cal_fluxx (state%es, state%es, tmp1)  
#else
      if (out_pt) call cal_fluxx (thermo%pt, thermo%pt, tmp1)  
#endif
      if (out_pt) call write_var ('ptv  ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
!
      if (out_qv) call cal_fluxx (thermo%qv, thermo%qv, tmp1)  
      if (out_qv) call write_var ('qvv  ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
    endif
    if (out_buoy) call write_var ('buoy ', nvar, ncid, varids, start, count, cnt, thermo%buoy, cal_mean, av, prof)
    if (out_beff) call write_var ('beff ', nvar, ncid, varids, start, count, cnt, diagnos%beff, cal_mean, av, prof)
    if (out_dp) call write_var ('dpb  ', nvar, ncid, varids, start, count, cnt, diagnos%pb, cal_mean, av, prof)
    if (out_dp) call write_var ('dpd  ', nvar, ncid, varids, start, count, cnt, diagnos%pd, cal_mean, av, prof)
    if (out_dp) call write_var ('dpnh ', nvar, ncid, varids, start, count, cnt, diagnos%pnh, cal_mean, av, prof)
    if (out_flut) then
      if (out_u) call cal_fluxu (pressure%dens, wind%w, wind%u, tmp1) 
      if (out_u) call write_var ('uw   ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
#ifdef MODEL_3D
      if (out_v) call cal_fluxv (pressure%dens, wind%w, wind%v, tmp1) 
      if (out_v) call write_var ('vw   ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
#endif
      if (out_buoy) call write_var ('bres ', nvar, ncid, varids, start, count, cnt, turbu_diag%bfres, cal_mean, av, prof)
#ifdef ISENTROPIC
      if (out_pt) call cal_flux (pressure%dens, wind%w, state%es, tmp1)
#else
      if (out_pt) call cal_flux (pressure%dens, wind%w, thermo%pt, tmp1)
#endif
      if (out_pt) call write_var ('ptw  ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
      if (out_qt) call cal_flux (pressure%dens, wind%w, state%qt, tmp1) 
      if (out_qt) call write_var ('qtw  ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
      if (out_ptv) call cal_flux (pressure%dens, wind%w, thermo%ptv, tmp1)
      if (out_ptv) call write_var ('ptvw ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
      if (out_sca.and.nscal>0) call cal_flux (pressure%dens, wind%w, state%scal(:,:,:,1), tmp1)
      if (out_sca.and.nscal>0) call write_var ('s1tw ', nvar, ncid, varids, start, count, cnt, tmp1, cal_mean, av, prof)
    endif
    if (out_fsgs) then
      if (out_u) call write_var ('usgs ', nvar, ncid, varids, start, count, cnt, turbu_diag%ufsgs, cal_mean, av, prof)
#ifdef MODEL_3D
      if (out_v) call write_var ('vsgs ', nvar, ncid, varids, start, count, cnt, turbu_diag%vfsgs, cal_mean, av, prof)
#endif
      if (out_buoy) call write_var ('bsgs ', nvar, ncid, varids, start, count, cnt, turbu_diag%bfsgs, cal_mean, av, prof)
      if (out_pt) call write_var ('ptsgs', nvar, ncid, varids, start, count, cnt, turbu_diag%ptfsgs, cal_mean, av, prof)
      if (out_qt) call write_var ('qtsgs', nvar, ncid, varids, start, count, cnt, turbu_diag%qtfsgs, cal_mean, av, prof)
    endif
    if (out_sca) then
      do h = 1, nscal
    	write(car1,'(i1)') h
    	call write_var ('sca'//car1//' ', nvar, ncid, varids, start, count, cnt, state%scal(:,:,:,h), cal_mean, av, prof)
      enddo
    endif
    if (out_dtnet) call write_var ('dtn  ', nvar, ncid, varids, start, count, cnt, rad%dtnet, cal_mean, av, prof)
#ifdef RAD_ENABLE
    if (out_dtnet) call write_var ('dtsw  ', nvar, ncid, varids, start, count, cnt, rad%dtsw, cal_mean, av, prof)
    if (out_dtnet) call write_var ('dtlw  ', nvar, ncid, varids, start, count, cnt, rad%dtlw, cal_mean, av, prof)
 
    if (out_frad) call write_var ('fnet ', nvar, ncid, varids, start, count, cnt, rad%frad, cal_mean, av, prof)
    if (out_frad) call write_var ('fsw  ', nvar, ncid, varids, start, count, cnt, rad%fluxs, cal_mean, av, prof)
    if (out_frad) call write_var ('flw  ', nvar, ncid, varids, start, count, cnt, rad%fluxir, cal_mean, av, prof)  
#endif 

#ifdef AERO_ENABLE
    if (out_aero) then
      do h = 1, nmode
        call write_var ('an   ', nvar, ncid, varids, start, count, cnt, aero3d(h)%n, cal_mean, av, prof)
        call write_var ('am   ', nvar, ncid, varids, start, count, cnt, aero3d(h)%m, cal_mean, av, prof)
        call write_var ('ama  ', nvar, ncid, varids, start, count, cnt, aero3d(h)%ma, cal_mean, av, prof)
      enddo
    endif
#endif  

#ifdef NUC_CNT
    if (out_in) then
      call write_var ('dif  ', nvar, ncid, varids, start, count, cnt, nucin2(1)%mode(2)%n, cal_mean, av, prof)
      call write_var ('dff  ', nvar, ncid, varids, start, count, cnt, nucin2(1)%mode(3)%n, cal_mean, av, prof)
      call write_var ('dtf  ', nvar, ncid, varids, start, count, cnt, nucin2(1)%mode(1)%n, cal_mean, av, prof)
      call write_var ('cif  ', nvar, ncid, varids, start, count, cnt, nucin2(2)%mode(2)%n, cal_mean, av, prof)
      call write_var ('cff  ', nvar, ncid, varids, start, count, cnt, nucin2(2)%mode(3)%n, cal_mean, av, prof)
      call write_var ('ctf  ', nvar, ncid, varids, start, count, cnt, nucin2(2)%mode(1)%n, cal_mean, av, prof)
      call write_var ('dic  ', nvar, ncid, varids, start, count, cnt, nucin2(1)%mode(2)%nc, cal_mean, av, prof)
      call write_var ('dfc  ', nvar, ncid, varids, start, count, cnt, nucin2(1)%mode(3)%nc, cal_mean, av, prof)
      call write_var ('dtc  ', nvar, ncid, varids, start, count, cnt, nucin2(1)%mode(1)%nc, cal_mean, av, prof)
      call write_var ('cic  ', nvar, ncid, varids, start, count, cnt, nucin2(2)%mode(2)%nc, cal_mean, av, prof)
      call write_var ('cfc  ', nvar, ncid, varids, start, count, cnt, nucin2(2)%mode(3)%nc, cal_mean, av, prof)
      call write_var ('ctc  ', nvar, ncid, varids, start, count, cnt, nucin2(2)%mode(1)%nc, cal_mean, av, prof)
    endif
#endif    
!
    if (ldiag) then
      do h = 1, ndiag
        if (out_diagu) call write_var ('udiag ', nvar, ncid, varids, start, count, cnt, diag(h)%u, cal_mean, av, prof)
        if (out_diagv) call write_var ('vdiag ', nvar, ncid, varids, start, count, cnt, diag(h)%v, cal_mean, av, prof)
        if (out_diagw) call write_var ('wdiag ', nvar, ncid, varids, start, count, cnt, diag(h)%w, cal_mean, av, prof)
        if (out_diagp) call write_var ('pdiag ', nvar, ncid, varids, start, count, cnt, diag(h)%p, cal_mean, av, prof)
        if (out_diagt) call write_var ('tdiag ', nvar, ncid, varids, start, count, cnt, diag(h)%pt, cal_mean, av, prof)
        if (out_diagtv) call write_var ('tvdiag', nvar, ncid, varids, start, count, cnt, diag(h)%ptv, cal_mean, av, prof)
        if (out_diagq) call write_var ('qdiag ', nvar, ncid, varids, start, count, cnt, diag(h)%qt, cal_mean, av, prof)
        if (out_diagk) call write_var ('kdiag ', nvar, ncid, varids, start, count, cnt, diag(h)%k, cal_mean, av, prof)
        if (lmicro>0.and.out_diagl) call write_var ('ldiag ', nvar, ncid, varids, start, count, cnt, diag(h)%qc, cal_mean, av, prof)
        if (lmicro>0.and.out_diagr) call write_var ('rdiag ', nvar, ncid, varids, start, count, cnt, diag(h)%qr, cal_mean, av, prof)
        if (lmicro>1.and.out_diagi) call write_var ('idiag ', nvar, ncid, varids, start, count, cnt, diag(h)%qi, cal_mean, av, prof)
        if (lmicro>0.and.out_micro) then
	  do k = 1, nhydro
            write(car1,'(i1)')k
	    call write_var ('qmic'//car1//' ', nvar, ncid, varids, start, count, cnt, diag(h)%micro(k)%q, cal_mean, av, prof)
	    if (moments==2) call write_var ('nmic'//car1//' ', nvar, ncid, varids, start, count, cnt, diag(h)%micro(k)%n, cal_mean, av, prof)
	  enddo
	endif
	do k = 1, nscal
          write(car1,'(i1)')k
          if (out_diags) call write_var ('sdiag'//car1, nvar, ncid, varids, start, count, cnt, diag(h)%sca(:,:,:,k), cal_mean, av, prof)
        enddo
#ifdef AERO_ENABLE
        if (out_diaga) then
	  do k = 1, nmode
            write(car1,'(i1)')k
	    call write_var ('ndiag'//car1, nvar, ncid, varids, start, count, cnt, diag(h)%aero(k)%n, cal_mean, av, prof)
	    call write_var ('mdiag'//car1, nvar, ncid, varids, start, count, cnt, diag(h)%aero(k)%m, cal_mean, av, prof)
	  enddo
	endif
#endif
      enddo
    endif
!
    return
    end subroutine

#ifdef LAGRANGE  
  subroutine write_lag_variables (ncid,varids,start,count)
  
    integer :: nvar, nv, nr, h
    integer, intent(in) :: ncid
    integer, dimension(:), intent(in) :: start, count
    integer, dimension(:), intent(in) :: varids
!
    real, dimension(:,:), allocatable :: tab
    character(len=1) :: car1
!        
    nvar = 2
    nv = 11 + nscal
    if (lmicro > 0) nv = nv + 1
    if (lmicro > 1) nv = nv + 1
    if (out_diagw) nv = nv + 9
    if (out_diagl) nv = nv + 9
!  
!  Allocate and fill array of parcels
!
    allocate( tab(1:nparc,1:nv) )
!
    call gather_parc ( nv, parcel_pos, parcel_sca, tab )
!  
!  Must be called in the correct order in which they have been registered
!   
    nr = 11
    call write_lag_var ('x    ', nvar, ncid, varids, start, count, tab(:,1))
    call write_lag_var ('y    ', nvar, ncid, varids, start, count, tab(:,2))
    call write_lag_var ('z    ', nvar, ncid, varids, start, count, tab(:,3))
    call write_lag_var ('u    ', nvar, ncid, varids, start, count, tab(:,4))
    call write_lag_var ('v    ', nvar, ncid, varids, start, count, tab(:,5))
    call write_lag_var ('w    ', nvar, ncid, varids, start, count, tab(:,6))
    call write_lag_var ('k    ', nvar, ncid, varids, start, count, tab(:,7))
    call write_lag_var ('t    ', nvar, ncid, varids, start, count, tab(:,8))
    call write_lag_var ('pt   ', nvar, ncid, varids, start, count, tab(:,9))
    call write_lag_var ('qt   ', nvar, ncid, varids, start, count, tab(:,10))
    call write_lag_var ('buoy ', nvar, ncid, varids, start, count, tab(:,11))
    if (lmicro > 0) then
      nr = nr + 1
      call write_lag_var ('ql   ', nvar, ncid, varids, start, count, tab(:,nr))
    endif
    if (lmicro > 1) then
      nr = nr + 1
      call write_lag_var ('qi   ', nvar, ncid, varids, start, count, tab(:,nr))
    endif
    if (nscal > 0) then
     do h = 1, nscal
        write(car1,'(i1)') h
        call write_lag_var ('sca'//car1//' ', nvar, ncid, varids, start, count, tab(:,nr+h))
      enddo
      nr = nr + nscal
    endif
    if (out_diagw) then
      do h = 1, 9
	  call write_lag_var ('wdg  ', nvar, ncid, varids, start, count, tab(:,nr+h))
      enddo
      nr = nr + 9
    endif
    if (out_diagl) then
      do h = 1, 9
        call write_lag_var ('ldg  ', nvar, ncid, varids, start, count, tab(:,nr+h))
      enddo
      nr = nr + 9
    endif
!
    deallocate(tab)
    return
    
    contains

    subroutine gather_parc (nv, pos, sca, tab)
   
    integer, intent(in) :: nv 
    type (position), dimension(nparc), intent(in) :: pos
    type (scalar), dimension(nparc), intent(in) :: sca  
    real, dimension(:,:), intent(inout) :: tab
    integer :: nr, i, h, ierr
!
    if (verbose > 2) call write_debug('Starting gather_parc')
!
    tab = 0.
!
    do i = 1, nparc
!
#ifdef SPMD
      if ( pos(i)%iproc == mypid .and. pos(i)%exist ) then
        nr = 11
        tab(i,1) = pos(i)%x
        tab(i,2) = pos(i)%y
        tab(i,3) = pos(i)%z
    
        tab(i,4) = sca(i)%u
        tab(i,5) = sca(i)%v
        tab(i,6) = sca(i)%w
        tab(i,7) = sca(i)%tke

        tab(i,8) = sca(i)%t
        tab(i,9) = sca(i)%pt
        tab(i,10) = sca(i)%qt
        tab(i,11) = sca(i)%buoy
    
        if (lmicro > 0) then
          nr = nr + 1
          tab(i,nr) = sca(i)%ql 
        endif
        if (lmicro > 1) then 
          nr = nr + 1
          tab(i,nr) = sca(i)%qi
        endif
        do h = 1, nscal
          tab(i,nr+h) = sca(i)%scal(h)
          nr = nr + 1
        enddo
        if (out_diagw) then 
          do h = 1, 9
            tab(i,nr+h) = sca(i)%dw(h)
          enddo
          nr = nr + 9
        endif
        if (out_diagl) then
          do h = 1, 9
            tab(i,nr+h) = sca(i)%dq(h)
          enddo
          nr = nr + 9
        endif
      endif
!
      if ( pos(i)%iproc /= 0 .and. pos(i)%iproc == mypid .and. pos(i)%exist ) then
        call MPI_SEND (tab(i,:), nv, REALTYPE, 0, i, MPI_COMM_WORLD, ierr)
      else if ( pos(i)%iproc /= 0 .and. mypid == 0 .and. pos(i)%exist ) then
        call MPI_RECV (tab(i,:), nv, REALTYPE, pos(i)%iproc, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      endif
!
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
#endif
!
    enddo
!
    if (verbose > 2) call write_debug('Terminating gather_parc')

    end subroutine

    subroutine write_lag_var (var, nvar, ncid, varids, start, count, tot)
    
    character(len=5), intent(in) :: var
    integer, intent(inout) :: nvar
    integer, intent(in) :: ncid
    integer, dimension(:), intent(in) :: varids
    integer, dimension(:), intent(in) :: start, count
    real, dimension(1:nparc), intent(in) :: tot
    
    nvar=nvar+1

    if (mypid == 0) call check( nf90_put_var(ncid, varids(nvar), tot, start=start, count=count), var )

    end
    
    end subroutine
#endif 
!    
    ! ================================================================
    !> Prepare and write data of output variables to netCDF file
    !! 
    !!  
    subroutine write_var (var,nvar,ncid,varids,start,count,cnt,tab,cal_mean,yav,filter)
    
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var
      character(len=*), intent(in) :: filter
      logical, intent(in) :: cal_mean, yav       !> switch for averaging on/off
      integer, intent(inout) :: nvar
      integer, dimension(:) :: varids
      integer, dimension(:), intent(in) :: start, count
      real, intent(in) :: cnt
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: tab !> Input data array
!      
      integer :: i, j, k, ierr, kz
      integer, dimension(3) :: index_max, index_min
      real, dimension(1:nx,1:ny,1:nz)  :: tab3d
      real, dimension(ip_start:ip_end,jp_start:jp_end)  :: tab2d
      real, dimension(nz) :: tabav
      real :: data_max, data_min, maxall, minall, zz, xx, yy
!      
      real, allocatable, dimension(:,:,:) :: tot, q
      real, allocatable, dimension(:,:) :: slice
!
      if (verbose > 1) call write_debug('Starting writing in output variable '//trim(var)) 
!    
!  Write profiles
!    
      nvar=nvar+1
!
      if (cal_mean) then
!
        if (trim(filter) == 'env' .or. trim(filter) == 'cl' .or. trim(filter) == 'int' .or.  &
            trim(filter) == 'out' .or. trim(filter) == 'cor' .or. trim(filter) == 'pig') then
          allocate( q(it_start:it_end,jt_start:jt_end,1:nz) )
          if (lmicro > 0) q = hydromtr(drop)%q + hydromtr(rain)%q
          if (lmicro > 1) q = q + hydromtr(ice)%q
          if (lmicro > 2) q = q + hydromtr(grau)%q + hydromtr(snow)%q 
        endif
!
        if (lmicro == 0) then
          call horav (tab, tabav)
        else if (trim(filter) == 'env') then
          call horav (tab, tabav, xq1=hydromtr(drop)%q, xc1=qthres, filter1='lower')	
        else if (trim(filter) == 'cl') then
          call horav (tab, tabav, xq1=hydromtr(drop)%q, xc1=qthres, filter1='larger')	
        else if (trim(filter) == 'cor') then
          call horav (tab, tabav, xq1=hydromtr(drop)%q, xc1=qthres, filter1='larger', xq2=wind%w, xc2=wthres, filter2='larger')
        else if (trim(filter) == 'upd') then
          call horav (tab, tabav, xq1=wind%w, xc1=wthres, filter1='larger')	
        else if (trim(filter) == 'dnd') then
          call horav (tab, tabav, xq1=wind%w, xc1=-wthres, filter1='lower')	
        else if (trim(filter) == 'int') then
          call horav (tab, tabav, xq1=hydromtr(drop)%q, xc1=qthres, filter1='larger', xq2=wind%w, xc2=wthres, filter2='larger', filter3='inside')
        else if (trim(filter) == 'out') then
          call horav (tab, tabav, xq1=hydromtr(drop)%q, xc1=qthres, filter1='larger', xq2=wind%w, xc2=wthres, filter2='larger', filter3='outside')
        else if (trim(filter) == 'ent') then
          call horav (tab, tabav, xq1=entrain%went, xc1=1.e-5, filter1='larger')	
        else if (trim(filter) == 'det') then
          call horav (tab, tabav, xq1=entrain%wdet, xc1=1.e-5, filter1='larger')	
        else if (trim(filter) == 'hal') then
          call horav (tab, tabav, xq1=entrain%activ, xc1=1.e-2, filter1='larger')	
	else
          call horav (tab, tabav)
        endif
!
        if ( mypid == 0 ) call check( nf90_put_var(ncid, varids(nvar), cnt*tabav(1:nz), start=start, count=count), var )
!
        if (trim(filter) == 'env' .or. trim(filter) == 'cl' .or. trim(filter) == 'int' .or.  &
            trim(filter) == 'out' .or. trim(filter) == 'cor' .or. trim(filter) == 'pig') deallocate(q)
!
!  Write mean y slices
!
      else if (yav) then
!
	allocate(tot(1:nx,1:ny,1:nz), slice(1:nx,1:nz))
!
#if ( defined SPMD )
        call collect ( tab,  tot )
#else
        tot = tab
#endif
!
        if ( mypid == 0 ) then
  	  slice = 0.0
  	  do j = j1, j2
    	    slice = slice + tot(:,j,:)
  	  enddo
  	  slice = slice/real(ny-5)
!
	  call check( nf90_put_var(ncid, varids(nvar), cnt*slice(i1:i2,1:NZO), start=start, count=count), var )
	endif
!
	deallocate(tot, slice)
!
      else
!
!  Write full data
!
        if ( minmax ) then
#if ( defined SPMD )
          data_max  = maxval( tab(it_start:it_end,jt_start:jt_end,1:nz) )
          data_min  = minval( tab(it_start:it_end,jt_start:jt_end,1:nz) )
          CALL MPI_ALLREDUCE (data_max, maxall, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
          CALL MPI_ALLREDUCE (data_min, minall, 1, REALTYPE, MPI_MIN, MPI_COMM_WORLD, ierr) 
#else
          maxall = maxval( tab(it_start:it_end,jt_start:jt_end,1:nz) )
          minall = minval( tab(it_start:it_end,jt_start:jt_end,1:nz) )
#endif
        endif
!
#ifndef PARALLEL_OUT
	allocate(tot(1:nx,1:ny,1:nz))
#else
	allocate(tot(ip_start:ip_end,jp_start:jp_end,1:nz))
#endif
!
!  Write full 3D
!
#if ( defined SPMD ) && ( !defined PARALLEL_OUT )
        call collect ( tab,  tot )
#else
        tot = tab
#endif
!   
#ifndef PARALLEL_OUT
        if ( mypid == 0 ) then 
#endif
!
#ifdef MODEL_3D
          call check( nf90_put_var(ncid, varids(nvar), cnt*tot(i1:i2,j1:j2,1:NZO), start=start, count=count), var )
#else
          call check( nf90_put_var(ncid, varids(nvar), cnt*tot(i1:i2,1,1:NZO), start=start, count=count), var )
#endif
!
	  deallocate(tot)
!
#ifndef PARALLEL_OUT
        endif
#endif
      endif
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
 
end module netcdfmod

