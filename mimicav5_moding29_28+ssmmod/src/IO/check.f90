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
!  CHECK.F90:
!	subroutine checkwind ()
!	subroutine checkx (min,max,x)
!
!  Purpose:
!	Searching wind speeds which potentially could violate the CFL
!	condition, searching for erroneous values for scalars.
!
!  Author:
!	Chien Wang
!	Center for Clouds, Chemistry, and Climate
!
! =========================================================

  module check

#ifdef SPMD
  USE mpi
#endif
  USE gridno
  USE shared_all
  USE typedef_svalue
  USE netcdfmod

  IMPLICIT NONE
  
  private
  
  public :: checkall
  
  CONTAINS

!   ===================================================
    subroutine checkall ( windl, statel, hydromtrl ) 
!   ===================================================

! ===
  type (atm_winds) :: windl
  type (atm_state) :: statel
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
!
  integer :: h
!
  if (verbose > 0) call write_debug('Starting checking')
!
!----------------------------------------------------------!
!                     Check interface                      !
!----------------------------------------------------------!
!
!  Check bounds of conserved variables
!
!   call checkx (200., 700., statel%es, 1)  
  call checkx (200., 4000., statel%es, 1)  
!
 ! call checkx (1.e-10, 0.1, statel%qt, 2)  
!
  call checkx(1.e-15,1.,statel%qt,2)
! 
!    do h = 1, nhydro
!    call checkx (-1.e-12, 1., hydromtrl(h)%q, 2+h)
!  enddo
!
!  Check CFL stability 
!
  call checkwind ( windl )
!
  if (mypid == 0) flush(7)
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating checking')
!
  return
  end
!
!   ===================================================
    subroutine checkwind ( windl ) 
!   ===================================================

! ===
! === CHECKWIND.F:  A program for checking if wind speed is too large
! ===                which might lead to the linear instability
! ===
    type(atm_winds) :: windl
!
    integer :: i, j, k
    integer :: stopflag, ierr, ntauold
    real    :: xx, fac
    real    :: cflx, cfly, cflz, cfl, cflmax, cflcmax
    real    :: dtref, dtold, dx1, dy1
    real    :: cflcx, cflcy, cflcz, cflc
!
    type (s_value) :: local_value
!
!----------------------------------------------------------!
!                 Define max allowed vel.                  !
!----------------------------------------------------------!
!
    stopflag = 0
    dx1 = 1./dx
#ifdef MODEL_3D
    dy1 = 1./dy
#else
    dy1 = 1.
#endif
!
    call init_svalue ( local_value )
!
!----------------------------------------------------------!
!                    Calculate max CFL                     !
!----------------------------------------------------------!
!
        dtref = dt0
	cfly = 0.
	cflcy = 0.
	cflmax = 0.
	cflcmax = 0.
	do k=1,nz
	do j=jt_start,jt_end
	do i=it_start,it_end
!
!  Calculate CFL
!
	  cflx = dt0*abs(windl%u(i,j,k))*dx1 
#ifdef MODEL_3D
	  cfly = dt0*abs(windl%v(i,j,k))*dy1 
#endif
	  cflz = dt0*abs(windl%w(i,j,k))*dz1w(k) 
!
	  cfl = max(cflx,cfly,cflz)
	  if(cfl > cflmax)then
	    cflmax = cfl
	    local_value%value = cfl
	    local_value%i     = i
	    local_value%j     = j
	    local_value%k     = k
	  endif
!
!  Calculate acoustic CFL if needed
!
#ifndef ANELASTIC
	  cflcx = dtau*(abs(windl%u(i,j,k))+sqrt(thermo_prop%cs2(i,j,k)))*dx1 
#ifdef MODEL_3D
	  cflcy = dtau*(abs(windl%v(i,j,k))+sqrt(thermo_prop%cs2(i,j,k)))*dy1 
#endif
	  cflcz = dtau*(abs(windl%w(i,j,k))+sqrt(thermo_prop%cs2(i,j,k)))*dz1w(k) 
!
          cflc  = max(cflcx,cflcy,cflcz) 
	  if (cflc > cflcmax) then
            cflcmax = cflc
          endif
#else
	  cflcmax = 0.
#endif
	end do
	end do
	end do
!
!----------------------------------------------------------!
!                    Find max CFL value                    !
!----------------------------------------------------------!
!
#if ( defined SPMD )
     CALL MPI_ALLREDUCE (cflmax, maxcfl, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
     CALL MPI_ALLREDUCE (cflcmax, maxcflc, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
#else
     maxcfl = cflmax
     maxcflc = cflcmax
#endif
!
!  Stop if needed
!
     if ( maxcfl > 2. .or. maxcfl < 1.e-7 ) then
#ifdef MODEL_3D
       call write_dim(3, 1., .false., .true.)
#else
       call write_dim(2, 1., .false., .true.)
#endif
!
       write(7,*)'before stop 1', maxcfl
       call stop_mimica ( 1 )
     endif
!
!----------------------------------------------------------!
!           Limit time step considering max cfl            !
!----------------------------------------------------------!
!
!  Advective time-step
!
	if (maxcfl >= cfl_max .and. mypid == 0 .and. ldtfix == 0) then
	  xx = maxcfl/dtref
	  fac = 0.95*cfl_max/xx
!
	  write(7,*)  'NTIME=', ntime, 'WARNING:'
	  write(7,*)  'CFL (', maxcfl, ') EXCEEDS MAX ALLOWED VALUE AT'
	  write(7,*)  local_value%i, local_value%j, local_value%k
	  dtold = dt0
	  dt0 = fac
!
	  write(7,*) 'NEW TIME STEP IS', dt0, ' (PREVIOUS: ', dtold, ')'
	  write(7,*) ''
!
	  if (limit_ts > 1. .and. dt0 <= dt0_i/max(1.,limit_ts)) then
	    write(7,*) 'TIME STEP IS TOO SMALL: STOPING THE RUN IN checkwind'
	    print*, 'ERROR: Time step too small, abort'
	    stopflag = 1
	  endif
	endif
!
!  Acoustic time-step (number of acoustic sub-steps)
!
#ifndef ANELASTIC
	if (maxcflc >= cfl_max .and. mypid == 0 .and. ldtfix == 0) then
	  xx = maxcflc*real(ntau)/dtref
	  fac = xx*dt0/(0.95*cfl_max)  ! the 0.9 factor is a security
!
	  write(7,*)  'NTIME=', ntime, 'WARNING:'
	  write(7,*)  'ACOUSTIC CFL (', maxcflc, ') EXCEEDS MAX ALLOWED VALUE'
	  ntauold = ntau
	  ntau = ceiling(fac)
!
	  write(7,*) 'NEW ACOUSTIC STEP # IS', ntau, ' (PREVIOUS: ', ntauold, ')'
	  write(7,*) ''
!
	  if (ntau >= 400) then
	    write(7,*) 'TOO MANY ACOUSTIC STEPS: STOPING THE RUN IN checkwind'
	    print*, 'ERROR: Too many acoustic steps, abort'
	    stopflag = 1
	  endif
	endif
#endif
!
!----------------------------------------------------------!
!          Increase time step considering min cfl          !
!----------------------------------------------------------!
!
!  Advective time-step
!
	if (maxcfl <= cfl_min .and. maxcfl > 0. .and. mypid == 0 .and. ldtfix == 0) then
	  xx = maxcfl/dtref
	  fac = 1.15*cfl_min/xx 
!
	  write(7,*)  'NTIME=', ntime, 'WARNING:'
	  write(7,*)  'CFL (', maxcfl, ') LOWER THAN MIN VALUE'
	  dtold = dt0
	  dt0 = min(fac,2.*dt0_i)
!
	  write(7,*) 'NEW TIME STEP IS', dt0, ' (PREVIOUS: ', dtold, ')'
	  write(7,*) ''
	endif
!
!  Acoustic time-step (number of acoustic sub-steps)
!
#ifndef ANELASTIC
	if (maxcflc <= cfl_min .and. maxcflc > 0. .and. mypid == 0 .and. ldtfix == 0) then
	  xx = maxcflc*real(ntau)/dtref
	  fac = xx*dt0/(1.15*cfl_min) 
!
	  write(7,*)  'NTIME=', ntime, 'WARNING:'
	  write(7,*)  'ACOUSTIC CFL (', maxcflc, ') LOWER THAN MIN VALUE'
	  ntauold = ntau
	  ntau = min(ceiling(fac),30)
!
	  write(7,*) 'NEW ACOUSTIC STEP # IS', ntau, ' (PREVIOUS: ', ntauold, ')'
	  write(7,*) ''
	endif
#endif
!
!----------------------------------------------------------!
!              Limit time step by saturation               !
!----------------------------------------------------------!
!
#ifndef SAT_ADJ
	if (ldtmic .and. mypid == 0 .and. ldtfix == 0) then
	  dtold = dt0
	  dt0 = max(min(dtold,0.5*dtmic),0.01)
	  if (dt0 < dtold) then
	    write(7,*) 'NTIME=', ntime, 'WARNING:'
	    write(7,*) 'MICROPHYSICS TIME SCALE TOO SMALL'
	    write(7,*) 'NEW TIME STEP IS', dt0
	    write(7,*) ''
	  endif
	endif
#endif
!
!----------------------------------------------------------!
!                Broadcast new time steps                  !
!----------------------------------------------------------!
!
#ifdef SPMD
	call MPI_BCAST (stopflag, 1, INTTYPE,  0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (dt0,      1, REALTYPE, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST (ntau,     1, INTTYPE, 0, MPI_COMM_WORLD, ierr)
#endif
!
!  Stop the run if necessary
!
       if (stopflag == 1) then
#ifdef MODEL_3D
         call write_dim(3, 1., .false., .true.)
#else
         call write_dim(2, 1., .false., .true.)
#endif
!
         write(7,*)'before stop 2', maxcfl
         call stop_mimica ( 1 ) 
       endif
!
!-----------------------------------------------------------!
!
  return
  end
!
!	===================================================
	subroutine checkx ( minx, maxx, x, flag )
!	===================================================

! ===
! === CHECKX.F:  A program for checking if Potential temperature
! ===            becomes too large or too small
! ===
!
	integer :: flag
	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x
!
    	integer :: ind(3), ierr, stopflag
	character(len=1) :: car1
	real  :: maxx, minx, maxq, minq, minqall, maxqall
	type(s_value) :: local_value
!
!-----------------------------------------------------------!
!                    Initializations                        !
!-----------------------------------------------------------!
!
	call init_svalue ( local_value )
!
	stopflag = 0
!
    	maxq = maxval(x(it_start:it_end,jt_start:jt_end,1:nz))
    	minq = minval(x(it_start:it_end,jt_start:jt_end,1:nz))
!
!  Find min/max
!
	if (maxq > maxx) then
	  ind = maxloc(x(it_start:it_end,jt_start:jt_end,1:nz))
	  local_value%i     = ind(1)
	  local_value%j     = ind(2)
	  local_value%k     = ind(3)
	endif
!
	if (minq < minx) then
	  ind = minloc(x(it_start:it_end,jt_start:jt_end,1:nz))
	  local_value%i     = ind(1)
	  local_value%j     = ind(2)
	  local_value%k     = ind(3) 
	endif
!
!  Stop the run if necessary
!
#if ( defined SPMD )
      call MPI_ALLREDUCE (maxq, maxqall, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
      call MPI_ALLREDUCE (minq, minqall, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
#else
      maxqall = maxq
      minqall = minq
#endif
!
!  Error
!
    	if ( flag == 1 .and. (maxqall > maxx .or. minqall < minx) .and. mypid == 0 ) then
	  stopflag = 1
	  write(7,*)  'NTIME=', ntime, 'ERROR: PT OUTSIDE OF ALLOWED BOUNDS: ', maxqall, minqall
	  write(7,*)  local_value%i, local_value%j, local_value%k
	endif
!
    	if ( flag == 2 .and. (maxqall > maxx .or. minqall < minx) .and. mypid == 0 ) then
	  stopflag = 1
	  write(7,*)  'NTIME=', ntime, 'ERROR: QT OUTSIDE OF ALLOWED BOUNDS: ', maxqall, minqall
	  write(7,*)  local_value%i, local_value%j, local_value%k
	endif
!
    	if ( flag > 2 .and. (maxqall > maxx .or. minqall < minx) .and. mypid == 0 ) then
	  write(car1,'(i1)') flag-2
	  write(7,*)  'NTIME=', ntime, 'ERROR: Q'//car1//' OUTSIDE OF ALLOWED BOUNDS: ', maxqall, minqall
	  write(7,*)  local_value%i, local_value%j, local_value%k
	endif
!
#if ( defined SPMD )
      call MPI_BCAST (stopflag,      1, INTTYPE, 0, MPI_COMM_WORLD, ierr)
#endif
!
!  Stop
!
       if (stopflag == 1) then
#ifdef MODEL_3D
         call write_dim(3, 1., .false., .true.)
#else
         call write_dim(2, 1., .false., .true.)
#endif
!
       write(7,*)'before stop 3', maxcfl
         call stop_mimica ( 1 )
       endif
!
 	return
      end	

      subroutine init_svalue (data1)

	USE shared_data
        USE typedef_svalue
        
        IMPLICIT NONE
        
        type(s_value) :: data1

        data1%value = 0. 
        data1%i     = 0
        data1%j     = 0
        data1%k     = 0
        return
      end subroutine init_svalue

  end module check
