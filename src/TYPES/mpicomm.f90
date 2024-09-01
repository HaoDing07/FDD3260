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
!  SPMD.F:                 
!
!  Purpose:
!	Collections of subroutines and functions for spmd runs.
!
!  Author
!	Julien Savre
!	MIM, Ludwig-Maximilian Universität, Munchen
!
! ================================================================

	module mpicomm

#if ( defined SPMD )
	USE mpi
	USE gridno
	USE shared_data	
	
	private
	
	interface collect
	  MODULE PROCEDURE collect_i_2d, collect_r_2d, collect_3d
	end interface collect
	
	interface distribute
	  MODULE PROCEDURE distribute_3d, distribute_2d
	end interface distribute
	
	interface exchange_x
	  MODULE PROCEDURE exchange_x_x, exchange_x_x_2d
	end interface exchange_x
	
	interface exchange_y
	  MODULE PROCEDURE exchange_y_x, exchange_y_x_2d
	end interface exchange_y
	
	public :: collect, distribute, exchange_x, exchange_y	

	CONTAINS

!	==========================================		
	subroutine exchange_x_x ( a, tag )
!	==========================================
!
! --- exchange data for u
!
	IMPLICIT NONE
	
	integer :: ierr, req(4),status(MPI_STATUS_SIZE,4)
	integer :: tag, dimy, dim
	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: a
	real, allocatable, dimension(:,:,:) :: recl, recr
	
! -------------------------------------------------------
!
! --- Allocate
!
	dimy = jp_end-jp_start+1
	dim = 3*dimy*nz
	req = MPI_REQUEST_NULL
	allocate( recl(1:3,1:dimy,1:nz), recr(1:3,1:dimy,1:nz) )

!
! --- first send to your left and right
!
    
	call MPI_ISEND ( a, 1, GHOSTX_send_s, left_p, tag, MPI_COMM_WORLD, req(1), ierr )

	call MPI_ISEND ( a, 1, GHOSTX_send_e, right_p, tag+1, MPI_COMM_WORLD, req(2), ierr )	

!
! --- then send to your right and receive from left
!
   
	call MPI_IRECV ( recr, dim, REALTYPE, right_p, MPI_ANY_TAG, MPI_COMM_WORLD, req(3), ierr )	
   
	call MPI_IRECV ( recl, dim, REALTYPE, left_p, MPI_ANY_TAG, MPI_COMM_WORLD, req(4), ierr )	

!
! --- Wait for all communications to end
!

	call MPI_WAITALL ( 4, req, status, ierr )
	
!
! --- Update ghost cells
!

        if (ip_start == 1) then
	  a(ip_end-1:ip_end,jp_start:jp_end,1:nz) = recr(1:2,1:dimy,1:nz)
	  if (.not.nest_run .and. (bcl(1) == 'per' .or. bcl(2) == 'per'))  &
	      a(ip_start:ip_start+2,jp_start:jp_end,1:nz) = recl(1:3,1:dimy,1:nz)
        else if (ip_end == nx) then
	  if (.not.nest_run .and. (bcl(1) == 'per' .or. bcl(2) == 'per'))  &
	      a(ip_end-1:ip_end,jp_start:jp_end,1:nz) = recr(1:2,1:dimy,1:nz)
	  a(ip_start:ip_start+1,jp_start:jp_end,1:nz) = recl(2:3,1:dimy,1:nz)
        else
	  a(ip_end-1:ip_end,jp_start:jp_end,1:nz) = recr(1:2,1:dimy,1:nz)
	  a(ip_start:ip_start+1,jp_start:jp_end,1:nz) = recl(2:3,1:dimy,1:nz)
	endif

	deallocate(recl, recr)
	
	return
	end subroutine exchange_x_x

!	==========================================		
	subroutine exchange_y_x ( a, tag )
!	==========================================
!
! --- exchange data for u
!

	IMPLICIT NONE
	
	integer :: ierr, req(4),status(MPI_STATUS_SIZE,4)
	integer :: tag, dim, dimx
	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: a
	real, allocatable, dimension(:,:,:) :: recl, recr
	
! -------------------------------------------------------
!
! --- Allocate
!
	dimx = ip_end-ip_start+1
	dim = 3*dimx*nz
	req = MPI_REQUEST_NULL
	allocate(recl(1:dimx,1:3,1:nz), recr(1:dimx,1:3,1:nz))

!
! --- first send to your left and right
!
    
	call MPI_ISEND ( a, 1, GHOSTY_send_s, back_p, tag, MPI_COMM_WORLD, req(1), ierr )	

	call MPI_ISEND ( a, 1, GHOSTY_send_e, front_p, tag+1, MPI_COMM_WORLD, req(2), ierr )	

!
! --- then send to your right and receive from left
!
   
	call MPI_IRECV ( recr, dim, REALTYPE, front_p, MPI_ANY_TAG, MPI_COMM_WORLD, req(3), ierr )	
   
	call MPI_IRECV ( recl, dim, REALTYPE, back_p, MPI_ANY_TAG, MPI_COMM_WORLD, req(4), ierr )	

!
! --- Wait for all communications to end
!

	call MPI_WAITALL ( 4, req, status, ierr )

!
! --- Update ghost cells
!

        if (jp_start == 1) then
  	  a(ip_start:ip_end,jp_end-1:jp_end,1:nz) = recr(1:dimx,1:2,1:nz)
  	  if (.not.nest_run .and. (bcl(3) == 'per' .or. bcl(4) == 'per'))  &
	      a(ip_start:ip_end,jp_start:jp_start+2,1:nz) = recl(1:dimx,1:3,1:nz)
        else if (jp_end == ny) then
	  if (.not.nest_run .and. (bcl(3) == 'per' .or. bcl(4) == 'per'))  &
	      a(ip_start:ip_end,jp_end-1:jp_end,1:nz) = recr(1:dimx,1:2,1:nz)
	  a(ip_start:ip_end,jp_start:jp_start+1,1:nz) = recl(1:dimx,2:3,1:nz)
        else
	  a(ip_start:ip_end,jp_end-1:jp_end,1:nz) = recr(1:dimx,1:2,1:nz)
	  a(ip_start:ip_end,jp_start:jp_start+1,1:nz) = recl(1:dimx,2:3,1:nz)
	endif
	
	deallocate(recl, recr)
	
	return
	end subroutine exchange_y_x
	 

!	==========================================		
	subroutine exchange_x_x_2d ( a, tag )
!	==========================================
!
! --- exchange data for u
!
	IMPLICIT NONE
	
	integer :: ierr, req(4), status(MPI_STATUS_SIZE,4)
	integer :: tag, dimy, dim
	real, dimension(ip_start:ip_end,jp_start:jp_end) :: a
	real, allocatable, dimension(:,:) :: recl, recr
	
! -------------------------------------------------------
!
! --- Allocate
!
	dimy = jp_end-jp_start+1
	dim = 3*dimy
	req = MPI_REQUEST_NULL
	allocate(recl(1:3,1:dimy), recr(1:3,1:dimy))

!
! --- first send to your left and right
!
    
	call MPI_ISEND ( a, 1, GHOSTX_send_s_2d, left_p, tag, MPI_COMM_WORLD, req(1), ierr )

	call MPI_ISEND ( a, 1, GHOSTX_send_e_2d, right_p, tag+1, MPI_COMM_WORLD, req(2), ierr )	

!
! --- then send to your right and receive from left
!
   
	call MPI_IRECV ( recr, dim, REALTYPE, right_p, MPI_ANY_TAG, MPI_COMM_WORLD, req(3), ierr )	
   
	call MPI_IRECV ( recl, dim, REALTYPE, left_p, MPI_ANY_TAG, MPI_COMM_WORLD, req(4), ierr )	

!
! --- Wait for all communications to end
!

	call MPI_WAITALL ( 4, req, status, ierr )
	
!
! --- Update ghost cells
!

        if (ip_start == 1) then
	  a(ip_end-1:ip_end,jp_start:jp_end) = recr(1:2,1:dimy)
	  if (.not.nest_run .and. (bcl(1) == 'per' .or. bcl(2) == 'per'))  &
	      a(ip_start:ip_start+2,jp_start:jp_end) = recl(1:3,1:dimy)
        else if (ip_end == nx) then
	  if (.not.nest_run .and. (bcl(1) == 'per' .or. bcl(2) == 'per'))  &
	      a(ip_end-1:ip_end,jp_start:jp_end) = recr(1:2,1:dimy)
	  a(ip_start:ip_start+1,jp_start:jp_end) = recl(2:3,1:dimy)
        else
	  a(ip_end-1:ip_end,jp_start:jp_end) = recr(1:2,1:dimy)
	  a(ip_start:ip_start+1,jp_start:jp_end) = recl(2:3,1:dimy)
	endif
	
	deallocate(recl, recr)
	
	return
	end subroutine exchange_x_x_2d

!	==========================================		
	subroutine exchange_y_x_2d ( a, tag )
!	==========================================
!
! --- exchange data for u
!

	IMPLICIT NONE
	
	integer :: ierr, req(4), status(MPI_STATUS_SIZE,4)
	integer :: tag, dim, dimx
	real, dimension(ip_start:ip_end,jp_start:jp_end) :: a
	real, allocatable, dimension(:,:) :: recl, recr
	
! -------------------------------------------------------
!
! --- Allocate
!
	dimx = ip_end-ip_start+1
	dim = 3*dimx
	req = MPI_REQUEST_NULL
	allocate(recl(1:dimx,1:3), recr(1:dimx,1:3))

!
! --- first send to your left and right
!
    
	call MPI_ISEND ( a, 1, GHOSTY_send_s_2d, back_p, tag, MPI_COMM_WORLD, req(1), ierr )	

	call MPI_ISEND ( a, 1, GHOSTY_send_e_2d, front_p, tag+1, MPI_COMM_WORLD, req(2), ierr )	

!
! --- then send to your right and receive from left
!
   
	call MPI_IRECV ( recr, dim, REALTYPE, front_p, MPI_ANY_TAG, MPI_COMM_WORLD, req(3), ierr )	
   
	call MPI_IRECV ( recl, dim, REALTYPE, back_p, MPI_ANY_TAG, MPI_COMM_WORLD, req(4), ierr )	

!
! --- Wait for all communications to end
!

	call MPI_WAITALL ( 4, req, status, ierr )

!
! --- Update ghost cells
!

        if (jp_start == 1) then
  	  a(ip_start:ip_end,jp_end-1:jp_end) = recr(1:dimx,1:2)
  	  if (.not.nest_run .and. (bcl(3) == 'per' .or. bcl(4) == 'per'))  &
	      a(ip_start:ip_end,jp_start:jp_start+2) = recl(1:dimx,1:3)
        else if (jp_end == ny) then
	  if (.not.nest_run .and. (bcl(3) == 'per' .or. bcl(4) == 'per'))  &
	      a(ip_start:ip_end,jp_end-1:jp_end) = recr(1:dimx,1:2)
	  a(ip_start:ip_end,jp_start:jp_start+1) = recl(1:dimx,2:3)
        else
	  a(ip_start:ip_end,jp_end-1:jp_end) = recr(1:dimx,1:2)
	  a(ip_start:ip_end,jp_start:jp_start+1) = recl(1:dimx,2:3)
	endif
	
	deallocate(recl, recr)
	
	return
	end subroutine exchange_y_x_2d


!	======================================================		
	subroutine collect_3d ( local,  global )
!	======================================================		

	IMPLICIT NONE

	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: local
	real, dimension(nx,ny,nz), intent(out) :: global

	integer :: ierr, status(MPI_STATUS_SIZE), send_req
	integer :: i, j, proc
	integer :: BLOCK3Dproc
	integer :: i1, j1, ni, nj

! -------------------------------------------------------
!
!  Master proc: receive data
!
    i1 = 1
    j1 = 1
    ni = it_end-it_start+1
    nj = jt_end-jt_start+1
!
    if(mypid == 0) then          
      global(ip_start:ip_end,jp_start:jp_end,1:nz) = local(ip_start:ip_end,jp_start:jp_end,1:nz) 
!
      do i=0,nprocx-1
        do j = 0, nprocy-1		
          proc = j+i*nprocy			     
 	  if (proc /= 0) then
	    i1 = 4+i*ni
#ifdef MODEL_3D
#ifdef DECOMP_2D
	    j1 = 4+j*nj
#else
	    j1 = jt_start
#endif
#endif
!
	    call MPI_TYPE_CREATE_SUBARRAY (3, (/nx,ny,nz/), (/ni,nj,nz/), (/i1-1,j1-1,0/), 	&
		    MPI_ORDER_FORTRAN, REALTYPE, BLOCK3Dproc, ierr)
!
	    call MPI_TYPE_COMMIT (BLOCK3Dproc, ierr)
!
            call MPI_Recv( global, 1, BLOCK3Dproc, proc, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr ) 
!
	    call MPI_TYPE_FREE (BLOCK3Dproc, ierr)
          endif	     
        enddo
      enddo
!
!  Periodic Boundary conditions
!
      if (bcl(1) == 'per' .or. bcl(2) == 'per') then
        global(nx,:,:) = global(5,:,:)
        global(nx-1,:,:) = global(4,:,:)
        global(1,:,:) = global(nx-4,:,:)
        global(2,:,:) = global(nx-3,:,:)
        global(3,:,:) = global(nx-2,:,:)
      endif
!
#ifdef MODEL_3D
      if (bcl(3) == 'per' .or. bcl(4) == 'per') then
        global(:,ny,:) = global(:,5,:)
        global(:,ny-1,:) = global(:,4,:)
        global(:,1,:) = global(:,ny-4,:)
        global(:,2,:) = global(:,ny-3,:)
        global(:,3,:) = global(:,ny-2,:)
      endif
#endif
!
!  Slave procs: send data
!
    else
      call MPI_Isend( local, 1, BLOCK3D, 0, mypid, MPI_COMM_WORLD, send_req, ierr)
!
      call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)   ! block until Isend is done
    endif
	return
	end subroutine collect_3d
!	 	 
!	===========================================	      
	subroutine collect_r_2d ( local, global )
!	===========================================	      

	IMPLICIT NONE

	real, dimension(ip_start:ip_end,jp_start:jp_end), intent(in) :: local
	real, dimension(nx,ny), intent(out) :: global

	integer :: ierr, send_req, status(MPI_STATUS_SIZE)
	integer :: i, j, proc
	integer :: BLOCK2Dproc
	integer :: i1, j1, ni, nj
!
! -------------------------------------------------------
!
!  Master proc: receive data
!
    i1 = 1
    j1 = 1
    ni = it_end-it_start+1
    nj = jt_end-jt_start+1
!
    if(mypid == 0) then
	
      global(ip_start:ip_end,jp_start:jp_end) = local(ip_start:ip_end,jp_start:jp_end)      
!
      do i=0,nprocx-1
        do j = 0, nprocy-1		
          proc = j+i*nprocy			     
 	  if (proc /= 0) then
	    i1 = 4+i*ni
#ifdef MODEL_3D
#ifdef DECOMP_2D
	    j1 = 4+j*nj
#else
	    j1 = jt_start
#endif
#endif
!
	    call MPI_TYPE_CREATE_SUBARRAY (2, (/nx,ny/), (/ni,nj/), (/i1-1,j1-1/), 	&
		    MPI_ORDER_FORTRAN, REALTYPE, BLOCK2Dproc, ierr)
!
	    call MPI_TYPE_COMMIT (BLOCK2Dproc, ierr)
            call MPI_Recv( global, 1, BLOCK2Dproc, proc, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr) 
	    call MPI_TYPE_FREE (BLOCK2Dproc, ierr)
          endif		     
        enddo
      enddo
!
!  Periodic Boundary conditions
!
      if (bcl(1) == 'per' .or. bcl(2) == 'per') then
        global(nx,:) = global(5,:)
        global(nx-1,:) = global(4,:)
        global(1,:) = global(nx-4,:)
        global(2,:) = global(nx-3,:)
        global(3,:) = global(nx-2,:)
      endif
!
#ifdef MODEL_3D
      if (bcl(3) == 'per' .or. bcl(4) == 'per') then
        global(:,ny) = global(:,5)
        global(:,ny-1) = global(:,4)
        global(:,1) = global(:,ny-4)
        global(:,2) = global(:,ny-3)
        global(:,3) = global(:,ny-2)
      endif
#endif
!
!  Slave procs: send data
!
    else
      call MPI_Isend( local, 1, BLOCK2D_r, 0, mypid, MPI_COMM_WORLD, send_req, ierr) 
!
      call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr) 
    endif
	return
	end subroutine collect_r_2d
!
!	===========================================	      
	subroutine collect_i_2d ( local, global )
!	===========================================	      

	IMPLICIT NONE

	integer, dimension(ip_start:ip_end,jp_start:jp_end), intent(in) :: local
	integer, dimension(nx,ny), intent(out) :: global

	integer :: ierr, send_req, status(MPI_STATUS_SIZE)
	integer :: i, j, proc
	integer :: BLOCK2Dproc
	integer :: i1, j1, ni, nj

! -------------------------------------------------------
!
!  Master proc: receive data
!
    i1 = 1
    j1 = 1
    ni = it_end-it_start+1
    nj = jt_end-jt_start+1
!
    if(mypid == 0) then          
      global(ip_start:ip_end,jp_start:jp_end) = local(ip_start:ip_end,jp_start:jp_end)      
!
      do i=0,nprocx-1
        do j = 0, nprocy-1		
          proc = j+i*nprocy			     
 	  if (proc /= 0) then
	    i1 = 4+i*ni
#ifdef MODEL_3D
#ifdef DECOMP_2D
	    j1 = 4*j*nj
#else
	    j1 = jt_start
#endif
#endif
!
	    call MPI_TYPE_CREATE_SUBARRAY (2, (/nx,ny/), (/ni,nj/), (/i1-1,j1-1/), 	&
		    MPI_ORDER_FORTRAN, INTTYPE, BLOCK2Dproc, ierr)
!
	    call MPI_TYPE_COMMIT (BLOCK2Dproc, ierr)
!
            call MPI_Recv( global, 1, BLOCK2Dproc, proc, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr) 
!
	    call MPI_TYPE_FREE (BLOCK2Dproc, ierr)
          endif		     
        enddo
      enddo
!
!  Periodic Boundary conditions
!
      if (bcl(1) == 'per' .or. bcl(2) == 'per') then
        global(nx,:) = global(5,:)
        global(nx-1,:) = global(4,:)
        global(1,:) = global(nx-4,:)
        global(2,:) = global(nx-3,:)
        global(3,:) = global(nx-2,:)
      endif
!
#ifdef MODEL_3D
      if (bcl(3) == 'per' .or. bcl(4) == 'per') then
        global(:,ny) = global(:,5)
        global(:,ny-1) = global(:,4)
        global(:,1) = global(:,ny-4)
        global(:,2) = global(:,ny-3)
        global(:,3) = global(:,ny-2)
      endif
#endif
!
!  Slave procs: send data
!
    else
      call MPI_Isend( local, 1, BLOCK2D_i, 0, mypid, MPI_COMM_WORLD, send_req, ierr)
!
      call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
    endif
	return
	end subroutine collect_i_2d
!
!	===========================================	     
  subroutine distribute_3d ( local, global )
!	===========================================  

	IMPLICIT NONE

	integer :: ierr, rec_req, send_req

	integer :: i, j, proc
	integer :: m, n, i1, i2, j1, j2
	integer :: BLOCK3Dproc

	real, dimension(nx,ny,nz) :: global
	real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: local
!
! -------------------------------------------------------
!
!  Master proc: send data
!
    m = it_end-it_start+1
    n = jt_end-jt_start+1
    j1 = jp_start
    j2 = jp_end
!
    if(mypid == 0) then          
      local(ip_start:ip_end,jp_start:jp_end,1:nz) = global(ip_start:ip_end,jp_start:jp_end,1:nz)
!
      do i = 0, nprocx-1
        do j = 0, nprocy-1		
          proc = j+i*nprocy			     
 	  if (proc /= 0) then
 	    i1 = 2+i*m
 	    i2 = i1+m+3
 	    if (i == 0) i1 = 1
 	    if (i == nprocx-1) i2 = nx
#if (defined MODEL_3D) && (defined DECOMP_2D)
 	    j1 = 2+j*n
 	    j2 = j1+n+3
 	    if ( j == 0 ) j1 = 1
 	    if ( j == nprocy-1 ) j2 = ny
#endif
!
	    call MPI_TYPE_CREATE_SUBARRAY (3, (/nx,ny,nz/), (/i2-i1+1,j2-j1+1,nz/), (/i1-1,j1-1,0/), 	&
		    MPI_ORDER_FORTRAN, REALTYPE, BLOCK3Dproc, ierr)
!
	    call MPI_TYPE_COMMIT (BLOCK3Dproc, ierr)
!
            call MPI_ISend( global, 1, BLOCK3Dproc, proc, proc+1, MPI_COMM_WORLD, send_req, ierr) 
!
	    call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
!
	    call MPI_TYPE_FREE (BLOCK3Dproc, ierr)
          endif		     
        enddo
      enddo
!
!  Slave procs: receive data (non-blocking instruction)
!
    else
      call MPI_IRecv( local, 1, BLOCK3Dp, 0, MPI_ANY_TAG, MPI_COMM_WORLD, rec_req, ierr )  
!
      call MPI_Wait( rec_req, MPI_STATUS_IGNORE, ierr )
    endif
!
    return
    end subroutine distribute_3d
!
!	===========================================	     
  subroutine distribute_2d ( local, global )
!	===========================================  

	IMPLICIT NONE

	integer :: ierr, rec_req, send_req

	integer :: i, j, proc
	integer :: tot, m, n, i1, i2, j1, j2
	integer :: BLOCK2Dproc

	real, dimension(nx,ny) :: global
	real, dimension(ip_start:ip_end,jp_start:jp_end) :: local
!
! -------------------------------------------------------
!
!  Master proc: send data
!
    m = it_end-it_start+1
    n = jt_end-jt_start+1
    j1 = jp_start
    j2 = jp_end
!
    if(mypid == 0) then          
      local(ip_start:ip_end,jp_start:jp_end) = global(ip_start:ip_end,jp_start:jp_end)
!
      do i = 0, nprocx-1
        do j = 0, nprocy-1		
          proc = j+i*nprocy			     
 	  if (proc /= 0) then
 	    i1 = 2+i*m
 	    i2 = i1+m+3
 	    if (i == 0) i1 = 1
 	    if (i == nprocx-1) i2 = nx
#if (defined MODEL_3D) && (defined DECOMP_2D)
 	    j1 = 2+j*n
 	    j2 = j1+n+3
 	    if ( j == 0 ) j1 = 1
 	    if ( j == nprocy-1 ) j2 = ny
#endif
	    tot = (i2-i1+1)*(j2-j1+1)
!
	    call MPI_TYPE_CREATE_SUBARRAY (2, (/nx,ny/), (/i2-i1+1,j2-j1+1/), (/i1-1,j1-1/), 	&
		    MPI_ORDER_FORTRAN, REALTYPE, BLOCK2Dproc, ierr)
!
	    call MPI_TYPE_COMMIT (BLOCK2Dproc, ierr)
!
            call MPI_ISend( global, 1, BLOCK2Dproc, proc, proc+1, MPI_COMM_WORLD, send_req, ierr) 
!
	    call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
!
	    call MPI_TYPE_FREE (BLOCK2Dproc, ierr)
          endif		     
        enddo
      enddo
!
!  Slave procs: receive data (non-blocking instruction)
!
    else
      call MPI_IRecv( local, 1, BLOCK2Dp, 0, MPI_ANY_TAG, MPI_COMM_WORLD, rec_req, ierr )  
!
      call MPI_Wait( rec_req, MPI_STATUS_IGNORE, ierr )
    endif
!
    return
    end subroutine distribute_2d
!  
#endif
!
  end module mpicomm
