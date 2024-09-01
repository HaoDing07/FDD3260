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
  
  module fft_solvers
        
!
!-----------------------------------------------------------!
!

  USE gridno
  USE shared_data
  USE shared_diag
  USE fft_tools
  USE, intrinsic :: iso_c_binding 

#ifdef SPMD
  USE mpi
  USE mpicomm
#endif  

  IMPLICIT NONE

  private

  public :: fft_solve_serial, fft_solve_parallel, fft_solve_parallel_2x1d

  contains

! ===============================================
  subroutine fft_solve_parallel ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: tend
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: output
  real, dimension(:,:,:), allocatable :: tend_shuff
  integer :: ierr
!
#if (defined SPMD) 
!
  integer(C_INTPTR_T) :: i, j, k, alloc_local, nxp, nyp, nzp, nxa, nya, i_offset
!
  TYPE(C_PTR) :: fft_plan_forward, fft_plan_backward, cdata
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: tendfft
#ifdef MODEL_3D
  complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: tendfft2d
#else
  complex(C_DOUBLE_COMPLEX), dimension(:), pointer :: tendfft1d
#endif
!
!-----------------------------------------------------------!
!
  include 'fftw3-mpi.f03'
!
  if (verbose > 0) call write_debug('Start fft_solve_parallel')
!
!  Initialize and allocate
!
!  nxp = int((nx-5)/(nprocx*nprocy),C_INTPTR_T)
  nxa = int(nx-5,C_INTPTR_T)
  nya = int(ny-5,C_INTPTR_T)
  nyp = int(ny-5,C_INTPTR_T)
  nzp = int(nz,C_INTPTR_T)
!
!-----------------------------------------------------------!
!                 Start with parallel FFT              	    !
!  In parallel, the 2D FFT is done as follows:		    !
!	(1) 1D FFT along y (non split direction)	    !
!	(2) Transpose the transformed plane	 	    !
!	(3) 1D FFT again along y (now transposed)	    !
!	(4) Solve tridiag system in the vertical	    !
!	(5) 1D backward FFT along y			    !
!	(6) Transpose again to retrieve original matrix     !
!	(7) 1D backward FFT along y			    !
!							    !
!  In parallel, the 1D FFT is simply done after gathering   !
!	all data on master proc				    !	
!-----------------------------------------------------------!
!
#ifdef MODEL_3D
!
!-----------------------------------------------------------!
!   Shuffle distributed data to get 1D structure along y    !
!-----------------------------------------------------------!
!
  alloc_local = fftw_mpi_local_size_2d( int(nxa,C_INTPTR_T), int(nya,C_INTPTR_T), MPI_COMM_WORLD, nxp, i_offset )
!
  allocate ( tend_shuff(1:nxp,1:nyp,1:nzp), tendfft(1:nxp,1:nyp,1:nzp) )
!
#ifdef DECOMP_2D
  call shuffle ( nxp, nyp, nzp, tend, tend_shuff, 0 )
#else
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        tend_shuff(i-it_start+1,j-jt_start+1,k) = tend(i,j,k)
      enddo
    enddo
  enddo
#endif
!
!-----------------------------------------------------------!
!     		        Forward FFT              	    !
!-----------------------------------------------------------!
!
!  Allocate FFTW pointers
! 
  cdata = fftw_alloc_complex ( alloc_local )
!
  call c_f_pointer ( cdata, tendfft2d, [int(nyp,C_INTPTR_T),int(nxp,C_INTPTR_T)] )
!
!  Create 1D plans for FFTs
!
  if (verbose > 2) call write_debug('Before create forward plan')
  fft_plan_forward = fftw_mpi_plan_dft_2d ( int(nxa,C_INTPTR_T), int(nya,C_INTPTR_T), tendfft2d, tendfft2d,     &
        MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE )
!
  if (verbose > 2) call write_debug('Before execute forward plan')
  do k = 1, nzp
    do j = 1, nyp
      do i = 1, nxp
        tendfft2d(j,i) = cmplx( tend_shuff(i,j,k), 0. )
      enddo
    enddo
!
!  Execute FFTs
!
    call fftw_mpi_execute_dft ( fft_plan_forward, tendfft2d, tendfft2d )
!
    do j = 1,nyp
      do i = 1,nxp
        tendfft(i,j,k) = tendfft2d(j,i)
      enddo
    enddo
  enddo
!
  call fftw_destroy_plan ( fft_plan_forward )
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Start tridiag')
  do j = 1,nyp
    call tridiag ( j, nxp, nzp, wv(:,j), tendfft(:,j,:) )  
  enddo
  if (verbose > 1) call write_debug('Terminate tridiag')
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
!  Create 1D plans for FFTs
!
  if (verbose > 2) call write_debug('Before create backward plan')
  fft_plan_backward = fftw_mpi_plan_dft_2d ( int(nxa,C_INTPTR_T), int(nya,C_INTPTR_T), tendfft2d, tendfft2d,    &
        MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE )
!
  if (verbose > 2) call write_debug('Before execute backward plan')
  do k = 1, nzp
    do j = 1,nyp
      do i = 1,nxp
        tendfft2d(j,i) = tendfft(i,j,k)
      enddo
    enddo
!
!  Execute (backward) FFTs
!
    call fftw_mpi_execute_dft ( fft_plan_backward, tendfft2d, tendfft2d )
!
    do j = 1, nyp
      do i = 1, nxp
        tend_shuff(i,j,k) = real( tendfft2d(j,i) )
      enddo
    enddo
  enddo
!
  call fftw_destroy_plan ( fft_plan_backward )
!
  call fftw_free ( cdata )
!
!-----------------------------------------------------------!
!  Reshuffle distributed data on original parallel layout   !
!-----------------------------------------------------------!
!
#ifdef DECOMP_2D
  call shuffle ( nxp, nyp, nzp, output, tend_shuff, 1 )
#else
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        output(i,j,k) = tend_shuff(i-it_start+1,j-jt_start+1,k)
      enddo
    enddo
  enddo
#endif
!
  deallocate ( tend_shuff, tendfft )
!
!-----------------------------------------------------------!
!  			 2D case			    !
!-----------------------------------------------------------!
!
#else
!
  allocate ( tend_shuff(1:nxp,1,1:nzp), tendfft(1:nxp,1,1:nzp) )
!
  do k = 1, nz
    do i = it_start, it_end
      tend_shuff(i-it_start+1,1,k) = tend(i,1,k)
    enddo
  enddo
!
!-----------------------------------------------------------!
!     		        Forward FFT              	    !
!-----------------------------------------------------------!
!
!  Allocate FFTW pointers
!
  cdata = fftw_alloc_complex ( int(nxp, C_SIZE_T) )
!
  call c_f_pointer ( cdata, tendfft1d, (/ int(nxp,C_INTPTR_T) /) )
!
!  Create 1D plans for FFTs
!
  fft_plan_forward = fftw_mpi_plan_dft_1d ( int(nxa,C_INTPTR_T), tendfft1d, tendfft1d,     &
        MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE )
!
  do k = 1, nzp
    do i = 1,nxp
      tendfft1d(i) = cmplx( tend_shuff(i,1,k),0. )
    enddo
!
!  Execute FFTs
!
    call fftw_mpi_execute_dft ( fft_plan_forward, tendfft1d, tendfft1d )
!
    tendfft(:,1,k) = tendfft1d
  enddo
!
  call fftw_destroy_plan ( fft_plan_forward )
!
  call fftw_free ( cdata )
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Start tridiag')
  call tridiag ( 1, nxp, nzp, wv(:,1), tendfft(:,1,:) )  
  if (verbose > 2) call write_debug('Start tridiag')
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
!  Allocate FFTW pointers
!
  cdata = fftw_alloc_complex ( int(nxp, C_SIZE_T) )
!
  call c_f_pointer ( cdata, tendfft1d, (/ int(nxp,C_INTPTR_T) /) )
!
!  Create 1D plans for FFTs
!
  fft_plan_backward = fftw_mpi_plan_dft_1d ( int(nxa,C_INTPTR_T), tendfft1d, tendfft1d,    &
        MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE )
!
  do k = 1, nzp
    tendfft1d = tendfft(:,1,k)
!
!  Execute (backward) FFTs
!
    call fftw_mpi_execute_dft ( fft_plan_backward, tendfft1d, tendfft1d )
!
    do i = 1, nxp
      tend_shuff(i,1,k) = real( tendfft1d(i) )
    enddo
  enddo
!
  call fftw_destroy_plan ( fft_plan_backward )
!
  call fftw_free ( cdata )
!
!  Output array
!
  do k = 1, nz
    do i = it_start, it_end
      output(i,1,k) = tend_shuff(i-it_start+1,1,k)
    enddo
  enddo
!
  deallocate ( tend_shuff, tendfft )
!
#endif
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate fft_solve_parallel')
!
#endif
!
  end subroutine fft_solve_parallel
!
! ===============================================
  subroutine fft_solve_parallel_2x1d ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: tend
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: output
!
#if (defined SPMD) && (defined MODEL_3D)
!
  integer :: i, j, k, nxp, nyp, nzp, nxa, nya, ierr
  real, dimension(:,:,:), allocatable :: tend_shuff
!
  TYPE(C_PTR) :: cdatax, cdatay
  TYPE(C_PTR) :: fft_plan_forward_x, fft_plan_backward_x, fft_plan_forward_y, fft_plan_backward_y
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: tendfft
  complex(C_DOUBLE_COMPLEX), dimension(:), pointer :: tendfft1dx, tendfft1dy
!
!-----------------------------------------------------------!
!
  include 'fftw3-mpi.f03'
!
  if (verbose > 0) call write_debug('Start fft_solve_parallel_2x1d')
!
!  Initialize and allocate
!
  nxp = (nx-5)/(nprocx*nprocy)
  nyp = ny-5
  nxa = nx-5
  nya = ny-5
  nzp = nz
!
!-----------------------------------------------------------!
!                 Start with parallel FFT              	    !
!  In parallel, the 2D FFT is done as follows:		    !
!	(1) 1D FFT along y (non split direction)	    !
!	(2) Transpose the transformed plane	 	    !
!	(3) 1D FFT again along y (now transposed)	    !
!	(4) Solve tridiag system in the vertical	    !
!	(5) 1D backward FFT along y			    !
!	(6) Transpose again to retrieve original matrix     !
!	(7) 1D backward FFT along y			    !
!							    !
!  In parallel, the 1D FFT is simply done after gathering   !
!	all data on master proc				    !	
!-----------------------------------------------------------!
!
!-----------------------------------------------------------!
!     		Initialize and define plans                 !
!-----------------------------------------------------------!
!
!  Shuffle distributed data to get 1D structure along y
!
  allocate ( tend_shuff(1:nxp,1:nyp,1:nzp), tendfft(1:nxp,1:nyp,1:nzp) )
!
#ifdef DECOMP_2D
  call shuffle ( nxp, nyp, nzp, tend, tend_shuff, 0 )
#else
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        tend_shuff(i-it_start+1,j-jt_start+1,k) = tend(i,j,k)
      enddo
    enddo
  enddo
#endif
!
!-----------------------------------------------------------!
!     		        Forward FFT              	    !
!-----------------------------------------------------------!
!
!  Allocate FFTW pointers
! 
  cdatax = fftw_alloc_complex ( int(nxp, C_SIZE_T) )
  cdatay = fftw_alloc_complex ( int(nyp, C_SIZE_T) )
!
  call c_f_pointer ( cdatax, tendfft1dx, (/ int(nxp,C_INTPTR_T) /) )
  call c_f_pointer ( cdatay, tendfft1dy, (/ int(nyp,C_INTPTR_T) /) )
!
!  Create 1D plans for FFTs
!
  if (verbose > 1) call write_debug('Before create forward plan x')
  fft_plan_forward_x = fftw_mpi_plan_dft_1d ( int(nxa,C_INTPTR_T), tendfft1dx, tendfft1dx,     &
          MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE )
!
  if (verbose > 1) call write_debug('Before create forward plan y')
  fft_plan_forward_y = fftw_plan_dft_1d ( int(nya,C_INT), tendfft1dy, tendfft1dy, FFTW_FORWARD, FFTW_ESTIMATE )
!
  if (verbose > 1) call write_debug('Before execute forward plan')
  do k = 1, nzp
!
!  Execute 1D serial FFTs along y
!
    do i = 1, nxp
      do j = 1, nyp
        tendfft1dy(j) = cmplx( tend_shuff(i,j,k), 0. )
      enddo
!
      call fftw_execute_dft ( fft_plan_forward_y, tendfft1dy, tendfft1dy )
!
      do j = 1, nyp
        tendfft(i,j,k) = tendfft1dy(j)
      enddo
    enddo
!
!  Execute 1D parallel FFTs along x
!
    do j = 1, nyp
      do i = 1, nxp
        tendfft1dx(i) = tendfft(i,j,k)
      enddo
!
      call fftw_mpi_execute_dft ( fft_plan_forward_x, tendfft1dx, tendfft1dx )
!
      do i = 1,nxp
        tendfft(i,j,k) = tendfft1dx(i)
      enddo
    enddo
  enddo
!
!  Deallocate plans and pointers
!
  call fftw_destroy_plan ( fft_plan_forward_x )
  call fftw_destroy_plan ( fft_plan_forward_y )
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Start tridiag')
  do j = 1,nyp
    call tridiag ( j, nxp, nzp, wv(:,j), tendfft(:,j,:) )  
  enddo
  if (verbose > 1) call write_debug('Terminate tridiag')
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
!  Create 1D plans for FFTs
!
  if (verbose > 2) call write_debug('Before create backward plan x')
  fft_plan_backward_x = fftw_mpi_plan_dft_1d ( int(nxa,C_INTPTR_T), tendfft1dx, tendfft1dx,     &
          MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE )
!
  if (verbose > 2) call write_debug('Before create backward plan y')
  fft_plan_backward_y = fftw_plan_dft_1d ( int(nya,C_INT), tendfft1dy, tendfft1dy, FFTW_BACKWARD, FFTW_ESTIMATE )
!
  if (verbose > 2) call write_debug('Before execute backward plan')
  do k = 1, nzp
!
!  Execute 1D serial FFTs along y
!
    do i = 1, nxp
      do j = 1, nyp
        tendfft1dy(j) = tendfft(i,j,k)
      enddo
!
      call fftw_execute_dft ( fft_plan_backward_y, tendfft1dy, tendfft1dy )
!
      do j = 1, nyp
        tendfft(i,j,k) = tendfft1dy(j)
      enddo
    enddo
!
!  Execute 1D parallel FFTs along x
!
    do j = 1, nyp
      do i = 1, nxp
        tendfft1dx(i) = tendfft(i,j,k)
      enddo
!
      call fftw_mpi_execute_dft ( fft_plan_backward_x, tendfft1dx, tendfft1dx )
!
      do i = 1,nxp
        tend_shuff(i,j,k) = real( tendfft1dx(i) )
      enddo
    enddo
  enddo
!
!  Deallocate plans and pointers
!
  call fftw_destroy_plan ( fft_plan_backward_x )
  call fftw_destroy_plan ( fft_plan_backward_y )
!
  call fftw_free ( cdatax )
  call fftw_free ( cdatay )
!
!-----------------------------------------------------------!
!     		        Terminating              	    !
!-----------------------------------------------------------!
!
!  Reshuffle distributed data on original parallel layout
!
#ifdef DECOMP_2D
  call shuffle ( nxp, nyp, nzp, output, tend_shuff, 1 )
#else
  do k = 1, nz
    do j = jt_start, jt_end
      do i = it_start, it_end
        output(i,j,k) = tend_shuff(i-it_start+1,j-jt_start+1,k)
      enddo
    enddo
  enddo
#endif
!
  deallocate ( tendfft, tend_shuff )
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate fft_solve_parallel_2x1d')
!
#endif
!
  end subroutine fft_solve_parallel_2x1d
!          
! ===============================================
  subroutine fft_solve_serial ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
!
  integer :: j, k, nxp, nyp, nzp
  real, dimension(ip_start:,jp_start:,:), intent(in) :: tend
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: output
!
  TYPE(c_ptr) :: fft_plan_forward, fft_plan_backward
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: tendfft
!
!-----------------------------------------------------------!
!
  include 'fftw3.f03'
!
  if (verbose > 0) call write_debug('Starting fft_solve_serial')
!
!  Initialize and allocate
!
  nxp = it_end-it_start+1
  nyp = jt_end-jt_start+1
  nzp = nz
!
  allocate ( tendfft(1:nxp,1:nyp,1:nzp) )
!
  tendfft = cmplx(tend(it_start:it_end,jt_start:jt_end,1:nzp),0.)
!
!-----------------------------------------------------------!
!     		        Forward FFT               	    !
!-----------------------------------------------------------!
!
#ifdef MODEL_3D
!
!  Create 2D plan for forward FFT
!
  fft_plan_forward = fftw_plan_dft_2d ( int(nyp,C_INT), int(nxp,C_INT), tendfft(:,:,1), tendfft(:,:,1),   &
       FFTW_FORWARD, FFTW_ESTIMATE )
!
!  Perform forward 2D FFT
!
  do k = 1, nzp
    call fftw_execute_dft ( fft_plan_forward, tendfft(:,:,k), tendfft(:,:,k) )
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_forward )
#else
!
!  Create 1D plan for forward FFT
!
  fft_plan_forward = fftw_plan_dft_1d ( int(nxp,C_INT), tendfft(:,1,1), tendfft(:,1,1),   &
       FFTW_FORWARD, FFTW_ESTIMATE )
!
!  Perform forward 1D FFT
!
  do k = 1, nzp
    call fftw_execute_dft ( fft_plan_forward, tendfft(:,1,k), tendfft(:,1,k) )
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_forward )
#endif
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  do j = 1, nyp
    call tridiag ( j, nxp, nzp, wv(:,j), tendfft(:,j,:) )  
  enddo
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
#ifdef MODEL_3D
!
!  Create 2D plan for backward FFT
!
  fft_plan_backward = fftw_plan_dft_2d ( int(nyp,C_INT), int(nxp,C_INT), tendfft(:,:,1), tendfft(:,:,1),   &
         FFTW_BACKWARD, FFTW_ESTIMATE )
!
!  Perform backward 2D FFT
!
  do k = 1, nzp
    call fftw_execute_dft ( fft_plan_backward, tendfft(:,:,k), tendfft(:,:,k) )
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_backward )
#else
!
!  Create 1D plan for backward FFT
!
  fft_plan_backward = fftw_plan_dft_1d ( int(nxp,C_INT), tendfft(:,1,1), tendfft(:,1,1),   &
         FFTW_BACKWARD, FFTW_ESTIMATE )
! 
!  Perform backward 1D FFT
!
  do k = 1, nzp
    call fftw_execute_dft ( fft_plan_backward, tendfft(:,1,k), tendfft(:,1,k) )
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_backward )
#endif
!
!  Store output and deallocate
!
  output(it_start:it_end,jt_start:jt_end,1:nzp) = real( tendfft )
!
  deallocate ( tendfft )
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating fft_solve_serial')
!
  end subroutine fft_solve_serial

end module fft_solvers
