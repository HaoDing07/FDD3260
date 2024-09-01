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

! ===============================================
  subroutine fft_solve_parallel_mpi ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
#ifdef SPMD
USE mpi
USE gridno
USE mpicomm
USE shared_data
USE shared_pressure
USE, intrinsic :: iso_c_binding 

IMPLICIT NONE
!
  integer :: i, j, k, nxp, nyp, nzp, ioff, joff
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: tend, output
  real :: out_tot
!
  integer(C_INTPTR_T) :: alloc_local, offset, offset2, local_ny, local_nx
  TYPE(C_PTR) :: fft_plan_forward, fft_plan_backward, fftw_plan_transpose, cdata
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
!  Initialize and allocate
!
  nxp = it_end-it_start+1
#ifdef MODEL_3D
  nyp = jt_end-jt_start+1
#else
  nyp = 1
#endif
  nzp = nz
!  
  if (mypid /= 0) ioff = it_start-4
  joff = 0
#if (defined DECOMP_2D)
  if (mypid /= 0) joff = jt_start-4
#endif
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
!     		Initialize and define plans                 !
!-----------------------------------------------------------!
!
  allocate ( tendfft(1:nyp,1:nxp,1:nzp) )
!
!  Allocate FFTW pointers
!
  alloc_local = fftw_mpi_local_size_2d ( int(ny-5,C_INTPTR_T), int(nx-5,C_INTPTR_T), MPI_COMM_WORLD, local_ny, offset )
!
  cdata = fftw_alloc_complex ( alloc_local )
!
  call c_f_pointer ( cdata, tendfft2d, [int(nx-5,C_INTPTR_T),local_ny] )
!
!  Create 1D plans for FFTs
!
  fft_plan_forward = fftw_mpi_plan_dft_2d ( int(ny-5,C_INTPTR_T), int(nx-5,C_INTPTR_T), tendfft2d, tendfft2d,     &
        MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE )
!
  fft_plan_backward = fftw_mpi_plan_dft_2d ( int(ny-5,C_INTPTR_T), int(nx-5,C_INTPTR_T), tendfft2d, tendfft2d,    &
        MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE )
!
!-----------------------------------------------------------!
!     		        Forward FFT              	    !
!-----------------------------------------------------------!
!
  do k = 1, nzp
    do i = it_start,it_end
      do j = jt_start,jt_end
        tendfft2d(j-jt_start+1,i-it_start+1) = cmplx( tend(i,j,k),0. )
      enddo
    enddo
!
    call fftw_mpi_execute_dft ( fft_plan_forward, tendfft2d, tendfft2d )
!
    tendfft(:,:,k) = tendfft2d
  enddo
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  call tridiag ( joff, ioff, nyp, nxp, nzp, tendfft )  
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
  do k = 1, nzp
    tendfft2d = tendfft(:,:,k)
!
    call fftw_mpi_execute_dft ( fft_plan_backward, tendfft2d, tendfft2d )
!
    do i = it_start,it_end
      do j = jt_start,jt_end
        output(i,j,k) = real( tendfft2d(j-jt_start+1,i-it_start+1) )
      enddo
    enddo
  enddo
!
!  Discard the information associated with plans.
!
  call fftw_destroy_plan ( fft_plan_forward )
!
  call fftw_destroy_plan ( fft_plan_backward )
!
  call fftw_free ( cdata )
!
#else
!
!-----------------------------------------------------------!
!  			 2D case			    !
!-----------------------------------------------------------!
!
  allocate ( tendfft(1:nxp,1,1:nzp) )
!
!  Allocate FFTW pointers
!
  alloc_local = fftw_mpi_local_size_1d ( int(nx-5,C_INTPTR_T), MPI_COMM_WORLD,	&
        FFTW_FORWARD, FFTW_ESTIMATE, local_nx, offset, local_ny, offset2 )
!
  cdata = fftw_alloc_complex ( alloc_local )
!
  call c_f_pointer ( cdata, tendfft1d, [local_nx] )
!
!  Create 1D plans for FFTs
!
  fft_plan_forward = fftw_mpi_plan_dft_1d ( int(nx-5,C_INTPTR_T), tendfft1d, tendfft1d,     &
        MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE )
!
  fft_plan_backward = fftw_mpi_plan_dft_1d ( int(nx-5,C_INTPTR_T), tendfft1d, tendfft1d,    &
        MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE )
!
!-----------------------------------------------------------!
!     		        Forward FFT              	    !
!-----------------------------------------------------------!
!
  do k = 1, nzp
    do i = it_start,it_end
      tendfft1d(i-it_start+1) = cmplx( tend(i,1,k),0. )
    enddo
!
    call fftw_mpi_execute_dft ( fft_plan_forward, tendfft1d, tendfft1d )
!
    tendfft(:,1,k) = tendfft1d
  enddo
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  call tridiag ( ioff, 0, nxp, nyp, nzp, tendfft )  
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
  do k = 1, nzp
    tendfft1d = tendfft(:,1,k)
!
    call fftw_mpi_execute_dft ( fft_plan_backward, tendfft1d, tendfft1d )
!
    do i = it_start,it_end
      output(i,1,k) = real( tendfft1d(i-it_start+1) )
    enddo
  enddo
!
!  Discard the information associated with plans.
!
  call fftw_destroy_plan ( fft_plan_forward )
!
  call fftw_destroy_plan ( fft_plan_backward )
!
  call fftw_free ( cdata )
!
#endif
!
  deallocate (tendfft)
#endif
!
!-----------------------------------------------------------!
!
  end subroutine fft_solve_parallel_mpi
!
! ===============================================
  subroutine fft_solve_serial ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
USE gridno
USE shared_data
USE shared_pressure
USE, intrinsic :: iso_c_binding 

IMPLICIT NONE
!
  integer :: i, j, k, i1, i2, j1, j2, nxp, nyp, nzp
  real :: out_tot
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: tend, output
!
  integer(C_INTPTR_T) :: alloc_local, local_nx, local_ny
  TYPE(c_ptr) :: fft_plan_forward, fft_plan_backward, cdatai, cdatao
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: tendfft, outputfft
!
!-----------------------------------------------------------!
!
  include 'fftw3.f03'
!
!  Initialize and allocate
!
  i1 = it_start
  i2 = it_end
  j1 = 1
  j2 = 1    
#ifdef MODEL_3D
  j1 = jt_start
  j2 = jt_end
#endif
!
  nxp = i2-i1+1
  nyp = j2-j1+1
  nzp = nz
!
  allocate ( tendfft(1:nxp,1:nyp,1:nzp), outputfft(1:nxp,1:nyp,1:nzp) )
!
!  Convert to complex
!
  tendfft = cmplx(tend(i1:i2,j1:j2,1:nzp),0.)
!
!-----------------------------------------------------------!
!     	            Initialize plans               	    !
!-----------------------------------------------------------!
!
#ifdef MODEL_3D
  fft_plan_forward = fftw_plan_dft_2d ( int(nyp,C_INT), int(nxp,C_INT), tendfft(:,:,1), outputfft(:,:,1),   &
       FFTW_FORWARD, FFTW_ESTIMATE )
!
  fft_plan_backward = fftw_plan_dft_2d ( int(nyp,C_INT), int(nxp,C_INT), outputfft(:,:,1), outputfft(:,:,1),   &
       FFTW_BACKWARD, FFTW_ESTIMATE )
#else
  fft_plan_forward = fftw_plan_dft_1d ( int(nxp,C_INT), tendfft(:,1,1), outputfft(:,1,1),   &
       FFTW_FORWARD, FFTW_ESTIMATE )
!
  fft_plan_backward = fftw_plan_dft_1d ( int(nxp,C_INT), outputfft(:,1,1), outputfft(:,1,1),   &
       FFTW_BACKWARD, FFTW_ESTIMATE )
#endif
!
!-----------------------------------------------------------!
!     		        Forward FFT               	    !
!-----------------------------------------------------------!
!
  do k = 1, nzp
#ifdef MODEL_3D
    call fftw_execute_dft ( fft_plan_forward, tendfft(:,:,k), outputfft(:,:,k) )
#else
    call fftw_execute_dft ( fft_plan_forward, tendfft(:,1,k), outputfft(:,1,k) )
#endif
  enddo
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  call tridiag ( 0, 0, nxp, nyp, nzp, outputfft )  
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
  do k = 1, nzp
#ifdef MODEL_3D
    call fftw_execute_dft ( fft_plan_backward, outputfft(:,:,k), outputfft(:,:,k) )
#else
    call fftw_execute_dft ( fft_plan_backward, outputfft(:,1,k), outputfft(:,1,k) )
#endif
  enddo
!
!-----------------------------------------------------------!
!     		         Terminate               	    !
!-----------------------------------------------------------!
!
  call fftw_destroy_plan ( fft_plan_forward )
!
  call fftw_destroy_plan ( fft_plan_backward )
!
!  Store output and deallocate
!
  output(i1:i2,j1:j2,1:nzp) = real(outputfft)
!
  deallocate (outputfft, tendfft)
!
!-----------------------------------------------------------!
!
  end subroutine fft_solve_serial

