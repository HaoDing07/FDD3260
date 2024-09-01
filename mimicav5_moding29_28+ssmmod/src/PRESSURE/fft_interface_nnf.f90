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
  
  module fft_solvers_nnf
        
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

  public :: fft_solve_nnf_serial, fft_solve_nnf_parallel

  contains

! ===============================================
  subroutine fft_solve_nnf_parallel ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: tend
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: output
!
#if (defined SPMD) 
!
  integer :: i, j, k, nxp, nyp, nzp, nxa
  real, dimension(:,:,:), allocatable :: tend_shuff
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
  include 'fftw3-mpi.f03'
!
  if (verbose > 0) call write_debug('Starting fft_solve_parallel')
!
!  Initialize and allocate
!
  nxp = (nx-5)/(nprocx*nprocy)
  nxa = nx-5
  nyp = 2*(ny-5)
  nzp = nz
!
!-----------------------------------------------------------!
!     		Initialize and define plans                 !
!-----------------------------------------------------------!
!
  allocate ( tend_shuff(1:nxp,1:nyp,1:nzp), tendfft(1:nxp,1:nyp,1:nzp) )
!
  !call duplicate (nxp, nyp, nzp, tend, tend_shuff, 0)
!
!  Allocate FFTW pointers
! 
  cdata = fftw_alloc_complex ( int(nyp * nxp, C_SIZE_T) )
!
  call c_f_pointer ( cdata, tendfft2d, [int(nyp,C_INTPTR_T),int(nxp,C_INTPTR_T)] )
!
!-----------------------------------------------------------!
!     		        Forward FFT              	    !
!-----------------------------------------------------------!
!
!  Create 1D plans for FFTs
!
  fft_plan_forward = fftw_mpi_plan_dft_2d ( int(nxa,C_INTPTR_T), int(nyp,C_INTPTR_T), tendfft2d, tendfft2d,     &
        MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE )
!
  do k = 1, nzp
!
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
!    
  enddo
!
  call fftw_destroy_plan ( fft_plan_forward )
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Start tridiag')
  do j = 1,nyp
    call tridiag ( j, nxp, nzp, wv(:,j), tendfft(:,j,:) )  
  enddo
  if (verbose > 2) call write_debug('Terminate tridiag')
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
!  Create 1D plans for FFTs
!
  fft_plan_backward = fftw_mpi_plan_dft_2d ( int(nxa,C_INTPTR_T), int(nyp,C_INTPTR_T), tendfft2d, tendfft2d,    &
        MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE )
!
  do k = 1, nzp
!
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
!
  enddo
!
  call fftw_destroy_plan ( fft_plan_backward )
!
!  Reshuffle distributed data on original parallel layout
!
#ifdef DECOMP_2D
  call shuffle ( nxp, ny-5, nzp, output, tend_shuff(:,1:ny-5,:), 1 )
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
!  Terminating
!
  call fftw_free ( cdata )
!
  deallocate ( tend_shuff, tendfft )
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating fft_solve_nnf_parallel')
!
#endif
!
  end subroutine fft_solve_nnf_parallel
!          
! ===============================================
  subroutine fft_solve_nnf_serial ( tend, output )
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
  real, dimension(:,:,:), allocatable :: tend_shuff
!
!-----------------------------------------------------------!
!
  include 'fftw3.f03'
!
  if (verbose > 0) call write_debug('Starting fft_solve_nnf_serial')
!
!  Initialize and allocate
!
  nxp = nx - 5 
  nyp = 2*(ny-5)
  nzp = nz
!
  allocate ( tendfft(1:nxp,1:nyp,1:nzp), tend_shuff(1:nxp,1:nyp,1:nzp) )
!
  !call duplicate (nxp, nyp, nzp, tend, tend_shuff, 0)
!
  tendfft = cmplx( tend_shuff, 0. )
!
!-----------------------------------------------------------!
!     		        Forward FFT               	    !
!-----------------------------------------------------------!
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
!
!  Store output and deallocate
!
  output(it_start:it_end,jt_start:jt_end,1:nzp) = real( tendfft(:,1:jt_end-jt_start+1,:) )
!
  deallocate ( tendfft, tend_shuff )
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating fft_solve_nnf_serial')
!
  end subroutine fft_solve_nnf_serial

end module fft_solvers_nnf
