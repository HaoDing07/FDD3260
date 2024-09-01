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
  
  module fft_solvers_open
        
!
!-----------------------------------------------------------!
!

  USE gridno
  USE shared_data
  USE shared_diag
  USE shared_pressure
  USE fft_tools
  USE, intrinsic :: iso_c_binding 

#ifdef SPMD
  USE mpi
  USE mpicomm
#endif  

  IMPLICIT NONE

  private

  public :: fft_solve_parallel_open, fft_solve_serial_open

  contains

! ===============================================
  subroutine fft_solve_parallel_open ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
  integer :: i, j, k, i1, i2, j1, j2, nxp, nxp2, nyp, nzp
  real :: out_tot
  real, dimension(ip_start:,jp_start:,:), intent(in) :: tend
  real, dimension(ip_start:,jp_start:,:), intent(out) :: output
!
#if (defined SPMD) && (defined NESTING)
!
  integer(C_INTPTR_T) :: alloc_local, local_nx, local_ny
  TYPE(c_ptr) :: fft_plan_forward, fft_plan_backward, cdatai, cdatao
  real(C_DOUBLE), dimension(:,:,:), allocatable :: tendfft, outputfft
!
!-----------------------------------------------------------!
!
  include 'fftw3.f03'
!
!  Initialize and allocate
!
#ifdef MODEL_3D
  i1 = max(it_start,3)
  i2 = min(it_end,nx-2)
  j1 = max(jt_start,3)
  j2 = min(jt_end,ny-2)
!
  nxp = i2-i1+1
  nyp = j2-j1+1
  nzp = nz
#else
  i1 = max(it_start,3)
  i2 = min(it_end,nx-2)
  j1 = 1
  j2 = 1    
!
  nxp2 = nx-4
  nxp = i2-i1+1
  nyp = 1
  nzp = nz
!
  allocate ( tendfft(1:nxp,1:nyp,1:nzp) )
!
  tendfft = tend(i1:i2,j1:j2,1:nzp)
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
!     		        Forward FFT              	    !
!-----------------------------------------------------------!
!
  allocate ( outputfft(1:nxp,1:nyp,1:nzp) )
!
!  Create 1D plan for forward FFT
!
  fft_plan_forward = fftw_plan_r2r_1d ( int(nyp,C_INT), tendfft(1,:,1), outputfft(1,:,1),  &
        FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform forward 2D FFT along y
!
  do k = 1, nzp
    do i = 1, nxp
      call fftw_execute_r2r ( fft_plan_forward, tendfft(i,:,k), outputfft(i,:,k) )
    enddo
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_forward )
!
!  Transpose output matrix
!
  call mpi_transpose_open (nxp,nyp,nzp,outputfft)
!
!  Create 1D plan for forward FFT on transpose
!
  fft_plan_forward = fftw_plan_r2r_1d ( int(nyp,C_INT), outputfft(1,:,1), outputfft(1,:,1),  &
        FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform forward 2D FFT along y
!
  do k = 1, nzp
    do i = 1, nxp
      call fftw_execute_r2r ( fft_plan_forward, outputfft(i,:,k), outputfft(i,:,k) )
    enddo
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
  call tridiag_open ( nxp, nyp, nzp, outputfft )  
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
!  Create 1D plan for forward FFT
!
  fft_plan_backward = fftw_plan_r2r_1d ( int(nyp,C_INT), outputfft(1,:,1), outputfft(1,:,1),  &
        FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform forward 2D FFT along y
!
  do k = 1, nzp
    do i = 1, nxp
      call fftw_execute_r2r ( fft_plan_backward, outputfft(i,:,k), outputfft(i,:,k) )
    enddo
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_backward )
!
!  Transpose output matrix
!
  call mpi_transpose_open (nxp,nyp,nzp,outputfft)
!
!  Create 1D plan for forward FFT on transpose
!
  fft_plan_backward = fftw_plan_r2r_1d ( int(nyp,C_INT), outputfft(1,:,1), outputfft(1,:,1),  &
      FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform forward 2D FFT along y
!
  do k = 1, nzp
    do i = 1, nxp
      call fftw_execute_r2r ( fft_plan_backward, outputfft(i,:,k), outputfft(i,:,k) )
    enddo
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_backward )
!
!  Store output and deallocate
!
  output(i1:i2,j1:j2,1:nzp) = outputfft
!
  deallocate (outputfft)
#else
!
!-----------------------------------------------------------!
!     		           2D case               	    !
!-----------------------------------------------------------!
!
  allocate ( outputfft(1:nxp2,1,1:nzp) )
!
!  Collect data on masterproc
!
  call block_sendrec_real (nxp,nzp,nxp2,tendfft,outputfft,0)
!
  if ( mypid == 0 ) then
!
!  Create 1D plan for backward FFT
!
    fft_plan_forward = fftw_plan_r2r_1d ( int(nxp2,C_INT), outputfft(:,1,1), outputfft(:,1,1),   &
        FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform backward 1D FFT
!
    do k = 1, nzp
      call fftw_execute_r2r ( fft_plan_forward, outputfft(:,1,k), outputfft(:,1,k) )
    enddo
!
!  Discard the information associated with the plans.
!
    call fftw_destroy_plan ( fft_plan_forward )
!
!  Solve tridiag system
!
    call tridiag_open ( nxp2, 1, nzp, outputfft )  
!
!  Create 1D plan for backward FFT
!
    fft_plan_backward = fftw_plan_r2r_1d ( int(nxp2,C_INT), outputfft(:,1,1), outputfft(:,1,1),   &
        FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform backward 1D FFT
!
    do k = 1, nzp
      call fftw_execute_r2r ( fft_plan_backward, outputfft(:,1,k), outputfft(:,1,k) )
    enddo
!
!  Discard the information associated with the plans.
!
    call fftw_destroy_plan ( fft_plan_backward )
  endif
!
!  Redistribute
!
  call block_sendrec_real (nxp,nzp,nxp2,tendfft,outputfft,1)
!
!  Store output and deallocate
!
  output(i1:i2,1,1:nzp) = tendfft(1:nxp,1,1:nzp)
!
  deallocate (outputfft)
#endif
!
  deallocate (tendfft)
#endif
!
!-----------------------------------------------------------!
!
#endif
!
  end subroutine fft_solve_parallel_open
        
! ===============================================
  subroutine fft_solve_serial_open ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
  integer :: i, j, k, i1, i2, j1, j2, nxp, nyp, nzp
  real :: out_tot
  real, dimension(ip_start:,jp_start:,:), intent(in) :: tend
  real, dimension(ip_start:,jp_start:,:), intent(out) :: output
!
#ifdef NESTING
!
  integer(C_INTPTR_T) :: alloc_local, local_nx, local_ny
  TYPE(c_ptr) :: fft_plan_forward, fft_plan_backward, cdatai, cdatao
  real(C_DOUBLE), dimension(:,:,:), allocatable :: tendfft, outputfft
!
!-----------------------------------------------------------!
!
  include 'fftw3.f03'
!
!  Initialize and allocate
!
  i1 = max(it_start,3)
  i2 = min(it_end,nx-2)
  j1 = 1
  j2 = 1    
#ifdef MODEL_3D
  j1 = max(jt_start,3)
  j2 = min(jt_end,ny-2)
#endif
!
  nxp = i2-i1+1
  nyp = j2-j1+1
  nzp = nz
!
  allocate ( tendfft(1:nxp,1:nyp,1:nzp), outputfft(1:nxp,1:nyp,1:nzp) )
  tendfft = tend(i1:i2,j1:j2,1:nzp)
!
!-----------------------------------------------------------!
!     		        Forward FFT               	    !
!-----------------------------------------------------------!
!
#ifdef MODEL_3D
!
!  Create 2D plan for forward FFT
!
  fft_plan_forward = fftw_plan_r2r_2d ( int(nyp,C_INT), int(nxp,C_INT), tendfft(:,:,1), outputfft(:,:,1),   &
       FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform forward 2D FFT
!
  do k = 1, nzp
    call fftw_execute_r2r ( fft_plan_forward, tendfft(:,:,k), outputfft(:,:,k) )
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_forward )
#else
!
!  Create 1D plan for forward FFT
!
  fft_plan_forward = fftw_plan_r2r_1d ( int(nxp,C_INT), tendfft(:,1,1), outputfft(:,1,1),   &
       FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform forward 1D FFT
!
  do k = 1, nzp
    call fftw_execute_r2r ( fft_plan_forward, tendfft(:,1,k), outputfft(:,1,k) )
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
  call tridiag_open ( nxp, nyp, nzp, outputfft )  
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
#ifdef MODEL_3D
!
!  Create 2D plan for backward FFT
!
    fft_plan_backward = fftw_plan_r2r_2d ( int(nyp,C_INT), int(nxp,C_INT), outputfft(:,:,1), outputfft(:,:,1),   &
         FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform backward 2D FFT
!
  do k = 1, nzp
    call fftw_execute_r2r ( fft_plan_backward, outputfft(:,:,k), outputfft(:,:,k) )
  enddo
!
!  Discard the information associated with the plans.
!
    call fftw_destroy_plan ( fft_plan_backward )
#else
!
!  Create 1D plan for backward FFT
!
  fft_plan_backward = fftw_plan_r2r_1d ( int(nxp,C_INT), outputfft(:,1,1), outputfft(:,1,1),   &
       FFTW_RODFT00, FFTW_ESTIMATE )
!
!  Perform backward 1D FFT
!
  do k = 1, nzp
    call fftw_execute_r2r ( fft_plan_backward, outputfft(:,1,k), outputfft(:,1,k) )
  enddo
!
!  Discard the information associated with the plans.
!
  call fftw_destroy_plan ( fft_plan_backward )
#endif
!
!  Store output and deallocate
!
  output(i1:i2,j1:j2,1:nzp) = outputfft
!
  deallocate (outputfft, tendfft)
!
!-----------------------------------------------------------!
!
#endif
!
  end subroutine fft_solve_serial_open
!
! ===============================================
  subroutine tridiag_open ( nxp, nyp, nzp, s1 )  
! ===============================================

  integer :: nxp, nyp, nzp
  real(C_DOUBLE) :: s1(1:nxp,1:nyp,1:nzp)

  integer :: i,j,k,l,m,ioff,joff,ierr
  real    :: dzz,wm1,nl,nm,xl,xm,af,cf
  real, parameter :: small=1.e-15

  real :: ak(1:nxp,1:nzp),dk(1:nxp,1:nzp),bk(1:nxp,1:nzp),ck(1:nxp,1:nzp)
  real :: xk(1:nxp,1:nzp),yk(1:nxp,1:nzp),wv(1:nxp,1:nyp)
!
!  Defined dimensions and offsets
!
    if (bcl(5) == 'ope') s1(:,:,1) = 0.
    if (bcl(6) == 'ope') s1(:,:,nzp) = 0.
!
    nl = 2.*(float(nx-4)+1.)
    nm = 2.*(float(max(ny-4,1))+1.)
!
    dzz  = 1./dz
    ioff = 0
    joff = 0
#if (defined SPMD)
    if (mypid /= 0) ioff = it_start-3
#endif
#if (defined SPMD) && (defined DECOMP_2D)
    if (mypid /= 0) joff = jt_start-3
#endif
!
!  Calculate horizontal wave number
!
    do l=1,nxp
       xl=float(l+ioff-1)
!
       do m=1,nyp
#ifdef MODEL_3D
         xm=float(m+joff-1)
#else
         xm=0.0
#endif
!
         wv(l,m) = -4.*(sin(pi*xl/nl)*sin(pi*xl/nl)/(dx*dx) + sin(pi*xm/nm)*sin(pi*xm/nm)/(dy*dy))
       enddo
    enddo
!
! Coefficients for tri-diagonal solver: a_k x_k-1 + d_k x_k + c_k x_k+1 = b_k
!
    do m=1,nyp
!
      do k=2,nzp-1 
        do l=1,nxp
          ak(l,k) = dz1(k)*dz1w(k-1) - stab0(k)/(dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
          bk(l,k) = s1(l,m,k)
          ck(l,k) = dz1(k)*dz1w(k) + stab0(k)/(dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
          dk(l,k) = wv(l,m) - dz1(k)*(dz1w(k) + dz1w(k-1))
        enddo
      enddo
!
      if (bcl(5) == 'nnf') then
        ak(1:nxp,1) = 0.
        bk(1:nxp,1) = s1(1:nxp,m,1)
        ck(1:nxp,1) = dz1(1)*dz1w(1) + 0.5*stab0(1)*dz1w(1)
        dk(1:nxp,1) = wv(1:nxp,m) - dz1(1)*dz1w(1) - 0.5*stab0(1)*dz1w(1)
      else if (bcl(5) == 'ope') then
        ak(1:nxp,1) = 0.
        bk(1:nxp,1) = s1(1:nxp,m,1)
        ck(1:nxp,1) = 0.
        dk(1:nxp,1) = 1.
      endif
!
      if (bcl(6) == 'nnf') then
        ak(1:nxp,nzp) = dz1(nzp)*dz1w(nzp-1) - 0.5*stab0(nzp)*dz1w(nzp)
        bk(1:nxp,nzp) = s1(1:nxp,m,nzp)
        ck(1:nxp,nzp) = 0.
        dk(1:nxp,nzp) = wv(1:nxp,m) - dz1(nzp)*dz1w(nzp-1) + 0.5*stab0(nzp)*dz1w(nzp)
      else if (bcl(6) == 'ope') then
        ak(1:nxp,nzp) = 0.
        bk(1:nxp,nzp) = s1(1:nxp,m,nzp)
        ck(1:nxp,nzp) = 0.
        dk(1:nxp,nzp) = 1.
      endif
!
! solve for fourier components, x_k, given a tri-diagonal matrnxp of the
! form a_k x_k-1 + d_k x_k + c_k x_k+1 = b_k.  y_k is a scratch array.
!
      do l=1,nxp
        call TRIDAG2 (ak(l,:),dk(l,:),ck(l,:),bk(l,:),xk(l,:),nzp,ierr)
      enddo
!
    enddo
  
  end subroutine tridiag_open

!  ===========================================
  subroutine block_sendrec_real (nxp,nzp,nxp2,inputfft,outputfft,flag)
!  ===========================================
!  
  integer :: i, nxp, nzp, nxp2, flag
  integer :: req(nprocx-1), ierr, req2
  real(C_DOUBLE), dimension(1:nxp2,1,1:nzp) :: outputfft
  real(C_DOUBLE), dimension(1:nxp,1,1:nzp) :: inputfft
  integer, dimension(nprocx-1) :: BLOCKr, BLOCKr_c

#if ( defined SPMD )
  real :: stat(MPI_STATUS_SIZE)
!
!-----------------------------------------------------------!
!           	      Collect on master          	    !
!-----------------------------------------------------------!
!
  if (flag == 0) then
!
    outputfft = 0.
!
    if (mypid == 0) then
      outputfft(1:nxp,1,1:nzp) = inputfft(1:nxp,1,1:nzp)
!
      do i = 1, nprocx-1
        call MPI_IRecv ( outputfft(:,1,:), 1, BLOCKr(i), i, MPI_ANY_TAG, MPI_COMM_WORLD, req(i), ierr )
      enddo
!
      call MPI_Waitall(nprocx-1, req, MPI_STATUSES_IGNORE, ierr) 
    else
      call MPI_ISend ( inputfft(:,1,:), nxp, REALTYPE, 0, mypid, MPI_COMM_WORLD, req2, ierr )
!
      call MPI_Wait(req2, MPI_STATUS_IGNORE, ierr) 
    endif
!
  else if (flag == 1) then
!
    inputfft = 0.
    if (mypid == 0) then
      inputfft(1:nxp,1,1:nzp) = outputfft(1:nxp,1,1:nzp)
!
      do i = 1, nprocx-1
        call MPI_ISend ( outputfft(:,1,:), 1, BLOCKr(i), i, i, MPI_COMM_WORLD, req(i), ierr )
      enddo
!
      call MPI_Waitall(nprocx-1, req, MPI_STATUSES_IGNORE, ierr) 
    else
      call MPI_IRecv ( inputfft(:,1,:), nxp, REALTYPE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, req2, ierr )
!
      call MPI_Wait(req2, MPI_STATUS_IGNORE, ierr) 
    endif
!
  endif
#endif
!
  return
  end subroutine block_sendrec_real

!  ===========================================
  subroutine mpi_transpose_open (nxp,nyp,nzp,outputfft)
!  ===========================================
!  
  integer :: i, j, k, ii, nxp, nyp, nzp, ierr
  integer, dimension(0:nprocx-1) :: sizes, sizer, disps, dispr
!
  real(C_DOUBLE), dimension(1:nxp,1:nyp,1:nzp) :: outputfft
  real(C_DOUBLE), dimension(1:nxt,1:nyt,1:nzp) :: transposes, transposer
!
  integer, dimension(0:nprocx-1) :: BLOCK_TRANSPOSEa, BLOCK_TRANSPOSEb
  integer, dimension(nprocx-1) :: BLOCKr, BLOCKr_c

#if ( defined SPMD )
!
!-----------------------------------------------------------!
!           	      Transpose matrix          	    !
!-----------------------------------------------------------!
!
!  Initialize MPI types
!
  do i = 0, nprocx-1
    call MPI_TYPE_CREATE_SUBARRAY (3, (/nxt,nyt,nz/), (/nxtp,nytp,nz/), (/mypid*nxtp,i*nytp,0/),   &
  	 MPI_ORDER_FORTRAN, REALTYPE, BLOCK_TRANSPOSEa(i), ierr)
    call MPI_TYPE_COMMIT (BLOCK_TRANSPOSEa(i), ierr)
!
    call MPI_TYPE_CREATE_SUBARRAY (3, (/nxt,nyt,nz/), (/nxtp,nytp,nz/), (/mypid*nxtp,i*nytp,0/),   &
  	 MPI_ORDER_FORTRAN, REALTYPE, BLOCK_TRANSPOSEb(i), ierr)
    call MPI_TYPE_COMMIT (BLOCK_TRANSPOSEb(i), ierr)
  enddo
!
  sizes = 1
  disps = 0
  sizer = 1
  dispr = 0
!
  transposes = 0.
  transposer = 0.
  transposes(max(it_start-2,1):it_end-2,max(jt_start-2,1):jt_end-2,1:nzp) = outputfft(1:nxp,1:nyp,1:nzp)
!
!  Transpose nprocx*nprocx blocks of dimension nxt*nyt
!
  call MPI_Alltoallw (transposes, sizes, disps, BLOCK_TRANSPOSEa, 	&
  		      transposer, sizer, dispr, BLOCK_TRANSPOSEb, MPI_COMM_WORLD, ierr)
!
!  Transpose inside each block on current proc
!
  do i = 0, nprocx-1
    do k = 1, nzp
      transposer(it_start-2:it_end-2,1+i*nytp:nytp+i*nytp,k) = 		&
      	  transpose( transposer(it_start-2:it_end-2,1+i*nytp:nytp+i*nytp,k) )
    enddo
  enddo
!
  outputfft(1:nxp,1:nyp,1:nzp) = transposer(max(it_start-2,1):min(it_end-2,nx-4),max(jt_start-2,1):min(jt_end-2,ny-4),1:nzp)
!
  do i = 0, nprocx-1
    call MPI_TYPE_FREE (BLOCK_TRANSPOSEa(i), ierr)
    call MPI_TYPE_FREE (BLOCK_TRANSPOSEb(i), ierr)
  enddo
#endif
!  
  return
  end subroutine mpi_transpose_open
!
!  ===========================================    
  SUBROUTINE TRIDAG2 (A,B,C,R,U,N,CODE)
!  ===========================================

  !*****************************************************************
  ! Solves for a vector U of length N the tridiagonal linear set
  ! M U = R, where A, B and C are the three main diagonals of matrix
  ! M(N,N), the other terms are 0. R is the right side vector.
  !*****************************************************************

  INTEGER :: N, J, CODE
  REAL :: BET,GAM(N),A(N),B(N),C(N),R(N),U(N)

  IF(B(1).EQ.0.D0) THEN
    CODE=1
    RETURN
  END IF

  BET=B(1)
  U(1)=R(1)/BET
  DO J=2,N                    !Decomposition and forward substitution
    GAM(J)=C(J-1)/BET
    BET=B(J)-A(J)*GAM(J)
    IF(BET.EQ.0.D0) THEN            !Algorithm fails
      CODE=2
      RETURN
    END IF
    U(J)=(R(J)-A(J)*U(J-1))/BET
  END DO

  DO J=N-1,1,-1                     !Back substitution
    U(J)=U(J)-GAM(J+1)*U(J+1)
  END DO
  
  CODE=0
  RETURN
  END SUBROUTINE TRIDAG2
  
end module fft_solvers_open
