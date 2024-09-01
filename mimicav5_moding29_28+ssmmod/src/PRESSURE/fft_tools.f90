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
  
  module fft_tools
        
!
!-----------------------------------------------------------!
!
#ifdef SPMD
  USE mpi
#endif
  USE gridno
  use shared_data
  USE shared_pressure
  USE, intrinsic :: iso_c_binding 

  IMPLICIT NONE

  private

  public :: get_wv, tridiag, mpi_transpose, block_sendrec, shuffle, duplicate

  contains

! ===============================================
  subroutine get_wv
! ===============================================

  integer :: l,m,nl,nm,ioff,joff,nxp,nyp
  real :: xl, xm
!
!  Initialize and allocate
!
#ifdef SPMD
    nxp = (nx-5)/(nprocx*nprocy)
#else
    nxp = nx-5
#endif
!
#ifdef MODEL_3D
#ifndef CHANNEL
    nyp = ny-5
#else
    nyp = 2*(ny-5)
#endif
#else
    nyp = 1
#endif
!
    nl = nx-5
    nm = max(nyp,1)
!    
    allocate( wv(1:nxp,1:nyp) )
!
!  Defined dimensions and offsets
!
    joff = 0
    ioff = mypid*nxp
!
!  Calculate horizontal wave number
!
    do l=1,nxp
       if( l+ioff <= nl/2 ) then
         xl=real(l+ioff-1)
       else
         xl=real(l+ioff-nl-1)
       endif
!
       do m=1,nyp
#ifdef MODEL_3D
         if ( m+joff <= nm/2 ) then
           xm=real(m+joff-1)
         else
           xm=real(m+joff-nm-1)
         endif
#else
         xm=0.0
#endif
!
!  Wave number
!
#ifdef MODEL_3D
         wv(l,m) = -4.*( sin(pi*xl/real(nl))*sin(pi*xl/real(nl))/(dx*dx) + sin(pi*xm/real(nm))*sin(pi*xm/real(nm))/(dy*dy) )
#else
         wv(l,m) = -4.*sin(pi*xl/real(nl))*sin(pi*xl/real(nl))/(dx*dx)
#endif
       enddo
    enddo
!
    if (mypid == 0) wv(1,1) = 1.
!
  end subroutine get_wv

! ===============================================
  subroutine tridiag ( m, nxp, nzp, wvl, s1 )  
! ===============================================

  IMPLICIT NONE

  integer, intent(in) :: m, nxp, nzp
  real(C_DOUBLE), intent(in) :: wvl(1:nxp)
  complex(C_DOUBLE_COMPLEX), intent(inout) :: s1(1:nxp,1:nzp)

  integer :: k,l,ioff,joff,ierr
  real(C_DOUBLE) :: ak(1:nxp,1:nzp),dk(1:nxp,1:nzp),bk(1:nxp,1:nzp),ck(1:nxp,1:nzp),xk(1:nxp,1:nzp)
!
!  Defined dimensions and offsets
!
      if (bcl(5) == 'ope') s1(:,1) = cmplx(0.,0.)
      if (bcl(6) == 'ope') s1(:,nzp) = cmplx(0.,0.)
!  
! Coefficients for tri-diagonal solver: a_k x_k-1 + d_k x_k + c_k x_k+1 = b_k
!
      if (fft_solv == 2) then
        ioff = it_start - 4 
        joff = jt_start - 4
      else
        ioff = mypid*(it_end-it_start+1)
        joff = 0
      endif
!
      do k=2,nzp-1 
        do l=1,nxp
          ak(l,k) = dz1(k-1)*dz1w(k) - stab0(k)/(dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
          bk(l,k) = s1(l,k)
          ck(l,k) = dz1(k)*dz1w(k) + stab0(k)/(dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
          dk(l,k) = wvl(l) - dz1w(k)*(dz1(k) + dz1(k-1))
        enddo
      enddo
!
      if (bcl(5) == 'nnf') then
        ak(1:nxp,1) = 0.
        bk(1:nxp,1) = s1(1:nxp,1)
        ck(1:nxp,1) = dz1(1)*dz1w(1) + stab0(1)*dz1(1)
        dk(1:nxp,1) = wvl(1:nxp) - dz1(1)*dz1w(1) - stab0(1)*dz1(1)
      else if (bcl(5) == 'ope') then
        ak(1:nxp,1) = 0.
        bk(1:nxp,1) = s1(1:nxp,1)
        ck(1:nxp,1) = 0.
        dk(1:nxp,1) = 1.
      endif
!
      if (bcl(6) == 'nnf') then
        ak(1:nxp,nzp) = dz1w(nzp)*dz1(nzp-1) - stab0(nzp)*dz1(nzp-1)
        bk(1:nxp,nzp) = s1(1:nxp,nzp)
        ck(1:nxp,nzp) = 0.
        dk(1:nxp,nzp) = wvl(1:nxp) - dz1w(nzp)*dz1(nzp-1) + stab0(nzp)*dz1(nzp-1)
      else if (bcl(6) == 'ope') then
        ak(1:nxp,nzp) = 0.
        bk(1:nxp,nzp) = s1(1:nxp,nzp)
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
      do k=1,nzp
        do l=1,nxp
          if (m+joff+l+ioff > 2) bk(l,k)=aimag(s1(l,k))
          if (m+joff+l+ioff > 2) s1(l,k)=xk(l,k)
        enddo
      enddo
!
      do l=1,nxp
        call TRIDAG2 (ak(l,:),dk(l,:),ck(l,:),bk(l,:),xk(l,:),nzp,ierr)
      enddo
!
      do k=1,nzp
        do l=1,nxp
          if (m+joff+l+ioff > 2) s1(l,k) = cmplx( real(s1(l,k)),xk(l,k) )
        enddo
      enddo

  end subroutine tridiag
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
  real(C_DOUBLE) :: BET,GAM(N),A(N),B(N),C(N),R(N),U(N)

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
!
!   ======================================================	    
    subroutine shuffle ( nnx, nny, nnz, xorig, xshuff, forw_back )
!   ======================================================	    

    integer :: nnx, nny, nnz, forw_back
    integer :: i, j, ii, dimx, dimy, dimz
    integer :: iproc

    real, dimension(nnx,nny,nnz) :: xshuff
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: xorig

    integer :: ierr, send_req(nprocy-1)
    integer, dimension(0:nprocy-1) :: BLOCK_ORIG, BLOCK_SHUFF
!
! -------------------------------------------------------
!
!  Initialize MPI types and dimensions
!	
#if ( defined SPMD )
    if (verbose > 1) call write_debug('Start shuffle')
!
    iproc = int(mypid/nprocy)
!
    dimx = ip_end-ip_start+1
    dimy = jp_end-jp_start+1
    dimz = nz
!
    do i = 0, nprocy-1
      call MPI_TYPE_CREATE_SUBARRAY (3, (/dimx,dimy,dimz/), (/nnx,nny/nprocy,nnz/), (/it_start-ip_start+i*nnx,jt_start-jp_start,0/),   &
             MPI_ORDER_FORTRAN, REALTYPE, BLOCK_ORIG(i), ierr)
      call MPI_TYPE_COMMIT (BLOCK_ORIG(i), ierr)
    enddo
!
    do j = 0, nprocy-1
      call MPI_TYPE_CREATE_SUBARRAY (3, (/nnx,nny,nnz/), (/nnx,nny/nprocy,nnz/), (/0,j*nny/nprocy,0/),   &
             MPI_ORDER_FORTRAN, REALTYPE, BLOCK_SHUFF(j), ierr)
      call MPI_TYPE_COMMIT (BLOCK_SHUFF(j), ierr)
    enddo
!
!  Send data
!
    ii = 0
    if (forw_back == 0) then
      do i = 0, nprocy-1
        if (nprocy*iproc+i /= mypid) then
          ii=ii+1
          call MPI_Isend( xorig, 1, BLOCK_ORIG(i), nprocy*iproc+i, 200+i, MPI_COMM_WORLD, send_req(ii), ierr )
        else
!          xshuff(1:nnx,jt_start-3:jt_end-3,1:nnz) = xorig(it_start+i*nnx:it_start+(i+1)*nnx-1,jt_start:jt_end,1:nz)
          xshuff(1:nnx,1+i*nny/nprocy:(i+1)*nny/nprocy,1:nnz) = xorig(it_start+i*nnx:it_start-1+(i+1)*nnx,jt_start:jt_end,1:nz)
        endif
      enddo
    else
      do i = 0, nprocy-1
        if (nprocy*iproc+i /= mypid) then
          ii=ii+1
          call MPI_Isend( xshuff, 1, BLOCK_SHUFF(i), nprocy*iproc+i, 200+i, MPI_COMM_WORLD, send_req(ii), ierr )
        else
!          xorig(it_start+i*nnx:it_start+(i+1)*nnx-1,jt_start:jt_end,1:nz) = xshuff(1:nnx,jt_start-3:jt_end-3,1:nnz)
          xorig(it_start+i*nnx:it_start-1+(i+1)*nnx,jt_start:jt_end,1:nz) = xshuff(1:nnx,1+i*nny/nprocy:(i+1)*nny/nprocy,1:nnz)
        endif
      enddo    
    endif
!
!  Wait for all send to complete
!
    call MPI_Waitall ( nprocy-1, send_req, MPI_STATUSES_IGNORE, ierr )
!
!  Receive data
!
    if (forw_back == 0) then
      do j = 0, nprocy-1
        if (nprocy*iproc+j /= mypid)  &
            call MPI_Recv( xshuff, 1, BLOCK_SHUFF(j), nprocy*iproc+j, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
      enddo
    else
      do j = 0, nprocy-1
        if (nprocy*iproc+j /= mypid)  &
            call MPI_Recv( xorig, 1, BLOCK_ORIG(j), nprocy*iproc+j, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
      enddo
    endif
!
!  Free types
!
    do i = 0, nprocy-1
      call MPI_TYPE_FREE (BLOCK_ORIG(i), ierr)
    enddo
!
    do j = 0, nprocy-1
      call MPI_TYPE_FREE (BLOCK_SHUFF(j), ierr)
    enddo
!
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
    if (verbose > 1) call write_debug('Terminate shuffle')
#endif
!
  return
  end subroutine shuffle
!
!  ===========================================
  subroutine mpi_transpose (nxp,nyp,nzp,outputfft)
!  ===========================================
!  
  integer :: nxp, nyp, nzp
  complex(C_DOUBLE_COMPLEX), dimension(1:nxp,1:nyp,1:nzp) :: outputfft
!
#if ( defined SPMD )
  integer :: i, k, nxpt, nypt, ierr, is
  complex(C_DOUBLE_COMPLEX), dimension(1:nx-5,1:ny-5) :: transposer, transposes
  integer, dimension(0:nprocx-1) :: BLOCK_TRANSPOSEs, BLOCK_TRANSPOSEr, sizes, sizer, disps, dispr
!
    if (verbose > 2) call write_debug('Start mpi_transpose')
!
!-----------------------------------------------------------!
!           	      Transpose matrix          	    !
!-----------------------------------------------------------!
!
!  Initializations
  sizes = 1
  sizer = 1
  disps = 0
  dispr = 0
!
  nxpt=(nx-5)/nprocx
  nypt=(ny-5)/nprocx
  is = mypid*nxpt+1
!
  do i = 0, nprocx-1
    call MPI_TYPE_CREATE_SUBARRAY (2, (/nx-5,ny-5/), (/nxpt,nypt/), (/mypid*nxpt,i*nypt/),   &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, BLOCK_TRANSPOSEs(i), ierr)
    call MPI_TYPE_COMMIT (BLOCK_TRANSPOSEs(i), ierr)
!
    call MPI_TYPE_CREATE_SUBARRAY (2, (/nx-5,ny-5/), (/nxpt,nypt/), (/mypid*nxpt,i*nypt/),   &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, BLOCK_TRANSPOSEr(i), ierr)
    call MPI_TYPE_COMMIT (BLOCK_TRANSPOSEr(i), ierr)
  enddo
!
!  Loop over k
!
  do k = 1, nzp
    transposer = 0.
    transposes = 0.
!
!  Transpose nprocx*nprocx blocks of dimension nxt*nyt
!
    transposes(is:is+nxpt-1,1:ny-5) = outputfft(1:nxp,1:nyp,k)
!
    call MPI_Alltoallw (transposes, sizes, disps, BLOCK_TRANSPOSEs, transposer, sizer, dispr, BLOCK_TRANSPOSEr, MPI_COMM_WORLD, ierr)
!
!  Transpose inside each block on current proc
!
    do i = 0, nprocx-1
      transposer(is:is+nxpt-1,i*nypt+1:(i+1)*nypt) = transpose( transposer(is:is+nxpt-1,i*nypt+1:(i+1)*nypt) )
    enddo
!
    outputfft(1:nxp,1:nyp,k) = transposer(is:is+nxpt-1,1:ny-5)
  enddo
!
  do i = 0, nprocx-1
    call MPI_TYPE_FREE (BLOCK_TRANSPOSEs(i), ierr)
    call MPI_TYPE_FREE (BLOCK_TRANSPOSEr(i), ierr)
  enddo
!
    if (verbose > 2) call write_debug('Terminate mpi_transpose')
#endif
!  
  return
  end subroutine mpi_transpose
!
!  ===========================================
  subroutine block_sendrec (nxp,nzp,inputfft,outputfft,flag)
!  ===========================================
  
  integer :: i, nxp, nzp, flag, nxpt
  integer :: req(nprocx-1), ierr, req2
  complex(C_DOUBLE_COMPLEX), dimension(1:nx-5,1,1:nzp) :: outputfft
  complex(C_DOUBLE_COMPLEX), dimension(1:nxp,1,1:nzp) :: inputfft
  integer, dimension(1:nprocx-1) :: BLOCKr_c
!
!-----------------------------------------------------------!
!           	      Collect on master          	    !
!-----------------------------------------------------------!
!
#if ( defined SPMD )
!
  nxpt=(nx-5)/nprocx
  do i = 1, nprocx-1
    call MPI_TYPE_CREATE_SUBARRAY (2, (/nx-5,nzp/), (/nxp,nzp/), (/i*nxpt,0/),   &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, BLOCKr_c(i), ierr)
    call MPI_TYPE_COMMIT ( BLOCKr_c(i), ierr )
  enddo
!
  if (flag == 0) then
!
    outputfft = 0.
!
    if (mypid == 0) then
      outputfft(1:nxp,1,1:nzp) = inputfft(1:nxp,1,1:nzp)
!
      do i = 1, nprocx-1
        call MPI_IRecv ( outputfft(:,1,:), 1, BLOCKr_c(i), i, MPI_ANY_TAG, MPI_COMM_WORLD, req(i), ierr )
      enddo
!
      call MPI_Waitall(nprocx-1, req, MPI_STATUSES_IGNORE, ierr) 
    else
      call MPI_ISend ( inputfft(:,1,:), nxp*nzp, MPI_DOUBLE_COMPLEX, 0, mypid, MPI_COMM_WORLD, req2, ierr )
!
      call MPI_Wait(req2, MPI_STATUS_IGNORE, ierr) 
    endif
!
  else if (flag == 1) then
!
    inputfft = 0.
!
    if (mypid == 0) then
      inputfft(1:nxp,1,1:nzp) = outputfft(1:nxp,1,1:nzp)
!
      do i = 1, nprocx-1
        call MPI_ISend ( outputfft(:,1,:), 1, BLOCKr_c(i), i, i, MPI_COMM_WORLD, req(i), ierr )
      enddo
!
      call MPI_Waitall(nprocx-1, req, MPI_STATUSES_IGNORE, ierr) 
    else
      call MPI_IRecv ( inputfft(:,1,:), nxp*nzp, MPI_DOUBLE_COMPLEX, 0, MPI_ANY_TAG, MPI_COMM_WORLD, req2, ierr )
!
      call MPI_Wait(req2, MPI_STATUS_IGNORE, ierr) 
    endif
!
  endif
!
  do i = 1, nprocx-1
    call MPI_TYPE_FREE ( BLOCKr_c(i), ierr )
  enddo
!
#endif
!
  return
  end subroutine block_sendrec
!
!   ======================================================	    
    subroutine duplicate ( nnx, nny, nnz, xorig, xcopy, forw_back )
!   ======================================================	    

    integer :: nnx, nny, nnz, nv, forw_back
    integer :: i, j
!
    real, dimension(nnx,nny,nnz) :: xcopy
    real, dimension(nnx,ny-5,nnz) :: xshuff
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: xorig
!
! -------------------------------------------------------
!
  if (verbose > 2) call write_debug('Start duplicate')
!
  nv = ny - 5
!
!  Rearranging data
!	
#ifdef SPMD
!
#ifdef DECOMP_2D
  call shuffle ( nnx, nv, nnz, xorig, xshuff, 0 )
#else
  do j = jt_start, jt_end
    do i = it_start, it_end
      xshuff(i-it_start+1,j-jt_start+1,:) = xorig(i,j,:)
    enddo
  enddo
#endif
!
#else
!
  xshuff = xorig(it_start:it_end,jt_start:jt_end,1:nz)
!
#endif
!
!  Duplication with symmetry
!
  xcopy(:,1:nv,:) = xshuff(:,1:nv,:)
!
  do j = nv, 1, -1
    xcopy(:,2*nv-j+1,:) = xshuff(:,j,:)
  enddo
!
  if (verbose > 2) call write_debug('Terminate duplicate')
!
  return
  end subroutine duplicate
  
end module fft_tools
