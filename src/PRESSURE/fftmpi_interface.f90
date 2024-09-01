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
  
module fftmpi_solver
        
!
!-----------------------------------------------------------!
!

  USE gridno
  USE shared_data
  USE allocation
  USE fft_tools
!
  use iso_c_binding
  use fft2d_wrap
!
#ifdef SPMD
  USE mpi
  USE mpicomm
#endif  

  IMPLICIT NONE

  private

  public :: fftmpi_solve

  contains

! ===============================================
  subroutine fftmpi_solve ( tend, output )
! ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: tend
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: output
  real, dimension(:), allocatable :: work
!
  integer :: i, j, k, l 
  integer :: nxp, nyp, nzp
  integer :: ierr, fftsize, prec=2
!
#if (defined SPMD) && (defined MODEL_3D) 
!
  TYPE(C_PTR) :: fft_plan
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: tendfft
  real(C_DOUBLE), dimension(:,:), allocatable :: wvl
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Start fftmpi_solve')
!
!-----------------------------------------------------------!
!     		        Allocations               	    !
!-----------------------------------------------------------!
!
  nxp = int((nx-5)/nprocx)
  nyp = int((ny-5)/nprocy)
  nzp = nz
!
  call fft2d_create(MPI_COMM_WORLD,prec,fft_plan)
!
  call plan ( nxp, nyp, fft_plan, fftsize ) 
!
  call alloc(work,2*fftsize) 
  allocate ( tendfft(1:nxp,1:nyp,1:nzp), wvl(1:nxp,1:nyp) )
!
  call calculate_wv ( wvl, nxp, nyp )
!
!-----------------------------------------------------------!
!     		        Forward FFT              	    !
!-----------------------------------------------------------!
!
  do k = 1, nzp
    l = 0
    do j = jt_start, jt_end
      do i = it_start, it_end
        work(2*l+1) = tend(i,j,k)
        work(2*l+2) = 0.
        l = l + 1
      enddo
    enddo
!
    call fft2d_compute(fft_plan,c_loc(work),c_loc(work),1)
!
    l = 0
    do j = 1, nyp 
      do i = 1, nxp 
        tendfft(i,j,k) = cmplx(work(2*l+1),work(2*l+2))
        l = l + 1
      enddo
    enddo
  enddo
!
!-----------------------------------------------------------!
!     		        Solve tridiag               	    !
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Start tridiag')
  do j = 1, nyp
    call tridiag ( j, nxp, nzp, wvl(:,j), tendfft(:,j,:) )  
  enddo
  if (verbose > 1) call write_debug('Terminate tridiag')
!
!-----------------------------------------------------------!
!     		        Backward FFT               	    !
!-----------------------------------------------------------!
!
  do k = 1, nzp
    l = 0
    do j = 1, nyp 
      do i = 1, nxp 
        work(2*l+1) = real(tendfft(i,j,k))
        work(2*l+2) = aimag(tendfft(i,j,k))
        l = l + 1
      enddo
    enddo
!
    call fft2d_compute(fft_plan,c_loc(work),c_loc(work),-1)
!
    l = 0
    do j = jt_start, jt_end
      do i = it_start, it_end
        output(i,j,k) = (nx-5)*(ny-5)*work(2*l+1)
        l = l + 1
      enddo
    enddo
  enddo
!
!-----------------------------------------------------------!
!     		        Deallocations               	    !
!-----------------------------------------------------------!
!
  deallocate( tendfft, wvl )
  call dealloc(work)
!
  call fft2d_destroy(fft_plan) 
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate fftmpi_solve')
!
#endif
!
  end subroutine fftmpi_solve

!
!-----------------------------------------------------------!
!

! ===============================================
  subroutine plan (nxp, nyp, fft_plan, fftsize)
! ===============================================

  implicit none

  type(c_ptr) :: fft_plan
  integer :: sendsize,recvsize,fftsize,flag,ierr
  integer :: inxlo, inxhi, inylo, inyhi
  integer :: nxp, nyp,  nxa, nya, ipx, ipy
  integer :: eflag, cflag = 2, pflag = 0, sflag = 1
!
#ifdef DECOMP_2D
  ! Flag for brick decomposition
  eflag = 1
#else
  ! Flag for pencil decomposition
  eflag = 0
#endif
!
  inxlo = it_start - 3
  inylo = jt_start - 3
  inxhi = inxlo + nxp - 1
  inyhi = inylo + nyp - 1
!  
  call fft2d_set(fft_plan,"collective",cflag)
  call fft2d_set(fft_plan,"exchange",eflag)
  call fft2d_set(fft_plan,"pack",pflag)
  call fft2d_set(fft_plan,"scaled",sflag)
!
  call fft2d_setup(fft_plan, nx-5, ny-5, inxlo, inxhi, inylo, inyhi, inxlo, inxhi, inylo, inyhi,      &
          0, fftsize, sendsize, recvsize)
!
end subroutine plan

!
!-----------------------------------------------------------!
!

! ===============================================
  subroutine calculate_wv(wvl,nxp,nyp)
! ===============================================

  implicit none

  integer :: nxp, nyp
  real, dimension(nxp,nyp) :: wvl
!
  integer :: l,m,nl,nm,ioff,joff,ipx,ipy
  real :: xl, xm
!
!  Defined dimensions and offsets
!
    nl = nx-5
    nm = ny-5 
!
    ioff = it_start - 4 
    joff = jt_start - 4
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
         if ( m+joff <= nm/2 ) then
           xm=real(m+joff-1)
         else
           xm=real(m+joff-nm-1)
         endif
!
!  Wave number
!
         wvl(l,m) = -4.*( sin(pi*xl/real(nl))*sin(pi*xl/real(nl))/(dx*dx) + sin(pi*xm/real(nm))*sin(pi*xm/real(nm))/(dy*dy) )
       enddo
    enddo
!
    if (mypid == 0) wvl(1,1) = 1.
!
  end subroutine calculate_wv

end module fftmpi_solver
