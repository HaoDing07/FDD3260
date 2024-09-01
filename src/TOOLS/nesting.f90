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
!  NESTING.F90                   
!
!  Purpose:
!	Sets boundary conditions for nested runs from external 
!	nesting.out file			  
!
!  Author
!	Julien Savre, LMU Munich
!
! ===========================================================

module nesting

#ifdef SPMD
  USE mpi
  USE mpicomm
#endif
  USE gridno
  USE shared_data
  USE shared_state
  USE shared_hydro
  USE shared_wind
  USE shared_thermo
  USE shared_pressure
!
  USE thermodynamics
  
  implicit none
  
  logical :: lfirst=.true.
  integer :: nx_nest, ny_nest, nz_nest
   
  real :: ftime, t_nest1, t_nest2, p00(nz)
  real, allocatable, dimension(:) :: xnest, ynest, znest
  real, allocatable, dimension(:,:,:) :: unest1, vnest1, wnest1
  real, allocatable, dimension(:,:,:) :: pnest1, ptnest1, qtnest1
  real, allocatable, dimension(:,:,:,:) :: qnest1, nnest1, snest1
  real, allocatable, dimension(:,:,:) :: unest2, vnest2, wnest2
  real, allocatable, dimension(:,:,:) :: pnest2, ptnest2, qtnest2
  real, allocatable, dimension(:,:,:,:) :: qnest2, nnest2, snest2
        
  private
    
  public :: savenest, readnest, load_nest, dealloc_nest
  public :: ftime, unest1, unest2, vnest1, vnest2, wnest1, wnest2
  public :: ptnest1, ptnest2, qtnest1, qtnest2, qnest1, nnest1, qnest2, nnest2
  
  save
  
  contains
  
!  ================================================
  subroutine savenest
!  ================================================
    
    implicit none
   
    integer :: i, j, h, ic, jc, i1, i2, j1, j2, nxn, nyn
    character(len = 100):: filename
    
    real, allocatable, dimension(:) :: xx, yy
    real, allocatable, dimension(:,:,:) :: q1      
!
    if (verbose > 2) call write_debug('Start savenest')
! 
!----------------------------------------------------------!
!
    allocate(q1(1:nx,1:ny,1:nz), xx(1:nx), yy(1:ny) )
!
!  Open nesting file and write dimensions
!
    filename = opdir(1:len_trim(opdir)) // '/nesting.dat'
    if (lfirst .and. mypid == 0) then
      open(108,file=filename(1:len_trim(filename)),form='unformatted',status='unknown',access='sequential')
    else if (mypid == 0) then
      open(108,file=filename(1:len_trim(filename)),form='unformatted',status='old',access='sequential',position='append')      
    endif
!
!  Find position of the nest
!
    i1 = 1
    i2 = 1
    j1 = 1
    j2 = 1
!
    nxn = ceiling(lx_nest/dx)+1
    nyn = ceiling(ly_nest/dy)+1
    if (mod(nxn,2) /= 0) nxn=nxn+1
    if (mod(nyn,2) /= 0) nyn=nyn+1
!
    do i = 1, nx
      xx(i) = real(i-1)*dx
    enddo
    do i = 2, nx
      if (xc_nest >= xx(i-1) .and. xc_nest < xx(i)) ic=i-1
    enddo
    i1 = ic-nxn/2
    i2 = ic+nxn/2+1
!
    do j = 1, ny
      yy(j) = real(j-1)*dy
    enddo
#ifdef MODEL_3D
    do j = 2, ny
      if (yc_nest >= yy(j-1) .and. yc_nest < yy(j)) jc=j-1
    enddo
    j1 = jc-nyn/2
    j2 = jc+nyn/2+1
#endif
!
!  Writing full data in nesting file
!
    if (lfirst .and. mypid == 0) then
      nx_nest = i2-i1+1
      ny_nest = j2-j1+1
      nz_nest = nz
      write(108)nx_nest,ny_nest,nz_nest
      write(108)xx(i1:i2)
      write(108)yy(j1:j2)
    endif
!
    if (mypid == 0) write(108)
    if (mypid == 0) write(108)time
!
#ifdef SPMD
    call collect ( wind2%u,  q1 )
    if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)
!
    call collect ( wind2%v,  q1 )
    if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)
!
    call collect ( wind2%w,  q1 )
    if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)
!
    call collect ( pressure2%p,  q1 )
    if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)
!
    call collect ( state2%es,  q1 )
    if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)
!
    call collect ( state2%qt,  q1 )
    if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)
!
    do h = 1, nhydro
      call collect ( hydromtr2(h)%q,  q1 )
      if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)
!
      call collect ( hydromtr2(h)%n,  q1 )
      if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)      
    enddo
!
    do h = 1, nscal
      call collect ( state2%scal(:,:,:,h),  q1 )
      if (mypid == 0) write(108)q1(i1:i2,j1:j2,1:nz)
    enddo
#else
    write(108)wind2%u(i1:i2,j1:j2,1:nz)
!
    write(108)wind2%v(i1:i2,j1:j2,1:nz)
!
    write(108)wind2%w(i1:i2,j1:j2,1:nz)
!
    write(108)pressure2%p(i1:i2,j1:j2,1:nz)
!
    write(108)state2%es(i1:i2,j1:j2,1:nz)
!
    write(108)state2%qt(i1:i2,j1:j2,1:nz)
!
    do h = 1, nhydro
      write(108)hydromtr2(h)%q(i1:i2,j1:j2,1:nz)
!
      write(108)hydromtr2(h)%n(i1:i2,j1:j2,1:nz)
    enddo
!
    do h = 1, nscal
      write(108)state2%scal(i1:i2,j1:j2,1:nz,h)
    enddo
#endif
!
!  Closing
!
    lfirst=.false.
!
    close(108)
!
    deallocate(q1, xx, yy)
! 
!----------------------------------------------------------!
!
    if (verbose > 2) call write_debug('Terminate savenest')
!
  end subroutine savenest
  
!  ================================================
  subroutine readnest
!  ================================================
    
    implicit none
    
    integer :: h, ierr
    character(len = 100):: filename
!
    if (verbose > 2) call write_debug('Start readnest')
! 
!----------------------------------------------------------!
!
!  Open nesting file, allocate and read dimensions
!
    filename = opdir(1:len_trim(opdir)) // '/nesting.dat'
!
    if (mypid == 0 .and. lfirst) then
      open(108,file=filename(1:len_trim(filename)),form='unformatted',status='old')
!
      read(108)nx_nest,ny_nest,nz_nest
!
      allocate(xnest(1:nx_nest), ynest(1:ny_nest), znest(1:nz_nest))
!
      read(108)xnest(1:nx_nest)
      read(108)ynest(1:ny_nest)
      read(108)znest(1:nz_nest)
!
      allocate(unest1(1:nx_nest,1:ny_nest,1:nz_nest), vnest1(1:nx_nest,1:ny_nest,1:nz_nest), wnest1(1:nx_nest,1:ny_nest,1:nz_nest))
      allocate(pnest1(1:nx_nest,1:ny_nest,1:nz_nest), ptnest1(1:nx_nest,1:ny_nest,1:nz_nest), qtnest1(1:nx_nest,1:ny_nest,1:nz_nest))
      allocate(qnest1(1:nx_nest,1:ny_nest,1:nz_nest,1:nhydro), nnest1(1:nx_nest,1:ny_nest,1:nz_nest,1:nhydro))
      allocate(snest1(1:nx_nest,1:ny_nest,1:nz_nest,1:nscal))
!
      allocate(unest2(1:nx_nest,1:ny_nest,1:nz_nest), vnest2(1:nx_nest,1:ny_nest,1:nz_nest), wnest2(1:nx_nest,1:ny_nest,1:nz_nest))
      allocate(pnest2(1:nx_nest,1:ny_nest,1:nz_nest), ptnest2(1:nx_nest,1:ny_nest,1:nz_nest), qtnest2(1:nx_nest,1:ny_nest,1:nz_nest))
      allocate(qnest2(1:nx_nest,1:ny_nest,1:nz_nest,1:nhydro), nnest2(1:nx_nest,1:ny_nest,1:nz_nest,1:nhydro))
      allocate(snest2(1:nx_nest,1:ny_nest,1:nz_nest,1:nscal))
!
!  Read main data: first tstep
!
      read(108)
      read(108)t_nest1
! 
      read(108)unest1(1:nx_nest,1:ny_nest,1:nz_nest)
! 
      read(108)vnest1(1:nx_nest,1:ny_nest,1:nz_nest)
!
      read(108)wnest1(1:nx_nest,1:ny_nest,1:nz_nest)
!
      read(108)pnest1(1:nx_nest,1:ny_nest,1:nz_nest)
!
      read(108)ptnest1(1:nx_nest,1:ny_nest,1:nz_nest)
!
      read(108)qtnest1(1:nx_nest,1:ny_nest,1:nz_nest)
!
      do h = 1, nhydro
        read(108)qnest1(1:nx_nest,1:ny_nest,1:nz_nest,h)
!
        read(108)nnest1(1:nx_nest,1:ny_nest,1:nz_nest,h)
      enddo
!
      do h = 1, nscal
        read(108)snest1(1:nx_nest,1:ny_nest,1:nz_nest,h)
      enddo
!
      p00 = p0
      t_nest2 = t_nest1
      unest2 = unest1
      vnest2 = vnest1
      wnest2 = wnest1
      pnest2 = pnest1
      ptnest2 = ptnest1
      qtnest2 = qtnest1
      qnest2 = qnest1
      nnest2 = nnest1
      snest2 = snest1
    endif
!
!  Communicate x, y, z dimensions
!
#if ( defined SPMD )
    if (lfirst) then
      call MPI_BCAST(nx_nest,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ny_nest,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nz_nest,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
!
      if (mypid /= 0) allocate(xnest(1:nx_nest), ynest(1:ny_nest), znest(1:nz_nest))
!
      call MPI_BCAST(xnest, nx_nest,  REALTYPE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ynest, ny_nest,  REALTYPE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(znest, nz_nest,  REALTYPE,0,MPI_COMM_WORLD,ierr)
    endif
#endif
!
!  Read main data: subsequent tsteps
!
    if (mypid == 0 .and. time >= t_nest2) then
      t_nest1 = t_nest2
      unest1 = unest2
      vnest1 = vnest2
      wnest1 = wnest2
      pnest1 = pnest2
      ptnest1 = ptnest2
      qtnest1 = qtnest2
      qnest1 = qnest2
      nnest1 = nnest2
      snest1 = snest2
!
      read(108)
      read(108)t_nest2
! 
      read(108)unest2(1:nx_nest,1:ny_nest,1:nz_nest)
! 
      read(108)vnest2(1:nx_nest,1:ny_nest,1:nz_nest)
!
      read(108)wnest2(1:nx_nest,1:ny_nest,1:nz_nest)
!
      read(108)pnest2(1:nx_nest,1:ny_nest,1:nz_nest)
!
      read(108)ptnest2(1:nx_nest,1:ny_nest,1:nz_nest)
!
      read(108)qtnest2(1:nx_nest,1:ny_nest,1:nz_nest)
!
      do h = 1, nhydro
        read(108)qnest2(1:nx_nest,1:ny_nest,1:nz_nest,h)
!
        read(108)nnest2(1:nx_nest,1:ny_nest,1:nz_nest,h)
      enddo
!
      do h = 1, nscal
        read(108)snest2(1:nx_nest,1:ny_nest,1:nz_nest,h)
      enddo
    endif
!
!  Set boundary values in ghost cells (these won't be updated by bc routines later on)
!
    if (mypid == 0) ftime = (time - t_nest1) / (t_nest2 - t_nest1)
!
    call load_nest (1, 2, nx-1, nx, 1, 2, ny-1, ny, ftime, unest1, unest2, wind2%u)
!
    call load_nest (1, 2, nx-1, nx, 1, 2, ny-1, ny, ftime, vnest1, vnest2, wind2%v)
!
    call load_nest (1, 2, nx-1, nx, 1, 2, ny-1, ny, ftime, wnest1, wnest2, wind2%w)
!
    call load_nest (1, 2, nx-1, nx, 1, 2, ny-1, ny, ftime, ptnest1, ptnest2, state2%es)
!
    call load_nest (1, 2, nx-1, nx, 1, 2, ny-1, ny, ftime, qtnest1, qtnest2, state2%qt)
!
    do h = 1, nhydro
      call load_nest (1, 2, nx-1, nx, 1, 2, ny-1, ny, ftime, qnest1(:,:,:,h), qnest2(:,:,:,h), hydromtr2(h)%q)
!
      call load_nest (1, 2, nx-1, nx, 1, 2, ny-1, ny, ftime, nnest1(:,:,:,h), nnest2(:,:,:,h), hydromtr2(h)%n)
    enddo
!
    do h = 1, nscal
      call load_nest (1, 2, nx-1, nx, 1, 2, ny-1, ny, ftime, snest1(:,:,:,h), snest2(:,:,:,h), state2%scal(:,:,:,h))
    enddo
!
    lfirst=.false.
!
!----------------------------------------------------------!
!
    if (verbose > 2) call write_debug('Terminate readnest')
!
  end subroutine readnest
!
!  ================================================
  subroutine load_nest ( i1, i2, i3, i4, j1, j2, j3, j4, ft, qnest1, qnest2, data )
!  ================================================

    implicit none
!
    integer :: i1, i2, i3, i4, j1, j2, j3, j4, ierr
    real, dimension(1:nx_nest,1:ny_nest,1:nz_nest) :: q1
!
    real, intent(in) :: ft
    real, dimension(ip_start:,jp_start:,:), intent(inout) :: data
    real, dimension(1:nx_nest,1:ny_nest,1:nz_nest), intent(in) :: qnest1, qnest2
!
!--------------------------------------------------------!
!
!  Time interpolation
!
    q1 = 0.
    if (mypid == 0) q1 = qnest1 + ft*(qnest2 - qnest1)
!
!  Communications
!
#ifdef SPMD
    call MPI_BCAST(q1, nx_nest*ny_nest*nz_nest, REALTYPE, 0, MPI_COMM_WORLD, ierr)
#endif
!
!  Interpolation
!
    call interp_nest ( i1, i2, i3, i4, j1, j2, j3, j4, q1, data )
  
  end subroutine load_nest

!  ================================================
  subroutine interp_nest ( i1, i2, i3, i4, j1, j2, j3, j4, nest, data )
!  ================================================

    implicit none
!
    integer :: i1, i2, i3, i4, j1, j2, j3, j4, i, j, k
    real :: xx, yy
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: data
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz_nest) :: data_bis
    real, dimension(1:nx_nest,1:ny_nest,1:nz_nest) :: nest
!
!--------------------------------------------------------!
!
!  Left boundary points in X
!
    if (ip_start == 1) then   
      do i = i1, i2
	xx = xc_nest - dx*real(nx/2-i+1)
        do j = jp_start, jp_end
	  yy = yc_nest - dy*real(ny/2-j+1)
!
          do k = 1, nz_nest
            call pwl_interp_2d ( nx_nest, ny_nest, xnest, ynest, nest(:,:,k), 1, xx, yy, data_bis(i,j,k) )
          enddo
!
          do k = 1, nz
            call pwl_interp_1d ( nz_nest, znest, data_bis(i,j,:), 1, z0(k), data(i,j,k) )
          enddo
        enddo
      enddo
    endif
!
!  Right boundary points in X
!
    if (ip_end == nx) then   
      do i = i3, i4
	xx = xc_nest - dx*real(nx/2-i+1)
        do j = jp_start, jp_end
	  yy = yc_nest - dy*real(ny/2-j+1)
!
          do k = 1, nz_nest
            call pwl_interp_2d ( nx_nest, ny_nest, xnest, ynest, nest(:,:,k), 1, xx, yy, data_bis(i,j,k) )
          enddo
!
          do k = 1, nz
            call pwl_interp_1d ( nz_nest, znest, data_bis(i,j,:), 1, z0(k), data(i,j,k) )
          enddo
        enddo
      enddo
    endif
!
!  Left boundary points in Y
!
#ifdef MODEL_3D
    if (jp_start == 1) then   
      do j = j1, j2
	yy = yc_nest - dy*real(ny/2-j+1)
        do i = ip_start, ip_end
	  xx = xc_nest - dx*real(nx/2-i+1)
!
          do k = 1, nz_nest
            call pwl_interp_2d ( nx_nest, ny_nest, xnest, ynest, nest(:,:,k), 1, xx, yy, data_bis(i,j,k) )
          enddo
!
          do k = 1, nz
            call pwl_interp_1d ( nz_nest, znest, data_bis(i,j,:), 1, z0(k), data(i,j,k) )
          enddo
        enddo
      enddo
    endif
!
!  Right boundary points in Y
!
    if (jp_end == nx) then   
      do j = j3, j4
	yy = yc_nest - dy*real(ny/2-j+1)
        do i = ip_start, ip_end
	  xx = xc_nest - dx*real(nx/2-i+1)
!
          do k = 1, nz_nest
            call pwl_interp_2d ( nx_nest, ny_nest, xnest, ynest, nest(:,:,k), 1, xx, yy, data_bis(i,j,k) )
          enddo
!
          do k = 1, nz
            call pwl_interp_1d ( nz_nest, znest, data_bis(i,j,:), 1, z0(k), data(i,j,k) )
          enddo
        enddo
      enddo
    endif
#endif
!  
  end subroutine interp_nest

  subroutine dealloc_nest
  
    USE gridno
    USE shared_data
  
    IMPLICIT NONE
!
    if (mypid == 0) then
      deallocate(xnest, ynest, znest)
!
      deallocate(unest1, vnest1, wnest1)
      deallocate(pnest1, ptnest1, qtnest1)
      deallocate(qnest1, nnest1)
      deallocate(snest1)
!
      deallocate(unest2, vnest2, wnest2)
      deallocate(pnest2, ptnest2, qtnest2)
      deallocate(qnest2, nnest2)
      deallocate(snest2)
    endif
!
  end subroutine dealloc_nest

end module nesting
