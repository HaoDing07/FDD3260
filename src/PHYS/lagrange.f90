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
!  LAGRANGE.F:                   
!
!  Purpose:
!       Lagrangian particle tracker for MIMICA                      
!	Note that so far, all parcels are treated on the master proc
!
!  Author
!       Julien Savre
!       MISU, Stockholm
!
! ================================================================
!
module lagrange
!
#ifdef SPMD
  USE mpi
#endif
  use gridno
  use shared_data
  use shared_lagrange
  use shared_all
  use averages
!
  IMPLICIT NONE
!
  private
!
#ifdef LAGRANGE
!
  public :: lagrange_step, lagrange_init, lagrange_interp, lagrange_proc
!
  CONTAINS
!
!  ===================
  subroutine lagrange_step
!  ===================
!
  if (verbose > 0) call write_debug('Starting lagrange_step')
!
!----------------------------------------------------------!
!      	     Evaluate turbulent fluctuations               !
!----------------------------------------------------------!
!
  call lagrange_turbu (parcel_pos, parcel_sca)
!
!----------------------------------------------------------!
!                  Aerosol particles                       !
!----------------------------------------------------------!
!
  if (aerosol_lag) call lagrange_aerosol
!
!----------------------------------------------------------!
!    	     Advance: find new parcel position	           !
!----------------------------------------------------------!
!
  call lagrange_adv
!
!----------------------------------------------------------!
!    	             Find new parcels proc   	           !
!----------------------------------------------------------!
!
  call lagrange_proc (.false.)
!
!----------------------------------------------------------!
!     	       Interpolate new parcel states               !
!----------------------------------------------------------!
!
  call lagrange_interp (.false., .false.)
!
!----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating lagrange_step')
!
end
!
!
!  ===================
  subroutine lagrange_init
!  ===================
!
  INTEGER :: i, ierr
  logical :: expart=.false.
!
  if (verbose > 0) call write_debug('Starting lagrange_init')
!
!----------------------------------------------------------!
!    		Initialize parcels position   	           !
!----------------------------------------------------------!
!
!  If parcel file exists:
!
  INQUIRE (FILE='./parcels.dat', EXIST=expart)
!
  if (expart) then
!
    if (mypid == 0) then
      open(900, FILE='./parcels.dat')
      read(900,*) nparc
    endif
!
#ifdef SPMD
    call MPI_BCAST(nparc, 1, INTTYPE,0,MPI_COMM_WORLD,ierr)
#endif
!
    allocate (parcel_loc(1:nparc), aer(1:nparc))
    call alloc_parc( parcel_pos )
    call alloc_parc( parcel_sca )
    parcel_loc = 0
!
    if (mypid == 0) then
      do i = 1, nparc
    	read(900,*) parcel_pos(i)%x, parcel_pos(i)%y, parcel_pos(i)%z
      enddo
      if (aerosol_lag) call reset_aerosols( aer )
      close (900)
    endif
!
!  ELSE: Generate initial parcel positions randomly
!
  else
    nparc = npart
!
    allocate (parcel_loc(1:nparc), aer(1:nparc))
    call alloc_parc( parcel_pos )
    call alloc_parc( parcel_sca )
    parcel_loc = 0
!      
    call reset_particles ( 0 )
!
    if (aerosol_lag) call reset_aerosols( aer )
  endif
!
  if (mypid==0) write(7,*) 'INITIALIZATION OF', nparc, 'PARCELS:'
!
!  Find parcel's proc location
!
  call lagrange_proc ( .true. )
!
!----------------------------------------------------------!
!    	         First scalar interpolations  	           !
!----------------------------------------------------------!
!
  call lagrange_interp ( .false., .true. )
!
!----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating lagrange_init')
!
  return
  end subroutine lagrange_init
!
!  ===================
  subroutine lagrange_adv
!  ===================
!
#if ( defined SPMD )
  integer :: ierr, status(MPI_STATUS_SIZE)
#endif
  integer :: i, n
  real    :: xx, yy, zz, lx, ly
!
  if (verbose > 1) call write_debug('Starting lagrange_adv')
!
!----------------------------------------------------------!
!    	       Advance parcels using wind2   	           !
!----------------------------------------------------------!
!
  lx = real(nx-5)*dx
#ifdef MODEL_3D
  ly = real(ny-5)*dy
#else
  ly = 0.
#endif
!
  do n = 1, nparc_loc
    i = parcel_loc(n)
!
    if ( parcel_pos(i)%exist ) then
!
!  Advance in time
!
      xx = parcel_sca(i)%u + parcel_sca(i)%uprim
      yy = parcel_sca(i)%v + parcel_sca(i)%vprim
      if (.not. aerosol_lag) then
        zz = parcel_sca(i)%w + parcel_sca(i)%wprim
      else
        zz = parcel_sca(i)%w + parcel_sca(i)%wprim + aer(i)%ws
      endif
!
      parcel_pos(i)%x = parcel_pos(i)%x + dt0*xx
      parcel_pos(i)%y = parcel_pos(i)%y + dt0*yy
      parcel_pos(i)%z = parcel_pos(i)%z + dt0*zz
!
!  Periodicity & Discard parcels leaving the domain
!
      if (bcl(1)=='per' .or. bcl(2)=='per') then
        if (parcel_pos(i)%x < 0.) parcel_pos(i)%x = lx + parcel_pos(i)%x
        if (parcel_pos(i)%x > lx) parcel_pos(i)%x = parcel_pos(i)%x - lx
      else if (bcl(1)/='per' .and. bcl(2)/='per') then
        if (parcel_pos(i)%x < 0. .or. parcel_pos(i)%x > lx) parcel_pos(i)%exist=.false.
      endif
!
#ifdef MODEL_3D
      if (bcl(3)=='per' .or. bcl(4)=='per') then
        if (parcel_pos(i)%y < 0.) parcel_pos(i)%y = ly + parcel_pos(i)%y
        if (parcel_pos(i)%y > ly) parcel_pos(i)%y = parcel_pos(i)%y - ly
      else if (bcl(3)/='per' .and. bcl(4)/='per') then
        if (parcel_pos(i)%y < 0. .or. parcel_pos(i)%y > ly) parcel_pos(i)%exist=.false.
      endif
#endif
!    
      if ( parcel_pos(i)%z <= 0. ) parcel_pos(i)%z = z0(1)
      if ( parcel_pos(i)%z >= zmax ) parcel_pos(i)%z = z0(nz)
!
!  Clean and reset parcels
!
      if ( .not.parcel_pos(i)%exist .and. res_lag ) then
        call reset_particles ( i )
        parcel_pos(i)%exist=.true.
      endif
!
    endif
!
  enddo
!
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating lagrange_adv')
!
  return
  end
!
!  ===================
  subroutine lagrange_turbu (parcel_pos, parcel_sca)
!  ===================
!
  integer :: i, j, n, SEED, T(8)
  real :: r(3)
  type (position), dimension(nparc) :: parcel_pos
  type (scalar), dimension(nparc) :: parcel_sca
!
!----------------------------------------------------------!
!    	    Find turbulent velocity fluctuations           !
!----------------------------------------------------------!
!
!  Initialize seed for RNG 
!
  CALL DATE_AND_TIME(VALUES = T)
  SEED = 56*mypid + ntime*T(1)+70*(T(2)/ntime+12*i*(T(3)+31*(T(5)+23*(T(6)/ntime+59*T(7)))))
  IF (MOD(SEED,2).EQ.0) SEED = SEED-1
!
  do n = 1, nparc_loc
    i = parcel_loc(n)
!
!  Calculate pseudo-random perturbations
!
    call r8vec_normal_01 ( 3, SEED, r ) 
!
    do j = 1, 3
      r(j) = sign(1.,r(j)) * min(max(abs(r(j))*sqrt(2./3.*max(parcel_sca(i)%tke,0.)),1.e-10),2.)
    enddo
!
    parcel_sca(i)%uprim = r(1)
#ifdef MODEL_3D
    parcel_sca(i)%vprim = r(2)
#else
    parcel_sca(i)%vprim = 0.
#endif
    parcel_sca(i)%wprim = r(3)
  enddo
!
  return
  end subroutine
!
!  ===================
  subroutine lagrange_interp ( output, init )
!  ===================
!
  logical :: output, init, lmask=.true.
  logical, save :: lfirst=.true.
  integer :: n, i, j, k, h, p, l, step
  integer :: nchange(6), nmask(3), nchangep(6), nmaskp(3)
  integer, save :: nchange_tot(6), nmask_tot(3)
  integer :: dim, nd(2), iproc, ierr
  integer, allocatable, dimension(:,:) :: ind
!
  real :: r_filter, ff
  real, allocatable, dimension(:) :: a, b
!
  real, dimension(1:3,1:3,1:3) :: filter
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: wav
  real, dimension(it_start-1:ip_end+1) :: xx
  real, dimension(jt_start-1:jp_end+1) :: yy
  real, dimension(1:nz) :: zz
!
  if (verbose > 2) call write_debug('Starting lagrange_interp')
!
!----------------------------------------------------------!
!    		    Initialize dimensions                  !
!----------------------------------------------------------!
!
!  Initialize auxiliary variables 
!
  wav = 0.
!
  do k = 1, nz-1
    wav(:,:,k) = 0.5*(wind%w(:,:,k+1) + wind%w(:,:,k))
  enddo
!
!  Dimensions and allocations
!
#ifdef MODEL_3D
  dim = 3
#else
  dim = 2
#endif
!
  l = lag_ord
!
  allocate (ind(1:2,1:dim),a(1:dim),b(1:dim))
!
!  Set global dimensions
!
  zz = 0.5*dz*fdz0(1)
  do i = it_start-1, it_end+1
    xx(i) = real(i-4)*dx
  enddo   
#ifdef MODEL_3D
  do j = jt_start-1, jt_end+1
    yy(j) = real(j-4)*dy
  enddo   
#endif
  do k = 2, nz
    zz(k) = zz(k-1) + 0.5*dz*(fdz0(k)+fdz0(k-1))
  enddo
!
!----------------------------------------------------------!
!    		      Start main loop    	           !
!----------------------------------------------------------!
!
  do n = 1, nparc_loc
    p = parcel_loc(n) 
!
! Calculate distance function or interpolation parameters 
!
    if (lag_ord == 0) then
!
      filter = 0.
      r_filter = sqrt(0.7)*min(dx,dy)
      do k = -1, 1
        do j = -1, 1
          do i = -1, 1
            ff = (parcel_pos(p)%x - xx(parcel_pos(p)%i+i))**2. + (parcel_pos(p)%y - yy(parcel_pos(p)%j+j))**2.  &
                    + (parcel_pos(p)%z - zz(parcel_pos(p)%k+k))**2.
            filter(i+2,j+2,k+2) = exp( -ff/r_filter**2. )
          enddo
        enddo
      enddo
      filter = filter / sum(filter)
!
    else
!
!  Initialize interpolations (use quadratic lagrange polynomials)
!
      ind(1,:) = 2  ! Linear interpolations
      ind(2,:) = 3  ! Quadratic interpolations
!
      call lagrange_interp_nd_size (dim, ind(1,:), nd(1))
      call lagrange_interp_nd_size (dim, ind(2,:), nd(2))
!
      if ( parcel_pos(p)%exist ) then
#ifdef MODEL_3D
        a(1) = xx(parcel_pos(p)%i-1)
        a(2) = yy(parcel_pos(p)%j-1)
        a(3) = zz(max(parcel_pos(p)%k-1,1))
        b(1) = xx(parcel_pos(p)%i+1)
        b(2) = yy(parcel_pos(p)%j+1)
        b(3) = zz(min(parcel_pos(p)%k+1,nz))
#else
        a(1) = xx(parcel_pos(p)%i-1)
        a(2) = zz(max(parcel_pos(p)%k-1,1))
        b(1) = xx(parcel_pos(p)%i+1)
        b(2) = zz(min(parcel_pos(p)%k+1,nz))
#endif
      endif
    endif
!
!----------------------------------------------------------!
!    		    Parcel interpolations   	           !
!----------------------------------------------------------!
!
!  Velocity vector and SGS TKE
!
    if (lag_ord /= 0) then
      call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, wind%u, parcel_pos(p), parcel_sca(p)%u )
    else
      call simple_interpolation ( filter, parcel_pos(p), wind%u, parcel_sca(p)%u )
    endif
!
#ifdef MODEL_3D
    if (lag_ord /= 0) then
      call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, wind%v, parcel_pos(p), parcel_sca(p)%v )
    else
      call simple_interpolation ( filter, parcel_pos(p), wind%v, parcel_sca(p)%v )
    endif
#endif
!
    if (lag_ord /= 0) then
      call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, wav, parcel_pos(p), parcel_sca(p)%w )
    else
      call simple_interpolation ( filter, parcel_pos(p), wind%w, parcel_sca(p)%w )
    endif
!
    if (lag_ord /= 0) then
      call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, turbu%ksgs, parcel_pos(p), parcel_sca(p)%tke )
    else
      call simple_interpolation ( filter, parcel_pos(p), turbu%ksgs, parcel_sca(p)%tke )
    endif
!
!  Scalars only interpolated if output or initialization
!
    if ( output .or. aerosol_lag ) then
!
      if (lag_ord /= 0) then
        call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, thermo%t, parcel_pos(p), parcel_sca(p)%t )
      else
        call simple_interpolation ( filter, parcel_pos(p), thermo%t, parcel_sca(p)%t )
      endif
!
      if (lag_ord /= 0) then
        call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, state%es, parcel_pos(p), parcel_sca(p)%pt )
      else
        call simple_interpolation ( filter, parcel_pos(p), state%es, parcel_sca(p)%pt )
      endif
!
      if (lag_ord /= 0) then
        call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, state%qt, parcel_pos(p), parcel_sca(p)%qt )
      else
        call simple_interpolation ( filter, parcel_pos(p), state%qt, parcel_sca(p)%qt )
      endif
!
      if (lag_ord /= 0) then
        call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, thermo%buoy, parcel_pos(p), parcel_sca(p)%buoy )
      else
        call simple_interpolation ( filter, parcel_pos(p), thermo%buoy, parcel_sca(p)%buoy )
      endif
!
      if (out_sat) then
        if (lag_ord /= 0) then
          call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, thermo%rh, parcel_pos(p), parcel_sca(p)%rh )
        else
          call simple_interpolation ( filter, parcel_pos(p), thermo%rh, parcel_sca(p)%rh )
        endif
      endif
!
      if (lmicro > 0) then
        if (lag_ord /= 0) then
          call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, hydromtr(drop)%q+hydromtr(rain)%q,    &
   	       parcel_pos(p), parcel_sca(p)%ql )
        else
          call simple_interpolation ( filter, parcel_pos(p), hydromtr(drop)%q+hydromtr(rain)%q, parcel_sca(p)%ql )
        endif
      endif
!
      if (lmicro > 1) then
        if (lag_ord /= 0) then
          call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, hydromtr(ice)%q, parcel_pos(p), parcel_sca(p)%qi )
        else
          call simple_interpolation ( filter, parcel_pos(p), hydromtr(ice)%q, parcel_sca(p)%qi )
        endif
      endif
!
!  passive scalars
!
      do k = 1, nscal
        if (lag_ord /= 0) then
          call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, state%scal(:,:,:,k), parcel_pos(p), parcel_sca(p)%scal(k) )
        else
          call simple_interpolation ( filter, parcel_pos(p), state%scal(:,:,:,k), parcel_sca(p)%scal(k) )
        endif
      enddo
!
!  Tendencies
!
      if (out_diagw) then
        do h = 1, 9
          if (lag_ord /= 0) then
            call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, diag(h)%w, parcel_pos(p), parcel_sca(p)%dw(h) )
          else
            call simple_interpolation ( filter, parcel_pos(p), diag(h)%w, parcel_sca(p)%dw(h) )
          endif
        enddo
      endif
!
      if (out_diagl) then
        do h = 1, 9
          if (lag_ord /= 0) then
            call lagrange_interpolation ( dim, nd(l), ind(l,:), a, b, diag(h)%qc+diag(h)%qr, parcel_pos(p), parcel_sca(p)%dq(h) )
          else
            call simple_interpolation ( filter, parcel_pos(p), diag(h)%qc+diag(h)%qr, parcel_sca(p)%dq(h) )
          endif
        enddo
      endif
!
    endif
!
!  End main loop
!
  enddo
!
!----------------------------------------------------------!
!    		    Set particle mask   	           !
!----------------------------------------------------------!
!
  if ( output .and. lmask ) then
    if (mypid==0.and.lfirst) nchange_tot = 0
    if (mypid==0.and.lfirst) nmask_tot = 0
    nchange = 0
    nmask = 0
    nchangep = 0
    nmaskp = 0
!
    do p = 1, nparc
      if ( mypid == parcel_pos(p)%iproc ) then
        call set_mask (p, parcel_pos(p), parcel_sca(p), nchange, nmask, lfirst)
      endif
    enddo
  endif
!
!  Deallocate
!
  deallocate (ind, a, b)
!
!----------------------------------------------------------!
!
  if (verbose > 2) call write_debug('Terminating lagrange_interp')
!
  return
  end
!
!  ===================
  subroutine simple_interpolation ( filter, pos, x, val )
!  ===================
!
  integer :: i, j, k
  type(position) :: pos
!
  real :: val
  real, dimension(1:3,1:3,1:3) :: filter
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x
!
!----------------------------------------------------------!
!    	     Interpolate x at parcel location    	   !
!----------------------------------------------------------!
!    
    val = 0.
    do k = -1, 1
      do j = -1, 1
        do i = -1, 1
          val = val + x(pos%i+i,pos%j+j,pos%k+k)*filter(i+2,j+2,k+2)
        enddo
      enddo
    enddo
!
  return
  end  
!
!  ===================
  subroutine lagrange_interpolation ( dim, nd, ind, a, b, sol3d, pos, val )
!  ===================
!
  integer :: n, i, j, k, ii, jj, kk, step, offset
  integer :: dim, nd, ind(dim)
!
  real :: val, a(dim), b(dim)
  real :: xd(dim,nd), xi(dim), zd(nd)
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: sol3d
  type(position) :: pos
!
!----------------------------------------------------------!
!    	    Interpolate sol at parcel location    	   !
!----------------------------------------------------------!
!
!  call lagrange_interp_nd_grid ( dim, ind, a, b, nd, xd )
!   
  xi(1) = pos%x
#ifdef MODEL_3D
  xi(2) = pos%y
  xi(3) = pos%z
#else
  xi(2) = pos%z
#endif
!
  if (ind(1) == 2) then
    step = 1
    offset = 2
  else if (ind(1) == 3) then
    step = 1
    offset = 1
  endif
!
  n = 0
#ifdef MODEL_3D
  do k = -step, step, offset
    do j = -step, step, offset
      do i = -step, step, offset
        n = n + 1
        ii = pos%i+i
        jj = pos%j+j
	kk = pos%k+k
        zd(n) = sol3d(ii,jj,kk)
      enddo
    enddo
  enddo
#else
  do k = -step, step, offset
    do i = -step, step, offset
      n = n + 1
      ii = pos%i+i
      kk = pos%k+k
      zd(n) = sol3d(ii,1,kk)
    enddo
  enddo
#endif
!
  call lagrange_interp_nd_value ( dim, ind, a, b, nd, zd, 1, xi, val )
!
!  val = sol3d(pos%i,1,pos%k)
!
  return
  end  
!
!  ===================
  subroutine lagrange_proc (init)
!  ===================
!
#if ( defined SPMD )
  integer, dimension(8) :: req 
  integer :: status(MPI_STATUS_SIZE)
  integer :: tag, ierr
#endif
!
  logical :: init
  logical :: exist_old(1:nparc)
  integer :: i, k, n, np, iproc, iprocy, ncellx, ncelly
  integer :: iproc_old(1:nparc)
!
  if (verbose > 1) call write_debug('Start lagrange_proc')
!
!----------------------------------------------------------!
!    	  	  Parcel location on grid   		   !
!----------------------------------------------------------!
!
  iproc = 0
  iproc_old = parcel_pos%iproc
  exist_old = parcel_pos%exist
!
  ncellx = (nx-5)/nprocx
  ncelly = max((ny-5)/nprocy,1)
!
  do n = 1, nparc_loc 
    i = parcel_loc(n)
!
    parcel_pos(i)%i = FLOOR(parcel_pos(i)%x/dx) + 4
    parcel_pos(i)%j = FLOOR(parcel_pos(i)%y/dy) + 4
    do k = 1, nz-1
      if (parcel_pos(i)%z >= z0(k) .and. parcel_pos(i)%z < z0(k+1)) parcel_pos(i)%k = k
    enddo
!
!  Find processor
!
#ifdef SPMD 
    iproc = FLOOR( real(parcel_pos(i)%i-4)/real(ncellx) )
    if (nprocy > 1) then
      iprocy = FLOOR( real(parcel_pos(i)%j-4)/real(ncelly) )
      iproc = iprocy + iproc*nprocy
    endif
    parcel_pos(i)%iproc = iproc
#else
    parcel_pos(i)%iproc = 0
#endif
!
    if (parcel_pos(i)%iproc < 0 .or. parcel_pos(i)%iproc >= nproc .or. .not.exist_old(i)) parcel_pos(i)%exist = .false.
  enddo
!
!----------------------------------------------------------!
!                    MPI communications                    !
!----------------------------------------------------------!
!
!  Communicate parcels absolute location
!
#ifdef SPMD  
  if (.not.init) then
    tag = 100
    do i = 1, nparc
      if ( exist_old(i) ) then
!
!  Broadcast parcels proc location
!
      call MPI_Bcast(parcel_pos(i)%iproc,1,INTTYPE,iproc_old(i),MPI_COMM_WORLD,ierr)
      call MPI_Bcast(parcel_pos(i)%exist,1,LOGTYPE,iproc_old(i),MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!  Send parcel positions
!
      if ( parcel_pos(i)%iproc /= iproc_old(i) .and. parcel_pos(i)%exist ) then
        if ( mypid == iproc_old(i) ) then
          call MPI_ISEND ( parcel_pos(i)%x, 1, REALTYPE, parcel_pos(i)%iproc,   tag, MPI_COMM_WORLD, req(1), ierr )	    
          call MPI_ISEND ( parcel_pos(i)%y, 1, REALTYPE, parcel_pos(i)%iproc, tag+1, MPI_COMM_WORLD, req(2), ierr )     
          call MPI_ISEND ( parcel_pos(i)%z, 1, REALTYPE, parcel_pos(i)%iproc, tag+2, MPI_COMM_WORLD, req(3), ierr )     
!
          call MPI_ISEND ( parcel_pos(i)%i, 1, INTTYPE, parcel_pos(i)%iproc, tag+3, MPI_COMM_WORLD, req(4), ierr )      
          call MPI_ISEND ( parcel_pos(i)%j, 1, INTTYPE, parcel_pos(i)%iproc, tag+4, MPI_COMM_WORLD, req(5), ierr )      
          call MPI_ISEND ( parcel_pos(i)%k, 1, INTTYPE, parcel_pos(i)%iproc, tag+5, MPI_COMM_WORLD, req(6), ierr )      
!
  	  if (aerosol_lag) then
	    call MPI_ISEND ( aer(i)%r, 1, REALTYPE, parcel_pos(i)%iproc, tag+6, MPI_COMM_WORLD, req(7), ierr )	  
	    call MPI_ISEND ( aer(i)%n_s, 1, REALTYPE, parcel_pos(i)%iproc, tag+7, MPI_COMM_WORLD, req(8), ierr )
	  endif
!
	  if (.not.aerosol_lag) then
            call MPI_WAITALL ( 6, req, MPI_STATUSES_IGNORE, ierr )
	  else
            call MPI_WAITALL ( 8, req, MPI_STATUSES_IGNORE, ierr )
	  endif
!
!  Receive parcel positions
!
        else if ( mypid == parcel_pos(i)%iproc ) then
          call MPI_RECV ( parcel_pos(i)%x, 1, REALTYPE, iproc_old(i),   tag, MPI_COMM_WORLD, status, ierr )      
          call MPI_RECV ( parcel_pos(i)%y, 1, REALTYPE, iproc_old(i), tag+1, MPI_COMM_WORLD, status, ierr )      
          call MPI_RECV ( parcel_pos(i)%z, 1, REALTYPE, iproc_old(i), tag+2, MPI_COMM_WORLD, status, ierr )      

          call MPI_RECV ( parcel_pos(i)%i, 1, INTTYPE, iproc_old(i), tag+3, MPI_COMM_WORLD, status, ierr )      
          call MPI_RECV ( parcel_pos(i)%j, 1, INTTYPE, iproc_old(i), tag+4, MPI_COMM_WORLD, status, ierr )      
          call MPI_RECV ( parcel_pos(i)%k, 1, INTTYPE, iproc_old(i), tag+5, MPI_COMM_WORLD, status, ierr )      
!
  	  if (aerosol_lag) then
            call MPI_RECV ( aer(i)%r, 1, REALTYPE, iproc_old(i), tag+6, MPI_COMM_WORLD, status, ierr )
	    call MPI_RECV ( aer(i)%n_s, 1, REALTYPE, iproc_old(i), tag+7, MPI_COMM_WORLD, status, ierr )
	  endif	   
        endif
      endif
!
      endif
    enddo
  endif
#endif
!
!----------------------------------------------------------!
!                 Count particles on procs                 !
!----------------------------------------------------------!
!
  nparc_loc = 0
  parcel_loc = 0
  do i = 1, nparc
    if ( parcel_pos(i)%iproc == mypid .and. parcel_pos(i)%exist ) then
      nparc_loc = nparc_loc + 1
      parcel_loc(nparc_loc) = i
    endif
  enddo
!
!----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminate lagrange_proc')
!
  return
  end  
!
!  ===================
  subroutine set_mask ( i, pos, scal, nchange, nmask, lfirst )
!  ===================
!
  integer :: i, mask_old, nmask(3), nchange(6)
  logical :: lfirst
  type(scalar) :: scal
  type(position) :: pos
!
!----------------------------------------------------------!
!    	    	   Set parcel's mask    		   !
!----------------------------------------------------------!
!
  if ( scal%ql > qthres .and. scal%w > wthres ) then
    scal%mask = 1
  else if ( scal%ql > qthres .and. scal%w <= wthres ) then
    scal%mask = 2
  else
    scal%mask = 0
  endif
!
  return
  end  
!
!  ===================
  subroutine reset_particles (I)
!  ===================
!
  INTEGER, intent(in) :: I
!
  INTEGER :: J, n, IS, IE, SEED
  INTEGER :: iproc, iprocy, err
  INTEGER, DIMENSION (8) :: T
!
  IF (I == 0) then
#ifdef SPMD
    IS = 1 + mypid*nparc/nproc
    IE = (mypid + 1)*nparc/nproc
#else
    IS = 1
    IE = nparc
#endif
    iproc = mypid    
    nparc_loc = 0
  ELSE
    IS = I
    IE = I
    iproc = FLOOR( nproc*(I - 1)/real(nparc) )
  ENDIF
!
  CALL DATE_AND_TIME(VALUES = T)
  SEED = 9*mypid+IS*T(1)+70*(T(2)+12*(T(3)+31*(IE*T(5)+23*(T(6)+59*T(7)))))
  IF (MOD(SEED,2).EQ.0) SEED = SEED-1
!
!----------------------------------------------------------!
!         Select random initialization method
!----------------------------------------------------------!
!
  if ( mypid == iproc ) then
    do J = IS, IE
      select case (lag_init)
        case (1)
          call bl_particles (parcel_pos(J), SEED, it_start, it_end, jt_start, jt_end)
        case (2)
          call layer_particles (parcel_pos(J), SEED, it_start, it_end, jt_start, jt_end)
        case (3)
          call line_particles (parcel_pos(J), SEED, it_start, it_end, jt_start, jt_end)
      end select
!
      nparc_loc = nparc_loc + 1
      parcel_loc(nparc_loc) = J
      parcel_pos(J)%iproc = iproc
      parcel_pos(J)%exist = .true.
    enddo
  endif
!
!----------------------------------------------------------!
!             Broadcast new parcels location               !
!----------------------------------------------------------!
!
#ifdef SPMD
  if (I == 0) then
    do J = 1, nparc
      if (mypid /= 0 .and. mypid == parcel_pos(J)%iproc) then
        call MPI_SEND(parcel_pos(J)%iproc, 1, INTTYPE, 0, 1, MPI_COMM_WORLD, err)
        call MPI_SEND(parcel_pos(J)%exist, 1, LOGTYPE, 0, 2, MPI_COMM_WORLD, err)
      else if (J > nparc_loc .and. mypid == 0) then
        CALL MPI_RECV(parcel_pos(J)%iproc, 1, INTTYPE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
        CALL MPI_RECV(parcel_pos(J)%exist, 1, LOGTYPE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      endif
!
      call MPI_BARRIER(MPI_COMM_WORLD, err)
      call MPI_BCAST(parcel_pos(J)%iproc, 1, INTTYPE, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(parcel_pos(J)%exist, 1, LOGTYPE, 0, MPI_COMM_WORLD, err)
    enddo
  else
    if (mypid /= 0 .and. mypid == parcel_pos(I)%iproc) then
      call MPI_SEND(parcel_pos(I)%iproc, 1, INTTYPE, 0, 1, MPI_COMM_WORLD, err)
      call MPI_SEND(parcel_pos(I)%exist, 1, LOGTYPE, 0, 2, MPI_COMM_WORLD, err)
    else if (I > nparc_loc .and. mypid == 0) then
      CALL MPI_RECV(parcel_pos(I)%iproc, 1, INTTYPE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      CALL MPI_RECV(parcel_pos(I)%exist, 1, LOGTYPE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
    endif
!
    call MPI_BARRIER(MPI_COMM_WORLD, err)
    call MPI_BCAST(parcel_pos(I)%iproc, 1, INTTYPE, 0, MPI_COMM_WORLD, err)
    call MPI_BCAST(parcel_pos(I)%exist, 1, LOGTYPE, 0, MPI_COMM_WORLD, err)
  endif
#endif
!
  return
!
  contains
!
  subroutine bl_particles (pos, seed, is, ie, js, je)
    
    INTEGER :: SEED, is, ie, js, je
    REAL :: X,Y,Z,RANF      
    TYPE(position) :: pos
      
    X=RANF(SEED)
    Y=RANF(SEED)
    Z=RANF(SEED)
!
    pos%x=real(is - 4)*dx + real(ie - is + 1)*dx*X
#ifdef MODEL_3D
    pos%y=real(js - 4)*dy + real(je - js + 1)*dy*Y
#endif
    pos%z=zbl*Z
      
    return
  end subroutine bl_particles
!
  subroutine line_particles (pos, seed, is, ie, js, je)
    
    INTEGER :: SEED, is, ie, js, je
    REAL :: X,Y,Z,RANF      
    TYPE(position) :: pos
      
    X=RANF(SEED)
    Y=RANF(SEED)
!
    pos%x=real(is - 4)*dx + real(ie - is + 1)*dx*X
#ifdef MODEL_3D
    pos%y=real(js - 4)*dy + real(je - js + 1)*dy*Y
#endif
    pos%z=zbl
      
    return
  end subroutine line_particles
!
  subroutine layer_particles (pos, seed, is, ie, js, je)
    
    INTEGER :: SEED, is, ie, js, je
    REAL :: X,Y,Z,RANF      
    TYPE(position) :: pos
      
    X=RANF(SEED)
    Y=RANF(SEED)
    Z=RANF(SEED)
!
    pos%x=real(is - 4)*dx + real(ie - is + 1)*dx*X
#ifdef MODEL_3D
    pos%y=real(js - 4)*dy + real(je - js + 1)*dy*Y
#endif
    pos%z=zl1 + (zl2-zl1)*Z

    return
  end subroutine layer_particles
!
end subroutine reset_particles
!
! =======================
  subroutine reset_aerosols( a )
! =======================
!
    integer :: i, n, SEED, T(8)
    real :: logn_value
    type(aerosol), dimension(1:nlag) :: a
    real, external :: r8_normal_ab
!
!  Chemical properties of aerosols
!
    real :: M_aer
    real :: rho_aer
!
!  Calculate parameters of dry aerosol
!
      if (compos_lag == 'NaCl') then
	M_aer	      = 58.4e-3 		      ! kg mol-1	      Molar weight if compound (NaCl)
	rho_aer       = 2165.0  		      ! kg m-3  	      Density of compound (NaCl)
      end if
!
!  Reinitialize seed 
!
      CALL DATE_AND_TIME(VALUES = T)
      SEED = 23*mypid + ntime*T(1)+70*(T(2)/ntime+12*i*(T(3)+31*(T(5)+23*(T(6)/ntime+59*T(7)))))
      IF (MOD(SEED,2).EQ.0) SEED = SEED-1
!
!  If aerosol_lag =true initialize log-normal distribution of particles
!
      do n = 1, nparc
	i = n
!
!  Draw from log-normal distribution with mean =mu_lag and standard deviation =sigma_lag     
!
	logn_value =  EXP( mu_lag + sigma_lag * r8_normal_ab ( mu_lag, sigma_lag, SEED ) )
	a(i)%r0 =     logn_value * 1.0e-6
	a(i)%r =	      a(i)%r0
!
	a(i)%m_s =    (4./3.) * pi * a(i)%r**3 * rho_aer
	a(i)%n_s =    a(i)%m_s / M_aer
      enddo
!
    return
  end subroutine reset_aerosols
!
! =======================
  subroutine lagrange_aerosol
! =======================
!
    integer :: i, n
    real, parameter :: aSigma_w = 7.28e-2     ! N m-1		      Droplet surface tension
    real, parameter :: aRho_w = 1000.0  	      ! kg m-3  	      Density of water
    real, parameter :: aM_w = 18.015e-3 	      ! kg mol-1	      Molecular weight of water
    real, parameter :: aRho_d = 1000.0  	      ! kg m-3  	      Density of droplet
    real, parameter :: aD_v = 0.282e-4  	      ! m2 s-1  	      Coefficient of diffusion for vapor
    real, parameter :: aL_v = 40660		      ! J mol-1 	      Enthalpy of vaporization
    real, parameter :: aK_t = 0.024		      ! W m-1 K-1	      Thermal conductivity for air
    real, parameter :: aR = 8.314		      ! J mol-1 K-1   Ideal gas constant
    real, parameter :: aR_v = 461.5		      ! K kg-1 K-1    Specific gas constant for water vapour
    real, parameter :: aLambda_0 = 6.62e-8    ! m			      Mean free path at, std atm
    real, parameter :: aP_0 = 101300	      ! Pa		      Pressure, std atm
    real, parameter :: aEtha_0 = 18.27e-6     ! Pa s-1  	      Viscosity for air, std atm
    
    real :: aS  							      ! 1			      Saturation of vapor
    real :: aT  							      ! K			      Temperature
    
    real,parameter :: E = 0.1			      ! 1		      Collision efficiency
    
    real :: a_aer, b_aer, c_aer, d_aer, c1_aer, c2_aer
    real :: collision_term, res
    
!
!  For each aerosol solv set of equations for condensational growth and velocity
!
      do n = 1, nparc_loc
	i = parcel_loc(n)
!
!  Calculate temperature and saturation
!
	aT	      = parcel_sca(i)%pt * (p0(parcel_pos(i)%k)/pref)**((cp_a-cv_a)/cp_a)
	aS	      = parcel_sca(i)%rh / 100.
!
!  Calculate parameters for the equations for aersosols
!
	a_aer	      = (2. * aM_w * aSigma_w) / (aR * aT * aRho_w)
	b_aer	      = (3. * aM_w * aer(i)%n_s) / (4. * pi * aRho_w)
	c_aer	      = aRho_d / (72. * aEtha_0)
	d_aer	      = 1.255 * aLambda_0
	c1_aer        = (aRho_w * aR * aT) / (aM_w * aD_v * cal_esw(aT))
	c2_aer        = (aRho_w * aL_v / (aM_w * aK_t * aT))*( (aL_v/(aR*aT)) - 1. )
!
!  Particle condensational growth
!
	res =	      0.
      res =   aer(i)%r + dt0 * ( (aS - 1.) - a_aer/aer(i)%r + b_aer/aer(i)%r**3 ) / ( aer(i)%r * (c1_aer + c2_aer) )
	if (res > aer(i)%r0) then
	  aer(i)%r = res
	  aer(i)%ws = -c_aer * (1 + d_aer/aer(i)%r) * aer(i)%r**2
	endif
      enddo
!
    return
  end subroutine lagrange_aerosol
!
#endif

end module lagrange
