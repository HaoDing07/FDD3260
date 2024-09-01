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
!  RESTART.F90                   
!
!  Purpose:
!	For restarting a run from a previous simulation			  
!
!  Author
!	Julien Savre, MISU
!
! ===========================================================

module restart

	USE shared_all
    	USE thermodynamics
	USE radiationmod
	USE boundary_conditions
	USE shared_lagrange
	USE allocation
	USE sources
	USE averages
	USE initialize
        USE lagrange
	USE aerosols
#ifdef SPMD
    	USE mpi
    	USE mpicomm
#endif

	IMPLICIT NONE
	
	integer :: ivers
	character(len = 3)  :: vers = 'V04'
        
  	private
    
  	public :: savedata, readdata
  
  	contains

!  ===========================================================
    Subroutine readdata
!  ===========================================================

	IMPLICIT NONE
!
	integer :: h, i, j, k, l, ierr, ininuc, nxp, nyp, nzp, ntaup, nscalp, nparcp, im
	logical :: interpolate, ex, expart
        character(len = 100):: filename
	character(len = 3)  :: vers
	real :: dtp, dxp, dyp, dzp, xn_in00, xn_ccn00
!
	real, dimension(:,:), allocatable :: scal0p
	real, dimension(:), allocatable :: xx, xxp, yy, yyp, z0p, fdz0p
	real, dimension(:), allocatable :: u00p,v00p,p0p,den0p,stab0p,n20p,g0p,avden0p, 	&
					   pt0p,mse0p,cv0p,t0p,qt0p,qv0p,ql0p,nl0p,w0p,k0p
	real, dimension(:,:,:), allocatable :: tmp
	character :: char_helper
!
        if (verbose > 0) call write_debug('Enter readdata')
!
!-----------------------------------------------------
!
!  Open restart files
!
	write(7,*)''
!
	if (mypid == 0) then
          filename = opdir(1:len_trim(opdir)) // '/nuclei.dat'
          if (lndrop /= 0) open(104,file=filename(1:len_trim(filename)), form='unformatted',status='unknown')
!
	  filename = opdir(1:len_trim(opdir)) // '/' // file_rest(1:len_trim(file_rest))
          open(106,file=filename(1:len_trim(filename)), form='unformatted',status='unknown')
!
	  read(106) vers
	  if ( .not.(vers(1:1) == 'V') ) then
	    ivers = 0
	    rewind(106)
	  else
	    read(vers(2:3),'(i2)') ivers
	  endif
! 
	  read(106)time,ntime,dtp,ntaup,nxp,nyp,nzp,dxp,dyp,dzp,nscalp,nparcp,xn_ccn00,xn_in00
	endif
!
!  Do we need to interpolate?
!
	interpolate=.false.
	if (mypid == 0 .and. (nxp /= nx .or. nyp /= ny)) then
	  write(7,*) 'WARNING: Non-conforming dimensions in the horizontal, will need to interpolate saved data on new grid'
	  interpolate=.true.
	endif
	if (mypid == 0 .and. nzp /= nz) then
	  write(7,*) 'WARNING: Non-conforming dimensions in the vertical, will need to interpolate saved data on new grid'
	  interpolate=.true.
	endif
!
#if ( defined SPMD )
        call MPI_BCAST(interpolate, 1,  LOGTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ivers, 1,  INTTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(ntime, 1,  INTTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(nscalp,1,  INTTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(nparcp,1,  INTTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(nxp,  1,  INTTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(nyp,  1,  INTTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(nzp,  1,  INTTYPE,0,MPI_COMM_WORLD,ierr)
  	call MPI_BCAST(time, 1,  REALTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(dxp,  1,  REALTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(dyp,  1,  REALTYPE,0,MPI_COMM_WORLD,ierr)
 	call MPI_BCAST(dzp,  1,  REALTYPE,0,MPI_COMM_WORLD,ierr)
#endif
!
	if (tstop <= time) then
          if (mypid == 0) then
	    write(7,*) 'ERROR RESTARTING: time >= tstop, (', time, '>=', tstop, ')'
	  endif
          call stop_mimica ( 1 ) 
	endif
	tstart = time
!
!  Define grid point locations in the horizontal
!
	allocate( xx(1:nx), xxp(1:nxp), yy(1:ny), yyp(1:nyp), z0p(1:nzp), fdz0p(1:nzp) )
	allocate( tmp(ip_start:ip_end,jp_start:jp_end,1:nz) )
!
	xxp(1) = 0.
	xx(1) = 0.
	if ( xc_nest > 0. .and. nest_run ) then
	  xx(1) = xc_nest - dx*real(nx/2)
	endif
	if ( nx /= nxp ) then
	  xx(1) = xc_nest - dx*real(nx/2)
	  xxp(1) = xc_nest - dxp*real(nxp/2)
	endif
!
	do i = 2, nx
	  xx(i) = xx(i-1) + dx
	enddo
	do i = 2, nxp
	  xxp(i) = xxp(i-1) + dxp
	enddo
!
#ifdef MODEL_3D
	yyp(1) = 0.
	yy(1) = 0.
	
	if ( yc_nest > 0. .and. nest_run ) then
	  yy(1) = yc_nest - dy*real(ny/2)
	endif
	if ( ny /= nyp ) then
	  yy(1) = yc_nest - dy*real(ny/2)
	  yyp(1) = yc_nest - dyp*real(nyp/2)
	endif
!
	do j = 2, ny
	  yy(j) = yy(j-1) + dy
	enddo
	do j = 2, nyp
	  yyp(j) = yyp(j-1) + dyp
	enddo
#endif
!
!  Define grid point locations in the vertical
!
        if ( mypid == 0 ) then
          inquire ( FILE=trim(gridfile), EXIST=ex )
	  call define_grid (nzp, dzp, fdz0p, z0p)
          if (.not.ex) then
            call define_grid (nz, dz, fdz0, z0)
	  else
	    dz = dzp; z0 = z0p; fdz0 = fdz0p
	  endif
	endif
!	  
#if ( defined SPMD )
        call MPI_BCAST(dzp,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)  
        call MPI_BCAST(z0p,   nzp, REALTYPE,0,MPI_COMM_WORLD,ierr)  
        call MPI_BCAST(fdz0p, nzp, REALTYPE,0,MPI_COMM_WORLD,ierr)  
!
        call MPI_BCAST(dz,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)  
        call MPI_BCAST(z0,     nz, REALTYPE,0,MPI_COMM_WORLD,ierr)  
        call MPI_BCAST(fdz0,   nz, REALTYPE,0,MPI_COMM_WORLD,ierr)  
#endif
!
!  Read and interpolate initial soundings
!
	allocate( u00p(1:nzp), v00p(1:nzp), p0p(1:nzp), den0p(1:nzp), stab0p(1:nzp), n20p(1:nz), g0p(1:nzp), 	 &
		  avden0p(1:nzp), pt0p(1:nzp), mse0p(1:nzp), cv0p(1:nzp), t0p(1:nzp), qt0p(1:nzp), 	 &
		  qv0p(1:nzp), ql0p(1:nzp), nl0p(1:nzp), w0p(1:nzp), k0p(1:nzp), scal0p(1:nzp,1:nscalp) )
!
	if (mypid == 0) then
	  if ( ivers > 1 .and. ivers < 4 ) then
	    read(106)u00p,v00p,p0p,den0p,stab0p,g0p,avden0p,pt0p,mse0p,t0p,qt0p,qv0p,ql0p,nl0p,w0p,k0p,scal0p
	  else if ( ivers >= 4 ) then
	    read(106)u00p,v00p,p0p,den0p,stab0p,n20p,g0p,avden0p,pt0p,mse0p,t0p,qt0p,qv0p,ql0p,nl0p,w0p,k0p,scal0p
	  else
	    read(106)u00p,v00p,p0p,den0p,stab0p,g0p,avden0p,pt0p,cv0p,t0p,qt0p,qv0p,ql0p,nl0p,w0p,k0p,scal0p
	  endif
!
	  call interp_1d (interpolate, nz, nzp, z0, z0p, u00, u00p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, v00, v00p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, p0, p0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, den0, den0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, g0, g0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, stab0, stab0p)
	  if ( ivers >= 4 ) call interp_1d (interpolate, nz, nzp, z0, z0p, n20, n20p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, avden0, avden0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, pt0, pt0p)
	  if ( ivers > 1 ) then
	    call interp_1d (interpolate, nz, nzp, z0, z0p, mse0, mse0p)
	  else
	    call interp_1d (interpolate, nz, nzp, z0, z0p, cv0, cv0p)
	  endif
	  call interp_1d (interpolate, nz, nzp, z0, z0p, t0, t0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, qt0, qt0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, qv0, qv0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, ql0, ql0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, nl0, nl0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, w0, w0p)
	  call interp_1d (interpolate, nz, nzp, z0, z0p, k0, k0p)
  	  do h = 1, nscal
	    if (h <= nscalp) then
	      call interp_1d (interpolate, nz, nzp, z0, z0p, scal0(:,h), scal0p(:,h))
	    else
	     scal0(:,h) = 0.
	    endif
	  enddo
	endif
!	
	deallocate(u00p,v00p,p0p,den0p,stab0p,n20p,g0p,avden0p,pt0p,mse0p,cv0p,t0p,qt0p,qv0p,ql0p,nl0p,w0p,k0p,scal0p)
!
	if (verbose > 2) call write_debug('Read dimensions and reference profiles')
!
! Read aerosol properties
!
#ifdef AERO_ENABLE
	if (mypid == 0) then
	  nmode = nmode0
	  if (reg_mode) nmode = nmode +1
	  allocate(aero(1:nmode), aero0(1:nmode))
	  do im = 1, nmode
		read(106)aero(im)%nelem
		do i=1,4
			read(106)aero(im)%init%present(i)
			read(106)aero(im)%init%frac(i)
		enddo
		read(106)aero(im)%init%n0
		if(aero_sfc_source) read(106)aero(im)%init%n_sfc_source
		read(106)aero(im)%size%rmean,aero(im)%size%mmean,aero(im)%size%sigma
		read(106)aero(im)%mix%kappa,aero(im)%mix%rho,aero(im)%mix%mw
		read(106)aero(im)%depv0
		read(106)aero0(im)%n
		read(106)aero0(im)%m
		read(106)aero0(im)%ma
	  enddo
	endif
!
!  Broadcast
!
#ifdef SPMD
	call MPI_BCAST(nmode,             1, INTTYPE,0,MPI_COMM_WORLD,ierr)
        if (mypid > 0) allocate(aero(1:nmode), aero0(1:nmode))
	do im = 1, nmode
		call MPI_BCAST(aero(im)%nelem,               1, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%init%present,        4, LOGTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%init%frac,           4, REALTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%init%n0,             1, REALTYPE,0,MPI_COMM_WORLD,ierr)
		if(aero_sfc_source) call MPI_BCAST(aero(im)%init%n_sfc_source,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%size%rmean,          1, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%size%mmean,          1, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%size%sigma,          1, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%mix%kappa,           1, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%mix%rho,             1, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%mix%mw,              1, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero(im)%depv0,               1, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero0(im)%n,                 nz, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero0(im)%m,                 nz, INTTYPE,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(aero0(im)%ma,                nz, INTTYPE,0,MPI_COMM_WORLD,ierr)
	enddo
#endif
#endif
!
!  Read and interpolate restart file: parallel case
!
        call load_data_u (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, wind2%u, wind2%v, wind2%w)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, pressure2%p)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, state2%es)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, state2%qt)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, pressure2%dens)
!
	do h = 1, nhydro
          call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, hydromtr2(h)%q)
!
          call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, hydromtr2(h)%n)
!
          if ( lmicro > 3 .and. h > ice ) call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, hydromtr2(h)%w)
	enddo
!
	do h = 1, nscal
	  state2%scal(ip_start:ip_end,jp_start:jp_end,1:nz,h) = 0.
!
	  if (h <= nscalp) then
            call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, state2%scal(ip_start:ip_end,jp_start:jp_end,1:nz,h))
	  else
	    call scalar_reset (h, state2%scal(ip_start:ip_end,jp_start:jp_end,1:nz,h), zbl, .true.)
	    !call horav (state2%scal(ip_start:ip_end,jp_start:jp_end,1:nz,h), scal)
	  endif
	enddo
!
	if (nscal < nscalp) then
	  do h = nscal+1, nscalp
            call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, tmp)
	  enddo
	endif
!
	if ( ivers > 0 ) then
          if ( ivers <= 2 ) call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, tmp)
!	
          call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, turbu%ksgs)
	endif
!
#ifdef SEIFERT
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, nuc2%ccn)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, nuc2%in)
#endif
!
        if (with_piggy) then
          call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, statep%es)
!
          call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, statep%qt)
!
	  do h = 1, nhydro
            call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, hydromtrp(h)%q)
!
            call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, hydromtrp(h)%n)
!
            if ( lmicro > 3 .and. h > ice ) call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, hydromtrp(h)%w)
	  enddo
        endif
!
#if (defined AERO_ENABLE)
	do h = 1, nmode
          call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aero3d2(h)%n)

          call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aero3d2(h)%m)

          call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aero3d2(h)%ma)
	enddo
#endif
!
#ifdef NUC_CNT
	do h = 1, 3
	  do l = 1, 3
            call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, nucin2%species(h)%mode(l)%n)
!     
            call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, nucin2%species(h)%mode(l)%nc)
!
#ifdef NUC_CNT1
            call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, nucin2%species(h)%mode(l)%m)
!     
            call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, nucin2%species(h)%mode(l)%mc)
#endif
	  enddo
	enddo
#endif
!
#ifdef CHEM_ENABLE
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%o3)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%co)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%ho)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%ho2)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%h2o2)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%xno)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%xno2)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%hno3)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%ch4)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%ch2o)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%ch3o2h)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%so2)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%h2so4)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, gas%dms)
#endif
!
#ifdef AQCHEM_ENABLE
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%o3)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%o3)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%civ)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%civ)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%h2o2)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%h2o2)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%nh4)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%nh4)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%xnv)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%xnv)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%siv)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%siv)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%svi)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%svi)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%ch2o)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%ch2o)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%ch3o2h)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%ch3o2h)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqc%hplus)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, aqr%hplus)
#endif
!
#ifdef SOLIDCHEM_ENABLE
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, solidi%o3)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, solidi%h2o2)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, solidi%nh4)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, solidi%xnv)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, solidi%ch2o)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, solidi%ch3o2h)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, solidi%siv)
!
        call load_data (interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, z0, z0p, solidi%svi)
#endif
!
!  Correct chemicals
!
	do k = 1, nz
	do j = jp_start,jp_end
	do i = ip_start,ip_end
#ifdef CHEM_ENABLE
	  if (gas(i,j,k)%o3    .lt.0.0) gas(i,j,k)%o3	  = 0.0
	  if (gas(i,j,k)%co    .lt.0.0) gas(i,j,k)%co	  = 0.0
	  if (gas(i,j,k)%zco2  .lt.0.0) gas(i,j,k)%zco2   = 0.0
	  if (gas(i,j,k)%ho    .lt.0.0) gas(i,j,k)%ho	  = 0.0
	  if (gas(i,j,k)%ho2   .lt.0.0) gas(i,j,k)%ho2    = 0.0
	  if (gas(i,j,k)%h2o2  .lt.0.0) gas(i,j,k)%h2o2   = 0.0
	  if (gas(i,j,k)%xno   .lt.0.0) gas(i,j,k)%xno    = 0.0
	  if (gas(i,j,k)%xno2  .lt.0.0) gas(i,j,k)%xno2   = 0.0
	  if (gas(i,j,k)%xno3  .lt.0.0) gas(i,j,k)%xno3   = 0.0
	  if (gas(i,j,k)%xn2o5 .lt.0.0) gas(i,j,k)%xn2o5  = 0.0
	  if (gas(i,j,k)%hno3  .lt.0.0) gas(i,j,k)%hno3   = 0.0
	  if (gas(i,j,k)%ch4   .lt.0.0) gas(i,j,k)%ch4    = 0.0
	  if (gas(i,j,k)%ch2o  .lt.0.0) gas(i,j,k)%ch2o   = 0.0
	  if (gas(i,j,k)%ch3o2h.lt.0.0) gas(i,j,k)%ch3o2h = 0.0
	  if (gas(i,j,k)%so2   .lt.0.0) gas(i,j,k)%so2    = 0.0
	  if (gas(i,j,k)%h2so4 .lt.0.0) gas(i,j,k)%h2so4  = 0.0
	  if (gas(i,j,k)%dms   .lt.0.0) gas(i,j,k)%dms    = 0.0
#endif

#ifdef AQCHEM_ENABLE
	  if (aqc(i,j,k)%o3    .lt.0.0) aqc(i,j,k)%o3    = 0.0
	  if (aqc(i,j,k)%civ   .lt.0.0) aqc(i,j,k)%civ   = 0.0
	  if (aqc(i,j,k)%h2o2  .lt.0.0) aqc(i,j,k)%h2o2  = 0.0
	  if (aqc(i,j,k)%xnv   .lt.0.0) aqc(i,j,k)%xnv   = 0.0
	  if (aqc(i,j,k)%siv   .lt.0.0) aqc(i,j,k)%siv   = 0.0
	  if (aqc(i,j,k)%svi   .lt.0.0) aqc(i,j,k)%svi   = 0.0
	  if (aqc(i,j,k)%ch2o  .lt.0.0) aqc(i,j,k)%ch2o  = 0.0
	  if (aqc(i,j,k)%ch3o2h.lt.0.0) aqc(i,j,k)%ch3o2h= 0.0
	  if (aqr(i,j,k)%o3    .lt.0.0) aqr(i,j,k)%o3    = 0.0
	  if (aqr(i,j,k)%civ   .lt.0.0) aqr(i,j,k)%civ   = 0.0
	  if (aqr(i,j,k)%h2o2  .lt.0.0) aqr(i,j,k)%h2o2  = 0.0
	  if (aqr(i,j,k)%xnv   .lt.0.0) aqr(i,j,k)%xnv   = 0.0
	  if (aqr(i,j,k)%siv   .lt.0.0) aqr(i,j,k)%siv   = 0.0
	  if (aqr(i,j,k)%svi   .lt.0.0) aqr(i,j,k)%svi   = 0.0
	  if (aqr(i,j,k)%ch2o  .lt.0.0) aqr(i,j,k)%ch2o  = 0.0
	  if (aqr(i,j,k)%ch3o2h.lt.0.0) aqr(i,j,k)%ch3o2h= 0.0
#endif

#ifdef SOLIDCHEM_ENABLE
	  if (solidi(i,j,k)%o3    .lt.0.0) solidi(i,j,k)%o3    = 0.0
	  if (solidi(i,j,k)%h2o2  .lt.0.0) solidi(i,j,k)%h2o2  = 0.0
	  if (solidi(i,j,k)%xnv   .lt.0.0) solidi(i,j,k)%xnv   = 0.0
	  if (solidi(i,j,k)%ch2o  .lt.0.0) solidi(i,j,k)%ch2o  = 0.0
	  if (solidi(i,j,k)%ch3o2h.lt.0.0) solidi(i,j,k)%ch3o2h= 0.0
	  if (solidi(i,j,k)%siv   .lt.0.0) solidi(i,j,k)%siv   = 0.0
	  if (solidi(i,j,k)%svi   .lt.0.0) solidi(i,j,k)%svi   = 0.0
#endif
!
	end do
	end do
	end do
!
!  Close restart files
!
	if (mypid == 0) then
	  close(106)
	  close(104)
	endif
!
	if (verbose > 2) call write_debug('Read all 3D fields')
!
	deallocate(xx, xxp, yy, yyp, z0p, fdz0p)
	deallocate(tmp)
!
!  Calculate 1D soundings
!
    	if (mypid == 0) then
          do h = 1, nhydro
            pxx(1:nz,h) = sqrt(den0(1)/den0(1:nz))
	  enddo
!
!  Time step / Dz
!
	  if(dt0.ne.dtp)then
	    write(7,*)'WARNING: TIME STEPS DO NOT MATCH'
	    write(7,*)'CONTINUING WITH DT=', dtp
	    dt0 = dtp 
	  endif
	  dt = dt0
!	
	  if(ntau.ne.ntaup)then
	    write(7,*)'WARNING: NUMBERS OF ACOUSTIC STEPS DO NOT MATCH'
	    write(7,*)'CONTINUING WITH NTAU=', ntaup
	    ntau = ntaup 
	  endif
!
!  IN/CCN
!
	  ininuc = 0
	  if ((xn_ccn0.ne.xn_ccn00) .or. (xn_in0.ne.xn_in00))then
	    write(7,*)'WARNING: REFERENCE IN/CCN NUMBER DO NOT MATCH'
	    write(7,*)'REINITIALIZING NUCLEI'
	    ininuc = 1
	  endif
!
!  Others
!
       	  call chem_m_init	   !note gas0 needed by radiation
        endif
!
!----------------------------------------------------------!
!                      Broadcasting                        ! 
!----------------------------------------------------------!
!
#if ( defined SPMD )
  if (verbose > 2) call write_debug('Broadcasting data')
  call MPI_BCAST(dt0,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(dt,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(ntau,  1, INTTYPE,0,MPI_COMM_WORLD,ierr)  

  call MPI_BCAST(kbl,    1, INTTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dz,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(zmax,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dz1,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dz1w, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dzw,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(u00,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(v00,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(qt0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(qv0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(pt0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(t0, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  if (lmicro > 1) then
    call MPI_BCAST(mse0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  else
    call MPI_BCAST(cv0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  endif
  call MPI_BCAST(ql0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(qi0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nl0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(w0, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(k0, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)      
!
  call MPI_BCAST(p0,     nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(den0,   nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(stab0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(n20,    nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(g0,     nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(avden0, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
!
  call MPI_BCAST(gas0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)    
  call MPI_BCAST(scal0, nz*nscal, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(pxx,	nz*nhydro, REALTYPE,0,MPI_COMM_WORLD,ierr)  
!
  call MPI_BCAST(ininuc,1, INTTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xn_ccn0,1,REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xn_in0,1, REALTYPE,0,MPI_COMM_WORLD,ierr)
#endif
!
!----------------------------------------------------------!
!                Initialize thermodynamics                 ! 
!----------------------------------------------------------!
!
  pressure%p = pressure2%p
!
#ifndef ANELASTIC
  pressure%dens = pressure2%dens
#else
  do k = 1, nz
    pressure%dens(:,:,k) = den0(k)
  enddo
#endif
!
  state = state2
!
  hydromtr = hydromtr2
!
  wind = wind2
!
#ifdef AERO_ENABLE
  aero3d = aero3d2
#else
  nuc = nuc2 
#endif
#ifdef NUC_CNT
  nucin = nucin2 
#endif
!
!  Initialize thermo
!
  call equation_of_state ( state, hydromtr, pressure_in=pressure%p, density_out=pressure%dens, thermo_out=thermo )
!
!----------------------------------------------------------!
!              Initialise lagrangian parcels		   !
!----------------------------------------------------------!
!
#ifdef LAGRANGE
      filename = opdir(1:len_trim(opdir)) // '/lagrange.dat'
      INQUIRE (FILE=filename, EXIST=expart)
!
      if ( expart ) then
        if (mypid == 0) then
          open( 107,file=filename(1:len_trim(filename)),         &
                form='unformatted',status='unknown' )
!
          read(107) parcel_pos%x
          read(107) parcel_pos%y         
          read(107) parcel_pos%z
          read(107) parcel_pos%exist
        endif
!
        call lagrange_proc ( .true. )
!
!  First interpolation, interpolate everything
!
        call lagrange_interp ( .false., .true. )
!
        if (mypid == 0) close (107)
!
      else if (.not.expart) then
        call lagrange_init
      endif
#endif
!
!  If new grid, write new restart file rightaway
!
    if ( interpolate ) call savedata      
!
!----------------------------------------------------------!
!                  Other initializations                   ! 
!----------------------------------------------------------!
!
!  Aerosols
!
  call init_aerosol
!
!  radiation
!
! #ifdef RAD_ENABLE
!     if (with_rad) call init_rad ( nz, p0 )
! 	if (with_rad .and. with_radmod .and. time==t_radmod) call init_rad(nz, p0)
! #endif
!
!  Diagnostics
!
    if (ldiag) call reset ( diag )
!
!----------------------------------------------------------!
!                   Output informations                    ! 
!----------------------------------------------------------!
!
      if (mypid == 0) then
        write(7,*)''
        write(7,1004)casename,time,ntime
        write(7,1005)dt0,ntau,dx,dy,dz,nx,ny,nz,lmicro
        write(7,1006)xnc0_d*1.e6,xnc0_s, xnc0_k
!
1004    format(1x/6x,'MIMICA-V5 ','  Case: ',a8,//                       &
        8x,'This is a run started with two-step-data (resume mode).',/     &
        10x,'The actual starting time is:',f15.4,' sec ',/                  &
        24x,'or, from step:',i8,/)
!
1005	format(8x,'Basic time step = ',f8.2,' second'/                 &
     	 '        Acoustic steps  = ',i3,   ' -    '/                  &
     	 '	  dx, dy, & dz    = ',3f7.1,' meter'/                  &
     	 '	  Gridpoint number= ',3i9,' in x, y, & z direction'/   &
	 '	  Microphysics level=',i2 )
!
1006   format(8x,'In calculating CN nucleation rate Nc = Cs^k (Curry & Khvorostyanov)'/   &
     	 '	  Da = ',f10.2,' micron;  sigma = ',f10.2, '  kappa = ', f10.2/)
	
	write(7,604)abs(iradx)*dt/60.0
604     format(8x,'Radiation is calculated at every ',f7.1,' minute'/)

#ifdef CHEM_ENABLE
        write(7,606)
606     format(8x,'Chemical Reaction are calculated '/)
#else
	write(7,605)
605     format(8x,'Run without Chemistry '/)
#endif
!
!  Printing initial header
!
        write(7,660)
        flush(7)
660     format(8x,60('=')/)
      endif
!
      if (verbose > 0) call write_debug('Terminate readdata')
!
    return
    end
!
!  ================================================
    Subroutine savedata
!  ================================================
!
	IMPLICIT NONE
!
	integer :: i, h, l, im
        character(len = 100):: filename
	character(len = 8)  :: car8
!
#if ( defined SPMD )
        integer, dimension(3) :: req 
        integer :: status(MPI_STATUS_SIZE)
        integer :: tag, ierr
#endif
!
        if (verbose > 0) call write_debug('Start savedata')
!
!--------------------------------------------------------!
!
!  Open restart files
!
	if (mypid == 0) then
!
	  if (all_rest) then
	    write(car8,'(i8)') int(time)
            filename = opdir(1:len_trim(opdir)) // '/restart.dat_'//adjustl(trim(car8))
	  else
            filename = opdir(1:len_trim(opdir)) // '/restart.dat'	  
	  endif
          open(106,file=filename(1:len_trim(filename)), form='unformatted',status='unknown')
!
#ifdef LAGRANGE
          filename = opdir(1:len_trim(opdir)) // '/lagrange.dat'
          open(107,file=filename(1:len_trim(filename)), form='unformatted',status='unknown')
#endif
	endif
!
!  Write initial soundings
!
	read(vers(2:3),'(i2)') ivers
	if (mypid == 0)	then
	  rewind (106)	
	  write(106)vers 
!
	  write(106)time,ntime,dt,ntau,nx,ny,nz,dx,dy,dz,nscal,nparc,xn_ccn0,xn_in0
	  write(106)u00,v00,p0,den0,stab0,n20,g0,avden0,pt0,mse0,t0,qt0,qv0,ql0,nl0,w0,k0,scal0
	  write(7,9000) 
	endif
!
! Write aerosol properties
!
#ifdef AERO_ENABLE
	if (mypid == 0) then
	  do im = 1, nmode
		write(106)aero(im)%nelem
		do i=1,4
			write(106)aero(im)%init%present(i)
			write(106)aero(im)%init%frac(i)
		enddo
		write(106)aero(im)%init%n0
		if(aero_sfc_source) write(106)aero(im)%init%n_sfc_source
		write(106)aero(im)%size%rmean,aero(im)%size%mmean,aero(im)%size%sigma
		write(106)aero(im)%mix%kappa,aero(im)%mix%rho,aero(im)%mix%mw
		write(106)aero(im)%depv0
		write(106)aero0(im)%n
		write(106)aero0(im)%m
		write(106)aero0(im)%ma
	  enddo
	endif
#endif
!
9000 format(11x,"Wrote restart data in OUTPUT/restart.dat and OUTPUT/REF0")
!
!  Writing file for parallel case
!
    	if (verbose > 2) call write_debug('Start saving full data')
!
        call save_data ( wind2%u )
!
        call save_data ( wind2%v )
!
        call save_data ( wind2%w )
!
        call save_data ( pressure2%p )
!
        call save_data ( state2%es )
!
        call save_data ( state2%qt )
!
        call save_data ( pressure2%dens )
!
	do h = 1, nhydro
          call save_data ( hydromtr2(h)%q )
!
          call save_data ( hydromtr2(h)%n )
!
	  if ( lmicro > 3 .and. h > ice ) then
	    call save_data ( hydromtr2(h)%w )
	  endif
	enddo
!
	do h = 1, nscal
	  call save_data ( state2%scal(:,:,:,h) )
	enddo
!
        if ( ivers > 0 ) then
	  if ( ivers <= 2 ) call save_data ( thermo%buoy )
!
	  if ( with_dif ) call save_data ( turbu%ksgs )
	endif
!	
#ifdef SEIFERT
	call save_data ( nuc2%ccn )
!
	call save_data ( nuc2%in )
#endif
!
        if (with_piggy) then
	  call save_data ( statep%es )
	  call save_data ( statep%qt )
	  do h = 1, nhydro
	    call save_data ( hydromtrp(h)%q )
	    call save_data ( hydromtrp(h)%n )
	    if ( lmicro > 3 .and. h > ice ) then
  	      call save_data ( hydromtrp(h)%w )
	    endif
	  enddo
        endif
!
#if (defined AERO_ENABLE)
        do h = 1, nmode
	  call save_data ( aero3d2(h)%n )
!         
	  call save_data ( aero3d2(h)%m )
!  
	  call save_data ( aero3d2(h)%ma )
        enddo
#endif
!
#ifdef NUC_CNT
	do h = 1, 3
	  do l = 1, 3
	    call save_data ( nucin2%species(h)%mode(l)%n )
! 
	    call save_data ( nucin2%species(h)%mode(l)%nc )
!
#ifdef NUC_CNT1
	    call save_data ( nucin2%species(h)%mode(l)%m )
! 
	    call save_data ( nucin2%species(h)%mode(l)%mc )
#endif
	  enddo
	enddo
#endif
!
#ifdef CHEM_ENABLE
	call save_data ( gas%o3 )
! 
	call save_data ( gas%co )
!	
	call save_data ( gas%ho )
! 
	call save_data ( gas%ho2 )
! 
	call save_data ( gas%h2o2 )
!	
	call save_data ( gas%nh3 )
! 
	call save_data ( gas%xno )
!	
	call save_data ( gas%xno2 )
! 
	call save_data ( gas%hno3 )
! 
	call save_data ( gas%ch4 )
!	
	call save_data ( gas%ch2o )
! 
	call save_data ( gas%ch3o2h )
! 
	call save_data ( gas%so2 )
! 
	call save_data ( gas%h2so4 )
! 
	call save_data ( gas%dms )
#endif		
!
#ifdef AQCHEM_ENABLE
	call save_data ( aqc%o3 )
! 
	call save_data ( aqr%o3 )
!	
	call save_data ( aqc%civ )
! 
	call save_data ( aqr%civ )
! 
	call save_data ( aqc%h2o2 )
!	
	call save_data ( aqr%h2o2 )
! 
	call save_data ( aqc%nh4 )
!	
	call save_data ( aqr%nh4 )
! 
	call save_data ( aqc%xnv )
! 
	call save_data ( aqr%xnv )
!
	call save_data ( aqc%siv )
! 
	call save_data ( aqr%siv )
!	
	call save_data ( aqc%svi )
! 
	call save_data ( aqr%svi )
! 
	call save_data ( aqc%ch2o )
!	
	call save_data ( aqr%ch2o )
! 
	call save_data ( aqc%ch3o2h )
!	
	call save_data ( aqr%ch3o2h )
! 
	call save_data ( aqc%hplus )
!	
	call save_data ( aqr%hplus )
#endif

#ifdef SOLIDCHEM_ENABLE
	call save_data ( solidi%o3 )
! 
	call save_data ( solidi%h2o2 )
!	
	call save_data ( solidi%nh4 )
! 
	call save_data ( solidi%xnv )
! 
	call save_data ( solidi%ch2o )
!	
	call save_data ( solidi%ch3o2h )
! 
	call save_data ( solidi%siv )
!	
	call save_data ( solidi%svi )
#endif
!
!  Lagrangian parcels
!
#ifdef LAGRANGE
#ifdef SPMD
!
!  Send parcel positions
!
	tag=100
        do i = 1, nparc
          if ( parcel_pos(i)%iproc == mypid .and. mypid /= 0 ) then
            call MPI_ISEND ( parcel_pos(i)%x, 1, REALTYPE, 0,   tag, MPI_COMM_WORLD, req(1), ierr )	    
            call MPI_ISEND ( parcel_pos(i)%y, 1, REALTYPE, 0, tag+1, MPI_COMM_WORLD, req(2), ierr )     
            call MPI_ISEND ( parcel_pos(i)%z, 1, REALTYPE, 0, tag+2, MPI_COMM_WORLD, req(3), ierr )     
            call MPI_WAITALL ( 3, req, MPI_STATUSES_IGNORE, ierr )
          else if ( parcel_pos(i)%iproc /= 0 .and. mypid == 0 ) then
            call MPI_RECV ( parcel_pos(i)%x, 1, REALTYPE, parcel_pos(i)%iproc,   tag, MPI_COMM_WORLD, status, ierr )      
            call MPI_RECV ( parcel_pos(i)%y, 1, REALTYPE, parcel_pos(i)%iproc, tag+1, MPI_COMM_WORLD, status, ierr )      
            call MPI_RECV ( parcel_pos(i)%z, 1, REALTYPE, parcel_pos(i)%iproc, tag+2, MPI_COMM_WORLD, status, ierr )      
          endif
	enddo
#endif
!
	if (mypid == 0) then
	  rewind(107)
	  write(107) parcel_pos%x
	  write(107) parcel_pos%y
	  write(107) parcel_pos%z
	  write(107) parcel_pos%iproc
	  write(107) parcel_pos%exist
	endif
!
        close(107)
#endif
!
!  Close and terminate
!
	if (mypid == 0) close(106)
!
!--------------------------------------------------------!
! 
    if (verbose > 0) call write_debug('Terminate savedata')
	
    return
  end
!
!  ================================================
  subroutine save_data ( data )
!  ================================================
!
  implicit none 
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: data
  real, dimension(1:nx,1:ny,1:nz) :: q1
!
#ifdef SPMD
	call collect ( data,  q1 )
	if (mypid == 0) write(106)q1
#else
       write(106)data
#endif  
!
!--------------------------------------------------------!
! 
    return
  end subroutine save_data
!
!  ================================================
  subroutine load_data ( interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, zz, zzp, data )
!  ================================================

    implicit none
!
    logical :: interpolate
    integer :: nxp, nyp, nzp
!
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: data
    real, dimension(nx) :: xx
    real, dimension(nxp) :: xxp
    real, dimension(ny) :: yy
    real, dimension(nyp) :: yyp
    real, dimension(nz) :: zz
    real, dimension(nzp) :: zzp
!
    real, allocatable, dimension(:,:,:) :: q2, q1
!
!--------------------------------------------------------!
! 
!  Interpolate
!
    if (interpolate) allocate( q1(1:nxp,1:nyp,1:nzp) )
    allocate( q2(1:nx,1:ny,1:nz) )
!  
    if (mypid == 0) then
      if (interpolate) then
        q2 = 0.
        read(106)q1
	if (.not.rep) then
    	  call interp (nxp, nyp, nzp, nx, ny, nz, xx, xxp, yy, yyp, zz, zzp, q1, q2)
	else
	  call repeated (nxp, nyp, nzp, nx, ny, nz, q1, q2)
	endif
      else
        read(106)q2
      endif
    endif
!
!  Distribute
!   
#ifdef SPMD
    call distribute ( data, q2)
#else
    data = q2
#endif
!
!  BCs
!
    call statebc ( data )
!
    deallocate (q2)
    if (interpolate) deallocate (q1)
!
!--------------------------------------------------------!
! 
    return
  end subroutine load_data
!
!  ================================================
  subroutine load_data_u ( interpolate, nxp, nyp, nzp, xx, xxp, yy, yyp, zz, zzp, data1, data2, data3 )
!  ================================================

    implicit none
!
    logical :: interpolate
    integer :: i, nxp, nyp, nzp
!
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: data1
#ifdef MODEL_3D
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: data2
#else
    real, dimension(:,:,:), allocatable :: data2
#endif
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: data3
!
    real, dimension(nxp) :: xxp
    real, dimension(nyp) :: yyp
    real, dimension(nzp) :: zzp
    real, dimension(nx) :: xx
    real, dimension(ny) :: yy
    real, dimension(nz) :: zz
!
    real, allocatable, dimension(:,:,:) :: q2, q1
!
!--------------------------------------------------------!
! 
!  Interpolate
!
    if (interpolate) allocate( q1(1:nxp,1:nyp,1:nzp) )
    allocate( q2(1:nx,1:ny,1:nz) )
!
    do i = 1, 3
      if (ivers > 2) then
#ifndef MODEL_3D
        if (i == 2) cycle
#endif
      endif
!
      if (mypid == 0) then
        if (interpolate) then
          q2 = 0.
          read(106)q1
	  if (.not.rep) then
    	    call interp (nxp, nyp, nzp, nx, ny, nz, xx, xxp, yy, yyp, zz, zzp, q1, q2)
	  else
    	    call repeated (nxp, nyp, nzp, nx, ny, nz, q1, q2)
	  endif
        else
    	  read(106)q2
        endif
      endif
!
!  Distribute
!   
#ifdef SPMD
      if (i == 1) call distribute ( data1, q2)
      if (i == 2) call distribute ( data2, q2)
      if (i == 3) call distribute ( data3, q2)
#else
      if (i == 1) data1 = q2
      if (i == 2) data2 = q2
      if (i == 3) data3 = q2
#endif
    enddo
!
!  BCs
!
    call windbc ( data1, data2, data3 )
!
    deallocate (q2)
    if (interpolate) deallocate (q1)
!
!--------------------------------------------------------!
!
    return
  end subroutine load_data_u
!
!  ================================================
  subroutine interp (nxp, nyp, nzp, nx, ny, nz, xx, xxp, yy, yyp, zz, zzp, q1, q2)
!  ================================================

    implicit none
!
    integer :: i, j, k, nxp, nyp, nzp, nx, ny, nz
    real, dimension(1) :: qinterp
!
    real, dimension(nxp,nyp,nzp) :: q1
    real, dimension(nx,ny,nz) :: q2
    real, dimension(nx) :: xx
    real, dimension(nxp) :: xxp
    real, dimension(ny) :: yy
    real, dimension(nyp) :: yyp
    real, dimension(nz) :: zz
    real, dimension(nzp) :: zzp

    real, allocatable, dimension(:,:,:) :: q1_bis, q1_ter
!
!--------------------------------------------------------!
!
  allocate ( q1_bis(nx,nyp,nzp), q1_ter(nx,ny,nzp) )
!
#ifdef MODEL_3D
!
!  Horizontal interpolations
!
    do k = 1, nzp
!
!  X interpolations
!
      do j = 1, nyp
        do i = 1, nx
          call pwl_interp_1d ( nxp, xxp, q1(:,j,k), 1, xx(i), qinterp )
	  q1_bis(i,j,k) = qinterp(1)
        enddo
      enddo
!
!  Y interpolations
!
      do j = 1, ny
        do i = 1, nx
          call pwl_interp_1d ( nyp, yyp, q1_bis(i,:,k), 1, yy(j), qinterp )
	  q1_ter(i,j,k) = qinterp(1)
        enddo
      enddo
!
    enddo
!
!  Vertical interpolations
!
    if ( nzp /= nz ) then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            call pwl_interp_1d ( nzp, zzp, q1_ter(i,j,:), 1, zz(k), q2(i,j,k) )
          enddo
        enddo
      enddo
    else
      q2 = q1_ter
    endif
#else
!
!  Horizontal interpolations
!
    do k = 1, nzp
      do i = 1, nx
        call pwl_interp_1d ( nxp, xxp, q1(:,1,k), 1, xx(i), q1_ter(i,1,k) )
      enddo
    enddo
!
!  Vertical interpolations
!
    if ( nzp /= nz ) then
      do k = 1, nz
        do i = 1, nx
          call pwl_interp_1d ( nzp, zzp, q1_ter(i,1,:), 1, zz(k), q2(i,1,k) )
        enddo
      enddo
    else
      q2 = q1_ter
    endif
#endif
!
  deallocate (q1_bis,q1_ter)
!  
!--------------------------------------------------------!
!
  return
  end subroutine interp
!
!  ================================================
  subroutine interp_1d (interpolate, nz, nzp, zz, zzp, q, qp)
!  ================================================

    implicit none
!
    logical :: interpolate
    integer :: k, nzp, nz
!
    real, dimension(nzp) :: qp
    real, dimension(nz) :: q
    real, dimension(nzp) :: zzp
    real, dimension(nz) :: zz
!
!--------------------------------------------------------!
!
!  No interpolation
!
  if ( .not.interpolate ) then
    q = qp
!
!  Interpolate
!
  else
    do k = 1, nz
      call pwl_interp_1d ( nzp, zzp, qp, 1, zz(k), q(k) )
    enddo  
  endif
!  
  return
  end subroutine interp_1d
!
!  ================================================
  subroutine repeated (nxp, nyp, nzp, nx, ny, nz, q1, q2)
!  ================================================

    implicit none
!
    integer :: i, j, nxp, nyp, nzp, nx, ny, nz
    integer :: l, m, fx, fy, is, ie, js, je
!
    real, dimension(nxp,nyp,nzp) :: q1
    real, dimension(nx,ny,nz) :: q2
!
!--------------------------------------------------------!
!
  if (nx /= nxp) then
    fx = (nx-5)/(nxp-5)
    do l = 1, fx
      if (l == 1) then
        is = 1
      else
        is = 4
      endif
      if (l == fx) then
        ie = nxp
      else
        ie = nxp - 2
      endif
!
#ifdef MODEL_3D
      if (ny /= nyp) then
        fy = (ny-5)/(nyp-5)
        do m = 1, fy
          if (m == 1) then
            js = 1
          else
            js = 4
          endif
          if (m == fy) then
            je = nyp
          else
            je = nyp - 2
          endif
!
          do j = js, je
            do i = is, ie
              q2(i+(l-1)*(nxp-5),j+(m-1)*(nyp-5),:) = q1(i,j,:)
            enddo
          enddo
        enddo
      else
        do i = is, ie
          q2(i+(l-1)*(nxp-5),:,:) = q1(i,:,:)
        enddo
      endif
#else
      do i = is, ie
        q2(i+(l-1)*(nxp-5),1,:) = q1(i,1,:)
      enddo
#endif
!
    enddo
  endif
!  
  return
  end subroutine repeated

end module restart
