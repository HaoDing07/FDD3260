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
!  diagnostics.f90                   
!
!  Purpose:
!      A package for calculating diagnostic variables before output                 
!
!  Author
!      Julien Savre
!      MIM, Ludwig Maximilian Universitat, Munich
!
! ================================================================
!

  module diagnostics

!
!-----------------------------------------------------------!
!
  USE shared_all
  USE allocation
  USE netcdfmod
  USE averages
  USE funcpack
  USE vof
  USE gradients
  USE advection, only : advection_ac
  USE thermodynamics
  USE pressure_solver
  USE boundary_conditions
  USE micropack
  USE sbpack
  USE kessler
  USE getcape
!
#ifdef SPMD
  USE mpi
  USE mpicomm
#endif
!
IMPLICIT NONE

  private

  real, dimension(:,:), allocatable :: zi0

  public :: thermo_diagnostics, special_diagnostics, time_series
  
  contains
!  
! ===============================================
  subroutine special_diagnostics
! ===============================================
!
! --- Calculate special diagnostics
!
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: pt
!
  if (verbose > 0) call write_debug('Starting special_diagnostics')
!
!-----------------------------------------------------------!
!
  if (spec_diag) then
!
!  General Diagnostics
!
  call get_diag_diag ( diag, thermo%pt, state%qt, hydromtr, wind )
!
!-----------------------------------------------------------!
!	      Diagnostic pressure and buoyancy	            !
!-----------------------------------------------------------!
!
    if (out_beff .or. out_dp) call pressure_diagnostics ( thermo%buoy, winds )
!
!-----------------------------------------------------------!
!		 SGS turbulence diagnostics	            !
!-----------------------------------------------------------!
!
    if (out_tke .or. out_buoy) call tke_diagnostics
!
!-----------------------------------------------------------!
!	   	      Cloud diagnostics	                    !
!-----------------------------------------------------------!
!
    if (lmicro > 0) call cloud_diagnostics 
!
!-----------------------------------------------------------!
!	            Dynamical diagnostics	            !
!-----------------------------------------------------------!
!
    if (out_div) then
      call caldiv ( wind, diagnos%div, horizontal=.true. )
      call calconv ( wind, state%qt, diagnos%qconv, horizontal=.true. )
    endif
!   
#ifdef MODEL_3D
    if (out_vort) call cal_rot ( wind, diagnos%vortx, diagnos%vorty, diagnos%vortz )
#endif
!  
  endif
!
!-----------------------------------------------------------!
!		   Entrainment/Detrainment	            !
!-----------------------------------------------------------!
!
!  if (spec_diag .and. out_ent) call level_set_entrainment
!
  if (spec_diag .and. out_ent) call romps_entrainment
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating special_diagnostics')
!
return
end subroutine special_diagnostics!
!
! ===============================================
  subroutine thermo_diagnostics
! ===============================================
!
! --- Calculate thermo diagnostics
!
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: tem, pt
!
  if (verbose > 1) call write_debug('Starting thermo_diagnostics')
!
!----------------------------------------------------------!
!	             Thermodynamics                        !
!----------------------------------------------------------!
!
  call equation_of_state ( state, hydromtr, pressure_in=pressure%p, thermo_out=thermo, conservative=.false. )
!
  call get_qv ( state%qt, hydromtr, thermo%qv )
! 
  if (out_ptv) then
    call get_ptv ( pressure, state, hydromtr, thermo%ptv )
  endif
! 
#ifdef ISENTROPIC
  if (out_mse) then
    call get_mse ( pressure, state, hydromtr, thermo%mse )
  endif
#endif
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating thermo_diagnostics')
!
return
end subroutine thermo_diagnostics

!
! ==================================================
  subroutine tke_diagnostics
! ==================================================
!
  integer :: i, j, k
  real, dimension(nz) :: uav2, vav2, wav2, bav
!
  if (verbose > 1) call write_debug('Starting tke_diagnostics')
!
!----------------------------------------------------------!
!
!  Horizontal averages
!
  call horav (wind%u, uav2)
#ifdef MODEL_3D
  call horav (wind%v, vav2)
#endif
  call horav (wind%w, wav2)
!
!  Turbulent kinetic energy
!
  if (out_tke) then
    turbu_diag%kres = 0.
    do k = 1, nz-1
      do j = jt_start, jt_end
        do i = it_start, it_end
          turbu_diag%kres(i,j,k) = 0.5*((0.5*(wind%u(i,j,k)+wind%u(i+1,j,k)) - uav2(k))**2.      &
#ifdef MODEL_3D
			         + (0.5*(wind%v(i,j,k)+wind%v(i,j+1,k)) - vav2(k))**2.  	 &
#endif
			         + (0.5*(wind%w(i,j,k)+wind%w(i,j,k+1)))**2.)
        enddo
      enddo
    enddo
!
!  W moments
!
    do k = 1, nz
      turbu_diag%wvar(:,:,k) = wind%w(:,:,k) - wav2(k)
      turbu_diag%wske(:,:,k) = wind%w(:,:,k) - wav2(k)
    enddo
!
    do k = 1, nz-1
      turbu_diag%wvar(:,:,k) = (0.5*(turbu_diag%wvar(:,:,k+1)+turbu_diag%wvar(:,:,k)))**2.
      turbu_diag%wske(:,:,k) = (0.5*(turbu_diag%wske(:,:,k+1)+turbu_diag%wske(:,:,k)))**3.
    enddo
  endif
!
!  Buoyancy fluxes
!
  if (out_flut.and.out_buoy) then
    call horav (thermo%buoy, bav)
    do k = 1, nz
      turbu_diag%bfres(:,:,k) = (wind%w(:,:,k) - wav2(k))*(thermo%buoy(:,:,k) - bav(k))
    enddo
  endif
!
  if (out_fsgs.and.out_buoy) then
    turbu_diag%bfsgs = -pressure%dens * turbu_diag%N2 * turbu%fkv / abs(pran)
  endif
!
  if (verbose > 1) call write_debug('Terminating tke_diagnostics')
!
return
end subroutine tke_diagnostics

!----------------------------------------------------------!
! Calculating vp of droplet in Stokes regime (Bergot,2016)!
!----------------------------------------------------------!
!
function cal_vdrop (h,q,n,pxx)

  integer :: h
  real :: q, n, vdrop_mean, lambda_drop, lambda_c, pxx, cal_vdrop
  real :: n0, q0, eta_vis, d2_expect

  cal_vdrop = 0.
  eta_vis = 1.827e-5

 
  if (h == drop .AND. with_dropsed) then 
  

    ! n0 = max(1.e6,n)
    n0 = n
    ! n0 = 300.e6
    q0 = q
    ! q0 = max(q,1.e-12)
    if (q0/=0 .and. n0/=0) then
    ! lambda_c = (pi/6*1000*gamma(nu_drop+3/alpha_drop)/gamma(nu_drop)*n0/1.29/q)**(1/3)
      lambda_c = (4.0/3.0*pi*1000.0/1.29*n0/q0*gamma(nu_drop+3.0/alpha_drop)/gamma(nu_drop))**(1.0/3.0)
    
      ! lambda_drop = min(max(lambda_c, 1.0),10.0)
      lambda_drop = lambda_c

      ! d2_expect = lambda_drop**((alpha_drop-1)/alpha_drop)/alpha_drop*(alpha_drop*nu_drop+1)/alpha_drop*gamma((alpha_drop*nu_drop+1)/alpha_drop)/gamma(nu_drop)

      !!!define alpha_drop and nu_drop in cm.nml and mimica.f90
      ! vdrop_mean = 1.0/18.0*9.8*d2_expect*(1000-1.29)/eta_vis
      vdrop_mean = 2.0/9.0*9.8*(1000.0-1.29)/eta_vis/(lambda_drop**2.0)*gamma(nu_drop+5.0/alpha_drop)/gamma(nu_drop+3.0/alpha_drop)

      ! cal_vdrop = max(vdrop_mean, 1.e-4)
      cal_vdrop = vdrop_mean
      
    end if
  
  end if

  return

  end function cal_vdrop
!
! ==================================================
  subroutine cloud_diagnostics
! ==================================================
!
  integer :: i, j, k, h, npts
  real    :: zz, bdiag
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: tmp
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ptv, suml, sumi, buoy
  real, dimension(1:nz) :: ptvav, ptav, qvav, grad_pt

!
  real, parameter :: eps = 1.e-13


!
  if (verbose > 0) call write_debug('Starting cloud_diagnostics')
!
!----------------------------------------------------------!
!                   Precip velocities                      !
!----------------------------------------------------------!
!
  do h = 2, nhydro
    hydromtr(h)%vp = 0.0 
    do k = 1, nz
      do j=jt_start,jt_end
        do i=it_start,it_end
#ifndef SEIFERT
          hydromtr(h)%vp(i,j,k) = cal_vp ( h, hydromtr(h)%q(i,j,k), hydromtr(h)%n(i,j,k), pxx(k,h) )
#else
          hydromtr(h)%vp(i,j,k) = calc_vel_k ( pressure%dens(i,j,k), thermo%T(i,j,k), hydromtr(h)%q(i,j,k) )
#endif
        enddo
      enddo
    enddo
  enddo

  do h = 1,1
    hydromtr(h)%vp = 0.0 
    do k = 1, nz
      do j=jt_start,jt_end
        do i=it_start,it_end
! #ifndef SEIFERT
          hydromtr(h)%vp(i,j,k) = cal_vdrop ( h, hydromtr(h)%q(i,j,k), hydromtr(h)%n(i,j,k), pxx(k,h) )
! #else
          ! hydromtr(h)%vp(i,j,k) = cal_vdrop ( h, hydromtr(h)%q(i,j,k), hydromtr(h)%n(i,j,k), pxx(k,h) )
! #endif
        enddo
      enddo
    enddo
  enddo
! !
!----------------------------------------------------------!
!	           Hydrometeor mean size                   !
!----------------------------------------------------------!
! 
  if (out_dc .or. out_dr .or. out_di .or. out_dg .or. out_ds) then
    call hydro_sizes ( hydromtr )
  endif
!
!----------------------------------------------------------!
!	       Equivalent Radar reflectivity               !
!----------------------------------------------------------!
! 
  if (out_z) then
    call reflectivity ( hydromtr, diagnos%Z )
  endif
!
!----------------------------------------------------------!
!	      Convective cloud diagnostics                 !
!----------------------------------------------------------!
!
  suml = hydromtr(drop)%q !+ hydromtr(rain)%q
  sumi = 0.
  if (lmicro > 1) sumi = hydromtr(ice)%q
  if (lmicro > 2) sumi = sumi + hydromtr(grau)%q + hydromtr(snow)%q
!  
  if (out_mf) then
    diagnos%mf = 0.
    diagnos%cc = 0.
    diagnos%mfm = 0.
    diagnos%ccm = 0.
    where ( suml+sumi >= qthres )
      diagnos%mf = pressure%dens*wind%w*dx*dy
      diagnos%cc = 1.
    end where
    where ( suml+sumi >= qthres .and. wind%w > wthres )
      diagnos%mfm = pressure%dens*wind%w*dx*dy
      diagnos%ccm = 1.
    end where
  endif
!
!----------------------------------------------------------!
!	        Find LWP, WVP and cloud top		   !
!----------------------------------------------------------!
!
!  LWP
!  
  if (out_lwp) then
    call vertint ( suml, column%lwp, weighted=.true. )
    if (lmicro > 1) call vertint ( sumi, column%iwp, weighted=.true. )
  endif
!
!  CWP
!  
  if (out_cwp) then
    call vertint ( hydromtr(drop)%q, column%cwp, weighted=.true. )  
    call vertint ( hydromtr(rain)%q, column%rwp, weighted=.true. )    
  endif
!  
!  WVP
!
  if (out_qv.and.out_wvp) then
    call vertint (thermo%qv, column%wvp, weighted=.true.)
    call vertint (thermo%qv, column%wvpb, zlim=zbl, weighted=.true.)
    call vertint (thermo%qv, column%wvpt, zbas=zbl, weighted=.true.)
  endif
!  
!  CMSE
!
  if (out_mse.and.out_cmse) then
    call vertint (thermo%mse, column%cmse, weighted=.true.)
    call vertint (thermo%mse, column%cmseb, zlim=zbl, weighted=.true.)
    call vertint (thermo%mse, column%cmset, zbas=zbl, weighted=.true.)
  endif
!  
!  CMFL
!
  if (out_cmfl) then
    call vertint (wind%w, column%cmfl, weighted=.true.)
  endif
!
!  Cloud top/base
!
  if (out_ctop) then
    do k = 2, nz-1
      where ( suml(:,:,k)+sumi(:,:,k) > qthres .and. wind%w(:,:,k) > wthres .and. column%cbas < dz/2. ) column%cbas = z0(k)
    enddo
!
    do k = nz, 2, -1
      where ( suml(:,:,k)+sumi(:,:,k) > qthres .and. wind%w(:,:,k) > wthres .and. column%ctop < dz/2. ) column%ctop = z0(k)
    enddo
  endif
!
!  BL/inversion top
!
  if (out_zinv) then
    do j = jt_start, jt_end
      do i = it_start, it_end
        call grad_1d ( state%es(i,j,:), grad_pt )
        column%zinv(i,j) = z0( SUM(maxloc(grad_pt,MASK=z0.lt.3000.)) )
      enddo
    enddo
!
    call get_ptv ( pressure, state, hydromtr, ptv )
    call horav ( ptv, ptvav )
    do j = jt_start, jt_end
      do i = it_start, it_end
        column%bltop(i,j) = z0( SUM(minloc(wind%w(i,j,:)*(ptv(i,j,:)-ptvav(:)),MASK=z0(:).lt.3000.)) )
      enddo
    enddo
  endif
!
!----------------------------------------------------------!
!              Precipitation fluxes mm/day	           !
!----------------------------------------------------------!
!
  if ( out_rrate ) then
    surf%precip = pressure%dens(:,:,1)*hydromtr(rain)%q(:,:,1)*hydromtr(rain)%vp(:,:,1)*3600.
    if (lmicro > 2) surf%precip = surf%precip + pressure%dens(:,:,1)*hydromtr(grau)%q(:,:,1)*hydromtr(grau)%vp(:,:,1)*3600. 
    surf%cumul = surf%cumul + surf%precip
  endif
!
!----------------------------------------------------------!
!		     Find CAPE and CIN		           !
!----------------------------------------------------------!
!  
  if (out_cape) then
!    call horav ( state%es, ptav )
!    call horav ( state%qt, qvav )
!
    do j = jt_start, jt_end
      do i = it_start, it_end
        !call get_cape ( nz, p0, ptav, qvav, column%cape(i,j), column%cin(i,j), th0=state%es(i,j,1), qv0=state%qt(i,j,1), zb=column%lfc(i,j), zc=column%lcl(i,j) )
        call get_cape ( nz, p0, state%es(i,j,:), state%qt(i,j,:), column%cape(i,j), column%cin(i,j), zb=column%lfc(i,j), zc=column%lcl(i,j), src=3 )
      enddo
    enddo
  endif
!
!----------------------------------------------------------!
!	     Find basic cold pool diagnostics              !
!----------------------------------------------------------!
!  
  if (out_cp) then
    call get_ptv ( pressure, state, hydromtr, ptv )
!
    call horav ( ptv, ptvav )
!
    buoy = 0.
    do k = 1, nz
      buoy(:,:,k) = g*(ptv(:,:,k) - ptvav(k)) / ptvav(k)
    enddo
!
    do j = jt_start, jt_end
      do i = it_start, it_end
	zz = 0.
        k = 1
	do while ( zz <= zcp )
          column%cpint(i,j) = column%cpint(i,j) - 2.*buoy(i,j,k)*dz*fdz0(k)
	  zz = zz + dz*fdz0(k)
          k = k + 1
	enddo
      enddo
    enddo
  endif
!
!----------------------------------------------------------!
!	       Vertically integrated scalar                !
!----------------------------------------------------------!
! 
  if (out_ints) then 
    call vertint (state%scal(:,:,:,1)*1.e-6, column%intsca, weighted=.true.) 
  endif
!
!----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating cloud_diagnostics')
!
return
end subroutine cloud_diagnostics
!
! ==================================================
  subroutine entrainment_diagnostics ( dens )
! ==================================================
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dens
  real :: thres=1.e-9
!
!----------------------------------------------------------!
!
!  Standard plume entrainment (1D only)
! 
!    if (out_plume) then
!      call plume_entrainment ( dens, wh, p, thermo%ptv, thermo%buoy, qt, hydromtr(drop)%q+hydromtr(rain)%q, 'wcore' )
!
!      call plume_entrainment ( dens, wh, p, thermo%ptv, thermo%buoy, qt, hydromtr(drop)%q+hydromtr(rain)%q, 'ccore' )
!    endif
!
!    if (pro_ent) call filtered_profiles ( entrain%went, entrain%wdet, hydromtr(drop)%q+hydromtr(rain)%q, sca1, thres, 1., 'wcore' )
!
return
end subroutine entrainment_diagnostics
!      
! === time_series:
!      calculates time series of global quantities
! 
!      =================================================   
  subroutine time_series (ts_tab)        
!      =================================================   
                         
      USE gridno
      USE shared_data

      IMPLICIT NONE
!
      real, dimension(ntsout) :: ts_tab
!
! ================================================================
!
  if (verbose > 0) call write_debug('Starting time_series')
!
  if (trim(ts_out) == 'stratus') then
    call stratus_timeser (ts_tab)
  else if (trim(ts_out) == 'cumulus') then
    call cumulus_timeser (ts_tab)
  else if (trim(ts_out) == 'deep') then
    call deep_timeser (ts_tab)
  else if (trim(ts_out) == 'conservation') then
    call conservation_timeser (ts_tab)
  endif
!
! ================================================================
!    
  if (verbose > 0) call write_debug('Terminating time_series')
!
  RETURN       
  END subroutine time_series    
!      
! === stratus_timeser:
!      Calculates output variables for time series specific to stratiform clouds
! 
!      =================================================                            
  subroutine stratus_timeser (ts_tab)        
!      =================================================                            
!
      integer        :: i, j, k, h, nn, nt, k1
!
      real, dimension(ntsout) :: ts_tab
      real, dimension(nz)    :: xav, xav2, x2all
!
      real           :: tmp, tmp2, z
      real           :: x1, x1all, x1all2
      real           :: cbh
 
#if ( defined SPMD )      
      integer :: ierr, status(MPI_STATUS_SIZE)
#endif
!
      real, dimension (:,:), allocatable :: zi, lwp, prec, zb, slice
      real, save :: time0
!
! ================================================================
!
  ts_tab = 0.
!
!  Allocate
!
  call alloc ( zi )
  call alloc ( zb )
  call alloc ( lwp )
  call alloc ( prec )
  call alloc ( slice )
!
  if ( .not.allocated(zi0) ) then
    allocate( zi0(ip_start:ip_end,jp_start:jp_end) )
    zi0 = 0.
    time0 = 0.
  endif
!
!  Mean LWP:
!
  if (lmicro > 0) then
    lwp = 0.0

    do i = it_start, it_end
      do j = jt_start, jt_end
        do k = 1, nz-1
          tmp = hydromtr(drop)%q(i,j,k) + hydromtr(rain)%q(i,j,k)
          lwp(i,j) = lwp(i,j) + den0(k)*max(tmp,0.0)*dz*fdz0(k)
        enddo
      enddo
    enddo
    call horav (lwp,x1all)

  else
    x1all = 0.0
  endif
  ts_tab(1) = x1all*1000.
!
!  Mean IWP:
!
  if (lmicro > 1) then
    slice = 0.0
    do i = it_start, it_end
      do j = jt_start, jt_end
        do k = 1, nz-1
	  tmp = 0.
          if (lmicro > 1) tmp = tmp + hydromtr(ice)%q(i,j,k) 
	  if (lmicro > 2) tmp = tmp + hydromtr(grau)%q(i,j,k) + hydromtr(snow)%q(i,j,k)
	  if (lmicro > 3) tmp = tmp + hydromtr(hail)%q(i,j,k)
	  slice(i,j) = slice(i,j) + den0(k)*max(tmp,0.0)*dz*fdz0(k)
        enddo
      enddo
    enddo
    call horav (slice,x1all)
  else
    x1all = 0.0
  endif
!
  ts_tab(2) = x1all*1000.
!
!  Mean cloud top height:
!
  if (lmicro > 0) then
    zi  = 0.0
    do i = it_start, it_end
      do j = jt_start, jt_end
        tmp2 = 0.0
        z = zmax
        do k = nz, 2, -1
          tmp = hydromtr2(drop)%q(i,j,k)
          if (tmp >= qthres .and. wind%w(i,j,k) >= wthres .and. tmp2 < qthres) then
            zi(i,j) = z
           exit
           else
             tmp2 = tmp
             z = z - dz*0.5*(fdz0(k)+fdz0(k-1))
            endif
        enddo
      enddo
    enddo
    if (maxval(zi0) == 0.) zi0 = zi
!
    call horav (zi,x1all)
    ts_tab(3) = x1all
  endif
!
!  Mean cloud base height:
!
  if (lmicro > 0) then
    zb = 0.
    do i = it_start, it_end
      do j = jt_start, jt_end
        tmp2 = 0.0
        z = 0.
        do k = 1, nz-1
          tmp = hydromtr2(drop)%q(i,j,k)
          if (tmp >= qthres .and. wind%w(i,j,k) > wthres) then
            zb(i,j) = z
            exit
          else
            z = z + dz*0.5*(fdz0(k)+fdz0(k+1))
          endif
        enddo
      enddo
    enddo
!
    call horav (zb,x1all)
    ts_tab(4) = x1all
    cbh = ts_tab(4)
  endif
!
!  Cloud fraction
!
  nn = 0      
  nt = 0
  x1all = 0.0
!
  if (lmicro > 0) then
  do i = it_start, it_end
    do j = jt_start, jt_end
      nt = nt + 1
      k = 1
      tmp = 0.
      do while (k < nz .and. tmp < qthres)
        tmp = 0.
	if (lmicro > 0) tmp = tmp + hydromtr(drop)%q(i,j,k) + hydromtr(rain)%q(i,j,k) 
	if (lmicro > 1) tmp = tmp + hydromtr(ice)%q(i,j,k) 
	if (lmicro > 2) tmp = tmp + hydromtr(grau)%q(i,j,k) + hydromtr(snow)%q(i,j,k)
	if (lmicro > 3) tmp = tmp + hydromtr(hail)%q(i,j,k)
        if (tmp >= qthres) nn = nn + 1
      k = k + 1
      enddo
    enddo
  enddo
  x1 = real(nn)/real(nt)
!
#ifdef SPMD
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  x1all = x1all/real(nproc)
#else
  x1all = x1
#endif 
!
  else
    x1all=0.0
  endif
!
  ts_tab(5) = x1all
!
!  Mean surface precipitation in mm/day
!
  if (lmicro > 0) then
    prec = 0.
    do h = 1, nhydro
      prec = prec + den0(1)*hydromtr(h)%q(:,:,1)*hydromtr(h)%vp(:,:,1)*86400.
    enddo
    call horav (prec, x1all)
  else
    x1all = 0.0
  endif
!
  ts_tab(6) = x1all     	
!
!  Max vertical wind 
!
  call max_all (wind%w, x1)
!
#ifdef SPMD
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
  x1all = x1
#endif 
!
  ts_tab(7) = x1all
!
!  Max vertical wind variance
!
  call max_all (wind%w*wind%w, x1)
!
#ifdef SPMD
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
  x1all = x1
#endif 
!
  ts_tab(8) = x1all
!
!  Deardorff convective velocity
!
!  BIR: Buoyancy Integral Ratio (decoupling index)
!
  if (out_flut) then
    xav = 0.
    xav2 = 0.
    if (out_buoy) call horav (turbu_diag%bfres, xav)
    if (out_fsgs) call horav (turbu_diag%bfsgs, xav2)
!
    k1 = 2
    xav = xav + xav2
    do k = 2, nzi
      if (xav(k) < 0. .and. xav(k+1) >= 0.) then
        k1 = k
        exit
      endif
    enddo
!
    z = 0.5*dz*fdz0(1)
    k = 2
    x1all = 0.
    x1all2 = 0.
    do while (z <= zbl .and. k < nz)
      if (k <= k1) then
        x1all2 = x1all2 + xav(k)*dz*0.5*(fdz0(k)+fdz0(k-1))
      else
        x1all = x1all + xav(k)*dz*0.5*(fdz0(k)+fdz0(k-1))
      endif
      z = z + dz*0.5*(fdz0(k)+fdz0(k-1))
      k = k + 1    
    enddo
!
    ts_tab(9) = -min(x1all2,-1.e-20)/max(x1all,1.e-20)
  endif
!
!  Variance of LWP
!
  if (lmicro > 0) then
    call horav (lwp, x1)
    lwp = (lwp - x1)**2.
    call horav (lwp, x1all)
  else
    x1all = 0.0
  endif
  ts_tab(10) = sqrt(x1all) 
!
!  Variance of Zi
!
  call horav (zi, x1)
  slice = (zi - x1)**2.
  call horav (slice, x1all)
!
  ts_tab(11) = sqrt(x1all) 
!
!  Variance of Zb
!
  call horav (zb, x1)
  slice = (zb - x1)**2.
  call horav (slice, x1all)
!
  ts_tab(12) = sqrt(x1all) 
!
!  Averaged Entrainment rate
!
  slice = 0.
  slice = (zi - zi0)/real(max(time,1.) - time0) + Ddiv*zi
  call horav (slice, x1)
  time0 = time
  zi0 = zi
!
  ts_tab(13) = 100.*x1   !cm/s
!
!  Max deposition rate
!
  if (out_diagi) then
    call max_all (diag(4)%qi, x1)
    ts_tab(14) = x1
  endif
!
#ifdef RAD_ENABLE
  ts_tab(15) = ssr
  ts_tab(16) = olr
  ts_tab(17) = olr_cs
  ts_tab(18) = toa
  ts_tab(19) = toa_cs
#endif
!
!  Mean NC 
!
  if (lmicro > 0) then
    call av_all (hydromtr(drop)%n, x1, xq1=hydromtr(drop)%q, xc1=qthres, filter1='larger')
    ts_tab(20) = x1
  endif
!
!  Allocate
!
  call dealloc ( zi )
  call dealloc ( zb )
  call dealloc ( lwp )
  call dealloc ( prec )
  call dealloc ( slice )
!
! ================================================================
!    
RETURN       
END subroutine stratus_timeser     
!      
! === cumulus_timeser:
!      Calculates output variables for time series specific to cumuliform clouds
! 
!      =================================================                            
  subroutine cumulus_timeser (ts_tab)      
!      =================================================                            
!
      integer :: i, j, k, h, nzt, nzb, nt, ntall, nzball, nztall, maxl(1)
!
      real, dimension(ntsout) :: ts_tab
      real, dimension(nz) :: a, xav, ptav, qtav
      real :: z, zm, val, val2, x1, x2, x3, x4, x1all, x2all, x3all, x4all
      real :: xav1, xav2, xav3, valall, denall
!      
#if ( defined SPMD )      
      integer :: ierr, status(MPI_STATUS_SIZE)
#endif
!
      real, dimension(3) :: alt
      real, allocatable, dimension(:,:,:) :: tem, rh, cmax
      real, allocatable, dimension (:,:) :: zb, mfb, mfbm, zi, prec, peff, cwp, pw, crh, cc, ccm
      real, parameter :: eps = 1.e-15
!
! ================================================================
!
  ts_tab = 0.
!
!  Allocate
!
  call alloc ( tem )
  call alloc ( rh )
  call alloc ( cmax )
  call alloc ( zi )
  call alloc ( zb )
  call alloc ( cc )
  call alloc ( mfb )
  call alloc ( ccm )
  call alloc ( mfbm )
  call alloc ( prec )
  call alloc ( peff )
  call alloc ( cwp )
  call alloc ( pw )
  call alloc ( crh )
!
!  Max vertical wind 
!
  call max_all (wind%w, x1)
!
#ifdef SPMD
  x1all = 0.
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
  x1all = x1
#endif 
  ts_tab(1) = x1all
!
!  Mean cloud top height:
!
  if (lmicro > 0) then 
!
  zi  = 0.0
  nzt = 0
  do i = it_start, it_end
    do j = jt_start, jt_end
      z = zmax
      do k = nz-1, 2, -1
	if ( hydromtr(drop)%q(i,j,k) >= qthres .and. hydromtr(drop)%q(i,j,k+1) < qthres) then
	  zi(i,j) = z
	  nzt=nzt+1
	  exit
	else
	  z = z - dz*0.5*(fdz0(k) + fdz0(k-1))
	endif
      enddo
    enddo
  enddo
  x1 = sum(zi)
!
#ifdef SPMD
  nztall = 0
  x1all = 0.
  call MPI_ALLReduce (nzt, nztall, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  nztall = nzt
  x1all = x1
#endif 
!
  ts_tab(2) = x1all/real(max(nztall,1))
!
!  Mean cloud base height:
!
  zb = 0.0
  nzb = 0
  do i = it_start, it_end
    do j = jt_start, jt_end
      k = 2
      z = dz*(fdz0(2) + 0.5*fdz0(1))
      do while (z < 2500.)
	if ( hydromtr(drop)%q(i,j,k) >= qthres .and. hydromtr(drop)%q(i,j,k-1) < qthres ) then
	  zb(i,j) = z
	  nzb=nzb+1
	  exit
	endif
        k = k + 1
	z = z + dz*0.5*(fdz0(k) + fdz0(k+1))
      enddo
    enddo
  enddo
  x1 = sum(zb)
!
#ifdef SPMD
  nzball = 0
  x1all = 0.
  call MPI_ALLReduce (nzb, nzball, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  nzball = nzb
  x1all = x1
#endif 
!
  ts_tab(3) = x1all/real(max(nzball,1))
!
!  Max cloud fraction height:
!
  cmax = 0.
  where ( hydromtr(drop)%q >= qthres ) cmax = 1.
  call horav (cmax, xav)
  zm = z0( SUM(maxloc(xav,MASK=z0.lt.2500.)) )
!
  ts_tab(4) = zm 
!
!  Cloud fraction and mass flux
!
  alt = (/zm,4000.,8000./)
  do h = 1, 3
    cc = 0.
    ccm = 0.
    mfb = 0.
    mfbm = 0.
    do k = 1, nz
      if (z0(k) >= alt(h)) then
        where ( hydromtr(drop)%q(:,:,k) >= qthres )
          cc(:,:) = 1
          mfb(:,:) = den0(k)*wind%w(:,:,k)*dx*dy
        end where
        where ( hydromtr(drop)%q(:,:,k) >= qthres .and. wind%w(:,:,k) >= wthres )
          ccm(:,:) = 1
          mfbm(:,:) = den0(k)*wind%w(:,:,k)*dx*dy
        end where
        exit
      endif
    enddo
    x1 = sum(cc(it_start:it_end,jt_start:jt_end))
    x2 = sum(mfb(it_start:it_end,jt_start:jt_end))
    x3 = sum(ccm(it_start:it_end,jt_start:jt_end))
    x4 = sum(mfbm(it_start:it_end,jt_start:jt_end))
!
#ifdef SPMD
    call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLReduce (x2, x2all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLReduce (x3, x3all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLReduce (x4, x4all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    x1all = x1
    x2all = x2
    x3all = x3
    x4all = x4
#endif 
!
    ts_tab(5+(h-1)*4) = x1all / (nx*ny)
    ts_tab(6+(h-1)*4) = x2all
    ts_tab(7+(h-1)*4) = x3all / (nx*ny)
    ts_tab(8+(h-1)*4) = x4all
  enddo
!
!  Mean CWP, PW and CRH:
!
  call get_temperature ( pressure, state, hydromtr, tem )
!
  cwp = 1.e-13
  pw = 1.e-13
  crh = 1.e-13
  do i = it_start, it_end
    do j = jt_start, jt_end
      do k = 1, nz-1
        val = hydromtr(drop)%q(i,j,k) + hydromtr(rain)%q(i,j,k)
	if (lmicro > 1) val = val + hydromtr(ice)%q(i,j,k)
	if (lmicro > 2) val = val + hydromtr(grau)%q(i,j,k) + hydromtr(snow)%q(i,j,k)
	cwp(i,j) = cwp(i,j) + den0(k)*val*dz*fdz0(k)
	pw(i,j) = pw(i,j) + den0(k)*state%qt(i,j,k)*dz*fdz0(k)
	crh(i,j) = crh(i,j) + den0(k)*cal_qsw(tem(i,j,k),p0(k))*dz*fdz0(k)
      enddo
    enddo
  enddo
  call horav (cwp, xav1)
  call horav (pw, xav2)
  call horav (pw/crh, xav3)
!
  ts_tab(17) = xav1
  ts_tab(18) = xav2
  ts_tab(19) = 100.*xav3
!
!  Variance of PW
!
  pw = (pw - xav2)**2.
  if (maxval(pw) > 0.1) then
    call horav (pw, xav2)
    ts_tab(20) = xav2
  endif
!
!  Mean surface precipitation in mm/day
!
  prec(it_start:it_end,jt_start:jt_end) = den0(1)*hydromtr(rain)%q(it_start:it_end,jt_start:jt_end,1)*hydromtr(rain)%vp(it_start:it_end,jt_start:jt_end,1)*86400.
!
  call horav (prec, x1all)
  call max_all (prec/24., x2all)
  call horav (surf%cond, x3all)
  call horav (surf%evap, x4all)
!
  ts_tab(21) = x1all     
  ts_tab(22) = x2all
  ts_tab(23) = 1. - abs(x4all)/(x3all+eps)
!
  endif
!
!  Sensible heat flux
!
  if (out_sfl) then
    call horav (thermo_prop%cp(:,:,1)*surf%esflux, x1all)

    ts_tab(12) = x1all
    ! ts_tab(24) = x1all
!
!  Latent heat flux
!
    call horav (flv00*surf%qvflux, x1all)

    ts_tab(13) = x1all
    ! ts_tab(25) = x1all
  endif
!
!  Radiative fluxes
!
#ifdef RAD_ENABLE
  ts_tab(26) = olr
  ts_tab(27) = olr_cs
  ts_tab(28) = toa
  ts_tab(29) = toa_cs
#endif
!
!  Averaged CAPE and CIN
!
  if (out_cape) then
    call horav ( state%es, ptav )
    call horav ( state%qt, qtav )
!
    call get_cape ( nz, p0, ptav, qtav, x1all, x2all, zb=x3all, src=3 )
!
    call horav (column%cape, x1all, xq=column%cape, xc=100., filter='larger')
    call horav (column%cin, x2all)
    call horav (column%lfc, x3all)
!
    ts_tab(30) = x1all
    ts_tab(31) = x2all
    ts_tab(32) = x3all
  endif
!
!  Deallocate
!
  call dealloc ( tem )
  call dealloc ( rh )
  call dealloc ( cmax )
  call dealloc ( zi )
  call dealloc ( zb )
  call dealloc ( cc )
  call dealloc ( mfb )
  call dealloc ( ccm )
  call dealloc ( mfbm )
  call dealloc ( prec )
  call dealloc ( peff )
  call dealloc ( cwp )
  call dealloc ( pw )
  call dealloc ( crh )
!
! ================================================================
!    
RETURN
END subroutine cumulus_timeser
!      
! === deep case: includes stratopause diagnostics
! 
!      =================================================                            
  subroutine deep_timeser (ts_tab)      
!      =================================================                            
!
      integer :: i, j, k, h, nzt, nzb, nt, ntall, nzball, nztall, mint(1)
!
      real, dimension(ntsout) :: ts_tab
      real, dimension(nz)    :: xav, xav2
      real :: z, val, val2, x1, x1all, valall, denall
      
#if ( defined SPMD )      
      integer :: ierr, status(MPI_STATUS_SIZE)
#endif
!
      real, allocatable, dimension (:,:,:) :: tmp
      real, allocatable, dimension (:,:) :: zb, zi, slice, prec
!
! ================================================================
!
  ts_tab = 0.
!
!  Allocate
!
  call alloc ( zi )
  call alloc ( zb )
  call alloc ( slice )
  call alloc ( prec )
  call alloc ( tmp )
!
!  Max vertical wind 
!
  call max_all (wind%w, x1)
!
#ifdef SPMD
  x1all = 0.
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
  x1all = x1
#endif 
  ts_tab(1) = x1all
!
!  Mean cloud top height:
!
  if (lmicro > 0) then 
!
  zi  = 0.0
  nzt = 0
  do i = it_start, it_end
    do j = jt_start, jt_end
      val2 = 0.0
      z = zmax
      do k = nz, 2, -1
        val = den0(k)*hydromtr(drop)%q(i,j,k)
	if ((val >= qthres .and. wind%w(i,j,k) >= wthres) .and. (val2 < qthres .or. wind%w(i,j,k+1) < wthres)) then
	  zi(i,j) = z
	  nzt=nzt+1
	  exit
	else
	  val2 = val
	  z = z - dz*0.5*(fdz0(k)+fdz0(k-1))
	endif
      enddo
    enddo
  enddo
  x1 = sum(zi)
!
#ifdef SPMD
  nztall = 0
  x1all = 0.
  call MPI_ALLReduce (nzt, nztall, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  nztall = nzt
  x1all = x1
#endif 
!
  ts_tab(2) = x1all/real(max(nztall,1))
!
!  Max cloud top height:
!
  call max_all (zi, x1)
!
#ifdef SPMD
  x1all = 0.
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
  x1all = x1
#endif 
!
  ts_tab(3) = x1all
!
!  Mean cloud base height:
!
  zb = 0.0
  nzb = 0
  nt = 0
  do i = it_start, it_end
    do j = jt_start, jt_end
      nt = nt + 1
      z = 0.5*dz*fdz0(1)
      do k = 1, nz-1
        val = den0(k)*hydromtr(drop)%q(i,j,k)
	if (val >= qthres .and. wind%w(i,j,k) > wthres) then
	  zb(i,j) = z
	  nzb=nzb+1
	  exit
	else
	  z = z + dz*0.5*(fdz0(k)+fdz0(k+1))
	endif
      enddo
    enddo
  enddo
  x1 = sum(zb)
!
#ifdef SPMD
  nzball = 0
  x1all = 0.
  call MPI_ALLReduce (nzb, nzball, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  nzball = nzb
  x1all = x1
#endif 
!
  ts_tab(4) = x1all/real(max(nzball,1))
!
!  Cloud fraction
!
  nt = 0
  nzb = 0
  do i = it_start, it_end
    do j = jt_start, jt_end
      nt = nt + 1
      val = 0.
      do k = 1, nz-1
  	if (hydromtr(drop)%q(i,j,k) >= qthres .and. wind%w(i,j,k) >= wthres) val = val + den0(k)*hydromtr(drop)%q(i,j,k)*dz*fdz0(k)
      enddo
      if (val > 1.e-2) nzb = nzb + 1
    enddo
  enddo
!
#ifdef SPMD
  ntall = 0
  nzball = 0
  call MPI_ALLReduce (nt, ntall, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLReduce (nzb, nzball, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  ntall = nt
  nzball = nzb
#endif 
!
  ts_tab(5) = real(nzball)/real(ntall)
!
!  Mean surface precipitation in mm/day
!
  prec = 0.
  do h = 1, nhydro
    prec = prec - den0(1)*hydromtr(h)%q(:,:,1)*hydromtr(h)%vp(:,:,1)*86400.
  enddo
!
  call horav (prec, x1all)
  ts_tab(6) = x1all     	
!
!  Max surface precipitation in mm/day
!
  call max_all (prec, ts_tab(7))
!
  endif
!
!  tropopause levels
!
  call horav (thermo%T, xav)
  mint = minloc(xav)
!
  xav2 = 0.
  if (sca_set == 2 .and. nscal > 0) then
    call horav (state%scal(:,:,:,1), xav2)
  endif
!
  ts_tab(8) = xav2(mint(1))
!
  ts_tab(9) = minval(xav)
!
  ts_tab(10) = z0(mint(1))
!
!  Integrated O3 (scalar 1)
!
  x1all = 0.
  if (sca_set == 2 .and. nscal > 0) then
    call vertint ( state%scal(:,:,:,1), slice, weighted=.true. )
    call horav (slice, x1all)
  endif
!
  ts_tab(11) = x1all
!
!  Sensible heat flux
!
  if (out_sfl) then
    call horav (thermo_prop%cp(:,:,1)*surf%esflux, x1all)
!
    ts_tab(12) = x1all
!
!  Latent heat flux
!
    call horav (flv00*surf%qvflux, x1all)
!
    ts_tab(13) = x1all
  endif
!
!  Radiative fluxes
!
#ifdef RAD_ENABLE
  call horav (rad%frad, xav)
  ts_tab(14) = sum(xav)
! #endif
! !
! #ifdef RAD_ENABLE
  ts_tab(15) = ssr
  ts_tab(16) = olr
  ts_tab(17) = olr_cs
  ts_tab(18) = toa
  ts_tab(19) = toa_cs
#endif
!
#ifdef RAD_ENABLE
  call horav (rad%frad(:,:,nz), x1all)
  ts_tab(20) = x1all
#endif
!
!  Deallocate
!
  call dealloc ( zi )
  call dealloc ( zb )
  call dealloc ( slice )
  call dealloc ( prec )
  call dealloc ( tmp )
!
! ================================================================
!    
RETURN
END subroutine deep_timeser
!      
! === conservation diagnostics 
! 
!      =================================================                            
  subroutine conservation_timeser (ts_tab)      
!      =================================================                            
!
      real, dimension(ntsout) :: ts_tab
      real :: x1, x1all, x2all
      
#if ( defined SPMD )      
      integer :: ierr, status(MPI_STATUS_SIZE)
#endif
!
! ================================================================
!
  ts_tab = 0.
!
!  Max vertical wind 
!
  call max_all (wind%w, x1)
!
#ifdef SPMD
  x1all = 0.
  call MPI_ALLReduce (x1, x1all, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
  x1all = x1
#endif 
  ts_tab(1) = x1all
!
!  Potential temperature, density and water conservation
!
  ts_tab(2) = des_tot
  ts_tab(3) = drh_tot
  ts_tab(4) = dqt_tot
!
!  total density
!
  call totalav (pressure%dens, x1all)
  ts_tab(5) = x1all
!
!  total moisture
!
  call totalav (pressure%dens*state%qt, x2all)
  ts_tab(6) = x2all
!
!  total mse & energy (or pt)
!
#if (defined CONSERVATIVE) || (!defined ANELASTIC)
    call totalav (pressure%dens*state%es, x2all)
    ts_tab(7) = x2all/x1all
#else
    call totalav (state%es, x2all)
    ts_tab(7) = x2all
#endif
!
!  Sensible heat flux
!
  if (out_sfl) then
    call horav (thermo_prop%cp(:,:,1)*surf%esflux, x1all)
!
    ts_tab(18) = x1all
!
!  Latent heat flux
!
    call horav (flv00*surf%qvflux, x1all)
!
    ts_tab(19) = x1all
  endif
!
! ================================================================
!    
RETURN
END subroutine conservation_timeser
! 
!      =================================================  
  subroutine layer_budget ( h1, h2, dens, w, sca1, sca2, q1, q2, thres, filter )
!      =================================================  
!
  IMPLICIT NONE
!
  integer :: k, count
  real :: thres, h1, h2, zz
  character(len=*) :: filter
  character(len=5) :: car5
  character(len=40) :: file_name
  logical :: per_tsstep, ex
!
  real :: mean, var, laytot, layadv, layoth
  real, dimension(1:nz) :: sav1, sav2, svar2, advav, othav
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dens, w, sca1, sca2, q1, q2
!  
  real, allocatable, dimension(:,:,:) :: tmp, tmp1, tmp2
!
! ================================================================
!
!  Allocate
!
  call alloc ( tmp )
  call alloc ( tmp1 )
  call alloc ( tmp2 )
!
!  Initialise profiles
!
!  means
!
  do k = 1, nz
    tmp1(ip_start:ip_end,jp_start:jp_end,k) = den0(k)*sca1(ip_start:ip_end,jp_start:jp_end,k)
  enddo
  call horav( tmp1, sav1, xq1=q1, xc1=thres, filter1=filter )
!
  do k = 1, nz
    tmp2(ip_start:ip_end,jp_start:jp_end,k) = dens(ip_start:ip_end,jp_start:jp_end,k)*sca2(ip_start:ip_end,jp_start:jp_end,k)
  enddo
  call horav( tmp2, sav2, xq1=q2, xc1=thres, filter1=filter )
!
!  Variance
!
  do k = 1, nz
    tmp(ip_start:ip_end,jp_start:jp_end,k) = (tmp2(ip_start:ip_end,jp_start:jp_end,k) - sav2(k))**2.
  enddo
  call horav( tmp, svar2, xq1=q2, xc1=thres, filter1=filter )
!
!  Tendencies (vertical advection and misc)
!
  call horav( diag(3)%qt, advav, xq1=q1, xc1=thres, filter1=filter )
  call horav( diag(8)%qt, othav, xq1=q1, xc1=thres, filter1=filter )
!
!  Layer budget
!
  count = 0
  laytot = 0.
  layadv = 0.
  layoth = 0.
  mean = 0.
  var = 0.
  zz = 0.5*dz*fdz0(1)
  do k = 1, nz-1
    if ( zz >= h1 .and. zz < h2 ) then
      count = count + 1
      mean = mean + sav2(k)
      var = var + svar2(k)
!
      laytot = laytot + 0.25*dz*(fdz0(k)+fdz0(k+1)) * ((sav2(k) - sav1(k)) + (sav2(k+1) - sav1(k+1))) / dt0
      layadv = layadv + 0.25*dz*(fdz0(k)+fdz0(k+1)) * (advav(k) + advav(k+1))
      layoth = layoth + 0.25*dz*(fdz0(k)+fdz0(k+1)) * (othav(k) + othav(k+1))
    endif
    zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
  enddo
!
!  Write data
!
  per_tsstep = (time-tstart >= real(n_ts)*its) .or. (time>=tstop-dt0*0.9)
  if ( mypid == 0 .and. per_tsstep ) then
    write(car5,'(i5)')floor(h1)
    file_name='./OUTPUT/T_S_'//trim(adjustl(car5))
!    
    inquire(FILE=trim(file_name),EXIST=ex)
    if ( ex .and. (ntime > 1) ) then
      open(100,FILE=trim(file_name),access='SEQUENTIAL',position='APPEND')
    else
      open(100,FILE=trim(file_name))
    endif
    write(100,*) time,mean/real(count),var/real(count),laytot,layadv,layoth
    close(100)
  endif
!
!  Deallocate
!
  call dealloc ( tmp )
  call dealloc ( tmp1 )
  call dealloc ( tmp2 )
!
! ================================================================
!    
  RETURN
  END subroutine layer_budget
! 
!      =================================================  
  subroutine romps_entrainment 
!      =================================================  
!
  integer :: i, j, k
  type(atm_winds) :: wind_tmp
  real, dimension(:,:,:), allocatable :: volume1, volume2, volume1n, volume2n, volume1a, volume2a
  real, dimension(:,:,:), allocatable :: s1, s2, s1n, s2n, s1a, s2a, sca, s1adv, s2adv, s1adv0, s2adv0, subg1, subg2
  real :: delta1, delta2
!
! ================================================================
!
  entrain%activ = 0.0
  entrain%dent = 0.0
  entrain%ddet = 0.0
  entrain%cent = 0.0
  entrain%cdet = 0.0
  entrain%went = 0.0
  entrain%wdet = 0.0
  entrain%went1 = 0.0
  entrain%wdet1 = 0.0
  entrain%went2 = 0.0
  entrain%wdet2 = 0.0
  entrain%went3 = 0.0
  entrain%wdet3 = 0.0
  entrain%went4 = 0.0
  entrain%wdet4 = 0.0
  entrain%went5 = 0.0
  entrain%wdet5 = 0.0
!
  call alloc ( volume1 )
  call alloc ( volume2 )
  call alloc ( volume1n )
  call alloc ( volume2n )
  call alloc ( volume1a )
  call alloc ( volume2a )
  call alloc ( s1 )
  call alloc ( s2 )
  call alloc ( s1n )
  call alloc ( s2n )
  call alloc ( s1a )
  call alloc ( s2a )
  call alloc ( sca )
  call alloc ( s1adv )
  call alloc ( s2adv )
  call alloc ( s1adv0 )
  call alloc ( s2adv0 )
  call alloc ( subg1 )
  call alloc ( subg2 )
  call alloc ( wind_tmp )
!
  do k = 1, nz-1
    do j = jp_start, jp_end
      do i = ip_start, ip_end
        s1(i,j,k) = den0(k)*(hydromtr(drop)%q(i,j,k) + hydromtr(rain)%q(i,j,k))
        s2(i,j,k) = 0.5*den0(k)*(wind%w(i,j,k) + wind%w(i,j,k+1))
      enddo
    enddo
  enddo
!
!  Advect tracers
!
  call advection_ac ( wind2, s1, s1adv0, ilim=.true. )
  call advection_ac ( wind2, s2, s2adv0, ilim=.true. )
!
  wind_tmp = wind2
!  call boxav ( wind2%u, wind_tmp%u, 1 )
!  call boxav ( wind2%v, wind_tmp%v, 1 )
!
  call advection_ac ( wind_tmp, s1, s1adv, ilim=.true. )
  call advection_ac ( wind_tmp, s2, s2adv, ilim=.true. )
!
!  Calculate VOF volumes
!
  call volume ( 1.e-5, s1, volume1 )
!
  call volume ( 1., s2, volume2 )
!
  s1n = s1 + dt0*(diag(9)%qc + diag(9)%qr)
  call volume ( 1.e-5, s1n, volume1n )
!
  s2n = s2 + dt0*(diag(9)%w)
  call volume ( 1., s2n, volume2n )
!
  s1a = s1 + dt0*s1adv
  call volume ( 1.e-5, s1a, volume1a )
!
  s2a = s2 + dt0*s2adv
  call volume ( 1., s2a, volume2a )
!
!  Cumulate entrainment
!
  do k = 1, nz-1
    do j = jt_start, jt_end
      do i = it_start, it_end
!
        s1(i,j,k) = den0(k)*(volume1n(i,j,k) - volume1a(i,j,k)) / dt0
        s2(i,j,k) = den0(k)*(volume2n(i,j,k) - volume2a(i,j,k)) / dt0
        sca(i,j,k) = den0(k)*(volume1n(i,j,k)*volume2n(i,j,k) - volume1a(i,j,k)*volume2a(i,j,k)) / dt0
        entrain%activ(i,j,k) = volume1(i,j,k)*volume2(i,j,k) 
!
        subg1(i,j,k) = s1adv0(i,j,k) - s1adv(i,j,k)
        subg2(i,j,k) = s2adv0(i,j,k) - s2adv(i,j,k)
!
        delta1 = 0.
        delta2 = 0.
	if ( abs(diag(4)%qc(i,j,k) + diag(4)%qr(i,j,k) + diag(6)%qc(i,j,k) + diag(6)%qr(i,j,k) + diag(7)%qc(i,j,k) + diag(7)%qr(i,j,k) + diag(8)%qc(i,j,k) + diag(8)%qr(i,j,k) + subg1(i,j,k)) > 1.e-10 ) &
	delta1 = s1(i,j,k) / (diag(4)%qc(i,j,k) + diag(4)%qr(i,j,k) + diag(6)%qc(i,j,k) + diag(6)%qr(i,j,k) + diag(7)%qc(i,j,k) + diag(7)%qr(i,j,k) + diag(8)%qc(i,j,k) + diag(8)%qr(i,j,k) + subg1(i,j,k))
	if ( abs(diag(6)%w(i,j,k) + diag(7)%w(i,j,k) + diag(8)%w(i,j,k) + subg2(i,j,k)) > 1.e-5 ) &
	delta2 = s2(i,j,k) / (diag(6)%w(i,j,k) + diag(7)%w(i,j,k) + diag(8)%w(i,j,k) + subg2(i,j,k))
!
	entrain%went1(i,j,k) = (diag(4)%qc(i,j,k) + diag(4)%qr(i,j,k))*delta1*volume2(i,j,k)
	entrain%wdet1(i,j,k) = (diag(6)%w(i,j,k) + diag(8)%w(i,j,k))*delta2*volume1(i,j,k)
!
	entrain%went2(i,j,k) = (diag(6)%qr(i,j,k) + diag(8)%qc(i,j,k))*delta1*volume2(i,j,k)
	entrain%wdet2(i,j,k) = (diag(7)%w(i,j,k))*delta2*volume1(i,j,k)
!
	entrain%went3(i,j,k) = (diag(9)%qc(i,j,k) + diag(9)%qr(i,j,k))*delta1*volume2(i,j,k)
	entrain%wdet3(i,j,k) = (diag(9)%w(i,j,k))*delta2*volume1(i,j,k)
!
	entrain%went4(i,j,k) = wind_tmp%u(i,j,k) !s1adv(i,j,k)*delta1*volume2(i,j,k)
	entrain%wdet4(i,j,k) = wind_tmp%v(i,j,k) !s2adv(i,j,k)*delta2*volume1(i,j,k)
!
	entrain%went5(i,j,k) = subg1(i,j,k)*delta1*volume2(i,j,k)
	entrain%wdet5(i,j,k) = subg2(i,j,k)*delta2*volume1(i,j,k)
!
        entrain%cent(i,j,k) = max(sca(i,j,k),0.)
        entrain%cdet(i,j,k) = max(-sca(i,j,k),0.)
!
        entrain%dent(i,j,k) = max(s1(i,j,k)*volume2(i,j,k) + s2(i,j,k)*volume1(i,j,k),0.)
        entrain%ddet(i,j,k) = max(-s1(i,j,k)*volume2(i,j,k) - s2(i,j,k)*volume1(i,j,k),0.)
!  
        entrain%went(i,j,k) = max(entrain%went1(i,j,k) + entrain%wdet1(i,j,k) + entrain%went2(i,j,k) + entrain%wdet2(i,j,k) + entrain%went5(i,j,k) + entrain%wdet5(i,j,k),0.)
        entrain%wdet(i,j,k) = max(-entrain%went1(i,j,k) - entrain%wdet1(i,j,k) - entrain%went2(i,j,k) - entrain%wdet2(i,j,k) - entrain%went5(i,j,k) - entrain%wdet5(i,j,k),0.)
!
      enddo
    enddo
  enddo
!
  call dealloc ( volume1 )
  call dealloc ( volume2 )
  call dealloc ( volume1n )
  call dealloc ( volume2n )
  call dealloc ( volume1a )
  call dealloc ( volume2a )
  call dealloc ( s1 )
  call dealloc ( s2 )
  call dealloc ( s1n )
  call dealloc ( s2n )
  call dealloc ( s1a )
  call dealloc ( s2a )
  call dealloc ( sca )
  call dealloc ( s1adv )
  call dealloc ( s2adv )
  call dealloc ( s1adv0 )
  call dealloc ( s2adv0 )
  call dealloc ( subg1 )
  call dealloc ( subg2 )
  call dealloc ( wind_tmp )
!
! ================================================================
!    
  RETURN
  END subroutine romps_entrainment
! 
!      =================================================  
  subroutine direct_entrainment ( thres, dens, sca1, sca2, object )
!      =================================================  
!
  integer :: i, j, k
  real :: thres, ent
  character(len=*) :: object
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dens, sca1, sca2
!
  real, allocatable, dimension(:,:,:) :: volume1, volume2
!
! ================================================================
!
!  Initialisations
!
  entrain%activ = 0.0
  entrain%cdet = 0.0
  entrain%cent = 0.0
!
  call alloc ( volume1 )
  call alloc ( volume2 )
!
!  Calculate cloudy volume at the end of time step
!
#ifdef MODEL_3D
  call volume ( thres, sca2, volume2 )
#else
  call volume ( thres, sca2(:,1,:), volume2(:,1,:) )
#endif
!
!  Calculate cloudy volume after advection only
!
#ifdef MODEL_3D
  call volume ( thres, sca1, volume1, wind%u, wind%v, wind%w )
#else
  call volume ( thres, sca1(:,1,:), volume1(:,1,:), wind%u(:,1,:), wind%w(:,1,:) )
#endif
!
!  Calculate entrainment/detrainment rates
!
  do k = 2, nz-2
    do j = jt_start, jt_end
      do i = it_start, it_end
        ent = dens(i,j,k)*(volume2(i,j,k) - volume1(i,j,k))/dt0
!
        entrain%activ(i,j,k) = volume1(i,j,k)
        entrain%cent(i,j,k) = max(ent,0.)
        entrain%cdet(i,j,k) = abs( min(ent,0.) )
      enddo
    enddo
  enddo
!
!  Deallocate
!
  call dealloc ( volume1 )
  call dealloc ( volume2 )
!
! ================================================================
!    
  RETURN
  END subroutine direct_entrainment
! 
!      =================================================  
  subroutine level_set_entrainment  
!      =================================================  
!
  integer :: i, j, k
  real :: sca1, sca2
  real, parameter :: d1e = 5.e-12, d2e = 5.e-3
  real, allocatable, dimension(:,:,:) :: volume1, volume2, delta1, delta2, sca, sadv
!
! ================================================================
!
!  Initialisations
!
  call alloc ( volume1 )
  call alloc ( volume2 )
  call alloc ( delta1 )
  call alloc ( delta2 )
  call alloc ( sca )
  call alloc ( sadv )
!
!  Regularized functions
!
  do k = 1, nz-1
    do j = jt_start, jt_end
      do i = it_start, it_end
        sca1 = hydromtr(drop)%q(i,j,k) + hydromtr(rain)%q(i,j,k)
        sca2 = 0.5*(wind%w(i,j,k) + wind%w(i,j,k+1))
!
        volume1(i,j,k) = 1. / (1. + exp(-(sca1 - 1.e-5)/sqrt(d1e))) 
!
        volume2(i,j,k) = 1. / (1. + exp(-(sca2 - 1.)/sqrt(d2e))) 
!
        delta1(i,j,k) = volume1(i,j,k)*(1. - volume1(i,j,k)) / sqrt(d1e)
!
        delta2(i,j,k) = volume2(i,j,k)*(1. - volume2(i,j,k)) / sqrt(d2e)
!
        sca(i,j,k) = den0(k) * volume1(i,j,k) * volume2(i,j,k)
      enddo
    enddo
  enddo
!
  call advection_ac ( wind, sca, sadv )
!
!  Calculate entrainment/detrainment rates
!
  do k = 1, nz-1
    do j = jt_start, jt_end
      do i = it_start, it_end
!	  
        if ( sqrt(d2e)*delta2(i,j,k)*volume1(i,j,k) > 1.e-7 .or. sqrt(d1e)*delta1(i,j,k)*volume2(i,j,k) > 1.e-7) then
	  entrain%went1(i,j,k) = (diag(4)%qc(i,j,k) + diag(4)%qr(i,j,k))*volume2(i,j,k)*delta1(i,j,k)
	  entrain%wdet1(i,j,k) = (diag(6)%w(i,j,k) + diag(8)%w(i,j,k))*volume1(i,j,k)*delta2(i,j,k)
!
	  entrain%went2(i,j,k) = (diag(6)%qc(i,j,k) + diag(6)%qr(i,j,k) + diag(7)%qc(i,j,k) + diag(7)%qr(i,j,k) + diag(8)%qc(i,j,k) + diag(8)%qr(i,j,k))*volume2(i,j,k)*delta1(i,j,k)
	  entrain%wdet2(i,j,k) = (diag(7)%w(i,j,k))*volume1(i,j,k)*delta2(i,j,k)
!
	  entrain%went3(i,j,k) = (diag(9)%qc(i,j,k) + diag(9)%qr(i,j,k))*volume2(i,j,k)*delta1(i,j,k)
	  entrain%wdet3(i,j,k) = (diag(9)%w(i,j,k))*volume1(i,j,k)*delta2(i,j,k)
!
	  entrain%went4(i,j,k) = (diag(1)%qc(i,j,k) + diag(1)%qr(i,j,k) + diag(2)%qc(i,j,k) + diag(2)%qr(i,j,k) + diag(3)%qc(i,j,k) + diag(3)%qr(i,j,k))*volume2(i,j,k)*delta1(i,j,k)
	  entrain%wdet4(i,j,k) = (diag(1)%w(i,j,k) + diag(2)%w(i,j,k) + diag(3)%w(i,j,k))*volume1(i,j,k)*delta2(i,j,k)
!  
!          entrain%went(i,j,k) = entrain%went(i,j,k) + max(entrain%went1(i,j,k) + entrain%went2(i,j,k) + entrain%wdet1(i,j,k) + entrain%wdet2(i,j,k),0.)
          entrain%went(i,j,k) = entrain%went(i,j,k) + max(entrain%went3(i,j,k) + entrain%wdet3(i,j,k) - sadv(i,j,k),0.)
!
!          entrain%wdet(i,j,k) = entrain%wdet(i,j,k) + max(-entrain%went1(i,j,k) - entrain%went2(i,j,k) - entrain%wdet1(i,j,k) - entrain%wdet2(i,j,k),0.) 
          entrain%wdet(i,j,k) = entrain%wdet(i,j,k) + max(-entrain%went3(i,j,k) - entrain%wdet3(i,j,k) + sadv(i,j,k),0.) 
        endif
!
      enddo
    enddo
  enddo
!
!  Deallocate
!
  call dealloc ( volume1 )
  call dealloc ( volume2 )
  call dealloc ( delta1 )
  call dealloc ( delta2 )
  call dealloc ( sca )
  call dealloc ( sadv )
!
! ================================================================
!    
  RETURN
  END subroutine level_set_entrainment
! 
!      =================================================  
  subroutine filtered_profiles ( ent, det, sca1, sca2, thres1, thres2, object )
!      =================================================  
!
  integer :: ierr, k, i, j, l, m, n, counte, countd, counta, nne, nnd, nna, nobj
  logical :: per_prstep
  real :: epse, qte, pte, ptve, buoye, beffe, we
  real :: epsd, qtd, ptd, ptvd, buoyd, beffd, wd
  real :: epsa, qta, pta, ptva, buoya, beffa, wa
  real, dimension(1:nz) :: bav
!
  character(len=*) :: object
  real :: thres1, thres2
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ent, det, sca1, sca2
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: ctop, cbas
  logical, dimension(ip_start:ip_end,jp_start:jp_end) :: lcloud
!
  type (filtered), dimension(nz) :: loc
!
! ================================================================
!
!  Initialize
!
  if ( trim(object) == 'wcore' ) then
    nobj = 1
  else if ( trim(object) == 'ccore' ) then
    nobj = 4
  endif
!
  per_prstep = (time-tstart >= real(n_pr)*ipro) .or. (time>=tstop-dt0*0.9)
!
  if ( per_prstep ) then
!
  ctop = 0.
  cbas = 0.
  do k = 1, nz
    where ( sca1(:,:,k) > qthres .and. cbas == 0. ) cbas = z0(k)
  enddo
  do k = nz,1,-1
    where ( sca1(:,:,k) > qthres .and. ctop == 0. ) ctop = z0(k)    
  enddo
!
  if ( out_buoy ) then
    call horav ( thermo%buoy, bav ) 
  else
    bav = 0.
  endif

  do l = 1, 3
!
  loc%epse = 0.
  loc%qte = 0.
  loc%pte = 0.
  loc%ptve = 0.
  loc%buoye = 0.
  loc%beffe = 0.
  loc%we = 0.
  loc%epsd = 0.
  loc%qtd = 0.
  loc%ptd = 0.
  loc%ptvd = 0.
  loc%buoyd = 0.
  loc%beffd = 0.
  loc%wd = 0.
  loc%epsa = 0.
  loc%qta = 0.
  loc%pta = 0.
  loc%ptva = 0.
  loc%buoya = 0.
  loc%beffa = 0.
  loc%wa = 0.
!
  lcloud = .false.
  if ( l == 1 ) then
    lcloud = .true.
  else if ( l == 2 ) then
    where ( cbas < 2000. .and. ctop <= 2500. ) lcloud = .true.
  else if ( l == 3 ) then
    where ( cbas < 2000. .and. ctop >= 3000. ) lcloud = .true.
  endif
!
!  Find filtered state
!
    do k = 1, nz-1
      epse = 0.
      qte = 0.
      pte = 0.
      ptve = 0.
      buoye = 0.
      beffe = 0.
      we = 0.
      counte = 0
      epsd = 0.
      qtd = 0.
      ptd = 0.
      ptvd = 0.
      buoyd = 0.
      beffd = 0.
      wd = 0.
      countd = 0
      epsa = 0.
      qta = 0.
      pta = 0.
      ptva = 0.
      buoya = 0.
      beffa = 0.
      wa = 0.
      counta = 0
!
      do j = jt_start, jt_end
        do i = it_start, it_end
!
!  Entrained points (outside cloud)
!
	  if ( ent(i,j,k) > thres1 ) then
!
#ifdef MODEL_3D
	    do n = j-1, j+1, 2
	      do m = i-1, i+1, 2
	        if ( sca2(m,n,k) < thres2 .and. lcloud(i,j) ) then
    	          counte = counte + 1
		  epse = epse + ent(i,j,k)
	          qte = qte + state%qt(m,n,k)
	          pte = pte + state%es(m,n,k)
	          ptve = ptve + thermo%ptv(m,n,k)
	          buoye = buoye + thermo%buoy(m,n,k) - bav(k)
	          beffe = beffe + diagnos%beff(m,n,k)
	          we = we + wind%w(m,n,k)
	        endif
	      enddo
	    enddo
#else
	    do m = i-1, i+1, 2
	      if ( sca2(m,1,k) < thres2 .and. lcloud(i,1) ) then
    	    	counte = counte + 1
		epse = epse + ent(i,1,k)
	    	qte = qte + state%qt(m,1,k)
	        pte = pte + state%es(m,1,k)
	        ptve = ptve + thermo%ptv(m,1,k)
	    	buoye = buoye + thermo%buoy(m,1,k) - bav(k)
	    	beffe = beffe + diagnos%beff(m,1,k)
	    	we = we + wind%w(m,1,k)
	      endif
	    enddo
#endif
!
	  endif
!
!  Detrained points (inside cloud)
!
	  if ( abs(det(i,j,k)) > thres1 ) then
!
#ifdef MODEL_3D
	    do n = j-1, j+1, 2
	      do m = i-1, i+1, 2
	        if ( sca2(m,n,k) >= thres2 .and. lcloud(m,n) ) then
    	          countd = countd + 1
		  epsd = epsd + det(i,j,k)
	          qtd = qtd + state%qt(m,n,k)
	          ptd = ptd + state%es(m,n,k)
	          ptvd = ptvd + thermo%ptv(m,n,k)
	          buoyd = buoyd + thermo%buoy(m,n,k) - bav(k)
	          beffd = beffd + diagnos%beff(m,n,k)
	          wd = wd + wind%w(m,n,k)
	        endif
	      enddo
	    enddo
#else
	    do m = i-1, i+1, 2
	      if ( sca2(m,1,k) >= thres2 .and. lcloud(m,1) ) then
    	    	countd = countd + 1
		epsd = epsd + det(i,1,k)
	    	qtd = qtd + state%qt(m,1,k)
	        ptd = ptd + state%es(m,1,k)
	        ptvd = ptvd + thermo%ptv(m,1,k)
	    	buoyd = buoyd + thermo%buoy(m,1,k) - bav(k)
	    	beffd = beffd + diagnos%beff(m,1,k)
	    	wd = wd + wind%w(m,1,k)
	      endif
	    enddo
#endif
!
	  endif
!
!  All boundary points
!
	  if ( (ent(i,j,k) > thres1 .and. lcloud(i,j)) .or. (abs(det(i,j,k)) > thres1 .and. 	&
	  	(lcloud(i+1,j) .or. lcloud(i,j) .or. lcloud(i-1,j) .or. lcloud(i+1,j+1) .or. lcloud(i,j+1) .or. lcloud(i-1,j+1) .or. lcloud(i+1,j-1) .or. lcloud(i,j-1) .or. lcloud(i-1,j-1))) ) then
!
#ifdef MODEL_3D
    	    counta = counta + 1
	    epsa = epsa + ent(i,j,k) + det(i,j,k)
	    qta = qta + state%qt(i,j,k)
	    pta = pta + state%es(i,j,k)
	    ptva = ptva + thermo%ptv(i,j,k)
	    buoya = buoya + thermo%buoy(i,j,k) - bav(k)
	    beffa = beffa + diagnos%beff(i,j,k)
	    wa = wa + wind%w(i,j,k)
#else
    	    counta = counta + 1
	    epsa = epsa + ent(i,1,k) + det(i,1,k)
	    qta = qta + state%qt(i,1,k)
	    pta = pta + state%es(i,1,k)
	    ptva = ptva + thermo%ptv(i,1,k)
	    buoya = buoya + thermo%buoy(i,1,k) - bav(k)
	    beffa = beffa + diagnos%beff(i,1,k)
	    wa = wa + wind%w(i,1,k)
#endif
!
	  endif
!
	enddo
      enddo
!
!  Store entrained properties  
!
#ifdef SPMD
      nne = 0
      nnd = 0
      nna = 0
!
      call MPI_ALLReduce (epse, loc(k)%epse, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (qte, loc(k)%qte, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (pte, loc(k)%pte, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (ptve, loc(k)%ptve, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (buoye, loc(k)%buoye, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (beffe, loc(k)%beffe, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (we, loc(k)%we, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (counte, nne, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
!      
      call MPI_ALLReduce (epsd, loc(k)%epsd, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (qtd, loc(k)%qtd, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (ptd, loc(k)%ptd, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (ptvd, loc(k)%ptvd, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (buoyd, loc(k)%buoyd, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (beffd, loc(k)%beffd, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (wd, loc(k)%wd, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (countd, nnd, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
!      
      call MPI_ALLReduce (epsa, loc(k)%epsa, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (qta, loc(k)%qta, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (pta, loc(k)%pta, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (ptva, loc(k)%ptva, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (buoya, loc(k)%buoya, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (beffa, loc(k)%beffa, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (wa, loc(k)%wa, 1, REALTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLReduce (counta, nna, 1, INTTYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
!
      loc(k)%epse=loc(k)%epse/max(real(nne),1.)
      loc(k)%qte=loc(k)%qte/max(real(nne),1.)
      loc(k)%pte=loc(k)%pte/max(real(nne),1.)
      loc(k)%ptve=loc(k)%ptve/max(real(nne),1.)
      loc(k)%buoye=loc(k)%buoye/max(real(nne),1.)
      loc(k)%beffe=loc(k)%beffe/max(real(nne),1.)
      loc(k)%we=loc(k)%we/max(real(nne),1.)
!
      loc(k)%epsd=loc(k)%epsd/max(real(nnd),1.)
      loc(k)%qtd=loc(k)%qtd/max(real(nnd),1.)
      loc(k)%ptd=loc(k)%ptd/max(real(nnd),1.)
      loc(k)%ptvd=loc(k)%ptvd/max(real(nnd),1.)
      loc(k)%buoyd=loc(k)%buoyd/max(real(nnd),1.)
      loc(k)%beffd=loc(k)%beffd/max(real(nnd),1.)
      loc(k)%wd=loc(k)%wd/max(real(nnd),1.)
!
      loc(k)%epsa=loc(k)%epsa/max(real(nna),1.)
      loc(k)%qta=loc(k)%qta/max(real(nna),1.)
      loc(k)%pta=loc(k)%pta/max(real(nna),1.)
      loc(k)%ptva=loc(k)%ptva/max(real(nna),1.)
      loc(k)%buoya=loc(k)%buoya/max(real(nna),1.)
      loc(k)%beffa=loc(k)%beffa/max(real(nna),1.)
      loc(k)%wa=loc(k)%wa/max(real(nna),1.)
#else
      loc(k)%epse=epse/max(real(counte),1.)
      loc(k)%qte=qte/max(real(counte),1.)
      loc(k)%pte=pte/max(real(counte),1.)
      loc(k)%ptve=ptve/max(real(counte),1.)
      loc(k)%buoye=buoye/max(real(counte),1.)
      loc(k)%beffe=beffe/max(real(counte),1.)
      loc(k)%we=we/max(real(counte),1.)
      
      loc(k)%epsd=epsd/max(real(countd),1.)
      loc(k)%qtd=qtd/max(real(countd),1.)
      loc(k)%ptd=ptd/max(real(countd),1.)
      loc(k)%ptvd=ptvd/max(real(countd),1.)
      loc(k)%buoyd=buoyd/max(real(countd),1.)
      loc(k)%beffd=beffd/max(real(countd),1.)
      loc(k)%wd=wd/max(real(countd),1.)
      
      loc(k)%epsa=epsa/max(real(counta),1.)
      loc(k)%qte=qte/max(real(counte),1.)
      loc(k)%pta=pta/max(real(counta),1.)
      loc(k)%ptva=ptva/max(real(counta),1.)
      loc(k)%buoya=buoya/max(real(counta),1.)
      loc(k)%beffa=beffa/max(real(counta),1.)
      loc(k)%wa=wa/max(real(counta),1.)
#endif
!
    enddo
!
!  Output if required
!    
    call write_1d_profile ( nobj+l-1, (n_pr == 0), loc=loc )    
!
  enddo
!
  endif
!
! ================================================================
!
  RETURN       
  END subroutine filtered_profiles   
! 
!      =================================================  
  subroutine plume_entrainment ( dens, w, p, ptv, b, qt, sca1, object )
!      =================================================  
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: sca1, sca2, dens, w, qt, ptv, p, b
!
  character(len=*) :: object
  logical :: per_prstep
  integer :: k, l, nobj
  real :: mflthres=1.e4
  logical, dimension(ip_start:ip_end,jp_start:jp_end) :: lcloud
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: ctop, cbas
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: unity, tmp
  real, dimension(1:nz) :: cover, qt_cor, ptv_cor, w_cor, gp_cor, b_cor
  real, dimension(1:nz) :: resq, resw, respt, qt_env, w_env, ptv_env 
  real, dimension(1:nz) :: mflc, qflc, tflc, wflc, gpflc, bflc, wvflc, wsflc, dqc_dt, dtc_dt, dwc_dt, dcov_dt
  real, dimension(1:nz) :: gq_cor, g_mfl, g_qfl, g_tfl, g_wfl, g_cov, gw_cor, gt_cor
!
  real, dimension(1:nz,3), save :: cover_w_old=0., qt_w_old=0., ptv_w_old=0., w_w_old=0., cover_b_old=0., qt_b_old=0., ptv_b_old=0., w_b_old=0.
!
  type (plumes), dimension(nz) :: plume
!
! ================================================================
!
!  Initialize
!
  if ( trim(object) == 'wcore' ) then
    nobj = 7
    sca2 = w
  else if ( trim(object) == 'bcore' ) then
    nobj = 10
    sca2 = b
  endif
!
  ctop = 0.
  cbas = 0.
  do k = 1, nz
    where ( sca1(:,:,k) > qthres .and. cbas == 0. ) cbas = z0(k)
  enddo
  do k = nz,1,-1
    where ( sca1(:,:,k) > qthres .and. ctop == 0. ) ctop = z0(k)    
  enddo
!
  per_prstep = (time-tstart >= real(n_pr)*ipro) .or. (time>=tstop-dt0*0.9)
  
  do l = 1, 3
!
  plume%qte = 0.
  plume%ptve = 0.
  plume%we = 0.
  plume%qt = 0.
  plume%ptv = 0.
  plume%w = 0.
  plume%qtb = 0.
  plume%ptvb = 0.
  plume%wb = 0.
  plume%dmfl = 0.
  plume%mfl = 0.
  plume%cov = 0.
  plume%dadt = 0.
  plume%mu = 0.
  plume%eq = 0.
  plume%et = 0.
  plume%ew = 0.
  plume%eqs = 0.
  plume%ets = 0.
  plume%ews = 0.
  plume%alphaq = 0.
  plume%alphat = 0.
  plume%alphaw = 0.
  plume%alphab = 0.
  plume%alphap = 0.
  plume%beta = 0.
  plume%gamma = 0.
  plume%delta1 = 0.
  plume%delta2 = 0.
  plume%wvar = 0.
  plume%wske = 0.
  plume%ri = 0.
  plume%mri = 0.
!
  lcloud = .false.
  if ( l == 1 ) then
    lcloud = .true.
  else if ( l == 2 ) then
    where ( cbas < 2000. .and. ctop <= 2500. ) lcloud = .true.
  else if ( l == 3 ) then
    where ( cbas < 2000. .and. ctop >= 3000. ) lcloud = .true.
  endif
!
!  Qt averages
!
  call horav ( qt, qt_env, xq1=sca1, xc1=qthres, filter1='lower' ) 
!
  call horav ( qt, qt_cor, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger' ) 
!
  call horav ( ptv, ptv_env, xq1=sca1, xc1=qthres, filter1='lower' ) 
!
  call horav ( ptv, ptv_cor, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger' ) 
!
  call horav ( w, w_env, xq1=sca1, xc1=qthres, filter1='lower' ) 
!
  call horav ( w, w_cor, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger' ) 
!
  call horav ( diag(7)%w+diag(8)%w, gp_cor, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger' )
!
  call horav ( diag(6)%w, b_cor, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger' )
!
!  Cloud cover
! 
  unity = 0.
  where ( sca1 >= qthres ) unity = dens
!
  call horsum ( unity, cover, scale=dx*dy )
!
  call grad_1d ( cover, g_cov ) 
!
!  Time variations
!
  if ( nobj == 7 ) then
!
    dqc_dt = (qt_cor - qt_w_old(:,l)) / dt0
!
    dtc_dt = (ptv_cor - ptv_w_old(:,l)) / dt0
!
    dwc_dt = (w_cor - w_w_old(:,l)) / dt0
!
    dcov_dt = (cover - cover_w_old(:,l)) / dt0
!
    qt_w_old(:,l) = qt_cor
!
    ptv_w_old(:,l) = ptv_cor
!
    w_w_old(:,l) = w_cor
!
    cover_w_old(:,l) = cover
! 
  else if ( nobj == 10 ) then
!
    dqc_dt = (qt_cor - qt_b_old(:,l)) / dt0
!
    dtc_dt = (ptv_cor - ptv_b_old(:,l)) / dt0
!
    dwc_dt = (w_cor - w_b_old(:,l)) / dt0
!
    dcov_dt = (cover - cover_b_old(:,l)) / dt0
!
    qt_b_old(:,l) = qt_cor
!
    ptv_b_old(:,l) = ptv_cor
!
    w_b_old(:,l) = w_cor
!
    cover_b_old(:,l) = cover
!
  endif
!
!  Calculate fluxes
!
    call horsum ( dens*w, mflc, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
    call horsum ( dens*w*qt, qflc, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
    call horsum ( dens*w*ptv, tflc, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
    call horsum ( dens*w*w, wflc, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
    call horsum ( w*diag(7)%w, gpflc, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
    call horsum ( w*diag(6)%w, bflc, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
!  W moments
!
  if ( per_prstep ) then
!
    do k = 1, nz
      tmp(:,:,k) = (w(:,:,k) - w_cor(k))**2.
    enddo
    call horav ( dens*w*tmp, wvflc, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger' ) 
!    
    do k = 1, nz
      tmp(:,:,k) = (w(:,:,k) - w_cor(k))**3.
    enddo
    call horav ( dens*w*tmp, wsflc, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger' ) 
!
!  Plume core gradients
!
    call grad_1d ( qt_cor, gq_cor ) 
!
    call grad_1d ( ptv_cor, gt_cor ) 
!
    call grad_1d ( w_cor, gw_cor ) 
!
    call grad_1d ( mflc, g_mfl ) 
!
    call grad_1d ( qflc, g_qfl ) 
!
    call grad_1d ( tflc, g_tfl ) 
!
    call grad_1d ( wflc, g_wfl ) 
!
!  Calculate residual profiles
!
    call horsum ( diag(5)%qt+diag(6)%qt+diag(7)%qt+diag(8)%qt,   			&
      resq, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
    call horsum ( diag(4)%ptv+diag(5)%ptv+diag(6)%ptv+diag(7)%ptv+diag(8)%ptv,		&
      respt, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
    call horsum ( diag(4)%w+diag(5)%w+diag(6)%w+diag(7)%w+diag(8)%w,      		&
      resw, xq1=sca1, xc1=qthres, filter1='larger', xflag1=lcloud, xq2=sca2, xc2=0., filter2='larger', scale=dx*dy )
!
!  Store plume diagnostics
!
    where ( abs(mflc) > mflthres )
      plume%qt  = qt_cor
!
      plume%ptv  = ptv_cor
!
      plume%w  = w_cor
!
      plume%qte  = qt_env
!
      plume%ptve  = ptv_env
!
      plume%we  = w_env
!
      plume%mfl  = mflc
!
      plume%cov  = cover
!
      plume%dmfl  = g_mfl
!
      plume%dadt  = dcov_dt
!
      plume%mu = (g_mfl + dcov_dt) / mflc
!
      plume%eq = ( resq - g_qfl + qt_cor*g_mfl - cover*dqc_dt ) / ( mflc*(qt_cor - qt_env) )
!
      plume%et = ( respt - g_tfl + ptv_cor*g_mfl - cover*dtc_dt ) / ( mflc*(ptv_cor - ptv_env) )
!
      plume%ew = ( resw - g_wfl + w_cor*g_mfl - cover*dwc_dt ) / ( mflc*(w_cor - w_env) )
!
      plume%eqs = ( resq/mflc - gq_cor ) / (qt_cor - qt_env)
!
      plume%ets = ( respt/mflc - gt_cor ) / (ptv_cor - ptv_env)
!
      plume%ews = ( resw/mflc - gw_cor ) / (w_cor - w_env)
!
      plume%qtb = qt_cor - ( resq - g_qfl + qt_cor*g_mfl - cover*dqc_dt ) / (g_mfl + dcov_dt)
!
      plume%ptvb = ptv_cor - ( respt - g_tfl + ptv_cor*g_mfl - cover*dtc_dt ) / (g_mfl + dcov_dt)
!
      plume%wb = w_cor - ( resw - g_wfl + w_cor*g_mfl - cover*dwc_dt ) / (g_mfl + dcov_dt)
!
      plume%alphaq = qflc / (mflc*qt_cor)
!
      plume%alphat = tflc / (mflc*ptv_cor)
!
      plume%alphaw = wflc / (mflc*w_cor)
!
      plume%alphab = bflc / (mflc*b_cor)
!
      plume%alphap = gpflc / (mflc*gp_cor)
!
      plume%beta = 1. + gp_cor / b_cor
!
      plume%wske = wsflc / mflc
    end where
!
!  Output  
!
    call write_1d_profile ( nobj+l-1, (n_pr == 0), plum=plume )   
!
  endif
!
  enddo 
!
! ================================================================
!
  RETURN       
  END subroutine plume_entrainment 
  
end module diagnostics  
