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
! ADVECTION.F90:
!	module advection.f90 
!
!  Purpose:
!	Provides an interface to calculating scalar advection.
!	Possibility to split advection along each direction, and
!	handles conservative and non-conservative advection.
!	Advection fluxes and tendencies are calculated in the
!	routines dedicated to each advection scheme available:
!	MUSCL, QUICK, PPM, LW, WENO
!
!  Author:
!	Julien Savre, Ludwig Maximilian Universitat, Munich
!
! ================================================================
!
  module advection

!  
!  ---------------------------------------------
!
!  General modules
!
  USE gridno
  USE shared_data
  USE shared_pressure
  USE shared_hydro
  USE shared_wind
  USE shared_thermo, only : thermo_prop
  USE shared_diag
!
!  Advection modules
!
  USE allocation
  USE boundary_conditions
  USE gradients
  USE funcpack
  USE advection_lw
  USE advection_lw_nl
  USE advection_ppm
  USE advection_muscl
  USE advection_quick
  USE advection_quick_nl
  USE advection_weno

  IMPLICIT NONE

  private 
 
  logical :: lim_def, with_dens, with_es, with_qt, with_hydro, with_state, lfirst=.true.
  
  save
  
  public :: advect, advection_x, advection_y, advection_z, advection_ac, advection_1d, advection_ac_1d 

  contains
!
!  =============================================
   subroutine advect ( dtl, dens, statel, hydromtrl,	&
                       densa, statea, hydromtra, limiter, full )
!  =============================================
!
  integer :: h
  logical :: lim, ful
  real :: dt1
!	 
  real, intent(in) :: dtl	
  real, dimension(ip_start:,jp_start:,:), optional, intent(in) :: dens
  type(atm_state), optional, intent(in) :: statel
  type(hydrometeor), dimension(:), optional, intent(in) :: hydromtrl
!
  real, dimension(ip_start:,jp_start:,:), optional, intent(inout) :: densa
  type(atm_state), optional, intent(inout) :: statea
  type(hydrometeor), dimension(:), optional, intent(inout) :: hydromtra
!
  logical, optional, intent(in) :: limiter, full
!
!  Local
!
  type(atm_state) :: state_rho
  type(hydrometeor), dimension(1:nhydro) :: hydromtr_rho  
  real, allocatable, dimension(:,:,:) :: dens_rho
!
#ifdef ADV_SPLIT
  type(atm_state) :: state_tmp
  type(hydrometeor), dimension(1:nhydro) :: hydromtr_tmp  
  real, allocatable, dimension(:,:,:) :: dens_tmp
#endif
!
  if (verbose > 1) call write_debug('Entering advection')
!
! -----------------------------------------------------------------
!
  dt1 = 1./dtl
!
!  Initialize advection
!
  if (lfirst) call advection_init 
  lfirst = .false.
!
!  Set switches
!
  if (present(limiter)) then
    lim = limiter
  else
    lim = lim_def
  endif
!
  if (present(full)) then
    ful = full
  else
    ful = .true.
  endif
!
  with_dens=.false.; with_es=.false.; with_qt=.false.; with_hydro=.false.; with_state=.false.
!
  if (present(densa).or.conservative) with_dens=.true.
  if (present(statel)) with_state=.true.
  if (present(statel)) with_es=.true.
  if (present(statel).and.ful) with_qt=.true.
  if (present(hydromtrl).and.present(hydromtra).and.nhydro>0) with_hydro=.true.
!
!-----------------------------------------------------------------!
!                 Initialize conserved variables                  !
!-----------------------------------------------------------------!
!
  call alloc ( dens_rho )
  call alloc ( state_rho )
  call alloc ( hydromtr_rho )
!
  if (with_dens) dens_rho = dens
  if (with_state) state_rho = statel
  if (with_hydro) then
    do h = 1, nhydro
      hydromtr_rho(h) = hydromtrl(h)
    enddo
  endif
!
!-----------------------------------------------------------------!
!   	   	Second step: Advection, no splitting 	          !
!   + divide by prognostic density at the end for conservation    !
!-----------------------------------------------------------------!
!
#ifndef ADV_SPLIT
!
!  Advance scalars along x
!
  call qnadv_scal_x ( dtl, wind_adv%u, dens_rho, state_rho, hydromtr_rho,	&
        	      dens_rho, state_rho, hydromtr_rho, 		&
		      densa, statea, hydromtra, lim )
!
!  Advance scalars along y
!
#ifdef MODEL_3D
  call qnadv_scal_y ( dtl, wind_adv%v, dens_rho, state_rho, hydromtr_rho,	&
   	   	      dens_rho, state_rho, hydromtr_rho, 		&
		      densa, statea, hydromtra, lim )
#endif
!
!  Advance scalars along z
!
  call qnadv_scal_z ( dtl, wind_adv%w, dens_rho, state_rho, hydromtr_rho,	&
   	 	      dens_rho, state_rho, hydromtr_rho, 		&
		      densa, statea, hydromtra, lim )
!
!-----------------------------------------------------------------!
!   	     Second step: Advection with Strang splitting         !
!	   alternate directions for odd and even time steps	  !
! + divide by prognostic density after each step for conservation !
!-----------------------------------------------------------------!
!
#else
!
!  Allocate temporary variables
!
  if (with_dens) call alloc ( dens_tmp )
  if (with_state) call alloc ( state_tmp )
  if (with_hydro) call alloc ( hydromtr_tmp )
!
  if (with_dens) dens_tmp = dens
  if (with_state) state_tmp = statel
  if (with_hydro) then
    do h = 1, nhydro
      hydromtr_tmp(h) = hydromtrl(h)
    enddo
  endif
!
  if (mod(ntime,2) == 0) then
!
!  Advance scalars along x
!
    call qnadv_scal_x ( 0.5*dtl, wind_adv%u, dens_tmp, state_tmp, hydromtr_tmp,	&
     	 	        dens_rho, state_rho, hydromtr_rho, 		&
			densa, statea, hydromtra, lim )
!
!  Advance scalars along y
!
#ifdef MODEL_3D
    call qnadv_scal_y ( dtl, wind_adv%v, dens_tmp, state_tmp, hydromtr_tmp,	&
     	 	        dens_rho, state_rho, hydromtr_rho, 		&
			densa, statea, hydromtra, lim )
#endif
!
!  Advance scalars along x
!
    call qnadv_scal_x ( 0.5*dtl, wind_adv%u, dens_tmp, state_tmp, hydromtr_tmp,	&
     	 	        dens_rho, state_rho, hydromtr_rho, 		&
			densa, statea, hydromtra, lim )
!
!  Advance scalars along z
!
    call qnadv_scal_z ( dtl, wind_adv%w, dens_tmp, state_tmp, hydromtr_tmp,	&
     	 	        dens_rho, state_rho, hydromtr_rho, 		&
			densa, statea, hydromtra, lim )
!
  else
!
!  Advance scalars along y
!
#ifdef MODEL_3D
    call qnadv_scal_y ( 0.5*dtl, wind_adv%v, dens_tmp, state_tmp, hydromtr_tmp,	&
     	 	        dens_rho, state_rho, hydromtr_rho, 		&
			densa, statea, hydromtra, lim )
#endif
!
!  Advance scalars along x
!
    call qnadv_scal_x ( dtl, wind_adv%u, dens_tmp, state_tmp, hydromtr_tmp,	&
     	 	        dens_rho, state_rho, hydromtr_rho, 		&
			densa, statea, hydromtra, lim )
!
!  Advance scalars along y
!
#ifdef MODEL_3D
    call qnadv_scal_y ( 0.5*dtl, wind_adv%v, dens_tmp, state_tmp, hydromtr_tmp,	&
     	 	        dens_rho, state_rho, hydromtr_rho, 		&
			densa, statea, hydromtra, lim )
#endif
!
!  Advance scalars along z
!
    call qnadv_scal_z ( dtl, wind_adv%w, dens_tmp, state_tmp, hydromtr_tmp,	&
     	 	        dens_rho, state_rho, hydromtr_rho, 		&
			densa, statea, hydromtra, lim )
!  
  endif
!
!  Deallocate temp variables
!
  if (with_dens) call dealloc ( dens_tmp )
  if (with_state) call dealloc ( state_tmp )
  if (with_hydro) call dealloc ( hydromtr_tmp )
!
!  Assemble advective tendencies
!
  if (with_dens) densa = densa + (dens_rho - dens)*dt1
!
  if (with_es) statea%es = statea%es + (state_rho%es - statel%es)*dt1
!
  if (with_qt) statea%qt = statea%qt + (state_rho%qt - statel%qt)*dt1
!
  if (with_hydro) then
    do h = 1, nhydro
      hydromtra(h)%q = hydromtra(h)%q + (hydromtr_rho(h)%q - hydromtrl(h)%q)*dt1
!
#ifdef SEIFERT
      if (moments == 2) hydromtra(h)%n = hydromtra(h)%n + (hydromtr_rho(h)%n - hydromtrl(h)%n)*dt1
#endif
!
#ifdef SEIFERT
      if (lmicro > 3 .and. h > ice) hydromtra(h)%w = hydromtra(h)%w + (hydromtr_rho(h)%w - hydromtrl(h)%w)*dt1
#endif
    enddo
  endif
!
#endif
!
!  Deallocate temp variables
!
  call dealloc ( dens_rho )
  call dealloc ( state_rho )
  call dealloc ( hydromtr_rho )
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating advection')
!
return                  
end subroutine advect                                                    
!
!  =============================================
   subroutine  qnadv_scal_x ( dtl, u, dens, statel, hydromtrl,		&
   	 	              dens_rho, state_rho, hydromtr_rho, 	&
			      denss, statea, hydromtra, lim )

!  =============================================

  real :: dtl
  integer :: h
  logical :: lim
  character(len=1) :: car1
!	 
! === Note: all parameters are exported with advection sign	
  type(atm_state), optional :: statel, state_rho, statea
  type(hydrometeor), dimension(1:nhydro), optional :: hydromtrl, hydromtr_rho, hydromtra
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), optional :: dens_rho, denss
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dens, u
!
!-----------------------------------------------------------------!
!  		  Advection of scalars along x    		  !
!-----------------------------------------------------------------!	
!
!  Advance density first for mass consistency
!
  if (with_dens) call advection_x ( 'de', dtl, dens_rho, dens, denss, u, lim, adv=.false. )
!
!  Advance all other scalars
!
  if (with_es) call advection_x ( 'es ', dtl, state_rho%es, statel%es, statea%es, u, lim )
!
  if (with_qt) call advection_x ( 'qt ', dtl, state_rho%qt, statel%qt, statea%qt, u, lim )
!
  if (with_hydro) then
    do h = 1, nhydro
      write(car1,'(i1)') h
!
      call advection_x ( 'q'//car1, dtl, hydromtr_rho(h)%q, hydromtrl(h)%q, hydromtra(h)%q, u, lim )
!
#ifdef SEIFERT
      if (moments == 2) call advection_x ( 'n'//car1, dtl, hydromtr_rho(h)%n, hydromtrl(h)%n, hydromtra(h)%n, u, lim )
#endif
!
#ifdef SEIFERT
      if ( lmicro > 3 .and. h > ice ) call advection_x ( 'w'//car1, dtl, hydromtr_rho(h)%w, hydromtrl(h)%w, hydromtra(h)%w, u, lim )
#endif
!
    enddo
  endif
!
!-----------------------------------------------------------------!
!  		   	 Density weighting  	  		  !
!-----------------------------------------------------------------!	
!
!  Update Boundary values
!
#ifdef ADV_SPLIT
  if (with_dens) call statebc ( dens_rho )
!
  if (with_es) call statebc ( state_rho%es )
!
  if (with_qt) call statebc ( state_rho%qt )
!
  if (with_hydro) call hydrobc ( hydromtr_rho )
!
  if (with_dens) dens = dens_rho
  if (with_es) statel%es = state_rho%es
  if (with_qt) statel%qt = state_rho%qt
  if (with_hydro) then
    do h = 1, nhydro
      hydromtrl(h) = hydromtr_rho(h) 
    enddo
  endif
#endif
!
!-----------------------------------------------------------------!
!
  return
!
  end subroutine qnadv_scal_x
!
!  =============================================
   subroutine  qnadv_scal_y ( dtl, v, dens, statel, hydromtrl,		&
   	 	              dens_rho, state_rho, hydromtr_rho,	&
			      denss, statea, hydromtra, lim )

!  =============================================
!
  real :: dtl
  integer :: h
  logical :: lim
  character(len=1) :: car1
!	 
! === Note: all parameters are exported with advection sign	
  type(atm_state), optional :: statel, state_rho, statea
  type(hydrometeor), dimension(1:nhydro), optional :: hydromtrl, hydromtr_rho, hydromtra
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), optional :: dens_rho, denss
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dens, v
!
!-----------------------------------------------------------------!
!  		    Advection of scalars along x    		  !
!-----------------------------------------------------------------!	
!
!  Advance density first for mass consistency
!
  if (with_dens) call advection_y ( 'de', dtl, dens_rho, dens, denss, v, lim, adv=.false. )
!
!  Advance all other scalars
!
  if (with_es) call advection_y ( 'es ', dtl, state_rho%es, statel%es, statea%es, v, lim )
!
  if (with_qt) call advection_y ( 'qt ', dtl, state_rho%qt, statel%qt, statea%qt, v, lim )
!
  if (with_hydro) then
    do h = 1, nhydro
      write(car1,'(i1)') h
!
      call advection_y ( 'q'//car1, dtl, hydromtr_rho(h)%q, hydromtrl(h)%q, hydromtra(h)%q, v, lim )
!
#ifdef SEIFERT
      if (moments == 2) call advection_y ( 'n'//car1, dtl, hydromtr_rho(h)%n, hydromtrl(h)%n, hydromtra(h)%n, v, lim )
#endif
!
#ifdef SEIFERT
      if ( lmicro > 3 .and. h > ice ) call advection_y ( 'w'//car1, dtl, hydromtr_rho(h)%w, hydromtrl(h)%w, hydromtra(h)%w, v, lim )
#endif
!
    enddo
  endif
!
!-----------------------------------------------------------------!
!  			 Density weighting   			  !
!-----------------------------------------------------------------!	
!
!  Update Boundary values
!
#ifdef ADV_SPLIT
  if (with_dens) call statebc ( dens_rho )
!
  if (with_es) call statebc ( state_rho%es )
!
  if (with_qt) call statebc ( state_rho%qt )
!
  if (with_hydro) call hydrobc ( hydromtr_rho )
!
  if (with_dens) dens = dens_rho
  if (with_es) statel%es = state_rho%es
  if (with_qt) statel%qt = state_rho%qt
  if (with_hydro) then
    do h = 1, nhydro
      hydromtrl(h) = hydromtr_rho(h) 
    enddo
  endif
#endif
!
!-----------------------------------------------------------------!
!
  return
!
  end subroutine qnadv_scal_y
!
!  =============================================
   subroutine  qnadv_scal_z ( dtl, w, dens, statel, hydromtrl,		&
   	 	              dens_rho, state_rho, hydromtr_rho, 	&
			      denss, statea, hydromtra, lim )

!  =============================================
!
  real :: dtl
  integer  :: h, k
  logical :: lim
  character(len=1) :: car1
!	 
! === Note: all parameters are exported with advection sign	
  type(atm_state), optional :: statel, state_rho, statea
  type(hydrometeor), dimension(1:nhydro), optional :: hydromtrl, hydromtr_rho, hydromtra
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), optional :: dens_rho, denss
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dens, w
!
!-----------------------------------------------------------------!
!  		    Advection of scalars along z    		  !
!-----------------------------------------------------------------!	
!
!  Advance density first for mass consistency
!
  if (with_dens) call advection_z ( 'de', dtl, dens_rho, dens, denss, w, lim, adv=.false. )
!
!  Advance all other scalars
!
  if (with_es) call advection_z ( 'es ', dtl, state_rho%es, statel%es, statea%es, w, lim )
!
  if (with_qt) call advection_z ( 'qt ', dtl, state_rho%qt, statel%qt, statea%qt, w, lim )
!
  if (with_hydro) then
    do h = 1, nhydro
      write(car1,'(i1)') h
!
      call advection_z ( 'q'//car1, dtl, hydromtr_rho(h)%q, hydromtrl(h)%q, hydromtra(h)%q, w, lim )
!
#ifdef SEIFERT
      if (moments == 2) call advection_z ( 'n'//car1, dtl, hydromtr_rho(h)%n, hydromtrl(h)%n, hydromtra(h)%n, w, lim )
#endif
!
#ifdef SEIFERT
      if ( lmicro > 3 .and. h > ice ) call advection_z ( 'w'//car1, dtl, hydromtr_rho(h)%w, hydromtrl(h)%w, hydromtra(h)%w, w, lim )
#endif
!
    enddo
  endif
!
!-----------------------------------------------------------------!
!  	    	         Density weighting    			  !
!-----------------------------------------------------------------!	
!
!  Update Boundary values
!
#ifdef ADV_SPLIT
  if (with_dens) call statebc ( dens_rho )
!
  if (with_es) call statebc ( state_rho%es )
!
  if (with_qt) call statebc ( state_rho%qt )
!
  if (with_hydro) call hydrobc ( hydromtr_rho )
!
  if (with_dens) dens = dens_rho
  if (with_es) statel%es = state_rho%es
  if (with_qt) statel%qt = state_rho%qt
  if (with_hydro) then
    do h = 1, nhydro
      hydromtrl(h) = hydromtr_rho(h) 
    enddo
  endif
#endif
!
!-----------------------------------------------------------------!
!
  return
  end subroutine qnadv_scal_z
!
!  =============================================
  subroutine advection_x (name, dtl, xrho, xtmp, xs, u, lim, adv)
!  =============================================
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: xtmp, xx, xs, xrho
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: u      
!
  integer :: j, k
  character(len=2) :: name
  logical, optional :: adv
  logical :: adv2, lim
  real :: dtl
!
!-----------------------------------------------------------------!
!
  xx(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
!
!  Advective form?
!
  if (present(adv)) then
    adv2 = adv
  else
    adv2 = .not.conservative
  endif
!  
!  Advective tendencies
!
  if (.not.lim) then
    do k = 1, nz
      do j = jt_start, jt_end
        select case ( trim(scal_adv) )
        case ('quick')
          call advsx_quick_nl (it_start,it_end, xtmp(:,j,k), u(:,j,k), xx(:,j,k), adv2)
        case default
          call advsx_lw_nl (it_start, it_end, xtmp(:,j,k), u(:,j,k), xx(:,j,k), adv2)
        end select
      enddo
    enddo
  else
!
!  Flux limited advection
!
  do k = 1, nz
    do j = jt_start, jt_end
!
      if (sorder == 1) then
        call advsx_lw (it_start,it_end, xtmp(:,j,k), u(:,j,k), xx(:,j,k), adv2)
      else
        select case( trim(scal_adv) )
        case ( 'weno' )
          call advsx_weno (it_start,it_end, xtmp(:,j,k), u(:,j,k), xx(:,j,k), adv2)
        case ( 'quick' )
          call advsx_quick (it_start,it_end, xtmp(:,j,k), u(:,j,k), xx(:,j,k), adv2)
        case ( 'ppm' )
          call advsx_ppm (it_start,it_end, xtmp(:,j,k), u(:,j,k), xx(:,j,k), adv2)
        case ( 'muscl' )
          call advsx_muscl (it_start,it_end, xtmp(:,j,k), u(:,j,k), xx(:,j,k), adv2)
        case default
          call advsx_lw (it_start,it_end, xtmp(:,j,k), u(:,j,k), xx(:,j,k), adv2)
        end select
      endif
!
    enddo
  enddo
!
  endif
!
#ifdef ADV_SPLIT
  xrho = xrho - xx*dtl 
#else
  xs = xs - xx
#endif
!
  if (out_diagt .and. trim(name)=='es') diag(1)%pt = diag(1)%pt - cdiag*xx
  if (out_diagq .and. trim(name)=='qt') diag(1)%qt = diag(1)%qt - cdiag*xx
  if (out_diagl .and. trim(name)=='q1') diag(1)%qc = diag(1)%qc - cdiag*xx
  if (out_diagr .and. trim(name)=='q2') diag(1)%qr = diag(1)%qr - cdiag*xx
  if (out_diagi .and. trim(name)=='q3') diag(1)%qi = diag(1)%qi - cdiag*xx
!
  end subroutine advection_x
!
!  =============================================
  subroutine advection_y (name, dtl, xrho, xtmp, xs, v, lim, adv)
!  =============================================
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: xtmp, xrho, xy, xs
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: v
!
  integer :: i, k
  character(len=2) :: name
  logical, optional :: adv
  logical :: adv2, lim
  real :: dtl
!
!-----------------------------------------------------------------!
!
  xy(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
!
!  Advective form?
!
  if (present(adv)) then
    adv2 = adv
  else
    adv2 = .not.conservative
  endif
!
!  Advective tendencies
!
  if (.not.lim) then
    do k = 1, nz
      do i = it_start, it_end
        select case ( trim(scal_adv) )
        case ('quick')
          call advsy_quick_nl (jt_start,jt_end, xtmp(i,:,k), v(i,:,k), xy(i,:,k), adv2)
        case default
          call advsy_lw_nl (jt_start,jt_end, xtmp(i,:,k), v(i,:,k), xy(i,:,k), adv2)
        end select
      enddo
    enddo
  else
!
!  Flux limited advection
!
  do k = 1, nz
    do i = it_start, it_end
!
      if (sorder == 1) then
        call advsy_lw (jt_start,jt_end, xtmp(i,:,k), v(i,:,k), xy(i,:,k), adv2)
      else
        select case( trim(scal_adv) )
        case ( 'weno' )
          call advsy_weno (jt_start,jt_end, xtmp(i,:,k), v(i,:,k), xy(i,:,k), adv2)
        case ( 'quick' )
          call advsy_quick (jt_start,jt_end, xtmp(i,:,k), v(i,:,k), xy(i,:,k), adv2)
        case ( 'ppm' )
          call advsy_ppm (jt_start,jt_end, xtmp(i,:,k), v(i,:,k), xy(i,:,k), adv2)
        case ( 'muscl' )
          call advsy_muscl (jt_start,jt_end, xtmp(i,:,k), v(i,:,k), xy(i,:,k), adv2)
        case default
          call advsy_lw (jt_start,jt_end, xtmp(i,:,k), v(i,:,k), xy(i,:,k), adv2)
        end select
      endif
!
    enddo
  enddo
!
  endif
!
#ifdef ADV_SPLIT
  xrho = xrho - xy*dtl
#else
  xs = xs - xy
#endif
!
  if (out_diagt .and. trim(name)=='es') diag(2)%pt = diag(2)%pt - cdiag*xy
  if (out_diagq .and. trim(name)=='qt') diag(2)%qt = diag(2)%qt - cdiag*xy
  if (out_diagl .and. trim(name)=='q1') diag(2)%qc = diag(2)%qc - cdiag*xy
  if (out_diagr .and. trim(name)=='q2') diag(2)%qr = diag(2)%qr - cdiag*xy
  if (out_diagi .and. trim(name)=='q3') diag(2)%qi = diag(2)%qi - cdiag*xy
  if (out_diags .and. trim(name)=='s1') diag(2)%sca(:,:,:,min(1,nscal)) = diag(2)%sca(:,:,:,min(1,nscal)) - cdiag*xy
!
  end subroutine advection_y
!
!  =============================================
  subroutine advection_z (name, dtl, xrho, xtmp, xs, w, lim, adv)
!  =============================================
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: xtmp, xz, xs, xrho
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: w
!
  integer :: i, j
  character(len=2) :: name
  logical, optional :: adv
  logical :: adv2, lim
  real :: dtl
!
!-----------------------------------------------------------------!
!
  xz(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
!
!  Advective form?
!
  if (present(adv)) then
    adv2 = adv
  else
    adv2 = .not.conservative
  endif
!
!  Advective tendencies
!
  if (.not.lim) then
    do j = jt_start, jt_end
      do i = it_start, it_end
        select case ( trim(scal_adv) )
        case ('quick')
          call advsz_quick_nl (1, nz, xtmp(i,j,:), w(i,j,:), xz(i,j,:), adv2)
        case default
          call advsz_lw_nl (1, nz, xtmp(i,j,:), w(i,j,:), xz(i,j,:), adv2)
        end select
      enddo
    enddo
  else
!
!  Flux limited advection
!
  do j = jt_start, jt_end
    do i = it_start, it_end
!
      if (sorder == 1) then
        call advsz_lw (1, nz, xtmp(i,j,:), w(i,j,:), xz(i,j,:), adv2)
      else
        select case( trim(scal_adv) )
        case ( 'weno' )
          call advsz_weno (1, nz, xtmp(i,j,:), w(i,j,:), xz(i,j,:), adv2)
        case ( 'quick' )
          call advsz_quick (1, nz, xtmp(i,j,:), w(i,j,:), xz(i,j,:), adv2)
        case ( 'ppm' )
          call advsz_ppm (1, nz, xtmp(i,j,:), w(i,j,:), xz(i,j,:), adv2)
        case ( 'muscl' )
          call advsz_muscl (1, nz, xtmp(i,j,:), w(i,j,:), xz(i,j,:), adv2)
        case default
          call advsz_lw (1, nz, xtmp(i,j,:), w(i,j,:), xz(i,j,:), adv2)
        end select
      endif
!
    enddo
  enddo
!
  endif
!
#ifdef ADV_SPLIT
  xrho = xrho - xz*dtl
#else 
  xs = xs - xz
#endif
!
  if (out_diagt .and. trim(name)=='es') diag(3)%pt = diag(3)%pt - cdiag*xz
  if (out_diagq .and. trim(name)=='qt') diag(3)%qt = diag(3)%qt - cdiag*xz
  if (out_diagl .and. trim(name)=='q1') diag(3)%qc = diag(3)%qc - cdiag*xz
  if (out_diagr .and. trim(name)=='q2') diag(3)%qr = diag(3)%qr - cdiag*xz
  if (out_diagi .and. trim(name)=='q3') diag(3)%qi = diag(3)%qi - cdiag*xz
  if (out_diags .and. trim(name)=='s1') diag(3)%sca(:,:,:,min(1,nscal)) = diag(3)%sca(:,:,:,min(1,nscal)) - cdiag*xz
!  
  end subroutine advection_z
!  
! ===================================================
  subroutine advection_ac ( windl, xx, xs, ilim, iscal, iaero )
! ===================================================
  
    Implicit none

    type(atm_winds), intent(in) :: windl
    real, dimension(ip_start:,jp_start:,:), intent(inout) :: xs
    real, dimension(ip_start:,jp_start:,:), intent(in) :: xx
    integer :: flag
    integer, optional :: iscal, iaero
    logical, optional :: ilim
!    
    integer :: k
    logical :: lim, aerod
    character(len=2) :: car2
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: xsx, xsy, xsz, xtmp, xso
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!
!
      car2 = 'xx'
!
      if (present(ilim)) then
        lim = ilim
      else
        lim = lim_def
      endif
!
      aerod = .false.
      if (present(iaero)) then
        if (iaero > 0) aerod = .true.
      endif 
!  
      xtmp(ip_start:ip_end,jp_start:jp_end,1:nz) = xx(ip_start:ip_end,jp_start:jp_end,1:nz)
      xsx(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.
      xsy(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.
      xsz(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.
!
!----------------------------------------------------------!
!                        Advection                         !
!----------------------------------------------------------!
!
#ifndef ADV_SPLIT
!
      call advection_x (car2, dt0, xtmp, xx, xsx, windl%u, lim)
!
      if (out_diaga .and. aerod .and. iaero<200 .and. iaero>100) diag(1)%aero(iaero-100)%n = diag(1)%aero(iaero-100)%n + cdiag*xsx
      if (out_diaga .and. aerod .and. iaero>200) diag(1)%aero(iaero-200)%m = diag(1)%aero(iaero-200)%m + cdiag*xsx
      if (out_diags .and. present(iscal)) diag(1)%sca(:,:,:,iscal) = diag(1)%sca(:,:,:,iscal) + cdiag*xsx
!
#ifdef MODEL_3D
      call advection_y (car2, dt0, xtmp, xx, xsy, windl%v, lim)
!
      if (out_diaga .and. aerod .and. iaero<200 .and. iaero>100) diag(2)%aero(iaero-100)%n = diag(2)%aero(iaero-100)%n + cdiag*xsy
      if (out_diaga .and. aerod .and. iaero>200) diag(2)%aero(iaero-200)%m = diag(2)%aero(iaero-200)%m + cdiag*xsy
      if (out_diags .and. present(iscal)) diag(2)%sca(:,:,:,iscal) = diag(2)%sca(:,:,:,iscal) + cdiag*xsy
#endif
!  
      call advection_z (car2, dt0, xtmp, xx, xsz, windl%w, lim)
!
      if (out_diaga .and. aerod .and. iaero<200 .and. iaero>100) diag(3)%aero(iaero-100)%n = diag(3)%aero(iaero-100)%n + cdiag*xsz
      if (out_diaga .and. aerod .and. iaero>200) diag(3)%aero(iaero-200)%m = diag(3)%aero(iaero-200)%m + cdiag*xsz
      if (out_diags .and. present(iscal)) diag(3)%sca(:,:,:,iscal) = diag(3)%sca(:,:,:,iscal) + cdiag*xsz
!
!  Get tendency
!
      xs(ip_start:ip_end,jp_start:jp_end,1:nz) = xs(ip_start:ip_end,jp_start:jp_end,1:nz) +  &
      	 xsx(ip_start:ip_end,jp_start:jp_end,1:nz) + xsy(ip_start:ip_end,jp_start:jp_end,1:nz) + xsz(ip_start:ip_end,jp_start:jp_end,1:nz)
!
!  Split advection
!
#else
!
      call advection_x (car2, dt0, xtmp, xtmp, xsx, windl%u, lim)
!
      call statebc ( xtmp )
!
      xsx = (xtmp - xx)/dt0
!
      if (out_diaga .and. aerod) diag(1)%aero(iaero)%n = diag(1)%aero(iaero)%n + cdiag*xsx
      if (out_diags .and. present(iscal)) diag(1)%sca(:,:,:,iscal) = diag(1)%sca(:,:,:,iscal) + cdiag*xsx
      xso = xsx
!
#ifdef MODEL_3D
      call advection_y (car2, dt0, xtmp, xtmp, xsy, windl%v, lim)
!
      call statebc ( xtmp )
!
      xsy = (xtmp - xx)/dt0      
!
      if (out_diaga .and. aerod) diag(2)%aero(iaero)%n = diag(2)%aero(iaero)%n + cdiag*(xsy - xso)
      if (out_diags .and. present(iscal)) diag(1)%sca(:,:,:,iscal) = diag(1)%sca(:,:,:,iscal) + cdiag*(xsy - xso)
      xso = xsy
#endif
!
      call advection_z (car2, dt0, xtmp, xtmp, xsz, windl%w, lim)
!
      call statebc ( xtmp )
!
      xsz = (xtmp - xx)/dt0      
!
      if (out_diaga .and. aerod) diag(3)%aero(iaero)%n = diag(3)%aero(iaero)%n + cdiag*(xsz - xso)
      if (out_diags .and. present(iscal)) diag(3)%sca(:,:,:,iscal) = diag(3)%sca(:,:,:,iscal) + cdiag*(xsz - sxo)
!
!  Get tendency
!
      xs(ip_start:ip_end,jp_start:jp_end,1:nz) = xs(ip_start:ip_end,jp_start:jp_end,1:nz) +  &
          (xtmp(ip_start:ip_end,jp_start:jp_end,1:nz) - xx(ip_start:ip_end,jp_start:jp_end,1:nz)) / dt0
#endif
!
!----------------------------------------------------------!
!
  return
  end subroutine advection_ac
!  
! ===================================================
  subroutine advection_1d ( w )
! ===================================================
  
  Implicit none

  real, dimension(:), intent(in) :: w
!  
  integer :: h, i, j, k
  character(len=1) :: car1
  real, dimension(1:nz) :: gxz
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dref
!
!----------------------------------------------------------!
!                        Advection                         !
!----------------------------------------------------------!
!
!  Reference density
!
#if (defined CONSERVATIVE)
#ifdef ANELASTIC
  do k = 1, nz
    dref(:,:,k) = den0(k)
  enddo
#else
  dref = pressure%dens
#endif
#else
  dref = 1.
#endif
!
!  Advective tendencies
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      gxz = 0.
      call grad_1d ( state%es(i,j,:)/dref(i,j,:), gxz )
      states%es(i,j,:) = states%es(i,j,:) - dref(i,j,:)*w*gxz
      if (out_diagt) diag(6)%pt(i,j,:) = diag(6)%pt(i,j,:) - cdiag*dref(i,j,:)*w*gxz
    enddo
  enddo
!
!  Optional
!
  if (with_lssub) then
!
    do j = jt_start, jt_end
      do i = it_start, it_end
        gxz = 0.
        call grad_1d ( state%qt(i,j,:)/dref(i,j,:), gxz ) 
        states%qt(i,j,:) = states%qt(i,j,:) - dref(i,j,:)*w*gxz
        if (out_diagq) diag(6)%qt(i,j,:) = diag(6)%qt(i,j,:) - cdiag*dref(i,j,:)*w*gxz
      enddo
    enddo
!
!    if (lmicro > 0) then
!      do h = 1, nhydro
!        do j = jt_start, jt_end
!          do i = it_start, it_end
!            gxz = 0.
!            call grad_1d ( hydromtr(h)%q(i,j,:)/dref(i,j,:), gxz ) 
!            hydromtrs(h)%q(i,j,:) = hydromtrs(h)%q(i,j,:) - dref(i,j,:)*w*gxz
!
!            if (moments > 1) then
!              gxz = 0.
!              call grad_1d ( hydromtr(h)%n(i,j,:)/dref(i,j,:), gxz ) 
!              hydromtrs(h)%n(i,j,:) = hydromtrs(h)%n(i,j,:) - dref(i,j,:)*w*gxz
!            endif
!          enddo
!        enddo
!      enddo
!    endif
!
  endif
!
!----------------------------------------------------------!
!
  return
  end subroutine advection_1d
!  
! ===================================================
  subroutine advection_ac_1d ( iscal, w, xx, xs )
! ===================================================
  
  Implicit none

  real, dimension(:), intent(in) :: w
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: xx, xs
!  
  integer :: iscal, i, j, k
  real, dimension(1:nz) :: xsz
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) ::  dref
!
!----------------------------------------------------------!
!                        Advection                         !
!----------------------------------------------------------!
!
!  Reference density
!
#if (defined CONSERVATIVE)
#ifdef ANELASTIC
  do k = 1, nz
    dref(:,:,k) = den0(k)
  enddo
#else
  dref = pressure%dens
#endif
#else
  dref = 1.
#endif
!
!  Advective scalar tendencies
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      xsz = 0.
      call grad_1d ( xx(i,j,:)/dref(i,j,:), xsz ) 
      xs(i,j,:) = xs(i,j,:) - dref(i,j,:)*w*xsz
      if (out_diags .and. iscal > 0) diag(6)%sca(i,j,:,iscal) = diag(6)%sca(i,j,:,iscal) - cdiag*dref(i,j,:)*w*xsz
    enddo
  enddo
!
!----------------------------------------------------------!
!
  return
  end subroutine advection_ac_1d
!
! ===================================================
  subroutine advection_init			     
! ===================================================
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!  
!
!  Scalar advection options
!
  lim_def = .true.
  if (.not.limit) lim_def = .false.
  if (.not.limit) ppmlimiter = 'no-limit'
!
#ifdef CONSERVATIVE
  conservative=.true.
#else
  conservative=.false.
#endif
!
!-----------------------------------------------------------!
!
  return                  
  end subroutine advection_init                                                    
  
  end module advection
