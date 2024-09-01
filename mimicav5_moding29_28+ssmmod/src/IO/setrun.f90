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
!  SETRUN.F:                   
!
!  Purpose:
!      To set up a model run.
!      
!      setrun_1d: using horizontally homogeneous initial data.
!      setrun_3d: using 2/3d initial data.
!                              
!  Special Notes:
!      The unit of initial P data should be in hPa (mb)
!                      T data should be in K
!                     qv data should be in kg/kg
!  
!      From the original setrun.F
!
!  Author
!      Chien Wang
!      MIT Joint Program on Science and Policy of Global Change
!            
! ================================================================

!      ===============================================
      subroutine setrun 
!      ===============================================

USE shared_all
USE thermodynamics
USE radiationmod
USE sources
USE subgrid
USE initialize
USE allocation
USE aerosols
USE lagrange
USE diagnostics
USE micro_diagnostics
USE boundary_conditions
!
#ifdef SPMD
USE mpi
#endif
!
IMPLICIT NONE
!
#if (defined CHEM_ENABLE)
#include "chemdef.h"
#endif
!
  integer :: k, h, ierr
!       
  save
!
!----------------------------------------------------------!
!                Get the initial values                    !
!----------------------------------------------------------!
!
 if (verbose > 0) call write_debug('Starting setrun')
!
  dt = dt0
!
  if (mypid == 0) then
    call define_grid (nz, dz, fdz0, z0)
!
    call getstatus
!
    call prt_conf
  end if
!
!  Transmitting initial data to all procs
!
#if ( defined SPMD )
  call MPI_BCAST(kbl,    1, INTTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dz,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)    
  call MPI_BCAST(zmax,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(fdz0, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dz1,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dz1w, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dzw,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
!
  call MPI_BCAST(u00,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(v00,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(qt0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(qv0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(rh0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(pt0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(t0,   nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(mse0, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ql0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(qi0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nl0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(w0,   nz, REALTYPE,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(k0,   nz, REALTYPE,0,MPI_COMM_WORLD,ierr)      
  call MPI_BCAST(z0,   nz, REALTYPE,0,MPI_COMM_WORLD,ierr)      
!
  call MPI_BCAST(s_scal0,nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(p0,     nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(p1,     nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(den0,   nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(stab0,  nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(n20,    nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(g0,     nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(avden0, nz, REALTYPE,0,MPI_COMM_WORLD,ierr)
!
  call MPI_BCAST(xrbubble,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(yrbubble,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(zrbubble,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xbubble,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ybubble,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(zbubble,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dtbubble,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
!
  call MPI_BCAST(ptop,          1, REALTYPE,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(scal0,  nz*nscal, REALTYPE,0,MPI_COMM_WORLD,ierr)      
  call MPI_BCAST(pxx,   nz*nhydro, REALTYPE,0,MPI_COMM_WORLD,ierr) 
!
  call MPI_BCAST_Forcing
#endif

!
!----------------------------------------------------------!
!                  Initializing winds                      !
!----------------------------------------------------------!
!
  do k=1,nz
    wind2%u(:,:,k) = u00(k)
#ifdef MODEL_3D
    wind2%v(:,:,k) = v00(k) 
#endif
    wind2%w(:,:,k) = 0.0 
  end do
! 
!----------------------------------------------------------!
!                  Initializing scalars                    !
!----------------------------------------------------------!
!
  do k=1,nz
    state2%qt(:,:,k) = qt0(k)
#ifdef ISENTROPIC
    state2%es(:,:,k) = pt0(k)
#else
    state2%es(:,:,k) = mse0(k)
#endif
    do h = 1, nscal
      state2%scal(:,:,k,h) = scal0(k,h)
    enddo
!
    pressure2%p(:,:,k) = 0.
    pressure2%dens(:,:,k) = den0(k)
  end do
!
  do h = 1, nscal
    if ( maxval(scal0(:,h)) == 0. ) then
      call scalar_reset ( h, state2%scal(:,:,:,h), zbl, .true. )
      scal0(:,h) = state2%scal(it_start,jt_start,:,h)
    endif
  enddo
!
!  TKE
!      
#ifdef TKE
  turbu%ksgs = 0.0
  if (maxval(k0) == 0.) then
    call tkeinit ( pressure2, wind2, state2, hydromtr2 )
  else
    do k=1,nz
      turbu%ksgs(:,:,k) = k0(k)      
    enddo
  endif
#endif
!
!----------------------------------------------------------!
!             Initializing hydrometeors to 0               !
!----------------------------------------------------------!
!   
  if (lmicro > 0) then         
    do k=1,nz
      hydromtr2(drop)%q(:,:,k) = ql0(k)
      hydromtr2(drop)%n(:,:,k) = 0.
      if (lmicro > 3) hydromtr2(drop)%w(:,:,k) = 0.0
      do h = rain, nhydro
        hydromtr2(h)%q(:,:,k) = 0.
        hydromtr2(h)%n(:,:,k) = 0.
        if (lmicro > 3) hydromtr2(h)%w(:,:,k) = 0.0
      enddo
    end do
  endif
!
!----------------------------------------------------------!
!                        Aerosols                          ! 
!----------------------------------------------------------!
!
  call init_aerosol
!
  call fill_aerosol
!
!----------------------------------------------------------!
!              Set Initial perturbations                   !
!----------------------------------------------------------!
!            
  call case_dependent
!
!----------------------------------------------------------!
!              Initial Boundary Conditions                 ! 
!----------------------------------------------------------!
!
  call windbc ( wind2 )
!
  call statebc ( state2 )      
!
  call pressurebc ( pressure2 )
!
  call hydrobc ( hydromtr2 )
!
!----------------------------------------------------------!
!                Initializing microphysics                 !
!----------------------------------------------------------!
!
  if ( with_mic .and. lmicro > 0 ) then
    if ( lndrop == 1 ) call init_activ ( pressure2, state2, hydromtr2, nuc2, aero3d2 ) 
!
    call diagnostic_micro ( pressure2, wind2, state2, hydromtr2, aero3d2 )
  endif
!
!----------------------------------------------------------!
!                 Initializing chemistry                   !
!----------------------------------------------------------!
!      
!  gas phase chemicals only
!
  if ( mypid .eq. 0 ) then 
    call chem_m_init         !note gas0 needed by radiation
  end if
!
#ifdef CHEM_ENABLE
  do k=1,nz
    do j=jp_start,jp_end
      do i=ip_start,ip_end
        gas2(i,j,k) = gas0 (k)
!
#ifdef AQCHEM_ENABLE
        aqc2(i,j,k) = 0.0
        aqr2(i,j,k) = 0.0
#endif
#ifdef SOLIDCHEM_ENABLE
        solidi2(i,j,k) = 0.0
#endif
      end do
    end do
  end do
!
!  No SO2, H2SO4 above 9th level
!
  do k=1,nz
    do j=jp_start,jp_end
      do i=ip_start,ip_end
        if (k >= 9) then
          gas2(i,j,k)%h2so4  = 0.0
          gas2(i,j,k)%so2    = 0.0
        endif
      end do
    end do 
  end do 
#endif
!
!----------------------------------------------------------!
!		   Other initializations		   !
!----------------------------------------------------------!
!
!  Fill data structures
!
  if (ldiag) call reset ( diag )
!
  if (ldiag) call reset ( column )
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
  nuc = nuc2 
#ifdef AERO_ENABLE
  aero3d = aero3d2
#endif
#ifdef NUC_CNT
  nucin = nucin2 
#endif
!
!  Initialization for piggy-backing
!
  if ( with_piggy ) then
    statep = state2
!
    hydromtrp = hydromtr2
  endif
!
!----------------------------------------------------------!
!                Initialize thermodynamics                 !
!----------------------------------------------------------!
!
#ifdef ANELASTIC
  call equation_of_state ( state, hydromtr, pressure_in=pressure%p, density_out=pressure%dens, thermo_out=thermo )
#else
  call equation_of_state ( state, hydromtr, density_in=pressure%dens, pressure_out=pressure%p, thermo_out=thermo )
#endif
!
!----------------------------------------------------------!
!                Set radiation properties                  !
!----------------------------------------------------------!
!            
  ! if (with_rad) call init_rad ( nz, p0 )
  ! if (with_rad .and. with_radmod .and. time==t_radmod) call init_rad(nz, p0)
!
!----------------------------------------------------------!
!             Initialize Lagrangian parcels                !
!----------------------------------------------------------!
!            
#ifdef LAGRANGE
  call lagrange_init
#endif
!
!----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating setrun')
!
return
end
!
subroutine MPI_BCAST_Forcing

use shared_data
use mpi

select case (trim(casename) )
  case('SE_ATLANTIC')
    call MPI_BCAST(lev_len,    1, INTTYPE,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(time_len,   1, INTTYPE,0,MPI_COMM_WORLD,ierr)   
    
    if (mypid /= 0) then
      !- Allocation of variables : 
      ALLOCATE(  Na_force(lev_len,time_len))  
      ALLOCATE(  qt_nud(lev_len,time_len)) 
      !ALLOCATE(  PT(lev_len,time_len)) 
      ALLOCATE(  ta_nud(lev_len,time_len))
      ALLOCATE(  pa_force(lev_len,time_len))  
      ALLOCATE(  u_nud(lev_len,time_len))
      ALLOCATE(  v_nud(lev_len,time_len))
      !ALLOCATE(  w_nud(lev_len,time_len))
      ALLOCATE(  ts_force(time_len))
      ALLOCATE(  ps_force(time_len))
    endif

    call MPI_BCAST(Na_force,   lev_len*time_len, REALTYPE,0,MPI_COMM_WORLD,ierr)            
    call MPI_BCAST(ta_nud,     lev_len*time_len, REALTYPE,0,MPI_COMM_WORLD,ierr) 
    call MPI_BCAST(qt_nud,     lev_len*time_len, REALTYPE,0,MPI_COMM_WORLD,ierr) 
    call MPI_BCAST(pa_force,   lev_len*time_len, REALTYPE,0,MPI_COMM_WORLD,ierr)                           
    call MPI_BCAST(u_nud,      lev_len*time_len, REALTYPE,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(v_nud,      lev_len*time_len, REALTYPE,0,MPI_COMM_WORLD,ierr)  
    call MPI_BCAST(ts_force,   time_len, REALTYPE,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ps_force,   time_len, REALTYPE,0,MPI_COMM_WORLD,ierr)  
  end select

end subroutine

