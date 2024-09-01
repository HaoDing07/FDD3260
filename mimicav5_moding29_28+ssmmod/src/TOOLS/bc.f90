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

!   =================================================
!
!   BC.F:  A package for calculating boundary
!          conditions of various quantities
!
!   Author:   Julien Savre, MISU
!
!   ==================================================

  module boundary_conditions

USE gridno
USE shared_data
USE shared_state
USE shared_pressure
USE shared_wind
USE shared_hydro
USE shared_aerosol_new
!
#ifdef SPMD
USE mpicomm
#endif
!
IMPLICIT NONE
!
  private
  
  interface windbc
    module procedure lbcu_p4, windbc_p4
  end interface windbc
  
  interface statebc
    module procedure qnbc_p4, statebc_p4
  end interface statebc
  
  interface hydrobc
    module procedure hydrobc_p4
  end interface hydrobc
  
  interface aerobc
    module procedure aerobc_p4
  end interface aerobc
  
  interface pressurebc
    module procedure pressurebc_p4
  end interface pressurebc
  
  interface gradbc
    module procedure lbcg_p4, gnbc_p4
  end interface gradbc

  public :: windbc, statebc, hydrobc, aerobc, pressurebc, gradbc
  public :: fluxbcu_x, fluxbcv_x, fluxbcw_x, fluxbcu_y, fluxbcv_y, fluxbcw_y

  CONTAINS
!
!  ====================================================
  subroutine lbcu_p4 ( u, v, w )
!  ====================================================
 
!
! --- Subroutine for calculating periodic-type of
! --- 	lateral boundary values of wind speeds.
! ---
! --- Last Revision:  Aug 3, 2011
!
!
!	Periodic boundary defined as:
!
!		     L.B.
!	+-----|-----+-----|-----+
!      u(1)  f(1)  u(2)  f(2) u(3)
!   u(nx-1) f(n-1) u(nx) f(n) u(nx+1) 	
!             	    |
!	   forced by average value		
!
!	-------------------------------------------------
!
!
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: u
#ifdef MODEL_3D
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: v
#else
  real, dimension(:,:,:), intent(inout) :: v
#endif
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: w
!
  if (verbose > 1) call write_debug('Starting windbc')
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
#ifndef SPMD
!
!  Single CPU mode
!
  if (.not.nest_run .and. (bcl(1) == 'per' .or. bcl(2) == 'per')) then
    u(1,jp_start:jp_end,1:nz) = u(nx-4,jp_start:jp_end,1:nz) 
    u(2,jp_start:jp_end,1:nz) = u(nx-3,jp_start:jp_end,1:nz)
    u(3,jp_start:jp_end,1:nz) = u(nx-2,jp_start:jp_end,1:nz)
!
    u(nx,jp_start:jp_end,1:nz) = u(5,jp_start:jp_end,1:nz)
    u(nx-1,jp_start:jp_end,1:nz) = u(4,jp_start:jp_end,1:nz)
!
#ifdef MODEL_3D
    v(1,jp_start:jp_end,1:nz) = v(nx-4,jp_start:jp_end,1:nz) 
    v(2,jp_start:jp_end,1:nz) = v(nx-3,jp_start:jp_end,1:nz)
    v(3,jp_start:jp_end,1:nz) = v(nx-2,jp_start:jp_end,1:nz)
!
    v(nx,jp_start:jp_end,1:nz) = v(5,jp_start:jp_end,1:nz)
    v(nx-1,jp_start:jp_end,1:nz) = v(4,jp_start:jp_end,1:nz)
#endif
!
    w(1,jp_start:jp_end,1:nz) = w(nx-4,jp_start:jp_end,1:nz) 
    w(2,jp_start:jp_end,1:nz) = w(nx-3,jp_start:jp_end,1:nz)
    w(3,jp_start:jp_end,1:nz) = w(nx-2,jp_start:jp_end,1:nz)
!        
    w(nx,jp_start:jp_end,1:nz) = w(5,jp_start:jp_end,1:nz)
    w(nx-1,jp_start:jp_end,1:nz) = w(4,jp_start:jp_end,1:nz)
  endif
!        
!----------------------------------------------------------!
!             For parallel run, periodic BCs               !
!              are treated via MPI exchange                !
!	 This also handles communication between	   !
!	      neighbouring internal domains	  	   !
!----------------------------------------------------------!
!
#else
    call exchange_x (u, 100)
!
#ifdef MODEL_3D     
    call exchange_x (v, 101)
#endif
!
    call exchange_x (w, 102)
#endif
!
!----------------------------------------------------------!
!                 Periodicity along y                      !
!        (assuming 1D domain decomposition along x)        !
!----------------------------------------------------------!
!
#ifndef CHANNEL
!
#ifdef MODEL_3D
!
#ifndef DECOMP_2D
  if (.not.nest_run .and. (bcl(3) == 'per' .or. bcl(4) == 'per')) then
    u(ip_start:ip_end,1,1:nz) = u(ip_start:ip_end,ny-4,1:nz) 
    u(ip_start:ip_end,2,1:nz) = u(ip_start:ip_end,ny-3,1:nz)
    u(ip_start:ip_end,3,1:nz) = u(ip_start:ip_end,ny-2,1:nz)
!        
    u(ip_start:ip_end,ny,1:nz) = u(ip_start:ip_end,5,1:nz)
    u(ip_start:ip_end,ny-1,1:nz) = u(ip_start:ip_end,4,1:nz)
!
    v(ip_start:ip_end,1,1:nz) = v(ip_start:ip_end,ny-4,1:nz) 
    v(ip_start:ip_end,2,1:nz) = v(ip_start:ip_end,ny-3,1:nz)
    v(ip_start:ip_end,3,1:nz) = v(ip_start:ip_end,ny-2,1:nz)
!        
    v(ip_start:ip_end,ny,1:nz) = v(ip_start:ip_end,5,1:nz)
    v(ip_start:ip_end,ny-1,1:nz) = v(ip_start:ip_end,4,1:nz)
!
    w(ip_start:ip_end,1,1:nz) = w(ip_start:ip_end,ny-4,1:nz) 
    w(ip_start:ip_end,2,1:nz) = w(ip_start:ip_end,ny-3,1:nz)
    w(ip_start:ip_end,3,1:nz) = w(ip_start:ip_end,ny-2,1:nz)
!        
    w(ip_start:ip_end,ny,1:nz) = w(ip_start:ip_end,5,1:nz)
    w(ip_start:ip_end,ny-1,1:nz) = w(ip_start:ip_end,4,1:nz)
  endif
!
!----------------------------------------------------------!
!             For parallel run, periodic BCs               !
!              are treated via MPI exchange                !
!	 This also handles communication between	   !
!	      neighbouring internal domains	  	   !
!----------------------------------------------------------!
!
#else
    call exchange_y (u, 103)
!
#ifdef MODEL_3D     
    call exchange_y (v, 104)
#endif
!
    call exchange_y (w, 105)
#endif
!
#endif
!
#endif
!
!----------------------------------------------------------!
!              Periodicity in the vertical                 !
!   Typically not used in atmosphere, but for evaluation   !
!----------------------------------------------------------!
!
  if (bcl(5) == 'per' .or. bcl(6) == 'per') then
    u(ip_start:ip_end,jp_start:jp_end,1) = u(ip_start:ip_end,jp_start:jp_end,nz-4) 
    u(ip_start:ip_end,jp_start:jp_end,2) = u(ip_start:ip_end,jp_start:jp_end,nz-3) 
    u(ip_start:ip_end,jp_start:jp_end,3) = u(ip_start:ip_end,jp_start:jp_end,nz-2) 
    u(ip_start:ip_end,jp_start:jp_end,nz) = u(ip_start:ip_end,jp_start:jp_end,5) 
    u(ip_start:ip_end,jp_start:jp_end,nz-1) = u(ip_start:ip_end,jp_start:jp_end,4) 
!
#ifdef MODEL_3D
    v(ip_start:ip_end,jp_start:jp_end,1) = v(ip_start:ip_end,jp_start:jp_end,nz-4) 
    v(ip_start:ip_end,jp_start:jp_end,2) = v(ip_start:ip_end,jp_start:jp_end,nz-3) 
    v(ip_start:ip_end,jp_start:jp_end,3) = v(ip_start:ip_end,jp_start:jp_end,nz-2) 
    v(ip_start:ip_end,jp_start:jp_end,nz) = v(ip_start:ip_end,jp_start:jp_end,5) 
    v(ip_start:ip_end,jp_start:jp_end,nz-1) = v(ip_start:ip_end,jp_start:jp_end,4) 
#endif
!
    w(ip_start:ip_end,jp_start:jp_end,1) = w(ip_start:ip_end,jp_start:jp_end,nz-4) 
    w(ip_start:ip_end,jp_start:jp_end,2) = w(ip_start:ip_end,jp_start:jp_end,nz-3) 
    w(ip_start:ip_end,jp_start:jp_end,3) = w(ip_start:ip_end,jp_start:jp_end,nz-2) 
    w(ip_start:ip_end,jp_start:jp_end,nz) = w(ip_start:ip_end,jp_start:jp_end,5) 
    w(ip_start:ip_end,jp_start:jp_end,nz-1) = w(ip_start:ip_end,jp_start:jp_end,4) 
  endif
!
!----------------------------------------------------------!
!    	    No flux/symmetry boundary conditions	   !
!----------------------------------------------------------!
!
!  U component along x: symmetry or no normal flow
!
  if (it_start == 4) then
    if (bcl(1) == 'nnf') then
      u(3,jp_start:jp_end,1:nz) = - u(4,jp_start:jp_end,1:nz)
      u(2,jp_start:jp_end,1:nz) = - u(5,jp_start:jp_end,1:nz)
      u(1,jp_start:jp_end,1:nz) = - u(6,jp_start:jp_end,1:nz)
    else if (bcl(1) == 'ope') then
      u(3,jp_start:jp_end,1:nz) = u(4,jp_start:jp_end,1:nz)
      u(2,jp_start:jp_end,1:nz) = u(5,jp_start:jp_end,1:nz)
      u(1,jp_start:jp_end,1:nz) = u(6,jp_start:jp_end,1:nz)
    endif
  endif
!
  if (it_end == nx-2) then
    if (bcl(2) == 'nnf') then
      u(nx,jp_start:jp_end,1:nz) = - u(nx-3,jp_start:jp_end,1:nz)
      u(nx-1,jp_start:jp_end,1:nz) = - u(nx-2,jp_start:jp_end,1:nz)
    else if (bcl(2) == 'ope') then
      u(nx,jp_start:jp_end,1:nz) = u(nx-3,jp_start:jp_end,1:nz)
      u(nx-1,jp_start:jp_end,1:nz) = u(nx-2,jp_start:jp_end,1:nz)
    endif
  endif
!
!  Other components along x: symmetry
!
  if ((bcl(1) == 'ope' .or. bcl(1) == 'nnf') .and. it_start == 3) then
#ifdef MODEL_3D
    v(2,jp_start:jp_end,1:nz) = v(4,jp_start:jp_end,1:nz) 
    v(1,jp_start:jp_end,1:nz) = v(5,jp_start:jp_end,1:nz) 
#endif
!
    w(2,jp_start:jp_end,1:nz) = w(4,jp_start:jp_end,1:nz) 
    w(1,jp_start:jp_end,1:nz) = w(5,jp_start:jp_end,1:nz) 
  endif
!
  if ((bcl(2) == 'ope' .or. bcl(2) == 'nnf') .and. it_end == nx-2) then
#ifdef MODEL_3D
    v(nx,jp_start:jp_end,1:nz) = v(nx-4,jp_start:jp_end,1:nz) 
    v(nx-1,jp_start:jp_end,1:nz) = v(nx-3,jp_start:jp_end,1:nz) 
#endif
!
    w(nx,jp_start:jp_end,1:nz) = w(nx-4,jp_start:jp_end,1:nz) 
    w(nx-1,jp_start:jp_end,1:nz) = w(nx-3,jp_start:jp_end,1:nz) 
  endif
!
!  V component along y: symmetry or no normal flow
!
#ifndef CHANNEL
#ifdef MODEL_3D
  if (jt_start == 4) then
    if (bcl(3) == 'nnf') then
      v(ip_start:ip_end,3,1:nz) = - v(ip_start:ip_end,4,1:nz)
      v(ip_start:ip_end,2,1:nz) = - v(ip_start:ip_end,5,1:nz)
      v(ip_start:ip_end,1,1:nz) = - v(ip_start:ip_end,6,1:nz)
    else if (bcl(3) == 'ope') then
      v(ip_start:ip_end,3,1:nz) = v(ip_start:ip_end,4,1:nz)
      v(ip_start:ip_end,2,1:nz) = v(ip_start:ip_end,5,1:nz)
      v(ip_start:ip_end,1,1:nz) = v(ip_start:ip_end,6,1:nz)
    endif
  endif
!
  if (jt_end == ny-2) then
    if (bcl(4) == 'nnf') then
      v(ip_start:ip_end,ny,1:nz) = - v(ip_start:ip_end,ny-3,1:nz)
      v(ip_start:ip_end,ny-1,1:nz) = - v(ip_start:ip_end,ny-2,1:nz)
    else if (bcl(4) == 'ope') then
      v(ip_start:ip_end,ny,1:nz) = v(ip_start:ip_end,ny-3,1:nz)
      v(ip_start:ip_end,ny-1,1:nz) = v(ip_start:ip_end,ny-2,1:nz)
    endif
  endif
!
  if ((bcl(3) == 'ope' .or. bcl(3) == 'nnf') .and. jt_start == 3) then
    u(ip_start:ip_end,2,1:nz) = u(ip_start:ip_end,4,1:nz) 
    u(ip_start:ip_end,1,1:nz) = u(ip_start:ip_end,5,1:nz) 
!
    w(ip_start:ip_end,2,1:nz) = w(ip_start:ip_end,4,1:nz) 
    w(ip_start:ip_end,1,1:nz) = w(ip_start:ip_end,5,1:nz)         
  endif
!
  if ((bcl(4) == 'ope' .or. bcl(4) == 'nnf') .and. jt_end == ny-2) then
    u(ip_start:ip_end,ny,1:nz) = u(ip_start:ip_end,ny-4,1:nz) 
    u(ip_start:ip_end,ny-1,1:nz) = u(ip_start:ip_end,ny-3,1:nz) 
!
    w(ip_start:ip_end,ny,1:nz) = w(ip_start:ip_end,ny-4,1:nz) 
    w(ip_start:ip_end,ny-1,1:nz) = w(ip_start:ip_end,ny-3,1:nz) 
  endif
#endif
#endif
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating windbc')
!
return
end
!
!  ====================================================
  subroutine lbcg_p4 ( grad )
!  ====================================================
 
  integer :: i
  real, dimension(ip_start:,jp_start:,:,:,:), intent(inout) :: grad
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
  call lbcu_p4 ( grad(:,:,:,1,1), grad(:,:,:,2,1), grad(:,:,:,3,1) )
!
#ifdef MODEL_3D
  call lbcu_p4 ( grad(:,:,:,1,2), grad(:,:,:,2,2), grad(:,:,:,3,2) )
#endif
!
  call lbcu_p4 ( grad(:,:,:,1,3), grad(:,:,:,2,3), grad(:,:,:,3,3) )
!
!----------------------------------------------------------!
!           BC for cell centered diagonal terms   	   !
!----------------------------------------------------------!
!
  do i = 1, 3
    call gnbc_p4 ( grad(:,:,:,i,i) )
  enddo
!
!-----------------------------------------------------------!
!
return
end
!
!  ====================================================
  subroutine qnbc_p4 ( x )
!  ====================================================
 
!
! --- Subroutine for calculating periodic-type of
! --- 	lateral boundary values for scalars.
! ---
! --- Last Revision:  Aug 4, 2011
!
!
!	Periodic boundary defined as:
!
!		     L.B.
!	+-----|-----+-----|-----+
!      u(1)  f(1)  u(2)  f(2) u(3)
!   u(nx-1) f(n-1) u(nx) f(n) u(nx+1) 	
!             	    |
!	   forced by average value		
!
!	-------------------------------------------------
!
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: x
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
!  Periodicity in x
!
#ifndef SPMD
  if (.not.nest_run .and. (bcl(1) == 'per' .or. bcl(2) == 'per')) then
    x(1,jp_start:jp_end,1:nz) = x(nx-4,jp_start:jp_end,1:nz)
    x(2,jp_start:jp_end,1:nz) = x(nx-3,jp_start:jp_end,1:nz)
    x(3,jp_start:jp_end,1:nz) = x(nx-2,jp_start:jp_end,1:nz)
!
    x(nx,jp_start:jp_end,1:nz) = x(5,jp_start:jp_end,1:nz)
    x(nx-1,jp_start:jp_end,1:nz) = x(4,jp_start:jp_end,1:nz)  
  endif
!
!  For parallel run, periodic BCs are treated via MPI exchange
!
#else
    call exchange_x (x, 500)
#endif
!
!  Periodicity in y
!
#ifndef CHANNEL
!
#ifdef MODEL_3D
!
#ifndef DECOMP_2D
  if (.not.nest_run .and. (bcl(3) == 'per' .or. bcl(4) == 'per')) then
    x(ip_start:ip_end,1,1:nz) = x(ip_start:ip_end,ny-4,1:nz)
    x(ip_start:ip_end,2,1:nz) = x(ip_start:ip_end,ny-3,1:nz)
    x(ip_start:ip_end,3,1:nz) = x(ip_start:ip_end,ny-2,1:nz)

    x(ip_start:ip_end,ny,1:nz) = x(ip_start:ip_end,5,1:nz)
    x(ip_start:ip_end,ny-1,1:nz) = x(ip_start:ip_end,4,1:nz)
  endif
!
!  For parallel run, periodic BCs are treated via MPI exchange
!
#else
  call exchange_y (x, 501)
#endif
!
#endif
!
#endif
!
!----------------------------------------------------------!
!              Periodicity in the vertical                 !
!   Typically not used in atmosphere, but for evaluation   !
!----------------------------------------------------------!
!
  if (bcl(5) == 'per' .or. bcl(6) == 'per') then
    x(ip_start:ip_end,jp_start:jp_end,1) = x(ip_start:ip_end,jp_start:jp_end,nz-4)
    x(ip_start:ip_end,jp_start:jp_end,2) = x(ip_start:ip_end,jp_start:jp_end,nz-3)
    x(ip_start:ip_end,jp_start:jp_end,3) = x(ip_start:ip_end,jp_start:jp_end,nz-2)
!
    x(ip_start:ip_end,jp_start:jp_end,nz) = x(ip_start:ip_end,jp_start:jp_end,5)
    x(ip_start:ip_end,jp_start:jp_end,nz-1) = x(ip_start:ip_end,jp_start:jp_end,4)
  endif
!
!----------------------------------------------------------!
!   	    No flux/symmetry boundary condition            !
!----------------------------------------------------------!
!
  if (it_start == 4) then
    if (bcl(1) == 'ope' .or. bcl(1) == 'nnf') then
      x(2,jp_start:jp_end,1:nz) = x(4,jp_start:jp_end,1:nz)
      x(1,jp_start:jp_end,1:nz) = x(5,jp_start:jp_end,1:nz)
    endif
  endif  
!
  if (it_end == nx-2) then
    if (bcl(2) == 'ope' .or. bcl(2) == 'nnf') then
      x(nx,jp_start:jp_end,1:nz) = x(nx-4,jp_start:jp_end,1:nz)
      x(nx-1,jp_start:jp_end,1:nz) = x(nx-3,jp_start:jp_end,1:nz)
    endif
  endif  
!
#ifndef CHANNEL
#ifdef MODEL_3D
  if (jt_start == 4) then
    if (bcl(3) == 'ope' .or. bcl(3) == 'nnf') then
      x(ip_start:ip_end,2,1:nz) = x(ip_start:ip_end,4,1:nz)
      x(ip_start:ip_end,1,1:nz) = x(ip_start:ip_end,5,1:nz)
    endif
  endif  
!
  if (jt_end == ny-2) then
    if (bcl(4) == 'ope' .or. bcl(4) == 'nnf') then
      x(ip_start:ip_end,ny,1:nz) = x(ip_start:ip_end,ny-4,1:nz)
      x(ip_start:ip_end,ny-1,1:nz) = x(ip_start:ip_end,ny-3,1:nz)
    endif
  endif  
#endif
#endif
!
!-----------------------------------------------------------!
!
return
end
!
!  ====================================================
  subroutine gnbc_p4 ( x )
!  ====================================================
 
!
! --- Subroutine for calculating periodic-type of
! --- 	lateral boundary values for scalars.
! ---
! --- Last Revision:  Aug 4, 2011
!
!
!	Periodic boundary defined as:
!
!		     L.B.
!	+-----|-----+-----|-----+
!      u(1)  f(1)  u(2)  f(2) u(3)
!   u(nx-1) f(n-1) u(nx) f(n) u(nx+1) 	
!             	    |
!	   forced by average value		
!
!	-------------------------------------------------
!
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: x
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
#ifndef SPMD
  if (.not.nest_run .and. (bcl(1) == 'per' .or. bcl(2) == 'per')) then
    x(1,jp_start:jp_end,1:nz) = x(nx-4,jp_start:jp_end,1:nz)
    x(2,jp_start:jp_end,1:nz) = x(nx-3,jp_start:jp_end,1:nz)
    x(3,jp_start:jp_end,1:nz) = x(nx-2,jp_start:jp_end,1:nz)
!
    x(nx,jp_start:jp_end,1:nz) = x(5,jp_start:jp_end,1:nz)
    x(nx-1,jp_start:jp_end,1:nz) = x(4,jp_start:jp_end,1:nz)  
  endif
!
!  For parallel run, periodic BCs are treated via MPI exchange
!
#else
    call exchange_x (x, 500)
#endif
!
!  Periodicity in y
!
#ifndef CHANNEL
!
#ifdef MODEL_3D
!
#ifndef DECOMP_2D
  if (.not.nest_run .and. (bcl(3) == 'per' .or. bcl(4) == 'per')) then
    x(ip_start:ip_end,1,1:nz) = x(ip_start:ip_end,ny-4,1:nz)
    x(ip_start:ip_end,2,1:nz) = x(ip_start:ip_end,ny-3,1:nz)
    x(ip_start:ip_end,3,1:nz) = x(ip_start:ip_end,ny-2,1:nz)

    x(ip_start:ip_end,ny,1:nz) = x(ip_start:ip_end,5,1:nz)
    x(ip_start:ip_end,ny-1,1:nz) = x(ip_start:ip_end,4,1:nz)
  endif
!
!  For parallel run, periodic BCs are treated via MPI exchange
!
#else
  call exchange_y (x, 501)
#endif
!
#endif
!
#endif
!
!----------------------------------------------------------!
!              Periodicity in the vertical                 !
!   Typically not used in atmosphere, but for evaluation   !
!----------------------------------------------------------!
!
  if (bcl(5) == 'per' .or. bcl(6) == 'per') then
    x(ip_start:ip_end,jp_start:jp_end,1) = x(ip_start:ip_end,jp_start:jp_end,nz-4)
    x(ip_start:ip_end,jp_start:jp_end,2) = x(ip_start:ip_end,jp_start:jp_end,nz-3)
    x(ip_start:ip_end,jp_start:jp_end,3) = x(ip_start:ip_end,jp_start:jp_end,nz-2)
!
    x(ip_start:ip_end,jp_start:jp_end,nz) = x(ip_start:ip_end,jp_start:jp_end,5)
    x(ip_start:ip_end,jp_start:jp_end,nz-1) = x(ip_start:ip_end,jp_start:jp_end,4)
  endif
!
!----------------------------------------------------------!
!	    No flux/symmetry boundary condition            !
!----------------------------------------------------------!
!
  if (it_start == 4) then
    if (bcl(1) == 'ope' .or. bcl(1) == 'nnf') then
      x(2,jp_start:jp_end,1:nz) = -x(4,jp_start:jp_end,1:nz)
      x(1,jp_start:jp_end,1:nz) = -x(5,jp_start:jp_end,1:nz)
    endif
  endif  
!
  if (it_end == nx-2) then
    if (bcl(2) == 'ope' .or. bcl(2) == 'nnf') then
      x(nx,jp_start:jp_end,1:nz) = -x(nx-4,jp_start:jp_end,1:nz)
      x(nx-1,jp_start:jp_end,1:nz) = -x(nx-3,jp_start:jp_end,1:nz)
    endif
  endif  
!
#ifndef CHANNEL
#ifdef MODEL_3D
  if (jt_start == 4) then
    if (bcl(3) == 'ope' .or. bcl(3) == 'nnf') then
      x(ip_start:ip_end,2,1:nz) = -x(ip_start:ip_end,4,1:nz)
      x(ip_start:ip_end,1,1:nz) = -x(ip_start:ip_end,5,1:nz)
    endif
  endif  
!
  if (jt_end == ny-2) then
    if (bcl(4) == 'ope' .or. bcl(4) == 'nnf') then
      x(ip_start:ip_end,ny,1:nz) = -x(ip_start:ip_end,ny-4,1:nz)
      x(ip_start:ip_end,ny-1,1:nz) = -x(ip_start:ip_end,ny-3,1:nz)
    endif
  endif  
#endif
#endif
!
!-----------------------------------------------------------!
!
return
end
!
!  ====================================================
  subroutine windbc_p4 ( x )
!  ====================================================
!
  type (atm_winds), intent(inout) :: x
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
  call lbcu_p4 ( x%u, x%v, x%w )
!
!-----------------------------------------------------------!
!
return
end
!
!  ====================================================
  subroutine statebc_p4 ( x )
!  ====================================================
!
  integer :: h
  type (atm_state), intent(inout) :: x
!
  if (verbose > 1) call write_debug('Starting statebc')
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
  call qnbc_p4 ( x%qt )
!
  call qnbc_p4 ( x%es )
!
  do h = 1, nscal
    call qnbc_p4 ( x%scal(:,:,:,h) )
  enddo
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating statebc')
!
return
end
!
!  ====================================================
  subroutine pressurebc_p4 ( x )
!  ====================================================
!
  type (atm_pressure), intent(inout) :: x
!
  if (verbose > 1) call write_debug('Starting pressurebc')
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
  call qnbc_p4 ( x%p )
!
  call qnbc_p4 ( x%dens )
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating pressurebc')
!
return
end
!
!  ====================================================
  subroutine hydrobc_p4 ( x )
!  ====================================================
!
  integer :: h
  type (hydrometeor), dimension(1:nhydro), intent(inout) :: x
!
  if (verbose > 1) call write_debug('Starting hydrobc')
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
  do h = 1, nhydro
    call qnbc_p4 ( x(h)%q )

#ifdef SEIFERT  
    if (moments == 2) call qnbc_p4 ( x(h)%n )
#endif

    if ( lmicro > 3 .and. h > ice ) call qnbc_p4 ( x(h)%w )
  enddo
!
!-----------------------------------------------------------!
!
  if (verbose > 1) call write_debug('Terminating hydrobc')
!
return
end
!
!  ====================================================
  subroutine aerobc_p4 ( x )
!  ====================================================
!
  integer :: h
  type (aero_3d), dimension(1:nmode), intent(inout) :: x
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
  do h = 1, nmode
    call qnbc_p4 ( x(h)%n )

    call qnbc_p4 ( x(h)%m )

    call qnbc_p4 ( x(h)%ma )
  enddo
!
!-----------------------------------------------------------!
!
return
end
!
!
  subroutine fluxbcu_x ( i1, i2, ut, u, Flux )
!
! --- Subroutine for calculating velocity no-flux BC
! --- Used in advection in the case of non-periodic BC		
!
!	-------------------------------------------------
!
  integer  :: i, i1, i2
!
  real, dimension (ip_start:ip_end) :: ut, u, Flux
!
!----------------------------------------------------------!
!            No flux lateral boundary conditions           !
!----------------------------------------------------------!
!
  if (i1-1 == 3) then
    i = i1
    if (bcl(1) == 'ope') then
      Flux(i) = ut(i+1)*u(i+1)
      Flux(i-1) = Flux(i+1)
    else if (bcl(1) == 'nnf') then
      Flux(i) = 0.0
      Flux(i-1) = -Flux(i+1)
    endif
  endif
!
  if (i2 == nx-2) then
    i = i2
    if (bcl(2) == 'ope') then
      Flux(i) = ut(i)*u(i)
      Flux(i+1) = Flux(i-1)
    else if (bcl(2) == 'nnf') then
      Flux(i) = 0.0
      Flux(i+1) = -Flux(i-1)
    endif
  endif
!
!-----------------------------------------------------------!
!
return
end
!
!
  subroutine fluxbcu_y ( j1, j2, u, Flux )
!
! --- Subroutine for calculating velocity no-flux BC
! --- Used in advection in the case of non-periodic BC		
!
!	-------------------------------------------------
!
  integer  :: j, j1, j2
  real, dimension (jp_start:jp_end) :: u, Flux
!
!----------------------------------------------------------!
!            No flux lateral boundary conditions           !
!----------------------------------------------------------!
!
  if (j1-1 == 3) then
    j = j1
    if (bcl(3) == 'ope') then
      Flux(j-1) = Flux(j)
    else if (bcl(3) == 'nnf') then
      Flux(j-1) = -Flux(j)
    endif
  endif
!
  if (j2 == ny-2) then
    j = min(j2-1,ny-3)
    if (bcl(4) == 'ope') then
      Flux(j+1) = Flux(j)
    else if (bcl(4) == 'nnf') then
      Flux(j+1) = -Flux(j)
    endif
  endif
!
!-----------------------------------------------------------!
!
return
end
!
!
  subroutine fluxbcv_x ( i1, i2, v, Flux )
!
! --- Subroutine for calculating velocity no-flux BC
! --- Used in advection in the case of non-periodic BC		
!
!	-------------------------------------------------
!
  integer  :: i, i1, i2
  real, dimension (ip_start:ip_end) :: v, Flux
!
!----------------------------------------------------------!
!            No flux lateral boundary conditions           !
!----------------------------------------------------------!
!
  if (i1-1 == 3) then
    i = i1
    if (bcl(1) == 'ope') then
      Flux(i-1) = Flux(i)
    else if (bcl(1) == 'nnf') then
      Flux(i-1) = -Flux(i)    
    endif
  endif
!
  if (i2 == nx-2) then
    i = min(i2-1,nx-3)
    if (bcl(2) == 'ope') then
      Flux(i+1) = Flux(i)
    else if (bcl(2) == 'nnf') then
      Flux(i+1) = -Flux(i)
    endif
  endif
!
!-----------------------------------------------------------!
!
return
end
!
!
  subroutine fluxbcv_y ( j1, j2, vt, v, Flux )
!
! --- Subroutine for calculating velocity no-flux BC
! --- Used in advection in the case of non-periodic BC		
!
!	-------------------------------------------------
!
  integer  :: j, j1, j2
  real, dimension (jp_start:jp_end) :: vt, v, Flux
!
!----------------------------------------------------------!
!            No flux lateral boundary conditions           !
!----------------------------------------------------------!
!
  if (j1-1 == 3) then
    j = j1
    if (bcl(3) == 'ope') then
       Flux(j) = vt(j+1)*v(j+1)
       Flux(j-1) = Flux(j+1) 
    else if (bcl(4) == 'nnf') then
       Flux(j) = 0.0
       Flux(j-1) = -Flux(j+1)  
    endif
  endif
!
  if (j2 == ny-2) then
    j = min(j2,ny-2)
    if (bcl(4) == 'ope') then
      Flux(j) = vt(j)*v(j)
      Flux(j+1) = Flux(j-1) 
    else if (bcl(4) == 'nnf') then
      Flux(j) = 0.0
      Flux(j+1) = -Flux(j-1) 
    endif
  endif
!
!-----------------------------------------------------------!
!
return
end
!
!
  subroutine fluxbcw_x ( i1, i2, w, Flux )
!
! --- Subroutine for calculating velocity no-flux BC
! --- Used in advection in the case of non-periodic BC		
!
!	-------------------------------------------------
!
  integer  :: i, i1, i2
  real, dimension (ip_start:ip_end) :: w, Flux
!
!----------------------------------------------------------!
!            No flux lateral boundary conditions           !
!----------------------------------------------------------!
!
  if (i1-1 == 3) then
    i = i1
    if (bcl(1) == 'ope') then
      Flux(i-1) = Flux(i)
    else if (bcl(1) == 'nnf') then
      Flux(i-1) = -Flux(i)
    endif
  endif
!
  if (i2 == nx-2) then
    i = min(i2-1,nx-3)
    if (bcl(2) == 'ope') then
      Flux(i+1) = Flux(i)
    else if (bcl(2) == 'nnf') then
      Flux(i+1) = -Flux(i)
    endif
  endif
!
!-----------------------------------------------------------!
!
return
end
!
!
  subroutine fluxbcw_y ( j1, j2, w, Flux )
!
! --- Subroutine for calculating velocity no-flux BC
! --- Used in advection in the case of non-periodic BC		
!
!	-------------------------------------------------
!
  integer  :: j, j1, j2
  real, dimension (jp_start:jp_end) :: w, Flux
!
!----------------------------------------------------------!
!            No flux lateral boundary conditions           !
!----------------------------------------------------------!
!
  if (j1-1 == 3) then
    j = j1
    if (bcl(3) == 'ope') then
      Flux(j-1) = Flux(j)
    else if (bcl(3) == 'nnf') then
      Flux(j-1) = -Flux(j)
    endif
  endif
!
  if (j2 == ny-2) then
    j = min(j2-1,ny-3)
    if (bcl(4) == 'ope') then
      Flux(j+1) = Flux(j)
    else if (bcl(4) == 'nnf') then
      Flux(j+1) = -Flux(j)
    endif
  endif
!
!-----------------------------------------------------------!
!
return
end

end module boundary_conditions
