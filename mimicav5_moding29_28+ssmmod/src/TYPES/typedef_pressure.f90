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
!  TYPEDEF_STATE.F:                   
!
!  Purpose:
!	Typedef for derived-type data: atm_pressure			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!		
! ================================================================
	
	module typedef_pressure
	
	USE gridno
	
	type :: atm_pressure
		! ---
		real, dimension(:,:,:), allocatable :: p
		real, dimension(:,:,:), allocatable :: dens
	end type atm_pressure
	
	interface assignment (=)
		module procedure const_to_pressure
	end interface

	interface operator (+)
		module procedure pressure_add
	end interface

	interface operator (-)
		module procedure pressure_sub
	end interface

	interface operator (*)
		module procedure pressure_times_real
		module procedure real_times_pressure
		module procedure pressure_times_int
		module procedure int_times_pressure
	end interface

	interface operator (/)
		module procedure pressure_div_real
		module procedure pressure_div_int
	end interface
	
	interface alloc
		module procedure alloc_pres
	end interface
	
	interface dealloc
		module procedure dealloc_pres
	end interface

	CONTAINS
	
	  subroutine alloc_pres (pres)
	  	type (atm_pressure), intent(inout) :: pres
		
		allocate ( pres%p(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate ( pres%dens(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		pres%p = 0.
		pres%dens = 0.
	  return
	  end subroutine alloc_pres
	
	  subroutine dealloc_pres (pres)
	  	type (atm_pressure), intent(inout) :: pres
		
		deallocate ( pres%p )
		deallocate ( pres%dens )
	  return
	  end subroutine dealloc_pres
	
	  subroutine const_to_pressure (press, const)
	  	type (atm_pressure), intent(inout) :: press
	  	real, intent(in)		:: const
		
		press%p    = const
		press%dens = const
	  return
	  end subroutine const_to_pressure

	  function pressure_add (data1, data2)
	  	type (atm_pressure) :: pressure_add
	  	type (atm_pressure), intent(in) :: data1, data2
		
		pressure_add%p    = data1%p    + data2%p
		pressure_add%dens = data1%dens + data2%dens
	  end function pressure_add
	  
	  function pressure_sub (data1, data2)
	  	type (atm_pressure) :: pressure_sub
	  	type (atm_pressure), intent(in) :: data1, data2
		
		pressure_sub%p    = data1%p    - data2%p
		pressure_sub%dens = data1%dens - data2%dens
	  end function pressure_sub
	  
	  function pressure_times_real (data1, const)
	  	type (atm_pressure) :: pressure_times_real
	  	type (atm_pressure), intent(in) :: data1
		real, intent(in)	     :: const
		
		pressure_times_real%p    = data1%p    *const
		pressure_times_real%dens = data1%dens *const
	  end function pressure_times_real

	  function real_times_pressure (const, data1)
	  	type (atm_pressure) :: real_times_pressure
	  	type (atm_pressure), intent(in) :: data1
		real, intent(in)	     :: const
		
		real_times_pressure%p    = const *data1%p
		real_times_pressure%dens = const *data1%dens
	  end function real_times_pressure
	  
	  function pressure_times_int (data1, const)
	  	type (atm_pressure) :: pressure_times_int
	  	type (atm_pressure), intent(in) :: data1
		integer, intent(in)	     :: const
		
		pressure_times_int%p    = data1%p    *real(const)
		pressure_times_int%dens = data1%dens *real(const)
	  end function pressure_times_int

	  function int_times_pressure (const, data1)
	  	type (atm_pressure) :: int_times_pressure
	  	type (atm_pressure), intent(in) :: data1
		integer, intent(in)	     :: const
		
		int_times_pressure%p    = real(const) *data1%p
		int_times_pressure%dens = real(const) *data1%dens
	  end function int_times_pressure	  		

	  function pressure_div_real (data1, const)
	  	type (atm_pressure) :: pressure_div_real
	  	type (atm_pressure), intent(in) :: data1
		real, intent(in)	     :: const
		
		pressure_div_real%p    = data1%p	 /const
		pressure_div_real%dens = data1%dens /const
	  end function pressure_div_real

	  function pressure_div_int (data1, const)
	  	type (atm_pressure) :: pressure_div_int
	  	type (atm_pressure), intent(in) :: data1
		integer, intent(in)	     :: const
		
		pressure_div_int%p    = data1%p    /real(const)
		pressure_div_int%dens = data1%dens /real(const)
	  end function pressure_div_int
	  	  	    			
	end module typedef_pressure

