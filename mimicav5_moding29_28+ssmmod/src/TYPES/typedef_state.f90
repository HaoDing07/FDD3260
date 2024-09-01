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
!	Typedef for derived-type data: atm_state			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!		
! ================================================================
	
	module typedef_state
	
	USE gridno
	
	type :: atm_state
		! ---
		real, dimension(:,:,:), allocatable :: es  		! Potential temperature
		real, dimension(:,:,:), allocatable :: qt      		! Total water mixing ratio (kg/kg)
		real, dimension(:,:,:,:), allocatable :: scal		! Additional tracers
	end type atm_state
	
	interface assignment (=)
		module procedure const_to_state
	end interface

	interface operator (+)
		module procedure state_add
	end interface

	interface operator (-)
		module procedure state_sub
	end interface

	interface operator (*)
		module procedure state_times_real
		module procedure real_times_state
		module procedure state_times_int
		module procedure int_times_state
	end interface

	interface operator (/)
		module procedure state_div_real
		module procedure state_div_int
	end interface
	
	interface alloc
		module procedure allocate_state
	end interface
	
	interface dealloc
		module procedure deallocate_state
	end interface

	CONTAINS
	
	  subroutine allocate_state (data)
	  	type(atm_state), intent(inout) :: data
		
		allocate ( data%qt(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate ( data%es(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate ( data%scal(ip_start:ip_end,jp_start:jp_end,1:nz,1:nscal) )
		
		data%qt = 0.
		data%es = 0.
		data%scal = 0.
	  return
	  end subroutine allocate_state
	
	  subroutine deallocate_state (data)
	  	type(atm_state), intent(inout) :: data
		
		deallocate ( data%qt )
		deallocate ( data%es )
		deallocate ( data%scal )
	  return
	  end subroutine deallocate_state
	
	  subroutine const_to_state (state, const)
	  	type (atm_state), intent(inout) :: state
	  	real, intent(in)		:: const
		
		state%es   = const
		state%qt   = const
		state%scal = const
	  return
	  end subroutine const_to_state

	  function state_add (data1, data2)
	  	type (atm_state) :: state_add
	  	type (atm_state), intent(in) :: data1, data2
		
		state_add%es   = data1%es   + data2%es
		state_add%qt   = data1%qt   + data2%qt
		state_add%scal = data1%scal + data2%scal
	  end function state_add
	  
	  function state_sub (data1, data2)
	  	type (atm_state) :: state_sub
	  	type (atm_state), intent(in) :: data1, data2
		
		state_sub%es   = data1%es   - data2%es
		state_sub%qt   = data1%qt   - data2%qt
		state_sub%scal = data1%scal - data2%scal
	  end function state_sub
	  
	  function state_times_real (data1, const)
	  	type (atm_state) :: state_times_real
	  	type (atm_state), intent(in) :: data1
		real, intent(in)	     :: const
		
		state_times_real%es   = data1%es   *const
		state_times_real%qt   = data1%qt   *const
		state_times_real%scal = data1%scal *const
	  end function state_times_real

	  function real_times_state (const, data1)
	  	type (atm_state) :: real_times_state
	  	type (atm_state), intent(in) :: data1
		real, intent(in)	     :: const
		
		real_times_state%es   = const *data1%es
		real_times_state%qt   = const *data1%qt
		real_times_state%scal = const *data1%scal
	  end function real_times_state
	  
	  function state_times_int (data1, const)
	  	type (atm_state) :: state_times_int
	  	type (atm_state), intent(in) :: data1
		integer, intent(in)	     :: const
		
		state_times_int%es   = data1%es   *real(const)
		state_times_int%qt   = data1%qt   *real(const)
		state_times_int%scal = data1%scal *real(const)
	  end function state_times_int

	  function int_times_state (const, data1)
	  	type (atm_state) :: int_times_state
	  	type (atm_state), intent(in) :: data1
		integer, intent(in)	     :: const
		
		int_times_state%es   = real(const) *data1%es
		int_times_state%qt   = real(const) *data1%qt
		int_times_state%scal = real(const) *data1%scal
	  end function int_times_state	  		

	  function state_div_real (data1, const)
	  	type (atm_state) :: state_div_real
	  	type (atm_state), intent(in) :: data1
		real, intent(in)	     :: const
		
		state_div_real%es   = data1%es   /const
		state_div_real%qt   = data1%qt   /const
		state_div_real%scal = data1%scal /const
	  end function state_div_real

	  function state_div_int (data1, const)
	  	type (atm_state) :: state_div_int
	  	type (atm_state), intent(in) :: data1
		integer, intent(in)	     :: const
		
		state_div_int%es   = data1%es   /real(const)
		state_div_int%qt   = data1%qt   /real(const)
		state_div_int%scal = data1%scal /real(const)
	  end function state_div_int
	  	  	    			
	end module typedef_state

