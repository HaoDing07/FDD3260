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
!  TYPEDEF_WIND.F:                   
!
!  Purpose:
!	Typedef for derived-type data: winds			  
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================
	
	module typedef_wind
	
	USE gridno

	type :: atm_winds
		! ---
		real, dimension(:,:,:), allocatable :: u	! Velocity in x 
		real, dimension(:,:,:), allocatable :: v	! Velocity in y
		real, dimension(:,:,:), allocatable :: w	! Velocity in z
	end type atm_winds
	
	interface assignment (=)
		module procedure const_to_wind
	end interface

	interface operator (+)
		module procedure wind_add
	end interface

	interface operator (-)
		module procedure wind_sub
	end interface

	interface operator (*)
		module procedure wind_times_real
		module procedure real_times_wind
		module procedure wind_times_int
		module procedure int_times_wind
	end interface

	interface operator (/)
		module procedure wind_div_real
		module procedure wind_div_int
	end interface
	
	interface alloc
		module procedure allocate_winds
	end interface
	
	interface dealloc
		module procedure deallocate_winds
	end interface

	CONTAINS
	
	  subroutine allocate_winds (windl)
	  	type(atm_winds), intent(inout) :: windl
		
		allocate ( windl%u(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate ( windl%v(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate ( windl%w(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		windl%u = 0.
		windl%v = 0.
		windl%w = 0.
	  return
	  end subroutine allocate_winds
	
	  subroutine deallocate_winds (windl)
	  	type(atm_winds), intent(inout) :: windl
		
		deallocate ( windl%u )
		deallocate ( windl%v )
		deallocate ( windl%w )
	  return
	  end subroutine deallocate_winds
	
	  subroutine const_to_wind (wind, const)
	  	type (atm_winds), intent(inout) :: wind
	  	real, intent(in)		:: const
		
		wind%u    = const
#ifdef MODEL_3D
		wind%v    = const
#endif
		wind%w    = const
	  return
	  end subroutine const_to_wind

	  function wind_add (data1, data2)
	  	type (atm_winds) :: wind_add
	  	type (atm_winds), intent(in) :: data1, data2
		
		wind_add%u    = data1%u    + data2%u
#ifdef MODEL_3D
		wind_add%v    = data1%v    + data2%v
#endif
		wind_add%w    = data1%w    + data2%w
	  end function wind_add
	  
	  function wind_sub (data1, data2)
	  	type (atm_winds) :: wind_sub
	  	type (atm_winds), intent(in) :: data1, data2
		
		wind_sub%u    = data1%u    - data2%u
#ifdef MODEL_3D
		wind_sub%v    = data1%v    - data2%v
#endif
		wind_sub%w    = data1%w - data2%w
	  end function wind_sub
	  
	  function wind_times_real (data1, const)
	  	type (atm_winds) :: wind_times_real
	  	type (atm_winds), intent(in) :: data1
		real, intent(in)	     :: const
		
		wind_times_real%u    = data1%u    *const
#ifdef MODEL_3D
		wind_times_real%v    = data1%v    *const
#endif
		wind_times_real%w    = data1%w *const
	  end function wind_times_real

	  function real_times_wind (const, data1)
	  	type (atm_winds) :: real_times_wind
	  	type (atm_winds), intent(in) :: data1
		real, intent(in)	     :: const
		
		real_times_wind%u    = const *data1%u
#ifdef MODEL_3D
		real_times_wind%v    = const *data1%v
#endif
		real_times_wind%w    = const *data1%w
	  end function real_times_wind
	  
	  function wind_times_int (data1, const)
	  	type (atm_winds) :: wind_times_int
	  	type (atm_winds), intent(in) :: data1
		integer, intent(in)	     :: const
		
		wind_times_int%u    = data1%u    *real(const)
		wind_times_int%v    = data1%v    *real(const)
		wind_times_int%w    = data1%w *real(const)
	  end function wind_times_int

	  function int_times_wind (const, data1)
	  	type (atm_winds) :: int_times_wind
	  	type (atm_winds), intent(in) :: data1
		integer, intent(in)	     :: const
		
		int_times_wind%u    = real(const) *data1%u
		int_times_wind%v    = real(const) *data1%v
		int_times_wind%w    = real(const) *data1%w
	  end function int_times_wind	  		

	  function wind_div_real (data1, const)
	  	type (atm_winds) :: wind_div_real
	  	type (atm_winds), intent(in) :: data1
		real, intent(in)	     :: const
		
		wind_div_real%u    = data1%u	 /const
		wind_div_real%v    = data1%v	 /const
		wind_div_real%w    = data1%w /const
	  end function wind_div_real

	  function wind_div_int (data1, const)
	  	type (atm_winds) :: wind_div_int
	  	type (atm_winds), intent(in) :: data1
		integer, intent(in)	     :: const
		
		wind_div_int%u    = data1%u    /real(const)
		wind_div_int%v    = data1%v    /real(const)
		wind_div_int%w    = data1%w /real(const)
	  end function wind_div_int
		  	  	    			
	end module typedef_wind

