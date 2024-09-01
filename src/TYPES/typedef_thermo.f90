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

	module typedef_thermo
	
	USE gridno
	USE shared_data
	
	type :: thermodyn
		! ---
		real, dimension(:,:,:), allocatable :: T   	! temperature
		real, dimension(:,:,:), allocatable :: exn      ! Exner function
!
		real, dimension(:,:,:), allocatable :: buoy     ! buoyancy 
		real, dimension(:,:,:), allocatable :: ptv  	! Equivalent Potential temperature
		real, dimension(:,:,:), allocatable :: mse	! Moist static energy 
		real, dimension(:,:,:), allocatable :: pt  	! Potential temperature
		real, dimension(:,:,:), allocatable :: qv  	! Water vapor mass mixing ratio
		real, dimension(:,:,:), allocatable :: rh, rhi  ! Relative humidities
	end type thermodyn
	
	type :: thermoprop
		! ---
		real, dimension(:,:,:), allocatable :: cs2	! sound speed
		real, dimension(:,:,:), allocatable :: cp	! heat capacity cst pressure
		real, dimension(:,:,:), allocatable :: cv	! heat capacity cst volume		
	end type thermoprop

	interface assignment (=)
		module procedure const_to_thermo
	end interface

	interface operator (+)
		module procedure thermo_add
	end interface

	interface operator (-)
		module procedure thermo_sub
	end interface

	interface operator (*)
		module procedure thermo_times_real
		module procedure real_times_thermo
		module procedure thermo_times_int
		module procedure int_times_thermo
	end interface

	interface operator (/)
		module procedure thermo_div_real
		module procedure thermo_div_int
	end interface
	
	interface alloc
		module procedure alloc_thermo, alloc_thermoprop
	end interface
	
	interface dealloc
		module procedure dealloc_thermo, dealloc_thermoprop
	end interface

	CONTAINS
	
	  subroutine alloc_thermo (thermol)
	  	type (thermodyn), intent(inout) :: thermol
		
		allocate( thermol%T(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate( thermol%exn(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate( thermol%buoy(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate( thermol%pt(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate( thermol%qv(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		if (out_ptv) allocate( thermol%ptv(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_mse) allocate( thermol%mse(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_sat) allocate ( thermol%rh(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        if (out_sat) allocate ( thermol%rhi(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		thermol%T = 0.
		thermol%exn = 0.
		thermol%buoy = 0.
		thermol%pt = 0.
		thermol%qv = 0.
		
		if (out_ptv) thermol%ptv = 0.
		if (out_mse) thermol%mse = 0.
		if (out_sat) thermol%rh = 0.
		if (out_sat) thermol%rhi = 0.
	  return
	  end subroutine alloc_thermo
	
	  subroutine dealloc_thermo (thermol)
	  	type (thermodyn), intent(inout) :: thermol
		
		deallocate( thermol%T )
		deallocate( thermol%exn )
		deallocate( thermol%buoy )
		deallocate( thermol%pt )
		deallocate( thermol%qv )
		
		if (out_ptv) deallocate( thermol%ptv )
		if (out_mse) deallocate( thermol%mse )
	        if (out_sat) deallocate ( thermol%rh )
	        if (out_sat) deallocate ( thermol%rhi )
	  return
	  end subroutine dealloc_thermo
	
	  subroutine alloc_thermoprop (thermol)
	  	type (thermoprop), intent(inout) :: thermol
		
		allocate( thermol%cp(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate( thermol%cv(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		thermol%cp = 0.
		thermol%cv = 0.

#if (!defined ANELASTIC) || (!defined ISENTROPIC)
		allocate( thermol%cs2(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		thermol%cs2 = 0.
#endif
	  return
	  end subroutine alloc_thermoprop
	
	  subroutine dealloc_thermoprop (thermol)
	  	type (thermoprop), intent(inout) :: thermol

		deallocate( thermol%cp )
		deallocate( thermol%cv )
		
#if (!defined ANELASTIC) || (!defined ISENTROPIC)
		deallocate( thermol%cs2 )
#endif
	  return
	  end subroutine dealloc_thermoprop
	
	  subroutine const_to_thermo (thermol, const)
	  	type (thermodyn), intent(inout) :: thermol
	  	real, intent(in)		:: const
		
		thermol%pt     = const
		thermol%mse    = const
		thermol%T      = const
		thermol%buoy   = const
		thermol%ptv    = const
		thermol%qv     = const
		thermol%exn    = const
	  return
	  end subroutine const_to_thermo

	  function thermo_add (data1, data2)
	  	type (thermodyn) :: thermo_add
	  	type (thermodyn), intent(in) :: data1, data2
		
		thermo_add%pt     = data1%pt     + data2%pt
		thermo_add%mse    = data1%mse    + data2%mse
		thermo_add%T      = data1%T      + data2%T
		thermo_add%buoy   = data1%buoy   + data2%buoy
		thermo_add%ptv    = data1%ptv    + data2%ptv
		thermo_add%qv     = data1%qv     + data2%qv
		thermo_add%exn    = data1%exn    + data2%exn
	  end function thermo_add
	  
	  function thermo_sub (data1, data2)
	  	type (thermodyn) :: thermo_sub
	  	type (thermodyn), intent(in) :: data1, data2
		
		thermo_sub%pt     = data1%pt     - data2%pt
		thermo_sub%mse    = data1%mse    - data2%mse
		thermo_sub%T      = data1%T      - data2%T
		thermo_sub%buoy   = data1%buoy   - data2%buoy
		thermo_sub%ptv    = data1%ptv    - data2%ptv
		thermo_sub%qv     = data1%qv     - data2%qv
		thermo_sub%exn    = data1%exn    - data2%exn
	  end function thermo_sub
	  
	  function thermo_times_real (data1, const)
	  	type (thermodyn) :: thermo_times_real
	  	type (thermodyn), intent(in) :: data1
		real, intent(in)	     :: const
		
		thermo_times_real%pt     = data1%pt   *const
		thermo_times_real%mse    = data1%mse  *const
		thermo_times_real%T      = data1%T    *const
		thermo_times_real%buoy   = data1%buoy *const
		thermo_times_real%ptv    = data1%ptv  *const
		thermo_times_real%qv     = data1%qv   *const
		thermo_times_real%exn    = data1%exn  *const
	  end function thermo_times_real

	  function real_times_thermo (const, data1)
	  	type (thermodyn) :: real_times_thermo
	  	type (thermodyn), intent(in) :: data1
		real, intent(in)	     :: const
		
		real_times_thermo%pt     = const *data1%pt
		real_times_thermo%mse    = const *data1%mse
		real_times_thermo%T      = const *data1%T
		real_times_thermo%buoy   = const *data1%buoy
		real_times_thermo%ptv    = const *data1%ptv
		real_times_thermo%qv     = const *data1%qv
		real_times_thermo%exn    = const *data1%exn
	  end function real_times_thermo
	  
	  function thermo_times_int (data1, const)
	  	type (thermodyn) :: thermo_times_int
	  	type (thermodyn), intent(in) :: data1
		integer, intent(in)	     :: const
		
		thermo_times_int%pt     = data1%pt   *real(const)
		thermo_times_int%mse    = data1%mse  *real(const)
		thermo_times_int%T      = data1%T    *real(const)
		thermo_times_int%buoy   = data1%buoy *real(const)
		thermo_times_int%ptv    = data1%ptv  *real(const)
		thermo_times_int%qv     = data1%qv   *real(const)
		thermo_times_int%exn    = data1%exn  *real(const)
	  end function thermo_times_int

	  function int_times_thermo (const, data1)
	  	type (thermodyn) :: int_times_thermo
	  	type (thermodyn), intent(in) :: data1
		integer, intent(in)	     :: const
		
		int_times_thermo%pt     = real(const) *data1%pt
		int_times_thermo%mse    = real(const) *data1%mse
		int_times_thermo%T      = real(const) *data1%T
		int_times_thermo%buoy   = real(const) *data1%buoy
		int_times_thermo%ptv    = real(const) *data1%ptv
		int_times_thermo%qv     = real(const) *data1%qv
		int_times_thermo%exn    = real(const) *data1%exn
	  end function int_times_thermo	  		

	  function thermo_div_real (data1, const)
	  	type (thermodyn) :: thermo_div_real
	  	type (thermodyn), intent(in) :: data1
		real, intent(in)	     :: const
		
		thermo_div_real%pt     = data1%pt   /const
		thermo_div_real%mse    = data1%mse  /const
		thermo_div_real%T      = data1%T    /const
		thermo_div_real%buoy   = data1%buoy /const
		thermo_div_real%ptv    = data1%ptv  /const
		thermo_div_real%qv     = data1%qv   /const 
		thermo_div_real%exn    = data1%exn  /const 
	  end function thermo_div_real

	  function thermo_div_int (data1, const)
	  	type (thermodyn) :: thermo_div_int
	  	type (thermodyn), intent(in) :: data1
		integer, intent(in)	     :: const
		
		thermo_div_int%pt     = data1%pt    /real(const)
		thermo_div_int%mse    = data1%mse   /real(const)
		thermo_div_int%T      = data1%T     /real(const)
		thermo_div_int%buoy   = data1%buoy  /real(const)
		thermo_div_int%ptv    = data1%ptv   /real(const)
		thermo_div_int%qv     = data1%qv    /real(const)
		thermo_div_int%exn    = data1%exn   /real(const)
	  end function thermo_div_int
	  	  	    			
	end module typedef_thermo

