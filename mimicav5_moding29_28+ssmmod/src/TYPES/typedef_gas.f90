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

	module typedef_gas

	type :: gas_chemical
		! ---		
		! --- Mixing ratios of gaseous chemical species
		! ---	Units are in ppb(m)
		! ---
		real :: o3	! ozone
		real :: co	! CO
		real :: zco2	! CO2
		real :: ho	! OH
		real :: ho2	! HO2
		real :: h2o2	! H2O2
		real :: xno	! NO
		real :: xno2	! NO2
		real :: xno3	! NO3
		real :: xn2o5	! N2O5
		real :: hno3	! HNO3
     		real :: ch4	! CH4
		real :: ch2o	! CH2O
		real :: ch3o2h	! CH3O2h
		real :: so2	! SO2
		real :: h2so4	! H2SO4
		real :: dms	! CH3SCH3
                real :: nh3     ! NH3
	end type gas_chemical
	
	interface assignment (=)
		module procedure const_to_gas
	end interface

	interface operator (+)
		module procedure gas_add
	end interface

	interface operator (-)
		module procedure gas_sub
	end interface

	interface operator (*)
		module procedure gas_times_real
		module procedure real_times_gas
		module procedure gas_times_int
		module procedure int_times_gas
	end interface

	interface operator (/)
		module procedure gas_div_real
		module procedure gas_div_int
	end interface

	CONTAINS
	
	  subroutine const_to_gas (gas, const)
	  	type (gas_chemical), intent(inout) :: gas
	  	real, intent(in)	  	   :: const
		
		gas%o3     = const
		gas%co     = const
		gas%zco2   = const
		gas%ho     = const
		gas%ho2    = const
		gas%h2o2   = const
		gas%xno    = const
		gas%xno2   = const
		gas%xno3   = const
		gas%xn2o5  = const
		gas%hno3   = const
		gas%ch4    = const
		gas%ch2o   = const
		gas%ch3o2h = const
		gas%so2    = const
		gas%h2so4  = const
		gas%dms    = const
                gas%nh3    = const
	  return
	  end subroutine const_to_gas

	  function gas_add (data1, data2)
	  	type (gas_chemical) :: gas_add
	  	type (gas_chemical), intent(in) :: data1, data2
		
		gas_add%o3     = data1%o3     + data2%o3
		gas_add%co     = data1%co     + data2%co
		gas_add%zco2   = data1%zco2   + data2%zco2
		gas_add%ho     = data1%ho     + data2%ho
		gas_add%ho2    = data1%ho2    + data2%ho2
		gas_add%h2o2   = data1%h2o2   + data2%h2o2
		gas_add%xno    = data1%xno    + data2%xno
		gas_add%xno2   = data1%xno2   + data2%xno2
		gas_add%xno3   = data1%xno3   + data2%xno3
		gas_add%xn2o5  = data1%xn2o5  + data2%xn2o5
		gas_add%hno3   = data1%hno3   + data2%hno3
		gas_add%ch4    = data1%ch4    + data2%ch4
		gas_add%ch2o   = data1%ch2o   + data2%ch2o
		gas_add%ch3o2h = data1%ch3o2h + data2%ch3o2h
		gas_add%so2    = data1%so2    + data2%so2
		gas_add%h2so4  = data1%h2so4  + data2%h2so4
		gas_add%dms    = data1%dms    + data2%dms
		gas_add%nh3    = data1%nh3    + data2%nh3
	  end function gas_add
	  
	  function gas_sub (data1, data2)
	  	type (gas_chemical) :: gas_sub
	  	type (gas_chemical), intent(in) :: data1, data2
		
		gas_sub%o3     = data1%o3     - data2%o3
		gas_sub%co     = data1%co     - data2%co
		gas_sub%zco2   = data1%zco2   - data2%zco2
		gas_sub%ho     = data1%ho     - data2%ho
		gas_sub%ho2    = data1%ho2    - data2%ho2
		gas_sub%h2o2   = data1%h2o2   - data2%h2o2
		gas_sub%xno    = data1%xno    - data2%xno
		gas_sub%xno2   = data1%xno2   - data2%xno2
		gas_sub%xno3   = data1%xno3   - data2%xno3
		gas_sub%xn2o5  = data1%xn2o5  - data2%xn2o5
		gas_sub%hno3   = data1%hno3   - data2%hno3
		gas_sub%ch4    = data1%ch4    - data2%ch4
		gas_sub%ch2o   = data1%ch2o   - data2%ch2o
		gas_sub%ch3o2h = data1%ch3o2h - data2%ch3o2h
		gas_sub%so2    = data1%so2    - data2%so2
		gas_sub%h2so4  = data1%h2so4  - data2%h2so4
		gas_sub%dms    = data1%dms    - data2%dms
		gas_sub%nh3    = data1%nh3    - data2%nh3
	  end function gas_sub
	  
	  function gas_times_real (data1, const)
	  	type (gas_chemical) :: gas_times_real
	  	type (gas_chemical), intent(in) :: data1
		real, intent(in)                :: const
		
		gas_times_real%o3     = data1%o3     *const
		gas_times_real%co     = data1%co     *const
		gas_times_real%zco2   = data1%zco2   *const
		gas_times_real%ho     = data1%ho     *const
		gas_times_real%ho2    = data1%ho2    *const
		gas_times_real%h2o2   = data1%h2o2   *const
		gas_times_real%xno    = data1%xno    *const
		gas_times_real%xno2   = data1%xno2   *const
		gas_times_real%xno3   = data1%xno3   *const
		gas_times_real%xn2o5  = data1%xn2o5  *const
		gas_times_real%hno3   = data1%hno3   *const
		gas_times_real%ch4    = data1%ch4    *const
		gas_times_real%ch2o   = data1%ch2o   *const
		gas_times_real%ch3o2h = data1%ch3o2h *const
		gas_times_real%so2    = data1%so2    *const
		gas_times_real%h2so4  = data1%h2so4  *const
		gas_times_real%dms    = data1%dms    *const
		gas_times_real%nh3    = data1%nh3    *const
	  end function gas_times_real

	  function real_times_gas (const, data1)
	  	type (gas_chemical) :: real_times_gas
	  	type (gas_chemical), intent(in) :: data1
		real, intent(in)                :: const
		
		real_times_gas%o3     = const *data1%o3
		real_times_gas%co     = const *data1%co
		real_times_gas%zco2   = const *data1%zco2
		real_times_gas%ho     = const *data1%ho
		real_times_gas%ho2    = const *data1%ho2
		real_times_gas%h2o2   = const *data1%h2o2
		real_times_gas%xno    = const *data1%xno
		real_times_gas%xno2   = const *data1%xno2
		real_times_gas%xno3   = const *data1%xno3
		real_times_gas%xn2o5  = const *data1%xn2o5
		real_times_gas%hno3   = const *data1%hno3
		real_times_gas%ch4    = const *data1%ch4
		real_times_gas%ch2o   = const *data1%ch2o
		real_times_gas%ch3o2h = const *data1%ch3o2h
		real_times_gas%so2    = const *data1%so2
		real_times_gas%h2so4  = const *data1%h2so4
		real_times_gas%dms    = const *data1%dms
		real_times_gas%nh3    = const *data1%nh3
	  end function real_times_gas
	  
	  function gas_times_int (data1, const)
	  	type (gas_chemical) :: gas_times_int
	  	type (gas_chemical), intent(in) :: data1
		integer, intent(in)             :: const
		
		gas_times_int%o3     = data1%o3     *real(const)
		gas_times_int%co     = data1%co     *real(const)
		gas_times_int%zco2   = data1%zco2   *real(const)
		gas_times_int%ho     = data1%ho     *real(const)
		gas_times_int%ho2    = data1%ho2    *real(const)
		gas_times_int%h2o2   = data1%h2o2   *real(const)
		gas_times_int%xno    = data1%xno    *real(const)
		gas_times_int%xno2   = data1%xno2   *real(const)
		gas_times_int%xno3   = data1%xno3   *real(const)
		gas_times_int%xn2o5  = data1%xn2o5  *real(const)
		gas_times_int%hno3   = data1%hno3   *real(const)
		gas_times_int%ch4    = data1%ch4    *real(const)
		gas_times_int%ch2o   = data1%ch2o   *real(const)
		gas_times_int%ch3o2h = data1%ch3o2h *real(const)
		gas_times_int%so2    = data1%so2    *real(const)
		gas_times_int%h2so4  = data1%h2so4  *real(const)
		gas_times_int%dms    = data1%dms    *real(const)
		gas_times_int%nh3    = data1%nh3    *real(const)
	  end function gas_times_int

	  function int_times_gas (const, data1)
	  	type (gas_chemical) :: int_times_gas
	  	type (gas_chemical), intent(in) :: data1
		integer, intent(in)             :: const
		
		int_times_gas%o3     = real(const) *data1%o3
		int_times_gas%co     = real(const) *data1%co
		int_times_gas%zco2   = real(const) *data1%zco2
		int_times_gas%ho     = real(const) *data1%ho
		int_times_gas%ho2    = real(const) *data1%ho2
		int_times_gas%h2o2   = real(const) *data1%h2o2
		int_times_gas%xno    = real(const) *data1%xno
		int_times_gas%xno2   = real(const) *data1%xno2
		int_times_gas%xno3   = real(const) *data1%xno3
		int_times_gas%xn2o5  = real(const) *data1%xn2o5
		int_times_gas%hno3   = real(const) *data1%hno3
		int_times_gas%ch4    = real(const) *data1%ch4
		int_times_gas%ch2o   = real(const) *data1%ch2o
		int_times_gas%ch3o2h = real(const) *data1%ch3o2h
		int_times_gas%so2    = real(const) *data1%so2
		int_times_gas%h2so4  = real(const) *data1%h2so4
		int_times_gas%dms    = real(const) *data1%dms
		int_times_gas%nh3    = real(const) *data1%nh3
	  end function int_times_gas	  		

	  function gas_div_real (data1, const)
	  	type (gas_chemical) :: gas_div_real
	  	type (gas_chemical), intent(in) :: data1
		real, intent(in)                :: const
		
		gas_div_real%o3     = data1%o3     /const
		gas_div_real%co     = data1%co     /const
		gas_div_real%zco2   = data1%zco2   /const
		gas_div_real%ho     = data1%ho     /const
		gas_div_real%ho2    = data1%ho2    /const
		gas_div_real%h2o2   = data1%h2o2   /const
		gas_div_real%xno    = data1%xno    /const
		gas_div_real%xno2   = data1%xno2   /const
		gas_div_real%xno3   = data1%xno3   /const
		gas_div_real%xn2o5  = data1%xn2o5  /const
		gas_div_real%hno3   = data1%hno3   /const
		gas_div_real%ch4    = data1%ch4    /const
		gas_div_real%ch2o   = data1%ch2o   /const
		gas_div_real%ch3o2h = data1%ch3o2h /const
		gas_div_real%so2    = data1%so2    /const
		gas_div_real%h2so4  = data1%h2so4  /const
		gas_div_real%dms    = data1%dms    /const
		gas_div_real%nh3    = data1%nh3    /const
	  end function gas_div_real

	  function gas_div_int (data1, const)
	  	type (gas_chemical) :: gas_div_int
	  	type (gas_chemical), intent(in) :: data1
		integer, intent(in)             :: const
		
		gas_div_int%o3     = data1%o3	  /real(const)
		gas_div_int%co     = data1%co	  /real(const)
		gas_div_int%zco2   = data1%zco2   /real(const)
		gas_div_int%ho     = data1%ho	  /real(const)
		gas_div_int%ho2    = data1%ho2    /real(const)
		gas_div_int%h2o2   = data1%h2o2   /real(const)
		gas_div_int%xno    = data1%xno    /real(const)
		gas_div_int%xno2   = data1%xno2   /real(const)
		gas_div_int%xno3   = data1%xno3   /real(const)
		gas_div_int%xn2o5  = data1%xn2o5  /real(const)
		gas_div_int%hno3   = data1%hno3   /real(const)
		gas_div_int%ch4    = data1%ch4    /real(const)
		gas_div_int%ch2o   = data1%ch2o   /real(const)
		gas_div_int%ch3o2h = data1%ch3o2h /real(const)
		gas_div_int%so2    = data1%so2    /real(const)
		gas_div_int%h2so4  = data1%h2so4  /real(const)
		gas_div_int%dms    = data1%dms    /real(const)
		gas_div_int%nh3    = data1%nh3    /real(const)
	  end function gas_div_int
	  	  	    			
	end module typedef_gas
