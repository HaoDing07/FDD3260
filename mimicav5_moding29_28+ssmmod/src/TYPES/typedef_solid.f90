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

	module typedef_solid
	
	type :: solid_chemical
		! ---
		! --- Mixing ratios of solid chemical species
		! ---	Units are in ppb(m) to air
		!	conversion can be ppb(m)*1.e-9*den/Mw -> mole/Lair
		!		or mole/Lair*Mw/den*1.e9 -> ppb(m)
		! ---
		real :: o3	! ozone
		real :: h2o2	! H2O2
		real :: xnv	! N(V) = HNO3 + N2O5
		real :: ch2o	! CH2O
		real :: ch3o2h	! CH3O2h
		real :: siv	! S(IV)
		real :: svi	! S(VI)
                real :: nh4     ! NH4
	end type solid_chemical

	interface assignment (=)
		module procedure const2solid
	end interface

	interface operator (+)
		module procedure solid_add
	end interface

	interface operator (-)
		module procedure solid_sub
	end interface

	interface operator (*)
		module procedure solid_times_real
		module procedure real_times_solid
		module procedure solid_times_int
		module procedure int_times_solid
	end interface

	interface operator (/)
		module procedure solid_div_real
		module procedure solid_div_int
	end interface

	CONTAINS
	
	  subroutine const2solid (solid, const)
	       type (solid_chemical), intent(inout) :: solid
	       real, intent(in) 	            :: const
	       
	       solid%o3     = const
	       solid%h2o2   = const
	       solid%xnv    = const
	       solid%ch2o   = const
	       solid%ch3o2h = const
	       solid%siv    = const
	       solid%svi    = const
	       solid%nh4   = const
	  return
	  end subroutine const2solid

	  function solid_add (data1, data2)
	  	type (solid_chemical) :: solid_add
	  	type (solid_chemical), intent(in) :: data1, data2
		
		solid_add%o3       = data1%o3	 + data2%o3
		solid_add%h2o2     = data1%h2o2  + data2%h2o2
		solid_add%xnv      = data1%xnv	 + data2%xnv
		solid_add%ch2o     = data1%ch2o  + data2%ch2o
		solid_add%ch3o2h   = data1%ch3o2h+ data2%ch3o2h
		solid_add%siv      = data1%siv	 + data2%siv
		solid_add%svi      = data1%svi	 + data2%svi
		solid_add%nh4      = data1%nh4    + data2%nh4
	  end function solid_add
	  
	  function solid_sub (data1, data2)
	  	type (solid_chemical) :: solid_sub
	  	type (solid_chemical), intent(in) :: data1, data2
		
		solid_sub%o3	= data1%o3	- data2%o3
		solid_sub%h2o2	= data1%h2o2	- data2%h2o2
		solid_sub%xnv	= data1%xnv	- data2%xnv
		solid_sub%ch2o = data1%ch2o	- data2%ch2o
		solid_sub%ch3o2h = data1%ch3o2h	- data2%ch3o2h
		solid_sub%siv	= data1%siv	- data2%siv
		solid_sub%svi	= data1%svi	- data2%svi
		solid_sub%nh4	= data1%nh4	- data2%nh4
	  end function solid_sub
	  
	  function solid_times_real (data1, const)
	  	type (solid_chemical) :: solid_times_real
	  	type (solid_chemical), intent(in) :: data1
		real, intent(in)                  :: const
		
		solid_times_real%o3      = data1%o3       *const
		solid_times_real%h2o2    = data1%h2o2     *const
		solid_times_real%xnv     = data1%xnv      *const
		solid_times_real%ch2o	 = data1%ch2o	  *const
		solid_times_real%ch3o2h  = data1%ch3o2h   *const
		solid_times_real%siv     = data1%siv      *const
		solid_times_real%svi     = data1%svi      *const
		solid_times_real%nh4     = data1%nh4     *const
	  end function solid_times_real

	  function real_times_solid (const, data1)
	  	type (solid_chemical) :: real_times_solid
	  	type (solid_chemical), intent(in) :: data1
		real, intent(in)                  :: const
		
		real_times_solid%o3       = const *data1%o3
		real_times_solid%h2o2     = const *data1%h2o2
		real_times_solid%xnv      = const *data1%xnv
		real_times_solid%ch2o	  = const *data1%ch2o
		real_times_solid%ch3o2h   = const *data1%ch3o2h
		real_times_solid%siv      = const *data1%siv
		real_times_solid%svi      = const *data1%svi
		real_times_solid%nh4      = const *data1%nh4
	  end function real_times_solid
	  
	  function solid_times_int (data1, const)
	  	type (solid_chemical) :: solid_times_int
	  	type (solid_chemical), intent(in) :: data1
		integer, intent(in)               :: const
		
		solid_times_int%o3       = data1%o3       *real(const)
		solid_times_int%h2o2     = data1%h2o2     *real(const)
		solid_times_int%xnv      = data1%xnv      *real(const)
		solid_times_int%ch2o	 = data1%ch2o	  *real(const)
		solid_times_int%ch3o2h   = data1%ch3o2h   *real(const)
		solid_times_int%siv      = data1%siv      *real(const)
		solid_times_int%svi      = data1%svi      *real(const)
		solid_times_int%nh4      = data1%nh4      *real(const)
	  end function solid_times_int

	  function int_times_solid (const, data1)
	  	type (solid_chemical) :: int_times_solid
	  	type (solid_chemical), intent(in) :: data1
		integer, intent(in)               :: const
		
		int_times_solid%o3       = real(const) *data1%o3
		int_times_solid%h2o2     = real(const) *data1%h2o2
		int_times_solid%xnv      = real(const) *data1%xnv
		int_times_solid%ch2o	 = real(const) *data1%ch2o
		int_times_solid%ch3o2h   = real(const) *data1%ch3o2h
		int_times_solid%siv      = real(const) *data1%siv
		int_times_solid%svi      = real(const) *data1%svi
		int_times_solid%nh4      = real(const) *data1%nh4
	  end function int_times_solid	  		

	  function solid_div_real (data1, const)
	  	type (solid_chemical) :: solid_div_real
	  	type (solid_chemical), intent(in) :: data1
		real, intent(in)                  :: const
		
		solid_div_real%o3       = data1%o3       /const
		solid_div_real%h2o2     = data1%h2o2     /const
		solid_div_real%xnv      = data1%xnv      /const
		solid_div_real%ch2o	= data1%ch2o	 /const
		solid_div_real%ch3o2h	= data1%ch3o2h   /const
		solid_div_real%siv      = data1%siv      /const
		solid_div_real%svi      = data1%svi      /const
		solid_div_real%nh4      = data1%nh4      /const
	  end function solid_div_real

	  function solid_div_int (data1, const)
	  	type (solid_chemical) :: solid_div_int
	  	type (solid_chemical), intent(in) :: data1
		integer, intent(in)               :: const
		
		solid_div_int%o3       = data1%o3       /real(const)
		solid_div_int%h2o2     = data1%h2o2     /real(const)
		solid_div_int%xnv      = data1%xnv      /real(const)
		solid_div_int%ch2o     = data1%ch2o	/real(const)
		solid_div_int%ch3o2h   = data1%ch3o2h	/real(const)
		solid_div_int%siv      = data1%siv      /real(const)
		solid_div_int%svi      = data1%svi      /real(const)
		solid_div_int%nh4      = data1%nh4      /real(const)
	  end function solid_div_int
	  	
	end module typedef_solid
