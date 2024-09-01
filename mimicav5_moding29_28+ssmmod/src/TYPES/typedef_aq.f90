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
	
	module typedef_aq
	
	type :: aq_chemical
		! ---
		! --- Mixing ratios of gaseous chemical species
		! ---	Units are in ppb(m)
		! ---
		real :: o3	! ozone
		real :: civ	! C(IV)
		real :: h2o2	! H2O2
		real :: xnv	! N(V) = HNO3 + N2O5
		real :: ch2o	! CH2O
		real :: ch3o2h	! CH3O2h
		real :: siv	! S(IV)
		real :: svi	! S(VI)
                real :: nh4     ! NH4
		real :: hplus	! [H+]
	end type aq_chemical

	interface assignment (=)
		module procedure const2aq
	end interface

	interface operator (+)
		module procedure aq_add
	end interface

	interface operator (-)
		module procedure aq_sub
	end interface

	interface operator (*)
		module procedure aq_times_real
		module procedure real_times_aq
		module procedure aq_times_int
		module procedure int_times_aq
	end interface

	interface operator (/)
		module procedure aq_div_real
		module procedure aq_div_int
	end interface

	CONTAINS
	
	  subroutine const2aq (aq, const)
	       type (aq_chemical), intent(inout) :: aq
	       real, intent(in) 	         :: const
	       
	       aq%o3     = const
	       aq%civ    = const
	       aq%h2o2   = const
	       aq%xnv    = const
	       aq%ch2o   = const
	       aq%ch3o2h = const
	       aq%siv    = const
	       aq%svi    = const
	       aq%nh4    = const
	       aq%hplus  = const
	  return
	  end subroutine const2aq

	  function aq_add (data1, data2)
	  	type (aq_chemical) :: aq_add
	  	type (aq_chemical), intent(in) :: data1, data2
		
		aq_add%o3       = data1%o3	 + data2%o3
		aq_add%civ      = data1%civ	 + data2%civ
		aq_add%h2o2     = data1%h2o2     + data2%h2o2
		aq_add%xnv      = data1%xnv	 + data2%xnv
		aq_add%ch2o     = data1%ch2o     + data2%ch2o
		aq_add%ch3o2h   = data1%ch3o2h   + data2%ch3o2h
		aq_add%siv      = data1%siv	 + data2%siv
		aq_add%svi      = data1%svi	 + data2%svi
		aq_add%nh4      = data1%nh4	 + data2%nh4
		aq_add%hplus    = data1%hplus    + data2%hplus
	  end function aq_add
	  
	  function aq_sub (data1, data2)
	  	type (aq_chemical) :: aq_sub
	  	type (aq_chemical), intent(in) :: data1, data2
		
		aq_sub%o3	= data1%o3	 - data2%o3
		aq_sub%civ	= data1%civ	 - data2%civ
		aq_sub%h2o2	= data1%h2o2	 - data2%h2o2
		aq_sub%xnv	= data1%xnv	 - data2%xnv
		aq_sub%ch2o	= data1%ch2o	 - data2%ch2o
		aq_sub%ch3o2h	= data1%ch3o2h   - data2%ch3o2h
		aq_sub%siv	= data1%siv	 - data2%siv
		aq_sub%svi	= data1%svi	 - data2%svi
		aq_sub%nh4	= data1%nh4	 - data2%nh4
		aq_sub%hplus	= data1%hplus	 - data2%hplus
	  end function aq_sub
	  
	  function aq_times_real (data1, const)
	  	type (aq_chemical) :: aq_times_real
	  	type (aq_chemical), intent(in) :: data1
		real, intent(in)               :: const
		
		aq_times_real%o3      = data1%o3       *const
		aq_times_real%civ     = data1%civ      *const
		aq_times_real%h2o2    = data1%h2o2     *const
		aq_times_real%xnv     = data1%xnv      *const
		aq_times_real%ch2o    = data1%ch2o     *const
		aq_times_real%ch3o2h  = data1%ch3o2h   *const
		aq_times_real%siv     = data1%siv      *const
		aq_times_real%svi     = data1%svi      *const
		aq_times_real%nh4     = data1%nh4      *const
		aq_times_real%hplus   = data1%hplus    *const
	  end function aq_times_real

	  function real_times_aq (const, data1)
	  	type (aq_chemical) :: real_times_aq
	  	type (aq_chemical), intent(in) :: data1
		real, intent(in)               :: const
		
		real_times_aq%o3       = const *data1%o3
		real_times_aq%civ      = const *data1%civ
		real_times_aq%h2o2     = const *data1%h2o2
		real_times_aq%xnv      = const *data1%xnv
		real_times_aq%ch2o     = const *data1%ch2o
		real_times_aq%ch3o2h   = const *data1%ch3o2h
		real_times_aq%siv      = const *data1%siv
		real_times_aq%svi      = const *data1%svi
		real_times_aq%nh4      = const *data1%nh4
		real_times_aq%hplus    = const *data1%hplus
	  end function real_times_aq
	  
	  function aq_times_int (data1, const)
	  	type (aq_chemical) :: aq_times_int
	  	type (aq_chemical), intent(in) :: data1
		integer, intent(in)            :: const
		
		aq_times_int%o3       = data1%o3       *real(const)
		aq_times_int%civ      = data1%civ      *real(const)
		aq_times_int%h2o2     = data1%h2o2     *real(const)
		aq_times_int%xnv      = data1%xnv      *real(const)
		aq_times_int%ch2o     = data1%ch2o     *real(const)
		aq_times_int%ch3o2h   = data1%ch3o2h   *real(const)
		aq_times_int%siv      = data1%siv      *real(const)
		aq_times_int%svi      = data1%svi      *real(const)
		aq_times_int%nh4      = data1%nh4      *real(const)
		aq_times_int%hplus    = data1%hplus    *real(const)
	  end function aq_times_int

	  function int_times_aq (const, data1)
	  	type (aq_chemical) :: int_times_aq
	  	type (aq_chemical), intent(in) :: data1
		integer, intent(in)            :: const
		
		int_times_aq%o3       = real(const) *data1%o3
		int_times_aq%civ      = real(const) *data1%civ
		int_times_aq%h2o2     = real(const) *data1%h2o2
		int_times_aq%xnv      = real(const) *data1%xnv
		int_times_aq%ch2o     = real(const) *data1%ch2o
		int_times_aq%ch3o2h   = real(const) *data1%ch3o2h
		int_times_aq%siv      = real(const) *data1%siv
		int_times_aq%svi      = real(const) *data1%svi
		int_times_aq%nh4      = real(const) *data1%nh4
		int_times_aq%hplus    = real(const) *data1%hplus
	  end function int_times_aq	  		

	  function aq_div_real (data1, const)
	  	type (aq_chemical) :: aq_div_real
	  	type (aq_chemical), intent(in) :: data1
		real, intent(in)               :: const
		
		aq_div_real%o3       = data1%o3       /const
		aq_div_real%civ      = data1%civ      /const
		aq_div_real%h2o2     = data1%h2o2     /const
		aq_div_real%xnv      = data1%xnv      /const
		aq_div_real%ch2o     = data1%ch2o     /const
		aq_div_real%ch3o2h   = data1%ch3o2h   /const
		aq_div_real%siv      = data1%siv      /const
		aq_div_real%svi      = data1%svi      /const
		aq_div_real%nh4      = data1%nh4      /const
		aq_div_real%hplus    = data1%hplus    /const
	  end function aq_div_real

	  function aq_div_int (data1, const)
	  	type (aq_chemical) :: aq_div_int
	  	type (aq_chemical), intent(in) :: data1
		integer, intent(in)            :: const
		
		aq_div_int%o3	    = data1%o3       /real(const)
		aq_div_int%civ      = data1%civ      /real(const)
		aq_div_int%h2o2     = data1%h2o2     /real(const)
		aq_div_int%xnv      = data1%xnv      /real(const)
		aq_div_int%ch2o     = data1%ch2o     /real(const)
		aq_div_int%ch3o2h   = data1%ch3o2h   /real(const)
		aq_div_int%siv      = data1%siv      /real(const)
		aq_div_int%svi      = data1%svi      /real(const)
		aq_div_int%nh4      = data1%nh4      /real(const)
		aq_div_int%hplus    = data1%hplus    /real(const)
	  end function aq_div_int
	  	
	end module typedef_aq
