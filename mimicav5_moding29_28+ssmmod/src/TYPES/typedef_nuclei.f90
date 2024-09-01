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
!  TYPEDEF_NUCLEI.F:                   
!
!  Purpose:
!	Typedef for nuclei			  
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================
	
	module typedef_nuclei
	
	USE gridno
	USE shared_data
	
	type :: nuclei
		! ---
		real, dimension(:,:,:), allocatable :: ccn		! Cloud condensation nuclei
		real, dimension(:,:,:), allocatable :: in		! Ice nuclei
	end type nuclei
!
	type :: nuclei_in
		real, dimension(:,:,:), allocatable :: n
		real, dimension(:,:,:), allocatable :: m
		real, dimension(:,:,:), allocatable :: nc
		real, dimension(:,:,:), allocatable :: mc
	end type
!
	type :: nuclei_in0
		real, dimension(nz) :: n
		real, dimension(nz) :: m
		real, dimension(nz) :: nc
		real, dimension(nz) :: mc
	end type
!
!  Classical Nucleation Theory block
!
        type :: nuclei_3d
	        type(nuclei_in), dimension(3) :: mode
        end type

        type :: nuclei0
	        type(nuclei_in0), dimension(3) :: mode
        end type

	type :: nuc_param
		! ---
		real :: n, d, s			! fine mode
		real :: nc, dc, sc		! Coarse mode
		real :: rhop
	end type

#ifdef NUC_CNT
	type :: nuc_prop
		! ---
		real :: frac, fracc		! Initial fractions
		real, dimension(3) :: theta, sigma
	end type
!
!  Diehl and Wurzler block
!
#else
	type :: nuc_prop
		! ---
		real :: frac, fracc	! Initial fractions
		real :: Beff, a		! Freezing efficiencies (Diehl & Wurzler)
		real :: ac		! Contact freezing parameter (Diehl et a. 2006)
		real :: bc		! Contact freezing parameter (Diehl et a. 2006)
		real :: theta		! Contact angle for deposition
	end type
#endif
	
	interface assignment (=)
		module procedure const_to_nuc
		module procedure nuclei_to_nuclei
	end interface

	interface operator (+)
		module procedure nuc_add
	end interface

	interface operator (-)
		module procedure nuc_sub
	end interface

	interface operator (*)
		module procedure nuc_times_real
		module procedure real_times_nuc
		module procedure nuc_times_int
		module procedure int_times_nuc
	end interface

	interface operator (/)
		module procedure nuc_div_real
		module procedure nuc_div_int
	end interface
	
	interface alloc
		module procedure alloc_nuc
		module procedure alloc_in
	end interface alloc
	
	interface dealloc
		module procedure dealloc_nuc
		module procedure dealloc_in
	end interface dealloc

	CONTAINS
	
	  subroutine alloc_nuc (nuc)
	  	type(nuclei), intent(inout) :: nuc

		allocate ( nuc%ccn(ip_start:ip_end,jp_start:jp_end,1:nz) )
		allocate ( nuc%in(ip_start:ip_end,jp_start:jp_end,1:nz) )
	  
	  	nuc%ccn = 0.0
		nuc%in = 0.0		
	  return
	  end subroutine alloc_nuc
	  
	  subroutine dealloc_nuc (nuc)
	  	type(nuclei), intent(inout) :: nuc
		
		deallocate (nuc%ccn)
		deallocate (nuc%in)
	  return
	  end subroutine dealloc_nuc
	
	  subroutine alloc_in (nucin)
	  	type(nuclei_3d), dimension(3), intent(inout) :: nucin
		integer :: h, l

		do h = 1, 3
		  do l = 1, 3
		    allocate ( nucin(h)%mode(l)%n(ip_start:ip_end,jp_start:jp_end,1:nz) )
		    allocate ( nucin(h)%mode(l)%nc(ip_start:ip_end,jp_start:jp_end,1:nz) )
		    allocate ( nucin(h)%mode(l)%m(ip_start:ip_end,jp_start:jp_end,1:nz) )
		    allocate ( nucin(h)%mode(l)%mc(ip_start:ip_end,jp_start:jp_end,1:nz) )
	  
	  	    nucin(h)%mode(l)%n = 0.0
		    nucin(h)%mode(l)%nc = 0.0	
	  	    nucin(h)%mode(l)%m = 0.0
		    nucin(h)%mode(l)%mc = 0.0	
		  enddo
		enddo	
	  return
	  end subroutine alloc_in
	  
	  subroutine dealloc_in (nucin)
	  	type(nuclei_3d), dimension(3), intent(inout) :: nucin
		integer :: h, l
		
		do h = 1, 3
		  do l = 1, 3
  		    deallocate (nucin(h)%mode(l)%n)
		    deallocate (nucin(h)%mode(l)%nc)
  		    deallocate (nucin(h)%mode(l)%m)
		    deallocate (nucin(h)%mode(l)%mc)
		  enddo
		enddo
	  return
	  end subroutine dealloc_in
	
	  subroutine const_to_nuc (nuc, const)
	  	type (nuclei), intent(inout) :: nuc
	  	real, intent(in)	          :: const
		
		nuc%ccn = const
		nuc%in = const
	  return
	  end subroutine const_to_nuc

	  subroutine nuclei_to_nuclei (nuclei1, nuclei2)
	  	type(nuclei_3d), dimension(3), intent(inout) :: nuclei1
	  	type(nuclei_3d), dimension(3), intent(in) :: nuclei2
		integer h, l
				
		do h = 1, 3
		  do l = 1, 3
		    nuclei1(h)%mode(l)%n=nuclei2(h)%mode(l)%n
		    nuclei1(h)%mode(l)%nc=nuclei2(h)%mode(l)%nc
#ifdef NUC_CNT1
		    nuclei1(h)%mode(l)%m=nuclei2(h)%mode(l)%m
		    nuclei1(h)%mode(l)%mc=nuclei2(h)%mode(l)%mc
#endif
		  enddo
		enddo
	  return
	  end subroutine nuclei_to_nuclei
	  
	  function nuc_add (data1, data2)
	  	type (nuclei) :: nuc_add
	  	type (nuclei), intent(in) :: data1, data2
		
		nuc_add%ccn = data1%ccn + data2%ccn
		nuc_add%in = data1%in + data2%in
	  end function nuc_add
	  
	  function nuc_sub (data1, data2)
	  	type (nuclei) :: nuc_sub
	  	type (nuclei), intent(in) :: data1, data2
		
		nuc_sub%ccn = data1%ccn - data2%ccn
		nuc_sub%in = data1%in - data2%in
	  end function nuc_sub
	  
	  function nuc_times_real (data1, const)
	  	type (nuclei) :: nuc_times_real
	  	type (nuclei), intent(in) :: data1
		real, intent(in)               :: const
		
		nuc_times_real%ccn = data1%ccn *const
		nuc_times_real%in = data1%in *const
	  end function nuc_times_real

	  function real_times_nuc (const, data1)
	  	type (nuclei) :: real_times_nuc
	  	type (nuclei), intent(in) :: data1
		real, intent(in)               :: const
		
		real_times_nuc%ccn = const *data1%ccn
		real_times_nuc%in = const *data1%in
	  end function real_times_nuc
	  
	  function nuc_times_int (data1, const)
	  	type (nuclei) :: nuc_times_int
	  	type (nuclei), intent(in) :: data1
		integer, intent(in)            :: const
		
		nuc_times_int%ccn = data1%ccn *real(const)
		nuc_times_int%in = data1%in *real(const)
	  end function nuc_times_int

	  function int_times_nuc (const, data1)
	  	type (nuclei) :: int_times_nuc
	  	type (nuclei), intent(in) :: data1
		integer, intent(in)            :: const
		
		int_times_nuc%ccn = real(const) *data1%ccn
		int_times_nuc%in = real(const) *data1%in
	  end function int_times_nuc	  		

	  function nuc_div_real (data1, const)
	  	type (nuclei) :: nuc_div_real
	  	type (nuclei), intent(in) :: data1
		real, intent(in)          :: const
		
		nuc_div_real%ccn = data1%ccn /const
		nuc_div_real%in = data1%in /const
	  end function nuc_div_real

	  function nuc_div_int (data1, const)
	  	type (nuclei) :: nuc_div_int
	  	type (nuclei), intent(in) :: data1
		integer, intent(in)            :: const
		
		nuc_div_int%ccn = data1%ccn /real(const)
		nuc_div_int%in = data1%in /real(const)
	  end function nuc_div_int
	  	  	    			
	end module typedef_nuclei

