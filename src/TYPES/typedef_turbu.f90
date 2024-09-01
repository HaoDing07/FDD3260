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
!  TYPEDEF_turbu.F:                   
!
!  Purpose:
!	Typedef for turbulence related variables		  
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================
	
	module typedef_turbu
	
	USE gridno
	USE shared_data
	
	type :: turbulence
		! ---
		real, dimension(:,:,:), allocatable :: ksgs	! SGS turbulent kinetic energy
		real, dimension(:,:,:), allocatable :: fkv  	! Vertical eddy diffusivity
	end type turbulence
	
	type :: turbudiag
		! ---
		real, dimension(:,:,:), allocatable :: kres	! Resolved turbulent kinetic energy
		real, dimension(:,:,:), allocatable :: bfsgs   ! SGS buoyancy fluxes
		real, dimension(:,:,:), allocatable :: bfres   ! resolved buoyancy fluxes
		real, dimension(:,:,:), allocatable :: ufsgs   ! U SGS turbulent flux
		real, dimension(:,:,:), allocatable :: vfsgs   ! V SGS turbulent flux
		real, dimension(:,:,:), allocatable :: ptfsgs   ! PT SGS turbulent flux
		real, dimension(:,:,:), allocatable :: qtfsgs   ! QT SGS turbulent flux
		real, dimension(:,:,:), allocatable :: wvar    ! W variance
		real, dimension(:,:,:), allocatable :: wske    ! W skewness
		real, dimension(:,:,:), allocatable :: N2      ! Brunt-Vaisala frequency
		real, dimension(:,:,:), allocatable :: S2      ! Deformation squared
	end type turbudiag

	interface assignment (=)
		module procedure const_to_turbu
	end interface

	interface operator (+)
		module procedure turbu_add
	end interface

	interface operator (-)
		module procedure turbu_sub
	end interface

	interface operator (*)
		module procedure turbu_times_real
		module procedure real_times_turbu
		module procedure turbu_times_int
		module procedure int_times_turbu
	end interface

	interface operator (/)
		module procedure turbu_div_real
		module procedure turbu_div_int
	end interface
	
	interface alloc
		module procedure alloc_turbu, alloc_turbudiag
	end interface
	
	interface dealloc
		module procedure dealloc_turbu, dealloc_turbudiag
	end interface

	CONTAINS
	
	  subroutine alloc_turbu (turbul)
	  	type (turbulence), intent(inout) :: turbul
		
		if ( with_dif ) then
  	  	  allocate( turbul%ksgs(ip_start:ip_end,jp_start:jp_end,1:nz) )
		  allocate( turbul%fkv(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		  turbul%ksgs = 0.
		  turbul%fkv = 0.
	  	endif
	  return
	  end subroutine alloc_turbu
	  
	  subroutine dealloc_turbu (turbul)
	  	type (turbulence), intent(inout) :: turbul
	  
		if ( with_dif ) then
	  	  deallocate (turbul%ksgs)
		  deallocate (turbul%fkv)
	  	endif
	  return
	  end subroutine dealloc_turbu
	
	  subroutine alloc_turbudiag (turbul)
	  	type (turbudiag), intent(inout) :: turbul
		
		if (out_tke) allocate( turbul%wvar(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_tke) allocate( turbul%wske(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_tke) allocate( turbul%kres(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_tke) allocate( turbul%S2(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_tke) allocate( turbul%N2(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_flut.and.out_buoy) allocate( turbul%bfres(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_fsgs.and.out_buoy) allocate( turbul%bfsgs(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_fsgs.and.out_u) allocate( turbul%ufsgs(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_fsgs.and.out_v) allocate( turbul%vfsgs(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_fsgs.and.out_pt) allocate( turbul%ptfsgs(ip_start:ip_end,jp_start:jp_end,1:nz) )
		if (out_fsgs.and.out_qt) allocate( turbul%qtfsgs(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		if (out_tke) turbul%wvar = 0.
		if (out_tke) turbul%wske = 0.
		if (out_tke) turbul%kres = 0.
		if (out_tke) turbul%S2 = 0.
		if (out_tke) turbul%N2 = 0.
		if (out_flut.and.out_buoy) turbul%bfres = 0.
		if (out_fsgs.and.out_buoy) turbul%bfsgs = 0.
		if (out_fsgs.and.out_u) turbul%ufsgs = 0.
		if (out_fsgs.and.out_v) turbul%vfsgs = 0.
		if (out_fsgs.and.out_pt) turbul%ptfsgs = 0.
		if (out_fsgs.and.out_qt) turbul%qtfsgs = 0.
	  return
	  end subroutine alloc_turbudiag
	  
	  subroutine dealloc_turbudiag (turbul)
	  	type (turbudiag), intent(inout) :: turbul
	  
		if (out_tke) deallocate (turbul%wvar)
		if (out_tke) deallocate (turbul%wske)
	  	if (out_tke) deallocate (turbul%kres)
		if (out_tke) deallocate (turbul%S2)
		if (out_tke) deallocate (turbul%N2)
		if (out_flut.and.out_buoy) deallocate (turbul%bfres)
		if (out_fsgs.and.out_buoy) deallocate (turbul%bfsgs)
		if (out_fsgs.and.out_u) deallocate (turbul%ufsgs)
		if (out_fsgs.and.out_v) deallocate (turbul%vfsgs)
		if (out_fsgs.and.out_pt) deallocate (turbul%ptfsgs)
		if (out_fsgs.and.out_qt) deallocate (turbul%qtfsgs)
	  return
	  end subroutine dealloc_turbudiag
	
	  subroutine const_to_turbu (turbul, const)
	  	type (turbulence), intent(inout) :: turbul
	  	real, intent(in)		:: const
		
		turbul%ksgs    = const
		turbul%fkv     = const
	  return
	  end subroutine const_to_turbu

	  function turbu_add (data1, data2)
	  	type (turbulence) :: turbu_add
	  	type (turbulence), intent(in) :: data1, data2
		
		turbu_add%ksgs    = data1%ksgs    + data2%ksgs
		turbu_add%fkv     = data1%fkv    + data2%fkv
	  end function turbu_add
	  
	  function turbu_sub (data1, data2)
	  	type (turbulence) :: turbu_sub
	  	type (turbulence), intent(in) :: data1, data2
		
		turbu_sub%ksgs    = data1%ksgs    - data2%ksgs
		turbu_sub%fkv     = data1%fkv    - data2%fkv
	  end function turbu_sub
	  
	  function turbu_times_real (data1, const)
	  	type (turbulence) :: turbu_times_real
	  	type (turbulence), intent(in) :: data1
		real, intent(in)	     :: const
		
		turbu_times_real%ksgs    = data1%ksgs    *const
		turbu_times_real%fkv     = data1%fkv    *const
	  end function turbu_times_real

	  function real_times_turbu (const, data1)
	  	type (turbulence) :: real_times_turbu
	  	type (turbulence), intent(in) :: data1
		real, intent(in)	     :: const
		
		real_times_turbu%ksgs    = const *data1%ksgs
		real_times_turbu%fkv     = const *data1%fkv
	  end function real_times_turbu
	  
	  function turbu_times_int (data1, const)
	  	type (turbulence) :: turbu_times_int
	  	type (turbulence), intent(in) :: data1
		integer, intent(in)	     :: const
		
		turbu_times_int%ksgs    = data1%ksgs    *real(const)
		turbu_times_int%fkv     = data1%fkv    *real(const)
	  end function turbu_times_int

	  function int_times_turbu (const, data1)
	  	type (turbulence) :: int_times_turbu
	  	type (turbulence), intent(in) :: data1
		integer, intent(in)	     :: const
		
		int_times_turbu%ksgs    = real(const) *data1%ksgs
		int_times_turbu%fkv     = real(const) *data1%fkv
	  end function int_times_turbu	  		

	  function turbu_div_real (data1, const)
	  	type (turbulence) :: turbu_div_real
	  	type (turbulence), intent(in) :: data1
		real, intent(in)	     :: const
		
		turbu_div_real%ksgs    = data1%ksgs	 /const
		turbu_div_real%fkv     = data1%fkv	 /const
	  end function turbu_div_real

	  function turbu_div_int (data1, const)
	  	type (turbulence) :: turbu_div_int
	  	type (turbulence), intent(in) :: data1
		integer, intent(in)	     :: const
		
		turbu_div_int%ksgs    = data1%ksgs    /real(const)
		turbu_div_int%fkv     = data1%fkv     /real(const)
	  end function turbu_div_int
	  	  	    			
	end module typedef_turbu

