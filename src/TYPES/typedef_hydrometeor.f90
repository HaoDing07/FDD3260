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
!  TYPEDEF_HYDROMETEOR.F:                   
!
!  Purpose:
!	Typedef for derived-type data: hydrometeor			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
! ================================================================
	
	module typedef_hydrometeor
	
	USE gridno
	USE shared_data
	
	type :: hydrometeor
		! ---		
		! --- Mixing ratios (kg/kg) and number concentrations (1/kg)
		! ---
		real, dimension(:,:,:), allocatable :: q	! mixing ratio of hydrometeors
		real, dimension(:,:,:), allocatable :: n	! number concentration of hydrometeors
		real, dimension(:,:,:), allocatable :: w	! liquid water mass (if applicable)
		real, dimension(:,:,:), allocatable :: d	! Mean size of hydrometeors
		real, dimension(:,:,:), allocatable :: vp	! Mass weighted precip velocity
	end type hydrometeor
	
	type :: hydro_param
		! ---		
		! --- Parameters for mass and terminal fall speed parameterizations
		! ---
		real :: cm	! Constant for mass power law
		real :: ctv	! Constant for terminal fall speed power law
		real :: am
		real :: btv
		real :: nu	! Parameter for size distributions
		real :: mu	! Parameter for size distributions
		real :: av, bv, cv  ! Parameters for ventilation
		real :: cap	! Capacitance
		real :: a0	! One-moment preexp factor
		real :: xmin, xmax ! Particle mass limits
	end type hydro_param	

	type :: hydro_const
		! ---		
		! --- Constants for Seifert-Beheng's scheme
		! ---
		real :: clam, clam1
		real :: cmm, cdm, cvm, czm, cem, cr6
		real :: cvpq, cvpn, cimp
		real :: aventn, bventn, cventn
		real :: aventq, bventq, cventq
		real :: deltaq, deltan
		real :: thetaq, thetan
		real :: deltaijq(6), deltaijn(6)
		real :: thetaijq(6), thetaijn(6)
	end type hydro_const	

	interface assignment (=)
		module procedure const_to_hydro
		module procedure hydro_to_hydro
	end interface

	interface operator (+)
		module procedure hydro_add
	end interface

	interface operator (-)
		module procedure hydro_sub
	end interface

	interface operator (*)
		module procedure hydro_times_real
		module procedure real_times_hydro
	end interface

	interface operator (/)
		module procedure hydro_div_real
	end interface
	
	interface alloc
		module procedure allocate_hydro
	end interface
	
	interface dealloc
		module procedure deallocate_hydro
	end interface

	CONTAINS
	
	  subroutine allocate_hydro (hydro)
	  	type(hydrometeor), dimension(1:nhydro), intent(inout) :: hydro
		integer :: h
		
		do h = 1, nhydro
		  allocate ( hydro(h)%q(ip_start:ip_end,jp_start:jp_end,1:nz) )
		  allocate ( hydro(h)%n(ip_start:ip_end,jp_start:jp_end,1:nz) )
		  allocate ( hydro(h)%vp(ip_start:ip_end,jp_start:jp_end,1:nz) )
		  if (lmicro > 3) allocate ( hydro(h)%w(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
		  hydro(h)%q = 0.
		  hydro(h)%n = 0.
		  hydro(h)%vp = 0.
		  if (lmicro > 3) hydro(h)%w = 0.
	        enddo
		
	  return
	  end subroutine allocate_hydro
	
	  subroutine deallocate_hydro (hydro)
	  	type(hydrometeor), dimension(1:nhydro), intent(inout) :: hydro
		integer :: h
		
		do h = 1, nhydro
		  deallocate ( hydro(h)%q )
		  deallocate ( hydro(h)%n )
		  deallocate ( hydro(h)%vp )
		  if (lmicro > 3) deallocate ( hydro(h)%w )
	        enddo
		
	  return
	  end subroutine deallocate_hydro
	
	  subroutine const_to_hydro (hydro, const)
	  	type (hydrometeor), dimension(nhydro), intent(inout) :: hydro
	  	real, intent(in) :: const
		integer :: h
		
		do h = 1, nhydro
  		  hydro(h)%q = const
		  hydro(h)%n = const
		  hydro(h)%vp = const
		  if (lmicro > 3) hydro(h)%w = const
	  	enddo
	  return
          end subroutine const_to_hydro
	
	  subroutine hydro_to_hydro (hydro1, hydro2)
	  	type (hydrometeor), dimension(nhydro), intent(inout) :: hydro1
	  	type (hydrometeor), dimension(nhydro), intent(in) :: hydro2
		integer :: h
		
		do h = 1, nhydro
  		  hydro1(h)%q = hydro2(h)%q
		  hydro1(h)%n = hydro2(h)%n
		  hydro1(h)%vp = hydro2(h)%vp
		  if (lmicro > 3) hydro1(h)%w = hydro2(h)%w
	  	enddo
	  return
          end subroutine hydro_to_hydro
	
	  function hydro_add (data1, data2)
	  	type (hydrometeor), dimension(nhydro) :: hydro_add
	  	type (hydrometeor), dimension(nhydro), intent(in) :: data1, data2
		integer :: h
		
		do h = 1, nhydro
		  hydro_add(h)%q = data1(h)%q + data2(h)%q
	 	  hydro_add(h)%n = data1(h)%n + data2(h)%n
	  	enddo
		  if (lmicro > 3) hydro_add(h)%w = data1(h)%w + data2(h)%w
	  end function hydro_add
	  
	  function hydro_sub (data1, data2)
	  	type (hydrometeor), dimension(nhydro) :: hydro_sub
	  	type (hydrometeor), dimension(nhydro), intent(in) :: data1, data2
		integer :: h
		
		do h = 1, nhydro
		  hydro_sub(h)%q = data1(h)%q - data2(h)%q
		  hydro_sub(h)%n = data1(h)%n - data2(h)%n
		  if (lmicro > 3) hydro_sub(h)%w = data1(h)%w - data2(h)%w
	  	enddo
	  end function hydro_sub

	  function real_times_hydro (const, data1) result(data2)
	  	type (hydrometeor), dimension(nhydro), intent(in) :: data1
	  	type (hydrometeor), dimension(nhydro) :: data2
		real, intent(in) :: const
		integer :: h
		
		do h = 1, nhydro
		  data2(h)%q = const * data1(h)%q
		  data2(h)%n = const * data1(h)%n
		  if (lmicro > 3) data2(h)%w = const * data1(h)%w
	  	enddo
	  end function real_times_hydro
	  
	  function hydro_times_real (data1, const) result(data2)
	  	type (hydrometeor), dimension(nhydro), intent(in) :: data1
	  	type (hydrometeor), dimension(nhydro) :: data2
		real, intent(in) :: const
		integer :: h
		
		do h = 1, nhydro
		  data2(h)%q = data1(h)%q * const
		  data2(h)%n = data1(h)%n * const
		  if (lmicro > 3) data2(h)%w = data1(h)%w * const
	  	enddo
	  end function hydro_times_real
	  
	  function hydro_div_real (data1, const) result(data2)
	  	type (hydrometeor), dimension(nhydro), intent(in) :: data1
	  	type (hydrometeor), dimension(nhydro) :: data2
		real, intent(in)               :: const
		integer :: h
		
		do h = 1, nhydro
		  data2(h)%q = data1(h)%q /const
		  data2(h)%n = data1(h)%n /const
		  if (lmicro > 3) data2(h)%w = data1(h)%w / const
	  	enddo
	  end function hydro_div_real
	  	  	    			
	end module typedef_hydrometeor

