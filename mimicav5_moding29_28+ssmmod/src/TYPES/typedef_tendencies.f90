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
!	Julien Savre, Cambridge
!		
! ================================================================
	
	module typedef_tend

	use gridno
	use shared_data
	use typedef_nuclei
	
	type :: htendencies	! Tendencies for hydrometeors
		! ---
		real, dimension(:,:,:), allocatable :: q	! Ice-liquid potential temperature
		real, dimension(:,:,:), allocatable :: n  
		real, dimension(:,:,:), allocatable :: w 
	end type htendencies
	
	type :: atendencies	! Tendencies for aerosols
		! ---
		real, dimension(:,:,:), allocatable :: n
		real, dimension(:,:,:), allocatable :: ma  
		real, dimension(:,:,:), allocatable :: m 
	end type atendencies
	
	type :: itendencies	! Tendencies for INs
		! ---
		type(nuclei_in), dimension(3) :: mode 
	end type itendencies
	
	type :: ntendencies	! Tendencies for CCNs
		! ---
		real, dimension(:,:,:), allocatable :: ccn
		real, dimension(:,:,:), allocatable :: in  
	end type ntendencies
	
	type :: tendencies
		! ---
		real, dimension(:,:,:), allocatable :: es  	! potential temperature
		real, dimension(:,:,:), allocatable :: qt      ! Total water mixing ratio (kg/kg)
		real, dimension(:,:,:), allocatable :: dens
		real, dimension(:,:,:), allocatable :: p
		real, dimension(:,:,:), allocatable :: u
		real, dimension(:,:,:), allocatable :: v
		real, dimension(:,:,:), allocatable :: w
		real, dimension(:,:,:,:), allocatable :: scal
	end type tendencies
	
	interface assignment (=)
		module procedure const_to_tend
		module procedure const_to_htend
		module procedure const_to_ntend
		module procedure const_to_atend
		module procedure const_to_itend
	end interface
	
	interface alloc
		module procedure alloc_ntend, alloc_htend, alloc_atend, alloc_itend, alloc_tend
	end interface alloc
	
	interface dealloc
		module procedure dealloc_ntend, dealloc_htend, dealloc_atend, dealloc_itend, dealloc_tend
	end interface dealloc

	CONTAINS
	
	  subroutine alloc_tend (tend)
	  	type (tendencies), intent(inout) :: tend
		
		allocate (tend%u(ip_start:ip_end,jp_start:jp_end,1:nz))
		allocate (tend%v(ip_start:ip_end,jp_start:jp_end,1:nz))
		allocate (tend%w(ip_start:ip_end,jp_start:jp_end,1:nz))
		allocate (tend%es(ip_start:ip_end,jp_start:jp_end,1:nz))
		allocate (tend%qt(ip_start:ip_end,jp_start:jp_end,1:nz))
		allocate (tend%dens(ip_start:ip_end,jp_start:jp_end,1:nz))
		allocate (tend%scal(ip_start:ip_end,jp_start:jp_end,1:nz,1:nscal))
	  return
	  end subroutine alloc_tend
	
	  subroutine dealloc_tend (tend)
	  	type (tendencies), intent(inout) :: tend
		
		deallocate (tend%u)
		deallocate (tend%v)
		deallocate (tend%w)
		deallocate (tend%es)
		deallocate (tend%qt)
		deallocate (tend%dens)
		deallocate (tend%scal)
	  return
	  end subroutine dealloc_tend
	
	  subroutine alloc_ntend (tend)
	  	type (ntendencies), intent(inout) :: tend
		
		allocate (tend%ccn(ip_start:ip_end,jp_start:jp_end,1:nz))
		allocate (tend%in(ip_start:ip_end,jp_start:jp_end,1:nz))
	  return
	  end subroutine alloc_ntend
	
	  subroutine dealloc_ntend (tend)
	  	type (ntendencies), intent(inout) :: tend
		
	 	deallocate (tend%ccn)
		deallocate (tend%in)
	  return
	  end subroutine dealloc_ntend
	
	  subroutine alloc_htend (tend)
	  	type (htendencies), dimension(nhydro), intent(inout) :: tend
		integer :: h
		
		do h = 1, nhydro
		  allocate (tend(h)%q(ip_start:ip_end,jp_start:jp_end,1:nz))
		  allocate (tend(h)%n(ip_start:ip_end,jp_start:jp_end,1:nz))
		  if ( lmicro > 3 ) allocate (tend(h)%w(ip_start:ip_end,jp_start:jp_end,1:nz))
		enddo
	  return
	  end subroutine alloc_htend
	
	  subroutine dealloc_htend (tend)
	  	type (htendencies), dimension(nhydro), intent(inout) :: tend
		integer :: h
		
		do h = 1, nhydro
		  deallocate (tend(h)%q)
		  deallocate (tend(h)%n)
		  if ( lmicro > 3 ) deallocate (tend(h)%w)
		enddo
	  return
	  end subroutine dealloc_htend
	
	  subroutine alloc_atend (tend)
	  	type (atendencies), dimension(nmode), intent(inout) :: tend
		integer :: h

                if (sorder > 1) then
		  do h = 1, nmode
		    allocate (tend(h)%n(ip_start:ip_end,jp_start:jp_end,1:nz))
		    allocate (tend(h)%ma(ip_start:ip_end,jp_start:jp_end,1:nz))
		    allocate (tend(h)%m(ip_start:ip_end,jp_start:jp_end,1:nz))
		  enddo
                else
                  do h = 1, nmode
                    allocate(tend(h)%n(1,1,1))
                    allocate(tend(h)%ma(1,1,1))
                    allocate(tend(h)%m(1,1,1))
                  enddo
                endif

          return
	  end subroutine alloc_atend
	
	  subroutine dealloc_atend (tend)
	  	type (atendencies), dimension(nmode), intent(inout) :: tend
		integer :: h
		
		do h = 1, nmode
		  deallocate (tend(h)%n)
		  deallocate (tend(h)%ma)
		  deallocate (tend(h)%m)
		enddo
	  return
	  end subroutine dealloc_atend
	
	  subroutine alloc_itend (tend)
	  	type (itendencies), dimension(3), intent(inout) :: tend
		integer :: h, l
		
		do h = 1, 3
		  do l = 1, 3
		    allocate (tend(h)%mode(l)%n(ip_start:ip_end,jp_start:jp_end,1:nz))
		    allocate (tend(h)%mode(l)%nc(ip_start:ip_end,jp_start:jp_end,1:nz))
		    allocate (tend(h)%mode(l)%m(ip_start:ip_end,jp_start:jp_end,1:nz))
		    allocate (tend(h)%mode(l)%mc(ip_start:ip_end,jp_start:jp_end,1:nz))
		  enddo
		enddo
	  return
	  end subroutine alloc_itend
	
	  subroutine dealloc_itend (tend)
	  	type (itendencies), dimension(3), intent(inout) :: tend
		integer :: h, l
		
		do h = 1, 3
		  do l = 1, 3
  		    deallocate (tend(h)%mode(l)%n)
		    deallocate (tend(h)%mode(l)%nc)
		    deallocate (tend(h)%mode(l)%m)
		    deallocate (tend(h)%mode(l)%mc)
		  enddo
		enddo
	  return
	  end subroutine dealloc_itend
	
	  subroutine const_to_tend (tend, const)
	  	type (tendencies), intent(inout) :: tend
	  	real, intent(in)	 	 :: const
		
		tend%es   = const
		tend%qt   = const
		tend%dens = const
		tend%u    = const
		tend%v    = const
		tend%w    = const
		tend%scal = const
	  return
	  end subroutine const_to_tend
	
	  subroutine const_to_ntend (tend, const)
	  	type (ntendencies), intent(inout) :: tend
	  	real, intent(in)		  :: const
		
 	  	tend%ccn = const
		tend%in = const
	  return
	  end subroutine const_to_ntend
	
	  subroutine const_to_htend (tend, const)
	  	type (htendencies), intent(inout) :: tend(:)
	  	real, intent(in)		  :: const
	  	integer				  :: h
		
		do h = 1, nhydro
 	  	  tend(h)%q = const
		  tend(h)%n = const
		  if ( lmicro > 3 ) tend(h)%w = const
	   	enddo
	  return
	  end subroutine const_to_htend
	
	  subroutine const_to_atend (tend, const)
	  	type (atendencies), intent(inout) :: tend(:)
	  	real, intent(in)		  :: const
	  	integer				  :: h
		
		do h = 1, nmode
 	  	  tend(h)%n = const
		  tend(h)%ma = const
		  tend(h)%m = const
	   	enddo
	  return
	  end subroutine const_to_atend
	
	  subroutine const_to_itend (tend, const)
	  	type (itendencies), intent(inout) :: tend(:)
	  	real, intent(in)		  :: const
	  	integer				  :: h, l
		
		do h = 1, 3
		  do l = 1, 3
 	  	    tend(h)%mode(l)%n = const
		    tend(h)%mode(l)%nc = const
		    tend(h)%mode(l)%m = const
		    tend(h)%mode(l)%mc = const
	   	  enddo
		enddo
	  return
	  end subroutine const_to_itend
	  	  	    			
	end module typedef_tend

