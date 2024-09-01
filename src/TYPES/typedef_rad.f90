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
!  TYPEDEF_RAD.F:                   
!
!  Purpose:
!	Typedef for radiation data			  
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================
	
	module typedef_rad
	
	USE gridno
	USE shared_data
	
	type :: radiation
		! ---
		real, dimension(:,:,:), allocatable :: dtnet   
		real, dimension(:,:,:), allocatable :: fluxs 
		real, dimension(:,:,:), allocatable :: fluxir  
		real, dimension(:,:,:), allocatable :: frad

		real, dimension(:,:,:), allocatable :: dtsw
		real, dimension(:,:,:), allocatable :: dtlw

		real, dimension(:,:), allocatable :: tau_aer_surf
		real, dimension(:,:), allocatable :: tau_cloud_surf
		real, dimension(:,:), allocatable :: tau_rain_surf
	end type radiation

	interface operator (+)
		module procedure rad_add
	end interface
	
	interface alloc
		module procedure alloc_radia
	end interface alloc
	
	interface dealloc
		module procedure dealloc_radia
	end interface dealloc

	CONTAINS
	
	  subroutine alloc_radia (rad)
	  	type (radiation), intent(inout) :: rad
	  
	        allocate ( rad%dtnet(ip_start:ip_end,jp_start:jp_end,1:nz) )
		
#ifdef RAD_ENABLE
	        allocate ( rad%fluxs(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        allocate ( rad%fluxir(ip_start:ip_end,jp_start:jp_end,1:nz) )
	        allocate ( rad%frad(ip_start:ip_end,jp_start:jp_end,1:nz) )		

		allocate ( rad%dtsw(ip_start:ip_end,jp_start:jp_end,1:nz))
		allocate ( rad%dtlw(ip_start:ip_end,jp_start:jp_end,1:nz))

		allocate ( rad%tau_aer_surf(ip_start:ip_end,jp_start:jp_end))
		allocate ( rad%tau_cloud_surf(ip_start:ip_end,jp_start:jp_end))
		allocate ( rad%tau_rain_surf(ip_start:ip_end,jp_start:jp_end))			
#endif

		rad%dtnet = 0.0
		
#ifdef RAD_ENABLE
		rad%fluxs = 0.0
		rad%fluxir = 0.0
		rad%frad = 0.0

		rad%dtsw = 0.0
		rad%dtlw = 0.0

		rad%tau_aer_surf = 0.0
		rad%tau_cloud_surf = 0.0
		rad%tau_rain_surf = 0.0
#endif
	  return
	  end subroutine alloc_radia
	  
	  subroutine dealloc_radia (rad)
	  	type (radiation), intent(inout) :: rad
	  	
		deallocate (rad%dtnet)
		
#ifdef RAD_ENABLE
		deallocate (rad%fluxs)
		deallocate (rad%fluxir)
		deallocate (rad%frad)

		deallocate (rad%dtsw)
		deallocate (rad%dtlw)

		deallocate (rad%tau_aer_surf)
		deallocate (rad%tau_cloud_surf)
		deallocate (rad%tau_rain_surf)
#endif
	  return
	  end subroutine dealloc_radia

	  function rad_add (data1, data2)
	  	type (radiation) :: rad_add
	  	type (radiation), intent(in) :: data1, data2
		
		rad_add%dtnet    = data1%dtnet + data2%dtnet
#ifdef RAD_ENABLE
		rad_add%fluxs    = data1%fluxs + data2%fluxs
		rad_add%fluxir   = data1%fluxir + data2%fluxir
		rad_add%frad     = data1%frad + data2%frad

		rad_add%dtsw     = data1%dtsw + data2%dtsw
		rad_add%dtlw     = data1%dtlw + data2%dtlw
#endif
	  end function rad_add
		
	end module typedef_rad

