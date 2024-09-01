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
!  TYPEDEF_surf.F:                   
!
!  Purpose:
!	Typedef for surface quantities			  
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================
	
	module typedef_surf
	
	USE gridno
	USE shared_data
	
	type :: surface
		! ---
		real, dimension(:,:), allocatable :: cumul	! Accumulated rain
		real, dimension(:,:), allocatable :: precip	! Instantaneous rain
		real, dimension(:,:), allocatable :: cond, evap    ! Total column condensation
		real, dimension(:,:), allocatable :: ustar	! Friction velocity
		real, dimension(:,:), allocatable :: mflux      ! Momentum flux
		real, dimension(:,:), allocatable :: esflux     ! Energy/sensible heat flux
		real, dimension(:,:), allocatable :: qvflux     ! Wapor/latent heat flux
		real, dimension(:,:), allocatable :: lwf, swf   ! Long wave and short wave radiative fluxes
		real, dimension(:,:), allocatable :: olr, toa	! Outgoing longwave radiation
		real, dimension(:,:), allocatable :: osr, isr	! Outgoing /incoming shortwave radiation
		real, dimension(:,:), allocatable :: olr_cs, osr_cs  ! Outgoing longwave/shortwave radiation clear sky		
	end type surface
!
	type :: columns
		real, dimension(:,:), allocatable :: lwp, iwp ! LWP, IWP 
		real, dimension(:,:), allocatable :: cwp, rwp ! CWP, RWP (cloud and rain water path)
		real, dimension(:,:), allocatable :: wvp     ! Water vapor path
		real, dimension(:,:), allocatable :: wvpt     ! Water vapor path in troposphere
		real, dimension(:,:), allocatable :: wvpb     ! Water vapor path in BL
		real, dimension(:,:), allocatable :: cmse    ! Column MSE
		real, dimension(:,:), allocatable :: cmset    ! Column MSE in troposphere
		real, dimension(:,:), allocatable :: cmseb    ! Column MSE in BL
		real, dimension(:,:), allocatable :: cmfl    ! Column mass flux
		real, dimension(:,:), allocatable :: cape    ! CAPE
		real, dimension(:,:), allocatable :: cin     ! CIN
		real, dimension(:,:), allocatable :: lcl, lfc  ! LCL and LFC
		real, dimension(:,:), allocatable :: ctop    ! Cloud top
		real, dimension(:,:), allocatable :: bltop   ! BL top
		real, dimension(:,:), allocatable :: cbas    ! Cloud base
		real, dimension(:,:), allocatable :: cpint0, cpint   ! Cold pool intensity
		real, dimension(:,:), allocatable :: intsca  ! Vertically integrated scalar (to be selected)
		real, dimension(:,:), allocatable :: zinv    ! Inversion height (SE Atlantic case)
	end type columns

	interface assignment (=)
		module procedure const_to_surf
	end interface

	interface operator (+)
		module procedure surf_add
	end interface
	
	interface alloc
		module procedure alloc_surface, alloc_columns
	end interface alloc
	
	interface dealloc
		module procedure dealloc_surface, dealloc_columns
	end interface dealloc

	CONTAINS
	
	  subroutine alloc_surface (surf)
	  	type (surface), intent(inout) :: surf
		
		allocate ( surf%ustar(ip_start:ip_end,jp_start:jp_end) )
		allocate ( surf%mflux(ip_start:ip_end,jp_start:jp_end) )
		allocate ( surf%esflux(ip_start:ip_end,jp_start:jp_end) )
		allocate ( surf%qvflux(ip_start:ip_end,jp_start:jp_end) )
		if (out_rrate) allocate ( surf%cumul(ip_start:ip_end,jp_start:jp_end) )
		if (out_rrate) allocate ( surf%precip(ip_start:ip_end,jp_start:jp_end) )
		if (out_rrate) allocate ( surf%cond(ip_start:ip_end,jp_start:jp_end) )
		if (out_rrate) allocate ( surf%evap(ip_start:ip_end,jp_start:jp_end) )
		if (out_srad) allocate ( surf%lwf(ip_start:ip_end,jp_start:jp_end) )
		if (out_srad) allocate ( surf%swf(ip_start:ip_end,jp_start:jp_end) )
		if (out_olr) allocate ( surf%olr(ip_start:ip_end,jp_start:jp_end) )
		if (out_olr) allocate ( surf%toa(ip_start:ip_end,jp_start:jp_end) )
		if (out_osr) allocate ( surf%osr(ip_start:ip_end,jp_start:jp_end) )
		if (out_osr) allocate ( surf%isr(ip_start:ip_end,jp_start:jp_end) )	
		if (out_ocs) allocate ( surf%olr_cs(ip_start:ip_end,jp_start:jp_end) )
		if (out_ocs) allocate ( surf%osr_cs(ip_start:ip_end,jp_start:jp_end) )			
		
		surf%ustar = 0.0
		surf%mflux = 0.0
		surf%esflux = 0.0
		surf%qvflux = 0.0
		if (out_rrate) surf%cumul = 0.0
		if (out_rrate) surf%precip = 0.0
		if (out_rrate) surf%cond = 0.0
		if (out_rrate) surf%evap = 0.0
		if (out_srad) surf%lwf = 0.0
		if (out_srad) surf%swf = 0.0
		if (out_olr) surf%olr = 0.0
		if (out_olr) surf%toa = 0.0		
		if (out_osr) surf%osr = 0.0		
		if (out_osr) surf%isr = 0.0
		if (out_ocs) surf%olr_cs = 0.0	
		if (out_ocs) surf%osr_cs = 0.0			
	  end subroutine alloc_surface
	
	  subroutine dealloc_surface (surf)
	  	type (surface), intent(inout) :: surf		
		deallocate ( surf%ustar )
		deallocate ( surf%mflux )
		deallocate ( surf%esflux )
		deallocate ( surf%qvflux )
		if (out_rrate) deallocate ( surf%cumul )
		if (out_rrate) deallocate ( surf%precip )
		if (out_rrate) deallocate ( surf%cond )
		if (out_rrate) deallocate ( surf%evap )
		if (out_srad) deallocate ( surf%lwf )
		if (out_srad) deallocate ( surf%swf )
		if (out_olr) deallocate ( surf%olr )
		if (out_olr) deallocate ( surf%toa )
		if (out_osr) deallocate ( surf%osr )
		if (out_osr) deallocate ( surf%isr )
		if (out_ocs) deallocate ( surf%olr_cs )	
		if (out_ocs) deallocate ( surf%osr_cs )			
	  end subroutine dealloc_surface
	
	  subroutine alloc_columns (surf)
	  	type (columns), intent(inout) :: surf
		
		if (out_lwp) allocate ( surf%lwp(ip_start:ip_end,jp_start:jp_end) )
		if (out_lwp) allocate ( surf%iwp(ip_start:ip_end,jp_start:jp_end) )
		if (out_cwp) allocate ( surf%cwp(ip_start:ip_end,jp_start:jp_end) )
		if (out_cwp) allocate ( surf%rwp(ip_start:ip_end,jp_start:jp_end) )		
		if (out_wvp) allocate ( surf%wvp(ip_start:ip_end,jp_start:jp_end) )
		if (out_wvp) allocate ( surf%wvpt(ip_start:ip_end,jp_start:jp_end) )
		if (out_wvp) allocate ( surf%wvpb(ip_start:ip_end,jp_start:jp_end) )
		if (out_cmse) allocate ( surf%cmse(ip_start:ip_end,jp_start:jp_end) )
		if (out_cmse) allocate ( surf%cmset(ip_start:ip_end,jp_start:jp_end) )
		if (out_cmse) allocate ( surf%cmseb(ip_start:ip_end,jp_start:jp_end) )
		if (out_cmfl) allocate ( surf%cmfl(ip_start:ip_end,jp_start:jp_end) )
		if (out_cape) allocate ( surf%cape(ip_start:ip_end,jp_start:jp_end) )
		if (out_cape) allocate ( surf%cin(ip_start:ip_end,jp_start:jp_end) )
		if (out_cape) allocate ( surf%lcl(ip_start:ip_end,jp_start:jp_end) )
		if (out_cape) allocate ( surf%lfc(ip_start:ip_end,jp_start:jp_end) )
		if (out_ctop) allocate ( surf%ctop(ip_start:ip_end,jp_start:jp_end) )
		if (out_ctop) allocate ( surf%cbas(ip_start:ip_end,jp_start:jp_end) )
		if (out_zinv) allocate ( surf%bltop(ip_start:ip_end,jp_start:jp_end) )
		if (out_zinv) allocate ( surf%zinv(ip_start:ip_end,jp_start:jp_end) )
		if (out_cp)  allocate ( surf%cpint0(ip_start:ip_end,jp_start:jp_end) )
		if (out_cp)  allocate ( surf%cpint(ip_start:ip_end,jp_start:jp_end) )
		if (out_ints) allocate ( surf%intsca(ip_start:ip_end,jp_start:jp_end) )
		
		if (out_lwp) surf%lwp = 0.0
		if (out_lwp) surf%iwp = 0.0
		if (out_cwp) surf%cwp = 0.0
		if (out_cwp) surf%rwp = 0.0		
		if (out_wvp) surf%wvp = 0.0
		if (out_wvp) surf%wvpt = 0.0
		if (out_wvp) surf%wvpb = 0.0
		if (out_cmse) surf%cmse = 0.0
		if (out_cmse) surf%cmset = 0.0
		if (out_cmse) surf%cmseb = 0.0
		if (out_cmfl) surf%cmfl = 0.0
		if (out_cape) surf%cape = 0.0
		if (out_cape) surf%cin = 0.0
		if (out_cape) surf%lcl = 0.0
		if (out_cape) surf%lfc = 0.0
		if (out_ctop) surf%ctop = 0.0
		if (out_ctop) surf%cbas = 0.0
		if (out_zinv) surf%bltop = 0.0
		if (out_zinv) surf%zinv=0.0
		if (out_cp)  surf%cpint0 = 0.0
		if (out_cp)  surf%cpint = 0.0
		if (out_ints) surf%intsca = 0.0
	  end subroutine alloc_columns
	
	  subroutine dealloc_columns (surf)
	  	type (columns), intent(inout) :: surf
		
		if (out_lwp) deallocate ( surf%lwp )
		if (out_lwp) deallocate ( surf%iwp )
		if (out_cwp) deallocate ( surf%cwp )
		if (out_cwp) deallocate ( surf%rwp )		
		if (out_wvp) deallocate ( surf%wvp )
		if (out_wvp) deallocate ( surf%wvpt )
		if (out_wvp) deallocate ( surf%wvpb )
		if (out_cmse) deallocate ( surf%cmse )
		if (out_cmse) deallocate ( surf%cmset )
		if (out_cmse) deallocate ( surf%cmseb )
		if (out_cmfl) deallocate ( surf%cmfl )
		if (out_cape) deallocate ( surf%cape )
		if (out_cape) deallocate ( surf%cin )
		if (out_cape) deallocate ( surf%lcl )
		if (out_cape) deallocate ( surf%lfc )
		if (out_ctop) deallocate ( surf%ctop )
		if (out_ctop) deallocate ( surf%cbas )
		if (out_zinv) deallocate ( surf%bltop )
		if (out_zinv) deallocate ( surf%zinv )
		if (out_cp)  deallocate ( surf%cpint0 )
		if (out_cp)  deallocate ( surf%cpint )
		if (out_ints) deallocate ( surf%intsca )
	  end subroutine dealloc_columns
	
	  subroutine const_to_surf (surf, const)
	  	type (surface), intent(inout) :: surf
	  	real, intent(in) :: const
		
		surf%mflux   = const
		surf%esflux   = const
		surf%qvflux   = const
		if (out_rrate) surf%cumul    = const
		if (out_rrate) surf%precip    = const
		if (out_rrate) surf%cond    = const
		if (out_rrate) surf%evap    = const
		if (out_srad) surf%lwf = const
		if (out_srad) surf%swf = const
		if (out_olr) surf%olr = const
		if (out_olr) surf%toa = const
		if (out_osr) surf%osr = const	
		if (out_osr) surf%isr = const
		if (out_ocs) surf%olr_cs = const	
		if (out_ocs) surf%osr_cs = const		
	  return
	  end subroutine const_to_surf

	  function surf_add (data1, data2)
	  	type (surface) :: surf_add
	  	type (surface), intent(in) :: data1, data2
		
		surf_add%mflux   = data1%mflux    + data2%mflux
		surf_add%esflux  = data1%esflux    + data2%esflux
		surf_add%qvflux  = data1%qvflux    + data2%qvflux
		if (out_rrate) surf_add%cumul   = data1%cumul    + data2%cumul
		if (out_rrate) surf_add%precip    = data1%precip      + data2%precip
		if (out_rrate) surf_add%cond    = data1%cond      + data2%cond
		if (out_rrate) surf_add%evap    = data1%evap      + data2%evap
		if (out_srad) surf_add%lwf     = data1%lwf    + data2%lwf
		if (out_srad) surf_add%swf     = data1%swf    + data2%swf
		if (out_olr) surf_add%olr     = data1%olr    + data2%olr
		if (out_olr) surf_add%toa     = data1%toa    + data2%toa
		if (out_osr) surf_add%osr     = data1%osr    + data2%osr		
		if (out_osr) surf_add%isr     = data1%isr    + data2%isr
		if (out_ocs) surf_add%olr_cs  = data1%olr_cs  + data2%olr_cs
		if (out_ocs) surf_add%osr_cs  = data1%osr_cs  + data2%osr_cs
	  end function surf_add
	  	  	    			
	end module typedef_surf

