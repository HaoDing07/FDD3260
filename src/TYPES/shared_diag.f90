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
!  SHARED_DIAG.F:                   
!
!  Purpose:
!	Module containing diagnostic variables
!
!	The dimension of diag is ndiag. Data are stored
!	as follows in diag:
!	- 1. Total tendency
!	- 2. Horizontal advection
!	- 3. Vertical advection
!	- 4. Horizontal mixing
!	- 5. Vertical mixing
!	- 6. Pressure term (winds) or radiation (scalars)
!	- 7. Buoyancy (wind) or precipitation (scalars)
!	- 8. Cond/evap or Dep/subl. (no wind component)
!	- from 10. Microphysics terms ....
!       - 15. CCN activation and ice nucleation
!
!  Author
!	Julien Savre, MISU
! ================================================================

	module shared_diag

	use gridno
	use shared_data
	USE shared_state
	USE shared_hydro
	USE shared_wind
	USE shared_thermo
	use typedef_diag
	
	public :: get_diag_diag
		
	SAVE
	
	type (tendencies), dimension(:), allocatable :: diag, diag_mean
	type (diagnosvar) :: diagnos
        type (entrainment) :: entrain

	CONTAINS
	
	subroutine get_diag_diag ( diag, pt, qt, hydromtrl, windl )
	      real, dimension(ip_start:,jp_start:,:), intent(in) :: pt, qt
	      type(hydrometeor), dimension(:), intent(in) :: hydromtrl
	      type(atm_winds), intent(in) :: windl
	      type(tendencies), dimension(:), intent(inout) :: diag
	      integer :: i, j, k, h, l
	      real :: sumh, qv, dqv
	
              if (out_diagtv) then
              do h = 1, ndiag
                do k = 1, nz
                  do j = jp_start, jp_end
                    do i = ip_start, ip_end
        	      qv = qt(i,j,k)
		      sumh = 0.
        	      do l = 1, nhydro
	  		qv = qv - hydromtrl(l)%q(i,j,k)
			sumh = sumh + hydromtrl(l)%q(i,j,k)
		      enddo
!
		      dqv = diag(h)%qt(i,j,k) - diag(h)%qc(i,j,k) - diag(h)%qr(i,j,k)
#ifdef SEIFERT
		      if (out_diagi) dqv = dqv - diag(h)%qi(i,j,k)
#endif
!
	              diag(h)%ptv(i,j,k) = (1. + epsm*qv - sumh)*diag(h)%pt(i,j,k) + epsm*pt(i,j,k)*dqv      &
		         - pt(i,j,k)*(diag(h)%qc(i,j,k) + diag(h)%qr(i,j,k))
#ifdef SEIFERT
		      if (out_diagi) diag(h)%ptv(i,j,k) = diag(h)%ptv(i,j,k) - pt(i,j,k)*diag(h)%qi(i,j,k)
#endif
	            enddo
		  enddo
	        enddo
	      enddo
	      endif

	      if (out_diagk) then
	        do h = 1, ndiag
	          diag(h)%k = 0.
	          if (out_diagu) diag(h)%k = diag(h)%k + windl%u*diag(h)%u
#ifdef MODEL_3D
	          if (out_diagv) diag(h)%k = diag(h)%k + windl%v*diag(h)%v
#endif
	          if (out_diagw) diag(h)%k = diag(h)%k + windl%w*diag(h)%w
	        enddo
	      endif
	        
	end subroutine
	
	end module shared_diag
