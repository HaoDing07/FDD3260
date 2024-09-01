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
!  SHARED_STATE.F:                   
!
!  Purpose:
!	Module containing pressure variables
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================

	module shared_pressure

	use gridno
	use shared_data
	use typedef_pressure
		
	SAVE
	
    	real, dimension(nz) :: p0              ! reference hydro pressure
    	real, dimension(nz) :: p1              ! hydro pressure correction
    	real, dimension(nz) :: den0            ! reference hydro density
    	real, dimension(nz) :: stab0, g0	   ! reference hydro stability and lapse rate
    	real, dimension(nz) :: avden0          ! den0 between 2 vert. layers
	
	type (atm_pressure) :: pressure, pressures, pressure2, pressure_mean
	
	end module shared_pressure
