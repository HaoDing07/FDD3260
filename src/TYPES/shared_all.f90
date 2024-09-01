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

	module shared_all

	use gridno
	use shared_data
  	USE shared_diag
	USE shared_tend
  	USE shared_state
  	USE shared_pressure
  	USE shared_wind
  	USE shared_hydro
  	USE shared_thermo
	USE shared_turbu
	USE shared_rad
	USE shared_surf
	
  	USE shared_aerosol_new
  	USE shared_gas
  	USE shared_aq
  	USE shared_solid
	USE shared_nuclei

	end module shared_all
