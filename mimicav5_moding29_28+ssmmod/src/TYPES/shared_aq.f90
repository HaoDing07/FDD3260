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
!  SHARED_HYDRO.F:                   
!
!  Purpose:
!	Module containing hydrometeor variables (replaces hydro12.h) 
!
!  Author
!	Julien Savre, MISU
!
! ================================================================

	module shared_aq

	use gridno
	use shared_data
	use typedef_aq
		
	SAVE
!
	type (aq_chemical), dimension(:,:,:), pointer :: aqc2, aqr2
     		
	end module shared_aq
