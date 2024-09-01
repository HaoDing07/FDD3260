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
!  SHARED_SECOND.F:                   
!
!  Purpose:
!	Module containing secondary variables		  
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================

	module shared_turbu

	use gridno
	use shared_data
	use typedef_turbu
		
	SAVE
	
	!
	!  Turbulence model parameters
	!
	real, parameter :: C_s=0.23
	real, parameter :: C_k=0.12
        real, parameter :: C_1=2.34
	real, parameter :: C_e=0.85
        real, parameter :: C_n=0.76
        real, parameter :: Ric=0.25
	real, parameter :: vk=0.41	
	real, parameter :: c_mix=0.1
	
	type (turbulence) :: turbu, turbu_mean
	type (turbudiag) :: turbu_diag
	
	end module shared_turbu
