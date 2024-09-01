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
!  TYPEDEF_SVALUE.F:                   
!
!  Purpose:
!	Typedef for derived-type data: svalue.			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
! ================================================================
	
	module typedef_svalue
	
	type :: s_value
		! ---
		! --- A scalar value with its 3D position
		! ---
		real 		    :: value	! value
		integer 	    :: i	! i-index
		integer 	    :: j	! j-index
		integer 	    :: k	! k-index
	end type s_value

	interface assignment (=)
		module procedure s2s
	end interface

	CONTAINS
	
	  subroutine s2s (s1, s2)
	       type (s_value), intent(inout) :: s1
	       type (s_value), intent(in)    :: s2
	       
	       s1%value = s2%value
	       s1%i	= s2%i
	       s1%j	= s2%j
	       s1%k	= s2%k
	  return
	  end subroutine s2s

	end module typedef_svalue


