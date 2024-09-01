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
!  TYPEDEF_LAGRANGE.F:                   
!
!  Purpose:
!	Type definition for parcel tracking in MIMICa	  
!
!  Author
!	Julien Savre, MISU
!		
! ================================================================
	
	module typedef_lagrange
	
	USE gridno
        USE shared_data	

	integer :: nparc
	integer :: nparc_loc
	integer, dimension(:), allocatable :: parcel_loc
	
	type :: position
	  real    :: x, y, z
	  integer :: i, j, k
	  integer :: iproc
	  logical :: exist
	end type position

	type :: scalar
	  real :: u, v, w
	  real :: uprim, vprim, wprim
	  real :: p, t, pt, qt
	  real :: rh, ql, qi
	  real :: buoy, tke
	  real, dimension(nscal) :: scal
          real, dimension(ndiag) :: dw, dq
          integer :: mask
	end type scalar
	
	type :: aerosol
	  real :: r0			! Initial aerosol radius
	  real :: r				! Radius of wet solute (aerosol / droplet)
	  real :: m_s			! Mass of solute (wet aerosol / droplet)
	  real :: n_s			! Amount of moles of solute dry aerosol
	  real :: ws			! Sedimentation velocity of aerosol / droplet
	end type

	interface alloc_parc
		module procedure allocate_parcsca, allocate_parcpos
	end interface alloc_parc

	CONTAINS

	subroutine allocate_parcsca ( parcsca )
	  type(scalar), dimension(:), allocatable, intent(inout) :: parcsca
	  integer :: h

	  allocate( parcsca(1:nparc) )
	  do h = 1, nparc
	    parcsca(h)%u = 0.
	    parcsca(h)%v = 0.
	    parcsca(h)%w = 0.
	    parcsca(h)%uprim = 0.
	    parcsca(h)%vprim = 0.
	    parcsca(h)%wprim = 0.
	    parcsca(h)%p = 0.
	    parcsca(h)%t = 0.
	    parcsca(h)%pt = 0.
	    parcsca(h)%qt = 0.
	    parcsca(h)%rh = 0.
	    parcsca(h)%ql = 0.
	    parcsca(h)%qi = 0.
	    parcsca(h)%buoy = 0.
	    parcsca(h)%tke = 0.
	    parcsca(h)%scal(1:nscal) = 0.
	    parcsca(h)%dw(1:ndiag) = 0.
	    parcsca(h)%dq(1:ndiag) = 0.
	    parcsca(h)%mask = 0
	  enddo  
	  return
	end subroutine allocate_parcsca

	subroutine allocate_parcpos ( parcpos )
	  type(position), dimension(:), allocatable, intent(inout) :: parcpos
	  integer :: h

	  allocate( parcpos(1:nparc) )
	  do h = 1, nparc
	    parcpos(h)%x = 0.
	    parcpos(h)%y = 0.
	    parcpos(h)%z = 0.
	    parcpos(h)%i = 0
	    parcpos(h)%j = 0
	    parcpos(h)%k = 0
	    parcpos(h)%iproc = 0
	    parcpos(h)%exist = .false.
	  enddo
	  return
	end subroutine allocate_parcpos

	end module typedef_lagrange

