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

      module radia
      
      save
      
      real, parameter :: cpr = 3.50174, rcp = 0.285572
      real, parameter :: mair = 28.967 , rowt = 1000., ep2 = 0.608
      real, parameter :: g=9.81, R=287.04
      real, parameter :: SolarConstant  = 1361.0
      real            :: totalpower
      integer  :: norig, nv, nv1, npts, mb, mbs, mbir
      	
      real, allocatable, dimension(:) ::  pp, pt, ph, po, plwc, piwc, pre, pde, 	&
           prwc, pgwc, fds, fus, fdir, fuir, taus, tauir

      end module radia
