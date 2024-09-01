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
!  SCALAR_SOURCES:
!	Subroutines calculating scalar sources in various different
!	ways: local sources, surface fluxes or other fluxes.
!
!  Author:
!	Julien Savre, Ludwig Maximilian Univeristat, Munich
!
! ================================================================	 

  module sources

  USE gridno
  USE shared_data
  USE shared_state
  USE shared_hydro
  USE shared_wind
  USE shared_pressure
  USE shared_thermo
  USE shared_diag
!
  IMPLICIT NONE
  
  private

  public :: scalar_reset, scalar_flux, scalar_source, other_source

  CONTAINS

subroutine scalar_reset ( h, scal, h_scal, reset )

! ----------------------------------------------
! --- Initialize a scalar by various means
! ----------------------------------------------
!
  integer :: h
  logical :: reset
  real :: h_scal
!
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: scal
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: scal_old
!
!-----------------------------------------------------------!
!
  if (reset) then
    scal_old = scal
!
    select case (sca_set)
!
    case (1)
      if (h == 3) then
        call scalar_reset_n (h, scal, h_scal, 3)
      else if (h == 4) then
        call scalar_reset_n (h, scal, h_scal, 4)
      else if (h == 5) then
        call scalar_reset_n (h, scal, h_scal, 5)
      endif
!
    case (2)
      if (h == 1) then
        call scalar_reset_n (h, scal, h_scal, 6)
      endif
!
    case (3)
      if (h == 1) then
        call scalar_reset_n (h, scal, h_scal, 1)
      else if (h == 2) then
        call scalar_reset_n (h, scal, h_scal, 1)
      endif
!
    case default
!
    end select
!
  endif
!
  return
!
contains
!
subroutine scalar_reset_n ( n, scal, h_scal, type )

! ----------------------------------------------
! --- Initialize a scalar by various means
! ----------------------------------------------
!
  integer :: n, type
  real :: h_scal
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: scal
!
  integer :: i, j, k
  real :: zz, epsw=0.1
!
!-----------------------------------------------------------!
!
  zz = 0.5*dz*fdz0(1)
  scal = 0.
!
  select case ( type )
!
!  Type 1: Initialize to cloud content (for entrainment)
!
  case (1)
!
    scal = hydromtr(drop)%q+hydromtr(rain)%q
!
!  Type 2: Initialize to unity below height h_scal
!
  case (2)
    do k = 1, nz-1
      if (zz <= h_scal) then	
        scal(ip_start:ip_end,jp_start:jp_end,k) = 1.
      endif
      zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
    enddo    
!
!  Type 3: Initialize to unity inside clouds (purity)
!
  case (3)
    do k = 1, nz-1
      do j=jt_start,jt_end
        do i=it_start,it_end
          if ( hydromtr(drop)%q(i,j,k) >= qthres .and. zz >= h_scal ) then	
            scal(i,j,k) = 1.
          endif
        enddo
      enddo
      zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
    enddo    
!
!  Type 4: Initialize to unity inside cloud cores
!
  case (4)
    do k = 1, nz-1
      do j=jt_start,jt_end
        do i=it_start,it_end
          if ( hydromtr(drop)%q(i,j,k) >= qthres .and. wind%w(i,j,k) >= wthres .and. zz >= h_scal ) then	
            scal(i,j,k) = 1.
          endif
        enddo
      enddo
      zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
    enddo    
!
!  Type 5: Initialize to unity inside cloud shell
!
  case (5)
    do k = 1, nz-1
      do j=jt_start,jt_end
        do i=it_start,it_end
          if ( hydromtr(drop)%q(i,j,k) >= qthres .and. wind%w(i,j,k) <= -wthres .and. zz >= h_scal ) then	
            scal(i,j,k) = 1.
          endif
        enddo
      enddo
      zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
    enddo    
!
!  Type 6: idealized Ozone mass fraction
!
  case (6)
    do k = 1, nz
      do j=jp_start,jp_end
        do i=ip_start,ip_end
          scal(i,j,k) = 3.6478 * (p0(k)/100.)**(0.83209) * exp(-p0(k)/1135.15)
        enddo
      enddo
    enddo
!
 end select
!
!-----------------------------------------------------------!
!
  return
  end
!
end subroutine scalar_reset
!
!
subroutine scalar_flux ( h, scal, flux )

! ----------------------------------------------
! --- Scalar source from surface fluxes
! ----------------------------------------------
!
  integer :: h
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: scal
  real, dimension(ip_start:,jp_start:), intent(inout) :: flux		
!
!-----------------------------------------------------------!
!
  flux = 0.
!
  if (h == 1 .or. h == 2) then
!    call scalar_flux_n (scal, flux, 1)
  elseif (h == 1000) then
    call scalar_flux_n (scal, flux, 2)
  endif
!
  return
!
contains
!
subroutine scalar_flux_n ( scal, flux, type )

! ----------------------------------------------
! --- Scalar source from surface fluxes
! ----------------------------------------------
!
  integer :: type
  integer :: i, j
  real :: windm
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: scal
  real, dimension(ip_start:,jp_start:), intent(inout) :: flux		
!
!-----------------------------------------------------------!
!
  select case ( type )
!
!  Type 1: Scalar flux from bulk surface model, with surface value of 1
!
  case (1)
    do j=jt_start,jt_end
      do i=it_start,it_end
        windm = ((wind%u(i+1,j,1)+wind%u(i,j,1))*0.5)**2.
#ifdef MODEL_3D
        windm = windm + ((wind%v(i,j+1,1)+wind%v(i,j,1))*0.5)**2.
#endif
!
        flux(i,j) = c_ds*sqrt(max(min_w,windm))*(scal(i,j,1) - 1.) / abs(pran)
      end do
    end do
!
!  Type 2: Fixed scalar flux
!
  case (2)
    flux(ip_start:ip_end,jp_start:jp_end) = scf
!
!  Type 3: Same a type 3 with surface value of 0
!
  case (3)
    do j=jt_start,jt_end
      do i=it_start,it_end
        windm = ((wind%u(i+1,j,1)+wind%u(i,j,1))*0.5)**2.
#ifdef MODEL_3D
        windm = windm + ((wind%v(i,j+1,1)+wind%v(i,j,1))*0.5)**2.
#endif
!
        flux(i,j) = c_ds*sqrt(max(min_w,windm))*(scal(i,j,1) - 0.) / abs(pran)
      end do
    end do
!
  end select
!
!-----------------------------------------------------------!
!
  return
  end
!
end subroutine scalar_flux
!
!
subroutine scalar_source ( h, h_scal, scal, scals )
!
  integer :: h
  real :: h_scal
!
  real, dimension(ip_start:,jp_start:,:), intent(in) :: scal
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: scals		
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: source		
!
!-----------------------------------------------------------!
!
  source = 0.0
!
  if ( h == 1 ) then
    call scalar_source_n (scal, h_scal, source, 3)
  else if ( h == 2 ) then
    call scalar_source_n (scal, h_scal, source, 2)
  endif
!
  scals = scals + source
!
!  Diagnostics
!
  if (out_diags) diag(7)%sca(:,:,:,h) = diag(7)%sca(:,:,:,h) + cdiag*source
!
  return
!
contains
!
subroutine scalar_source_n ( scal, h_scal, source, type )

! ----------------------------------------------
! --- Scalar sources 
! ----------------------------------------------
!
  integer :: type
  integer :: i, j, k
  real :: h_scal, zz
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: source		
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: scal
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: flu3d
  real, dimension(1:nz) :: flu, gflu
!
  real :: k4p
  real :: tau_dec1=1800., tau_dec2=120.
  real :: k2=3.91e5, k3=0.001, k4=1.e8, fo2=0.21, e=2060.
!
!-----------------------------------------------------------!
!
  select case ( type )
!
!  Type 1: Radioactive tracer
!
  case (1)
    source(ip_start:ip_end,jp_start:jp_end,1:nz) = source(ip_start:ip_end,jp_start:jp_end,1:nz)  &
    	- scal(ip_start:ip_end,jp_start:jp_end,1:nz) / tau_dec1
!
!  Type 2: Radioactive tracer in environment only
!
  case (2)
    zz = 0.5*dz*fdz0(1)
    do k = 1, nz-1
      do j = jt_start, jt_end
        do i = it_start, it_end
          if ( hydromtr(drop)%q(i,j,k) < qthres .or. wind%w(i,j,k) < wthres .and. zz >= h_scal ) 	&	
          	source(i,j,k) = source(i,j,k) - scal(i,j,k) / tau_dec2
        enddo
      enddo
      zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
    enddo    
!
!  Type 3: Simple Ozone model
!
  case (3)
    !zz = 0.5*dz*fdz0(1)
    !do k = 1, nz
    !  if ( zz <= 15000. ) then
    !    tau_dec2 = 10. + 90.*zz/15000.	
    !  else
    !    tau_dec2 = 100. - 99*(zz - 15000.)/30000.
    !  endif
    !  source(ip_start:ip_end,jp_start:jp_end,k) = source(ip_start:ip_end,jp_start:jp_end,k) -   &
    !  		(scal(ip_start:ip_end,jp_start:jp_end,k) - scal0(k,1)) / (86400.*tau_dec2)
    !  zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
    !enddo    
  end select
!
!-----------------------------------------------------------!
!
  return
  end
!
end subroutine scalar_source
!
!  ===================================================
   subroutine other_source ( type, dens, xs )
!  ===================================================
!
  integer :: type
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: xs, dens
!
  character(len=3) :: c3
  integer :: i, j, k
  real :: xx, yy, lx, ly, d
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: dref, source
!
  real :: h = 1000., r = 1500., dqv_surf = 2.e-6
!
!----------------------------------------------------------!
!              Time varying scalar sources                 !
!----------------------------------------------------------!
!
  source = 0.
!
#if (defined CONSERVATIVE)
  dref = dens
#else
  dref = 1.
#endif
!
  select case (casename)
!
    case ('KM-Sc')
      select case (type)
        case (1)
          do k = 1, nz
            source(:,:,k) = -3. *(pref/p0(k))**((cp_a-cv_a)/cp_a) / cp_a / 1000.
          enddo
        case (2)
          source = 3. / flv00 / 1000.
      end select
!
    case ('CP_IDEAL')
!
      if ( time < 1800. ) then
        lx = real(nx-5)*dx
        ly = real(ny-5)*dy
!
        select case (type)
!
!  Temperature
!	
	case (1)
	  c3 = 'pt '
          do k = 1, nz
            if ( z0(k) <= h ) then
              do j = jt_start, jt_end
	        do i = it_start, it_end
	          xx = real(i-4)*dx
	          yy = real(j-4)*dy
	          d = sqrt( (xx - (lx-dx)/2.)**2. + (yy - (ly-dy)/2.)**2. )
	          if ( d <= r ) source(i,j,k) = -dqv_surf*flv00/cp_a * dref(i,j,k) * (h - z0(k)) / h
	        enddo
	      enddo
	    endif
          enddo
!
!  Total moisture
!
	case (2)
	  c3 = 'qt '
          do k = 1, nz
            if ( z0(k) <= h ) then
              do j = jt_start, jt_end
	        do i = it_start, it_end
	          xx = real(i-4)*dx
	          yy = real(j-4)*dy
	          d = sqrt( (xx - (lx-dx)/2.)**2. + (yy - (ly-dy)/2.)**2. )
	          if ( d <= r ) source(i,j,k) = dqv_surf * dref(i,j,k) * (h - z0(k)) / h
	        enddo
	      enddo
	    endif
          enddo
!
!  Default
!
	case default
	  c3 = '   '
!
        end select
      endif
!  
  end select
!
  xs(it_start:it_end,jt_start:jt_end,1:nz) = xs(it_start:it_end,jt_start:jt_end,1:nz) + source(it_start:it_end,jt_start:jt_end,1:nz)
!
  if (out_diagt .and. type == 1) diag(8)%pt = diag(8)%pt + cdiag*source
  if (out_diagq .and. type == 2) diag(8)%qt = diag(8)%qt + cdiag*source
  if (out_diags .and. c3=='sca') diag(7)%sca(:,:,:,1) = diag(7)%sca(:,:,:,1) + cdiag*source
!
return
end subroutine other_source
!
end module sources
