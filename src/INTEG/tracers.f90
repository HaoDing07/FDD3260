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
!  Purpose:
!
!  Author:
!	Julien Savre, LMU
!
! ================================================================

  module tracersmod
  
  USE shared_all
  USE advection
  USE sources
  USE subgrid
  USE funcpack
!
  IMPLICIT NONE
!
  private
!
  public :: tracers

  CONTAINS

!	===================================================
  subroutine tracers			     
!	===================================================
!    
  integer :: k, h
  real :: dref0(nz)
  real :: flux_s(ip_start:ip_end,jp_start:jp_end)
!
  if (verbose > 0) call write_debug('Entering tracers')
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!
!
!  Reference soundings
!
#if (defined CONSERVATIVE)
  dref0 = den0
#else
  dref0 = 1.
#endif
!
!  Loop over scalars 
!
  do h = 1, nscal
!
!----------------------------------------------------------!
!                        Sources                           !
!----------------------------------------------------------!
!
!  Turbulent diffusion
!
    if ( with_dif ) then 
      call scalar_flux ( h, state2%scal(:,:,:,h), flux_s )
!
      call effs ( 2+2*nhydro+h, pressure2%dens, state2%scal(:,:,:,h), states%scal(:,:,:,h), flux_s )
    endif
!
!  Nudging
!
    if ( with_nudg ) &
        call nudge_x ( 'sca', dref0*scal0(:,h), state%scal(:,:,:,h), states%scal(:,:,:,h) )
!
!  Large-scale subsidence and upwelling
!
    if ( with_lssub ) &
    	call advection_ac_1d ( h, w0, state%scal(:,:,:,h), states%scal(:,:,:,h) )
!
!  Other sources
!
    call scalar_source ( h, zbl, state%scal(:,:,:,h), states%scal(:,:,:,h) )  
!
!----------------------------------------------------------!
!    		      Scalar advection	  		   !
!----------------------------------------------------------!
!
    if ( with_adv ) &
        call advection_ac ( wind_adv, state%scal(:,:,:,h), states%scal(:,:,:,h), iscal=h )
!
!----------------------------------------------------------!
!               	  Terminate   	                   !
!----------------------------------------------------------!
!
    if (ldiag.and.out_diags) then
      diag(9)%sca(:,:,:,h) = diag(9)%sca(:,:,:,h) + cdiag*states%scal(:,:,:,h)
    endif
!
!  End loop
!
  enddo
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating tracers')
!
return                  
end

end module tracersmod
