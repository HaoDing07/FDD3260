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
!  UD.F:
!
!  Purpose:
!	Scalars pt, qv, q, n: All scalar tendencies are calculated first
!	including advection, sgs mixing, damping and microphysics,
!	then scalars are integrated using either forward euler or RK.
!	As an input, all tendency terms already contain SGS mixing tendencies.
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================

  module columnmod
  
  USE shared_all
  USE advection
  USE micro
  USE sources
  USE funcpack
  USE columnphys
  USE radiationmod
!
  IMPLICIT NONE
!
  private
!
  public :: advance_column

  CONTAINS

!	===================================================
  subroutine advance_column			     
!	===================================================
!    
  integer :: k, h
  character(len=1) :: car1
  real :: dref0(nz), zero(nz)
!
  if (verbose > 0) call write_debug('Entering scalars')
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!
!
#if (defined CONSERVATIVE)
  dref0 = den0
#else
  dref0 = 1.
#endif
!
  zero = 0.
!
!----------------------------------------------------------!
!                        Sources                           !
!----------------------------------------------------------!
!
!  Nudging
!
  if (with_nudg) then
    call nudge_x ( 'pt ', dref0*pt0, state%es, states%es )
!
    call nudge_x ( 'qt ', dref0*qt0, state%qt, states%qt )
!
    do h = 1, nscal
        call nudge_x ( 'sca', dref0*scal0(:,h), state%scal(:,:,:,h), states%scal(:,:,:,h) )
    enddo
  endif
!
!  LS advection
!
  if (with_lsadv) then
    call ls_advec ( 1, state%es, states%es )
!
    call ls_advec ( 2, state%qt, states%qt )
  endif
!
!  Other sources
!
  if (with_lssrc) then
    call other_source ( 1, pressure%dens, states%es )
!
    call other_source ( 2, pressure%dens, states%qt )
!
    do h = 1, nscal
      call scalar_source ( h, zbl, state%scal(:,:,:,h), states%scal(:,:,:,h) )  
    enddo
  endif
!
!----------------------------------------------------------!
!                       Radiation                          !
!----------------------------------------------------------!
!
  if (with_rad) call rad_source ( pressure%dens, states%es )
!
!----------------------------------------------------------!
!                 Calculate PBL mixing                     !
!----------------------------------------------------------!
!
  if (with_dif) call pbl_mixing ( pressure%dens, zbl )
!
!----------------------------------------------------------!
!    		         Convection	  		   !
!----------------------------------------------------------!
!
  call convection ( pressure2, state2, hydromtr2 )
!
!----------------------------------------------------------!
!    		      Scalar advection	  		   !
!----------------------------------------------------------!
!
  if (with_adv) call advection_1d ( w0 )
!
!----------------------------------------------------------!
!	     Additional sources and relaxations            !
!----------------------------------------------------------!
!
  call damps ( 'pt', state%es, dref0*pt0, states%es )
!
  call damps ( 'qt', state%qt, dref0*qt0, states%qt )
!
  do h = 1, nhydro
    write(car1,'(i1)') h 
!  
    call damps('q'//car1, hydromtr(h)%q, zero, hydromtrs(h)%q)
!
    if (moments == 2) call damps('n'//car1, hydromtr(h)%n, zero, hydromtrs(h)%n)
  enddo
!
!----------------------------------------------------------!
!              	        Diagnostics   	                   !
!----------------------------------------------------------!
!  
  if (ldiag) then
    if (out_diagt) diag(9)%pt = states%es
    if (out_diagq) diag(9)%qt = states%qt
    do h = 1, nscal
      diag(9)%sca(:,:,:,h) = states%scal(:,:,:,h)
    enddo
!
    if (out_diagl.and.lmicro>0) diag(9)%qc = hydromtrs(drop)%q
    if (out_diagr.and.lmicro>0) diag(9)%qr = hydromtrs(rain)%q
    if (out_diagi.and.lmicro>1) diag(9)%qi = hydromtrs(ice)%q
  endif
!
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminating scalars')
!
return                  
end subroutine advance_column

end module columnmod
