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
!  INTEGPACK:
!	Module containing routines used for time integration
!
!  Author:
!	Julien Savre, MISU, Stockholm
!
! ================================================================
!
  module integpack
!
  USE gridno
  USE shared_data
  USE shared_wind
  USE shared_state
  USE shared_hydro
  USE shared_thermo
  USE shared_pressure
  USE shared_tend
  USE averages

  IMPLICIT NONE
  
  private
  
  public :: integ_mom, integ_scal, integ_other, tridiff
  
  contains

  subroutine integ_mom ( im, dtl, windold, windl, winda )
!
!----------------------------------------------------------!
!                                                          !
!  Integrate winds using either Euler explicit (default)   !
!  or 3rd order Adams-Bashforth-Moulton predictor-corrector!
!  scheme (predictor: 2nd order explicit Adams-Bashforth,  !
!  corrector: 3rd order implicit Adams-Moulton)            !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  integer  :: im
  real :: dtl, k, H, tauw(nz), taud=0.25
!
  type (atm_winds), intent(in) :: windold, winda
  type (atm_winds), intent(inout) :: windl
!
#include "integ.h"
!
!----------------------------------------------------------!
!              Integrate velocities: Leapfrog              !
!----------------------------------------------------------!
!
  if (ntime == 1) then
    k = 1.
  else
    k = 2.
  endif
!
!  Implicit sponge layer
!
  tauw = 1.
  if (sponge == 2) then
    H = zmax - zdamp
    do k = 1, nz
      if (z0(k) >= zdamp) then
        tauw = taud * max(min(cos(0.5*pi*(zmax - z0(k))/H)**2.,1.),0.)
      endif
    enddo
    tauw = 1. / (1. + a_lm*dtl*tauw)
  endif
!
!  Integration
!
  if ( sorder > 1 ) then
    tend%u = winda%u + b_lm*tend%u
#ifdef MODEL_3D
    tend%v = winda%v + b_lm*tend%v
#endif
    tend%w = winda%w + b_lm*tend%w
  endif
!
!  Integrate u
!
  if ( sorder > 1 ) then
    windl%u = windold%u + a_lm*dtl*tend%u
!
#ifdef MODEL_3D
    windl%v = windold%v + a_lm*dtl*tend%v
#endif
!
    do k = 1, nz
      windl%w(:,:,k) = tauw(k) * (windold%w(:,:,k) + a_lm*dtl*tend%w(:,:,k))
    enddo
  else
    windl%u = windold%u + dtl*winda%u  
!
#ifdef MODEL_3D
     windl%v = windold%v + dtl*winda%v  
#endif
!
    do k = 1, nz
      windl%w(:,:,k) = tauw(k) * (windold%w(:,:,k) + dtl*winda%w(:,:,k))
    enddo 
  endif
!  
!----------------------------------------------------------!
!
  return
  end subroutine integ_mom
!  
  subroutine integ_scal (im, dtl, pressureold, stateold, hydromtrold,     &
 			 pressurel, statel, hydromtrl, 			  &
			 pressurea, statea, hydromtra )
!
!----------------------------------------------------------!
!                                                          !
!  Integrate scalars using either Euler explicit (default) !
!  or 3rd order Adams-Bashforth-Moulton predictor-corrector!
!  scheme (predictor: 2nd order explicit Adams-Bashforth,  !
!  corrector: 3rd order implicit Adams-Moulton)            !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  integer  :: im, h
  real  :: k, dtl, tot, tot_old
!
  type(atm_state), optional, intent(in) :: stateold, statea
  type(atm_state), optional, intent(inout) :: statel
  type(atm_pressure), optional, intent(in) :: pressureold, pressurea
  type(atm_pressure), optional, intent(inout) :: pressurel
  type(hydrometeor), optional, dimension(:), intent(in) :: hydromtrold, hydromtra
  type(hydrometeor), optional, dimension(:), intent(inout) :: hydromtrl
!
#include "integ.h"
!
!----------------------------------------------------------!
!              Integrate velocities: Leapfrog              !
!----------------------------------------------------------!
!
  if (ntime == 1) then
    k = 1.
  else
    k = 2.
  endif
!
  if ( sorder > 1 ) then
    if ( present(pressurea) ) tend%dens = pressurea%dens + b_lm*tend%dens
    if ( present(statea) ) tend%es = statea%es + b_lm*tend%es
    if ( present(statea) ) tend%qt = statea%qt + b_lm*tend%qt
    if ( present(statea) ) tend%scal = statea%scal + b_lm*tend%scal
    if ( present(hydromtra) ) then
      do h=1,nhydro
        htend(h)%q = hydromtra(h)%q + b_lm*htend(h)%q
        htend(h)%n = hydromtra(h)%n + b_lm*htend(h)%n
      enddo
    endif
  endif
!
!  Density (will be overwritten after pressure step)
!
  if (present(pressureold) .and. present(pressurea) .and. present(pressurel)) then
    if ( sorder > 1 ) then
      pressurel%dens = pressureold%dens + a_lm*dtl*tend%dens
    else
      pressurel%dens = pressureold%dens + dtl*pressurea%dens
    endif
  endif
!
!  Sensible energy
!
  if (present(stateold) .and. present(statea) .and. present(statel)) then
!
    if ( sorder > 1 ) then
      statel%es = stateold%es + a_lm*dtl*tend%es
    else
      statel%es = stateold%es + dtl*statea%es
    endif
!
!  Water vapor qt 
!
    if ( sorder > 1 ) then
      statel%qt = stateold%qt + a_lm*dtl*tend%qt
    else
      statel%qt = stateold%qt + dtl*statea%qt
    endif
!                                           
!  Passive scalars
!
    if ( sorder > 1 ) then
      do h=1,nscal 
        statel%scal(:,:,:,h) = stateold%scal(:,:,:,h) + a_lm*dtl*tend%scal(:,:,:,h)
      enddo
    else
      do h=1,nscal 
        statel%scal(:,:,:,h) = stateold%scal(:,:,:,h) + dtl*statea%scal(:,:,:,h)
      enddo
    endif
!
  endif
!                                           
!  Hydrometeors
!
  if (present(hydromtrold) .and. present(hydromtra) .and. present(hydromtrl)) then
    if ( sorder > 1 ) then
      do h=1,nhydro
        hydromtrl(h)%q = hydromtrold(h)%q + a_lm*dtl*htend(h)%q
!
#ifdef SEIFERT
        if (moments == 2) then
          hydromtrl(h)%n = hydromtrold(h)%n + a_lm*dtl*htend(h)%n
        endif
!
        if (h > ice .and. lmicro > 3) then
          hydromtrl(h)%w = hydromtrold(h)%w + a_lm*dtl*htend(h)%w
        endif
#endif
      end do
    else
      do h=1,nhydro
        hydromtrl(h)%q = hydromtrold(h)%q + dtl*hydromtra(h)%q
!
#ifdef SEIFERT
        if (moments == 2) then
          hydromtrl(h)%n = hydromtrold(h)%n + dtl*hydromtra(h)%n
        endif
!
        if (h > ice .and. lmicro > 3) then
          hydromtrl(h)%w = hydromtrold(h)%w + dtl*hydromtra(h)%w
        endif
#endif
      end do
    endif
  endif
!
!  Store total variations
!
  if (present(pressureold) .and. present(pressurel)) then
    call totalav (pressurel%dens, tot)
    call totalav (pressureold%dens, tot_old)
    drh_tot = (tot - tot_old)/dtl
  endif
!
  if (present(stateold) .and. present(statel)) then
    call totalav (statel%es, tot)
    call totalav (stateold%es, tot_old)
    des_tot = (tot - tot_old)/dtl
!
    call totalav (statel%qt, tot)
    call totalav (stateold%qt, tot_old)
    dqt_tot = (tot - tot_old)/dtl
  endif
!  
!----------------------------------------------------------!
!
  return
  end subroutine integ_scal
!
  subroutine integ_other ( im, dtl, xold, xl, xa, tendx )
!
!----------------------------------------------------------!
!                                                          !
!  Integrate winds using either Euler explicit (default)   !
!  or 3rd order Adams-Bashforth-Moulton predictor-corrector!
!  scheme (predictor: 2nd order explicit Adams-Bashforth,  !
!  corrector: 3rd order implicit Adams-Moulton)            !
!                                                          !
!----------------------------------------------------------!
!
!  Input/Output
!
  real :: dtl
  integer :: im
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(in) :: xold, xa
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz), intent(inout) :: xl, tendx
!
#include "integ.h"
!
!----------------------------------------------------------!
!              Integrate velocities: Leapfrog              !
!----------------------------------------------------------!
!
  if ( sorder > 1 ) then
    tendx = xa + b_lm*tendx
!
    xl = xold + a_lm*dtl*tendx
  else
    xl = xold + dtl*xa
  endif
!
!----------------------------------------------------------!
!
  return
  end subroutine integ_other
!
  subroutine tridiff(nh1,nh2,nz,nhdo,cix,ci,cip1,rhs,cj,cjp1)

    integer, intent(in) :: nh1,nh2,nz,nhdo
    real, intent(in)    :: cix(nh1:nh2,nz),ci(nh1:nh2,nz),cip1(nh1:nh2,nz)
    real, intent(inout) :: rhs(nh1:nh2,nz),cj(nh1:nh2,nz),cjp1(nh1:nh2,nz)

    integer k,i
    real eps

    eps=sqrt(tiny(1.))

    do i=nh1,nh2
       cjp1(i,2)=cip1(i,2) / ci(i,2)
       rhs(i,2)=rhs(i,2) / ci(i,2)
       do k=3,nz-1
          cj(i,k)=ci(i,k)-cix(i,k)*cjp1(i,k-1)
          cjp1(i,k)=cip1(i,k) / cj(i,k)
          rhs(i,k)=(rhs(i,k)-cix(i,k)*rhs(i,k-1)) / cj(i,k)
       enddo

       cj(i,nz-1) = rhs(i,nz-1)
       do k=nz-2,2,-1
          cj(i,k) = rhs(i,k) - sign(1.,cj(i,k+1)*cjp1(i,k))*max(abs(cj(i,k+1)*cjp1(i,k)),eps)
       enddo
    end do

  end subroutine tridiff
  
end module integpack
