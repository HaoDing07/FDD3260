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
!  ADVCM.F90
!
!  Purpose:
!	A package for calculating advection and turbulent mixing of
!       chemicals and aerosols
!
!  Units:
!	Gaseous chemicals - ppb(m) -> 1.e-9 kg/m^3air -> ppb(m)
!	Aqueous chemicals - moles/Lair			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
! ================================================================
!
!  ==========================================
    subroutine advc_m ( im )
!  ==========================================
!
! --- Note the units of aq chemicals are mole/Lair,
!	the units of gases and solid species are ppb(m)
!	the conversion is mole/Lair*Mw/den*1.e9 -> ppb(m)
!
USE shared_all
USE aerosols
!
IMPLICIT NONE
!
  integer :: k, h, l, im, iel, flag
  real :: flux_0, flux_n, flux_m          ! turbulent surface flux
  real, dimension(1:nz) :: zero
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz)  :: xs, scal
!
#ifdef NUC_CNT
  type(nuclei_3d), dimension(3) :: nucin	
#endif
!
!----------------------------------------------------------!
!                     Initializations		           !
!----------------------------------------------------------!
!
  flag = 0
  zero = 0.0
  xs(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
  flux_0 = 0.   ! turbulent surface flux; 2d initialization of zeros
!
!----------------------------------------------------------!
!           	      Solve for aerosols                   !
!----------------------------------------------------------!
!
#ifdef AERO_ENABLE
!
!  First calculate aerosol physics
!
    if ( aero_flg%chem) call aerosol_calc ( aero3d )
!
    do l = 1, nmode
      ! Define aerosol surface source  
      if (aero_sfc_source) then
        flux_n = -aero(l)%init%n_sfc_source
        flux_m = flux_n * aero(l)%size%mmean
      else
        flux_n = flux_0
        flux_m = flux_0
      end if
!
      call advance_ac ( im, l, aero3d2(l)%n, aero3d(l)%n, aero3ds(l)%n, atend(l)%n, flux_n )
!
      call advance_ac ( im, 0, aero3d2(l)%m, aero3d(l)%m, aero3ds(l)%m, atend(l)%m, flux_m )
!
      call advance_ac ( im, 100+l, aero3d2(l)%ma, aero3d(l)%ma, aero3ds(l)%ma, atend(l)%ma, flux_0)
    enddo
#else
    call advance_ac ( im, 0, nuc2%ccn, nuc%ccn, xs, ntend%ccn, flux_0)
!
    call advance_ac ( im, 0, nuc2%in, nuc%in, xs, ntend%in, flux_0)
#endif
!
!  INs
!
#ifdef NUC_CNT
    do h = 1, 3
       do l = 1, 3
	  call advance_ac ( im, flag, nucin2(h)%mode(l)%n, nucin(h)%mode(l)%n, nucins(h)%mode(l)%n, itend(h)%mode(l)%n, flux_0 )
!
	  call advance_ac ( im, flag, nucin2(h)%mode(l)%nc, nucin(h)%mode(l)%nc, nucins(h)%mode(l)%nc, itend(h)%mode(l)%nc, flux_0 )
!
#ifdef NUC_CNT1
	  call advance_ac ( im, flag, nucin2(h)%mode(l)%m, nucin(h)%mode(l)%m, nucins(h)%mode(l)%m, itend(h)%mode(l)%m, flux_0 )
!
	  call advance_ac ( im, flag, nucin2(h)%mode(l)%mc, nucin(h)%mode(l)%mc, nucins(h)%mode(l)%mc, itend(h)%mode(l)%mc, flux_0 )
#endif
       enddo
    enddo
#endif
!
!----------------------------------------------------------!
!           	      Solve for chemicals                  !
!----------------------------------------------------------!
!
#ifdef CHEM_ENABLE
!
!  Gaseous chemistry
!
    call advance_ac ( im, flag, gas2%o3, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%co, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%h2o2, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%xno, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%xno2, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%xn2o5, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%hno3, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%ch2o, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%ch3o2h, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%so2, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%h2so4, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%dms, xs, flux_0)
!    
    call advance_ac ( im, flag, gas2%nh3, xs, flux_0)
#endif
!
#ifdef AQCHEM_ENABLE
!
!  Cloud droplets
!    
    call advance_ac ( im, flag, aqc2%o3, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqc2%civ, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqc2%h2o2, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqc2%xnv, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqc2%ch2o, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqc2%ch3o2h, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqc2%siv, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqc2%svi, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqc2%nh4, xs, flux_0 )
!
!  Rain drops
!  
    call advance_ac ( im, flag, aqr2%o3, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqr2%civ, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqr2%h2o2, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqr2%xnv, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqr2%ch2o, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqr2%ch3o2h, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqr2%siv, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqr2%svi, xs, flux_0 )
!    
    call advance_ac ( im, flag, aqr2%nh4, xs, flux_0 )
#endif
!
#ifdef SOLIDCHEM_ENABLE
!
!  Ice particles
!
    call advance_ac ( im, flag, solidi2%o3, xs, flux_0 )
!
    call advance_ac ( im, flag, solidi2%h2o2, xs, flux_0 )
!
    call advance_ac ( im, flag, solidi2%xnv, xs, flux_0 )
!
    call advance_ac ( im, flag, solidi2%ch2o, xs, flux_0 )
!
    call advance_ac ( im, flag, solidi2%ch3o2h, xs, flux_0 )
!
    call advance_ac ( im, flag, solidi2%siv, xs, flux_0 )
!
    call advance_ac ( im, flag, solidi2%svi, xs, flux_0 )
!
    call advance_ac ( im, flag, solidi2%nh4, xs, flux_0 )
#endif
!
!----------------------------------------------------------!
!
    return
    end 
  
!  ==========================================
    subroutine advance_ac ( im, flag, x2, xx, xs, tendx, flux_s )
!  ==========================================
  
    USE gridno
    USE shared_data
    USE shared_pressure
    USE shared_diag
    USE advection
    USE subgrid
    USE funcpack
    USE boundary_conditions
    USE integpack
  
    IMPLICIT NONE
    
    logical :: adv
    integer :: flag, im
    real, intent(in) :: flux_s
    real, dimension(ip_start:ip_end,jp_start:jp_end) :: flux
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x2, xx, xs, tendx  
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: u, v, w
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: xtmp, dref
!  
!----------------------------------------------------------!
!
    xtmp(ip_start:ip_end,jp_start:jp_end,1:nz) = xx(ip_start:ip_end,jp_start:jp_end,1:nz)
!
#ifdef CONSERVATIVE
    dref = pressure2%dens
#else
    dref = 1.
#endif
!
    flux(ip_start:ip_end,jp_start:jp_end) = flux_s
!
!  Turbulent diffusion
!
    if ( with_dif ) call effs ( 100+flag, pressure2%dens, x2, xs, flux )
!
!  Advection
!
    if (with_adv) call advection_ac ( wind_adv, xtmp, xs, iaero=100+flag )
!
! Diagnostics
! 
    if ( out_diaga .and. flag>0 .and. flag<100 ) diag(9)%aero(flag)%n = diag(9)%aero(flag)%n + cdiag*xs
    if ( out_diaga .and. flag>100 ) diag(9)%aero(flag-100)%m = diag(9)%aero(flag-100)%m + cdiag*xs
!
!  Update scalar
!
    call integ_other ( im, dt0, dref*x2, xtmp, xs, tendx )
!
!  Boundary Conditions
!
    call statebc ( xtmp )     
!  
    xx(ip_start:ip_end,jp_start:jp_end,1:nz) = xtmp(ip_start:ip_end,jp_start:jp_end,1:nz)
!
    xs(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.0
!
!----------------------------------------------------------!
!
    return
    end subroutine advance_ac     
