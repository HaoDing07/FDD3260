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
!  CHEMSCAVENGE.F                   
!
!  Purpose:
!	Subroutines for calculating scavenging of chemicals
!       aqscaveng is for aqueous chemicals during coagulation/conversion
!       solscaveng is for solid chemicals during coagulation/conversion
!       evscaveng is for scavenging during evaporation (from liquid/ice to gas)
!       frscaveng is for scavenging during freezing (from liquid to ice)
!       aqprepscav is for precipitation scavenging by rain
!       solprepscav is for precipitation scavenging by ice
!
!  Author
!	Julien Savre
!       MISU
!
! ================================================================
!
!
!   ==========================================
    subroutine aqscaveng(i, j, k, hydromtrl,	    &
 			 auq, clr, cli, clg,        &
 			 xmir, cirr, crg,           &
			 xmsr, csrr, cls            &
#if ( defined AQCHEM_ENABLE ) || ( defined SOLIDCHEM_ENABLE )			 
     		     	 , gas                      &
#endif				 
#ifdef AQCHEM_ENABLE
			 , aqc, aqr                 &
#endif			 
#ifdef SOLIDCHEM_ENABLE
     		     	 , solidi                   &
#endif
                                               )
!				       
!    ==========================================
!
USE gridno
USE typedef_gas
USE typedef_aq
USE typedef_solid
USE shared_hydro
USE shared_data
!
IMPLICIT NONE 
!
  integer  :: i, j, k, h
!
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
!
#include "chemdef.h"
!
  real :: auq, clr, cli, clg, xmir, xmsr, cirr, csrr, cls, crg
  real :: qconv
!
#if ( defined AQCHEM_ENABLE ) && ( defined SOLIDCHEM_ENABLE )
!
!----------------------------------------------------------!
!           Scavenging of cloud drop chemicals             !
!----------------------------------------------------------!
!
  if (lmicro > 0) then
  if (hydromtrl(drop)%q(i,j,k) > qmin(drop)) then
!
!  Collection of cloud droplets 
!
    qconv = 1.0 - dt*(auq + clr + cli + clg + cls)/ hydromtr(drop)%q(i,j,k)
!
    if(qconv <= 0.0) qconv = 0.0
    aqc(i,j,k)%o3     = aqc%o3    *qconv
    aqc(i,j,k)%civ    = aqc(i,j,k)%civ   *qconv
    aqc(i,j,k)%h2o2   = aqc(i,j,k)%h2o2  *qconv
    aqc(i,j,k)%xnv    = aqc(i,j,k)%xnv   *qconv
    aqc(i,j,k)%ch2o   = aqc(i,j,k)%ch2o  *qconv
    aqc(i,j,k)%ch3o2h = aqc(i,j,k)%ch3o2h*qconv
    aqc(i,j,k)%siv    = aqc(i,j,k)%siv   *qconv
    aqc(i,j,k)%svi    = aqc(i,j,k)%svi   *qconv
    aqc(i,j,k)%nh4    = aqc(i,j,k)%nh4   *qconv
!
!  Evaporation of ice phase during rimming
!
    qconv = dt*(cli + clg + cls)/hydromtr(drop)%q(i,j,k)
!
    gas(i,j,k)%o3   = gas(i,j,k)%o3 + aqc(i,j,k)%o3  *qconv*escape_o3
    gas(i,j,k)%zco2 = gas(i,j,k)%zco2 + aqc(i,j,k)%civ *qconv*escape_civ
    gas(i,j,k)%h2o2 = gas(i,j,k)%h2o2 + aqc(i,j,k)%h2o2*qconv*escape_h2o2
    gas(i,j,k)%hno3 = gas(i,j,k)%hno3 + aqc(i,j,k)%xnv *qconv*escape_xnv
    gas(i,j,k)%ch2o = gas(i,j,k)%ch2o + aqc(i,j,k)%ch2o*qconv*escape_ch2o
    gas(i,j,k)%ch3o2h = gas(i,j,k)%ch3o2h + aqc(i,j,k)%ch3o2h*qconv*escape_ch3o2h
    gas(i,j,k)%so2  = gas(i,j,k)%so2 + aqc(i,j,k)%siv *qconv*escape_siv
    gas(i,j,k)%h2so4= gas(i,j,k)%h2so4 + aqc(i,j,k)%svi *qconv*escape_svi
    gas(i,j,k)%nh3 = gas(i,j,k)%nh3 + aqc(i,j,k)%nh4*qconv*escape_nh4
  end if
!
!----------------------------------------------------------!
!       Chemical transfer from ice and cloud to rain       !
!----------------------------------------------------------!
!
  qconv = 0.0
  if (hydromtrl(rain)%q(i,j,k) > qmin(rain) .and.       &
      hydromtrl(rain)%n(i,j,k) > xnmin(rain)) then
!
!  rain/cloud transformations 
!
    if (hydromtrl(drop)%q(i,j,k) > qmin(drop) .and.	&
        hydromtrl(drop)%n(i,j,k) > xnmin(drop)) then
      qconv = dt*(auq + clr)/hydromtrl(drop)%q(i,j,k)
!      
      if(qconv >= 1.0) qconv = 1.0
      aqr(i,j,k)%o3   = aqr(i,j,k)%o3 + aqc(i,j,k)%o3  *qconv
      aqr(i,j,k)%civ  = aqr(i,j,k)%civ + aqc(i,j,k)%civ *qconv
      aqr(i,j,k)%h2o2 = aqr(i,j,k)%h2o2 + aqc(i,j,k)%h2o2*qconv
      aqr(i,j,k)%xnv  = aqr(i,j,k)%xnv + aqc(i,j,k)%xnv *qconv
      aqr(i,j,k)%ch2o = aqr(i,j,k)%ch2o + aqc(i,j,k)%ch2o *qconv
      aqr(i,j,k)%ch3o2h = aqr(i,j,k)%ch3o2h + aqc(i,j,k)%ch3o2h *qconv
      aqr(i,j,k)%siv  = aqr(i,j,k)%siv + aqc(i,j,k)%siv *qconv
      aqr(i,j,k)%svi  = aqr(i,j,k)%svi + aqc(i,j,k)%svi *qconv
      aqr(i,j,k)%nh4  = aqr(i,j,k)%nh4 + aqc(i,j,k)%nh4*qconv
    endif
!
!  rain/ice transformations
!
    if (lmicro > 1) then
!
#ifdef SOLIDCHEM_ENABLE
    if (hydromtrl(ice)%q(i,j,k) > qmin(ice) .and.    &
        hydromtrl(ice)%n(i,j,k) > xnmin(ice)) then
      qconv = dt*xmir/hydromtrl(ice)%q(i,j,k)
!      
      if(qconv >= 1.0) qconv = 1.0
      aqr(i,j,k)%o3   = aqr(i,j,k)%o3 + solidi(i,j,k)%o3  *qconv
      aqr(i,j,k)%h2o2 = aqr(i,j,k)%h2o2 + solidi(i,j,k)%h2o2*qconv
      aqr(i,j,k)%xnv  = aqr(i,j,k)%xnv + solidi(i,j,k)%xnv *qconv
      aqr(i,j,k)%ch2o = aqr(i,j,k)%ch2o + solidi(i,j,k)%ch2o*qconv
      aqr(i,j,k)%ch3o2h = aqr(i,j,k)%ch3o2h + solidi(i,j,k)%ch3o2h *qconv
      aqr(i,j,k)%siv  = aqr(i,j,k)%siv + solidi(i,j,k)%siv *qconv
      aqr(i,j,k)%svi  = aqr(i,j,k)%svi + solidi(i,j,k)%svi *qconv
      aqr(i,j,k)%nh4  = aqr(i,j,k)%nh4 + solidi(i,j,k)%nh4 *qconv
    endif
    endif
#endif
!
!----------------------------------------------------------!
!    Scavenging of rain drops chemicals through rimming    !
!----------------------------------------------------------!
!
!  Evaporation during rimming: sinks of rain drops
!
    qconv = dt*(cirr + crg + csrr)/hydromtrl(rain)%q(i,j,k)
!
    gas(i,j,k)%o3   = gas(i,j,k)%o3 + aqr(i,j,k)%o3  *qconv*escape_o3
    gas(i,j,k)%zco2 = gas(i,j,k)%zco2 + aqr(i,j,k)%civ *qconv*escape_civ
    gas(i,j,k)%h2o2 = gas(i,j,k)%h2o2 + aqr(i,j,k)%h2o2*qconv*escape_h2o2
    gas(i,j,k)%hno3 = gas(i,j,k)%hno3 + aqr(i,j,k)%xnv *qconv*escape_xnv
    gas(i,j,k)%ch2o = gas(i,j,k)%ch2o + aqr(i,j,k)%ch2o*qconv*escape_ch2o
    gas(i,j,k)%ch3o2h = gas(i,j,k)%ch3o2h + aqr(i,j,k)%ch3o2h*qconv*escape_ch3o2h
    gas(i,j,k)%so2  = gas(i,j,k)%so2 + aqr(i,j,k)%siv *qconv*escape_siv
    gas(i,j,k)%h2so4= gas(i,j,k)%h2so4 + aqr(i,j,k)%svi *qconv*escape_svi
    gas(i,j,k)%nh3  = gas(i,j,k)%nh3 + aqr(i,j,k)%nh4*qconv*escape_nh4
    endif
!
!  Liquid phase (rain)
!
    qconv = 1.0 - qconv
    if(qconv <= 0.0) qconv = 0.0
    aqr(i,j,k)%o3     = aqr(i,j,k)%o3    *qconv
    aqr(i,j,k)%civ    = aqr(i,j,k)%civ   *qconv
    aqr(i,j,k)%h2o2   = aqr(i,j,k)%h2o2  *qconv
    aqr(i,j,k)%xnv    = aqr(i,j,k)%xnv   *qconv
    aqr(i,j,k)%ch2o   = aqr(i,j,k)%ch2o  *qconv
    aqr(i,j,k)%ch3o2h = aqr(i,j,k)%ch3o2h*qconv     
    aqr(i,j,k)%siv    = aqr(i,j,k)%siv   *qconv
    aqr(i,j,k)%svi    = aqr(i,j,k)%svi   *qconv     
    aqr(i,j,k)%nh4    = aqr(i,j,k)%nh4   *qconv
  endif
!  
  endif
#endif
!
!----------------------------------------------------------!
!
return
end
!  
!	  
!	==========================================
	subroutine solscaveng(i, j, k, hydromtrl,      &
	                      cli, aig, ais,          &
			      cir, cis, xmir          &
#if ( defined AQCHEM_ENABLE ) || ( defined SOLIDCHEM_ENABLE )
     		     	 , gas                        &
#endif
#ifdef AQCHEM_ENABLE
     		     	 , aqc, aqr                   &
#endif
#ifdef SOLIDCHEM_ENABLE
     		     	 , solidi                     &
#endif			 
			               )
!
!     	==========================================
!
USE gridno
USE typedef_aq
USE typedef_solid
USE typedef_gas
USE shared_hydro
USE shared_data
!
IMPLICIT NONE
!
  integer  :: i, j, k, h
!  
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
!
#include "chemdef.h"
!
  real :: cli, aig, ais, cir, cis, xmir
  real :: qconv
!
#if ( defined AQCHEM_ENABLE ) && ( defined SOLIDCHEM_ENABLE )
!
!----------------------------------------------------------!
!           Transfer of cloud chemicals to ice             !
!----------------------------------------------------------!
!
  if (lmicro > 1) then
!
  qconv = 0.0
  if (hydromtrl(ice)%q(i,j,k) > qmin(ice) .and.     &
      hydromtrl(ice)%n(i,j,k) > xnmin(ice)) then
!
!  ice/cloud interactions
!
#ifdef AQCHEM_ENABLE
    if (hydromtrl(drop)%q(i,j,k) > qmin(drop) .and.    &
        hydromtrl(drop)%n(i,j,k) > xnmin(drop)) then
      qconv = dt*cli/hydromtrl(drop)%q(i,j,k)
!
      if(qconv >= 1.0) qconv = 1.0
      solidi(i,j,k)%o3    = solidi(i,j,k)%o3 + aqc(i,j,k)%o3*qconv*retent_o3
      solidi(i,j,k)%h2o2  = solidi(i,j,k)%h2o2 + aqc(i,j,k)%h2o2*qconv*retent_h2o2
      solidi(i,j,k)%xnv   = solidi(i,j,k)%xnv + aqc(i,j,k)%xnv*qconv*retent_xnv
      solidi(i,j,k)%ch2o  = solidi(i,j,k)%ch2o + aqc(i,j,k)%ch2o*qconv*retent_ch2o
      solidi(i,j,k)%ch3o2h= solidi(i,j,k)%ch3o2h + aqc(i,j,k)%ch3o2h*qconv*retent_ch3o2h
      solidi(i,j,k)%siv   = solidi(i,j,k)%siv + aqc(i,j,k)%siv*qconv*retent_siv
      solidi(i,j,k)%svi   = solidi(i,j,k)%svi + aqc(i,j,k)%svi*qconv*retent_svi
      solidi(i,j,k)%nh4   = solidi(i,j,k)%nh4 + aqc(i,j,k)%nh4*qconv*retent_nh4
    endif
#endif
!
!----------------------------------------------------------!
!  Ice chemicals scavenging via melting/rimming/conversion !
!----------------------------------------------------------!
!
    qconv = 1.0 - dt*(aig + ais + cir + cis + xmir)/hydromtrl(ice)%q(i,j,k)
!
    if(qconv <= 0.0) qconv = 0.0
    solidi(i,j,k)%o3	= solidi(i,j,k)%o3	*qconv
    solidi(i,j,k)%h2o2  = solidi(i,j,k)%h2o2  *qconv
    solidi(i,j,k)%xnv   = solidi(i,j,k)%xnv	*qconv
    solidi(i,j,k)%ch2o  = solidi(i,j,k)%ch2o  *qconv
    solidi(i,j,k)%ch3o2h= solidi(i,j,k)%ch3o2h*qconv
    solidi(i,j,k)%siv   = solidi(i,j,k)%siv	*qconv
    solidi(i,j,k)%svi   = solidi(i,j,k)%svi	*qconv
    solidi(i,j,k)%nh4   = solidi(i,j,k)%nh4	*qconv
  endif
!  
  endif
#endif
!
!----------------------------------------------------------!
!
return
end
!
!
!   ==========================================
    subroutine evscaveng(i, j, k, hydromtrl, 	    &
 			 dqc, dqr, dqi              &
#if ( defined AQCHEM_ENABLE ) || ( defined SOLIDCHEM_ENABLE )			 
     		     	 , gas                      &
#endif				 
#ifdef AQCHEM_ENABLE
			 , aqc, aqr                 &
#endif			 
#ifdef SOLIDCHEM_ENABLE
     		     	 , solidi                   &
#endif
                                               )
!				       
!    ==========================================
!
USE gridno
USE typedef_gas
USE typedef_aq
USE typedef_solid
USE shared_hydro
USE shared_data
!
IMPLICIT NONE
!
  integer :: i, j, k, h
!
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
!
#include "chemdef.h"
!
  real :: dqc, dqr, dqi, qconv
!
#if ( defined AQCHEM_ENABLE ) && ( defined SOLIDCHEM_ENABLE )
!
!----------------------------------------------------------!
!            Scavenging of rain drop chemicals             !
!----------------------------------------------------------!
!
  if (lmicro > 0) then
!
  if(hydromtrl(rain)%q(i,j,k) > qmin(rain) .and. hydromtrl(rain)%n(i,j,k) > xnmin(rain))then
    qconv = - dqr/hydromtrl(rain)%q(i,j,k)*0.5
  else
    qconv = 1.0
  endif
!
  if(qconv >= 1.0) qconv = 1.0

  gas(i,j,k)%o3    = gas(i,j,k)%o3 + qconv*aqr(i,j,k)%o3
  gas(i,j,k)%zco2  = gas(i,j,k)%zco2 + qconv*aqr(i,j,k)%civ
  gas(i,j,k)%h2o2  = gas(i,j,k)%h2o2 + qconv*aqr(i,j,k)%h2o2
  gas(i,j,k)%hno3  = gas(i,j,k)%hno3 + qconv*aqr(i,j,k)%xnv
  gas(i,j,k)%ch2o  = gas(i,j,k)%ch2o + qconv*aqr(i,j,k)%ch2o
  gas(i,j,k)%ch3o2h= gas(i,j,k)%ch3o2h + qconv*aqr(i,j,k)%ch3o2h
  gas(i,j,k)%so2   = gas(i,j,k)%so2 + qconv*aqr(i,j,k)%siv
  gas(i,j,k)%h2so4 = gas(i,j,k)%h2so4 + qconv*aqr(i,j,k)%svi
  gas(i,j,k)%nh3  = gas(i,j,k)%nh3 + qconv*aqr(i,j,k)%nh4
!
  qconv = 1.0 - qconv
  aqr(i,j,k)%o3    = aqr(i,j,k)%o3    *qconv
  aqr(i,j,k)%civ   = aqr(i,j,k)%civ   *qconv
  aqr(i,j,k)%h2o2  = aqr(i,j,k)%h2o2  *qconv
  aqr(i,j,k)%xnv   = aqr(i,j,k)%xnv   *qconv
  aqr(i,j,k)%ch2o  = aqr(i,j,k)%ch2o  *qconv
  aqr(i,j,k)%ch3o2h= aqr(i,j,k)%ch3o2h*qconv
  aqr(i,j,k)%siv   = aqr(i,j,k)%siv   *qconv
  aqr(i,j,k)%svi   = aqr(i,j,k)%svi   *qconv
  aqr(i,j,k)%nh4   = aqr(i,j,k)%nh4  *qconv
!
!----------------------------------------------------------!
!            Scavenging of cloud drop chemicals            !
!----------------------------------------------------------!
!
  if(hydromtrl(drop)%q(i,j,k) > qmin(drop) .and. hydromtrl(drop)%n(i,j,k) > xnmin(drop))then
    qconv = - dqc/hydromtrl(drop)%q(i,j,k)*0.5
  else
    qconv = 1.0
  endif
!
  if (qconv >= 1.0) qconv = 1.0
!
  gas(i,j,k)%o3    = gas(i,j,k)%o3 + qconv*aqc(i,j,k)%o3
  gas(i,j,k)%zco2  = gas(i,j,k)%zco2 + qconv*aqc(i,j,k)%civ
  gas(i,j,k)%h2o2  = gas(i,j,k)%h2o2 + qconv*aqc(i,j,k)%h2o2
  gas(i,j,k)%hno3  = gas(i,j,k)%hno3 + qconv*aqc(i,j,k)%xnv
  gas(i,j,k)%ch2o  = gas(i,j,k)%ch2o + qconv*aqc(i,j,k)%ch2o
  gas(i,j,k)%ch3o2h= gas(i,j,k)%ch3o2h + qconv*aqc(i,j,k)%ch3o2h
  gas(i,j,k)%so2   = gas(i,j,k)%so2 + qconv*aqc(i,j,k)%siv
  gas(i,j,k)%h2so4 = gas(i,j,k)%h2so4 + qconv*aqc(i,j,k)%svi
  gas(i,j,k)%nh3  = gas(i,j,k)%nh3 + qconv*aqc(i,j,k)%nh4
!
  qconv = 1.0 - qconv
  aqc(i,j,k)%o3	   = aqc(i,j,k)%o3    *qconv
  aqc(i,j,k)%civ   = aqc(i,j,k)%civ   *qconv
  aqc(i,j,k)%h2o2  = aqc(i,j,k)%h2o2  *qconv
  aqc(i,j,k)%xnv   = aqc(i,j,k)%xnv   *qconv
  aqc(i,j,k)%ch2o  = aqc(i,j,k)%ch2o  *qconv
  aqc(i,j,k)%ch3o2h= aqc(i,j,k)%ch3o2h*qconv
  aqc(i,j,k)%siv   = aqc(i,j,k)%siv   *qconv
  aqc(i,j,k)%svi   = aqc(i,j,k)%svi   *qconv
  aqc(i,j,k)%nh4   = aqc(i,j,k)%nh4  *qconv
!  
  endif
!
!----------------------------------------------------------!
!               Scavenging of ice chemicals                !
!----------------------------------------------------------!
!
  if (lmicro > 1) then
!
  if(hydromtrl(ice)%q(i,j,k) > qmin(ice) .and. hydromtrl(ice)%n(i,j,k) > xnmin(ice))then
    qconv = - dqi/hydromtrl(ice)%q(i,j,k)*0.5
  else
    qconv = 1.0
  endif
  if(qconv.ge.1.0) qconv = 1.0

  gas(i,j,k)%o3    = gas(i,j,k)%o3 + qconv*solidi(i,j,k)%o3
  gas(i,j,k)%h2o2  = gas(i,j,k)%h2o2 + qconv*solidi(i,j,k)%h2o2
  gas(i,j,k)%hno3  = gas(i,j,k)%hno3 + qconv*solidi(i,j,k)%xnv
  gas(i,j,k)%ch2o  = gas(i,j,k)%ch2o + qconv*solidi(i,j,k)%ch2o
  gas(i,j,k)%ch3o2h= gas(i,j,k)%ch3o2h + qconv*solidi(i,j,k)%ch3o2h
  gas(i,j,k)%so2   = gas(i,j,k)%so2 + qconv*solidi(i,j,k)%siv
  gas(i,j,k)%h2so4 = gas(i,j,k)%h2so4 + qconv*solidi(i,j,k)%svi
  gas(i,j,k)%nh3   = gas(i,j,k)%nh3 + qconv*solidi(i,j,k)%nh4
!
  qconv = 1.0 - qconv
  solidi(i,j,k)%o3    = solidi(i,j,k)%o3    *qconv
  solidi(i,j,k)%h2o2  = solidi(i,j,k)%h2o2  *qconv
  solidi(i,j,k)%xnv   = solidi(i,j,k)%xnv   *qconv
  solidi(i,j,k)%ch2o  = solidi(i,j,k)%ch2o  *qconv
  solidi(i,j,k)%ch3o2h= solidi(i,j,k)%ch3o2h*qconv
  solidi(i,j,k)%siv   = solidi(i,j,k)%siv   *qconv
  solidi(i,j,k)%svi   = solidi(i,j,k)%svi   *qconv
  solidi(i,j,k)%nh4   = solidi(i,j,k)%nh4   *qconv
  endif
!
#endif
!
return
end
!
!
!   ==========================================
    subroutine frscaveng(i, j, k, hydromtrl,         &
 			 frg, fci                   &
#if ( defined AQCHEM_ENABLE ) || ( defined SOLIDCHEM_ENABLE )			 
     		     	 , gas                      &
#endif				 
#ifdef AQCHEM_ENABLE
			 , aqc, aqr                 &
#endif			 
#ifdef SOLIDCHEM_ENABLE
     		     	 , solidi                   &
#endif
                                               )
!				       
!    ==========================================
!
USE gridno
USE typedef_gas
USE typedef_aq
USE typedef_solid
USE shared_hydro
USE shared_data
!
IMPLICIT NONE
!
  integer :: i, j, k, h
!
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
!
#include "chemdef.h"
!
  real :: frg, fci, qconv
!
#if ( defined AQCHEM_ENABLE ) && ( defined SOLIDCHEM_ENABLE )
!
!----------------------------------------------------------!
!            Scavenging of rain drop chemicals             !
!----------------------------------------------------------!
!
  if (lmicro > 1) then
!
  if (hydromtrl(rain)%q(i,j,k) > qmin(rain) .and.     &
      hydromtrl(rain)%n(i,j,k) > xnmin(rain)) then
!
    qconv = dt2*frg/hydromtrl(rain)%q(i,j,k)
!
!  Solid (ice)
!
#ifdef SOLIDCHEM_ENABLE
    if(qconv <= 0.0) qconv = 0.0
    solidi(i,j,k)%o3	  = solidi(i,j,k)%o3 + aqr(i,j,k)%o3* retent_o3*qconv
    solidi(i,j,k)%h2o2   = solidi(i,j,k)%h2o2 + aqr(i,j,k)%h2o2*retent_h2o2*qconv
    solidi(i,j,k)%xnv    = solidi(i,j,k)%xnv + aqr(i,j,k)%xnv *retent_xnv*qconv
    solidi(i,j,k)%ch2o   = solidi(i,j,k)%ch2o + aqr(i,j,k)%ch2o *retent_ch2o*qconv
    solidi(i,j,k)%ch3o2h = solidi(i,j,k)%ch3o2h + aqr(i,j,k)%ch3o2h*retent_ch3o2h*qconv
    solidi(i,j,k)%siv    = solidi(i,j,k)%siv + aqr(i,j,k)%siv *retent_siv*qconv
    solidi(i,j,k)%svi    = solidi(i,j,k)%svi + aqr(i,j,k)%svi *retent_svi*qconv
    solidi(i,j,k)%nh4    = solidi(i,j,k)%nh4 + aqr(i,j,k)%nh4*retent_nh4*qconv
#endif
!
!  Liquid
!
#ifdef AQCHEM_ENABLE
    gas(i,j,k)%o3   = gas(i,j,k)%o3 + aqr(i,j,k)%o3  *qconv*escape_o3
    gas(i,j,k)%zco2 = gas(i,j,k)%zco2 + aqr(i,j,k)%civ *qconv*escape_civ
    gas(i,j,k)%h2o2 = gas(i,j,k)%h2o2 + aqr(i,j,k)%h2o2*qconv*escape_h2o2
    gas(i,j,k)%hno3 = gas(i,j,k)%hno3 + aqr(i,j,k)%xnv *qconv*escape_xnv
    gas(i,j,k)%ch2o = gas(i,j,k)%ch2o + aqr(i,j,k)%ch2o*qconv*escape_ch2o
    gas(i,j,k)%ch3o2h = gas(i,j,k)%ch3o2h + aqr(i,j,k)%ch3o2h*qconv*escape_ch3o2h
    gas(i,j,k)%so2  = gas(i,j,k)%so2 + aqr(i,j,k)%siv *qconv*escape_siv
    gas(i,j,k)%h2so4= gas(i,j,k)%h2so4 + aqr(i,j,k)%svi *qconv*escape_svi
    gas(i,j,k)%nh3  = gas(i,j,k)%nh3 + aqr(i,j,k)%nh4*qconv*escape_nh4
!
    if(qconv <= 0.0) qconv = 0.0
    aqr(i,j,k)%o3     = aqr(i,j,k)%o3    *(1.0 - qconv)
    aqr(i,j,k)%civ    = aqr(i,j,k)%civ   *(1.0 - qconv)
    aqr(i,j,k)%h2o2   = aqr(i,j,k)%h2o2  *(1.0 - qconv)
    aqr(i,j,k)%xnv    = aqr(i,j,k)%xnv   *(1.0 - qconv)
    aqr(i,j,k)%ch2o   = aqr(i,j,k)%ch2o  *(1.0 - qconv)
    aqr(i,j,k)%ch3o2h = aqr(i,j,k)%ch3o2h*(1.0 - qconv)     
    aqr(i,j,k)%siv    = aqr(i,j,k)%siv   *(1.0 - qconv)
    aqr(i,j,k)%svi    = aqr(i,j,k)%svi   *(1.0 - qconv)     
    aqr(i,j,k)%nh4    = aqr(i,j,k)%nh4   *(1.0 - qconv)
  endif
#endif
  endif
!
!----------------------------------------------------------!
!            Scavenging of cloud drop chemicals            !
!----------------------------------------------------------!
!
  if (hydromtrl(drop)%q(i,j,k) > qmin(drop) .and.     &
      hydromtrl(drop)%n(i,j,k) > xnmin(drop)) then
!
  qconv = dt2*fci/hydromtrl(drop)%q(i,j,k)
!
!  Solid
!
#ifdef SOLIDCHEM_ENABLE
  solidi(i,j,k)%o3	= solidi(i,j,k)%o3 + aqc(i,j,k)%o3* retent_o3*qconv
  solidi(i,j,k)%h2o2	= solidi(i,j,k)%h2o2 + aqc(i,j,k)%h2o2*retent_h2o2*qconv
  solidi(i,j,k)%xnv	= solidi(i,j,k)%xnv + aqc(i,j,k)%xnv *retent_xnv*qconv
  solidi(i,j,k)%ch2o	= solidi(i,j,k)%ch2o + aqc(i,j,k)%ch2o *retent_ch2o*qconv
  solidi(i,j,k)%ch3o2h = solidi(i,j,k)%ch3o2h + aqc(i,j,k)%ch3o2h*retent_ch3o2h*qconv
  solidi(i,j,k)%siv	= solidi(i,j,k)%siv + aqc(i,j,k)%siv *retent_siv*qconv
  solidi(i,j,k)%svi	= solidi(i,j,k)%svi + aqc(i,j,k)%svi *retent_svi*qconv
  solidi(i,j,k)%nh4	= solidi(i,j,k)%nh4 + aqc(i,j,k)%nh4*retent_nh4*qconv
#endif  
!
!  Liquid
!
#ifdef AQCHEM_ENABLE
  gas(i,j,k)%o3     = gas(i,j,k)%o3 + aqc(i,j,k)%o3    *escape_o3*qconv
  gas(i,j,k)%zco2   = gas(i,j,k)%zco2 + aqc(i,j,k)%civ   *escape_civ*qconv
  gas(i,j,k)%h2o2   = gas(i,j,k)%h2o2 + aqc(i,j,k)%h2o2  *escape_h2o2*qconv
  gas(i,j,k)%hno3   = gas(i,j,k)%hno3 + aqc(i,j,k)%xnv   *escape_xnv*qconv
  gas(i,j,k)%ch2o   = gas(i,j,k)%ch2o + aqc(i,j,k)%ch2o  *escape_ch2o*qconv
  gas(i,j,k)%ch3o2h = gas(i,j,k)%ch3o2h + aqc(i,j,k)%ch3o2h*escape_ch3o2h*qconv
  gas(i,j,k)%so2    = gas(i,j,k)%so2 + aqc(i,j,k)%siv  *escape_siv*qconv
  gas(i,j,k)%h2so4  = gas(i,j,k)%h2so4 + aqc(i,j,k)%svi   *escape_svi*qconv
  gas(i,j,k)%nh3    = gas(i,j,k)%nh3 + aqc(i,j,k)%nh4  *escape_nh4*qconv
!
  aqc(i,j,k)%o3     = aqc(i,j,k)%o3	 *(1.0 - qconv)
  aqc(i,j,k)%zco2   = aqc(i,j,k)%zco2  *(1.0 - qconv)
  aqc(i,j,k)%h2o2   = aqc(i,j,k)%h2o2  *(1.0 - qconv)
  aqc(i,j,k)%hno3   = aqc(i,j,k)%hno3  *(1.0 - qconv)
  aqc(i,j,k)%ch2o   = aqc(i,j,k)%ch2o  *(1.0 - qconv)
  aqc(i,j,k)%ch3o2h = aqc(i,j,k)%ch3o2h*(1.0 - qconv)
  aqc(i,j,k)%so2    = aqc(i,j,k)%so2   *(1.0 - qconv)
  aqc(i,j,k)%h2so4  = aqc(i,j,k)%h2so4 *(1.0 - qconv)
  aqc(i,j,k)%nh3    = aqc(i,j,k)%nh3   *(1.0 - qconv)
#endif
  endif
!  
  endif
#endif
!
return
end
!
!  ==================================================
!
   subroutine aqprecip   (pup, pdn, gas                 &
#ifdef SOLIDCHEM_ENABLE
                          , solidi                      &
#endif			  
			  , aqr, aqc)   
!
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Subroutines calculating precipitation scavenging of aqueous
!       phase chemicals
!
! ================================================================
!
USE gridno
USE typedef_aq
USE typedef_gas
USE typedef_solid
USE shared_hydro
USE shared_data
!
  integer :: i, j, k, h
  real  :: qmean_up, qmean_dn, cmean_up, cmean_dn,         &
           qconv_up, qconv_dn 
!	   
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz,1:nhydro)  :: pup, pdn
!
#include "chemdef.h"
!
#if ( defined AQCHEM_ENABLE )
!
!----------------------------------------------------------!
!      	
  if (lmicro > 0) then
!
  do k = 2, nzi
    do j = js, je
      do i = is, ie
!
	qmean_up = 0.5*(hydromtr2(rain)%q(i,j,k+1)+hydromtr2(rain)%q(i,j,k))
        if (qmean_up > qmin(rain)) then
          qconv_up = pup(i,j,k,rain)*dt/qmean_up
          if(qconv_up >= 1.0) qconv_up = 1.0
        else
          qconv_up = 0.0
        endif
!
	qmean_dn = 0.5*(hydromtr2(rain)%q(i,j,k-1)+hydromtr2(rain)%q(i,j,k))
        if (qmean_dn > qmin(rain)) then
          qconv_dn = pdn(i,j,k,rain)*dt/qmean_dn
          if(qconv_dn >= 1.0) qconv_dn = 1.0
        else
          qconv_dn = 0.0
        endif
!
        cmean_up = (aqr(i,j,k)%o3 + aqr(i,j,k)%o3)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%o3 + aqr(i,j,k)%o3)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%o3 = aqr(i,j,k)%o3 + qconv
        if(aqr(i,j,k)%o3.le.0.0) aqr(i,j,k)%o3 = 0.0
!
        cmean_up = (aqr(i,j,k)%civ + aqr(i,j,k)%civ)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%civ + aqr(i,j,k)%civ)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%civ = aqr(i,j,k)%civ + qconv
        if(aqr(i,j,k)%civ.le.0.0) aqr(i,j,k)%civ = 0.0
!
        cmean_up = (aqr(i,j,k)%h2o2 + aqr(i,j,k)%h2o2)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%h2o2 + aqr(i,j,k)%h2o2)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%h2o2 = aqr(i,j,k)%h2o2 + qconv
        if(aqr(i,j,k)%h2o2.le.0.0) aqr(i,j,k)%h2o2 = 0.0
!
        cmean_up = (aqr(i,j,k)%nh4 + aqr(i,j,k)%nh4)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%nh4 + aqr(i,j,k)%nh4)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%nh4 = aqr(i,j,k)%nh4 + qconv
        if(aqr(i,j,k)%nh4.le.0.0) aqr(i,j,k)%nh4 = 0.0
!
        cmean_up = (aqr(i,j,k)%xnv + aqr(i,j,k)%xnv)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%xnv + aqr(i,j,k)%xnv)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%xnv = aqr(i,j,k)%xnv + qconv
        if(aqr(i,j,k)%xnv.le.0.0) aqr(i,j,k)%xnv = 0.0
!
        cmean_up = (aqr(i,j,k)%ch2o + aqr(i,j,k)%ch2o)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%ch2o + aqr(i,j,k)%ch2o)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%ch2o = aqr(i,j,k)%ch2o + qconv
        if(aqr(i,j,k)%ch2o.le.0.0) aqr(i,j,k)%ch2o = 0.0
!
        cmean_up = (aqr(i,j,k)%ch3o2h + aqr(i,j,k)%ch3o2h)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%ch3o2h + aqr(i,j,k)%ch3o2h)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%ch3o2h = aqr(i,j,k)%ch3o2h + qconv
        if(aqr(i,j,k)%ch3o2h.le.0.0) aqr(i,j,k)%ch3o2h = 0.0
!
        cmean_up = (aqr(i,j,k)%siv + aqr(i,j,k)%siv)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%siv + aqr(i,j,k)%siv)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%siv = aqr(i,j,k)%siv + qconv
        if(aqr(i,j,k)%siv.le.0.0) aqr(i,j,k)%siv = 0.0
!
        cmean_up = (aqr(i,j,k)%svi + aqr(i,j,k)%svi)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%svi + aqr(i,j,k)%svi)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%svi = aqr(i,j,k)%svi + qconv
        if(aqr(i,j,k)%svi.le.0.0) aqr(i,j,k)%svi = 0.0
!
        cmean_up = (aqr(i,j,k)%nh4 + aqr(i,j,k)%nh4)*0.5 * qconv_up
        cmean_dn = (aqr(i,j,k)%nh4 + aqr(i,j,k)%nh4)*0.5 * qconv_dn
        qconv	 = cmean_up - cmean_dn
        aqr(i,j,k)%nh4 = aqr(i,j,k)%nh4 + qconv
        if(aqr(i,j,k)%nh4.le.0.0) aqr(i,j,k)%nh4 = 0.0
!
      enddo
    enddo
  enddo
!    
  endif
#endif
!  
return
end 
!
!
!
!  ==================================================
   subroutine solprecip   (pup, pdn, gas                   &
#ifdef AQCHEM_ENABLE
                           , aqr, aqc                      &
#endif
                           ,solidi)
!
!  ==================================================
!
! ================================================================
!
!  Purpose:
!	Subroutines calculating precipitation scavenging of solid
!       phase chemicals
!
! ================================================================
!
USE gridno
USE typedef_solid
USE typedef_aq
USE typedef_gas
USE shared_hydro
USE shared_data
!
  integer :: i, j, k
  real  :: qmean_up, qmean_dn, cmean_up, cmean_dn,    &
           qconv_up, qconv_dn
!	   
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz,1:nhydro)  :: pup, pdn
!
#include "chemdef.h"
!
#if ( defined SOLIDCHEM_ENABLE )
!
!----------------------------------------------------------!
!   
  if (lmicro > 1) then
!
  do k = 2, nzi
    do j = js, je
      do i = is, ie
!
	qmean_up = 0.5*(hydromtr2(i,j,k+1,rain)%q+hydromtr2(i,j,k,rain)%q)
  	if (qmean_up > qmin(ice)) then
  	  qconv_up = pup(i,j,k+1,ice)*dt / qmean_up
  	  if(qconv_up >= 1.0) qconv_up = 1.0
  	else
  	  qconv_up = 0.0
  	endif
! 	
	qmean_dn = 0.5*(hydromtr2(i,j,k-1,rain)%q+hydromtr2(i,j,k,rain)%q)
  	if (qmean_dn > qmin(ice)) then
  	  qconv_dn = pdn(i,j,k,ice)*dt/qmean_dn
  	  if(qconv_dn >= 1.0) qconv_dn = 1.0
  	else
  	  qconv_dn = 0.0
  	endif
! 	
  	cmean_up = (solidi(i,j,k)%o3 + solidi(i,j,k)%o3)*0.5 * qconv_up
  	cmean_dn = (solidi(i,j,k)%o3 + solidi(i,j,k)%o3)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%o3 = solidi(i,j,k)%o3 + qconv
  	if(solidi(i,j,k)%o3.le.0.0) solidi(i,j,k)%o3 = 0.0
! 	
  	cmean_up = (solidi(i,j,k)%h2o2 + solidi(i,j,k)%h2o2)*0.5 * qconv_up
  	cmean_dn = (solidi(i,j,k)%h2o2 + solidi(i,j,k)%h2o2)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%h2o2 = solidi(i,j,k)%h2o2 + qconv
  	if(solidi(i,j,k)%h2o2.le.0.0) solidi(i,j,k)%h2o2 = 0.0
! 	
  	cmean_up = (solidi(i,j,k)%nh4 + solidi(i,j,k)%nh4)*0.5 * qconv_up
  	cmean_dn = (solidi(i,j,k)%nh4 + solidi(i,j,k)%nh4)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%nh4 = solidi(i,j,k)%nh4 + qconv
  	if(solidi(i,j,k)%nh4.le.0.0) solidi(i,j,k)%nh4 = 0.0
! 	
  	cmean_up = (solidi(i,j,k)%xnv + solidi(i,j,k)%xnv)*0.5 * qconv_up
  	cmean_dn = (solidi(i,j,k)%xnv + solidi(i,j,k)%xnv)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%xnv = solidi(i,j,k)%xnv + qconv
  	if(solidi(i,j,k)%xnv.le.0.0) solidi(i,j,k)%xnv = 0.0
! 	
  	cmean_up = (solidi(i,j,k)%ch2o + solidi(i,j,k)%ch2o)*0.5 * qconv_up
  	cmean_dn = (solidi(i,j,k)%ch2o + solidi(i,j,k)%ch2o)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%ch2o = solidi(i,j,k)%ch2o + qconv
  	if(solidi(i,j,k)%ch2o.le.0.0) solidi(i,j,k)%ch2o = 0.0
! 	
  	cmean_up = (solidi(i,j,k)%ch3o2h + solidi(i,j,k)%ch3o2h)*0.5 * qconv_up
  	cmean_dn = (solidi(i,j,k)%ch3o2h + solidi(i,j,k)%ch3o2h)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%ch3o2h = solidi(i,j,k)%ch3o2h + qconv
  	if(solidi(i,j,k)%ch3o2h.le.0.0) solidi(i,j,k)%ch3o2h = 0.0
! 	
  	cmean_up = (solidi(i,j,k)%siv + solidi(i,j,k)%siv)*0.5* qconv_up
  	cmean_dn = (solidi(i,j,k)%siv + solidi(i,j,k)%siv)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%siv = solidi(i,j,k)%siv + qconv
  	if(solidi(i,j,k)%siv.le.0.0) solidi(i,j,k)%siv = 0.0
! 	
  	cmean_up = (solidi(i,j,k)%svi + solidi(i,j,k)%svi)*0.5 * qconv_up
  	cmean_dn = (solidi(i,j,k)%svi + solidi(i,j,k)%svi)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%svi = solidi(i,j,k)%svi + qconv
  	if(solidi(i,j,k)%svi.le.0.0) solidi(i,j,k)%svi = 0.0
! 	
  	cmean_up = (solidi(i,j,k)%nh4 + solidi(i,j,k)%nh4)*0.5 * qconv_up
  	cmean_dn = (solidi(i,j,k)%nh4 + solidi(i,j,k)%nh4)*0.5 * qconv_dn
  	qconv	 = cmean_up - cmean_dn
  	solidi(i,j,k)%nh4 = solidi(i,j,k)%nh4 + qconv
  	if(solidi(i,j,k)%nh4.le.0.0) solidi(i,j,k)%nh4 = 0.0
! 	
      enddo
    enddo
  enddo
!  
  endif
#endif
!	
return
end 
