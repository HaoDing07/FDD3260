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
!  AQ.F:
!	subroutine aqtrop (dtr, ifss, ktrop)
!	subroutine tropaq (tout, neq, y)
!	subroutine fex_aq (neq, t, y, ydot)
!	subroutine jex_aq (neq, t, y, ml, mu, pd, nrpd)
!	subroutine rateaq (T, patm, Rflux)                   
!	function cal_eta (xkn, xalfa)
!	subroutine calh( i, j, k, aqc, aqr )
!	function zbrent( nnn, i, j, k, aqx )
!	function func_c( x, i, j, k, aqc )
!	function func_r( x, i, j, k, aqr )
!
!  Purpose:
!	A package for calculating aqueous reactions in troposphere.			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
! ================================================================


!	================================================
	subroutine aqtrop ( ifss    &
#if (defined AQCHEM_ENABLE)
      		, gas, aqc, aqr     &
#endif
				)
!	================================================

	USE gridno
	USE typedef_gas
	USE typedef_aq
	USE typedef_solid
	USE shared_data
	USE shared_state
	USE shared_rad
	USE shared_hydro
	USE shared_thermo
	USE shared_pressure

	IMPLICIT NONE
	
! --------------------------------------------------
!  Brief Description of Dummy Variables:
!
!	dtr:	time step of integration in second
!	ifss:	1 for using steady-state assumption
!		0 for not using
!       ktrop:  layer label of tropopause
!
! --------------------------------------------------
!  Units of Variables:
!
!	Aqueous phase:
!		inside reaction:   mole/Lwater
!		outside reaction:  ppb(m)
!	Gaseous phase:
!		inside reaction:   mole/Lair
!		outside reaction:  ppb(m)
!
! --------------------------------------------------
!  Table of Variables:
!
!	ychem( 1) = o3(g)
!	ychem( 2) = zco2(g)
!	ychem( 3) = h2o2(g)
!	ychem( 4) = hno3(g) + xn2o5(g)
!	ychem( 5) = so2(g)
!	ychem( 6) = h2so4(g)
!
!	ychem( 7) = o3_c
!	ychem( 8) = civ_c
!	ychem( 9) = h2o2_c
!	ychem(10) = xnv_c
!	ychem(11) = siv_c
!	ychem(12) = svi_c
!
!	ychem(13) = ho_c
!	ychem(14) = ho2_c
!	ychem(15) = ch2(oh)2_c
!	ychem(16) = ch3o2h_c
!
!	ychem(17) = o3_r
!	ychem(18) = civ_r
!	ychem(19) = h2o2_r
!	ychem(20) = xnv_r
!	ychem(21) = siv_r
!	ychem(22) = svi_r
!
!	ychem(23) = ho_r
!	ychem(24) = ho2_r
!	ychem(25) = ch2(oh)2_r
!	ychem(26) = ch3o2h_r
!
!	ychem(27) = ho(g)
!	ychem(28) = ho2(g)
!	ychem(29) = ch2o(g)
!	ychem(30) = ch3o2h(g)
!
!       ychem(31) = nh3(g)
!       ychem(32) = nh4_c
!       ychem(33) = nh4_r
! --------------------------------------------------
	
	integer :: iii, ifss, i, j, k, ind
	real    :: ptam, Rflux,                                    &
     		   xxx, yyy, qconv, wgt_t, wgt_qc, wgt_qr, T298, patm

	real    :: pres  = 0.0
	real    :: Temp  = 0.0
	real    :: den   = 0.0
	real    :: r_den = 0.0

	real, dimension(naqvaria) :: ychem 	!= (/ (0.0, i=1,naqvaria) /)
	real, dimension(9)        :: awchem    	!= (/ (0.0, i=1,9) /)
	real, dimension(9)        :: revawchem 	!= (/ (0.0, i=1,9) /)

	logical :: l_qc, l_qr
	
	! --- tmp space saving parameters:
	! --- hplus_t 1 and 2 for [H+] of cloud and rain,
	! --- c_kn and r_kn for Knudson number of cloud and rain,
	! --- c_ar and r_ar for mean r of cloud and rain,
	! --- c_l and r_l for mole fraction of lwc of cloud and rain
	! --- RtimesT = 0.08205*T
	real, dimension(2) :: hplus_t
	real		   :: c_kn, r_kn, c_ar, r_ar, c_l, r_l
	real		   :: RtimesT

	common /hplustmp/hplus_t,c_kn,r_kn,c_ar,r_ar,c_l,r_l,RtimesT

#include "chemdef.h"

#ifdef AQCHEM_ENABLE

	! --- note: currently 19 diagnostic variables
	integer				:: mnv
	real, dimension (nz,19)		:: tdiag, t2diag
	real, dimension (nz,19)		:: gdiag
	real, dimension (nz)		:: wdiag
	
	save

! ===============================================================

	ychem    (:) = 0.0
	awchem   (:) = 0.0
	revawchem(:) = 0.0
	
	tdiag(:, :)     = 0.0

	! 
	! --- Note: all units involved in budget are in ppb(m)
	!
	do k=1,ktrop
	do j=jt_start,jt_end
	do i=it_start,it_end
	  tdiag (k, 1) = tdiag (k, 1) + gas(i,j,k)%o3*16.0/awO3
	  tdiag (k, 2) = tdiag (k, 2) + gas(i,j,k)%zco2*12.0/awCO2
	  tdiag (k, 3) = tdiag (k, 3) + gas(i,j,k)%h2o2/awH2O2
	  tdiag (k, 4) = tdiag (k, 4) + gas(i,j,k)%ho /awHO                &
     				     + gas(i,j,k)%ho2/awHO2
         tdiag (k, 5) = tdiag (k, 5) + gas(i,j,k)%ch2o*12.0/awCH2O
         tdiag (k, 6) = tdiag (k, 6) + gas(i,j,k)%ch3o2h*12.0/awCH3O2H
         tdiag (k, 7) = tdiag (k, 7) + gas(i,j,k)%so2*32.0/awSO2
         tdiag (k, 8) = tdiag (k, 8) + gas(i,j,k)%h2so4*32.0/awH2SO4
         tdiag (k, 9) = tdiag (k, 9) +(gas(i,j,k)%hno3 /awHNO3		   &
     				     + gas(i,j,k)%xn2o5/awN2O5)*14.0

         ! --- all aq diagnostic variables are in mole/Lair
         tdiag (k,10) = tdiag (k,10)					   &
     		      + (aqc(i,j,k)%o3   + aqr(i,j,k)%o3)*16.0/awO3
         tdiag (k,11) = tdiag (k,11)    				   &
     		      + (aqc(i,j,k)%civ  + aqr(i,j,k)%civ)*12.0/awCO2
         tdiag (k,12) = tdiag (k,12)  					   &
     		      + (aqc(i,j,k)%h2o2 + aqr(i,j,k)%h2o2)/awH2O2
         tdiag (k,13) = tdiag (k,13)  					   &
     		      + (aqc(i,j,k)%ch2o + aqr(i,j,k)%ch2o)*12.0/awCH2O
         tdiag (k,14) = tdiag (k,14)					   &
     		      + (aqc(i,j,k)%ch3o2h				   &
     		       + aqr(i,j,k)%ch3o2h)*12.0/awCH3O2H
         tdiag (k,15) = tdiag (k,15)    				   &
     		      + (aqc(i,j,k)%siv  + aqr(i,j,k)%siv)*32.0/awSO2
         tdiag (k,16) = tdiag (k,16)    				   &
     		      + (aqc(i,j,k)%svi  + aqr(i,j,k)%svi)*32.0/awH2SO4
         tdiag (k,17) = tdiag (k,17)					   &
     		      + (aqc(i,j,k)%xnv  + aqr(i,j,k)%xnv)*14.0/awHNO3
									    
!    	   --- nh3 and nh4
         tdiag (k,18) = tdiag (k,18) + gas(i,j,k)%nh3/awNH3
         tdiag (k,19) = tdiag (k,19)  					   &
     		      + (aqc(i,j,k)%nh4 + aqr(i,j,k)%nh4)/awNH3

       end do
       end do
       end do
            

       awchem(1) = awO3*1.e9
       awchem(2) = awCO2*1.e9
       awchem(3) = awH2O2*1.e9
       awchem(4) = awHNO3*1.e9
       awchem(5) = awSO2*1.e9
       awchem(6) = awH2SO4*1.e9
       awchem(7) = awCH2O*1.e9
       awchem(8) = awCH3O2H*1.e9
       awchem(9) = awNH3*1.e9

       do ind = 1,9
         revawchem(ind) = 1.0/awchem(ind)
       end do

       main_loop: do k=1,ktrop

       do j=jt_start,jt_end
       do i=it_start,it_end

         l_qc = hydromtr2(i,j,k,drop)%q.gt.qmin(drop).and.		   &
     		hydromtr2(i,j,k,drop)%n.gt.xnmin(drop)
         l_qr = hydromtr2(i,j,k,rain)%q.gt.qmin(rain).and.		   &
     		hydromtr2(i,j,k,rain)%n.gt.xnmin(rain)
	  
!	if( T1(i,j,k).gt.255.0.and.(l_qc.or.l_qr) )then
	if ( l_qc.or.l_qr ) then

	  do iii = 1,naqvaria
		ychem(iii) = 0.0
	  enddo
	  hplus_t(1) = 0.0
	  hplus_t(2) = 0.0

	  pres = pressure%p(i,j,k)+p0(k)
	  den  = pres/(2.87e2*thermo%T(i,j,k))
	  Temp = thermo%T(i,j,k)
	  r_den= 1./den

	  ! --- Pressure in Pa
	  patm = pres/101327.0

#ifdef RAD_ENABLE
	  ! --- Incoming solar flux W/m^2
	  Rflux= (rad%fluxsd(i,j,k+1) + rad%fluxsu(i,j,k)) * 0.5678
	  if(Rflux.le.50.0) Rflux = 0.0	!no photochemistry
 					!for nighttime
#else
	  Rflux= 1080.0*uuu0
#endif

	  !
	  ! --- Convert all aq species from ppb(m) to mole/Lair
	  !	ppb(m)*1.e-9*den/aw => mol/Lair
	  !
	  aqc(i,j,k)%o3     = aqc(i,j,k)%o3    *revawchem(1)*den
	  aqc(i,j,k)%civ    = aqc(i,j,k)%civ   *revawchem(2)*den
	  aqc(i,j,k)%h2o2   = aqc(i,j,k)%h2o2  *revawchem(3)*den
	  aqc(i,j,k)%xnv    = aqc(i,j,k)%xnv   *revawchem(4)*den
	  aqc(i,j,k)%siv    = aqc(i,j,k)%siv   *revawchem(5)*den
	  aqc(i,j,k)%svi    = aqc(i,j,k)%svi   *revawchem(6)*den
	  aqc(i,j,k)%ch2o   = aqc(i,j,k)%ch2o  *revawchem(7)*den
	  aqc(i,j,k)%ch3o2h = aqc(i,j,k)%ch3o2h*revawchem(8)*den
	  aqc(i,j,k)%nh4    = aqc(i,j,k)%nh4   *revawchem(9)*den
	
	  aqr(i,j,k)%o3     = aqr(i,j,k)%o3    *revawchem(1)*den
	  aqr(i,j,k)%civ    = aqr(i,j,k)%civ   *revawchem(2)*den
	  aqr(i,j,k)%h2o2   = aqr(i,j,k)%h2o2  *revawchem(3)*den
	  aqr(i,j,k)%xnv    = aqr(i,j,k)%xnv   *revawchem(4)*den
	  aqr(i,j,k)%siv    = aqr(i,j,k)%siv   *revawchem(5)*den
	  aqr(i,j,k)%svi    = aqr(i,j,k)%svi   *revawchem(6)*den
	  aqr(i,j,k)%ch2o   = aqr(i,j,k)%ch2o  *revawchem(7)*den
	  aqr(i,j,k)%ch3o2h = aqr(i,j,k)%ch3o2h*revawchem(8)*den
	  aqr(i,j,k)%nh4    = aqr(i,j,k)%nh4   *revawchem(9)*den

	  ! --- Calculate [H+]
	  T298 = (298.0 - thermo%T(i,j,k))/(298.0*thermo%T(i,j,k))
	  call calh( i, j, k, hydromtr2(i,j,k,1:nhydro),              &
     		     aqc(i,j,k), aqr(i,j,k), T298, pres, den )	  
	  hplus_t(1) = aqc(i,j,k)%hplus
	  hplus_t(2) = aqr(i,j,k)%hplus

	  call rateaq (Temp, patm, Rflux)
	
	  ! --- Convert mixing ratio from ppb(m) to mol/Lair:
	  ! --- Note:	N2O5 + H2O -> 2 HNO3 has been implemented
	  ! ---		here by combining it into HNO3(aq) before
	  ! ---		calculating aqueous reactions.
	  ychem( 1) = gas(i,j,k)%o3    *revawchem(1)*den
	  ychem( 2) = gas(i,j,k)%zco2  *revawchem(2)*den
	  ychem( 3) = gas(i,j,k)%h2o2  *revawchem(3)*den
	  ychem( 4) = gas(i,j,k)%hno3  *revawchem(4)*den
	  ychem( 5) = gas(i,j,k)%so2   *revawchem(5)*den
	  ychem( 6) = gas(i,j,k)%h2so4 *revawchem(6)*den
	
	  ychem(27) = gas(i,j,k)%ho    /awHO*1.e-9*den
	  ychem(28) = gas(i,j,k)%ho2   /awHO2*1.e-9*den
	  ychem(29) = gas(i,j,k)%ch2o  *revawchem(7)*den
	  ychem(30) = gas(i,j,k)%ch3o2h*revawchem(8)*den

	  ychem(31) = gas(i,j,k)%nh3   *revawchem(9)*den

	  wgt_t  = 1./(hydromtr2(i,j,k,drop)%q + hydromtr2(i,j,k,rain)%q)
	  wgt_qc = hydromtr2(i,j,k,drop)%q*wgt_t
	  wgt_qr = hydromtr2(i,j,k,rain)%q*wgt_t

	  wgt_t  = 2.0*gas(i,j,k)%xn2o5*1.e-9/awN2O5*den	  
	  aqc(i,j,k)%xnv = aqc(i,j,k)%xnv + wgt_t*wgt_qc	! Mair
      	  aqr(i,j,k)%xnv = aqr(i,j,k)%xnv + wgt_t*wgt_qr	! Mair
          gas(i,j,k)%xn2o5 = 0.0
	  
	  RtimesT = 0.08205*Temp	!RT
	  	
	  if (l_qc) then
     		! --- Convert units of aqueous chemicals
     	  	! --- from mole/Lair to mole/Lwater
     		qconv = 1.e3/(den*hydromtr2(i,j,k,drop)%q)
	  	ychem( 7) = aqc(i,j,k)%o3    *qconv
	  	ychem( 8) = aqc(i,j,k)%civ   *qconv
	  	ychem( 9) = aqc(i,j,k)%h2o2  *qconv
	  	ychem(10) = aqc(i,j,k)%xnv   *qconv
	  	ychem(11) = aqc(i,j,k)%siv   *qconv
	  	ychem(12) = aqc(i,j,k)%svi   *qconv
		ychem(15) = aqc(i,j,k)%ch2o  *qconv
		ychem(16) = aqc(i,j,k)%ch3o2h*qconv
	  	ychem(32) = aqc(i,j,k)%nh4  *qconv
		
	  	! === initializing short-lived species with
		! === equilibrium aqueous concentrations
		! === [C(aq)] = Hc*Pc/(1 + Hc*L*R*T)
		if (ifss.eq.1) then
		  xxx = hydromtr2(i,j,k,drop)%q*den*1.e-3*RtimesT	! L
		  ychem(13) = rke_c(6)*gas(i,j,k)%ho                      &
     			   *28.97296245e-9/awHO 			  &
     			   /(1. + rke_c(6)*xxx)
        	 ychem(14) = rke_c(8)*gas(i,j,k)%ho2			  &
     			   *28.97296245e-9/awHO2			  &
     			   /(1. + rke_c(8)*xxx) 		 
     
        	 ! --- speedup calculation by giving an initial
        	 ! --- pssa solution to [OH] and [HO2]
        	 xxx = rka_c(1)*ychem(9) 				  &
     		     + rka_c(10)*ychem(7)*ychem(14)
        	 yyy =(rka_c(3) + rka_c(5))*ychem(14)			  &
     		     + rka_c(4) *ychem(9)  + rka_c(6) *ychem(7) 	  &
     		     + rka_c(7) *ychem(13) + rka_c(11)*ychem(16)	  &
     		     + rka_c(12)*ychem(15) + rka_c(15)*ychem(8) 	  &
     		     +(rka_c(21) + rka_c(22))*ychem(11)
        	 if (yyy.gt.0.0)					  &
     		       ychem(13) = xxx/yyy
     
        	 xxx = rka_c(16)*ychem(14) + rka_c(17)*ychem(9)
        	 if (xxx.gt.0.0) then
        	       xxx  = rka_c(17)*ychem(9)/xxx
        	 else
        	       xxx  = 0.0
        	 endif
        	 xxx = (rka_c(14)*ychem(14) 				  &
     		      + rka_c(15)*ychem(13))*ychem(8)*xxx		  &
     		     + (rka_c(4) *ychem(9) 				  &
     		      + rka_c(6) *ychem(7) 				  &
     		      + 2.0*rka_c(12)*ychem(15))*ychem(13)
        	 yyy = (rka_c(3) + rka_c(5))*ychem(13)			  &
     		     +  rka_c(8)*ychem(14) + rka_c(10)*ychem(7)
        	 if(yyy.gt.0.0) 					  &
     		       ychem(14) = xxx/yyy
     
               endif

               !mean radius	       
               c_ar = 0.047538						  &
     		    *(hydromtr2(i,j,k,drop)%q/hydromtr2(i,j,k,drop)%n)**(1.0/3.0)
     		!knudson number
     		c_kn = 6.6e-8*(Temp/273.15*patm)/c_ar
		!mole fraction of water vapor
		c_l  = hydromtr2(i,j,k,drop)%q*den*1.e-3
	  else
	  	aqc(i,j,k)%o3	  = 0.0
	  	aqc(i,j,k)%civ	  = 0.0
	  	aqc(i,j,k)%h2o2	  = 0.0
	  	aqc(i,j,k)%xnv	  = 0.0
	  	aqc(i,j,k)%siv	  = 0.0
	  	aqc(i,j,k)%svi	  = 0.0
		aqc(i,j,k)%ch2o	  = 0.0
		aqc(i,j,k)%ch3o2h = 0.0
		aqc(i,j,k)%nh4    = 0.0
		aqc(i,j,k)%hplus  = 0.0

	  	ychem( 7) = 0.0
	  	ychem( 8) = 0.0
	  	ychem( 9) = 0.0
	  	ychem(10) = 0.0
	  	ychem(11) = 0.0
	  	ychem(12) = 0.0
		ychem(13) = 0.0
	  	ychem(14) = 0.0
	  	ychem(15) = 0.0
	  	ychem(16) = 0.0
	  	ychem(32) = 0.0
		
		c_ar = 0.0
		c_kn = 0.0
		c_l  = 0.0
	  endif
	
	  if (l_qr) then
     		! --- Convert units of aqueous chemicals
     	  	! --- from mole/Lair to moel/Lwater
     		qconv = 1.e3/(den*hydromtr2(i,j,k,rain)%q)
	  	ychem(17) = aqr(i,j,k)%o3    *qconv
	  	ychem(18) = aqr(i,j,k)%civ   *qconv
	  	ychem(19) = aqr(i,j,k)%h2o2  *qconv
	  	ychem(20) = aqr(i,j,k)%xnv   *qconv
	  	ychem(21) = aqr(i,j,k)%siv   *qconv
	  	ychem(22) = aqr(i,j,k)%svi   *qconv
		ychem(25) = aqr(i,j,k)%ch2o  *qconv
		ychem(26) = aqr(i,j,k)%ch3o2h*qconv
		ychem(33) = aqr(i,j,k)%nh4   *qconv

	  	! === initializing short-lived species with
		! === equilibrium aqueous concentrations
		! === [C(aq)] = Hc*Pc/(1 + Hc*L*R*T)
		if (ifss.eq.1) then
		  xxx = hydromtr2(i,j,k,rain)%q*den*1.e-3*RtimesT	! L
	  	  ychem(23) = rke_r(6)*gas(i,j,k)%ho                   &		 
     			   *28.97296245e-9/awHO 		       &
     			   /(1. + rke_r(6)*xxx)
        	 ychem(24) = rke_r(8)*gas(i,j,k)%ho2		       &
     			   *28.97296245e-9/awHO2		       &
     			   /(1. + rke_r(8)*xxx)

        	 ! --- speedup calculation by giving an initial
        	 ! --- pssa solution to [OH] and [HO2]
        	 xxx = rka_r(1)*ychem(19) 			       &
     		     + rka_r(10)*ychem(17)*ychem(24)
        	 yyy =(rka_r(3) + rka_r(5))*ychem(24)		       &
     		     + rka_r(4) *ychem(19) + rka_r(16)*ychem(17)       &
     		     + rka_r(7) *ychem(23) + rka_r(21)*ychem(26)       &
     		     + rka_r(12)*ychem(25) + rka_r(25)*ychem(18)       &
     		     +(rka_r(21) + rka_r(22))*ychem(21)
        	 if (yyy.gt.0.0)				       &
     		       ychem(23) = xxx/yyy
     
        	 xxx = rka_r(16)*ychem(24) + rka_r(17)*ychem(19)
        	 if (xxx.gt.0.0) then
        	       xxx  = rka_r(17)*ychem(19)/xxx
        	 else
        	       xxx  = 0.0
        	 endif
        	 xxx = (rka_r(14)*ychem(24) 			       &
     		      + rka_r(15)*ychem(23))*ychem(18)*xxx	       &
     		     + (rka_r(4) *ychem(19) 			       &
     		      + rka_r(6) *ychem(17) 			       &
     		      + 2.0*rka_r(12)*ychem(25))*ychem(23)
        	 yyy = (rka_r(3) + rka_r(5))*ychem(23)		       &
     		     +  rka_r(8)*ychem(24) + rka_r(10)*ychem(17)
        	 if(yyy.gt.0.0) 				       &
     		       ychem(24) = xxx/yyy			        
     	 
               endif						        
     		       
               !mean radius	       
               r_ar = 0.034139					       &
     		    *(hydromtr2(i,j,k,rain)%q/                          &
		    hydromtr2(i,j,k,rain)%n)**(1.0/3.0)
     		!knudson number
     		r_kn = 6.6e-8*(Temp/273.15*patm)/r_ar
		!mole fraction of water vapor
		r_l  = hydromtr2(i,j,k,rain)%q*den*1.e-3
	  else
	  	aqr(i,j,k)%o3	  = 0.0
	  	aqr(i,j,k)%civ	  = 0.0
	  	aqr(i,j,k)%h2o2	  = 0.0
	  	aqr(i,j,k)%xnv	  = 0.0
	  	aqr(i,j,k)%siv	  = 0.0
	  	aqr(i,j,k)%svi	  = 0.0
		aqr(i,j,k)%ch2o	  = 0.0
		aqr(i,j,k)%ch3o2h = 0.0
		aqr(i,j,k)%hplus  = 0.0
		aqr(i,j,k)%nh4  = 0.0

	  	ychem(17) = 0.0
	  	ychem(18) = 0.0
	  	ychem(19) = 0.0
	  	ychem(20) = 0.0
	  	ychem(21) = 0.0
	  	ychem(22) = 0.0
	  	ychem(23) = 0.0
	  	ychem(24) = 0.0
	  	ychem(25) = 0.0
	  	ychem(26) = 0.0
	  	ychem(33) = 0.0
				
		r_ar = 0.0
		r_kn = 0.0
		r_l  = 0.0
	  end if

	  do iii = 1,naqvaria
		if(ychem(iii).lt.1.e-24) ychem(iii) = 0.0
	  enddo
	  
! ===
	  call tropaq(naqvaria, ychem)
! ===

	  do iii = 1,naqvaria
	  	if(ychem(iii).lt.1.e-24) ychem(iii) = 0.0
	  end do

	  ! --- Give value back to the common block:
	  gas(i,j,k)%o3	    = ychem( 1)*awchem(1)*r_den
	  gas(i,j,k)%zco2   = ychem( 2)*awchem(2)*r_den
	  gas(i,j,k)%h2o2   = ychem( 3)*awchem(3)*r_den
	  gas(i,j,k)%hno3   = ychem( 4)*awchem(4)*r_den
	  gas(i,j,k)%so2    = ychem( 5)*awchem(5)*r_den
	  gas(i,j,k)%h2so4  = ychem( 6)*awchem(6)*r_den
          gas(i,j,k)%nh3    = ychem(31)*awchem(9)*r_den

	  ! 
	  ! --- convert from mol/Lwater to ppb(m)
	  !	mol/Lwater*q*aw*1.e9/1.e3(denw) => ppb(m)
	  !
	  qconv = 1.e-3*hydromtr2(i,j,k,drop)%q
	  aqc(i,j,k)%o3	    = ychem( 7)*awchem(1)*qconv
	  aqc(i,j,k)%civ    = ychem( 8)*awchem(2)*qconv
	  aqc(i,j,k)%h2o2   = ychem( 9)*awchem(3)*qconv
	  aqc(i,j,k)%xnv    = ychem(10)*awchem(4)*qconv
	  aqc(i,j,k)%siv    = ychem(11)*awchem(5)*qconv
	  aqc(i,j,k)%svi    = ychem(12)*awchem(6)*qconv
	  aqc(i,j,k)%ch2o   = ychem(15)*awchem(7)*qconv
	  aqc(i,j,k)%ch3o2h = ychem(16)*awchem(8)*qconv
	  aqc(i,j,k)%nh4    = ychem(32)*awchem(9)*qconv

	  qconv = 1.e-3*hydromtr2(i,j,k,rain)%q
	  aqr(i,j,k)%o3	    = ychem(17)*awchem(1)*qconv
	  aqr(i,j,k)%civ    = ychem(18)*awchem(2)*qconv
	  aqr(i,j,k)%h2o2   = ychem(19)*awchem(3)*qconv
	  aqr(i,j,k)%xnv    = ychem(20)*awchem(4)*qconv
	  aqr(i,j,k)%siv    = ychem(21)*awchem(5)*qconv
	  aqr(i,j,k)%svi    = ychem(22)*awchem(6)*qconv
	  aqr(i,j,k)%ch2o   = ychem(25)*awchem(7)*qconv
	  aqr(i,j,k)%ch3o2h = ychem(26)*awchem(8)*qconv
	  aqr(i,j,k)%nh4    = ychem(33)*awchem(9)*qconv
	  
	  gas(i,j,k)%ho     = ychem(27)*awHO*1.e9*r_den
	  gas(i,j,k)%ho2    = ychem(28)*awHO2*1.e9*r_den
	  gas(i,j,k)%ch2o   = ychem(29)*awCH2O*1.e9*r_den
	  gas(i,j,k)%ch3o2h = ychem(30)*awCH3O2H*1.e9*r_den

	else
	  hplus_t(1) = 0.0
	  hplus_t(2) = 0.0
	  aqc(i,j,k)%hplus = 0.0
	  aqr(i,j,k)%hplus = 0.0
	endif
	
	end do	!i loop
	end do	!j loop

	end do main_loop
!	----------------

	t2diag(:,:) = 0.0 
	! 
	! --- Note: all units involved in budget are in ppb(m)
	!
	do k=1,ktrop
	do j=jt_start,jt_end
	do i=it_start,it_end
	  t2diag (k, 1) = t2diag (k, 1) + gas(i,j,k)%o3*16.0/awO3
	  t2diag (k, 2) = t2diag (k, 2) + gas(i,j,k)%zco2*12.0/awCO2
	  t2diag (k, 3) = t2diag (k, 3) + gas(i,j,k)%h2o2/awH2O2
	  t2diag (k, 4) = t2diag (k, 4) + gas(i,j,k)%ho /awHO             &
     				       + gas(i,j,k)%ho2/awHO2
         t2diag (k, 5) = t2diag (k, 5) + gas(i,j,k)%ch2o*12.0/awCH2O
         t2diag (k, 6) = t2diag (k, 6) + gas(i,j,k)%ch3o2h*12.0/awCH3O2H
         t2diag (k, 7) = t2diag (k, 7) + gas(i,j,k)%so2*32.0/awSO2
         t2diag (k, 8) = t2diag (k, 8) + gas(i,j,k)%h2so4*32.0/awH2SO4
         t2diag (k, 9) = t2diag (k, 9) +(gas(i,j,k)%hno3 /awHNO3	  &
     				       + gas(i,j,k)%xn2o5/awN2O5)*14.0

         ! --- all aq diagnostic variables are in mole/Lair
         t2diag (k,10) = t2diag (k,10)    				  &
     		      + (aqc(i,j,k)%o3   + aqr(i,j,k)%o3)*16.0/awO3
         t2diag (k,11) = t2diag (k,11)   				  &
     		      + (aqc(i,j,k)%civ  + aqr(i,j,k)%civ)*12.0/awCO2
         t2diag (k,12) = t2diag (k,12)  				  &
     		      + (aqc(i,j,k)%h2o2 + aqr(i,j,k)%h2o2)/awH2O2
         t2diag (k,13) = t2diag (k,13)  				  &
     		      + (aqc(i,j,k)%ch2o + aqr(i,j,k)%ch2o)*12.0/awCH2O
         t2diag (k,14) = t2diag (k,14)					  &
     		      + (aqc(i,j,k)%ch3o2h				  &
     		       + aqr(i,j,k)%ch3o2h)*12.0/awCH3O2H
         t2diag (k,15) = t2diag (k,15)   				  &
     		      + (aqc(i,j,k)%siv  + aqr(i,j,k)%siv)*32.0/awSO2
         t2diag (k,16) = t2diag (k,16)   				  &
     		      + (aqc(i,j,k)%svi  + aqr(i,j,k)%svi)*32.0/awH2SO4
         t2diag (k,17) = t2diag (k,17)    				  &
     		      + (aqc(i,j,k)%xnv  + aqr(i,j,k)%xnv)*14.0/awHNO3
									   
!    	   --- add nh3 and nh4
         t2diag (k,18) = t2diag (k,18) + gas(i,j,k)%nh3/awNH3
         t2diag (k,19) = t2diag (k,19)  				  &
     		      + (aqc(i,j,k)%nh4 + aqr(i,j,k)%nh4)/awNH3
	end do
	end do
	end do

	tdiag(:,:) = t2diag(:,:) - tdiag(:,:)
	
#endif

	return
	 end

 
!	================================
	subroutine tropaq (neq, y)
!	================================

! --------------------------------------------------------
! --- TROPREACT1.F:  Subroutine for calculating chemical
! --- 			reactions in troposphere
! ---
! --- Last Revised on:    July 19, 1998
! --------------------------------------------------------

	external fex_aq, jex_aq

	real :: y(neq), rwork(1408),iwork(53), atol(neq)

	integer :: neq, itol, itask, istate, iopt, lrw, liw, mf
	real    :: t, rtol

#ifdef AQCHEM_ENABLE
	t = 0.
	itol = 2	!1
	rtol = 1.e-5
!	atol = 1.e-6	! mols/Lwater

	atol( 1) = 1.e-16	!o3(g)
	atol( 2) = 1.e-16	!zco2(g)
	atol( 3) = 1.e-22	!h2o2(g)
	atol( 4) = 1.e-22	!hno3(g) + xn2o5(g)
	atol( 5) = 1.e-20	!so2(g)
	atol( 6) = 1.e-20	!h2so4(g)
	atol( 7) = 1.e-12	!o3_c
	atol( 8) = 1.e-12	!civ_c
	atol( 9) = 1.e-12	!h2o2_c
	atol(10) = 1.e-12	!xnv_c
	atol(11) = 1.e-12	!siv_c
	atol(12) = 1.e-12	!svi_c
	atol(13) = 1.e-20	!ho_c
	atol(14) = 1.e-20	!ho2_c
	atol(15) = 1.e-18	!ch2(oh)2_c
	atol(16) = 1.e-18	!ch3o2h_c
	atol(17) = 1.e-12	!o3_r
	atol(18) = 1.e-12	!civ_r
	atol(19) = 1.e-12	!h2o2_r
	atol(20) = 1.e-12	!xnv_r
	atol(21) = 1.e-12	!siv_r
	atol(22) = 1.e-12	!svi_r
	atol(23) = 1.e-20	!ho_r
	atol(24) = 1.e-20	!ho2_r
	atol(25) = 1.e-18	!ch2(oh)2_r
	atol(26) = 1.e-18	!ch3o2h_r
	atol(27) = 1.e-24	!ho(g)
	atol(28) = 1.e-24	!ho2(g)
	atol(29) = 1.e-22	!ch2o(g)
	atol(30) = 1.e-22	!ch3o2h(g)
	atol(31) = 1.e-22	!nh3(g)
	atol(32) = 1.e-12	!nh3_c
	atol(33) = 1.e-12	!nh3_r
		
	itask  = 1
	istate = 1
	iopt   = 0
	lrw    = 1408	!= 22 + 9*neq + neq**2 for mf = 21
	liw    = 53	!= 20 + neq for mf = 21
	mf     = 21	!stiff (bdf) and user supply full jacobian

	call lsodenew(fex_aq,neq,y,t,                   &
     		     time,itol,rtol,atol,itask,istate,	&
     		     iopt,rwork,lrw,iwork,liw,jex_aq,mf)
#endif

	return
	 end

!	===================================
	subroutine fex_aq (neq, t, y, ydot)
!	===================================
!
!  Subroutine for defining of ydot
!
     
	USE gridno
	USE shared_data

	integer :: neq
	real    :: t, y(neq), ydot(neq)     

#ifdef AQCHEM_ENABLE
	real    :: so5m, so5m2, xxx, yyy, fkmt,                       &
     		  rr1, rr2, rr3, rr4, rr5, rr6, rr7, rr8,	      &
     		  rr10,rr11,rr12,rr13,rr14,rr15,rr19,rr20,	      &
     		  rr21,rr22,rr24,rr25,rr26,rr27,rr30,		      &
     		  e11y1y7,e3y2y8,e9y3y9,e13y4y10,e15y5y11,e18y6y12,   &
     		  e11y1y17,e3y2y18,e9y3y19,e13y4y20,e15y5y21,	      &
     		  e18y6y22,e6y27y23,e8y28y24,e12y29y25,e21y30y26,     &
     		  e22y31y32,e22y31y33
     
       ! --- tmp space saving parameters:
       ! --- hplus_t 1 and 2 for [H+] of cloud and rain,
       ! --- c_kn and r_kn for Knudson number of cloud and rain,
       ! --- c_ar and r_ar for mean r of cloud and rain,
       ! --- c_l and r_l for mole fraction of lwc of cloud and rain
       ! --- RtimesT = 0.08205*T
       real, dimension(2) :: hplus_t   
       common /hplustmp/hplus_t,c_kn,r_kn,c_ar,r_ar,c_l,r_l,RtimesT
		
! === Equations Applied with PSSA Solution:
! -------------------------------------------------------------------
! --- PSSA Solutions:
! --- [HCOOH]ss = (RA12/RA13)[CH2(OH)2]
! --- [CO3-]ss	= {RA14[HO2]+RA15[OH]}[C(IV)]/{RA16[HO2]+RA17[H2O2]}
! --- [SO5-]ss	= SQRT{(RA21+RA22)[S(IV)][OH]/(2RA31)}
! --- [HSO5-]ss	= (RA24 + RA25)[ SO5-]/(RA32 [H+])
! --- [SO4-]ss	= {(RA26+RA27)[S(IV)][SO5-]
! ---            +2RA30[SO5-]**2}/{(RA28+RA29)[S(IV)]}
! --- [SO3-]ss	= {(RA24+RA25+RA26+RA27)[S(IV)][SO5-]
! ---            +2(RA30+RA31)[SO5-]**2}/(RA23[O2])
! -------------------------------------------------------------------
! --------------------------------------------------
!  Table of Variables:
!	ychem( 1) = o3(g)
!	ychem( 2) = zco2(g)
!	ychem( 3) = h2o2(g)
!	ychem( 4) = hno3(g) + xn2o5(g)
!	ychem( 5) = so2(g)
!	ychem( 6) = h2so4(g)
!
!	ychem( 7) = o3_c
!	ychem( 8) = civ_c
!	ychem( 9) = h2o2_c
!	ychem(10) = xnv_c
!	ychem(11) = siv_c
!	ychem(12) = svi_c
!
!	ychem(13) = ho_c
!	ychem(14) = ho2_c
!	ychem(15) = ch2(oh)2_c
!	ychem(16) = ch3o2h_c
!
!	ychem(17) = o3_r
!	ychem(18) = civ_r
!	ychem(19) = h2o2_r
!	ychem(20) = xnv_r
!	ychem(21) = siv_r
!	ychem(22) = svi_r

!	ychem(23) = ho_r
!	ychem(24) = ho2_r
!	ychem(25) = ch2(oh)2_r
!	ychem(26) = ch3o2h_r
!
!	ychem(27) = ho(g)
!	ychem(28) = ho2(g)
!	ychem(29) = ch2o(g)
!	ychem(30) = ch3o2h(g)
!
!       ychem(31) = nh3(g)
!       ychem(32) = nh4_c
!       ychem(33) = nh4_r
! --------------------------------------------------

	do iii=1,neq
	  ydot(iii) = 0.0
!	  if(y(iii).le.0.0) y(iii) = 0.0
	end do
	
! --- For cloud water

	if(hplus_t(1).gt.0.0.and.c_l.gt.0.0)then
	  fkmt= 3.0e-5/c_ar**2
	  	  	  
	  ! === [SO5-]ss**2 and [SO5-]ss
	  if(rka_c(31).gt.0.0                                         &
     	   .and.y(11).gt.1.e-24.and.y(13).gt.1.e-24)then
           so5m2= (rka_c(21) + rka_c(22))			      &
     		*y(11)*y(13)/(2*rka_c(31))
           so5m = sqrt(so5m2)
         else
           so5m2= 0.0
           so5m = 0.0
         endif

         ! === [CO3-]ss left RK17[H2O2]/(RK16[HO2] + RK17[H2O2])
         xxx = rka_c(16)*y(14) + rka_c(17)*y(9)
         if(xxx.gt.1.e-26)then
           xxx  = rka_c(17)*y(9)/xxx
         else
           xxx  = 0.0
         endif

         ! === reduce number of multiplications        
         rr1  = rka_c(1) *y(9)         !rk1[H2O2]
         rr2  = rka_c(2) *y(7)         !rk2[O3]
         rr3  = rka_c(3) *y(13)*y(14)  !rk3[OH][HO2]
         rr4  = rka_c(4) *y(13)*y(9)   !rk4[OH][H2O2]
         rr5  = rka_c(5) *y(13)*y(14)  !rk5[OH][HO2]
         rr6  = rka_c(6) *y(13)*y(7)   !rk6[OH][O3]
         rr7  = rka_c(7) *y(13)**2     !rk7[OH]^2
         rr8  = rka_c(8) *y(14)**2     !rk8[HO2]^2
         rr10 = rka_c(10)*y(14)*y(7)   !rk10[HO2][O3]
         rr11 = rka_c(11)*y(13)*y(16)  !rk11[OH][CH3O2H]
         rr12 = rka_c(12)*y(13)*y(15)  !rk12[OH][CH2(OH)2]
         rr13 = rka_c(12)*y(13)*y(15)  !rk12[OH][CH2(OH)2] (HCOOHss)
         rr14 = rka_c(14)*y(14)*y(8)   !rk14[HO2][C(IV)]
         rr15 = rka_c(15)*y(13)*y(8)   !rk15[OH][C(IV)]
         rr19 = rka_c(19)*y(11)*y(7)   !rk19[S(IV)][O3]
         rr20 = rka_c(20)*y(11)*y(9)   !rk20[S(IV)][H2O2]
         rr21 = rka_c(21)*y(13)*y(11)  !rk21[S(IV)][OH]
         rr22 = rka_c(22)*y(13)*y(11)  !rk22[S(IV)][OH]
         rr24 = rka_c(24)*y(11)*so5m   !rk24[S(IV)][SO5-]
         rr25 = rka_c(25)*y(11)*so5m   !rk25[S(IV)][SO5-]
         rr26 = rka_c(26)*y(11)*so5m   !rk26[S(IV)][SO5-]
         rr27 = rka_c(27)*y(11)*so5m   !rk27[S(IV)][SO5-]
         rr30 = rka_c(30)*so5m2        !rk30[SO5-]^2

         e11y1y7  = fkmt*cal_eta (c_kn, alphaO3)		      &
     		  *( y(1) - y( 7)/(rke_c(11)*RtimesT) )
         e3y2y8   = fkmt*cal_eta (c_kn, alphaCO2)		      &
     		  *( y(2) - y( 8)/(rke_c( 3)*RtimesT) )
         e9y3y9   = fkmt*cal_eta (c_kn, alphaH2O2)		      &
     		  *( y(3) - y( 9)/(rke_c( 9)*RtimesT) )
         e13y4y10 = fkmt*cal_eta (c_kn, alphaHNO3)		      &
     		  *( y(4) - y(10)/(rke_c(13)*RtimesT) )
         e15y5y11 = fkmt*cal_eta (c_kn, alphaSO2)		      &
     		  *( y(5) - y(11)/(rke_c(15)*RtimesT) )
         e18y6y12 = fkmt*cal_eta (c_kn, alphaH2SO4)		      &
     		  *( y(6) - y(12)/(rke_c(18)*RtimesT) )
         e6y27y13 = fkmt*cal_eta (c_kn, alphaEST)		      &
     		  *( y(27)- y(13)/(rke_c( 6)*RtimesT) )
         e8y28y14 = fkmt*cal_eta (c_kn, alphaEST)		      &
     		  *( y(28)- y(14)/(rke_c( 8)*RtimesT) )
         e12y29y15= fkmt*cal_eta (c_kn, alphaEST)		      &
     		  *( y(29)- y(15)/(rke_c(12)*RtimesT) )
         e21y30y16= fkmt*cal_eta (c_kn, alphaEST)		      &
     		  *( y(30)- y(16)/(rke_c(21)*RtimesT) )
         e22y31y32= fkmt*cal_eta (c_kn, alphaNH3)		      &
     		  *( y(31) - y(32)/(rke_c(22)*RtimesT) )	       
     	 
         ydot( 1) = - c_l*e11y1y7 
         ydot( 2) = - c_l*e3y2y8
         ydot( 3) = - c_l*e9y3y9
         ydot( 4) = - c_l*e13y4y10
         ydot( 5) = - c_l*e15y5y11
         ydot( 6) = - c_l*e18y6y12
               
         ydot( 7) = e11y1y7 - rr2 - rr6 - rr10 - rr19
         ydot( 8) = e3y2y8 + rr13
         ydot( 9) = e9y3y9 + rr2 + rr7 + rr8			      &
     		  - rr1 - rr4 - rr20				      &
     		  - (rr14 + rr15)*xxx
         ydot(10) = e13y4y10
         ydot(11) = e15y5y11 - rr19 - rr20 - rr21 - rr22	      &
     		  - 2.0*(rr24 + rr25 + rr26 + rr27 + rr30)
         ydot(12) = e18y6y12 + rr19 + rr20 			      &
     		  + rr24 + rr25 + 2.0*(rr26 + rr27 + rr30)
         ydot(13) = e6y27y13 + 2.0*rr1 + rr10			      &
     		  - rr3 - rr4 - rr5 - rr6 - rr7 		      &
     		  - rr11 - rr12 - rr13 - rr15 - rr21 - rr22
         ydot(14) = e8y28y14 + (rr14 + rr15)*xxx  		      &
     		  + rr4 + rr6 + rr12 + rr13			      &
     		  - rr3 - rr5 - rr8 - rr10
	  ydot(15) = e12y29y15 + rr11 - rr12
	  ydot(16) = e21y30y16 - rr11
	  
	  ydot(27) = - c_l*e6y27y13
	  ydot(28) = - c_l*e8y28y14
	  ydot(29) = - c_l*e12y29y15
	  ydot(30) = - c_l*e21y30y16
          ydot(31) = - c_l*e22y31y32
          ydot(32) = e22y31y32
	endif

! --- For rain water

	if(hplus_t(2).gt.0.0.and.r_l.gt.0.0)then
	  fkmt= 3.0e-5/r_ar**2

	  ! === [SO5-]ss**2 and [SO5-]ss
	  if(rka_r(31).gt.0.0                                        &
     	   .and.y(21).gt.1.e-24.and.y(23).gt.1.e-24)then
           so5m2= (rka_r(21) + rka_r(22)) 			     &
     		*y(21)*y(23)/(2*rka_r(31))
           so5m = sqrt(so5m2)
         else
           so5m2= 0.0
           so5m = 0.0
         endif

         ! === [CO3-]ss left RK17[H2O2]/(RK16[HO2] + RK17[H2O2])
         xxx = rka_r(16)*y(24) + rka_r(17)*y(19)
         if(xxx.gt.1.e-26)then
           xxx  = rka_r(17)*y(19)/xxx
         else
           xxx  = 0.0
         endif

         ! === reduce number of multiplications        
         rr1  = rka_r(1) *y(19)        !rk1[H2O2]
         rr2  = rka_r(2) *y(17)        !rk2[O3]
         rr3  = rka_r(3) *y(23)*y(24)  !rk3[OH][HO2]
         rr4  = rka_r(4) *y(23)*y(19)  !rk4[OH][H2O2]
         rr5  = rka_r(5) *y(23)*y(24)  !rk5[OH][HO2]
         rr6  = rka_r(6) *y(23)*y(17)  !rk6[OH][O3]
         rr7  = rka_r(7) *y(23)**2     !rk7[OH]^2
         rr8  = rka_r(8) *y(24)**2     !rk8[HO2]^2
         rr10 = rka_r(10)*y(24)*y(17)  !rk10[HO2][O3]
         rr11 = rka_r(11)*y(23)*y(26)  !rk11[OH][CH3O2H]
         rr12 = rka_r(12)*y(23)*y(25)  !rk12[OH][CH2(OH)2]
         rr13 = rka_r(12)*y(25)*y(23)  !rk13[OH][CH2(OH)2] (HCOOHss)
         rr14 = rka_r(14)*y(24)*y(18)  !rk14[HO2][C(IV)]
         rr15 = rka_r(15)*y(23)*y(18)  !rk15[OH][C(IV)]
         rr19 = rka_r(19)*y(21)*y(17)  !rk19[S(IV)][O3]
         rr20 = rka_r(20)*y(21)*y(19)  !rk20[S(IV)][H2O2]
         rr21 = rka_r(21)*y(21)*y(23)  !rk21[S(IV)][OH]
         rr22 = rka_r(22)*y(21)*y(23)  !rk22[S(IV)][OH]
         rr24 = rka_r(24)*y(21)*so5m   !rk24[S(IV)][SO5-]
         rr25 = rka_r(25)*y(21)*so5m   !rk25[S(IV)][SO5-]
         rr26 = rka_r(26)*y(21)*so5m   !rk26[S(IV)][SO5-]
         rr27 = rka_r(27)*y(21)*so5m   !rk27[S(IV)][SO5-]
         rr30 = rka_r(30)*so5m2        !rk30[SO5-]^2

         e11y1y17 = fkmt*cal_eta (r_kn, alphaO3)		     &
     		  *( y(1) - y(17)/(rke_r(11)*RtimesT) )
         e3y2y18  = fkmt*cal_eta (r_kn, alphaCO2)		     &
     		  *( y(2) - y(18)/(rke_r( 3)*RtimesT) )
         e9y3y19  = fkmt*cal_eta (r_kn, alphaH2O2)		     &
     		  *( y(3) - y(19)/(rke_r( 9)*RtimesT) )
         e13y4y20 = fkmt*cal_eta (r_kn, alphaHNO3)		     &
     		  *( y(4) - y(20)/(rke_r(13)*RtimesT) )
         e15y5y21 = fkmt*cal_eta (r_kn, alphaSO2)		     &
     		  *( y(5) - y(21)/(rke_r(15)*RtimesT) )
         e18y6y22 = fkmt*cal_eta (r_kn, alphaH2SO4)		     &
     		  *( y(6) - y(22)/(rke_r(18)*RtimesT) )
         e6y27y23 = fkmt*cal_eta (r_kn, alphaEST)		     &
     		  *( y(27)- y(23)/(rke_r( 6)*RtimesT) )
         e8y28y24 = fkmt*cal_eta (r_kn, alphaEST)		     &
     		  *( y(28)- y(24)/(rke_r( 8)*RtimesT) )
         e12y29y25= fkmt*cal_eta (r_kn, alphaEST)		     &
     		  *( y(29)- y(25)/(rke_r(12)*RtimesT) )
         e21y30y26= fkmt*cal_eta (r_kn, alphaEST)		     &
     		  *( y(30)- y(26)/(rke_r(21)*RtimesT) )
         e22y31y33= fkmt*cal_eta (c_kn, alphaNH3)		     &
     		  *( y(31)- y(33)/(rke_r(22)*RtimesT) )
     
         ydot( 1) = ydot( 1) - r_l*e11y1y17
         ydot( 2) = ydot( 2) - r_l*e3y2y18
         ydot( 3) = ydot( 3) - r_l*e9y3y19
         ydot( 4) = ydot( 4) - r_l*e13y4y20
         ydot( 5) = ydot( 5) - r_l*e15y5y21
         ydot( 6) = ydot( 6) - r_l*e18y6y22
               
         ydot(17) = e11y1y17 - rr2 - rr6 - rr10 - rr19
         ydot(18) = e3y2y18 + rr13
         ydot(19) = e9y3y19 + rr2 + rr7 + rr8			     &
     		  - rr1 - rr4 - rr20				     &
     		  - (rr14 + rr15)*xxx
         ydot(20) = e13y4y20
         ydot(21) = e15y5y21 - rr19 - rr20 - rr21 - rr22	     &
     		  - 2.0*(rr24 + rr25 + rr26 + rr27 + rr30)
         ydot(22) = e18y6y22 + rr19 + rr20 			     &
     		  + rr24 + rr25 + 2.0*(rr26 + rr27 + rr30)
         ydot(23) = e6y27y23 + 2.0*rr1 + rr10			     &
     		  - rr3 - rr4 - rr5 - rr6 - rr7 		     &
     		  - rr11 - rr12 - rr13 - rr15 - rr21 - rr22
         ydot(24) = e8y28y24 + (rr14 + rr15)*xxx  		     &
     		  + rr4 + rr6 + rr12 + rr13			     &
     		  - rr3 - rr5 - rr8 - rr10
	  ydot(25) = e12y29y25 + rr11 - rr12
	  ydot(26) = e21y30y26 - rr11

	  ydot(27) = ydot(27) - r_l*e6y27y23
	  ydot(28) = ydot(28) - r_l*e8y28y24
	  ydot(29) = ydot(29) - r_l*e12y29y25
	  ydot(30) = ydot(30) - r_l*e21y30y26	  
	  ydot(31) = ydot(31) - r_l*e22y31y33
          ydot(33) = e22y31y33
	endif	  	
#endif

	return
	 end

!	===============================================
	subroutine jex_aq (neq, t, y, ml, mu, pd, nrpd)
!	===============================================

! -----------------------------------------------------
! --- Subroutine for calculating the Jacobian for f(i)
! ---
! ---  P(m,n) here represents df(m)/dn
! --- 		or d ydot(m)/dn
! -----------------------------------------------------

	USE gridno
	USE shared_data

	integer :: neq, ml, mu, nrpd
	real    :: t, y(neq), pd(nrpd,neq)

#ifdef AQCHEM_ENABLE
	real    :: xxx, yyy, zzz, ddd,so5m2,so5m,fkmt,fff,              &
     		  r3y13,r3y14,r3y23,r3y24,r4y9,r4y13,r4y19,r4y23,	&
     		  r5y13,r5y14,r5y23,r5y24,r6y7,r6y13,r6y17,r6y23,	&
     		  r7y13,r7y23,r8y14,r8y24,r10y7,r10y14,r10y17,r10y24,	&
     		  r11y13,r11y16,r11y23,r11y26,r12y13,r12y15, 		&
     		  r12y23,r12y25,r14y14,r1417y9,r14y24,r1417y19, 	&
     		  r15y8,r15y13,r15y18,r15y23,r16y148,r16y2418,  	&
     		  r19y7,r19y11,r19y17,r19y21,r20y9,r20y11,		&
     		  r20y19,r20y21,r21y11,r21y13,r21y21,r21y23,		&
     		  r22y11,r22y13,r22y21,r22y23,r30y11,r30y13,		&
     		  r30y21,r30y23,	        			&
     		  e3y2,e3y8,e3y18,e9y3,e9y9,e9y19,			&
     		  e11y1,e11y7,e11y17,e13y4,e13y10,e13y20, 		&
     		  e15y5,e15y11,e15y21,e18y6,e18y12,e18y22,		&
     		  e22y31,e22y32,e22y33  			  
     
	! --- tmp space saving parameters:
	! --- hplus_t 1 and 2 for [H+] of cloud and rain,
	! --- c_kn and r_kn for Knudson number of cloud and rain,
	! --- c_ar and r_ar for mean r of cloud and rain,
	! --- c_l and r_l for mole fraction of lwc of cloud and rain
	! --- RtimesT = 0.08205*T
	real, dimension(2) :: hplus_t	
	common /hplustmp/hplus_t,c_kn,r_kn,c_ar,r_ar,c_l,r_l,RtimesT
	
! --------------------------------------------------
!  Table of Variables:
!	ychem( 1) = o3(g)
!	ychem( 2) = zco2(g)
!	ychem( 3) = h2o2(g)
!	ychem( 4) = hno3(g) + xn2o5(g)
!	ychem( 5) = so2(g)
!	ychem( 6) = h2so4(g)
!
!	ychem( 7) = o3_c
!	ychem( 8) = civ_c
!	ychem( 9) = h2o2_c
!	ychem(10) = xnv_c
!	ychem(11) = siv_c
!	ychem(12) = svi_c
!
!	ychem(13) = ho_c
!	ychem(14) = ho2_c
!	ychem(15) = ch2(oh)2_c
!	ychem(16) = ch3o2h_c
!
!	ychem(17) = o3_r
!	ychem(18) = civ_r
!	ychem(19) = h2o2_r
!	ychem(20) = xnv_r
!	ychem(21) = siv_r
!	ychem(22) = svi_r

!	ychem(23) = ho_r
!	ychem(24) = ho2_r
!	ychem(25) = ch2(oh)2_r
!	ychem(26) = ch3o2h_r
!
!	ychem(27) = ho(g)
!	ychem(28) = ho2(g)
!	ychem(29) = ch2o(g)
!	ychem(30) = ch3o2h(g)
!
!       ychem(31) = nh3(g)
!       ychem(32) = nh4_c
!       ychem(33) = nh4_r
! --------------------------------------------------

	do j=1,neq
	  if(y(j).le.0.0) y(j) = 0.0  
	do i=1,neq
	  pd(i,j) = 0.0
 	end do
 	end do

! ===
! === pd(i,j) = d ydot(i)/dy(j):
! ===

! --- For cloud water

	if(hplus_t(1).gt.1.e-14)then
	  fkmt = 3.0e-5/c_ar**2
	  	
	  ! === [SO5-]ss**2 and [SO5-]ss
	  if(rka_c(31).gt.0.0.and.y(11).gt.1.e-24.and.y(13).gt.1.e-24)then
	  	ddd  = (rka_c(21) + rka_c(22))/(2*rka_c(31))
	  	so5m2= ddd*y(11)*y(13)
	  	so5m = sqrt(so5m2)
	  else
	  	so5m2= 0.0
	  	so5m = 0.0
	  endif
	
	  ! === [CO3-]ss left RK17[H2O2]/(RK16[HO2] + RK17[H2O2])
	  zzz  = rka_c(16)*y(14) + rka_c(17)*y(9)
	  if(zzz.gt.1.e-16)then
	    zzz = 1./zzz
	  else
	    zzz = 1.e16
	  endif
	  xxx  = rka_c(17)*y(9)*zzz
	  yyy  = rka_c(17)*zzz**2
	  if(yyy.ge.1.e30) yyy = 1.e30

	  r3y13   = rka_c(3) *y(13)
	  r3y14   = rka_c(3) *y(14)
	  r4y9    = rka_c(4) *y(9)
	  r4y13   = rka_c(4) *y(13)
	  r5y13   = rka_c(5) *y(13)
	  r5y14   = rka_c(5) *y(14)
	  r6y7    = rka_c(6) *y(7)
	  r6y13   = rka_c(6) *y(13)
	  r7y13   = rka_c(7) *y(13)
	  r8y14   = rka_c(8) *y(14)
	  r10y7   = rka_c(10)*y(7)
	  r10y14  = rka_c(10)*y(14)
	  r11y13  = rka_c(11)*y(13)
	  r11y16  = rka_c(11)*y(16)
	  r12y13  = rka_c(12)*y(13)
	  r12y15  = rka_c(12)*y(15)
	  r14y14  = rka_c(14)*y(14)
	  r1417y9 = rka_c(14)*rka_c(17)*y(9)
	  r15y8   = rka_c(15)*y(8)
	  r15y13  = rka_c(15)*y(13)
	  r16y148 = rka_c(16)*y(14)*y(8)
	  r19y7   = rka_c(19)*y(7)
	  r19y11  = rka_c(19)*y(11)
	  r20y9   = rka_c(20)*y(9)
	  r20y11  = rka_c(20)*y(11)
	  r21y11  = rka_c(21)*y(11)
	  r21y13  = rka_c(21)*y(13)
	  r22y11  = rka_c(22)*y(11)
	  r22y13  = rka_c(22)*y(13)
	  r30y11  = rka_c(30)*y(11)
	  r30y13  = rka_c(30)*y(13)	

	  e3y2    = fkmt*cal_eta (c_kn, alphaCO2)
	  e3y8    = e3y2/(rke_c( 3)*RtimesT)
	  e9y3    = fkmt*cal_eta (c_kn, alphaH2O2)
	  e9y9    = e9y3/(rke_c( 9)*RtimesT)
	  e11y1   = fkmt*cal_eta (c_kn, alphaO3)
	  e11y7   = e11y1/(rke_c(11)*RtimesT)
	  e13y4   = fkmt*cal_eta (c_kn, alphaHNO3)
	  e13y10  = e13y4/(rke_c(13)*RtimesT)
	  e15y5   = fkmt*cal_eta (c_kn, alphaSO2)
	  e15y11  = e15y5/(rke_c(15)*RtimesT)
	  e18y6   = fkmt*cal_eta (c_kn, alphaH2SO4)
	  e18y12  = e18y6/(rke_c(18)*RtimesT)
	  e6y27   = fkmt*cal_eta (c_kn, alphaEST)
	  e6y13   = e6y27/(rke_c( 6)*RtimesT)
     	  e8y28   = fkmt*cal_eta (c_kn, alphaEST)
	  e8y14   = e8y28/(rke_c( 8)*RtimesT)
     	  e12y29  = fkmt*cal_eta (c_kn, alphaEST)
	  e12y15  = e12y29/(rke_c(12)*RtimesT)
     	  e21y30  = fkmt*cal_eta (c_kn, alphaEST)
	  e21y16  = e21y30/(rke_c(21)*RtimesT)
	  e22y31  = fkmt*cal_eta (c_kn, alphaNH3)
	  e22y32  = e22y31/(rke_c(22)*RtimesT)


       	  pd(1,1)  = - c_l*e11y1
       	  pd(1,7)  =   c_l*e11y7
       	  pd(2,2)  = - c_l*e3y2
       	  pd(2,8)  =   c_l*e3y8
       	  pd(3,3)  = - c_l*e9y3
       	  pd(3,9)  =   c_l*e9y9
       	  pd(4,4)  = - c_l*e13y4
       	  pd(4,10) =   c_l*e13y10
       	  pd(5,5)  = - c_l*e15y5
       	  pd(5,11) =   c_l*e15y11
	  pd(6,6)  = - c_l*e18y6
          pd(6,12) =   c_l*e18y12
     		
	  pd(7,1)  = e11y1
	  pd(7,7)  = - (rka_c(2) + e11y7 + r6y13 + r10y14 + r19y11)
	  pd(7,11) = - r19y7
	  pd(7,13) = - r6y7
	  pd(7,14) = - r10y7
	  	
	  pd(8,2)  = e3y2
	  pd(8,8)  = - e3y8
	  pd(8,13) = r12y15
	  pd(8,15) = r12y13
	
	  pd(9,3)  = e9y3
	  pd(9,7)  = rka_c(2)
	  pd(9,8)  = - xxx*(r14y14 + r15y13)
	  pd(9,9)  = - e9y9 - rka_c(1)                           &
     		    - r4y13 - r20y11				 &
     		    - yyy*r16y148*(r14y14 + r15y13)
         pd(9,11) = - r20y9
         pd(9,13) = 2.0*r7y13 - r4y9 - xxx*r15y8
         pd(9,14) = 2.0*r8y14 					 &
     		  - yyy*(r1417y9 + rka_c(16)*r15y13)*y(8)*y(9)
              
         pd(10,4) = e13y4
         pd(10,10)= - e13y10

         pd(11,5) = e15y5
         pd(11,7) = - r19y11
         pd(11,9) = - r20y11
         pd(11,11)= - e15y11 - r19y7 - r20y9 - r21y13 - r22y13
         pd(11,13)= - r21y11 - r22y11
         if(y(13).gt.0.0)then
           pd(11,11) = pd(11,11)				 &
     		     - 3.0*(rka_c(24) + rka_c(25) 		 &
     			  + rka_c(26) + rka_c(27))*so5m 	 &
     		     - 2.0*ddd*r30y13				  
     	   pd(11,13) = pd(11,13)				 &
     		     - y(11)/y(13)*(rka_c(24) + rka_c(25) 	 &
     				  + rka_c(26) + rka_c(27))*so5m  &
     		     - 2.0*ddd*r30y11
         endif
         
         pd(12,6) = e18y6
         pd(12,7) = r19y11
         pd(12,9) = r20y11
         pd(12,11)= r19y7 + r20y9
         pd(12,12)= - e18y12
         if(y(13).gt.0.0)then
           pd(12,11) = pd(12,11)				 &
     		     + 1.5*(rka_c(24) + rka_c(25))*so5m 	 &
     		     + 3.0*(rka_c(26) + rka_c(27))*so5m 	 &
     		     + 2.0*ddd*r30y13				  
     	   pd(12,13) = y(11)/y(13)				 &
     		      *(0.5*(rka_c(24) + rka_c(25))		 &
     			  +  rka_c(26) + rka_c(27) )*so5m	 &
     		     + 2.0*ddd*r30y11
         endif
        	 
         pd(13,7) = r10y14 - r6y13
         pd(13,8) = - r15y13
         pd(13,9) = 2.0*rka_c(1) - r4y13
         pd(13,11)= - r21y13 - r22y13
         pd(13,13)= - r3y14 - r5y14 - r4y9 - r6y7 - 2.0*r7y13	 &
     		    - r11y16 - r12y15 - rka_c(13)*y(15) 	 &
     		    - r15y8 - r21y11 - r22y11			 &
     		    - e6y13
         pd(13,14)= r10y7 - r3y13 - r5y13			  
     	 pd(13,15)= - r12y13 - rka_c(13)*y(13)
         pd(13,16)= - r11y13
         pd(13,27)= e6y27

         pd(14,7) = r6y13 - r10y14
         pd(14,8) = xxx*(r14y14 + r15y13)
         pd(14,9) = yyy*r16y148*(r14y14 + r15y13) + r4y13
         pd(14,13)= xxx*r15y8 + r4y9 + r6y7 + 2.0*r12y15	 &
     		  - r3y14 - r5y14
         pd(14,14)= yyy*(r1417y9 + rka_c(16)*r15y13)*y(8)*y(9)	 &
     		  - r3y13 - r5y13 - 2.0*r8y14 - r10y7		 &
     		  - e8y14
	  pd(14,15)= 2.0*r12y13
	  pd(14,28)= e8y28
	
	  pd(15,13)= r11y16 - r12y15
	  pd(15,15)= - e12y15 - r12y13
	  pd(15,16)= r11y13
	  pd(15,29)= e12y29
	
	  pd(16,13)= - r11y16
	  pd(16,16)= - e21y16 - r11y13
	  pd(16,30)= e21y30
	  
	  pd(27,27)= - c_l*e6y27
       	  pd(27,13)=   c_l*e6y13
       	  pd(28,28)= - c_l*e8y28
       	  pd(28,14)=   c_l*e8y14
       	  pd(29,29)= - c_l*e12y29
       	  pd(29,15)=   c_l*e12y15
	  pd(30,30)= - c_l*e21y30
          pd(30,16)=   c_l*e21y16

          pd(31,31)= - c_l*e22y31
          pd(31,32)=   c_l*e22y32

          pd(32,31)=   e22y31
          pd(32,32)=  -e22y32
	endif

! --- For rain water

	if(hplus_t(2).gt.1.e-14)then
	  fkmt = 3.0e-5/r_ar**2
	  
	  ! === [SO5-]ss**2 and [SO5-]ss
	  if(rka_r(31).gt.0.0.and.y(21).gt.1.e-24.and.y(23).gt.1.e-24)then
	  	ddd  = (rka_r(21) + rka_r(22))/(2*rka_r(31))
	  	so5m2= ddd*y(21)*y(23)
	  	so5m = sqrt(so5m2)
	  else
	  	so5m2= 0.0
	  	so5m = 0.0
	  endif

	  ! === [CO3-]ss left RK17[H2O2]/(RK16[HO2] + RK17[H2O2])
	  zzz  = rka_r(16)*y(24) + rka_r(17)*y(19)
	  if(zzz.gt.1.e-16)then
	    zzz = 1./zzz
	  else
	    zzz = 1.e16
	  endif
	  xxx  = rka_r(17)*y(19)*zzz
	  yyy  = rka_r(17)*zzz**2
	  if(yyy.ge.1.e30) yyy = 1.e30

	  r3y23   = rka_r(3) *y(23)
	  r3y24   = rka_r(3) *y(24)
	  r4y19   = rka_r(4) *y(19)
	  r4y23   = rka_r(4) *y(23)
	  r5y23   = rka_r(5) *y(23)
	  r5y24   = rka_r(5) *y(24)
	  r6y17   = rka_r(6) *y(17)
	  r6y23   = rka_r(6) *y(23)
	  r7y23   = rka_r(7) *y(23)
	  r8y24   = rka_r(8) *y(24)
	  r10y17  = rka_r(10)*y(17)
	  r10y24  = rka_r(10)*y(24)
	  r11y23  = rka_r(11)*y(23)
	  r11y26  = rka_r(11)*y(26)
	  r12y23  = rka_r(12)*y(23)
	  r12y25  = rka_r(12)*y(25)
	  r14y24  = rka_r(14)*y(24)
	  r1417y19= rka_r(14)*rka_r(17)*y(19)
	  r15y18  = rka_r(15)*y(18)
	  r15y23  = rka_r(15)*y(23)
	  r16y2418= rka_r(16)*y(24)*y(18)
	  r19y17  = rka_r(19)*y(17)
	  r19y21  = rka_r(19)*y(21)
	  r20y19  = rka_r(20)*y(19)
	  r20y21  = rka_r(20)*y(21)
	  r21y21  = rka_r(21)*y(21)
	  r21y23  = rka_r(21)*y(23)
	  r22y21  = rka_r(22)*y(21)
	  r22y23  = rka_r(22)*y(23)
	  r30y21  = rka_r(30)*y(21)
	  r30y23  = rka_r(30)*y(23)

	  e3y2    = fkmt*cal_eta (r_kn, alphaCO2)
	  e3y18    = e3y2/(rke_r( 3)*RtimesT)
	  e9y3    = fkmt*cal_eta (r_kn, alphaH2O2)
	  e9y19    = e9y3/(rke_r( 9)*RtimesT)
	  e11y1   = fkmt*cal_eta (r_kn, alphaO3)
	  e11y17   = e11y1/(rke_r(11)*RtimesT)
	  e13y4   = fkmt*cal_eta (r_kn, alphaHNO3)
	  e13y20  = e13y4/(rke_r(13)*RtimesT)
	  e15y5   = fkmt*cal_eta (r_kn, alphaSO2)
	  e15y21  = e15y5/(rke_r(15)*RtimesT)
	  e18y6   = fkmt*cal_eta (r_kn, alphaH2SO4)
	  e18y22  = e18y6/(rke_r(18)*RtimesT)
	  e6y27   = fkmt*cal_eta (r_kn, alphaEST)
	  e6y23   = e6y27/(rke_r( 6)*RtimesT)
     	  e8y28   = fkmt*cal_eta (r_kn, alphaEST)
	  e8y24   = e8y28/(rke_r( 8)*RtimesT)
     	  e12y29  = fkmt*cal_eta (r_kn, alphaEST)
	  e12y25  = e12y29/(rke_r(12)*RtimesT)
     	  e21y30  = fkmt*cal_eta (r_kn, alphaEST)
	  e21y26  = e21y30/(rke_r(21)*RtimesT)
	  e22y31   = fkmt*cal_eta (r_kn, alphaCO2)
	  e22y33   = e22y31/(rke_r(22)*RtimesT)


       	  pd(1,1)  = pd(1,1) - r_l*e11y1
       	  pd(1,17) = r_l*e11y17
       	  pd(2,2)  = pd(2,2) - r_l*e3y2
       	  pd(2,18) = r_l*e3y18
       	  pd(3,3)  = pd(3,3) - r_l*e9y3
       	  pd(3,19) = r_l*e9y19
       	  pd(4,4)  = pd(4,4) - r_l*e13y4
       	  pd(4,20) = r_l*e13y20
       	  pd(5,5)  = pd(5,5) - r_l*e15y5
       	  pd(5,21) = r_l*e15y21
	  pd(6,6)  = pd(6,6) - r_l*e18y6
          pd(6,22) = r_l*e18y22
	  	
	  pd(17,1)  = e11y1
	  pd(17,17) = - (rka_r(2) + e11y17 + r6y23 + r10y24 + r19y21) 
	  pd(17,21) = - r19y17
	  pd(17,23) = - r6y17
	  pd(17,24) = - r10y17
	  	
	  pd(18,2)  = e3y2
	  pd(18,18) = - e3y18
	  pd(18,23) = r12y25
	  pd(18,25) = r12y23
	  	
	  pd(19,3)  = e9y3
	  pd(19,17) = rka_r(2)
	  pd(19,18) = - xxx*(r14y24 + r15y23)
	  pd(19,19) = - e9y19 - rka_r(1) - r4y23 - r20y21            &
     		     - yyy*r16y2418*(r14y24 + r15y23)
         pd(19,23) = 2.0*r7y23 - r4y19 - xxx*r15y18
         pd(19,21) = - r20y19
         pd(19,24) = 2.0*r8y24  				     &
     		   - yyy*(r1417y19 + rka_r(16)*r15y23)*y(18)*y(19)

         pd(20,4)  = e13y4
         pd(20,20) = - e13y20

         pd(21,5)  = e15y5
         pd(21,17) = - r19y21
         pd(21,19) = - r20y21
         pd(21,21) = - e15y21 - r19y17 - r20y19 - r21y23 - r22y23
         pd(21,23) = - r21y21 - r22y21
         if(y(23).gt.0.0)then
           pd(21,21) = pd(21,21)				     &
     		     - 3.0*(rka_r(24) + rka_r(25)		     &
     			  + rka_r(26) + rka_r(27))*so5m 	     &
     		     - 2.0*ddd*r30y23				      
     	   pd(21,23) = pd(21,23)				     &
     		     - y(21)/y(23)*(rka_r(24) + rka_r(25)	     &
     				  + rka_r(26) + rka_r(27))*so5m      &
     		     - 2.0*ddd*r30y21
         endif
               
         pd(22,6)  = e18y6
         pd(22,17) = r19y21
         pd(22,19) = r20y21
         pd(22,21) = r19y17 + r20y19
         pd(22,22) = - e18y22
         if(y(23).gt.0.0)then
           pd(22,21) = pd(22,21)				     &
     		     + 1.5*(rka_c(24) + rka_c(25))*so5m 	     &
     		     + 3.0*(rka_r(26) + rka_r(27))*so5m 	     &
     		     + 2.0*ddd*r30y23				      
     	   pd(22,23) = y(21)/y(23)				     &
     		      *(0.5*(rka_c(24) + rka_c(25))		     &
     			  +  rka_c(26) + rka_c(27) )*so5m	     &
     		     + 2.0*ddd*r30y21
         endif

         pd(23,17) = r10y24 - r6y23
         pd(23,18) = - r15y23
         pd(23,19) = 2.0*rka_r(1) - r4y23
         pd(23,21) = - r21y23 - r22y23
         pd(23,23) = - r3y24 - r5y24 - r4y19 - r6y17 - 2.0*r7y23     &
     		     - r11y26 - r12y25 - rka_r(13)*y(25)	     &
     		     - r15y18 - r21y21 - r22y21 		     &
     		     - e6y23
	  pd(23,24) = r10y17 - r3y23 - r5y23
          pd(23,25) = - r12y23 - rka_r(13)*y(23)
	  pd(23,26) = - r11y23
	  pd(23,27) = e6y27
	
	  pd(24,17) = r6y23 - r10y24
	  pd(24,18) = xxx*(r14y24 + r15y23)
	  pd(24,19) = yyy*r16y2418*(r14y24 + r15y23) + r4y23
	  pd(24,23) = xxx*r15y18 + r4y19 + r6y17 + 2.0*r12y25      &
     		   - r3y24 - r5y24
         pd(24,24) = yyy*(r1417y19 + rka_r(16)*r15y23)*y(18)*y(19) &
     		   - r3y23 - r5y23 - 2.0*r8y24 - r10y17 	   &
     		   - e8y24
	  pd(24,25) = 2.0*r12y23
	  pd(24,28) = e8y28
	
	  pd(25,23) = r11y26 - r12y25
	  pd(25,25) = - e12y25 - r12y23
	  pd(25,26) = r11y23
	  pd(25,29) = e12y29
	
	  pd(26,23) = - r11y26
	  pd(26,26) = - e21y26 - r11y23
	  pd(26,30) = e21y30
	  
       	  pd(27,27) = pd(27,27) - r_l*e6y27
       	  pd(27,23) =   r_l*e6y23
       	  pd(28,28) = pd(28,28) - r_l*e8y28
       	  pd(28,24) =   r_l*e8y24
       	  pd(29,29) = pd(29,29) - r_l*e12y29
       	  pd(29,25) =   r_l*e12y25
	  pd(30,30) = pd(30,30) - r_l*e21y30
          pd(30,26) =   r_l*e21y26

          pd(31,31) = pd(31,31) - r_l*e22y31
          pd(31,33) =   r_l*e22y33

          pd(33,31) = e22y31
          pd(33,33) = - e22y33 
	endif
#endif 

	return
	 end
      	 	 
!	==================================
	subroutine rateaq (T, patm, Rflux)
!	==================================

	USE gridno
	USE shared_data

	real    :: T, patm, Rflux
	
#ifdef AQCHEM_ENABLE
	real, dimension(nequil)  :: rke
	real, dimension(naqreact):: rka

	real    :: f1, f2, f3, f4, f5, f6, f7, f8, f9, f10
        real    :: f11, fsiv, T298, est1
	
	! --- tmp space saving parameters:
	! --- hplus_t 1 and 2 for [H+] of cloud and rain,
	! --- c_kn and r_kn for Knudson number of cloud and rain,
	! --- c_ar and r_ar for mean r of cloud and rain,
	! --- c_l and r_l for mole fraction of lwc of cloud and rain
	! --- RtimesT = 0.08205*T
	real, dimension(2) :: hplus_t	
	common /hplustmp/hplus_t,c_kn,r_kn,c_ar,r_ar,c_l,r_l,RtimesT

	do iii=1,nequil
	  rke  (iii) = 0.0
	  rke_c(iii) = 0.0
	  rke_r(iii) = 0.0
	end do
	
	do iii=1,naqreact
	  rka  (iii) = 0.0
	  rka_c(iii) = 0.0
	  rka_r(iii) = 0.0
	end do 
			
! === 1/T - 1/298 = (298 - T)/(298T)

	T298 = (298.0 - T)/(298.0*T)
	  
! === Equilibrium Reactions

	! --- RE1	O2(g) <=> O2(aq)
	rke(1)  = 1.30e-3*exp(1500.0*T298)

	! --- RE2	H2O <=> H+ + OH-
	rke(2)  = 1.00e-14*exp(-6716.0*T298) 

	! --- RE3	CO2(g) <=> CO2(aq)
	rke(3)  = 3.40e-2*exp(2420.0*T298)

	! --- RE4	CO2(aq) <=> HCO3- + H+ 
	rke(4)  = 4.46e-7*exp(-1000.0*T298) 

	! --- RE5	HCO3- <=> CO32- + H+ 
	rke(5)  = 4.48e-11*exp(-1760.0*T298)

	! --- RE6	OH(g) <=> OH(aq)
	rke(6)  = 2.50e1*exp(5280.0*T298)

	! --- RE7	HO2(g) <=> HO2(aq)
	rke(7)  = 2.30e3*exp(6640.0*T298)

	! --- RE8	HO2(aq) <=> O2- + H+
	rke(8)  = 3.50e-5

	! --- RE9	H2O2(g) <=> H2O2 (aq)
	rke(9)  = 7.40e4*exp(6615.0*T298)
	
	! --- RE10	H2O2(aq) <=> HO2- + H+
	rke(10) = 2.20e-12*exp(-3730.0*T298)

	! --- RE11	O3(g) <=> O3(aq)
	rke(11) = 1.13e-2*exp(2300.0*T298)

	! --- RE12	CH2O(g) <=> CH2(OH)2
	rke(12) = 3.10e3*exp(6500.0*T298)

	! --- RE13	HNO3(g) <=> HNO3 (aq)
	rke(13) = 2.10e5*exp(8700.0*T298)

	! --- RE14	HNO3(aq) <=> NO3- + H+
	rke(14) = 1.54e1*exp(8700.0*T298)

	! --- RE15	SO2(g) <=> SO2(aq)
	rke(15) = 1.20*exp(3100.0*T298)

	! --- RE16	SO2(aq) <=> HSO3- + H+
	rke(16) = 1.23e-2*exp(1960.0*T298)

	! --- RE17	HSO3- <=> SO32- + H+
	rke(17) = 6.61e-8*exp(1500.0*T298)

	! --- RE18	H2SO4(g) <=> H2SO4(aq) set to be large enough
	rke(18) = 1.0e8

	! --- RE19	H2SO4(aq) <=> HSO4- + H+
	rke(19) = 1.00e3

	! --- RE20	HSO4- <=> SO42- + H+
	rke(20) = 1.02e-2*exp(1960.0*T298)

	! --- RE21	CH3O2H(g) <=> CH3O2H(aq)
	rke(21) = 3.10e2*exp(5200.0*T298)

	! --- RE22	NH3(g) <=> NH3 (aq)
	rke(22)  = 62.0*exp(4119.0*T298)
	
	! --- RE23	NH3(aq) <=> NH4+ + OH-
	rke(23) = 1.7e-5*exp(-450.0*T298)

! === exp(-1500.0*T298) is a common factor based estimate
	est1 = exp(-1500.0*T298)

! === Aqueous Phase Reactions

	! --- RA1	H2O2 + hv -> 2OH
	rka(1)  = 0.3*rk(17)	! 30% of gaseous reaction
	
	! --- RA2	O3 + hv (+H2O) -> H2O2
	rka(2)  = 0.8*rk(1)	! 80% of gaseous reaction
		
	! --- RA3	OH + HO2 -> H2O + O2
	rka(3)  = 7.0e9*est1
	
	! --- RA4	OH + H2O2 -> H2O + HO2
	rka(4)  = 2.7e7*exp(-1700.0*T298)

	! --- RA5	OH + O2- -> OH- + O2
	rka(5)  = 1.0e10*est1

	! --- RA6	OH + O3 -> HO2 + O2
	rka(6)  = 2.0e9

	! --- RA7	OH + OH -> H2O2
	rka(7)  = 5.2e9*est1

	! --- RA8	HO2 + O2- -> HO2- + O2
	rka(8)  = 1.0e8*est1

	! --- RA9	HO2- + H+ -> H2O2
	rka(9)  = 5.0e10*est1

	! --- RA10	O2- + O3 (+ H2O) -> OH + 2O2 + OH-
	rka(10) = 1.5e9*est1

	! --- RA11	CH3O2H + OH -> CH2(OH)2 + OH
	rka(11) = 2.7e7*exp(-1700.0*T298)

	! --- RA12	CH2(OH)2 + OH (+ O2) -> H2O + HCOOH + HO2
	rka(12) = 2.0e9*est1

	! --- RA13	HCOOH + OH (+ O2) -> CO2 + H2O + HO2
	rka(13) = 2.0e8*est1

	! --- RA14	HCO3- + O2- -> HO2- + CO3-
	rka(14) = 1.5e6

	! --- RA15	HCO3- + OH -> H2O + CO3-
	rka(15) = 1.5e7*exp(-1910.0*T298)

	! --- RA16	CO3- + O2- (+ H2O) -> HCO3- + O2 + OH-
	rka(16) = 4.0e8*est1

	! --- RA17	CO3- + H2O2 -> HO2 + HCO3-
	rka(17) = 8.0e5*exp(-2820.0*T298)

	! --- RA18	N2O5(g) + H2O(aq) -> 2HNO3
	rka(18) = 0.0

	! --- RA19	S(IV) + O3 -> S(VI) + O2
	rka(19) = 0.0	! to be set later

	! --- RA20	S(IV) + H2O2 -> S(VI) + H2O
	rka(20) = 7.5e7*exp(-4430.0*T298)

	! --- RA21	HSO3- + OH -> SO3- + H2O
	rka(21) = 4.5e9*est1

	! --- RA22	SO32- + OH -> SO3- + OH-
	rka(22) = 5.2e9*est1

	! --- RA23	SO3- + O2 -> SO5-
	rka(23) = 1.5e9*est1

	! --- RA24	SO5- + HSO3- -> HSO5- + SO3-
	rka(24) = 2.5e4*exp(-3100.0*T298)

	! --- RA25	SO5- + SO32- (+ H2O) -> HSO5- + SO3- + OH-
	rka(25) = 2.5e4*exp(-2000.0*T298)

	! --- RA26	SO5- + HSO3- -> SO4- + SO42- + H+
	rka(26) = 7.5e4*exp(-3100.0*T298)

	! --- RA27	SO5- + SO32- -> SO4-  + SO42-
	rka(27) = 7.5e4*exp(-2000.0*T298)

	! --- RA28	SO4- + HSO3- -> SO42- + SO3- + H+
	rka(28) = 7.5e8*est1

	! --- RA29	SO4- + SO32- -> SO42- + SO3-
	rka(29) = 5.5e8*est1

	! --- RA30	SO5- + SO5- -> 2SO4- + O2
	rka(30) = 6.0e8*est1

	! --- RA31	SO5- + SO5- -> S2O82- + O2
	rka(31) = 1.4e8*est1

	! --- RA32	HSO5- + HSO3- (+ H+) -> 2 SO42- + 3H+
	rka(32) = 7.1e6*exp(-3100.0*T298)


! === For cloud droplets:
! ---	
! --- conversion factors f = 1/f
! --- fsiv= 1./([H+]*([H+] + re16) + re16*re17)	
! --- f6  = 1./(1.0 + [H+]/re8)
! --- f7  = 1./(1.0 + [H+]/re10)
! --- f8  = 1./(1.0 + [H+]/re4  + re5/[H+])
! --- f9  = 1./(1.0 + [H+]/re16 + re17/[H+])
! --- f10 = 1./(1.0 + [H+]/re17 + [H+]**2/(re16*re17))
! ---
	rke_c( 1) = rke( 1)
	rke_c( 2) = rke( 2)
	rke_c( 4) = rke( 4)
	rke_c( 5) = rke( 5)
	rke_c( 6) = rke( 6)
	rke_c( 8) = rke( 8)
	rke_c(10) = rke(10)
	rke_c(11) = rke(11)
	rke_c(12) = rke(12)
	rke_c(14) = rke(14)
	rke_c(16) = rke(16)
	rke_c(17) = rke(17)
	rke_c(18) = rke(18)
	rke_c(19) = rke(19)
	rke_c(20) = rke(20)
	rke_c(21) = rke(21)
	rke_c(22) = rke(22)
	rke_c(23) = rke(23)
	
	rka_c( 1) = rka( 1)
	rka_c( 2) = rka( 2)
	rka_c( 3) = rka( 3)
	rka_c( 4) = rka( 4)
	rka_c( 6) = rka( 6)
	rka_c( 7) = rka( 7)
	rka_c(11) = rka(11)
	rka_c(12) = rka(12)
	rka_c(13) = rka(13)
	rka_c(17) = rka(17)
	rka_c(18) = rka(18)
	rka_c(23) = rka(23)
	rka_c(30) = rka(30)
	rka_c(31) = rka(31)
	
	if(hplus_t(1).gt.1.e-14)then	                    
	  fsiv= 1.0/( hplus_t(1)*(hplus_t(1) + rke_c(16))    &
     				  + rke_c(16)*rke_c(17) )   
         f1  = (1.0 + rke_c(4)/hplus_t(1) 		     &
     		    + rke_c(4)*rke_c(5)/hplus_t(1)**2)
         f2  = (1.0 + rke_c(8)/hplus_t(1))
         f3  = (1.0 + rke_c(10)/hplus_t(1))
         f4  = (1.0 + rke_c(14)/hplus_t(1))		    
         f5  = (1.0 + rke_c(16)/hplus_t(1) 		     &
     		    + rke_c(16)*rke_c(17)/hplus_t(1)**2)        
         f6  = rke_c(8)/(rke_c(8) + hplus_t(1))               
         f7  = rke_c(10)/(rke_c(10) + hplus_t(1))	        
         f8  = rke_c(4)*hplus_t(1)			     &
     	      /( hplus_t(1)*(hplus_t(1) + rke_c(4)) 	     &
     				+rke_c(4)*rke_c(5) )	        			
         f9  = rke_c(16)*hplus_t(1)*fsiv	            
         f10 = rke_c(16)*rke_c(17)*fsiv 		      
     	 f11 = rke_c(22)*(1+rke_c(23)*hplus_t(1)/rke_c(2))

         rke_c( 3) = rke( 3)*f1
         rke_c( 7) = rke( 7)*f2
         rke_c( 9) = rke( 9)*f3
         rke_c(13) = rke(13)*f4 			    
         rke_c(15) = rke(15)*f5 			      
     	 rke_c(22) = rke(22)*f11
        	 
         rka_c( 5) = rka( 5)*f6
         rka_c( 8) = rka( 8)*f6
         rka_c( 9) = rka( 9)*f7        
         rka_c(10) = rka(10)*f6
         rka_c(14) = rka(14)*f8*f6
         rka_c(15) = rka(15)*f8
         rka_c(16) = rka(16)*f6 			    
         rka_c(19) = ( hplus_t(1)**2*2.4e4		     &
     		     + rke_c(16)*hplus_t(1)		     &
     		      *3.7e5*exp(-5530.0*T)		     &
     		     + rke_c(16)*rke_c(17)		     &
     		      *1.5e9*exp(-5280.0*T) )*fsiv	    
         rka_c(20) = rka(20)*rke_c(16)*hplus_t(1)**2	     &
     		 /(1.0 + 13.0*hplus_t(1))*fsiv
	  rka_c(21) = rka(21)*f9
	  rka_c(22) = rka(22)*f10
	  rka_c(24) = rka(24)*f9
	  rka_c(25) = rka(25)*f10
	  rka_c(26) = rka(26)*f9
	  rka_c(27) = rka(27)*f10
	  rka_c(28) = rka(28)*f9
	  rka_c(29) = rka(29)*f10
	  rka_c(32) = rka(32)*f9
	else
	  rke_c( 3) = rke( 3)
	  rke_c( 7) = rke( 7)
	  rke_c( 9) = rke( 9)
	  rke_c(13) = rke(13)
	  rke_c(15) = rke(15)
          rke_c(22) = rke(22)
	  
	  rka_c( 5) = 0.0
	  rka_c( 8) = 0.0
	  rka_c( 9) = 0.0	
	  rka_c(10) = 0.0
	  rka_c(14) = 0.0
	  rka_c(15) = 0.0
	  rka_c(16) = 0.0
	  rka_c(19) = 0.0
	  rka_c(20) = 0.0
	  rka_c(21) = 0.0
	  rka_c(22) = 0.0
	  rka_c(24) = 0.0
	  rka_c(25) = 0.0
	  rka_c(26) = 0.0
	  rka_c(27) = 0.0
	  rka_c(28) = 0.0
	  rka_c(29) = 0.0
	  rka_c(32) = 0.0
	endif
	
! === For rain drops:
! ---	
! --- conversion factors f = 1/f
! --- fsiv= 1./([H+]*([H+] + re16) + re16*re17)	
! --- f6  = 1./(1.0 + [H+]/re8)
! --- f7  = 1./(1.0 + [H+]/re10)
! --- f8  = 1./(1.0 + [H+]/re4  + re5/[H+])
! --- f9  = 1./(1.0 + [H+]/re16 + re17/[H+])
! --- f10 = 1./(1.0 + [H+]/re17 + [H+]**2/(re16*re17))
! ---
	rke_r( 1) = rke( 1)
	rke_r( 2) = rke( 2)
	rke_r( 4) = rke( 4)
	rke_r( 5) = rke( 5)
	rke_r( 6) = rke( 6)
	rke_r( 8) = rke( 8)
	rke_r(10) = rke(10)
	rke_r(11) = rke(11)
	rke_r(12) = rke(12)
	rke_r(14) = rke(14)
	rke_r(16) = rke(16)
	rke_r(17) = rke(17)
	rke_r(18) = rke(18)
	rke_r(19) = rke(19)
	rke_r(20) = rke(20)
	rke_r(21) = rke(21)
	rke_r(22) = rke(22)
	rke_r(23) = rke(23)
	
	rka_r( 1) = rka( 1)
	rka_r( 2) = rka( 2)
	rka_r( 3) = rka( 3)
	rka_r( 4) = rka( 4)
	rka_r( 6) = rka( 6)
	rka_r( 7) = rka( 7)
	rka_r(11) = rka(11)
	rka_r(12) = rka(12)
	rka_r(13) = rka(13)
	rka_r(17) = rka(17)
	rka_r(18) = rka(18)
	rka_r(23) = rka(23)
	rka_r(30) = rka(30)
	rka_r(31) = rka(31)
	
	if(hplus_t(2).gt.1.e-14)then                         
	  fsiv= 1.0/( hplus_t(2)*(hplus_t(2) + rke_r(16))    &
     				  + rke_r(16)*rke_r(17) )    
         f1  = (1.0 + rke_r(4)/hplus_t(2) 		     &
     		    + rke_r(4)*rke_r(5)/hplus_t(2)**2)
         f2  = (1.0 + rke_r(8)/hplus_t(2))
         f3  = (1.0 + rke_r(10)/hplus_t(2))
         f4  = (1.0 + rke_r(14)/hplus_t(2))		     
         f5  = (1.0 + rke_r(16)/hplus_t(2) 		     &
     		    + rke_r(16)*rke_r(17)/hplus_t(2)**2)        
         f6  = rke_r(8)/(rke_r(8) + hplus_t(2))               
         f7  = rke_r(10)/(rke_r(10) + hplus_t(2))	        
         f8  = rke_r(4)*hplus_t(2)			     &
     	      /( hplus_t(2)*(hplus_t(2) + rke_r(4)) 	     &
     				+rke_r(4)*rke_r(5) )	        			
         f9  = rke_r(16)*hplus_t(2)*fsiv	             
         f10 = rke_r(16)*rke_r(17)*fsiv 		      
     	 f11 = rke_r(22)*(1+rke_r(23)*hplus_t(1)/rke_r(2))

         rke_r( 3) = rke( 3)*f1
         rke_r( 7) = rke( 7)*f2
         rke_r( 9) = rke( 9)*f3
         rke_r(13) = rke(13)*f4
         rke_r(15) = rke(15)*f5
         rke_r(22) = rke(22)*f11
        	 
         rka_r( 5) = rka( 5)*f6
         rka_r( 8) = rka( 8)*f6
         rka_r( 9) = rka( 9)*f7        
         rka_r(10) = rka(10)*f6
         rka_r(14) = rka(14)*f8*f6
         rka_r(15) = rka(15)*f8
         rka_r(16) = rka(16)*f6 			     
         rka_r(19) = ( hplus_t(2)**2*2.4e4		     &
     		     + rke_r(16)*hplus_t(2)		     &
     		      *3.7e5*exp(-5530.0*T)		     &
     		     + rke_r(16)*rke_r(17)		     &
     		      *1.5e9*exp(-5280.0*T) )*fsiv	     
         rka_r(20) = rka(20)*rke_r(16)*hplus_t(2)**2	     &
     		 /(1.0 + 13.0*hplus_t(2))*fsiv
	  rka_r(21) = rka(21)*f9
	  rka_r(22) = rka(22)*f10
	  rka_r(24) = rka(24)*f9
	  rka_r(25) = rka(25)*f10
	  rka_r(26) = rka(26)*f9
	  rka_r(27) = rka(27)*f10
	  rka_r(28) = rka(28)*f9
	  rka_r(29) = rka(29)*f10
	  rka_r(32) = rka(32)*f9
	else
	  rke_r( 3) = rke( 3)
	  rke_r( 7) = rke( 7)
	  rke_r( 9) = rke( 9)
	  rke_r(13) = rke(13)
	  rke_r(15) = rke(15)
	  rke_r(22) = rke(22)
	  
	  rka_r( 5) = 0.0
	  rka_r( 8) = 0.0
	  rka_r( 9) = 0.0	
	  rka_r(10) = 0.0
	  rka_r(14) = 0.0
	  rka_r(15) = 0.0
	  rka_r(16) = 0.0
	  rka_r(19) = 0.0
	  rka_r(20) = 0.0
	  rka_r(21) = 0.0
	  rka_r(22) = 0.0
	  rka_r(24) = 0.0
	  rka_r(25) = 0.0
	  rka_r(26) = 0.0
	  rka_r(27) = 0.0
	  rka_r(28) = 0.0
	  rka_r(29) = 0.0
	  rka_r(32) = 0.0
	endif
#endif 
	
	return
	 end

!	=============================
	real function cal_eta (xkn, xalfa)
!	=============================

!
! === A function for calculating ETA for gas-aqueous
! === 	mass transfer
!
	real :: xkn, xalfa
	real :: xkn1, yyy
	
!	xkn1  = 1.0/xkn
!	yyy   = 1.0 + ( (1.33 + 0.71*xkn1)/(1.0 + xkn1)
!     &                + 4.0*(1.0 - xalfa)/(3.0*xalfa) )*xkn
!
!	cal_eta = 1.0/yyy
	
	xkn1    = 1.0 + xkn

	cal_eta = 3.0*xkn1*xalfa                             &
     	       /( xkn1*(3.0*xalfa + 4.0*xkn*(1.0 - xalfa))   &
     		+ xkn*xalfa*(3.99*xkn + 2.13) ) 

       return
        end
        
!      =================================================
       subroutine calh( i, j, k, hydro, aqc, aqr, 	     &
     			T298, pres, den )
!	=================================================

! ================================================================
!
!  CALH
!
!  Purpose:
!	Calculate [H+] in M.
!	
!!  FUNCTION zbrent(func,x1,x2,tol)
!	van Wijingaarden-Dekker-Brent roots finding package,
!	  based on Press et al., Numerical Recipes, 1992. !
!
!  Author:
!	Chien Wang
!	MIT
!
!  Revision:
!	Date	By		Brief Description
!	----	--		-----------------	
!	093098	Chien Wang	created
!	110698	Chien Wang	use q_min & n_min
!
! ================================================================

	USE gridno
	USE shared_data
	USE typedef_aq
	USE typedef_hydrometeor
		
	integer :: i, j, k
	real    :: T298, pres, den
	type (hydrometeor) :: hydro(nhydro)
	type (aq_chemical) :: aqc, aqr

	
! ==============================================================
			
#ifdef AQCHEM_ENABLE
	rke_c(2)  = 1.00e-14*exp(-6716.0*T298)
	rke_c(4)  = 4.46e-7*exp(-1000.0*T298)
	rke_c(5)  = 4.48e-11*exp(-1760.0*T298)
	rke_c(14) = 1.54e1*exp(8700.0*T298)
	rke_c(16) = 1.23e-2*exp(1960.0*T298)
	rke_c(17) = 6.61e-8*exp(1500.0*T298)
	rke_c(20) = 1.02e-2*exp(1960.0*T298)
        rke_c(23) = 1.7e-5*exp(-450.0*T298)
	
	rke_r(2)  = rke_c(2)
	rke_r(4)  = rke_c(4)
	rke_r(5)  = rke_c(5)
	rke_r(14) = rke_c(14)
	rke_r(16) = rke_c(16)
	rke_r(17) = rke_c(17)
	rke_r(20) = rke_c(20)
	rke_r(23) = rke_c(23)
		
	if (lmicro > 0) then	
	if(hydro(drop)%q.gt.qmin(drop).and.           &
     	  hydro(drop)%n.gt.xnmin(drop))then
         aqc%hplus = zbrent( 1, i, j, k, hydro, aqc,  &
     			     pres, den )
       else
         aqc%hplus = 0.0
       endif

       if(hydro(rain)%q.gt.qmin(rain).and.	      &
     	  hydro(rain)%n.gt.xnmin(rain))then
         aqr%hplus = zbrent( 2, i, j, k, hydro, aqr,  &
     			     pres, den )
	else
	  aqr%hplus = 0.0
	endif	
	endif
	
#endif
	return
	 end

!	===========================================
	real function zbrent( nnn, i, j, k, hydro, aqx,   &
     			 pres, den )
!	===========================================

	USE gridno
	USE typedef_aq
	USE typedef_hydrometeor

	! --- nnn = 1 for cloud drop
	! --- nnn = 2 for rain drop      
	
	integer :: nnn, i, j, k
	real	:: pres, den
	type (hydrometeor) :: hydro(nhydro)
	type (aq_chemical) :: aqx
	

#ifdef AQCHEM_ENABLE
	integer :: ITMAX
	real    :: tol,x1,x2,func_c,func_r,EPS

	external func_c,func_r

	parameter (ITMAX=100,EPS=3.e-8)
	parameter (x1=1.e-10,x2=5.e-1,tol=1.e-18)

	integer :: iter
	real    :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

!	----------------------------------------------

	a  = x1
	b  = x2
	
	if(nnn.eq.1)then	! cloud droplet
	  fa = func_c(a,i,j,k, hydro, aqx, pres, den)
	  fb = func_c(b,i,j,k, hydro, aqx, pres, den)
	else			! rain drop
	  fa = func_r(a,i,j,k, hydro, aqx, pres, den)
	  fb = func_r(b,i,j,k, hydro, aqx, pres, den)
	endif

!      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause
!     *'root must be bracketed for zbrent'

      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))then
	zbrent = 0.0
	return
      endif

	c  = b
	fc = fb
	do 11 iter=1,ITMAX
	  if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
		c  = a
		fc = fa
		d  = b-a
		e  = d
	  endif
	  if(abs(fc).lt.abs(fb)) then
		a  = b
		b  = c
		c  = a
		fa = fb
		fb = fc
		fc = fa
	  endif
	  tol1 = 2.*EPS*abs(b)+0.5*tol
	  xm   = 0.5*(c-b)
	  if(abs(xm).le.tol1 .or. fb.eq.0.)then
		zbrent = b
		return
	  endif
	  if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
		s = fb/fa
		if(a.eq.c) then
		  p = 2.*xm*s
		  q = 1.-s
		else
		  q = fa/fc
		  r = fb/fc
		  p = s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
		  q = (q-1.)*(r-1.)*(s-1.)
		endif
		
		if(p.gt.0.) q = -q		
		p = abs(p)
		
		if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
		  e = d
		  d = p/q
		else
		  d = xm
		  e = d
		endif
	  else
		d = xm
		e = d
	  endif
	  a  = b
	  fa = fb
	  if(abs(d) .gt. tol1) then
		b = b+d
	  else
		b = b+sign(tol1,xm)
	  endif
	  
	  if(nnn.eq.1)then
	    fb = func_c( b, i, j, k, hydro, aqx, pres, den )
	  else
	    fb = func_r( b, i, j, k, hydro, aqx, pres, den )
	  endif
	  
11	continue


        !pause 'zbrent exceeding maximum iterations'
	zbrent = b
#else
	zbrent = 0.
#endif	
      
	return
	 end

!	====================================================
	real function func_c( x, i, j, k, hydro, aqc, pres, den )
!	====================================================
!
! For return electroneutrality equation in cloud droplet
!
	USE gridno
	USE shared_data
	USE typedef_aq
	USE typedef_hydrometeor
	
	integer :: i, j, k
	real    :: x, pres, den

	type (hydrometeor) :: hydro(nhydro)
	type (aq_chemical) :: aqc
	
#ifdef AQCHEM_ENABLE
        real :: qconv
	
! =================================================================

	qconv   = 1.e3/(hydro(drop)%q*den)
	aqc%civ = aqc%civ*qconv
	aqc%siv = aqc%siv*qconv
	aqc%svi = aqc%svi*qconv
	aqc%xnv = aqc%xnv*qconv
        aqc%nh4 = aqc%nh4*qconv
		
	func_c = rke_c(2)/x                                      &
     	      + (rke_c(4)*x)					 &
     	    /(x**2 + rke_c(4)*x + rke_c(4)*rke_c(5))*aqc%civ	 &
     	      + (rke_c(16)*x + 2.0*rke_c(16)*rke_c(17)) 	 &
     	    /(x**2 + rke_c(16)*x + rke_c(16)*rke_c(17))*aqc%siv  &
     	      + (x + 2.0*rke_c(20))/(x + rke_c(20))*aqc%svi	 &
     	      + rke_c(14)/(x + rke_c(14))*aqc%xnv		 &
     	      - rke_c(23)*x/(rke_c(2) + x*rke_c(14))*aqc%nh4	 &
     	      - x 
     
	qconv   = 1.e-3*hydro(drop)%q*den
	aqc%civ = aqc%civ*qconv
	aqc%siv = aqc%siv*qconv
	aqc%svi = aqc%svi*qconv
	aqc%xnv = aqc%xnv*qconv
	aqc%nh4 = aqc%nh4*qconv
#else
	func_c = 0.0
#endif
		
	return
	 end
	 
!	====================================================
	real function func_r( x, i, j, k, hydro, aqr, pres, den )
!	====================================================
!
! For return electroneutrality equation in rain drop
!

	USE gridno
	USE shared_data
	USE typedef_gas
	USE typedef_aq
	USE typedef_hydrometeor
	
	integer :: i, j, k
	real    :: x, pres, den
	type (hydrometeor) :: hydro(nhydro)
	type (aq_chemical) :: aqr
			
#ifdef AQCHEM_ENABLE
	real :: qconv

! =================================================================

	qconv   = 1.e3/(hydro(rain)%q*den)
	aqr%civ = aqr%civ*qconv
	aqr%siv = aqr%siv*qconv
	aqr%svi = aqr%svi*qconv
	aqr%xnv = aqr%xnv*qconv
	aqr%nh4 = aqr%nh4*qconv
		
	func_r = rke_r(2)/x                                      &
     	      + (rke_r(4)*x)					 &
     	    /(x**2 + rke_r(4)*x + rke_r(4)*rke_r(5))*aqr%civ	 &
     	      + (rke_r(16)*x + 2.0*rke_r(16)*rke_r(17)) 	 &
     	    /(x**2 + rke_r(16)*x + rke_r(16)*rke_r(17))*aqr%siv  &
     	      + (x + 2.0*rke_r(20))/(x + rke_r(20))*aqr%svi	 &
     	      + rke_r(14)/(x + rke_r(14))*aqr%xnv		 &
     	      - rke_r(23)*x/(rke_r(2) + x*rke_r(14))*aqr%nh4	 &
     	      - x
     
	qconv   = 1.e-3*hydro(rain)%q*den
	aqr%civ = aqr%civ*qconv
	aqr%siv = aqr%siv*qconv
	aqr%svi = aqr%svi*qconv
	aqr%xnv = aqr%xnv*qconv
	aqr%nh4 = aqr%nh4*qconv
#else
	func_r = 0.0
#endif
		
	return
	 end		
	
