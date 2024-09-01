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
!  NUCLEATE.F:
!	subroutine nucleate
!
!  Purpose:
!	Calculate nuclation rate of aerosols (bnrate, cm-3 s-1) according 
!       to Kulmala et al. (1998)'s formula for binary nucleation.
!
!  Author:
!	Annica Ekman
!	MIT Joint Program on Science and Policy of Global Change
!       After work done by Julian Wilson, Environment Institute, Ispra, Italy.
!
! ================================================================

!	===================================================
	subroutine nucleate                                        &
#ifdef AERO_ENABLE
     		       ( temp, dens, rh, gas, gas0 )
!	===================================================

! -------------------------------------------------
!  Brief Description of Dummy Variables:
!
!	dt:	time step of integration in second
!       gas4l:  gas phase h2so4 concentration (molec cm-3)
!       g4nuc:  integral of h2so4 gas loss due to nucleation
!               over one timestep (ppbm)
!	ngain:	mass gain for nucleation mode (ppbm)
!       aerop%critn:  critical number of sulfate molecules per nucleation cluster
!       rh:     relative humidity
!       bnrate: nucleation rate (cm-3 s-1)
!
! -------------------------------------------------	

	USE gridno
	USE typedef_hydrometeor
	USE typedef_gas
	USE typedef_aq
	USE typedef_solid
	USE shared_data
	USE shared_state
	USE shared_thermo
	USE shared_aerosol_new
	
	IMPLICIT NONE

	real    :: temp0, boltzmann, f1
	real    :: peh2o, peh2so4, eta, nac, alpha, nwv, rac, xal
	real    :: temp, dens, rh, delta
	real    :: gas4t, gas4l, gas4l0
	real    :: eps_a, nuccond, masscond, massc1, massc2
	real    :: bnrate, g4nuc, ngain

        type(gas_chemical) :: gas
        type(gas_chemical)  :: gas0

        temp0     = 273.15
        boltzmann = 1.38e-16
        eps_a     = 1.e-15

        if(temp.ge.220.0.and.rh.ge.0.10) then
!        if(rh.ge.0.10) then

!       === convert gas%h2so4 to molec cm-3

#ifdef CHEM_ENABLE
        gas4l = gas%h2so4
#else
	gas4l = gas0%h2so4
#endif
	gas4l0 = gas4l
        if(gas4l.le.1.e-15) gas4l=eps_a
        call ppbm2mcm(1,gas4l,dens,1./awh2so4)

!       === Calculate equilibrium vapor pressures 
!       === (Jaeker-Mirabel, Atmos. Environ., 1988)

!       === H2O saturation vapor pressure (Tabata)

        peh2o = 0.750064*(10.**( 8.42926609-(1827.17843       &
             +71208.271/temp)/temp ))*1333./(boltzmann*temp)

!       === h2so4 saturation vapor pressure at (Ayers, GRL, 1980)

        peh2so4 = exp(-10156./temp+16.259)*7.6e2*1333./(boltzmann*temp)

!       === Relative acidity ( percentage from 0 to 1)

        rac = gas4l/peh2so4

!       === Water vapor molecule concentration (molecules cm-3)

        nwv = rh * peh2o

!
        delta = temp/temp0

!       === Molefraction of h2so4 in the critical cluster minus the h2so4g term
!       ===              (Eq. 17 in Kulmala et al.)
 
        xal = 1.2233-0.0154*rac/(rac+rh)-0.0415*log(nwv)                  &
     	   + 0.0016*temp
									 
!      === Exponent of the critical cluster (Eq. 18 in Kulmala et al.)
									 
       nac = -14.5125+0.1335*temp-(10.5462-1958.4/temp)*rh
									 
!      === Sum of all terms in Eq.20 containing h2so4g
									 
       eta = 25.1289-(4890.8-76.96166)/temp				&
     	    -(2.2479*rh + 0.02010624/rh)*delta
									 
!      === Sum all terms not containing h2so4g 
									 
       alpha = nac*(-25.1289+4890.8/temp+2.2479*delta*rh)		&
     	    - 1743.3/temp+xal*(7643.4/temp-1.9712*delta/rh) 

!       === Nucleation rate (cm-3 s-1) (eq.19)
!       === The nucleation parameterization is only valid if nuc. rate < 10e5
!       === The nucleation rate cannot be larger than the supply of h2so4

        nuccond  = log(gas4l)*eta+alpha
        if(nuccond.gt.11.5) eta = (11.5-alpha)/log(gas4l)
        massc1   = exp(alpha+log(gas4l)*eta)
        massc2   = (gas4l/(aerop%critn*dt))
        masscond = (log(gas4l/(aerop%critn*dt))-alpha)/log(gas4l)
        if(massc1.gt.massc2) then
          eta = masscond
        endif
        bnrate = exp(alpha+log(gas4l)*eta)


!       === Integral of Jnuc over timestep assuming no new h2so4g production

        if(massc1.gt.0.0) then
          f1 = gas4l**(1.0-eta)-aerop%critn*exp(alpha)*(1.0-eta)*dt
          if(f1.gt.0.0) then
            gas4t = exp(log(f1)/(1.0-eta))
            gas4t = max(gas4t, 0.0)
          else
            gas4t = gas4l
            bnrate=0.
          endif
        else
          gas4t = gas4l
          bnrate = 0.
        endif
        gas4t = min(gas4t,gas4l)

!       === Net h2so4g loss to nucleation


!       === Convert gas loss from molec cm-3 to ppb(m)
!       === Calculate change in nuclei aerosol concentration and mass
!       === Calculate change in h2so4 concentration

        g4nuc = gas4l - gas4t
        ngain = g4nuc/(aerop%critn*dt)

        call mcm2ppbm (1,gas4t,dens,awh2so4)

        g4nuc = (gas4l0 - gas4t)
        g4nuc = max(0.,g4nuc)

#ifdef CHEM_ENABLE
        gas%h2so4 = gas4t
#else
        gas0%h2so4 = gas4t
#endif

        endif     !    rh and temp conditions


#else  
        ()
#endif

	return
	 end	


! ================================================================
!
!  NUCLEATE.F:
!	subroutine bin_nuc
!
!  Purpose:
!	Calculate nuclation rate of aerosols (bnrate, cm-3 s-1) according 
!       to Vehkamaki et al. (2002)'s formula for binary nucleation.
!
!  Author:
!	Annica Ekman
!	MIT Joint Program on Science and Policy of Global Change
!       After work done by Don Lucas, MIT
!
!  Revision:
!	Date    By		Brief Description
!	----    --		-----------------	
!       121703  Annica Ekman	created
!
! ================================================================

!	===================================================
	subroutine bin_nuc                                       &
#ifdef AERO_ENABLE
     			( temp, dens, rh, gas, gas0  )
!	===================================================

! -------------------------------------------------
!  Brief Description of Dummy Variables:
!
!	dt:	time step of integration in second
!       gas4l:  gas phase h2so4 concentration (molec cm-3)
!       g4nuc:  integral of h2so4 gas loss due to nucleation
!               over one timestep (ppbm)
!	ngain:	mass gain for nucleation mode (ppbm)
!       aerop%critn:  critical number of sulfate molecules per nucleation cluster
!       rh:     relative humidity
!       bnrate: nucleation rate (cm-3 s-1)
!
! -------------------------------------------------	

	USE gridno
	USE typedef_hydrometeor
	USE typedef_gas
	USE typedef_aq
	USE typedef_solid
	USE shared_data
	USE shared_state
	USE shared_thermo
	USE shared_aerosol_new
	
	IMPLICIT NONE

        real    :: temp0, boltzmann, f1
        real    :: temp, dens, rh
        real    :: gas4t, gas4l, gas4l0
        real    :: eps_a 
        real    :: xstar, aax, bbx, ccx, ddx, eex, ffx
        real    :: ggx, hhx, iix, jjx, ntot
        real    :: ax,bx,cxx,dxx,ex,fx,gx,hx,ix,jx
	real    :: bnrate, g4nuc, ngain
	
        type(gas_chemical) :: gas
        type(gas_chemical)  :: gas0

        temp0     = 273.15
        boltzmann = 1.38e-16
        eps_a     = 1.e-15

        if(temp.ge.190.0.and.temp.le.305.15) then
        if(rh.ge.0.0001) then

!       === convert gas%h2so4 to molec cm-3

#ifdef CHEM_ENABLE
        gas4l = gas%h2so4
#else
	gas4l = gas0%h2so4
#endif
	gas4l0 = gas4l
	
        if(gas4l.le.1.e-15) gas4l=eps_a
        call ppbm2mcm(1,gas4l,dens,1./awh2so4)

!        print*,i,j,k,gas4l,gas%h2so4,dens,awh2so4
        if(gas4l.ge.1.0e4.and.gas4l.le.1.0e11) then

!   Calculate sulfuric acid mole fraction in critical cluster
        xstar = 0.740997 - 0.00266379*temp - 0.00349998*log(gas4l) +     &
     	     0.0000504022*temp*log(gas4l) + 0.00201048*log(rh) - 	 &
     	     0.000183289*temp*log(rh) + 0.00157407*log(rh)**2. -	 &
     	     0.0000179059*temp*log(rh)**2. + 0.000184403*log(rh)**3. -   &
     	     1.50345e-6*temp*log(rh)**3.
        xstar = max(xstar,eps_a)
!-------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------!
!   Check the number of molecules in the critical cluster
        aax = - 0.00295413  - 0.0976834000*temp + 0.0010248500*      &
     	     temp**2.-2.18646e-06 *temp**3. - 0.101717000/xstar       
       bbx = - 0.00205064  - 0.0075850400*temp + 0.0001926540*	     &
     	     temp**2.- 6.70430e-07*temp**3. - 0.255774000/xstar       
       ccx =   0.00322308  + 0.0008526370*temp - 0.0000154757*	     &
     	     temp**2.+ 5.66661e-08*temp**3. + 0.033844400/xstar       
       ddx =   0.0474323   - 0.0006251040*temp + 2.6506600e-6*	     &
     	     temp**2. -3.67471e-09*temp**3. - 0.000267251/xstar       
       eex = - 0.0125211   + 0.0058065500*temp - 0.0001016740*	     &
     	     temp**2. +2.88195e-07*temp**3. + 0.094224300/xstar       
       ffx = - 0.038546    - 0.0006723160*temp + 2.6028800e-6*	     &
     	     temp**2. +1.19416e-08*temp**3. - 0.008515150/xstar       
       ggx = - 0.0183749   + 0.0001720720*temp - 3.7176600e-7*	     &
     	     temp**2. -5.14875e-10*temp**3. + 0.000268660/xstar       
       hhx = - 0.0619974   + 0.0009069580*temp - 9.1172800e-7*	     &
     	     temp**2. -5.36796e-09*temp**3. - 0.007742340/xstar       
       iix = 0.0121827   - 0.0001066500*temp + 2.5346000e-7*	     &
     	     temp**2. -3.63519e-10*temp**3. + 0.000610065/xstar       
       jjx =  0.000320184 - 0.0000174762*temp + 6.0650400e-8*	     &
     	     temp**2. -1.42177e-11*temp**3. + 0.000135751/xstar
								      
       ntot = exp(aax + bbx*log(rh) + ccx*log(rh)**2. + 	     &
     	      ddx*log(rh)**3. + eex*log(gas4l) +     		     &
     	      ffx*log(rh)*log(gas4l) + ggx*log(gas4l)*log(rh)**2. +  &
     	      hhx*log(gas4l)**2. + iix*log(rh)*log(gas4l)**2. +      &
     	      jjx*log(gas4l)**3.)

        if (ntot < 4) then
          bnrate = eps_a
        else

!------------------------------------------------------------------------!

!------------------------------------------------------------------------!
!   Calculate nucleation rate
        ax =  0.14309  + 2.21956*temp - 0.0273911*temp**2.    +        &
     	     0.0000722811*temp**3. + 5.91822/xstar		        
       bx =  0.117489 + 0.462532*temp - 0.0118059*temp**2.    +        &
     	     0.0000404196*temp**3. + 15.7963/xstar		        
       cx = -0.215554 - 0.0810269*temp + 0.00143581*temp**2.   -       &
     	     4.7758e-6*temp**3. - 2.91297/xstar 		        
       dxx = -3.58856  + 0.049508*temp   - 0.00021382*temp**2.   +     &
     	     3.10801e-7*temp**3.   - 0.0293333/xstar		        
       ex =  1.14598  - 0.600796*temp	+ 0.00864245*temp**2.	-      &
     	     0.0000228947*temp**3. - 8.44985/xstar		        
       fx =  2.15855  + 0.0808121*temp  - 0.000407382*temp**2.  -      &
     	     4.01957e-7*temp**3.   + 0.721326/xstar		        
       gx =  1.6241   - 0.0160106*temp  + 0.0000377124*temp**2. +      &
     	     3.21794e-8*temp**3.   - 0.0113255/xstar		        
       hx =  9.71682  - 0.115048*temp	+ 0.000157098*temp**2.  +      &
     	     4.00914e-7*temp**3.   + 0.71186/xstar		        
       ix = -1.05611  + 0.00903378*temp - 0.0000198417*temp**2. +      &
     	     2.46048e-8*temp**3.   - 0.0579087/xstar		        
       jx = -0.148712 + 0.00283508*temp - 9.24619e-6*temp**2.	+      &
     	     5.00427e-9*temp**3. - 0.0127081/xstar
								        
       bnrate = exp(ax + bx*log(rh) + cx*log(rh)**2. + 	       &
     	     dxx*log(rh)**3. + ex*log(gas4l) + fx*log(rh)*log(gas4l) + &
     	     gx*log(gas4l)*log(rh)**2. + hx*log(gas4l)**2. + 	       &
     	     ix*log(rh)*log(gas4l)**2. + jx*log(gas4l)**3.)
        endif
!--------------------------------------------------------------------------!

        g4nuc=bnrate*ntot*dt
        g4nuc=max(0.,g4nuc)

        gas4t=gas4l - g4nuc
!        print*,i,j,k,xstar,ntot,aerop%critn
        call mcm2ppbm (1,gas4t,dens,awh2so4)
        call mcm2ppbm (1,g4nuc,dens,awh2so4)

#ifdef CHEM_ENABLE
        gas%h2so4=max(gas4t,eps_a)
#else
	gas0%h2so4=max(gas4t,eps_a)
#endif


        endif     !    sulfur conditions
        endif     !    rh and temp conditions
        endif     !    rh and temp conditions

#else
        ()
#endif

	return
	 end	

!
! ================================================================
!
!  NUCLEATE.F:
!	subroutine tern_nuc
!
!  Purpose:
!	Calculate nuclation rate of aerosols (bnrate, cm-3 s-1) according 
!       to Noppel et al. (2002)'s formula for ternary nucleation.
!       If temperature is less than 240K, or more than 305K, 
!       nucleation rates are calculated
!       according to Vehkamaki et al. (2002) binary nucleation param.
!
!  Author:
!	Annica Ekman
!	MIT Joint Program on Science and Policy of Global Change
!
!  Revision:
!	Date    By		Brief Description
!	----    --		-----------------	
!       121703  Annica Ekman	created
!
! ================================================================

!	===================================================
	subroutine tern_nuc                                       &
#ifdef AERO_ENABLE
     			( temp, dens, rh, gas, gas0 )
!	===================================================

! -------------------------------------------------
!  Brief Description of Dummy Variables:
!
!	dt:	time step of integration in second
!       gas4l:  gas phase h2so4 concentration (molec cm-3)
!       g4nuc:  integral of h2so4 gas loss due to nucleation
!               over one timestep (ppbm)
!	ngain:	mass gain for nucleation mode (ppbm)
!       aerop%critn:  critical number of sulfate molecules per nucleation cluster
!       rh:     relative humidity
!       bnrate: nucleation rate (cm-3 s-1)
!
! -------------------------------------------------	

	USE gridno
	USE typedef_hydrometeor
	USE typedef_gas
	USE typedef_aq
	USE typedef_solid
	USE shared_data
	USE shared_state
	USE shared_thermo
	USE shared_aerosol_new
	
	IMPLICIT NONE

        real    :: temp0, boltzmann, f1
        real    :: temp, dens, rh
        real    :: gas4t, gas4l, gas4l0, gasnhl, gasnht, nhtemp
        real    :: eps_a 
        real    :: xstar, aax, bbx, ccx, ddx, eex, ffx
        real    :: ggx, hhx, iix, jjx, ntot
        real    :: ncrit3, ncrit4, nhnuc, rstar
        real    :: ax,bx,cxx,dxx,ex,fx,gx,hx,ix,jx
        real    :: kx,lx,mxx,nxx,ox,px,qx,rx,sx,tx
	real    :: bnrate, g4nuc, ngain, bn

        type(gas_chemical) :: gas
        type(gas_chemical)  :: gas0

        temp0     = 273.15
        boltzmann = 1.38e-16
        eps_a     = 1.e-15

!       === convert gas%h2so4 to molec cm-3

#ifdef CHEM_ENABLE
        gas4l = gas%h2so4
#else
	gas4l = gas0%h2so4
#endif
	gas4l0 = gas4l

        if(gas4l.le.1.e-15) gas4l=eps_a
        call ppbm2mcm(1,gas4l,dens,1./awh2so4)

        if(temp.ge.240.0.and.temp.le.300.0.and.    &
     	  rh.ge.0.05.and.rh.le.0.95.and.	   &
     	  gas4l.ge.1.0e4.and.gas4l.le.1.0e9) then

!       === condition for ammonia

#ifdef CHEM_ENABLE
        gasnhl = gas%nh3*1.e3 ! ppt
#else
	gasnhl = gas0%nh3*1.e3
#endif
!         gasnhl = 10.0

! parameterization only valid 0.1ppt<nh3<100ppt
        if(gasnhl.ge.90.0) gasnhl = 90.0  

        if(gasnhl.ge.0.1.and.gasnhl.le.100) then

!       === coefficients for polynomial

        ax=-0.355297 - 33.8449*temp + 0.34536*temp**2. -                &    
     	   0.000824007*temp**3. 					 
       bx= 3.13735 - 0.772861*temp + 0.00561204*temp**2. -		&
     	   9.74576e-6*temp**3.						 
       cx= 19.0359 - 0.170957*temp + 0.000479808*temp**2. -		&
     	   4.14699e-7*temp**3.						 
       dxx= 1.07605 + 1.48932*temp - 0.00796052*temp**2. +		&
     	    7.61229e-6*temp**3. 					 
       ex= 6.0916 - 1.25378*temp + 0.00939836*temp**2. -		&
     	   0.0000174927*temp**3.					 
       fx= 0.31176 + 1.64009*temp - 0.00343852*temp**2. -		&
     	   0.0000109753*temp**3.					 
       gx=-0.0200738 - 0.752115*temp + 0.00525814*temp**2. -		&
     	   8.98038e-6*temp**3.						 
       hx= 0.165536 + 3.26623*temp - 0.0489703*temp**2. +		&
     	   0.000146967*temp**3. 					 
       ix= 6.52645 - 0.258002*temp + 0.00143456*temp**2. -		&
     	   2.02036e-6*temp**3.						 
       jx= 3.68024 - 0.204098*temp + 0.00106259*temp**2. -		&
     	   1.2656e-6*temp**3.						 
       kx=-0.066514 - 7.82382*temp + 0.0122938*temp**2. +		&
     	   0.0000618554*temp**3.					 
       lx= 0.65874 + 0.190542*temp - 0.00165718*temp**2. +		&
     	   3.41744e-6*temp**3.						 
       mxx= 0.0599321 + 5.96475*temp - 0.0362432*temp**2. +		&
     	   0.0000493337*temp**3.					 
       nxx=-0.732731 - 0.0184179*temp + 0.000147186*temp**2. -		&
     	   2.37711e-7*temp**3.						 
       ox= 0.728429 + 3.64736*temp - 0.027422*temp**2. +		&
     	   0.0000493478*temp**3.					 
       px= 41.3016 - 0.35752*temp + 0.000904383*temp**2. -		&
     	   5.73788e-7*temp**3.						 
       qx=-0.160336 + 0.00889881*temp - 0.0000539514*temp**2. + 	&
     	   8.39522e-8*temp**3.						 
       rx= 8.57868 - 0.112358*temp + 0.000472626*temp**2. -		&
     	   6.48365e-7*temp**3.						 
       sx= 0.0530167 - 1.98815*temp + 0.0157827*temp**2. -		&
     	   0.0000293564*temp**3.					 
       tx=-2.32736 + 0.0234646*temp - 0.000076519*temp**2. +		&
     	   8.0459e-8*temp**3.						 
        								 
!     === nucleation rate
 									 
       bnrate = -84.7551+ax/log(gas4l)+bx*log(gas4l)+		&
     			cx*(log(gas4l))**2.+dxx*log(gasnhl) +		&
     			ex*(log(gasnhl))**2.+fx*rh+gx*log(rh) + 	&
     			hx*log(gasnhl)/log(gas4l)+ix*log(gasnhl)*	&
     			log(gas4l)+jx*rh*log(gas4l)+kx*rh/		&
     			log(gas4l)+lx*rh*log(gasnhl)+mxx*log(rh)/	&
     			log(gas4l)+nxx*log(rh)*log(gasnhl)+		&
     			ox*(log(gasnhl))**2./log(gas4l)+px*		&
     			log(gas4l)*(log(gasnhl))**2.+qx*		&
     			(log(gas4l))**2.*log(gasnhl)+rx*rh*		&
     			(log(gasnhl))**2.+sx*rh*			&
     			log(gasnhl)/log(gas4l)+tx*(log(gas4l))**	&
     			2.*(log(gasnhl))**2.


       bnrate = exp(bnrate)
       bn=bnrate

!      === number of particles of each compound in crit. nuc.

        ncrit4 = 38.1645 + 0.774106*log(bn) + 0.00298879*        &
     		log(bn)**2. - 0.357605*temp - 0.00366358*	 &
     		temp * log(bn)+0.0008553*temp**2.
								  
       ncrit3 = 26.8982 + 0.682905*log(bnrate) + 	 &
     		0.00357521*log(bn)**2.-0.265748*temp -		 &
     		0.00341895*temp*log(bn) + 0.000673454*		 &
     		temp**2.
								  
       ntot = 79.3484 + 1.7384*log(bn) + 0.00711403 *		 &
     	      log(bn)**2. - 0.744993*temp - 0.00820608 * 	 &
     	      temp*log(bn) + 0.0017855*temp**2.
								  
       rstar = 0.141027 - 0.00122625*log(bn) - 7.82211e-6*	 &
     	       log(bn)**2. - 0.00156727*temp - 0.000003076*	 &
     	       temp*log(bn) + 0.0000108375*temp**2.


        g4nuc=bnrate*ncrit4*dt
        g4nuc=max(0.,g4nuc)

        nhtemp=gasnhl*1.e-3 !from ppt to ppb
        call ppbm2mcm(1,nhtemp,dens,1./awnh3)

        nhnuc=bnrate*ncrit3*dt

        nhtemp=nhtemp-nhnuc
        call mcm2ppbm(1,nhtemp,dens,awnh3)
        gasnht=nhtemp
        
        gas4t=gas4l-g4nuc
        call mcm2ppbm (1,gas4t,dens,awh2so4)

        call mcm2ppbm(1,g4nuc,dens,awh2so4)
        call mcm2ppbm(1,nhnuc,dens,awnh3)

! --- At the moment we assume that the nucleated aerosol is sulfur....

! it's the next row......
#ifdef CHEM_ENABLE
        gas%h2so4=max(gas4t,eps_a)
        gas%nh3=max(gas%nh3-nhnuc,eps_a)
#else
        gas0%h2so4=max(gas4t,eps_a)
        gas0%nh3=max(gas0%nh3-nhnuc,eps_a)
#endif
        endif 
!
        elseif(temp.ge.190.0.and.temp.le.305.15.and.   &
     	      rh.ge.0.0001.and. 		       &
     	      gas4l.ge.1.0e4.and.gas4l.le.1.0e11) then   ! cond. for ternary nuc.

!   Calculate sulfuric acid mole fraction in critical cluster 
        xstar = 0.740997 - 0.00266379*temp - 0.00349998*log(gas4l) +     &
     	     0.0000504022*temp*log(gas4l) + 0.00201048*log(rh) - 	 &
     	     0.000183289*temp*log(rh) + 0.00157407*log(rh)**2. -	 &
     	     0.0000179059*temp*log(rh)**2. + 0.000184403*log(rh)**3. -   &
     	     1.50345e-6*temp*log(rh)**3.
        xstar = max(xstar,eps_a)
!-------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------!
!   Check the number of molecules in the critical cluster
        aax = - 0.00295413  - 0.0976834000*temp + 0.0010248500*       &
     	     temp**2.-2.18646e-06 *temp**3. - 0.101717000/xstar        
       bbx = - 0.00205064  - 0.0075850400*temp + 0.0001926540*	      &
     	     temp**2.- 6.70430e-07*temp**3. - 0.255774000/xstar        
       ccx =   0.00322308  + 0.0008526370*temp - 0.0000154757*	      &
     	     temp**2.+ 5.66661e-08*temp**3. + 0.033844400/xstar        
       ddx =   0.0474323   - 0.0006251040*temp + 2.6506600e-6*	      &
     	     temp**2. -3.67471e-09*temp**3. - 0.000267251/xstar        
       eex = - 0.0125211   + 0.0058065500*temp - 0.0001016740*	      &
     	     temp**2. +2.88195e-07*temp**3. + 0.094224300/xstar        
       ffx = - 0.038546    - 0.0006723160*temp + 2.6028800e-6*	      &
     	     temp**2. +1.19416e-08*temp**3. - 0.008515150/xstar        
       ggx = - 0.0183749   + 0.0001720720*temp - 3.7176600e-7*	      &
     	     temp**2. -5.14875e-10*temp**3. + 0.000268660/xstar        
       hhx = - 0.0619974   + 0.0009069580*temp - 9.1172800e-7*	      &
     	     temp**2. -5.36796e-09*temp**3. - 0.007742340/xstar        
       iix = 0.0121827   - 0.0001066500*temp + 2.5346000e-7*	      &
     	     temp**2. -3.63519e-10*temp**3. + 0.000610065/xstar        
       jjx =  0.000320184 - 0.0000174762*temp + 6.0650400e-8*	      &
     	     temp**2. -1.42177e-11*temp**3. + 0.000135751/xstar
								       
       ntot = exp(aax + bbx*log(rh) + ccx*log(rh)**2. + 	      &
     	      ddx*log(rh)**3. + eex*log(gas4l) +     		      &
     	      ffx*log(rh)*log(gas4l) + ggx*log(gas4l)*log(rh)**2. +   &
     	      hhx*log(gas4l)**2. + iix*log(rh)*log(gas4l)**2. +       &
     	      jjx*log(gas4l)**3.)

        if (ntot < 4) then
          bnrate = eps_a
        else

!        aerop%critn=xstar*ntot
!------------------------------------------------------------------------!

!------------------------------------------------------------------------!
!   Calculate nucleation rate
        ax =  0.14309  + 2.21956*temp - 0.0273911*temp**2.    +         &
     	     0.0000722811*temp**3. + 5.91822/xstar 			 
       bx =  0.117489 + 0.462532*temp - 0.0118059*temp**2.    + 	&
     	     0.0000404196*temp**3. + 15.7963/xstar			 
       cx = -0.215554 - 0.0810269*temp + 0.00143581*temp**2.   - 	&
     	     4.7758e-6*temp**3. - 2.91297/xstar 			 
       dxx = -3.58856  + 0.049508*temp   - 0.00021382*temp**2.   + 	&
     	     3.10801e-7*temp**3.   - 0.0293333/xstar			 
       ex =  1.14598  - 0.600796*temp	+ 0.00864245*temp**2.	- 	&
     	     0.0000228947*temp**3. - 8.44985/xstar			 
       fx =  2.15855  + 0.0808121*temp  - 0.000407382*temp**2.  - 	&
     	     4.01957e-7*temp**3.   + 0.721326/xstar			 
       gx =  1.6241   - 0.0160106*temp  + 0.0000377124*temp**2. + 	&
     	     3.21794e-8*temp**3.   - 0.0113255/xstar			 
       hx =  9.71682  - 0.115048*temp	+ 0.000157098*temp**2.  + 	&
     	     4.00914e-7*temp**3.   + 0.71186/xstar			 
       ix = -1.05611  + 0.00903378*temp - 0.0000198417*temp**2. + 	&
     	     2.46048e-8*temp**3.   - 0.0579087/xstar			 
       jx = -0.148712 + 0.00283508*temp - 9.24619e-6*temp**2.	+ 	&
     	     5.00427e-9*temp**3. - 0.0127081/xstar
									 
       bnrate = exp(ax + bx*log(rh) + cx*log(rh)**2. + 		&
     	     dxx*log(rh)**3. + ex*log(gas4l) + fx*log(rh)*log(gas4l) +  &
     	     gx*log(gas4l)*log(rh)**2. + hx*log(gas4l)**2. + 		&
     	     ix*log(rh)*log(gas4l)**2. + jx*log(gas4l)**3.)
        endif
!--------------------------------------------------------------------------!

        g4nuc=bnrate*ntot*dt
        g4nuc=max(0.,g4nuc)

        gas4t=gas4l-g4nuc

        call mcm2ppbm (1,gas4t,dens,awh2so4)
        call mcm2ppbm(1,g4nuc,dens,awh2so4)

#ifdef CHEM_ENABLE
        gas%h2so4=max(gas4t,eps_a)
#else
        gas0%h2so4=max(gas4t,eps_a)
#endif
        endif     !    cond. for ternary nuc.

#else
        ()
#endif

	return
	 end	
      
!	========================================
	subroutine ppbm2mcm (np, xx, den, revaw) 
!	========================================

! ------------------------------------------------
! --- A program for convert mixing ratio in
! --- 	ppb(m) to concentration in molecules/cm^3
! ------------------------------------------------

	integer :: np, i
	real    :: xx(np), revaw(np), den, ddd 

	ddd = 6.02217e+11*den
	do i=1,np
	  xx (i) = xx (i)*ddd*revaw(i)
	end do

	return
	 end

!	=====================================
	subroutine mcm2ppbm (np, xx, den, aw)
!	=====================================

! -----------------------------------------------
! --- A program for convert concentration in
! --- 	molecules/cm^3 to mixing ratio in ppb(m)
! -----------------------------------------------

	integer :: np, i
	real    :: xx(np),aw(np),den,ddd

	ddd = 1./(6.02217e+11*den)
	do i=1,np
	  xx (i) = xx (i)*ddd*aw(i)
	end do

	return
	 end
       

