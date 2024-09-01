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
!  AEROSOL.F:
!      subroutine aerosol_init,
!      subroutine aerosol_init_0,
!       subroutine aerosol_calc
!
!  Purpose:
!      For handling aerosol physics and chemistry in CM5 model.
!
!  Author:
!      Chien Wang
!      MIT Joint Program on the Science and Policy of Global Change
!
! =====================================================
!
  module aerosols
!
  USE gridno
  USE shared_data
  USE shared_all
  USE aeroactivate  
  USE thermodynamics

  IMPLICIT NONE

  private
  
  public :: init_aerosol, fill_aerosol, regeneration, aerosol_calc

  contains

! =====================================================
  subroutine init_aerosol
! =====================================================

#ifdef AERO_ENABLE
      real, dimension(nz) :: nuc_profile, ait_profile
      real, dimension(nz) :: acc_profile, coa_profile
        
      character(len = 200):: file_acc,file_nuc,file_ait, file_coa
      character(len = 200):: filename

      namelist /cm_aero/file_acc,file_ait, file_coa,file_nuc
#endif
!
      integer :: i,j,k,l,h,im
!
!----------------------------------------------------------!
!         Bulk CCN and IN initialized in any case 	   !
!----------------------------------------------------------!
!
  nl0 = 0.
  ccn0 = xn_ccn0 / den0
  in0  = xn_in0 / den0 
!
!----------------------------------------------------------!
!		     Initialize aerosols      		   !
!----------------------------------------------------------!
!
!  Initial aerosol properties
!
#ifdef AERO_ENABLE
  allocate(aero(1:nmode), aero0(1:nmode))
!
  aero(1:nmode0)%nelem = aeroi(1:nmode0)%nelem
  aero(1:nmode0)%init  = aeroi(1:nmode0)%init
  aero(1:nmode0)%size  = aeroi(1:nmode0)%size
!
  if (reg_mode) then
    aero(nmode)%nelem = aeror%nelem
    aero(nmode)%init  = aeror%init
    aero(nmode)%size  = aeror%size
  endif
!    
  do im = 1, nmode
    call initialize_mix (im)
  enddo
  !
!Mean mass of aerosol size distribution [kg]: 
  aero(1:nmode)%size%mmean = 4./3.*pi*aero(1:nmode)%mix%rho *	&
  	(1.e-6*aero(1:nmode)%size%rmean)**3. * exp(4.5*log(aero(1:nmode)%size%sigma)*log(aero(1:nmode)%size%sigma))
!
  aero%depv0 = 0.1
!  
!  allocate (tab_coag(1:nmode,1:nmode,1:nmode))
!  call coag_init 
!
!  Initialize aerosols

if (casename=='SE_ATLANTIC') then
   
    call forcing_SE_Atlantic_aerosols
  
  else
!
  do im = 1, nmode
    do k=1,nz
      aero0(im)%n(k)  = aero(im)%init%n0/ den0(k)                          !Number concentration (1D) [#/kg]
      aero0(im)%m(k)  = aero(im)%size%mmean * aero(im)%init%n0 /den0(k)    !Total mass (1D) [kg/kg] = Mean mass [kg] * Number concentration [#/kg]
      aero0(im)%ma(k) = 0.                                                 !Mass of activated aerosols (1D)
    end do
  end do
endif

#endif
!
!  Cloud droplet number concentration
!
  do k = 1, nz
    if ( lndrop == 1 .and. ql0(k) > qmin(drop) ) then
#ifdef AERO_ENABLE
      nl0(k) = sum( aero0%n(k) )
#else
      nl0(k) = ccn0(k) 
#endif
    else if ( lndrop == 0 ) then
#ifdef AERO_ENABLE
      nl0(k) = sum( aero0%n(k) )
#else
      nl0(k) = ccn0(k) 
#endif
    endif
  enddo
!
!----------------------------------------------------------!
!                  Initialize ice nuclei                   !
!----------------------------------------------------------!
!
#ifdef NUC_CNT
        call cnt_init
!
	do k = 1, nz
	  do h = 1, 3
!
!  Total aero populations
!
	    nucin0(h)%mode(1)%n(k) = nucp(h)%n * nucprop(h)%frac
	    nucin0(h)%mode(1)%nc(k) = nucp(h)%nc * nucprop(h)%fracc 
#ifdef NUC_CNT1
	    nucin0(h)%mode(1)%m(k) = pi/6.*nucp(h)%rhop * cal_lnmom(3.,nucp(h)%s,nucp(h)%d) * nucin0(h)%mode(1)%n(k)
	    nucin0(h)%mode(1)%mc(k) = pi/6.*nucp(h)%rhop * cal_lnmom(3.,nucp(h)%sc,nucp(h)%dc) * nucin0(h)%mode(1)%nc(k)
#endif
!
!  Immersed aerosols
!
	    if ( lndrop == 1 ) then
  	      nucin0(h)%mode(2)%n(k) = nucin0(h)%mode(1)%n(k) * nl0(k) / ccn0(k)
	      nucin0(h)%mode(2)%nc(k) = nucin0(h)%mode(1)%nc(k) * nl0(k) / ccn0(k)
#ifdef NUC_CNT1
	      nucin0(h)%mode(2)%m(k) = nucin0(h)%mode(1)%m(k) * nl0(k) / ccn0(k)
	      nucin0(h)%mode(2)%mc(k) = nucin0(h)%mode(1)%mc(k) * nl0(k) / ccn0(k)
#endif
	    else
  	      nucin0(h)%mode(2)%n(k) = 0.
	      nucin0(h)%mode(2)%nc(k) = 0.
#ifdef NUC_CNT1
	      nucin0(h)%mode(2)%m(k) = 0.
	      nucin0(h)%mode(2)%mc(k) = 0.    
#endif
	    endif
!
!  Frozen aero populations (0.)
!
	    nucin0(h)%mode(3)%n(k) = 0.
	    nucin0(h)%mode(3)%nc(k) = 0.
#ifdef NUC_CNT1
	    nucin0(h)%mode(3)%m(k) = 0.
	    nucin0(h)%mode(3)%mc(k) = 0.
#endif
	  enddo
	enddo
#endif
!
!----------------------------------------------------!
!
  return
  end
!
!  ===================================================
  subroutine fill_aerosol 
!  ===================================================
!
      integer :: i,j,k,l,im
!
!----------------------------------------------------------!
!         Bulk CCN and IN initialized in any case 	   !
!----------------------------------------------------------!
!
#ifdef SEIFERT
  do k = 1, nz
    do j = jp_start, jp_end
      do i = ip_start, ip_end
  	nuc2%ccn(i,j,k) = ccn0(k)
        nuc2%in(i,j,k)  = in0(k)
      enddo
    enddo   
  end do
#endif
!
!----------------------------------------------------------!
!		     Initialize aerosols      		   !
!----------------------------------------------------------!
!
#ifdef AERO_ENABLE
  do k = 1, nz
    do im = 1, nmode
      if ( lndrop == 1 .and. aero_flg%any .and. nl0(k) > xnmin(drop) ) then
        aero3d2(im)%n(:,:,k) = 0.
        aero3d2(im)%m(:,:,k) = 0.
        aero3d2(im)%ma(:,:,k) = aero0(im)%m(k) 
      else
        aero3d2(im)%n(:,:,k)  = aero0(im)%n(k)    !Number concentration (3D) [#/kg]  
        aero3d2(im)%m(:,:,k)  = aero0(im)%m(k)    !Total mass (3D) [kg/kg] = Mean mass [kg] * Number concentration [#/kg]
        aero3d2(im)%ma(:,:,k) = 0.                !Mass of activated aerosols (3D)
      endif
    end do
  end do
#endif
!
!  Cloud droplet number concentration
!
  do k = 1, nz
    hydromtr2(drop)%n(:,:,k) = nl0(k)
  enddo
!
!----------------------------------------------------------!
!                  Initialize ice nuclei                   !
!----------------------------------------------------------!
!
#ifdef NUC_CNT
	do h = 1, 3
	  do l = 1, 3
	    do k = 1, nz
  	      nucin2(h)%mode(l)%n(:,:,k) = nucin0(h)%mode(l)%n(k)
  	      nucin2(h)%mode(l)%nc(:,:,k) = nucin0(h)%mode(l)%nc(k)
#ifdef NUC_CNT1
  	      nucin2(h)%mode(l)%m(:,:,k) = nucin0(h)%mode(l)%m(k)
  	      nucin2(h)%mode(l)%mc(:,:,k) = nucin0(h)%mode(l)%mc(k)
#endif
	    enddo
	  enddo
	enddo
#endif
!
!----------------------------------------------------!
!
  return
  end
!
!  ===================================================
  subroutine regeneration (i, j, k, dref, nall, aerol, dnall, aero3dm )
!  ===================================================
!
  type (hydrometeor), dimension(1:nhydro) :: hydromtrl
  type (aero_3d), dimension(1:nmode) :: aerol
  type (aero_3d), dimension(1:nmode), optional :: aero3dm
!
  integer :: i, j, k, im
  real    :: renc, nall, dnall, d0, dref
  real, dimension(1:nmode) :: naer, dn, dm
!
!----------------------------------------------------------!
!            Calculate aerosol release due to              !
!                 evaporation/sublimation        	   !
!----------------------------------------------------------!
!
#ifdef CONSERVATIVE
  d0 = dref
#else
  d0 = 1.
#endif
!
!  Mean mass and number
!
  renc = abs( min(dnall,0.) )
  naer = aerol%ma(i,j,k)/aero%size%mmean
!
!  Get tendencies
!
  if ( renc > 0. .and. maxval(naer) > 1.e-6 ) then
    naer = naer / sum(naer)
    do im = 1, nmode
      dm(im) = min(aerol(im)%ma(i,j,k)*renc/nall, aerol(im)%ma(i,j,k)/dt0)
      dn(im) = naer(im)*renc
    enddo
!
    do im = 1, nmode
      if ( present(aero3dm) ) then
        if (reg_mode) then
          aero3dm(nmode)%m(i,j,k) = aero3dm(nmode)%m(i,j,k) + dm(im) 
          aero3dm(nmode)%n(i,j,k) = aero3dm(nmode)%n(i,j,k) + dn(im)
        else
          aero3dm(im)%m(i,j,k) = aero3dm(im)%m(i,j,k) + dm(im)
          aero3dm(im)%n(i,j,k) = aero3dm(im)%n(i,j,k) + dn(im)
        endif
        aero3dm(im)%ma(i,j,k) = aero3dm(im)%ma(i,j,k) - dm(im) 
      else
        if (reg_mode) then
          aerol(nmode)%m(i,j,k) = aerol(nmode)%m(i,j,k) + dm(im)*dt0 
          aerol(nmode)%n(i,j,k) = aerol(nmode)%n(i,j,k) + dn(im)*dt0
        else
          aerol(im)%m(i,j,k) = aerol(im)%m(i,j,k) + dm(im)*dt0
          aerol(im)%n(i,j,k) = aerol(im)%n(i,j,k) + dn(im)*dt0
        endif
        aerol(im)%ma(i,j,k) = aerol(im)%ma(i,j,k) - dm(im)*dt0
      endif
!
      if (ldiag .and. out_diaga) then
        if (reg_mode) then    
          diag(6)%aero(nmode)%n(i,j,k) = diag(6)%aero(nmode)%n(i,j,k) + d0*cint*dn(im)
          diag(9)%aero(nmode)%n(i,j,k) = diag(9)%aero(nmode)%n(i,j,k) + d0*cint*dn(im)
        else
          diag(6)%aero(im)%n(i,j,k) = diag(6)%aero(im)%n(i,j,k) + d0*cint*dn(im)
          diag(9)%aero(im)%n(i,j,k) = diag(9)%aero(im)%n(i,j,k) + d0*cint*dn(im)
        endif
        diag(6)%aero(im)%m(i,j,k) = diag(6)%aero(im)%m(i,j,k) - d0*cint*dm(im)
        diag(9)%aero(im)%m(i,j,k) = diag(9)%aero(im)%m(i,j,k) - d0*cint*dm(im)
      endif
!
    enddo
  endif
!
  return
  end
!
! ===========================================
  subroutine aerosol_calc ( aerol )
! ==========================================

! ------------------------------------------------------------------
!  Brief description of dummy variables:
!
!      ttn:      aerosol mass for aerosol mode n=1,3 (molecules cm-3)
!      c2:      condensation coefficient for mode n=1,3
!      rr:      equillibrium radius (wet ambient) for mode n=1,3
!      com:      coagulation coefficient for mode n=1,3
!       tri:      index for temperature
!       rhi:      index for RH
!      pri      index for pressure
!       gas4l      h2so4 concentration (molecules cm-3) 
!       sulfem: direct sulfate emission source
!      foscarbem: fossil fuel carbon emissions
!      biocarbem: biomass burning carbon emissions
!       dia:       diameter of average mass
!       nmd:      number median diameter
!--------------------------------------------------------------------
      
#include "chemdef.h"

       real, parameter :: bk   = 1.38e-16 
       real, parameter :: pdfj = 0.3
                                                                                          
       integer :: i, j, k, im, im1, iel
       logical :: present, presentbc, presentsu    
       real    :: tr, tkmax, rhr, pr, tk         
       real    :: gas4l, sumcond, sulfem    
       real    :: lair, Kn, alpha, awall, xxx, delv, dmax, movma, movnr, maxpara      
       real    :: dens,xfmu,xfkd,xib,velb,relhum,Tem,ppp,qv,ak
       
       type (aero_3d), dimension(1:nmode)  :: aerol
       real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: temperature
                    
       real, dimension (nmode,nmode) :: coag, com 
       real, dimension (nmode)       :: conf, rad, rho, c2,     &
                                        rr, raa, condterm,    	&
                                        densae, denswe, nmd,    &
					fcv2x, fpd2x, 	   	&
					fh2x, fvs2x

       real, dimension (2,nmode,nelem) :: epsl
       real, dimension (2,nmode)       :: ttn, num
       real, dimension(nmode,15)       :: diag

! ================================================================================
!
        call get_temperature( pressure, state, hydromtr, temperature )
!
        do k = 1, nz
          do j = jt_start,jt_end
            do i = it_start,it_end
!
!----------------------------------------------------------!
!                     Initializations                      !
!----------------------------------------------------------!
!
	     diag = 0.0
     	     condterm = 0.0 
             conf     = 0.0
	     rr       = 0.0
	     fvs2x    = 0.0
	     fcv2x    = 0.0
	     fpd2x    = 0.0
	     fh2x     = 0.0
!
             ppp  = pressure%p(i,j,k)+p0(k)
             dens = pressure%dens(i,j,k)
	     Tem  = temperature(i,j,k)
	     ak   = 2.*cal_sigmawv(Tem) / (Rw*Tem*rhow) 
!
	     qv  = state%qt(i,j,k)
	     if (lmicro > 0) qv = qv - hydromtr(drop)%q(i,j,k) - hydromtr(rain)%q(i,j,k)
	     if (lmicro > 1) qv = qv - hydromtr(ice)%q(i,j,k)
	     if (lmicro > 2) qv = qv - hydromtr(grau)%q(i,j,k) - hydromtr(snow)%q(i,j,k)
	     if (lmicro > 3) qv = qv - hydromtr(hail)%q(i,j,k)
             relhum = qv / cal_qsw (Tem, ppp)
!
	     xfkd = cal_xfkd (Tem,ppp)
             xfmu = cal_xfmu (Tem)
             lair = 21.55e-6*xfmu*sqrt(Tem)/ppp
!
!  Find BC & aged BC (if they exist)
!	
	     bc  = 0
	     bcs = 0
	     presentbc = .false.
	     presentsu = .false.
	     do im = 1, nmode
	       call have_element (im, 'BC  ', presentbc, im1)
	       if (presentbc .and. aero(im)%nelem == 1) bc = im
	       call have_element (im, 'SULF', presentsu, im1)
	       if (presentbc .and. presentsu .and. aero(im)%nelem == 2) bcs = im
	     enddo
!
!  Find sulfate modes (if they exist)
!	
	     anu = 0
	     ait = 0
	     acc = 0
	     presentsu = .false.
	     do im = 1, nmode
	       call have_element (im, 'SULF', presentsu, im1)
	       if (presentsu .and. aero(im)%nelem == 1) then
	         if (aero(im)%size%rmean <= 0.02) anu = im
	         if (aero(im)%size%rmean > 0.02 .and. aero(im)%size%rmean <= 0.2) ait = im
	         if (aero(im)%size%rmean > 0.2  .and. aero(im)%size%rmean <= 2.)  acc = im
	       endif
	     enddo
!
!----------------------------------------------------------!
!               Initialize temporary arrays                !
!----------------------------------------------------------!
!
	    do im = 1, nmode
	      epsl(1:2,im,:) = 0.0
              ttn(1,im) = aerol(im)%m(i,j,k)
              num(1,im) = aerol(im)%n(i,j,k)
	      do iel = 1, aero(im)%nelem
	        epsl(1,im,iel) = aero(im)%init%frac(iel)
	      enddo
!
	      ttn(2,im) = ttn(1,im)
	      num(2,im) = num(1,im)    
	      epsl(2,im,:) = epsl(1,im,:)
	    enddo
!
!  Calculate properties of mixed aerosols
!
!            call calc_mix_prop (epsl(1,:,:))	! Needed to comment the line as calc_mix_prop no longer exists
!
!----------------------------------------------------------!
!         Condensation and coagulation parameters          !
!----------------------------------------------------------!
!
	  do im = 1, nmode
	    if (ttn(1,im) > 1.e-12 .and. num(1,im) > 1.e-6) then
	      raa (im) = cal_rgeo (aero(im)%size%sigma, aero(im)%mix%rho, ttn(1,im), num(1,im))
	      !call activate ( relhum, 0., raa(im), aero(im)%size%sigma, aero(im)%mix%rho, ak, aero(im)%mix%kappa, rwet=rr(im) )
!
              nmd (im) = 2. * raa (im)
	      densae(im) = aero(im)%mix%rho
	      denswe(im) = rhow + (densae(im) - rhow)*(raa(im)/rr(im))**3.
!  
  	      awall = ((rr(im)**3. - raa(im)**3.)*awh2o + raa(im)**3.*aero(im)%mix%mw) / rr(im)**3.      
	      Kn   = 0.5*lair / rr(im)
              alpha = 2.514 + 0.8*exp(-0.55/Kn)
              fpd2x (im) =  dens*kb*Tem*(1. + alpha*Kn) / (6.*pi*xfmu*rr(im))
!
	      fvs2x (im) = fpd2x(im)/(kb*Tem) * 4./3.*pi*rr(im)**3. * g * (denswe(im) - dens)
!
              velb   = sqrt(8.*kb*Tem / (pi*awH2SO4*1.e-3))
              xib    = 4.*xfkd / (alpha*velb)
              delv   = ( (rr(im) + xib)**3 - (rr(im)**2 + xib**2)**1.5 ) 			     &
	             / (3. * rr(im) * xib) - rr(im)		       
!	      
              fcv2x (im) = sqrt(8.*kb*Tem / (pi*awall*1.e-3))			       
              xxx        = 4.*fpd2x(im) / (alpha*fcv2x(im))
              fh2x (im)  = ( ( (2.*rr(im) + 2.*xxx)**3. - (4.*rr(im)**2. + 4.*xxx**2.)**1.5 )	     &
             		 / ( 12.*rr(im) * xxx) - 2.*rr(im) ) * sqrt(2.0)                                                      
!
              conf(im) = 4.*pi*rr(im)**2. * xfkd / (xib/rr(im) + 1.0/(1.0 + delv / rr(im)))
	    endif
          end do
!
!----------------------------------------------------------!
!             	      Condensation                         !
!----------------------------------------------------------!   
!
#ifdef CHEM_ENABLE
	    gas4l = gas2(i,j,k)%h2so4
#else
            gas4l = gas0(k)%h2so4
#endif
            call ppbm2mcm(1,gas4l,dens,1./awh2so4)
!
!  Deposition of h2so4 on sulfate containing aerosols
!
#ifdef AERO_COND
	    call condensation (gas4l, conf, condterm, ttn, eps, diag)
#endif 
!
!----------------------------------------------------------!
!           		Coagulation                        !
!----------------------------------------------------------!
!
#ifdef AERO_COAG
!
!  Coagulation kernels
!
	    call coag_kernel ( fcv2x, fpd2x, fh2x, fvs2x, rr, coag )
!
!  Calculate coagulation: particles collecting H2SO4
!
            call coagulation ( nmd, ttn, num, eps, coag, diag )
!
#endif
!
!----------------------------------------------------------!
!           	  	  BC aging                         !
!----------------------------------------------------------!
!
#ifdef AERO_AGE
!
            call bc_ageing (condterm, nmd, ttn, num, eps, diag)
!
#endif
!
!----------------------------------------------------------!
!           		Aero check                         !
!----------------------------------------------------------!
!
#ifdef AERO_CHECK                                                                        
!
!  Move from nucleation to aitken mode
!                                      
             if (anu /= 0) then    
	       dmax    = 2.*cal_inclogn (1.e-6*aero(anu)%size%rmean, aero(anu)%size%sigma, 0.99) 	     
               movnr   = num(2,anu)*0.05  			 
               maxpara = dmax * exp(1.5*(log(aero(anu)%size%sigma))**2)		    
               maxpara = max(maxpara,nmd(anu))  			 
               movma   = movnr * (aero(anu)%mix%rho*pi/6.*maxpara**3)		     
!            						
               ttn(2,anu) = ttn(2,anu) - movma			   
               ttn(2,ait) = ttn(2,ait) + movma			   
               num(2,anu) = num(2,anu) - movnr			   
               num(2,ait) = num(2,ait) + movnr			   
             endif					       
!
!  Move from aitken to accumulation mode
!                                                         
             if (ait /= 0) then			   
	       dmax    = 2.*cal_inclogn (1.e-6*aero(ait)%size%rmean, aero(ait)%size%sigma, 0.99)   
               movnr   = num(2,ait)*0.05  			 
               maxpara = dmax * exp(1.5*(log(aero(ait)%size%sigma))**2)		   
               maxpara = max(maxpara,nmd(ait))  			 
               movma   = movnr * (aero(ait)%mix%rho*pi/6.*maxpara**3)   

               ttn(2,ait) = ttn(2,ait) - movma
               ttn(2,acc) = ttn(2,acc) + movma
               num(2,ait) = num(2,ait) - movnr
               num(2,acc) = num(2,acc) + movnr
             endif
#endif
!
!----------------------------------------------------------!
!                     Dry deposition                       !
!----------------------------------------------------------!
!      
#ifdef AERO_DDEP
             if(k == 1) call drydep (aerol, i, j)
#endif
!
!----------------------------------------------------------!
!           		     End                           !
!----------------------------------------------------------!
!
	     do im = 1, nmode
               aero3ds(im)%n(i,j,k) = aero3ds(im)%n(i,j,k) + 		&
               		(max(0.0,num(2,im)) - aerol(im)%n(i,j,k))/dt
               aero3ds(im)%m(i,j,k) = aero3ds(im)%m(i,j,k) + 		&
               		(max(0.0,ttn(2,im)) - aerol(im)%m(i,j,k))/dt
	     enddo
!
             call mcm2ppbm (1,gas4l,dens,awh2so4)

#ifdef CHEM_ENABLE
	     gas2(i,j,k)%h2so4 = gas4l
#else
             gas0(k)%h2so4 = gas4l
#endif
!
!----------------------------------------------------------!
!               Nucleation of new aerosols                 !
!----------------------------------------------------------!
!      
#ifdef AERO_NUC
!
#ifdef AERO_NUCBIN
        call bin_nuc (Tem, ppp, dens, relhum, 	&
	              aero3d, gas1(i,j,k), gas0(k) )
#elif (defined AERO_NUCTERN)
        call tern_nuc(Tem, ppp, dens, relhum, 	&
	              aero3d, gas1(i,j,k), gas0(k) )
#else
        call nucleate(Tem, ppp, dens, relhum, 	&
	              aero3d, gas1(i,j,k), gas0(k) )
#endif
!
#endif
!
            enddo
          enddo  
        enddo 
!
!----------------------------------------------------------!
!
return
end
!
!   ===================================================
   subroutine coag_kernel ( fcv2x, fpd2x, fh2x, fvs2x, rr, coag )
!   ===================================================

   integer :: im, im1
   
   real  :: rpvm1, cv21, pd21, h21
   real  :: rpvm2, cv22, pd22, h22
   real  :: rpav, h2av, pd2av, cv2av
   real  :: coc, grc
   real, dimension(nmode)  :: fcv2x, fpd2x, fh2x, fvs2x, rr
   real, dimension(nmode,nmode) :: coag
!
!-----------------------------------------------------------!
!
  coag = 0.0
  do im = 1,nmode
    rpvm1= rr(im)
    cv21 = fcv2x(im)
    pd21 = fpd2x(im)
    h21  = fh2x (im)
    do im1 = im,nmode
!
!  averaged properties between particles im and im1
!
      rpvm2 = rr(im1)
      rpav  = (rpvm1 + rpvm2) * 0.5
!
      cv22  = fcv2x(im1)
      cv2av = sqrt(cv21 + cv22)
!
      pd22  = fpd2x(im1)
      pd2av = (pd21 + pd22) * 0.5
!
      h22   = fh2x (im1)
      h2av  = sqrt((h21**2 + h22**2) * 0.5)
!
!  Brownian coagulation kernel coc
!
      coc = 16.*pi * pd2av * rpav
!
!  Gravitational coagulation kernel coc
!
      grc = 0.5*pi*rr(im1)**2. * (fvs2x(im1) - fvs2x(im))
!
!  Total kernel
!
      coag(im1,im) = coc + grc
    enddo
  enddo
!
end subroutine
!
!   ===================================================
   subroutine coag_init
!   ===================================================

   integer :: im, im1, iel, ie1, im2
   logical  :: present
   logical, dimension(nelem) 	   :: pres
   logical, dimension(nmode,nelem) :: pres2
!
!-----------------------------------------------------------!
!
!  Create tab_coag of coagulation possibilities: 
!
  pres2=.false.
  do im2 = 1, nmode
    do iel = 1, aero(im2)%nelem
      do ie1 = 1, nelem
  	if ( aero(im2)%init%present(iel) ) pres2(im2,ie1) = .true.
      enddo  
    enddo
  enddo
!
  tab_coag = 0.0
  do im = 1, nmode
    do im1 = im, nmode        
      pres = .false.
      do iel = 1, aero(im)%nelem
        do ie1 = 1, nelem
          if ( aero(im)%init%present(ie1) ) pres(ie1) = .true.
        enddo
      enddo
      do iel = 1, aero(im1)%nelem
        do ie1 = 1, nelem
          if ( aero(im1)%init%present(ie1) ) pres(ie1) = .true.
        enddo
      enddo
      do im2 = 1, nmode
        present=.true.
        do iel = 1, nelem
          present = present .and. ((pres(iel).and.pres2(im2,iel)).or.(.not.pres(iel).and..not.pres2(im2,iel)))
        enddo
        if (present) tab_coag(im,im1,im2) = 1.
        if (max(aero(im)%size%rmean,aero(im1)%size%rmean) > aero(im2)%size%rmean) tab_coag(im,im1,im2) = 0.
	if (aero(im2)%size%rmean > max(aero(im)%size%rmean,aero(im1)%size%rmean) .and. 	&
	    aero(im)%size%rmean == aero(im1)%size%rmean) tab_coag(im,im1,im2) = 0.
      enddo
    enddo
  enddo
!
  end subroutine
!
!   ===================================================
   subroutine coagulation ( nmd, ttn, num, eps, coag, diag )
!   ===================================================

   integer :: im, im1, iel, ie1, ie2, im2
   real  :: coamass, coamass1, coamass2, coamass3, epstot
   logical  :: present
!
   real, dimension(nmode)   :: nmd
   real, dimension(2,nmode) :: num, ttn
   real, dimension(nmode,nmode)   :: coag, com
   real, dimension(2,nmode,nelem) :: eps
   real, dimension(nmode,15)      :: diag
!
!-----------------------------------------------------------!
!
!  Compute coagulation rates
!
            do im = 1, nmode
              do im1 = im, nmode
            	if(num(1,im) > 1.e-6 .and. num(1,im1) > 1.e-6)then  
	          com(im1,im) = coag(im1,im) * num(1,im1) * num(1,im) * dt
!
!  Take care of number
!	  
              	  num(2,im)   = num(2,im)  - com(im1,im) 
              	  num(2,im1)  = num(2,im1) - com(im1,im) 
		  do im2 = 1, nmode
		    num(2,im2) = num(2,im2) + tab_coag(im,im1,im2) * com(im1,im) 
		  enddo
!
!  Take care of mass
!
		  coamass1  = com(im1,im)*pi/6.*aero(im)%mix%rho*nmd(im)**3. 
		  coamass2  = com(im1,im)*pi/6.*aero(im1)%mix%rho*nmd(im1)**3.
		  coamass   = coamass1 + coamass2 
                  ttn(2,im)  = ttn(2,im)  - coamass1
                  ttn(2,im1) = ttn(2,im1) - coamass2
		  diag(im,4)  = diag(im,4) - coamass1/dt
		  diag(im1,4) = diag(im1,4) - coamass2/dt
		  do im2 = 1, nmode
		    ttn(2,im2) = ttn(2,im2) + tab_coag(im,im1,im2) * coamass
		    diag(im2,4) = diag(im2,4) + tab_coag(im,im1,im2) * coamass/dt
		  enddo
!
!  Take care of element fractions
!
		  do im2 = 1, nmode
		    if (tab_coag(im,im1,im2) /= 0. .and. coamass /= 0.) then
  		      epstot = 0.
		      coamass3 = 0.
		      do iel = 1, aero(im2)%nelem
			if ( aero(im2)%init%present(iel) ) coamass3 = coamass3 + eps(1,im,ie1)*coamass1
			if ( aero(im2)%init%present(ie2) ) coamass3 = coamass3 + eps(1,im,ie2)*coamass2
! 
		        eps(2,im2,iel) = (eps(1,im2,iel)*ttn(1,im2) + coamass3) / ttn(2,im2)
		        epstot = epstot + eps(2,im2,iel)
		      enddo
		      eps(2,im2,:) = eps(2,im2,:) / epstot
		    endif
		  enddo
!
		endif
              enddo
	    enddo
!
end subroutine
!
!   ===================================================
   subroutine condensation ( gascon, conf, condterm, ttn, eps, diag )
!   ===================================================

   integer :: im, iel, iname
   real  :: condmass, fs, gascon
   logical  :: present
!
   real, dimension(nmode)   :: conf, condterm
   real, dimension(2,nmode) :: ttn
   real, dimension(2,nmode,nelem) :: eps
   real, dimension(nmode,15)      :: diag
!
!-----------------------------------------------------------!
!
   condmass = 0.0
   do im = 1, nmode
     call have_element (im, 'SULF', present, iname)
!
     if(present .and. conf(im) > 1.e-20) then
       condterm(im) = conf(im) * gascon * dt
       if(condterm(im) < 1.e-20) condterm(im)=0.0
       ttn(2,im) = ttn(2,im) + condterm(im)
       condmass = condmass + condterm(im)
       diag(im,5) = condterm(im)/dt
!
       do iel = 1, aero(im)%nelem
         if ( iel == 1 .and. aero(im)%init%present(1) ) then
           eps(2,im,1) = (eps(1,im,1)*ttn(1,im) + condterm(im)) / ttn(2,im)
         else  
           eps(2,im,iel) = eps(1,im,iel)*ttn(1,im) / ttn(2,im)
         endif
       enddo
     endif
   enddo
! 
   gascon = gascon - condmass
   gascon = max(gascon,0.)
!
end subroutine 
!
!   ===================================================
   subroutine bc_ageing ( condterm, nmd, ttn, num, eps, diag )
!   ===================================================
!
   integer :: im, iel
   real    :: agetermf
!
   real, dimension(nmode)  :: nmd, condterm
   real, dimension(2,nmode) :: num, ttn
   real, dimension(2,nmode,nelem) :: eps
   real, dimension(nmode,15)      :: diag
!
!-----------------------------------------------------------!
!
!  Fossil carbon (bc is the index for pure BC)
!        
            if(ttn(1,bc) > 1.e-12 .and. num(1,bc) > 1.e-6 .and. condterm(bc) > 1.e-20) then
              agetermf = condterm(bc)*(ttn(1,bc)/num(1,bc)) / aerop%bcagemass  
!
!  Update aerosol mass & number
!                       
               num(2,bc)  = num(2,bc) - agetermf / (pi/6.*aero(bc)%mix%rho*nmd(bc)**3.)		     
               num(2,bcs) = num(2,bcs) + agetermf / (pi/6.*aero(bcs)%mix%rho*nmd(bcs)**3.)  
               ttn(2,bc)  = ttn(2,bc) - agetermf  			
               ttn(2,bcs) = ttn(2,bcs) + agetermf 
	       diag(bc,6) = -agetermf
	       diag(bcs,6) = agetermf		     
	       
       	       do iel = 1, aero(bcs)%nelem
	         if ( iel == 2 .and. aero(bcs)%init%present(2) ) then
	           eps(2,bcs,2) = (eps(1,bcs,2)*(ttn(1,bcs) - agetermf) + eps(1,bc,2)*agetermf) / ttn(2,bcs)
	         else
	           eps(2,bcs,iel) = eps(1,bcs,iel)*(ttn(1,bcs) - agetermf) / ttn(2,bcs)
	         endif
               enddo
	     endif
!
  end subroutine 
!
! ===================================================
    subroutine drydep (aerol,i,j)
! ===================================================

! -------------------------------------------------
!  Brief Description of Dummy Variables:
!
!  aeroc%depv0      dry deposition velocity  at 1m (cm s-1)
!  depvl      dry deposition velocity at lowest 
!            model level (cm s-1)
!  ddflux       dry deposition flux (kg m-2 s-1)
! -------------------------------------------------      

      integer :: i, j, it
      integer :: ipstart, ipend, jpstart, jpend
      real    :: aeromass, dens, dp , ddflux
      real    :: depvl, zeta
      real, dimension (nmode) :: ddep

      type (aero_3d), dimension(1:nmode)  :: aerol

! =================================================================

      if (turbu%fkv(i,j,1) /= 0.) then
!      
        dens = pressure%dens(i,j,1)
        dp = pressure%p(i,j,1) - pressure%p(i,j,2)
        zeta = dp/(9.8185*dens)

        do it=1,nmode
          depvl = aero(it)%depv0 * (1./(1+aero(it)%depv0*zeta/turbu%fkv(i,j,1)))
          aeromass = aerol(it)%m(i,j,1)*1.e-9 * dens
          ddflux = depvl * aeromass
          ddep(it) = ddflux * 1.e9 * 9.8185 / dp

          aero3ds(it)%m(i,j,1) = aero3ds(it)%m(i,j,1) - ddep(it)
        enddo
!
      endif
!
      return
    end

    
      !======================================================
        subroutine forcing_SE_Atlantic_aerosols
      !=======================================================
      USE shared_data 
      USE shared_all   
      
      integer :: im, time_simulation,k 

      if (verbose > 2) call write_debug('Starting subroutine forcing_SE_Atlantic_aerosols_time1')
      
      time_simulation=1
      
      do im=1,nmode
        !Aerosol mass concentration is computed directly from the aerosol number concentration      
        FORALL(k=1:nz) aero0(im)%n(k)=1e6*Na_force(lev_len-k+1,time_simulation) 
        FORALL(k=1:nz) aero0(im)%m(k)= aero(im)%size%mmean*aero0(im)%n(k)
        FORALL(k=1:nz) aero0(im)%ma(k) = 0.
      end do
      
      if (verbose > 2) call write_debug('Finishing Forcing Aerosols SE Atlantic time 1')
      
      return
      end        


end module aerosols
