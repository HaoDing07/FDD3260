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
!  CHEMTROP.F:
!	subroutine chemtrop (dtr, ifss, ktrop, uuu0_useless)
!	subroutine tropreact1 (tout, neq, c00, y)
!	subroutine fex (neq, t, y, ydot)
!	subroutine jex (neq, t, y, ml, mu, pd, nrpd)
!	subroutine ppbm2mcm (np, xx, den, revaw)
!	subroutine mcm2ppbm (np, xx, den, aw)
!	subroutine rateph (T, ntemp, Rflux, coszagcm)
!	subroutine ratebm (T, ntemp, patm)
!	subroutine ratetm (T, ntemp, ym)
!	function calkmt (rkt0, rkt00)
!	Block data rktable3dat                    
!
!  Purpose:
!	A package for calculating chemical reactions in troposphere.			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
! ================================================================
!

!	====================================================
	subroutine chemtrop ( ifss, uuu0_useless, gas0, gas )
!	====================================================

! ----------------------------------------------------------
!  Brief Description of Dummy Variables:
!
!	dtr:	time step of integration in second
!	ifss:	1 for using steady-state assumption
!		0 for not using
!       ktrop:  layer label of tropopause
!	Temp:	temperature in Kelvin
!	p:	air pressure in hPa
!	qv:	water vapor density ratio in kg/m^3
!	den:	air density in kg/m^3
!	uuu0:	cos(zenith angle)
!	uuu0_useless:
!		as uuu0 defined in com.h here the dummy uuu0
!		becomes useless
!
! ----------------------------------------------------------

#ifdef SPMD
    USE mpi
#endif
	USE gridno
	USE typedef_gas
	USE typedef_aq
	USE shared_data
	USE shared_state
	USE shared_thermo
	USE shared_rad

	IMPLICIT NONE

	integer :: ipstart, ipend, jpstart, jpend
       integer :: ifss, i, j, k, ind, ntemp, iii
       real    :: pres = 0.0
       real    :: Temp = 0.0
       real    :: den  = 0.0
       real    :: ptam, Rflux, patm,				   &
     		  uuu0_useless, yyy, yyy1, yyy2, 		   &
     		  xxxn2o5, xxxhno3, xxxtotal, ynox0, ynox1, yso20, &
     		  dnox, dso2, oh_k, ddms, dtmp

	real, dimension(nvaria) :: ychem     = (/ (0.0, i=1,nvaria) /)
	real, dimension(nvaria) :: awchem    = (/ (0.0, i=1,nvaria) /)
	real, dimension(nvaria) :: revawchem = (/ (0.0, i=1,nvaria) /)
	real, dimension(2)      :: ychem0    = (/ 0.0, 0.0 /)

	type (gas_chemical), dimension(nz) :: gas0
#if  ( ! defined CHEM_ENABLE )
        type (gas_chemical) :: gas, tmp_gas
#else
        type (gas_chemical),                              &
       dimension(ip_start:ip_end,jp_start:jp_end,1:nz)	  &
     			       :: gas, tmp_gas

#if ( defined SPMD )
	integer :: ierr, status(MPI_STATUS_SIZE), bufsize, nn
     	integer, save :: lastp
#endif

	! --- note: currently 10 diagnostic variables
	integer				:: mnv
	real, dimension (nz,10)		:: tdiag, t2diag
	real, dimension (nz,10)		:: gdiag, g2diag
	real, dimension (nz)		:: wdiag

! ========================================================================

	awchem(1) = awO3
	awchem(2) = awCO
	awchem(3) = awCO2
	awchem(4) = awHO
	awchem(5) = awHO2
	awchem(6) = awH2O2
	awchem(7) = awNO
	awchem(8) = awNO2
	awchem(9) = awCH4
	awchem(10)= awCH2O
	awchem(11)= awCH3O2H
	awchem(12)= awSO2
	awchem(13)= awN2O5
	awchem(14)= awHNO3
	awchem(15)= awH2SO4

	do ind = 1,15
	  revawchem(ind) = 1./awchem(ind)
	end do
       
	yyy       = Avogadro*1.e-3
	yyy1      = yyy/28.964 		!dry air
	yyy2      = yyy/18.0152		!water vapor

	! --- photochemical budget
	do k=1,nz
	do j=jt_start,jt_end
	do i=it_start,it_end
	  tmp_gas(i,j,k)%o3     = gas(i,j,k)%o3*16.0/awO3
	  tmp_gas(i,j,k)%co     = gas(i,j,k)%co*12.0/awCO
	  tmp_gas(i,j,k)%ho     = gas(i,j,k)%ho/awHO                         &
     			       + gas(i,j,k)%ho2/awHO2	       !using ho
         tmp_gas(i,j,k)%xno    = (gas(i,j,k)%xno/awNO 			     &
     				+ gas(i,j,k)%xno2/awNO2)*14.0  !using no
         tmp_gas(i,j,k)%ch2o   = gas(i,j,k)%ch2o*12.0/awCH2O
         tmp_gas(i,j,k)%ch3o2h = gas(i,j,k)%ch3o2h*12.0/awCH3O2H
         tmp_gas(i,j,k)%so2    = gas(i,j,k)%so2*32.0/awSO2
         tmp_gas(i,j,k)%h2so4  = gas(i,j,k)%h2so4*32.0/awH2SO4
         tmp_gas(i,j,k)%hno3   = (gas(i,j,k)%hno3/awHNO3 		     &
     				+ gas(i,j,k)%xn2o5/awN2O5)*14.0 !using hno3
         tmp_gas(i,j,k)%ch4    = gas(i,j,k)%ch4*12.0/awCH4
         tmp_gas(i,j,k)%h2o2   = gas(i,j,k)%h2o2/awH2O2
       end do
       end do
       end do

       t2diag (:,:) = 0.0      
       do k=1,ktrop
       do j=jt_start,jt_end
       do i=it_start,it_end
         t2diag(k, 1) = t2diag(k, 1) + gas(i,j,k)%o3*16.0/awO3
         t2diag(k, 2) = t2diag(k, 2) + gas(i,j,k)%co*12.0/awCO
         t2diag(k, 3) = t2diag(k, 3) + gas(i,j,k)%ho/awHO 		     &
     				     + gas(i,j,k)%ho2/awHO2
         t2diag(k, 4) = t2diag(k, 4) + (gas(i,j,k)%xno /awNO		     &
     				     +  gas(i,j,k)%xno2/awNO2)*14.0
         t2diag(k, 5) = t2diag(k, 5) + gas(i,j,k)%ch2o*12.0/awCH2O
         t2diag(k, 6) = t2diag(k, 6) + gas(i,j,k)%ch3o2h*12.0/awCH3O2H
         t2diag(k, 7) = t2diag(k, 7) + gas(i,j,k)%so2*32.0/awSO2
         t2diag(k, 8) = t2diag(k, 8) + gas(i,j,k)%h2so4*32.0/awH2SO4
         t2diag(k, 9) = t2diag(k, 9) + (gas(i,j,k)%hno3 /awHNO3 	     &
     				     +  gas(i,j,k)%xn2o5/awN2O5)*14.0
	  t2diag(k,10) = t2diag(k,10) + gas(i,j,k)%h2o2/awH2O2
	end do
	end do
	end do

	do k = 1,ktrop
	do j = jt_start,jt_end
	do i = it_start,it_end

	  ! === Patch for CM5:
	  pres = pressure2%p(i,j,k)+p0(k)
	  den  = pres/(2.87e2*thermo%T(i,j,k))
	  Temp = thermo%T(i,j,k)
	  ! ==================

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

	  ! --- Convert air density from kg/m^3 to mols/cm^3:
	  ychem0(1) = den*yyy1		!dry air
	  ychem0(2) = qv*den *yyy2	!water vapor

	  ! --- Convert mixing ratio from ppb(m) to mols/cm^3:
	  ychem( 1) = gas(i,j,k)%o3
	  ychem( 2) = gas(i,j,k)%co
	  ychem( 3) = gas(i,j,k)%zco2
	  ychem( 4) = gas(i,j,k)%ho
	  ychem( 5) = gas(i,j,k)%ho2
	  ychem( 6) = gas(i,j,k)%h2o2
	  ychem( 7) = gas(i,j,k)%xno
	  ychem( 8) = gas(i,j,k)%xno2
	  ychem( 9) = gas(i,j,k)%ch4
	  ychem(10) = gas(i,j,k)%ch2o
	  ychem(11) = gas(i,j,k)%ch3o2h
	  ychem(12) = gas(i,j,k)%so2

	  ychem(13) = gas(i,j,k)%xn2o5
	  ychem(14) = gas(i,j,k)%hno3
	  ychem(15) = gas(i,j,k)%h2so4

	  do iii = 1,nvaria
	    if(ychem(iii).lt.0.0) ychem(iii) = 0.0
	  enddo

	  call ppbm2mcm (nvaria, ychem, den, revawchem)

	  ! --- Calculate reaction rates:
	  ntemp = nint((Temp - 200.0)*2.0)
	  if(ntemp.le.1)   ntemp = 1
	  if(ntemp.ge.300) ntemp = 300

	  ! ===
	  ! 042596 add zenith angle function:
	  call rateph (Temp, ntemp, Rflux, uuu0)
	  call ratebm (Temp, ntemp, patm)
	  call ratetm (Temp, ntemp, ychem0)

!	  if(ifss.eq.1) call chemsteady(nvaria, ychem0,ychem)

	  ! ===
	  ! === New scheme for calculating NOy and S(VI)
	  ! ===     Chien Wang, 092195
	  ! ===
	  ! --- weighting d[N2O5] and d[HNO3]:
	  xxxn2o5 = rk(14)*ychem(1)	!  k14[O3]
	  xxxhno3 = rk(13)*ychem(4)	!  k13[HO]
	  xxxtotal= 2*xxxn2o5 + xxxhno3	! [NO2] = 1./
					! (2xxxn2o5 + xxxhno3)
	  ! --- save initial NOx and SO2:
	  if ( ychem(7) .lt. 0.0 ) ychem(7) = 0.0
	  if ( ychem(8) .lt. 0.0 ) ychem(8) = 0.0
	  ynox0 = ychem(7) + ychem(8)
	  yso20 = ychem(12)

! ===
	  call tropreact1(12, ychem0, ychem)
! ===

	  ! --- Check if NOx budget has error:
	  if ( ychem(7) .lt. 0.0 ) ychem(7) = 0.0
	  if ( ychem(8) .lt. 0.0 ) ychem(8) = 0.0
	  ynox1= ychem(7) + ychem(8) 
	  dnox = ynox1 - ynox0
	  if ( dnox .gt. 0.0 .and. ynox1 .gt. 0.0 ) then
		dtmp = (1.0 - dnox/ynox1)
		ychem(7) = ychem(7)*dtmp
		ychem(8) = ychem(8)*dtmp
		if ( ychem(7) .lt. 0.0 ) ychem(7) = 0.0
		if ( ychem(8) .lt. 0.0 ) ychem(8) = 0.0
	  end if
	  dnox = min(0.0, dnox)
	  ! --- Derive [N2O5], [HNO3], and [H2SO4]; 092195:
	  if(xxxtotal.gt.1.e-16)then
	    ychem(13) = ychem(13) - xxxn2o5/xxxtotal * dnox 
	    ychem(14) = ychem(14) - xxxhno3/xxxtotal * dnox
	  endif

	  dso2 = ychem(12) - yso20
	  if(dso2.gt.0.0) dso2 = 0.0
	  ychem(15) = ychem(15) - dso2

	  ! --- 112498, Calculate DMS + OH -> product 
	  ! ---   then SO2 production = yield*product
	  ! --- Note: [OH]*k is dimensionless.
	  oh_k = exp(-1.2e-11*exp(-260.0/Temp)*ychem(4)*time)
	  if(oh_k.ge.0.999999) oh_k = 0.999999	!0 <= oh_k< 1
	  if(oh_k.le.0.0)      oh_k = 0.0

	  if(Temp.lt.200.0)then
	    ddms = 200.0
	  else if(Temp.gt.300.0)then
	    ddms = 300.0
	  else
	    ddms = Temp
	  endif

	  ddms = (0.8446*ddms - 168.85)*0.01		!yield (0-1)
	  ddms = gas(i,j,k)%dms*(1.0 - oh_k)*ddms * awSO2/awDMS

	  gas(i,j,k)%dms = gas(i,j,k)%dms*oh_k

	  ! --- Convert concentration to ppb(m):
	  call mcm2ppbm (nvaria, ychem, den, awchem)

	  do iii = 1,nvaria
	    if(ychem(iii).lt.0.0) ychem(iii) = 0.0
	  end do

	  ! --- Give value back to the common block:
	  gas(i,j,k)%o3     = ychem(01)
	  gas(i,j,k)%co     = ychem(02)
	  gas(i,j,k)%zco2   = ychem(03)
	  gas(i,j,k)%ho     = ychem(04)
	  gas(i,j,k)%ho2    = ychem(05)
	  gas(i,j,k)%h2o2   = ychem(06)
	  gas(i,j,k)%xno    = ychem(07)
	  gas(i,j,k)%xno2   = ychem(08)
	  gas(i,j,k)%ch4    = ychem(09)
	  gas(i,j,k)%ch2o   = ychem(10)
	  gas(i,j,k)%ch3o2h = ychem(11)
	  gas(i,j,k)%so2    = ychem(12)                             &
     			   + ddms	   ! 112498
         gas(i,j,k)%xn2o5  = ychem(13)
         gas(i,j,k)%hno3   = ychem(14)
         gas(i,j,k)%h2so4  = ychem(15)

       end do
       end do
       end do
!      ----------------

       ! --- photochemical budget
       tdiag (:,:) = 0.0
       do k=1,ktrop
       do j=jt_start,jt_end
       do i=it_start,it_end
         tdiag(k, 1) = tdiag(k, 1)				    &
     		     + gas(i,j,k)%o3*16.0/awO3  		    &
     		     - tmp_gas(i,j,k)%o3
         tdiag(k, 2) = tdiag(k, 2)				    &
     		     + gas(i,j,k)%co*12.0/awCO  		    &
     		     - tmp_gas(i,j,k)%co			    &
     		     + (gas(i,j,k)%ch4*12.0/awCH4 		    &
     		     - tmp_gas(i,j,k)%ch4)
         tdiag(k, 3) = tdiag(k, 3)				    &
     		     + gas(i,j,k)%ho/awHO + gas(i,j,k)%ho2/awHO2    &
     		     - tmp_gas(i,j,k)%ho
         tdiag(k, 4) = tdiag(k, 4)				    &
     		     + (gas(i,j,k)%xno/awNO 			    &
     		      + gas(i,j,k)%xno2/awNO2)*14.0 		    &
     		     - tmp_gas(i,j,k)%xno
         tdiag(k, 5) = tdiag(k, 5)				    &
     		     + gas(i,j,k)%ch2o*12.0/awCH2O 		    &
     		     - tmp_gas(i,j,k)%ch2o
         tdiag(k, 6) = tdiag(k, 6)				    &
     		     + gas(i,j,k)%ch3o2h*12.0/awCH3O2H  	    &
     		     - tmp_gas(i,j,k)%ch3o2h
         tdiag(k, 7) = tdiag(k, 7)				    &
     		     + gas(i,j,k)%so2*32.0/awSO2 		    &
     		     - tmp_gas(i,j,k)%so2
         tdiag(k, 8) = tdiag(k, 8)				    &
     		     + gas(i,j,k)%h2so4*32.0/awH2SO4 		    &
     		     - tmp_gas(i,j,k)%h2so4
         tdiag(k, 9) = tdiag(k, 9)				    &
     		     + (gas(i,j,k)%hno3/awHNO3  		    &
     		     + gas(i,j,k)%xn2o5/awN2O5)*14.0 		    &
     		     - tmp_gas(i,j,k)%hno3
         tdiag(k,10) = tdiag(k,10)				    &
     		     + gas(i,j,k)%h2o2/awH2O2 			    &
     		     - tmp_gas(i,j,k)%h2o2
	end do
	end do
	end do

#if ( defined SPMD )
	 call MPI_BARRIER ( MPI_COMM_WORLD, ierr )
	 bufsize = nz*10
	 call MPI_REDUCE (tdiag(1,1),gdiag(1,1),bufsize,MPI_REAL,   &
     		          MPI_SUM,0, MPI_COMM_WORLD, ierr)

	 call MPI_REDUCE (t2diag(1,1),g2diag(1,1),bufsize,MPI_REAL, &
     		          MPI_SUM,0, MPI_COMM_WORLD, ierr)
#else
	gdiag (:,:) = tdiag (:,:)
	g2diag(:,:) = t2diag(:,:)
#endif

	! --- write gaseous diag
	if ( mypid .eq. 0 ) then	
	do mnv = 1,10
	  wdiag(1:nz) = gdiag (1:nz,mnv)	  
	  write(150) wdiag
	  wdiag(1:nz) = g2diag(1:nz,mnv)	  
	  write(150) wdiag
	end do
	end if

#endif
	return
	 end
 
!	=========================================
	subroutine tropreact1 (neq, c00, y)
!	=========================================

! --------------------------------------------------------
! --- TROPREACT1.F:  Subroutine for calculating chemical
! --- 			reactions in troposphere
! ---
! --- Last Revised on:    July 19, 1998
! --------------------------------------------------------

	external fex, jex

	real :: y(neq), rwork(872), iwork(45), c00(2), c00tmp(2), time
	common /c00tmp/c00tmp

	integer :: neq, itol, itask, istate, iopt, lrw, liw, mf
	real    :: t, rtol, atol 

	c00tmp(1) = c00(1)
	c00tmp(2) = c00(2)

	t = 0.
	itol = 1
	rtol = 1.e-3
	atol = 1.e+3	!molecules/cm3

!	atol(1) = 1.e-6
!	atol(2) = 1.e-10
!	atol(3) = 1.e-6

	itask  = 1
	istate = 1
	iopt   = 0
	lrw    = 872
	liw    = 45
	mf     = 21

!	call lsodenew(fex,neq,y,t,time,itol,rtol,atol,itask,istate,  &
!          iopt,rwork,lrw,iwork,liw,jex,mf)

	return
	 end

!	================================
	subroutine fex (neq, t, y, ydot)
!	================================

!
! --- Subroutine for defining of ydot
!
     
	USE gridno
	USE shared_data

	integer :: neq
	real    :: t, c00(2), y(neq), ydot(neq)      
	common /c00tmp/c00

	real    :: bbb, ccc, xbeta1, xbeta2, xbeta,             &
     		  rr1, rr5, rr7, rr8, rr10,rr11,rr12,rr13,	&
     		  rr14,rr16,rr17,rr18,rr19,rr21,rr23,rr24,	&
     		  rr25,rr26,rr27,rr29,rr32,rr33,rr34,		&
     		  xxx, gamax,gamay

        bbb = 0.7809*c00(1)
        ccc = 0.2095*c00(1)

! === redesigned scheme, 092195
!
	xbeta1 = rk(3)*bbb + rk(4)*ccc	!k3[N2] + k4[O2]
	xbeta2 = rk(2)*c00(2)		!k2[H2O]
	xbeta  = xbeta2/(xbeta1 + xbeta2)

	rr1 = rk(1) *y(1)		!j1[O3]
	rr5 = rk(5) *y(2) *y(4)		!k5[CO][HO]
	rr7 = rk(7) *y(5) *y(7)		!k7[HO2][NO]
	rr8 = rk(8) *y(8)		!k8[NO2]
	rr10= rk(10)*y(1) *y(5)		!k10[O3][HO2]
	rr11= rk(11)*y(1) *y(4)		!k11[O3][HO]
	rr12= rk(12)*y(1) *y(7)		!k12[O3][NO]
	rr13= rk(13)*y(4) *y(8)		!k13[HO][NO2]
	rr14= rk(14)*y(1) *y(8)		!k14[O3][NO2]
	rr16= rk(16)*y(5) *y(5)		!k16[HO2][HO2]
	rr17= rk(17)*y(6)		!j17[H2O2]
	rr18= rk(18)*y(4) *y(6)		!k18[HO][H2O2]
	rr19= rk(19)*y(4) *y(9)		!k19[HO][CH4]
	rr21= rk(21)*y(7)		!k21[NO] - for ss only
	rr23= rk(23)*y(5)		!k23[HO2]- for ss only
	rr24= rk(24)*y(11)		!j24[CH3O2H]
	rr25= rk(25)*y(4) *y(11)	!k25[HO][CH3O2H]
	rr26= rk(26)*y(10)		!j26[CH2O]
	rr27= rk(27)*y(4) *y(10)	!k27[HO][CH2O]
	rr29= rk(29)*y(4) *y(12)	!k29[HO][SO2]

! --- New reactions - 062995:

	rr32= rk(32)*y(4)*y(5)		!k32[HO][HO2]
	rr33= rk(33)*y(4)*y(4)		!k33[HO][HO]
	rr34= rk(34)*y(4)*y(4)		!k34[HO][HO]

	xxx   = (rr21 + rr23)

	gamax = 0.0
	if(rr21.gt.1.e-12) gamax = rr21/xxx

	gamay = 0.0
	if(xxx.gt.0.0)gamay = rr23/xxx

	ydot(1)  = rr8 + rr33		 		      &
     		-(xbeta*rr1 + rr10 + rr11 + rr12 + rr14)	   	

       ydot(2)  = rr26 + rr27			       	  &
     		- rr5

       ydot(3)  = rr5				       

       ydot(4)  = 2.0*xbeta*rr1 + rr7 + rr10	       	  &
     		+ 2.0*rr17+ rr24				  &
     		-(rr5 + rr11 + rr13 + rr18 + rr19 		  &
     		 +rr25+ rr27 + rr29				  &
     		 +rr32+ rr33 + rr34 )

       ydot(5)  = rr5 + rr11 + rr18 + rr27 + rr29      	  &
     		+ rr24+ 2.0*rr26				  &
     		+(2.0*gamax - 1.0)*(rr19 + rr25)		  &
     		-(rr7 + rr10 + rr16 + rr32)

       ydot(6)  = rr16 + rr34			       	  &
     		-(rr17+ rr18)

       ydot(7)  = rr8				       	  &
     		-(rr7 + rr12 + gamax*(rr19 + rr25))

       ydot(8)  = rr7 + rr12 + gamax*(rr19 + rr25)     	  &
     		-(rr8 + rr13+ rr14+ rr14)

       ydot(9)  =-rr19  			       

       ydot(10) = rr24 + gamax*(rr19 + rr25)	       	  &
     		-(rr26 + rr27)

       ydot(11) = gamay*(rr19 + rr25)		       	  &
     		-(rr24 + rr25) 

	ydot(12) =-rr29 				

	return
	 end

!	============================================
	subroutine jex (neq, t, y, ml, mu, pd, nrpd)
!	============================================

! -----------------------------------------------------
! --- Subroutine for calculating the Jacobian for f(i)
! ---
! ---  P(m,n) here represents df(m)/dn
! --- 		or d ydot(m)/dn
! -----------------------------------------------------

	USE gridno
	USE shared_data

	integer :: neq, ml, mu, nrpd, i, j
	real    :: t, c00(2), y(neq), pd(nrpd,neq)
	common /c00tmp/c00

	real    :: bbb, ccc, xbeta1, xbeta2, xbeta, xxx,          &
     		  rr21,rr23,gamax,  gamay,			  &
     		  rk5y2, rk5y4, rk7y5, rk7y7, rk10y1,rk10y5,	  &
     		  rk11y1,rk11y4,rk12y1,rk12y7,rk13y4,rk13y8,	  &
     		  rk14y1,rk14y8,rk18y4,rk18y6,rk19y4,rk19y9,	  &
     		  rk25y4,rk25y11,rk27y4,rk27y10,rk29y4,rk29y12,   &
     		  rk32y4,rk32y5

	bbb = 0.7809*c00(1)
	ccc = 0.2095*c00(1)

        xbeta1 = rk(3)*bbb + rk(4)*ccc   !k3[N2] + k4[O2]
        xbeta2 = rk(2)*c00(2)            !k2[H2O]
        xbeta  = xbeta2/(xbeta1 + xbeta2)

        rr21= rk(21)*y(7)               !k21[NO] - for ss only
        rr23= rk(23)*y(5)               !k23[HO2]- for ss only
        xxx   = (rr21 + rr23)

        gamax = 0.0
        if(rr21.gt.1.e-12) gamax = rr21/xxx

        gamay = 0.0
        if(xxx.gt.0.0) gamay = rr23/xxx

	do j=1,neq
	do i=1,neq
	  pd(i,j) = 0.0
	end do
	end do

! ===
! === pd(i,j) = d ydot(i)/dy(j):
! ===	
	rk5y2    = rk(5)*y(2)
	rk5y4    = rk(5)*y(4)
	rk7y5    = rk(7)*y(5)
	rk7y7    = rk(7)*y(7)
	rk10y1   = rk(10)*y(1)
	rk10y5   = rk(10)*y(5)
	rk11y1   = rk(11)*y(1)
	rk11y4   = rk(11)*y(4)
	rk12y1   = rk(12)*y(1)
	rk12y7   = rk(12)*y(7)
	rk13y4   = rk(13)*y(4)
	rk13y8   = rk(13)*y(8)
	rk14y1   = rk(14)*y(1)
	rk14y8   = rk(14)*y(8)
	rk18y4   = rk(18)*y(4)
	rk18y6   = rk(18)*y(6)
	rk19y4   = rk(19)*y(4)
	rk19y9   = rk(19)*y(9)
	rk25y4   = rk(25)*y(4)
	rk25y11  = rk(25)*y(11)
	rk27y4   = rk(27)*y(4)
	rk27y10  = rk(27)*y(10)
	rk29y4   = rk(29)*y(4)
	rk29y12  = rk(29)*y(12)
	rk32y4   = rk(32)*y(4)
	rk32y5   = rk(32)*y(5)

	pd(1,1)  = -(xbeta*rk(1)   + rk10y5 + rk11y4       &
     			 +rk12y7  + rk14y8)
       pd(4,1)  = 2.0*xbeta*rk(1) + rk10y5 - rk11y4
       pd(5,1)  = rk11y4 - rk10y5
       pd(8,1)  = rk12y7 - 2.*rk14y8

       pd(2,2)  =-rk5y4
       pd(3,2)  = rk5y4
       pd(4,2)  =-rk5y4
       pd(5,2)  = rk5y4

       pd(1,4)  = 2.0*rk(33)*y(4) - rk11y1
       pd(2,4)  = rk27y10    - rk5y2
       pd(3,4)  = rk5y2
       pd(4,4)  =-( 					   &
     		  rk5y2  + rk11y1 + rk(13)*y(8) 	   &
     		+ rk18y6 + rk19y9 + rk25y11		   &
     		+ rk27y10+ rk29y12+ rk32y5		   &
     		+(rk(33) + rk(34))*y(4) 		   &
     		   )
       pd(5,4)  =(rk5y2 + rk11y1 + rk18y6		   &
     		+ rk27y10+ rk29y12)			   &
     		+(2.0*gamax - 1.0)*(rk19y9 + rk25y11)	   &
     		- rk32y5
       pd(6,4)  =-rk18y6
       pd(7,4)  =-gamax*(rk19y9 + rk25y11)
       pd(8,4)  = gamax*(rk19y9 + rk25y11)		   &
     		- rk13y8
       pd(9,4)  =-rk19y9
       pd(10,4) = gamax*(rk19y9 + rk25y11)		   &
     		- rk27y10
       pd(11,4) = gamay*(rk19y9 + rk25y11)		   &
     		- rk25y11
       pd(12,4) =-rk29y12

       pd(1,5)  = rk11y1
       pd(4,5)  = rk7y7 + rk10y1			   &
     		- rk32y4
       pd(5,5)  =					   &
     		-(rk7y7 + rk10y1 + rk(16)*y(5) + rk32y4)
       pd(6,5)  = 2.0*   (rk(16) + rk(34))*y(5)
       pd(7,5)  =-rk7y7
       pd(8,5)  = rk7y7
!      pd(11,5)

       pd(4,6)  = 2.0*rk(17)				   &
     		- rk18y4
	pd(5,6)  = rk18y4
	pd(6,6)  =-(rk(17) + rk18y4)

	pd(1,7)  =-rk12y1
	pd(4,7)  = rk7y5
	pd(5,7)  =-rk7y5
	pd(7,7)  =-(rk7y5 + rk12y1)
	pd(8,7)  = (rk7y5 + rk12y1)
	
	pd(1,8)  = rk(8) - rk14y1
	pd(4,8)  =-rk13y4
	pd(7,8)  = rk(8)
	pd(8,8)  =-(rk(8) + rk13y4 + 2.0*rk14y1)

	pd(4,9)  =-rk19y4
	pd(5,9)  =(2.0*gamax - 1.0)*rk19y4
	pd(7,9)  =-gamax*rk19y4
	pd(8,9)  = gamax*rk19y4
	pd(9,9)  =-rk19y4
	pd(10,9) = gamax*rk19y4
	pd(11,9) = gamay*rk19y4

	pd(2,10) = rk(26) + rk27y4
	pd(4,10) =-rk27y4
	pd(5,10) = rk27y4 + 2.0*rk(26)
	pd(10,10)=-(rk(26)+ rk27y4)

	pd(4,11) = rk(24) - rk25y4
	pd(5,11) = rk(24) + (2.0*gamax - 1.0)*rk25y4
	pd(7,11) =-gamax*rk25y4
	pd(8,11) = gamax*rk25y4
	pd(10,11)= rk(24) + gamax*rk25y4
	pd(11,11)=-(rk(24) + rk25y4*max(0.0,1.0 - gamay))

	pd(4,12) =-rk(29)*y(4)
	pd(5,12) = rk29y4
	pd(12,12)=-rk29y4

      return
      end

!	=============================================
	subroutine rateph (T, ntemp, Rflux, coszagcm)
!	=============================================

! --------------------------------------------------
! --- A program for calculating rate constants
! --- 	of photochemical reactions
! ---
! --- The time unit for all reaction rates is second
! ---
! --- Last revised: July 19, 1998
! ---
! --------------------------------------------------

	USE gridno
	USE shared_data

	integer :: ntemp, nband
	real    :: T, Rflux, coszagcm, wgtxx, wgt1, wgt2
#include "rktable.h"
       
	if(coszagcm.lt.cosza4rk(10)) coszagcm = cosza4rk(10)
	if(coszagcm.gt.cosza4rk(1))  coszagcm = cosza4rk(1)

! --- For znith angle.lt.86:
	if(coszagcm.ne.cosza4rk(10))then

	do 1 nband=1,9

	if(coszagcm.le.cosza4rk(nband).and.   &
          coszagcm.gt.cosza4rk(nband+1))then
	  wgtxx = 1.0/(cosza4rk(nband) - cosza4rk(nband+1))
	  wgt1  = (cosza4rk(nband) - coszagcm)  *wgtxx
	  wgt2  = (coszagcm - cosza4rk(nband+1))*wgtxx
!
!	+--- wgt2 ---+--- wgt1 ---+
!	nband --- coszagcm --- nband+1

! --- R8: NO2 + hv -> NO + O
!	Atkinson et al., 1992:
!
!	  rk(8) = (18096.7314929 - 3.9704*(T-273.))
!     &        *5.0e-10
!     &        *Rflux

!	  rk(8) = rktable1(08,ntemp)*Rflux

	  rk(8) = ((rk08gama(nband)*wgt2 + rk08gama(nband+1)*wgt1)  &
     	       +  (rk08aaa (nband)*wgt2 + rk08aaa (nband+1)*wgt1)   &
     	       *(T-273.))*Rflux

! --- R1: O3 + hv -> O(1D) + O2
!	Atkinson et al., 1992:
!
!	  rk(1)  = 0.0028*rk(8)
!!	  rk(1) = 3.467625419939578E-08*Rflux

!	  rk(1) = (rk01table(nband)*wgt2 
!     &          + rk01table(nband+1)*wgt1)*Rflux
	  rk(1) = 0.0028*rk(8)

! --- R17: H2O2 + hv -> 2OH
!	Atkinson et al., 1992:
!
!	  rk(17)= 6.263464249748237E-09*Rflux

	  rk(17) = (rk17table(nband)*wgt2 +  rk17table(nband+1)*wgt1)*Rflux

! --- R24: CH3O2H + hv -> CH3O + HO
!	Atkinson et al., 1992:
!
!	  rk(24)= 4.474208962739173E-09*Rflux

	  rk(24) = (rk24table(nband)*wgt2 +  rk24table(nband+1)*wgt1)*Rflux

! --- R26: CH2O + hv -> CHO + H
!	Atkinson et al., 1992:
!
!	  rk(26)= 9.410107812688823E-08*Rflux

	  rk(26) = (rk26table(nband)*wgt2 +  rk26table(nband+1)*wgt1)*Rflux

	endif

1	continue

	else
! --- For znith angle.ge.86:

	  rk(8) = (rk08gama(10)    &
     	       +  rk08aaa (10)	   &
     	       *(T-273.))	   &
     	       *Rflux

!        rk(1) = rk01table(10)	    
!    &  	* Rflux 	    
     	 rk(1)  = 0.0028*rk(8)

         rk(17) = rk17table(10)    &
     		* Rflux

         rk(24) = rk24table(10)    &
     		* Rflux

         rk(26) = rk26table(10)    &
     		*Rflux

	endif

	return
	 end

!	==================================
	subroutine ratebm (T, ntemp, patm)
!	==================================

! --------------------------------------------------
! --- A program for calculating rate constants
! --- 	of bimolecular reactions
! ---
! --- The time unit for all reaction rates is second
! ---
! --------------------------------------------------

	USE gridno
	USE shared_data

	integer :: ntemp
	real    :: T, patm, Trev
#include "rktable.h"

	Trev = 1./T

! --- R2: O(1D) + H2O -> 2OH

	rk(2)  = 2.2e-10

! --- R3: O(1D) + N2 -> O + N2

!	rk(3)  = 1.8e-11 * exp( 107.0*Trev)
	rk(3)  = rktable1(03,ntemp)

! --- R4: O(1D) + O2 -> O + O2

!	rk(4)  = 3.2e-11 * exp(  67.0*Trev)
	rk(4)  = rktable1(04,ntemp)

! --- R5: CO + HO -> H + CO2

	rk(5)  = 1.5e-13 * (1.0 + 0.6*patm)

! --- R7: HO2 + NO -> HO + NO2

!	rk(7)  = 3.7e-12 * exp( 240.0*Trev)
	rk(7)  = rktable1(07,ntemp)

! --- R10: HO2 + O3 -> HO + 2O2

!	rk(10) = 1.1e-14 * exp(-500.0*Trev)
	rk(10) = rktable1(10,ntemp)

! --- R11: HO + O3 -> HO2 + O2

!	rk(11) = 1.6e-12 * exp(-940.0*Trev)
	rk(11) = rktable1(11,ntemp)
!	rk(11) = 0.0

! --- R12: NO + O3 -> NO2 + O2

!	rk(12) = 2.0e-12 * exp(-1400.0*Trev)
	rk(12) = rktable1(12,ntemp)

! --- R14: NO2 + O3 -> NO3 + O2

!	rk(14) = 1.2e-13 * exp(-2450.0*Trev)
	rk(14) = rktable1(14,ntemp)

! --- R16: HO2 + HO2 -> H2O2 + O2

!	rk(16) = 2.3e-13 * exp( 600.0*Trev)

! --- R18: H2O2 + HO -> HO2 + H2O

!	rk(18) = 2.9e-12 * exp(-160.0*Trev)
	rk(18) = rktable1(18,ntemp)

! --- R19: CH4 + HO -> CH3 + H2O

!	rk(19) = 2.65e-12 * exp(-1800.0*Trev)
	rk(19) = rktable1(19,ntemp)

! --- R21: CH3O2 + NO -> CH3O + NO2

!	rk(21) = 4.2e-12 * exp( 180.0*Trev)
	rk(21) = rktable1(21,ntemp)

! --- R22: CH3O + O2 -> CH2O + HO2

!	rk(22) = 3.9e-14 * exp(-900.0*Trev)
	rk(22) = rktable1(22,ntemp)

! --- R23: CH3O2 + HO2 -> CH3O2H + O2

!	rk(23) = 3.8e-13 * exp( 780.0*Trev)
	rk(23) = rktable1(23,ntemp)

! --- R25: CH3O2H + HO -> CH3O2 + H2O

!	rk(25) = 1.9e-12 * exp( 190.0*Trev)
	rk(25) = rktable1(25,ntemp)

! --- R27: CH2O + HO -> CHO + H2O

	rk(27) = 1.0e-11

! --- R28: CHO + O2 -> CO + HO2

!	rk(28) = 3.5e-12 * exp( 140.0*Trev)
	rk(28) = rktable1(28,ntemp)

! --- R30: HOSO2 + O2 -> HO2 + SO3

!	rk(30) = 1.3e-12 * exp(-330.0*Trev)
	rk(30) = rktable1(30,ntemp)

! --- R31: SO3 + H2O -> H2SO4

	rk(31) = 2.4e-15

! === New reactions - 062995:

! --- R32: HO + HO2 -> H2O + O2

!	rk(32) = 4.8e-11 * exp( 250.0*Trev)
	rk(32) = rktable1(32,ntemp)

! --- R33: HO + HO -> H2O + O

!	rk(33) = 4.2e-12 * exp(-240.0*Trev)
	rk(33) = rktable1(33,ntemp)

	return
	 end

!	================================
	subroutine ratetm (T, ntemp, ym)
!	================================

! --------------------------------------------------
! --- A program for calculating rate constants
! --- 	of termolecular reactions
! ---
! --- The time unit for all reaction rates is second
! --------------------------------------------------

	USE gridno
	USE shared_data

	integer :: ntemp
	real    :: T, ym(2), T300, xm, xn2, h2o, rkt0, rkt00
	real, external :: calkmt
#include "rktable.h"

        T300 = T/300.
        xm   = ym(1)
        xn2  = 0.7809*xm
        h2o  = ym(2)

! --- R6:  H + O2 + M -> HO2 + M

!	rkt0   = 6.2e-32 * xm * T300**(-1.6)
	rkt0   = rktable1(06,ntemp) * xm
	rkt00  = 7.5e-11
	rk(6)  = calkmt(rkt0, rkt00)

! --- R9:  O + O2 + M -> O3 + M  (no high pres. limit)

!	rkt0   = 5.6e-34 * xm * T300**(-2.8)
	rkt0   = rktable1(09,ntemp) * xm
	rk(9)  = rkt0

! --- R13: NO2 + HO + M -> HNO3 + M

!	rkt0   = 2.6e-30 * xm * T300**(-3.2)
!	rkt00  = 2.4e-11 * T300**(-1.3)
	rkt0   = rktable1(13,ntemp) * xm
	rkt00  = rktable2(1, ntemp)
	rk(13) = calkmt(rkt0, rkt00)

! --- R15: NO3 + NO2 + M -> N2O5 + M

!	rkt0   = 2.2e-30 * xm * T300**(-3.9)
!	rkt00  = 1.5e-12 * T300**(-0.7)
	rkt0   = rktable1(15,ntemp) * xm
	rkt00  = rktable2(2, ntemp)
	rk(15) = calkmt(rkt0, rkt00)

! --- R16: HO2 + HO2 + M -> H2O2 + O2
! --- 062895:

!	rk(16) = 1.7e-33 * xn2 * exp(1000.0/T)
	rk(16) = rktable1(16,ntemp) * xn2 * (1.+ 1.4e-21*h2o*exp(2200.0/T))

! --- R20: CH3 + O2 + M -> CH3O2 + M

!	rkt0   = 4.5e-31 * xm * T300**(-3.0)
!	rkt00  = 1.8e-12 * T300**(-1.7)
	rkt0   = rktable1(20,ntemp) * xm
	rkt00  = rktable2(3, ntemp)
	rk(20) = calkmt(rkt0, rkt00)

! --- R29: SO2 + HO + M -> HOSO2 + M

!	rkt0   = 3.0e-31 * xm * T300**(-3.3)
	rkt0   = rktable1(29,ntemp) * xm
	rkt00  = 1.5e-12
	rk(29) = calkmt(rkt0, rkt00)

! === New reactions, 063095:

! --- R34: HO + HO + M -> H2O2 + M

!        rkt0   = 6.9e-31 * xm * T300**(-0.8)
        rkt0   = rktable1(34,ntemp) * xm
        rkt00  = 1.5e-11
        rk(34) = calkmt(rkt0, rkt00)


	return
	 end

!	=============================
	real function calkmt (rkt0, rkt00)
!	=============================

!
!  A function for calculate k(M,T) of termolecular
!    reaction rate from low and high limit
!
	real :: rkt0, rkt00, aaa, bbb, ccc

	aaa = rkt0/rkt00
	bbb = log10(aaa)**2
	ccc = 1./(1.+bbb)

	calkmt = rkt0/(1.+aaa) * 0.6**ccc

	return
	 end


