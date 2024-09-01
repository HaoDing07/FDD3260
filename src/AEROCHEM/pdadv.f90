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
!  PDADV.F90:
!	subroutine pdadv1(c,c_w4,c_w2,c_w1,N)
!	subroutine pdadv2(c,q,c_w4,c_w2,c_w1,c_ww,N,NOOS)                
!
!  Purpose:
!	A package of the improved Bott's advection scheme based on
!		Chien Wang and Xingzhong Liang's work.			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
! ================================================================

	module pdadv_index
	
	SAVE
	
     	integer :: index_p(5000), index_n(5000), nip, nin
	
	end module pdadv_index
	
!	=====================================
 	subroutine pdadv1(c,c_w4,c_w2,c_w1,n)
!	=====================================

! -------------------------------------------------------------
! --- This is a subroutine for the first part 
! --- 	of Bott's advection scheme.
! ---
! --- Andreas Bott  1989:  A Positive Definite 
! --- 		Advection scheme obtained by Nonlinear
! --- 		Renormalization of the advective fluxes
! --- 		Mon. Wea. Rev. 117 1006-15
! ---
! --- Fourth Order: with coefficients from    
! --- 		Mon. Wea. Rev. 117 2633-36
! ---
! --- Input: C=U*DT/DX[N+1]  Output: W4[3:N1,5],W2[2;3;n1;n,3]
! ---				 and   W1[1;2;n;n+1,2]
! ---   On the Staggered Grid: C(i')----Q(i)----C(i'+1)
! -------------------------------------------------------------

	USE pdadv_index
	
	real, parameter :: C0=1.0/1920.0,C1=1.0/384.00,C2=1.0/384.0
     	real, parameter :: C3=1.0/768.00,C4=1.0/3840.0,EP=1.0E-15
	real, parameter :: cc0=-1./24.,cc1=1./16.,cc2=1./16.

	integer :: N, n1, n2, n3, i 
	real    :: c(n+1),c_w4(3:n-1,5),c_w2(4,3),c_w1(4,2)
     	real    :: rr1, rr2, r1, r2, r3, r4
     

	n1 = n-1
	n2 = n-2
	n3 = n-3

	index_p(1:n+1) = 0
	index_n(1:n+1) = 0
	nip            = 0
	nin            = 0	
	do i=1,n+1
	  if(c(i).gt.0.0)then
	    nip = nip + 1
	    index_p(nip) = i
	  else
	    nin = nin + 1
	    index_n(nin) = i
	  endif
	end do
	    

	! --- GET THE COEFFICIENTS DEPENDENT ON C ONLY
	c_w1(1,1) = abs(c(1))
	c_w1(1,2) = 0.0
	c_w1(2,1) = abs(c(2))
	c_w1(2,2) = 0.5*(1.-(1.-2.*c_w1(2,1))**2)
	c_w1(3,1) = abs(c(n))
	c_w1(3,2) = 0.5*(1.-(1.-2.*c_w1(3,1))**2)
	c_w1(4,1) = abs(c(n+1))
	c_w1(4,2) = 0.0

	rr1 = abs(c(2))
	rr2 = 1.-(rr1+rr1)
	r1  = rr2**2
	r2  = r1*rr2
	c_w2(1,1) = rr1*cc0
	c_w2(1,2) = (1.-r1)*cc1
	c_w2(1,3) = (1.-r2)*cc2

	rr1 = abs(c(3))
	rr2 = 1.-(rr1+rr1)
	r1  = rr2**2
	r2  = r1*rr2
	c_w2(2,1) = rr1*cc0
	c_w2(2,2) = (1.-r1)*cc1
	c_w2(2,3) = (1.-r2)*cc2

	rr1 = abs(c(n1))
	rr2 = 1.-(rr1+rr1)
	r1  = rr2**2
	r2  = r1*rr2
	c_w2(3,1) = rr1*cc0
	c_w2(3,2) = (1.-r1)*cc1
	c_w2(3,3) = (1.-r2)*cc2

	rr1 = abs(c(n))
	rr2 = 1.-(rr1+rr1)
	r1  = rr2**2
	r2  = r1*rr2
	c_w2(4,1) = rr1*cc0
	c_w2(4,2) = (1.-r1)*cc1
	c_w2(4,3) = (1.-r2)*cc2

	do i = 3 ,n1
	  rr1  = abs( c(i) )
	  rr2  = 1.0 - (rr1+rr1)
	  r1   = rr2*rr2
	  r2   = r1*rr2
	  r3   = r2*rr2
	  r4   = r3*rr2
 
	  c_w4(i,1)  = rr1 *c0
	  c_w4(i,2)  = (1.0-r1)*c1
	  c_w4(i,3)  = (1.0-r2)*c2
	  c_w4(i,4)  = (1.0-r3)*c3
	  c_w4(i,5)  = (1.0-r4)*c4
	end do

	return
	 end


!	=================================================
	subroutine pdadv2(c,q,c_w4,c_w2,c_w1,N,NOOS)
!	=================================================

! ---------------------------------------------------------
! --- This is a subroutine for the second part of Bott's 
! --- 		advection scheme.
! ---
! --- Andreas Bott  1989:  A Positive Definite Advection 
! --- 		scheme obtained by Nonlinear Renormalization
! --- 		of the advective fluxes
! --- 		Mon. Wea. Rev. 117 1006-15
! ---
! --- Fourth Order: with coefficients from    
! --- 		Mon. Wea. Rev. 117 2633-36
! ---
! --- Input: C=U*DT/DX[N+1]  &  Q[N]   Output: Q[2 N-1]
! ---   On the Staggered Grid:    C(i')----Q(i)----C(i'+1)
! ---
! ---   NOSS     = 1: Perform non-oscillatory option
! ---------------------------------------------------------

	USE pdadv_index
	
	real, parameter :: c0=1.0/1920.0,c1=1.0/384.00,c2=1.0/384.0
	real, parameter :: c3=1.0/768.00,c4=1.0/3840.0,ep=1.0e-15
	real, parameter :: cc0=-1./24.,cc1=1./16.,cc2=1./16.

	integer :: n, noos, i, j, k, n1, n2, n3
	real    :: c(n+1),q(n),c_w4(3:n-1,5),c_w2(4,3),c_w1(4,2),c_ww(n+1,5)
     	real    :: a0, a1, a2, a3, a4, q_l1, q_l2, q00, q_r1, q_r2
     	real    :: q_p1, q_p2, q_m1, q_m2, xxx1, xxx2
     
     	real	:: tmp6(6), tmp4(4)
        
	n1       = n-1
	n2       = n-2
	n3       = n-3
!
!     FOR ANY POSITIVE-DEFINITE Q ADVECTION
!

! --- 1. First order scheme for i=2 and n:
	a0         = q(1)
	a1         = q(2) - q(1)
	c_ww(1,1)  = a0
	c_ww(1,2)  = a0*c_w1(1,1)
	c_ww(2,3)  = a0*c_w1(2,1) + a1*c_w1(2,2)

! --- 2. Second order scheme for i=2,3,n1,n:
	a0         = q(3) - 26.*q(2) + q(1)
	a1         = q(3) - q(1)
	a2         = q(3) - 2.*q(2)  + q(1)
	c_ww(2,1)  = cc0*a0 + cc2*a2
	c_ww(2,2)  = a0*c_w2(1,1)-a1*c_w2(1,2)+a2*c_w2(1,3)
	c_ww(3,3)  = a0*c_w2(2,1)+a1*c_w2(2,2)+a2*c_w2(2,3) 
       
! --- 3. Fourth order scheme for i=3,n1:
	do i = 3 ,n2
	  q_l2 = q(i-2)
	  q_l1 = q(i-1)
	  q00  = q(i)
	  q_r1 = q(i+1)
	  q_r2 = q(i+2)
	  q_p1 = q_r1 + q_l1
	  q_p2 = q_r2 + q_l2
	  q_m1 = q_r1 - q_l1
	  q_m2 = q_r2 - q_l2
	  
	! --- coefficients: area preserving flux form
	  a0   = 9.0*q_p2 - 116.0*q_p1 + 2134.0*q00
	  a1   =-5.0*q_m2 +  34.0*q_m1
	  a2   =    -q_p2 +  12.0*q_p1 -   22.0*q00
	  a3   =     q_m2 -   2.0*q_m1
	  a4   =     q_p2 -   4.0*q_p1 +    6.0*q00
	  
	! --- integrals: for the use of in/out flux of the grid
!	  c_ww(i,1)   = c0*(a0+10.0*a2+a4)
	  c_ww(i,1)   = q00
	  c_ww(i,2)   = a0*c_w4(i,1)-a1*c_w4(i,2)+a2*c_w4(i,3)	&
                      - a3*c_w4(i,4)+a4*c_w4(i,5)
	  c_ww(i+1,3) = a0*c_w4(i+1,1)+a1*c_w4(i+1,2)	&
        	      + a2*c_w4(i+1,3)  		&
        	      + a3*c_w4(i+1,4)+a4*c_w4(i+1,5)
	end do

! --- 4. Second order scheme for i=2,3,n1,n:
	a0         = q(n) - 26.*q(n1) + q(n2)
	a1         = q(n) - q(n2) 
	a2         = q(n) - 2.*q(n1)  + q(n2)
	c_ww(n1,1) = cc0*a0 + cc2*a2
	c_ww(n1,2) = a0*c_w2(3,1)-a1*c_w2(3,2)+a2*c_w2(3,3)
	c_ww(n,3)  = a0*c_w2(4,1)+a1*c_w2(4,2)+a2*c_w2(4,3)

! --- 5. First order scheme for i=2 and n:
	a0         = q(n)
	a1         = q(n)-q(n1)
	c_ww(n,1)  = a0
	c_ww(n,2)  = a0*c_w1(3,1) - a1*c_w1(3,2)
	c_ww(n+1,3)= a0*c_w1(4,1)

! --- restrict the integrals to preserve the sign
	do i=1,nip
	  c_ww(index_p(i),2) = 0.0
	  if(c_ww(index_p(i),3).lt.0.0) c_ww(index_p(i),3) = 0.0
	end do
	do i=1,nin
	  if(c_ww(index_n(i),2).lt.0.0) c_ww(index_n(i),2) = 0.0
	  c_ww(index_n(i),3) = 0.0
	end do

	do i = 1 ,n
	  xxx1  = c_ww(i,2)+c_ww(i+1,3)+ep
	  if(c_ww(i,1).lt.xxx1) c_ww(i,1) = xxx1
	end do

! --- get the weighting factor
	do i = 1 ,n
	  c_ww(i,1)   = q(i) / c_ww(i,1)
	end do

!                                         <= c_ww(i,2)
! --- get the in/out flux of the grid  i --- i+1/2
!                                            c_ww(i,3) =>
	do i = 1 ,n
	  c_ww(i,2)   = c_ww(i,2)*c_ww(i,1)
	  c_ww(i+1,3) = c_ww(i+1,3)*c_ww(i,1)
	end do

! ===
! === 050198 add boundary value c_ww(1,3) and c_ww(n+1,2)
! ===	in-flow boundary set gradient of fluxes = 0
! ===	out-flow boundary no inflow flux
! ===
	if(c(1).gt.0.0)then
	  c_ww(1,3) = c_ww(2,3) - c_ww(2,2) + c_ww(1,2)
	else
	  c_ww(1,3) = 0.0
	endif

	if(c(n+1).lt.0.0)then
	  c_ww(n+1,2) = c_ww(n,2) - c_ww(n,3) + c_ww(n+1,3)
	else
	  c_ww(n+1,2) = 0.0
	endif
	
      if( noos.ne.1 ) then
! --- compute the total advection tendency

	do i = 2 ,n1
	  q(i)     = c_ww(i+1,2)-c_ww(i+1,3) - c_ww(i,2)+c_ww(i,3)+q(i)
	end do
!
      else
! ---
! --- NON-OSCILLATORY OPTION: Zalesak FCT LIMITER
! --- P.K.Smolarkiewicz  &  W.W.Grabowski, 1990: The multidimensional
! --- positive definite advection transport algorithm: Nonoscillatory
! --- option,  J. Comput. Phys., 86, 355-375
! ---

	! --- Get the donor-cell fluxes (low-order)
	do i = 2 ,n
	  if( c(i).gt.0.0 ) then
		c_ww(i,1) = q(i-1)
	  else
		c_ww(i,1) =-q(i)
	  endif
	end do

	c_ww(1,1) = abs(q(1)*c(1))	
	do i = 2 ,n
	  c_ww(i,1)  = c_ww(i,1) * c(i)
	end do
	c_ww(n+1,1) = abs(q(n)*c(n+1))
	
	do i=1,nip
	  c_ww(index_p(i),4) = 0.0
	  c_ww(index_p(i),5) = c_ww(index_p(i),1)
	end do
	do i=1,nin
	  c_ww(index_n(i),4) = c_ww(index_n(i),1)
	  c_ww(index_n(i),5) = 0.0
	end do

	do i = 1 ,n
	  c_ww(i,1) = c_ww(i+1,4) - c_ww(i+1,5)		&
        	    - c_ww(i,4) + c_ww(i,5)		&
        	    + q(i)
	end do

	! --- get the a-flux = f(high-order)-f(low-order)
	do i = 1 ,n + 1
	  c_ww(i,4)   = c_ww(i,2) - c_ww(i,4)                     
	  c_ww(i,5)   = c_ww(i,3) - c_ww(i,5)
	  
	  xxx1        = c_ww(i,4)
	  if(xxx1.lt.0.0) xxx1 = 0.0
	  xxx2        = c_ww(i,5)
	  if(xxx2.gt.0.0) xxx2 = 0.0
	  c_ww(i,2)   = xxx1 - xxx2
	  
	  xxx1        = c_ww(i,5)
	  if(xxx1.lt.0.0) xxx1 = 0.0
	  xxx2        = c_ww(i,4)
	  if(xxx2.gt.0.0) xxx2 = 0.0
	  c_ww(i,3)   = xxx1 - xxx2	  
	end do

	c_ww(1,4)= min(c_ww(1,1),c_ww(2,1),q(1),q(2))
	c_ww(1,5)= max(c_ww(1,1),c_ww(2,1),q(1),q(2))
	do i = 2 ,n1
	   j = i-1
	   k = i+1
	  c_ww(i,4)= min(c_ww(j,1),c_ww(i,1),c_ww(k,1),q(j),q(i),q(k)) 
	  c_ww(i,5)= max(c_ww(j,1),c_ww(i,1),c_ww(k,1),q(j),q(i),q(k))
	end do
	c_ww(n,4) = min(c_ww(n1,1),c_ww(n,1),q(n1),q(n))
	c_ww(n,5) = max(c_ww(n1,1),c_ww(n,1),q(n1),q(n))

	do i = 1 ,n
	  c_ww(i,4) = (c_ww(i,1)-c_ww(i,4)) / (c_ww(i,2)+c_ww(i+1,3)+ep) 
	  c_ww(i,5) = (c_ww(i,5)-c_ww(i,1)) / (c_ww(i,3)+c_ww(i+1,2)+ep)  
	  q(i)      =  c_ww(i,1)
	end do

	do i = 2 ,n
	  c_ww(i,1) = min(1.0,c_ww(i-1,5),c_ww(i,4) )
	  c_ww(i,2) = c_ww(i,2) * c_ww(i,1)                     
	  c_ww(i,1) = min(1.0,c_ww(i-1,4),c_ww(i,5) )
	  c_ww(i,3) = c_ww(i,3) * c_ww(i,1)
	end do

! ---
! --- compute the high-order advection tendency
! --- then compute the total advection tendency
! ---
	do i = 2 ,n1
	  c_ww(i,1)   = c_ww(i+1,2) - c_ww(i+1,3)	&
      	              - c_ww(i,2)   + c_ww(i,3)
     	  q(i) = c_ww(i,1) + q(i)
	end do

      endif
 
	return
	 end


! ================================================================
!
!  PDADV.F90:
!	subroutine pdadv1(c,c_w4,c_w2,c_w1,N)
!	subroutine pdadv2(c,q,c_w4,c_w2,c_w1,c_ww,N,NOOS)                
!
!  Purpose:
!	A package of the improved Bott's advection scheme based on
!		Chien Wang and Xingzhong Liang's work.
!	Modified by Julien Savre to account for non constant grid size
!	based on the modified flux form proposed by Easter (MWR, 1993).			  
!
!  Author
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
!  Revision:
!	Date	By		Brief Description
!	----	--		-----------------	
!	070698	Chien Wang	cpp
!	072098	Chien Wang	F90/F95 version
!	072098	Chien Wang	revised
!	012400	Chien Wang	added index_p,index_n
!	021903	Chien Wang	rev. for -r8
!	081607	Chien Wang	f90 free format
!       310113  Julien Savre    Copied from CRM-MIT-ice (non optimized)
!       050213  Julien Savre    New version entirely rewritten
!
! ================================================================
! $Id: pdadv.f90 561 2013-04-15 15:52:36Z x_julsa $
! ================================================================

!	=====================================
 	subroutine pdadv1_mod (c,c_w1,c_w2,c_w3,c_w4,dtx,n)
!	=====================================

! -------------------------------------------------------------
! --- This is a subroutine for the first part 
! --- 	of Bott's advection scheme.
! ---
! --- Andreas Bott  1989:  A Positive Definite 
! --- 		Advection scheme obtained by Nonlinear
! --- 		Renormalization of the advective fluxes
! --- 		Mon. Wea. Rev. 117 1006-15
! ---
! --- Fourth Order: with coefficients from    
! --- 		Mon. Wea. Rev. 117 2633-36
! ---
! --- Input: C=U*DT/DX[N+1]  Output: W4[3:N1,5],W2[2;3;n1;n,3]
! ---				 and   W1[1;2;n;n+1,2]
! ---   On the Staggered Grid: C(i')----Q(i)----C(i'+1)
! -------------------------------------------------------------
	
	IMPLICIT NONE

	integer :: N, n1, n2, n3, i 
	real    :: dtx(n), c(n+1)
	real    :: c_w1(n,5), c_w2(n,5), c_w3(n,5), c_w4(n,5)
     	real    :: rr1, rr2, r1, r2, r3, r4
!
	n1 = n-1
	n2 = n-2
	n3 = n-3
!
!  Interior points
!
	do i = 2 ,n-1
	! cp
	  rr1  = 0.5*(c(i) + abs(c(i)))*dtx(i)
	  rr2  = 1.0 - (rr1+rr1)
	  r1   = rr2*rr2
	  r2   = r1*rr2
	  r3   = r2*rr2
	  r4   = r3*rr2
 
	  c_w1(i,1)  = rr1+rr1 
	  c_w1(i,2)  = (1.0-r1)
	  c_w1(i,3)  = (1.0-r2)
	  c_w1(i,4)  = (1.0-r3)
	  c_w1(i,5)  = (1.0-r4)

	! cp0
	  rr1  = 0.5*(c(i-1) + abs(c(i-1)))*dtx(i)
	  rr2  = 1.0 - (rr1+rr1)
	  r1   = rr2*rr2
	  r2   = r1*rr2
	  r3   = r2*rr2
	  r4   = r3*rr2
 
	  c_w2(i,1)  = rr1+rr1
	  c_w2(i,2)  = (1.0-r1)
	  c_w2(i,3)  = (1.0-r2)
	  c_w2(i,4)  = (1.0-r3)
	  c_w2(i,5)  = (1.0-r4)

	! cm
	  rr1  = -0.5*(c(i) - abs(c(i)))*dtx(i)
	  rr2  = 1.0 - (rr1+rr1)
	  r1   = rr2*rr2
	  r2   = r1*rr2
	  r3   = r2*rr2
	  r4   = r3*rr2
 
	  c_w3(i,1)  = rr1+rr1
	  c_w3(i,2)  = (1.0-r1)
	  c_w3(i,3)  = (1.0-r2)
	  c_w3(i,4)  = (1.0-r3)
	  c_w3(i,5)  = (1.0-r4)

	! cm0
	  rr1  = -0.5*(c(i-1) - abs(c(i-1)))*dtx(i)
	  rr2  = 1.0 - (rr1+rr1)
	  r1   = rr2*rr2
	  r2   = r1*rr2
	  r3   = r2*rr2
	  r4   = r3*rr2
 
	  c_w4(i,1)  = rr1+rr1
	  c_w4(i,2)  = (1.0-r1)
	  c_w4(i,3)  = (1.0-r2)
	  c_w4(i,4)  = (1.0-r3)
	  c_w4(i,5)  = (1.0-r4)
	end do
			
	return
	 end


!	=================================================
	subroutine pdadv2_mod(c,q,dtx,c_w1,c_w2,c_w3,c_w4,N,NOOS)
!	=================================================

! ---------------------------------------------------------
! --- This is a subroutine for the second part of Bott's 
! --- 		advection scheme.
! ---
! --- Andreas Bott  1989:  A Positive Definite Advection 
! --- 		scheme obtained by Nonlinear Renormalization
! --- 		of the advective fluxes
! --- 		Mon. Wea. Rev. 117 1006-15
! ---
! --- Fourth Order: with coefficients from    
! --- 		Mon. Wea. Rev. 117 2633-36
! ---
! --- Input: C=U*DT/DX[N+1]  &  Q[N]   Output: Q[2 N-1]
! ---   On the Staggered Grid:    C(i')----Q(i)----C(i'+1)
! ---
! ---   NOSS     = 1: Perform non-oscillatory option
! ---------------------------------------------------------
	
	IMPLICIT NONE

	real, parameter :: c0=1.0/3840.0,c1=1.0/384.0,c2=1.0/1152.0	
	real, parameter :: c3=1.0/768.0,c4=1.0/3840.0,ep=1.0e-12
	real, parameter :: cc0=1./48.,cc1=1./16.,cc2=1./48.
	real, parameter :: ccc0=1./2., ccc1=1./8.
	real, parameter :: i0=1./1920., i2=1./576., i4=1./1920.
	real, parameter :: ii0=1./24., ii2=1./24.

	integer :: n, noos, i, n1, n2, n3, j, k, h
	real    :: dtx(n), c(n+1), q(n)
	real    :: c_w1(n,5), c_w2(n,5), c_w3(n,5), c_w4(n,5)
	real    :: a(n,5), am(n,5), ap(n,5), I00(n)
	real    :: frp(n), frm(n), flp(n), flm(n), fr(n), fl(n)
	real    :: frup(n), flup(n), qup(n), qmax(n), qmin(n) 
	real    :: ar(n), al(n), rp(n), rm(n)
	real    :: pp, pm, cc, cco, a0, qim1, cim1, iim1
!
!-----------------------------------------------------------------!
!    		      Advection of scalars                        !
!-----------------------------------------------------------------!
!
	frp = 0.0
	flp = 0.0
	frm = 0.0
	flm = 0.0
	fr  = 0.0
	fl  = 0.0
	I00 = 0.0
!
!-----------------------------------------------------------------!
!    		     Calculate coefficients                       !
!-----------------------------------------------------------------!
!
	a(1,1) = q(1)
	a(1,2) = q(2) - q(1)
	a(1,3) = 0.
	a(1,4) = 0.
	a(1,5) = 0.
	I00(1) = a(1,1)
!
	a(2,1) = -q(3) + 26.*q(2) - q(1)
	a(2,2) = q(3) - q(1)
	a(2,3) = q(3) - 2.*q(2) + q(1)
	a(2,4) = 0.
	a(2,5) = 0.
	I00(2) = ii0*a(2,1) + ii2*a(2,3)
!
	do i = 3, n-2
	  a(i,1) = 9.*q(i-2) - 116.*q(i-1) + 2134.*q(i) - 116.*q(i+1) + 9.*q(i+2)
	  a(i,2) = 5.*q(i-2) - 34.*q(i-1) + 34.*q(i+1) - 5.*q(i+2)
	  a(i,3) = -3.*q(i-2) + 36.*q(i-1) - 66.*q(i) + 36.*q(i+1) - 3.*q(i+2)
	  a(i,4) = -q(i-2) + 2.*q(i-1) - 2.*q(i+1) + q(i+2)
	  a(i,5) = q(i-2) - 4.*q(i-1) + 6.*q(i) - 4.*q(i+1) + q(i+2)
	  I00(i) = i0*a(i,1) + i2*a(i,3) + i4*a(i,5)
	enddo
!
	a(n-1,1) = -q(n) + 26.*q(n-1) - q(n-2)
	a(n-1,2) = q(n) - q(n-2)
	a(n-1,3) = q(n) - 2.*q(n-1) + q(n-2)
	a(n-1,4) = 0.
	a(n-1,5) = 0.
	I00(n-1) = ii0*a(n-1,1) + ii2*a(n-1,3)
!
	a(n,1) = q(n)
	a(n,2) = q(n) - q(n-1)
	a(n,3) = 0.
	a(n,4) = 0.
	a(n,5) = 0.
	I00(n) = a(n,1)
!
	do h = 1, 5
	  a0 = a(1,h)
	  do i = 2, n-1
	    am(i,h)= a0
	    a0 = a(i,h)
	  enddo
!
	  a0 = a(n,h)
	  do i = n-1, 2, -1
	    ap(i,h)= a0
	    a0 = a(i,h)
	  enddo
	enddo
!
!-----------------------------------------------------------------!
!    		        Boundary fluxes                           !
!-----------------------------------------------------------------!
!
!  Point 2: needs second order adn 1st order
!
	frp(2) = CC0*c_w1(2,1)*a(2,1) + CC1*c_w1(2,2)*a(2,2) + CC2*c_w1(2,3)*a(2,3)
	flp(2) = CCC0*c_w2(2,1)*a(1,1) + CCC1*c_w2(2,2)*a(1,2)
!
	frm(2) = C0*c_w3(2,1)*a(3,1) - C1*c_w3(2,2)*a(3,2) + C2*c_w3(2,3)*a(3,3)			   &
	       - C3*c_w3(2,4)*a(3,4) + C2*c_w3(2,5)*a(3,5)
	flm(2) = CC0*c_w4(2,1)*a(2,1) - CC1*c_w4(2,2)*a(2,2) + CC2*c_w4(2,3)*a(2,3)
!
!  Point 3: needs second order
!
	frp(3) = C0*c_w1(3,1)*a(3,1) + C1*c_w1(3,2)*a(3,2) + C2*c_w1(3,3)*a(3,3)	       		   &
	       + C3*c_w1(3,4)*a(3,4) + C4*c_w1(3,5)*a(3,5)
	flp(3) = CC0*c_w2(3,1)*a(2,1) + CC1*c_w2(3,2)*a(2,2) + CC2*c_w2(3,3)*a(2,3)
! 
	frm(3) = C0*c_w3(3,1)*a(4,1) - C1*c_w3(3,2)*a(4,2) + C2*c_w3(3,3)*a(4,3)	       	  	   &
	       - C3*c_w3(3,4)*a(4,4) + C4*c_w3(3,5)*a(4,5)
	flm(3) = C0*c_w4(3,1)*a(3,1) - C1*c_w4(3,2)*a(3,2) + C2*c_w4(3,3)*a(3,3)			   &
	       - C3*c_w4(3,4)*a(3,4) + C4*c_w4(3,5)*a(3,5)
!
!  point n-2: needs 2nd order
!
	frp(n-2) = C0*c_w1(n-2,1)*a(n-2,1) + C1*c_w1(n-2,2)*a(n-2,2) + C2*c_w1(n-2,3)*a(n-2,3)		   &
	         + C3*c_w1(n-2,4)*a(n-2,4) + C4*c_w1(n-2,5)*a(n-2,5)
	flp(n-2) = C0*c_w2(n-2,1)*a(n-3,1) + C1*c_w2(n-2,2)*a(n-3,2) + C2*c_w2(n-2,3)*a(n-3,3) 		   &
	         + C3*c_w2(n-2,4)*a(n-3,4) + C4*c_w2(n-2,5)*a(n-3,5)
!
	frm(n-2) = CC0*c_w3(n-2,1)*a(n-1,1) - CC1*c_w3(n-2,2)*a(n-1,2) + CC2*c_w3(n-2,3)*a(n-1,3)
	flm(n-2) = C0*c_w4(n-2,1)*a(n-2,1) - C1*c_w4(n-2,2)*a(n-2,2) + C2*c_w4(n-2,3)*a(n-2,3)		   &
		 - C3*c_w4(n-2,4)*a(n-2,4) + C4*c_w4(n-2,5)*a(n-2,5)
!
!  point n-1: needs 2nd order and 1st order
!
	frp(n-1) = CC0*c_w1(n-1,1)*a(n-1,1) + CC1*c_w1(n-1,2)*a(n-1,2) + CC2*c_w1(n-1,3)*a(n-1,3)
	flp(n-1) = C0*c_w2(n-1,1)*a(n-2,1) + C1*c_w2(n-1,2)*a(n-2,2) + C2*c_w2(n-1,3)*a(n-2,3) 		   &
	         + C3*c_w2(n-1,4)*a(n-2,4) + C4*c_w2(n-1,5)*a(n-2,5)
!
	frm(n-1) = CCC0*c_w3(n-1,1)*a(n,1) - CCC1*c_w3(n-1,2)*a(n,2) 
	flm(n-1) = CC0*c_w4(n-1,1)*a(n-1,1) - CC1*c_w4(n-1,2)*a(n-1,2) + CC2*c_w4(n-1,3)*a(n-1,3)
!
!-----------------------------------------------------------------!
!               Calculate other fluxes : 4th order                !
!-----------------------------------------------------------------!
!
	do i = 4, n-3
	  frp(i) = C0*c_w1(i,1)*a(i,1) + C1*c_w1(i,2)*a(i,2) + C2*c_w1(i,3)*a(i,3)		&
	         + C3*c_w1(i,4)*a(i,4) + C4*c_w1(i,5)*a(i,5)
	  flp(i) = C0*c_w2(i,1)*am(i,1) + C1*c_w2(i,2)*am(i,2) + C2*c_w2(i,3)*am(i,3)		&
	         + C3*c_w2(i,4)*am(i,4) + C4*c_w2(i,5)*am(i,5)
!
	  frm(i) = C0*c_w3(i,1)*ap(i,1) - C1*c_w3(i,2)*ap(i,2) + C2*c_w3(i,3)*ap(i,3)		&
	         - C3*c_w3(i,4)*ap(i,4) + C4*c_w3(i,5)*ap(i,5)
	  flm(i) = C0*c_w4(i,1)*a(i,1) - C1*c_w4(i,2)*a(i,2) + C2*c_w4(i,3)*a(i,3)		&
	         - C3*c_w4(i,4)*a(i,4) + C4*c_w4(i,5)*a(i,5)
	enddo
!
!  Limitations
!
	flp(2) = 0.
	flm(2) = 0.
	do i = 2, n-1
	  frp(i) = max(frp(i),0.0)
	  frm(i) = max(frm(i),0.0)
	  flp(i) = max(flp(i),0.0)
	  flm(i) = max(flm(i),0.0)
	  I00(i) = max(I00(i),frp(i)+flm(i),ep)
	enddo
!
!-----------------------------------------------------------------!
!    		    Cumulate high order fluxes                    !
!-----------------------------------------------------------------!
!
	fr(2) = max(c(2),0.)*frp(2)*q(2)/I00(2) - max(-c(2),0.)*frm(2)*q(3)/I00(3)
	fl(2) = max(c(1),0.)*flp(2) - max(-c(1),0.)*flm(2)*q(2)/I00(2)
	cim1 = c(2)
	qim1 = q(2)
	iim1 = I00(2)
	do i = 3, n-2
	  fr(i) = max(c(i),0.)*frp(i)*q(i)/I00(i) - max(-c(i),0.)*frm(i)*q(i+1)/I00(i+1)
	  fl(i) = max(cim1,0.)*flp(i)*qim1/iim1 - max(-cim1,0.)*flm(i)*q(i)/I00(i)
	  cim1 = c(i)
	  qim1 = q(i)
	  iim1 = I00(i)
	enddo
	fr(n-1) = max(c(n-1),0.)*frp(n-1)*q(n-1)/I00(n-1) - max(-c(n-1),0.)*frm(n-1)
	fl(n-1) = max(c(n-2),0.)*flp(n-1)*q(n-2)/I00(n-2) - max(-c(n-2),0.)*flm(n-1)*q(n-1)/I00(n-1)
!
!-----------------------------------------------------------------!
!    		 Update q without flux limitation                 !
!-----------------------------------------------------------------!
!
      if( NOOS /= 1 ) then
	do i = 2, n-1
	  q(i) = q(i) - (fr(i) - fl(i)) * dtx(i)
	end do
!
!-----------------------------------------------------------------!
!    	    Flux limited tendencies: FCT method, Zalesak          !
!-----------------------------------------------------------------!
!
      else
        ar = 0.0
	al = 0.0
        pp = 0.0
	pm = 0.0
        rp = 0.0
	rm = 0.0
	qup  = q
	qmin = q
	qmax = q
!
!-----------------------------------------------------------------!
!          Calculate low order fluxes : upstream method           !
!-----------------------------------------------------------------!
!
	cim1 = c(1)
	qim1 = q(1) 
	do i = 2, n-1
	  frup(i) = 0.5*(c(i)*(q(i+1) + q(i)) - abs(c(i))*(q(i+1) - q(i)))
	  flup(i) = 0.5*(cim1*(q(i) + qim1) - abs(cim1)*(q(i) - qim1))
!
	  qup(i) = q(i) - (frup(i) - flup(i)) * dtx(i)
	  cim1 = c(i)
	  qim1 = q(i)
	enddo	
!
!-----------------------------------------------------------------!
!                   Get antidiffusive fluxes                      !
!-----------------------------------------------------------------!
!
	qmax(2) = max(q(2),q(3),qup(3),qup(2)) 
	qmin(2) = min(q(2),q(3),qup(3),qup(2)) 
	do i = 3, n-2
	  j = i+1
	  k = i-1
	  qmax(i) = max(q(i),q(j),q(k),qup(i),qup(j),qup(k)) 
	  qmin(i) = min(q(i),q(j),q(k),qup(i),qup(j),qup(k)) 
	enddo
	qmax(n-1) = max(q(n-2),q(n-1),qup(n-2),qup(n-1)) 
	qmin(n-1) = min(q(n-2),q(n-1),qup(n-2),qup(n-1)) 
!
	do i = 2, n-1
	  ar(i) = fr(i) - frup(i)
	  al(i) = fl(i) - flup(i)
!
	  pp = max(0.0,al(i)) - min(0.0,ar(i))
	  pm = max(0.0,ar(i)) - min(0.0,al(i))
!
	  if (pp > 0.0) rp(i) = min(1.0,(qmax(i) - qup(i)) / pp)
	  if (pm > 0.0) rm(i) = min(1.0,(qup(i) - qmin(i)) / pm)
	enddo
!
!-----------------------------------------------------------------!
!                      Get correction terms                       !
!-----------------------------------------------------------------!
!	
	cco = 0.
	do i = 2, n-2
	  if (ar(i) >= 0.0) then
	    cc = min(rp(i+1), rm(i))
	  else
	    cc = min(rp(i), rm(i+1))
	  endif
	  ar(i) = cc*ar(i)
	  al(i) = cco*al(i)
	  cco = cc
	enddo
!
!-----------------------------------------------------------------!
!                           Update q                              !
!-----------------------------------------------------------------!
!	
	do i = 3, n-2
	  q(i) = q(i) - (ar(i) - al(i)) * dtx(i)
	enddo
      endif
!  
!-----------------------------------------------------------------!
! 
return
end
