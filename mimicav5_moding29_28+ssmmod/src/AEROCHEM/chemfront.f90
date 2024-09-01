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
!  CHEM_M_FRONT.F:
!	subroutine chem_m_front(nod, per_creact, dtr)
!
!  Purpose:
!	Front gate for chemistry modules.
!
!  Author:
!	Chien Wang
!	MIT Joint Program on Science and Policy of Global Change
!
! ================================================================

!	===================================================
	subroutine chem_m_front                     &
     		       ( ipstart, ipend, jpstart, jpend,  &
		         nod, per_creact,	    &
     			 gas0                       &
#ifdef CHEM_ENABLE
     			, gas	                    &
#endif
#ifdef AQCHEM_ENABLE 
     			, aqc,  aqr                 &
#endif
#ifdef SOLIDCHEM_ENABLE
     			, solidi                    &
#endif
					)
!	===================================================
!
! -------------------------------------------------
!  Brief Description of Dummy Variables:
!
!	dt:	time step of integration in second
!	ifss:	1 for using steady-state assumption
!		0 for not using
!	ktrop:  layer label of tropopause
!	Temp:	temperature in Kelvin
!	p:	air pressure in hPa
!	qv:	water vapor density ratio in kg/m^3
!	den:	air density in kg/m^3
!	uuu0:	cos(zenith angle)
!
! -------------------------------------------------	
!
#ifdef SPMD
USE mpi
#endif
USE gridno
USE typedef_gas
USE typedef_aq
USE typedef_solid
USE shared_data
USE shared_state
USE shared_wind
USE shared_turbu
USE shared_nuclei
USE shared_thermo
USE shared_rad
USE shared_hydro

	integer :: nod, ifss, i, j, k, ierr, ipstart, ipend, jpstart, jpend
	real    :: qc_max, qr_max, qi_max
	logical :: per_creact

#ifdef CHEM_ENABLE
	real, dimension(nz)   :: ad_o3,ad_co,ad_nox,ad_h2o2,ad_ch2o,  &
     				ad_ch3o2h,ad_so2,ad_svi,ad_nv,	      &
     				ad_nh4     
       real, dimension(nz)   :: ed_o3,ed_co,ed_nox,ed_h2o2,ed_ch2o,   &
     				ed_ch3o2h,ed_nv,ed_so2,ed_svi,	      &
     				ed_nh4
       common /diagtmp/ad_o3,ad_co,ad_nox,ad_h2o2,ad_ch2o,	      &
     		       ad_ch3o2h,ad_so2,ad_svi,ad_nv,		      &
     		       ad_nh4,					      &
     		       ed_o3,ed_co,ed_nox,ed_h2o2,ed_ch2o,	      &
     		       ed_ch3o2h,ed_nv,ed_so2,ed_svi,		      &
     		       ed_nh4
#endif

	type (gas_chemical), dimension(nz) :: gas0
	
#include "chemdef.h"

	! --- note: currently 10 diagnostic variables
	integer				:: mnv
	real, dimension (nz,10)		:: tdiag, t2diag
	real, dimension (nz,10)		:: gdiag, g2diag
	real, dimension (nz)		:: wdiag

	ifss = 0
!
! =================================================================         	
!
#ifdef CHEM_ENABLE
!
!----------------------------------------------------------!
!   			 Diagnostics			   !
!----------------------------------------------------------!
!
  if (per_creact) then
    tdiag(1:nz, 1) = ad_o3     (1:nz)
    tdiag(1:nz, 2) = ad_co     (1:nz)
    tdiag(1:nz, 3) = ad_h2o2   (1:nz)
    tdiag(1:nz, 4) = ad_nox    (1:nz)
    tdiag(1:nz, 5) = ad_ch2o   (1:nz)
    tdiag(1:nz, 6) = ad_ch3o2h (1:nz)
    tdiag(1:nz, 7) = ad_so2    (1:nz)
    tdiag(1:nz, 8) = ad_svi    (1:nz)
    tdiag(1:nz, 9) = ad_nv     (1:nz)
    tdiag(1:nz,10) = ad_nh4    (1:nz)
!
    t2diag(1:nz, 1) = ed_o3	(1:nz)
    t2diag(1:nz, 2) = ed_co	(1:nz)
    t2diag(1:nz, 3) = ed_h2o2	(1:nz)
    t2diag(1:nz, 4) = ed_nox	(1:nz)
    t2diag(1:nz, 5) = ed_ch2o	(1:nz)
    t2diag(1:nz, 6) = ed_ch3o2h (1:nz)
    t2diag(1:nz, 7) = ed_so2	(1:nz)
    t2diag(1:nz, 8) = ed_svi	(1:nz)
    t2diag(1:nz, 9) = ed_nv	(1:nz)
    t2diag(1:nz,10) = ed_nh4	(1:nz)
!
!  Reduce all procs.
!
#if ( defined SPMD )
    call MPI_BARRIER ( MPI_COMM_WORLD, ierr )
    bufsize = nz*10
    call MPI_REDUCE ( tdiag(1,1),gdiag(1,1),bufsize,MPI_REAL,	  &
    		      MPI_SUM,0, MPI_COMM_WORLD, ierr )
!
    call MPI_REDUCE ( t2diag(1,1),g2diag(1,1),bufsize,MPI_REAL,   &
    		      MPI_SUM,0, MPI_COMM_WORLD, ierr )
#else
    gdiag (:,:) = tdiag (:,:)
    g2diag(:,:) = t2diag(:,:)
#endif
!
!  Gaseous diagnostics
!
    if ( mypid .eq. 0 ) then	    
      do mnv = 1,10
    	wdiag(1:nz) = gdiag (1:nz,mnv)        
    	write(151) wdiag
    	wdiag(1:nz) = g2diag(1:nz,mnv)        
    	write(152) wdiag
      end do
    end if
!		    
    ad_o3(:)	 = 0.0
    ad_co(:)	 = 0.0
    ad_h2o2(:)   = 0.0
    ad_nox(:)	 = 0.0
    ad_ch2o(:)   = 0.0
    ad_ch3o2h(:) = 0.0
    ad_so2(:)	 = 0.0
    ad_svi(:)	 = 0.0
    ad_nv(:)	 = 0.0
    ad_nh4(:)	 = 0.0
!
    ed_o3(:)	 = 0.0
    ed_co(:)	 = 0.0
    ed_h2o2(:)   = 0.0
    ed_nox(:)	 = 0.0
    ed_ch2o(:)   = 0.0
    ed_ch3o2h(:) = 0.0
    ed_so2(:)	 = 0.0
    ed_svi(:)	 = 0.0
    ed_nv(:)	 = 0.0
    ed_nh4(:)	 = 0.0
!
!----------------------------------------------------------!
!   		  Tropospheric chemistry		   !
!----------------------------------------------------------!
!
#ifdef CHEM_REACT
    call chemtrop( ifss, uuu0, gas0, gas )
#endif
!
!----------------------------------------------------------!
!   		         Limitations		           !
!----------------------------------------------------------!
!
#ifdef AQCHEM_ENABLE
   if (lmicro > 0) then
   qc_max = 0.0
   qr_max = 0.0
    do k=1,nz
      do j=jt_start,jt_end
        do i=it_start,it_end
 	  if ( hydromtr2(i,j,k,drop)%q > qmax(drop) ) then
 	    qmax(drop) = hydromtr2(i,j,k,drop)%q
	  endif
	  if ( hydromtr2(i,j,k,rain)%q > qmax(rain) ) then
 	    qmax(rain) = hydromtr2(i,j,k,rain)%q
          endif
        end do
      end do
    end do
!
    if (qmax(drop).gt.qmin(drop) .or. qmax(rain).gt.qmin(rain)) then
      call aqtrop( 1		    &
#if (defined AQCHEM_ENABLE)
      		, gas, aqc, aqr     &
#endif
				  )
    endif
#endif
!
#ifdef SOLIDCHEM_ENABLE
    if (lmicro > 1) then
    qi_max = 0.0
    do k=1,nz
      do j=jt_start,jt_end
        do i=it_start,it_end
          if ( hydromtr2(i,j,k,ice)%q > qmax(ice) ) then
            qmax(ice) = hydromtr2(i,j,k,ice)%q
          endif
        end do
      end do
    end do

    if (qmax(ice) > qmin(ice)) then
      call solidtrop( 0, gas, solidi )
    endif
#endif
!
  endif	
!
#endif  	
!
!----------------------------------------------------------!
!	
  return
  end
!
!
!    ==========================================                                                     
	subroutine chem_m_init 
!    ==========================================
!
USE gridno
USE shared_data
USE shared_gas
USE typedef_aq
USE typedef_solid
!
! ---
! --- Units         :  1.e-9 kg/kg
! ---
!
	integer :: i, j, k
	real    :: z, hgt, dens
!        
	real, dimension(nz,17) :: ptr
!
#ifdef CHEM_ENABLE	
       real, dimension(nz)   :: ad_o3,ad_co,ad_nox,ad_h2o2,ad_ch2o,    &
     				ad_ch3o2h,ad_so2,ad_svi,ad_nv
!
       real, dimension(nz)   :: ed_o3,ed_co,ed_nox,ed_h2o2,ed_ch2o,    &
     				ed_ch3o2h,ed_nv,ed_so2,ed_svi
!
       common /diagtmp/ad_o3,ad_co,ad_nox,ad_h2o2,ad_ch2o,	       &
     		       ad_ch3o2h,ad_so2,ad_svi,ad_nv,		       &
     		       ed_o3,ed_co,ed_nox,ed_h2o2,ed_ch2o,	       &
     		       ed_ch3o2h,ed_nv,ed_so2,ed_svi
!
!----------------------------------------------------------!
!   	        Initialize gas0 from file	           !
!----------------------------------------------------------!
!
  do k = 1, nz
    read(105) ptr(k, 1),ptr(k, 2),ptr(k, 3),ptr(k, 4),ptr(k, 5),    & 
	      ptr(k, 6),ptr(k, 7),ptr(k, 8),ptr(k, 9),ptr(k,10),    &
	      ptr(k,11),ptr(k,12),ptr(k,13),ptr(k,14),ptr(k,15),    &
	      ptr(k,16),ptr(k,17)
  enddo
!
  gas0(:)%o3     = ptr(:,1)
  gas0(:)%co     = ptr(:,2)
  gas0(:)%zco2   = ptr(:,3)
  gas0(:)%ho     = ptr(:,4)
  gas0(:)%ho2    = ptr(:,5)
  gas0(:)%h2o2   = ptr(:,6)
  gas0(:)%xno    = ptr(:,7)
  gas0(:)%xno2   = ptr(:,8)
  gas0(:)%xno3   = ptr(:,9)
  gas0(:)%xn2o5  = ptr(:,10)
  gas0(:)%hno3   = ptr(:,11)
  gas0(:)%ch4    = ptr(:,12)
  gas0(:)%ch2o   = ptr(:,13)
  gas0(:)%ch3o2h = ptr(:,14)
  gas0(:)%so2    = ptr(:,15)
  gas0(:)%h2so4  = ptr(:,16)
  gas0(:)%dms    = ptr(:,17)
!
  close(105)
!
#ifdef AMMONIA
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        hgt = (den0(k)/den0(1))**3.
        if(k.le.ktrop) then
          gas0(k)%nh3 = 0.15/den0(1) * hgt  ! ppb
        else
          gas0(k)%nh3 = 0.0
        endif
      enddo
    enddo
  enddo
#endif  
!
!----------------------------------------------------------!
!   		 Output diagnostics header	           !
!----------------------------------------------------------!
!
  write(7,660)
660	format(1x/,8x,60('=')/)
  write(7,2022)
2022	format(                                                     &
     	    6x,2x,'Z(M)',5x,'O3(ppb)',  		 	    &
     	    4x,'CO(ppb)',2X,'CH2O(ppt)',1x,'SO2(ppt)',		    &
     	    2x,'SUL(ppt)',2x,'NOx(ppt)',2x,'H2O2(ppt)'/)	    	  
!
  z = 0.
  do k = 1, nz
    z = z + dz*fdz0(k)
  enddo
!
  do k = nz,1,-1						
    z  = z - dz*fdz0(k) 					
    write(7,2026)k,z,gas0(k)%o3,gas0(k)%co, 	       &
		     gas0(k)%ch2o*1.e3,gas0(k)%so2*1.e3,   &
		     gas0(k)%h2so4*1.e3,		       &
		     (gas0(k)%xno+gas0(k)%xno2)*1.e3,      &
		     gas0(k)%h2o2*1.e3
  end do
!
2026   format(2x,i3,2x,f8.2,2x,f8.2,2x,f8.2,		       &
     		2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.2,2x,f8.2)
  write(7,2027)
2027	format(1x/)
!
!----------------------------------------------------------!
!   		     	 Initialize to 0		         	   !
!----------------------------------------------------------!
!
  ad_o3(:)     = 0.0
  ad_co(:)     = 0.0
  ad_h2o2(:)   = 0.0
  ad_nox(:)    = 0.0
  ad_ch2o(:)   = 0.0
  ad_ch3o2h(:) = 0.0
  ad_so2(:)    = 0.0
  ad_svi(:)    = 0.0
  ad_nv(:)     = 0.0
! 
  ed_o3(:)     = 0.0
  ed_co(:)     = 0.0
  ed_h2o2(:)   = 0.0
  ed_nox(:)    = 0.0
  ed_ch2o(:)   = 0.0
  ed_ch3o2h(:) = 0.0
  ed_so2(:)    = 0.0
  ed_svi(:)    = 0.0
  ed_nv(:)     = 0.0
!
!----------------------------------------------------------!
!   			Default case, no chemicals        		   !
!----------------------------------------------------------!
!
#else
  gas0(:)%o3     = 0.
  gas0(:)%co     = 0.
  gas0(:)%zco2   = 0.
  gas0(:)%ho     = 0.
  gas0(:)%ho2    = 0.
  gas0(:)%h2o2   = 0.
  gas0(:)%xno    = 0.
  gas0(:)%xno2   = 0.
  gas0(:)%xno3   = 0.
  gas0(:)%xn2o5  = 0.
  gas0(:)%hno3   = 0.
  gas0(:)%ch4    = 0.
  gas0(:)%ch2o   = 0.
  gas0(:)%ch3o2h = 0.
  gas0(:)%so2    = 0.
  gas0(:)%h2so4  = 0.
  gas0(:)%dms    = 0.
#endif
!
!----------------------------------------------------------!
!		   
  return
  end
