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
!  MPI_PREPARE.F:                   
!
!  Purpose:
!	Prepares MPI simulation:
!	  - domain decomposition (1D-2D)
!	  - MPI type definitions
!	  - Array allocations
!	  - Exchange of shared parameters		  
!
!  Author
!	Julien Savre, MISU
!
! ================================================================
!
!	=================
	subroutine decomp_1D
!	=================

#ifdef SPMD
	USE mpi
#endif
	USE gridno
	USE shared_data

	IMPLICIT NONE
	
#if (defined SPMD)
	integer, save :: m, lastp, ierr

	
! ----------------------------------------------------------------------

	lastp = nproc - 1
!
! 1D decomposition:
!
	if ( mypid .ne. 0 ) then 
	  left_p = mypid - 1 
	else	
	  left_p = lastp
	end if
	
	if ( mypid .ne. lastp ) then
	  right_p = mypid + 1
	else
	  right_p = 0
	end if

	m  = (nx-5)/nproc
	
	it_start = mypid*m + 4
	it_end   = it_start + m - 1

	! --- find borders for the patch
	ip_start = it_start - 2
	ip_end   = it_end   + 2
	if (it_start == 4) ip_start = 1
#endif

	return
	 end

!	=================
	subroutine decomp_2D
!	=================

#ifdef SPMD
	USE mpi
#endif
	USE gridno
	USE shared_data

	IMPLICIT NONE

	integer, save :: m, n, i, j, lastp, flag

#if (defined SPMD)
	integer :: ierr	
	
! ----------------------------------------------------------------------

	lastp = nproc - 1
!
!  Decomposition in x
!
	if ( mypid >= nprocy ) then 
          left_p  = mypid - nprocy
	else	
          left_p  = mypid + nproc - nprocy 
	end if
!
	if ( mypid <= lastp-nprocy ) then
          right_p = mypid + nprocy
	else
          right_p = mypid - nproc + nprocy
	end if
!
!  Decomposition in y
!
	if ( mod(mypid,nprocy) > 0 ) then 
          back_p  = mypid - 1
	else	
          back_p  = mypid + nprocy - 1 
	end if
!
	if ( mod(mypid,nprocy) < nprocy - 1 ) then
          front_p = mypid + 1
	else
          front_p = mypid - nprocy + 1
	end if	
!
! --- 2 D decomposition:
!
	m  = (nx-5)/nprocx
	n  = (ny-5)/nprocy

! --- x direction
	flag = 0
	do i = 0, lastp-nprocy+1, nprocy
	  if (mypid >= i .and. mypid < i+nprocy) then
	    it_start = flag*m + 4
	    it_end   = it_start + m - 1
	  endif
	  flag = flag + 1
	enddo
	if (it_start == 3) it_start = 4

! --- y direction
	do i = 0, nprocy-1
	  do j = 0, nprocx-1
	    if (mypid == i+j*nprocy) then
	      jt_start = i*n + 4
	      jt_end   = jt_start + n - 1
	    endif
	  enddo
	enddo	
	if (jt_start == 3) jt_start = 4
	
! --- find borders for the patch
	ip_start = it_start - 2
	ip_end   = it_end   + 2
        jp_start = jt_start - 2
	jp_end   = jt_end   + 2
	if (it_start == 4) ip_start = 1
	if (jt_start == 4) jp_start = 1
#endif

	return
	 end

!	===================
	subroutine allocations ( )
!	===================

#ifdef SPMD
	USE mpi
#endif
	USE gridno
	USE typedef_hydrometeor
	USE typedef_gas
	USE typedef_aq
	USE typedef_solid
	USE typedef_svalue
	USE shared_state
	USE shared_pressure
	USE shared_turbu
	USE shared_wind
	USE shared_surf
	USE shared_nuclei
	USE shared_thermo
	USE shared_rad
	USE shared_data
	USE shared_hydro
	USE shared_aerosol_new
	USE shared_tend
	USE shared_diag

	IMPLICIT NONE
!
! --- define global variables:
!
#ifdef MEANDATA
	integer :: mean_cnt	! step counter for averaging
#endif   

!	------------------------------------

#ifdef CHEM_ENABLE     
        type (gas_chemical), dimension(:,:,:), pointer :: gas
#ifdef MEANDATA
        type (gas_chemical), dimension(:,:,:), pointer :: gas_mean
#endif
#endif

!	-------------------------------------
#ifdef AQCHEM_ENABLE    
        type (aq_chemical),  dimension(:,:,:), pointer :: aqc, aqr
#ifdef MEANDATA
        type (aq_chemical),  dimension(:,:,:), pointer :: aqc_mean, aqr_mean
#endif
#endif

!	-------------------------------------
#ifdef SOLIDCHEM_ENABLE    
        type (solid_chemical),  dimension(:,:,:),  pointer :: solidi
#ifdef MEANDATA
        type (solid_chemical),  dimension(:,:,:),  pointer :: solidi_mean
#endif
#endif
!
	save
!	    	    
! -----------------------------------------------------------------
!
        if (verbose > 0) call write_debug('Starting allocations ')
!
! --- allocate variables:
!
       call alloc ( wind_adv )
       call alloc ( wind )
       call alloc ( wind2 )
       call alloc ( winds )
       
       call alloc ( state )
       call alloc ( state2 )
       call alloc ( states )
       if (with_piggy) call alloc ( statep )
       
       call alloc ( pressure )
       call alloc ( pressure2 )
       call alloc ( pressures )
       
       allocate ( htend(1:nhydro) )
       allocate ( hydromtr(1:nhydro) )
       allocate ( hydromtr2(1:nhydro) )
       allocate ( hydromtrs(1:nhydro) )
       if (with_piggy) allocate ( hydromtrp(1:nhydro) )
       call alloc ( hydromtr )
       call alloc ( hydromtr2 )
       call alloc ( hydromtrs )
       if (with_piggy) call alloc ( hydromtrp )

       call alloc ( turbu )
       call alloc ( turbu_diag )
       
       call alloc ( thermo )
       call alloc ( thermo_prop )
       
       call alloc ( surf )
       call alloc ( column)
       
       call alloc ( rad )
       
#ifdef SEIFERT		
       call alloc ( nuc )
       call alloc ( nuc2 )
#endif
			    
#ifdef NUC_CNT
       allocate ( itend(1:3) )
       call alloc ( nucin )	
       call alloc ( nucin2 )	
       call alloc ( nucins )	
       
#ifdef MEANDATA
       call alloc ( nucin_mean )
#endif
#endif

#ifdef AERO_ENABLE
       allocate ( aero3d(1:nmode) )
       allocate ( aero3d2(1:nmode) )
       allocate ( aero3ds(1:nmode) )
       allocate ( atend(1:nmode) )
       call alloc ( aero3d )
       call alloc ( aero3d2 )
       call alloc ( aero3ds )
       call alloc ( atend )
       
#ifdef MEANDATA
       allocate ( aero3d_mean(1:nmode) )
       call alloc ( aero3d_mean )
#endif
#endif
       
       allocate ( diag(1:ndiag) )
       call alloc ( diag )
       call alloc ( diagnos )
       call alloc ( entrain )

#ifdef MEANDATA

       call alloc ( wind_mean )

       call alloc ( state_mean )
       
       call alloc ( pressure_mean )
      
       allocate ( hydromtr_mean(1:nhydro) )
       call alloc ( hydromtr_mean )
        
       call alloc ( thermo_mean )
       
       call alloc ( turbu_mean )
       
       call alloc ( surf_mean )
       
       call alloc ( rad_mean )
       
       call alloc ( nuc_mean )

       allocate ( diag_mean(1:ndiag) )
       call alloc (diag)

#endif

       if ( sorder > 1 ) then
         call alloc ( htend )
         call alloc ( tend )
#ifndef AERO_ENABLE
         call alloc ( ntend )
#endif
#ifdef NUC_CNT
         call alloc ( itend )
#endif
       endif

#ifdef CHEM_ENABLE
	allocate ( gas(ip_start:ip_end,jp_start:jp_end,1:nz) )
#ifdef MEANDATA
	allocate ( gas_mean(ip_start:ip_end,jp_start:jp_end,1:nz) )
#endif
#endif

#ifdef AQCHEM_ENABLE
	allocate ( aqc(ip_start:ip_end,jp_start:jp_end,1:nz) )
	allocate ( aqr(ip_start:ip_end,jp_start:jp_end,1:nz) )
#ifdef MEANDATA
	allocate ( aqc_mean(ip_start:ip_end,jp_start:jp_end,1:nz) )
	allocate ( aqr_mean(ip_start:ip_end,jp_start:jp_end,1:nz) )
#endif
#endif

#ifdef SOLIDCHEM_ENABLE
	allocate ( solidi(ip_start:ip_end,jp_start:jp_end,1:nz) )
#ifdef MEANDATA
	allocate ( solidi_mean(ip_start:ip_end,jp_start:jp_end,1:nz) )
#endif
#endif
!
        if (verbose > 0) call write_debug('Terminating allocations ')
!
return
end
!
!	===================
	subroutine deallocations ()
!	===================

#ifdef SPMD
	USE mpi
#endif
	USE gridno
	USE typedef_hydrometeor
	USE typedef_gas
	USE typedef_aq
	USE typedef_solid
	USE typedef_svalue
	USE shared_state
	USE shared_pressure
	USE shared_turbu
	USE shared_wind
	USE shared_surf
	USE shared_nuclei
	USE shared_thermo
	USE shared_rad
	USE shared_data
	USE shared_hydro
	USE shared_aerosol_new
	USE shared_tend
	USE shared_diag

	IMPLICIT NONE
!
! --- define global variables:
!
#ifdef MEANDATA
	integer :: mean_cnt	! step counter for averaging
#endif   

!	------------------------------------

#ifdef CHEM_ENABLE     
        type (gas_chemical), dimension(:,:,:), pointer :: gas
#ifdef MEANDATA
        type (gas_chemical), dimension(:,:,:), pointer :: gas_mean
#endif
#endif

!	-------------------------------------
#ifdef AQCHEM_ENABLE    
        type (aq_chemical),  dimension(:,:,:), pointer :: aqc, aqr
#ifdef MEANDATA
        type (aq_chemical),  dimension(:,:,:), pointer :: aqc_mean, aqr_mean
#endif
#endif

!	-------------------------------------
#ifdef SOLIDCHEM_ENABLE    
        type (solid_chemical),  dimension(:,:,:),  pointer :: solidi
#ifdef MEANDATA
        type (solid_chemical),  dimension(:,:,:),  pointer :: solidi_mean
#endif
#endif
!
	save
!
! ===	Free defined type and delocate memory
!	
	call dealloc ( wind_adv )
	call dealloc ( wind )
	call dealloc ( wind2 )
	call dealloc ( winds )
	
	call dealloc ( state )
	call dealloc ( state2 )
	if (with_piggy) call dealloc ( statep )
	call dealloc ( states )
	
	call dealloc ( pressure )
	call dealloc ( pressure2 )
	call dealloc ( pressures )
	
	call dealloc ( hydromtr )
	call dealloc ( hydromtr2 )
	call dealloc ( hydromtrs )
	if (with_piggy) call dealloc ( hydromtrp )
	deallocate ( hydromtr )
	deallocate ( hydromtr2 )
	deallocate ( hydromtrs )
	if (with_piggy) deallocate ( hydromtrp )
	
	call dealloc ( thermo )
	call dealloc ( thermo_prop )
	
	call dealloc ( turbu )
	call dealloc ( turbu_diag )
	
	call dealloc ( surf )
	call dealloc ( column )
	
	call dealloc ( rad )
	
#ifdef SEIFERT		
	call dealloc ( nuc )
	call dealloc ( nuc2 )
#endif
	
#ifdef NUC_CNT
	call dealloc ( nucin )
	call dealloc ( nucin2 )
	call dealloc ( nucins )

#ifdef MEANDATA
	call dealloc ( nucin_mean )
#endif
#endif

#ifdef AERO_ENABLE
        call dealloc (aero3d)
        call dealloc (aero3d2)
        call dealloc (aero3ds)
	deallocate ( aero3d )
	deallocate ( aero3d2 )
	deallocate ( aero3ds )
	
#ifdef MEANDATA
        call dealloc (aero3d_mean)
	deallocate ( aero3d_mean )
#endif
#endif

	call dealloc (diag)
	deallocate ( diag )
	call dealloc ( diagnos )
	call dealloc ( entrain )

#ifdef MEANDATA
        call dealloc ( wind_mean )
	
        call dealloc ( state_mean )
	
	call dealloc ( pressure_mean )
	
	call dealloc ( hydromtr_mean )
	deallocate ( hydromtr_mean )
	
	call dealloc ( turbu_mean )
	
	call dealloc ( thermo_mean )
	
	call dealloc ( surf_mean )
	
	call dealloc ( rad_mean )

     	call dealloc ( nuc_mean )

	call dealloc (diag_mean)
	deallocate ( diag_mean )
#endif

        if ( sorder > 1 ) then 
	  call dealloc ( tend )
	  call dealloc ( htend )
	  deallocate ( htend )
#ifdef AERO_ENABLE
	  call dealloc ( atend )
	  deallocate ( atend )
#else
         call dealloc ( ntend )
#endif
#ifdef NUC_CNT
         call dealloc ( itend )
         deallocate ( itend )
#endif
       endif

#ifdef CHEM_ENABLE
	deallocate ( gas )
#ifdef MEANDATA
	deallocate ( gas_mean )
#endif
#endif

#ifdef AQCHEM_ENABLE
	deallocate ( aqc )
	deallocate ( aqr )
#ifdef MEANDATA
	deallocate ( aqc_mean )
	deallocate ( aqr_mean )
#endif
#endif

#ifdef SOLIDCHEM_ENABLE
	deallocate ( solidi )
#ifdef MEANDATA
	deallocate ( solidi_mean )
#endif
#endif
!
return
end

!	====================
	subroutine isendcom 
!	====================

#ifdef SPMD
	USE mpi
#endif
	USE gridno
	USE shared_data
	USE shared_aerosol_new
	IMPLICIT NONE
		
#if ( defined SPMD )
	integer :: ierr, im

! ================================================================
!
!  cm_run parameters
!
	call MPI_BCAST(ldtfix,    1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ntau,      1, INTTYPE,0,MPI_COMM_WORLD,ierr)
  	call MPI_BCAST(verbose,   1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(dt0,       1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(tstart,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(tstop,	  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(limit_ts,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(new_run,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(nest_run,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
!
!  cm_init parameters
!
    	call MPI_BCAST(j_day,	   1, INTTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(kpert,	   1, INTTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(sca_set,    1, INTTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(psurf,	   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(dpt,	   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(dqv,	   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(dw,	   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(ctr_lat,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(t_local,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(tpert,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(rep,        1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(casename,  len(casename), CARTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(datadir,   len(datadir), CARTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(incdir,    len(incdir), CARTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(file_init, len(file_init), CARTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(file_rest, len(file_rest), CARTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(forcing_file,  len(forcing_file), CARTYPE,0,MPI_COMM_WORLD,ierr)
!	
!  cm_grid parameters
!	
	call MPI_BCAST(dx,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(dy,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(dz,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(xc_nest, 1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(yc_nest, 1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lx_nest, 1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ly_nest, 1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ztop,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(z1,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)	
	call MPI_BCAST(z2,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(sratio,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(zdamp,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)		
	call MPI_BCAST(dxdamp,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)			
        call MPI_BCAST(dydamp,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)			
        call MPI_BCAST(tdamp,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)		
        call MPI_BCAST(tnudg,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)		
   	call MPI_BCAST(rep     ,1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
   	call MPI_BCAST(gridfile, len(gridfile), CARTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(sponge,	1, INTTYPE,0,MPI_COMM_WORLD,ierr)
!
!  cm_out parameters
!
	call MPI_BCAST(iax,       1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iav,       1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ires,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)	
	call MPI_BCAST(inest,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)	
	call MPI_BCAST(its,       1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ipro,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(isli,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(zbl,       1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ztrop,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qthres,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(wthres,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(slicex_io, 16, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(slicey_io, 16, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(slicez_io, 16, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nslicex,   1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nslicey,   1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nslicez,   1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(kout,      1, INTTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(ts_out,    len(ts_out), CARTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(no_out,    1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(spec_diag, 1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(out_surf,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(out_hov,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(out_yav,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(all_rest,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
!
!  cm_num parameters
!
    	call MPI_BCAST(scal_adv,  len(scal_adv), CARTYPE,0,MPI_COMM_WORLD,ierr)
   	call MPI_BCAST(limit,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(p_mcons,	1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(imp_buoy,1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(split_mic,1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lavisc,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lddamp,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(sorder,	1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(mom_ord,	1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(diff_ord,1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nsubp,	1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(fft_solv,1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(cfl_max, 1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(cfl_min, 1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lim_tol,	1, REALTYPE,0,MPI_COMM_WORLD,ierr)
!
!  cm_phys parameters
!
    	call MPI_BCAST(with_mom,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_scal, 1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_adv,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_dif,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_buoy, 1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_mic,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(anis_k,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(nl_sgs,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_nudg, 1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_lsadv, 1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_lssrc, 1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_lssub, 1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_cor,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_rad,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_piggy,1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(with_tvar, 1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(cst_cp,    1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rad_sw,    1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rad_o3,    1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iradx,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(zdec,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(pran,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(diff,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(Ddiv,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)			
	call MPI_BCAST(w_up,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)			
	call MPI_BCAST(u0shift,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(v0shift,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	
!
!  cm_bc parameters
!
    	call MPI_BCAST(bcl,     len(bcl), CARTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(isurf,   1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(momsf,   1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ust,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)	
	call MPI_BCAST(shf0,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)	
	call MPI_BCAST(lhf0,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(scf,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)	
	call MPI_BCAST(sst,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ssm,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(c_dm,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(c_ds,    1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(zrough,  1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(alb,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(emi,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
!
!  cm_micro parameters
!	
	call MPI_BCAST(ice_habit,   len(ice_habit), CARTYPE,0,MPI_COMM_WORLD,ierr)	
    	call MPI_BCAST(micro_dif,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(lvent,       1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(ldrizz,      1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(lrime,       1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(ldtmic,      1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lmicro,      1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(laero,       1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aero_flg%tran,     1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aero_flg%act,      1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aero_flg%act_scv,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aero_flg%imp_scv,  1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aero_flg%reg,      1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aero_flg%chem,     1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lndrop,      1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(lfreeze,     1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(auto,        1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(moments,     1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ice_delay,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(xauto,       1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qauto,       1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(dtcon,       1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(xnc0_d,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(xnc0_s,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(xnc0_k,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(xn_ccn0,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(xn_in0,      1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nlsod,       1, INTTYPE,0,MPI_COMM_WORLD,ierr)		
!
!  cm_aero parameters
!	
#ifdef AERO_ENABLE
	call MPI_BCAST(nmode0,            1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aero_sfc_source,   1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
    	call MPI_BCAST(lkin,              1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	do im = 1, nmode0
	  call MPI_BCAST(aeroi(im)%nelem,               1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeroi(im)%init%frac,           4, REALTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeroi(im)%init%n0,             1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeroi(im)%size%rmean,          1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeroi(im)%size%sigma,          1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeroi(im)%init%present,        4, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeroi(im)%init%n_sfc_source,   1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	enddo
!
	if (reg_mode) then
	  call MPI_BCAST(aeror%nelem,          1, INTTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeror%init%frac,      4, REALTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeror%init%n0,        1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeror%size%rmean,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeror%size%sigma,     1, REALTYPE,0,MPI_COMM_WORLD,ierr)
	  call MPI_BCAST(aeror%init%present,   4, LOGTYPE,0,MPI_COMM_WORLD,ierr)
	endif
#endif
!
!  cm_lag parameters
!
#ifdef LAGRANGE
        call MPI_BCAST(aerosol_lag,	1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(compos_lag,	len(compos_lag), CARTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(mu_lag,  	1, REALTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(sigma_lag,	1, REALTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ilag,		1, REALTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(zl1,		1, REALTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(zl2,		1, REALTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nlag,		1, INTTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(res_lag, 	1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(lag_init,	1, INTTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(lag_ord, 	1, INTTYPE,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(lag_mix, 	1, LOGTYPE,0,MPI_COMM_WORLD,ierr)
#endif
!	
#ifdef CHEM_ENABLE			       
	call MPI_BCAST(rktable1,10200,REALTYPE,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rktable2,  900,REALTYPE,0,MPI_COMM_WORLD,ierr)
#endif

#endif

	return
	end
!
!	====================
	subroutine create_mpi_types 
!	====================
!
!  Creates MPI derived data types to handle array communication
!

#ifdef SPMD
	USE mpi
	USE gridno
	USE shared_data
!
	IMPLICIT NONE
!	
	integer :: ierr
	integer :: dimx, dimy, dimz, dimxt, dimyt
	
	dimx = ip_end-ip_start+1
	dimy = jp_end-jp_start+1
	dimz = nz
!
	dimxt = it_end-it_start+1
	dimyt = jt_end-jt_start+1
!
!  Types for full domain collection
!	
	call MPI_TYPE_CREATE_SUBARRAY (3, (/dimx,dimy,dimz/), (/dimxt,dimyt,dimz/), (/it_start-ip_start,jt_start-jp_start,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, BLOCK3D, ierr)
	call MPI_TYPE_COMMIT (BLOCK3D, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (2, (/dimx,dimy/), (/dimxt,dimyt/), (/it_start-ip_start,jt_start-jp_start/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, BLOCK2D_r, ierr)
	call MPI_TYPE_COMMIT (BLOCK2D_r, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (2, (/dimx,dimy/), (/dimxt,dimyt/), (/it_start-ip_start,jt_start-jp_start/), 	&
		MPI_ORDER_FORTRAN, INTTYPE, BLOCK2D_i, ierr)
	call MPI_TYPE_COMMIT (BLOCK2D_i, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (3, (/dimx,dimy,dimz/), (/dimx,dimy,dimz/), (/0,0,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, BLOCK3Dp, ierr)
	call MPI_TYPE_COMMIT (BLOCK3Dp, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (2, (/dimx,dimy/), (/dimx,dimy/), (/0,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, BLOCK2Dp, ierr)
	call MPI_TYPE_COMMIT (BLOCK2Dp, ierr)
!
!  Types for boundary ghost points
!
	call MPI_TYPE_CREATE_SUBARRAY (3, (/dimx,dimy,dimz/), (/3,dimy,dimz/), (/it_start-ip_start,0,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, GHOSTX_send_s, ierr)
	call MPI_TYPE_COMMIT (GHOSTX_send_s, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (3, (/dimx,dimy,dimz/), (/3,dimy,dimz/), (/ip_end-ip_start-4,0,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, GHOSTX_send_e, ierr)
	call MPI_TYPE_COMMIT (GHOSTX_send_e, ierr)

#if (defined MODEL_3D) && (defined DECOMP_2D)
	call MPI_TYPE_CREATE_SUBARRAY (3, (/dimx,dimy,dimz/), (/dimx,3,dimz/), (/0,jt_start-jp_start,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, GHOSTY_send_s, ierr)
	call MPI_TYPE_COMMIT (GHOSTY_send_s, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (3, (/dimx,dimy,dimz/), (/dimx,3,dimz/), (/0,jp_end-jp_start-4,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, GHOSTY_send_e, ierr)
	call MPI_TYPE_COMMIT (GHOSTY_send_e, ierr)
#endif
!
!  Types for boundary ghost points for slices
!
#ifdef MODEL_3D	
	call MPI_TYPE_CREATE_SUBARRAY (2, (/dimx,dimy/), (/3,dimy/), (/it_start-ip_start,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, GHOSTX_send_s_2d, ierr)
	call MPI_TYPE_COMMIT (GHOSTX_send_s_2d, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (2, (/dimx,dimy/), (/3,dimy/), (/ip_end-ip_start-4,0/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, GHOSTX_send_e_2d, ierr)
	call MPI_TYPE_COMMIT (GHOSTX_send_e_2d, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (2, (/dimx,dimy/), (/dimx,3/), (/0,jt_start-jp_start/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, GHOSTY_send_s_2d, ierr)
	call MPI_TYPE_COMMIT (GHOSTY_send_s_2d, ierr)
!
	call MPI_TYPE_CREATE_SUBARRAY (2, (/dimx,dimy/), (/dimx,3/), (/0,jp_end-jp_start-4/), 	&
		MPI_ORDER_FORTRAN, REALTYPE, GHOSTY_send_e_2d, ierr)
	call MPI_TYPE_COMMIT (GHOSTY_send_e_2d, ierr)
#endif
!
#endif
!
    return
    end	
!
!	====================
	subroutine free_mpi_types 
!	====================
!
!  Creates MPI derived data types to handle array communication
!
#ifdef SPMD
	USE mpi
	USE gridno
	USE shared_data
!
	IMPLICIT NONE
	
	integer :: ierr
!
	call MPI_TYPE_FREE (BLOCK3D, ierr)
	call MPI_TYPE_FREE (BLOCK2D_r, ierr)
	call MPI_TYPE_FREE (BLOCK2D_i, ierr)
	call MPI_TYPE_FREE (BLOCK3Dp, ierr)
	call MPI_TYPE_FREE (BLOCK2Dp, ierr)
!
	call MPI_TYPE_FREE (GHOSTX_send_s, ierr)
	call MPI_TYPE_FREE (GHOSTX_send_e, ierr)

#if (defined MODEL_3D) && (defined DECOMP_2D)
	call MPI_TYPE_FREE (GHOSTY_send_s, ierr)
	call MPI_TYPE_FREE (GHOSTY_send_e, ierr)
#endif
!
#ifdef MODEL_3D	
	call MPI_TYPE_FREE (GHOSTX_send_s_2d, ierr)
	call MPI_TYPE_FREE (GHOSTX_send_e_2d, ierr)

	call MPI_TYPE_FREE (GHOSTY_send_s_2d, ierr)
	call MPI_TYPE_FREE (GHOSTY_send_e_2d, ierr)
#endif
!
#endif
!
    return
    end	
