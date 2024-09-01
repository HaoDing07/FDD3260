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
!  output.f90                   
!
!  Purpose:
!      A package for handling data output of MIMICA                    
!
!  Author
!      Chien Wang
!      MIT Joint Program on Science and Policy of Global Change
!
! ================================================================

!      =================================================
      subroutine output_main 
!      =================================================                            

USE shared_all
USE diagnostics
USE netcdfmod
USE netcdfslice
USE lagrange
#ifdef LAGRANGE
USE shared_lagrange
#endif
  
  IMPLICIT NONE
!
! --- define global variables:
!
  integer  :: i, h, mean_cnt
  real, dimension(ntsout)   :: ts_tab
  real :: cnt_step
!
  logical  :: lnew, per_iostep, per_tsstep, per_avstep, per_prstep, per_slstep, per_lagstep
  character(len = 100):: filename
!
  save
!
  if (verbose > 0) call write_debug('Start output_main')
!
!----------------------------------------------------------!
!                      Set logicals                        !
!----------------------------------------------------------!
!
  lnew = new_run.and.ntime==0
  per_iostep = (time-tstart >= real(n_io)*iax) .or. (time>=tstop) .or. (.not.new_run.and.nite==1)
  per_tsstep = (time-tstart >= real(n_ts)*its) .or. (time>=tstop) .or. (.not.new_run.and.nite==1)
  per_prstep = (time-tstart >= real(n_pr)*ipro) .or. (time>=tstop) .or. (.not.new_run.and.nite==1)
  per_slstep = (time-tstart >= real(n_sl)*isli) .or. (time>=tstop) .or. (.not.new_run.and.nite==1)
  per_lagstep = (time-tstart >= real(n_lag)*ilag) .or. (time>=tstop) .or. (.not.new_run.and.nite==1)
  per_avstep = (time-tstart >= real(n_av)*iav) .or. (time==tstart) .or. (.not.new_run.and.nite==1)
!
  if (time-tstart >= real(n_io)*iax) n_io=n_io+1 
  if (time-tstart >= real(n_av)*iav) n_av=n_av+1 
  if (time-tstart >= real(n_ts)*its) n_ts=n_ts+1 
  if (time-tstart >= real(n_pr)*ipro) n_pr=n_pr+1 
  if (time-tstart >= real(n_sl)*isli) n_sl=n_sl+1 
  if (time-tstart >= real(n_lag)*ilag) n_lag=n_lag+1 
!
!----------------------------------------------------------!
!		  Cumulate time averages                   !
!----------------------------------------------------------!
!
#ifdef MEANDATA
  call calmean_cum  ( mean_cnt )
#else
  mean_cnt=1
#endif
!
  cnt_step = 1./real(mean_cnt)
!
!----------------------------------------------------------!
!	      Calculate diagnostic variables               !
!----------------------------------------------------------!
!
  if ( lnew .or. per_tsstep .or. per_iostep .or. per_prstep .or. per_slstep ) then
    call thermo_diagnostics
!
    call special_diagnostics
  endif
!
!----------------------------------------------------------!
!                      Time series                         !
!----------------------------------------------------------!
!
  if (per_tsstep) then
    call time_series (ts_tab)
!
    if (mypid == 0) then
      opdir = './OUTPUT'
      filename = trim(opdir) // '/T_S'   
      open(202,file=filename(1:len_trim(filename)),     &
  	     form='formatted',			        &
  	     access='SEQUENTIAL', 		        &
  	     position=how_access, 		        &
  	     status='unknown')
!
      if (ntime > 0) call repositionts (202)
!
      call output_ts ( 202, ts_tab )
!
      close(202)
!
      write(7,901)
      flush(7)
    endif
  endif
!
!----------------------------------------------------------!
!                     Output 3D data                       !
!----------------------------------------------------------!
!
  if (per_iostep .and. .not.no_out) then
!      
!  Call netcdf interface for outputs
!
#ifdef MODEL_3D
    call write_dim ( 3, cnt_step, lnew, .false. )
#else
    call write_dim ( 2, cnt_step, lnew, .false. )
#endif
!
! Write chemicals and aerosols in formatted output
! 
#ifdef CHEM_ENABLE
    call output_gas ( gas_tmp, cnt_step )
#endif

#ifdef AQCHEM_ENABLE
    call output_aq ( aqc_tmp, aqr_tmp, cnt_step )
#endif

#ifdef SOLIDCHEM_ENABLE
    call output_solid ( solidi_tmp, cnt_step )
#endif
!
    if (mypid==0) then
      write(7,902)
      flush(7)
    endif
!
  endif
!
!----------------------------------------------------------!
!                     Output slices                        !
!----------------------------------------------------------!
!
#ifdef MODEL_3D
  if ( per_slstep ) then
!
    if ( max(nslicex,nslicey,nslicez) > 0 ) call write_2d_slice ( cnt_step, lnew )
!
    if (mypid==0) then
      write(7,903)
      flush(7)
    endif
!
    if ( out_surf ) call write_2d_surf ( cnt_step, lnew )
!
    if (out_surf.and.mypid==0) then
      write(7,904)
      flush(7)
    endif
!
    if ( out_yav ) call write_yav ( cnt_step )
!
  endif
#endif
!
!----------------------------------------------------------!
!                     Output profiles                      !
!----------------------------------------------------------!
!
  if (per_prstep .and. nvarpro > 0) then
!
    do h = 1, nvarpro
      call write_2d ( cnt_step, lnew, prname(h) )
    enddo
!
    if ( out_hov .and. out_surf ) call write_hov ( cnt_step )
!
    if (mypid==0) then
      write(7,905)
      flush(7)
    endif
!
  endif
!
!----------------------------------------------------------!
!                   Lagrangian parcels                     !
!----------------------------------------------------------!
!
#ifdef LAGRANGE
  if ( per_lagstep ) then
    call lagrange_interp ( .true., .false. )
!
    call write_lagrange
!
    if (mypid==0) then
      write(7,906)
      flush(7)
    endif
  endif
#endif
!
!----------------------------------------------------------!
!              Terminating and closing files               !
!----------------------------------------------------------!
!
!  Reset means if needed
!
#ifdef MEANDATA
  if (per_avstep) call calmean_reset  ( cnt_step )
#endif
!
!  Closing files
!  
  if (verbose > 0) call write_debug('Terminate output_main')
!
!----------------------------------------------------------!
!
901 format(11x,"Finished writing time series in OUTPUT/T_S")
902 format(11x,"Finished writing 2D/3D outputs in OUTPUT/output.nc")
903 format(11x,"Finished writing slices in OUTPUT/slice_xxx.nc")
904 format(11x,"Finished writing surface slice in OUTPUT/slice_surf.nc")
905 format(11x,"Finished writing profiles in OUTPUT/profiles_tot.nc")
906 format(11x,"Finished writing lagrangian parcels in OUTPUT/lagrange.nc")
!
return
end
!
!      =================================================
      subroutine output_pig 
!      =================================================                            

USE shared_all
USE diagnostics
USE netcdfmod
USE netcdfslice
  
  IMPLICIT NONE
!
! --- define global variables:
!
  logical  :: lnew, per_iostep, per_prstep, per_slstep
!
  save
!
!----------------------------------------------------------!
!                      Set logicals                        !
!----------------------------------------------------------!
!
  lnew = new_run.and.(ntime==0) 
  per_prstep = (time-tstart >= real(n_pr)*ipro) .or. (time>=tstop-0.5*dt0) .or. (.not.new_run.and.nite==1)
  per_iostep = (time-tstart >= real(n_io)*iax) .or. (time>=tstop-0.5*dt0) .or. (.not.new_run.and.nite==1)
  per_slstep = (time-tstart >= real(n_sl)*isli) .or. (time>=tstop-0.5*dt0) .or. (.not.new_run.and.nite==1)
!
!----------------------------------------------------------!
!                    Extra diagnostics                     !
!----------------------------------------------------------!
!
  if ( per_prstep .or. per_iostep .or. per_slstep ) then
    call thermo_diagnostics
!
    call special_diagnostics
  endif
!
!----------------------------------------------------------!
!                       Output data                        !
!----------------------------------------------------------!
!
!  Output profiles
!
  if ( per_prstep ) call write_2d ( 1., lnew, prname(1), pig=.true. )
!
!  Output 3D
!
  if ( per_iostep .and. .not.no_out ) then
#ifdef MODEL_3D
    call write_dim ( 3, 1., lnew, .false., pig=.true. )
#else
    call write_dim ( 2, 1., lnew, .false., pig=.true. )
#endif
  endif
!
!  Output slices   
!
#ifdef MODEL_3D
    if ( per_slstep .and. max(nslicex,nslicey,nslicez) > 0 ) call write_2d_slice ( 1., lnew )
#endif
!
return
end
!
!      =============================
      subroutine output_gas ( gas, cnt )
!      =============================

      USE gridno
      USE typedef_gas
      USE shared_data

      IMPLICIT NONE
      
      real    :: cnt
#ifdef CHEM_ENABLE
      type (gas_chemical), dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: gas
#else
        type (gas_chemical) :: gas
#endif

#ifdef CHEM_ENABLE
      real    :: data_max, data_min, aw_x
      integer :: i, j, k
      
      integer, dimension(3) :: index_max, index_min
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: x

! ============================================================

      !
      ! --- O3: in ppbv
      !
      aw_x = 28.97296245/awO3
      call outx ( 'GASO3 ', gas%o3*aw_x, cnt, 1)

      !
      ! --- CO: in ppbv
      !
      aw_x = 28.97296245/awCO
      call outx ( 'GASCO ', gas%co*aw_x, cnt, 1)

       !
       ! --- OH: 10^3 radicals/cm^3
       !

       ! === ppb(m) to 10^5 radicals/cm^3
       do k = 1, nz
       do j = jp_start,jp_end
       do i = ip_start,ip_end
         x(i,j,k) = gas(i,j,k)%ho*den0(k)*6.022e6/17.0
       end do
       end do
       end do

       aw_x = 100.0
       call outx ( 'GASHO ', x*aw_x, cnt, 1)


       if ( masterp ) then             
         write(7,9424)data_max,index_max,                        &
                       data_min,index_min
       end if

       !
       ! --- HO2:
       !
       aw_x = 28.97296245/awHO2
      call outx ( 'GASHO2', gas%ho2*aw_x, cnt, 1)


       !
       ! --- H2O2: pptv
       !
       aw_x = 28.97296245/awH2O2*1.e3
             call outx ( 'GASCO ', gas%co*aw_x, cnt, 1)

       !
       ! --- NH3: pptv
       !
       aw_x = 28.97296245/awNH3*1.e3
      call outx ( 'GASNH3', gas%co*aw_x, cnt, 1)


       !
       ! --- NO:
       !
       aw_x = 28.97296245/awNO
      call outx ( 'GASNO ', gas%xno*aw_x, cnt, 1)


       !
       ! --- NO2:
       !
       aw_x = 28.97296245/awNO2
      call outx ( 'GASNO2', gas%xno2*aw_x, cnt, 1)


       !
       ! --- HNO3: pptv
       !
       aw_x = 28.97296245/awHNO3*1.e3
       call outx ( 'GASHNO3', gas%hno3*aw_x, cnt, 1)


       !
       ! --- CH4: ppbv
       !
       aw_x = 28.97296245/awCH4
       call outx ( 'GASCH4', gas%ch4*aw_x, cnt, 1)


       !
       ! --- CH2O: pptv
       !
       aw_x = 28.97296245/awCH2O*1.e3
       call outx ( 'GASCH2O', gas%ch2o*aw_x, cnt, 1)


       !
       ! --- CH3O2H: pptv
       !
       aw_x = 28.97296245/awCH3O2H*1.e3
       call outx ( 'GASCH3O2H', gas%ch3o2h*aw_x, cnt, 1)


       !
       ! --- SO2: pptv
       !
       aw_x = 28.97296245/awSO2*1.e3
       call outx ( 'GASSO2', gas%so2*aw_x, cnt, 1)


       !
       ! --- H2SO4: pptv
       !
       aw_x = 28.97296245/awH2SO4*1.e3
       call outx ( 'GASH2SO4', gas%h2so4*aw_x, cnt, 1)


       !
       ! --- DMS: 0.1 pptv
       !
       aw_x = 10.*28.97296245/awDMS*1.e3
       call outx ( 'GASDMS', gas%dm2*aw_x, cnt, 1)


#endif
      return
        end subroutine output_gas
!
!      ==================================
      subroutine output_aq ( aqc, aqr, cnt )
!      ==================================

      USE gridno
      USE typedef_aq
      USE shared_data

      IMPLICIT NONE
      
!
! === Subroutine to write aq values into files
!

      real    :: cnt
#ifdef AQCHEM_ENABLE
      type (aq_chemical), dimension(ip_start:ip_end,jp_start:jp_end,1:nz)       &
			:: aqc, aqr
#else
      type (aq_chemical) :: aqc, aqr
#endif

#ifdef AQCHEM_ENABLE
      real    :: data_max, data_min
      integer,dimension(3) :: index_max, index_min

! ============================================================

      !
      ! --- O3:
      !
      call outx ( 'AQCO3  ', aqc%o3, cnt, 0)
      call outx ( 'AQRO3  ', aqr%o3, cnt, 0)

       !
       ! --- C(IV):
       !
      call outx ( 'AQCCIV ', aqc%civ, cnt, 0)
      call outx ( 'AQRCIV ', aqr%civ, cnt, 0)

       !       
       ! --- H2O2:
       !
      call outx ( 'AQCH2O2', aqc%h2o2, cnt, 0)
      call outx ( 'AQRH2O2', aqr%h2o2, cnt, 0)

       !       
       ! --- NH4:
       !
      call outx ( 'AQCNH4 ', aqc%nh4, cnt, 0)
      call outx ( 'AQRNH4 ', aqr%nh4, cnt, 0)

       !
       ! --- N(V):
       !
      call outx ( 'AQCNV  ', aqc%xnv, cnt, 0)
      call outx ( 'AQRNV  ', aqr%xnv, cnt, 0)

       !       
       ! --- S(IV):
       !
      call outx ( 'AQCSIV ', aqc%siv, cnt, 0)
      call outx ( 'AQRSIV ', aqr%siv, cnt, 0)

       !
       ! --- S(VI):
       !
      call outx ( 'AQCSVI ', aqc%svi, cnt, 0)
      call outx ( 'AQRSVI ', aqr%svi, cnt, 0)

       !
       ! --- CH2O:
       !
      call outx ( 'AQCCH2O', aqc%ch2o, cnt, 0)
      call outx ( 'AQRCH2O', aqr%ch2o, cnt, 0)

       !
       ! --- CH3O2H:
       !
      call outx ( 'AQCCH3O2H', aqc%ch3o2h, cnt, 0)
      call outx ( 'AQRCH3O2H', aqr%ch3o2h, cnt, 0)

       !
       ! --- pH:
       !
      call outx ( 'AQCPH  ', 1000.*aqc%hplus, cnt, 0)
      call outx ( 'AQRPH  ', 1000.*aqr%hplus, cnt, 0)
#endif

      return
        end subroutine output_aq
!
!      ==================================
      subroutine output_solid ( solidi, cnt )
!      ==================================

      USE gridno
      USE typedef_solid
      USE shared_data

      IMPLICIT NONE

      real :: cnt
#ifdef SOLIDCHEM_ENABLE
      type (solid_chemical), dimension(ip_start:ip_end,jp_start:jp_end,1:nz)      &
			:: solidi
#else
      type (solid_chemical) :: solidi
#endif

#ifdef SOLIDCHEM_ENABLE
      real    :: data_max, data_min
      integer, dimension(3) :: index_max, index_min

! ============================================================
      
      !
      ! --- O3:
      !
      call outx ( 'SOLO3 ', solidi%o3, cnt, 0)

       !       
       ! --- H2O2:
       !
      call outx ( 'SOLH2O2', solidi%h2o2, cnt, 0)

       !       
       ! --- NH4:
       !
      call outx ( 'SOLNH4', solidi%nh4, cnt, 0)

       !
       ! --- N(V):
       !
      call outx ( 'SOLNV ', solidi%xnv, cnt, 0)

       !
       ! --- CH2O, note = use data_min to store max in rain:
       !
      call outx ( 'SOLCH2O', solidi%ch2o, cnt, 0)

       !
       ! --- CH3O2H, note = use data_min to store max in rain:
       !
      call outx ( 'SOLCH3O2H', solidi%ch3o2h, cnt, 0)

       !       
       ! --- S(IV), note = use data_min to store max in rain:
       !
      call outx ( 'SOLSIV', solidi%siv, cnt, 0)

       !
       ! --- S(VI):
       !
      call outx ( 'SOLSVI', solidi%svi, cnt, 0)
#endif

      return
        end subroutine output_solid
!
!      =================================================                            
      subroutine calmean_cum ( mean_cnt )
!      =================================================                            

    USE shared_all
    USE averages

    IMPLICIT NONE

    integer :: mean_cnt      ! step counter for averaging
    integer :: i, j, k, h
    real    :: wav(nz)

! ================================================================
	
	!
	! === accumulation:
	!
    	call horav (wind%w, wav)

	do k=1,nz
	  wind_mean%u(:,:,k) = wind_mean%u(:,:,k) + wind%u(:,:,k)
	  wind_mean%v(:,:,k) = wind_mean%v(:,:,k) + wind%v(:,:,k)
	  wind_mean%w(:,:,k) = wind_mean%w(:,:,k) + (wind%w(:,:,k) - wav(k))
	end do

	state_mean = state_mean + state
	
	pressure_mean = pressure_mean + pressure
	
	hydromtr_mean = hydromtr_mean + hydromtr

	thermo_mean = thermo_mean + thermo

	turbu_mean = turbu_mean + turbu
	
	nuc_mean = nuc_mean + nuc
	
	surf_mean = surf_mean + surf
	
	rad_mean = rad_mean + rad

#ifdef AERO_ENABLE
        aero3d_mean = aero3d_mean + aero3d
#endif
		
	do k=1,nz
	do j=jp_start,jp_end
	do i=ip_start,ip_end


#ifdef CHEM_ENABLE
        gas_mean(i,j,k) = gas_mean(i,j,k) + gas2(i,j,k)
#ifdef AQCHEM_ENABLE
        aqc_mean(i,j,k) = aqc_mean(i,j,k) + aqc2(i,j,k)
        aqr_mean(i,j,k) = aqr_mean(i,j,k) + aqr2(i,j,k)
#endif
#ifdef SOLIDCHEM_ENABLE
        solidi_mean(i,j,k) = solidi_mean(i,j,k) + solidi2(i,j,k)
#endif
#endif
	end do
	end do
	end do
!
	mean_cnt = mean_cnt + 1

      RETURN       
      END             
!  
!      =================================================                            
      subroutine calmean_reset ( mean_cnt )
!      =================================================                            

USE shared_all

IMPLICIT NONE

      integer :: mean_cnt      ! step counter for averaging
      integer :: i, j, k, h

! ================================================================

        write(7,9100)
9100 format(11x,"Reset calculation of time averages")
      
      !
      ! === reset all:
      !
	wind_mean = wind

	surf_mean = surf

	state_mean = state
	
	pressure_mean = pressure
	
	do h=1,nhydro
	  hydromtr_mean(h) = hydromtr(h)
	enddo
	
	thermo_mean = thermo
	
	turbu_mean = turbu
	
	nuc_mean = nuc
	
	rad_mean = rad
	
	surf_mean = surf

#ifdef AERO_ENABLE
        aero3d_mean = aero3d
#endif
			
	do k=1,nz
	do j=jp_start,jp_end
	do i=ip_start,ip_end
#ifdef CHEM_ENABLE
        gas_mean(i,j,k) = gas2(i,j,k)
#ifdef AQCHEM_ENABLE
        aqc_mean(i,j,k) = aqc2(i,j,k)
        aqr_mean(i,j,k) = aqr2(i,j,k)
#endif
#ifdef SOLIDCHEM_ENABLE
        solidi_mean(i,j,k) = solidi2(i,j,k)
#endif
#endif
	end do
	end do
	end do
					 
	mean_cnt = 1
	
	if (mypid == 0) write(7,9111) 

9111 format(11x,"Reset average quantities for output")

      RETURN       
      END      
!      
!      =================================================                            
      subroutine count_outputs
!      =================================================                            
     
      USE gridno
      USE shared_data
      
      IMPLICIT NONE
!
      if (verbose > 0) call write_debug('Start count_outputs')

      nvarout=0
      nvarsurf=0
      nvarpro=0
      ldiag=.false.
!
!  Full output variables
!
      if(out_u) nvarout=nvarout+1
#ifdef MODEL_3D
      if(out_v) nvarout=nvarout+1	
#endif
      if(out_w) nvarout=nvarout+1	
      if(out_p) nvarout=nvarout+2
      if(out_t) nvarout=nvarout+1	
      if(out_pt) nvarout=nvarout+1	
      if(out_ptv) nvarout=nvarout+1
      if(out_mse) nvarout=nvarout+1
      if(out_rho) nvarout=nvarout+2
      if(out_sca) nvarout=nvarout+nscal	
      if(out_qv) nvarout=nvarout+1	
      if(lmicro>0 .and. out_qc) nvarout=nvarout+1	
      if(lmicro>0 .and. out_nc) nvarout=nvarout+1	
      if(lmicro>0 .and. out_dc) nvarout=nvarout+1	
      if(lmicro>0 .and. out_qr) nvarout=nvarout+1	
      if(lmicro>0 .and. out_nr) nvarout=nvarout+1	
      if(lmicro>0 .and. out_dr) nvarout=nvarout+1	
      if(lmicro>0 .and. out_vp) nvarout=nvarout+1	
      if(lmicro>1 .and. out_qi) nvarout=nvarout+1	
      if(lmicro>1 .and. out_ni) nvarout=nvarout+1	
      if(lmicro>1 .and. out_di) nvarout=nvarout+1	
      if(lmicro>2 .and. out_qg) nvarout=nvarout+1	
      if(lmicro>2 .and. out_ng) nvarout=nvarout+1	
      if(lmicro>2 .and. out_dg) nvarout=nvarout+1	
      if(lmicro>3 .and. out_wg) nvarout=nvarout+1	
      if(lmicro>2 .and. out_qs) nvarout=nvarout+1	
      if(lmicro>2 .and. out_ns) nvarout=nvarout+1	
      if(lmicro>2 .and. out_ds) nvarout=nvarout+1	
      if(lmicro>3 .and. out_ws) nvarout=nvarout+1	
      if(lmicro>3 .and. out_qh) nvarout=nvarout+1	
      if(lmicro>3 .and. out_nh) nvarout=nvarout+1	
      if(lmicro>3 .and. out_dh) nvarout=nvarout+1	
      if(lmicro>3 .and. out_wh) nvarout=nvarout+1	
      if(lmicro>0 .and. out_prec) nvarout=nvarout+1
      if(lmicro>1 .and. out_prec) nvarout=nvarout+1
      if(out_qt) nvarout=nvarout+1
      if(out_sat) nvarout=nvarout+2		 
      if(out_mf) nvarout=nvarout+4
      if(out_buoy) nvarout=nvarout+1    
      if(out_k) nvarout=nvarout+1
      if(out_ccn) nvarout=nvarout+1
      if(out_in) nvarout=nvarout+1
      if(out_dtnet) nvarout=nvarout+1
#ifdef RAD_ENABLE
      if(out_dtnet) nvarout=nvarout+2
      if(out_frad) nvarout=nvarout+3
#else
      if(out_frad) out_frad=.false.
#endif    

#ifdef NUC_CNT
      if(out_in) nvarout=nvarout+12
#endif    
#ifdef AERO_ENABLE
      if(out_aero) nvarout=nvarout+nmode*3
#endif 

!
!  Special diagnostics
!
      if (spec_diag) then
        if(out_z) nvarout=nvarout+1
        if(out_beff) nvarout=nvarout+1    
        if(out_tke) nvarout=nvarout+6
        if(out_grad) nvarout=nvarout+15
        if(out_div) nvarout=nvarout+2
        if(out_vort) nvarout=nvarout+3
        if(out_ent) nvarout=nvarout+17
        if(out_dp) nvarout=nvarout+3
        if(out_var.and.out_pt) nvarout=nvarout+1
        if(out_var.and.out_qv) nvarout=nvarout+1
        if(out_flut.and.out_u) nvarout=nvarout+1
        if(out_flut.and.out_v) nvarout=nvarout+1
        if(out_flut.and.out_buoy) nvarout=nvarout+1
        if(out_flut.and.out_pt) nvarout=nvarout+1 
        if(out_flut.and.out_qt) nvarout=nvarout+1
        if(out_flut.and.out_ptv) nvarout=nvarout+1
        if(out_flut.and.out_sca.and.nscal>0) nvarout=nvarout+1
        if(out_fsgs.and.out_u) nvarout=nvarout+1
        if(out_fsgs.and.out_v) nvarout=nvarout+1
        if(out_fsgs.and.out_buoy) nvarout=nvarout+1
        if(out_fsgs.and.out_pt) nvarout=nvarout+1 
        if(out_fsgs.and.out_qt) nvarout=nvarout+1
      else
  	out_z = .false.
	out_beff = .false.
	out_tke = .false.
	out_grad = .false.
	out_div = .false.
	out_vort = .false.
	out_ent = .false.
	out_dp = .false.
	out_var = .false.
	out_flut = .false.
	out_fsgs = .false.
      endif
!
!  Diagnostics
!
      if( out_diagu.or.out_diagv.or.out_diagw.or.out_diagp.or.out_diagt     &
	  .or.out_diagq.or.out_diagl.or.out_diagr.or.out_diagi.or.out_diags.or.out_micro.or.out_diaga ) ldiag=.true.
!
      if(out_diagu) nvarout=nvarout+ndiag	       
      if(out_diagv) nvarout=nvarout+ndiag
      if(out_diagw) nvarout=nvarout+ndiag
      if(out_diagp) nvarout=nvarout+ndiag
      if(out_diagtv) out_diagt=.true.
      if(out_diagtv) out_diagq=.true.
      if(out_diagtv) out_diagl=.true.
      if(out_diagtv) out_diagr=.true.
      if(out_diagt) nvarout=nvarout+ndiag
      if(out_diagtv) nvarout=nvarout+ndiag
      if(out_diagk) nvarout=nvarout+ndiag
      if(out_diagq) nvarout=nvarout+ndiag
      if(out_diagl) nvarout=nvarout+ndiag
      if(out_diagr) nvarout=nvarout+ndiag
      if(out_diagi) nvarout=nvarout+ndiag
      if(out_diags) nvarout=nvarout+ndiag*nscal
      if(out_micro) nvarout=nvarout+ndiag*nhydro
      if(out_micro.and.moments==2) nvarout=nvarout+ndiag*nhydro
#ifdef AERO_ENABLE
      if(out_diaga) nvarout=nvarout+2*ndiag*nmode
#endif    
!
!  2D (surface) variables
!
      if (out_surf) then
        if(out_lwp) nvarsurf=nvarsurf+2
        if(out_cwp) nvarsurf=nvarsurf+2
        if(out_wvp) nvarsurf=nvarsurf+3
	if(out_cmse) nvarsurf=nvarsurf+3
	if(out_cmfl) nvarsurf=nvarsurf+1
        if(out_ctop) nvarsurf=nvarsurf+2
        if(out_zinv) nvarsurf=nvarsurf+2
        if(out_sfl) nvarsurf=nvarsurf+2
        if(out_sst) nvarsurf=nvarsurf+2
        if(out_rrate) nvarsurf=nvarsurf+2
        if(out_srad) nvarsurf=nvarsurf+2
        if(out_olr) nvarsurf=nvarsurf+2
        if(out_osr) nvarsurf=nvarsurf+2
        if(out_ocs) nvarsurf=nvarsurf+2
        if(out_cape) nvarsurf=nvarsurf+4
        if(out_cp) nvarsurf=nvarsurf+2
        if(out_ints) nvarsurf=nvarsurf+1
!#ifdef AERO_RADIA
        if(out_opthic) nvarsurf=nvarsurf+3
!#endif 
      endif
!
!  Profile variables
!
      if (pro_tot) then
        nvarpro=nvarpro+1
	prname(nvarpro) = 'tot'
      endif
      if (pro_env) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'env'
      endif
      if (pro_cl) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'cl'
      endif
      if (pro_cor) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'cor'
      endif    
      if (pro_int) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'int'
      endif 
      if (pro_out) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'out'
      endif 
      if (pro_cfev) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'cfev'
      endif
      if (pro_ent) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'ent'
      endif
      if (pro_det) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'det'
      endif
      if (pro_halo) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'hal'
      endif 
      if (pro_up) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'upd'
      endif 
      if (pro_dn) then
        nvarpro=nvarpro+1
        prname(nvarpro) = 'dnd'
      endif 

      if (verbose > 0) call write_debug('Terminate count_outputs')

      return
      end
