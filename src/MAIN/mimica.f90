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
!  MIMICA.F:                   
!
!  Purpose:
!      MIMICA's Main Program  - LES Model: Version 5.0              
!
!  Author
!      Chien Wang
!      MIT Joint Program on Science and Policy of Global Change
!
! ================================================================

!      =============
      PROGRAM mimica
!      =============

#if (defined SPMD)
      USE mpi
#endif
      USE gridno
      USE shared_data
      USE shared_aerosol_new
      USE, intrinsic :: iso_c_binding 
      
      IMPLICIT NONE
      
      integer :: ierr
      real, dimension(34,300)  :: r_rktable1
      real, dimension( 3,300)  :: r_rktable2
#include "rktable.h"

      character(len = 100):: filename
      
      namelist /cm_run/dt0,ntau,new_run,nest_run,verbose,ldtfix,limit_ts,tstart,tstop

      namelist /cm_out/no_out,iax,iav,its,ipro,isli,ires,inest,minmax,ts_out,file_output,zbl,ztrop,nslicex,nslicey,nslicez,slicex_io,slicey_io,  &
                       slicez_io,spec_diag,out_surf,kout,all_rest,out_hov,out_yav,qthres,wthres

      namelist /cm_grid/dx,dy,dz,xc_nest,yc_nest,lx_nest,ly_nest,ztop,z1,z2,sratio,sponge,                     &
                        zdamp,dxdamp,dydamp,tdamp,tnudg,gridfile,rep

      namelist /cm_init/psurf,dpt,dqv,dw,kpert,j_day,sca_set,ctr_lat,t_local,           &
                        casename,file_init,file_rest,tpert,rep,forcing_file,forcing_file_ssm,forcing_file_shf,forcing_file_lhf,forcing_surface
                        !,moist_var

      namelist /cm_phys/with_mom,with_scal,with_adv,with_dif,with_mic,with_buoy,with_nudg,with_lsadv,with_lssrc,with_lssub,with_cor,with_rad,u0shift,v0shift,Ddiv,w_up,   &
                        pran,diff,anis_k,nl_sgs,cst_cp,iradx,rad_sw,rad_o3,zdec,with_piggy,with_tvar, with_radmod,t_radmod, with_warmadv, tstart1_warmadv, tstop1_warmadv, tstart2_warmadv, tstop2_warmadv, tstart3_warmadv, tstop3_warmadv, &
                        warm_adv_H1,warm_adv_dT1,q_adv_H1,q_adv_dQ1,warm_adv_H2,warm_adv_dT2,q_adv_H2,q_adv_dQ2,warm_adv_H3,warm_adv_dT3,q_adv_H3,q_adv_dQ3
                        
      namelist /cm_num/sorder,scal_adv,cfl_max,cfl_min,limit,lavisc,lim_tol,mom_ord,nsubp,fft_solv,p_mcons,imp_buoy,split_mic,diff_ord,lddamp

      namelist /cm_bc/bcl,sst,ssm,c_dm,c_ds,isurf,momsf,ust,shf0,lhf0,scf,alb,emi,zrough,min_w
      
      namelist /cm_micro/lmicro,lndrop,lvent,ldrizz,ice_habit,micro_dif,lfreeze,qauto,with_dropsed,with_collision,xauto,with_aero_swelling,hgf,     	 &
                         xnc0_d,xnc0_s,xnc0_k,xn_ccn0,xn_in0,ice_delay,auto,moments,lrime,dtcon,ldtmic,alpha_drop,nu_drop

      namelist /cm_aero/laero,aero_flg,nmode0,lkin,reg_mode,aero_sfc_source,aeroi,aeror
      
      namelist /cm_lag/aerosol_lag,compos_lag,mu_lag,sigma_lag,ilag,nlag,res_lag,lag_init,lag_ord,lag_mix,zl1,zl2
!
      namelist /ou_in/out_u, out_v, out_w,             &
                  out_p, out_sca,     		       &
                  out_t, out_pt, out_ptv, out_mse,     &
                  out_rho, out_qv, out_qc, out_nc,     &
                  out_qr, out_nr, out_qi, out_ni,      &
                  out_qg, out_ng, out_qs, out_ns,      &
                  out_dc, out_dr, out_di, out_dg,      &
		  out_qh, out_nh, out_dh, out_wh,      &
		  out_wg, out_ws, out_ds, out_vp,      &
                  out_qt, out_ccn, 	    	       &
                  out_in, out_k, out_tke, out_beff,    &
                  out_dp, out_dtnet, out_frad,         &
                  out_buoy, out_mf,                    &
                  out_sat, out_flut, out_fsgs,         &
                  out_grad, out_ent, out_div, 	       &
		  out_vort, out_z, out_prec,	       & 
		  out_var, out_aero, out_plume
  
      namelist /ou2d_in/out_lwp,out_cwp,out_wvp,out_cmse,out_cmfl,out_ctop,out_sfl,out_rrate,out_srad,out_cape,out_cp,out_ints,&
                   out_olr,out_osr,out_ocs,out_opthic,out_zinv,out_sst

      namelist /oud_in/out_diagu, out_diagv, out_diagw, out_diagp,      &
                   out_diagt, out_diagtv, out_diagq, out_diagl,         &
                   out_diagr, out_diagi, out_diaga, out_diags,          &
                   out_diagk, out_micro

      namelist /pro_in/pro_tot,pro_env,pro_cl,pro_cor,pro_int,pro_out,pro_cfev,pro_ent,pro_det,pro_halo,pro_up,pro_dn
!
! --------------------------------------------------------
!
#include "default.h"

!
! ===  Read Namelists:
! 
      filename = 'cm.nml'
      open(3,file=filename(1:len_trim(filename)),status='old')
      read(3,nml=cm_run)
      read(3,nml=cm_init)
      read(3,nml=cm_grid)
      read(3,nml=cm_out)
      read(3,nml=cm_num)
      read(3,nml=cm_phys)
      read(3,nml=cm_bc)
      read(3,nml=cm_micro)
#ifdef AERO_ENABLE
      read(3,nml=cm_aero)
#endif
#ifdef LAGRANGE
      read(3,nml=cm_lag)
#endif
      close(3)
!
! ===  Open Files:
!
      opdir = './OUTPUT'

      filename = 'out.nml'
      open(3,file=filename(1:len_trim(filename)),status='old')
      read(3,nml=ou_in)
      read(3,nml=ou2d_in)
      read(3,nml=oud_in)
      read(3,nml=pro_in)
      close(3)
      
      ! === adjust time steps
      dt0_i = dt0
!
#ifdef AQCHEM_ENABLE
      if (lmicro < 1) then
        print*, 'ERROR: FOR LIQUID PHASE CHEMISTRY, lmicro must be >= 1'
        stop
      endif
#endif
!
#ifdef SOLIDCHEM_ENABLE
      if (lmicro < 2) then
        print*, 'ERROR: FOR SOLID PHASE CHEMISTRY, lmicro must be >= 2'
        stop
      endif
#endif
!
! ===  Initial Data Input:
!
#ifdef CHEM_ENABLE
      read(12)r_rktable1      !reaction rates
      read(12)r_rktable2
      rktable1(:,:) = r_rktable1(:,:)
      rktable2(:,:) = r_rktable2(:,:)
#endif

!
! === Initialize mimica log
!
      mypid = 0
      if (new_run) then
        how_access = 'rewind'
      else
        how_access = 'append'
      endif
      prtname = opdir(1:len_trim(opdir)) // '/cm.prt'
      call write_to_prt(7,'STARTING MIMICA RUN',how_access)
!
! ===  Start MPI:
! 

#if ( defined SPMD )
!
!  Initialize MPI
!
!      call MPI_INITIALIZED(mpi_running,ierr)
      call MPI_INIT( ierr )

      call MPI_COMM_RANK( MPI_COMM_WORLD, mypid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )
!
!  Check number of processes
!
      if ( mod(nx-5, nprocx) .ne. 0 ) then
        print *,"WRONG NUMBER OF PROCESSES IN X:",nx-5,nprocx,"STOPPING THE RUN"
	call MPI_FINALIZE ( ierr )
        stop
      end if
!
#ifdef DECOMP_2D
      if ( mod(ny-5, nprocy) .ne. 0 ) then
        print *,"WRONG NUMBER OF PROCESSES IN Y:",ny-5,nprocy,"STOPPING THE RUN"
	call MPI_FINALIZE ( ierr )
        stop
      end if
#endif
!      
#else
       nproc=1
       mypid=0
#endif
 
!
! ===  Processing Model Run:
!
      call modelctl     

      end       

