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
!  netcdf_interface.f90                   
!
!  Purpose:
!      Module file containing routines to output netcdf files of various dimensions                 
!
!  Author
!      Julien Savre
!      MIM, Ludwig Maximilian Universitat, Munich
!
! ================================================================

module netcdfdef

!
!-----------------------------------------------------------!
!

  USE gridno
  USE shared_data
  USE netcdf
  
  implicit none
        
  private
  
  character (len = *), parameter :: UNITS = "units", INFO = "description"
  
  public :: add_variables, add_2d_variables, add_lag_variables, register_var, define_indices
  
  contains
  !
  ! ================================================================
  !> Prepare output variables for adding to netCDF file
  !! 
  !! Checks output flags and prepares output
  !! calls "register_var(...)" for each output variable
  !! 
  subroutine add_variables (ndims,ivar,ncid,dimids,varids)
    integer, intent(in) :: ncid               !< ID of netCDF file
  integer, intent(in) :: ndims, ivar          !< Dimentionality of run
  integer, dimension(:), intent(in) :: dimids !< IDs of netCDF dimensions
  integer, dimension(:), intent(inout) :: varids !< IDs of netCDF variables
  !
  integer :: h, k, nvar
  character(len=1) :: car1
  
  nvar = ivar

  if (out_u) call register_var ('u    ', nvar, ncid, dimids, varids)
#ifdef MODEL_3D
  if (out_v) call register_var ('v    ', nvar, ncid, dimids, varids)
#endif  
  if (out_w) call register_var ('w    ', nvar, ncid, dimids, varids)
  if (out_p) call register_var ('p    ', nvar, ncid, dimids, varids)
  if (out_p) call register_var ('ptot ', nvar, ncid, dimids, varids)
  if (out_t) call register_var ('t    ', nvar, ncid, dimids, varids)
  if (out_pt) call register_var ('pt   ', nvar, ncid, dimids, varids)
  if (out_mse) call register_var ('mse  ', nvar, ncid, dimids, varids)
  if (out_ptv) call register_var ('pte  ', nvar, ncid, dimids, varids)
  if (out_rho) call register_var ('rho  ', nvar, ncid, dimids, varids)
  if (out_rho) call register_var ('drho ', nvar, ncid, dimids, varids)
  if (out_qv) call register_var ('qv   ', nvar, ncid, dimids, varids)
  if (out_qt) call register_var ('qt   ', nvar, ncid, dimids, varids)
  if (out_z) call register_var ('ref  ', nvar, ncid, dimids, varids)
  if (lmicro>0.and.out_qc) call register_var ('qc   ', nvar, ncid, dimids, varids)
  if (lmicro>0.and.out_qr) call register_var ('qr   ', nvar, ncid, dimids, varids)
  if (lmicro>1.and.out_qi) call register_var ('qi   ', nvar, ncid, dimids, varids)
  if (lmicro>2.and.out_qg) call register_var ('qg   ', nvar, ncid, dimids, varids)
  if (lmicro>2.and.out_qs) call register_var ('qs   ', nvar, ncid, dimids, varids)
  if (lmicro>3.and.out_qh) call register_var ('qh   ', nvar, ncid, dimids, varids)
  if (lmicro>0.and.out_nc) call register_var ('nc   ', nvar, ncid, dimids, varids)
  if (lmicro>0.and.out_nr) call register_var ('nr   ', nvar, ncid, dimids, varids)
  if (lmicro>1.and.out_ni) call register_var ('ni   ', nvar, ncid, dimids, varids)
  if (lmicro>2.and.out_ng) call register_var ('ng   ', nvar, ncid, dimids, varids)
  if (lmicro>2.and.out_ns) call register_var ('ns   ', nvar, ncid, dimids, varids)
  if (lmicro>3.and.out_nh) call register_var ('nh   ', nvar, ncid, dimids, varids)
  if (lmicro>0.and.out_dc) call register_var ('dc   ', nvar, ncid, dimids, varids)
  if (lmicro>0.and.out_dr) call register_var ('dr   ', nvar, ncid, dimids, varids)
  ! if (lmicro>0.and.out_vp) call register_var ('vp   ', nvar, ncid, dimids, varids)
  if (lmicro>0.and.out_vp) call register_var ('vdrop   ', nvar, ncid, dimids, varids)
  if (lmicro>1.and.out_di) call register_var ('di   ', nvar, ncid, dimids, varids)
  if (lmicro>2.and.out_dg) call register_var ('dg   ', nvar, ncid, dimids, varids)
  if (lmicro>2.and.out_ds) call register_var ('ds   ', nvar, ncid, dimids, varids)
  if (lmicro>3.and.out_dh) call register_var ('dh   ', nvar, ncid, dimids, varids)
  if (lmicro>3.and.out_wg) call register_var ('wg   ', nvar, ncid, dimids, varids)
  if (lmicro>3.and.out_ws) call register_var ('ws   ', nvar, ncid, dimids, varids)
  if (lmicro>3.and.out_wh) call register_var ('wh   ', nvar, ncid, dimids, varids)
  if (lmicro>0.and.out_prec) call register_var ('prel ', nvar, ncid, dimids, varids)
  if (lmicro>1.and.out_prec) call register_var ('prei ', nvar, ncid, dimids, varids)
  if (out_ccn) call register_var ('ccn  ', nvar, ncid, dimids, varids)
  if (out_in) call register_var ('in   ', nvar, ncid, dimids, varids)
  if (out_k) call register_var ('k    ', nvar, ncid, dimids, varids)
  if (out_tke) call register_var ('ke   ', nvar, ncid, dimids, varids)
  if (out_tke) call register_var ('ke2  ', nvar, ncid, dimids, varids)
  if (out_tke) call register_var ('wp2  ', nvar, ncid, dimids, varids)
  if (out_tke) call register_var ('wp3  ', nvar, ncid, dimids, varids)
  if (out_tke) call register_var ('s2   ', nvar, ncid, dimids, varids)
  if (out_tke) call register_var ('n2   ', nvar, ncid, dimids, varids)
  if (out_sat) call register_var ('rh   ', nvar, ncid, dimids, varids)
  if (out_sat) call register_var ('rhi  ', nvar, ncid, dimids, varids)
  if (out_div) call register_var ('div  ', nvar, ncid, dimids, varids)
  if (out_div) call register_var ('qcon ', nvar, ncid, dimids, varids)
  if (out_vort) call register_var ('vox  ', nvar, ncid, dimids, varids)
  if (out_vort) call register_var ('voy  ', nvar, ncid, dimids, varids)
  if (out_vort) call register_var ('voz  ', nvar, ncid, dimids, varids)
  if (out_mf) call register_var ('mf   ', nvar, ncid, dimids, varids)
  if (out_mf) call register_var ('cc   ', nvar, ncid, dimids, varids)
  if (out_mf) call register_var ('mfm  ', nvar, ncid, dimids, varids)
  if (out_mf) call register_var ('ccm  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gu1  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gu2  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gu3  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gv1  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gv2  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gv3  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gw1  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gw2  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gw3  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gt1  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gt2  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gt3  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gq1  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gq2  ', nvar, ncid, dimids, varids)
  if (out_grad) call register_var ('gq3  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('activ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('en1  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('en2  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('en3  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('en4  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('en5  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('en6  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('en7  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('en8  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('de1  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('de2  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('de3  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('de4  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('de5  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('de6  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('de7  ', nvar, ncid, dimids, varids)
  if (out_ent) call register_var ('de8  ', nvar, ncid, dimids, varids)
  if (out_var.and.out_pt) call register_var ('ptv  ', nvar, ncid, dimids, varids)
  if (out_var.and.out_qv) call register_var ('qvv  ', nvar, ncid, dimids, varids)
  if (out_buoy) call register_var ('buoy ', nvar, ncid, dimids, varids)
  if (out_beff) call register_var ('beff ', nvar, ncid, dimids, varids)
  if (out_dp) call register_var ('dpb  ', nvar, ncid, dimids, varids)
  if (out_dp) call register_var ('dpd  ', nvar, ncid, dimids, varids)
  if (out_dp) call register_var ('dpnh ', nvar, ncid, dimids, varids)
  if (out_flut) then
    if (out_u) call register_var ('uw   ', nvar, ncid, dimids, varids)
#ifdef MODEL_3D
    if (out_v) call register_var ('vw   ', nvar, ncid, dimids, varids)  
#endif
    if (out_buoy) call register_var ('bres ', nvar, ncid, dimids, varids)
    if (out_pt) call register_var ('ptw  ', nvar, ncid, dimids, varids)
    if (out_qt) call register_var ('qtw  ', nvar, ncid, dimids, varids)
    if (out_ptv) call register_var ('ptvw ', nvar, ncid, dimids, varids)
    if (out_sca.and.nscal>0) call register_var ('s1tw ', nvar, ncid, dimids, varids)
  endif
  if (out_fsgs) then
    if (out_u) call register_var ('usgs ', nvar, ncid, dimids, varids)  
#ifdef MODEL_3D
    if (out_v) call register_var ('vsgs ', nvar, ncid, dimids, varids)  
#endif
    if (out_buoy) call register_var ('bsgs ', nvar, ncid, dimids, varids)  
    if (out_pt) call register_var ('ptsgs', nvar, ncid, dimids, varids)  
    if (out_qt) call register_var ('qtsgs', nvar, ncid, dimids, varids)  
  endif
  if (out_sca) then
    do h = 1, nscal
      write(car1,'(i1)') h
      call register_var ('sca  ', nvar, ncid, dimids, varids, h=h)
    enddo
  endif
  if (out_dtnet) call register_var ('dtn  ', nvar, ncid, dimids, varids)
#ifdef RAD_ENABLE
  if (out_dtnet) call register_var ('dtsw  ', nvar, ncid, dimids, varids)
  if (out_dtnet) call register_var ('dtlw  ', nvar, ncid, dimids, varids)

  if (out_frad) call register_var ('fnet ', nvar, ncid, dimids, varids)
  if (out_frad) call register_var ('fsw  ', nvar, ncid, dimids, varids)
  if (out_frad) call register_var ('flw  ', nvar, ncid, dimids, varids)  
#endif    

#ifdef AERO_ENABLE
  if (out_aero) then  
    do h = 1, nmode
      call register_var ('an   ', nvar, ncid, dimids, varids, h=h)
      call register_var ('am   ', nvar, ncid, dimids, varids, h=h)
      call register_var ('ama  ', nvar, ncid, dimids, varids, h=h)
    enddo
  endif
#endif

#ifdef NUC_CNT
  if (out_in) then  
    call register_var ('dif  ', nvar, ncid, dimids, varids)
    call register_var ('dff  ', nvar, ncid, dimids, varids)
    call register_var ('dtf  ', nvar, ncid, dimids, varids)
    call register_var ('cif  ', nvar, ncid, dimids, varids)
    call register_var ('cff  ', nvar, ncid, dimids, varids)
    call register_var ('ctf  ', nvar, ncid, dimids, varids)
    call register_var ('dic  ', nvar, ncid, dimids, varids)
    call register_var ('dfc  ', nvar, ncid, dimids, varids)
    call register_var ('dtc  ', nvar, ncid, dimids, varids)
    call register_var ('cic  ', nvar, ncid, dimids, varids)
    call register_var ('cfc  ', nvar, ncid, dimids, varids)
    call register_var ('ctc  ', nvar, ncid, dimids, varids)
  endif
#endif
!
  if (ldiag) then
    do h = 1, ndiag
      if (out_diagu) call register_var ('udg  ', nvar, ncid, dimids, varids, h)
      if (out_diagv) call register_var ('vdg  ', nvar, ncid, dimids, varids, h)
      if (out_diagw) call register_var ('wdg  ', nvar, ncid, dimids, varids, h)
      if (out_diagp) call register_var ('pdg  ', nvar, ncid, dimids, varids, h)
      if (out_diagt) call register_var ('tdg  ', nvar, ncid, dimids, varids, h)
      if (out_diagtv) call register_var ('tvdg ', nvar, ncid, dimids, varids, h)
      if (out_diagq) call register_var ('qdg  ', nvar, ncid, dimids, varids, h)
      if (out_diagk) call register_var ('kdg  ', nvar, ncid, dimids, varids, h)
      if (lmicro>0.and.out_diagl) call register_var ('ldg  ', nvar, ncid, dimids, varids, h)
      if (lmicro>0.and.out_diagr) call register_var ('rdg  ', nvar, ncid, dimids, varids, h)
      if (lmicro>1.and.out_diagi) call register_var ('idg  ', nvar, ncid, dimids, varids, h)
      if (lmicro>0.and.out_micro) then
        do k = 1, nhydro
          write(car1,'(i1)')k
          call register_var ('mdq'//car1//' ', nvar, ncid, dimids, varids, h)
          if(moments==2) call register_var ('mdn'//car1//' ', nvar, ncid, dimids, varids, h)
        enddo
      endif
      do k = 1, nscal
        write(car1,'(i1)')k
        if (out_diags) call register_var ('sdg'//car1//' ', nvar, ncid, dimids, varids, h)
      enddo
#ifdef AERO_ENABLE
      do k = 1, nmode
        write(car1,'(i1)')k
        if (out_diaga) call register_var ('ndg'//car1//' ', nvar, ncid, dimids, varids, h)
        if (out_diaga) call register_var ('mdg'//car1//' ', nvar, ncid, dimids, varids, h)
      enddo
#endif
    enddo
  endif
  
  end subroutine
      
  subroutine add_2d_variables (ndims,ivar,ncid,dimids,varids)
    
  integer :: nvar,h
  integer, intent(in) :: ncid,ndims,ivar
  integer, dimension(:), intent(in) :: dimids
  integer, dimension(:), intent(inout) :: varids
  character(len=1) :: car1

  nvar = ivar
  
  if (out_lwp) call register_var ('lwp  ', nvar, ncid, dimids, varids)
  if (out_lwp) call register_var ('iwp  ', nvar, ncid, dimids, varids)
  if (out_cwp) call register_var ('cwp  ', nvar, ncid, dimids, varids)
  if (out_cwp) call register_var ('rwp  ', nvar, ncid, dimids, varids)  
  if (out_wvp) call register_var ('wvp  ', nvar, ncid, dimids, varids)
  if (out_wvp) call register_var ('wvpt ', nvar, ncid, dimids, varids)
  if (out_wvp) call register_var ('wvpb ', nvar, ncid, dimids, varids)
  if (out_cmse) call register_var ('cmse ', nvar, ncid, dimids, varids)
  if (out_cmse) call register_var ('cmset', nvar, ncid, dimids, varids)
  if (out_cmse) call register_var ('cmseb', nvar, ncid, dimids, varids)
  if (out_cmfl) call register_var ('cmfl ', nvar, ncid, dimids, varids)
  if (out_cape) call register_var ('cape ', nvar, ncid, dimids, varids)
  if (out_cape) call register_var ('cin  ', nvar, ncid, dimids, varids)
  if (out_cape) call register_var ('lcl  ', nvar, ncid, dimids, varids)
  if (out_cape) call register_var ('lfc  ', nvar, ncid, dimids, varids)
  if (out_cp) call register_var ('cpi0 ', nvar, ncid, dimids, varids)
  if (out_cp) call register_var ('cpi  ', nvar, ncid, dimids, varids)
  if (out_ctop) call register_var ('ctop ', nvar, ncid, dimids, varids)
  if (out_ctop) call register_var ('cbas ', nvar, ncid, dimids, varids)
  if (out_zinv) call register_var ('bltop', nvar, ncid, dimids, varids)
  if (out_zinv) call register_var('zinv ', nvar, ncid, dimids, varids)
  if (out_ints) call register_var ('isca ', nvar, ncid, dimids, varids)
  if (out_rrate) call register_var ('rrate', nvar, ncid, dimids, varids)
  if (out_rrate) call register_var ('cumul', nvar, ncid, dimids, varids)
  if (out_sfl) call register_var ('shf  ', nvar, ncid, dimids, varids)
  if (out_sfl) call register_var ('lhf  ', nvar, ncid, dimids, varids)
  if (out_sst) call register_var ('sst  ', nvar, ncid, dimids, varids)
  if (out_sst) call register_var ('ssq  ', nvar, ncid, dimids, varids)
  if (out_srad) call register_var ('slw  ', nvar, ncid, dimids, varids)
  if (out_srad) call register_var ('ssw  ', nvar, ncid, dimids, varids)
  if (out_olr) call register_var ('olr  ', nvar, ncid, dimids, varids)
  if (out_olr) call register_var ('toa  ', nvar, ncid, dimids, varids)
  if (out_osr) call register_var ('osr  ', nvar, ncid, dimids, varids)
  if (out_osr) call register_var ('isr  ', nvar, ncid, dimids, varids)  
  if (out_ocs) call register_var ('olrcs', nvar, ncid, dimids, varids)
  if (out_ocs) call register_var ('osrcs', nvar, ncid, dimids, varids)

  #ifdef AERO_RADIA
  if (out_opthic) call register_var('opta2', nvar, ncid, dimids, varids)
  if (out_opthic) call register_var('optc2', nvar, ncid, dimids, varids)  
  if (out_opthic) call register_var('optr2', nvar, ncid, dimids, varids)
  #endif
  end subroutine
      
  subroutine add_lag_variables (ncid,dimids,varids)
  
  integer :: h, nvar
  integer, intent(in) :: ncid
  integer, dimension(:), intent(in) :: dimids
  integer, dimension(:), intent(inout) :: varids
  character(len=1) :: car1

  nvar = 2  

  call register_var('x    ',nvar, ncid, dimids, varids)
  call register_var('y    ',nvar, ncid, dimids, varids)
  call register_var('z    ',nvar, ncid, dimids, varids)
  call register_var('u    ',nvar, ncid, dimids, varids)
  call register_var('v    ',nvar, ncid, dimids, varids)
  call register_var('w    ',nvar, ncid, dimids, varids)
  call register_var('k    ',nvar, ncid, dimids, varids)
  call register_var('t    ',nvar, ncid, dimids, varids)
  call register_var('pt   ',nvar, ncid, dimids, varids)
  call register_var('qt   ',nvar, ncid, dimids, varids)
  call register_var('buoy ',nvar, ncid, dimids, varids)
  if (lmicro > 0) call register_var('ql   ',nvar, ncid, dimids, varids)
  if (lmicro > 1) call register_var('qi   ',nvar, ncid, dimids, varids)
  if (nscal > 0) then
    do h = 1, nscal
      write(car1,'(i1)') h
      call register_var ('sca  ', nvar, ncid, dimids, varids, h=h)
    enddo
  endif
  if (out_diagw) then
    do h = 1, 9
      call register_var ('wdg  ', nvar, ncid, dimids, varids, h)
    enddo
  endif
  if (out_diagl) then
    do h = 1, 9
      call register_var ('ldg  ', nvar, ncid, dimids, varids, h)
    enddo
  endif

  end subroutine

  ! ================================================================
  !> Adding variable to netCDF file
  !! 
  !! Adds variables to netcdf and addes attribuites according to information
  !! from "return_infor(...)"
  subroutine register_var (var,nvar,ncid,dimids,varids,h)
  
    character(len=5) :: var
    integer, intent(inout) :: nvar
    integer, intent(in) :: ncid
    integer, optional, intent(in) :: h
    integer, dimension(:), intent(in) :: dimids
    integer, dimension(:), intent(inout) :: varids
    integer :: l
    character(len=100) :: varunit,varname,varinfo

    if (present(h)) then
      l = h
    else
      l = 0
    endif

    nvar = nvar+1

    call return_info (var,varunit,varname,varinfo,l)
    call check( nf90_def_var(ncid, trim(varname), nf90_real, dimids, varids(nvar)), varname )
    call check( nf90_put_att(ncid, varids(nvar), UNITS, trim(varunit)), varname )
    call check( nf90_put_att(ncid, varids(nvar), INFO, trim(varinfo)), varname )
  
  contains
  
    subroutine return_info (var,units,name,info,l)
    
    implicit none
    
    integer :: l, h
    character(len=*), intent(out), optional :: units,name,info
    character(len=1) :: car1
    character(len=2) :: car2
    character(len=5) :: var
    
    write(car1,'(a1)') var
    if (l /= 0) then
      if ( l < 10 ) then
        write(car1,'(i1)') l
	car2 = car1//' '
      else
        write(car2,'(i2)') l
      endif
    endif
    
    select case ( var )
      case('x    ')
        if (present(name)) name = 'X'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Parcel position in X'
      case('y    ')
        if (present(name)) name = 'Y'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Parcel position in Y'
      case('z    ')
        if (present(name)) name = 'Z'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Parcel position in Z'
      case('u    ')
        if (present(name)) name = 'U'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Zonal velocity U'
      case('v    ')
        if (present(name)) name = 'V'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Meridional velocity V'
      case('w    ')
        if (present(name)) name = 'W'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Vertical velocity W'
      case('p    ')
        if (present(name)) name = 'P_pert'
        if (present(units)) units = 'Pa'
	if (present(info)) info = 'Perturbation pressure'
      case('ptot ')
        if (present(name)) name = 'P_tot'
        if (present(units)) units = 'Pa'
	if (present(info)) info = 'Total pressure (hydrostatic + perturbation)'
      case('t    ')
        if (present(name)) name = 'T'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Absolute temperature'
      case('pt   ')
        if (present(name)) name = 'PT'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Potential temperature'
      case('pte  ')
        if (present(name)) name = 'PTv'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Virtual potential temperature'
      case('mse  ')
        if (present(name)) name = 'MSE'
        if (present(units)) units = 'J'
	if (present(info)) info = 'Moist static energy'
      case('drho ')
        if (present(name)) name = 'RHO_pert'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Perturbation density'
      case('rho  ')
        if (present(name)) name = 'RHO'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Full density'
      case('qv   ')
        if (present(name)) name = 'Qv'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Water vapour mixing fraction'
      case('qt   ')
        if (present(name)) name = 'Qt'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Total water mixing fraction'
      case('ref  ')
        if (present(name)) name = 'Ze'
        if (present(units)) units = 'dBz'
	if (present(info)) info = 'Equivalent radar reflectivity'
      case('ql   ')
        if (present(name)) name = 'Ql'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Liquid water mixing fraction'
      case('qc   ')
        if (present(name)) name = 'Qc'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Cloud droplet mixing fraction'
      case('qr   ')
        if (present(name)) name = 'Qr'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Rain drops mixing fraction'
      case('qi   ')
        if (present(name)) name = 'Qi'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Ice mixing fraction'
      case('qg   ')
        if (present(name)) name = 'Qg'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Graupel mixing fraction'
      case('qs   ')
        if (present(name)) name = 'Qs'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Snow mixing fraction'
      case('qh   ')
        if (present(name)) name = 'Qh'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Hail mixing fraction'
      case('nc   ')
        if (present(name)) name = 'Nc'
        if (present(units)) units = '1/kg'
	if (present(info)) info = 'Cloud droplet number density'
      case('nr   ')
        if (present(name)) name = 'Nr'
        if (present(units)) units = '1/kg'
	if (present(info)) info = 'Rain drops number density'
      case('ni   ')
        if (present(name)) name = 'Ni'
        if (present(units)) units = '1/kg'
	if (present(info)) info = 'Ice number density'
      case('ng   ')
        if (present(name)) name = 'Ng'
        if (present(units)) units = '1/kg'
	if (present(info)) info = 'Graupel number density'
      case('ns   ')
        if (present(name)) name = 'Ns'
        if (present(units)) units = '1/kg'
	if (present(info)) info = 'Snow number density'
      case('nh   ')
        if (present(name)) name = 'Nh'
        if (present(units)) units = '1/kg'
	if (present(info)) info = 'Hail number density'
      case('dc   ')
        if (present(name)) name = 'Dc'
        if (present(units)) units = 'mm'
	if (present(info)) info = 'Cloud droplet mean size'
      case('dr   ')
        if (present(name)) name = 'Dr'
        if (present(units)) units = 'mm'
	if (present(info)) info = 'Rain drops mean size'
      case('di   ')
        if (present(name)) name = 'Di'
        if (present(units)) units = 'mm'
	if (present(info)) info = 'Ice mean size'
      case('dg   ')
        if (present(name)) name = 'Dg'
        if (present(units)) units = 'mm'
	if (present(info)) info = 'Graupel mean size'
      case('ds   ')
        if (present(name)) name = 'Ds'
        if (present(units)) units = 'mm'
	if (present(info)) info = 'Snow mean size'
      case('dh   ')
        if (present(name)) name = 'Dh'
        if (present(units)) units = 'mm'
	if (present(info)) info = 'Hail mean size'
      case('wg   ')
        if (present(name)) name = 'Wg'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Rimed mass of water in graupel'
      case('ws   ')
        if (present(name)) name = 'Ws'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Rimed mass of water in graupel'
      case('wh   ')
        if (present(name)) name = 'Wh'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Rimed mass of water in hail'
  !     case('vp   ')
  !       if (present(name)) name = 'Vp'
  !       if (present(units)) units = 'm/s'
	! if (present(info)) info = 'Rain precipitation velocity'
      case('vdrop ')
        if (present(name)) name = 'Vdrop'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Droplet sedimentation velocity'
      case('prel ')
        if (present(name)) name = 'PREC_RAIN'
        if (present(units)) units = 'kg/m2/s'
	if (present(info)) info = 'Liquid (rain) precipitation flux'
      case('prei ')
        if (present(name)) name = 'PREC_ICE'
        if (present(units)) units = 'kg/m2/s'
	if (present(info)) info = 'Ice precipitation flux'
      case('rh   ')
        if (present(name)) name = 'RH'
        if (present(units)) units = '%'
	if (present(info)) info = 'Relative humidity'
      case('rhi  ')
        if (present(name)) name = 'RHi'
        if (present(units)) units = '%'
	if (present(info)) info = 'Relative humidity (over ice)'
      case('ccn  ')
        if (present(name)) name = 'CCN'
        if (present(units)) units = '1/m3'
	if (present(info)) info = 'Bulk total CCN number density'
      case('in   ')
        if (present(name)) name = 'IN'
        if (present(units)) units = '1/m3'
	if (present(info)) info = 'Bulk IN number density'
      case('k    ')
        if (present(name)) name = 'Ksgs'
        if (present(units)) units = 'm2/s'
	if (present(info)) info = 'Horizontal SGS eddy viscosity'
      case('ke   ')
        if (present(name)) name = 'TKE'
        if (present(units)) units = 'm2/s2'
	if (present(info)) info = 'Resolved turbulent kinetic energy'
      case('ke2  ')
        if (present(name)) name = 'TKE_sgs'
        if (present(units)) units = 'm2/s2'
	if (present(info)) info = 'SGS turbulent kinetic energy'
      case('div  ')
        if (present(name)) name = 'Divergence'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'Horizontal divergence'
      case('qcon ')
        if (present(name)) name = 'Q_conv'
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Horizontal moisture convergence'
      case('vox  ')
        if (present(name)) name = 'Vorticity_x'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'Vorticity component x'
      case('voy  ')
        if (present(name)) name = 'Vorticity_y'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'Vorticity component y'
      case('voz  ')
        if (present(name)) name = 'Vorticity_z'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'Vorticity component z'
      case('gu1  ')
        if (present(name)) name = 'Grad_U_x'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'U velocity gradient in the x direction'
      case('gu2  ')
        if (present(name)) name = 'Grad_U_y'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'U velocity gradient in the y direction'
      case('gu3  ')
        if (present(name)) name = 'Grad_U_z'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'U velocity gradient in the z direction'
      case('gv1  ')
        if (present(name)) name = 'Grad_V_x'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'V velocity gradient in the x direction'
      case('gv2  ')
        if (present(name)) name = 'Grad_V_y'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'V velocity gradient in the y direction'
      case('gv3  ')
        if (present(name)) name = 'Grad_V_z'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'V velocity gradient in the z direction'
      case('gw1  ')
        if (present(name)) name = 'Grad_W_x'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'W velocity gradient in the x direction'
      case('gw2  ')
        if (present(name)) name = 'Grad_W_y'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'W velocity gradient in the y direction'
      case('gw3  ')
        if (present(name)) name = 'Grad_W_z'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'W velocity gradient in the z direction'
      case('gt1  ')
        if (present(name)) name = 'Grad_PT_x'
        if (present(units)) units = 'K/s'
	if (present(info)) info = 'PT gradient in the x direction'
      case('gt2  ')
        if (present(name)) name = 'Grad_PT_y'
        if (present(units)) units = 'K/s'
	if (present(info)) info = 'PT gradient in the y direction'
      case('gt3  ')
        if (present(name)) name = 'Grad_PT_z'
        if (present(units)) units = 'K/s'
	if (present(info)) info = 'PT gradient in the z direction'
      case('gq1  ')
        if (present(name)) name = 'Grad_Qt_x'
        if (present(units)) units = 'kg/kg/s'
	if (present(info)) info = 'QT gradient in the x direction'
      case('gq2  ')
        if (present(name)) name = 'Grad_Qt_y'
        if (present(units)) units = 'kg/kg/s'
	if (present(info)) info = 'QT gradient in the y direction'
      case('gq3  ')
        if (present(name)) name = 'Grad_Qt_z'
        if (present(units)) units = 'kg/kg/s'
	if (present(info)) info = 'QT gradient in the z direction'
      case('s2   ')
        if (present(name)) name = 'S2'
        if (present(units)) units = 'm2/s2'
	if (present(info)) info = 'Squared rate of strain tensor'
      case('activ')
        if (present(name)) name = 'Activity'
        if (present(units)) units = '-'
	if (present(info)) info = 'Activity tracer'
      case('de1  ')
        if (present(name)) name = 'Det_tend'
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Detrainment rate from tendencies'
      case('de2  ')
        if (present(name)) name = 'Ent_buoy'
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Buoyancy contribution to detrainment rate'
      case('de3  ')
        if (present(name)) name = 'Ent_pres'
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Pressure gradient contribution to detrainment rate'
      case('de4  ')
        if (present(name)) name = 'Ent_subw'
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Subgrid transport contribution to detrainment rate'
      case('de5  ')
        if (present(name)) name = 'Ent_totalw'
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Total w tendency contribution to detrainmnent rate'
      case('de6  ')
        if (present(name)) name = 'Ent_advw'
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'W advection contribution to detrainment rate'
      case('de7  ')
        if (present(name)) name = 'Det_separate'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Detrainment rate from separate w and q contributions'
      case('de8  ')
        if (present(name)) name = 'Det_activity'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Detrainment rate from activity tracer'
      case('en1  ')
        if (present(name)) name = 'Ent_tend'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Entrainment rate from tendencies'
      case('en2  ')
        if (present(name)) name = 'Ent_evap'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Evaporation/condensation contribution to entrainment/detrainment'
      case('en3  ')
        if (present(name)) name = 'Ent_micro'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Microphysics contribution to entrainment/detrainment'
      case('en4  ')
        if (present(name)) name = 'Ent_subq'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Subgrid transport contribution to entrainment/detrainment'
      case('en5  ')
        if (present(name)) name = 'Ent_totalq'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Total q tendency contribution to entrainment/detrainment'
      case('en6  ')
        if (present(name)) name = 'Ent_advq'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Q advection contribution to entrainment/detrainment'
      case('en7  ')
        if (present(name)) name = 'Ent_separate'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Entrainment rate from separate w and q contributions'
      case('en8  ')
        if (present(name)) name = 'Ent_activity'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Entrainment rate from activity tracer'
      case('wp2  ')
        if (present(name)) name = 'W_var'
        if (present(units)) units = 'm2/s2'
	if (present(info)) info = 'Vertical velocity variance'
      case('wp3  ')
        if (present(name)) name = 'W_skew'
        if (present(units)) units = 'm3/s3'
	if (present(info)) info = 'Vertical velocity skewness'
      case('ptv  ')
        if (present(name)) name = 'Pt_var'
        if (present(units)) units = 'K2'
	if (present(info)) info = 'Potential temperature variance'
      case('qvv  ')
        if (present(name)) name = 'Qv_var'
        if (present(units)) units = '(kg/m3)2'
	if (present(info)) info = 'Water vapour variance'
      case('ptw  ')
        if (present(name)) name = 'WPt_flux'
        if (present(units)) units = 'kg.K/m2/s'
	if (present(info)) info = 'Potential temperature vertical flux'
      case('ptvw ')
        if (present(name)) name = 'WPtv_flux'
        if (present(units)) units = 'kg.K/m2/s'
	if (present(info)) info = 'Virtual potential temperature vertical flux'
      case('qtw  ')
        if (present(name)) name = 'WQt_flux'
        if (present(units)) units = 'kg/m2/s'
	if (present(info)) info = 'Total water vertical flux'
      case('s1tw ')
        if (present(name)) name = 'WSCA1_flux'
        if (present(units)) units = 'kg/m2/s'
	if (present(info)) info = 'Scalar vertical flux'
      case('buoy ')
        if (present(name)) name = 'Buoy'
        if (present(units)) units = 'm2/s3'
	if (present(info)) info = 'Resolved buoyancy'
      case('eva  ')
        if (present(name)) name = 'PT_tend_4'
        if (present(units)) units = 'K/s'
	if (present(info)) info = 'Evaporation cooling'
      case('uw   ')
        if (present(name)) name = 'UW_flux'
        if (present(units)) units = 'kg/m/s2'
	if (present(info)) info = 'Horizontal momentum vertical flux'
      case('vw   ')
        if (present(name)) name = 'VW_flux'
        if (present(units)) units = 'kg/m/s2'
	if (present(info)) info = 'Horizontal momentum vertical flux'
      case('bres ')
        if (present(name)) name = 'WB_flux'
        if (present(units)) units = 'kg/m/s3'
	if (present(info)) info = 'Resolved buoyancy vertical flux'
      case('usgs ')
        if (present(name)) name = 'UWsgs_flux'
        if (present(units)) units = 'kg/m/s2'
	if (present(info)) info = 'Sub-grid vertical U flux (parameterized)'
      case('vsgs ')
        if (present(name)) name = 'VWsgs_flux'
        if (present(units)) units = 'kg/m/s2'
	if (present(info)) info = 'Sub-grid vertical V flux (parameterized)'
      case('bsgs ')
        if (present(name)) name = 'WBsgs_flux'
        if (present(units)) units = 'kg/m/s3'
	if (present(info)) info = 'Sub-grid vertical buoyancy flux (parameterized)'
      case('ptsgs')
        if (present(name)) name = 'PTsgs_flux'
        if (present(units)) units = 'kg K/m2/s'
	if (present(info)) info = 'Sub-grid vertical potential temperature flux (parameterized)'
      case('qtsgs')
        if (present(name)) name = 'QTsgs_flux'
        if (present(units)) units = 'kg/m2/s'
	if (present(info)) info = 'Sub-grid vertical total water flux (parameterized)'
      case('beff ')
        if (present(name)) name = 'Beff'
        if (present(units)) units = 'm2/s3'
	if (present(info)) info = 'Effective buoyancy'
      case('dpb  ')
        if (present(name)) name = 'P_b'
        if (present(units)) units = 'Pa'
	if (present(info)) info = 'Buoyant pressure perturbation'
      case('dpd  ')
        if (present(name)) name = 'P_d'
        if (present(units)) units = 'Pa'
	if (present(info)) info = 'Dynamic pressure perturbation'
      case('dpnh ')
        if (present(name)) name = 'P_nh'
        if (present(units)) units = 'Pa'
	if (present(info)) info = 'Non-hydrostatic pressure perturbation'
      case('n2   ')
        if (present(name)) name = 'N2'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'Buoyancy frequency'
      case('ri   ')
        if (present(name)) name = 'Ri'
        if (present(units)) units = '-'
	if (present(info)) info = 'Richardson number'
      case('sca  ')
        if (present(name)) name = 'SCA'//car1
        if (present(units)) units = '-'
	if (present(info)) info = 'Passive scalar number '//car1
      case('dtn  ')
        if (present(name)) name = 'DT_RAD'
        if (present(units)) units = 'K/s'
	if (present(info)) info = 'Net radiative heating'
      case('fnet ')
        if (present(name)) name = 'FRAD_net'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Net radiative flux'
#ifdef RAD_ENABLE
      case('fsw  ')
        if (present(name)) name = 'FRAD_sw'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Short wave radiative flux'
      case('flw  ')
        if (present(name)) name = 'FRAD_lw'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Long wave radiative flux'

      case('dtsw  ')
        if (present(name)) name = 'DT_SW'
        if (present(units)) units = 'K/s'
  if (present(info)) info = 'Short wave radiative heating'
  
      case('dtlw  ')
        if (present(name)) name = 'DT_LW'
        if (present(units)) units = 'K/s'
  if (present(info)) info = 'Long wave radiative heating'

#endif

!#ifdef AERO_RADIA
    case('opta3')
      if (present(name)) name = 'OPTAER3D'
      if (present(units)) units = '-'
  if (present(info)) info = 'Aerosol optical depth'


    case('optc3')
      if (present(name)) name = 'OPTC3D'
      if (present(units)) units = '-'
  if (present(info)) info = 'Cloud optical depth'
  
  case('optr3')
    if (present(name)) name = 'OPTR3D'
    if (present(units)) units = '-'
  if (present(info)) info = 'Rain optical depth' 
!#endif   

#ifdef AERO_ENABLE
      case('an   ')
        if (present(name)) name = 'N_AERO'//car1
        if (present(units)) units = '1/kg'
	if (present(info)) info = 'Number concentration of aerosols in mode '//car1
      case('am   ')
        if (present(name)) name = 'M_AERO'//car1
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Mass concentration of aerosols in mode '//car1
      case('ama  ')
        if (present(name)) name = 'MA_AERO'//car1
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Mass concentration of activated aerosols in mode '//car1
#endif
#ifdef NUC_CNT
      case('dif  ')
        if (present(name)) name = 'IN_fine_dust_imm'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Immersed fine dust concentration'
      case('dff  ')
        if (present(name)) name = 'IN_fine dust_froz'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Frozen fine dust concentration'
      case('dtf  ')
        if (present(name)) name = 'IN_fine_dust_tot'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Total fine dust concentration'
      case('cif  ')
        if (present(name)) name = 'IN_fine_bc_imm'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Immersed fine black carbon concentration'
      case('cff  ')
        if (present(name)) name = 'IN_fine_bc_froz'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Frozen fine black carbon concentration'
      case('ctf  ')
        if (present(name)) name = 'IN_fine_bc_tot'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Total fine black carbon concentration'
      case('dic  ')
        if (present(name)) name = 'IN_coarse_dust_imm'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Immersed coarse dust concentration'
      case('dfc  ')
        if (present(name)) name = 'IN_coarse_dust_froz'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Frozen coarse dust concentration'
      case('dtc  ')
        if (present(name)) name = 'IN_coarse_dust_tot'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Total coarse dust concentration'
      case('cic  ')
        if (present(name)) name = 'IN_coarse_bc_imm'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Immersed coarse black carbon concentration'
      case('cfc  ')
        if (present(name)) name = 'IN_coarse_bc_froz'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Frozen coarse black carbon concentration'
      case('ctc  ')
        if (present(name)) name = 'IN_coarse_bc_tot'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Total coarse black carbon concentration'
#endif
      case('udg  ')
        if (present(name)) name = 'U_tend_'//trim(car2)
        if (present(units)) units = 'm/s2'
	if (present(info)) info = 'Zonal velocity tendency '//trim(car2)
      case('vdg  ')
        if (present(name)) name = 'V_tend_'//trim(car2)
        if (present(units)) units = 'm/s2'
	if (present(info)) info = 'Meridional velocity tendency '//trim(car2)
      case('wdg  ')
        if (present(name)) name = 'W_tend_'//trim(car2)
        if (present(units)) units = 'm/s2'
	if (present(info)) info = 'Vertical velocity tendency '//trim(car2)
      case('pdg  ')
        if (present(name)) name = 'P_tend_'//trim(car2)
        if (present(units)) units = 'Pa/s'
	if (present(info)) info = 'Pressure tendency '//trim(car2)
      case('tdg  ')
        if (present(name)) name = 'PT_tend_'//trim(car2)
        if (present(units)) units = 'K/s'
	if (present(info)) info = 'Potential temperature tendency '//trim(car2)
      case('tvdg ')
        if (present(name)) name = 'PTV_tend_'//trim(car2)
        if (present(units)) units = 'K/s'
	if (present(info)) info = 'Virtual potential temperature tendency '//trim(car2)
      case('qdg  ')
        if (present(name)) name = 'QT_tend_'//trim(car2)
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Total water mixing fraction tendency '//trim(car2)
      case('kdg  ')
        if (present(name)) name = 'K_tend_'//trim(car2)
        if (present(units)) units = 'm2/s3'
	if (present(info)) info = 'Kinetic energy tendency '//trim(car2)
      case('sdg1 ','sdg2 ','sdg3 ','sdg4 ','sdg5 ','sdg6 ','sdg7 ','sdg8 ','sdg9 ')
        if (present(name)) name = 'SCAL'//var(4:4)//'_tend_'//trim(car2)
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Tracer tendency '//trim(car2)
      case('ndg1 ','ndg2 ','ndg3 ','ndg4 ','ndg5 ','ndg6 ','ndg7 ','ndg8 ','ndg9 ')
        if (present(name)) name = 'NAERO'//var(4:4)//'_tend_'//trim(car2)
        if (present(units)) units = '#/kg/s'
	if (present(info)) info = 'Aerosol tendencies mode '//trim(car2)
      case('mdg1 ','mdg2 ','mdg3 ','mdg4 ','mdg5 ','mdg6 ','mdg7 ','mdg8 ','mdg9 ')
        if (present(name)) name = 'MAAERO'//var(4:4)//'_tend_'//trim(car2)
        if (present(units)) units = 'kg/kg/s'
	if (present(info)) info = 'Aerosol tendencies mode '//trim(car2)
      case('ldg  ')
        if (present(name)) name = 'QC_tend_'//trim(car2)
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Cloud droplet mixing fraction tendency '//trim(car2)
      case('rdg  ')
        if (present(name)) name = 'QR_tend_'//trim(car2)
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Rain drop mixing fraction tendency '//trim(car2)
      case('idg  ')
        if (present(name)) name = 'QI_tend_'//trim(car2)
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Ice crystal mixing fraction tendency '//trim(car2)
      case('mdq1 ','mdq2 ','mdq3 ','mdq4 ','mdq5 ','mdq6 ','mdq7 ','mdq8 ','mdq9 ')
        if (present(name)) name = 'Q'//var(4:4)//'_micro_'//trim(car2)
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Mass microphysical tendency '//trim(car2)//' for hydrometeor type '//var(4:4)
      case('mdn1 ','mdn2 ','mdn3 ','mdn4 ','mdn5 ','mdn6 ','mdn7 ','mdn8 ','mdn9 ')
        if (present(name)) name = 'N'//var(4:4)//'_micro_'//trim(car2)
        if (present(units)) units = 'kg/m3/s'
	if (present(info)) info = 'Number microphysical tendency '//trim(car2)//' for hydrometeor type '//var(4:4)
      case('cid  ')
        if (present(name)) name = 'CLOUD_ID'
        if (present(units)) units = '#'
	if (present(info)) info = 'Cloud ID (integer)'
      case('uid  ')
        if (present(name)) name = 'UPDRAFT_ID'
        if (present(units)) units = '#'
	if (present(info)) info = 'Updraft ID (integer)'
      case('lwp  ')
        if (present(name)) name = 'LWP'
        if (present(units)) units = 'kg/m2'
	if (present(info)) info = 'Liquid Water Path'
      case('iwp  ')
        if (present(name)) name = 'IWP'
        if (present(units)) units = 'kg/m2'
	if (present(info)) info = 'Ice Water Path'
      case('cwp  ')
  if (present(name)) name = 'CWP'
  if (present(units)) units = 'kg/m2'
if (present(info)) info = 'Cloud Water Path'
      case('rwp  ')
  if (present(name)) name = 'RWP'
  if (present(units)) units = 'kg/m2'
if (present(info)) info = 'Rain Water Path'
      case('wvp  ')
        if (present(name)) name = 'WVP'
        if (present(units)) units = 'kg/m2'
	if (present(info)) info = 'Column integrated water vapor'
      case('wvpt ')
        if (present(name)) name = 'WVP_TROP'
        if (present(units)) units = 'kg/m2'
	if (present(info)) info = 'Column integrated water vapor'
      case('wvpb ')
        if (present(name)) name = 'WVP_BL'
        if (present(units)) units = 'kg/m2'
	if (present(info)) info = 'Column integrated water vapor'
      case('cmse ')
        if (present(name)) name = 'CMSE'
        if (present(units)) units = 'J/m2'
	if (present(info)) info = 'Column integrated moist static energy'
      case('cmset')
        if (present(name)) name = 'CMSE_TROP'
        if (present(units)) units = 'J/m2'
	if (present(info)) info = 'Column integrated moist static energy'
      case('cmseb')
        if (present(name)) name = 'CMSE_BL'
        if (present(units)) units = 'J/m2'
	if (present(info)) info = 'Column integrated moist static energy'
      case('cmfl ')
        if (present(name)) name = 'CMFL'
        if (present(units)) units = 'kg m/s'
	if (present(info)) info = 'Column integrated cloud mass flux'
      case('cape ')
        if (present(name)) name = 'CAPE'
        if (present(units)) units = 'J/kg'
	if (present(info)) info = 'Convective Available Potential Energy' 
      case('cin  ')
        if (present(name)) name = 'CIN'
        if (present(units)) units = 'J/kg'
	if (present(info)) info = 'Convective Inhibition'
      case('lcl  ')
        if (present(name)) name = 'LCL'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Lifted Condensation Level'
      case('lfc  ')
        if (present(name)) name = 'LFC'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Level of Free Convection'
      case('cpi0 ')
        if (present(name)) name = 'Intensity_t0'
        if (present(units)) units = 'm2�/2s'
	if (present(info)) info = 'Cold pool intensity at beginning of time step' 
      case('cpi  ')
        if (present(name)) name = 'Intensity'
        if (present(units)) units = 'm2/s2'
	if (present(info)) info = 'Cold pool intensity' 
      case('trig ')
        if (present(name)) name = 'Triggering'
        if (present(units)) units = '-'
	if (present(info)) info = 'Triggering index (0 or 1)' 
      case('tint ')
        if (present(name)) name = 'dLWP_dt'
        if (present(units)) units = '1/s'
	if (present(info)) info = 'Triggering intensity (dLWP_dt)'
      case('cpb  ')
        if (present(name)) name = 'Out_intensity'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Intensity outside cold pools'
      case('ctop ')
        if (present(name)) name = 'CT_height'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Cloud top height'
      case('bltop')
        if (present(name)) name = 'BL_height'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Boundary layer top height'
      case('cbas ')
        if (present(name)) name = 'CB_height'
        if (present(units)) units = 'm'
	if (present(info)) info = 'Cloud base height'
      case('isca ')
        if (present(name)) name = 'SCA_integral'
        if (present(units)) units = 'kg/m2'
	if (present(info)) info = 'Vertically integrated scalar 1'
      case('shf')
        if (present(name)) name = 'SHF'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Sensible Heat Fluxes'
      case('lhf')
        if (present(name)) name = 'LHF'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Latent Heat Fluxes'
      case('sst')
        if (present(name)) name = 'SST'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Surface (potential) temperature'
      case('ssq')
        if (present(name)) name = 'SSQ'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Surface vapor mixing ratio'
      case('rrate')
        if (present(name)) name = 'R_rate'
        if (present(units)) units = 'mm/d'
	if (present(info)) info = 'Surface rain rate'
      case('cumul')
        if (present(name)) name = 'R_cumul'
        if (present(units)) units = 'mm'
	if (present(info)) info = 'Accumulated surface precipitation'
      case('slw  ')
        if (present(name)) name = 'Surface_LW'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Surface long wave net radiative flux'
      case('ssw  ')
        if (present(name)) name = 'Surface_SW'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Surface short wave net radiative flux'
      case('olr  ')
        if (present(name)) name = 'OLR'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Outgoing Long-wave Radiation'
      case('toa  ')
        if (present(name)) name = 'TOA'
        if (present(units)) units = 'W/m2'
	if (present(info)) info = 'Net flux at top of the atmosphere'
      case('osr  ')
        if (present(name)) name = 'OSR'
        if (present(units)) units = 'W/m2'
  if (present(info)) info = 'Outgoing Short-wave Radiation'
      case('isr  ')
        if (present(name)) name = 'ISR'
        if (present(units)) units = 'W/m2'
  if (present(info)) info = 'Incoming Short-wave Radiation'
      case('olrcs')
        if (present(name)) name = 'OLRCS'
        if (present(units)) units = 'W/m2'
  if (present(info)) info = 'Outgoing Long-wave Radiation Clear Sky'
      case('osrcs')
        if (present(name)) name = 'OSRCS'
        if (present(units)) units = 'W/m2'
  if (present(info)) info = 'Outgoing Short-wave Radiation Clear Sky'
      case('mf  ')
        if (present(name)) name = 'MFL'
        if (present(units)) units = 'kg/s'
	if (present(info)) info = 'Total mass flux'
      case('cc  ')
        if (present(name)) name = 'CC'
        if (present(units)) units = '-'
	if (present(info)) info = 'Cloud cover'
      case('mfm ')
        if (present(name)) name = 'MFLup'
        if (present(units)) units = 'kg/s'
	if (present(info)) info = 'Convective mass flux'
      case('ccm ')
        if (present(name)) name = 'CCup'
        if (present(units)) units = '-'
	if (present(info)) info = 'Convective cloud cover'
      case('dmf  ')
        if (present(name)) name = 'dMFL_dz'
        if (present(units)) units = 'kg/m.s'
	if (present(info)) info = 'Vertical derivative of integrated convective mass flux'
      case('cov  ')
        if (present(name)) name = 'Cover'
        if (present(units)) units = '-'
	if (present(info)) info = 'Cloud cover'
      case('dad  ')
        if (present(name)) name = 'dAdt'
        if (present(units)) units = 'kg/m/s'
	if (present(info)) info = 'Cloud cover time rate of change'
      case('mu   ')
        if (present(name)) name = 'MFL_ent'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Total entrainment rate from mass flux'
      case('eq   ')
        if (present(name)) name = 'Qt_ent'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Total entrainment rate from qt'
      case('et   ')
        if (present(name)) name = 'Pt_ent'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Total entrainment rate from pt'
      case('ew   ')
        if (present(name)) name = 'W_ent'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Total entrainment rate from w'
      case('eqs  ')
        if (present(name)) name = 'Qt_ent_s'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Simplified entrainment rate from qt'
      case('ets  ')
        if (present(name)) name = 'Pt_ent_s'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Simplified entrainment rate from pt'
      case('ews  ')
        if (present(name)) name = 'W_ent_s'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Simplified entrainment rate from w'
      case('alq  ')
        if (present(name)) name = 'Alpha_qt'
        if (present(units)) units = '-'
	if (present(info)) info = 'Ratio of total water flux over mass flux times qt'
      case('alt  ')
        if (present(name)) name = 'Alpha_pt'
        if (present(units)) units = '-'
	if (present(info)) info = 'Ratio of potential temperature flux over mass flux times pt'
      case('alw  ')
        if (present(name)) name = 'Alpha_w'
        if (present(units)) units = '-'
	if (present(info)) info = 'Ratio of momentum flux over mass flux times w'
      case('alb  ')
        if (present(name)) name = 'Alpha_b'
        if (present(units)) units = '-'
	if (present(info)) info = 'Ratio of buoyancy flux over mass flux times buoyancy'
      case('alp  ')
        if (present(name)) name = 'Alpha_p'
        if (present(units)) units = '-'
	if (present(info)) info = 'Ratio of pressure gradient flux over mass flux times pressure gradient'
      case('mri  ')
        if (present(name)) name = 'MRi'
        if (present(units)) units = '-'
	if (present(info)) info = 'Mass flux weighted Richardson number'
      case('dl1  ')
        if (present(name)) name = 'Delta1'
        if (present(units)) units = '-'
	if (present(info)) info = 'Delta 1 coefficient'
      case('dl2  ')
        if (present(name)) name = 'Delta2'
        if (present(units)) units = '-'
	if (present(info)) info = 'Delta 2 coefficient'
      case('bet  ')
        if (present(name)) name = 'Alpha'
        if (present(units)) units = '-'
	if (present(info)) info = 'Alpha coefficient'
      case('gam  ')
        if (present(name)) name = 'Beta'
        if (present(units)) units = '-'
	if (present(info)) info = 'Beta coefficient'
      case('ene  ')
        if (present(name)) name = 'ENT_all'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Entrainment rate net'
      case('qte  ')
        if (present(name)) name = 'Qt_env'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Total water mixing ratio averaged over the environment'
      case('tee  ')
        if (present(name)) name = 'PTV_env'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Virtual potential temperature averaged over the environment'
      case('we  ')
        if (present(name)) name = 'W_env'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Vertical velocity averaged over the environment'
      case('enb  ')
        if (present(name)) name = 'EPS_boundary'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Entrainment rate at the plume boundary'
      case('qtb  ')
        if (present(name)) name = 'Qt_boundary'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Total water mixing ratio along the plume boundary'
      case('tb  ')
        if (present(name)) name = 'PT_boundary'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Potential temperature averaged over plumes boundary'
      case('teb  ')
        if (present(name)) name = 'PTV_boundary'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Potential temperature along the plume boundary'
      case('wb  ')
        if (present(name)) name = 'W_boundary'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Vertical velocity along the plume boundary'
      case('bub  ')
        if (present(name)) name = 'Buoy_boundary'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Buoyancy averaged over plumes boundary'
      case('beb  ')
        if (present(name)) name = 'Beff_boundary'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Effective buoyancy averaged over plumes boundary'
      case('een  ')
        if (present(name)) name = 'EPS_entrained'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Entrainment rate entrained'
      case('qen  ')
        if (present(name)) name = 'Qt_entrained'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Total water mixing ratio entrained'
      case('ten  ')
        if (present(name)) name = 'PT_entrained'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Potential temperature entrained'
      case('ven  ')
        if (present(name)) name = 'PTV_entrained'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Virtual potential temperature entrained'
      case('buen ')
        if (present(name)) name = 'Buoy_entrained'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Buoyancy entrained'
      case('been ')
        if (present(name)) name = 'Beff_entrained'
        if (present(units)) units = 'kg/m2/s2'
	if (present(info)) info = 'Effective buoyancy entrained'
      case('wen  ')
        if (present(name)) name = 'W_entrained'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Vertical velocity entrained'
      case('ede  ')
        if (present(name)) name = 'DELT_detrained'
        if (present(units)) units = '1/m'
	if (present(info)) info = 'Detrainment rate detrained'
      case('qde  ')
        if (present(name)) name = 'Qt_detrained'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Total water mixing ratio detrained'
      case('tde  ')
        if (present(name)) name = 'PT_detrained'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Potential temperature detrained'
      case('vde  ')
        if (present(name)) name = 'PTV_detrained'
        if (present(units)) units = 'kg/kg'
	if (present(info)) info = 'Virtual potential temperature detrained'
      case('bude ')
        if (present(name)) name = 'Buoy_detrained'
        if (present(units)) units = 'K'
	if (present(info)) info = 'Buoyancy detrained'
      case('bede ')
        if (present(name)) name = 'Beff_detrained'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Effective buoyancy detrained'
      case('wde  ')
        if (present(name)) name = 'W_detrained'
        if (present(units)) units = 'm/s'
	if (present(info)) info = 'Vertical velocity detrained'
      case('opta2')
        if (present(name)) name = 'OPTAER'
        if (present(units)) units = '-'
  if (present(info)) info = 'Aerosol optical depth'
      case('optc2')
        if (present(name)) name = 'OPTC'
        if (present(units)) units = '-'
  if (present(info)) info = 'Cloud optical depth'  
      case('optr2')
        if (present(name)) name = 'OPTR'
        if (present(units)) units = '-'
  if (present(info)) info = 'Rain optical depth' 
    case('zinv ')
      if (present(name)) name = 'Zi'
      if (present(units)) units = 'm'
  if (present(info)) info = 'Inversion height'
  
      case default
        write(7,*) 'Unrecognized variable name in netcdf output: ', trim(var) 
    end select
    
    return
    end subroutine
    
  end subroutine  register_var
  
  subroutine define_indices ( i1, i2, j1, j2, is, js, NXO, NYO, NZO, NVTO )
  
  integer :: i1, i2, j1, j2, is, js
  integer :: NXO, NYO, NZO, NVTO
!
#ifndef PARALLEL_OUT
      i1 = 4
      i2 = nx-2
#ifdef MODEL_3D
      j1 = 4
      j2 = ny-2
#else
      j1 = 1
      j2 = 1
#endif
#else
      i1 = it_start
      i2 = it_end
#ifdef MODEL_3D
      j1 = jt_start
      j2 = jt_end
#else
      j1 = 1
      j2 = 1
#endif
#endif
!
  is = it_start-4
  js = max(0,jt_start-4)
!
  NXO = i2-i1+1
  NYO = j2-j1+1 
  NZO = nz
  NVTO = 2
  if (kout > 0) NZO = kout
!
  end subroutine define_indices
  
  subroutine check(status,name)
    integer, intent (in) :: status
    character(len=*), intent (in) :: name

    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      print *, 'AN ERROR WAS DETECTED IN NETCDF INTERFACE, FILE/VARIABLE: ', trim(name)
      print *, 'FORCING MIMICA TO STOP (WITHOUT FINALIZING MPI)'
      stop
    end if
  end subroutine check  
  
end module netcdfdef
