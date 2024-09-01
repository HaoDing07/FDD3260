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
!  INITIAL.F                   
!
!  Purpose:
!      A package containing initialization subroutines.                    
!
!  Author
!      Chien Wang
!      MIT Joint Program on Science and Policy of Global Change
!
! ================================================================

  module initialize
  
  USE shared_all
  USE micro_diagnostics
  USE thermodynamics
!
  IMPLICIT NONE 
!  
  private
!  
  public :: getstatus, prt_conf, case_dependent, init_other
! 
  CONTAINS
!
!      ===============================
      Subroutine getstatus
!      ===============================
!
      integer :: k, h, l
      logical :: ideal, with_t=.false., with_rh=.false.
!      
      real :: pp00, zz, dsurf, ptsurf
      real, dimension(nz) :: cp0
      real, dimension(nz,nhydro) :: hydromtr0
!
      if (verbose > 0) call write_debug('Starting getstatus'   )
!      
! ===============================================
!
      if ( trim(casename) == '1D_ADV' .or. trim(casename) == 'SPIRAL' .or.	    &
           trim(casename) == 'CONE' .or. trim(casename) == 'TG_VORTEX' .or.	    &
           trim(casename) == 'VORTEX' .or. trim(casename) == 'SINE' .or.            &
	   trim(casename) == '2D_ADV' .or. trim(casename) == '1D_DIFF' .or.         &
	   trim(casename) == '2D_DIFF' .or. trim(casename) == 'ADJUST_1D' .or.      & 
	   trim(casename) == '2D_NL' .or. trim(casename) == 'IGW' .or.              &
           trim(casename) == 'KM-Sc' ) then
        ideal = .true.
      else
        ideal = .false.
      endif
!
      if (.not. ideal) then
!
      select case ( trim(casename) )
!
!  Idealized cases
!
      case ('DYCOMS','ISDAC31','ISDAC','ISDAC16','ISDAC25','ISDAC34',   &
            'BUBBLE','BUBBLE3D','BUBBLE_PREC','BUBBLE_PREC3D',          &
	    'MOIST_BUBBLE','CURRENT','COLD_POOL','ADJUST_2D',  	        &
	    'RICO', 'BOMEX', 'GERMANY', 'MIDLAT', 'RCERH',     & 
	    'RCEMIP','SE_ATLANTIC','SE_ATLANTIC2')
!
	call init_soundings
!
!  Interpolate initial data on grid
!
      case default
!    
        call read_soundings 
!
      end select
!
      nl0 = 0.
      hydromtr0 = 0.
      if(trim(casename) == 'SE_ATLANTIC' .or. trim(casename) == 'SE_ATLANTIC2') then 
         pt0 = t0*(pref/p0)**((cp_a-cv_a)/cp_a)            
      else
         t0 = pt0*(p0/pref)**((cp_a-cv_a)/cp_a)
      endif
!
      do k = 1, nz
        call get_cp_cv (qt0(k), ql0(k), qi0(k), cp0(k), cv0(k))
        mse0(k) = cp0(k)*(t0(k) - 273.15) + flv00*qv0(k) + (flv00-fls00)*qi0(k) + g*z0(k)
      enddo
!
!  Initial test cases
!
      else
!
        select case (trim(casename))
!
          case ('1D_ADV')
#include "1d_adv.h"
!
          case ('SINE')
#include "sine.h"
!
          case ('2D_NL')
#include "2d_nl.h"
!
          case ('2D_ADV')
#include "2d_adv.h"
!
          case ('SPIRAL')
#include "spiral.h"
!
          case ('CONE')
#include "cone.h"
!
          case ('TG_VORTEX')
#include "tgvortex.h"
!
          case ('VORTEX')
#include "vortex.h"
!
          case ('1D_DIFF')
#include "1d_diff.h"
!
          case ('2D_DIFF')
#include "2d_diff.h"
!
          case ('ADJUST_1D')
#include "adjust1d.h"
!
          case ('IGW')
#include "igw.h"
!
          case ('KM-Sc')
#include "kmsc.h"
!
       end select
!
!  Complementary initialisations
!
        do k = 1, nz
          pp00 = (p0(k)/pref)**((cp_a-cv_a)/cp_a)
	  if (with_t) then
	    pt0(k) = t0(k)/pp00
	  else
	    t0(k) = pt0(k)*pp00
	  endif
	  qv0(k) = qt0(k)
	  ql0(k) = 0.
	  qi0(k) = 0.
	  call get_cp_cv (qt0(k), ql0(k), qi0(k), cp0(k), cv0(k))
          mse0(k) = cp0(k)*(t0(k) - 273.15) + flv00*qv0(k) + (flv00-fls00)*qi0(k) + g*z0(k)
	  den0(k) = p0(k) / ((cp0(k)-cv0(k))*t0(k))
          if (k > 1) avden0(k) = (fdz0(k)*den0(k) + fdz0(k-1)*den0(k-1)) / (fdz0(k)+fdz0(k-1))
        enddo
	avden0(1) = den0(1)
!
      endif
!
!  Precipitation correction
!
      if ( lmicro > 0 ) then
        do h = 1, nhydro
          pxx(1:nz,h) = sqrt(den0(1)/den0(1:nz)) 
        enddo
      endif
!
      if (verbose > 0) call write_debug('Terminating getstatus')
!
      return
      end subroutine getstatus
!
!      ===============================
      Subroutine read_soundings
!      ===============================
!
      integer :: k, l, nl, stat, iscal
      character(len=1) :: car1
      character(len=9) :: line(9)
      character(len=100) :: filename
      logical :: header=.false., rescale=.true., with_rh=.false., with_t=.false., with_ke=.false., with_ce=.false., with_ql=.false., with_w=.false., with_o3=.false.
!      
      real  :: xx, zz
      real, allocatable, dimension(:) :: r_h, r_u, r_v, r_w, r_q, r_ql, r_pt, r_o3, r_so3
!
      if (verbose > 0) call write_debug('Starting read_soundings')
! 
! ===============================================
!
    iscal = min(1,nscal)
!
!  Reading the initial profiles in UNIT 10 (assign.h)
!
      if (file_init(1:len_trim(file_init)) /= 'none') then
        filename = incdir(1:len_trim(incdir)) //'/'// file_init(1:len_trim(file_init))
        open( 10,file=filename(1:len_trim(filename)), form='formatted', status='old' )
      else
        if (mypid == 0) write(7,*) 'ERROR: input file name incorrect, stopping'
	call stop_mimica (0) 
      endif
!
      read(10,*) car1
      if (car1 == '#') header=.true.
!
      stat = 0
      nl = 1
      if (header) nl = 0
      do while (stat >= 0)
        read(10,*,iostat=stat) xx
        nl = nl+1
      enddo
      rewind(10)
      nl=nl-1
!
      allocate(r_h(1:nl), r_u(1:nl), r_v(1:nl), r_w(1:nl), r_pt(1:nl), r_q(1:nl), r_ql(1:nl), r_o3(1:nl), r_so3(1:nl))
!
      if (header) then
        with_t=.true.
        read(10,*) line
	do l = 1, 9
	  if (trim(line(l)) == 'RH') with_rh=.true.
	  if (trim(line(l)) == 'PT') with_t=.false.
	  if (trim(line(l)) == 'QT(*)') rescale=.false.
	  if (trim(line(l)) == 'T(C)') with_ce=.true.
	  if (trim(line(l)) == 'T(K)') with_ke=.true.
	  if (trim(line(l)) == 'QL') with_ql=.true.
	  if (trim(line(l)) == 'W') with_w=.true.
	  if (trim(line(l)) == 'O3') with_o3=.true.
	enddo
      endif
!
      rewind(10)
      if (header) read(10,*) car1
!
      if(with_ql) then
        do l = 1, nl
          read(10,*) r_h(l),r_pt(l),r_q(l),r_ql(l),r_u(l),r_v(l)
        enddo
      else if(with_w) then 
        do l = 1, nl
          read(10,*) r_h(l),r_pt(l),r_q(l),r_u(l),r_v(l),r_w(l),r_o3(l) 
        enddo     
      else if(with_o3) then
        do l = 1, nl
          read(10,*) r_h(l),r_pt(l),r_q(l),r_u(l),r_v(l),r_o3(l)
        enddo
      else
        do l = 1, nl
          read(10,*) r_h(l),r_pt(l),r_q(l),r_u(l),r_v(l)
        enddo
      endif

      if (with_rh) then
        r_q = r_q / 100.		! % --> no unit
      else if (rescale) then
        r_q = r_q / 1000.		! g/kg --> kg/kg
      endif
      if (with_ce) then
	r_pt = r_pt+273.15
      endif
      if (with_o3) then
	r_o3 = r_o3*1.e6		! --> mg/m3
      endif
!
!  Interpolate
!
      k0  = 1.0
      ql0 = 0.0
      qi0 = 0.0
      scal0 = 0.0
      s_scal0 = 0.0
      zz = 0.5*dz*fdz0(1)
      do k = 1, nz
	call pwl_interp_1d ( nl, r_h, r_u, 1, zz, u00(k) )
	call pwl_interp_1d ( nl, r_h, r_v, 1, zz, v00(k) )
        if (.not.with_rh) then
	  call pwl_interp_1d ( nl, r_h, r_q, 1, zz, qv0(k) )
        else
	  call pwl_interp_1d ( nl, r_h, r_q, 1, zz, rh0(k) )
        endif
        if (.not.with_ql) then
	  ql0(k) = 0.
        else
	  call pwl_interp_1d ( nl, r_h, r_ql, 1, zz, ql0(k) )
        endif
        if (.not.with_t) then
	  call pwl_interp_1d ( nl, r_h, r_pt, 1, zz, pt0(k) )
        else
	  call pwl_interp_1d ( nl, r_h, r_pt, 1, zz, t0(k) )
        endif
	if (with_w) call pwl_interp_1d ( nl, r_h, r_w, 1, zz, w0(k) )
	if (with_o3) call pwl_interp_1d ( nl, r_h, r_o3, 1, zz, scal0(k,iscal) )
!        
        if (k < nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))
      enddo
!
      deallocate(r_h, r_u, r_v, r_w, r_pt, r_q, r_ql, r_o3, r_so3)
!
!  Initialise pressure and density
!
      call get_pressure_density (with_rh,with_t)
!
!  Complementary initialisations
!
      if (.not.with_w) then
        if (w_up /= 0.) then
          w0 = w_up
        else
          zz = 0.
          do k = 1, nz
            w0(k) = -Ddiv * zz
            if (k < nz) zz = zz + dz*fdz0(k)
          enddo
        endif
      endif
!
      close (10)
!
      if (verbose > 0) call write_debug('Terminating read_soundings')
!      
      return
      end subroutine read_soundings
!       
!  ===============================
      Subroutine init_soundings
!  ===============================
!
!-----------------------------------------------------------!
!   Enables the construction of initial profiles for        !
!   u, v, p, T and qv using prescribed functions.           !
!   Hard-coded but could be modified easily                 !
!-----------------------------------------------------------!
!
      logical :: set_P, with_t, with_rh
      integer :: k
      real    :: zz, zw, thetae0, rh00, rh01, rh02
!
      if (verbose > 0) call write_debug('Starting init_soundings')
!
! ===============================================
!
    set_p = .true.
    with_t = .false.
    with_rh = .false.
!
!  Default values
!
    u00(:) = 0.
    v00(:) = 0.
    w0(:)  = 0.
    p0(:)  = 100000.
    pt0(:) = 300.
    qt0(:) = 0.
    qv0(:) = 0.
    ql0(:) = 0.
    qi0(:) = 0.
    k0(:)  = 1.
    den0(:)= p0(:) / ((cp_a-cv_a)*pt0(:))
    scal0(:,:) = 0.  
    usurf = u00(1)
    vsurf = v00(1)
!
!  Case dependent initialization
!
    select case (casename)
!
      case ('DYCOMS')
#include "dycoms_ii.h"
!
      case ('ISDAC31','ISDAC')
#include "isdac_f31.h"
!
      case ('ISDAC16')
#include "isdac_f16.h"
!
      case ('ISDAC25')
#include "isdac_f25.h"
!
      case ('ISDAC34')
#include "isdac_f34.h"
!
          case ('BUBBLE')
#include "bubble.h"
!
          case ('BUBBLE3D')
#include "bubble3d.h"
!
          case ('BUBBLE_PREC')
#include "bubble_precip.h"
!
          case ('BUBBLE_PREC3D')
#include "bubble_precip3d.h"
!
          case ('MOIST_BUBBLE')
#include "moist_bubble.h"
!
          case ('CURRENT')
#include "current.h"
!
          case ('COLD_POOL')
#include "cold_pool.h"
!
          case ('ADJUST_2D')
#include "adjust2d.h"
!
          case ('RCEMIP')
#include "rcemip.h"
!
          case ('RICO')
#include "rico.h"
!
          case ('BOMEX')
#include "bomex.h"
!
          case ('GERMANY')
#include "icon.h"
!
          case ('MIDLAT')
#include "midlat.h"
!
          case ('RCERH')
!#include "rcerh55.h"
!#include "rcerh65.h"
!#include "rcerh75.h"
#include "rcerh85.h"
!
          case ('SE_ATLANTIC')
#include "se_atlantic.h"
          with_t = .true.
!
          call read_forcing_SE_Atlantic
!
          call forcing_SE_Atlantic_time1
!
!          case ('SE_ATLANTIC2')
!#include "se_atlantic.h"

!call forcing_SE_Atlantic_time1
      
      end select
!
!  Initialise pressure and density
!
      if(set_p) call get_pressure_density (with_rh,with_t)
!
      if (verbose > 0) call write_debug('Terminating init_soundings')
!
      return
      end subroutine init_soundings
!    
!  ===============================
  subroutine get_pressure_density (with_rh,with_t)
!  ===============================
!
  integer :: k, l, h
  logical :: with_t, with_rh
  real  :: dsurf, cpm, cvm
!
!  Initialize pressure and gas density: first level
!
      p0(1) = psurf
      do l = 1, 10
        if (with_t) then
          pt0(1) = t0(1) / (p0(1)/pref)**((cp_a-cv_a)/cp_a)
        else
          t0(1) = pt0(1) * (p0(1)/pref)**((cp_a-cv_a)/cp_a)
        endif
!
        if (with_rh) then
	  qv0(1) = rh0(1)*cal_qsw(t0(1),p0(1))
	else
	  rh0(1) = qv0(1)/cal_qsw(t0(1),p0(1))
	endif
	qt0(1) = qv0(1) + ql0(1)
!
        call get_cp_cv ( qt0(1), 0., 0., cpm, cvm )
!
        if (with_t) then
          den0(1) = p0(1) / ( (cpm-cvm)*t0(1) )
          dsurf = psurf / ( (cpm-cvm)*t0(1) )
        else
          den0(1) = p0(1) / ((cpm-cvm)*pt0(1)*(p0(1)/pref)**((cp_a-cv_a)/cp_a))
          dsurf = psurf / ((cpm-cvm)*pt0(1)*(psurf/pref)**((cp_a-cv_a)/cp_a))
        endif
!
        p0(1) = psurf - g*0.5*(dsurf+den0(1))*0.5*dz*fdz0(1)
      enddo 
!
!  Initialize pressure and gas density: other level
!
      do k = 2, nz
        p0(k) = p0(k-1) - g*den0(k-1)/dz1(k-1)
        do l = 1, 10
          if (with_t) then
            pt0(k) = t0(k) / (p0(k)/pref)**((cp_a-cv_a)/cp_a)
          else
            t0(k) = pt0(k) * (p0(k)/pref)**((cp_a-cv_a)/cp_a)
          endif
!
          if (with_rh) then
    	    qv0(k) = rh0(k)*cal_qsw(t0(k),p0(k))
	  else
	    rh0(k) = qv0(k)/cal_qsw(t0(k),p0(k))
    	  endif
          qt0(k) = qv0(k) + ql0(k)
!
          call get_cp_cv ( qt0(k), 0., 0., cpm, cvm )
!
          if (with_t) then
            den0(k) = p0(k) / ( (cpm-cvm)*t0(k) )
          else
            den0(k) = p0(k) / ((cpm-cvm)*pt0(k)*(p0(k)/pref)**((cp_a-cv_a)/cp_a))
          endif
!
          p0(k) = p0(k-1) - g*0.5*(den0(k-1)+den0(k))/dz1(k-1)
        enddo 
      enddo 
!
!  Mid-point density
!
      do k = 2, nz
        avden0(k) = (fdz0(k)*den0(k) + fdz0(k-1)*den0(k-1)) / (fdz0(k)+fdz0(k-1))
      enddo 
      avden0(1) = den0(1)
!
!  Reference stability parameter
!
      do k = 2, nz-1
        stab0(k) = (den0(k+1)-den0(k-1)) / (den0(k)*dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
	g0(k) = (p0(k+1)-p0(k-1)) / (p0(k)*dz*(0.5*fdz0(k+1)+fdz0(k)+0.5*fdz0(k-1)))
	g0(k) = g0(k) / stab0(k)
      enddo
      stab0(1) = (den0(2)-den0(1)) / (den0(1)*0.5*dz*(fdz0(2)+fdz0(1)))
      stab0(nz) = (den0(nz)-den0(nz-1)) / (den0(nz)*0.5*dz*(fdz0(nz)+fdz0(nz-1)))
      g0(1) = g0(2)
      g0(nz) = g0(nz-1)
      n20 = -g*stab0
!
!  Pressure correction
!
      p1 = 0.0
!
    end subroutine get_pressure_density
!    
!  ===============================
  subroutine case_dependent
!  ===============================

  integer :: i, k, l, iscal
  integer :: ii, kk
  real :: xx, zz, r0, rr, rx, ru, rz, T, Lx, Hz, pt
!
! =============================================== 
! 
!  Case specific init.
!
    iscal = min(1,nscal)
! 
    select case (casename)
!
      case ("1D_ADV")
        ii = (ip_end-ip_start-mod(ip_end-ip_start,2))/2
        do i = ip_start, ii
   	      state2%scal(i,:,:,iscal)  = 1.
	    enddo
!
      case ("1D_ADV_VERT")
        kk = nz/2+1
        do k = 1, kk
   	      state2%scal(i,:,:,iscal)  = 1.
	    enddo
!
      case ("SINE")
	xx = 0.0
        do i = it_start, it_end
   	  state2%scal(i,:,:,iscal) = 1.0 + sin(2.0*pi*xx)
   	  xx = xx + dx
	enddo
!
      case ("1D_DIFF")
		xx = 0.0
        do i = it_start, it_end
   	      if (xx >= 4.0 .and. xx <= 6.0) state2%scal(i,:,:,iscal) = 1.0
   	      xx = xx + dx
	    enddo
!
      case ("2D_ADV")
        xx = (it_start-4)*dx
        do i = it_start, it_end
          zz = 0.5*dz
          do k = 1, nz
            rr = min(1.,4.*sqrt((xx - .5)**2. + (zz - .5)**2.))
   	        state2%scal(i,:,k,iscal)  = 0.5*(1. + cos(pi*rr))
            zz = zz + dz
	      enddo
          xx = xx + dx
	    enddo
!
      case ("2D_NL")
        xx = (it_start-4)*dx
        do i = it_start, it_end
          zz = 0.5*dz
          do k = 1, nz
	    if ( xx >= 0.5 .and. xx <= 1. .and. zz >= 0.5 .and. zz <= 1. ) then
              wind2%u(i,:,k) = 2.
              wind2%w(i,:,k) = 2.
            endif 
            zz = zz + dz
	  enddo
          xx = xx + dx
	enddo
!
      case ("2D_DIFF")
        xx = (it_start-4)*dx
        do i = it_start, it_end
          zz = 0.5*dz
          do k = 1, nz
            rr = sqrt((xx - 2.5)**2. + (zz - 2.5)**2.)
   	        if (rr <= 1.) state2%scal(i,:,k,iscal) = 1.
            zz = zz + dz
	      enddo
          xx = xx + dx
	    enddo
!
      case ("SPIRAL")
        xx = 0.0
        state2%scal(:,:,:,iscal) = 0.0
        do i = it_start, it_end
          zz = dz*fdz0(1)
          do k = 4, nz-2
    	    wind2%u(i,:,k) = sin(pi*(xx-dx/2.))*sin(pi*(xx-dx/2.))*sin(2.*pi*zz)
    	    wind2%w(i,:,k) = -sin(pi*(zz-dz/2.))*sin(pi*(zz-dz/2.))*sin(2.*pi*xx)     
           
            rr = min(1.,4.*sqrt((xx - .25)**2. + (zz - .25)**2.))
            state2%scal(i,:,k,iscal) = 0.5*(1. + cos(pi*rr))
            
            if (k < nz) zz = zz + 0.5*dz*(fdz0(k+1)+fdz0(k))
          enddo
          xx = xx + dx
        enddo
!
      case ("CONE")
        xx = (it_start-5)*dx
        state2%scal(:,:,:,iscal) = 0.0
        do i = it_start, it_end
          zz = dz*fdz0(1)
          do k = 1, nz
            wind2%u(i,:,k) = -(zz - 0.5)
            wind2%w(i,:,k) = (xx - 0.5)
            
            rr = min(sqrt((xx - 0.5)**2. + (zz - 0.25)**2.),0.2)/0.2
            state2%scal(i,:,k,iscal) = 0.
            if (rr < 1.) state2%scal(i,:,k,iscal) = 1. - rr
            
            if (k < nz) zz = zz + 0.5*dz*(fdz0(k+1)+fdz0(k))
          enddo
          xx = xx + dx
        enddo
!
      case ("TG_VORTEX")
        xx = 0.0
        r0 = 0.05
        do i = it_start, it_end
          zz = 0.5*dz
          do k = 4, nz-2
            wind2%u(i,:,k) = 0.0
            wind2%w(i,:,k) = 0.0 
            
            wind2%u(i,:,k) = wind2%u(i,:,k) + cos(pi*(xx-dx/2.))*sin(pi*zz)
            wind2%w(i,:,k) = wind2%w(i,:,k) - sin(pi*xx)*cos(pi*(zz-dz/2.))
            do l = 1, 10
              pressure2%p(i,:,k) = -0.25*pressure2%dens(i,:,k)*(cos(2.*pi*xx) + cos(2.*pi*zz))
              pressure2%dens(i,:,k) = (pressure2%p(i,:,k)+p0(k)) / ((cp_a-cv_a)*Pt0(1))
            enddo            
            zz = zz + dz
          enddo
          xx = xx + dx
        enddo
!
      case ("VORTEX")
        xx = (it_start-4)*dx
        do i = it_start, it_end
          zz = 0.5*dz
          do k = 1, nz
            wind2%u(i,:,k) = 0.0 
            wind2%w(i,:,k) = 0.0 
            
            rr = 0.05
            rx = sqrt((xx - 0.5)**2. + (zz - 0.5)**2.) / rr
            ru = sqrt((xx-dx/2. - 0.5)**2. + (zz - 0.5)**2.) / rr
            rz = sqrt((xx - 0.5)**2. + (zz-dz/2. - 0.5)**2.) / rr
            if (k <= nz) then
              wind2%u(i,:,k) = u00(k) - 0.1*u00(k)*exp(-ru*ru/2.) * (zz - 0.5)/rr
              wind2%w(i,:,k) = 0.1*u00(k)*exp(-rz*rz/2.) * (xx - 0.5)/rr
              T = Pt0(1) - 0.5*(0.1*u00(k)*exp(-rx*rx/2.))**2. / cp_a
              pressure2%dens(i,:,k) = den0(k)*(T/Pt0(1))**(cp_a/(cp_a-cv_a))
              pressure2%p(i,:,k) = pressure2%dens(i,:,k)*(cp_a-cv_a)*T - p0(k)
              if (k <= nz) state2%es(i,:,k) = T*((pressure2%p(i,:,k)+p0(k))/1.E+05)**(-(cp_a-cv_a)/cp_a)
            endif
            
            zz = zz + dz
          enddo
          xx = xx + dx
        enddo
!
      case ("IGW")
        xx = 0.0
        do i = it_start, it_end
          zz = 0.5*dz
          do k = 1, nz
            wind2%u(i,:,k) = 20.
            wind2%w(i,:,k) = 0.
            T = 1.e-2*sin(pi*zz/10000.) / (1. + (xx - 100000.)**2./(5000.)**2.) 
            state2%es(i,:,k) = pt0(k) + T
            pressure2%dens(i,:,k) = den0(k)*(1. - T/pt0(k))
            zz = zz + dz
          enddo
          xx = xx + dx
        enddo
!
      case ("KM-Sc")
	Lx = 0.5*real(it_end-it_start+1)*dx
	Hz = real(nz)*dz
        xx = 0.
        do i = it_start, it_end
          zz = 0.
          do k = 1, nz
    	    wind2%u(i,:,k) = 1.7*Lx/Hz*cos(pi*xx/Lx)*cos(pi*(zz+0.5*dz)/Hz)
    	    wind2%w(i,:,k) = 1.7*sin(pi*(xx+0.5*dx)/Lx)*sin(pi*zz/Hz)
            zz = zz + dz
          enddo
          xx = xx + dx
        enddo
!
      case default
!
!  Perturbations and bubbles
!
      if (casename == 'BUBBLE' .or. casename == 'BUBBLE3D' .or. casename == 'CURRENT' .or. casename == 'COLD_POOL') then
        call init_bubble (state2%es, pressure2%dens)
      else if (casename == 'MOIST_BUBBLE') then
        call init_moist_bubble (state2%es, pressure2%dens)
      else if (casename == 'BUBBLE_PREC' .or. casename == 'BUBBLE_PREC3D') then
        call init_prec_bubble (state2%es, pressure2%dens, state2%qt)
      else
        call init_random (pressure2%dens, state2%es, turbu%ksgs)
      endif
!
!  Perturb vertical wind
!
      if (dw > 0.) call init_random_wind ( wind2%w )
!  
    end select
!
contains

!      ====================================
      subroutine init_bubble (ptd, den)
!      ====================================
!
      INTEGER, DIMENSION (8) :: T
      REAL :: RANF, tmp
      INTEGER :: SEED1, SEED2   
!
      integer  :: k, i, j
      real     :: zz, xx, yy, rr, pp0, cpm, cvm
 !
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: den
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ptd
!            
! -----------------------------------------------------------------------
!
      if (casename == 'SQUALL_LINE') then
        xrbubble = 4000.
        zrbubble = 1500.
        xbubble = 500000.
        zbubble = 1500.
        dtbubble = 3.
       else if (casename == 'DEEP_CONVECTION') then
        xrbubble = 10000.
	yrbubble = 10000.
        zrbubble = 1500.
	xbubble = 65000.
	ybubble = 65000.
	zbubble = 1500.
	dtbubble = 2.
      endif
!
!  Initial Random number generator
!
      CALL DATE_AND_TIME(VALUES = T)
      SEED1 = 43*mypid+7*T(1)+70*(T(2)+12*(T(3)+31*(T(5)+23*(T(6)+59*T(7)))))
      SEED2 = 78*mypid+7*T(1)+89*(T(2)+3*(14*T(3)+(T(5)+9*(5*T(6)+59*T(7)))))
      IF (MOD(SEED1,2).EQ.0) SEED1 = SEED1-1
      IF (MOD(SEED2,2).EQ.0) SEED2 = SEED2-1
!
!  Initialize dry bubble
!
      zz = 0.5*dz*fdz0(1)
      do k = 1, nz
        call get_cp_cv (qt0(k), ql0(k), qi0(k), cpm, cvm)
!
	yy = max(jt_start-4,0)*dy
        do j = jt_start, jt_end
	  xx = (it_start-4)*dx
	  do i = it_start, it_end
	    pp0 = (p0(k)/pref)**((cp_a-cv_a)/cp_a)
 	    rr = ((xx - xbubble)/xrbubble)**2. + ((zz - zbubble)/zrbubble)**2.
#ifdef MODEL_3D
            rr = rr + ((yy - ybubble)/yrbubble)**2.
#endif
!
	    rr = sqrt(rr)
	    tmp = max(min(sqrt(-2.*log(max(RANF(SEED1),1.e-10)))*sin(2.*pi*RANF(SEED2))/sqrt(20.),1.),-1.)
	    if (rr < 1.) ptd(i,j,k) = ptd(i,j,k) + 0.5*(dtbubble+dpt*tmp)*(1. + cos(pi*rr))
!
	    den(i,j,k) = p0(k) / ((cpm-cvm)*ptd(i,j,k)*pp0)
	    xx = xx + dx
	  enddo
	yy = yy + dy
	enddo
!
        if (k/=nz) zz = zz + dz*0.5*(fdz0(k)+fdz0(k+1))  
      enddo
!     
      return
     end

!      ====================================
      subroutine init_prec_bubble (ptd, den, qt)
!      ====================================
!
      integer  :: k, i, j
      real     :: zz, xx, yy, rr
      real     :: pp0, cpm, cvm, qsw
 !
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: den, qt,ptd
!            
! -----------------------------------------------------------------------
!
!  Initialize dry bubble
!
      zz = 0.5*dz*fdz0(1)
      do k = 1, nz
      yy = max(jt_start-4,0)*dy
      do j = jt_start, jt_end
	  xx = (it_start-4)*dx
	  do i = it_start, it_end
	  pp0 = (p0(k)/pref)**((cp_a-cv_a)/cp_a) 
	  qsw = cal_qsw(pt0(k)*pp0, p0(k))
!
 	  rr = ((xx - xbubble))**2. + ((zz - zbubble))**2.
#ifdef MODEL_3D
	  rr = rr + ((yy - ybubble))**2.
#endif
!
	  rr = sqrt(rr)
	  if (rr <= 200.) then
	    qt(i,j,k) = qsw
	  else if (rr > 200. .and. rr <= 300.) then
	    qt(i,j,k) = (0.2 + 0.8*cos(0.005*pi*(rr - 200.))*cos(0.005*pi*(rr - 200.)))*qsw
	  else
	    qt(i,j,k) = 0.2*qsw
	  endif
!
          call get_cp_cv (qt(i,j,k), 0., 0., cpm, cvm)
!
	  den(i,j,k) = p0(k) / ((cpm-cvm)*ptd(i,j,k)*pp0)
	  xx = xx + dx
	  enddo
	  yy = yy + dy
	  enddo
!
        if (k/=nz) zz = zz + dz*0.5*(fdz0(k)+fdz0(k+1))   
      enddo
!     
      return
      end

!      ====================================
      subroutine init_moist_bubble ( ptd, den )
!      ====================================
!
      integer  :: k, i, j, l
      real     :: zz, xx, rr, pp0, qv, ev, tt(nz), cpm, cvm, thetae0, ttold, tol
 !
 	  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: den
	  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ptd
!
#include "moist_bubble.h"
!            
! -----------------------------------------------------------------------
!
      j = ny
      zz = 0.5*dz*fdz0(1)
      qv0 = 0.0
      do k = 1, nz
	xx = (it_start-4)*dx
!
!  Internal sub-iterations 1: find hydrostatic state
!
	ev = 0.0
        qv = 0.0
	tt(k) = pt0(k)
	ttold=tt(k)
	tol = 1.
	do while (tol > 1.e-12)
	  if (qt0(k) > 0.) ev = cal_esw(tt(k))
	  if (qt0(k) > 0.) qv0(k) = cal_qsw(tt(k), p0(k))
          ql0(k) = max(qt0(k) - qv0(k),0.)
          call get_cp_cv (qt0(k), ql0(k), qi0(k), cpm, cvm)
!
 	  tt(k) = thetae0*((p0(k)-ev)/pref)**((cp_a-cv_a)/(cp_a+cp_l*qt0(k)))*exp(-cal_flv(tt(k))*qv0(k)/(tt(k)*(cp_a+cp_l*qt0(k))))
	  pp0 = (p0(k)/pref)**((cp_a-cv_a)/cp_a)
	  pt0(k) = tt(k)/pp0
	  den0(k) = p0(k) / ((cpm-cvm)*tt(k))
	  if (k > 1) p0(k) = p0(k-1) - g*0.5*(den0(k-1)+den0(k))/dz1(k-1)
!
	  tol = abs(tt(k)-ttold)/tt(k)
	  ttold=tt(k)
	enddo
	ptd(:,:,k) = pt0(k)
	hydromtr2(drop)%q(:,:,k) = ql0(k)
!
	do i = it_start, it_end
          rr = sqrt(((xx - xbubble)/xrbubble)**2. + ((zz - zbubble)/zrbubble)**2.)
!
!  Internal sub-iterations 2: find perturbed state
!
	  do l = 1, 15
	    if (qt0(k) > 0.) ev = cal_esw(tt(k))
	    if (qt0(k) > 0.) qv = cal_qsw(tt(k), p0(k))
	    hydromtr2(drop)%q(i,j,k) = max(qt0(k) - qv,0.)
            call get_cp_cv (qt0(k), hydromtr2(drop)%q(i,j,k), qi0(k), cpm, cvm)
!  
            if (rr <= 1.) ptd(i,j,k) = pt0(k)*(1. + (cp_v-cv_v)/(cp_a-cv_a)*qv0(k)) / (1. + (cp_v-cv_v)/(cp_a-cv_a)*qv) *     &
  	  		    (1. + dtbubble/300.*cos(0.5*pi*rr)*cos(0.5*pi*rr))
!
            pp0 = (p0(k)/pref)**((cp_a-cv_a)/cp_a)
  	    tt(k) = ptd(i,j,k)*pp0
            den(i,j,k) = p0(k) / ((cpm-cvm)*tt(k))
          enddo
!
	  xx = xx + dx
	enddo
!
        if (k/=nz) zz = zz + dz*0.5*(fdz0(k)+fdz0(k+1))  
      enddo
!     
      return
      end
!
!      ====================================
      subroutine init_random (den, ptd, ksgs)
!      ====================================
!
      integer  :: i,j,k
      
      INTEGER, DIMENSION (8) :: T
      REAL :: X,Y,RANF      
      INTEGER :: SEED      
      
      real  :: fac, tmp, cpm, cvm
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: den
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ptd
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ksgs
!
!----------------------------------------------
!
!  Initial Random number generator
!
      CALL DATE_AND_TIME(VALUES = T)
      SEED = 43*mypid+7*T(1)+70*(T(2)+12*(T(3)+31*(T(5)+23*(T(6)+59*T(7)))))
      IF (MOD(SEED,2).EQ.0) SEED = SEED-1
!
!  Applying random perturbations
!
      do k = 1, kpert
        fac = 1. !real(kpert-k)/real(kpert-1)
        call get_cp_cv (qt0(k), ql0(k), qi0(k), cpm, cvm)
!
        do j = jp_start, jp_end
          do i = ip_start, ip_end
            X = RANF(SEED)
            Y = RANF(SEED)
            tmp = max(min(sqrt(-2.*log(max(X,1.e-10)))*sin(2.*pi*Y)/2.,1.),-1.)
!
            ptd(i,j,k) = fac*tmp*dpt + ptd(i,j,k)
            den(i,j,k) = p0(k)**(cv_a/cp_a) / ((cpm-cvm)*ptd(i,j,k)*pref**(cv_a/cp_a-1.))
!
#ifdef TKE
            ksgs(i,j,k) = fac*tmp*k0(k) + ksgs(i,j,k)
#endif	    
          enddo
        enddo
      enddo
!
      return
      end
!
!      ====================================
      subroutine init_random_wind (w)
!      ====================================

      integer  :: i,j,k
      
      INTEGER, DIMENSION (8) :: T
      REAL :: X,Y,RANF      
      INTEGER :: SEED      
      real  :: fac, tmp
      real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: w
!
!----------------------------------------------
!
!  Initial Random number generator
!
      CALL DATE_AND_TIME(VALUES = T)
      SEED = 13*mypid+7*T(1)+70*(T(2)+12*(T(3)+31*(T(5)+23*(T(6)+59*T(7)))))
      IF (MOD(SEED,2).EQ.0) SEED = SEED-1
!
!  Applying random perturbations
!
      do k = 1, kpert
        fac = 1. !real(kpert-k)/real(kpert-1)
        do j = jp_start, jp_end
          do i = ip_start, ip_end
            X = RANF(SEED)
            Y = RANF(SEED)
            tmp = max(min(sqrt(-2.*log(max(X,1.e-10)))*sin(2.*pi*Y)/2.,1.),-1.)
            w(i,j,k) = fac*tmp*dw
          enddo
        enddo
      enddo
!
      return
      end
! 
  end subroutine case_dependent
!
  subroutine prt_conf

!    =================================================
!    Subroutine for printing out initial state and options
!    =================================================
      
      integer :: k, mmx
      
      character (len = 12) :: yes_no
      
      write(7,1004)version(1:8),casename
1004      format(1x/6x,' MIMICA-V5, Version: ', a8,', Case: ',a20,/)

#ifdef MEANDATA
          yes_no = 'Time Means'
#else
          yes_no = 'Instant Data'
#endif
      write(7,1204)yes_no
1204  format(8x,'Output Fields are Derived based on ',a12/,           &
             8x,'================================================'/)

#ifndef FINE
      write(7,1005)dt,dx,dy,dz,nx,ny,nz
#else
      write(7,1005)dt,dx,dy,dz,dz*sratio,nx,ny,nz
#endif
1005  format(8x,'Basic time step = ',f8.2,' second'/                 &
#ifndef FINE
       '        dx, dy, & dz    = ',3f7.1,' meter'/                  &
#else
       '        dx, dy, & dz (max,min)    = ',4f7.1,' meter'/        &
#endif
       '        Gridpoint number= ',3i7,' in x, y, & z direction'/   &
       '        Lateral boundary: ',a12,' type condition'/           &
       '                          ',a12,' type condition'/)

#if ( defined SPMD )
      write(7,1010)nproc
1010      format(8x,'SPMD run using    ',i6,' processors'/)
#else
      write(7,1110)
1110      format(8x,'Single CPU run'/)
#endif

#if ( defined CORIOLIS )
      write(7,1014)fcor(1,1),fcor(nx,ny)
1014  format(8x,'Coriolis coefficient (1,1) and (nx,ny) = ',2f12.5/)
#endif

#ifdef AERO_ENABLE

       write(7,2006)xn_ccn0*1.e-6,xn_in0*1.e-3
2006   format(8x,'Initial Accuml-mode CN = ',                      &
                 f10.2,' 1/cc; IN = ',                            &
                 f10.2,' 1/L'/)
#else
       write(7,1006)xnc0_d*1.e-6,xnc0_k
1006   format(8x,'In calculating CN nucleation rate Nc = Cs^k'/     &
       '      C = ',f10.2,' 1/cc;  k = ',f10.2/)

       write(7,2006)xn_ccn0*1.e-6,xn_in0*1.e-3
2006   format(8x,'Reference concentrations are '/                &
       '      CCN = ',f10.2,' 1/cc; IN = ',                      &
                 f10.2,' 1/L'/)
#endif

#ifdef RAD_ENABLE
      write(7,604)abs(iradx)*dt/60.0
604      format(8x,'Radiation is calculated at every ',f7.1, ' minute'/)
#else
      write(7,6041)
6041      format(8x,'Run without Radiation '/)
#endif

#ifdef CHEM_ENABLE

#ifdef CHEM_REACT
606      format(8x,'Chemical Reactions are calculated '/)
#else
      write(7,6061)
6061    format(8x,'Transport-only Chemistry '/)
#endif
#else
      write(7,605)
605     format(8x,'Run without Chemistry '/)
#endif


#ifdef AQCHEM_ENABLE
608      format(8x,'Aqueous Reactions are calculated '/)
#else
      write(7,607)
607     format(8x,'Run without Aqueous Chemistry '/)
#endif

#ifdef SOLIDCHEM_ENABLE
609      format(8x,'Heterogeneous Reactions are calculated '/)
#else
      write(7,610)
610     format(8x,'Run without Heterogeneous Chemistry '/)
#endif

#ifdef AERO_ENABLE
      write(7,611)
611  format(8x,'Run with multimode and multimoment aerosol scheme ')
      write(7,612) , aero_flg%tran, aero_flg%act,  aero_flg%act_scv, aero_flg%reg, aero_flg%imp_scv, aero_flg%chem
612  format(10x,'simulate aerosol processes: '/ ,&
                '            transport:             ', l, /&
                '            activation:            ', l,/&
                '            activation scavenging: ', l,/&
                '            regeneration:          ', l,/&
                '            impaction scavenging:  ', l,/&
                '            chemistry:             ',l )
#else
      write(7,613)
613     format(8x,'Run with single mode and moment aerosol scheme '/)
#endif

      write(7,660)
660      format(8x,60('=')/)

      write(7,2022)
2022      format(                               &
               6x,3x,'Z(M)',8x,'PI0',           &
               4x,'PT0(K)',1X,'QT0(G/KG)',      &
               3x,'U0(m/s)',2x,'V0(m/s)',2x,'W0(m/s)',5x,'DEN0(KG/M3)'/)        

       mmx = nx/2
! 
       do k = nz,1,-1                                 
            write(7,2026)k,z0(k),p0(k),pt0(k),1000.*qt0(k),u00(k),v00(k),w0(k),den0(k)
       end do

       flush(7)

2026   format(1x,1x,i3,2x,f8.1,2x,f8.1,2x,f8.2,2x,f8.2,  &             
        2x,f8.2,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,3x,f8.4)


      return
      end subroutine prt_conf
!    
     !======================================================
      subroutine forcing_SE_Atlantic_time1
      !=====================================================
      USE shared_data 
      USE shared_all   
      
      integer :: time_simulation, k, t 
      
      time_simulation=1

      !nz-k+1 is used to reverse the order of the arrays that are
      FORALL (k=1:nz) t0(k)=ta_nud(lev_len-k+1,time_simulation)
      FORALL (k=1:nz) qv0(k)=qt_nud(lev_len-k+1,time_simulation)
      FORALL (k=1:nz) u00(k)=u_nud(lev_len-k+1,time_simulation)
      FORALL (k=1:nz) v00(k)=v_nud(lev_len-k+1,time_simulation) 
      
      !ta_nud (absolute temperature) is converted to potential temperature
      FORALL(k=1:lev_len,t=1:time_len) ta_nud(k,t) = ta_nud(k,t)*((pref/pa_force(k,t))&
                                                      **((cp_a-cv_a)/cp_a))
      
      return
      
      end

      !======================================================
      subroutine read_forcing_SE_Atlantic
      !=======================================================
      USE netcdf
      USE shared_data 
      USE shared_all   
      USE shared_aerosol_new
      USE shared_wind
      USE shared_state
      USE gridno
      USE shared_thermo
      USE mpi
      
      IMPLICIT NONE
      
      integer :: ii, ncid, status, time_simulation, k, im, ierr
      integer :: time_id,lat_id,lon_id,lev_id
      integer :: lat_len,lon_len 
      integer :: tsec_id, Na_accum_id, q_id !,PT_id 
      integer :: ta_id, pa_id, ts_id, ps_id
      integer :: wind_u_id, wind_v_id !, wind_w_id
      
      if (verbose > 2) call write_debug ('Starting read_forcing_SE_Atlantic')
      
      !Opening file with atmospheric conditions for force the model.
      call check(NF90_OPEN(forcing_file,NF90_NOWRITE,ncid))
      
      !getting the dimension ids 
      call check(NF90_INQ_DIMID(ncid,"time",time_id))  
      call check(NF90_INQ_DIMID(ncid,"zh",lev_id))    
      call check(NF90_INQ_DIMID(ncid,"lat",lat_id))
      call check(NF90_INQ_DIMID(ncid,"lon",lon_id))
      
      !getting the dimension lengths
      call check(NF90_INQUIRE_DIMENSION(ncid,time_id,len=time_len))
      call check(NF90_INQUIRE_DIMENSION(ncid,lev_id,len=lev_len))  
      call check(NF90_INQUIRE_DIMENSION(ncid,lat_id,len=lat_len))
      call check(NF90_INQUIRE_DIMENSION(ncid,lon_id,len=lon_len)) 
      
      !- Allocation of variables : 
      ALLOCATE(  Na_force(lev_len,time_len))  
      ALLOCATE(  qt_nud(lev_len,time_len)) 
      !ALLOCATE(  PT(lev_len,time_len)) 
      ALLOCATE(  ta_nud(lev_len,time_len))
      ALLOCATE(  pa_force(lev_len,time_len))  
      ALLOCATE(  u_nud(lev_len,time_len))
      ALLOCATE(  v_nud(lev_len,time_len))
      !ALLOCATE(  w_nud(lev_len,time_len))
      ALLOCATE(  ts_force(time_len))
      ALLOCATE(  ps_force(time_len))
            
      !Getting the variables ids 
      call check(NF90_INQ_VARID(ncid,"na_nud",Na_accum_id))
      call check(NF90_INQ_VARID(ncid,"qt_nud",q_id))
      call check(NF90_INQ_VARID(ncid,"ta_nud",ta_id))  
      call check(NF90_INQ_VARID(ncid,"pa_force",pa_id)) 
      !call check(NF90_INQ_VARID(ncid,"thetal_nud",PT_id))
      call check(NF90_INQ_VARID(ncid,"ua_nud",wind_u_id))
      call check(NF90_INQ_VARID(ncid,"va_nud",wind_v_id))
      !call check(NF90_INQ_VARID(ncid,"wa",wind_w_id))
      call check(NF90_INQ_VARID(ncid,"ts_force",ts_id))
      call check(NF90_INQ_VARID(ncid,"ps_force",ps_id))  
            
      !Reading the data
      call check(NF90_GET_VAR(ncid,Na_accum_id,Na_force))  
      call check(NF90_GET_VAR(ncid,q_id,qt_nud))  
      call check(NF90_GET_VAR(ncid,ta_id,ta_nud))
      call check(NF90_GET_VAR(ncid,pa_id,pa_force))  
      call check(NF90_GET_VAR(ncid,wind_u_id,u_nud))
      call check(NF90_GET_VAR(ncid,wind_v_id,v_nud))
      !call check(NF90_GET_VAR(ncid,wind_w_id,w_nud))
      call check(NF90_GET_VAR(ncid,ts_id,ts_force))
      call check(NF90_GET_VAR(ncid,ps_id,ps_force))  
      
      !closing file
      CALL check(nf90_close(ncid))

      if (verbose > 2) call write_debug ('Terminating reaad_forcing_SE_Atlantic')

      return
      end

      !===============================================

      SUBROUTINE check(istatus)
      !!Checks the return status of a nf90 function. If an error
      !!occured, the according error message is written to the
      !!screen and the program exits
    
      USE netcdf
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: istatus
      IF (istatus /= nf90_noerr) THEN
         WRITE(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
         WRITE(*,*) 'exiting'
         STOP 
      END IF
    END SUBROUTINE check
    
    !===============================================

    subroutine init_other

    IMPLICIT NONE
    real :: dlat
!
!  Coriolis
!
    dlat = abs(ctr_lat)
    fcor = f2omega*sin(pi*dlat/180.)
!
!  Surface fluxes
!
    shf = shf0
    lhf = lhf0
!
    return
    end subroutine 

  end module initialize
