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
! ================================================================

  module columnphys
  
  USE shared_all
  USE sources
  USE surfacemod
  USE gradients
  USE thermodynamics
  USE getcape
!
  IMPLICIT NONE
!
  private
!
  integer, parameter :: adjust = 1  ! =0 moist adiabat, =1 fixed lapse rate, =2 adjustment to rce
!
  public :: pbl_mixing, convection

  CONTAINS
!
! ===============================================
  subroutine convection ( pressurel, statel, hydromtrl )
! ===============================================
!
  type(atm_state), intent(inout) :: statel
  type(atm_pressure) :: pressurel
  type(hydrometeor), dimension(1:nhydro) :: hydromtrl
! 
  logical :: cyc
  integer :: i, j, k, l, nl, stat
  real :: tn, cape, cin, zb, zt, dt, dq, qm, tau=1800.
  real :: rhi, rhc, qsc, qsi
  real, dimension(:), allocatable :: zref, ptref, qtref
  real, dimension(1:nz) :: dth, dqv, exn, tref, rhref, ptadj, qtadj
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: tem, tdew
!
!-----------------------------------------------------------!
!              	        Calculate CAPE                      !
!-----------------------------------------------------------!
!
  call get_temperature ( pressurel, statel, hydromtrl, tem )
!
  call get_tdew ( tem, pressurel, statel, hydromtrl, tdew )
!
!  Various initializations
!
  exn = (p0/pref)**((cp_a-cv_a)/cp_a)
!
  tref(1) = sst - 0.5*0.0065/dz1(1)
  rhi = 0.8
  rhc = 0.5
  do k = 2, nz
    tref(k) = tref(k-1) - 0.5*0.0065*(1./dz1(k-1) + 1./dz1(k))
  enddo
!
  w0 = w_up
!
!  Read external profiles
!
  call read_file ( zref, ptref, qtref )
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      cyc = .false.
      dth = 0.
      dqv = 0.
      ptadj = statel%es(i,j,:)
      qtadj = statel%qt(i,j,:)
      if (adjust == 0) then
        !call get_cape ( nz, p0, tem(i,j,:), tdew(i,j,:), cape, cin, dth=dth, dqv=dqv, zt=zt )
        ptadj = statel%es(i,j,:) + dth
        do k = 1, nz
          if (dth(k) /= 0.) qtadj(k) = max(0.5*cal_qsw(tem(i,j,k)+dth(k)*exn(k),p0(k)),qm)
        enddo
      else if (adjust == 1) then
        do k = 2, nz
          qsc = rhc*cal_qsw(tem(i,j,k),p0(k))
          qsi = rhi*cal_qsi(tem(i,j,k),p0(k))
          qtadj(k) = max(min(qsc,qsi),0.)
          if (qtadj(k) > qtadj(k-1)) qtadj(k) = qtadj(k-1)
          if (rad%dtnet(i,j,k) > 0. .and. z0(k) > 8000.) cyc = .true.
          if (cyc) cycle
          ptadj(k) = tref(k)/exn(k)
          w0(k) = 0.
          zt = z0(k)
        enddo
      else
        !call get_cape( nz, p0, tem(i,j,:), tdew(i,j,:), cape, cin, zt=zt )
        do k = 1, nz
          if (zref(k) < zt .and. k < nl) then
            ptadj(k) = ptref(k)
            qtadj(k) = qtref(k)
          endif
        enddo
      endif
!
      dt = 0.; dq = 0.
      exn = (p0/pref)**((cp_a-cv_a)/cp_a)
      do l = 1, 2
        do k = 2, nz
          if (l == 1) then
            if (z0(k) <= zt) then
              dt = dt + 0.5*((ptadj(k)*exn(k) + ptadj(k+1)*exn(k+1) - tem(i,j,k) - tem(i,j,k+1)) + flv00/cp_a*(qtadj(k) + qtadj(k+1) - statel%qt(i,j,k) - statel%qt(i,j,k+1)))/dz1(k)/zt
            endif
          else
            if (z0(k) <= zt) then
              dth(k) = ptadj(k) + dt/exn(k) - statel%es(i,j,k)
            endif
            dqv(k) = qtadj(k) - statel%qt(i,j,k)
          endif
        enddo
      enddo
!
      states%es(i,j,:) = states%es(i,j,:) + dth(:)/tau
      states%qt(i,j,:) = states%qt(i,j,:) + dqv(:)/tau
      if (out_diagt) diag(4)%pt(i,j,:) = diag(4)%pt(i,j,:) + dth(:)/tau
      if (out_diagq) diag(4)%qt(i,j,:) = diag(4)%qt(i,j,:) + dqv(:)/tau
    enddo
  enddo
!
  if (allocated(zref)) deallocate(zref,ptref,qtref)
!
!-----------------------------------------------------------!
!
return
!
contains
!
  subroutine read_file(zref,ptref,qtref)
        logical :: ex
        integer :: nl, stat
        real, dimension(:), allocatable :: zref, ptref, qtref

        inquire(FILE="./rce_reference", EXIST=ex)
        if (ex) then
          open( 10,file='./rce_reference', form='formatted', status='old' )
          read(10,*)
          stat=0; nl=0
          do while (stat == 0)
            read(10,*,iostat=stat)
            nl = nl + 1
          enddo 
          nl = nl - 1
          rewind(10)
          read(10,*)
          allocate(zref(1:nl), ptref(1:nl), qtref(1:nl))
          do l = 1, nl
            read(10,*) zref(l), ptref(l), qtref(l)
          enddo
          close(10)
        endif
  end 
!
end
!
! ===============================================
  subroutine pbl_mixing ( dens, height )
! ===============================================
!
  real, intent(in) :: height
  real, dimension(ip_start:,jp_start:,:), intent(in) :: dens
!
  integer :: k, h
  real :: zmin=5.
  real, dimension(1:nz) :: zero
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: flux
!
!-----------------------------------------------------------!
!              	        Surface fluxes                      !
!-----------------------------------------------------------!
!
  if ( isurf >= 0 ) call calc_surface ( pressure2, wind2, state2, hydromtr2 )
!
!----------------------------------------------------------!
!                Calculate eddy viscosity                  !
!----------------------------------------------------------!
!
  turbu%fkv = 0.
  do k = 2, nz
    if (z0(k) <= height) turbu%fkv(:,:,k) = vk*surf%ustar(:,:)*z0(k)*(1. - z0(k)/height)**2.
  enddo
  turbu%fkv(:,:,1) = vk*surf%ustar(:,:)*zmin*(1. - zmin/height)**2.
!
!----------------------------------------------------------!
!                       SGS Fluxes                         !
!----------------------------------------------------------!
!
  call eff1d (1, -surf%esflux, pt0, states%es, state2%es, turbu%fkv)
!
  call eff1d (2, -surf%qvflux, qt0, states%qt, state2%qt, turbu%fkv)
!
  zero = 0.
  flux = 0.
!
  if (micro_dif) then
    do h = 1, nhydro
!
      call eff1d (2+h, flux, zero, hydromtrs(h)%q, hydromtr2(h)%q, turbu%fkv)
!
#ifdef SEIFERT
      if (moments == 2) call eff1d (2+nhydro+h, flux, zero, hydromtrs(h)%n, hydromtr2(h)%n, turbu%fkv)
#endif
    enddo 
  endif
!
  if (nscal > 0) then
    do h = 1, 1 !nscal
      call eff1d (2+2*nhydro+h, flux, zero, states%scal(:,:,:,h), state2%scal(:,:,:,h), turbu%fkv)
    enddo
  endif
!
!-----------------------------------------------------------!
!
return
end
!
!   ================================================
    subroutine eff1d ( flag, flux_s, x0, ec, x, fkv )
!   ================================================

! ----------------------------------------------
! --- Subroutine calculating eddy terms of 
! ---            scalar variables
! ----------------------------------------------
!
  integer :: i, j, k, flag
  real, dimension(:), intent(in) :: x0
  real, dimension(ip_start:,jp_start:), intent(in) :: flux_s
  real, dimension(ip_start:,jp_start:,:), intent(in) :: x, fkv
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: ec
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: ez, gradz
!
!----------------------------------------------------------!
!          	  Calculate scalar gradient                !
!----------------------------------------------------------!
!
  ez = 0.
!
  do j = jt_start, jt_end
    do i = it_start, it_end
      call grad_1d ( x(i,j,:), gradz(i,j,:) )
    enddo
  enddo
!
!  Diagnostics
!
  if (out_diagt .and. out_grad .and. flag==1) diag(3)%grad%pt = gradz
  if (out_diagq .and. out_grad .and. flag==2) diag(3)%grad%qt = gradz
 
!
!----------------------------------------------------------!
!              Calculate scalar SGS fluxes                 !
!----------------------------------------------------------!
!
  do k = 1, nz-1
    gradz(:,:,k) = 0.5*(den0(k+1)*fkv(:,:,k+1) + den0(k)*fkv(:,:,k)) * gradz(:,:,k+1)
  enddo
!
  do k = 2, nz-1
    ez(:,:,k) = (gradz(:,:,k) - gradz(:,:,k-1))*dz1w(k)
  enddo
  ez(:,:,1) = (gradz(:,:,1) - flux_s)*dz1w(1)
!
!  Non-conservative case: divide by density afterwards
!
#ifndef CONSERVATIVE
  do k = 1, nz
    ez(:,:,k) = ez(:,:,k) / den0(k)
  enddo
#endif
!
!----------------------------------------------------------!
!          	         Assemble    		           !
!----------------------------------------------------------!
!
!  Diagnostics
!
  if (out_diagt .and. flag==1) diag(5)%pt = diag(5)%pt + ez
  if (out_diagq .and. flag==2) diag(5)%qt = diag(5)%qt + ez
!
  if (out_diagl .and. flag==3) diag(5)%qc = diag(5)%qc + ez
  if (out_diagr .and. flag==4) diag(5)%qr = diag(5)%qr + ez
!
  if (out_diags .and. flag > 2+2*nhydro) then
    diag(5)%sca(:,:,:,flag-2*(nhydro+1)) = diag(5)%sca(:,:,:,flag-2*(nhydro+1)) + ez
  endif
!
  ec = ec + ez
!
!----------------------------------------------------------!
!
return                                                                   
end
!
end module columnphys
