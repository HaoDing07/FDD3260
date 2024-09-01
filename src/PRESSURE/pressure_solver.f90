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
!  PRESSURE_SOLVER.F:                   
!
!  Purpose:
!      Diagnostic evaluation of pressure: solving Poisson (anelastic
!	   case) or a generalized Helmholtz equation (compressible case)
!      using a multigrid method from the MUDPACK open-source library. 
!
!      - "tend" represents the gradient of momentum tendency      
!      - "pp" or "pp1" is a correction pressure which enforces div(u) = 0. 
!
!  Author
!      Julien Savre, MISU
!
! ================================================================
  
  module pressure_solver
        
!
!-----------------------------------------------------------!
!
USE shared_all
USE advection
USE gradients
USE averages
USE allocation
USE funcpack
USE fftmpi_solver
USE fft_solvers
USE fft_solvers_nnf
USE fft_solvers_open
USE boundary_conditions
USE, intrinsic :: iso_c_binding 
!
#ifdef SPMD
USE mpi
USE mpicomm
#endif
!
IMPLICIT NONE

  private

  public :: anelastic_pressure, pressure_diagnostics

  contains
!
!      ===============================================
       subroutine anelastic_pressure ( wold ) 
!      ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
    integer :: ip
    real :: m
    real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: wold
    real, allocatable, dimension(:,:,:) :: deltap, div
!
#if ( defined SPMD )      
    integer :: ierr, status(MPI_STATUS_SIZE)
#endif
!
  if (verbose > 0) call write_debug('Start anelastic_pressure')
!
!  Allocate and initialize
!
    maxdiv = 0.
!
    call alloc ( deltap )
    call alloc ( div )
!
!-----------------------------------------------------------!
!          	     Solve for pressure		            !
!-----------------------------------------------------------!
!
   if ( nsubp > 0 ) then
     do ip = 1, nsubp
!
!  Calculate non-hydrostatic pressure
!
       call solve_pressure ( deltap )
!
!  Apply pressure gradients
!
       call correction ( deltap )
     enddo
   endif
!
!  Hydrostatic pressure correction (so that <w> = 0)
!
  call cal_p1 ( wind%w, wold )
!
!----------------------------------------------------------!
!       Recalculate residual for conservation check        !
!----------------------------------------------------------!
!
   call caldiv ( wind, div )
!
   if ( out_diagp ) diag(7)%p = div
!
   m = maxval( div(it_start:it_end,jt_start:jt_end,2:nz-1) )
#ifdef SPMD
   call MPI_ALLReduce (m, maxdiv, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
   maxdiv = m
#endif
!
   call dealloc ( deltap )
   call dealloc ( div )
!        
!-----------------------------------------------------------!
!
  if (verbose > 0) call write_debug('Terminate anelastic_pressure')
!
return
end subroutine anelastic_pressure
!
!      ===============================================
       subroutine solve_pressure ( deltap )
!      ===============================================

  real, dimension (ip_start:,jp_start:,:), intent(out) :: deltap
!
  integer :: k
  real :: dt1
  real, allocatable, dimension(:,:,:) :: tend
!
!  Allocate and initialize
!
  call alloc ( tend )
!
  dt1 = 1./dt0
!
!-----------------------------------------------------------!
!           Compute divergence for pressure solver          !
!-----------------------------------------------------------!
!
  call caldiv ( wind, tend )
!
  if (out_diagp) diag(1)%p = tend
!
  do k=1,nz    
    tend(:,:,k) = dt1 * tend(:,:,k) / den0(k)
  enddo
!
  if ( out_diagp ) diag(2)%p = tend
!
!-----------------------------------------------------------!
!           	   Call pressure solvers	            !
!-----------------------------------------------------------!
!
  call fft_solve ( tend, deltap )
!
!  Correct dynamical pressure (constant of integration)
!
  call rescale_pressure ( deltap )
!
  if ( out_diagp ) diag(3)%p = deltap
!
  call dealloc ( tend )
!        
!-----------------------------------------------------------!
!
return
end subroutine solve_pressure
!
!   ===============================================
  subroutine pressure_diagnostics ( buoy, windl )
!   ===============================================
!
! --- Solving Poisson equation to evaluate pressure diagnostically
!
  integer :: k
  real, dimension (ip_start:,jp_start:,:), intent(in) :: buoy
  type(atm_winds), intent(in) :: windl
!
  real, allocatable, dimension (:,:,:) :: divb, beff
  real, allocatable, dimension (:,:,:) :: deltap, deltap1, deltap2, div
  type(atm_winds) :: wind1
!
#ifdef SPMD
  include 'fftw3-mpi.f03'
#else
  include 'fftw3.f03'
#endif
!
!-----------------------------------------------------------!
!                Conmpute effective buoyancy                !
!-----------------------------------------------------------!
!
  if (out_beff) then
!
    call alloc ( wind1 )
    call alloc ( divb )
    call alloc ( beff )
!
    call grad_quick ( buoy, wind1%u, wind1%v, wind1%w )
!
    wind1%w = 0.
!
    call caldiv ( wind1, divb )
!
!  Call pressure solvers
!
    call fft_solve ( divb, beff )
!
  endif
!
!-----------------------------------------------------------!
!           Conmpute non-hydrostatic pressure pert.         !
!-----------------------------------------------------------!
!
  if (out_dp) then
!
    call alloc ( div )
    call alloc ( deltap )
    call alloc ( deltap1 )
    call alloc ( deltap2 )
!
    deltap = 0.0
    wind1%u = windl%u
    wind1%v = windl%v
    wind1%w = windl%w - buoy
!
    call caldiv ( wind1, div )
!
    do k=1,nz    
      div(:,:,k) = (div(:,:,k) + divb(:,:,k)) / den0(k)
    enddo
!
!  Call pressure solvers
!
    call fft_solve ( div, deltap )
!
!-----------------------------------------------------------!
!              Compute buoyant pressure pert.               !
!-----------------------------------------------------------!
!
    deltap1 = 0.0
    wind1%u = 0.0
    wind1%v = 0.0
    wind1%w = buoy
!
    call caldiv ( wind1, div )
!
    do k=1,nz    
      div(:,:,k) = div(:,:,k) / den0(k)
    enddo
!
!  Call pressure solvers
!
    call fft_solve ( div, deltap1 )
!
!-----------------------------------------------------------!
!             Conmpute dynamic pressure pert.               !
!-----------------------------------------------------------!
!
    deltap2 = 0.0
    wind1%u = windl%u
    wind1%v = windl%v
    wind1%w = windl%w - buoy
!
    call caldiv ( wind1, div )
!
    do k=1,nz    
      div(:,:,k) = div(:,:,k) / den0(k)
    enddo
!
!  Call pressure solvers
!
    call fft_solve ( div, deltap2 )
!
    call dealloc ( div )
    call dealloc ( divb )
    call dealloc ( wind1 )
  endif
!
!-----------------------------------------------------------!
!                     Store diagnostics                     !
!-----------------------------------------------------------!
!
  if (out_beff) then
    call rescale_pressure ( beff )
!
    diagnos%beff = beff
  endif
!
  if (out_dp) then
    call rescale_pressure ( deltap )
!
    call rescale_pressure ( deltap1 )
!
    call rescale_pressure ( deltap2 )
!
    do k=1,nz    
      diagnos%pnh(:,:,k) = den0(k) * deltap(:,:,k)
      diagnos%pb(:,:,k) = den0(k) * deltap1(:,:,k)
      diagnos%pd(:,:,k) = den0(k) * deltap2(:,:,k)
    enddo
  endif
!
!  Deallocate
!
  if (out_beff) call dealloc ( beff )  
!
  if (out_dp) then
    call dealloc ( deltap )
    call dealloc ( deltap1 )
    call dealloc ( deltap2 )
  endif
!        
!-----------------------------------------------------------!
!
return
end subroutine pressure_diagnostics
!
!      ===============================================
       subroutine correction ( deltap )
!      ===============================================
!
! --- Update velotities using new pressure
!
  integer :: k
  real, dimension (ip_start:,jp_start:,:), intent(in) :: deltap
  real, dimension (ip_start:ip_end,jp_start:jp_end,1:nz) :: gpx, gpy, gpz
!
!----------------------------------------------------------!
!               Calculate pressure gradient                !
!----------------------------------------------------------!
!
  do k = 1, nz
    pressure%p(:,:,k) = pressure%p(:,:,k) + den0(k)*deltap(:,:,k)
  enddo
!
!----------------------------------------------------------!
!                    Add pressure term                     !
!        The pressure terms are directly integrated        !
!           to give the final velocities at n+1            !
!     Pressure terms already multiplied by the time step   !
!----------------------------------------------------------!
!
!  Find pressure gradient
!
  call gradp ( deltap, wind%w, gpx, gpy, gpz )
!
!  Update velocity
!
  do k=1,nz
    wind%u(:,:,k) = wind%u(:,:,k) - dt0*den0(k)*gpx(:,:,k)
!
#ifdef MODEL_3D
    wind%v(:,:,k) = wind%v(:,:,k) - dt0*den0(k)*gpy(:,:,k)
#endif
!
    wind%w(:,:,k) = wind%w(:,:,k) - dt0*avden0(k)*gpz(:,:,k)
  end do
!
!  BCs 
!
  call windbc ( wind )
!
!  Diagnostics
!
  if (ldiag) then
    if (out_diagp) then
      diag(4)%p = gpx
#ifdef MODEL_3D
      diag(5)%p = gpy
#endif
      diag(6)%p = gpz
    endif
!
    do k = 1, nz  
      if (out_diagu) diag(7)%u(:,:,k) = diag(7)%u(:,:,k) - cdiag*den0(k)*gpx(:,:,k)
      if (out_diagu) diag(9)%u(:,:,k) = diag(9)%u(:,:,k) - cdiag*den0(k)*gpx(:,:,k)
!
#ifdef MODEL_3D
      if (out_diagv) diag(7)%v(:,:,k) = diag(7)%v(:,:,k) - cdiag*den0(k)*gpy(:,:,k)
      if (out_diagv) diag(9)%v(:,:,k) = diag(9)%v(:,:,k) - cdiag*den0(k)*gpy(:,:,k)
#endif
!
      if (out_diagw) diag(7)%w(:,:,k) = diag(7)%w(:,:,k) - cdiag*avden0(k)*gpz(:,:,k)
      if (out_diagw) diag(9)%w(:,:,k) = diag(9)%w(:,:,k) - cdiag*avden0(k)*gpz(:,:,k)
    enddo
  endif
!
return
end subroutine correction
!
!  ===========================================
   subroutine rescale_pressure ( pp )
!  ===========================================
!
! ---------------------------------
! --- Correct pressure obtained from FFTW solver: using FFTW, the
!	pressure is calculated to a constant of integration and is 
!	not normalized. This is handled here by ensuring total energy
!	conservation in the domain (to get the constant of integration)
!	and renormalizing the solution.
! ---------------------------------
!
  real, dimension(ip_start:,jp_start:,:), intent(inout) :: pp
!
  integer :: k
  real :: pav(nz)
!
!----------------------------------------------------------!
!            Rescale pressure from FFT solver   	   !
!----------------------------------------------------------!
!
!  Divide by total number of points in plane as FFTW does not normalize:
!  Careful, the 2D FFT is performed differently in serial and parallel
!
#ifdef MODEL_3D
#ifndef CHANNEL
  pp = pp / real( (nx-5)*(ny-5) )
#else
  pp = pp / real( 2*(nx-5)*(ny-5) )
#endif
#else
  pp = pp / real( nx-5 )
#endif
!
! Make sure perturbation averages to 0
!
  call horav(pp, pav)
!
  do k = 1, nz
    pp(:,:,k) = pp(:,:,k) - pav(k)
  enddo
!
! ---------------------------------------------------------!
!
return
end subroutine rescale_pressure
!       
! ===============================================
  subroutine fft_solve ( tend, output )
! ===============================================
!
  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: tend, output
  logical :: fft1d, fftmpi
!
#ifdef SPMD
  include 'fftw3-mpi.f03'
#else
  include 'fftw3.f03'
#endif
!
  fft1d = .false.
  fftmpi = .false.
  if (fft_solv == 1) then
    fft1d = .true.
  else if (fft_solv == 2) then
    fftmpi = .true.
  endif
!
!  Initialize memory
!
#ifdef SPMD
  call fftw_mpi_init ()
#endif
!  
!  Call appropriate fft_solve routine
!
#ifdef SPMD
!
  if (.not.nest_run) then
#if (defined MODEL_3D)
#if (defined CHANNEL)
    call fft_solve_nnf_parallel ( tend, output )
#else
    if (fftmpi) then
      call fftmpi_solve ( tend, output )
    else
      if ( mod(nx-5,nprocx*nprocy) /= 0 .or. fft1d ) then
        call fft_solve_parallel_2x1d ( tend, output )
      else
        call fft_solve_parallel ( tend, output )
      endif
    endif
#endif
#else
    call fft_solve_parallel ( tend, output )
#endif
  else
    call fft_solve_parallel_open ( tend, output )
  endif
!
#else
!
#if (defined CHANNEL)
  call fft_solve_nnf_serial ( tend, output )
#else
  if (.not.nest_run) then
    call fft_solve_serial ( tend, output )
  else
    call fft_solve_serial_open ( tend, output )
  endif
#endif
!
#endif
!
!  Clean planes and memory
!
#ifdef SPMD
  call fftw_mpi_cleanup ()
#else
  call fftw_cleanup ()
#endif
!
!-----------------------------------------------------------!
!
  end subroutine fft_solve
  
end module pressure_solver
