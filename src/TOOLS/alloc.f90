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
!  ALLOC:
!	Package containing allocation/deallocation functions
!	for clean definition of local arrays.
!
!  Author:
!	Julien Savre, Ludwig-Maximilian-Universitat, Munich
!
! ================================================================

module allocation

USE gridno
USE shared_data
USE shared_diag
USE shared_surf
USE shared_nuclei
!
IMPLICIT NONE
!
  private
!
  interface alloc
      MODULE PROCEDURE allocate5d, allocate4d, allocate3d, allocate2d, allocate1d
  end interface alloc
!
  interface dealloc
      MODULE PROCEDURE deallocate5d, deallocate4d, deallocate3d, deallocate2d, deallocate1d
  end interface dealloc
!
  interface reset
      MODULE PROCEDURE reset3d, reset4d, resetmicro, resetin, resetaero, resetdiag, resetcol
  end interface reset
!
  public :: alloc, dealloc, reset
!
CONTAINS
!
  subroutine allocate5d ( data, dim1, dim2 )
  
    integer, intent(in) :: dim1, dim2
    real, allocatable, dimension(:,:,:,:,:), intent(inout) :: data
    
    if ( .not.allocated(data) ) then
      allocate( data(ip_start:ip_end,jp_start:jp_end,1:nz,1:dim1,1:dim2) )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to allocate a 5d array already allocated'
    endif

    data(ip_start:ip_end,jp_start:jp_end,1:nz,1:dim1,1:dim2) = 0.
    
  end subroutine allocate5d
!
  subroutine deallocate5d ( data )
  
    real, allocatable, dimension(:,:,:,:,:), intent(inout) :: data
    
    if ( allocated(data) ) then
      deallocate( data )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to deallocate unallocated 5d array'
    endif
    
  end subroutine deallocate5d
!
  subroutine allocate4d ( data, dim )
  
    integer, intent(in) :: dim
    real, allocatable, dimension(:,:,:,:), intent(inout) :: data
    
    if ( .not.allocated(data) ) then
      allocate( data(ip_start:ip_end,jp_start:jp_end,1:nz,1:dim) )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to allocate a 4d array already allocated'
    endif

    data(ip_start:ip_end,jp_start:jp_end,1:nz,1:dim) = 0.
    
  end subroutine allocate4d
!
  subroutine deallocate4d ( data )
  
    real, allocatable, dimension(:,:,:,:), intent(inout) :: data
    
    if ( allocated(data) ) then
      deallocate( data )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to deallocate unallocated 4d array'
    endif
    
  end subroutine deallocate4d
!
  subroutine allocate3d ( data )
  
    real, allocatable, dimension(:,:,:), intent(inout) :: data
    
    if ( .not.allocated(data) ) then
      allocate( data(ip_start:ip_end,jp_start:jp_end,1:nz) )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to allocate a 3d array already allocated'
    endif

    data(ip_start:ip_end,jp_start:jp_end,1:nz) = 0.
    
  end subroutine allocate3d
!
  subroutine deallocate3d ( data )
  
    real, allocatable, dimension(:,:,:), intent(inout) :: data
    
    if ( allocated(data) ) then
      deallocate( data )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to deallocate unallocated 3d array'
    endif
    
  end subroutine deallocate3d
!
  subroutine allocate2d ( data )
  
    real, allocatable, dimension(:,:), intent(inout) :: data
    
    if ( .not.allocated(data) ) then
      allocate( data(ip_start:ip_end,jp_start:jp_end) )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to allocate a 2d array already allocated'
    endif

    data(ip_start:ip_end,jp_start:jp_end) = 0.
    
  end subroutine allocate2d
!
  subroutine deallocate2d ( data )
  
    real, allocatable, dimension(:,:), intent(inout) :: data
    
    if ( allocated(data) ) then
      deallocate( data )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to deallocate unallocated 2d array'
    endif
    
  end subroutine deallocate2d
!
  subroutine allocate1d ( data, dim )
  
    real, allocatable, dimension(:), intent(inout) :: data
    integer, intent(in) :: dim

    if ( .not.allocated(data) ) then
      allocate( data(1:dim) )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to allocate 1d array already allocated'
    endif

    data(:) = 0.
    
  end subroutine allocate1d
!
  subroutine deallocate1d ( data )
  
    real, allocatable, dimension(:), intent(inout) :: data
    
    if ( allocated(data) ) then
      deallocate( data )
    else
      if (mypid.eq.0) write(7,*) 'WARNING: Trying to deallocate unallocated 2d array'
    endif
    
  end subroutine deallocate1d
!
  subroutine reset3d ( data )
  
    real, dimension(:,:,:), intent(inout) :: data
    
    data = 0.
    
  end subroutine reset3d
!
  subroutine reset4d ( data )
  
    real, dimension(:,:,:,:), intent(inout) :: data
    
    data = 0.
    
  end subroutine reset4d
!
  subroutine resetin ( data )
  
    type (nuclei_3d), dimension(:), intent(inout) :: data
    integer :: h, l
    
    do h = 1, 3
      do l = 1, 3
        data(h)%mode(l)%n = 0.
        data(h)%mode(l)%nc = 0.
#ifdef NUC_CNT1  
        data(h)%mode(l)%m = 0.
        data(h)%mode(l)%mc = 0.
#endif
      enddo
    enddo
    
  end subroutine resetin
!
  subroutine resetmicro ( data )
  
    type (microdiag), dimension(:), intent(inout) :: data
    integer :: h
    
    do h = 1, nhydro
      data(h)%q = 0.
      data(h)%n = 0.
    enddo
    
  end subroutine resetmicro
!
  subroutine resetaero ( data )
  
    type (aerodiag), dimension(:), intent(inout) :: data
    integer :: h
    
    do h = 1, nmode
      data(h)%n = 0.
      data(h)%m = 0.
    enddo
    
  end subroutine resetaero
!
  subroutine resetdiag ( data )
  
    type (tendencies), dimension(:), intent(inout) :: data
    integer :: h

    do h = 1, ndiag
      if (out_diagu) call reset3d ( data(h)%u )
      if (out_diagv) call reset3d ( data(h)%v )
      if (out_diagw) call reset3d ( data(h)%w )
      if (out_diagq) call reset3d ( data(h)%qt )
      if (out_diagt) call reset3d ( data(h)%pt )
      if (out_diagtv) call reset3d ( data(h)%ptv )
      if (out_diagp) call reset3d ( data(h)%p )
      if (out_diagk) call reset3d ( data(h)%k )
      if (out_diags.and.nscal > 0) call reset4d ( data(h)%sca )
      if (out_diagl.and.lmicro > 0) call reset3d ( data(h)%qc )
      if (out_diagr.and.lmicro > 0) call reset3d ( data(h)%qr )
      if (out_diagi.and.lmicro > 1) call reset3d ( data(h)%qi )
      if (out_micro.and.lmicro > 0) call resetmicro ( data(h)%micro )
      if (out_diaga) call resetaero ( data(h)%aero )
    enddo

  end subroutine resetdiag
!
  subroutine resetcol ( column )
 
    type (columns), intent(inout) :: column
  
    if (out_lwp) column%lwp = 0. 
    if (out_lwp) column%iwp = 0. 
    if (out_cwp) column%cwp = 0. 
    if (out_cwp) column%rwp = 0. 
    if (out_wvp) column%wvp = 0. 
    if (out_wvp) column%wvpb = 0. 
    if (out_wvp) column%wvpt = 0. 
    if (out_cmse) column%cmse = 0.
    if (out_cmse) column%cmseb = 0.
    if (out_cmse) column%cmset = 0.
    if (out_cmfl) column%cmfl = 0.
    if (out_ctop) column%ctop = 0.
    if (out_ctop) column%cbas = 0.
    if (out_zinv) column%bltop = 0.
    if (out_zinv) column%zinv = 0.
    if (out_cape) column%cape = 0.
    if (out_cape) column%cin = 0.
    if (out_cape) column%lfc = 0.
    if (out_cape) column%lcl = 0.
    if (out_cp) column%cpint = 0.
    if (out_cp) column%cpint0 = 0.
    if (out_ints) column%intsca = 0.

  end subroutine resetcol
!    
end module allocation
