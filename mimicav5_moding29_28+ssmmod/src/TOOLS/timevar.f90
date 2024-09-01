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
!  TIMEVAR:
!	Updates potential time varying quantities
!
!  Author:
!	Julien Savre, MISU
!
! ================================================================	 
!
!  ===================================================
   subroutine timevar
!  ===================================================
!
  USE gridno
  USE shared_data
  USE surface_var
  USE hygro_growth
  ! USE nudge_adv
!
IMPLICIT NONE
!
  ! if (verbose > 0) call write_debug('Starting timevar')
  if (verbose > 0) write(7,*) 'Starting timevar'
!
!----------------------------------------------------------!
!    Call functions to calculate time varying variables    !
!----------------------------------------------------------!
!
  if ( with_dif ) call time_surf 
!
  if ( with_lsadv ) call time_lsadv
!
  if (forcing_surface) call surface_forcing

  if (with_aero_swelling) call aerosol_swelling

  ! if (with_warmadv) call nudge_warm_adv


!----------------------------------------------------------!
!
  ! if (verbose > 0) call write_debug('Terminating timevar')
  if (verbose > 0) write(7,*) 'Terminating timevar'
!
return
end
!

!  ===================================================
subroutine nudge_warm_adv
  !  ===================================================
  
    USE gridno
    USE shared_data
    USE shared_state
    USE shared_all
    integer ::  l, m, n, i, j, k
  
  
  
    if (with_warmadv .and. time == tstart1_warmadv) then

      write(7,*) 'Reading warm advection files 1:', warm_adv_H1, warm_adv_dT1
  
      open(333, file=warm_adv_H1, status='old')
      n = 0
      do
          read(333,*,iostat=m)  ! read the data of next line
          if (m /= 0) exit     ! if no data, quit the loop
          n = n + 1            
      end do
      if (.not. allocated(H_warm_adv1)) then 
        allocate(H_warm_adv1(n))       ! allocate the array
      end if
      rewind(333)              ! let the pointer be at the start of the file
      do m = 1, n
          read(333,*) H_warm_adv1(m)   ! store the data in ts_force array
      end do
      close(333)
    
      open(444, file=warm_adv_dT1, status='old')
      n = 0
      do
          read(444,*,iostat=m)  ! read the data of next line
          if (m /= 0) exit     ! if no data, quit the loop
          n = n + 1            
      end do
      if (.not. allocated(dT_warm_adv1)) then 
        allocate(dT_warm_adv1(n))       ! allocate the array
      end if
      rewind(444)              ! let the pointer be at the start of the file
      do m = 1, n
          read(444,*) dT_warm_adv1(m)   ! store the data in ts_force array
      end do
      close(444)
  
      write(7,*) 'Starting nudging warm advection 1'
    endif

    if (with_warmadv .and. time >= tstart1_warmadv .and. time <= tstop1_warmadv) then
 
      do k = int(H_warm_adv1(1))/dz,  int(H_warm_adv1(size(H_warm_adv1)))/dz
        do j = jp_start, jp_end
          do i = ip_start, ip_end
            state2%es(i,j,k) = state2%es(i,j,k)+ dT_warm_adv1((int(k*dz-H_warm_adv1(1)))/dz+1)
          enddo
        enddo
      enddo
      write(7,*) 'nudging warm advection 1 is working at this timestep'
    endif 

    !!! Period2
    if (with_warmadv .and. time == tstart2_warmadv) then
      write(7,*) 'Reading warm advection 2 files:', warm_adv_H2, warm_adv_dT2
  
      open(344, file=warm_adv_H2, status='old')
      n = 0
      do
          read(344,*,iostat=m)  ! read the data of next line
          if (m /= 0) exit     ! if no data, quit the loop
          n = n + 1            
      end do
      if (.not. allocated(H_warm_adv2)) then 
        allocate(H_warm_adv2(n))       ! allocate the array
      end if
      rewind(344)              ! let the pointer be at the start of the file
      do m = 1, n
          read(344,*) H_warm_adv2(m)   ! store the data in ts_force array
      end do
      close(344)
    
      open(455, file=warm_adv_dT2, status='old')
      n = 0
      do
          read(455,*,iostat=m)  ! read the data of next line
          if (m /= 0) exit     ! if no data, quit the loop
          n = n + 1            
      end do
      if (.not. allocated(dT_warm_adv2)) then 
        allocate(dT_warm_adv2(n))       ! allocate the array
      end if
      rewind(455)              ! let the pointer be at the start of the file
      do m = 1, n
          read(455,*) dT_warm_adv2(m)   ! store the data in ts_force array
      end do
      close(455)
  
      write(7,*) 'Starting nudging warm advection 2'
    endif

    if (with_warmadv .and. time >= tstart2_warmadv .and. time <= tstop2_warmadv) then
      ! write(7,*) 'thermo%T(1,1,90) =', thermo%T(1,1,90)
      do k = int(H_warm_adv2(1))/dz,  int(H_warm_adv2(size(H_warm_adv2)))/dz
        do j = jp_start, jp_end
          do i = ip_start, ip_end
            state2%es(i,j,k) = state2%es(i,j,k)+ dT_warm_adv2((int(k*dz-H_warm_adv2(1)))/dz+1)
          enddo
        enddo
      enddo
      write(7,*) 'nudging warm advection 2 is working at this timestep'
    endif

    !!! Period3
    if (with_warmadv .and. time == tstart3_warmadv) then
      write(7,*) 'Reading warm advection 3 files:', warm_adv_H3, warm_adv_dT3
  
      open(355, file=warm_adv_H3, status='old')
      n = 0
      do
          read(355,*,iostat=m)  ! read the data of next line
          if (m /= 0) exit     ! if no data, quit the loop
          n = n + 1            
      end do
      if (.not. allocated(H_warm_adv3)) then 
        allocate(H_warm_adv3(n))       ! allocate the array
      end if
      rewind(355)              ! let the pointer be at the start of the file
      do m = 1, n
          read(355,*) H_warm_adv3(m)   ! store the data in ts_force array
      end do
      close(355)
    
      open(466, file=warm_adv_dT3, status='old')
      n = 0
      do
          read(466,*,iostat=m)  ! read the data of next line
          if (m /= 0) exit     ! if no data, quit the loop
          n = n + 1            
      end do
      if (.not. allocated(dT_warm_adv3)) then 
        allocate(dT_warm_adv3(n))       ! allocate the array
      end if
      rewind(466)              ! let the pointer be at the start of the file
      do m = 1, n
          read(466,*) dT_warm_adv3(m)   ! store the data in ts_force array
      end do
      close(466)
  
      write(7,*) 'Starting nudging warm advection 3'
    endif

    if (with_warmadv .and. time >= tstart3_warmadv .and. time <= tstop3_warmadv) then
      ! write(7,*) 'thermo%T(1,1,90) =', thermo%T(1,1,90)
      do k = int(H_warm_adv3(1))/dz,  int(H_warm_adv3(size(H_warm_adv3)))/dz
        do j = jp_start, jp_end
          do i = ip_start, ip_end
            state2%es(i,j,k) = state2%es(i,j,k)+ dT_warm_adv3((int(k*dz-H_warm_adv3(1)))/dz+1)
          enddo
        enddo
      enddo
      write(7,*) 'nudging warm advection 3 is working at this timestep'
    endif 


    if (with_warmadv .and. time == tstop3_warmadv) write(7,*) 'Terminating nudging warm advection'
    return
    end subroutine nudge_warm_adv
  
  !! Nudging q advection (po valley case as a example)
  !  ===================================================
    subroutine nudge_q_adv
      !  ===================================================
      
        USE gridno
        USE shared_data
        USE shared_state
        USE shared_all
        integer ::  l, m, n, i, j, k
      
      
      
        if (with_warmadv .and. time == tstart1_warmadv) then
          write(7,*) 'Reading q advection files 1:', q_adv_H1, q_adv_dQ1
      
          open(555, file=q_adv_H1, status='old')
          n = 0
          do
              read(555,*,iostat=m)  ! read the data of next line
              if (m /= 0) exit     ! if no data, quit the loop
              n = n + 1            
          end do
          if (.not. allocated(H_q_adv1)) then 
            allocate(H_q_adv1(n))       ! allocate the array
          end if
          rewind(555)              ! let the pointer be at the start of the file
          do m = 1, n
              read(555,*) H_q_adv1(m)   ! store the data in ts_force array
          end do
          close(555)
        
          open(666, file=q_adv_dQ1, status='old')
          n = 0
          do
              read(666,*,iostat=m)  ! read the data of next line
              if (m /= 0) exit     ! if no data, quit the loop
              n = n + 1            
          end do
          if (.not. allocated(dQ_q_adv1)) then 
            allocate(dQ_q_adv1(n))       ! allocate the array
          end if
          rewind(666)              ! let the pointer be at the start of the file
          do m = 1, n
              read(666,*) dQ_q_adv1(m)   ! store the data in ts_force array
          end do
          close(666)
      
          write(7,*) 'Starting nudging q advection 1'
      
        endif
        if (with_warmadv .and. time >= tstart1_warmadv .and. time <= tstop1_warmadv) then
          ! write(7,*) 'thermo%T(1,1,90) =', thermo%T(1,1,90)
          do k = int(H_q_adv1(1))/dz,  int(H_q_adv1(size(H_q_adv1)))/dz
            do j = jp_start, jp_end
              do i = ip_start, ip_end
                ! write(7,*) 'check', k, dz,H_q_adv1(1),int(k*dz-H_q_adv1(1)), (int(k*dz-H_q_adv1(1)))/dz
                state2%qt(i,j,k) = state2%qt(i,j,k)+ dQ_q_adv1((int(k*dz-H_q_adv1(1)))/dz+1)
              enddo
            enddo
          enddo
          write(7,*) 'nudging q advection 1 is working at this timestep'
        endif 

        !!! Period 2
        if (with_warmadv .and. time == tstart2_warmadv) then
          write(7,*) 'Reading q advection files 2:', q_adv_H2, q_adv_dQ2
      
          open(566, file=q_adv_H2, status='old')
          n = 0
          do
              read(566,*,iostat=m)  ! read the data of next line
              if (m /= 0) exit     ! if no data, quit the loop
              n = n + 1            
          end do
          if (.not. allocated(H_q_adv2)) then 
            allocate(H_q_adv2(n))       ! allocate the array
          end if
          rewind(566)              ! let the pointer be at the start of the file
          do m = 1, n
              read(566,*) H_q_adv2(m)   ! store the data in ts_force array
          end do
          close(566)
        
          open(677, file=q_adv_dQ2, status='old')
          n = 0
          do
              read(677,*,iostat=m)  ! read the data of next line
              if (m /= 0) exit     ! if no data, quit the loop
              n = n + 1            
          end do
          if (.not. allocated(dQ_q_adv2)) then 
            allocate(dQ_q_adv2(n))       ! allocate the array
          end if
          rewind(677)              ! let the pointer be at the start of the file
          do m = 1, n
              read(677,*) dQ_q_adv2(m)   ! store the data in ts_force array
          end do
          close(677)
      
          write(7,*) 'Starting nudging q advection 2 '
      
        endif

        if (with_warmadv .and. time >= tstart2_warmadv .and. time <= tstop2_warmadv) then
          ! write(7,*) 'thermo%T(1,1,90) =', thermo%T(1,1,90)
          do k = int(H_q_adv2(1))/dz,  int(H_q_adv2(size(H_q_adv2)))/dz
            do j = jp_start, jp_end
              do i = ip_start, ip_end
                state2%qt(i,j,k) = state2%qt(i,j,k)+ dQ_q_adv2((int(k*dz-H_q_adv2(1)))/dz+1)
              enddo
            enddo
          enddo
          write(7,*) 'nudging q advection 2 is working at this timestep'
        endif 

      !!! Period 3
      if (with_warmadv .and. time == tstart3_warmadv) then
        write(7,*) 'Reading q advection files 2:', q_adv_H3, q_adv_dQ3
    
        open(577, file=q_adv_H3, status='old')
        n = 0
        do
            read(577,*,iostat=m)  ! read the data of next line
            if (m /= 0) exit     ! if no data, quit the loop
            n = n + 1            
        end do
        if (.not. allocated(H_q_adv3)) then 
          allocate(H_q_adv3(n))       ! allocate the array
        end if
        rewind(577)              ! let the pointer be at the start of the file
        do m = 1, n
            read(577,*) H_q_adv3(m)   ! store the data in ts_force array
        end do
        close(577)
      
        open(688, file=q_adv_dQ3, status='old')
        n = 0
        do
            read(688,*,iostat=m)  ! read the data of next line
            if (m /= 0) exit     ! if no data, quit the loop
            n = n + 1            
        end do
        if (.not. allocated(dQ_q_adv3)) then 
          allocate(dQ_q_adv3(n))       ! allocate the array
        end if
        rewind(688)              ! let the pointer be at the start of the file
        do m = 1, n
            read(688,*) dQ_q_adv3(m)   ! store the data in ts_force array
        end do
        close(688)
    
        write(7,*) 'Starting nudging q advection 3 '

      endif
      if (with_warmadv .and. time >= tstart3_warmadv .and. time <= tstop3_warmadv) then
        ! write(7,*) 'thermo%T(1,1,90) =', thermo%T(1,1,90)
        do k = int(H_q_adv3(1))/dz,  int(H_q_adv3(size(H_q_adv3)))/dz
          do j = jp_start, jp_end
            do i = ip_start, ip_end
              state2%qt(i,j,k) = state2%qt(i,j,k)+ dQ_q_adv3((int(k*dz-H_q_adv3(1)))/dz+1)
            enddo
          enddo
        enddo
        write(7,*) 'nudging q advection 3 is working at this timestep'
        
      endif 


      if (with_warmadv .and. time == tstop3_warmadv) write(7,*) 'Terminating nudging q advection'
      return
      end subroutine nudge_q_adv

!  ===================================================
   subroutine time_surf
!  ===================================================

USE gridno
USE shared_data
!
IMPLICIT NONE
!
  real :: ft
  logical, save :: lfirst=.true.
!
!----------------------------------------------------------!
!               Time varying surface fluxes:               !
!         Expressed as a linear function of time           !
!----------------------------------------------------------!
!
  select case (trim(casename))
!
    case ('TRANSIT')
      ft = max(0.,cos(0.5*pi*(5.25 - time/3600. + t_local)/5.25))
      shf = shf0 * (ft)**1.5
      lhf = lhf0 * (ft)**1.3
!
    case ('GERMANY')
      ft = max(0.,sin(pi*(time/3600. + t_local - 6.)/14.))
      shf = shf0*ft**6.5 - 10.
      lhf = lhf0*ft**6.5 + 10.
!
    case ('MIDLAT')
      ft = max(0.,sin(pi*(time/3600. + t_local - 6.)/12.))
      shf = shf0*ft
      lhf = lhf0*ft
!
    case default
      shf = shf0
      lhf = lhf0
!
  end select   
!
return
end
!
!
!  ===================================================
   subroutine time_lsadv
!  ===================================================

#ifdef SPMD
USE mpi
#endif
USE gridno
USE shared_data
!
IMPLICIT NONE
!
  logical  :: ladv=.false.
  integer  :: ierr, it
!
!----------------------------------------------------------!
!           Time varying large scale adveection            !
!----------------------------------------------------------!
!
  lsp = 0.0
  lsqt = 0.0
  lstem = 0.0
!
  select case (casename)
!
  case ('ISDAC')
  ladv = .true.
  it = time/10800.
!
  if (mypid == 0) then
    lsp(1) = 1.e+3*psurf**(cp/Ra)
    lsp(2) = 975.
    lsp(3) = 950.
    lsp(4) = 925.
    lsp(5) = 900.
    lsp(6) = 875.
    lsp(7) = 850.
    lsp(8) = 825.
!
!    if (it == 0) then
      lsqt0(1)   = 0.2e-7
      lsqt0(2)   = 0.2e-7
      lsqt0(3)   = 0.2e-7
      lsqt0(4)   = 0.2e-7
      lsqt0(5)   = 0.2e-7
      lsqt0(6)   = 0.2e-7
      lsqt0(7)   = 0.2e-7
      lsqt0(8)   = 0.2e-7
      lstem0(1) = 0.
      lstem0(2) = 0.
      lstem0(3) = 0.
      lstem0(4) = 0.
      lstem0(5) = 0.
      lstem0(6) = 0.
      lstem0(7) = 0.
      lstem0(8) = 0.   
!
      lsqt(1)   = 0.2e-7
      lsqt(2)   = 0.2e-7
      lsqt(3)   = 0.2e-7
      lsqt(4)   = 0.2e-7
      lsqt(5)   = 0.2e-7
      lsqt(6)   = 0.2e-7
      lsqt(7)   = 0.2e-7
      lsqt(8)   = 0.2e-7
      lstem(1) = 0.
      lstem(2) = 0.
      lstem(3) = 0.
      lstem(4) = 0.
      lstem(5) = 0.
      lstem(6) = 0.
      lstem(7) = 0.
      lstem(8) = 0.	
!    else if (it == 1) then
!      lsqt0  = lsqt
!      lstem0 = lstem
!      lsqt(1)   = 0.03e-7
!      lsqt(2)   = 0.06e-7
!      lsqt(3)   = 0.1e-7
!      lsqt(4)   = 0.13e-7
!      lsqt(5)   = 0.17e-7
!      lsqt(6)   = 0.17e-7
!      lsqt(7)   = 0.17e-7
!      lsqt(8)   = 0.17e-7
!      lstem(1) = 0.
!      lstem(2) = 0.
!      lstem(3) = 0.
!      lstem(4) = 0.
!      lstem(5) = 0.
!      lstem(6) = 0.
!      lstem(7) = 0.
!      lstem(8) = 0.	   
!    else if (it == 2) then
!      lsqt0  = lsqt
!      lstem0 = lstem
!      lsqt(1)   = 0.03e-7
!      lsqt(2)   = 0.06e-7
!      lsqt(3)   = 0.1e-7
!      lsqt(4)   = 0.13e-7
!      lsqt(5)   = 0.17e-7
!      lsqt(6)   = 0.17e-7
!      lsqt(7)   = 0.17e-7
!      lsqt(8)   = 0.17e-7
!      lstem(1) = 0.
!      lstem(2) = 0.
!      lstem(3) = 0.
!      lstem(4) = 0.
!      lstem(5) = 0.
!      lstem(6) = 0.
!      lstem(7) = 0.
!      lstem(8) = 0.	   
!    endif 
  endif
!  
  end select
!
!  Broadcast new large scale tendencies
!
#ifdef SPMD
  if (ladv) then
    call MPI_BCAST(lsp,      8, REALTYPE,0,MPI_COMM_WORLD,ierr)  
    call MPI_BCAST(lsqt,     8, REALTYPE,0,MPI_COMM_WORLD,ierr)      
    call MPI_BCAST(lsqt0,    8, REALTYPE,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(lstem,    8, REALTYPE,0,MPI_COMM_WORLD,ierr)      
    call MPI_BCAST(lstem0,   8, REALTYPE,0,MPI_COMM_WORLD,ierr)
  endif
#endif    
!
return
end
!
!
!  ===================================================
   subroutine add_pert 
!  ===================================================

USE gridno
USE shared_data
USE shared_state
#ifdef SPMD
USE mpicomm
#endif
!
IMPLICIT NONE
!
  integer :: i, j, ii, jj, ii1, jj1, k, ni, nj
  INTEGER, DIMENSION (8) :: T
  REAL :: rr,p,x,RANF 
  REAL :: xx(nx), yy(ny)
  INTEGER :: SEED      
!
  real, dimension(nx,ny) :: pos
  real, dimension(ip_start:ip_end,jp_start:jp_end) :: pos_l
!
  real, parameter :: d_pt = 2.5, xr = 1500., yr = 1500.
!
!----------------------------------------------------------!
!               Add pseudo-random perturbations            !
!----------------------------------------------------------!
!
  select case (trim(casename))
!
    case ('MIDLAT')
!
!  Initialisations
!
      if (mypid == 0) then
        CALL DATE_AND_TIME(VALUES = T)
        SEED = 43*mypid+7*T(1)+70*(T(2)+12*(T(3)+31*(T(5)+23*(T(6)+59*T(7)))))
        IF (MOD(SEED,2).EQ.0) SEED = SEED-1
!
        do i = 1, nx
          xx(i) = real(i-1)*dx
        enddo
        do j = 1, ny
          yy(j) = real(j-1)*dy
        enddo
        ni = int(xr/dx)+1
        nj = int(yr/dy)+1
!
!  Define randomly distributed bubbles in plane
!
        pos = 0.
        p = 800. / real((nx-5)*(ny-5))
        do j = 4, ny-2
          do i = 4, nx-2
            x = RANF(SEED)
            if (x < p) then
              do jj = j-nj, j+nj
                jj1 = jj
                if (jj1 > ny-2) jj1 = jj1 - (ny-5)
                if (jj1 < 4) jj1 = jj1 + (ny-5)
                do ii = i-ni, i+ni
                  ii1 = ii
                  if (ii1 > nx-2) ii1 = ii1 - (nx-5)
                  if (ii1 < 4) ii1 = ii1 + (nx-5)
 	          rr = sqrt( ((xx(ii1) - xx(i))/xr)**2. + ((yy(jj1) - yy(j))/yr)**2. )
	          if (rr < 1.) pos(ii1,jj1) = pos(ii1,jj1) + 0.5*(1. + cos(pi*rr))
                enddo
              enddo
            endif
          enddo
        enddo
      endif
!
#ifdef SPMD
      call distribute(pos_l, pos)
#else
      pos_l = pos
#endif
!
!  Perturb potential temperature field
!
      do k = 11, 27	!between 300m and 800m
        do j = jp_start, jp_end
          do i = ip_start, ip_end
            state2%es(i,j,k) = state2%es(i,j,k) + d_pt*pos_l(i,j) 
          enddo
        enddo
      enddo
!
  end select
!
return 
end

      !



