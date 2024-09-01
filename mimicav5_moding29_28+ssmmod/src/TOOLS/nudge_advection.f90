#include "ctrparam.h"

module nudge_adv

 !General modules

    USE gridno 
    USE shared_data    
    USE shared_all
    USE shared_state
    USE averages
    USE gradients
  
    IMPLICIT NONE

    public :: nudge_warm_adv
    public :: nudge_q_adv
    
    contains

   !! Nudging warm advection (po valley case as a example)
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
      ! write(7,*) 'thermo%T(1,1,90) =', thermo%T(1,1,90)
      do k = int(H_warm_adv1(1))/5,  int(H_warm_adv1(size(H_warm_adv1)))/5
        do j = jp_start, jp_end
          do i = ip_start, ip_end
            state2%es(i,j,k) = state2%es(i,j,k)+ dT_warm_adv1((int(k*5-H_warm_adv1(1)))/5+1)
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
      do k = int(H_warm_adv2(1))/5,  int(H_warm_adv2(size(H_warm_adv2)))/5
        do j = jp_start, jp_end
          do i = ip_start, ip_end
            state2%es(i,j,k) = state2%es(i,j,k)+ dT_warm_adv2((int(k*5-H_warm_adv2(1)))/5+1)
          enddo
        enddo
      enddo
      write(7,*) 'nudging warm advection 2 is working at this timestep'
    endif 

    if (with_warmadv .and. time == tstop2_warmadv) write(7,*) 'Terminating nudging warm advection'
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
          do k = int(H_q_adv1(1))/5,  int(H_q_adv1(size(H_q_adv1)))/5
            do j = jp_start, jp_end
              do i = ip_start, ip_end
                state2%qt(i,j,k) = state2%qt(i,j,k)+ dQ_q_adv1((int(k*5-H_q_adv1(1)))/5+1)
              enddo
            enddo
          enddo
          write(7,*) 'nudging q advection 1 is working at this timestep'
        
        !!! Period 2
        endif 
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
          do k = int(H_q_adv2(1))/5,  int(H_q_adv2(size(H_q_adv2)))/5
            do j = jp_start, jp_end
              do i = ip_start, ip_end
                state2%qt(i,j,k) = state2%qt(i,j,k)+ dQ_q_adv2((int(k*5-H_q_adv2(1)))/5+1)
              enddo
            enddo
          enddo
          write(7,*) 'nudging q advection 2 is working at this timestep'
          
        endif 

        if (with_warmadv .and. time == tstop2_warmadv) write(7,*) 'Terminating nudging q advection'
        return
        end subroutine nudge_q_adv
      !
      

      end module nudge_adv