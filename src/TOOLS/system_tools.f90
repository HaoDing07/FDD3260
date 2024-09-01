#include "ctrparam.h"

!	====================
	subroutine stop_mimica (stopflag) 
!	====================

#ifdef SPMD
USE mpi
#endif
USE gridno
USE shared_data
USE, intrinsic :: iso_c_binding 

#ifdef LAGRANGE
USE shared_lagrange
#endif
#ifdef NESTING
USE nesting
#endif

  IMPLICIT NONE

  integer :: ierr, stopflag
  logical :: op

!
! === Synchronize MPI before terminating
!
#ifdef SPMD
  call MPI_BARRIER ( MPI_COMM_WORLD, ierr )
#endif

!
! === Write info
!
  if ( mypid .eq. 0 ) then
    write(7,*) ''
    if (stopflag == 0) then
      write(7,*) 'MIMICA ENDED NORMALLY'  
    else
      write(7,*) 'AN ERROR WAS DETECTED: ABORTING MIMICA'  
    endif
!
    close(7)
!
    inquire(UNIT=108, OPENED=op)
    if (op) close (108)
  endif

!
! === Deallocate general
!
  call deallocations
  
  call dealloc_micro
  
#ifdef RAD_ENABLE
  call dealloc_rad
#endif

!
! === Deallocate specific
!
#ifdef LAGRANGE
  if (allocated(parcel_pos)) deallocate (parcel_pos)
  if (allocated(parcel_loc)) deallocate (parcel_loc)
  if (allocated(parcel_sca)) deallocate (parcel_sca) 
  if (allocated(aer)) deallocate (aer)
#endif
!  
#ifdef NESTING
  if (nest_run) call dealloc_nest
#endif
!
#ifdef SPMD
  call free_mpi_types 
#endif
!
#if (defined ANELASTIC) && (!defined COLUMN)
  deallocate( wv )
#endif
!
! ===  Finalize MPI:
!      

#ifdef SPMD
  call MPI_FINALIZE ( ierr )
#endif

!
! === Abort
!  
  if (stopflag /= 0) stop  
!
return
end

! ==========================================
  subroutine dealloc_micro
! ==========================================
!
  USE shared_data
  USE shared_hydro

  IMPLICIT NONE 
!
  integer :: h
!
  deallocate ( pxx, qmin, lmax, xnmin )
  deallocate ( hydrop, hydroc ) 
!  
!  deallocate ( mk, mk_tmp, rwork, iwork )
!  
  return
  end subroutine

!  ===============================
  subroutine dealloc_rad
! ================================

  USE radia

#ifdef RAD_ENABLE
  if (allocated(pp)) deallocate(pp)
  if (allocated(pt)) deallocate(pt)
  if (allocated(ph)) deallocate(ph)
  if (allocated(po)) deallocate(po)
  if (allocated(plwc)) deallocate(plwc)
  if (allocated(piwc)) deallocate(piwc)
  if (allocated(pre)) deallocate(pre)
  if (allocated(pde)) deallocate(pde)
  if (allocated(prwc)) deallocate(prwc)
  if (allocated(pgwc)) deallocate(pgwc)
  if (allocated(fds)) deallocate(fds)
  if (allocated(fus)) deallocate(fus)
  if (allocated(fdir)) deallocate(fdir)
  if (allocated(fuir)) deallocate(fuir)
  if (allocated(taus)) deallocate(taus)
  if (allocated(tauir)) deallocate(tauir)
#endif
!
  return
  end subroutine dealloc_rad

!	====================
subroutine systMemUsage()
!	====================

    USE ifport
    IMPLICIT NONE
    character(len=200):: filename=' '
    character(len=80) :: line
    character(len=8)  :: pid_char=' '
    integer ::  pid, ios, fu, valueRSS

    valueRSS = 0
    pid=getpid()
    write(pid_char,'(I8)') pid
    filename='/proc/'//trim(adjustl(pid_char))//'/status'
    open(newunit=fu, file=trim(filename), action='read')
    do
        read(fu,'(A)',iostat=ios) line
        if(ios/=0) exit
        if(line(1:6) == 'VmRSS:') then
            read (line(7:),*) valueRSS
            write(7,101) valueRSS/1000.
            exit
        endif
    enddo
    close(fu)

101 format(11x, "Memory usage = ",f10.3," MB ")
return
end subroutine

!	====================
Subroutine dtime_std(time)
!	====================

  IMPLICIT NONE
  Real(4) time(2)
  Real,Save :: last_time = 0
  Real this_time
  Intrinsic Cpu_Time
  Call Cpu_Time(this_time)
  time(1) = this_time - last_time
  time(2) = 0
  last_time = this_time
  return
End Subroutine

!	====================
subroutine write_debug(fname) 
!	====================

  USE shared_data
  IMPLICIT NONE
  character(*) :: fname

  if ( mypid == 0 ) then
      write(7,*) trim(fname) 
      flush(7)
  endif

  return
end subroutine write_debug

subroutine debug_minmax ( data, name )

#ifdef SPMD
  USE mpi
#endif
  USE gridno
  USE shared_data
        
  IMPLICIT NONE

  real, dimension(ip_start:ip_end,jp_start:jp_end,1:nz) :: data
  real :: maxl, minl, maxg, ming
  character(len=*) :: name
  integer :: ierr
!
!----------------------------------------------------------!
!                 Find local min/max values                !
!----------------------------------------------------------!
!
  maxl = maxval(data(it_start:it_end,jt_start:jt_end,1:nz))
  minl = minval(data(it_start:it_end,jt_start:jt_end,1:nz))
!
!  Find global min/max values
!
#if ( defined SPMD )
  CALL MPI_ALLREDUCE (maxl, maxg, 1, REALTYPE, MPI_MAX, MPI_COMM_WORLD, ierr) 
  CALL MPI_ALLREDUCE (minl, ming, 1, REALTYPE, MPI_MIN, MPI_COMM_WORLD, ierr) 
#else
  maxg = maxl
  ming = minl
#endif
!
!  Print out result
!
  if (  mypid == 0 ) then
    write(7,*)  '    VARIABLE ', trim(name), ', MAX/MIN VALUES: ', maxg, ming
    flush(7)
  endif
!
  return
end subroutine debug_minmax

