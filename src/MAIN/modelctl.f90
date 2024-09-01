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
!  MODELCTL.F:                   
!
!  Purpose:
!       A center piece of the model that actually processes the
!              model setup and integration, and produces output.                       
!
!  Author
!       Chien Wang
!       MIT Joint Program on Science and Policy of Global Change
!
! ================================================================

!       ===================
       subroutine modelctl
!       ===================

#ifdef SPMD
       USE mpi
#endif
!
!  Data types modules
!
       USE shared_all
       USE typedef_svalue
!
!  Physics modules
!
       USE time_step
       USE column_step
       USE piggyback
       USE radiationmod
       USE nesting
       USE restart
       USE diagnostics
       USE micropack
       USE fft_tools
       USE check
       USE lagrange
       USE initialize

       IMPLICIT NONE

       integer            :: h, k, m, nod
       integer            :: my_val(8)
#ifdef SPMD
       integer            :: ierr, status(MPI_STATUS_SIZE) 
#endif
!
       real(4)            :: tarry(2)
       !real, allocatable  :: tcpu(:)
!
       logical            :: per_radstep, per_creact, per_ires, per_inest, per_pert 
!
       character (len = 100) :: car
       character (len = 8) :: my_date
       character (len =10) :: my_time
       character (len = 5) :: my_zone
       character (len = 3) :: my_month(12) =                             &
                                (/'JAN','FEB','MAR','APR','MAY','JUN',   &
                                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
!
       save
!
!----------------------------------------------------------!
!                  Synchronize and start                   !
!----------------------------------------------------------!
!
#if ( defined SPMD )
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif
!
  CALL DTIME_STD(tarry)
!
!  Open log file
!
  if ( mypid == 0 ) then
    open(7,file=prtname(1:len_trim(prtname)), access='SEQUENTIAL', position='append', status='unknown')
!
    call date_and_time (my_date, my_time, my_zone, my_val)
    car = 'Start Time: '//my_month(my_val(2))//' '//my_date(7:8)//', '//my_date(1:4)//'   '      &
         //my_time(1:2)//':'//my_time(3:4)//':'//my_time(5:10)
    write(7,*) car
  endif
!
!----------------------------------------------------------!
!        Splitting, allocations and initialisations        !
!----------------------------------------------------------!
!
!  Setting flags
!
  if (lmicro == 0) then
    nhydro = 0
  else if (lmicro == 1) then
    nhydro = 2
  else if (lmicro == 2) then
    nhydro = 3
  else if (lmicro == 3) then
    nhydro = 5
  else if (lmicro == 4) then
    nhydro = 6
  endif
!
#ifdef CHEM_REACT
  per_creact = .true.
#else
  per_creact = .false.
#endif
!
!  Initialise dims
!
#ifndef COLUMN
  ip_start  = 1
  ip_end    = nx

  it_start  = ip_start+3
  it_end    = ip_end-2

  jp_start  = 1
  jp_end    = ny

#ifdef MODEL_3D
  jt_start  = jp_start+3
  jt_end    = jp_end-2
#else
  jt_start = 1
  jt_end = jp_end
#endif
!
#ifdef SPMD
#ifdef DECOMP_2D
  call decomp_2D
#else
  call decomp_1D
#endif
#endif
!
!  Broadcast namelist data
!
       call isendcom
!
#ifdef SPMD
       call create_mpi_types
#endif
!
!  For pressure solver
!
#ifdef ANELASTIC
       call get_wv
#endif
!
#else
      ip_start = 1
      it_start = 1
      ip_end = 1
      it_end = 1
!
      jp_start = 1
      jt_start = 1
      jp_end = 1
      jt_end = 1
#endif
!
#ifdef AERO_ENABLE
      nmode = nmode0
      if (reg_mode) nmode = nmode + 1
#else
      nmode = 0
#endif
!
!  Count outputs
!
      call count_outputs
!
!  Allocate work arrays
!
      call allocations
!
!----------------------------------------------------------!
!                       Initialize                         !
!----------------------------------------------------------!
!
#if ( defined SPMD )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
!  Set following flags to 1 if you dont want first passage
!
      n_io = 0
      n_ts = 0
      n_pr = 0
      n_sl = 0
      n_av = 1
      n_lag = 0
      n_res = 1
      n_rad = 1
      n_pe = 1
      n_ne = 0
!
      nite=0
      ntime=0
      time=0.
!
!  Initialise microphysics
!
      call init_micro
!
!  setup up aerosol flags
!
#ifdef AERO_ENABLE
      call setup_aero_flags
#endif
!
!  Other initialisations
!
      call init_other
!
!  Main nitialization
!
      if (new_run) then 
        call setrun
      else
        call readdata
      endif
!
!  First nest
!
#ifdef NESTING
      if ( .not.nest_run ) call savenest
#endif
!
      if (mypid == 0) then
        flush(7)
        !allocate(tcpu(1))
        !tcpu(1) = 0.
      endif
!
!  First outputs
!
      if (new_run) then
        if (with_piggy) call output_pig
!
        call output_main
      endif
!
!-----------------------------------------------------------!
!                                                           !
!                    Starting time loop                     !
!                                                           !
!-----------------------------------------------------------!
!
      do
!
        if (verbose > 0) call write_debug('Starting time loop')
!
        nite  = nite + 1
        ntime = ntime + 1
        time  = time + dt0
        nod = int(2.*(real(ntime-1)/2. - real(ntime-1)/2.))
!
!  Logicals
!
        per_ires = (time-tstart >= real(n_res)*ires)
        per_radstep = (time-tstart >= real(n_rad)*abs(iradx) .or. nite == 1)
        per_inest = (time-tstart >= real(n_ne)*inest) .or. (time.eq.tstart)
        per_pert = (time-tstart >= real(n_pe)*tpert)
!
        if (time-tstart >= real(n_res)*ires) n_res=n_res+1 
        if (time-tstart >= real(n_rad)*abs(iradx)) n_rad=n_rad+1 
        if (time-tstart >= real(n_ne)*inest) n_ne=n_ne+1 
        if (time-tstart >= real(n_pe)*tpert) n_pe=n_pe+1000000 
!
!----------------------------------------------------------!
!             Update time varying variables                !
!----------------------------------------------------------!
!
!  Lateral boundaries, nest
!
#ifdef NESTING
        if (nest_run) call readnest
#endif
!
!  Surface fluxes and large scale forcings
!
        if (with_tvar) call timevar
!
!  Global perturbations
!
        if (per_pert) call add_pert
        
!  Nudging warm advection
        if (with_warmadv) call nudge_warm_adv
        if (with_warmadv) call nudge_q_adv
!
!  Piggy-backing
!
        if (with_piggy) call piggy (per_radstep)


!
!  Radiation
!
! #ifdef RAD_ENABLE
        if(with_rad .and. with_radmod .and. time==dt0) call init_rad(nz, p0)
        if(with_rad .and. with_radmod .and. time==t_radmod) call init_rad(nz, p0)
! #endif 
!----------------------------------------------------------!
!                      time stepping                       !
!----------------------------------------------------------!
!
#ifdef COLUMN
	call columnstep (per_radstep)
#else
        call stepping (per_radstep) 
#endif
!
!----------------------------------------------------------!
!          Calculate and integrate chemicals               !
!----------------------------------------------------------!
!
#ifdef CHEM_ENABLE
        call chem_m_front( ip_start, ip_end, jp_start, jp_end,	&
                           nod, per_creact )
#endif
!
!----------------------------------------------------------!
!		Track Lagrangian Parcels		   !
!----------------------------------------------------------!
!
#ifdef LAGRANGE
        call lagrange_step
#endif
!
!----------------------------------------------------------!
!                        Output                            !
!----------------------------------------------------------!
!
        call output_main
!
!----------------------------------------------------------!
!                 Check scalars & winds                    !
!----------------------------------------------------------!
!
#ifndef COLUMN
        call checkall ( wind2, state2, hydromtr2 )
#endif
!
!----------------------------------------------------------!
!                     Restart data                         !
!----------------------------------------------------------!
!
!  Global restart file
!
        if ( per_ires ) call savedata
!
!  Write nest restart file
!
#ifdef NESTING
        if ( .not.nest_run .and. per_inest ) call savenest
#endif
!
!  Print summary and terminate time loop
!
        if ( mypid == 0 ) then
          write(7,1004)ntime,time,dt0,dtau,maxcfl,maxcflc,maxdiv
1004      format(' ',5x,'  Step =',I6, '  Time=', f10.2, ' sec.',      &
                 ' time step=', e8.2, ' acoustic tstep=', e8.2,      &
                 ' Max. CFL=', f6.3, ' Max. CFL_acoustic=', f6.3,    &
                 ' Max. Div=', e9.2)
!
          CALL DTIME_STD(tarry)
          !call store_time(ntime,tarry(1),tcpu)
          write(7,100) tarry(1)
          if (verbose > 0) call systMemUsage()
          write(7,*)
          flush(7)
        end if
!
!  Ending time loop
!
        if ( time >= tstop) exit
!
      enddo
!
!----------------------------------------------------------!
!	                  Finishing                        !
!----------------------------------------------------------!
!
!  Get end time before finalising
!
    if ( mypid == 0 ) then
      call date_and_time (my_date, my_time, my_zone, my_val)
      write(7,*)' '
      write(7,*)'Ending Time: ',                                    &
                   my_month(my_val(2)),' ',                         &
                   my_date(7:8),', ',my_date(1:4),'	',          &
                   my_time(1:2),':',my_time(3:4),':',my_time(5:10)
      flush(7)
      !print*, tcpu
      !deallocate(tcpu)
    endif
!
!  Save final state for restart
!
    call savedata 
      
!
! === Stopping MIMICA cleanly
!
    call stop_mimica ( 0 ) 
!  
100 format(11x, "Elapsed time = ",f10.3," sec ")
!
!-----------------------------------------------------------!
!
return

contains

!	====================
Subroutine store_time(nt,ts,tcpu)
!	====================

   IMPLICIT NONE
   integer :: nt
   real(4) :: ts
   real, allocatable :: tcpu(:), tcpu2(:)

   allocate(tcpu2(1:nt))
   tcpu2=tcpu
   deallocate(tcpu)
   allocate(tcpu(1:nt+1))
   tcpu(1:nt)=tcpu2
   tcpu(nt+1)=ts
   deallocate(tcpu2)
  return
end Subroutine

end
