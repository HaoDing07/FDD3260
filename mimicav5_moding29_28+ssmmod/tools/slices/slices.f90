
!      =================================================                            
  program slices
!      =================================================                            

  USE sharedmod
  USE netcdfmod
  USE diagnostics

  IMPLICIT NONE

  integer :: unit1, unit2, n, i, j, k, t
  real :: time, time1, time2, tstep, dt
!
  logical :: ex, scale_t=.false.
  character(len = 7) :: car71, car72
  character(len = 100) :: file_output, file_input, time_series, filename
!
  real, dimension(:,:,:,:), allocatable :: slice
!
  namelist /pro/file_output,file_input,time1,time2,tstep,scale_t,dx
!
!----------------------------------------------------------!
!                 First, read input file                   !
!----------------------------------------------------------!
!
!  Option namelist
!
  filename = 'namelist'
  open(3,file=filename(1:len_trim(filename)),status='old')
  read(3,nml=pro)
  close(3)
!
  call read_dim (file_input)
!
  if ( scale_t ) tt = tt*86400.
!
  ntime = 1
  if (time2 > time1) then
    ntime = int((time2-time1) / tstep) + 1
  endif
!
  allocate( slice(1:NXO,1:NYO,1:nvarslic,1:ntime) )
!
!----------------------------------------------------------!
!                   Post-process slices                    !
!----------------------------------------------------------!
! 
  k = 1
  time = time1
  dt = 60.
  slice = 0.
!
  do t = 1, NTO
    if ( tt(t) > time-dt .and. tt(t) < time+dt ) then
!
    do j = 1, NYO
      do i = 1, NXO
        if ( state(i,j,1,t)%pt > 1. ) then 
!
!  slices
!
          slice(i,j,1,k) = wind(i,j,1,t)%u
!
          slice(i,j,2,k) = wind(i,j,1,t)%v
!
          slice(i,j,3,k) = wind(i,j,1,t)%w
!
          slice(i,j,4,k) = state(i,j,1,t)%pt
! 
          slice(i,j,5,k) = state(i,j,1,t)%ptv
!
          slice(i,j,6,k) = state(i,j,1,t)%qv
!
          slice(i,j,7,k) = state(i,j,1,t)%qc
!
          slice(i,j,8,k) = state(i,j,1,t)%qr
!
          slice(i,j,9,k) = state(i,j,1,t)%rh
!  
          slice(i,j,10,k) = state(i,j,1,t)%buoy
!
          slice(i,j,11,k) = state(i,j,1,t)%beff
!
          slice(i,j,12,k) = state(i,j,1,t)%mse
!
!          if ( state(i,j,1,t)%qc + state(i,j,1,t)%qr > 1.e-6 ) slice(i,j,13,k) = 1.
          if ( state(i,j,1,t)%qc + state(i,j,1,t)%qr > 1.e-6 .and. wind(i,j,1,t)%w > 0.1 ) slice(i,j,13,k) = 1.
!
          slice(i,j,14,k) = state(i,j,1,t)%wvp
!
          slice(i,j,15,k) = state(i,j,1,t)%cth
!
          slice(i,j,16,k) = state(i,j,1,t)%ptvar
!
          slice(i,j,17,k) = state(i,j,1,t)%qvvar
!
          slice(i,j,18,k) = tend(i,j,1,t)%w(1) + tend(i,j,1,t)%w(2) + tend(i,j,1,t)%w(3)
!
          slice(i,j,19,k) = tend(i,j,1,t)%w(3)
!
          slice(i,j,20,k) = tend(i,j,1,t)%w(6) + tend(i,j,1,t)%w(7) + tend(i,j,1,t)%w(8)
!
!         Elements 21 and 22 will contain cloud labels and cloud sizes resp.
!
          slice(i,j,23,k) = tend(i,j,1,t)%w(9)
!
          slice(i,j,24,k) = tend(i,j,1,t)%qt(1) + tend(i,j,1,t)%qt(2) + tend(i,j,1,t)%qt(3)
!
          slice(i,j,25,k) = tend(i,j,1,t)%qt(3)
!
          slice(i,j,26,k) = tend(i,j,1,t)%qt(9)
!
          slice(i,j,27,k) = tend(i,j,1,t)%pt(1) + tend(i,j,1,t)%pt(2) + tend(i,j,1,t)%pt(3)
!
          slice(i,j,28,k) = tend(i,j,1,t)%pt(3)
!
          slice(i,j,29,k) = tend(i,j,1,t)%pt(4) !* (-1004./2501000.)
!
          slice(i,j,30,k) = tend(i,j,1,t)%pt(9)
!
        endif
      enddo
    enddo
!
!  Renormalization
!
    slice(:,:,10,k) = slice(:,:,10,k) - sum( slice(:,:,10,k) ) / real(NXO*NYO)
!
    if ( time2 > time1 .and. time < time2 ) then
      k = k + 1
      time = time + tstep
      goto 101
    else
      exit
    endif
!
    endif
101 continue 
  enddo
!
!----------------------------------------------------------!
!		    Calculate diagnostics                  !
!----------------------------------------------------------!
!
!  Cloud labeling
!
  call cloud_ids ( 0.5, slice(:,:,:,1) )
!
  do k = 2, ntime
    call cloud_ids ( 0.5, slice(:,:,:,k) )
  enddo
!
!  Frequency distributions
!
!  call distributions ( 3, 0.5, slice(:,:,:,1), file_output )	! Vertical velocity
!  call distributions ( 5, 0.5, slice(:,:,:,1), file_output )	! Virtual potential temperature
!  call distributions ( 6, 0.5, slice(:,:,:,1), file_output )	! Water vapour
!
!  Empty space
!
!  call empty_space ( NXO*NYO, 0.5, slice(:,:,:,1), file_output, .true. )
!
!  time correlations
!
!  call time_correlations ( 11, 14, 22, 0.5, slice, tstep, file_output, .true. )
!
!  call time_corr_2D ( 21, 4, 6, 22, 0.5, slice, tstep, file_output, .true. )
!
!  space correlations
!
!  call space_correlations ( 60, 5, 5, 0.5, slice(:,:,:,1), file_output, .true. )
!
!  Radial distributions
!
  call distance_stats ( 50, 27, 28, 29, 30, 0.5, slice(:,:,:,1), file_output )
!
!----------------------------------------------------------!
!		     Deallocate and end                    !
!----------------------------------------------------------!
!
  deallocate( xx, yy, zz, tt )
  deallocate( wind, state, tend )
!
!----------------------------------------------------------!
!  
  end program slices
