
!      =================================================                            
  program profiles
!      =================================================                            

  USE sharedmod
  USE netcdfmod

  IMPLICIT NONE

  integer :: unit1, unit2, n, k, t
  real :: tmin, tmax
!
  logical :: ex, scale_t=.false.
  character(len = 7) :: car71, car72
  character(len = 100) :: file_output, file_input, time_series, filename
!
  integer :: k1, k2, k3, k4, k5
  integer, dimension(:), allocatable :: count
  real, dimension(:,:), allocatable :: prof, ts_tab
!
  namelist /pro/file_output,file_input,time_series,tmin,tmax,scale_t
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
  if ( scale_t ) time = time*86400.
!
  allocate( prof(1:NXO,1:nvarin), ts_tab(1:NTO,nts), count(1:NXO) )
!
!----------------------------------------------------------!
!                   Post-process profiles                  !
!----------------------------------------------------------!
!
  k1 = 0; k2 = 0; k3 = 0; k4 = 0; k5 = 0
  do k = 1, NXO-1
    if ( xx(k) >= 750. .and. k1 == 0 ) k1 = k 
    if ( xx(k) >= 1500. .and. k2 == 0 ) k2 = k 
    if ( xx(k) >= 4000. .and. k3 == 0 ) k3 = k 
    if ( xx(k) >= 7500. .and. k4 == 0 ) k4 = k 
    if ( xx(k) >= 11000. .and. k5 == 0 ) k5 = k 
  enddo
!
  count = 0
  prof = 0.
  ts_tab = 0.
  do t = 1, NTO
    if ( time(t) > tmin .and. time(t) < tmax ) then
!
    do k = 1, NXO
      if ( state(k,1,1,t)%pt > 1. ) then 
        count(k) = count(k) + 1
!
!  Profiles
!
        prof(k,1) = prof(k,1) + wind(k,1,1,t)%u
!
        prof(k,2) = prof(k,2) + wind(k,1,1,t)%v
!
        prof(k,3) = prof(k,3) + wind(k,1,1,t)%w
!
        prof(k,4) = prof(k,4) + state(k,1,1,t)%pt
! 
        prof(k,5) = prof(k,5) + state(k,1,1,t)%mse
!
        prof(k,6) = prof(k,6) + state(k,1,1,t)%qv + state(k,1,1,t)%qc + state(k,1,1,t)%qr
!
        prof(k,7) = prof(k,7) + state(k,1,1,t)%qc
!
        prof(k,8) = prof(k,8) + state(k,1,1,t)%qr
!
        prof(k,9) = prof(k,9) + state(k,1,1,t)%rh
!  
        prof(k,10) = prof(k,10) + state(k,1,1,t)%buoy
!
        prof(k,11) = prof(k,11) + state(k,1,1,t)%beff
!
        prof(k,12) = prof(k,12) + state(k,1,1,t)%evap
!
        prof(k,13) = prof(k,13) + state(k,1,1,t)%cent
!
        prof(k,14) = prof(k,14) + state(k,1,1,t)%cdet
!
        prof(k,15) = prof(k,15) + state(k,1,1,t)%sent
!
        prof(k,16) = prof(k,16) + state(k,1,1,t)%sdet
!
      endif
    enddo
!
!  Time series
!
!        ts_tab(t,1) = sum(wind(k1:k2,1,1,t)%w) / real(k2-k1+1)
!
!        ts_tab(t,2) = sum(wind(k2:k3,1,1,t)%w) / real(k3-k2+1)
!
        ts_tab(t,1) = sum(state(k1:k2,1,1,t)%buoy) / real(k2-k1+1)
!
        ts_tab(t,2) = sum(state(k2:k3,1,1,t)%buoy) / real(k3-k2+1)
!
        ts_tab(t,3) = sum(state(k3:k4,1,1,t)%buoy) / real(k4-k3+1)
!
        ts_tab(t,4) = sum(state(k1:k2,1,1,t)%beff) / real(k2-k1+1)
!
        ts_tab(t,5) = sum(state(k2:k3,1,1,t)%beff) / real(k3-k2+1)
!
        ts_tab(t,6) = sum(state(k3:k4,1,1,t)%beff) / real(k4-k3+1)
!
    endif
  enddo
!
  do n = 1, nvarin
    prof(:,n) = prof(:,n) / real(count)
  enddo
!
!----------------------------------------------------------!
!		   Diagnostics and outputs                 !
!----------------------------------------------------------!
!
  write(car71,'(i7)') int(tmin)
  write(car72,'(i7)') int(tmax)
!
!  Open file
!
  unit1=203
  unit2=204
!
  inquire (FILE=trim(file_output), OPENED=ex)
  if (.not.ex) then
    open(unit1,file=file_output(1:len_trim(file_output)),	 &
           form='formatted', status='unknown')
  endif
!
  inquire (FILE=trim(time_series), OPENED=ex)
  if (.not.ex) then
    open(unit2,file=time_series(1:len_trim(time_series)),	&
    	   form='formatted', status='unknown')
  endif
!
!  Write in file
!
  do k = 1, NXO
    write(unit1,100) xx(k), prof(k,:)
  enddo
!
  do t = 1, NTO
    if ( time(t) > tmin .and. time(t) < tmax ) write(unit2,*) time(t), ts_tab(t,:)
  enddo
!
  close (unit1)
  close (unit2)
!
100 format (18(f13.6,1x))
!
!----------------------------------------------------------!
!		     Deallocate and end                    !
!----------------------------------------------------------!
!
  deallocate( xx, yy, zz, time )
  deallocate( wind, state )
!
!----------------------------------------------------------!
!  
  end program profiles
