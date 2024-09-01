
!      =================================================                            
  program slices
!      =================================================                            

  USE sharedmod
  USE netcdfmod
  USE diagnostics

  IMPLICIT NONE

  integer :: nslic, ntim, unit1, unit2, n, i, j, k, t, t1
  real :: time, time1, tstop, dt, dtim, dslic, lslic
  real :: exn
!
  logical :: ex, scale_t=.false.
  character(len = 3) :: car3
  character(len = 4) :: car4
  character(len = 5) :: car5
  character(len = 7) :: car71, car72
  character(len = 100) :: file_output, file_name, file_input, time_series, filename
!
  real, dimension(:,:,:,:,:), allocatable :: slice
!
  namelist /pro/file_output,file_input,time1,tstop,scale_t,dx,nslic,ntim,dslic,dtim,norm
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
  ntime = 1
  time = time1
  dt = 2.
!
!----------------------------------------------------------!
!                   Post-process slices                    !
!----------------------------------------------------------!
!
  lslic = 0.
  do n = 1, nslic
!
    if ( n == nslic/2+1 ) lslic = 0.
    lslic = lslic + dslic
!
    if ( n < nslic/2+1 ) then
      file_name = trim(file_input)//'_x_'
    else
      file_name = trim(file_input)//'_y_'
    endif
!
    if ( lslic < 1000. ) then
      write(car3,'(i3)')int(lslic)
      call read_dim ( trim(file_name)//car3//'.nc' )
    else if ( lslic < 10000. ) then
      write(car4,'(i4)')int(lslic)
      call read_dim ( trim(file_name)//car4//'.nc' )
    else if ( lslic < 100000. ) then
      write(car5,'(i5)')int(lslic)
      call read_dim ( trim(file_name)//car5//'.nc' )
    endif
!
    if (.not.allocated(slice)) then
      allocate( slice(1:NYO,1:NZO,1:nvarslic,1:nslic,1:ntim) )
      slice = 0.
    endif
!
    if ( scale_t ) tt = tt*86400.
!
    t1 = 1
    do t = 1, NTO
!
    if ( (tt(t) > time+real(t1-1)*dtim-dt .and. tt(t) < time+real(t1-1)*dtim+dt) ) then
!
    do i = 1, NYO
      do j = 1, NZO
        if ( state(1,i,j,t)%pt > 1. ) then 
!
!  slices
!
          slice(i,j,1,n,t1) = wind(1,i,j,t)%u
!
          slice(i,j,2,n,t1) = wind(1,i,j,t)%v
!
          slice(i,j,3,n,t1) = wind(1,i,j,t)%w
!
          slice(i,j,4,n,t1) = state(1,i,j,t)%pt
! 
          slice(i,j,5,n,t1) = state(1,i,j,t)%ptv
!
          slice(i,j,6,n,t1) = state(1,i,j,t)%qv
!
          slice(i,j,7,n,t1) = state(1,i,j,t)%qc
!
          slice(i,j,8,n,t1) = state(1,i,j,t)%qr + state(1,i,j,t)%qc
!
          slice(i,j,9,n,t1) = state(1,i,j,t)%rh
!  
          slice(i,j,10,n,t1) = state(1,i,j,t)%buoy
!
          slice(i,j,11,n,t1) = state(1,i,j,t)%beff
!
          slice(i,j,12,n,t1) = state(1,i,j,t)%mse
!
          if ( state(1,i,j,t)%vof > 0.1 ) slice(i,j,13,n,t1) = 1.
!          if ( state(1,i,j,t)%qc + state(1,i,j,t)%qr > 1.e-5 ) slice(i,j,13,n,t1) = 1.
!
          slice(i,j,14,n,t1) = state(1,i,j,t)%ent
!
          slice(i,j,15,n,t1) = state(1,i,j,t)%det
!
          slice(i,j,16,n,t1) = state(1,i,j,t)%ptvar
!
          slice(i,j,17,n,t1) = state(1,i,j,t)%qvvar
!
          slice(i,j,18,n,t1) = tend(1,i,j,t)%w(1) + tend(1,i,j,t)%w(2) + tend(1,i,j,t)%w(3)
!
          slice(i,j,19,n,t1) = tend(1,i,j,t)%w(6)
!
          slice(i,j,20,n,t1) = tend(1,i,j,t)%w(7) + tend(1,i,j,t)%w(8)
!
!         Elements 21, 22, 23 and 24 will contain cloud labels, cloud height, cloud base and max velocity
!
          slice(i,j,25,n,t1) = tend(1,i,j,t)%w(9)
!
          slice(i,j,26,n,t1) = tend(1,i,j,t)%qc(1) + tend(1,i,j,t)%qr(2) + tend(1,i,j,t)%qc(1) + tend(1,i,j,t)%qr(2)
!
          slice(i,j,27,n,t1) = tend(1,i,j,t)%qc(3) + tend(1,i,j,t)%qr(3)
!
          slice(i,j,28,n,t1) = tend(1,i,j,t)%qc(4) + tend(1,i,j,t)%qr(4)
!
          slice(i,j,29,n,t1) = tend(1,i,j,t)%qc(9) + tend(1,i,j,t)%qr(9)
!
          slice(i,j,30,n,t1) = tend(1,i,j,t)%ptv(1) + tend(1,i,j,t)%ptv(2) !+ tend(1,i,j,t)%ptv(3)
!
          slice(i,j,31,n,t1) = tend(1,i,j,t)%ptv(3)
!
          slice(i,j,32,n,t1) = tend(1,i,j,t)%ptv(4) !* (-1004./2501000.)
!
          slice(i,j,33,n,t1) = tend(1,i,j,t)%ptv(5)
!
          slice(i,j,34,n,t1) = tend(1,i,j,t)%ptv(9)
! 
          slice(i,j,35,n,t1) = tend(1,i,j,t)%pt(1) + tend(1,i,j,t)%pt(2) !+ tend(1,i,j,t)%pt(3)
!
          slice(i,j,36,n,t1) = tend(1,i,j,t)%pt(3)
!
          slice(i,j,37,n,t1) = tend(1,i,j,t)%pt(4) !* (-1004./2501000.)
!
          slice(i,j,38,n,t1) = tend(1,i,j,t)%pt(5)
!
          slice(i,j,39,n,t1) = tend(1,i,j,t)%pt(9)
!
          slice(i,j,40,n,t1) = state(1,i,j,t)%p
!
          slice(i,j,41,n,t1) = wind(1,i,j,t)%div
!
        endif
      enddo
    enddo
!
!  Renormalization
!
    do j = 1, NZO
      slice(:,j,10,n,t1) = 1000.*(slice(:,j,10,n,t1) - sum( slice(:,j,10,n,t1) ) / real(NYO))
!      slice(:,j,18,n,t1) = (slice(:,j,18,n,t1) - sum( slice(:,j,18,n,t1) ) / real(NYO))
!      slice(:,j,19,n,t1) = (slice(:,j,19,n,t1) - sum( slice(:,j,19,n,t1) ) / real(NYO))
!      slice(:,j,20,n,t1) = (slice(:,j,20,n,t1) - sum( slice(:,j,20,n,t1) ) / real(NYO))
!      slice(:,j,25,n,t1) = slice(:,j,25,n,t1)
    enddo
!
!  Diagnose contributions to buoyancy tendencies: first 4: PT contribution, last 4: hydrometeor loading
!
!    slice(:,:,30,n,t1) = (1. + 0.602*(slice(:,:,6,n,t1)+slice(:,:,8,n,t1)) - 1.602*slice(:,:,8,n,t1)) * slice(:,:,35,n,t1)
!    slice(:,:,31,n,t1) = (1. + 0.602*(slice(:,:,6,n,t1)+slice(:,:,8,n,t1)) - 1.602*slice(:,:,8,n,t1)) * slice(:,:,36,n,t1)
!    slice(:,:,32,n,t1) = (1. + 0.602*(slice(:,:,6,n,t1)+slice(:,:,8,n,t1)) - 1.602*slice(:,:,8,n,t1)) * slice(:,:,37,n,t1)
!    slice(:,:,34,n,t1) = (1. + 0.602*(slice(:,:,6,n,t1)+slice(:,:,8,n,t1)) - 1.602*slice(:,:,8,n,t1)) * slice(:,:,39,n,t1)
!
!    slice(:,:,26,n,t1) = - 1.602*slice(:,:,4,n,t1)*slice(:,:,26,n,t1)
!    slice(:,:,27,n,t1) = - 1.602*slice(:,:,4,n,t1)*slice(:,:,27,n,t1)
!    slice(:,:,28,n,t1) = - 1.602*slice(:,:,4,n,t1)*slice(:,:,28,n,t1)
!    slice(:,:,29,n,t1) = - 1.602*slice(:,:,4,n,t1)*slice(:,:,29,n,t1)
!
!  Diagnose relative humidity tendencies: first 4: log(RH) tendencies, last 4: qv tendencies
!
!    slice(:,:,26,n,t1) = slice(:,:,30,n,t1) - (1. + 0.602*(slice(:,:,6,n,t1)+slice(:,:,8,n,t1)) - 1.602*slice(:,:,8,n,t1))*slice(:,:,35,n,t1) + slice(:,:,4,n,t1)*slice(:,:,26,n,t1)
!    slice(:,:,27,n,t1) = slice(:,:,31,n,t1) - (1. + 0.602*(slice(:,:,6,n,t1)+slice(:,:,8,n,t1)) - 1.602*slice(:,:,8,n,t1))*slice(:,:,36,n,t1) + slice(:,:,4,n,t1)*slice(:,:,27,n,t1)
!    slice(:,:,28,n,t1) = slice(:,:,32,n,t1) - (1. + 0.602*(slice(:,:,6,n,t1)+slice(:,:,8,n,t1)) - 1.602*slice(:,:,8,n,t1))*slice(:,:,37,n,t1) + slice(:,:,4,n,t1)*slice(:,:,28,n,t1)
!    slice(:,:,29,n,t1) = slice(:,:,33,n,t1) - (1. + 0.602*(slice(:,:,6,n,t1)+slice(:,:,8,n,t1)) - 1.602*slice(:,:,8,n,t1))*slice(:,:,38,n,t1) 
!
!    slice(:,:,26,n,t1) = slice(:,:,26,n,t1) / (0.602*slice(:,:,4,n,t1))/slice(:,:,6,n,t1)
!    slice(:,:,27,n,t1) = slice(:,:,27,n,t1) / (0.602*slice(:,:,4,n,t1))/slice(:,:,6,n,t1)
!    slice(:,:,28,n,t1) = slice(:,:,28,n,t1) / (0.602*slice(:,:,4,n,t1))/slice(:,:,6,n,t1)
!    slice(:,:,29,n,t1) = slice(:,:,29,n,t1) / (0.602*slice(:,:,4,n,t1))/slice(:,:,6,n,t1)
!
!    slice(:,:,5,n,t1) = slice(:,:,4,n,t1) * (slice(:,:,40,n,t1)/100000.)**(288./1004.)
!    slice(:,:,14,n,t1) = min(max( (slice(:,:,5,n,t1) - 253.)/20., 0.), 1.)
!    slice(:,:,15,n,t1) = 4283.7744 / (slice(:,:,5,n,t1) - 30.03)**2.
!    slice(:,:,16,n,t1) = 5807.8125 / (slice(:,:,5,n,t1) -  7.65)**2.
!
!    slice(:,:,30,n,t1) = slice(:,:,26,n,t1) - slice(:,:,5,n,t1)/slice(:,:,4,n,t1) * (slice(:,:,14,n,t1)*slice(:,:,15,n,t1) + (1. - slice(:,:,14,n,t1))*slice(:,:,16,n,t1)) * slice(:,:,35,n,t1)
!    slice(:,:,31,n,t1) = slice(:,:,27,n,t1) - slice(:,:,5,n,t1)/slice(:,:,4,n,t1) * (slice(:,:,14,n,t1)*slice(:,:,15,n,t1) + (1. - slice(:,:,14,n,t1))*slice(:,:,16,n,t1)) * slice(:,:,36,n,t1)
!    slice(:,:,32,n,t1) = slice(:,:,28,n,t1) - slice(:,:,5,n,t1)/slice(:,:,4,n,t1) * (slice(:,:,14,n,t1)*slice(:,:,15,n,t1) + (1. - slice(:,:,14,n,t1))*slice(:,:,16,n,t1)) * slice(:,:,37,n,t1)
!    slice(:,:,34,n,t1) = slice(:,:,29,n,t1) - slice(:,:,5,n,t1)/slice(:,:,4,n,t1) * (slice(:,:,14,n,t1)*slice(:,:,15,n,t1) + (1. - slice(:,:,14,n,t1))*slice(:,:,16,n,t1)) * slice(:,:,38,n,t1)
!
!  Cloud labeling
!
    call cloud_ids ( 0.5, zz, maxval(slice(:,:,21,:,:)), slice(:,:,:,n,t1) )
!
    if ( t1 < ntim ) then
      t1 = t1 + 1
      if ( (tt(t) > tstop-dt .and. tt(t) < tstop+dt) ) exit
    else
      exit
    endif
!
    endif
!   
    enddo
!
  enddo
!
!----------------------------------------------------------!
!		    Calculate diagnostics                  !
!----------------------------------------------------------!
!
!  Frequency distributions
!
!  call distributions ( 0.5, nslic, ntim, slice )	! Vertical velocity
!
!  Radial distributions
!
  call distance_stats ( 101, zz, (/3, 4, 10, 41, 14, 15, 32, 34, 26, 27, 28, 29/), 0.5, nslic, ntim, slice, file_output )
!
!  Analyze thermals
!
!  call thermals ( zz, (/3, 8, 10, 41, 30, 32, 33, 34, 18, 19, 20, 25/), 0.5, nslic, ntim, slice, './thermals_rh80_dbtv.dat                     ' )
!
!----------------------------------------------------------!
!		     Deallocate and end                    !
!----------------------------------------------------------!
!
  deallocate( xx, yy, zz, tt )
  deallocate( wind, state, tend )
  deallocate( slice )
!
!----------------------------------------------------------!
!  
  end program slices
