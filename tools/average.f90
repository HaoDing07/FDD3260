program average
!
!
!  Calculates time averages of MIMICA time-series
!
!
implicit none

  character(len=20) :: file_name
  character(len=15) :: car15
  character(len=4)  :: car41, car42, car43
  character(len=15), allocatable, dimension(:)  :: namevar
! 
  integer  :: i, j, k, l, nvar, ii, jj, kk, nall
  integer  :: nt, nv, nl, stat, err
  integer  :: iqi13, iqi12, iqi10, iqi8
!
  real  :: stime, etime, x(3)
  real, allocatable, dimension(:)     :: all
  real, allocatable, dimension(:,:)   :: xx, zz
  real, allocatable, dimension(:,:,:) :: var
  real, allocatable, dimension(:,:)   :: mean
  real, allocatable, dimension(:,:)   :: tab
  real, allocatable, dimension(:)     :: averages, ntime
!
!  First diagnostics or profiles time series
!
  print*, 'File name:'
  read*, file_name
  print*, 'Start time for averaging:'
  read*, stime
  print*, 'Final time for averaging:'
  read*, etime
!  
  open(UNIT=10, FILE=file_name, STATUS='OLD')
!  
  read(10,*)
  read(10,*) car15
  k = 1
  do while (car15(1:4) /= 'ZONE')
    read(10,*) car15
    k = k + 1
  enddo
  nvar = k - 4
  allocate(namevar(1:nvar))
  rewind(10)
!
  read(10,*)
  read(10,*) 
  read(10,*) 
  read(10,*) car15
  do k = 1, nvar
    read(10,*) car15
    i = index(car15(1:15),'.')
    j = index(car15(1:i+2),' ')
    l = index(car15(1:15),'pro')
    if (l /= 0) then
      namevar(k) = car15(1:i-1)  
    else
      if (j /= 0) then
        namevar(k) = car15(1:i-1)//car15(j+1:j+1)
      else
        namevar(k) = car15(1:i-1)//car15(i+1:i+2)
      endif
    endif
  enddo
!
  read(10,*)
  read(10,*) car41, ii, car42, jj, car43, kk
  read(10,*)
  nall = ii*kk*(nvar+3)
  allocate(xx(1:ii,1:kk), zz(1:ii,1:kk))
  allocate(var(1:ii,1:kk,1:nvar), all(1:nall))
!  
  read(10,*,iostat=err) all(1:nall)
!
  k = 0
  do j = 1, kk
    do i = 1, ii
      k = k + 1
      xx(i,j) = all(k)
    enddo
  enddo
!
  do j = 1, kk
    do i = 1, ii
      k = k + 1
    enddo
  enddo  
!
  do j = 1, kk
    do i = 1, ii
      k = k + 1
      zz(i,j) = all(k)
    enddo
  enddo
!
  do l = 1, nvar
    do j = 1, kk
      do i = 1, ii
        k = k + 1
        var(i,j,l) = all(k)      
      enddo
    enddo
  enddo
!
  close(10)
!
  allocate (mean(1:kk,1:nvar))
  mean = 0.
  do j = 1, kk
    do i = 1, ii
      if (xx(i,j) >= stime .and. xx(i,j) <= etime .and. abs(var(i,j,iqi8)) > 1.e-10) then
        var(i,j,iqi13) = (var(i,j,iqi8) - var(i-1,j,iqi8)) / (xx(i,j) - xx(i-1,j))
	var(i,j,iqi13) = var(i,j,iqi13) / var(i,j,iqi8)
      endif
    enddo  
!
    k = 0
    do i = 1, ii
      if (xx(i,j) >= stime .and. xx(i,j) <= etime) then
        k = k + 1
	mean(j,:) = mean(j,:) + var(i,j,:)
      endif
    enddo
    mean(j,:) = mean(j,:) / real(k)
  enddo
!  
  open(11,FILE='average_'//file_name(1:4)//'.out')
  write(11,102) 'ZZ', namevar(:)
  do j = 1, kk
    write(11,103) zz(1,j), mean(j,:)
  enddo
  close(11)
!
!
!  Now scalar time series
!
!
  nv = 16
  allocate (ntime(1:10000), averages(1:nv), tab(1:10000,1:nv))
!
!  Open file time
!
    file_name = 'OUTPUT/T_S'   
    open(UNIT=100,FILE=file_name,STATUS='OLD',IOSTAT=stat)
    if (stat > 0) then
      print*, 'The file ', file_name, ' does not exist, skipping instruction'
      goto 1000
    endif
!
!  Read data
!
    read(100,*)
!
    err = 0
    j = 1
    do while (err == 0)
      read(100,*,IOSTAT=err) ntime(j), tab(j,1:nv)
      j = j + 1
    enddo
    nl = j-1
    close(100)
!
!  calculate averages
!
  k = 0
  l = 0
  averages = 0.
  do j = 1, nl
    if (ntime(j) >= stime .and. ntime(j) <= etime) then
      averages(1:nv-1) = averages(1:nv-1) + tab(j,1:nv-1)
      k = k + 1
!
      if (tab(j,nv) < 0.02) then
        averages(nv) = averages(nv) + tab(j,nv)
        l = l + 1
      endif
    endif
  enddo
  averages(1:nv-1) = averages(1:nv-1) / real(k)
  averages(nv) = averages(nv) / real(l)
!
!  Write result
!
  open(UNIT=12, FILE='average_ts.out')
    write(12,*) averages
  close(12)
!
1000 continue
!
101 format(3(f15.7,3x))
102 format(3x,a2,13x,120(a18))
103 format(120(f15.7,3x))

end program
