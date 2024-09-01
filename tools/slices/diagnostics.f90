
module diagnostics

  USE sharedmod
  
  implicit none
        
  private
  
  public :: empty_space, space_correlations, time_correlations, time_corr_2D, distance_stats, distributions, cloud_ids
  
  contains
!
  SUBROUTINE cloud_ids ( thres, slice )
    
    real, intent(in) :: thres
    real, dimension(NXO,NYO,nvarslic), intent(inout) :: slice
!    
    real, dimension(NXO,NYO) :: old
    real :: ref1, ref2, size(50000)
    integer :: i, j, k, m, npts, nclouds, minm
!
!  Label clouds
!
    nclouds = 0
    slice(:,:,21) = 0.0
    do j = 2, NYO
      do i = 2, NXO
        if ( slice(i,j,13) > thres ) then
  	  if ( slice(max(i-1,2),j,21) > 0.5 .or. slice(i,max(j-1,2),21) > 0.5 ) then
  	    slice(i,j,21) = max(slice(max(i-1,2),j,21),slice(i,max(j-1,2),21))
	    if (slice(i,j,21) == slice(i,max(j-1,2),21) .and. slice(max(i-1,2),j,21) > 0.5)    &
	    	where (slice(:,:,21) == slice(max(i-1,2),j,21)) slice(:,:,21) = slice(i,j,21)
	    if (slice(i,j,21) == slice(max(i-1,2),j,21) .and. slice(i,max(j-1,2),21) > 0.5)    &
	    	where (slice(:,:,21) == slice(i,max(j-1,2),21)) slice(:,:,21) = slice(i,j,21)
  	  else
	    nclouds = nclouds + 1
	    slice(i,j,21) = real(nclouds)
	  endif
	else
	  slice(i,j,21) = 0.
        endif    
      enddo
    enddo
!  
!  i=1 boundary
!
    do j = 2, NYO
      if ( slice(1,j,13) > thres ) then
        if ( slice(2,j,21) > 0.5 .or. slice(1,max(j-1,2),21) > 0.5 ) then
    	  slice(1,j,21) = min(slice(2,j,21),slice(1,max(j-1,2),21))
	  if (slice(1,j,21) == slice(1,max(j-1,2),21) .and. slice(2,j,21) > 0.5)    &
	  	where (slice(:,:,21) == slice(2,j,21)) slice(:,:,21) = slice(1,j,21)
	  if (slice(1,j,21) == slice(2,j,21) .and. slice(1,max(j-1,2),21) > 0.5)    &
	  	where (slice(:,:,21) == slice(1,max(j-1,2),21)) slice(:,:,21) = slice(1,j,21)
        else
    	  nclouds = nclouds + 1
    	  slice(1,j,21) = real(nclouds)
        endif
      endif  
! 
      if ( slice(NXO,j,21) > 0. ) then
        ref1 = slice(1,j,21)
	ref2 = slice(NXO,j,21)
        where (slice(:,:,21) == ref2) slice(:,:,21) = ref1
      endif  
    enddo
!  
!  j = 1, boundary
!
    do i = 2, NXO
      if ( slice(i,1,13) > thres ) then
        if ( slice(i-1,1,21) > 0.5 .or. slice(i,2,21) > 0.5 ) then
          slice(i,1,21) = min(slice(max(i-1,2),1,21),slice(i,2,21))
	  if (slice(i,1,21) == slice(max(i-1,2),1,21) .and. slice(i,2,21) > 0.5)    &
	  	where (slice(:,:,21) == slice(i,2,21)) slice(:,:,21) = slice(i,1,21)
	  if (slice(i,1,21) == slice(i,2,21) .and. slice(max(i-1,2),1,21) > 0.5)    &
	  	where (slice(:,:,21) == slice(max(i-1,2),1,21)) slice(:,:,21) = slice(i,1,21)
        else
          nclouds = nclouds + 1
          slice(i,1,21) = real(nclouds)
        endif
      endif    
! 
      if ( slice(i,NYO,21) > 0. ) then
        ref1 = slice(i,1,21)
	ref2 = slice(i,NYO,21)
        where (slice(:,:,21) == ref2) slice(:,:,21) = ref1
      endif  
    enddo
!
!  Get cloud sizes and relabel
!
    size = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( slice(i,j,21) > 0. ) then
	  size(int(slice(i,j,21))) = size(int(slice(i,j,21))) + 1.
	endif  
      enddo
    enddo
!      
    npts = 0
    slice(:,:,22) = 0.0
    do m = 1, nclouds
      if ( size(m) > 0. ) then
        where (slice(:,:,21) == real(m)) slice(:,:,22) = size(m)
	npts = npts + 1 
      endif
    enddo
!
    call sort(size,50000)
    
    minm = 1
    do while (size(minm) == 0)
      minm = minm + 1
    enddo
!
    old = slice(:,:,21)
    slice(:,:,21) = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( slice(i,j,22) > 0. .and. slice(i,j,21) == 0. ) then
          do k = minm, 50000
	    if (slice(i,j,22) == size(k) .and. size(k) /= 0.) then
	      where (old == old(i,j)) slice(:,:,21) = real(k - minm + 1)
	      size(k) = 0.
	      exit
	    endif 
	  enddo
	endif
      enddo
    enddo
!    
    print*, 'Number of clouds:', maxval(slice(:,:,21)), maxval(slice(:,:,22))
    
!    where (slice(:,:,21) > 0.) slice(:,:,21) = 1.
!  
  END SUBROUTINE cloud_ids

  SUBROUTINE empty_space ( npts, thres, slice, file_output, lout )
    
    integer, intent(in) :: npts
    character(len = 100), intent(in) :: file_output
    real, intent(in) :: thres
    real, dimension(NXO,NYO,nvarslic), intent(inout) :: slice
    logical, intent(in) :: lout
!    
    real :: pi=3.141592653589793
    real :: X1, X2, x, y, R, dist, maxdist
    real, dimension(npts) :: mindist
    real, dimension(npts,2) :: cumul
    integer :: values(1:8), k, n, unit1, ndraw
    integer :: i, j, inew, jnew, i1, j1
    integer, allocatable :: seed(:)
    logical :: ex
!    
    call date_and_time(values=values)
!
    call random_seed(size=k)
    allocate( seed(1:k) )
    seed(:) = values(8)
    call random_seed(put=seed)
!    
    ndraw = 0
    maxdist = 1.e6
    mindist = maxdist
!
    do n = 1, npts
!    
      CALL RANDOM_NUMBER( X1 )
      CALL RANDOM_NUMBER( X2 )
!      
      x = X1*xx(NXO)
      y = X2*yy(NYO)
      i = int(X1*real(NXO)) + 1
      j = int(X2*real(NYO)) + 1
!      
      R = dx
      if ( slice(i,j,13) < thres ) then
!      
        ndraw = ndraw + 1
        do jnew = j-100, j+100        
          do inew = i-100, i+100
      	    i1 = inew
      	    j1 = jnew
      	    if ( inew < 1 ) i1 = i1 + NXO
      	    if ( inew > NXO ) i1 = i1 - NXO
      	    if ( jnew < 1 ) j1 = j1 + NYO
      	    if ( jnew > NYO ) j1 = j1 - NYO
!	    
	    if ( slice(i1,j1,13) >= thres ) then
	      dist = R*sqrt( real(i - inew)**2. + real(j - jnew)**2. )
	      mindist(n) = min(mindist(n),dist)
	    endif
          enddo
        enddo
!
      endif
!  
      if ( mindist(n) < maxdist ) slice(i,j,14) = mindist(n)
!  
    enddo
!
    where ( mindist == maxdist ) mindist = 0.
    mindist = pi*mindist**2.
!
    call Sort ( mindist, npts )
!
    k = 0
    cumul = 0.
    do n = npts-ndraw, npts
      if ( mindist(n) > 0. .and. mindist(n) /= cumul(k,1) ) then
        k = k + 1
        cumul(k,1) = mindist(n)
        cumul(k,2) = real( (npts - n) + 1 ) / real( ndraw + 1 )
      endif
    enddo
!
!  Open file and write output
!
    if ( lout ) then
      unit1=201
      inquire (FILE=trim(file_output), OPENED=ex)
      if (.not.ex) then
        open(unit1,file=file_output(1:len_trim(file_output)),	 &
               form='formatted', status='unknown')
      endif
!
      do i = 1, k
        write(unit1,100) cumul(i,1), cumul(i,2)
      enddo
!
      close (unit1)
    endif
!
100 format (2(f12.4,1x))
!    
    deallocate ( seed )
!  
  END SUBROUTINE empty_space

  SUBROUTINE time_correlations ( nd, n1, n2, thres, slice, tstep, file_output, lout )
    
    real, intent(in) :: thres, tstep
    real, dimension(NXO,NYO,nvarslic,ntime), intent(inout) :: slice
    integer, intent(in) :: n1, n2, nd
    character(len = 100), intent(in) :: file_output
    logical, intent(in) :: lout
!    
    real :: mean1, mean2, mean3, mm, mm1, vv1, ss1, var1, var2, var3
    real :: minv, maxv, loc, dv, R, dist, maxdist 
    real, dimension(NXO,NYO) :: cov, mindist, flag
    real, dimension(nd,3) :: stat
!    
    logical :: ex
    integer :: k, unit1, nn(nd), nn0(nd), nn1(nd)
    integer :: i, j, t, inew, jnew, i1, j1, npts1, npts2, npts3, nt
    character(len=3) :: car3
!
!  Point distance to clouds
!
    mindist = 0.
    maxdist = 1.e6
    do j = 1, NYO
      do i = 1, NXO
!    
      if ( slice(i,j,13,1) < thres ) then
        dist = maxdist
        do jnew = j-50, j+50       
          do inew = i-50, i+50
      	    i1 = inew
      	    j1 = jnew
      	    if ( inew < 1 ) i1 = i1 + NXO
      	    if ( inew > NXO ) i1 = i1 - NXO
      	    if ( jnew < 1 ) j1 = j1 + NYO
      	    if ( jnew > NYO ) j1 = j1 - NYO
!	    
	    if ( slice(i1,j1,13,1) >= thres ) then
	      R = dx*sqrt( real(i - inew)**2. + real(j - jnew)**2. )
	      dist = min(dist,R)
	    endif
          enddo
        enddo
      endif
!  
      if ( dist < maxdist ) mindist(i,j) = dist
      enddo
    enddo
!
!  Initiliaze stats
!
    flag = 0.
    npts1 = 0
    mean1 = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres ) then
          flag(i,j) = 1.
	  npts1 = npts1 + 1
	  mean1 = mean1 + slice(i,j,n1,1)
	endif
      enddo
    enddo
    mean1 = mean1 / real(npts1)
!
    var1 = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres ) then
	  var1 = var1 + (slice(i,j,n1,1) - mean1)**2.
        endif
      enddo
    enddo
    var1 = var1 / real(npts1)
!	  
    minv =  1000000.
    maxv = -1000000.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres ) then
	  minv = min(slice(i,j,n1,1),minv)
	  maxv = max(slice(i,j,n1,1),maxv)
	endif
      enddo
    enddo
    minv = minv - 1.e-6
    maxv = maxv + 1.e-6
!
!  Loop over all time steps
!
    do t = 1, ntime-1
    nt = t + 1
!
!  Initiliaze stats
!
    npts2 = 0
    mean2 = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres .and. flag(i,j) > 0.5 ) then
  	  npts2 = npts2 + 1
	  mean2 = mean2 + slice(i,j,21,nt)	
	endif
      enddo
    enddo
    mean2 = mean2 / real(npts2)
!
    var2 = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres .and. flag(i,j) > 0.5 ) then
	  var2 = var2 + (slice(i,j,21,nt) - mean2)**2.
        endif
      enddo
    enddo
    var2 = var2 / real(npts2)
!
    npts3 = 0
    mean3 = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres .and. flag(i,j) > 0.5 .and. slice(i,j,n2,nt) > 0. ) then
  	  npts3 = npts3 + 1
	  mean3 = mean3 + slice(i,j,n2,nt)	
	endif
      enddo
    enddo
    mean3 = mean3 / real(npts3)
!
    var3 = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres .and. flag(i,j) > 0.5 .and. slice(i,j,n2,nt) > 0. ) then
	  var3 = var3 + (slice(i,j,n2,nt) - mean3)**2.
        endif
      enddo
    enddo
    var3 = var3 / real(npts3)
!
    print*, 'Step',t,'Mean & variance for scalars 1 and 2:', mean1, sqrt(var1), mean2, sqrt(var2), mean3, sqrt(var3)
!
    nn = 0
    mm = 0.
    nn0 = 0
    nn1 = 0
    mm1 = 0.
    stat = 0.
    dv = abs(maxv - minv) / real(nd-1)
    do j = 1, NYO
      do i = 1, NXO
        do k = 1, nd
          loc = minv + real(k-1)*dv
	  if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres .and. (slice(i,j,n1,1) >= loc) .and. (slice(i,j,n1,1) < loc+dv) ) then
            nn0(k) = nn0(k) + 1
	    if ( maxval(slice(i,j,13,1:nt)) > thres ) mm1 = mm1 + slice(i,j,n1,1)
!
            if ( flag(i,j) > 0.5 ) then
              nn(k) = nn(k) + 1
	      mm = mm + slice(i,j,n1,1)
              stat(k,1) = stat(k,1) + slice(i,j,n1,1)
              stat(k,2) = stat(k,2) + slice(i,j,21,nt)
!
              if ( slice(i,j,n2,nt) > 0. ) then
                nn1(k) = nn1(k) + 1
                stat(k,3) = stat(k,3) + slice(i,j,n2,nt)
	        flag(i,j) = 0.
              endif
	    endif
	  endif
        enddo
      enddo
    enddo
!
    stat(:,1) = stat(:,1)/real(max(nn,1))
    stat(:,2) = (stat(:,2)/real(max(nn,1)) - mean2) / sqrt( var2 )
    stat(:,3) = (stat(:,3)/real(max(nn1,1)) - mean3) / sqrt( var3 )
!
!  Open file and write output
!
    print*, 'Mean & std dev. of triggered points:', mm/real(sum(nn)), mm1/real(sum(nn1))
!
    if ( lout ) then
      unit1=201
      inquire (FILE=trim(file_output), OPENED=ex)
!
      write(car3,'(i3)') t*int(tstep/60.) 
      if (.not.ex) then
        open(unit1,file=file_output(1:len_trim(file_output))//'_'//trim(adjustl(car3)),	 &
               form='formatted', status='unknown')
      endif
! 
      do k = 1, nd
        if ( nn1(k) > 2 ) write(unit1,*) stat(k,1), stat(k,2), stat(k,3), real(nn0(k))/real(sum(nn0)), real(nn1(k))/real(sum(nn0)), real(nn(k))/real(sum(nn0))
      enddo
!
      close (unit1)
    endif
!
!  End loop over all time steps
!
    enddo
!  
  END SUBROUTINE time_correlations

  SUBROUTINE time_corr_2D ( nd, n1, n2, n3, thres, slice, tstep, file_output, lout )
    
    real, intent(in) :: thres, tstep
    real, dimension(NXO,NYO,nvarslic,ntime), intent(inout) :: slice
    integer, intent(in) :: n1, n2, n3, nd
    character(len = 100), intent(in) :: file_output
    logical, intent(in) :: lout
!    
    real :: minv1, maxv1, minv2, maxv2, dv1, dv2, loc1, loc2
    real :: dist, maxdist, R
    real, dimension(nd,nd,4) :: stat
    real, dimension(NXO,NYO) :: mindist, flag
!    
    logical :: ex
    integer :: unit1, nn(nd,nd), nn0(nd,nd)
    integer :: i, j, k, l, t, nt, i1, j1, inew, jnew
    character(len=3) :: car3
!
!  Point distance to clouds
!
    flag = 0.
    mindist = 0.
    maxdist = 1.e6
    do j = 1, NYO
      do i = 1, NXO
!    
      if ( slice(i,j,13,1) < thres ) then
        flag(i,j) = 1.
        dist = maxdist
        do jnew = j-50, j+50       
          do inew = i-50, i+50
      	    i1 = inew
      	    j1 = jnew
      	    if ( inew < 1 ) i1 = i1 + NXO
      	    if ( inew > NXO ) i1 = i1 - NXO
      	    if ( jnew < 1 ) j1 = j1 + NYO
      	    if ( jnew > NYO ) j1 = j1 - NYO
!	    
	    if ( slice(i1,j1,13,1) >= thres ) then
	      R = dx*sqrt( real(i - inew)**2. + real(j - jnew)**2. )
	      dist = min(dist,R)
	    endif
          enddo
        enddo
      endif
!  
      if ( dist < maxdist ) mindist(i,j) = dist
      enddo
    enddo
!
!  Loop over all time steps
!
    do t = 1, ntime-1
    nt = t + 1
!
!  Initiliaze stats
!	  
    minv1 =  1000000.
    maxv1 = -1000000.
    minv2 =  1000000.
    maxv2 = -1000000.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres ) then
!	  if ( slice(i,j,13,nt) > thres .and. flag(i,j) > 0.5 ) then
	    minv1 = min(slice(i,j,n1,1),minv1)
	    maxv1 = max(slice(i,j,n1,1),maxv1)
	    minv2 = min(slice(i,j,n2,1),minv2)
	    maxv2 = max(slice(i,j,n2,1),maxv2)
!          endif
	endif
      enddo
    enddo
    maxv1 = maxv1 + 1.e-12
    maxv2 = maxv2 + 1.e-12
!
    dv1 = abs(maxv1 - minv1) / real(nd-1)
    dv2 = abs(maxv2 - minv2) / real(nd-1)
!
    nn = 0
    nn0 = 0
    stat = 0.
    do j = 1, NYO
      do i = 1, NXO
        do k = 1, nd
          loc1 = minv1 + real(k-1)*dv1
          do l = 1, nd
            loc2 = minv2 + real(l-1)*dv2
	    if ( mindist(i,j) > 1000. .and. slice(i,j,13,1) < thres .and.  		&
	        (slice(i,j,n1,1) >= loc1) .and. (slice(i,j,n1,1) < loc1+dv1) .and.  	&
		(slice(i,j,n2,1) >= loc2) .and. (slice(i,j,n2,1) < loc2+dv2) ) then
!
              nn0(k,l) = nn0(k,l) + 1
              stat(k,l,1) = stat(k,l,1) + slice(i,j,n1,1)
              stat(k,l,2) = stat(k,l,2) + slice(i,j,n2,1)
!
              if ( slice(i,j,13,nt) > thres .and. flag(i,j) > 0.5 ) then
	        nn(k,l) = nn(k,l) + 1
          	stat(k,l,3) = stat(k,l,3) + slice(i,j,n3,nt)
		flag(i,j) = 0.
              endif
!
	    endif
          enddo
        enddo
      enddo
    enddo
!
    print*, 'Average scalars 1 and 2 at time', t, ':', sum(stat(:,:,1))/real(sum(nn0)), sum(stat(:,:,2))/real(sum(nn0))
!
    where (nn0 > 0) stat(:,:,1) = stat(:,:,1) / real(nn0)
    where (nn0 > 0) stat(:,:,2) = stat(:,:,2) / real(nn0)
    where (nn0 > 0) stat(:,:,4) = real(nn) / real(nn0)
    where (nn > 0) stat(:,:,3) = stat(:,:,3) / real(nn)
!
    print*, 'Total number of environmental and cloudy grid points at time', t, ':', sum(nn0), sum(nn)
!
!  Open file and write output
!
    if ( lout ) then
      unit1=201
      inquire (FILE=trim(file_output), OPENED=ex)
!
      write(car3,'(i3)') t*int(tstep/60.) 
      if (.not.ex) then
        open(unit1,file=file_output(1:len_trim(file_output))//'_'//trim(adjustl(car3)),	 &
               form='formatted', status='unknown')
      endif
! 
      do k = 1, nd
        loc1 = minv1 + real(k-1)*dv1
        do l = 1, nd
          loc2 = minv2 + real(l-1)*dv2
	  if (nn(k,l) > 5) write(unit1,*) loc1+dv1/2., loc2+dv2/2., stat(k,l,3), stat(k,l,4), real(nn0(k,l)), real(nn(k,l))
	enddo
      enddo
!
      close (unit1)
    endif
!
!  End loop over all time steps
!
    enddo
!  
  END SUBROUTINE time_corr_2D

  SUBROUTINE space_correlations ( nd, n1, n2, thres, slice, file_output, lout )
    
    real, intent(in) :: thres
    real, dimension(NXO,NYO,nvarslic), intent(inout) :: slice
    integer, intent(in) :: n1, n2, nd
    character(len = 100), intent(in) :: file_output
    logical, intent(in) :: lout
!    
    real :: rr, mean1, mean2, std1, std2
    real, dimension(NXO,NYO,(2*nd+1)*(2*nd+1)) :: cov
    real, dimension(2*nd+1,2*nd+1) :: dist
    real, dimension((2*nd+1)*(2*nd+1)) :: dist1d, rcov, mcov
!    
    logical :: ex
    integer, dimension(2*nd+1,2*nd+1) :: index
    integer, dimension((2*nd+1)*(2*nd+1)) :: index1d
    integer :: values(1:8), k, unit1, nt, npts, maxk
    integer :: i, j, inew, jnew, i1, j1
    integer, allocatable :: seed(:)
!
!  Initialize distances and indices
!
    k = 0
    dist = 0.
    do jnew = -nd, nd	
      do inew = -nd, nd
        k = k+1
        rr = dx*sqrt( real(inew)**2. + real(jnew)**2. )
	dist(inew+nd+1,jnew+nd+1) = rr
	dist1d(k) = rr
      enddo
    enddo
!
    k = 1
    nt = 2*nd+1
    rcov(1) = dist1d(1)
    where ( dist == dist1d(1) ) index = 1
    do i = 2, nt*nt
      if ( .not.any( dist1d(i) == rcov(1:k) ) ) then
        k = k + 1
	rcov(k) = dist1d(i)
        where ( dist == dist1d(i) ) index = k
      endif
    enddo
    maxk = k
!
!  Initiliaze stats
!
    npts = 0
    mean1 = 0.
    mean2 = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( slice(i,j,13) < thres ) then
	  npts = npts + 1
	  mean1 = mean1 + slice(i,j,n1)
	  mean2 = mean2 + slice(i,j,n2)	
	endif
      enddo
    enddo
    mean1 = mean1 / real(npts)
    mean2 = mean2 / real(npts)
!    
    std1 = 0.
    std2 = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( slice(i,j,13) < thres ) then
	  std1 = std1 + (slice(i,j,n1) - mean1)**2.
	  std2 = std2 + (slice(i,j,n2) - mean2)**2.	
	endif
      enddo
    enddo
    std1 = sqrt( std1 / real(npts) )
    std2 = sqrt( std2 / real(npts) )
!
!  Calculate covariances
!	  
    cov = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( slice(i,j,13) < thres ) then
          index1d = 0
          do jnew = j-nd, j+nd        
            do inew = i-nd, i+nd
              i1 = inew
              j1 = jnew
              if ( inew < 1 ) i1 = i1 + NXO
              if ( inew > NXO ) i1 = i1 - NXO
              if ( jnew < 1 ) j1 = j1 + NYO
              if ( jnew > NYO ) j1 = j1 - NYO
!     	      
      	      k = index(inew-i+nd+1,jnew-j+nd+1)
     	      index1d(k) = index1d(k) + 1
      	      cov(i,j,k) = cov(i,j,k) + (slice(i,j,n1) - mean1) * (slice(i1,j1,n2) - mean2)
      	    enddo
          enddo
!
      	  do k = 1, maxk
      	    if (index1d(k) > 0) cov(i,j,k) = cov(i,j,k) / real(index1d(k))
      	  enddo
        endif
      enddo
    enddo
!
!  Normalization of cross-correlation coefficients
!
    mcov = 0.
    npts = 0
    do j = 1, NYO
      do i = 1, NXO
        if ( slice(i,j,13) < thres ) then
	  npts = npts + 1
          mcov = mcov + cov(i,j,:)
        endif
      enddo
    enddo
    mcov = mcov / real(npts) / std1 / std2
!
!  Open file and write output
!
    if ( lout ) then
      unit1=201
      inquire (FILE=trim(file_output), OPENED=ex)
      if (.not.ex) then
        open(unit1,file=file_output(1:len_trim(file_output)),	 &
               form='formatted', status='unknown')
      endif
! 
      do k = 1, maxk
        write(unit1,*) rcov(k), mcov(k)
      enddo
!
      close (unit1)
    endif
!  
  END SUBROUTINE space_correlations

  SUBROUTINE distributions ( index, thres, slice, file_output )
    
    character(len = 100), intent(in) :: file_output
    real, intent(in) :: thres
    real, dimension(NXO,NYO,nvarslic), intent(inout) :: slice
    integer, intent(in) :: index
!    
    real :: maxv, minv, dv, pdfv(2,100)
    integer :: i, j, k, unit1, npts=100
    logical :: ex
!
    maxv = 0.
    minv = 1.e6
    do j = 1, NYO
      do i = 1, NXO
        if ( slice(i,j,13) < thres ) then
          if (slice(i,j,index) > maxv) maxv = slice(i,j,index)
	  if (slice(i,j,index) < minv) minv = slice(i,j,index)
        endif
      enddo
    enddo
!    
    dv = (maxv-minv)/real(npts-1)
!    
    pdfv = 0.
    do j = 1, npts
      pdfv(1,j) = minv+real(j-1)*dv
    enddo   
!
    do j = 1, NYO
      do i = 1, NXO
        if ( slice(i,j,13) < thres ) then
	  do k = 1, npts
            if (slice(i,j,index) >= pdfv(1,k)-0.5*dv .and. slice(i,j,index) < pdfv(1,k)+0.5*dv) &
	       pdfv(2,k) = pdfv(2,k) + 1.
          enddo
	endif
      enddo
    enddo
    pdfv(2,:) = pdfv(2,:)/sum(pdfv(2,:))
!
!  Open file and write output
!
    unit1=201
    inquire (FILE=trim(file_output), OPENED=ex)
    if (.not.ex) then
      open(unit1,file=file_output(1:len_trim(file_output)), form='formatted', status='unknown')
    endif
!
    do k = 1, npts
      write(unit1,*) pdfv(1,k), pdfv(2,k)
    enddo
!
    close (unit1)
!  
    return
  END SUBROUTINE distributions
!
  SUBROUTINE distance_stats ( nd, index1, index2, index3, index4, thres, slice, file_output )
    
    character(len = 100), intent(in) :: file_output
    real, intent(in) :: thres
    real, dimension(NXO,NYO,nvarslic), intent(inout) :: slice
    integer, intent(in) :: index1, index2, index3, index4, nd
!    
    real :: R, maxdist, dist, scl, ncl, maxncl, maxscl
    real, dimension(NXO,NYO) :: mindist, scloud, ncloud
    real, dimension(NXO*NYO,4) :: store
    integer :: i, j, k, inew, jnew, i1, j1, icl
!    
    scloud = 0.
    mindist = 0.
    maxdist = 1.e6
!
    do j = 1, NYO
      do i = 1, NXO
!    
      if ( slice(i,j,13) < thres ) then
        dist = maxdist
	scl = 0.
	ncl = 0.
        do jnew = j-nd, j+nd        
          do inew = i-nd, i+nd
      	    i1 = inew
      	    j1 = jnew
      	    if ( inew < 1 ) i1 = i1 + NXO
      	    if ( inew > NXO ) i1 = i1 - NXO
      	    if ( jnew < 1 ) j1 = j1 + NYO
      	    if ( jnew > NYO ) j1 = j1 - NYO
!	    
	    if ( slice(i1,j1,13) >= thres ) then
	      R = dx*sqrt( real(i - inew)**2. + real(j - jnew)**2. )
	      if ( min(dist,R) == R ) then
	        ncl = slice(i1,j1,21)
	        scl = slice(i1,j1,22)
	      endif
	      dist = min(dist,R)
	    endif
          enddo
        enddo
      endif
!  
      if ( dist < maxdist ) mindist(i,j) = dist
      if ( dist < maxdist ) scloud(i,j) = scl
      if ( dist < maxdist ) ncloud(i,j) = ncl
!  
      enddo
    enddo
    maxncl = maxval(ncloud)
!
!  Write stats for all clouds
!
    call write_distributions (0., (/index1, index2, index3, index4/), file_output, slice, mindist, scloud)
!
!  Write stats for 20% largest clouds
!
    icl = floor(0.8*maxncl)
    do j = 1, NYO
      do i = 1, NXO      
        if (slice(i,j,21) == real(icl)) then 
	  maxscl = slice(i,j,22)
          exit
	endif
      enddo
    enddo
!
    call write_distributions (maxscl, (/index1, index2, index3, index4/), file_output, slice, mindist, scloud)
!
!  Write stats for 20% smallest clouds
!
    icl = floor(0.2*maxncl)
    do j = 1, NYO
      do i = 1, NXO      
        if (slice(i,j,21) == real(icl)) then 
	  maxscl = slice(i,j,22)
          exit
	endif
      enddo
    enddo
!
    call write_distributions (-maxscl, (/index1, index2, index3, index4/), file_output, slice, mindist, scloud)
!    
    return
!    
    contains
!
    subroutine write_distributions (limsize, index, file_output, slice, mindist, scloud)
    
    integer :: index(4)
    character(len = 100) :: file_output
    real, dimension(NXO,NYO,nvarslic) :: slice
    real, dimension(NXO,NYO) :: mindist, scloud
    real :: limsize
!    
    logical :: ex
    integer :: i, j, k, l, n, unit1, nn
    character(len = 20) :: suffix
    character(len = 1) :: car1
    real, dimension(NXO*NYO,4) :: store
    real :: mean
!
    if (limsize == 0.) then
      suffix = 'tot'
    else if (limsize > 0.) then
      suffix = 'large'
    else if (limsize < 0.) then
      suffix = 'small'
      scloud = -scloud
    endif
!
    do n = 1, 4
    write(car1,'(i1)') n
!
    nn = 0
    mean = 0.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 0. .and. mindist(i,j) < 2500. .and. slice(i,j,21) < thres ) then
	  nn = nn + 1
	  mean = mean + slice(i,j,index(n))
	endif
      enddo
    enddo
    mean = mean / real(nn)
!
    print*, 'Mean scalar value', n, ':', mean
!
    k = 0
    store = 0.
    store(:,4) = 1.
    do j = 1, NYO
      do i = 1, NXO
        if ( mindist(i,j) > 0. .and. mindist(i,j) < 2500. .and. slice(i,j,21) < thres .and. .not.any(mindist(i,j) == store(:,1)) .and. scloud(i,j) >= limsize ) then
	  k = k+1
	  store(k,1) = mindist(i,j)
	  store(k,2) = slice(i,j,index(n)) !- mean
	  store(k,3) = (slice(i,j,index(n)) - mean)**2.
	else if ( mindist(i,j) > 0. .and. slice(i,j,21) < thres .and. scloud(i,j) >= limsize ) then 
	  do l = 1, k
	    if ( mindist(i,j) == store(l,1) ) then
	      store(l,2) = store(l,2) + slice(i,j,index(n)) !- mean
	      store(l,3) = store(l,3) + (slice(i,j,index(n)) - mean)**2.
	      store(l,4) = store(l,4) + 1.
	      exit
	    endif
	  enddo
	endif
      enddo
    enddo
    store(:,2) = store(:,2) / store(:,4)
    store(:,3) = store(:,3) / store(:,4) - store(:,2)**2.
!
!  Open file and write output
!
    unit1=201
    inquire (FILE=trim(file_output), OPENED=ex)
    if (.not.ex) then
      open(unit1,file=file_output(1:len_trim(file_output))//'_'//car1//'_'//trim(suffix),      &
    	     form='formatted', status='unknown')
    endif
!
    do j = 1, NXO*NYO
      if ( store(j,1) > 0. ) write(unit1,*) store(j,1), store(j,2), store(j,3)
    enddo
!
    close (unit1)
!
    enddo
!    
    print*, 'Written distributions for category:', suffix
!
    end subroutine write_distributions
!  
  END SUBROUTINE distance_stats
!  
  SUBROUTINE  Sort(x, Size)

  IMPLICIT  NONE
  INTEGER		 :: Size
  INTEGER		 :: i
  INTEGER		 :: Location
  REAL, DIMENSION(1:Size) :: x

  DO i = 1, Size-1		    ! except for the last
   CALL FindMinimum (x, i, Size, Location)	    ! find min from this to last
   CALL Swap (x(i), x(Location))  ! swap this and the minimum
  END DO

  END SUBROUTINE  Sort

  SUBROUTINE  Swap(a, b)

  IMPLICIT  NONE
  REAL :: a, b
  REAL :: Temp

  Temp = a
  a    = b
  b    = Temp

  END SUBROUTINE  Swap

  SUBROUTINE  FindMinimum(x, Start, End, Location)

  IMPLICIT  NONE
  INTEGER		 :: Start, End
  INTEGER		 :: Location
  INTEGER		 :: i
  REAL  		 :: Minimum
  REAL, DIMENSION(1:End) :: x

  Minimum  = x(Start)		    ! assume the first is the min
  Location = Start		    ! record its position
  DO i = Start+1, End		    ! start with next elements
    IF (x(i) < Minimum) THEN	    !	if x(i) less than the min?
      Minimum  = x(i) 	   	    !	   Yes, a new minimum found
      Location = i		    !	   record its position
    END IF
  END DO

  END SUBROUTINE  FindMinimum

end module diagnostics

