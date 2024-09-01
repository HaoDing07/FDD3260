
module diagnostics

  USE sharedmod
  
  implicit none
        
  private
  
  public :: distance_stats, distributions, cloud_ids, thermals
  
  contains
!
  SUBROUTINE cloud_ids ( thres, z, maxl, slice )
    
    real, intent(in) :: thres, maxl, z(:)
    real, dimension(NYO,NZO,nvarslic), intent(inout) :: slice
!    
    real, dimension(NYO,NZO) :: old
    real :: ref1, ref2, size(20000), base(20000), top(20000), wmax(20000)
    integer :: i, j, k, m, npts, nclouds, minm
!
!  Label clouds
!
    nclouds = 0
    slice(:,:,21) = 0.0
    do j = 5, NZO
      do i = 2, NYO
        if ( (slice(i,j,13) > thres) .and. z(j) >= base_min .or. (slice(i,j-4,13) > thres .and. slice(i,j-1,21) > 0.5 .and. z(j) > top_max .and. .not.norm) ) then
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
        endif    
      enddo
    enddo
!  
!  i=1 boundary
!
    do j = 5, NZO
      if ( (slice(1,j,13) > thres) .and. z(j) >= base_min .or. (slice(1,j-4,13) > thres .and. slice(1,j-1,21) > 0.5 .and. z(j) > top_max .and. .not.norm) ) then
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
      if ( slice(NYO,j,21) > 0. .and. z(j) >= base_min ) then
        ref1 = slice(1,j,21)
	ref2 = slice(NYO,j,21)
        where (slice(:,:,21) == ref2) slice(:,:,21) = ref1
      endif  
    enddo
!
!  Get cloud sizes and relabel
!
    size = 0.
    wmax = 0.
    base = 1.e6
    top = 0.
    do i = 1, NYO
      do j = 1, NZO
        if (slice(i,j,21) > 0. .and. z(j) >= base_min ) then
	  top(int(slice(i,j,21))) = max(top(int(slice(i,j,21))),z(j))
	  base(int(slice(i,j,21))) = min(base(int(slice(i,j,21))),z(j))
	  wmax(int(slice(i,j,21))) = max(wmax(int(slice(i,j,21))), slice(i,j,3))
	  size(int(slice(i,j,21))) = top(int(slice(i,j,21))) - base(int(slice(i,j,21)))
	endif
      enddo
    enddo
!      
    npts = 0
    slice(:,:,22) = 0.0
    slice(:,:,23) = 0.0
    slice(:,:,24) = 0.0
    do m = 1, nclouds
      if ( base(m) > base_min .and. base(m) <= base_max ) then
        where (slice(:,:,21) == real(m)) 
	  slice(:,:,22) = size(m)
	  slice(:,:,23) = base(m)
	  slice(:,:,24) = wmax(m)
	end where
	npts = npts + 1 
      endif
    enddo
!
    call sort(base,20000)
!    
    minm = 1
    do while (base(minm) == 0)
      minm = minm + 1
    enddo
!
    old = slice(:,:,21)
    slice(:,:,21) = 0.
    do j = 1, NZO
      do i = 1, NYO
        if ( slice(i,j,23) > 0. .and. slice(i,j,21) == 0. ) then
          do k = minm, 20000
	    if (slice(i,j,23) == base(k) .and. base(k) < 1.e6) then
	      where (old == old(i,j)) slice(:,:,21) = maxl + real(k - minm + 1)
	      size(k) = 0.
    	      base(k) = 1.e6
	      exit
	    endif 
	  enddo
	endif
!
!	if (slice(i,j,21) > 0.5) slice(i,j,13) = 1.
      enddo
    enddo
!    
    print*, 'Number of clouds (cumulated):', int(maxval(slice(:,:,21))), 'and properties (not cumulated):', maxval(slice(:,:,22)), maxval(slice(:,:,23)), maxval(slice(:,:,24))
!  
  END SUBROUTINE cloud_ids

  SUBROUTINE distance_stats ( nd, z, index, thres, nslic, ntim, slice, file_output )
    
    character(len = 100), intent(in) :: file_output
    integer, intent(in) :: index(12), nd, ntim, nslic
    real, intent(in) :: thres, z(:)
    real, dimension(NYO,NZO,nvarslic,nslic,ntim), intent(inout) :: slice
!    
    real :: R, maxdist, dist, scl, ncl, wcl, bcl, maxncl, maxscl
    real, dimension(NYO,NZO,nslic,ntim) :: mindist, scloud, bcloud, ncloud, wcloud
    integer :: i, j, k, t, n, inew, i1, icl
!    
    scloud = 0.
    bcloud = 0.
    wcloud = 0.
    ncloud = 0.
    maxdist = 1.e6
    mindist = -maxdist
!
!  Limit the analysis to clouds with bases below 2500 m, above 750 m and depth larger than 500 m
!
    do n = 1, nslic
!
    do t = 1, ntim
!
    do j = 1, NZO
      do i = 1, NYO
!    
      scl = 0.
      wcl = 0.
      scl = 0.
      bcl = 0.
      if ( slice(i,j,13,n,t) < thres ) then
        dist = maxdist
        do inew = i-nd, i+nd
      	  i1 = inew
      	  if ( inew < 1 ) i1 = i1 + NYO
      	  if ( inew > NYO ) i1 = i1 - NYO
!	  
	  if ( slice(i1,j,13,n,t) > thres .and. slice(i1,j,22,n,t) >= depth_min .and. slice(i1,j,23,n,t) < base_max .and. slice(i1,j,23,n,t) > base_min ) then
!	  if ( slice(i1,j,3,n,t) >= 0. .and. slice(i1,j,22,n,t) >= depth_min .and. slice(i1,j,23,n,t) < base_max .and. slice(i1,j,23,n,t) > base_min ) then
	    R = dx*sqrt( real(i - inew)**2. )
	    if ( min(dist,R) == R ) then
	      ncl = slice(i1,j,21,n,t)
	      scl = slice(i1,j,22,n,t)
	      bcl = slice(i1,j,23,n,t)
	      wcl = slice(i1,j,24,n,t)
	      dist = min(dist,R)
	    endif
	  endif
        enddo
      else
        dist = -maxdist
        do inew = i-nd, i+nd
      	  i1 = inew
      	  if ( inew < 1 ) i1 = i1 + NYO
      	  if ( inew > NYO ) i1 = i1 - NYO
!	  
	  if ( slice(i1,j,13,n,t) < thres .and. slice(i,j,22,n,t) >= depth_min .and. slice(i,j,23,n,t) < base_max .and. slice(i,j,23,n,t) > base_min ) then
!	  if ( slice(i1,j,3,n,t) < 0. .and. slice(i,j,22,n,t) >= depth_min .and. slice(i,j,23,n,t) < base_max .and. slice(i,j,23,n,t) > base_min ) then
	    R = - dx*sqrt( real(i - inew)**2. )
	    if ( max(dist,R) == R ) then
	      ncl = slice(i,j,21,n,t)
	      scl = slice(i,j,22,n,t)
	      bcl = slice(i,j,23,n,t)
	      wcl = slice(i,j,24,n,t)
	      dist = max(dist,R)
	    endif
	  endif
        enddo
      endif
!  
      if ( dist < maxdist .and. dist > -maxdist ) mindist(i,j,n,t) = dist
      if ( dist < maxdist .and. dist > -maxdist ) ncloud(i,j,n,t) = ncl
      if ( dist < maxdist .and. dist > -maxdist ) scloud(i,j,n,t) = scl
      if ( dist < maxdist .and. dist > -maxdist ) bcloud(i,j,n,t) = bcl
      if ( dist < maxdist .and. dist > -maxdist ) wcloud(i,j,n,t) = wcl
!  
      enddo
    enddo
!
    enddo
!
    enddo
!
!    call write_prof (0., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_prof (2500., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_prof (5., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_prof (-1250., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_corr (0., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_corr (2500., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_corr (5., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_corr (-1250., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_mean (0., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_mean (3000., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_mean (6., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!    call write_mean (-1500., z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!  Write stats for all clouds
!
    call write_distributions (0., z, index, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!  Write stats for larger clouds
!
    call write_distributions (2200., z, index, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!  Write stats for activ
!
    call write_distributions (5., z, index, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!
!  Write stats for smaller clouds
!
    call write_distributions (-1200., z, index, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
!    
    return
!    
    contains
!
    subroutine write_prof ( lsize, z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud )
    
    integer :: nslic,ntim
    character(len = 100) :: file_output
    real :: thres, lsize
    real, dimension(NZO) :: z
    real, dimension(NYO,NZO,nvarslic,nslic,ntim) :: slice
    real, dimension(NYO,NZO,nslic,ntim) :: mindist, ncloud, scloud, bcloud, wcloud
!    
    logical :: ex
    character(len = 20) :: suffix
    integer :: i, j, l, t, m, h, maxnc, nh
    real :: dh, htot, height
    real, allocatable, dimension(:) :: alt
    real, allocatable, dimension(:,:) :: propcm, propem
    real, allocatable, dimension(:,:,:) :: propc, prope, count
    logical, allocatable, dimension(:,:) :: take
    logical, dimension(NYO,NZO,nslic,ntim) :: local
!
    local = .false.
    if (lsize == 0.) then
      suffix = 'tot'
      local = .true.
    else if (lsize == 5.) then
      suffix = 'activ'
      where ( scloud >= 2200. .and. wcloud > lsize ) local = .true.
    else if (lsize > 0.) then
      suffix = 'large'
      where ( scloud > lsize ) local = .true.
    else if (lsize < 0.) then
      suffix = 'small'
      where ( scloud < abs(lsize) ) local = .true.
    endif
!
    if (norm) then
      dh = 0.05
      htot = 1.
    else
      dh = 200.
      htot = 10000.
    endif
    nh = 1+int(htot/dh)
!
    maxnc = maxval(ncloud)
    allocate( propc(1:5,1:nh,1:maxnc), prope(1:5,1:nh,1:maxnc), propcm(1:5,1:maxnc), propem(1:5,1:maxnc) )
    allocate( alt(1:nh), count(2,1:nh,1:maxnc), take(1:nh,1:maxnc) )
!
    propc = 0.
    prope = 0.
    propcm = 0.
    propem = 0.
    count = 0.
    take = .true.
    alt = 0.
    do h = 2, nh
      alt(h) = alt(h-1) + dh
    enddo
!
    do j = 2, NZO
!
      do m = 1, nslic
        do t = 1, ntim
          do i = 1, NYO
!
            if ( local(i,j,m,t) .and. mindist(i,j,m,t) >= -dist_max .and. mindist(i,j,m,t) <= dist_max .and. z(j) <= top_max ) then 
              if (norm) then
	        height = (z(j) - bcloud(i,j,m,t)) / abs(scloud(i,j,m,t))
	      else
		height = z(j)
	      endif
	      do h = 1, nh
	        if ( height >= dh*real(h-1) .and. height < dh*real(h) ) exit
	      enddo
!
	      n = ncloud(i,j,m,t)
	      if ( mindist(i,j,m,t) < 0. ) then
		count(1,h,n) = count(1,h,n) + 1.
		propc(1,h,n) = propc(1,h,n) + slice(i,j,3,m,t)
		propc(2,h,n) = propc(2,h,n) + slice(i,j,10,m,t)
		propc(3,h,n) = propc(3,h,n) + slice(i,j,41,m,t)
		propc(4,h,n) = propc(4,h,n) + slice(i,j,14,m,t)
		propc(5,h,n) = propc(5,h,n) + slice(i,j,15,m,t)
	      else if ( slice(i,j,3,m,t) < 0. .and. take(h,n) ) then
		count(2,h,n) = count(2,h,n) + 1.
		prope(1,h,n) = prope(1,h,n) + slice(i,j,3,m,t)
		prope(2,h,n) = prope(2,h,n) + slice(i,j,10,m,t)
		prope(3,h,n) = prope(3,h,n) + slice(i,j,41,m,t)
		prope(4,h,n) = prope(4,h,n) + slice(i,j,14,m,t)
		prope(5,h,n) = prope(5,h,n) + slice(i,j,15,m,t) 
	      else
	        take(h,n) = .false.	     
	      endif
	    endif
!
	  enddo
        enddo
      enddo
!
    enddo
!
    do h = 1, nh
      do i = 1, 5
        if(sum(count(1,h,:)) > 0.) propcm(i,h) = sum(propc(i,h,:)) / sum(count(1,h,:))
        if(sum(count(2,h,:)) > 0.) propem(i,h) = sum(prope(i,h,:)) / sum(count(2,h,:))
        where(count(1,h,:) > 0.) propc(i,h,:) = propc(i,h,:) / count(1,h,:)
        where(count(2,h,:) > 0.) prope(i,h,:) = prope(i,h,:) / count(2,h,:)
      enddo
    enddo
!
    open(201,file=file_output(1:len_trim(file_output))//'_'//trim(suffix), form='formatted', status='unknown')
    open(202,file=file_output(1:len_trim(file_output))//'_scatter_'//trim(suffix), form='formatted', status='unknown')
    do h = 1, nh
      write(201,*) alt(h), propcm(:,h), propem(:,h)
      do n = 1, maxnc
        write(201,*) alt(h), propc(:,h,n), prope(:,h,n)
      enddo
    enddo
    close(201)
    close(202)
!
    deallocate( propc, prope, propcm, propem, count, alt, take )
!
    end subroutine write_prof
!
    subroutine write_mean ( lsize, z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud )
    
    integer :: nslic,ntim
    character(len = 100) :: file_output
    real :: thres, lsize
    real, dimension(NZO) :: z
    real, dimension(NYO,NZO,nvarslic,nslic,ntim) :: slice
    real, dimension(NYO,NZO,nslic,ntim) :: mindist, ncloud, scloud, bcloud, wcloud
!    
    logical :: ex, add
    character(len = 20) :: suffix
    integer :: i, i1, j, t, l, m, n, h, c, nh, unit1, count, count2, start, labels(4000)
    integer, dimension(:), allocatable :: npts
    real, dimension(:,:), allocatable :: pdf1, pdf2, pdf3
    real, dimension(NYO,NZO,nslic,ntim) :: tcloud
    real, dimension(:,:), allocatable :: store
    real, dimension(:,:,:), allocatable :: prop
    real :: height, dist, dh, htot, minw, minb, deltaw, deltab
    logical, dimension(NYO,NZO,nslic,ntim) :: local
!
    if (norm) then
      dh = 0.05
      htot = 1.
    else
      dh = 200.
      htot = 10000.
    endif
    nh = 1+int(htot/dh)
!    
    allocate( store(1:nh,5), prop(1:nh,1:20000,1:4), npts(1:nh) )
    allocate( pdf1(1:5000,1:4), pdf2(1:5000,1:4), pdf3(1:5000,1:4) )
!
    local = .false.
    if (lsize == 0.) then
      suffix = 'tot'
      local = .true.
    else if (lsize == 6.) then
      suffix = 'activ'
      where ( scloud >= 3000. .and. wcloud > lsize ) local = .true.
    else if (lsize > 0.) then
      suffix = 'large'
      where ( scloud > lsize ) local = .true.
    else if (lsize < 0.) then
      suffix = 'small'
      where ( scloud < abs(lsize) ) local = .true.
    endif
!
    pdf1 = 0.; pdf2 = 0.; pdf3 = 0.
    count = 0
    labels = 0
    npts = 0
    store = 0.
    prop = 0.
    do h = 2, nh
      store(h,1) = store(h-1,1) + dh
    enddo
!
    do j = 2, NZO
!
      c = 0
      do m = 1, nslic
        do t = 1, ntim
          do i = 1, NYO
!
            if ( local(i,j,m,t) .and. mindist(i,j,m,t) == -dx .and. z(j) <= top_max ) then 
	      if ( .not.any(labels == int(ncloud(i,j,m,t))) ) then
	        count = count + 1
		labels(count) = int(ncloud(i,j,m,t))
	      endif	    
!
              if (norm) then
	        height = (z(j) - bcloud(i,j,m,t)) / abs(scloud(i,j,m,t))
	      else
		height = z(j)
	      endif
!
	      deltaw = 0.
	      deltab = 0.
	      minw = 0.
	      minb = 0.
	      add = .true.
	      if ( mindist(min(i+1,NYO),j,m,t) >= 0. ) then
	        do l = 1, 100
      	          i1 = i+l
      	          if ( i+l > NYO ) i1 = i1 - NYO
		  if ( slice(i1,j,3,m,t) < 0. .and. mindist(i1,j,m,t) > 0. .and. ncloud(i1,j,m,t) == ncloud(i,j,m,t) ) then
		    deltaw = deltaw + dx
		    if ( slice(i1,j,10,m,t) < 0. .and. add ) then
		      deltab = deltab + dx
		    else
		      add = .false.
		    endif
		    minw = min(minw,slice(i1,j,3,m,t))
		    minb = min(minb,slice(i1,j,10,m,t))		  
		  else
		    exit
		  endif
		enddo
	      else 
	        do l = 1, 100
      	          i1 = i-l
      	          if ( i-l < 1 ) i1 = i1 + NYO
		  if ( slice(i1,j,3,m,t) < 0. .and. mindist(i1,j,m,t) > 0. .and. ncloud(i1,j,m,t) == ncloud(i,j,m,t) ) then
		    deltaw = deltaw + dx
		    if ( slice(i1,j,10,m,t) < 0. .and. add ) then
		      deltab = deltab + dx
		    else
		      add = .false.
		    endif
		    minw = min(minw,slice(i1,j,3,m,t))
		    minb = min(minb,slice(i1,j,10,m,t))
		  else
		    exit
		  endif
		enddo	      
	      endif
!
	      do h = 1, nh
	        if ( height >= dh*real(h-1) .and. height < dh*real(h) ) exit
	      enddo
!
	      if ( minw <= 0. .and. deltaw < 3000. ) then
  	        npts(h) = npts(h) + 1
	        prop(h,npts(h),1) = minw
	        prop(h,npts(h),2) = minb
	        prop(h,npts(h),3) = deltaw
	        prop(h,npts(h),4) = deltab
	        store(h,2) = store(h,2) + minw
	        store(h,3) = store(h,3) + minb
	        store(h,4) = store(h,4) + deltaw
	        store(h,5) = store(h,5) + deltab
		if (z(j-1) < 2000. .and. z(j) >= 2000.) then
		  c = c + 1
		  pdf1(c,:) = (/minw, minb, deltaw, deltab/)
		endif
		if (z(j-1) < 4000. .and. z(j) >= 4000.) then
		  c = c + 1
		  pdf2(c,:) = (/minw, minb, deltaw, deltab/)
		endif
		if (z(j-1) < 6000. .and. z(j) >= 6000.) then
		  c = c + 1
		  pdf3(c,:) = (/minw, minb, deltaw, deltab/)
		endif
              endif
	    endif
!
	  enddo
        enddo
      enddo
!
    enddo
!
    do j = 2, 5
      where ( npts > 0 ) store(:,j) = store(:,j) / real(npts)
    enddo
!
    print*, 'Sampled', count, 'in ', trim(suffix), ' analysis'
!
!  Open file and write output
!
    inquire (FILE=trim(file_output), OPENED=ex)
    if (.not.ex) then
      open(200,file=file_output(1:len_trim(file_output))//'_'//trim(suffix),	&
      		form='formatted', status='unknown')
    endif
!
    do j = 1, nh
      if ( npts(j) > 1 ) then    
        call sort(prop(j,:,1),20000)
        call sort(prop(j,:,2),20000)
        call sort(prop(j,:,3),20000)
        call sort(prop(j,:,4),20000)
!
        start=20000-npts(j)+1
        write(200,*) store(j,1), store(j,2), store(j,3), store(j,4), store(j,5), prop(j,floor(0.95*real(npts(j)))+1,1), prop(j,ceiling(0.05*real(npts(j))),1), prop(j,floor(0.95*real(npts(j)))+1,2), prop(j,ceiling(0.05*real(npts(j))),2), prop(j,start+floor(0.95*real(npts(j))),3), prop(j,start+ceiling(0.05*real(npts(j)))-1,3), prop(j,start+floor(0.95*real(npts(j))),4), prop(j,start+ceiling(0.05*real(npts(j)))-1,4)
      endif
    enddo
!
!  Output PDFs
!
    open(201,file='pdf_shell_2000.dat', form='formatted', status='unknown')
    do j = 1, 4000
      if (pdf1(j,1) < 0.) write(201,*) pdf1(j,:)
    enddo
    close(201)
!
    open(202,file='pdf_shell_4000.dat', form='formatted', status='unknown')
    do j = 1, 4000
      if (pdf2(j,1) < 0.) write(202,*) pdf2(j,:)
    enddo
    close(202)
!
    open(203,file='pdf_shell_6000.dat', form='formatted', status='unknown')
    do j = 1, 4000
      if (pdf3(j,1) < 0.) write(203,*) pdf3(j,:)
    enddo
    close(203)
!
    deallocate( store, prop, npts )
!
    end subroutine write_mean
!
    subroutine write_corr ( lsize, z, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud )
    
    integer :: nslic,ntim
    character(len = 100) :: file_output
    real :: thres, lsize
    real, dimension(NZO) :: z
    real, dimension(NYO,NZO,nvarslic,nslic,ntim) :: slice
    real, dimension(NYO,NZO,nslic,ntim) :: mindist, ncloud, scloud, bcloud, wcloud
!    
    logical :: ex
    real :: height, deltaw
    character(len = 20) :: suffix
    integer :: c, i, j, l, t, m, nc, maxnc
    real, allocatable, dimension(:) :: minw, minb, mindp, minbuo, mindpre, mina, posw, posb, posdp, depth
    logical, dimension(NYO,NZO,nslic,ntim) :: local
!
    maxnc = maxval(ncloud)
!
    allocate( minw(1:maxnc), minb(1:maxnc), mindp(1:maxnc), minbuo(1:maxnc), mindpre(1:maxnc), mina(1:maxnc),  &
    	posw(1:maxnc), posb(1:maxnc), posdp(1:maxnc), depth(1:maxnc) )
!
    local = .false.
    if (lsize == 0.) then
      suffix = 'tot'
      local = .true.
    else if (lsize == 6.) then
      suffix = 'activ'
      where ( scloud >= 3000. .and. wcloud > lsize ) local = .true.
    else if (lsize > 0.) then
      suffix = 'large'
      where ( scloud > lsize ) local = .true.
    else if (lsize < 0.) then
      suffix = 'small'
      where ( scloud < abs(lsize) ) local = .true.
    endif
!
    minw = 0.
    minb = 0.
    mindp = 0.
    minbuo = 0.
    mindpre = 0.
    mina = 0.
    posw = 0.
    posb = 0.
    posdp = 0.
    depth = 0.
!
    do j = 2, NZO
!
      height = z(j)
      do m = 1, nslic
        do t = 1, ntim
          do i = 1, NYO
!
            if ( local(i,j,m,t) .and. mindist(i,j,m,t) >= 0. .and. mindist(i,j,m,t) <= 2500. .and. z(j) <= top_max ) then 
!
	      deltaw = 0.
	      nc = int(ncloud(i,j,m,t))
	      if ( mindist(min(i+1,NYO),j,m,t) > 0. ) then
	        do l = 1, 50
      	          i1 = i+l
      	          if ( i+l > NYO ) i1 = i1 - NYO
		  if ( slice(i1,j,3,m,t) < 0. .and. mindist(i1,j,m,t) > 0. .and. ncloud(i1,j,m,t) == ncloud(i,j,m,t) ) then
	      	    minw(nc) = min(minw(nc),slice(i1,j,3,m,t))
	      	    minbuo(nc) = min(minbuo(nc),slice(i1,j,19,m,t))
	      	    mindpre(nc) = min(mindpre(nc),slice(i1,j,20,m,t))
	      	    if ( minw(nc) == slice(i1,j,3,m,t) ) then
			posw(nc) = height
	        	mina(nc) = slice(i1,j,18,m,t)
	        	minb(nc) = slice(i1,j,19,m,t)
	       		mindp(nc) = slice(i1,j,20,m,t)
			depth(nc) = bcloud(i1,j,m,t) + scloud(i1,j,m,t)
	      	    endif
	      	    if ( minbuo(nc) == slice(i1,j,19,m,t) ) then
	        	posb(nc) = height
	      	    endif
	      	    if ( mindpre(nc) == slice(i1,j,20,m,t) ) then
	        	posdp(nc) = height
	      	    endif
		  else
		    exit
		  endif
		enddo
	      else 
	        do l = 1, 50
      	          i1 = i-l
      	          if ( i-l < 1 ) i1 = i1 + NYO
		  if ( slice(i1,j,3,m,t) < 0. .and. mindist(i1,j,m,t) > 0. .and. ncloud(i1,j,m,t) == ncloud(i,j,m,t) ) then
	      	    minw(nc) = min(minw(nc),slice(i1,j,3,m,t))
	      	    minbuo(nc) = min(minbuo(nc),slice(i1,j,19,m,t))
	      	    mindpre(nc) = min(mindpre(nc),slice(i1,j,20,m,t))
	      	    if ( minw(nc) == slice(i1,j,3,m,t) ) then
			posw(nc) = height
	        	mina(nc) = slice(i1,j,18,m,t)
	        	minb(nc) = slice(i1,j,19,m,t)
	       		mindp(nc) = slice(i1,j,20,m,t)
			depth(nc) = bcloud(i1,j,m,t) + scloud(i1,j,m,t)
	      	    endif
	      	    if ( minbuo(nc) == slice(i1,j,19,m,t) ) then
	        	posb(nc) = height
	      	    endif
	      	    if ( mindpre(nc) == slice(i1,j,20,m,t) ) then
	        	posdp(nc) = height
	      	    endif
		  else
		    exit
		  endif
		enddo	      
	      endif
!
	    endif
!
	  enddo
        enddo
      enddo
!
    enddo
!
    c = 0
    open(201,file='corr_shell_all.dat_'//suffix, form='formatted', status='unknown')
    do j = 1, maxnc
      if (posw(j) > 0.) then
        c = c + 1
        write(201,*) c, depth(j), posw(j), posb(j), posdp(j), minw(j), minb(j), mindp(j), mina(j)
      endif
    enddo
    close(201)
!
    deallocate( posw, posb, posdp, depth, minw, minb, mindp, minbuo, mindpre, mina )
!
    end subroutine write_corr
!
    subroutine write_distributions (lsize, z, index, thres, file_output, nslic, ntim, slice, mindist, ncloud, scloud, bcloud, wcloud)
    
    character(len = 100) :: file_output
    integer :: index(12),nslic,ntim,ny
    real :: thres, lsize
    real, dimension(NZO) :: z
    real, dimension(NYO,NZO,nvarslic,nslic,ntim) :: slice
    real, dimension(NYO,NZO,nslic,ntim) :: mindist, ncloud, scloud, bcloud, wcloud
    real, dimension(NYO,NZO,nslic,ntim) :: ncloudnew
    logical, dimension(NYO,NZO,nslic,ntim) :: local
!    
    logical :: ex
    character(len = 1) :: car1
    character(len = 2) :: car2
    character(len = 20) :: suffix
    integer :: i, j, k, t, l, m, n, h, nh, unit1, count, countall, labels(4000)
    integer, dimension(NZO) :: nn
    real, dimension(NZO) :: mean
    real, dimension(NYO,NZO,nslic,ntim) :: tcloud
    real, dimension(:,:,:), allocatable :: store
    real :: height, dist, dh, htot
!
    if (norm) then
      dh = 0.05
      htot = 1.
    else
      dh = 200.
      htot = 10000.
    endif
    nh = 1+int(htot/dh)
    ny = int(2.*dist_max/dx) + 1
    allocate( store(ny,nh,4) )
!
    local = .false.
    if (lsize == 0.) then
      suffix = 'tot'
      local = .true.
    else if (lsize == 5.) then
      suffix = 'activ'
      where ( scloud >= 2200. .and. wcloud > lsize ) local = .true.
    else if (lsize > 0.) then
      suffix = 'large'
      where ( scloud > lsize ) local = .true.
    else if (lsize < 0.) then
      suffix = 'small'
      where ( scloud < abs(lsize) ) local = .true.
    endif
!
    do n = 1, size(index)
!
    nn = 0
    mean = 0.
    countall = 0
    labels = 0
    ncloudnew = 0.
    do j = 1, NZO
      do m = 1, nslic
        do t = 1, ntim
          do i = 1, NYO
            if ( slice(i,j,21,m,t) < thres .and. mindist(i,j,m,t) > 0. .and. mindist(i,j,m,t) <= dist_max ) then
	      if ( .not.any(labels == int(slice(i,j,21,m,t))) ) then
	        countall = countall + 1
		labels(countall) = int(slice(i,j,21,m,t))
	      endif	    
	      nn(j) = nn(j) + 1
	      mean(j) = mean(j) + slice(i,j,index(n),m,t)
	      ncloudnew(i,j,m,t) = 1.
	    endif
          enddo
        enddo
      enddo
      if (nn(j) > 0) mean(j) = mean(j) / real(nn(j))
    enddo
!
    count = 0
    labels = 0
    store = 0.
    store(:,:,4) = 1.
    store(1,:,1) = -dist_max
    do h = 2, nh
      store(:,h,2) = store(:,h-1,2) + dh
    enddo
    do k = 2, ny
      store(k,:,1) = store(k-1,:,1) + dx
    enddo
!
    do j = 1, NZO
      do m = 1, nslic
        do t = 1, ntim
          do i = 1, NYO
!
            if ( mindist(i,j,m,t) > -dist_max .and. mindist(i,j,m,t) <= dist_max .and. local(i,j,m,t) .and. z(j) <= top_max ) then 
	      if ( .not.any(labels == int(ncloud(i,j,m,t))) ) then
	        count = count + 1
		labels(count) = int(ncloud(i,j,m,t))
	      endif	    
!
              if (norm) then
	        height = (z(j) - bcloud(i,j,m,t)) / abs(scloud(i,j,m,t))
	      else
		height = z(j)
	      endif
!
	      do h = 1, nh
	        if ( height >= dh*real(h-1) .and. height < dh*real(h) ) exit
	      enddo
!
	      dist = mindist(i,j,m,t)
	      do k = 0, ny-1
	        if ( dist >= dx*real(k-ny/2) .and. dist < dx*real(k-ny/2+1) ) exit
	      enddo
	      if (dist > 0.) k = k-1
	      if (k == ny) k = ny
!
  	      store(k,h,3) = store(k,h,3) + slice(i,j,index(n),m,t)
	      store(k,h,4) = store(k,h,4) + 1.
            endif
!
	  enddo
        enddo
      enddo
    enddo
!
    store(:,:,3) = store(:,:,3) / store(:,:,4)
!
    print*, 'Sampled', count, 'in ', trim(suffix), ' analysis'
!
!  Open file and write output
!
    unit1=201
    if (n < 10) then
      write(car1,'(i1)') n
      open(unit1,file=file_output(1:len_trim(file_output))//'_'//car1//'_'//trim(suffix),	&
      		form='formatted', status='unknown')
    else
      write(car2,'(i2)') n
      open(unit1,file=file_output(1:len_trim(file_output))//'_'//car2//'_'//trim(suffix),	&
      		form='formatted', status='unknown')
    endif
!
    do j = 1, nh
      do i = 1, ny
        if (store(i,j,1) >= -dist_max .and. store(i,j,1) <= dist_max) write(unit1,*) store(i,j,1), store(i,j,2), store(i,j,3), store(i,j,4)
      enddo
    enddo
!
    close (unit1)
!
    enddo
!
    deallocate( store )
!
    end subroutine write_distributions
!  
  END SUBROUTINE distance_stats

  SUBROUTINE distributions ( thres, nslic, ntim, slice )
    
    integer, intent(in) :: nslic, ntim
    real, intent(in) :: thres
    real, dimension(NYO,NZO,nvarslic,nslic,ntim) :: slice
!    
    real :: size(10000), top(10000), wmax(10000), labels(10000)
    integer :: k, i, t, m, j, count, unit1, unit2, unit3
!
    count = 0
    labels = 0.
    size = 0.
    do j = 1, NZO
      do m = 1, nslic
        do t = 1, ntim
          do i = 1, NYO
	    if ( slice(i,j,22,m,t) >= 500. .and. slice(i,j,23,m,t) < 2500. .and. slice(i,j,23,m,t) > 750.  .and. .not.any(slice(i,j,21,m,t) == labels(1:max(count,1))) ) then
	      count = count + 1
	      labels(count) = slice(i,j,21,m,t)
	      size(count) = slice(i,j,22,m,t)
	      top(count) = slice(i,j,23,m,t)+slice(i,j,22,m,t)
	      wmax(count) = slice(i,j,24,m,t)
	    endif
	  enddo
	enddo
      enddo
    enddo
!
    call sort(size,10000)
    call sort(top,10000)
    call sort(wmax,10000)
!
!  Open file and write output
!
    unit1=201
    unit2=202
    unit3=203
    open(unit1,file='distribution_size.dat', form='formatted', status='unknown')
    open(unit2,file='distribution_top.dat', form='formatted', status='unknown')
    open(unit3,file='distribution_wmax.dat', form='formatted', status='unknown')
!
    do k = 10000-count, 10000
      write(unit1,*) size(k), k-10000
      write(unit2,*) top(k), k-10000
      write(unit3,*) wmax(k), k-10000
    enddo
!
    close (unit1)
    close (unit2)
    close (unit3)
!  
    return
  END SUBROUTINE distributions
!
  SUBROUTINE thermals ( z, index, thres, nslic, ntim, slice, file_output )
    
    character(len = 40), intent(in) :: file_output
    integer, intent(in) :: index(12), ntim, nslic
    real, intent(in) :: thres, z(:)
    real, dimension(NYO,NZO,nvarslic,nslic,ntim), intent(inout) :: slice
!    
    real :: bmax(4000), xx(81), zz(21), zmid
    integer :: i, j, k, t, m, n, i1, j1, jbase, jtop
    integer :: unit1, count, count2, labels(4000), labels2(4000), bmloc(4000,2), ind(4000,2)
    real, dimension(81,31,size(index)+1) :: store
!
    zmid = 4750.
    jbase = 0
    jtop = 0
    do j = 1, NZO
      if ( z(j) >= zmid-750. .and. jbase == 0 ) jbase = j
      if ( z(j) >  zmid+750. .and. jtop  == 0 ) jtop  = j-1
    enddo
!
    do i = 0, 80
      xx(i+1) = dx*real(i-40)
    enddo
    do j = 0, 20
      zz(j+1) = 200.*real(j-10)
    enddo
!
!  Limit the analysis to a layer between 4000 and 6000 m
!
    count = 0
    labels = 0
    bmax = 0.
!
    do n = 1, nslic
    do t = 1, ntim
    do j = jbase, jtop
      do i = 1, NYO
!    
        if ( slice(i,j,13,n,t) > thres ) then
	  if ( .not.any(labels == int(slice(i,j,21,n,t))) ) then
	    count = count + 1
	    k = count
	    ind(k,:) = (/n,t/)
	    labels(k) = int(slice(i,j,21,n,t))
	    bmax(k) = max(bmax(k), slice(i,j,3,n,t))
	    if (bmax(k) == slice(i,j,3,n,t)) bmloc(k,:) = (/i,j/)
	  else
	    do k = 1, count
	      if ( labels(k) == int(slice(i,j,21,n,t)) ) then
	        bmax(k) = max(bmax(k),slice(i,j,3,n,t))
	        if (bmax(k) == slice(i,j,3,n,t)) bmloc(k,:) = (/i,j/)
		exit	      
	      endif
	    enddo
	  endif
	endif
!  
      enddo
    enddo
    enddo
    enddo
!   
    store = 0.
    count2 = 0
    labels2 = 0
    do n = 1, count
      if ( bmax(n) > 5. ) then
        do j = 0, 30
          do i = 0, 80
      	     i1 = bmloc(n,1) + real(i-40)
	     j1 = bmloc(n,2) + real(j-15) 
      	     if ( i1 < 1 ) i1 = i1 + NYO
      	     if ( i1 > NYO ) i1 = i1 - NYO
	     do k = 1, 20
	       if ( z(j1) >= zmid + zz(k) .and. z(j1) < zmid + zz(k+1) ) exit
	     enddo
	     if ( .not.any(labels2 == labels(n)) ) then
	       count2 = count2 + 1
	       labels2(count2) = labels(n)
	     endif
!
	     store(i+1,k,1) = store(i+1,k,1) + 1.
    	     do m = 1, size(index)
	       store(i+1,k,1+m) = store(i+1,k,1+m) + slice(i1,j1,index(m),ind(n,1),ind(n,2))
    	     enddo
!
	  enddo
	enddo
      endif
    enddo
!
    print*, 'Number of thermals identified:', count2
!
    open(201,file='thermals_rh80_cloud_drh.dat', form='formatted', status='unknown')
!
    do j = 1, 21
      do i = 1, 81
        if ( store(i,j,1) > 1. ) then
	  write(201,*) xx(i), zz(j), store(i,j,2:)/store(i,j,1), store(i,j,1)
	else
	  write(201,*) xx(i), zz(j), store(i,j,2:), store(i,j,1)
	endif
      enddo
    enddo
!
    close(201)
!    
    return
  end subroutine thermals
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

