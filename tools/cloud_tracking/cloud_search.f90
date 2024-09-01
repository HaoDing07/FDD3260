  
  module searchmod

  USE cloudsmod
    
  public :: object_search, subobject_search, object_split, find_halos, time_tracking, cloud_diagnostics,  &
  	    cloud_stats, cloud_output
  
  CONTAINS

!      =================================================                            
  subroutine object_search ( ifdrop, id_cl, r_cl, clouds )
!      =================================================                            
!
  IMPLICIT NONE
!
  logical :: lmerge
  integer :: icloud, i, j, k, l, m, n
!     
  type (cloud_type), dimension(ncl) :: clouds
!
  integer, dimension (NXO,NYO) :: ifdrop
  real, dimension (NXO,NYO) :: id_cl, r_cl
!
!----------------------------------------------------------!
!	       Cloud search on master proc	      	   !
!----------------------------------------------------------!
!
  nclouds=0
  do n = 1, ncl
    clouds(n)%npts = 0
    clouds(n)%size = 0
    clouds(n)%id   = 0
    clouds(n)%i(:) = 0
    clouds(n)%j(:) = 0
    clouds(n)%ic   = 0
    clouds(n)%jc   = 0
!
    clouds(n)%mfl   = 0.0
    clouds(n)%dist  = 0.0
  enddo
!
!  Starting cloud search
!
  do j = min(js+1,je), je
    do i = is+1, ie
      call search ( ncl, i, j, nclouds, ifdrop, ifdrop, clouds )
    enddo
  enddo
!
  call qnbc_clouds ( ncl, nclouds, ifdrop, ifdrop, clouds )
!
!----------------------------------------------------------!
!                  Merge and filter clouds     	  	   !
!----------------------------------------------------------!
!
!  Merge clouds
!
  call object_merge ( ncl, nclouds, clouds )
!
!  Remove small clouds
!
  do n = 1, ncl
    if ( clouds(n)%size < nc_min .and. clouds(n)%size > 0 ) then  
      clouds(n)%i(:) = 0
      clouds(n)%j(:) = 0
      clouds(n)%id   = 0
      clouds(n)%npts = 0
      clouds(n)%size = 0
    endif
  enddo 
!
!  Recount clouds and reset IDs
!
  m=0
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
      m=m+1
      clouds(n)%id = m
    endif
  enddo
  nclouds=m
!
!  Assign IDs in 2D
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
      do l = 1, clouds(n)%npts
        id_cl(clouds(n)%i(l),clouds(n)%j(l)) = clouds(n)%id
        r_cl(clouds(n)%i(l),clouds(n)%j(l)) = sqrt(real(clouds(n)%size)*dxl*dyl)
      enddo
    endif
  enddo
!
  call qnbc_id ( id_cl )
  call qnbc_id ( r_cl )
!
!----------------------------------------------------------!
!
!  print*, "Terminating cloud_search"
!
  return
  
  end subroutine object_search
!
!      =================================================                            
  subroutine subobject_search ( ifcl, ifup, id_up, clouds, updrafts )
!      =================================================                            
!
  IMPLICIT NONE
!
  integer :: i, j, k, m, n
!     
  type(cloud_type), dimension(ncl) :: clouds
  type(cloud_type), dimension(ncl,nup) :: updrafts
!
  integer, dimension (NXO,NYO) :: ifcl, ifup
  real, dimension (NXO,NYO) :: id_up
!
!----------------------------------------------------------!
!	        Updraft search on master proc	      	   !
!----------------------------------------------------------!
!
  do n = 1, ncl
    nupdrafts(n)=0
    do m = 1, nup
      updrafts(n,m)%npts = 0
      updrafts(n,m)%size = 0
      updrafts(n,m)%id   = 0
      updrafts(n,m)%i(:) = 0
      updrafts(n,m)%j(:) = 0
      updrafts(n,m)%ic   = 0.0
      updrafts(n,m)%jc   = 0.0
!
      updrafts(n,m)%mfl  = 0.0
      updrafts(n,m)%dist = 0.0
    enddo
  enddo
!
!  Starting updraft search
!
  do n = 1, ncl
    if (clouds(n)%size >= nc_min) then
!
      do k = 1, clouds(n)%npts
        i = clouds(n)%i(k)
        j = clouds(n)%j(k)
!    
	call search ( nup, i, j, nupdrafts(n), ifup, ifcl, updrafts(n,:) )
      enddo
    endif
!
    call qnbc_clouds ( nup, nupdrafts(n), ifup, ifcl, updrafts(n,:) )
  enddo
!
!----------------------------------------------------------!
!                Merge and filter updrafts     	  	   !
!----------------------------------------------------------!
!
!  Merge detected updrafts
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) call object_merge ( nup, nupdrafts(n), updrafts(n,:) )
  enddo
!
!  Remove small updrafts
!
  do n = 1, ncl
    do m = 1, nup 
      if ( updrafts(n,m)%size < nu_min .and. updrafts(n,m)%size > 0 ) then  
    	updrafts(n,m)%i(:) = 0
    	updrafts(n,m)%j(:) = 0
    	updrafts(n,m)%id   = 0
    	updrafts(n,m)%npts = 0
    	updrafts(n,m)%size = 0
      endif
    enddo 
  enddo
!
!  Recount updrafts
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
      m=0
      do k = 1, nup
        if (updrafts(n,k)%size >= nu_min) m=m+1
      enddo
      nupdrafts(n)=m
    endif
  enddo
!
!  Assign IDs in 2D
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
      do m = 1, nup
        if ( updrafts(n,m)%size >= nu_min ) then
          do k = 1, updrafts(n,m)%npts
	    id_up(updrafts(n,m)%i(k),updrafts(n,m)%j(k)) = updrafts(n,m)%id
          enddo
	endif
      enddo
    endif
  enddo
!
  call qnbc_id ( id_up )
!
!----------------------------------------------------------!
!
!  print*, "Terminating updraft_search"
!
  return
!
  end subroutine subobject_search
!
!      =================================================                            
  subroutine search ( ntot, i, j, npts, if1, if2, clouds )
!      =================================================                            
!
  IMPLICIT NONE
!  
  integer :: i, j, npts, ntot
  integer, dimension (NXO,NYO) :: if1, if2
  type (cloud_type), dimension (ntot) :: clouds
!
  integer :: m, l, ii
  logical :: ex
!
!----------------------------------------------------------!
!
  if ( if1(i,j) > 0 ) then
!  
    ex=.false.
    if ( if1(i-1,j) == 1 .and. if2(i-1,j) == 1 ) then
      do m = 1, npts
  	if ( any(i-1 == clouds(m)%i) .and. any(j == clouds(m)%j) ) then
  	  clouds(m)%npts=clouds(m)%npts+1
  	  ii=m
  	  ex=.true.
  	  exit 
  	endif
      enddo
    else if ( NYO > 1 ) then
      if (if1(i,j-1) == 1 .and. if2(i,j-1) == 1 ) then
        do m = 1, npts
  	  if ( any(i == clouds(m)%i) .and. any(j-1 == clouds(m)%j) ) then
  	    clouds(m)%npts=clouds(m)%npts+1
  	    ii=m
  	    ex=.true.
  	    exit
  	  endif
        enddo
      endif
    endif    
!
    if ( .not.ex ) then
      npts=npts+1
      clouds(npts)%npts=1
      clouds(npts)%id=npts
      ii=npts
    endif
!
    clouds(ii)%i(clouds(ii)%npts)=i
    clouds(ii)%j(clouds(ii)%npts)=j
    clouds(ii)%size=clouds(ii)%npts
!
  endif 
!
!----------------------------------------------------------!
!  
  return
  end subroutine search
!
!      =================================================                            
  subroutine object_split ( id_cl, clouds, updrafts )
!      =================================================                            
!
  IMPLICIT NONE
!     
  type(cloud_type), dimension(ncl) :: clouds
  type(cloud_type), dimension(ncl,nup) :: updrafts
  real, dimension (NXO,NYO) :: id_cl
!
  integer :: i, j, k, l, m, n
  integer :: ncloud, nscl, nspl, minup
  real :: dist, mindist
!
  integer, dimension(nup) :: ic, jc
  logical, dimension(200) :: lsplit
  type(cloud_type), dimension(nup) :: split_clouds
!
!----------------------------------------------------------!
!	    Split clouds with multiple updrafts	      	   !
!----------------------------------------------------------!
!
!  Starting cloud splits search
!
  ncloud = 0
  do n = 1, ncl
!
    nscl = 0
    do m = 1, nup
      split_clouds(m)%mfl   = 0.
      split_clouds(m)%dist  = 0.
      split_clouds(m)%angle = 0.
      split_clouds(m)%npts  = 0
      split_clouds(m)%size  = 0
      split_clouds(m)%stat  = 0
      split_clouds(m)%id    = 0
      split_clouds(m)%i     = 0
      split_clouds(m)%j     = 0
      split_clouds(m)%ic    = 0
      split_clouds(m)%jc    = 0
!
      if (clouds(n)%npts >= nc_min .and. updrafts(n,m)%npts >= nu_min) nscl = nscl + 1
!
      call find_center (updrafts(n,m)%npts, updrafts(n,m)%i, updrafts(n,m)%j, ic(m), jc(m))
    enddo
!
!  Cloud status
!
    if (nscl == 0) then
      clouds(n)%stat = 0
    else if (nscl == 1) then
      clouds(n)%stat = 1
    else
      clouds(n)%stat = 2
    endif
!
    if (clouds(n)%npts >= nc_min) then
      split_clouds(1:nup)%npts = updrafts(n,1:nup)%npts
      split_clouds(1:nup)%size = updrafts(n,1:nup)%size
      do k = 1, clouds(n)%npts
        split_clouds(1:nup)%i(k) = updrafts(n,1:nup)%i(k)
        split_clouds(1:nup)%j(k) = updrafts(n,1:nup)%j(k)
      enddo
!
      lsplit = .false.
      nspl = 0
      do k = 1, clouds(n)%npts
        if ( clouds(n)%i(k) > 0 .and. clouds(n)%j(k) > 0 ) then
	  lsplit(k) = .true.
          nspl = nspl + 1
	endif
!
	do m = 1, nup
          if ( split_clouds(m)%npts >= nu_min ) then
            do l = 1, split_clouds(m)%npts
	      if ( clouds(n)%i(k) == split_clouds(m)%i(l) .and. clouds(n)%j(k) == split_clouds(m)%j(l) ) then
	        lsplit(k) = .false.
                nspl = nspl - 1
	      endif
	    enddo
	  endif
	  if ( .not.lsplit(k) ) exit
	enddo
      enddo
!
      do k = 1, clouds(n)%npts
	if ( lsplit(k) ) then
	  mindist = 1.E+12
	  minup = 0
!
          do m = 1, nup
	    if ( updrafts(n,m)%npts >= nu_min ) then
	      if ( clouds(n)%i(k)-NXO/2 > is ) then
   	        dist = (ic(m) - (clouds(n)%i(k) - NXO + is))**2.
	      else
	        dist = (ic(m) - clouds(n)%i(k))**2.
	      endif
	      if ( clouds(n)%j(k)-NYO/2 > js ) then
   	        dist = dist + (jc(m) - (clouds(n)%j(k) - NYO + js))**2.
              else
	        dist = dist + (jc(m) - clouds(n)%j(k))**2.
	      endif	      
	      dist = sqrt(dist)
!
	      if (dist < mindist) then
	        mindist = dist
		minup = m
	      endif
	    endif
	  enddo
!
	  if (minup > 0) then
	    split_clouds(minup)%npts = split_clouds(minup)%npts + 1
	    split_clouds(minup)%size = split_clouds(minup)%size + 1
	    split_clouds(minup)%i(split_clouds(minup)%npts) = clouds(n)%i(k)
	    split_clouds(minup)%j(split_clouds(minup)%npts) = clouds(n)%j(k)
	  endif
	endif
      enddo
!
!  Update clouds ID and properties
!
      do m = 1, nup
        if (split_clouds(m)%npts >= nu_min) then
	  ncloud = ncloud + 1
    	  clouds(ncloud)%i = split_clouds(m)%i
    	  clouds(ncloud)%j = split_clouds(m)%j
    	  clouds(ncloud)%npts = split_clouds(m)%npts
    	  clouds(ncloud)%size = split_clouds(m)%size
    	  clouds(ncloud)%id = ncloud
	  clouds(ncloud)%stat = 1
	  do l = 1, split_clouds(m)%npts
	    id_cl(split_clouds(m)%i(l),split_clouds(m)%j(l)) = ncloud	
	  enddo
	endif
      enddo
!
    endif
  enddo
!
!----------------------------------------------------------!
!
!  print*, "Terminating cloud_split"
!
  return
!
  end subroutine object_split
!
!      =================================================                            
  subroutine object_merge ( ntot, npts, clouds )  
!      =================================================                            
!
  IMPLICIT NONE
!  
  integer :: ntot, npts
  type (cloud_type), dimension(ntot) :: clouds
!
  integer :: n, m, k, l, j
  logical :: lmerge
!
!----------------------------------------------------------!
!
  n_loop: do n = 1, npts
    m = 1
    m_loop: do
      lmerge=.false.
!
!  Searching for contiguous clouds (including through boundaries)
!
      if ( m /= n .and. clouds(m)%size /= 0 ) then
  	do l = 1, clouds(n)%npts
	  do j = 1, clouds(m)%npts
  	    if ( (clouds(n)%i(l)+1 == clouds(m)%i(j) .and. clouds(n)%j(l) == clouds(m)%j(j)) .or.       &
	         (clouds(n)%i(l)-1 == clouds(m)%i(j) .and. clouds(n)%j(l) == clouds(m)%j(j)) .or.	&
  	         (clouds(n)%j(l)+1 == clouds(m)%j(j) .and. clouds(n)%i(l) == clouds(m)%i(j)) .or.	&
  	         (clouds(n)%j(l)-1 == clouds(m)%j(j) .and. clouds(n)%i(l) == clouds(m)%i(j)) ) then   
  	      lmerge=.true.
  	      exit
  	    endif
	  enddo
  	enddo
!
!  Merging 2 clouds
!
  	if (lmerge) then
  	  clouds(n)%i(clouds(n)%npts+1:clouds(n)%npts+clouds(m)%npts)=clouds(m)%i(1:clouds(m)%npts)
  	  clouds(n)%j(clouds(n)%npts+1:clouds(n)%npts+clouds(m)%npts)=clouds(m)%j(1:clouds(m)%npts)
  	  clouds(n)%npts=clouds(n)%npts+clouds(m)%npts
  	  clouds(n)%size=clouds(n)%size+clouds(m)%size
!
    	  clouds(m)%i(:) = 0
    	  clouds(m)%j(:) = 0
    	  clouds(m)%npts = 0
    	  clouds(m)%size = 0
    	  clouds(m)%id   = 0
	  m = 1
  	endif
      endif
!
      m = m+1
      if (m > npts) exit m_loop
!
    enddo m_loop
  enddo n_loop
!
!----------------------------------------------------------!
!  
  return
  end subroutine object_merge
!
!      =================================================                            
  subroutine find_halos ( id_cl, clouds )  
!      =================================================                            
!
  IMPLICIT NONE
!  
  real, dimension (NXO,NYO) :: id_cl
  type (cloud_type), dimension(ncl) :: clouds
!
  logical :: lres
  integer :: n, k, l, i, j, ip, im, jp, jm
  real, allocatable, dimension (:,:) :: id_ha
!
!  print*, "Starting find_halos"
!
!----------------------------------------------------------!
!
  allocate( id_ha(1:NXO,1:NYO) )
!
  do n = 1, ncl
    do l = 1, 6
      clouds(n)%halo(l)%i = 0
      clouds(n)%halo(l)%j = 0
      clouds(n)%halo(l)%npts = 0
    enddo
  enddo
!
  do n = 1, ncl
    id_ha = 0.
    if ( clouds(n)%npts >= nc_min ) then
!
!  First layer
!
    where ( id_cl == clouds(n)%id ) id_ha = clouds(n)%id
!
    do k = 1, clouds(n)%npts
      i = clouds(n)%i(k)
      j = clouds(n)%j(k)
      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      if (ip > NXO) ip = ip - NXO
      if (im < 1) im = im + NXO
      if (jp > NYO) jp = jp - NYO
      if (jm < 1) jm = jm + NYO
!
100   continue
      lres = .false.
      if ( id_ha(ip,j) /= clouds(n)%id ) then
        clouds(n)%halo(1)%npts = clouds(n)%halo(1)%npts + 1
	clouds(n)%halo(1)%i(clouds(n)%halo(1)%npts) = ip
	clouds(n)%halo(1)%j(clouds(n)%halo(1)%npts) = j
	id_ha(ip,j) = clouds(n)%id
	lres = .true.
      else if ( id_ha(i,jp) /= clouds(n)%id ) then
        clouds(n)%halo(1)%npts = clouds(n)%halo(1)%npts + 1
	clouds(n)%halo(1)%i(clouds(n)%halo(1)%npts) = i
	clouds(n)%halo(1)%j(clouds(n)%halo(1)%npts) = jp
	id_ha(i,jp) = clouds(n)%id
	lres = .true.
      else if ( id_ha(im,j) /= clouds(n)%id ) then
        clouds(n)%halo(1)%npts = clouds(n)%halo(1)%npts + 1
	clouds(n)%halo(1)%i(clouds(n)%halo(1)%npts) = im
	clouds(n)%halo(1)%j(clouds(n)%halo(1)%npts) = j
	id_ha(im,j) = clouds(n)%id
	lres = .true.
      else if ( id_ha(i,jm) /= clouds(n)%id ) then
        clouds(n)%halo(1)%npts = clouds(n)%halo(1)%npts + 1
	clouds(n)%halo(1)%i(clouds(n)%halo(1)%npts) = i
	clouds(n)%halo(1)%j(clouds(n)%halo(1)%npts) = jm
	id_ha(i,jm) = clouds(n)%id
	lres = .true.
      endif
      if (lres) goto 100
    enddo
!
!  Subsequent layers
!
    do l = 2, 8
      do k = 1, clouds(n)%halo(l-1)%npts
        i = clouds(n)%halo(l-1)%i(k)
	j = clouds(n)%halo(l-1)%j(k)
        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        if (ip > NXO) ip = ip - NXO
        if (im < 1) im = im + NXO
        if (jp > NYO) jp = jp - NYO
        if (jm < 1) jm = jm + NYO
!
101     continue
        lres = .false.
        if ( id_ha(ip,j) /= clouds(n)%id ) then
          clouds(n)%halo(l)%npts = clouds(n)%halo(l)%npts + 1
	  clouds(n)%halo(l)%i(clouds(n)%halo(l)%npts) = ip
	  clouds(n)%halo(l)%j(clouds(n)%halo(l)%npts) = j
	  id_ha(ip,j) = clouds(n)%id
	  lres = .true.
	else if ( id_ha(i,jp) /= clouds(n)%id ) then
          clouds(n)%halo(l)%npts = clouds(n)%halo(l)%npts + 1
	  clouds(n)%halo(l)%i(clouds(n)%halo(l)%npts) = i
	  clouds(n)%halo(l)%j(clouds(n)%halo(l)%npts) = jp
	  id_ha(i,jp) = clouds(n)%id
	  lres = .true.
	else if ( id_ha(im,j) /= clouds(n)%id ) then
          clouds(n)%halo(l)%npts = clouds(n)%halo(l)%npts + 1
	  clouds(n)%halo(l)%i(clouds(n)%halo(l)%npts) = im
	  clouds(n)%halo(l)%j(clouds(n)%halo(l)%npts) = j
	  id_ha(im,j) = clouds(n)%id
	  lres = .true.
	else if ( id_ha(i,jm) /= clouds(n)%id ) then
          clouds(n)%halo(l)%npts = clouds(n)%halo(l)%npts + 1
	  clouds(n)%halo(l)%i(clouds(n)%halo(l)%npts) = i
	  clouds(n)%halo(l)%j(clouds(n)%halo(l)%npts) = jm
	  id_ha(i,jm) = clouds(n)%id
	  lres = .true.
        endif
        if (lres) goto 101
      enddo
    enddo
!
!  Cloud shape factor
!
    clouds(n)%sf = real(clouds(n)%halo(1)%npts) / (2.*sqrt(pi)*sqrt(real(clouds(n)%npts)))
!
    endif
  enddo
!
  deallocate( id_ha )
!
!----------------------------------------------------------!
!
!  print*, "Terminating find_halos"
!  
  return
  end subroutine find_halos
!
!      =================================================                            
  subroutine time_tracking ( t, timel, clouds, clouds_old, id_cl )  
!      =================================================                            
!
  IMPLICIT NONE
!  
  real, dimension (:,:) :: id_cl
  real, dimension (:) :: timel
  type (cloud_type), dimension(ncl) :: clouds, clouds_old
!
  logical :: newcl(1:ncl)
  integer :: t, n, m, l, k, j, max_id, min_n, mm
  integer :: overlap(1:ncl,1:ncl), maxn(1:ncl), maxm(1:ncl)
!
!----------------------------------------------------------!
!
  id_cl = 0.
  clouds%id = 0
  overlap = 0
  maxn = 0
  maxm = 0
  newcl = .false.
!
!  Find overlapping clouds and assign IDs
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
!
      do m = 1, ncl
	if ( clouds_old(m)%size >= nc_min ) then
!
	  do j = 1, clouds(n)%npts
            do l = 1, clouds_old(m)%npts
	      if ( clouds(n)%i(j) == clouds_old(m)%i(l) .and. clouds(n)%j(j) == clouds_old(m)%j(l) ) &
	      	overlap(n,m) = overlap(n,m) + 1
	    enddo
          enddo
!
          maxm(m) = maxloc(overlap(:,m),dim=1)
	endif
      enddo
!
      maxn(n) = maxloc(overlap(n,:),dim=1)     
    endif
  enddo
!  
!  IDs of largest overlapping clouds
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
      if ( .not.any(overlap(n,:) > 0) ) then
        newcl(n) = .true.
      else 
        do m = 1, ncl
	  if ( maxm(m) == n .and. maxn(n) == m ) then
	    clouds(n)%id = clouds_old(m)%id
            exit
	  endif
	enddo
	if ( m > ncl ) newcl(n) = .true.
      endif
    endif
  enddo
!
!  Set new IDs for new clouds
!
  max_id = maxval(clouds%id)
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
      if ( newcl(n) ) then
        max_id = max_id + 1
        clouds(n)%id = max_id
      endif
!
      do l = 1, clouds(n)%npts
        id_cl(clouds(n)%i(l),clouds(n)%j(l)) = clouds(n)%id
      enddo
    endif
  enddo
!
!  Clouds age
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min .and. .not.newcl(n) .and. t > 1 ) then
       clouds(n)%age = clouds(n)%age + (timel(t) - timel(t-1))
    else 
      clouds(n)%age = 0
    endif
  enddo   
!
  call qnbc_id ( id_cl )
!
!----------------------------------------------------------!
!  
  return
  end subroutine time_tracking
!
!      =================================================        
  subroutine cloud_diagnostics ( mom, qt, mse, buo, bef, ent, det, lwp, ctop, dq, qe, clouds, updrafts )   		          
!      =================================================                            
!
  IMPLICIT NONE
!  
  real, dimension (NXO,NYO) :: mom, qt, mse, buo, bef, ent, det, lwp, ctop
  real, dimension (NXO,NYO,5) :: dq
  type(cloud_type), dimension(ncl) :: clouds
  type(cloud_type), dimension(ncl,nup) :: updrafts
!
  real, dimension(8) :: mom_hal, qt_hal, buo_hal, nhal
  integer :: i, j, k, l, n, m, noff1, noff2
  real :: di, dj, didj, qe
!
!----------------------------------------------------------!
!              Find clouds & updrafts centres     	   !
!----------------------------------------------------------!
!
  do n = 1, ncl
!
    call find_center ( clouds(n)%npts, clouds(n)%i, clouds(n)%j, clouds(n)%ic, clouds(n)%jc )
!
    if ( with_up ) then
      do m = 1, nup
        call find_center ( updrafts(n,m)%npts, updrafts(n,m)%i, updrafts(n,m)%j, updrafts(n,m)%ic, updrafts(n,m)%jc )
      enddo
    endif
!
  enddo
!
!----------------------------------------------------------!
!	        Calculate cloud diagnostics		   !
!----------------------------------------------------------!
!
!  Calculate distance and angle between clouds
!
  clouds%mfl = 0.001
  clouds%qfl = 0.001
  clouds%dmdz = 0.
  clouds%qt  = 0.
  clouds%mse = 0.
  clouds%buo = 0.
  clouds%bef = 0.
  clouds%ent = 0.
  clouds%det = 0.
  clouds%lwp = 0.
  clouds%ctop = 0.
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
!
      noff1 = min(n-1,ndist/2)
      noff2 = min(ncl-n+1,ndist/2)
      k = 0
      do m = n-noff1, n+noff2
	k = k + 1
        if ( clouds(m)%size >= nc_min .and. m /= n ) then
!
	  if (from_centre) then		!Distance between clouds calculated from centres
            di = real(clouds(n)%ic - clouds(m)%ic)
	    if ( di > real(NXO/2) ) di = di - real(NXO)
	    if ( di < -real(NXO/2) ) di = di + real(NXO)
!
            dj = real(clouds(n)%jc - clouds(m)%jc)
	    if ( dj > real(NYO/2) ) dj = dj - real(NYO)
	    if ( dj < -real(NYO/2) ) dj = dj + real(NYO)
!
            clouds(n)%dist(k) = sqrt((dxl*di)**2. + (dyl*dj)**2.)
            if (clouds(n)%dist(k) > 0.) clouds(n)%angle(k) = (180./pi)*acos(dxl*di / clouds(n)%dist(k))
	  else				!Minimum distance between clouds
	    clouds(n)%dist(k) = 1000000.
	    do l = 1, clouds(n)%npts
	      do j = 1, clouds(m)%npts
                di = real(clouds(n)%i(l) - clouds(m)%i(j))
	        if ( di > real(NXO/2) ) di = di - real(NXO)
	        if ( di < -real(NXO/2) ) di = di + real(NXO)
!
                dj = real(clouds(n)%j(l) - clouds(m)%j(j))
	        if ( dj > real(NYO/2) ) dj = dj - real(NYO)
	        if ( dj < -real(NYO/2) ) dj = dj + real(NYO)
!
	        didj = sqrt( (dxl*di)**2. + (dyl*dj)**2. )
	        if ( didj < clouds(n)%dist(k) ) clouds(n)%dist(k) = didj
	      enddo
	    enddo
	  endif
!
	endif
      enddo
!
!  Total mass-flux and mean MSE per cloud
!
      do l = 1, clouds(n)%npts
	clouds(n)%dmdz= clouds(n)%dmdz + dq(clouds(n)%i(l),clouds(n)%j(l),1)
	clouds(n)%mfl = clouds(n)%mfl  + mom(clouds(n)%i(l),clouds(n)%j(l))
	clouds(n)%qfl = clouds(n)%qfl  + mom(clouds(n)%i(l),clouds(n)%j(l))*(qt(clouds(n)%i(l),clouds(n)%j(l)) - qe)
	
	clouds(n)%qt  = clouds(n)%qt   + (qt(clouds(n)%i(l),clouds(n)%j(l)) - qe) / real(clouds(n)%npts)
	clouds(n)%mse = clouds(n)%mse  + mse(clouds(n)%i(l),clouds(n)%j(l)) / real(clouds(n)%npts)
	clouds(n)%buo = clouds(n)%buo  + buo(clouds(n)%i(l),clouds(n)%j(l)) / real(clouds(n)%npts)
	clouds(n)%bef = clouds(n)%bef  + bef(clouds(n)%i(l),clouds(n)%j(l)) / real(clouds(n)%npts)
	clouds(n)%lwp = clouds(n)%lwp  + lwp(clouds(n)%i(l),clouds(n)%j(l)) / real(clouds(n)%npts)
	clouds(n)%ctop= clouds(n)%ctop + ctop(clouds(n)%i(l),clouds(n)%j(l)) / real(clouds(n)%npts)
	clouds(n)%mom = clouds(n)%mom  + mom(clouds(n)%i(l),clouds(n)%j(l)) / real(clouds(n)%npts)
	clouds(n)%ent = clouds(n)%ent  + ent(clouds(n)%i(l),clouds(n)%j(l)) / real(clouds(n)%npts)
	clouds(n)%det = clouds(n)%det  + det(clouds(n)%i(l),clouds(n)%j(l)) / real(clouds(n)%npts)
      enddo
!
    endif
  enddo
  mclouds = sum(clouds%mfl)
!
  do n = 1, ncl
    clouds(n)%mfl  = dxl*dyl*clouds(n)%mfl
    clouds(n)%qfl  = dxl*dyl*clouds(n)%qfl
    clouds(n)%dmdz = dxl*dyl*clouds(n)%dmdz
    clouds(n)%ent  = dxl*dyl*clouds(n)%ent
    clouds(n)%det  = dxl*dyl*clouds(n)%det
  enddo
!
!----------------------------------------------------------!
!	        Calculate halo diagnostics		   !
!----------------------------------------------------------!
!
  nhal = 0.
  mom_hal = 0.
  buo_hal = 0.
  qt_hal = 0.
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
!
      do l = 1, 8
        do k = 1, clouds(n)%halo(l)%npts
          nhal(l)    = nhal(l) + 1.
          mom_hal(l) = mom_hal(l) + mom(clouds(n)%halo(l)%i(k),clouds(n)%halo(l)%j(k))
          qt_hal(l)  = qt_hal(l) + qt(clouds(n)%halo(l)%i(k),clouds(n)%halo(l)%j(k))
          buo_hal(l) = buo_hal(l) + buo(clouds(n)%halo(l)%i(k),clouds(n)%halo(l)%j(k))
        enddo
      enddo
!
    endif
  enddo
!  
  do l = 1, 8
    print*, l, mom_hal(l)/nhal(l), qt_hal(l)/nhal(l), buo_hal(l)/nhal(l)
  enddo
!
!----------------------------------------------------------!
!               Calculate updraft diagnostics              !
!----------------------------------------------------------!
!
  if ( with_up ) then
!
  do n = 1, ncl
    if ( clouds(n)%size >= nc_min ) then
!
      do m = 1, nup
        if ( updrafts(n,m)%size >= nu_min ) then
!
          k = 0
          do l = 1, ncl
            if ( clouds(l)%size >= nc_min ) then
!
              do j = 1, nup
                if ( updrafts(l,j)%size >= nu_min .and. (l /= n .and. j /= m) ) then
                  di = real(updrafts(n,m)%ic - updrafts(l,j)%ic)
	          if ( di > real(NXO/2) ) di = di - real(NXO)
	          if ( di < -real(NXO/2) ) di = di + real(NXO)
!
                  dj = real(updrafts(n,m)%jc - updrafts(l,j)%jc)
	          if ( dj > real(NYO/2) ) dj = dj - real(NYO)
	          if ( dj < -real(NYO/2) ) dj = dj + real(NYO)
!
		  k = k+1
                  updrafts(n,m)%dist(k) = sqrt((dxl*di)**2. + (dyl*dj)**2.)
                  updrafts(n,m)%angle(k) = (180./pi)*acos(dxl*di / updrafts(n,m)%dist(k))
		endif
              enddo
!
            endif
          enddo
!
!  Mean mass-flux per updraft
!
          updrafts(n,m)%mfl = 0.
          do l = 1, updrafts(n,m)%npts
	    updrafts(n,m)%mfl = updrafts(n,m)%mfl + dxl*dyl*mom(updrafts(n,m)%i(l),updrafts(n,m)%j(l))
	  enddo
!
	endif
      enddo
!
    endif
  enddo
  uclouds = sum(updrafts%mfl)
!
  endif
!
!----------------------------------------------------------!
!
  return
  end subroutine cloud_diagnostics
!
!      =================================================                  		          
  subroutine cloud_stats ( clouds, clouds_old, updrafts, tav, nav, navb, nupav, nupavb,	 	&
  			   nvar, nvarb, nupvar, nupvarb, pdf, mav, mavb, mvar, mvarb, pdfm, 	&
  			   mupav, mupavb, mupvar, mupvarb, msummax, mupsummax, pdfupm, 		&
			   dav, dvar, pdfd, dupav, dupvar, pdfupd, pdft, pdfdq, pdfda, pdfdd )
!      =================================================                            
!
  IMPLICIT NONE
!
  integer :: k, l, m, n, ndist, nupdist, iblock, jblock, lblock
  real :: tav, nav, navb, nupav, nupavb, nvar, nvarb, nupvar, nupvarb
  real :: mav, mavb, mvar, mvarb, mupav, mupavb, mupvar, mupvarb, msummax, mupsummax
  real :: dav, dvar, dupav, dupvar
  real :: ff, fft, fft1, fft2, ffupt
!
  type (cloud_type), dimension(ncl) :: clouds, clouds_old
  type (cloud_type), dimension(ncl,nup) :: updrafts
!
  real, dimension(ncl,2) :: quart
  real, dimension(ncl) :: siz
  real, dimension(NXO/nblock*NYO/nblock) :: navbl, mavbl, nupavbl, mupavbl
!
  real, dimension(npdf) :: pdf, pdfm, pdfupm, pdfd, pdfupd, pdft
  real, dimension(npdf/2,npdf/2) :: pdfda, pdfdd
  real, dimension(npdf,2) :: pdfdq
!
!----------------------------------------------------------!
!                Mean mass flux and number                 !
!----------------------------------------------------------!
!
  if (nclouds > 0) then
!
!  Basic cloud stats
!
    nav  = sum(real(clouds%size))/real(nclouds)
    nvar = sum((real(clouds%size) - nav)**2.) / real(nclouds)
!
    mav  = sum(clouds%mfl)/real(nclouds)
    mvar = sum((clouds%mfl - mav)**2.) / real(nclouds)
!    
    msummax = sum(clouds%mfl) / maxval(max(clouds%mfl,1.))
!
!  Basic updraft stats
!
    nupav  = sum(updrafts%size) / real(sum(nupdrafts))
    nupvar = sum((updrafts%size - nupav)**2.) / real(sum(nupdrafts))
!
    mupav  = sum(updrafts%mfl) / real(sum(nupdrafts))
    mupvar = sum((updrafts%mfl - mupav)**2.) / real(sum(nupdrafts))
!    
    mupsummax = sum(updrafts%mfl) / maxval(max(updrafts%mfl,1.))
!
!  Cloud lifetime stats
!
    tav = 0.
    k = 0
    if (time_clt) then
      do n = 1, ncl
        if ( clouds(n)%age < clouds_old(n)%age .and. clouds_old(n)%age >= tmin ) then
          tav = tav + clouds_old(n)%age
          k = k + 1
        endif
      enddo
      tav = tav/max(real(k),1.)
    endif
!
!  Stats per blocks
!
    navbl = 0.
    mavbl = 0.
    if (nblock > 1) then
      do n = 1, ncl
        if ( clouds(n)%size >= nc_min ) then
          iblock = (clouds(n)%ic-is)/nblock + 1
	  jblock = (clouds(n)%jc-js)/nblock + 1
	  lblock = iblock + (jblock - 1)*(NXO/nblock)
	  navbl(lblock) = navbl(lblock) + 1.
	  mavbl(lblock) = mavbl(lblock) + clouds(n)%mfl
        endif
      enddo
    endif
!
    navb = sum(navbl)/real(size(navbl))
    nvarb = sum((navbl - navb)**2.) / real(size(navbl))
!
    mavb = sum(mavbl)/real(size(mavbl))
    mvarb = sum((mavbl - mavb)**2.) / real(size(mavbl))
!
    nupavbl = 0.
    mupavbl = 0.
    if (nblock > 1) then
      do n = 1, ncl
        do m = 1, nup
          if ( updrafts(n,m)%size >= nu_min ) then
            iblock = (updrafts(n,m)%ic-is)/nblock + 1
	    jblock = (updrafts(n,m)%jc-js)/nblock + 1
	    lblock = iblock + (jblock - 1)*(NXO/nblock)
	    nupavbl(lblock) = nupavbl(lblock) + 1.
	    mupavbl(lblock) = mupavbl(lblock) + updrafts(n,m)%mfl
          endif
        enddo
      enddo
    endif
!
    nupavb = sum(nupavbl)/real(size(nupavbl))
    nupvarb = sum((nupavbl - nupavb)**2.) / real(size(nupavbl))
!
    mupavb = sum(mupavbl)/real(size(mupavbl))
    mupvarb = sum((mupavbl - mupavb)**2.) / real(size(mupavbl))
!
!----------------------------------------------------------!
!                Size/Mass flux distributions              !
!----------------------------------------------------------!
!
    ndist=0
    nupdist=0
    fft = 0.
    fft1 = 0.
    fft2 = 0.
    ffupt = 0.
    siz=0.0
    pdf=0.0
    pdfm=0.0
    pdfupm=0.0
    pdfd=0.0
    pdfupd=0.0
    pdft=0.0
    pdfdq=0.0
    pdfda=0.0
    pdfdd=0.0
    dav=0.0
    dvar=0.0
    dupav=0.0
    dupvar=0.0
!
    k = 0
    do n = 1, ncl
      if (clouds(n)%size >= nc_min) then
!    
!  Cloud PDF
!  
        do m = 1, size(pdf)
          if ( clouds(n)%mfl >= (real(m)-0.5)*dm .and. clouds(n)%mfl < (real(m)+0.5)*dm) then
            pdfm(m)=pdfm(m)+1.
          endif
!
          if ( clouds(n)%size >= (real(m)-0.5) .and. clouds(n)%size < (real(m)+0.5) ) then
            pdf(m)=pdf(m)+1.
          endif
        enddo
!    
!  Cloud lifetime PDF
!  
        do m = 1, size(pdft)
          if ( clouds(n)%age > real(m-1)*30.*tmin/real(size(pdft)) .and. clouds(n)%age <= real(m)*30.*tmin/real(size(pdft))) then
            pdft(m)=pdft(m)+1.
          endif
        enddo
!
!  Updraft PDF
!
        do m = 1, maxval(updrafts(n,:)%id)
	  if ( updrafts(n,m)%size >= nu_min ) then
!
            do l = 1, size(pdfupm)
              if ( updrafts(n,m)%mfl > real(l-1)*dm .and. updrafts(n,m)%mfl <= real(l)*dm ) then
                pdfupm(l)=pdfupm(l)+1.
                exit
              endif
            enddo
!
	  endif
        enddo
!
	k = k + 1
	siz(k) = real(clouds(n)%size)*dxl*dyl
      endif
    enddo
!
!  Size quartiles
!
    call sort( siz, nclouds )
!
    quart = 0.
    do n = 1, ncl
      if ( clouds(n)%size >= nc_min ) then
!    
        do l = 1, nclouds
	  if ( siz(l) >= real(nc_min)*dxl*dyl .and. clouds(n)%size*dxl*dyl == siz(l) ) then
	    if ( l <= k/4 ) then
              quart(n,1) = 1.
	    else if ( l >= 3*k/4 ) then
              quart(n,2) = 1.
            endif
          endif
	enddo
!
      endif
    enddo
!
!----------------------------------------------------------!
!                 Cloud distance statistics                !
!----------------------------------------------------------!
!
    do n = 1, ncl
      if ( clouds(n)%size >= nc_min ) then
!
        do m = 1, ncl
          if ( clouds(m)%size >= nc_min .and. m /= n ) then
!
!  Cloud distance PDF (total and by quartiles)
!
            do l = 1, size(pdfd)
	      ff = ((dxl*real(NXO)/2.)**2.) / (((real(l)+0.5)*dmin)**2. - ((real(l)-0.5)*dmin)**2.)
              if ( clouds(n)%dist(m) >= (real(l)-0.5)*dmin .and. clouds(n)%dist(m) < (real(l)+0.5)*dmin .and. 	&
	      	   clouds(n)%dist(m) < min(dxl*real(NXO)/2.,dyl*real(NYO)/2.) ) then
                pdfd(l) = pdfd(l) + ff
		fft = fft + 1.
!
		if ( quart(n,1) == 1. ) then
                  pdfdq(l,1) = pdfdq(l,1) + ff
		  fft1 = fft1 + 1.
		else if ( quart(n,2) == 1. ) then
                  pdfdq(l,2) = pdfdq(l,2) + ff
		  fft2 = fft2 + 1.
		endif
                exit
              endif
            enddo
            dav = dav + clouds(n)%dist(m)
            dvar = dvar + clouds(n)%dist(m)*clouds(n)%dist(m)
            ndist = ndist+1
!
!  Joint distance-angle PDF
!
            do l = 1, size(pdfda,1)
	      ff = ((dxl*real(NXO)/2.)**2.) / (((real(l)+0.5)*dmin)**2. - ((real(l)-0.5)*dmin)**2.)
	      do k = 1, size(pdfda,2)
                if ( clouds(n)%dist(m) >= (real(l)-0.5)*dmin*real(size(pdfd))/real(size(pdfda,1)) .and.	&
	      	     clouds(n)%dist(m) < (real(l)+0.5)*dmin*real(size(pdfd))/real(size(pdfda,1)) .and.	&
		     clouds(n)%angle(m) >= (real(k)-0.5)*amin/real(size(pdfda,2)) .and.			&
	      	     clouds(n)%angle(m) < (real(k)+0.5)*amin/real(size(pdfda,2)) .and.			&
		     clouds(n)%dist(m) < min(dxl*real(NXO)/2.,dyl*real(NYO)/2.) ) then
                  pdfda(l,k) = pdfda(l,k) + ff
                  exit
                endif
              enddo
	    enddo
!
!  Joint distance-size PDF
!
            do l = 1, size(pdfdd,1)
	      ff = ((dxl*real(NXO)/2.)**2.) / (((real(l)+0.5)*dmin)**2. - ((real(l)-0.5)*dmin)**2.)
	      do k = 1, size(pdfdd,2)
                if ( clouds(n)%dist(m) >= (real(l)-0.5)*dmin*real(size(pdfd))/real(size(pdfdd,1)) .and.	&
	      	     clouds(n)%dist(m) < (real(l)+0.5)*dmin*real(size(pdfd))/real(size(pdfdd,1)) .and.	&
		     clouds(n)%size*dxl*dyl >= (real(k)-0.5)*smin/real(size(pdfdd,2)) .and.		&
	      	     clouds(n)%size*dxl*dyl < (real(k)+0.5)*smin/real(size(pdfdd,2)) .and.		&
		     clouds(n)%dist(m) < min(dxl*real(NXO)/2.,dyl*real(NYO)/2.) ) then
                  pdfdd(l,k) = pdfdd(l,k) + ff
                  exit
                endif
              enddo
	    enddo
!
          endif
        enddo
!
!  Updraft distance PDF
!
        do m = 1, nup
	  if ( updrafts(n,m)%size >= nu_min ) then
!
            do l = 1, sum(nupdrafts)
              do k = 1, size(pdfupd)
	        ff = ((dxl*real(NXO)/2.)**2.) / ((real(k)*dmin)**2. - (real(k-1)*dmin)**2.)
                if ( updrafts(n,m)%dist(l) > real(k-1)*dmin .and. updrafts(n,m)%dist(l) <= real(k)*dmin .and.	&
	    	     updrafts(n,m)%dist(l) < min(dxl*real(NXO)/2.,dyl*real(NYO)/2.) .and. 			&
		     updrafts(n,m)%dist(l) /= 0. ) then
            	  pdfupd(k) = pdfupd(k) + ff
		  ffupt = ffupt + 1.
            	  exit
                endif
              enddo
              dupav = dupav + updrafts(n,m)%dist(l)
              dupvar = dupvar + updrafts(n,m)%dist(l)*updrafts(n,m)%dist(l)
              nupdist=nupdist+1
	    enddo
!
          endif
        enddo
!	
      endif
    enddo
!
!  Normalization
!
    dav = dav/real(max(ndist,1))
!    dvar = sum((clouds%dist - dav)**2.)/real(max(ndist,1))
    dupav = dupav/real(max(nupdist,1))
!    dupvar = sum((updrafts%dist - dupav)**2.)/real(max(nupdist,1))
!
    pdfd = pdfd/max(fft,1.)
    pdfupd = pdfupd/max(ffupt,1.)
    pdfdq(:,1) = pdfdq(:,1)/max(fft1,1.)
    pdfdq(:,2) = pdfdq(:,2)/max(fft2,1.)
    pdfda = pdfda/max(fft,1.)
    pdfdd = pdfdd/max(fft,1.)
  endif
!
!----------------------------------------------------------!
!
  return
  end subroutine cloud_stats
!
!      =================================================                  		          
  subroutine cloud_output ( unit1, unit2, unit3, t, tav, nmax, nav, nvar, navb, nvarb, 	&
			    nupav, nupvar, nupavb, nupvarb, mav, mvar, mavb, mvarb,	&
			    mupav, mupvar, mupavb, mupvarb, msummax, mupsummax,		&
			    dav, dvar, dupav, dupvar, pdf, pdfm, pdfupm, pdfd, pdfupd,  &
			    pdft, pdfdq, pdfda, pdfdd )
!      =================================================                            
!
  IMPLICIT NONE
!
  integer :: k, l
  integer :: unit1, unit2, unit3, nmax, nclouds_tot
  real :: t, tav, nav, nvar, navb, nvarb, nupav, nupvar, nupavb, nupvarb
  real :: mav, mvar, mavb, mvarb, mupav, mupvar, mupavb, mupvarb, msummax, mupsummax
  real :: dav, dvar, dupav, dupvar
!
  real, dimension(npdf) :: pdf, pdfm, pdfupm, pdfd, pdfupd, pdft
  real, dimension(npdf/2,npdf/2) :: pdfda, pdfdd
  real, dimension(npdf,2) :: pdfdq
!
  real, dimension(27) :: ts_tab
  real, dimension(npdf) :: pdf_mean, pdfm_mean, pdfupm_mean, pdfd_mean, pdfupd_mean, pdft_mean
  real, dimension(npdf/2,npdf/2) :: pdfda_mean, pdfdd_mean
  real, dimension(npdf,2) :: pdfdq_mean
!
  integer :: count = 0
  real, dimension(npdf/2) :: dd, ss, aa
!
  save
!
!----------------------------------------------------------!
!             Output time-series basic stats               !
!----------------------------------------------------------!
!
  if (count == 0) then
    nclouds_tot = 0
    pdf_mean = 0.0
    pdfm_mean = 0.0
    pdfd_mean = 0.0
    pdfupm_mean = 0.0
    pdfupd_mean = 0.0
    pdft_mean = 0.0
    pdfdq_mean = 0.0
    pdfda_mean = 0.0
    pdfdd_mean = 0.0
  endif
!
  if (nclouds > 0) then
    nclouds_tot = nclouds_tot + nclouds
!
    ts_tab(1) = real(nclouds_tot)	! Cumulated number of clouds
    ts_tab(2) = real(nclouds)		! Instantaneous number of clouds
    ts_tab(3) = real(sum(nupdrafts))	! Instantaneous number of updrafts
    ts_tab(4) = tav			! Average cloud age
    ts_tab(5) = nmax			! Max cloud size
    ts_tab(6) = navb			! Average cloud number per block
    ts_tab(7) = nvarb			! Variance of cloud number per block
    ts_tab(8) = nav			! Average cloud size
    ts_tab(9) = nvar			! Variance of cloud size
    ts_tab(10) = mav			! Average mass flux per cloud
    ts_tab(11) = mvar			! Variance of mass flux per cloud
    ts_tab(12) = mavb			! Average mass flux per block
    ts_tab(13) = mvarb			! Variance of mass flux per block
    ts_tab(14) = nupavb			! Number of updrafts per block
    ts_tab(15) = nupvarb		! Variance of updrafts per block
    ts_tab(16) = nupav			! Average updraft size
    ts_tab(17) = nupvar			! Variance of updraft size
    ts_tab(18) = mupav			! Average mass flux per updraft
    ts_tab(19) = mupvar			! Variance of mass flux per updraft
    ts_tab(20) = mupavb			! Average updraft mass flux per block
    ts_tab(21) = mupvarb		! Variance of updraft mass flux per block
    ts_tab(22) = dav			! Average distance between clouds
    ts_tab(23) = dvar			! Variance of distance between clouds
    ts_tab(24) = dupav			! Average distance between updrafts
    ts_tab(25) = dupvar			! Variance of distance between updrafts
    ts_tab(26) = msummax		! Cloud mass flux sum over max
    ts_tab(27) = mupsummax		! Updraft mass flux sum over max
!
    call output_clts (unit1, t, ts_tab)
!
!----------------------------------------------------------!
!             Average and output distributions     	   !
!----------------------------------------------------------!
!
    count=count+1
    pdf_mean=pdf_mean+pdf
    pdfm_mean=pdfm_mean+pdfm
    pdfupm_mean=pdfupm_mean+pdfupm
    pdfd_mean=pdfd_mean+pdfd
    pdfupd_mean=pdfupd_mean+pdfupd
    pdft_mean=pdft_mean+pdft
    pdfdq_mean=pdfdq_mean+pdfdq
    pdfda_mean=pdfda_mean+pdfda
    pdfdd_mean=pdfdd_mean+pdfdd
!
!  1D PDFs
!
    call output_pdf ( unit2, pdf_mean/real(count), pdfm_mean/real(count), 	&
    		      pdfupm_mean/real(count), pdfd_mean/real(count), 		&
		      pdfupd_mean/real(count), pdft_mean/real(count), 		&
		      pdfdq_mean(:,1)/real(count), pdfdq_mean(:,2)/real(count) )
!
!  2D PDFs
!
    do k = 1, npdf/2
      dd(k) = (real(k)-0.5)*dmin*2.
    enddo
    do k = 1, npdf/2
      aa(k) = (real(k)-0.5)*amin/real(size(pdfda,2))
    enddo
    do k = 1, npdf/2
      ss(k) = (real(k)-0.5)*smin/real(size(pdfdd,2))
    enddo
!
    do k = 1, npdf/2
      do l = 1, npdf/2
!        write(unit3,100)dd(k),aa(l),pdfda_mean(k,l)/max(sum(pdfda_mean(:,l)),1.)
        write(unit3,100)dd(k),ss(l),pdfdd_mean(k,l)/max(sum(pdfdd_mean(:,l)),1.)
      enddo
      write(unit3,*)
    enddo
!
100  format (2x, 3e15.7)
!
  endif 
!
!----------------------------------------------------------!
!
!    print*, "SUCCESS writing cloud tracking outputs"
!
  return
  end subroutine cloud_output
!
!  ====================================================
  subroutine find_center ( npts, i, j, ic, jc )
!  ====================================================
!
  IMPLICIT NONE
!
  integer :: npts, ic, jc
  integer, dimension(200) :: i, j
!
  integer :: l, sizei, sizej
!
!----------------------------------------------------------!
!
    if ( npts >= nc_min ) then
!
      sizei = 0
      do l = 1, npts
        if ( i(l)-NXO/2 > is ) then
	  ic = ic + (i(l)-NXO)
        else 
          ic = ic + i(l)
	endif
	sizei = sizei+1
      enddo
!
      sizej = 0
      do l = 1, npts
        if ( j(l)-NYO/2 > js ) then
	  jc = jc + (j(l)-NYO)
        else 
          jc = jc + j(l)
	endif
	sizej = sizej+1
      enddo
!
      ic = ic / sizei
      jc = jc / sizej
!
      if (ic < is) ic = ic + NXO
      if (jc < js) jc = jc + NYO
!
    endif
!
!----------------------------------------------------------!
!
  return
  end subroutine find_center
!
!  ====================================================
  subroutine qnbc_clouds ( ntot, npts, if1, if2, clouds )
!  ====================================================
!
  IMPLICIT NONE
!
  logical :: ex
  integer :: ntot, npts
  integer :: m, i, j, ii
  integer, dimension (:,:) :: if1, if2
  type (cloud_type), dimension(1:ntot) :: clouds
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
!  Periodicity in x
!
  do j = min(js+1,je), je
    if ( if1(is,j) == 1 ) then
!  
    ex=.false.
    if ( if1(ie,j) == 1 .and. if2(ie,j) == 1 ) then
      do m = 1, npts
  	if ( any(ie == clouds(m)%i) .and. any(j == clouds(m)%j) ) then
  	  clouds(m)%npts=clouds(m)%npts+1
  	  ii=m
  	  ex=.true.
  	  exit 
  	endif
      enddo
    else if (NYO > 1) then
      if ( if1(is,j-1) == 1 .and. if2(is,j-1) == 1 ) then
        do m = 1, npts
  	  if ( any(is == clouds(m)%i) .and. any(j-1 == clouds(m)%j) ) then
  	    clouds(m)%npts=clouds(m)%npts+1
  	    ii=m
  	    ex=.true.
  	    exit 
  	  endif
        enddo
      endif
    endif 
!
    if ( .not.ex ) then
      npts=npts+1
      clouds(npts)%npts=1
      clouds(npts)%id=npts
      ii=npts
    endif
!
    clouds(ii)%i(clouds(ii)%npts)=is
    clouds(ii)%j(clouds(ii)%npts)=j
    clouds(ii)%size=clouds(ii)%npts
!
    endif 
  enddo
!
!  Periodicity in y
!    
  if (NYO > 1) then
    do i = is, ie-1
      if ( if1(i,js) == 1 ) then
!  
      ex=.false.
      if ( if1(i,je) == 1 .and. if2(i,je) == 1 ) then
        do m = 1, npts
  	  if ( any(i == clouds(m)%i) .and. any(je == clouds(m)%j) ) then
  	    clouds(m)%npts=clouds(m)%npts+1
  	    ii=m
  	    ex=.true.
  	    exit
  	  endif
        enddo
      else if ( if1(i+1,js) == 1 .and. if2(i+1,js) == 1 ) then
        do m = 1, npts
  	  if ( any(js == clouds(m)%j) .and. any(i-1 == clouds(m)%i) ) then
  	    clouds(m)%npts=clouds(m)%npts+1
  	    ii=m
  	    ex=.true.
  	    exit 
  	  endif
        enddo
      endif 
!
      if ( .not.ex ) then
        npts=npts+1
        clouds(npts)%npts=1
        clouds(npts)%id=npts
        ii=npts
      endif
!
      clouds(ii)%i(clouds(ii)%npts)=i
      clouds(ii)%j(clouds(ii)%npts)=js
      clouds(ii)%size=clouds(ii)%npts
!
      endif 
    enddo
  endif
!
!-----------------------------------------------------------!
!
return
end
!
!  ====================================================
  subroutine qnbc_id ( x )
!  ====================================================
!
  IMPLICIT NONE
!
  integer :: i, j
  real, dimension(NXO,NYO) :: x
!
!----------------------------------------------------------!
!                   Periodic lateral BCs                   !
!----------------------------------------------------------!
!
!  Periodicity in x
!
  do j = js, je
    if ( x(ie,j) /= 0. .and. x(is,j) /= 0. ) then
      where ( x == x(is,j) ) x = x(ie,j)
    endif
  enddo
!
!  Periodicity in y
!
  if (NYO > 1) then
    do i = is, ie
      if ( x(i,je) /= 0. .and. x(i,js) /= 0. ) then
        where ( x == x(i,js) ) x = x(i,je)
      endif
    enddo
  endif
!
!-----------------------------------------------------------!
!
return
end
!
!    =================================================
  subroutine output_clts ( unit_number, t, ts_tab )	
!    =================================================
!
      integer :: unit_number
      real, dimension(27) :: ts_tab
      real :: t

! ============================================================
!
      write(unit_number,101) t,ts_tab
!
101  format (2x, 28e15.7)
!	
RETURN
END

!    =================================================
  subroutine output_pdf ( unit_number, pdf1, pdf2, pdf3, pdf4, pdf5, pdf6, pdf7, pdf8 )	
!    =================================================
!
	integer :: unit_number, i
	real, dimension(npdf) :: pdf1,pdf2,pdf3,pdf4,pdf5,pdf6,pdf7,pdf8

! ============================================================
!
    write(unit_number,*)
    do i = 1, npdf
      write(unit_number,100)real(i),pdf1(i),pdf2(i),pdf3(i),pdf4(i),pdf5(i),pdf6(i),pdf7(i),pdf8(i)
    enddo
!
100  format (2x, 9e15.7)
!	
RETURN
END

end module searchmod
