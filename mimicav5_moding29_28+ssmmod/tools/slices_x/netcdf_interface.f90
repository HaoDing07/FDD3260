
module netcdfmod

  USE netcdf
  USE sharedmod
  
  implicit none
        
  private
  
  character(len=10), dimension(nvarin) :: varlist=(/ 'U         ', 'V         ', 'W         ', 'PT        ', 'PTV       ', 'Qv        ', 'Qc        ', 'Qr        ', 'RH        ', 'Buoy      ', 'Beff      ', 'MSE       ', 'WVP       ',  'CT_height ', 'Pt_var    ', 'Qv_var    ', 'W_tend_1  ', 'W_tend_2  ', 'W_tend_3  ', 'W_tend_4  ', 'W_tend_5  ', 'W_tend_6  ', 'W_tend_7  ', 'W_tend_8  ', 'W_tend_9  ', 'QC_tend_1 ', 'QC_tend_2 ', 'QC_tend_3 ', 'QC_tend_4 ', 'QC_tend_5 ', 'QC_tend_6 ', 'QC_tend_7 ', 'QC_tend_8 ', 'QC_tend_9 ', 'QR_tend_1 ', 'QR_tend_2 ', 'QR_tend_3 ', 'QR_tend_4 ', 'QR_tend_5 ', 'QR_tend_6 ', 'QR_tend_7 ', 'QR_tend_8 ', 'QR_tend_9 ', 'PT_tend_1 ', 'PT_tend_2 ', 'PT_tend_3 ', 'PT_tend_4 ', 'PT_tend_5 ', 'PT_tend_6 ', 'PT_tend_7 ', 'PT_tend_8 ', 'PT_tend_9 ', 'PTV_tend_1', 'PTV_tend_2', 'PTV_tend_3', 'PTV_tend_4', 'PTV_tend_5', 'PTV_tend_6', 'PTV_tend_7', 'PTV_tend_8', 'PTV_tend_9', 'P_tot     ', 'Qt        ', 'Divergence', 'Ent_wcore ', 'Det_wcore ', 'VOF       ' /)
  character(len=*), parameter :: UNITS = "units", INFO = "description"
  character(len=*), parameter :: REC_NAME = "time", REC_UNITS = "s"
  character(len=*), parameter :: X_NAME = "X", X_UNITS = "m", Y_NAME = "Y", Y_UNITS = "m"

  real, dimension(:,:,:,:), allocatable :: var
  
  public :: read_dim, write_slice
  
  contains

  subroutine read_dim ( file_output )

  implicit none

  logical :: ex
  integer :: i

  integer :: ncid, include_parents=0
  integer, dimension(:), allocatable :: dimids, varids, count
  character(len = *) :: file_output
  character(len = 100) :: FILE_NAME, name
  character(len = 100), dimension(:), allocatable :: varname
  
  ! Create the file(s).
  FILE_NAME = trim(file_output)
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  if (.not.ex) then
    print*, 'File does not exist, exiting'
    stop
  else
    print*, ''
    print*, 'Reading file: ', trim(FILE_NAME)
  endif

  !
  ! Open new file and get dimensions
  !  
  call check( nf90_open(FILE_NAME, nf90_write, ncid) )
      
  call check( nf90_inquire(ncid, ndims, nvar) )
  
  allocate ( dimids(1:ndims), varids(1:nvar) )

  call check( nf90_inq_dimids(ncid, ndims, dimids, include_parents) )
    
  call check( nf90_inq_varids(ncid, nvar, varids) )
  
  !
  !  Inquire and get dimensions
  !
  if ( allocated(xx) ) deallocate ( xx )
  if ( allocated(yy) ) deallocate ( yy )
  if ( allocated(zz) ) deallocate ( zz )
  if ( allocated(tt) ) deallocate ( tt )
  if ( allocated(state) ) deallocate ( state )
  if ( allocated(wind) ) deallocate ( wind )
  if ( allocated(tend) ) deallocate ( tend )
  
  if (ndims == 4) then
  
    call check( nf90_inquire_dimension(ncid, dimids(1), name=name, len=NXO) )
  
    allocate( xx(1:NXO) )

    call check( nf90_inquire_dimension(ncid, dimids(2), name=name, len=NYO) )
  
    allocate( yy(1:NYO) )
    
    call check( nf90_get_var(ncid, varids(2), yy, start=(/ 1 /), count=(/ NYO /)) )
      
    call check( nf90_inquire_dimension(ncid, dimids(3), name=name, len=NZO) )
  
    allocate( zz(1:NZO) )
    
    call check( nf90_get_var(ncid, varids(3), zz, start=(/ 1 /), count=(/ NZO /)) )
      
  else if (ndims == 3) then
  
    NXO = 1
  
    allocate( xx(1:NXO) )
  
    call check( nf90_inquire_dimension(ncid, dimids(1), name=name, len=NYO) )
  
    allocate( yy(1:NYO) )
    
    call check( nf90_get_var(ncid, varids(1), yy, start=(/ 1 /), count=(/ NYO /)) )
      
    call check( nf90_inquire_dimension(ncid, dimids(2), name=name, len=NZO) )
  
    allocate( zz(1:NZO) )
    
    call check( nf90_get_var(ncid, varids(2), zz, start=(/ 1 /), count=(/ NZO /)) )
      
  else if (ndims == 2) then
  
    NYO = 1
    NZO = 1
  
    allocate( yy(1:NYO) )
    allocate( zz(1:NZO) )
    
  else
    print*, 'Wrong number of dimensions, exiting'
    stop    
  endif    
      
  call check( nf90_inquire_dimension(ncid, dimids(ndims), name=name, len=NTO) )
  
  allocate( tt(1:NTO) )
    
  call check( nf90_get_var(ncid, varids(ndims), tt, start=(/ 1 /), count=(/ NTO /)) )
      
  print*, 'Output times found in file: ', NTO, tt
      
  !
  !  Inquire and get variables
  !
  allocate ( count(1:ndims-1), varname(1:nvar), var(1:NXO,1:NYO,1:NZO,1:NTO) )
  
  allocate ( state(1:NXO,1:NYO,1:NZO,1:NTO), wind(1:NXO,1:NYO,1:NZO,1:NTO), tend(1:NXO,1:NYO,1:NZO,1:NTO) )
  
  ! These settings tell netcdf to write first timestep of data.
  if (ndims == 4) then
    count = (/ NXO, NYO, NZO, NTO /)
  else if (ndims == 3) then
    count = (/ NYO, NZO, NTO /)
  else if (ndims == 2) then
    count = (/ NXO, NTO /)
  endif

  do i = ndims+1, nvar
  
    call check( nf90_inquire_variable(ncid, varids(i), name=varname(i)) )
    
    if ( any(varname(i) == varlist) ) then
    
      call check( nf90_get_var(ncid, varids(i), var, count=count) )
      
      call store( varname(i), var )
      
      print*, 'Found variable in file: ', trim(varname(i))
  
    endif
  
  enddo
  
  deallocate( dimids, varids, count, varname, var )
  
  print*, "SUCCESS reading netcdf file ", trim(FILE_NAME)
  print*, 
  
  end subroutine read_dim

  subroutine write_slice (t,height,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
  
  implicit none
  
  integer, intent(in) :: t
  real, intent(in) :: height
  real, dimension(1:NXO,1:NYO), intent(in) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
  
  logical :: ex, new
  integer :: ncid, ntmp, ndim, include_parents=0
  integer :: y_dimid, x_dimid, t_dimid, rec, i, n
  integer, dimension(:), allocatable :: start, count
  integer, dimension(:), allocatable :: dimids, varids

  character (len = 100) :: FILE_NAME
  character(len=10) :: name
  character(len=4) :: car4

  ! Create the file.
  write(car4,'(i4)')floor(height)
  FILE_NAME = './slice_'//trim(adjustl(car4))//'.nc'  
  
  inquire(FILE=trim(FILE_NAME), EXIST=ex)
  new = .not.ex
  
  ntmp = 8
  ndim = 3
  
  allocate( start(1:ndim), count(1:ndim-1) )
  allocate( dimids(1:ndim), varids(1:ntmp) )
  
  !
  ! New file: define all variables
  !
  if (new) then
    call check( nf90_create(FILE_NAME, nf90_clobber, ncid) )
  
    ! Define the dimensions. The record dimension is defined to have
    ! unlimited length - it can grow as needed. In this example it is
    ! the time dimension.
    call check( nf90_def_dim(ncid, X_NAME, NXO, x_dimid) )
    call check( nf90_def_dim(ncid, Y_NAME, NYO, y_dimid) )
    call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, t_dimid) )

    ! Define the coordinate variables. We will only define coordinate
    ! variables for iy and ix.  Ordinarily we would need to provide
    ! an array of dimension IDs for each variable's dimensions, but
    ! since coordinate variables only have one dimension, we can
    ! simply provide the address of that dimension ID (x_dimid) and
    ! similarly for (y_dimid).
    call check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, varids(1)) )
    call check( nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, varids(2)) )
    call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, t_dimid, varids(3)) )

    ! Assign units attributes to coordinate variables.
    call check( nf90_put_att(ncid, varids(1), UNITS, X_UNITS) )
    call check( nf90_put_att(ncid, varids(2), UNITS, Y_UNITS) )
    call check( nf90_put_att(ncid, varids(3), UNITS, REC_UNITS) )

    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables. Both of the netCDF variables we are creating
    ! share the same four dimensions. In Fortran, the unlimited
    ! dimension must come last on the list of dimids.
    dimids = (/ x_dimid, y_dimid, t_dimid /)
  
    n = ndim
    call register_var ('wmo  ', n, ncid, dimids, varids)
    call register_var ('bre  ', n, ncid, dimids, varids)
    call register_var ('rh   ', n, ncid, dimids, varids)
    call register_var ('ql   ', n, ncid, dimids, varids)
    call register_var ('cid  ', n, ncid, dimids, varids)
    call register_var ('uid  ', n, ncid, dimids, varids)

  ! End define mode.
    call check( nf90_enddef(ncid) )
  
  ! Write the coordinate variable data.   
    call check( nf90_put_var(ncid, varids(1), xx) )
    call check( nf90_put_var(ncid, varids(2), yy) )
    call check( nf90_put_var(ncid, varids(3), tt(t)) )
    
    ! These settings tell netcdf to write first timestep of data.
    count = (/ NXO, NYO /)
    start = (/ 1, 1, 1 /)
  else
  
    !
    ! Existing file: appending new dataset
    !  
    call check( nf90_open(FILE_NAME, nf90_write, ncid) )
  	
    call check( nf90_inq_dimids(ncid, ndim, dimids, include_parents) )
  	
    call check( nf90_inq_varids(ncid, ntmp, varids) )
  	
    call check( nf90_inquire_dimension(ncid, dimids(3), name=name, len=NTO) )

    call check( nf90_put_var(ncid, varids(ndim), tt(t), start=(/ NTO+1 /)) )
    
    ! These settings tell netcdf to write first timestep of data.
    count = (/ NXO, NYO /)
    start = (/ 1, 1, NTO+1 /)
  endif
        
  ! Write the pretend data. 
  n = ndim
  call check( nf90_put_var(ncid, varids(n+1), tmp1, start=start, count=count) )
  call check( nf90_put_var(ncid, varids(n+2), tmp2, start=start, count=count) )
  call check( nf90_put_var(ncid, varids(n+3), tmp3, start=start, count=count) )
  call check( nf90_put_var(ncid, varids(n+4), tmp4, start=start, count=count) )
  call check( nf90_put_var(ncid, varids(n+5), tmp5, start=start, count=count) )
  call check( nf90_put_var(ncid, varids(n+6), tmp6, start=start, count=count) )
  
  call check( nf90_close(ncid) )
  
  deallocate (start, count, dimids, varids)
  
  print*,  "SUCCESS writing netcdf file ", trim(FILE_NAME)
  
  end subroutine write_slice
  
  subroutine store( varname, var )
  
  implicit none
  
  character(len=*) :: varname
  real, dimension(:,:,:,:) :: var
  
  if (trim(varname) == 'U') then
    wind%u = var
  else if (trim(varname) == 'V') then
    wind%v = var    
  else if (trim(varname) == 'W') then
    wind%w = var   
  else if (trim(varname) == 'Divergence') then
    wind%div = var    
  else if (trim(varname) == 'P_tot') then
    state%p = var   
  else if (trim(varname) == 'Qc') then
    state%qc = var  
  else if (trim(varname) == 'Qr') then
    state%qr = var  
  else if (trim(varname) == 'Qv') then
    state%qv = var  
  else if (trim(varname) == 'Qt') then
    state%qt = var  
  else if (trim(varname) == 'RH') then
    state%rh = var  
  else if (trim(varname) == 'PT') then
    state%pt = var  
  else if (trim(varname) == 'PTV') then
    state%ptv = var  
  else if (trim(varname) == 'Buoy') then
    state%buoy = var  
  else if (trim(varname) == 'Beff') then
    state%beff = var  
  else if (trim(varname) == 'MSE') then
    state%mse = var  
  else if (trim(varname) == 'VOF') then
    state%vof = var  
  else if (trim(varname) == 'Ent_wcore') then
    state%ent = var  
  else if (trim(varname) == 'Det_wcore') then
    state%det = var  
  else if (trim(varname) == 'Pt_var') then
    state%ptvar = var  
  else if (trim(varname) == 'Qv_var') then
    state%qvvar = var  
  else if (trim(varname) == 'W_tend_1') then
    tend%w(1) = var  
  else if (trim(varname) == 'W_tend_2') then
    tend%w(2) = var  
  else if (trim(varname) == 'W_tend_3') then
    tend%w(3) = var  
  else if (trim(varname) == 'W_tend_4') then
    tend%w(4) = var  
  else if (trim(varname) == 'W_tend_5') then
    tend%w(5) = var  
  else if (trim(varname) == 'W_tend_6') then
    tend%w(6) = var  
  else if (trim(varname) == 'W_tend_7') then
    tend%w(7) = var  
  else if (trim(varname) == 'W_tend_8') then
    tend%w(8) = var  
  else if (trim(varname) == 'W_tend_9') then
    tend%w(9) = var  
  else if (trim(varname) == 'QT_tend_1') then
    tend%qt(1) = var  
  else if (trim(varname) == 'QT_tend_2') then
    tend%qt(2) = var  
  else if (trim(varname) == 'QT_tend_3') then
    tend%qt(3) = var  
  else if (trim(varname) == 'QT_tend_4') then
    tend%qt(4) = var  
  else if (trim(varname) == 'QT_tend_5') then
    tend%qt(5) = var  
  else if (trim(varname) == 'QT_tend_6') then
    tend%qt(6) = var  
  else if (trim(varname) == 'QT_tend_7') then
    tend%qt(7) = var  
  else if (trim(varname) == 'QT_tend_8') then
    tend%qt(8) = var  
  else if (trim(varname) == 'QT_tend_9') then
    tend%qt(9) = var  
  else if (trim(varname) == 'QC_tend_1') then
    tend%qc(1) = var  
  else if (trim(varname) == 'QC_tend_2') then
    tend%qc(2) = var  
  else if (trim(varname) == 'QC_tend_3') then
    tend%qc(3) = var  
  else if (trim(varname) == 'QC_tend_4') then
    tend%qc(4) = var  
  else if (trim(varname) == 'QC_tend_5') then
    tend%qc(5) = var  
  else if (trim(varname) == 'QC_tend_6') then
    tend%qc(6) = var  
  else if (trim(varname) == 'QC_tend_7') then
    tend%qc(7) = var  
  else if (trim(varname) == 'QC_tend_8') then
    tend%qc(8) = var  
  else if (trim(varname) == 'QC_tend_9') then
    tend%qc(9) = var  
  else if (trim(varname) == 'QR_tend_1') then
    tend%qr(1) = var  
  else if (trim(varname) == 'QR_tend_2') then
    tend%qr(2) = var  
  else if (trim(varname) == 'QR_tend_3') then
    tend%qr(3) = var  
  else if (trim(varname) == 'QR_tend_4') then
    tend%qr(4) = var  
  else if (trim(varname) == 'QR_tend_5') then
    tend%qr(5) = var  
  else if (trim(varname) == 'QR_tend_6') then
    tend%qr(6) = var  
  else if (trim(varname) == 'QR_tend_7') then
    tend%qr(7) = var  
  else if (trim(varname) == 'QR_tend_8') then
    tend%qt(8) = var  
  else if (trim(varname) == 'QR_tend_9') then
    tend%qr(9) = var  
  else if (trim(varname) == 'PT_tend_1') then
    tend%pt(1) = var  
  else if (trim(varname) == 'PT_tend_2') then
    tend%pt(2) = var  
  else if (trim(varname) == 'PT_tend_3') then
    tend%pt(3) = var  
  else if (trim(varname) == 'PT_tend_4') then
    tend%pt(4) = var  
  else if (trim(varname) == 'PT_tend_5') then
    tend%pt(5) = var  
  else if (trim(varname) == 'PT_tend_6') then
    tend%pt(6) = var  
  else if (trim(varname) == 'PT_tend_7') then
    tend%pt(7) = var  
  else if (trim(varname) == 'PT_tend_8') then
    tend%pt(8) = var  
  else if (trim(varname) == 'PT_tend_9') then
    tend%pt(9) = var  
  else if (trim(varname) == 'PTV_tend_1') then
    tend%ptv(1) = var  
  else if (trim(varname) == 'PTV_tend_2') then
    tend%ptv(2) = var  
  else if (trim(varname) == 'PTV_tend_3') then
    tend%ptv(3) = var  
  else if (trim(varname) == 'PTV_tend_4') then
    tend%ptv(4) = var  
  else if (trim(varname) == 'PTV_tend_5') then
    tend%ptv(5) = var  
  else if (trim(varname) == 'PTV_tend_6') then
    tend%ptv(6) = var  
  else if (trim(varname) == 'PTV_tend_7') then
    tend%ptv(7) = var  
  else if (trim(varname) == 'PTV_tend_8') then
    tend%ptv(8) = var  
  else if (trim(varname) == 'PTV_tend_9') then
    tend%ptv(9) = var  
  else
    print*, 'Variable does not exist'
  endif
  
  end subroutine store
  
  subroutine register_var (var,nloc,ncid,dimids,varids)
  
    character(len=5) :: var
    integer, intent(inout) :: nloc
    integer, intent(in) :: ncid
    integer, dimension(:), intent(in) :: dimids
    integer, dimension(:), intent(inout) :: varids
    
    character(len=100) :: varunit,varname,varinfo

    nloc = nloc+1
    call return_info (var,varunit,varname,varinfo)
    call check( nf90_def_var(ncid, trim(varname), NF90_REAL, dimids, varids(nloc)) )
    call check( nf90_put_att(ncid, varids(nloc), UNITS, trim(varunit)) )
    call check( nf90_put_att(ncid, varids(nloc), INFO, trim(varinfo)) )
  
  contains
  
    subroutine return_info (var,units,name,info)
    
    implicit none
    
    character(len=*), intent(out), optional :: units,name,info
    character(len=5) :: var
    character(len=1) :: car1
    character(len=3) :: car3
    
    write(car3,'(a3)') var
    
    select case (trim(car3))
      case('wmo')
        if (present(name)) name = 'W momentum'
        if (present(units)) units = 'kg/m2/s'
	if (present(info)) info = 'Vertical momentum'
      case('ql')
        if (present(name)) name = 'Ql'
        if (present(units)) units = 'kg/m3'
	if (present(info)) info = 'Liquid water mixing fraction'
      case('bre')
        if (present(name)) name = 'Buoy'
        if (present(units)) units = 'm2/s3'
	if (present(info)) info = 'Resolved buoyancy'
      case('bef')
        if (present(name)) name = 'Beff'
        if (present(units)) units = 'm2/s3'
	if (present(info)) info = 'Effective buoyancy'
      case('rh')
        if (present(name)) name = 'RH'
        if (present(units)) units = '%'
	if (present(info)) info = 'Relative humidity'
      case('cid')
        if (present(name)) name = 'CLOUD_ID'
        if (present(units)) units = '#'
	if (present(info)) info = 'Cloud ID (integer)'
      case('uid')
        if (present(name)) name = 'UPDRAFT_ID'
        if (present(units)) units = '#'
	if (present(info)) info = 'Updraft ID (integer)'
    end select
    
    return
    end subroutine
    
  end subroutine  register_var
  
  subroutine check(status)
    integer, intent (in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      print *, 'AN ERROR WAS DETECTED IN NETCDF INTERFACE'
      stop
    end if
  end subroutine check  
  
end module netcdfmod
