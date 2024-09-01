
	module cloudsmod
	
	SAVE
	
	public

        integer :: is, ie, js, je
	integer :: NXO, NYO, NZO, NTO

	character(len=10) :: object
	integer :: nc_min, nu_min, nblock
	logical :: time_clt, split_clt, stat_clt, with_up, with_halo
  	real :: dxl, dyl, h_clt, qcmin, wmin, bmin, lwpmin, ctmin

  	logical :: from_centre=.false.
	integer, parameter :: ncl=2000, ndist=100, nup=20, npdf=50
  	real, parameter :: tmin=300., dmin=2000., amin=180., smin=120.E+6, dm=2.E+6
   	real, parameter :: pi=3.14159265

	integer :: nclouds
	integer, dimension(:), allocatable :: nupdrafts
	real :: mclouds, uclouds
	
	type :: halo_type
		integer 	    :: npts		! Number of points
		integer 	    :: i(100)		! i-index
		integer 	    :: j(100)		! j-index
	end type halo_type
	
	type :: cloud_type
		real 		    :: mfl, qfl		! Mass flux and water mass flux
		real		    :: mom              ! Cloud averaged momentum
		real		    :: dmdz             ! Cloud averaged convective tendency
		real 		    :: mse, qt		! moist static energy & total water
		real		    :: buo, bef		! cloud buoyancy 
		real		    :: det,ent		! cloud entrainment rates
		real		    :: lwp,ctop		! LWP and cloud top
		real		    :: age		! Life time
		real		    :: sf		! Shape factor (perimeter/4*radius)
		real 		    :: dist(100)	! Distance between cloud centres
		real 		    :: angle(100)	! Angle between cloud centres
		integer             :: stat		! Cloud status (passive, single updraft, multiple updrafts)
		integer 	    :: id		! Cloud id
		integer 	    :: npts		! Number of points
		integer 	    :: size		! Size in grid cells
		integer 	    :: i(100)		! i-index
		integer 	    :: j(100)		! j-index
		integer		    :: ic, jc		! Cloud centre index
		type(halo_type)     :: halo(8)          ! halos
	end type cloud_type

	interface assignment (=)
		module procedure const2c
		module procedure c2c
	end interface

	CONTAINS
	
	  subroutine const2c (c1, const)
	       type (cloud_type), intent(inout) :: c1
	       real, intent(in) :: const
	       integer :: n
	       
	       c1%mfl  	= const
	       c1%qfl  	= const
	       c1%dmdz 	= const
	       c1%mom 	= const
	       c1%qt  	= const
	       c1%mse  	= const
	       c1%buo  	= const
	       c1%bef  	= const
	       c1%det  	= const
	       c1%ent  	= const
	       c1%lwp  	= const
	       c1%ctop 	= const
	       c1%age   = const
	       c1%dist 	= const
	       c1%angle = const
	       c1%npts 	= int(const)
	       c1%size 	= int(const)
	       c1%stat	= int(const)
	       c1%id	= int(const)
	       c1%i	= int(const)
	       c1%j	= int(const)
	       c1%ic	= int(const)
	       c1%jc	= int(const)
	  return
	  end subroutine const2c
	
	  subroutine c2c (c1, c2)
	       type (cloud_type), intent(inout) :: c1
	       type (cloud_type), intent(in)    :: c2
	       
	       c1%mfl   = c2%mfl
	       c1%qfl   = c2%qfl
	       c1%dmdz  = c2%dmdz
	       c1%mom   = c2%mom
	       c1%qt    = c2%qt
	       c1%mse   = c2%mse
	       c1%buo   = c2%buo
	       c1%bef   = c2%bef
	       c1%det   = c2%det
	       c1%ent   = c2%ent
	       c1%lwp   = c2%lwp
	       c1%ctop  = c2%ctop
	       c1%age   = c2%age
	       c1%dist  = c2%dist
	       c1%angle = c2%angle
	       c1%npts  = c2%npts
	       c1%size  = c2%size
	       c1%stat  = c2%stat
	       c1%id    = c2%id
	       c1%i     = c2%i
	       c1%j     = c2%j
	       c1%ic    = c2%ic
	       c1%jc    = c2%jc
	  return
	  end subroutine c2c
   
! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

    SUBROUTINE  FindMinimum(x, Start, End, Location)
      IMPLICIT  NONE
      INTEGER                :: Start, End
      INTEGER                :: Location
      INTEGER                :: i
      REAL                   :: Minimum
      REAL, DIMENSION(1:End) :: x

      Minimum  = x(Start)               ! assume the first is the min
      Location = Start                  ! record its position
      DO i = Start+1, End               ! start with next elements
         IF (x(i) < Minimum) THEN       !   if x(i) less than the min?
            Minimum  = x(i)             !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO

    END SUBROUTINE  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

    SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      REAL :: a, b
      REAL :: Temp

      Temp = a
      a    = b
      b    = Temp

    END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

    SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      INTEGER                :: Size
      INTEGER                :: i
      INTEGER                :: Location
      REAL, DIMENSION(1:Size) :: x

      DO i = 1, Size-1                  ! except for the last
         CALL FindMinimum (x, i, Size, Location)        ! find min from this to last
         CALL Swap (x(i), x(Location))  ! swap this and the minimum
      END DO

    END SUBROUTINE  Sort

	end module cloudsmod


