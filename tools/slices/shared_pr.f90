
  module sharedmod

    type :: winds
   	    real :: u		    ! Velocity in x 
   	    real :: v		    ! Velocity in y
   	    real :: w		    ! Velocity in z
    end type winds

    type :: states
   	    real :: qc, qr	    ! Droplet and rain mixing ratios
   	    real :: qv, rh	    ! Water vapor mixing ratio
   	    real :: pt, ptv 	    ! Virtual potential temperature
   	    real :: buoy, beff	    ! Buoyancy
	    real :: ptvar, qvvar    ! variances
	    real :: mse		    ! MSE
	    real :: wvp		    ! Water Vapour Path
	    real :: cth		    ! Cloud top height
    end type states

    type :: tendencies
   	    real, dimension(9) :: qt, pt, w
    end type tendencies

    integer, parameter :: nvarin=43
    integer, parameter :: nvarslic=30
    integer :: NXO, NYO, NZO, NTO, NDIMS, NVAR, ntime
    
    real :: dx
    real, dimension(:), allocatable :: xx, yy, zz, tt
    
    type(winds), dimension(:,:,:,:), allocatable :: wind
    type(states), dimension(:,:,:,:), allocatable :: state
    type(tendencies), dimension(:,:,:,:), allocatable :: tend

    save wind, state, tend, xx, yy, zz, tt

  end module sharedmod
