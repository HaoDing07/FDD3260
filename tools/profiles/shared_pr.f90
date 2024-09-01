
  module sharedmod

    type :: winds
   	    real :: u		    ! Velocity in x 
   	    real :: v		    ! Velocity in y
   	    real :: w		    ! Velocity in z
    end type winds

    type :: states
   	    real :: qc, qr	    ! Droplet and rain mixing ratios
   	    real :: qv, rh	    ! Water vapor mixing ratio
   	    real :: pt, mse 	    ! potential temperature and MSE
   	    real :: buoy, beff	    ! Buoyancy
	    real :: evap	    ! Evaporation
	    real :: cent, cdet	    ! Entrainment
	    real :: sent, sdet	    ! Entrainment
    end type states

    integer, parameter :: nvarin=16, nts=6
    integer :: NXO, NYO, NZO, NTO, NDIMS, NVAR
    
    real, dimension(:), allocatable :: xx, yy, zz, time
    type(winds), dimension(:,:,:,:), allocatable :: wind
    type(states), dimension(:,:,:,:), allocatable :: state

    save wind, state, xx, yy, zz, time

  end module sharedmod
