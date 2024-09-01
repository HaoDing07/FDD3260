
  module sharedmod

    type :: winds
   	    real :: u		    ! Velocity in x 
   	    real :: v		    ! Velocity in y
   	    real :: w		    ! Velocity in z
   	    real :: div		    ! Horizontal divergence
    end type winds

    type :: states
   	    real :: qc, qr	    ! Droplet and rain mixing ratios
   	    real :: qv, qt, rh	    ! Water vapor mixing ratio
   	    real :: pt, ptv 	    ! Virtual potential temperature
   	    real :: buoy, beff	    ! Buoyancy
	    real :: ptvar, qvvar    ! variances
	    real :: p, mse	    ! MSE
	    real :: wvp		    ! Water Vapour Path
	    real :: cth		    ! Cloud top height
	    real :: ent, det, vof   ! Entrainment/detrainment
    end type states

    type :: tendencies
   	    real, dimension(9) :: qt, qc, qr, pt, ptv, w
    end type tendencies

    integer, parameter :: nvarin=67
    integer, parameter :: nvarslic=41
    integer :: NXO, NYO, NZO, NTO, NDIMS, NVAR, ntime
    
    logical :: norm
    
    real :: dx
    real, dimension(:), allocatable :: xx, yy, zz, tt
    real :: base_min=1000., base_max=2500., depth_min=250., top_max=10000., dist_max=3000.
    
    type(winds), dimension(:,:,:,:), allocatable :: wind
    type(states), dimension(:,:,:,:), allocatable :: state
    type(tendencies), dimension(:,:,:,:), allocatable :: tend

    save wind, state, tend, xx, yy, zz, tt

  end module sharedmod
