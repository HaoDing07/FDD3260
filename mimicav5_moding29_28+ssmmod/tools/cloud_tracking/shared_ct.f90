
  module sharedmod
  
    SAVE

    type :: winds
   	    real :: u		    ! Velocity in x 
   	    real :: v		    ! Velocity in y
   	    real :: w		    ! Velocity in z
    end type winds

    type :: states
   	    real :: qc, qr	    ! Droplet and rain mixing ratios
   	    real :: qv, rh	    ! Water vapor mixing ratio
   	    real :: buoy, beff	    ! Buoyancies
	    real :: dens	    ! Density
   	    real :: pt  	    ! Potential temperature
   	    real :: ptv 	    ! Virtual potential temperature
   	    real :: mse 	    ! Moist static energy
   	    real :: lwp	   	    ! Liquid water path
   	    real :: ctop	    ! Cloud top height
   	    real :: det,ent	    ! entrainment/detrainment rates
    end type states

    type :: tendencies
	    real :: dqcon             ! qt convection (to diagnose entrainment)
	    real :: dqadv             ! qt convection (to diagnose entrainment)
	    real :: dqdif             ! qt convection (to diagnose entrainment)
	    real :: dqfor             ! qt convection (to diagnose entrainment)
	    real :: dqdt              ! qt convection (to diagnose entrainment)
    end type tendencies

    real, dimension(:), allocatable :: xx, yy, zz, time
    type(winds), dimension(:,:,:,:), allocatable :: wind
    type(states), dimension(:,:,:,:), allocatable :: state
    type(tendencies), dimension(:,:,:,:), allocatable :: tend

  end module sharedmod
