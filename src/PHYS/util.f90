Module util

    Contains  

!Taken from "util"

FUNCTION closest(array,val)
    ! Find the index of "array" with value closest to "val"
    IMPLICIT NONE
    INTEGER :: closest
    REAL, INTENT(in) :: array(:)
    REAL, INTENT(in) :: val
    
    REAL, ALLOCATABLE :: diff(:)
    REAL :: mindiff
    
    INTEGER :: NN, i
    
    NN = SIZE(array)
    ALLOCATE(diff(NN))
    
    diff = ABS(array-val)
    mindiff = MINVAL(diff)
    
    DO i = 1,NN
       IF (diff(i) == mindiff) THEN
          closest = i
          EXIT
       END IF
    END DO
    
    closest = MAX(MIN(closest,NN),1)
    DEALLOCATE(diff)
    
  END FUNCTION closest
 
 ! --------------------------------------------------------------------------- 
  ! For Level >= 4: Returns the index for mass mixing ratio in the prognostic
  ! tracer arrays for a given SALSA size bin and a given aerosol species (or 
  ! water)
  ! Input arguments: nbtot = total number of bins
  !                  nb    = number of the bin for which the index is fetched
  !                  nm    = index of the aerosol species (use classSpecies for this)
  !
  INTEGER FUNCTION getMassIndex(nbtot,nb,nm)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nbtot, nb, nm
  
  getMassIndex = (nm-1)*nbtot+nb
  
 END FUNCTION getMassIndex

end module util
 