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
      REAL   	             :: Minimum
      REAL, DIMENSION(1:End) :: x

      Minimum  = x(Start)		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, End		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
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

      DO i = 1, Size-1			! except for the last
         CALL FindMinimum (x, i, Size, Location)	! find min from this to last
         CALL Swap (x(i), x(Location))	! swap this and the minimum
      END DO
      
   END SUBROUTINE  Sort

