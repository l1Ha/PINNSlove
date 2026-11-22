SUBROUTINE FndInt (n, array, value, start, FindLoc)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
! >>> FndInt -- finds the interval where value lies in the array
!
! --- on entry the calling routine supplies --
!
!     n     -- the number of elements in the array
!     array -- a monotonically nondecreasing array of floating point
!                numbers (i.e., increasing array )
!     value -- test value to find in the array
!     start -- INTEGER initial guess where value is in the array
!
! --- on exit the routine returns --
!
!     FindLoc  -- INTEGER value locating where value is in the array, so
!                that array(FindLoc) <= value < array(FindLoc+1), and
!                if value < array(1), then FindLoc = 0, and
!                if value >= array(n), then FindLoc = n
!
! --- Method:  routine first checks to see if the test value is outside
!                the array bounds (FindLoc=0 or n), then checks
!                the current interval and adjacent intervals,
!                where FindLoc=start, start-1, or start+1.  If correct
!                interval is not located by this point, routine
!                branches into a binary search.  Because of the
!                interval end conditions, it takes
!                   1 test to locate FindLoc=0
!                   2 tests to locate FindLoc=n
!                   4 tests to locate FindLoc=start, or start-1, and
!                   5 tests to locate FindLoc=start+1
INTEGER :: high, low, FindLoc, start, Look, N
REAL(Kind=WP_Kind) :: array (n), value
!
high = n
low = 1
!
!     CHECK SPECIAL CASES
!
!     -- check out of bounds conditions
!
!      ntests = 1
IF(value<array (1) )THEN
   FindLoc = 0
   RETURN
ELSEIF (value>=array (n) )THEN
   FindLoc = n
!         ntests = 2
   RETURN
ENDIF
!
!    -- check intervals next to the start search
!
look = min0 (n, max0 (1, start) )
!      ntests = 4
IF(value<array (look) )THEN
   high = look
   IF(value<array (high - 1) )THEN
      high = high - 1
   ELSE
      FindLoc = high - 1
      RETURN
   ENDIF
ELSE
   low = look
   IF(value<array (low + 1) )THEN
      FindLoc = low
      RETURN
   ELSE
      low = low + 1
!            ntests = 5
      IF(value<array (low + 1) )THEN
         FindLoc = low
         RETURN
      ELSE
         low = low + 1
      ENDIF
   ENDIF
ENDIF
!
!     NOT A SIMPLE ANSWER, BEGIN BINARY SEARCH
!
look = (high + low) / 2
look = min0 (high - 1, look)
look = max0 (low + 1, look)
    1 CONTINUE
IF(high - low==1)THEN
   FindLoc = low
   RETURN
ELSE
!         ntests = ntests + 1
   IF(value<array (look) )THEN
      high = look
   ELSE
      low = look
   ENDIF
!
   look = (high + low) / 2
   look = min0 (high - 1, look)
   look = max0 (low + 1, look)
ENDIF
goto 1
ENDSUBROUTINE FndInt
