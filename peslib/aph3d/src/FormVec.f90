SUBROUTINE FormVec(BeginVec,EndVec,NVecs,VecList,NWords)
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE
!         P U R P O S E   O F    S U B R O U T I N E
!  This Routine stores the initial and final positions along with
!  the initial and final widths into an array. A C routine is called
!  to WRITE out the vector list for plotting.
 
!         I N P U T    A R G U M E N T S
!  BeginVec   An array of points in 3 space for the starting
!               position of the vectors.
!  EndVec  An array of points for the ending positions of the
!               vectors.
!  NVecs      Number of vectors.
 
!         O U T P U T    A R G U M E N T S
!  VecList    An array with the (x,y,z) starting position and
!               width (w=1.0) followed by the ending position
!               and width.
!  NWords     Number of words in this array.
 
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 
INTEGER I, NVecs, NWords
REAL(KIND=WP_Kind) BeginVec(3,0:NVecs), EndVec(3,0:Nvecs)
REAL(KIND=WP_Kind) Veclist(NWords)
 
!
!  Store the vectors in VecList.
!
DO i = 0, NVecs - 1
   VecList(8*i+1) = BeginVec(1,i)
   VecList(8*i+2) = BeginVec(2,i)
   VecList(8*i+3) = BeginVec(3,i)
   VecList(8*i+4) = 1.0
   VecList(8*i+5) = EndVec(1,i)
   VecList(8*i+6) = EndVec(2,i)
   VecList(8*i+7) = EndVec(3,i)
   VecList(8*i+8) = 1.0
ENDDO
 
!
!  CALL C routine WriteVecs to WRITE VecList for subsequent
!    plotting. The loop increment parameter of 4200 is
!    used to fill the buffer in one CALL.
!
NWords = 8*NVecs
DO i = 1, NWords, 4200
!cccc         CALL WriteVecs (VecList, min(4200,NWords+1-i) )
ENDDO
 
ENDSUBROUTINE FORMVEC
