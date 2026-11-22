SUBROUTINE Figure (j_minimum,jmax,mvib,jrot,lorb,megax,negvl,minvib ,maxvib)
USE Numeric_Kinds_Module
!-----------------------------------------------------------------------
!this routine is called by:
!         qlevel
!this routine calls
!         popt
!-------------------------------------------------------------------
! j_minimum(i)/JMAX(i) have the min/max rotational quantum number for
! the nu vibrational state, indexed by i=nu+1; not altered on return
! MINVIB is the minimum value of nu+1; not altered
! MAXVIB is the maximum value of nu+1; not altered
! JTOT is the total angular momentum; not altered
!
! The following variables are loaded in this routine:
!
! NEGVL is the number of vibrational-rotational states
! in this channel
!
! The following variables are indexed from 1 to negvl:
!
! JROT has the rotational quantum number
! MVIB has the vibrational quantum number, nu
! LORB has the orbital quantum number
! MEGA is the total angualr momentum projection quantum number
!-------------------------------------------------------------------
USE FileUnits_Module
USE Narran_Module
USE melm_Module
USE totj_Module
USE elemnt_Module
USE approx_Module
USE NState_Module
INTEGER I, IthCll, IthSub, JHigh, Jlitl, JLow, KVib, MaxVib, Meg
INTEGER MegMax, MinVib, Negvl
INTEGER j_minimum(0:*), jmax(0:*), mvib(*), lorb(*), jrot(*), megax(*)
LOGICAL brnch1, brnch2
LOGICAL little, medium, full
!-------------------------------------------------------------------
! determines the body-fixed quantum numbers.
!-------------------------------------------------------------------
!      DATA mbar/0/
DATA ithcll/0/,ithsub/0/
DATA little/.false./,medium/.false./,full/.false./
CALL popt('figure  ',little,medium,full,ithcll,ithsub)
!-----------------------------------------------------------------------
mbar=0
IF(little)THEN
WRITE(Out_Unit,*)'symmetry=',symmetry,' jeven=',jeven,' npar=',npar
WRITE(Out_Unit,*)'pseudo=',pseudo,' censud=',censud,' engsud=',engsud
WRITE(Out_Unit,*)'mbar=',mbar,' melmnt=',melmnt
WRITE(Out_Unit,*)'minvib=',minvib,' maxvib=',maxvib
ENDIF
!-----------------------------------------------------------------------
! loop over the vibrations.
!-----------------------------------------------------------------------
negvl=0
DO 12 kvib=minvib,maxvib
   jhigh=jmax(kvib)
   jlow=j_minimum(kvib)
   IF(melmnt=='molbf')THEN
      jhigh=jtot+jmax(kvib)
      jlow=0
   ENDIF
   IF(medium)WRITE(Out_Unit,*)'jlow=',jlow,' jhigh=',jhigh
!-----------------------------------------------------------------------
! Loop over the rotations.
! determine rotational quantum numbers for each vibrational state
!-----------------------------------------------------------------------
   DO 11 jlitl=jlow,jhigh
      megmax=min0(jlitl,jtot)
      DO 10 meg=megmax,-megmax,-1
         IF(full) WRITE(Out_Unit,*)'meg=',meg,' negvl=',negvl
         IF(.NOT.pseudo.AND.censud.AND.(meg/=mbar)) GOTO 10
         IF(.NOT.pseudo.AND.engsud.AND.(meg/=mbar)) GOTO 10
!-----------------------------------------------------------------------
! Check for even or odd momentum
!-----------------------------------------------------------------------
         IF(((jlitl+1)/2)/=(jlitl/2))THEN
            brnch1=.false.
         ELSE
            brnch1=.true.
         ENDIF
!-----------------------------------------------------------------------
! Check for positive or negative projections.
!-----------------------------------------------------------------------
         IF(meg>=0)THEN
            brnch2=.true.
         ELSE
            brnch2=.false.
         ENDIF
         IF(symmetry.AND.jeven.AND..NOT.brnch1) GOTO 10
         IF(symmetry.AND..NOT.jeven.AND.brnch1) GOTO 10
!-----------------------------------------------------------------------
! Check to see if the Total angular momentum is even or odd.
!-----------------------------------------------------------------------
         IF(2*(jtot/2)==jtot)THEN
            IF(npar.AND.parity==0.AND..NOT.brnch2) GOTO 10
            IF(npar.AND.parity==1.AND.brnch2) GOTO 10
         ELSE
            IF(npar.AND.parity==0.AND.brnch2) GOTO 10
            IF(npar.AND.parity==1.AND..NOT.brnch2) GOTO 10
         ENDIF

         negvl=negvl+1
         IF(negvl>nstate)THEN
           WRITE(Out_Unit,*)' negvl>nstate:',negvl,nstate
           WRITE(Msg_Unit,*)' negvl>nstate:',negvl,nstate
           STOP 'Figure: negvl>nstate:'
         ENDIF
         mvib(negvl)=kvib
         megax(negvl)=meg
         IF(melmnt/='molbf')THEN
            jrot(negvl)=jlitl
            IF(meg>=0) lorb(negvl)=jtot+jrot(negvl)-2*meg
            IF(meg<0) lorb(negvl)=jtot+jrot(negvl)+2*meg+1
         ELSE
            lorb(negvl)=jlitl
            IF(meg>=0) jrot(negvl)=jtot+lorb(negvl)-2*meg
            IF(meg<0) jrot(negvl)=jtot+lorb(negvl)+2*meg+1
         ENDIF
   10       CONTINUE
   11    CONTINUE
   12 CONTINUE
!-----------------------------------------------------------------------
! If there are zero states STOP the calculation.
!-----------------------------------------------------------------------
IF(negvl==0)THEN
   WRITE(Out_Unit,30)
   STOP 'Figure: negvl==0'
ENDIF
!-----------------------------------------------------------------------
! If desired print out quantum numbers.
!-----------------------------------------------------------------------
IF(little.or.medium.or.full)THEN
   WRITE(Out_Unit,'(5A7)')'i','mvib','jrot','lorb','megax'
   DO i=1,negvl
      WRITE(Out_Unit,'(5I7)') i,mvib(i),jrot(i),lorb(i),megax(i)
   ENDDO
ENDIF
!-----------------------------------------------------------------------
!  return to routine qlevel.
!-----------------------------------------------------------------------
RETURN

   30 FORMAT(/,' negvl is equal to zero in figure')
END SUBROUTINE Figure
