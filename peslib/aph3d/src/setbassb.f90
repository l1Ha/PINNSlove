SUBROUTINE setbassb (nbasis, ngleg, nherm, nbasiss)

!            P U R P O S E   O F    S U B R O U T I N E
! This routine determines the number of angles and
! sets up the quantum numbers for the basis.
!  This version is for the AphDel_Old or AphDel_New transformation.


!            O U T P U T    A R G U M E N T S
! nang     number of APH theta and chi angles.
! nbasis   number of basis functions.
! maxn     maximum vibrational quantum number.
! maxj     maximum rotational quantum number.
!this routine is called by:
!         intbas
!this routine calls
!         popt
!-----------------------------------------------------------------------

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>



USE FileUnits_Module

!            I N C L U D E    S T A T E M E N T S
USE Narran_Module
USE gaussb_Module
USE quantb_Module
USE totj_Module
IMPLICIT NONE

!            I N T E G E R S
INTEGER ibasis, rot, vib, maxn, maxj, karran, nang, jbasis,nangch(narran+1), &
        nherm(narran), ngleg(narran), nbasis(0:mxmega), nbasiss(0:mxmega)


!            E X T E R N A L S
EXTERNAL popt

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('setbassb', little, medium, full, ithcall, ithsub)

!---------------------------------------------------------------------
!  set some initial parameters.
!---------------------------------------------------------------------
WRITE(Out_Unit,*)'begining of routine setbassb'

maxn = 0
maxj = 0
nangch(1) = 1
megamin=0
IF((-1)**(jtot+parity)/=1) megamin=1
DO  mega = megamin, megamax
ibasis = 0
jbasis = 0
DO  karran = 1, narran
!-----------------------------------------------------------------------
!Check data to make sure that dimension will not be exceeded.
!-----------------------------------------------------------------------
   IF(nherm(karran)>maxhermt)THEN  
      WRITE(*,*)'Error: nherm(karran)>maxhermt', nherm(karran),maxhermt
      WRITE(Out_Unit,*)'Error: nherm(karran)>maxhermt', nherm(karran),maxhermt
      STOP 'setbassb loc1'
   ENDIF
   IF(ngleg(karran)>mxglegn)THEN
      WRITE(*,*)'Error: ngleg(karran)>mxglegn', ngleg(karran),mxglegn
      WRITE(Out_Unit,*)'Error: ngleg(karran)>mxglegn', ngleg(karran),mxglegn
      STOP 'setbassb loc2'
   ENDIF
   IF(noscil(karran)>nherm(karran))THEN
      WRITE(*,*)'Warning: ', 'noscil(karran)>nherm(karran) ', noscil(karran),nherm(karran)
      WRITE(Out_Unit,*)'Warning: ', 'noscil(karran)>nherm(karran) ', noscil(karran),nherm(karran)
!            STOP 'setbassb loc3'
   ENDIF
   IF(nlegndre(karran)>ngleg(karran))THEN
      WRITE(*,*)'Warning: ', 'nlegndre(karran)>ngleg(karran) ', nlegndre(karran),ngleg(karran)
      WRITE(Out_Unit,*)'Warning: ', 'nlegndre(karran)>ngleg(karran) ', nlegndre(karran),ngleg(karran)
!            STOP 'setbassb loc4'
   ENDIF

!-----------------------------------------------------------------------
! Determine basis
!  note.  AphDel_Old or AphDel_New does not use symmetry and uses all angles, so this must
!  construct functions at all angles.
!-----------------------------------------------------------------------
   nangch(karran+1) = nangch(karran)+ nherm(karran)*ngleg(karran)

   DO  vib = minvib(karran), maxvib(karran)
!  ---------------------------------------------------------------
!  determine jmin for arbitrary mega with account of symmetry
!  ---------------------------------------------------------------
      jmin(vib,karran,mega)=max(jmin(vib,karran,0), mega)
      IF(symmetry.AND.karran==1)THEN
         jskip(1) = 2
         IF(jeven)THEN
            jmin(vib,karran,mega) = jmin(vib,karran,mega) + mod(jmin(vib,karran,mega),2)
         ELSE
            jmin(vib,karran,mega) = jmin(vib,karran,mega) - mod(jmin(vib,karran,mega),2) + 1
         ENDIF
      ELSE
         jskip(karran) = 1
      ENDIF
!  ---------------------------------------------------------------
      DO  rot = jmin(vib,karran,mega), jmax(vib,karran), jskip(karran)
         ibasis = ibasis +1
     IF(karran/=3)jbasis = jbasis + 1
         chanl(ibasis, mega) = karran
         nvib(ibasis)  = vib
         jrot(ibasis,mega)  = rot
         IF(vib>maxn) maxn = vib
         IF(rot>maxj) maxj = rot
    ENDDO
  ENDDO
 ENDDO
  nang = nangch(narran+1) - 1

!  --------------------------------------------------------------
!  nbasis is the total no. of raw primitive basis functions
!  nbasiss is the number of symmetry-adapted primitive basis fcns.
!  nbasis-nbasiss is the number of basis fcns in channel 2, if
!  symmetry.
!  --------------------------------------------------------------

nbasis(mega) = ibasis
nbasiss(mega) = ibasis
IF(symmetry) nbasiss(mega)= jbasis

!------------------------------------------------------------------
!  check to see that current parameters do not exceed dimensions.
!------------------------------------------------------------------

IF(nbasis(mega)>mxbasis)THEN
   WRITE(*,*)'Error: nbasis>mxbasis ', nbasis(mega),mxbasis
   WRITE(Out_Unit,*)'Error: nbasis>mxbasis ', nbasis(mega),mxbasis
   STOP 'setbassb  loc5'
ENDIF
IF(nbasiss(mega)>mxbasiss)THEN
   WRITE(*,*)'Error: nbasiss>mxbasiss', nbasiss(mega),mxbasiss
   WRITE(Out_Unit,*)'Error: nbasiss>mxbasiss', nbasiss(mega),mxbasiss
   STOP 'setbassb  loc6'
ENDIF


IF(maxn>mxn)THEN
   WRITE(*,*)'Error: maxn>mxn ',maxn,mxn
   WRITE(Out_Unit,*)'Error: maxn>mxn ',maxn,mxn
   STOP 'setbassb  loc7'
ENDIF
IF(maxj>mxl)THEN
   WRITE(*,*)'Error: maxj>mxl ',maxj,mxl
   WRITE(Out_Unit,*)'Error: maxj>mxl ',maxj,mxl
   STOP 'setbassb  loc8'
ENDIF
IF(nang>mxang)THEN
   WRITE(*,*)'Error: nang>mxang ',nang,mxang
   WRITE(Out_Unit,*)'Error: nang>mxang ',nang,mxang
   STOP 'setbassb loc9'
ENDIF
IF(megamax>mxmega)THEN
   WRITE(*,*)'Error: megamax>mxmega',megamax,mxmega
   WRITE(Out_Unit,*)'Error: megamax>mxmega',megamax,mxmega
   STOP 'setbassb loc10'
ENDIF

IF(little)THEN
   WRITE(Out_Unit,*)'routine setbassb'
   WRITE(Out_Unit,*)'nang,nbasis(mega),maxn,maxj,nangch= ', nang,nbasis(mega),maxn,maxj, (nangch(karran), karran = 1, narran +1)
   WRITE(Out_Unit,*)'megamin, megamax=',megamin, megamax
   IF(symmetry) WRITE(Out_Unit,*) 'nbasiss(mega)= ',nbasiss(mega)
ENDIF

IF(medium)THEN
   WRITE(Out_Unit,*)' ibasis, chanl, nvib, jrot'
   DO  ibasis = 1, nbasis(mega)
      WRITE(Out_Unit,*)mega,ibasis,chanl(ibasis, mega), nvib(ibasis), jrot(ibasis,mega)
   ENDDO
ENDIF
    ENDDO
 WRITE(Out_Unit,*)'END of routine setbassb'

RETURN
END
