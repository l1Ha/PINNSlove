SUBROUTINE SetBasis (nang, nbasis, maxn, maxl, nangch, nbasiss)
!

!            P U R P O S E   O F    S U B R O U T I N E
! This routine determines the number of APH theta and chi points and
! set up the quantum numbers for the basis.


!            O U T P U T    A R G U M E N T S
! nang     number of APH theta and chi angles.
! nbasis   number of basis functions.
! maxn     maximum vibrational quantum number.
! maxl     maximum rotational quantum number.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

USE Parms_Module
USE TotJ_Module
USE FileUnits_Module
USE Quantb_Module
USE QDiv_Module
USE QCase_Module
USE Gaussb_Module
IMPLICIT NONE

!            I N T E G E R S
INTEGER ibasis, rot, vib, karran, nang, jbasis, nth, nch, ithdiv
INTEGER ichdiv, maxn, maxl
INTEGER nangch(narran+1),nbasis(0:mxmega), nbasiss(0:mxmega)
INTEGER i,maxtemp
!            E X T E R N A L S
EXTERNAL popt

!            I N T R I N S I C S
INTRINSIC mod

!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('setbasis', little, medium, full, ithcall, ithsub)

ibasis = 0
jbasis = 0
maxn = 0
maxl = 0
nangch(1) = 1
DO karran = 1, narran
!-----------------------------------------------------------------------
! Check DATA to make sure that dimension will not be exceeded.
!-----------------------------------------------------------------------
   IF(nhermt(karran)>mxhermt)THEN
      WRITE(Out_Unit,*)'Error: nhermt(karran)>mxhermt', nhermt(karran),mxhermt
      STOP 'setbasis1'
   ENDIF
   IF(nglegn(karran)>mxglegn)THEN
      WRITE(Out_Unit,*)'Error: nglegn(karran)>mxglegn',  nglegn(karran),mxglegn
      STOP 'setbasis2'
   ENDIF
   IF(noscil(karran)>nhermt(karran))THEN
      WRITE(Out_Unit,*)'Error: noscil(karran)>nhermt(karran) ', noscil(karran),nhermt(karran)
      STOP 'setbasis3'
   ENDIF
   IF(nlegndre(karran)>nglegn(karran))THEN
      WRITE(Out_Unit,*)'Error: ', 'nlegndre(karran)>nglegn(karran) ',  nlegndre(karran),nglegn(karran)
      STOP 'setbasis4'
   ENDIF

!-----------------------------------------------------------------------
! Determine the number of points for numberical quadrature
!-----------------------------------------------------------------------
   IF(.not.qcase)THEN
      IF(symmetry)THEN
         IF(karran==1)THEN
            nangch(karran+1) = nangch(karran) + nhermt(karran)*(nglegn(karran)/2)
         ELSEIF(karran==2)THEN
            nangch(karran+1) = nangch(karran) + nhermt(karran)*nglegn(karran) !GregParker TEMPMOD ?
            !nangch(karran+1) = nangch(karran) + nhermt(karran)*(nglegn(karran)/2)
         ELSE
           nangch(karran+1) = nangch(karran) ! GregParker  TEMPMOD ?
           !nangch(karran+1) = nangch(karran) + nhermt(karran)*(nglegn(karran)/2)
         ENDIF
      ELSE
         nangch(karran+1) = nangch(karran) + nhermt(karran)*nglegn(karran)
      ENDIF
   ENDIF

   DO vib = minvib(karran), maxvib(karran)
!  ---------------------------------------------------------------
! Determine basis
!  determine jmin for arbitrary mega with account of symmetry
!  ---------------------------------------------------------------
      jmin(vib,karran,mega)=max(jmin(vib,karran,0), mega)
      IF(symmetry.AND.karran==1)THEN
         jskip(1) = 2
         IF(jeven)THEN
            jmin(vib,karran,mega) = jmin(vib,karran,mega)  + mod(jmin(vib,karran,mega),2)
         ELSE
            jmin(vib,karran,mega) = jmin(vib,karran,mega)  - mod(jmin(vib,karran,mega),2) + 1
         ENDIF
      ELSE
         jskip(karran) = 1
      ENDIF
!  ---------------------------------------------------------------
      DO rot = jmin(vib,karran,mega), jmax(vib,karran), jskip(karran)
         ibasis = ibasis +1
         IF(karran/=3)jbasis = jbasis + 1
         !write(Msg_Unit,*)ibasis, karran, vib, rot, rot, jbasis
         chanl(ibasis,mega) = karran
         nvib(ibasis)  = vib
         jrot(ibasis,mega)  = rot
         IF(vib+1>maxn) maxn = vib+1
         IF(rot+1>maxl) maxl = rot+1
      ENDDO
   ENDDO
ENDDO
!
!  determine nang.
IF(qcase)THEN
   nth = 0
   DO ithdiv = 1, nthdiv
      nth = nth + nqth(ithdiv)
   ENDDO
   nch = 0
   DO ichdiv = 1, nchdiv
      nch = nch + nqch(ichdiv)
   ENDDO

   nang = nth*nch

   IF(symmetry)THEN
      IF(mod(nch,2)/=0)THEN
         WRITE(Out_Unit,*)'nch must be even for symmetry true'
         STOP 'setbasis5'
      ENDIF
      nang=nang/2
   ENDIF
   nangch(1) = 1
   nangch(2) = 1 + nang
   nangch(3) = 1 + nang
   nangch(4) = 1 + nang
   WRITE(Out_Unit,*)'nang,nth,nch=',nang,nth,nch
ELSE
   nang = nangch(narran+1) - 1
ENDIF
WRITE(Out_Unit,*)'nangch=',nangch
mxang=nang

maxtemp=0
DO i=1,narran
   IF(i==1.AND.symmetry)THEN
      maxtemp=maxtemp+nhermt(i)*nglegn(i)/2
   ELSE
      maxtemp=maxtemp+nhermt(i)*nglegn(i)
   ENDIF
ENDDO

mxang=max(nang,maxtemp)

!  --------------------------------------------------------------
!  nbasis is the total no. of raw primitive basis functions
!  nbasiss is the number of symmetry-adapted primitive basis fcns.
!  nbasis-nbasiss is the number of basis fcns in channel 2, IF
!  symmetry.
!  --------------------------------------------------------------

nbasis(mega) = ibasis
nbasiss(mega) = ibasis
IF(symmetry) nbasiss(mega)= jbasis

!  ------------------------------------------------------------------
!  check to see that current parameters DO not exceed dimensions.
!  ------------------------------------------------------------------

IF(nbasis(mega)>mxbasis)THEN
   WRITE(Out_Unit,*)'Error: ', 'nbasis>mxbasis ',nbasis(mega),mxbasis
   STOP 'setbasis6'
ENDIF
IF(nbasiss(mega)>mxbasiss)THEN
   WRITE(Out_Unit,*)'Error: ', 'nbasiss>mxbasiss',nbasiss(mega),mxbasiss
   STOP 'setbasis7'
ENDIF
IF(maxn>mxn)THEN
   WRITE(Out_Unit,*)'Error: maxn>mxn ',maxn,mxn
   STOP 'setbasis8'
ENDIF
IF(maxl>mxl)THEN
   WRITE(Out_Unit,*)'Error: maxl>mxl ',maxl,mxl
   STOP 'setbasis9'
ENDIF
IF(nang>mxang)THEN
   WRITE(Out_Unit,*)'Error: nang>mxang ',nang,mxang
   STOP 'setbasis10'
ENDIF
IF(megamax>mxmega)THEN
   WRITE(Out_Unit,*)'Error: megamax>mxmega',megamax,mxmega
   STOP 'setbasis11'
ENDIF


IF(little)THEN
   WRITE(Out_Unit,*)'routine setbasis'
   WRITE(Out_Unit,*)'nang,nbasis(mega),maxn,maxl,nangch= ', nang,nbasis(mega),maxn,maxl, (nangch(karran), karran = 1, narran +1)
   WRITE(Out_Unit,*)'megamin, megamax=',megamin, megamax
   IF(symmetry) WRITE(Out_Unit,*) 'nbasiss(mega)= ',nbasiss(mega)
ENDIF
IF(medium)THEN
   WRITE(Out_Unit,*)' ibasis, chanl, nvib, jrot'
   DO ibasis = 1, nbasis(mega)
      WRITE(Out_Unit,*)ibasis, chanl(ibasis,mega),nvib(ibasis),  jrot(ibasis,mega)
   ENDDO
ENDIF

RETURN
ENDSUBROUTINE SetBasis
