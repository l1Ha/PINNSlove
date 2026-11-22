SUBROUTINE OvrBas (rholast, rho, nang, nbasis, xchanl, nangch, philast, phinow, vectlast, vect, s,    &
                   amatrx, intwt, ithrho, nbasiss,  neigmin, iunit1, jrot)
!

!            P U R P O S E   O F    S U B R O U T I N E
! This routine determines the overlap matrix elements between surface
! functions at different rho values.


!            I N P U T    A R G U M E N T S
!   rholast  previous hyperradius.
!   rho      current hyperradius
!   nang     number of quadrature angles.
!   nbasis   number of raw primitive basis functions.
!   nbasiss  number of symmetry-adapted primitive basis functions.
!   narran   number of arrangement channels.
!   xchanl    an array giving the channel number of each basis
!            function.
!   nangch   an array containing the starting position of the
!            angles for each of the arrangement channels.
!            nangch(narran+1) is equal to nang+1.
!   philast  primitive basis of previous rho at present quadrature pts.
!   phinow   primitive basis of present rho at present quadrature pts.
!   vectlast coeffs of surface fcns of previous rho.
!   vect     coeffs of present surface functions.

!            O U T P U T    A R G U M E N T S
!   s        overlap matrix elements.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE Parms_Module
USE TotJ_Module
USE FileUnits_Module
USE QCase_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER nang, ibasis, jbasis, nbasis, iang, iangst1, iangend1
INTEGER iangst2, iangend2, jpar
INTEGER ithrho, neigmin, nbasiss, nbasis2
INTEGER iunit1, kbasis, lbasis
INTEGER xchanl(nbasis),nangch(narran+1)
INTEGER intwt(narran), jrot(nbasis)
REAL(Kind=WP_Kind) overlap, rholast, rho, overlp1, overlp2, xfac, xnorm
REAL(Kind=WP_Kind) philast(nang, nbasis), s(nbasiss, nbasiss)
REAL(Kind=WP_Kind) phinow(nang, nbasis)
REAL(Kind=WP_Kind) vectlast(nbasiss, nbasiss), vect(nbasiss, nbasiss)
REAL(Kind=WP_Kind) amatrx(nbasiss,nbasiss)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('ovrbas  ', little, medium, full, ithcall, ithsub)

IF(qcase)THEN
   CALL ovrbas2(rholast, rho, nang, nbasis, xchanl, philast, phinow, &
    vectlast, vect, s,  amatrx, ithrho, nbasiss, neigmin, iunit1, jrot)
   RETURN
ENDIF
!
!  --------------------------------------------------------------
!  evaluate overlaps of primitive bases of previous and present rho.
!  --------------------------------------------------------------
xnorm = sqrt(half)
xfac = sqrt(two)

! PORT (2026-10-04): OpenMP over the outer basis loop.  Each iteration writes
! only s(ibasis,jbasis) (its own row) and accumulates into private scalars, so
! the loop is race-free and results are bit-identical to the serial version.
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP& PRIVATE(ibasis, jbasis, iang, iangst1, iangend1, iangst2, iangend2, &
!$OMP&         kbasis, lbasis, jpar, nbasis2, overlap, overlp1, overlp2)
DO ibasis = 1, nbasiss
      iangst1 = nangch(xchanl(ibasis))
      iangend1 = nangch(xchanl(ibasis)+1)-1
!  this overlap matrix should not be symmetric, so full range of
!  jbasis is necessary.
   DO jbasis = 1, nbasiss
      iangst2 = nangch(xchanl(jbasis))
      iangend2 = nangch(xchanl(jbasis)+1)-1
      overlap = zero

     IF(symmetry)THEN
            nbasis2=nbasis-nbasiss
            kbasis = jbasis + nbasis2
            lbasis = ibasis + nbasis2
          IF(jeven)THEN
            jpar = 0
          ELSE
            jpar = 1
          ENDIF
       IF(xchanl(ibasis)==1.AND.xchanl(jbasis)==1)THEN
            DO iang = iangst1, iangend1
            overlap = overlap + philast(iang,ibasis) *phinow(iang,jbasis)
            ENDDO
       ELSEIF(xchanl(ibasis)==1.AND.xchanl(jbasis)==2)THEN
         IF(intwt(xchanl(ibasis))==intwt(xchanl(jbasis)))THEN
           overlp1=zero
           overlp2=zero
           IF(mod(jrot(jbasis),2)==jpar)THEN
             DO iang = iangst1, iangend1
               overlp1 = overlp1 + philast(iang, ibasis) *(phinow(iang,jbasis)+phinow(iang,kbasis))
            ENDDO
             overlp1 = xnorm*overlp1
           ELSE
             DO iang = iangst1, iangend1
               overlp1 = overlp1 + philast(iang, ibasis)*(phinow(iang,jbasis)-phinow(iang,kbasis))
            ENDDO
             overlp1 = xnorm*overlp1
           ENDIF

           DO iang = iangst2, iangend2
             overlp2 = overlp2+philast(iang,ibasis)*phinow(iang,jbasis)
          ENDDO
           overlp2 = xfac*overlp2
           overlap = half*(overlp1+overlp2)

         ELSEIF(intwt(xchanl(ibasis))>intwt(xchanl(jbasis)))THEN
           IF(mod(jrot(jbasis),2)==jpar)THEN
             DO iang = iangst1, iangend1
               overlap = overlap + philast(iang, ibasis)*(phinow(iang,jbasis)+phinow(iang,kbasis))
            ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = iangst1, iangend1
               overlap = overlap + philast(iang, ibasis)*(phinow(iang,jbasis)-phinow(iang,kbasis))
            ENDDO
             overlap = xnorm*overlap
           ENDIF

         ELSE
           DO iang = iangst2, iangend2
             overlap = overlap+philast(iang,ibasis)*phinow(iang,jbasis)
          ENDDO
           overlap = xfac*overlap
         ENDIF
        ELSEIF(xchanl(ibasis)==2.AND.xchanl(jbasis)==1)THEN
         IF(intwt(xchanl(ibasis))==intwt(xchanl(jbasis)))THEN
           overlp1=zero
           overlp2=zero
           IF(mod(jrot(ibasis),2)==jpar)THEN
             DO iang = iangst2, iangend2
               overlp1 = overlp1 + phinow(iang, jbasis) *(philast(iang,ibasis)+philast(iang,lbasis))
             ENDDO
             overlp1 = xnorm*overlp1
           ELSE
             DO iang = iangst2, iangend2
               overlp1 = overlp1 + phinow(iang, jbasis) *(philast(iang,ibasis)-philast(iang,lbasis))
             ENDDO
             overlp1 = xnorm*overlp1
           ENDIF

           DO iang = iangst1, iangend1
             overlp2 = overlp2+philast(iang,ibasis)*phinow(iang,jbasis)
           ENDDO
           overlp2 = xfac*overlp2
           overlap = half*(overlp1+overlp2)

         ELSEIF(intwt(xchanl(jbasis))>intwt(xchanl(ibasis)))THEN
           IF(mod(jrot(ibasis),2)==jpar)THEN
             DO iang = iangst2, iangend2
               overlap = overlap + phinow(iang, jbasis)*(philast(iang,ibasis)+philast(iang,lbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = iangst2, iangend2
               overlap = overlap + phinow(iang, jbasis)*(philast(iang,ibasis)-philast(iang,lbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF

         ELSE
           DO iang = iangst1, iangend1
             overlap = overlap+philast(iang,ibasis)*phinow(iang,jbasis)
           ENDDO
           overlap = xfac*overlap
         ENDIF

       ELSE
         IF(mod(jrot(jbasis),2)==jpar)THEN
           DO iang = iangst1, iangend1
             overlap = overlap + philast(iang, ibasis)*(phinow(iang, jbasis) + phinow(iang, kbasis))
           ENDDO
         ELSE
           DO iang = iangst1, iangend1
             overlap = overlap + philast(iang, ibasis)*(phinow(iang, jbasis) - phinow(iang, kbasis))
           ENDDO
         ENDIF
       ENDIF
     ELSE

      IF(xchanl(ibasis)==xchanl(jbasis))THEN
         DO iang = iangst1, iangend1
            overlap = overlap + philast(iang, ibasis)* phinow(iang,jbasis)
         ENDDO
      ELSEIF(intwt(xchanl(ibasis))==intwt(xchanl(jbasis)))THEN
!  these overlaps should be same for quadrature in either arrangement.
         overlp1=zero
         overlp2=zero
         DO iang = iangst1, iangend1
            overlp1 = overlp1 + philast(iang, ibasis)* phinow(iang,jbasis)
         ENDDO
         DO iang = iangst2, iangend2
            overlp2 = overlp2 + philast(iang, ibasis)* phinow(iang,jbasis)
         ENDDO
         overlap = half*(overlp1+overlp2)
      ELSEIF(intwt(xchanl(ibasis))>intwt(xchanl(jbasis)))THEN
         DO iang = iangst1, iangend1
            overlap = overlap + philast(iang, ibasis)* phinow(iang,jbasis)
         ENDDO
      ELSE
         DO iang = iangst2, iangend2
            overlap = overlap + philast(iang, ibasis)* phinow(iang,jbasis)
         ENDDO
      ENDIF
     ENDIF
      s(ibasis, jbasis) = overlap
   ENDDO
ENDDO
!$OMP END PARALLEL DO

!  ---------------------------------------------------------------
!  s now contains overlaps of primitives of previous and present rho.
!  now READ coeffs of surface fcns at previous rho, vectlast, from
!  iunit1.  This is kase=1 (READ) for type=1 (Ovrlp) integrals.
!  ---------------------------------------------------------------
IF(mega==megamin) REWIND iunit1
CALL rdwrsurb(1, iunit1, mega, megamin, nang, nbasis,neigmin, philast, philast, vectlast, 1, nbasis,nbasiss, nbasiss)

!  --------------------------------------------------------------
!  DO similarity transform to get overlaps of surface functions.
!  snew=vectlast(trans)*s*vect contains overlaps of surface fcns.
!  -----------------------------------------------------------------
CALL dgemm('N','N',nbasiss,neigmin,nbasiss, 1.d0, s,       nbasiss,vect,  nbasiss, 0.d0, amatrx, nbasiss)
CALL dgemm('T','N',neigmin,neigmin,nbasiss, 1.d0, vectlast,nbasiss,amatrx,nbasiss, 0.d0, s,      nbasiss)
!  --------------------------------------------------------------
!  WRITE overlap matrix on unit Ovr_Unit 'Ovrlp'
!  --------------------------------------------------------------
   WRITE(Ovr_Unit) ithrho-1, rholast, rho
   WRITE(Ovr_Unit) neigmin, rholast, rho
   DO jbasis=1,neigmin
      WRITE(Ovr_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
   ENDDO
!  construct s*s(trans) to test completeness of new basis.
   CALL dgemm('N','T',neigmin,neigmin,neigmin, 1.d0, s,nbasiss,s,nbasiss, 0.d0, amatrx,nbasiss)
!-----------------------------------------------------------------------
! IF desired WRITE out the overlap matrix.
!-----------------------------------------------------------------------
   IF(little)THEN
      WRITE(Out_Unit,*)'routine ovrbas'
      WRITE(Out_Unit,*)'nbasiss = ', nbasiss
      WRITE(Out_Unit,'(A,2ES15.7)') 'new basis completeness test diagonal elements for rholast, rho= ',rholast, rho
      WRITE(Out_Unit,'(10ES11.3)')(amatrx(ibasis,ibasis),ibasis=1,neigmin)
   ENDIF
   IF(medium.AND..NOT.full)THEN
      WRITE(Out_Unit,'(A,2ES15.7)') 'Part of the overlap matrix at rholast, rho = ',rholast, rho
      DO ibasis = 1, min(nbasiss,10)
         WRITE(Out_Unit,'(10ES11.3)')(s(ibasis,jbasis),jbasis = 1, min(nbasiss,10))
      ENDDO
      WRITE(Out_Unit,'(A,2ES15.7)') 'new basis completeness test at rholast, rho= ',  rholast, rho
      DO ibasis = 1, min(nbasiss,10)
         WRITE(Out_Unit,'(10ES11.3)')(amatrx(ibasis,jbasis),jbasis = 1, min(nbasiss,10))
      ENDDO
   ENDIF
   IF(full)THEN
      WRITE(Out_Unit,'(A,2ES15.7)')' Overlap matrix (OvrBas) at rholast, rho = ', rholast, rho, " NBasiss=", NBasiss
      CALL MxOut(s, nbasiss, nbasiss)
      WRITE(Out_Unit,'(A,2ES15.7)') 'new basis completeness test OvrBas) rholast, rho= ', rholast, rho, " NBasiss=", NBasiss
      CALL MxOut(amatrx, nbasiss, nbasiss)
!  note.  IF(neigmin )<nbasiss, the last rows and columns
!  of s and amatrx as just written will contain junk.
ENDIF

RETURN
ENDSUBROUTINE OvrBas
