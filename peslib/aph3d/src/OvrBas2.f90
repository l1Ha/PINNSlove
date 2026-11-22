SUBROUTINE OvrBas2(rholast, rho, nang, nbasis, xchanl, philast, phinow, vectlast, vect, s, amatrx, ithrho, nbasiss, &
                   neigmin, iunit1, jrot)
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
!   philast  primitive basis of previous rho at present quadrature pts.
!   phinow   primitive basis of present rho at present quadrature pts.
!   vectlast coeffs of surface fcns of previous rho.
!   vect     coeffs of present surface functions.

!            O U T P U T    A R G U M E N T S
!   s        overlap matrix elements.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE TotJ_Module
USE FileUnits_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER nang, ibasis, jbasis, nbasis, iang, jpar
INTEGER ithrho, neigmin, nbasiss, nbasis2
INTEGER iunit1, kbasis, lbasis
INTEGER xchanl(nbasis), jrot(nbasis)
REAL(Kind=WP_Kind) overlap, rholast, rho, xfac, xnorm
REAL(Kind=WP_Kind) philast(nang, nbasis), s(nbasiss, nbasiss)
REAL(Kind=WP_Kind) phinow(nang, nbasis)
REAL(Kind=WP_Kind) vectlast(nbasiss, nbasiss), vect(nbasiss, nbasiss)
REAL(Kind=WP_Kind) amatrx(nbasiss,nbasiss)
EXTERNAL rdwrsurb
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('ovrbas  ', little, medium, full, ithcall, ithsub)

!
!  --------------------------------------------------------------
!  evaluate overlaps of primitive bases of previous and present rho.
!  --------------------------------------------------------------
xnorm = sqrt(half)
xfac = sqrt(two)

DO ibasis = 1, nbasiss
   DO jbasis = 1, nbasiss
      overlap = zero
!  case of symmetry=true follows

     IF(symmetry)THEN
            nbasis2=nbasis-nbasiss
            lbasis = ibasis + nbasis2
            kbasis = jbasis + nbasis2
          IF(jeven)THEN
            jpar = 0
          ELSE
            jpar = 1
          ENDIF
! 1-1 integrals
       IF(xchanl(ibasis)==1.AND.xchanl(jbasis)==1)THEN
            DO iang = 1, nang
            overlap = overlap + philast(iang,ibasis)*phinow(iang,jbasis)
            ENDDO
!  1-2 integrals
       ELSEIF(xchanl(ibasis)==1.AND.xchanl(jbasis)==2)THEN
!  proper linear combinations of fcns in channels 2 and 3
           IF(mod(jrot(jbasis),2)==jpar)THEN
             DO iang = 1, nang
               overlap = overlap + philast(iang, ibasis) *(phinow(iang, jbasis) + phinow(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = 1, nang
               overlap = overlap + philast(iang, ibasis)*(phinow(iang, jbasis) - phinow(iang, kbasis))
            ENDDO
             overlap = xnorm*overlap
           ENDIF
!  2-1 integrals

       ELSEIF(xchanl(ibasis)==2.AND.xchanl(jbasis)==1)THEN
!  proper linear combinations of fcns in channels 2 and 3
           IF(mod(jrot(ibasis),2)==jpar)THEN
             DO iang = 1, nang
               overlap = overlap + phinow(iang, jbasis)*(philast(iang,ibasis) + philast(iang,lbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = 1, nang
               overlap = overlap + phinow(iang, jbasis)*(philast(iang,ibasis) - philast(iang,lbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF
!  2-2 and 2-3 integrals together, for proper symmetry and norm
       ELSE
         IF(mod(jrot(ibasis),2)==jpar.AND.  mod(jrot(jbasis),2)==jpar)THEN
           DO iang = 1, nang
             overlap = overlap +(philast(iang,ibasis) + philast(iang,lbasis))*(phinow(iang, jbasis) + phinow(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
         IF(mod(jrot(ibasis),2)==jpar.AND. mod(jrot(jbasis),2)/=jpar)THEN
           DO iang = 1, nang
             overlap = overlap +(philast(iang,ibasis) + philast(iang,lbasis))*(phinow(iang, jbasis) - phinow(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF

         IF(mod(jrot(ibasis),2)/=jpar.AND. mod(jrot(jbasis),2)==jpar)THEN
           DO iang = 1, nang
             overlap = overlap +(philast(iang,ibasis) - philast(iang,lbasis))*(phinow(iang, jbasis) + phinow(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
         IF(mod(jrot(ibasis),2)/=jpar.AND. mod(jrot(jbasis),2)/=jpar)THEN
           DO iang = 1, nang
             overlap = overlap +(philast(iang,ibasis) - philast(iang,lbasis))*(phinow(iang, jbasis) - phinow(iang, kbasis))
            ENDDO
           overlap = half*overlap
         ENDIF
!  ENDIF on arrangement channels of basis fcns
       ENDIF
!  case of symmetry=false follows
     ELSE
      DO  iang = 1, nang
         overlap=overlap+philast(iang,ibasis)*phinow(iang,jbasis)
      ENDDO
!  END of symmetry IF follows
     ENDIF
      s(ibasis, jbasis) = overlap
   ENDDO
ENDDO

!  ---------------------------------------------------------------
!  s now contains overlaps of primitives of previous and present rho.
!  now READ coeffs of surface fcns at previous rho, vectlast, from
!  iunit1.  This is kase=1 (READ) for type=1 (Ovrlp) integrals.
!  ---------------------------------------------------------------
IF(mega==megamin) REWIND iunit1
CALL rdwrsurb(1, iunit1, mega, megamin, nang, nbasis, neigmin, philast, philast, vectlast, 1, nbasis, nbasiss, nbasiss)

!  --------------------------------------------------------------
!  DO similarity transform to get overlaps of surface functions.
!  snew=vectlast(trans)*s*vect contains overlaps of surface fcns.
!  -----------------------------------------------------------------
CALL dgemm('N','N',nbasiss,neigmin,nbasiss, 1.d0 ,s,       nbasiss,vect,  nbasiss, 0.d0, amatrx,nbasiss)
CALL dgemm('T','N',neigmin,neigmin,nbasiss, 1.d0, vectlast,nbasiss,amatrx,nbasiss, 0.d0, s,     nbasiss)

!  --------------------------------------------------------------
!  WRITE overlap matrix on unit Ovr_Unit 'Ovrlp'
!  --------------------------------------------------------------
   WRITE(Ovr_Unit) ithrho-1, rholast, rho
   WRITE(Ovr_Unit) neigmin, rholast, rho
   DO jbasis=1,neigmin
      WRITE(Ovr_Unit)(s(ibasis,jbasis),ibasis=1,neigmin)
   ENDDO
!-----------------------------------------------------------------------
! IF desired WRITE out the overlap matrix.
!-----------------------------------------------------------------------
   IF(little)THEN
      WRITE(Out_Unit,*)'routine ovrbas2'
      WRITE(Out_Unit,*)'nbasiss = ', nbasiss
!  construct s*s(trans) to test completeness of new basis.
      CALL dgemm('N','T',neigmin,neigmin,neigmin, 1.d0, s,nbasiss,s,nbasiss,0.d0, amatrx,nbasiss)
      WRITE(Out_Unit,*) 'new basis completeness test.  ','rholast, rho= ',rholast, rho
      WRITE(Out_Unit,*)(amatrx(ibasis,ibasis),ibasis=1,neigmin)
   ENDIF
   IF(medium.AND..NOT.full)THEN
      WRITE(Out_Unit,*) 'Part of the overlap matrix at ', 'rholast, rho = ', rholast, rho
      CALL MxOut(s,nbasiss,nbasiss)
      WRITE(Out_Unit,*) 'new basis completeness test.  ','rholast, rho= ', rholast, rho
      CALL MxOut(amatrx,nbasiss,nbasiss)
   ENDIF
   IF(full)THEN
      WRITE(Out_Unit,*)' Overlap matrix at ', 'rholast, rho = ', rholast,rho
      CALL MxOut(s, nbasiss, nbasiss)
      WRITE(Out_Unit,*) 'new basis completeness test.  ', 'rholast, rho= ', rholast, rho
      CALL MxOut(amatrx, nbasiss, nbasiss)
!  note.  IF(neigmin )<nbasiss, the last rows and columns
!  of s and amatrx as just written will contain junk.
ENDIF

RETURN
ENDSUBROUTINE OvrBas2
