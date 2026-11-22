SUBROUTINE PotElmnt (nang, nbasis, xchanl, nangch, phi, v, s, rho, intwt, nbasisdi, nbasiss,   nbasssdi, jrot)
!

!            P U R P O S E   O F    S U B R O U T I N E
! This routine calculates potential energy matrix elements.
!  For more explanation, see ovlaps.F

!            I N P U T    A R G U M E N T S
!   nang     number of quadrature angles.
!   nbasis   number of basis functions.
!   nbasiss  no. of symmetry adapted basis functions.
!   nbasisdi dims of matrices. diff for calls by sfunbas and matbas
!   narran   number of arrangement channels.
!   xchanl    an array giving the channel number of each basis
!            function.
!   nangch   an array containing the starting position of the
!            angles for each of the arrangement channels.
!            nangch(narran+1) is equal to nang+1.
!   phi      primitive basis functions evaluated at each of the
!            angles.
!   v        an array of the potential energy evaluated at each
!            of the angles.
!   rho      current hyperradius.

!            O U T P U T    A R G U M E N T S
!   s        potential matrix elements.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE Parms_Module
USE TotJ_Module
USE FileUnits_Module
USE QCase_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER nang, ibasis, jbasis, nbasis, iang, iangst1, iangend1
INTEGER iangst2, iangend2, kbasis
INTEGER jpar, nbasisdi, nbasiss, nbasis2, nbasssdi
INTEGER jrot(nbasisdi), nangch(narran+1), intwt(narran), xchanl(nbasisdi)
REAL(Kind=WP_Kind) overlap, rho, overlp1, overlp2, xnorm, xfac
REAL(Kind=WP_Kind) phi(nang,nbasisdi), s(nbasssdi,nbasssdi),  v(nang)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('potelmnt', little, medium, full, ithcall, ithsub)

IF(qcase)THEN
   CALL potelmnt2(nang, nbasis, xchanl, phi, v, s, rho, nbasisdi, nbasiss,   nbasssdi, jrot)
   RETURN
ENDIF
!
xnorm=sqrt(half)
xfac = sqrt(two)
! PORT (2026-10-03): OpenMP over the outer basis loop. Iterations are
! independent: iteration i writes row i (columns >= i) and column i
! (rows >= i) only, so no element is written twice and the arithmetic
! order per element is unchanged (results are bit-identical to serial).
! The loop over ibasis is the dominant cost at large nbasiss
! (~nbasiss^2 x nang) and is the practical bottleneck of the ABM stage.
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP& PRIVATE(ibasis, jbasis, iang, iangst1, iangend1, iangst2, iangend2, &
!$OMP&         kbasis, nbasis2, jpar, overlap, overlp1, overlp2)
DO ibasis = 1, nbasiss
      iangst1 = nangch(xchanl(ibasis))
      iangend1 = nangch(xchanl(ibasis)+1)-1
   DO jbasis = ibasis, nbasiss
      iangst2 = nangch(xchanl(jbasis))
      iangend2 = nangch(xchanl(jbasis)+1)-1
      overlap = zero

     IF(symmetry)THEN
            nbasis2 = nbasis-nbasiss
            kbasis = jbasis + nbasis2
            IF(jeven)THEN
               jpar = 0
            ELSE
               jpar = 1
            ENDIF
       IF(xchanl(ibasis)==1.AND.xchanl(jbasis)==1)THEN
            DO iang = iangst1, iangend1
            overlap = overlap + phi(iang, ibasis)*phi(iang,jbasis)*v(iang)
            ENDDO
       ELSEIF(xchanl(ibasis)==1.AND.xchanl(jbasis)==2)THEN
         IF(intwt(xchanl(ibasis))==intwt(xchanl(jbasis)))THEN
           overlp1=zero
           overlp2=zero
           IF(mod(jrot(jbasis),2)==jpar)THEN
             DO iang = iangst1, iangend1
               overlp1 = overlp1 + phi(iang, ibasis)*v(iang)*(phi(iang, jbasis) + phi(iang, kbasis))
             ENDDO
             overlp1 = xnorm*overlp1
           ELSE
             DO iang = iangst1, iangend1
               overlp1 = overlp1 + phi(iang, ibasis)*v(iang)*(phi(iang, jbasis) - phi(iang, kbasis))
             ENDDO
             overlp1 = xnorm*overlp1
           ENDIF

           DO iang = iangst2, iangend2
             overlp2 = overlp2+phi(iang,ibasis)*phi(iang,jbasis)*v(iang)
           ENDDO
           overlp2 = xfac*overlp2
           overlap = half*(overlp1+overlp2)

         ELSEIF(intwt(xchanl(ibasis))>intwt(xchanl(jbasis)))THEN
           IF(mod(jrot(jbasis),2)==jpar)THEN
             DO iang = iangst1, iangend1
               overlap = overlap + phi(iang, ibasis)*v(iang)*(phi(iang, jbasis) + phi(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = iangst1, iangend1
               overlap = overlap + phi(iang, ibasis)*v(iang)*(phi(iang, jbasis) - phi(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF

         ELSE
           DO iang = iangst2, iangend2
             overlap = overlap+phi(iang,ibasis)*phi(iang,jbasis)*v(iang)
           ENDDO
           overlap = xfac*overlap
         ENDIF
       ELSE
         IF(mod(jrot(jbasis),2)==jpar)THEN
           DO iang = iangst1, iangend1
             overlap = overlap + phi(iang, ibasis)*v(iang)*(phi(iang, jbasis) + phi(iang, kbasis))
           ENDDO
         ELSE
           DO iang = iangst1, iangend1
             overlap = overlap + phi(iang, ibasis)*v(iang)*(phi(iang, jbasis) - phi(iang, kbasis))
           ENDDO

         ENDIF

       ENDIF
     ELSE
      IF(xchanl(ibasis)==xchanl(jbasis))THEN
         DO iang = iangst1, iangend1
            overlap = overlap + phi(iang, ibasis)*phi(iang,jbasis)* v(iang)
         ENDDO
      ELSEIF(intwt(xchanl(ibasis))==intwt(xchanl(jbasis)))THEN
         overlp1=zero
         overlp2=zero
         DO iang = iangst1, iangend1
            overlp1 = overlp1 + phi(iang, ibasis)*phi(iang,jbasis) * v(iang)
         ENDDO
         DO iang = iangst2, iangend2
            overlp2 = overlp2 + phi(iang, ibasis)*phi(iang,jbasis) * v(iang)
         ENDDO
         overlap = half*(overlp1+overlp2)
      ELSEIF(intwt(xchanl(ibasis))>intwt(xchanl(jbasis)))THEN
         DO iang = iangst1, iangend1
            overlap = overlap + phi(iang, ibasis)*phi(iang,jbasis) * v(iang)
         ENDDO
      ELSE
         DO iang = iangst2, iangend2
            overlap = overlap + phi(iang, ibasis)*phi(iang,jbasis) * v(iang)
         ENDDO
      ENDIF
     ENDIF
      IF(ABS(overlap).le.100*Epsilon(One))overlap=Zero      !TMPmod GregParker
      s(ibasis, jbasis) = overlap
      s(jbasis, ibasis) = s(ibasis, jbasis)
   ENDDO
ENDDO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
! IF desired WRITE out the potential matrix.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine potelmnt'
   WRITE(Out_Unit,*)'nbasiss = ', nbasiss
ENDIF
IF(medium.AND..NOT.full)THEN
   WRITE(Out_Unit,'(A,ES20.12)') 'Part of the Potential matrix (PotElmnt) at rho = ', rho
   DO ibasis = 1, min(nbasiss,10)
      WRITE(Out_Unit,'(10ES15.7)')(s(ibasis,jbasis),jbasis = 1, min(nbasiss,10))
   ENDDO
ENDIF
IF(full)THEN
   WRITE(Out_Unit,'(A,ES20.12,A,I5)')' Potential matrix (PotElmnt) at rho = ', rho, " NBasisdi=", NBasisdi
   CALL MxOut(s, nbasiss, nbasiss)   ! CALL MxOut(s, nbasisdi, nbasisdi) GregParker TEMPMOD
ENDIF

RETURN
ENDSUBROUTINE PotElmnt
