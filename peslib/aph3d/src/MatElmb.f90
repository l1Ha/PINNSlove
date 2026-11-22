SUBROUTINE MatElmb (nang, nbasi1, nbasi2, nangch, phi1, v, phi2, s, rho, intwt, chanl1, chanl2,   &
                    nbasiss1, nbasiss2, jrot1, jrot2)

!            P U R P O S E   O F    S U B R O U T I N E
!  This routine calculates coupling matrix elements
!  over primitive basis functions.

!            I N P U T    A R G U M E N T S
!   nang     number of quadrature angles.
!   nbasi   number of basis functions.
!   narran   number of arrangement channels.
!   chanl    an array giving the channel number of each basis
!            function.
!   nangch   an array containing the starting position of the
!            angles for each of the arrangement channels.
!            nangch(narran+1) is equal to nang+1.
!   phi      primitive basis functions evaluated at each of the
!            angles.
!   v        an array of the coupling potential evaluated at each
!            of the angles.
!   rho      current hyperradius.
!   jeven    true  gives jrot even in the incident arrangement  for
!            symmetric (A+B2) systems.  This
!            is the conserved quantity, not reflection in chi.  This
!            routine does Coriolis matrix elements which
!            couple functions
!            with different parities under reflection in chi.  rtp 92/8.

!            O U T P U T    A R G U M E N T S
!   s        coupling matrix elements.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE Parms_Module
USE TotJ_Module
USE FileUnits_Module
USE QCase_Module
IMPLICIT NONE

!            I N T E G E R S
INTEGER nang, ibasis, jbasis, nbasi1, iang, iangst1, iangend1
INTEGER iangst2, iangend2
INTEGER nbasi2, jpar, nbasiss1, nbasiss2
INTEGER nbasis21, nbasis22, kbasis, lbasis
INTEGER chanl1(nbasi1), chanl2(nbasi2), jrot1(nbasi1)
INTEGER nangch(narran+1), intwt(narran), jrot2(nbasi2)

!            R E A L S
REAL(Kind=WP_Kind) overlap, rho, overlp1, overlp2, xfac, xnorm

!            D I M E N S I O N S
REAL(Kind=WP_Kind) phi1(nang,nbasi1), s(mxbasis,mxbasis), v(nang)
REAL(Kind=WP_Kind) phi2(nang, nbasi2)

!              P A R A M E T E R S
REAL(Kind=WP_Kind),  PARAMETER:: zero=0.0d0, half=0.5d0, two=2.0d0
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('matelmb ', little, medium, full, ithcall, ithsub)

IF(qcase)THEN
   CALL matelmb2(nang, nbasi1, nbasi2, phi1, v, phi2, s, rho, nbasiss1, nbasiss2, chanl1, chanl2, jrot1, jrot2)
   RETURN
ENDIF

!  ------------------------------------------------------------------
!  start loops over basis functions.
!  ------------------------------------------------------------------
xnorm = sqrt(half)
xfac = sqrt(two)

DO ibasis = 1, nbasiss1
      iangst1 = nangch(chanl1(ibasis))
      iangend1 = nangch(chanl1(ibasis)+1)-1
!  these matrix elements should not be symmetric, so full range of
!  jbasis is necessary.
   DO jbasis = 1, nbasiss2
      iangst2 = nangch(chanl2(jbasis))
      iangend2 = nangch(chanl2(jbasis)+1)-1
      overlap = zero

     IF(symmetry)THEN
            nbasis21=nbasi1-nbasiss1
            nbasis22=nbasi2-nbasiss2
            kbasis = jbasis + nbasis22
            lbasis = ibasis + nbasis21
          IF(jeven)THEN
            jpar = 0
          ELSE
            jpar = 1
          ENDIF
       IF(chanl1(ibasis)==1.AND.chanl2(jbasis)==1)THEN
            DO iang = iangst1, iangend1
            overlap = overlap + phi1(iang,ibasis)*phi2(iang,jbasis)*v(iang)
            ENDDO
       ELSEIF(chanl1(ibasis)==1.AND.chanl2(jbasis)==2)THEN
         IF(intwt(chanl1(ibasis))==intwt(chanl2(jbasis)))THEN
           overlp1=zero
           overlp2=zero
           IF(mod(jrot2(jbasis),2)==jpar)THEN
             DO iang = iangst1, iangend1
               overlp1 = overlp1 + phi1(iang,ibasis)*v(iang) *(phi2(iang,jbasis)+phi2(iang,kbasis))
             ENDDO
             overlp1 = xnorm*overlp1
           ELSE
             DO iang = iangst1, iangend1
               overlp1 = overlp1 + phi1(iang,ibasis)*v(iang)*(phi2(iang,jbasis)-phi2(iang,kbasis))
             ENDDO
             overlp1 = xnorm*overlp1
           ENDIF

           DO iang = iangst2, iangend2
             overlp2 = overlp2+phi1(iang,ibasis)*phi2(iang,jbasis)*v(iang)
           ENDDO
           overlp2 = xfac*overlp2
           overlap = half*(overlp1+overlp2)

         ELSEIF(intwt(chanl1(ibasis))>intwt(chanl2(jbasis)))THEN
           IF(mod(jrot2(jbasis),2)==jpar)THEN
             DO iang = iangst1, iangend1
               overlap = overlap + phi1(iang,ibasis)*v(iang)*(phi2(iang,jbasis)+phi2(iang,kbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = iangst1, iangend1
               overlap = overlap + phi1(iang,ibasis)*v(iang)*(phi2(iang,jbasis)-phi2(iang,kbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF

         ELSE
           DO iang = iangst2, iangend2
             overlap = overlap+phi1(iang,ibasis)*phi2(iang,jbasis)*v(iang)
           ENDDO
           overlap = xfac*overlap
         ENDIF
        ELSEIF(chanl1(ibasis)==2.AND.chanl2(jbasis)==1)THEN
         IF(intwt(chanl1(ibasis))==intwt(chanl2(jbasis)))THEN
           overlp1=zero
           overlp2=zero
           IF(mod(jrot1(ibasis),2)==jpar)THEN
             DO iang = iangst2, iangend2
               overlp1 = overlp1 + phi2(iang,jbasis)*v(iang)*(phi1(iang,ibasis)+phi1(iang,lbasis))
             ENDDO
             overlp1 = xnorm*overlp1
           ELSE
             DO iang = iangst2, iangend2
               overlp1 = overlp1 + phi2(iang,jbasis)*v(iang)*(phi1(iang,ibasis)-phi1(iang,lbasis))
             ENDDO
             overlp1 = xnorm*overlp1
           ENDIF

           DO iang = iangst1, iangend1
             overlp2 = overlp2+phi1(iang,ibasis)*phi2(iang,jbasis)*v(iang)
           ENDDO
           overlp2 = xfac*overlp2
           overlap = half*(overlp1+overlp2)

         ELSEIF(intwt(chanl2(jbasis))>intwt(chanl1(ibasis)))THEN
           IF(mod(jrot1(ibasis),2)==jpar)THEN
             DO iang = iangst2, iangend2
               overlap = overlap + phi2(iang,jbasis)*v(iang)*(phi1(iang,ibasis)+phi1(iang,lbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = iangst2, iangend2
               overlap = overlap + phi2(iang,jbasis)*v(iang)*(phi1(iang,ibasis)-phi1(iang,lbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF

         ELSE
           DO iang = iangst1, iangend1
             overlap = overlap+phi1(iang,ibasis)*phi2(iang,jbasis)*v(iang)
           ENDDO
           overlap = xfac*overlap
         ENDIF

       ELSE
         IF(mod(jrot2(jbasis),2)==jpar)THEN
           DO iang = iangst1, iangend1
             overlap = overlap + phi1(iang,ibasis)*v(iang)*(phi2(iang, jbasis) + phi2(iang, kbasis))
           ENDDO
         ELSE
           DO iang = iangst1, iangend1
             overlap = overlap + phi1(iang,ibasis)*v(iang)*(phi2(iang, jbasis) - phi2(iang, kbasis))
           ENDDO

         ENDIF

       ENDIF

     ELSE

      IF(chanl1(ibasis)==chanl2(jbasis))THEN
         DO iang = iangst1, iangend1
            overlap = overlap+phi1(iang,ibasis)*phi2(iang,jbasis)* v(iang)
         ENDDO
      ELSEIF(intwt(chanl1(ibasis))==intwt(chanl2(jbasis)))THEN
         overlp1=zero
         overlp2=zero
         DO iang = iangst1, iangend1
            overlp1 = overlp1+phi1(iang,ibasis)*phi2(iang,jbasis) * v(iang)
         ENDDO
         DO iang = iangst2, iangend2
            overlp2 = overlp2+phi1(iang,ibasis)*phi2(iang,jbasis) * v(iang)
         ENDDO
         overlap = half*(overlp1+overlp2)
      ELSEIF(intwt(chanl1(ibasis))>intwt(chanl2(jbasis)))THEN
         DO iang = iangst1, iangend1
            overlap = overlap+phi1(iang,ibasis)*phi2(iang,jbasis)* v(iang)
         ENDDO
      ELSE
         DO iang = iangst2, iangend2
            overlap = overlap+phi1(iang,ibasis)*phi2(iang,jbasis) * v(iang)
         ENDDO
      ENDIF
     ENDIF

      s(ibasis, jbasis) = overlap
   ENDDO
ENDDO

!-----------------------------------------------------------------------
! IF desired WRITE out the primitive coupling matrix.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine matelmb'
   WRITE(Out_Unit,*)'nbasiss1, nbasiss2 = ', nbasiss1, nbasiss2
ENDIF
IF(medium.AND..NOT.full)THEN
   WRITE(Out_Unit,*) 'Part of primitive coupling ', 'matrix at rho = ', rho
   CALL MxOut(s, nbasiss1, nbasiss2)
ENDIF
IF(full)THEN
   WRITE(Out_Unit,*)' Primitive coupling matrix at rho = ', rho
   CALL MxOut(s, nbasiss1, nbasiss2)
!  note.  this matrix may not be very readable.  it is dimensioned
!  mxbasis*mxbasis.
ENDIF

RETURN
ENDSUBROUTINE MatElmb
