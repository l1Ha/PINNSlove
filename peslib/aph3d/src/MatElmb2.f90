SUBROUTINE MatElmb2(nang, nbasi1, nbasi2, phi1, v, phi2, s, rho, nbasiss1, nbasiss2, chanl1, chanl2, jrot1, jrot2)

!            P U R P O S E   O F    S U B R O U T I N E
!  This routine calculates coupling matrix elements
!  over primitive basis functions.

!            I N P U T    A R G U M E N T S
!   nang     number of quadrature angles.
!   nbasi   number of basis functions.
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
USE Numbers_Module
IMPLICIT NONE
INTEGER nang, ibasis, jbasis, nbasi1, iang
INTEGER nbasi2, nbasiss1, nbasiss2, jpar
INTEGER nbasis21, nbasis22, kbasis, lbasis
INTEGER chanl1(nbasi1), chanl2(nbasi2)
INTEGER jrot1(nbasi1), jrot2(nbasi2)
REAL(Kind=WP_Kind) overlap, rho, xfac, xnorm
REAL(Kind=WP_Kind) phi1(nang,nbasi1), s(mxbasis,mxbasis), v(nang)
REAL(Kind=WP_Kind) phi2(nang, nbasi2)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('matelmb ', little, medium, full, ithcall, ithsub)

!  ------------------------------------------------------------------
!  start loops over basis functions.
!  ------------------------------------------------------------------
xnorm = sqrt(half)
xfac = sqrt(two)

DO ibasis = 1, nbasiss1
   DO jbasis = 1, nbasiss2
      overlap = zero

!  case of symmetry=true follows

     IF(symmetry)THEN
            nbasis21=nbasi1-nbasiss1
            nbasis22=nbasi2-nbasiss2
            lbasis = ibasis + nbasis21
            kbasis = jbasis + nbasis22
          IF(jeven)THEN
            jpar = 0
          ELSE
            jpar = 1
          ENDIF
! 1-1 integrals

       IF(chanl1(ibasis)==1.AND.chanl2(jbasis)==1)THEN
            DO iang = 1, nang
            overlap = overlap + phi1(iang,ibasis)*phi2(iang,jbasis)*v(iang)
            ENDDO

!  1-2 integrals

       ELSEIF(chanl1(ibasis)==1.AND.chanl2(jbasis)==2)THEN
!  proper linear combinations of fcns in channels 2 and 3
           IF(mod(jrot2(jbasis),2)==jpar)THEN
             DO iang = 1, nang
               overlap = overlap + phi1(iang, ibasis)*v(iang) *(phi2(iang, jbasis) + phi2(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = 1, nang
               overlap = overlap + phi1(iang, ibasis)*v(iang)*(phi2(iang, jbasis) - phi2(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF

!  2-1 integrals

       ELSEIF(chanl1(ibasis)==2.AND.chanl2(jbasis)==1)THEN
!  proper linear combinations of fcns in channels 2 and 3
           IF(mod(jrot1(ibasis),2)==jpar)THEN
             DO iang = 1, nang
               overlap = overlap + phi2(iang, jbasis)*v(iang) *(phi1(iang,ibasis) + phi1(iang,lbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = 1, nang
               overlap = overlap + phi2(iang, jbasis)*v(iang)*(phi1(iang,ibasis) - phi1(iang,lbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF


!  2-2 and 2-3 integrals together, for proper symmetry and norm

       ELSE
         IF(mod(jrot1(ibasis),2)==jpar.AND. mod(jrot2(jbasis),2)==jpar)THEN
           DO iang = 1, nang
             overlap = overlap + v(iang)*(phi1(iang,ibasis) + phi1(iang,lbasis))*(phi2(iang, jbasis) + phi2(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
         IF(mod(jrot1(ibasis),2)==jpar.AND. mod(jrot2(jbasis),2)/=jpar)THEN
           DO iang = 1, nang
             overlap = overlap + v(iang)* (phi1(iang,ibasis) + phi1(iang,lbasis))*(phi2(iang, jbasis) - phi2(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF

         IF(mod(jrot1(ibasis),2)/=jpar.AND. mod(jrot2(jbasis),2)==jpar)THEN
           DO iang = 1, nang
             overlap = overlap + v(iang)*(phi1(iang,ibasis) - phi1(iang,lbasis))*(phi2(iang, jbasis) + phi2(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
         IF(mod(jrot1(ibasis),2)/=jpar.AND. mod(jrot2(jbasis),2)/=jpar)THEN
           DO iang = 1, nang
             overlap = overlap + v(iang)*(phi1(iang,ibasis) - phi1(iang,lbasis))*(phi2(iang, jbasis) - phi2(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
!  ENDIF on arrangement channels of basis fcns
       ENDIF
!  case of symmetry=false follows
     ELSE
      DO iang = 1, nang
         overlap = overlap+ phi1(iang,ibasis)*v(iang)*phi2(iang,jbasis)
      ENDDO
!  END of symmetry IF follows
     ENDIF
     s(ibasis, jbasis) = overlap
   ENDDO
ENDDO

!-----------------------------------------------------------------------
! IF desired WRITE out the primitive coupling matrix.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine matelmb2'
   WRITE(Out_Unit,*)'nbasiss1, nbasiss2 = ', nbasiss1, nbasiss2
ENDIF

IF(medium.AND..NOT.full)THEN
   WRITE(Out_Unit,*) 'Part of primitive coupling ', 'matrix at rho = ', rho
   CALL MxOut(s, nbasiss1, nbasiss2)
ENDIF

IF(full)THEN
   WRITE(Out_Unit,*)' Primitive coupling matrix at rho = ', rho
   CALL MxOut(s, nbasiss1, nbasiss2)
ENDIF

RETURN
ENDSUBROUTINE MatElmb2
