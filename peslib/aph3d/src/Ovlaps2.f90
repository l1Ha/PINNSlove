SUBROUTINE Ovlaps2(nang, nbasis, xchanl, phi, s, rho, nbasiss, jrot)
!

!            P U R P O S E   O F    S U B R O U T I N E
! This routine determines the overlap matrix elements since basis
! set is a nonorthogonal basis.
!  accessed IF qcase=true.  Quadrature in APH coordinates.

!            I N P U T    A R G U M E N T S
!   nang     number of quadrature angles.
!   nbasis   number of basis functions.
!   nbasiss  no. of symmetry adapted basis functions.
!   phi      primitive basis functions evaluated at each of the
!            angles.

!            O U T P U T    A R G U M E N T S
!   s        overlap matrix elements.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE Parms_Module
USE TotJ_Module
USE FileUnits_Module
USE Converge_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER nang, ibasis, jbasis, nbasis, iang, nbasiss
INTEGER kbasis, lbasis, jpar, nbasis2
INTEGER xchanl(nbasis), jrot(nbasis)
REAL(Kind=WP_Kind) overlap, rho, xnorm, xfac, serr, terr, abs
REAL(Kind=WP_Kind) phi(nang,nbasis), s(nbasiss,nbasiss)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('ovlaps  ', little, medium, full, ithcall, ithsub)
!

xnorm = sqrt(half)
xfac = sqrt(two)
DO ibasis = 1, nbasiss
   DO jbasis = ibasis, nbasiss
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
            overlap = overlap + phi(iang,ibasis)*phi(iang,jbasis)
            ENDDO

!  1-2 integrals

       ELSEIF(xchanl(ibasis)==1.AND.xchanl(jbasis)==2)THEN
!  proper linear combinations of fcns in channels 2 and 3
           IF(mod(jrot(jbasis),2)==jpar)THEN
             DO iang = 1, nang
               overlap = overlap + phi(iang, ibasis) *(phi(iang, jbasis) + phi(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = 1, nang
               overlap = overlap + phi(iang, ibasis)*(phi(iang, jbasis) - phi(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF

!  2-2 and 2-3 integrals together, for proper symmetry and norm

       ELSE
         IF(mod(jrot(ibasis),2)==jpar.AND. mod(jrot(jbasis),2)==jpar)THEN
           DO iang = 1, nang
             overlap = overlap + (phi(iang, ibasis) + phi(iang, lbasis))*(phi(iang, jbasis) + phi(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
         IF(mod(jrot(ibasis),2)==jpar.AND. mod(jrot(jbasis),2)/=jpar)THEN
           DO iang = 1, nang
             overlap = overlap + (phi(iang, ibasis) + phi(iang, lbasis))*(phi(iang, jbasis) - phi(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF

         IF(mod(jrot(ibasis),2)/=jpar.AND. mod(jrot(jbasis),2)==jpar)THEN
           DO iang = 1, nang
             overlap = overlap +(phi(iang, ibasis) - phi(iang, lbasis))*(phi(iang, jbasis) + phi(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
         IF(mod(jrot(ibasis),2)/=jpar.AND. mod(jrot(jbasis),2)/=jpar)THEN
           DO iang = 1, nang
             overlap = overlap +(phi(iang, ibasis) - phi(iang, lbasis))*(phi(iang, jbasis) - phi(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF

!  ENDIF on arrangement channels of basis fcns

       ENDIF

!  case of symmetry=false follows

     ELSE

      DO iang = 1, nang
         overlap = overlap + phi(iang,ibasis)*phi(iang,jbasis)
      ENDDO

!  END of symmetry IF follows

     ENDIF

      s(ibasis, jbasis) = overlap
      s(jbasis, ibasis) = s(ibasis, jbasis)
   ENDDO
ENDDO


!-----------------------------------------------------------------------
! IF desired WRITE out the overlap matrix.
!-----------------------------------------------------------------------
   IF(little)THEN
      WRITE(Out_Unit,*)'routine ovlaps2'
      WRITE(Out_Unit,*)'nbasiss = ', nbasiss
      WRITE(Out_Unit,10)(s(ibasis,ibasis),ibasis=1,nbasiss)
  10        FORMAT(1x,'diag_ovlap:',6f10.5)
   ENDIF
   serr = 0.d0
   DO ibasis = 1, nbasiss
      terr = ABS(s(ibasis,ibasis)-1.d0)
      IF(terr>serr)serr=terr
   ENDDO
   WRITE(Out_Unit,*)'Largest Percent Error in ','Overlap=',serr*100.d0
   IF(serr>ovrerr)THEN
      WRITE(Out_Unit,*)'Warning serr>ovrerr'
      WRITE(Out_Unit,*)'ovrerr=',ovrerr,'   serr=',serr
   ENDIF
   IF(medium.AND..NOT.full)THEN
      WRITE(Out_Unit,*) 'Part of the overlap matrix at rho = ', rho
      DO ibasis = 1, min(nbasiss,10)
         WRITE(Out_Unit,*)(s(ibasis,jbasis),jbasis = 1, min(nbasiss,10))
      ENDDO
   ENDIF
   
   IF(full)THEN
      WRITE(Out_Unit,*)'Overlap matrix at rho (Ovlaps2)= ', rho," NBasiss=", NBasiss
      CALL MxOut(s, nbasiss, nbasiss)
   ENDIF

RETURN
ENDSUBROUTINE Ovlaps2
