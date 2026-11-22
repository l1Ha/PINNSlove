SUBROUTINE PotElmnt2(nang, nbasis, xchanl, phi, v, s, rho, nbasisdi, nbasiss,  nbasssdi, jrot)
!

!            P U R P O S E   O F    S U B R O U T I N E
! This routine calculates potential energy matrix elements.

!            I N P U T    A R G U M E N T S
!   nang     number of quadrature angles.
!   nbasis   number of basis functions.
!   nbasiss  no. of symmetry adapted basis functions.
!   nbasisdi dims of matrices. diff for calls by sfunbas and matbas
!   narran   number of arrangement channels.
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
USE Numbers_Module
IMPLICIT NONE
INTEGER nang, ibasis, jbasis, nbasis, iang, nbasis2
INTEGER nbasisdi, nbasiss, nbasssdi, kbasis, lbasis, jpar
INTEGER xchanl(nbasisdi), jrot(nbasisdi)
REAL(Kind=WP_Kind) overlap, rho, xnorm, xfac
REAL(Kind=WP_Kind) phi(nang,2*nbasisdi), s(nbasssdi,nbasssdi), v(nang)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('potelmnt', little, medium, full, ithcall, ithsub)

!
xnorm=sqrt(half)
xfac = sqrt(two)
! PORT (2026-10-03): OpenMP over the outer basis loop (same argument as
! in PotElmnt: each iteration writes only the L-shape of index ibasis,
! so iterations are independent and results stay bit-identical).
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) &
!$OMP& PRIVATE(ibasis, jbasis, iang, kbasis, lbasis, nbasis2, jpar, overlap)
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
               overlap = overlap + v(iang)*phi(iang,ibasis)*phi(iang,jbasis)
            ENDDO

!  1-2 integrals

       ELSEIF(xchanl(ibasis)==1.AND.xchanl(jbasis)==2)THEN
!  proper linear combinations of fcns in channels 2 and 3
           IF(mod(jrot(jbasis),2)==jpar)THEN
             DO iang = 1, nang
               overlap = overlap + phi(iang, ibasis)*v(iang)*(phi(iang, jbasis) + phi(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ELSE
             DO iang = 1, nang
               overlap = overlap + phi(iang, ibasis)*v(iang)*(phi(iang, jbasis) - phi(iang, kbasis))
             ENDDO
             overlap = xnorm*overlap
           ENDIF

!  2-2 and 2-3 integrals together, for proper symmetry and norm

       ELSE
         IF(mod(jrot(ibasis),2)==jpar.AND. mod(jrot(jbasis),2)==jpar)THEN
           DO iang = 1, nang
             overlap = overlap + v(iang)*(phi(iang, ibasis) + phi(iang, lbasis))*(phi(iang, jbasis) + phi(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
         IF(mod(jrot(ibasis),2)==jpar.AND.  mod(jrot(jbasis),2)/=jpar)THEN
           DO iang = 1, nang
             overlap = overlap + v(iang)*(phi(iang, ibasis) + phi(iang, lbasis))*(phi(iang, jbasis) - phi(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF

         IF(mod(jrot(ibasis),2)/=jpar.AND. mod(jrot(jbasis),2)==jpar)THEN
           DO iang = 1, nang
             overlap = overlap + v(iang)*(phi(iang, ibasis) - phi(iang, lbasis))*(phi(iang, jbasis) + phi(iang, kbasis))
           ENDDO
           overlap = half*overlap
         ENDIF
         IF(mod(jrot(ibasis),2)/=jpar.AND.  mod(jrot(jbasis),2)/=jpar)THEN
           DO iang = 1, nang
             overlap = overlap + v(iang)*(phi(iang, ibasis) - phi(iang, lbasis))*(phi(iang, jbasis) - phi(iang, kbasis))
          ENDDO
           overlap = half*overlap
         ENDIF

!  ENDIF on arrangement channels of basis fcns

       ENDIF

!  case of symmetry=false follows

     ELSE


      DO iang = 1, nang
         overlap = overlap + phi(iang,ibasis)*v(iang)*phi(iang,jbasis)
      ENDDO

!  END of symmetry IF follows

     ENDIF

      s(ibasis, jbasis) = overlap
      s(jbasis, ibasis) = s(ibasis, jbasis)
   ENDDO
ENDDO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
! IF desired WRITE out the potential matrix.
!-----------------------------------------------------------------------
IF(little)THEN
   WRITE(Out_Unit,*)'routine potelmnt2'
   WRITE(Out_Unit,*)'nbasiss = ', nbasiss
ENDIF
IF(medium.AND..NOT.full)THEN
   WRITE(Out_Unit,'(A,ES20.12)') 'Part of the Potential matrix (PotElmnt2) at rho = ', rho
   DO 6 ibasis = 1, min(nbasiss,10)
      WRITE(Out_Unit,*)(s(ibasis,jbasis),jbasis = 1, min(nbasiss,10))
    6    ENDDO
ENDIF
IF(full)THEN
   WRITE(Out_Unit,'(A,ES20.12,A,I5)')' Potential matrix (PotElmnt2) at rho = ', rho, " NBasisdi=", NBasisdi
   CALL MxOut(s, nbasisdi, nbasisdi)
ENDIF

RETURN
ENDSUBROUTINE PotElmnt2
