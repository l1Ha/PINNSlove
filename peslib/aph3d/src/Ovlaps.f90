SUBROUTINE Ovlaps (nang, nbasis, xchanl, nangch, phi, s, rho, intwt, nbasiss, jrot)
!
!            P U R P O S E   O F    S U B R O U T I N E
! This routine determines the overlap matrix elements since basis
! set is a nonorthogonal basis.
!  quadratures in channel delves coordinates.

!            I N P U T    A R G U M E N T S
!   nang     number of quadrature angles.
!   nbasis   number of basis functions.
!   nbasiss  no. of symmetry adapted basis functions.
!   narran   number of arrangement channels.
!   xchanl    an array giving the channel number of each basis
!            function.
!   nangch   an array containing the starting position of the
!            angles for each of the arrangement channels.
!            nangch(narran+1) is equal to nang+1.
!   phi      primitive basis functions evaluated at each of the
!            angles.
!   rho1     previous hyperradius.
!   rho2     current hyperradius.

!            O U T P U T    A R G U M E N T S
!   s        overlap matrix elements.

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
USE Numeric_Kinds_Module
USE Parms_Module
USE TotJ_Module
USE FileUnits_Module
USE QCase_Module
USE Converge_Module
USE Numbers_Module
IMPLICIT NONE
INTEGER nang, ibasis, jbasis, nbasis, iang, iangst1, iangend1
INTEGER iangst2, iangend2, kbasis
INTEGER jpar, novrerr, nbasiss, nbasis2
INTEGER xchanl(nbasis),nangch(narran+1),intwt(narran),jrot(nbasis)
REAL(Kind=WP_Kind) overlap, rho, overlp1, overlp2, err, xnorm, xfac
REAL(Kind=WP_Kind) phi(nang,nbasis), s(nbasiss,nbasiss)
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub
DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./
CALL popt ('ovlaps  ', little, medium, full, ithcall, ithsub)

IF(qcase)THEN
   CALL ovlaps2(nang, nbasis, xchanl, phi, s, rho, nbasiss, jrot)
   RETURN
ENDIF
!
IF(full)THEN
   WRITE(Out_Unit,*)'chanl,iang,phi'
   DO ibasis = 1, 15, 7
      iangst1 = nangch(xchanl(ibasis))
      iangend1 = nangch(xchanl(ibasis)+1)-1
      DO iang = iangst1, iangend1
         !IF(ABS(phi(iang,1))>0.05)THEN
            WRITE(Out_Unit,'(2I5,3es21.13)')ibasis, iang, phi(iang,1), phi(iang,8), phi(iang,15)
         !ENDIF
      ENDDO
   ENDDO
ENDIF

novrerr=0
xnorm = sqrt(half)
xfac = sqrt(two)
DO ibasis = 1, nbasiss
   iangst1 = nangch(xchanl(ibasis))
   iangend1 = nangch(xchanl(ibasis)+1)-1
   DO jbasis = ibasis, nbasiss
      iangst2 = nangch(xchanl(jbasis))
      iangend2 = nangch(xchanl(jbasis)+1)-1
      overlap = zero
      IF(symmetry)THEN  !  case of symmetry=true follows
         nbasis2=nbasis-nbasiss
         kbasis = jbasis + nbasis2
         IF(jeven)THEN
            jpar = 0
         ELSE
            jpar = 1
         ENDIF
         IF(xchanl(ibasis)==1.AND.xchanl(jbasis)==1)THEN        ! 1-1 integrals
            DO iang = iangst1, iangend1
               overlap = overlap + phi(iang,ibasis)*phi(iang,jbasis)
            ENDDO
         ELSEIF(xchanl(ibasis)==1.AND.xchanl(jbasis)==2)THEN     !  1-2 integrals
            !  DO in both channels, first in channel 1
            IF(intwt(xchanl(ibasis))==intwt(xchanl(jbasis)))THEN
               overlp1=zero
               overlp2=zero
               !  proper linear combinations of fcns in channels 2 and 3
               IF(mod(jrot(jbasis),2)==jpar)THEN
                  DO iang = iangst1, iangend1
                     overlp1=overlp1+phi(iang,ibasis)*(phi(iang,jbasis)+phi(iang,kbasis))
                  ENDDO
                  overlp1 = xnorm*overlp1
               ELSE
                  DO iang = iangst1, iangend1
                     overlp1=overlp1+phi(iang,ibasis)*(phi(iang,jbasis)-phi(iang,kbasis))
                  ENDDO
                  overlp1 = xnorm*overlp1
               ENDIF
               !  now DO in channel 2 and average
               DO iang = iangst2, iangend2
                  overlp2 = overlp2+phi(iang,ibasis)*phi(iang,jbasis)
               ENDDO
               overlp2 = xfac*overlp2
               !IF(mod(jrot(jbasis),2)==jpar)THEN
               !   DO iang = iangst2, iangend2
               !      overlp2=overlp2+phi(iang,ibasis)*(phi(iang,jbasis)+phi(iang,kbasis))
               !   ENDDO
               !   overlp2 = xnorm*overlp2
               !ELSE
               !   DO iang = iangst2, iangend2
               !      overlp2=overlp2+phi(iang,ibasis)*(phi(iang,jbasis)-phi(iang,kbasis))
               !   ENDDO
               !   overlp2 = xnorm*overlp2
               !ENDIF
               overlap = half*(overlp1+overlp2)
               err=ABS(overlp2-overlp1)
               IF(ABS(err)>ovrerr)THEN
                  WRITE(Out_Unit,FMT=66)ibasis,jbasis,overlp1,overlp2,err
                  WRITE(Msg_Unit,FMT=66)ibasis,jbasis,overlp1,overlp2,err,overlp1/overlp2
                  novrerr=novrerr+1
               ENDIF
               !  only DO in channel with greater intwt
            ELSEIF(intwt(xchanl(ibasis))>intwt(xchanl(jbasis)))THEN
               !  DO only in channel 1 with proper linear
               !  combination of fcns in 2 and 3
               IF(mod(jrot(jbasis),2)==jpar)THEN
                  DO iang = iangst1, iangend1
                     overlap=overlap+phi(iang,ibasis)*(phi(iang,jbasis)+phi(iang,kbasis))
                  ENDDO
                  overlap = xnorm*overlap
               ELSE
                  DO iang = iangst1, iangend1
                     overlap=overlap+phi(iang,ibasis)*(phi(iang,jbasis)-phi(iang,kbasis))
                  ENDDO
                  overlap = xnorm*overlap
               ENDIF
              !  DO only in channel 2
           ELSE
              DO iang = iangst2, iangend2
                 overlap = overlap+phi(iang,ibasis)*phi(iang,jbasis)
              ENDDO
              overlap = xfac*overlap
           ENDIF
           !  2-2 and 2-3 integrals together, for proper symmetry and norm
         ELSE
            IF(mod(jrot(jbasis),2)==jpar)THEN
               DO iang = iangst1, iangend1
                  overlap=overlap+phi(iang,ibasis)*(phi(iang,jbasis)+phi(iang,kbasis))
               ENDDO
            ELSE
               DO iang = iangst1, iangend1
                  overlap=overlap+phi(iang,ibasis)*(phi(iang,jbasis)-phi(iang,kbasis))
               ENDDO
            ENDIF
            !  ENDIF on arrangement channels of basis fcns
         ENDIF
      ELSE  !  case of symmetry=false follows
         IF(xchanl(ibasis)==xchanl(jbasis))THEN
            DO iang = iangst1, iangend1
               overlap = overlap + phi(iang,ibasis)*phi(iang,jbasis)
            ENDDO
         ELSEIF(intwt(xchanl(ibasis))==intwt(xchanl(jbasis)))THEN
            overlp1=zero
            overlp2=zero
            DO iang = iangst1, iangend1
               overlp1 = overlp1 + phi(iang,ibasis)*phi(iang,jbasis)
            ENDDO
            DO iang = iangst2, iangend2
               overlp2 = overlp2 + phi(iang,ibasis)*phi(iang,jbasis)
            ENDDO
            overlap = half*(overlp1+overlp2)
            err=(overlp2-overlp1)
            IF(little.AND.ABS(err)>ovrerr)THEN
               WRITE(Out_Unit,FMT=66)ibasis,jbasis,overlp1,overlp2,err
               novrerr=novrerr+1
            ENDIF
         ELSEIF(intwt(xchanl(ibasis))>intwt(xchanl(jbasis)))THEN
            DO iang = iangst1, iangend1
               overlap = overlap + phi(iang,ibasis)*phi(iang,jbasis)
            ENDDO
         ELSE
            DO iang = iangst2, iangend2
               overlap = overlap + phi(iang,ibasis)*phi(iang,jbasis)
            ENDDO
         ENDIF
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
   WRITE(Out_Unit,*)'routine ovlaps'
   WRITE(Out_Unit,*)'nbasiss = ', nbasiss
   WRITE(Out_Unit,*)'novrerr = ', novrerr
   WRITE(Msg_Unit,*)'novrerr = ', novrerr
ENDIF

IF(medium.AND..NOT.full)THEN
   WRITE(Out_Unit,*) 'Part of the overlap matrix at rho = ',rho
   DO ibasis = 1, min(nbasiss,10)
         WRITE(Out_Unit,*)(s(ibasis,jbasis),jbasis = 1, min(nbasiss,10))
   ENDDO
ENDIF

IF(full)THEN
   WRITE(Out_Unit,*)'Overlap matrix at rho (Ovlaps)= ', rho," NBasiss=", NBasiss
   CALL MxOut(s, nbasiss, nbasiss)
ENDIF


66 FORMAT(1x,' ovrerr ',2i5,4e15.7)
RETURN
ENDSUBROUTINE Ovlaps
