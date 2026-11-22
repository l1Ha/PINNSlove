SUBROUTINE addban (a, maxa, s, lm, nid, nidm, idi, beta, neq)
USE Numeric_Kinds_Module
USE FileUnits_Module
!
!   $RCSfile: addban.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:02:48 $
!   $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E
!----------------------------------------------------------------------
!     .     a = global stiffness                                      .
!     .     s = element stiffness                                     .
!     .                                                               .
!     .               s(1)        s(2)        s(3)        . . .       .
!     .     s   =                 s(18+1)     s(18+2)     . . .       .
!     .                                       s(2*18)     . . .       .
!     .                                                   . . .       .
!     .                                                               .
!     .                                                               .
!     .               a(1)        a(3)        a(6)        . . .       .
!     .     a   =                 a(2)        a(5)        . . .       .
!     .                                       a(4)        . . .       .
!----------------------------------------------------------------------


!     I N P U T    A R G U M E N T S
!  a
!  maxa
!  s
!  lm
!  nid
!  nidm
!  idi
!  beta
!  neq


IMPLICIT NONE
INTEGER i, ii, ij, j, jj, k, kk, ks, kss, l, mce, mi, nce, ndi, neq, nidm
INTEGER maxa(1), lm(1), nid(1), idi(nidm,1)
REAL(Kind=WP_Kind)   fac, fact
REAL(Kind=WP_Kind) a(1), s(1), beta(nidm,1)

 10   ndi=0
DO i=1,18
   ii=lm(i)
   !IF(ii) 5,160,90
   IF(ii==0)goto 160
   IF(ii>0)goto 90
 5       CONTINUE
   nce=-ii
   DO k=1,nid(nce)
      ii=idi(k,nce)
      IF(ii<1.or.ii>neq)THEN
         WRITE(Out_unit,*)'error ii'
         GOTO 80
      ENDIF
      fac=beta(k,nce)
      mi=maxa(ii)
      ks=i
      DO j=1,18
         jj=lm(j)
         !IF(jj) 20,70,50
         IF(jj==0)goto 70
         IF(jj>0)goto 50
 20            mce=-jj
         kss=ks
         IF(j>=i) kss=j+ndi
         DO 40 l=1,nid(mce)
            jj=idi(l,mce)
            ij=ii-jj
            !IF(ij) 40,30,30
            IF(ij<0)goto 40
 30               kk=mi+ij
            fact=beta(l,mce)
            a(kk)=a(kk)+fac*s(kss)*fact
 40            CONTINUE
         GOTO 70
 50            ij=ii-jj
         !IF(ij) 70,60,60
         IF(ij<0)goto 70
 60            kk=mi+ij
         kss=ks
         IF(j>=i) kss=j+ndi
         a(kk)=a(kk)+s(kss)*fac
 70            ks=ks+18-j
    ENDDO
 80         CONTINUE
 ENDDO
      GOTO 160
 90         IF(ii<1.or.ii>neq)THEN
         WRITE(Out_unit,*)' error in addban: ii<1.or.ii>neq ',ii, neq
         STOP 'addban'
      ENDIF
      mi=maxa(ii)
      ks=i
      DO j=1,18
         jj=lm(j)
         !IF(jj) 100,150,130
         IF(jj==0)goto 150
         IF(jj>0)goto 130
 100           nce=-jj
         kss=ks
         IF(j>=i) kss=j+ndi
         DO 120 k=1,nid(nce)
            jj=idi(k,nce)
            fac=beta(k,nce)
            ij=ii-jj
            !IF(ij) 120,110,110
            IF(ij<0)goto 120
 110              kk=mi+ij
            a(kk)=a(kk)+s(kss)*fac
 120           CONTINUE
         GOTO 150
 130           ij=ii-jj
         !IF(ij) 150,140,140
         IF(ij<0)goto 150
 140           kk=mi+ij
         kss=ks
         IF(j>=i) kss=j+ndi
         a(kk)=a(kk)+s(kss)
 150           ks=ks+18-j
     ENDDO
 160           ndi=ndi+18-i
     ENDDO
         RETURN
         END
