SUBROUTINE NEWHam_3d(t_rho,rho_val,nrho,ntheta,nchi,lam,jtot,parity,a,b,d,t_chir1,npro,npro1,npro2,isurface,ecut3,kindex)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE DVR_Module
USE DVRMasses_Module
USE OneD_Module
USE DVR2_Module
USE Parms_Module
USE TwoD_Module
USE Masses_Module
!
!     this SUBROUTINE takes the 2-d hamiltonians from NEWham_2d and the
!     coupling matrices from angmom and turns them into a full 3-d
!     hamiltonian matrix.
!     the resulting matrix is stored in h3d which
!     is THEN diagonalized to give eigenvectors (in h3d) and
!     eigenvalues (in e3d).
!
IMPLICIT NONE
INTEGER nrho, ntheta, nchi, ntot, ntvmaxx, ngri, jrho, i2, k, i, j, lam, jtot
INTEGER itheta, ichi, indexk, indexi, nchir, ngrj, npro, indexj, npro1
INTEGER ntempi, item4, i4, i1, item1, npro2, parity
REAL(Kind=WP_Kind) t_rho(nrho,nrho),a(ntheta),b(ntheta),d(ntheta)
REAL(Kind=WP_Kind) t_chir1(nchi,nchi),rho_val(nrho),ecut3,etemp1,etemp2
REAL(Kind=WP_Kind),ALLOCATABLE::h3d(:,:),sd(:),subd(:)
REAL(Kind=WP_Kind) waa, htemp, cfactor, cfac, temp, afactor, afac
INTEGER irho, isurface, ii, jj, kindex, itemp1, itemp2
REAL(Kind=WP_Kind) vtemp1(nrho),vtemp2(nrho)
LOGICAL Gtemp
EXTERNAL mprint
CHARACTER(LEN=1)lamlabel1
CHARACTER(LEN=4)Glabel
CHARACTER(LEN=22)label1,label2
CHARACTER(LEN=19)Glabel2


ntot=nrho

ntvmaxx=(ntot)*(ntot+1)/2
!
!     h3d--3D hamiltonian

IF(Ntheta>NthetMax.or.Nchi>NchiMax)THEN
   WRITE(*,*)"Ntheta or Nchi exceeds maximum allowed values"
   WRITE(*,*)"Ntheta=",Ntheta,"  NthetMax=",NthetMax
   WRITE(*,*)"Nchi=",Nchi,"  NchiMax=",NchiMax
   STOP "NewHam3D"
ENDIF
ngri=0
ALLOCATE(h3d(ntot,ntot))
h3d=0.0d0

!OPEN(Unit=sfunfbr_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//'BinOut/h3d_temp',form='unformatted')

!
!     cc2dtot is in the truncated surface basis: D=dC
!





!
!     begin outer loop over rho
!
DO irho=1,nrho
!
!     ngrj and ngri will help us keep track of where we are in
!     the full hamiltonian matrix. basically, they are partial
!     sums of the number of 2-d eigenvectors

!
!     begin inner loop over rho
!
   DO jrho=irho,nrho
      waa=t_rho(irho,jrho)

!
!     now transform the (diagonal) block of the coupling matrix
!

       htemp=0.0d0
!
!     because waa is diagonal , we only need to DO one loop
!     over angles. this is the heart of the fast transformation
!
       DO i2=1,ntheta*nchi
          htemp=htemp+cc2dtot(i2,isurface,irho)*cc2dtot(i2,isurface,jrho)*waa
       ENDDO

!
!     when done, place terms in the correct place in h3d
!
            h3d(irho,jrho)=htemp
            h3d(jrho,irho)=htemp

   ENDDO
!
!     another time saving maneuver. don't add in full 2-d hamiltonian
!     just add in 2-d eigenvalues along the diagonal
!

      h3d(irho,irho)=h3d(irho,irho)+e2d_new(isurface,irho)
ENDDO
!
!     done constructing h3d
!


ngri=nrho
k=0

ALLOCATE(sd(ngri*(ngri+1)/2),subd(ngri))
DO 9 i=1,ngri
   DO 29 j=1,i
      k=k+1
      sd(k)=h3d(i,j)
 29      CONTINUE
 9    CONTINUE
h3d=0.0d0

CALL tdiagrw(sd,e3d_pro(1,isurface),cc3dtot(1,1,isurface),subd,ngri,ngri,ngri*(ngri+1)/2)

DEALLOCATE(sd,subd)
DEALLOCATE(h3d)

k=kindex
IF(e3d_pro(1,isurface)<ecut3)THEN

DO irho=1,nrho
k=k+1
IF(e3d_pro(irho,isurface)<ecut3)THEN
IF(k<10)THEN
   WRITE(GLabel,'(4I1)')0,0,0,k
ELSEIF(k<100)THEN
   WRITE(GLabel,'(2I1,I2)')0,0,k
ELSEIF(k<1000)THEN
   WRITE(GLabel,'(I1,I3)')0,k
ELSEIF(k<10000)THEN
   WRITE(GLabel,'(I4)')k
ELSE
   WRITE(Out_Unit,*)'state number is greater than 9999:'
   STOP 'Diagonalization'
ENDIF
GLabel2='BinOut/Gfunction'//Glabel
IF(.true.)THEN
OPEN(Unit=GLabel_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//Glabel2,form='unformatted')
REWIND(GLabel_Unit)
!--------------------------------------------------------------------------------------------------------
!Start to WRITE the G basis on the disk
!--------------------------------------------------------------------------------------------------------
ii=0
jj=1
DO j=1,nrho*nchi*ntheta
ii=ii+1
IF(ii>nchi*ntheta)THEN
ii=1
jj=jj+1
ENDIF
WRITE(GLabel_Unit)cc2dtot(ii,isurface,jj)*cc3dtot(jj,irho,isurface)
ENDDO
CLOSE(GLabel_Unit)
ENDIF
Gindex(irho,isurface)=k
!--------------------------------------------------------------------------------------------------------
!writing is finished
!--------------------------------------------------------------------------------------------------------
ELSE
kindex=k-1
GOTO 123
ENDIF

ENDDO
kindex=k
123 CONTINUE

ELSE
kindex=k
ENDIF

!--------------------------------------------------------------------------------------------------------
!checking basis oliver need to be deleted
!--------------------------------------------------------------------------------------------------------
IF(.false.)THEN
DO i=1,nrho
   DO j=1,nrho
      htemp=0.0d0
      DO k=1,nrho
      htemp=htemp+cc3dtot(k,i,isurface)*cc3dtot(k,j,isurface)
      ENDDO
      IF(i==j)THEN
           IF(ABS(htemp-1.0d0)>1.0d-6)THEN
         WRITE(Out_Unit, * )'wrong 3dbasis=',i,j,htemp
         STOP 'error 3d'
         ENDIF
      ELSE
           IF(ABS(htemp)>1.0d-6)THEN
           WRITE(Out_Unit, * )'wrong 3dbasis=',i,j,htemp
         STOP 'error 3d'
         ENDIF
      ENDIF
   ENDDO
ENDDO
ENDIF






ii=1
IF(isurface/=1)THEN

!--------------------------------------------------------------------------------------------------------
!Start the energy sorting process
!--------------------------------------------------------------------------------------------------------
DO i=1,irho-1
DO j=ii,nstate
   IF(e3d_pro(i,isurface)<Genergy(j))THEN
   DO jj=j+1,nstate
   IF(jj==(j+1))THEN
   etemp1=Genergy(jj-1)
   itemp1=Gindex2(jj-1)
   DO jrho=1,nrho
      vtemp1(jrho)=v_pro(jrho,jj-1)
   ENDDO
   ENDIF
   etemp2=Genergy(jj)
   itemp2=Gindex2(jj)
   DO jrho=1,nrho
   vtemp2(jrho)=v_pro(jrho,jj)
   ENDDO

   Genergy(jj)=etemp1
   Gindex2(jj)=itemp1
   DO jrho=1,nrho
   v_pro(jrho,jj)=vtemp1(jrho)
   ENDDO

   etemp1=etemp2
   itemp1=itemp2
   vtemp1=vtemp2
   ENDDO

   GENERGY(j)=e3d_pro(i,isurface)
   Gindex2(j)=Gindex(i,isurface)
   DO jrho=1,nrho
   v_pro(jrho,j)=e2d_new(isurface,jrho)
   ENDDO

   Gtemp=.true.
   GOTO 124
   ENDIF
ENDDO
124 CONTINUE
IF(.NOT.Gtemp)GOTO 125
Gtemp=.false.
ii=j
ENDDO

125 CONTINUE
!--------------------------------------------------------------------------------------------------------
!finish the energy sorting process
!--------------------------------------------------------------------------------------------------------
ELSE
DO i=1,irho-1
IF(i<=nstate)THEN
Gindex2(i)=Gindex(i,1)
Genergy(i)=e3d_pro(i,1)
   DO j=1,nrho
   v_pro(j,i)=e2d_new(isurface,j)
   ENDDO
ENDIF
ENDDO
ENDIF


!temporary RETURN
RETURN





!--------------------------------------------------------------------------------------------------------
!New strategy construct the lower half of the Hamiltonian row by row
!--------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------
!1.The diagonal part
!--------------------------------------------------------------------------------------------------------
!DO i=1,ngri
!DO j=1,i
!indexi=i+npro
!indexj=j+npro
!WRITE(GFunction_Unit)indexi,indexj,h3d(i,j)
!Hamil(indexi,indexj)=h3d(i,j)
!Hamil(indexj,indexi)=Hamil(indexi,indexj)
!ENDDO
!ENDDO



!-------------------------------------------------------------------
! 2 construct Coriolis coupling terms----off diagonal for the same row
!-------------------------------------------------------------------
IF(lam>=1)THEN

!BETTER WAY TO TRANSFORM [lam,lam-1] block
vectortemp=0.0d0
cfactor = cfac (jtot, lam, lam - 1, parity)

WRITE(lamlabel1,'(I1)')lam-1
Label1='BinOut/transformation'//lamlabel1
OPEN(Unit=Label1_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//Label1,form='unformatted')
REWIND(Label1_Unit)
DO irho=1,nrho
!WRITE(Msg_Unit,*)'here irho=',irho,lam-1
!read in cc2dtot1 for irho
DO i=1,ntheta*nchi
DO j=1,i
READ(Label1_Unit)cc2dtemp(i,j)
cc2dtemp(j,i)=cc2dtemp(i,j)
ENDDO
ENDDO

  DO i=1,ntheta*nchi
    DO j=1,nev_3d1(irho)
       itheta=(i-1)/nchi+1
       ichi=i-(itheta-1)*nchi
       DO k=1,nchi
       temp= - 1.0d0/rho_val (irho)**2 * d (itheta) * t_chir1 (ichi, k) * cfactor
       indexk=(itheta-1)*nchi+k
!      vectortemp(i,j)=vectortemp(i,j)+temp*cc2dtot1(indexk,j,irho)
       vectortemp(i,j)=vectortemp(i,j)+temp*cc2dtemp(indexk,j)
       ENDDO
    ENDDO
 ENDDO


  DO i=1,nev_3d(irho)
    DO j=1,nev_3d1(irho)
        temp=0.0d0
      DO k=1,ntheta*nchi
        temp=temp+cc2dtot(k,i,irho)*vectortemp(k,j)
      ENDDO
      indexi=npro+i
      indexj=npro1+j
      WRITE(GFunction_Unit)indexi,indexj,temp
    ENDDO
  ENDDO

ENDDO
!CLOSE(CC2dtemp_Unit)
CLOSE(Label1_Unit)
ENDIF
!-------------------------------------------------------------------
!3 Construct Asymmetric Top Coupling terms.
!-------------------------------------------------------------------
IF(lam>=2)THEN
vectortemp=0.0d0
afactor = afac (jtot, lam, lam - 2, parity)
   DO irho = 1, nrho
      DO itheta = 1, ntheta
         DO ichi = 1, nchir
            i =(irho-1)*ntheta*nchi+(itheta-1)*nchi+ichi
            vectortemp(i, i) = (a (itheta) - b (itheta) ) / 2.d0 /rho_val (irho)**2  * afactor / usys2
         ENDDO
      ENDDO
   ENDDO
ngrj=0


!DO the transformation faster


WRITE(lamlabel1,'(I1)')lam-2
Label2='BinOut/transformation'//lamlabel1
OPEN(Unit=Label2_Unit,File=OutDIR(1:LEN(TRIM(OutDIR)))//Label2,form='unformatted')

DO  irho=1,nrho

!read in cc2dtot2 for irho
DO i=1,ntheta*nchi
DO j=1,i
READ(Label2_Unit)cc2dtemp(i,j)
cc2dtemp(j,i)=cc2dtemp(i,j)
ENDDO
ENDDO

   ntempi=nchi*ntheta
      DO  i4=1,ntempi
         item4=ngrj+i4
         DO i1=1,ntempi
            item1=ngrj+i1
            htemp=0.0d0
            DO  i2=1,ntempi
          indexi=i2+(irho-1)*ntempi
               htemp=htemp+cc2dtot(i2,i4,irho)*cc2dtemp(i2,i1)*vectortemp(indexi,indexi)
            ENDDO
         indexi=npro+item4
         indexj=npro2+item1
         IF(item4<=ngri.AND.item1<=nev_3d2(irho))THEN
         WRITE(GFunction_Unit)indexi,indexj,htemp
         ENDIF
!            h3d(item4,item1)=htemp
!            h3d(item1,item4)=htemp
         ENDDO
     ENDDO
   ngrj=ngrj+ntempi
ENDDO

!CLOSE(CC2dtemp_Unit)
CLOSE(Label2_Unit)
ENDIF

IF(jtot/=0)THEN
   cc2dtot1=cc2dtot
   cc2dtot2=cc2dtot1
   nev_3d1=nev_3d
   nev_3d2=nev_3d1
ENDIF



!DEALLOCATE(h3d)

RETURN

END
