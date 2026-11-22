SUBROUTINE RdWrSurb(kase, iunit, mega, megamin, nang, nbasis, neigmin, phipr, phider, vect, ktype, nbasisdi, nbasiss, nbasssdi)

!                 P U R P O S E   O F   S U B R O U T I N E

!  reads or WRITEs primitive basis, its derivative w.o. chi, and the
!  surface function coefficients of the primitive basis.
!  Used by matbas and ovrbas.

!                 A R G U M E N T S

!  kase      kase=1 reads. kase=2 WRITEs.
!  iunit     unit number sf33_unit or sf34_unit) written to or READ from.
!  phipr     primitive basis at the quadrature angles.
!  phider    derivative of primitive basis w.o. chi at quad. angles.
!  vect      vect(nbasis,neigmin).  coeffs of surface functions relative
!            to the primitive basis.
!  nbasis    number of raw primitive basis functions.
!  nbasiss   number of symmetry-adapted primitive basis functions.
!  nbasisdi  dimension of matrices using nbasis.
!  nbasssdi  dimension of matrices using nbasiss.
USE Numeric_Kinds_Module
USE FileUnits_Module
IMPLICIT NONE

!                   I N T E G E R S

INTEGER ibasis, jbasis, iang, kase, iunit, mega, megamin
INTEGER nang, nbasis, ktype, megar, nbasisr, nangr
INTEGER neigmin, neigminr, nbasisdi, nbasiss, nbasssdi

!                   R E A L S

REAL(Kind=WP_Kind) phipr(nang, nbasisdi), phider(nang, nbasisdi)
REAL(Kind=WP_Kind) vect(nbasssdi, nbasssdi)

IF(kase==1)THEN
!  --------------------------------------------------------------
!  READ surface functions.
!  position head correctly first.
!  --------------------------------------------------------------
1  READ(iunit) megar, nbasisr, neigminr
   IF(megar>mega)THEN
      WRITE(Out_Unit,*)'megar, mega=', megar, mega
      STOP 'rdwrsurb1'
   ENDIF
   IF(megar<mega)THEN
      DO jbasis = 1, neigminr
         READ(iunit)
      ENDDO
      GOTO 1
   ELSE
      ! PORT (2026-10-12): Fix historic bug that blocked large-basis convergence.
      ! Across sectors, the number of good eigenvalues neigen can vary by a few
      ! modes due to threshold filtering. The stored record has neigminr vectors,
      ! while the current sector only requests the leading neigmin vectors.
      ! Requiring neigminr <= neigmin was backwards: if the file has MORE vectors
      ! than requested (neigminr > neigmin), we simply read the leading neigmin
      ! and skip the rest to maintain exact tape alignment. An error only occurs
      ! if the file has FEWER vectors than requested (neigminr < neigmin).
      IF(nbasisr/=nbasiss .or. neigminr<neigmin)THEN
         WRITE(Out_Unit,*)'nbasisr,nbasiss,neigminr,neigmin=', nbasisr,nbasiss,neigminr,neigmin
         WRITE(Msg_Unit,*)'nbasisr,nbasiss=', nbasisr, nbasiss
         WRITE(Msg_Unit,*)'neigminr,neigmin=', neigminr,neigmin
         STOP 'RdWrSurb2'
      ENDIF
!  -----------------------------------------------------------------
!  READ vector of surface fcn coeffs for mega value called.
!  ------------------------------------------------------------------
      DO jbasis = 1, neigmin
         READ(iunit)(vect(ibasis,jbasis), ibasis=1,nbasiss)
      ENDDO
      IF(neigminr > neigmin)THEN
         DO jbasis = neigmin + 1, neigminr
            READ(iunit)
         ENDDO
      ENDIF
   ENDIF

   IF(ktype==2)THEN
!  ------------------------------------------------------------------
!  potential matrix els.  READ phipr for mega and advance phider tape.
!  ------------------------------------------------------------------
2     READ(phip_unit) megar, nangr, nbasisr
      READ(phid_unit)
      IF(nangr/=nang.or.megar>mega)THEN
         WRITE(Out_Unit,*)'nangr, nang, megar, mega=', nangr, nang, megar, mega
         STOP 'rdwrsurb3'
      ENDIF
      IF(megar<mega)THEN
         DO jbasis=1, nbasisr
            READ(phip_unit)
            READ(phid_unit)
         ENDDO
         GOTO 2
      ELSE
         IF(nbasisr/=nbasis)THEN
            WRITE(Out_Unit,*)'nbasisr, nbasis=',nbasisr, nbasis
            STOP 'rdwrsurb4'
         ENDIF
         DO jbasis=1, nbasis
            READ(phip_unit)(phipr(iang,jbasis), iang=1, nang)
            READ(phid_unit)
         ENDDO
      ENDIF

   ELSEIF(ktype==3)THEN
!  ------------------------------------------------------------------
!  coriolis matrix elements.  READ phider for megap=mega+1 and
!  advance phipr tape.
!  ------------------------------------------------------------------
      READ(phip_unit)
      READ(phid_unit) megar, nangr, nbasisr
         IF(megar/=mega.or.nangr/=nang.or.nbasisr/=nbasis)THEN
            WRITE(Out_Unit,*)'megar,mega,nangr,', 'nang,nbasisr,nbasis=',  megar,mega,nangr,nang,nbasisr,nbasis
            STOP 'rdwrsurb5'
         ENDIF
      DO jbasis=1, nbasis
         READ(phip_unit)
         READ(phid_unit)(phider(iang,jbasis), iang=1, nang)
      ENDDO

   ELSEIF(ktype==4)THEN
!  ------------------------------------------------------------------
!  asym top matrix elements.  READ phipr for megap=mega+2 and put
!  it into phider, which is phi2 of matbas.
!  ------------------------------------------------------------------
      READ(phip_unit)megar, nangr, nbasisr
         IF(megar/=mega.or.nangr/=nang.or.nbasisr/=nbasis)THEN
            WRITE(Out_Unit,*)'megar,mega,nangr,', 'nang,nbasisr,nbasis=', megar,mega,nangr,nang,nbasisr,nbasis
            STOP 'rdwrsurb6'
         ENDIF
      DO jbasis=1, nbasis
         READ(phip_unit)(phider(iang,jbasis), iang=1, nang)
      ENDDO

   ELSEIF(ktype<1.or.ktype>4)THEN
      WRITE(Out_Unit,*)'ktype=', ktype
      STOP 'rdwrsurb7'
   ENDIF



ELSE
!  --------------------------------------------------------------
!  WRITE surface functions.
!  --------------------------------------------------------------
   IF(mega==megamin)THEN
      REWIND iunit
      REWIND phip_unit
      REWIND phid_unit
   ENDIF
   WRITE(phip_unit) mega, nang, nbasis
   DO jbasis=1, nbasis
      WRITE(phip_unit) (phipr(iang,jbasis), iang=1, nang)
   ENDDO
   WRITE(phid_unit) mega, nang, nbasis
   DO jbasis=1, nbasis
      WRITE(phid_unit) (phider(iang,jbasis), iang=1, nang)
   ENDDO
   WRITE(iunit) mega, nbasiss, neigmin
   DO jbasis=1, neigmin
      WRITE(iunit) (vect(ibasis,jbasis), ibasis=1, nbasiss)
   ENDDO
ENDIF
RETURN
ENDSUBROUTINE RdWrSurb
