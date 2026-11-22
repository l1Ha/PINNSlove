SUBROUTINE diagbasb (rho, nbasis, nhermt, noscil, oscil, nglegn, nlegndre, pn, parityfn, &
                     phi, tee, lam, vect,  bfaci, nxdim, nydim, nbasiss)

!            P U R P O S E   O F    S U B R O U T I N E
!   On the first CALL, this routine sets up the primitive basis set at
!   the quadrature points.  On all calls it transforms to produce the
!   tee'th surface function phi at the quadrature points.

!            I N P U T    A R G U M E N T S
!   nbasis   number of basis functions.
!   narran   number of arrangement channels.
!   nhermt   an array containing the number of Gauss_Hermite
!            points used in each arrangement channel.
!   noscil   an array containing the max order of oscillators used in
!            each arrangement channel.
!   nmax     maximum number of oscillators used in any arrangement
!            channel.
!   oscil    an array containing the oscillators evaluated at each
!             of the angles.
!   nglegn   an array containing the number of Gauss_Legendre
!            points used in each arrangement channel.
!   nlegndre an array containing the order of the Legendre polynomial
!            used in each arrangement channel.
!   lmax     maximum order of the Legendre polynomial in any channel.
!   pn       an array containing the Legendre polynomials evaluated
!            at each of the angles.
!   parityfn an array containing the cosine of the APH chi angle
!            for odd parity functions and 1 for even parity
!            functions.
!   bfaci    a factor that goes with the Jacobian and sho functions.

!            O U T P U T    A R G U M E N T S
!   phipr    an array of the primitive basis evaluated at each of
!            the quadrature points.
!   phi      the aph surface functions at the quadrature points.
!
!this routine is called by:
!         sfunbasb
!this routine calls
!         popt
!-----------------------------------------------------------------------
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!               I N C L U D E   F I L E S
USE Numeric_Kinds_Module
USE fileunits_Module
USE Narran_Module
USE quantb_Module
USE totj_Module
USE Numbers_Module
SAVE
INTEGER ibasis, jbasis, nbasisb, nbasis, karran, marran, khermt, n
INTEGER kglegn, megas, l, tee, lam, neigmin, nxdim, nydim
INTEGER nbasiss, nbasis2, nbasis1, jpar, nthcall, lambda
INTEGER nhermt(narran), nglegn(narran), noscil(narran), nlegndre(narran)
REAL(Kind=WP_Kind), ALLOCATABLE:: phipr(:,:,:,:)
REAL(Kind=WP_Kind) xnorm, rho, rhob, delrho
REAL(Kind=WP_Kind) phi(nxdim,nydim,narran), pn(0:mxl,mxglegn,narran)
REAL(Kind=WP_Kind) oscil(0:mxn, maxhermt, narran)
REAL(Kind=WP_Kind) parityfn(mxglegn,maxhermt,narran), vect(nbasiss,nbasiss)
REAL(Kind=WP_Kind) bfaci(maxhermt,narran)

DATA nthcall /0/, megas /-1/
!-----------------------------------------------------------------------
! Determine printing options.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub

DATA ithcall /0/, ithsub /0/
DATA little /.false./, medium /.false./, full /.false./

CALL popt ('diagbasb', little, medium, full, ithcall, ithsub)
IF(.NOT.allocated(phipr))ALLOCATE(phipr(mxbasis,mxglegn,maxhermt,narran))
!  ----------------------------------------------------------------
!  build total primitive as a product of a legendre poly, a sho,
!  a parity function, and a jacobian factor.
!  karran runs over the quadrature arrangement channel.
!  loops 4 and 5 run over the quadrature points
!  marran runs over the function arrangement channels.
!  loops 1 and 2 run over the functions.
!  at this large rho, neglect overlaps between diff arrangements,
!  but run over all raw primitives in all arrangements to keep
!  ibasis a good index for use later in this SUBROUTINE.
!  ----------------------------------------------------------------
IF(Medium)WRITE(Out_Unit,*)'begining of routine diagbasb'
IF(nthcall==0.or.megas/=lam)THEN
   IF(symmetry)THEN
      nbasis2=nbasis-nbasiss
      nbasis1=nbasiss-nbasis2
      jpar=1
      IF(jeven) jpar=0
   ENDIF
!         WRITE(Msg_Unit,*) 'construct prim for nthcall,megas,lam=',
!                     nthcall,
!     >               megas, lam
   megas=lam
   xnorm=sqrt(half)
   nthcall=1
DO  karran = 1, narran
   DO  khermt = 1, nhermt(karran)
     DO  kglegn = 1, nglegn(karran)
         ibasis = 0
         DO  marran = 1, narran
          IF(marran==karran)THEN
           IF(symmetry.AND.marran/=1)THEN
            IF(marran==2)THEN
             DO  n = minvib(marran), maxvib(marran)
               DO  l = jmin(n,marran,lam), jmax(n,marran), jskip(marran)
                  ibasis = ibasis + 1
                  phipr(ibasis,kglegn,khermt,karran) = xnorm*oscil(n, khermt, karran) *parityfn(kglegn,khermt,karran)&
                             *bfaci(khermt,karran)*pn(l,kglegn,karran)
                ENDDO
             ENDDO
            ELSE
             DO  n = minvib(marran), maxvib(marran)
               DO  l = jmin(n,marran,lam), jmax(n,marran), jskip(marran)
                ibasis = ibasis + 1
                  phipr(ibasis,kglegn,khermt,karran) = xnorm*oscil(n, khermt, karran) *parityfn(kglegn,khermt,karran)&
                             *bfaci(khermt,karran)*pn(l,kglegn,karran) *(-1)**(l+jpar)
               ENDDO
             ENDDO

            ENDIF

           ELSE
             DO  n = minvib(marran), maxvib(marran)
               DO  l = jmin(n,marran,lam), jmax(n,marran), jskip(marran)
                  ibasis = ibasis + 1
                  phipr(ibasis,kglegn,khermt,karran) = oscil(n, khermt, karran) *parityfn(kglegn,khermt,karran)&
                             *bfaci(khermt,karran)*pn(l,kglegn,karran)
               ENDDO
             ENDDO
           ENDIF
          ELSE
            DO  n = minvib(marran), maxvib(marran)
               DO  l = jmin(n,marran,lam), jmax(n,marran), jskip(marran)
                  ibasis = ibasis + 1
                  phipr(ibasis,kglegn,khermt,karran) = zero
               ENDDO
            ENDDO
          ENDIF
      ENDDO
     ENDDO
    ENDDO
   ENDDO

IF(ibasis/=nbasis)THEN
   WRITE(Out_Unit,*)'ibasis/=nbasis',ibasis,nbasis
   WRITE(Msg_Unit,*)'ibasis/=nbasis',ibasis,nbasis
   STOP 'diagbasb: ibasis/=nbasis'
ENDIF

!  --------------------------------------------------------------
!  READ vect(nbasis,neigmin) from unit vec_unit='vector'
!  and check to be sure correct DATA is obtained.
!  --------------------------------------------------------------

14    READ(vec_unit,END=100) nbasisb, neigmin, rhob, lambda
WRITE(Out_Unit,*) 'nbasisb=',nbasisb
WRITE(Out_Unit,*) 'neigmin=',neigmin
WRITE(Out_Unit,*) 'rhob=',rhob
WRITE(Out_Unit,*) 'lambda=',lambda

IF(lambda<lam)THEN
   DO  jbasis=1,neigmin
      READ(vec_unit)
   ENDDO
   GOTO 14
ELSEIF(lambda>lam)THEN
   WRITE(Out_Unit,*)' Error, lambda>lam: lambda,lam= ', lambda,lam
   WRITE(Msg_Unit,*)' Error, lambda>lam: lambda,lam= ', lambda,lam
   STOP 'diagbasb: lambda>lam'
ENDIF
delrho=ABS((rho-rhob)/rhob)
IF(nbasiss/=nbasisb.or.delrho>1.e-08)THEN
   WRITE(Out_Unit,*) 'ERROR: nbasiss,nbasisb=',nbasiss,nbasisb
   WRITE(Out_Unit,*) 'ERROR: rho, rhob, delrho=', rho, rhob, delrho
   WRITE(Msg_Unit,*) 'ERROR: nbasiss,nbasisb=',nbasiss,nbasisb
   WRITE(Msg_Unit,*) 'ERROR: rho, rhob, delrho=', rho, rhob, delrho
   !STOP 'diagbasb nbasiss/=nbasisb.or.delrho>1.e-08' TMPMODGregParker
ENDIF
DO  jbasis=1,neigmin
   READ(vec_unit)(vect(ibasis,jbasis), ibasis=1,nbasiss)
ENDDO
!  --------------------------------------------------------------
!  on subsequent calls, make sure too many eigenfunctions are not
!  requested.
!  --------------------------------------------------------------
ENDIF
IF(tee>neigmin)THEN
   WRITE(Out_Unit,*) 'tee gt neigmin', tee, neigmin
   WRITE(Msg_Unit,*) 'tee gt neigmin', tee, neigmin
   STOP 'diagbasb: tee>neigmin'
ENDIF
!  --------------------------------------------------------------
!  transform the primitive basis to produce the surface function
!  at each of the quadrature points.
!  the special form for symmetry is necessary because the index of
!  vect THEN only runs to nbasiss.
!  --------------------------------------------------------------
DO  karran = 1, narran
   DO  khermt = 1, nhermt(karran)
      DO  kglegn = 1, nglegn(karran)
        phi(kglegn,khermt,karran)=zero
        IF(symmetry)THEN
         DO  ibasis=1, nbasis1
            phi(kglegn,khermt,karran)=phi(kglegn,khermt,karran)+phipr(ibasis,kglegn,khermt,karran)*vect(ibasis,tee)
         ENDDO
!
         DO  ibasis=nbasis1+1, nbasiss
            phi(kglegn,khermt,karran)=phi(kglegn,khermt,karran)+(phipr(ibasis,kglegn,khermt,karran)&
            +phipr(ibasis+nbasis2,kglegn,khermt,karran))*vect(ibasis,tee)
         ENDDO
!
        ELSE
         DO  ibasis=1, nbasis
            phi(kglegn,khermt,karran)=phi(kglegn,khermt,karran)+phipr(ibasis,kglegn,khermt,karran)*vect(ibasis,tee)
         ENDDO
        ENDIF
     ENDDO
   ENDDO
    ENDDO
IF(little)THEN
   WRITE(Out_Unit,*)' nbasis,nbasiss = ', nbasis,nbasiss
ENDIF

IF(Medium)WRITE(Out_Unit,*)'END of routine diagbasb'
RETURN
100   WRITE(Msg_Unit,*)' Error, EOF detected in file vector in routine DiagBasb'

STOP 'diagbasb: EOF when reading Vec_Unit'
END
