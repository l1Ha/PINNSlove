SUBROUTINE sfunc(id, theta, chi, phi, yz, lm, nid, idi, beta, nel, nodj, maxa, mht, numnp, ndisce, nidm, neq, nfreq,   &
           nwk, nnc, nc, hamil, sOvrlp, ar, br, vec, eigv,  d, tt, w, nloc, rtolv, evc1, evc2, freq, ww,               &
           nsit, nlancz, nblock, maxj, nwork, mxnode, pmax)
USE Numeric_Kinds_Module
USE FileUnits_Module
USE var_Module
USE mshdat_Module,only:aterm,bterm
USE Numbers_Module
IMPLICIT NONE
SAVE
!
!     $RCSfile: sfunc.f,v $ $Revision: 1.15 $
!     $Date: 89/11/16 11:03:31 $
!     $State: Stable $
!
!     P U R P O S E   O F    S U B R O U T I N E
!     I N P U T    A R G U M E N T S

!-----------------------------------------------------------------------
!     This routine calculates the surface functions.
!     theta     theta values of the nodal points
!     chi       chi values of the nodal points.
!     phi       surface function values at the nodal points.
!     beta      constraint coefficients.
!     mht       column heights.
!     maxa      addresses of the diagonal elements.
!     numnp     number of nodal points.
!     ndisce    number of displacement constraints.
!     neq       number of independent equations.
!-----------------------------------------------------------------------
LOGICAL little, medium, full
INTEGER ithcall, ithsub, lm, nid, nel, nodj, maxa, mht, numnp, ndisce, nidm, neq, nfreq, nwk, nnc, nc ,tconv
INTEGER maxw, maxiw, errno, id, nsit, nlancz, nblock, maxj, pmax,  nwork, mxnode, numest, ite, isub, nstape
INTEGER maxest,  nume, i, j, idi, nd, nod, nloc, ma, dummy, ar, nextra
REAL(Kind=WP_Kind) theta, chi, phi, yz, beta, a, b, second, tmstlm, t
REAL(Kind=WP_Kind) hamil, sOvrlp, br, eigv, d, tt, w, rtolv,  evc1, evc2, ww, vec, freq, s, xx, ctime

DIMENSION id(numnp), theta(numnp), chi(numnp), phi(numnp,nc), yz(18,nel), lm(18,nel), nid(ndisce)
DIMENSION idi(nidm,2*ndisce), beta(nidm,2*ndisce), nodj(nidm),  maxa(neq+1), mht(neq) , nod(9), s(378), xx(18)
DIMENSION hamil(nwk), sOvrlp(nwk), ar(nnc), br(nnc), vec(nc,nc), eigv(nc), d(nc), tt(neq), w(neq), nloc(neq)
DIMENSION rtolv(nc), evc1(nc), evc2(mxnode), freq(nc), ww(neq),  nsit(nc)
!     CPU1 The temporary variable used to hold the accumulated CPU time
!     at the start of a procedure to be timed.
REAL(Kind=WP_Kind) cpu1

!     CPU2 The temporary variable used to hold the accumulated CPU time
!     at the END of a procedure to be timed.
REAL(Kind=WP_Kind) cpu2

!     OVERHEAD The overhead (in cpu seconds) of the timing procedure
!     CPUTIME.
REAL(Kind=WP_Kind) overhead

DATA ithcall/0/, ithsub/0/
DATA little/.false./, medium/.false./, full/.false./
CALL popt ('sfunc   ', little, medium, full, ithcall, ithsub)

!-----------------------------------------------------------------------
!     IF nwk is zero the surface function DATA needs to be READ in and
!     nwk needs to be calculated.
!-----------------------------------------------------------------------
WRITE(Msg_Unit,*)'start sfunc'
IF(nwk==0)THEN
!     Zero the pointer arrays so that values stored previously at those
!     addresses DO not destroy the logic in SUBROUTINE colht.
   CALL nsets(neq+1, 0, 0, maxa, 1)
   CALL nsets(neq, 0, 0, mht, 1)
 10      numest=0
   maxest=0
   ite=0
   isub=0
   nstape=14
!-----------------------------------------------------------------------
!     READ finite element mesh DATA.
!-----------------------------------------------------------------------
   CALL adini (id, theta, chi, nid, idi, beta, nodj, nume)
   Full=.true.
   IF(full)THEN
      WRITE(Out_unit,*) 'after adini'
      WRITE(Out_unit,'(A,/,(10I11))') 'id =',(id(i),i=1,numnp)
      WRITE(Out_unit,'(A,/,(10ES11.3))') 'theta =',(theta(i),i=1,numnp)
      WRITE(Out_unit,'(A,/,(10ES11.3))') 'chi =',(chi(i),i=1,numnp)
      WRITE(Out_unit,'(A,/,(10I11))') 'nid =',(nid(i),i=1,ndisce)
      WRITE(Out_unit,'(A,/,(10I11))') 'idi =',((idi(i,j),j=1,ndisce),i=1,nidm)
      WRITE(Out_unit,'(A,/,(10ES11.3))') 'beta =',((beta(i,j),j=1,ndisce),i=1,nidm)
      WRITE(Out_unit,'(A,/,(10I11))') 'nodj =',(nodj(i),i=1,nidm) !tempmod should this be nidm or numnp?
   ENDIF
!-----------------------------------------------------------------------
!     obtain element information.
!-----------------------------------------------------------------------
   WRITE(Msg_Unit,*)'tdfe pos1'
   CALL tdfe (id, theta, chi, mht, lm, yz, nd, neq, nume, nid, nidm, idi, 0, hamil, nod, s, xx, beta, maxa, nwk )

   IF(full)THEN
      WRITE(Out_unit,*) 'after tdfe'
      WRITE(Out_unit,'(A,/,(10ES11.3))') 'theta =',(theta(i),i=1,numnp)
      WRITE(Out_unit,'(A,/,(10ES11.3))') 'chi =',(chi(i),i=1,numnp)
      WRITE(Out_unit,'(A,/,(10I11))') 'mht =',(mht(i),i=1,neq)
      WRITE(Out_unit,'(A,/,(10I11))') 'lm =',((lm(i,j),j=1,nel),i=1,18)
      WRITE(Out_unit,'(A,/,(10ES11.3))') 'yz =',((yz(i,j),j=1,nel),i=1,18)
      WRITE(Out_unit,'(A,/,(10I11))') 'nid =',(nid(i),i=1,ndisce)
      WRITE(Out_unit,'(A,/,(10I11))') 'idi =',((idi(i,j),j=1,ndisce),i=1,nidm)
      WRITE(Out_unit,'(A,/,(10ES11.3))') 'beta =',((beta(i,j),j=1,ndisce),i=1,nidm)
   ENDIF
!-----------------------------------------------------------------------
!     Compute the addresses of the diagonal elements.
!     These addresses are stored in the array maxa.
!     The number of elements below skyline is nwk.
!-----------------------------------------------------------------------
   CALL addres (maxa, mht, neq, nwk, ma)
   IF(full)THEN
      WRITE(Out_unit,*) 'after addres'
      WRITE(Out_unit,*) 'nwk =',nwk
      WRITE(Out_unit,'(A,/,(20I6))') 'mht =',(mht(i),i=1,neq)
      WRITE(Out_unit,'(A,/,(20I6))') 'maxa =',(maxa(i),i=1,neq+1)
   ENDIF
ELSE
!-----------------------------------------------------------------------
!     Surface function DATA has been READ in and storage has been
!     allocated.
!     This section produces matrix elements and produces
!     the surface functions and eigenenergies.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Calculate the Hamiltonian matrix elements.
!-----------------------------------------------------------------------

   DO i=1,nwk
      hamil(i)=zero
   ENDDO
   WRITE(Msg_Unit,*)'tdfe pos2'
   CALL tdfe (id, theta, chi, mht, lm, yz, nd, neq, nume, nid, nidm, idi, 1, hamil, nod, s, xx, beta, maxa, nwk )
   IF(full)THEN
      WRITE(Out_unit,*)' Hamiltonian matrix'
      WRITE(Out_unit,'(10ES11.3)')(hamil(i),i=1,nwk)
   ENDIF

!-----------------------------------------------------------------------
!     Calculate the Overlap matrix elements.
!-----------------------------------------------------------------------
   DO i=1,nwk
      sOvrlp(i)=zero
   ENDDO
   WRITE(Msg_Unit,*)'tdfe pos3'
   CALL tdfe (id, theta, chi, mht, lm, yz, nd, neq, nume, nid, nidm, idi, 2, sOvrlp, nod, s, xx, beta, maxa, nwk )
   IF(full)THEN
      WRITE(Out_unit,*)' Overlap matrix'
      WRITE(Out_unit,'(10ES11.3)')(sOvrlp(i),i=1,nwk)
   ENDIF
!-----------------------------------------------------------------------
!     Diagonalize the hamiltonian using the subspace iteration
!     technique.
!     The eigenvectors are the values of the surface functions at the
!     nodal points.
!-----------------------------------------------------------------------
      IF(little.AND.nlancz==0)THEN
         WRITE(Out_unit,*)'calling sspace'
         WRITE(Out_unit,*)'sfunc: neq, nc, nfreq, nwk = ',neq, nc, nfreq, nwk
      ENDIF
      IF(little.AND.nlancz==1)THEN
         WRITE(Out_unit,*)'calling lanczo'
         WRITE(Out_unit,*)'sfunc: neq, nc, nfreq, nwk, nblock, maxj = ', neq, nc, nfreq, nwk, nblock, maxj
      ENDIF
   REWIND Hamil_FEM_Unit
   WRITE(Hamil_FEM_Unit)hamil
   REWIND SOvrlp_Unit
   WRITE(SOvrlp_Unit)sOvrlp

!     -------------------------
!     Get overhead for CPUTIME.
!     -------------------------

   CALL cputime (cpu1)
   CALL cputime (cpu2)
   overhead = cpu2 - cpu1

   CALL cputime (cpu1)
   REWIND Tidi_Unit
!     Use idi as a temporary storage vector in sspace.
   WRITE(Out_unit,*)'nidm,ndisce,numnp=',nidm,ndisce,numnp
   WRITE(Tidi_Unit) ((idi(i,j),i=1,nidm),j=1,numnp)
         IF(nlancz==0)THEN
!         CALL sspace (maxa, hamil, sOvrlp, ar, br, vec, eigv, d, tt,  w, nloc, rtolv, evc1, evc2, phi, freq, ww, nsit, neq,   &
!              neq, nc, nfreq, nwk, numnp, mxnode, id, idi)

     ELSEIF(nlancz==1)THEN
!     tt, eigv, vec, ar, br, d, and nloc are all work vectors
!     and the names have no significance.  See lanczo.

     nextra = nc - nfreq
     STOP "lanczo not working"
     !CALL lanczo (maxa, hamil, sOvrlp, phi, freq, neq, nfreq, nwk, nblock, maxj, nwork,t, ar, br, d, nloc, numnp, mxnode,    &
     !              id, idi, evc2, nc)
         ELSEIF(nlancz==2)THEN

!---------------------------------------------------------------------
!  ******************set parameters for running STLM******************
!---------------------------------------------------------------------

!  look for eigenvalues on the interval (a,b) where
      a=5.d+2
      b=1.d+5
      a=aterm
      b=bterm
!  ***memstore and mem_sfunc must have values consistent with this.***
!  stlm will not run IF maxw is too small.
      maxw=max(3*nwk,2*nwk+(2*pmax+5)*neq+pmax*(pmax+5))
      maxiw=max(2*neq,2*neq+6*pmax +4)
      WRITE(Out_unit,*) 'maxw, maxiw,nwk,pmax,neq=', maxw, maxiw,nwk,pmax,neq

      tmstlm=second()
!  -------------------------------------------------------------------
!  ****CALL STLM driver*****
!  -------------------------------------------------------------------

      CALL stlmdr (neq, maxa, hamil, sOvrlp, ar, a, b, tconv, freq, phi, nc, maxw, maxiw, errno, pmax)

      tmstlm=second()-tmstlm
      WRITE(Out_unit,*) 'Time in STLM=',tmstlm
         ENDIF
   REWIND Tidi_Unit
   READ(Tidi_Unit) ((idi(i,j),i=1,nidm),j=1,numnp)

   CALL cputime (cpu2)
   ctime = cpu2 - cpu1 - overhead
   IF(little.AND.nlancz==0)THEN
      WRITE(Out_unit, '(''total time in sspace = '', f10.3, '' seconds'')') ctime
   ENDIF
   IF(little.AND.nlancz==1)THEN
      WRITE(Out_unit, '(''total time in lanczo = '', f10.3, '' seconds'')') ctime
   ENDIF

!-----------------------------------------------------------------------
!     WRITE out the surface functions.
!-----------------------------------------------------------------------
   CALL wrmod (freq, hamil, idi, numnp, neq, nc, nc, phi, beta, nid, nidm, ndisce, id, nlancz, nel)
ENDIF

RETURN
END
