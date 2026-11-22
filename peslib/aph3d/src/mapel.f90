SUBROUTINE mapel(id, lcor, mxelmnt, idiv, rho, nelem, symmetry, &
                 ndiv, numnp, maxelp, mega, idatwr, ido, idivo, &
                 lcoro, mapnod)
USE restar_Module
USE FileUnits_Module
!
!     $RCSfile: mapel.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:07 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E

!     this routine maps the old mesh into the new mesh when the  mesher
!     is run at least twice for different rho values


!     I N P U T    A R G U M E N T S
!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     lcor
!     mxelmnt
!     idiv
!     rho
!     nelem
!     symmetry
!     ndiv
!     numcp
!     maxelp
!     mega
!     idatwr
!     ido
!     idivo
!     lcoro
!     mapnod


IMPLICIT NONE
LOGICAL first, symmetry
INTEGER mxelmnt, nelem, ndiv, numnp, maxelp,  mega, idatwr
INTEGER iel, k, i, nrd, nelemo, ndivo, numnpo
INTEGER ielo, ielr, ielor, nnnn,  icor, jcor, icoro, jcoro, n
INTEGER id(9,nelem), idiv(nelem), lcor(2,nelem), ido(9,nelem)
INTEGER idivo(nelem), lcoro(2,nelem), mapnod(numnp)
REAL(Kind=WP_Kind) rho

SAVE first
DATA first /.true./

IF(lrestar)first=.false.
IF(first)THEN
   first=.false.
   REWIND Zromsh_Bin_Unit
   WRITE(Zromsh_Bin_Unit) nelem, ndiv, numnp, ((id(k,iel),k=1,9),iel=1,nelem),   &
              ((lcor(i,iel),i=1,2),iel=1,nelem),(idiv(iel),iel=1,nelem)
!     WRITE(Out_unit,*)'nelem=',nelem
!     WRITE(Out_unit,*)'to Zromsh_Bin_Unit; idiv=',(idiv(iel), iel=1,nelem)

ELSE

!     map the current mesh to the last mesh

!     IF(idatwr<1)THEN
   IF(idatwr/=0)THEN
      WRITE(Out_unit,*)'id='
      WRITE(Out_unit,*)((id(k,iel), k=1,9), iel=1,nelem)
   ENDIF

!     IF mega=0, THEN READ the old mesh from the zero mesh
!     file, not from the last mesh file

   IF(mega==0)THEN
      nrd=Zromsh_Bin_Unit
   ELSE
      nrd=OldMsh_Unit
   ENDIF

!     READ in the old mesh DATA

   REWIND nrd
   READ(nrd) nelemo,ndivo,numnpo,((ido(k,iel),k=1,9),iel=1,nelemo),  &
              ((lcoro(i,iel),i=1,2),iel=1,nelemo),(idivo(iel),iel=1,nelemo)

!     WRITE(Out_unit,*)'nelemo=',nelemo
!     WRITE(Out_unit,*)'READ from nrd=',nrd,
!     &          ' idivo=',(idivo(iel), iel=1,nelemo)

!     save the current mega=0 mesh DATA so we can map
!     it to the next mega=0 mesh

   IF(mega==0)THEN
      REWIND Zromsh_Bin_Unit
      WRITE(Zromsh_Bin_Unit) nelem,ndiv,numnp,((id(k,iel),k=1,9),iel=1,nelem),  &
                 ((lcor(i,iel),i=1,2),iel=1,nelem),(idiv(iel),iel=1,nelem)
!     WRITE(Out_unit,*)'nelem=',nelem
!     WRITE(Out_unit,*)'to Zromsh_Bin_Unit: idiv=',(idiv(iel), iel=1,nelem)
   ENDIF

!     IF the range of chi was doubled, THEN we only have to
!     loop thru the first half of the elements

!     SUBROUTINE double flips the last element and the element
!     with the largest chi (in the array id only); we have to
!     flip them back to be consistent

   IF(symmetry)CALL swap (nelem, maxelp, lcor, idiv, id)
   nelem=nelem/2
   nelemo=nelemo/2

   iel=1
   ielo=1
   ielr=1+nelem
   ielor=1+nelemo
   DO 10 nnnn=1,mxelmnt

!     find the corner (upper left) and the
!     length in i&j to the nodes on the side  of the two elements

      icor=lcor(1,iel)
      jcor=lcor(2,iel)
      icoro=lcoro(1,ielo)
      jcoro=lcoro(2,ielo)

!     IF we changed the number of subdividions allowed , take care of it

      IF(ndiv==ndivo+1)THEN
         icor=2*icor-1
         jcor=2*jcor-1
      ELSEIF (ndiv==ndivo-1)THEN
         icoro=2*icoro-1
         jcoro=2*jcoro-1
      ELSEIF (ndiv/=ndivo)THEN
         WRITE(Out_unit,*) '***error***: ndiv changed by more than one'
         WRITE(Out_unit,*) 'ndiv=',ndiv,' ndivo=',ndivo
         STOP 'mapel'
      ENDIF

      IF(icor/=icoro.or.jcor/=jcoro)THEN
         WRITE(Out_unit,*) '***error***: the base elements DO not have the same corner'
         WRITE(Out_unit,*) 'icor=',icor,' jcor=',jcor,' icoro=', icoro,' jcoro=',jcoro
         STOP 'mapel'
      ENDIF

!     now sort the elements according to which is bigger

      CALL map (iel, ielo, ielr, ielor, id, ido, mapnod, idiv, idivo, mxelmnt, symmetry)

      IF(iel>nelem)THEN
         IF(ielo>nelemo)THEN
            GOTO 20
         ELSE
            WRITE(Out_unit,*) '***error***: elements out of phase'
            WRITE(Out_unit,*) 'nelem=',nelem,' nelemo=',nelemo
            STOP 'mapel'
         ENDIF
      ENDIF
 10      CONTINUE
 20      CONTINUE

!     put the id arrarys back how we got them

   nelem=nelem*2
   IF(symmetry)CALL swap (nelem, maxelp, lcor, idiv, id)

!     WRITE out the DATA for the next rho value and for sfunc

   REWIND Mapmsh_Bin_Unit
   WRITE(Mapmsh_Bin_Unit) rho,mega,numnp,numnpo,(mapnod(n),n=1,numnp)
!     IF(idatwr<1)THEN
   IF(idatwr/=0)THEN
      WRITE(Out_unit,*)'rho=',rho,' numnp=',numnp,'map is'
      WRITE(Out_unit,*)(mapnod(n), n=1,numnp)
   ENDIF
ENDIF

REWIND OldMsh_Unit
WRITE(OldMsh_Unit)nelem, ndiv, numnp, ((id(k,iel),k=1,9),iel=1,nelem),  &
        ((lcor(i,iel),i=1,2),iel=1,nelem),(idiv(iel), iel=1,nelem)
!     WRITE(Out_unit,*)'nelem=',nelem
!     WRITE(Out_unit,*)'to OldMsh_Unit; idiv=',(idiv(iel), iel=1,nelem)
RETURN
ENDSUBROUTINE Mapel
