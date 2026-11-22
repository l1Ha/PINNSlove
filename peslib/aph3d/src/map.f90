SUBROUTINE map (iel, ielo, ielr, ielor, id, ido, mapnod, idiv, idivo, mxelmnt, symmetry)
USE FileUnits_Module
!
!     $RCSfile: map.f,v $ $Revision: 1.14 $
!     $Date: 89/10/18 14:18:06 $
!     $State: Stable $
!

!     P U R P O S E   O F    S U B R O U T I N E

!     this routine actually maps the old mesh onto the new mesh; that is IF
!     the new node number is put into mapnod, THEN you get out the
!     corresponding old node number. the mapping is done by cycling
!     through the elements.  note that there are three main cases:
!     1) old and new elements are the same size, 2) old element is smaller
!     than the new element, and 3) old element is bigger than the new elemen
!     we have allowed the smaller element to subdivided up to three more
!     times than the larger element.


!     I N P U T    A R G U M E N T S
!
!     iel
!     ielo
!     eilr
!     ielor
!     id(numnp)   id=-nce if the is constrained and positive otherwise
!     ido
!     mapnod
!     idiv
!     idivo
!     mxelmnt
!     symmetry


IMPLICIT NONE

!     L O G I C A L S

LOGICAL symmetry

!     I N T E G E R S



INTEGER iel, ielo, ielr, ielor, ido, idiv, idivo, mapnod, mxelmnt, n1map, n1mapr, n2mapr, n2map, n2mapo, n2mapu, &
           kold, irefl, kold2, koldr2, k, nquad, kstar, kstarr, id,  mel, iiel

!     D I M E N S I O N S

DIMENSION n1map(4), n1mapr(4), n2mapr(4), n2map(4), n2mapo(4), n2mapu(4), kold(9,4), idiv(mxelmnt), idivo(mxelmnt), &
           id(9,mxelmnt), irefl(4), mapnod(*), ido(9,mxelmnt), kold2(4,4), koldr2(4,4)



DATA n1map /1,2,4,5/
DATA n1mapr /4,5,1,2/, n2mapr /7,9,1,3/
DATA n2map /0,1,3,4/, n2mapo /1,3,7,9/, n2mapu /3,4,0,1/
DATA kold /1,1,2,1,5,5,4,5,5,2,3,3,5,5,3,5,5,6,4,5,5,7,5,5,7,7,8,5,5,6,5,5,9,8,9,9/
DATA irefl /3,4,1,2/
DATA kold2 /1,2,4,5,2,3,5,6,4,5,7,8,5,6,8,9/
DATA koldr2 /7,8,4,5,8,9,5,6,4,5,1,2,5,6,2,3/

!     WRITE(Out_unit,*)'idiv=',idiv(iel),' for iel=',iel
!     WRITE(Out_unit,*)'idivo=',idivo(ielo),' for ielo=',ielo
IF(idiv(iel)==idivo(ielo))THEN

!     the elements are the same size

   DO 10 k=1,9
      mapnod(id(k,iel))=ido(k,ielo)
      mapnod(id(k,ielr))=ido(k,ielor)
 10      CONTINUE
   iel=iel+1
   ielo=ielo+1
   ielr=ielr+1
   ielor=ielor+1

ELSEIF (idiv(iel)<idivo(ielo))THEN

!     the old element is smaller than the new element
!     loop thru the four quadrants of the new element

   DO 40 nquad=1,4

!     find the starting k for this quadrant of the new element; each
!     quadrant is defind by four nodes, the starting k tells us which four.
      kstar=n1map(nquad)
      kstarr=n1mapr(nquad)

      IF(idivo(ielo)-idiv(iel)==1)THEN

!     there is only one old element in this quadrant of the new element

         mapnod(id(kstar,iel))=ido(1,ielo)
         mapnod(id(kstar+1,iel))=ido(3,ielo)
         mapnod(id(kstar+3,iel))=ido(7,ielo)
         mapnod(id(kstar+4,iel))=ido(9,ielo)
         IF(symmetry)THEN
!     notice that the local node numbers stay the same, but the elements are
!     numbered in reverse order upon reflection
            mapnod(id(kstarr,ielr))=ido(1,ielor)
            mapnod(id(kstarr+1,ielr))=ido(3,ielor)
            mapnod(id(kstarr+3,ielr))=ido(7,ielor)
            mapnod(id(kstarr+4,ielr))=ido(9,ielor)
         ELSEIF(.NOT.symmetry)THEN
!     both node numbers and element numbers are unchanged upon adding pi
            mapnod(id(kstar,ielr))=ido(1,ielor)
            mapnod(id(kstar+1,ielr))=ido(3,ielor)
            mapnod(id(kstar+3,ielr))=ido(7,ielor)
            mapnod(id(kstar+4,ielr))=ido(9,ielor)
         ENDIF
         ielo=ielo+1
         ielor=ielor+1

      ELSE

!     there are at least 4 elements, maybe 16, in this quadrant

         DO 30 mel=1,4

!     loop thru the 4 parts of this quadrant

            IF(idivo(ielo)-idiv(iel)==2)THEN

!     this element has only been subdivided twice

               mapnod(id(kstar+n2map(mel),iel))=ido(n2mapo(mel),ielo)
               IF(symmetry)THEN
                  mapnod(id(kstarr+n2mapu(mel),ielr))=ido(n2mapr(mel),ielor)
               ELSEIF(.NOT.symmetry)THEN
                  mapnod(id(kstar+n2map(mel),ielr))=ido(n2mapo(mel),ielor)
               ENDIF
               ielo=ielo+1
               ielor=ielor+1
            ELSE

!     this element has been subdivided three times, so we have 4 elements
!     to map here. Since only one old element touches an new node
!     we only have to map one node; we DO need to check to make sure
!     there are no further subdivisions

               mapnod(id(kstar+n2map(mel),iel))=ido(n2mapo(mel),ielo+mel-1)
               IF(symmetry)THEN
                  mapnod(id(kstarr+n2mapu(mel),ielr))=ido(n2mapr(mel),ielor+mel-1)
               ELSEIF(.NOT.symmetry)THEN
                  mapnod(id(kstar+n2map(mel),ielr))=ido(n2mapo(mel),ielor+mel-1)
               ENDIF
               DO 20 iiel=0,3
                  IF(idivo(ielo+iiel)-idiv(iel)/=3)THEN
                     WRITE(Out_unit,*) '***error***: too much change',' in subdivisions to map'
                     WRITE(Out_unit,*) 'idivo=',idivo(ielo+iiel),' for element',ielo+iiel
                     WRITE(Out_unit,*) 'idiv=',idiv(iel),' for element 1',iel
                     STOP 'map'
                  ENDIF
!     IF(symmetry.AND.idivo(ielor+iiel)-idiv(ielr)/=3)THEN
!     WRITE(Out_unit,*) '***error***: too much change in subdivisions to map'
!     WRITE(Out_unit,*) 'idivo=',idivo(ielor+iiel),' for element',ielor+iiel
!     WRITE(Out_unit,*) 'idiv=',idiv(ielr),' for element',ielr,' upper half 1
!     STOP 'map'
!     ENDIF
 20                  CONTINUE
               ielo=ielo+4
               ielor=ielor+4
            ENDIF
 30            CONTINUE
      ENDIF
 40      CONTINUE
   iel=iel+1
   ielr=ielr+1
ELSE

!     the new elements are smaller than the old element

   DO 100 nquad=1,4

!     loop thru the four quadrants of the old element

      kstar=n1map(nquad)
      kstarr=n1mapr(nquad)

      IF(idiv(iel)-idivo(ielo)==1)THEN

!     there is only one new element in this quadrant of the old element

         DO 50 k=1,9
            mapnod(id(k,iel))=ido(kold(k,nquad),ielo)
            IF(symmetry)THEN
               mapnod(id(k,ielr))=ido(kold(k,irefl(nquad)),ielor)
            ELSEIF(.NOT.symmetry)THEN
               mapnod(id(k,ielr))=ido(kold(k,nquad),ielor)
            ENDIF
 50            CONTINUE
         iel=iel+1
         ielr=ielr+1
      ELSE

!     there are at least 4 elements, maybe 16, in this quadrant

         DO 90 mel=1,4
            IF(idiv(iel)-idivo(ielo)==2)THEN

!     this element has only been subdivided twice

               DO 60 k=1,9
                  mapnod(id(k,iel))=ido(kold2(mel,nquad),ielo)
                  IF(symmetry)THEN
                     mapnod(id(k,ielr))= ido(koldr2(mel,nquad),ielor)
                  ELSEIF(.NOT.symmetry)THEN
                     mapnod(id(k,ielr))=ido(kold2(mel,nquad),ielor)
                  ENDIF
 60                  CONTINUE
               iel=iel+1
               ielr=ielr+1
            ELSE

!     this element has been subdivided at least three times;
!     make sure it has not been divided more than that

               DO 80 iiel=1,4
                  IF(idiv(iel)-idivo(ielo)/=3)THEN
                     WRITE(Out_unit,*) '***error***: too much change',' in subdivisions to map'
                     WRITE(Out_unit,*) 'idivo=',idivo(ielo), ' for element',ielo
                     WRITE(Out_unit,*) 'idiv=',idiv(iel), ' for element 2',iel
                     STOP 'map'
                  ENDIF
!     IF(symmetry.AND.idiv(ielr)-idivo(ielor)/=3)THEN
!     WRITE(Out_unit,*) '***error***: too much change in subdivisions to map'
!     WRITE(Out_unit,*) 'idivo=',idivo(ielor),' for element',ielor
!     WRITE(Out_unit,*) 'idiv=',idiv(ielr),' for element',ielr,' upper half 2
!     STOP 'map'
!     ENDIF
                  DO 70 k=1,9
                     mapnod(id(k,iel))=ido(kold2(mel,nquad),ielo)
                     IF(symmetry)THEN
                        mapnod(id(k,ielr))=ido(koldr2(mel,nquad),ielor)
                     ELSEIF(.NOT.symmetry)THEN
                        mapnod(id(k,ielr))=ido(kold2(mel,nquad),ielor)
                     ENDIF
 70                     CONTINUE
                  ielr=ielr+1
                  iel=iel+1
 80                  CONTINUE
            ENDIF
 90            CONTINUE
      ENDIF
 100     CONTINUE
   ielo=ielo+1
   ielor=ielor+1
ENDIF
RETURN
END
