SUBROUTINE STLMDR (NEQ, MAXA, HAMIL, SOVL, IW, A, B, TCONV,  FREQ, PHI, MAXL, MAXW, MAXIW, ERRNO, PMAX)
USE FileUnits_Module
USE Numeric_Kinds_Module
!
!   $RCSfile: stlmdr.f,v $ $Revision: 1.3 $
!   $Date: 89/07/28 09:57:44 $
!   $State: Stable $
!
!
! **********************************************************************
!
!     PURPOSE -
!
!         THE FOLLOWING SUBROUTINE IS A DRIVER FOR THE STLM-PACKAGE
!         OF RUHE AND ERICHSON. IT TAKE THE K AND M MATRICES STORED
!         IN SKYLINE FORM AND TRANSFORM THEM TO PROFILESTORAGE IN ORDER
!         TO GET THE RIGHT INPUT TO STLM.
!
!     INPUT
!
!         NEQ     DIMENSION OF THE K AND M MATRICES
!         MAXA    INTEGERARRAY CONTAINING THE POSITIONS OF DIAGONAL-
!                 ELEMENTS IN HAMIL AND SOVL
!         HAMIL   THE HAMILTON-MATRIX PACKED IN SKYLINE FORM
!         SOVL    THE OVERLAB-MATRIX PACKED IN SKYLINE FORM
!         PHI     THE EIGENVECTORS
!         W       WORKARRAY NEEDED FOR STLM. IN THIS VERSION W IS
!                 SIMPLY THE COMBINED SPACES (HAMIL,SOVL,AR)
!         MAXW    DIMENSION OF W MUST BE AT LEAST
!                 2*(MAXA(NEQ+1)-1)+2*NEQ+PMAX*(5+PMAX)+MAXL
!         IW      INTEGER WORKARRAY NEEDED FOR STLM
!         MAXIW   DIMENSION OF IW MUST BE AT LEAST
!                 NEQ+2*PMAX+MAXL+2
!         A,B     STLM WILL FIND ALL EIGENVALUES IN [A,B]
!         MAXL    MAXIMUM EIGENVALUES WANTED
!         PMAX    MAXIMUM NUMBER OF LANCZOS ITERATIONS ALLOWED.  MUST
!                 BE AT LEAST 15, BUT SHOULD BE 50, 100, OR MORE FOR
!                 LARGE NEQ.
!         PROFILE = 1 FOR NON-DIAGONAL SOVL.
!         MXREST  LIMITS NO. OF EIGENVALUES REMAINING TO BE COMPUTED IN
!                 AN INTERVAL.  INITIALLY SET NEGATIVE.  PROGRAM RESETS.
!         MSGLVL  MESSAGE LEVEL.  CONTROLS AMOUNT OF OUTPUT WRITTEN.
!
!     OUTPUT -
!
!         TCONV   NUMBER OF CONVERGED EIGENVALUES
!         W       (1:TCONV) CONTAINS THE EIGENVALUES
!         PHI     (1:TCONV,1:NEQ) CONTAINS THE EIGENVECTORS
!         ERRNO   TELL IF THE RUN HAS BEEN COMPLEETED SUCCESSFULLY
!                 ==0 YES ; /=0 NO AND TELL WHICH ERROR HAS OCCURED.
!
! **********************************************************************
!
SAVE
INTEGER, PARAMETER:: NORST = 2
INTEGER   NEQ, MAXA(400), IW(400), DAFILE, KFILE, HSFILE,   MAXIW,  MAXL,  MAXREC,  MAXW,  MSGLVL,  MXREST,    & !TMPModGregParker
          NERR, NLEFT, NOUT, PMAX, PROFIL, TCONV, X,   ERRNO , I, J, IDUMMY, IST,LAST, INORST
REAL(Kind=WP_Kind)    A, B, HAMIL(2*7318), SOVL(2*7318), PHI(NEQ,MAXL), BG, FREQ(2*7318)   !TMPModGregParker

COMMON   /STLMPR/ IDUMMY, NERR, NOUT
!
!     *********************************************
!     TRANSFORM HAMIL OG SOVL INTO PROFILE-STORAGE.
!     *********************************************
!
!      Write(Msg_Unit,*) neq,maxl, maxw,maxiw
INORST = 0
3     CONTINUE
REWIND (Hamil_FEM_Unit)
READ(Hamil_FEM_Unit) ((HAMIL(J),J=MAXA(I+1)-1,MAXA(I),-1),I=1,NEQ)
REWIND (SOvrlp_Unit)
READ(SOvrlp_Unit) ((SOVL(J),J=MAXA(I+1)-1,MAXA(I),-1),I=1,NEQ)
DO I=1,NEQ
  IW(I)=MAXA(I+1)-1
ENDDO
!
!     ****************************************
!     INITIALIZE THE INPUT PARAMETERS TO STLM.
!     ****************************************
!
PROFIL = 1
!     PMAX   = 100
MXREST = -1
MSGLVL = 0
MAXREC = 5*PMAX + MAXL + 2
!
!     *****************************************
!     THE USER MAY HAVE TO CHANGE THESE VALUES.
!     PLEASE SEE THE USER GUIDE.
!     *****************************************
!
DAFILE = 91
KFILE  = 92
!
!     ************************************************************
!     OPEN DAFILE AND KFILE (FORTRAN77). PLEASE SEE THE USER GUIDE
!     AND THE LOCAL INSTALLATION DOCUMENT FOR MORE DETAILS.
!     NOTE THAT RECL IS COMPUTER AND INSTALLATION DEPENDENT.
!     ************************************************************
!
OPEN(KFILE,File= OutDIR(1:LEN(TRIM(OutDIR)))//'KSCR_STLM.bin', FORM = 'UNFORMATTED',  ACCESS = 'SEQUENTIAL',STATUS = 'UNKNOWN')
!
OPEN(DAFILE,File= OutDIR(1:LEN(TRIM(OutDIR)))//'DASCR_STLM.bin', FORM = 'UNFORMATTED',  ACCESS = 'DIRECT',STATUS = 'UNKNOWN', RECL = 8*NEQ)
!
!     **********
!     CALL STLM.
!     **********
!
!      Write(Msg_Unit,*) neq,maxl,maxw,maxiw
CALL STLM(NEQ, A, B, MAXL, PROFIL, PMAX, MXREST, MSGLVL, MAXW,  MAXIW, DAFILE, MAXREC, KFILE, X, BG, TCONV, NLEFT,   &
                ERRNO, HAMIL, IW)
!     **************************
!     CHECK IF UNEXPECTED ERROR.
!     **************************
!
!       WRITE(0, *) 'IN STLMDR, ERRNO=', ERRNO

IF(.NOT. ((ERRNO/=1701.AND.ERRNO/=0).OR. TCONV == 0)) GOTO 20
      INORST = INORST + 1
  PMAX = PMAX + 10
  WRITE(NERR, 200)
  WRITE(NERR, *) ERRNO
  WRITE(0, *) 'IN STLMDR, ERRNO=', ERRNO
  WRITE(0, 200)
200     FORMAT(//1X, "***** UNEXPECTED ERROR.")
!
  CLOSE(KFILE, STATUS='DELETE')
  CLOSE(DAFILE, STATUS='DELETE')
  IF(INORST > NORST ) STOP 'ERROR in STLM'
  GOTO 3
!
!
!     *****************
!     WRITE EIGENPAIRS.
!     *****************
!
20    CONTINUE
REWIND (Temp2_FEM_Unit)
DO 30 I = 1, TCONV
!
  FREQ(I)=HAMIL(I)
!       WRITE(NOUT, *) I, HAMIL(I),HAMIL(I)/7468.125998d0
!       **********************************************************
!       THE FOLLOWING FORTRAN77 STATEMENT MAY HAVE TO BE REPLACED.
!       THE STATEMENT SHOULD READ THE EIGENVECTOR Z, OF LENGTH N,
!       WITH INDEX = IW(I), FROM THE DIRECT ACCESS FILE.
!       **********************************************************
!
  READ(DAFILE, REC = IW(I)) (PHI(J,I),J=1,NEQ)
  WRITE(Temp2_FEM_Unit) (PHI(J,I),J=1,NEQ)
!
!       WRITE(NOUT, 300) I, HAMIL(I)
!300    FORMAT(1X//' EIGENVALUE NUMBER ', I5, ' = ', 1P1D20.7)
!
!       WRITE(NOUT, 400) Z
!400    FORMAT(/1X, 14HEIGENVECTOR =  / (1P3D20.7))
!
30    CONTINUE
!
!
9999  CONTINUE
CLOSE(KFILE, STATUS='DELETE')
CLOSE(DAFILE, STATUS='DELETE')
RETURN
END
