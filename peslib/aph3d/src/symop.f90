SUBROUTINE symop (irrep, j, i, group, gamma, n, isymop, n1, nirrep, nsymop, kase, system, fac, par)
USE FileUnits_Module
IMPLICIT NONE
CHARACTER(LEN=3) :: group, system
INTEGER :: n, irrep, j, i, gamma, gammaC6v (12, 8), isymop, &
 gammaC1 (1, 1), gammaC2 (2, 2), gammaC1v (2, 2), gammaC2v (4, 4), &
 gammaC3v (6, 4), gammaC3 (3, 3), gammaC6 (6, 6), n1 (8), nirrep, &
 nsymop, kase, fac, par (8)
DATA gammaC1 / 1 /
DATA gammaC2 / 1, 1, 1,-1 /
DATA gammaC1v / 1, 1, 1,-1 /
DATA gammaC2v / 1, 1, 1, 1, 1, 1,-1,-1, 1,-1, 1,-1, 1,-1,-1, 1 /
DATA gammaC3 / 1, 1, 1, 2,-1,-1, 2,-1,-1 /
DATA gammaC3v / 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1, 2, &
-1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0 /
DATA gammaC6 / 1, 1, 1, 1, 1, 1, 1,-1, 1,-1, 1,-1, 2, 1, &
-1,-2,-1, 1, 0, 1, 1, 0,-1,-1, 2,-1,-1, 2,-1, &
-1, 0, 1,-1, 0, 1,-1 /
DATA gammaC6v / 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
 1,-1,-1,-1,-1,-1,-1, 1,-1,-1, 1, 1,-1, 1, 1, 1, &
-1,-1,-1, 1,-1,-1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 2, &
 1, 1,-1,-1,-2, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, &
 0, 0, 0, 0, 2,-1,-1,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 1, &
-1,-1, 1, 0, 0, 0, 0, 0, 0, 0 /
fac = 0
IF(group=='C1 ')THEN
   IF(n<0)THEN
      WRITE(Out_Unit, * ) 'Error: n=', n
      STOP 'symop'
   ENDIF
   IF(kase==0)THEN
      WRITE(Out_Unit,  * ) 'Warning NOT exploiting full symmetry for an ABC system'
      nirrep = 1
      nsymop = 1
      n1 (1) = n
      par (1) = 0
      IF(system=='ABC'.or.system=='ABB'.or.system=='AAA')THEN
         WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for system=', system
      ENDIF
      RETURN
   ENDIF
   IF(irrep<1.or.irrep>1)THEN
      WRITE(Out_Unit, * ) 'irrep=', irrep
      STOP 'symop'
   ENDIF
   IF(isymop<1.or.isymop>1)THEN
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
   gamma = gammaC1 (isymop, irrep)
   IF(isymop==1)THEN
      i = j
   ELSE
      WRITE(Out_Unit, * ) 'Error: isymop=', isymop
   ENDIF
ELSEIF (group=='C2 ')THEN
   IF(n<0.or. (n / 2) * 2/=n)THEN
      WRITE(Out_Unit, * ) 'Error: n=', n
      STOP 'symop'
   ENDIF
   IF(kase==0)THEN
      WRITE(Out_Unit, * ) 'Exploiting full symmetry for an ABC system'
      nirrep = 2
      nsymop = 2
      n1 (1) = n / 2
      n1 (2) = n / 2
      par (1) = 0
      par (2) = 1
      IF(system=='ABC'.or.system=='ABB'.or.system=='AAA')THEN
         WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for system=', system
      ENDIF
      RETURN
   ENDIF
   IF(irrep<1.or.irrep>2)THEN
      WRITE(Out_Unit, * ) 'irrep=', irrep
      STOP 'symop'
   ENDIF
   IF(isymop<1.or.isymop>2)THEN
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
   gamma = gammaC2 (isymop, irrep)
   IF(isymop==1)THEN
      i = j
   ELSEIF (isymop==2)THEN
      i = mod (j-1 + 6 * n / 12, n) + 1
   ELSE
      WRITE(Out_Unit, * ) 'Error: isymop=', isymop
   ENDIF
ELSEIF (group=='C1v')THEN
   IF(system=='ABC')THEN
      WRITE(Out_Unit, * ) 'system=', system
      STOP 'symop'
   ENDIF
   IF(n<0.or. (n / 2) * 2/=n)THEN
      WRITE(Out_Unit, * ) 'Error: n=', n
      STOP 'symop'
   ENDIF
   IF(kase==0)THEN
      WRITE(Out_Unit,  * ) 'Warning NOT exploiting full symmetry for an AB2 system'
      nirrep = 2
      nsymop = 2
      n1 (1) = n / 2 + 1
      n1 (2) = n / 2-1
      par (1) = 0
      par (1) = 0
      IF(system=='ABB'.or.system=='AAA')THEN
         WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for an system=', system
      ENDIF
      RETURN
   ENDIF
   IF(irrep<1.or.irrep>2)THEN
      WRITE(Out_Unit, * ) 'irrep=', irrep
      STOP 'symop'
   ENDIF
   IF(isymop<1.or.isymop>2)THEN
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
   gamma = gammaC1v (isymop, irrep)
   IF(isymop==1)THEN
      i = j
   ELSEIF (isymop==2)THEN
      i = mod (n + 1-j, n) + 1
   ELSE
      WRITE(Out_Unit, * ) 'Error: isymop=', isymop
   ENDIF
ELSEIF (group=='C2v')THEN
   IF(system=='ABC')THEN
      WRITE(Out_Unit, * ) 'system=', system
      STOP 'symop'
   ENDIF
   IF(n<0.or. (n / 4) * 4/=n)THEN
      WRITE(Out_Unit, * ) 'Error: n=', n
      STOP 'symop'
   ENDIF
   IF(kase==0)THEN
      IF(system=='AB2')THEN
         WRITE(Out_Unit, * ) 'Exploiting full symmetry for an AB2 system'
      ENDIF
      nirrep = 4
      nsymop = 4
      n1 (1) = n / 4 + 1
      n1 (2) = n / 4-1
      n1 (3) = n / 4
      n1 (4) = n / 4
      par (1) = 0
      par (2) = 1
      par (3) = 0
      par (4) = 1
      IF(system=='AAA')THEN
         WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for system=', system
      ENDIF
      RETURN
   ENDIF
   IF(irrep<1.or.irrep>4)THEN
      WRITE(Out_Unit, * ) 'irrep=', irrep
      STOP 'symop'
   ENDIF
   IF(isymop<1.or.isymop>4)THEN
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
   gamma = gammaC2v (isymop, irrep)
   IF(isymop==1)THEN
      i = j
   ELSEIF (isymop==2)THEN
      i = mod (j-1 + 6 * n / 12, n) + 1
   ELSEIF (isymop==3)THEN
      i = mod (n + 1-j, n) + 1
   ELSEIF (isymop==4)THEN
      i = mod (n + 1-j + 6 * n / 12, n) + 1
   ELSE
      WRITE(Out_Unit, * ) 'Error: isymop=', isymop
   ENDIF
ELSEIF (group=='C3 ')THEN
   IF(system=='ABC'.or.system=='ABB')THEN
      WRITE(Out_Unit, * ) 'Error: system=', system
      STOP 'symop'
   ENDIF
   IF(n<0.or. (n / 6) * 6/=n)THEN
      WRITE(Out_Unit, * ) 'Error: n=', n
      STOP 'symop'
   ENDIF
   IF(kase==0)THEN
      WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for an A3 system'
      nirrep = 3
      nsymop = 3
      n1 (1) = n / 3
      n1 (2) = n / 3
      n1 (3) = n / 3
      par (1) = 0
      par (1) = 0
      par (1) = 0
      IF(system=='AAA')THEN
         WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for system=', system
      ENDIF
      RETURN
   ENDIF
   IF(irrep<1.or.irrep>3)THEN
      WRITE(Out_Unit, * ) 'irrep=', irrep
      STOP 'symop'
   ENDIF
   IF(isymop<1.or.isymop>3)THEN
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
   IF(irrep==1)THEN
      fac = 0
   ELSEIF (irrep==2)THEN
      fac = 1
   ELSE
      fac =-1
   ENDIF
   gamma = gammaC3 (isymop, irrep)
   IF(isymop==1)THEN
      i = j
   ELSEIF (isymop==2)THEN
      i = mod (j-1 + 4 * n / 12, n) + 1
   ELSEIF (isymop==3)THEN
      i = mod (j-1 + 8 * n / 12, n) + 1
   ELSE
      WRITE(Out_Unit, * ) 'Error: isymop=', isymop
   ENDIF
ELSEIF (group=='C3v')THEN
   IF(system=='ABC'.or.system=='ABB')THEN
      WRITE(Out_Unit, * ) 'Error: system=', system
      STOP 'symop'
   ENDIF
   IF(n<0.or. (n / 6) * 6/=n)THEN
      WRITE(Out_Unit, * ) 'Error: n=', n
      STOP 'symop'
   ENDIF
   IF(kase==0)THEN
      WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for an A3 system'
      nirrep = 4
      nsymop = 6
      n1 (1) = n / 6 + 1
      n1 (2) = n / 6-1
      n1 (3) = 2 * n / 6
      n1 (4) = 2 * n / 6
      par (1) = 0
      par (2) = 0
      par (3) = 0
      par (4) = 0
      IF(system=='AAA')THEN
         WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for system=', system
      ENDIF
      RETURN
   ENDIF
   IF(irrep<1.or.irrep>4)THEN
      WRITE(Out_Unit, * ) 'irrep=', irrep
      STOP 'symop'
   ENDIF
   IF(isymop<1.or.isymop>6)THEN
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
   IF(irrep<=2)THEN
      fac = 0
   ELSEIF (irrep==3)THEN
      fac = 1
   ELSE
      fac =-1
   ENDIF
   gamma = gammaC3v (isymop, irrep)
   IF(isymop==1)THEN
      i = j
   ELSEIF (isymop==2)THEN
      i = mod (j-1 + 4 * n / 12, n) + 1
   ELSEIF (isymop==3)THEN
      i = mod (j-1 + 8 * n / 12, n) + 1
   ELSEIF (isymop==4)THEN
      i = mod (n + 1-j, n) + 1
   ELSEIF (isymop==5)THEN
      i = mod (n + 1-j + 4 * n / 12, n) + 1
   ELSEIF (isymop==6)THEN
      i = mod (n + 1-j + 8 * n / 12, n) + 1
   ELSE
      WRITE(Out_Unit, * ) 'Error: isymop=', isymop
   ENDIF
ELSEIF (group=='C6 ')THEN
   IF(system=='ABC'.or.system=='ABB')THEN
      WRITE(Out_Unit, * ) 'Error: system=', system
      STOP 'symop'
   ENDIF
   IF(n<0.or. (n / 6) * 6/=n)THEN
      WRITE(Out_Unit, * ) 'Error: n=', n
      STOP 'symop'
   ENDIF
   IF(kase==0)THEN
      WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for an A3 system'
      nirrep = 6
      nsymop = 6
      n1 (1) = n / 6
      n1 (2) = n / 6
      n1 (3) = n / 6
      n1 (4) = n / 6
      n1 (5) = n / 6
      n1 (6) = n / 6
      par (1) = 0
      par (2) = 0
      par (3) = 0
      par (4) = 0
      par (5) = 0
      par (6) = 0
      IF(system=='AAA')THEN
         WRITE(Out_Unit, * ) 'Warning NOT exploiting full symmetry for system=', system
      ENDIF
      RETURN
   ENDIF
   IF(irrep<1.or.irrep>6)THEN
      WRITE(Out_Unit, * ) 'irrep=', irrep
      STOP 'symop'
   ENDIF
   IF(isymop<1.or.isymop>6)THEN
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
   IF(irrep<=2)THEN
      fac = 0
   ELSEIF (irrep==3.or.irrep==5)THEN
      fac = 1
   ELSE
      fac =-1
   ENDIF
   gamma = gammaC6 (isymop, irrep)
   IF(isymop==1)THEN
      i = j
   ELSEIF (isymop==2)THEN
      i = mod (j-1 + 2 * n / 12, n) + 1
   ELSEIF (isymop==3)THEN
      i = mod (j-1 + 4 * n / 12, n) + 1
   ELSEIF (isymop==4)THEN
      i = mod (j-1 + 6 * n / 12, n) + 1
   ELSEIF (isymop==5)THEN
      i = mod (j-1 + 8 * n / 12, n) + 1
   ELSEIF (isymop==6)THEN
      i = mod (j-1 + 10 * n / 12, n) + 1
   ELSE
      WRITE(Out_Unit, * ) 'Error: isymop=', isymop
   ENDIF
ELSEIF (group=='C6v')THEN
   IF(system=='ABC'.or.system=='ABB')THEN
      WRITE(Out_Unit, * ) 'Error: system=', system
      STOP 'symop'
   ENDIF
   IF(n<0.or. (n / 12) * 12/=n)THEN
      WRITE(Out_Unit, * ) 'Error: n=', n
      STOP 'symop'
   ENDIF
   IF(kase==0)THEN
      WRITE(Out_Unit, * ) 'Exploiting full symmetry for A3'
      nirrep = 8
      nsymop = 12
      n1 (1) = n / 12 + 1
      n1 (2) = n / 12-1
      n1 (3) = n / 12
      n1 (4) = n / 12
      n1 (5) = 2 * n / 12
      n1 (6) = 2 * n / 12
      n1 (7) = 2 * n / 12
      n1 (8) = 2 * n / 12
      par (1) = 0
      par (2) = 0
      par (3) = 1
      par (4) = 1
      par (5) = 1
      par (6) = 1
      par (7) = 0
      par (8) = 0
      RETURN
   ENDIF
   IF(irrep<1.or.irrep>8)THEN
      WRITE(Out_Unit, * ) 'irrep=', irrep
      STOP 'symop'
   ENDIF
   IF(isymop<1.or.isymop>12)THEN
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
   IF(irrep<=4)THEN
      fac = 0
   ELSEIF (irrep==5.or.irrep==7)THEN
      fac = 1
   ELSE
      fac =-1
   ENDIF
   gamma = gammaC6v (isymop, irrep)
   IF(isymop==1)THEN
      i = j
   ELSEIF (isymop==2)THEN
      i = mod (j-1 + 2 * n / 12, n) + 1
   ELSEIF (isymop==3)THEN
      i = mod (j-1 + 10 * n / 12, n) + 1
   ELSEIF (isymop==4)THEN
      i = mod (j-1 + 4 * n / 12, n) + 1
   ELSEIF (isymop==5)THEN
      i = mod (j-1 + 8 * n / 12, n) + 1
   ELSEIF (isymop==6)THEN
      i = mod (j-1 + 6 * n / 12, n) + 1
   ELSEIF (isymop==7)THEN
      i = mod (n + 1-j, n) + 1
   ELSEIF (isymop==8)THEN
      i = mod (n + 1-j + 4 * n / 12, n) + 1
   ELSEIF (isymop==9)THEN
      i = mod (n + 1-j + 8 * n / 12, n) + 1
   ELSEIF (isymop==10)THEN
      i = mod (n + 1-j + 2 * n / 12, n) + 1
   ELSEIF (isymop==11)THEN
      i = mod (n + 1-j + 6 * n / 12, n) + 1
   ELSEIF (isymop==12)THEN
      i = mod (n + 1-j + 10 * n / 12, n) + 1
   ELSE
      WRITE(Out_Unit, * ) 'isymop=', isymop
      STOP 'symop'
   ENDIF
ELSE
   WRITE(Out_Unit, * ) 'Incorrect group: group=', group
   STOP 'symop'
ENDIF
IF(i<1.or.i>n)THEN
WRITE(Out_Unit,  * ) 'i=', i, '   n=', n
   STOP 'symop'
ENDIF
RETURN
ENDSUBROUTINE symop
