SUBROUTINE Morse_Pot(r, n, mu, v, vp, vpp)
USE Numeric_Kinds_Module
USE Numbers_Module
USE FileUnits_Module
USE FileUnits_OneDim_Module
USE CriticalDist_Module
USE DiatomicPot_MODULE
USE Convrsns_Module
IMPLICIT NONE
INTEGER :: n, i, nu, numax

REAL(KIND=WP_Kind) :: mu, w0, wnu, Eigens, arg, a

REAL(KIND=WP_Kind) :: r(0:n)
REAL(KIND=WP_Kind) :: v(0:n)
REAL(KIND=WP_Kind) :: vp(0:n)
REAL(KIND=WP_Kind) :: vpp(0:n)

w0=Sqrt(Two*k/mu)

a=Beta
WRITE(Out_Unit,*)
WRITE(Out_Unit,*)"Morse Parameters"
WRITE(Out_Unit,*)"  De=",De
WRITE(Out_Unit,*)"  re=",re
WRITE(Out_Unit,*)"  k=",k
WRITE(Out_Unit,*)"  Beta=",Beta
WRITE(Out_Unit,*)"  w0=",w0

OPEN(unit=Morse_Unit,file=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"__MorseEigen.txt")
WRITE(Morse_Unit,'(A5,2A15)') "nu", "Eigens(Ha)","Eigens(cm^-1)"
numax=Two*De/w0-One/Two
DO nu=0,numax
    wnu=w0*(One*nu+Half)
    Eigens=wnu*(One-wnu/(Four*De))+vasy-De
    WRITE(Morse_Unit,'(I5,4ES15.7)') nu, Eigens, Eigens/cmm1toau
ENDDO
CLOSE(unit=Morse_Unit)

OPEN(Unit=Morse_Unit,File=TRIM(OutDIR)//"DiatomicOut/"//TRIM(DiatomicPot)//"__MorsePot.txt")
WRITE(Morse_Unit,'(4A15)')"r(i)", "v(i)", "vp(i)", "vpp(i)"
DO i=0,n
    arg=Beta*(re-r(i))
    v(i)=De*(One-exp(arg))**2-De+vasy
    vp(i)=Two*Beta*De*exp(arg)*(One-exp(arg))
    vpp(i)=Two*Beta**2*De*exp(arg)*(Two*exp(arg)-One)
    WRITE(Morse_Unit,'(4ES15.7)')r(i),v(i),vp(i),vpp(i)
ENDDO
CLOSE(unit=Morse_Unit)

RETURN
END SUBROUTINE Morse_Pot

