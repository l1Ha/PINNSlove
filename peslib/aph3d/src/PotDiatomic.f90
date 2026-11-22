SUBROUTINE PotDiatomic(rval, hstep, v0, v1, v2, jvalue, mu)
USE Numeric_Kinds_Module       ! Set numeric kinds
USE Numbers_Module
USE FileUnits_Module           ! Fileunit numbers
USE OneDim_Parms_Module
IMPLICIT NONE
SAVE
INTEGER :: jvalue
REAL(KIND=WP_Kind) r(3), rval, vr, vrmh, vrph
REAL(KIND=WP_Kind) h, hstep
REAL(KIND=WP_Kind) v0, v1, v2, vh, v2h, v3h, v4h
REAL(KIND=WP_Kind) mu, cent0, cent1, cent2

mu=reducedmass
h=1.0d-3
r(1)=2.154434690031884d12  ! PORT FIX: was Ten**(RANGE(One)/3)=1e102 in double precision -> PES overflow (NaN)
r(2)=2.154434690031884d12  ! PORT FIX: was Ten**(RANGE(One)/3)=1e102 in double precision -> PES overflow (NaN)
r(3)=2.154434690031884d12  ! PORT FIX: was Ten**(RANGE(One)/3)=1e102 in double precision -> PES overflow (NaN)
r(IArran)=rval
IF(rval>=four*h)THEN
   CALL Surface(vr,r)
   r(IArran)=rval-h
   CALL Surface(vrmh,r)
   r(IArran)=rval+h
   CALL Surface(vrph,r)
   v0=vr
   v1=(vrph-vrmh)/h
   v2=(vrph-Two*vr+vrmh)/(Two*h**2)
ELSE
    r(IArran)=h
    call Surface(vh,r)
    r(IArran)=h*two
    call Surface(v2h,r)
    r(IArran)=h*three
    call Surface(v3h,r)
    r(IArran)=h*four
    call Surface(v4h,r)
    V0=0
    v0=V0+(rval-h*two)*(rval-h*three)*(rval-h*four)*vh/(-6*h**3)
    v0=v0+(rval-h)*(rval-h*three)*(rval-h*four)*v2h/(2*h**3)
    v0=v0+(rval-h)*(rval-h*two)*(rval-h*four)*v3h/(-2*h**3)
    v0=v0+(rval-h)*(rval-h*two)*(rval-h*three)*v4h/(6*h**3)

    v1=0
    v1=v1+(rval-h*three)*(rval-h*four)*vh/(-6*h**3)
    v1=v1+(rval-h*three)*(rval-h*four)*v2h/(2*h**3)
    v1=v1+(rval-h*two)*(rval-h*four)*v3h/(-2*h**3)
    v1=v1+(rval-h*two)*(rval-h*three)*v4h/(6*h**3)
    v1=v1+(rval-h*two)*(rval-h*four)*vh/(-6*h**3)
    v1=v1+(rval-h)*(rval-h*four)*v2h/(2*h**3)
    v1=v1+(rval-h)*(rval-h*four)*v3h/(-2*h**3)
    v1=v1+(rval-h)*(rval-h*three)*v4h/(6*h**3)
    v1=v1+(rval-h*two)*(rval-h*three)*vh/(-6*h**3)
    v1=v1+(rval-h)*(rval-h*three)*v2h/(2*h**3)
    v1=v1+(rval-h)*(rval-h*two)*v3h/(-2*h**3)
    v1=v1+(rval-h)*(rval-h*two)*v4h/(6*h**3)
    
    V2=0
    v2=v2+(rval-h*two)*vh/(-6*h**3)
    v2=v2+(rval-h)*v2h/(2*h**3)
    v2=v2+(rval-h)*v3h/(-2*h**3)
    v2=v2+(rval-h)*v4h/(6*h**3)
    v2=v2+(rval-h*three)*vh/(-6*h**3)
    v2=v2+(rval-h*three)*v2h/(2*h**3)
    v2=v2+(rval-h*two)*v3h/(-2*h**3)
    v2=v2+(rval-h*two)*v4h/(6*h**3)
    v2=v2+(rval-h*four)*vh/(-6*h**3)
    v2=v2+(rval-h*four)*v2h/(2*h**3)
    v2=v2+(rval-h*four)*v3h/(-2*h**3)
    v2=v2+(rval-h*three)*v4h/(6*h**3)
    v2=2*v2
ENDIF

CALL CenTerm(rval, jvalue, mu, cent0, cent1, cent2, hstep)
v0=v0+cent0
v1=v1+cent1
v2=v2+cent2

RETURN
ENDSUBROUTINE PotDiatomic