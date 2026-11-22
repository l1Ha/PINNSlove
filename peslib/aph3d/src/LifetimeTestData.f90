SUBROUTINE LifetimeTestData(NState, S_Real, S_Imag, SE_Real, SE_Imag)
USE Numeric_Kinds_Module       ! Set numeric kinds
INTEGER, INTENT(IN):: NState
REAL(KIND=WP_Kind)              :: Energy
REAL(KIND=WP_Kind)              :: Theta(2,2), Umat(2,2)
REAL(KIND=WP_Kind), INTENT(OUT) :: S_Real(NState,NState),  S_Imag(NState,NState)
REAL(KIND=WP_Kind), INTENT(OUT) :: SE_Real(NState,NState),  SE_Imag(NState,NState)

Energy=1.d0
Efun=Cosh(Energy)
Theta=5.d0
Theta(1,1)=0.5
Theta(2,2)=0.75
Umat(1,1)= 1.d0/Sqrt(2.d0)
Umat(2,1)= 1.d0/Sqrt(2.d0)
Umat(1,2)= 1.d0/Sqrt(2.d0)
Umat(2,2)=-1.d0/Sqrt(2.d0)
S_Real=MatMul(Matmul(Transpose(Umat),Sin(Theta)),Umat)
S_Imag=MatMul(Matmul(Transpose(Umat),Cos(Theta)),Umat)
SE_Real=Efun*S_Real
SE_Imag=Efun*S_Imag

RETURN   
END SUBROUTINE LifetimeTestData