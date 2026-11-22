REAL(Kind=WP_Kind) FUNCTION second(dummy)
USE Numeric_Kinds_Module
! PORT NOTE (2026-10-01): original used the VMS/Sun etime extension.
INTEGER :: Dummy
second = 0.d0
RETURN
END FUNCTION second
