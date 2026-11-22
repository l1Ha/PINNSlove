SUBROUTINE Manolopous(u, z, w, n, rstart, rend, nstpl, zed, zp,      &
                      usys2, w2, zedz, EDeriv, zeta, omega, new,   &
                      g0, g1, g2, dmu, psi, gamma,                  &
                      rhoval, xksq, energy, w0, w1, w2b,            &
                      diag, cor, asym, ndim, ipiv, ui, work)

!            P U R P O S E   O F    S U B R O U T I N E
!  This routine propagates a set of coupled second-order differential
!  equations using Bob Johnson's Log Derivative method.

!            I N P U T    A R G U M E N T S
!  u       An temporary array for storing the interaction matrix
!  z       On entering contains the Log Derivative matrix at rstart.
!  n       Number of coupled channels.
!  rstart  Starting distance.
!  rend    Ending distance.
!  nstpl   Number of steps per asymptotic wavelength.
!
!  zeta    gamma vector at rstart
!  omega   temporary vector used in gamma propagation
!  new     temporary vector used in gamma propagation
!

!            O U T P U T    A R G U M E N T S
! z        On return z contains the Wiger R-matrix at rend.
!
! zeta     gamma vector at rend
!

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!                               NOTE

!                  Derivatives Not Working Presently
!
!trp>>
!  Modified 9/12/97:  Changed calls to dsysv to correct LDA and LDB.
!    9/16/97  TRP  Modified calls to addhone() in accordance with changes
!                  in that routine.
!                  Added ipiv, ui, work to argument list.  Changed their
!                  dimensions to n and nXn.  This also meant that two
!                  of the calls to dsysv which were modified on 9/12
!                  had to be set back to their original state.
!<<trp
!
USE Numbers_Module
USE dip_MODULE
USE engpro_MODULE
USE FileUnits_MODULE
USE TotalEng_Module
IMPLICIT NONE
!     ----------------------
INTEGER, PARAMETER:: lwork = 300
!     ----------------------
LOGICAL EDeriv
INTEGER n, nstpl, nstep, istep, nread, ncall, mcall, ndim
INTEGER ilo, ii, il , ir, ic, i, id, jj, info
REAL(Kind=WP_Kind) drad, f3, f6, fh, sq, ag, tn, th, sca, xktot
REAL(Kind=WP_Kind) rstart, rend, r, d1, d2, d3, d4, drnow, wavel, h
REAL(Kind=WP_Kind) usys2, derinit, hstep
REAL(Kind=WP_Kind) dae(lwork,lwork), dal(lwork,lwork), ui(n,n)
REAL(Kind=WP_Kind) us(lwork,lwork),ql(lwork,lwork), qm(lwork,lwork)
REAL(Kind=WP_Kind) qr(lwork,lwork)
REAL(Kind=WP_Kind) z1(lwork), z2(lwork), z1d(lwork), z2d(lwork)
REAL(Kind=WP_Kind) z4(lwork), z3(lwork)
!            D I M E N S I O N S
REAL(Kind=WP_Kind) u(1), z(n,n), w(n,n), zed(n,n), zp(n,n), w2(n,n)
REAL(Kind=WP_Kind) zedz(n,n)
REAL(Kind=WP_Kind) g0(n,n), g1(n,n), g2(n,n), dmu(n,n), psi(n), gamma(n)
REAL(Kind=WP_Kind) rhoval(*), xksq(n), energy(n), w0(n,n), w1(n,n), w2b(n,n)
REAL(Kind=WP_Kind) zeta(n), omega(n), new(n), diag(n,n), cor(n,n), asym(n,n)
REAL(Kind=WP_Kind) uu(lwork,lwork)
REAL(Kind=WP_Kind) work(n), wr(lwork),wl(lwork)
INTEGER ipiv(n)
SAVE
!
!     ---------------------
COMMON /utotal/ xktot
!     ---------------------
!
DATA derinit / 1.d0 /


!            E X T E R N A L S
EXTERNAL pot_shift, smxinv
!
EXTERNAL phot_read
!

!
! Read in energy independent.
!
    drad = 180d0 / pi
!
!         READ(frst_unit)wavel,nstep,drnow,h,d1,d2,d4

! Calculate energy independent.
    wavel=2.d0*pi/sqrt(usys2*engmax)
    IF(wavel>2.d0) wavel=2.d0
    nstep=nstpl*(rend-rstart)/wavel+.9d0
    IF(nstep==0) nstep=1
    drnow=(rend-rstart)/nstep
    h=drnow*0.5d0
    d1=h*h/3.d0
    d2=2.d0*d1
    d4=-d1/16.d0
    d3=2.d0*d2


!
!        ================================
!        Define Constants For Quadratures
!        ================================
!
    f3 = h*h / 3d0
    f6 = h*h / 6d0
    fh = 4d0
!

!
!
!     =======================
!     Calculate -W at Point A
!     =======================
!
!     ==================================================================
!     Potmatn + Pot_shift  Return The Negative Of The Interaction Matrix
!     ==================================================================
!
!
!     ------------------------------------------------------------
r=rstart
!      READ(frst_unit)nread,r,hstep
CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval,     &
                xksq, energy, w0, w1, w2b, zp, diag, cor, asym, ndim)
CALL pot_shift (n, r, u)
!     ------------------------------------------------------------
!
!
!     ======================================================
!     Define The Reference Potential as [ E  -  eigenvalue ]
!     ======================================================
!
!
!     =====================================
!     Reference Potential For Sector = xksq
!     =====================================
!
!     ====================================================
!     xksq =  [ 2 * mu [ E - eigenvalue At Sector Center ]
!     ====================================================
!
!     ============================================
!     Determine Homogeneous Half-Sector Propagator
!     ============================================
!
!     ====================================================================
!     This Portion For Defining The Trigonometric And Hyperbolic Functions
!     That Define z1 Through z4 Taken From M. Alexander.
!     ====================================================================
!
    ilo = ilo + 1
!
!
!     ==========================
!     Calculate h * Y(0) = Z (0)
!     ==========================
!
DO ii = 1,n
!        -----------------------------
    CALL dscal (n, h, z(1,ii), 1)
!        -----------------------------
ENDDO
!
!     ==================================================================
!     Fill In The Two Dimensional Array UU From U And Subtract Reference
!     ==================================================================
!
!     ==================================================================
!     UU Is The Negative Of The Interaction Matrix As Defined By Johnson
!     ==================================================================
!
!     ==================================================
!     UU  =  [2 * mu ] [ E  -  eigenvalue  - V (resid) ]
!     ==================================================
!
!     ------------------------------------------------------------
il = 0
DO ir = 1,n
    DO ic = 1,ir
    il = il + 1
    uu(ir,ic) = u(il)
    uu(ic,ir) = u(il)
    ENDDO
ENDDO
!     ------------------------------------------------------------
!
!     ======================================================
!     U Is The Manalopoulos U  Which In This Code Is Also UU
!     ======================================================
!
!     ==================================================
!     determine z(0)  = h * y (0) + (h * h / 3) *  U (a)
!     ==================================================
!
DO ii = 1,n
!     --------------------------------
CALL dscal ( n, f3, uu(1,ii), 1)
!     --------------------------------
ENDDO
!
!     ================================================
!     Copy [ (h * h ) / 3 ] * U(Left Of Sector) To QL
!     ================================================
!
!
DO ii = 1,n
!     -----------------------------------------
CALL dcopy ( n, uu(1,ii), 1, ql(1,ii), 1)
!     -----------------------------------------
ENDDO
!
!
!
!
!      =================
!      Start Integration
!      =================
!
!
!     modification
!
DO ii = 1,n
    wl(ii) = 0d0
ENDDO
!
DO 100 istep = 1,nstep
!
!
r = rstart + h * (2 * istep - 1)
!
!
!     ======================================================
!     If Only Two Subsectors Here Then We Are At Midsector c
!     ======================================================
!
!     ======================================================
!     Need To Get Interaction Matrix At c To Determine The
!     Quadrature Contribution At Midsector Q(c) Which enters
!     In Forming Script Z4 Which In Its Turn Is Needed for
!     A --------> c Propagation
!     ======================================================
!
!      READ(frst_unit)nread,r,hstep
!     ------------------------------------------------------------
CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval,   &
                xksq, energy, w0, w1, w2b , zp, diag, cor, asym, ndim)
CALL pot_shift (n, r, u)
!     ------------------------------------------------------------
!
!
!     DO ii = 1,n
!        wr(ii) = -xksq(ii)
!     ENDDO
!
!
!
!
    do  i = 1,n
    id = (i * (i + 1)) / 2
    wr(i) = u(id)
    ENDDO
!
!
!

DO i = 1,n
    sq = sqrt(abs(wr(i)))
    ag = (h/2) * sq
    IF(wr(i)<0)THEN
        tn = tan(ag)
        z1(i) =   (ag / tn - ag*tn)
        z2(i) =   (ag / tn + ag*tn)
    ELSE
        th = tanh(ag)
        z1(i) =   (ag / th + ag * th)
        z2(i) =   (ag / th - ag * th)
    ENDIF
ENDDO
!
!     modification
!
DO ii = 1,n
    ql(ii,ii) = ql(ii,ii) - f3 * ( wr(ii) - wl(ii))
ENDDO
!
!     =================================
!     Fill In UU And Subtract Reference
!     =================================
!
!     ----------------------------------------
il = 0
DO ii = 1,n
    DO jj = 1,ii
        il = il + 1
        uu(jj,ii) = u(il)
        uu(ii,jj) = u(il)
    ENDDO
    uu(ii,ii) = uu(ii,ii) - wr(ii)
ENDDO
!     ---------------------------------------
!
!     =======================================
!     Form [ I - ( ( h * h ) / 6 ) * UU (C) ]
!     =======================================
!
DO ii = 1,n
!     --------------------------------
CALL dscal (n, -f6, uu(1,ii), 1)
!     --------------------------------
ENDDO
!
!     -------------------------
CALL addhone (n, uu, lwork, One)
!     -------------------------
!
!     UU = [ I - ( ( h * h ) / 6 ) * UU (c) ]
!     =======================================
!
!
!     =============================================
!     To Invert UU Form The Identity Matrix UI = I
!     Get The Inverse By solving        UU * X = UI
!     =============================================
!
!     -----------------
CALL unit (n, ui)
!     -----------------
!
!     DO 30 ii = 2,n
!        im1 = ii - 1
!        DO 30 jj = 1,im1
!        uu(jj,ii) = uu(ii,jj)
!30    CONTINUE
!     ------------------------------------------------------------
!     CALL dsysv ('l', n, n, uu, n, ipiv, ui, n, work, lwork, info)
CALL dsysv ('l', n, n, uu, lwork, ipiv,  ui, n, work, lwork, info)
!     ------------------------------------------------------------
!
!
!     UI = [ I - ( ( h * h ) / 6 ) * UU (C) ]^{-1}
!     ============================================
!
!
!     ================================================
!     Form [ I - ( ( h * h ) / 6 ) * UU (c) ]^{-1} - I
!     ================================================
!
!     --------------------------
CALL addhone (n, ui, n, -One)
!     --------------------------
!
!     =================================================
!     Determine Quadrature Contribution UI At Midsector
!     =================================================
!
DO ii = 1,n
!     -------------------------------
CALL dscal (n, fh, ui(1,ii), 1)
!     -------------------------------
ENDDO
!
!     UI = 4 * { [ I - ( ( h * h ) / 6 ) * UU (c) ]^{-1} - I }
!     ========================================================
!
!
!     =============
!     Copy UI TO QM
!     =============
!
DO ii = 1,n
!     ----------------------------------------
CALL dcopy (n, ui(1,ii), 1, qm(1,ii), 1)
!     ----------------------------------------
ENDDO
!
!     UI  = Q(c) = (4 / h) [ I - h * h / 6 UU(c) ] ^{-1} - (4 / h) * I
!     ================================================================
!
!
!     =============================
!     Propagate Z To C (Mid-Sector)
!     =============================
!
!     ================================================================
!     Script Z1 = Q(A) + z1(diagonal)
!     where Q(a) = ql = Quadrature Contribution At Left Side Of Sector
!     ================================================================
!
DO ii = 1,n
    DO jj = 1,n
    z(jj,ii) = z(jj,ii) + ql(jj,ii)
    ENDDO
    z(ii,ii) = z(ii,ii) + z1(ii)
ENDDO
!
!     Z = Z(a) + Z1(a,c)
!     ==================
!
!     ==================================
!     Determine [ Z(a) + Z1(a,c) ] ^{-1}
!     ==================================
!
!     -----------------
CALL unit (n, ui)
!     -----------------
!
!     DO 10 ii = 2,n
!        im1 = ii - 1
!        DO 10 jj = 1,im1
!        z(jj,ii) = z(ii,jj)
!10    CONTINUE
!
!
!
!     ------------------------------------------------------------
CALL dsysv ('l', n, n, z, n, ipiv, ui, n, work, lwork, info)
!     ------------------------------------------------------------
!
!
!
!     UI = [ Z(a) + Z1(1,c) ] ^{-1}
!     =============================
!
!     =======================
!     Save UI For Derivatives
!     =======================
!
!     DO ii = 1,n
!     ----------------------------------------
!     CALL dcopy (n, ui(1,ii), 1, us(1,ii), 1)
!     ----------------------------------------
!     ENDDO
!
!     US = [ Z(a) + Z1(1,c) ] ^{-1}
!     =============================
!
!
!     =============================================
!     Form UI = [ Z(a) + Z1(1,c) ] ^{-1} * Z2 (a,c)
!     =============================================
!
!
DO jj = 1,n
sca = z2(jj)
!     --------------------------------
CALL dscal (n, sca, ui(1,jj), 1)
!     --------------------------------
ENDDO

!
!
!     UI = [ Z(a) + Z1(1,c) ] ^{-1} * Z2 (a,c)
!     ========================================
!
!     ====================================================================
!     Calculate -Z3(a,c) * [ Z(a) + Z1(a,c) ] ^{-1} * Z2 (a,c) ] Transpose
!     ====================================================================
!
DO jj = 1,n
    DO ii = 1,n
    sca = -z2(ii)
    ui(ii,jj) = sca * ui(ii,jj)
    ENDDO
ENDDO
!
!
!    ==========================================================
!    Calculate Z4 - Z3E * [ Z(a) + Z1(1,c) ] ^{-1} * Z2 (a,c) ]
!    ==========================================================
!
!     =====================================================================
!     U(ii,jj) Is Added To Z(jj,ii) Since UI Transpose Was Calculated Above
!     =====================================================================
!
DO ii = 1,n
DO jj = 1,n
    z(jj,ii) = ui(jj,ii) + qm(jj,ii)
ENDDO
    z(ii,ii) = z(ii,ii) + z1(ii)
ENDDO
!
!
!
!    Z = Z(a,c)
!    ==========
!
!
!
!
!
!     IF(EDeriv)THEN
!     kkkkkkkkkkkkkkkkk
!     ||||||||||||||||||||||||||||||||||||||
!     Begin Derivative Calculation dZ_n / dE
!     ||||||||||||||||||||||||||||||||||||||
!
!     ||||||||||||||||||||||||||||||||||||
!     End Derivative Calculation dZ_n / dE
!     ||||||||||||||||||||||||||||||||||||
!
!     ENDIF
!     kkkkk
!
!
!     ==============================
!     Propagate Z Matrix From C to B
!     ==============================
!
!
DO ii = 1,n
!     ----------------------------------------
CALL dcopy (n, qm(1,ii), 1, ql(1,ii), 1)
!     ----------------------------------------
ENDDO
!
r = rstart + h*(2*istep)
!      READ(frst_unit)nread,r,hstep
!     ------------------------------------------------------------
CALL potmatn ( n, r, u, g0, g1, g2, dmu, psi, gamma, rhoval,  &
              xksq, energy, w0, w1, w2b , zp, diag, cor, asym, ndim)
!
CALL pot_shift (n, r, u)
!     ------------------------------------------------------------
il = 0
DO ii = 1,n
    DO jj = 1,ii
        il = il + 1
        uu(jj,ii) = u(il)
        uu(ii,jj) = u(il)
    ENDDO
    uu(ii,ii) = uu(ii,ii) - wr(ii)
ENDDO
!     -------------------------------------------------------------
!
!     ====================================
!     Get Quadrature Contribution At Right
!     ====================================
DO ii = 1,n
!     --------------------------------
CALL dscal ( n, f3, uu(1,ii), 1)
!     --------------------------------
ENDDO
!
!     modification
!
DO ii = 1,n
    wl(ii) = wr(ii)
ENDDO
!
!     modification
!
DO ii = 1,n
!     -----------------------------------------
CALL dcopy ( n, uu(1,ii), 1, qr(1,ii), 1)
!     -----------------------------------------
ENDDO
!
DO ii = 1,n
    DO jj = 1,n
    z(jj,ii) = z(jj,ii) + ql(jj,ii)
ENDDO
    z(ii,ii) = z(ii,ii) + z1(ii)
ENDDO
!
DO ii = 1,n
!     -----------------------------------------
CALL dcopy ( n, qr(1,ii), 1, ql(1,ii), 1)
!     -----------------------------------------
ENDDO
!
!     Z  =  [ Z(c) + Z1(c,b) ]
!     ========================
!
!     -----------------
CALL unit (n, ui)
!     -----------------
!
!     DO 20 ii = 2,n
!        im1 = ii - 1
!        DO 20 jj = 1,im1
!        z(jj,ii) = z(ii,jj)
!20    CONTINUE
!
!     ------------------------------------------------------------
CALL dsysv ('l', n, n, z, n, ipiv, ui, n, work, lwork, info)
!     ------------------------------------------------------------
!
!     UI = [ Z(c) + Z1(c,b) ] ^{-1}
!     =============================
!
!
!     =============================
!     Save UI In US For Derivatives
!     =============================
!
!     DO ii = 1,n
!     ----------------------------------------
!     CALL dcopy (n, ui(1,ii), 1, us(1,ii), 1)
!     ----------------------------------------
!     ENDDO
!
!     US = [ Z(a) + Z1(1,c) ] ^{-1}
!     =============================
!
DO jj = 1,n
sca = z2(jj)
!     --------------------------------
CALL dscal (n, sca, ui(1,jj), 1)
!     --------------------------------
ENDDO
!
!     UI = [ Z(c) + Z1(c,b) ] ^{-1} * Z2 (c,b)
!     ========================================
!
!
!
!     IF(EDeriv)THEN
!     bbbbbbbbbbbbbbbbb
!     ||||||||||||||||||||||||||||||||||||||
!     Begin Derivative Calculation dZ_n / dE
!     ||||||||||||||||||||||||||||||||||||||
!
!     |||||||||||||||||||||||||||
!     End Derivative Calculations
!     |||||||||||||||||||||||||||
!     ENDIF
!     bbbbb
!
DO jj = 1,n
    DO ii = 1,n
    sca = -z2(ii)
    ui(ii,jj) = sca * ui(ii,jj)
    ENDDO
ENDDO
!
!
!
!     UI =  - Z3 (a,c) * [ Z(a) + Z1(1,c) ] ^{-1} * Z2 (a,c)
!     ======================================================
!
!
DO ii = 1,n
    DO jj = 1,n
    z(jj,ii) =    ui(jj,ii) + qr(jj,ii)
    ENDDO
    z(ii,ii)    = z(ii,ii) + z1(ii)
ENDDO
!
!
!     DO ii = 1,n
!     ----------------------------------------
!     CALL dcopy (n, qr(1,ii), 1, ql(1,ii), 1)
!     ----------------------------------------
!     ENDDO
!
!     Z  =  Z4(c,b)  - Z3 (c,b) * [ Z(c) + Z1(c,b) ] ^{-1} * Z2 (c,b)
!     ===============================================================
!
!
!
100 CONTINUE
!
!
!     =========================
!     Reads To Stay In Register
!     =========================
!
!
!     =====================================
!     Determine y(2*nstep)=[z(2*n)+I]/h-I/h
!     =====================================
!
z=z/h
!     CALL addhone (n, z, n, -1.d0/h)
!
!
! detrmine gamma_2 at end of sector
!
IF(lphoto)THEN
zeta=zeta/h
ENDIF
!

IF(EDeriv)THEN
   zed=zed/h
ENDIF
!
!
205   FORMAT(5(d14.6,1x))
!***********************************************************************
!     now z = Y, zed = YE
!***********************************************************************
RETURN
ENDSUBROUTINE Manolopous
