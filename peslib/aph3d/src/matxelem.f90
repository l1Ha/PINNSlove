SUBROUTINE matxelem (nelem, naph1, phi1, th, ch, naph2, phi2, numnp, id, fbar, fbarp, thetaval, chivals,         &
                           type, ntheta, nchi, hamil, rho, rhocent,   xpt, wht, value, nquad, h, p, hcoef)
USE Numeric_Kinds_Module
USE Masses_Module
USE csbasis_Module
USE totj_Module
Use FileUnits_Module
!
!   $RCSfile: matxelem.f,v $ $Revision: 1.15 $
!   $Date: 89/11/16 11:03:06 $
!   $State: Stable $
!

!              P U R P O S E   O F   S U B R O U T I N E

!        Calculates matrix elements  < Phi|V|Phi  >

!                  I N P U T   A R G U M E N T S

!    nelem      Number of elements.
!    naph1      Number of surface functions in array phi1.
!    phi1       Surface functions at rho 1.
!    th         List of theta indices.
!    thetaval   List of theta values.
!    ntheta     Number of theta angles.
!    ch         List of chi indicies.
!    chivals    List of chi values.
!    nchi       Number of chi angles.
!    naph2      Number of surface functions in array phi2.
!    phi2       Surface functions at rho 2. [Note: rho 1 is equal to
!               rho 2 for matrix elements and they are not equal for
!               overlap matrix elements.
!    numnp      Number of nodal points.
!    id         List of node numbers.
!    type       Used to specify which type of integral to evaluate.
!               type = 1  Overlap matrix elements.
!               type = 2  Potential matrix elements.
!               type = 3  Coriolis coupling matrix elements.
!               type = 4  Asymetric Top coupling matrix elements.
!    fbar       Temporary array used to evaluate the phi1 surface
!               functions at the Gauss_Legendre quadrature points.
!    fbarp      Temporary array used to evaluate the phi2 surface
!               functions at the Gauss_Legendre quadrature points.
!    xpt        Abscissas for the Gauss_Legendre integration.
!    wht        Weights for the Gauss_Legendre integration.
!    nquad      Number of Gauss_Legendre quadrature points.
!    ch
!    rho
!    rhocent
!    value
!    h
!    p
!    hcoef


!               O U T P U T   A R G U M E N T S

!    hamil      Resulting matrix elements.

IMPLICIT NONE

!      I N T E G E R S
INTEGER th, ch, numnp, naph1, ntheta, nchi, id, nelem, ielem, type, i, j, k, naph2, imap, nquad

!      R E A L S
REAL(Kind=WP_Kind) phi1, phi2, hamil, thetaval, chivals, h, p, fbar, fbarp, hcoef, value, pot, potcent, rho, rhocent
REAL(Kind=WP_Kind) thetas, chis, stheta, ctheta, theta, chi, tmursi, xtra, xtrafac,  xpt, wht, zero, one, two

!      L O G I C A L S
LOGICAL first

!      D I M E N S I O N S
DIMENSION phi1(numnp, naph1), hamil(naph1, naph2), phi2(numnp, naph2), thetaval(ntheta), chivals(nchi)
DIMENSION th(numnp), ch(numnp), id(9,nelem), h(9,nquad,nquad), p(2, 9, nquad, nquad), fbar(naph1, nquad*nquad)
DIMENSION fbarp(naph2, nquad*nquad), xpt(nquad), wht(nquad), hcoef(9,nquad, nquad), value(nquad, nquad), thetas(9)
DIMENSION chis(9), imap(9)

!      I N T R I N S I C S  F U N C T I O N S
INTRINSIC sin, cos
EXTERNAL potaph, interpfn, vsets, surfunct, quadsfun, thetachi
DATA imap /2, 6, 3, 5, 9, 7, 1, 8, 4/
DATA first /.true./

!      P A R A M E T E R S
PARAMETER (zero=0.d0, one=1.d0, two=2.d0)
IF(first)THEN
first=.false.
!-----------------------------------------------------------------------
! calculate interpolation functions.
!-----------------------------------------------------------------------
CALL interpfn (h, p, xpt, wht, nquad)
ENDIF

!-----------------------------------------------------------------------
! Zero the matrix.
!-----------------------------------------------------------------------
CALL vsets(naph1*naph2, hamil, 1, zero)

IF(type>4) RETURN
!----------------------------------------------------------------------
! calculate matrix elements.
!----------------------------------------------------------------------
tmursi=1./(usys2*rho**2)
xtrafac=-0.25*tmursi*(jtot*(jtot+1))*(-1)**(jtot+parity)
!  start loop on the elements

DO 106 ielem = 1, nelem/2

!----------------------------------------------------------------------
! Determine theta and chi values at each of the nodal points.
!----------------------------------------------------------------------
CALL thetachi (th, ch, thetaval, chivals, ntheta, nchi, numnp, id(1,ielem), thetas, chis)

!----------------------------------------------------------------------
! Overlap matrix elements have the integrand equal to 1.
!  k runs over the nodal points of the element for interpolation.
!  j and i run over the quadrature points of the element.
!----------------------------------------------------------------------
   IF(type==1)THEN
      DO 105 j = 1,nquad
         DO 102 i = 1,nquad
            value(i,j) = one
            DO 101 k = 1,9
               hcoef(k, i, j) = h(k, i, j)
  101             CONTINUE
  102          CONTINUE
  105       CONTINUE

!----------------------------------------------------------------------
! Potential matrix coupling elements.
!----------------------------------------------------------------------
   ELSEIF(type==2)THEN
      DO 205 j = 1,nquad
         DO 202 i = 1,nquad
            chi = zero
            theta = zero
            DO 201 k = 1,9
               hcoef(k, i, j) = h(k, i, j)
               theta = theta + thetas(k)*h(imap(k), i, j)
               chi = chi + chis(k)*h(imap(k), i, j)
  201             CONTINUE
            CALL potaph (rhocent, theta, chi, potcent)
            CALL potaph (rho    , theta, chi, pot)
            value(i,j) =  pot - (rhocent/rho)**2*potcent
            IF(mega==1.AND..NOT.cstest)THEN
            stheta=sin(theta)
            xtra=xtrafac*(two/(1.+stheta)-one/stheta**2)
!  xtra is -(-1)**(j+p)*j(j+1)*(a-b)/4
            value(i,j)=value(i,j)+xtra
            ENDIF
  202          CONTINUE
  205       CONTINUE

!----------------------------------------------------------------------
! Coriolis coupling matrix elements.
!----------------------------------------------------------------------
   ELSEIF(type==3)THEN
      DO 305 j = 1,nquad
         DO 302 i = 1,nquad
            chi = zero
            theta = zero
            DO 301 k = 1,9
               hcoef(k, i, j) = p(2, k, i, j)
               theta = theta + thetas(k)*h(imap(k), i, j)
               chi = chi + chis(k)*h(imap(k), i, j)
  301             CONTINUE
            stheta = sin(theta)
            ctheta = cos(theta)
            value(i,j) = -tmursi*ctheta/stheta**2
  302          CONTINUE
  305       CONTINUE

!----------------------------------------------------------------------
! Asymetric top coupling matrix elements.
!----------------------------------------------------------------------
   ELSEIF(type==4)THEN
      DO 405 j = 1,nquad
         DO 402 i = 1,nquad
            chi = zero
            theta = zero
            DO 401 k = 1,9
               hcoef(k, i, j) = h( k, i, j)
               theta = theta + thetas(k)*h(imap(k), i, j)
               chi = chi + chis(k)*h(imap(k), i, j)
  401             CONTINUE
            stheta = sin(theta)
            value(i,j) = tmursi*(2./(1.+stheta)-1./stheta**2)
!  value here is a-b.
  402          CONTINUE
  405       CONTINUE

!----------------------------------------------------------------------
! Illegal type.
!----------------------------------------------------------------------
   ELSE
      WRITE(Out_unit,*)'matxelem: Wrong type...error , type = ', type
      STOP 'matxelem'
   ENDIF

!----------------------------------------------------------------------
! evaluate the functions at the Gaussian integration points using
! the interpolation coefficients that were calculated earlier.
!----------------------------------------------------------------------
   CALL surfunct (naph1,  fbar, h    , id(1,ielem), phi1, numnp, nquad)
   CALL surfunct (naph2, fbarp, hcoef, id(1,ielem), phi2, numnp, nquad)

!----------------------------------------------------------------------
!  Evaluate the matrix element contribution for this element.
!----------------------------------------------------------------------
   CALL quadsfun (thetas, chis, hamil, fbar, fbarp, value, naph1, naph2, h, p, xpt, wht, nquad)
  106 CONTINUE
RETURN
!
END
