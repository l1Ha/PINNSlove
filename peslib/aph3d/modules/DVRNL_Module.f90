MODULE DVRNL_Module
   USE Numeric_Kinds_Module
   IMPLICIT NONE
   SAVE
   CHARACTER(LEN=3) :: pointu = 'c2v'       ! point group of upper PES
   CHARACTER(LEN=3) :: pointl = 'c2v'         ! point group of the lower PES
   LOGICAL :: cstestu =.False.               ! Coupled-States on upper surface
   LOGICAL :: cstestl = .False.               ! Coupled-States on lower surface
   LOGICAL :: restart = .False.     ! LOGICAL switch IF true will enable restart option
   LOGICAL :: jevenu  = .True.                ! Rotation even states on upper surface
   INTEGER :: nrhos1u = 1               ! the no. of points where the dvr energy parameters are set.
   INTEGER :: iparu = 0                 ! parity for states on upper PES.
   INTEGER :: jtotu = 0                 ! total J for states on the upper surface.
   INTEGER :: megau = 0                 ! Projection quantum number for states on the upper surface.
   INTEGER :: nbreak=999998        ! after each NBREAK rho WRITE out a separate matrix element, overlap, dipole moment and
                                  ! wave function file, to disable make nbreak very large nbreak must be an even number !!!
   INTEGER :: irhost=3            ! number of rho where you want to restart irhost must be odd !!
   LOGICAL :: jevenl = .True.                 ! Rotation even states on lower surface
   INTEGER :: nrhos1l = 1                ! the no. of points where the dvr energy parameters are set.
   INTEGER :: iparl = 0                 ! parity for states on the lower surface
   INTEGER :: jtotl = 0                 ! total J for states on the lower surface
   INTEGER :: megal = 0                 ! Projection quantum number for states on the lower surface.
   REAL(Kind=WP_Kind) :: pevms1u(2) = 4.d0 ! max energy for which 1-d eigenvectors on upper PES are kept at rhos1u(i).
   REAL(Kind=WP_Kind) :: rhos1u(2)  = 0.d0  ! rho values where the dvr energy parameters are set.
   REAL(Kind=WP_Kind) :: pevms1l(2) = 0.d0  ! max energy for which 1-d eigenvectors on lower PES are kept at rhos1l(i).
   REAL(Kind=WP_Kind) :: rhos1l(2)  = 0.d0 ! rho values where the dvr energy parameters are set.
   REAL(Kind=WP_Kind) :: rhoendl    = 10.d0  ! the point after which the lower state surface functions are no longer calculated and the dipole
                                  !     is assumed to be zero
   NAMELIST/sysu/jevenu,nrhos1u,pevms1u,iparu,jtotu,megau, cstestu,pointu,rhos1u

!  --------------------------------------------------------------------
!     sysl contains DATA needed for lower PES
!  --------------------------------------------------------------------
!     IF nsurf=1, sysl is not used.
!     most parameters in sysl are analogous to those in sysu except that
!     rhoendl is the point after which the lower state surface
!     functions are no longer calculated and the dipole matrix element
!     is assumed to be zero

   NAMELIST/sysl/jevenl,nrhos1l,pevms1l,iparl,jtotl,megal, cstestl,pointl,rhos1l,rhoendl

!  --------------------------------------------------------------------
!     rstart contains DATA needed for restart option
!  --------------------------------------------------------------------
!    restart - LOGICAL switch IF true will enable restart option
!    irhost - number of rho where you want to restart
!              irhost must be odd !!

   NAMELIST/rstart/restart,nbreak,irhost
END MODULE DVRNL_Module
