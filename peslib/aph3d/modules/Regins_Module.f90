MODULE Regins_Module
!  NOTE the following inequalities must be satisfied:
!           START <= ENDAPH <= ENDDELVE <= FINISH
!  and
!           START <= SwitchToVIVS <= FINISH
!
!  The following schematic may be useful.
!
!             Log-derivative              VIVS
!          --------------------- -------------------------
!         |         |           |            |            |
!      start      endaph     SwitchToVIVS       enddelve      finish
!          --------- ------------------------ ------------
!             APH            Delves             Jacobi
!
!                      distance ----->

   USE Numeric_Kinds_Module
   SAVE
   REAL(Kind=WP_Kind) :: start=2.1d0         ! starting distance used in the propagation.
   REAL(Kind=WP_Kind) :: endaph=8.1d0        ! distance which signals the end of the Delves region.
   REAL(Kind=WP_Kind) :: SwitchToVIVS=8.1d0  ! distance to switch from the Log-derivative propagator to the VIVS propagator.
   REAL(Kind=WP_Kind) :: enddelve=8.1d0      ! distance which signals the end of the Delves region.
   REAL(Kind=WP_Kind) :: finish=8.1d0        ! maximum propagation distance. Distance at which asymptotic functions are used to obtain the S-matr
   NAMELIST/regins/ start, SwitchToVIVS, finish, endaph, enddelve
END MODULE Regins_Module