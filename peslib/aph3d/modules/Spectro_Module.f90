MODULE Spectro_Module
   USE Narran_Module
   USE Numeric_Kinds_Module
   SAVE
   LOGICAL nrigid
   CHARACTER(LEN=8) :: diatom(narran) = ['H-H    ','H-H    ','H-H    ']  !
   
   ! Spectroscopic Constants in Atomic Units
   REAL(Kind=WP_Kind) wzero(narran)
   REAL(Kind=WP_Kind) we(narran)       ! vibrational constant in atomic units
   REAL(Kind=WP_Kind) wexe(narran)     ! anharmonicity constant in atomic units
   REAL(Kind=WP_Kind) weye(narran)     !
   REAL(Kind=WP_Kind) be(narran)       ! rotational constant in atomic units
   REAL(Kind=WP_Kind) de(narran)       ! rotational anharmonicity in atomic units
   REAL(Kind=WP_Kind) he(narran)
   REAL(Kind=WP_Kind) alfe(narran)     ! vibrational-rotational coupling term.
   REAL(Kind=WP_Kind) bete(narran)
   REAL(Kind=WP_Kind) game(narran)
   REAL(Kind=WP_Kind) :: reau(narran) = [0.1401081d+01, 0.1401081d+01, 0.1401081d+01]    ! equilibrium distance
   REAL(Kind=WP_Kind) :: vshift(3)    = [0.0000000d+00, 0.0000000d+00, 0.0000000d+00]    ! potential shift
   
   ! Spectroscopic Constants in Inverse Centemeters
   REAL(Kind=WP_Kind) wzerocm(narran)
   REAL(Kind=WP_Kind) :: wecm(narran)   = [0.4401213d+04, 0.4401213d+04, 0.4401213d+04]  ! vibrational constant in inverse centimeters.
   REAL(Kind=WP_Kind) :: wexecm(narran) = [0.1213360d+03, 0.1213360d+03, 0.1213360d+03]  ! anharmonicity constant in inverse centimeters.
   REAL(Kind=WP_Kind) :: weyecm(narran) = [0.8129000d+00, 0.8129000d+00, 0.8129000d+00]  !
   REAL(Kind=WP_Kind) :: becm(narran)   = [0.6085300d+02, 0.6085300d+02, 0.6085300d+02]  ! rotational constant in inverse centimeters.
   REAL(Kind=WP_Kind) :: decm(narran)   = [0.4710000d-01, 0.4710000d-01, 0.4710000d-01]  ! rotational anharmonicity in inverse centimeters.
   REAL(Kind=WP_Kind) :: alfecm(narran) = [0.3062200d+01, 0.3062200d+01, 0.3062200d+01]  ! vibrational-rotational coupling term.
   REAL(Kind=WP_Kind) hecm(narran)
   REAL(Kind=WP_Kind) betecm(narran)
   REAL(Kind=WP_Kind) gamecm(narran)
   
   NAMELIST/Spectro/ wecm, wexecm, weyecm, becm, decm, alfecm, reau, diatom, vshift
END MODULE Spectro_Module
