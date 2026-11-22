MODULE Opts_Module
   USE sfntyp_Module
   USE CSBasis_Module
   USE Schem_Module
   USE Crays_Module
   USE DoEnlvls_Module
   IMPLICIT NONE
   SAVE
   LOGICAL :: lsfunc=.True.    !  Set to true when calculating surface functions.
   LOGICAL :: lmatelm=.True.   !  Set to true when calculating matrix elements.
   LOGICAL :: loverlap=.True.  !  Set to true when calculating overlap matrix elements.
   LOGICAL :: lmesher=.True.   !
   LOGICAL :: laph3d=.False.   !
   LOGICAL :: lxsec=.False.    !
   LOGICAL :: lplot=.False.    !
   LOGICAL :: lmatch=.False.   !
   LOGICAL :: lsample=.False.  !
   NAMELIST/options/  lsfunc, lmatelm, loverlap, lmesher, laph3d, lxsec, lplot, lmatch, lsample,  &
                          cstest, lenlvls, cray, scheme, sfuntype
END MODULE Opts_Module