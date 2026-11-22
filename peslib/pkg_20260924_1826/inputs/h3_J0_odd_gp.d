# H3 rate, J=0, odd exchange, GP (two-state diabatic)
# Kendrick 2018 Fig.2 comparison: ultracold 1 K -> 1.4 K, emax=4.75 (full double-state basis)
&input
  mass = 1, 1, 1
  jtot = 0
  ipar = 1
  jpar = -1
  jmax = 25
  kmax = 0
  rmax = 12.0
  mtr = 250
  emax = 4.75
  enrg = 2.1610001
  dnrg = 0.0000001
  nnrg = 300
  nout = 0
  jout = 0
  nmax = 30000
  NPES = 2
  PESMODE = 'H3'
  rmin_override = 0.65
&end
