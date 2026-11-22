# H3 rate, J=0, even exchange, gp (two-state diabatic), muK grid  (Kendrick 2018 Fig.2 1muK anchor)
# E_coll = 1.16 muK .. 348 muK  (300 pts, linear 1.16 muK step)
&input
  mass = 1, 1, 1
  jtot = 0
  ipar = 1
  jpar = 1
  jmax = 25
  kmax = 0
  rmax = 12.0
  mtr = 250
  emax = 4.75
  enrg = 2.1609116233
  dnrg = 0.00000000001
  nnrg = 300
  nout = 0
  jout = 0
  nmax = 30000
  NPES = 2
  PESMODE = 'H3'
  rmin_override = 0.65
&end
