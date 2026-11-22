# H3 rate, J=1, physical parity ipar=-1, even exchange (para H2 j=0), NGP
# (this block carries the l=1 partial wave -> Kendrick's 1 K shape resonance)
&input
  mass = 1, 1, 1
  jtot = 1
  ipar = -1
  jpar = 1
  jmax = 25
  kmax = 1
  rmax = 12.0
  mtr = 250
  emax = 4.75
  enrg = 2.1610001
  dnrg = 0.0000001
  nnrg = 300
  nout = 0
  jout = 0
  nmax = 30000
  NPES = 1
  PESMODE = 'H3'
  rmin_override = 0.65
&end
