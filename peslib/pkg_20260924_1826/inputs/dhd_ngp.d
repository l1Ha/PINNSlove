# D+HD(v=4,j=0), J=0, NGP (single surface), fine energy grid (mass=2,1,2 = Kendrick 2016 reactant D+HD)
# HD(v=4) threshold 1.9051 eV (code) = 1.91 eV (Kendrick); E_coll 1.03-10.1 K at 0.023 K resolution
&input
  mass = 2, 1, 2
  jtot = 0
  ipar = 1
  jpar = 0
  jmax = 15
  kmax = 0
  rmax = 14.0
  mtr = 200
  emax = 2.50
  enrg = 1.9051972
  dnrg = 0.000002
  nnrg = 435
  nout = 0
  jout = 0
  nmax = 30000
  NPES = 1
  PESMODE = 'H3'
  rmin_override = 0.8
&end
