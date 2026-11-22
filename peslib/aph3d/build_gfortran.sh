#!/bin/bash
#-----------------------------------------------------------------------
# build_gfortran.sh -- rebuild the whole APH3D port with gfortran.
#   usage:  ./build_gfortran.sh [-j N]
# Ported to gfortran 2026-10; sources live in modules/ and src/, the
# system potentials in potentials/ (e.g. potentials/H3_DMBE), and the
# drivers in Main_Programs/.
#-----------------------------------------------------------------------
set -u
cd "$(dirname "$0")"

FC=${FC:-gfortran}
FFLAGS="-g -O2 -std=legacy -fallow-argument-mismatch -ffree-line-length-none -fmax-errors=5 -mcmodel=medium"
# BLAS/LAPACK link flags; override on machines without system libraries, e.g.
#   export EXTRALIBS="-L$HOME/abc_test/lib -lopenblas -Wl,-rpath,$HOME/abc_test/lib"
if [ -z "${EXTRALIBS:-}" ] && [ -f "/storage/home/lih/abc_test/lib/libopenblas.so" ]; then
    EXTRALIBS="-L/storage/home/lih/abc_test/lib -lopenblas -Wl,-rpath,/storage/home/lih/abc_test/lib"
fi
LIBS=${EXTRALIBS:--lblas -llapack}
OBJ=objs/gfortran
mkdir -p "$OBJ" lib bin

# Files that contain OpenMP directives get -fopenmp; -fopenmp implies
# -frecursive in gfortran, which would put the (very large) automatic arrays
# of several legacy routines on the stack and crash them, so it must NOT be
# applied globally.
OMP_FILES=" PotElmnt PotElmnt2 CalcPot CalcPot_2x2 OvrBas "
compile() {  # compile $1 -> $OBJ/<base>.o ; echo error output if any
  local f="$1" b o out extra=""
  b=$(basename "${f%.*}")
  case "$OMP_FILES" in *" $b "*) extra="-fopenmp" ;; esac
  o="$OBJ/$b.o"
  out=$("$FC" $FFLAGS $extra -J "$OBJ" -I "$OBJ" -c "$f" -o "$o" 2>&1)
  if echo "$out" | grep -qE "^(Error|Fatal)"; then
      echo "$out" | grep -E "^(Error|Fatal)" | head -3 | sed "s|^|  [$f] |"
      return 1
  fi
  return 0
}

#-----------------------------------------------------------------------
# 1) modules: iterate until all compile (module interdependencies)
#-----------------------------------------------------------------------
echo "== modules =="
pend=$(ls modules/*.f90 2>/dev/null)
pass=0
while [ -n "$pend" ] && [ $pass -lt 8 ]; do
   pass=$((pass+1)); next=""; nfail=0; ndone=0
   for f in $pend; do
      if compile "$f"; then ndone=$((ndone+1)); else next="$next $f"; nfail=$((nfail+1)); fi
   done
   echo "  pass $pass: compiled $ndone, pending $nfail"
   # always run a second pass: re-export hubs (NameList_Module, FileUnits_Module,
   # ...) snapshot the properties of the modules they USE, so they must be
   # recompiled after any edited dependency.
   if [ $nfail -eq 0 ] && [ $pass -ge 2 ]; then break; fi
   pend="$next"
   [ -z "$pend" ] && pend=$(ls modules/*.f90 2>/dev/null)
   [ $ndone -eq 0 ] && [ $nfail -gt 0 ] && { echo "  STUCK on:$next"; break; }
done

#-----------------------------------------------------------------------
# 2) sources
#-----------------------------------------------------------------------
echo "== src =="
nfail=0
for f in src/*.f90 src/*.f; do
   [ -e "$f" ] || continue
   compile "$f" || nfail=$((nfail+1))
done
echo "  src failures: $nfail"

#-----------------------------------------------------------------------
# 3) potentials (system surface routines) and PES support libraries
#-----------------------------------------------------------------------
echo "== potentials =="
nfail=0
for f in potentials/*/*.f90 potentials/*/*.f KRb2_PESs/RKHS.f90 KRb2_PESs/krb2_pes.f90; do
   [ -e "$f" ] || continue
   compile "$f" || nfail=$((nfail+1))
done
echo "  potential failures: $nfail"

#-----------------------------------------------------------------------
# 4) archive (all local objects except the drivers)
#-----------------------------------------------------------------------
echo "== archive =="
rm -f lib/libaph3d_gfortran.a
find "$OBJ" -maxdepth 1 -name '*.o' ! -name 'Program_*' -print0 \
   | xargs -0 ar rcs lib/libaph3d_gfortran.a
ls -la lib/libaph3d_gfortran.a

#-----------------------------------------------------------------------
# 5) drivers
#-----------------------------------------------------------------------
echo "== drivers =="
for f in Main_Programs/*.f90; do
   b=$(basename "${f%.*}")
   out=$("$FC" $FFLAGS -J "$OBJ" -I "$OBJ" -c "$f" -o "$OBJ/$b.o" 2>&1)
   if echo "$out" | grep -qE "^(Error|Fatal)"; then
      echo "  SKIP $b (compile error)"
      continue
   fi
   case "$b" in
     Program_H3_DMBE)        exe=xeq_H3_DMBE.x ;;
     Program_H3_PKH3)        exe=xeq_H3_PKH3.x ;;
     Program_H3_LSTH)        exe=xeq_H3_LSTH.x ;;
     Program_H3_CCI)         exe=xeq_H3_CCI.x ;;
     Program_H3_PAN)         exe=xeq_H3_PAN.x ;;
     Program_HeH2)           exe=xeq_HeH2.x ;;
     Program_KRb2)           exe=xeq_KRb2.x ;;
     *)                      exe="" ;;
   esac
   [ -n "$exe" ] || continue
   if "$FC" -g -mcmodel=medium -fopenmp -o "bin/$exe" "$OBJ/$b.o" lib/libaph3d_gfortran.a $LIBS 2> /tmp/link_$b.err; then
      echo "  linked bin/$exe"
   else
      echo "  LINK FAIL $b:"; head -3 /tmp/link_$b.err | sed 's/^/    /'
   fi
done
echo "== done =="
