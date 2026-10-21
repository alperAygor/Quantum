#!/bin/bash
# Kullanim: bash calistir.sh [PARALEL_IS]   (varsayilan 1 = en guvenli, bellek tasmaz)
N=${1:-1}; T=$(( $(nproc) / N )); export OMP_NUM_THREADS=$T
MOLS="porfin klorin bakteriyoklorin tetrametilporfin porfirazin Zn-porfin Mg-porfin Pd-porfin"
EKSIK=""; for m in $MOLS; do [ -f $m.xyz ] || EKSIK="$EKSIK $m"; done
[ -n "$EKSIK" ] && { echo "Geometriler:$EKSIK"; python3 geo.py $EKSIK; }
KALAN=""; for m in $MOLS; do [ -f ${m}_v3_sonuc.json ] || KALAN="$KALAN $m"; done
echo "Kalan molekuller:$KALAN  ($N paralel is, her biri $T thread)"
echo $KALAN | tr ' ' '\n' | xargs -P $N -I{} sh -c 'python3 tarama_qm.py {} >> {}.log 2>&1; echo "{} bitti: $(tail -c 120 {}.log | tail -1 | cut -c1-80)"'
python3 ozet.py
