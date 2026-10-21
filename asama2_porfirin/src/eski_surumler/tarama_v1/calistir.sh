#!/bin/bash
# Kullanim: bash calistir.sh  PARALEL_IS_SAYISI   (ornek: bash calistir.sh 4)
N=${1:-2}; T=$(( $(nproc) / N )); export OMP_NUM_THREADS=$T
MOLS="porfin klorin bakteriyoklorin tetrametilporfin porfirazin Zn-porfin Mg-porfin Pd-porfin"
echo "Geometriler..."; python3 geo.py $MOLS
echo "Kuantum triyaj: $N paralel is, her biri $T thread"
echo $MOLS | tr ' ' '\n' | xargs -P $N -I{} sh -c 'python3 tarama_qm.py {} > {}.log 2>&1; echo "{} bitti"'
python3 ozet.py
