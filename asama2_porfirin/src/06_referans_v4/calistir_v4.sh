#!/bin/bash
# v4: klasik referans (TD-DFT) + buyuk aktif uzay (16 kubit esdegeri) + final figurler.  ~/tarama icinde calistir.
export OMP_NUM_THREADS=$(nproc)
MOLS="porfin klorin Zn-porfin Mg-porfin porfirazin bakteriyoklorin tetrametilporfin Pd-porfin"
echo "== 1) TD-DFT referansi (B3LYP + CAM-B3LYP), molekul basina ~5-15 dk =="
for m in $MOLS; do
  if [ -f ${m}_tddft.json ] && grep -q camb3lyp ${m}_tddft.json; then echo "$m: TD-DFT zaten var"; continue; fi
  python3 tddft_referans.py $m >> ${m}_tddft.log 2>&1; tail -3 ${m}_tddft.log | head -2
done
python3 final_figurler.py
if [ -n "$KOPYA_DIR" ]; then mkdir -p "$KOPYA_DIR"; cp F*_*.png final_tablo.csv final_ozet.txt "$KOPYA_DIR"/ 2>/dev/null; fi
[ -n "$KOPYA_DIR" ] && echo "== Figurler $KOPYA_DIR klasorune kopyalandi =="
echo "== 2) Buyuk aktif uzay CAS(8,8)+NEVPT2: bakteriyoklorin, porfin (her biri ~1-3 saat) =="
for m in bakteriyoklorin porfin; do
  if [ -f ${m}_cas88.json ] && grep -q NEVPT2 ${m}_cas88.json; then echo "$m: CAS(8,8) zaten var"; continue; fi
  python3 cas_buyuk.py $m 8 8 >> ${m}_cas88.log 2>&1; tail -1 ${m}_cas88.log | cut -c1-200
  python3 final_figurler.py > /dev/null; cp F*_*.png final_tablo.csv final_ozet.txt $D/ 2>/dev/null
done
echo "== HEPSI BITTI =="; cat final_ozet.txt
