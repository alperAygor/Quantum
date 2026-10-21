#!/usr/bin/env bash
# Tüm notebook'ları sırayla çalıştırır ve çıktılarıyla birlikte kaydeder.
# Kullanım:  QFE_PRESET=demo bash tools/run_all.sh
set -e
cd "$(dirname "$0")/../notebooks"
export QFE_PRESET=${QFE_PRESET:-demo}
echo "Ön ayar: $QFE_PRESET"
for nb in 0*.ipynb; do
  echo "=== $nb  ($(date +%H:%M:%S))"
  jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 "$nb"
done
echo "Bitti ($(date +%H:%M:%S)). Sonuçlar: results/$QFE_PRESET/"
