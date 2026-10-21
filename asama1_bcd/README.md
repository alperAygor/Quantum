# Aşama 1: β-siklodekstrin + siklopentanol

Projenin doğrulama aşaması. Küçük, deneysel verisi bilinen bir konak–misafir sisteminde
(deney: ΔG = −3,05 kcal/mol, Rekharsky 1997) klasik + kuantum iş akışının doğru çalıştığı gösterildi.
Sonuç: bağlanma dispersiyon güdümlü, kuantum düzeltmesi pratikte sıfır, yani triyaj "klasik yeterli" der.

Bu klasörde aynı sistem için iki ayrı çalışma var:

| Klasör | Ne | Araçlar |
|---|---|---|
| [`gromacs_vqe/`](gromacs_vqe/) | Asıl araştırma: 10 ns MD, 37 pencereli umbrella sampling + WHAM, GFN2-xTB etkileşimi, MFCC parçalama, misafir/bölge VQE | GROMACS 2026.1 (GROMOS 54A7/ATB), tblite, PySCF, PennyLane |
| [`qfe_demo/`](qfe_demo/) | Eğitim simülasyonu: Günther et al. (JCTC 2026) FreeQuantum pipeline'ının laptop ölçekli kopyası, ekibe adım adım anlatım için 8 notebook | OpenMM (GAFF/q4md-CD), pymbar, PySCF, scikit-learn, OpenFermion |

İkisi farklı kuvvet alanları ve serbest enerji yöntemleri (PMF ile alkimik FEP) kullanıyor, bu yüzden
MM sayıları birebir aynı çıkmaz. `gromacs_vqe` sonuçlarının yorumu için
[`../docs/01_BCD_PMF_sonuc_ozeti.md`](../docs/01_BCD_PMF_sonuc_ozeti.md) dosyasına bakın.
