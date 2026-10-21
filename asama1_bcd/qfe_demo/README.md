# FreeQuantum pipeline'ı: laptop ölçeğinde simülasyon

Bu proje, **Günther et al., "How to Use Quantum Computers for Biomolecular Free Energies",
J. Chem. Theory Comput. 2026, 22, 4329–4345** makalesindeki FreeQuantum pipeline'ının küçültülmüş,
eğitim amaçlı bir kopyası. Her aşamayı kendi bilgisayarınızda çalıştırıp ekibinize adım adım
anlatabilmeniz için hazırlandı.

- **Sistem:** β-siklodekstrin (konakçı) + siklopentanol (misafir), suda
- **Gerçek veri:** Girdi dosyaları MobleyLab *benchmarksets* (CD Set 2, guest 3) setinden geliyor.
  Deneysel ΔG_bağlanma = −3.050 kcal/mol = **−12.76 kJ/mol** (Rekharsky et al., *J. Phys. Chem. B* 1997, 101, 87).
- **Katmanlar:** MM-FEP → QM/MM (DFT) + ML1 → QM/QM/MM (CASCI+NEVPT2) + ML2 → kuantum algoritmaları
  (VQE, QPE, Trotter) → kaynak tahmini

> **Not:** Notebook'larda şu an görünen çıktılar, paketi doğrulamak için `test` ön ayarıyla 2 çekirdekli
> bir CPU'da alındı. Tüm akışın hatasız çalıştığını gösteriyorlar ama sayılar kısa örnekleme yüzünden gürültülü.
> Kendi GPU'nuzda `demo` ile yeniden çalıştırdığınızda çıktılar güncellenir.

## Klasörler

```
qfe_demo/
├── README.md              ← bu dosya
├── environment.yml        ← conda ortamı
├── qfe/                   ← pipeline kodu (notebook'lar bunu çağırıyor)
│   ├── config.py          ← ön ayarlar (test / demo / standart) ve sabitler
│   ├── mmsys.py           ← veri indirme, MM sistemleri, dengeleme
│   ├── fep.py             ← alkimik FEP + MBAR
│   ├── qmmm.py            ← QM bölgesi, DFT ve CASCI/NEVPT2 motorları
│   ├── ml.py              ← tanımlayıcılar, Gauss süreci, aktif öğrenme, Zwanzig
│   ├── quantum.py         ← qubit Hamiltoniyeni, VQE, QPE, Trotter, kaynak tahmini
│   ├── store.py           ← aşamalar arası veri (makaledeki veritabanının basit karşılığı)
│   └── plotting.py
├── notebooks/             ← 00 … 07: ekibe anlatım sırası
├── tools/run_all.sh       ← bütün notebook'ları sırayla çalıştırır
├── tests/                 ← hızlı birim testleri (python -m pytest tests)
├── data/                  ← indirilen girdiler (ilk çalıştırmada oluşur)
└── results/<ön_ayar>/     ← sonuçlar, süreler, Hamiltoniyen dosyaları
```

## Kurulum: Windows + NVIDIA GPU (WSL2 ile)

PySCF Windows'ta doğrudan çalışmıyor. Bu yüzden **WSL2 (Windows içinde Linux)** kullanıyoruz.
OpenMM, WSL2 içinden NVIDIA GPU'nuzu CUDA ile kullanabiliyor.

1. **WSL2 + Ubuntu kurun.** PowerShell'i *yönetici olarak* açın:
   ```powershell
   wsl --install -d Ubuntu
   ```
   Bilgisayarı yeniden başlatın ve Ubuntu kullanıcı adı/şifresini belirleyin.

2. **NVIDIA sürücüsü:** Windows tarafında güncel NVIDIA sürücüsü kurulu olsun. WSL'nin CUDA desteği
   bu sürücüyle birlikte geliyor. **WSL içine ayrıca Linux NVIDIA sürücüsü kurmayın.**
   Ubuntu terminalinde `nvidia-smi` komutu GPU'nuzu göstermeli.

3. **Miniforge (conda/mamba) kurun.** Ubuntu terminalinde:
   ```bash
   curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
   bash Miniforge3-Linux-x86_64.sh -b
   ~/miniforge3/bin/conda init bash && exec bash
   ```

4. **Projeyi WSL'ye kopyalayın.** Linux dosya sistemi çok daha hızlı. Zip'i Windows'ta `İndirilenler`
   klasörüne kaydettiyseniz:
   ```bash
   cd ~ && cp /mnt/c/Users/<KULLANICI_ADINIZ>/Downloads/qfe_demo.zip . && unzip qfe_demo.zip && cd qfe_demo
   ```

5. **Ortamı oluşturun** (5–15 dakika sürebilir):
   ```bash
   mamba env create -f environment.yml
   conda activate qfe
   python -m openmm.testInstallation      # "CUDA" platformu listede görünmeli
   ```

6. **Çalıştırın:**
   ```bash
   jupyter lab
   ```
   Tarayıcıda açılan adresten `notebooks/` klasörüne girin ve 00'dan 07'ye kadar sırayla açıp çalıştırın.
   Ya da hepsini arka arkaya çalıştırın:
   ```bash
   QFE_PRESET=demo bash tools/run_all.sh
   ```

> Linux/macOS kullanıyorsanız 3. adımdan başlayın. macOS'ta CUDA olmadığı için MD CPU'da (veya OpenCL ile) çalışır.

## Ön ayarlar ve tahmini süreler

Ön ayarı `QFE_PRESET` ortam değişkeniyle seçiyorsunuz (varsayılan `demo`). Ayarlar `qfe/config.py` dosyasında.

| Ön ayar | Amaç | GPU'lu laptop (tahmini) | Sadece CPU |
|---|---|---|---|
| `test` | Kodun çalıştığını görmek; sonuçlar gürültülü | ~15–25 dk | ~40–60 dk |
| `demo` | Ekibe gösterim; kabaca anlamlı sonuçlar | ~1–2 saat | önerilmez |
| `standart` | Daha güvenilir sayılar | birkaç saat + QM için gece boyu | önerilmez |

Süreyi en çok etkileyen iki şey var: MD (02, GPU'da hızlı) ve DFT hesapları (03–05, CPU'da çalışır ve
çekirdek sayısıyla hızlanır). Her aşama ölçtüğü süreyi `results/<ön_ayar>/timings.json` dosyasına yazıyor.
07. notebook bu sürelerden "gerçek projede nerede HPC gerekir" tablosunu üretiyor.

### Ne tür sonuçlar beklemeli?

Bu paketi hazırlarken FEP kodunu 2 çekirdekli bir CPU'da, `demo` λ şemasıyla (13 durum × 70 ps) ayrıca
doğruladık. Sonuç **ΔG_MM = −16.2 ± 2.8 kJ/mol**, deney ise −12.8 kJ/mol. GAFF/TIP3P'nin bu tür
sistemlerde bağlanmayı biraz abartması literatürle uyumlu. `test` ön ayarı çok kısa örnekleme yaptığı
için ΔG'yi ±10 kJ/mol veya daha fazla saptırabilir; o ön ayar yalnızca kodun çalıştığını görmek için.

## Notebook'lar (ekibe anlatım sırası)

| # | İçerik | Makaledeki karşılığı |
|---|---|---|
| 00 | Genel bakış, eşleştirme tablosu, ortam kontrolü | Bölüm 1, Şekil 1–2 |
| 01 | Gerçek veri, kompleks ve çözücü bacakları, NPT dengeleme | Şekil 2 (gri) |
| 02 | Alkimik FEP + MBAR → ΔG(MM), anlık görüntü havuzu | Bölüm 2.1, Denklem 3–5 |
| 03 | QM bölgesi, elektrostatik gömme, DFT maliyeti, ΔU dağılımı | Bölüm 2.3, Denklem 6 |
| 04 | ML1 + aktif öğrenme → ΔG(MM+ML1) | Bölüm 2.2, Şekil 6–7 |
| 05 | Kuantum çekirdeği CASCI+NEVPT2, ML2 → ΔG(MM+ML2) | Bölüm 2.3–2.4, Denklem 7 |
| 06 | Qubit Hamiltoniyeni, VQE, QPE devre simülasyonu, Trotter, HF örtüşmesi | Bölüm 3, Şekil 4d |
| 07 | Kaynak tahmini, hata yayılımı, sonuç özeti, HPC/kuantum yol haritası | Bölüm 3, Şekil 4a–c |

## Makaleye göre bilinçli basitleştirmeler

Ekibe bunları açıkça söyleyin. Bu bir **eğitim simülasyonu**, makalenin yeniden üretimi değil.

1. **Sistem:** Protein + Ru ilacı yerine küçük bir konakçı–misafir sistemi kullanıyoruz. Siklopentanol kapalı
   kabuklu ve zayıf korelasyonlu (makale Şekil 1'in sol alt köşesi). Kuantum avantajı burada beklenmiyor;
   amaç iş akışını göstermek.
2. **Kuantum çekirdeği:** Makalede çekirdek, DFT bölgesinin içine Huzinaga gömmesiyle yerleştirilmiş Ru d-orbitalleri.
   Burada misafir zaten küçük olduğu için çekirdek bütün misafir ve aktif uzay CAS(4,4).
3. **ML potansiyeli:** Makale kuvvet veren bir ML potansiyeliyle MD ve 10 ps NEQ anahtarlama yapıyor. Biz ML'i
   yalnızca enerji farkı (ΔU) tahmincisi olarak kullanıp tek adımlı pertürbasyon (Zwanzig + kümülant) uyguluyoruz.
4. **QM seviyesi:** D3 dispersiyon düzeltmesi yok. Misafir–çevre van der Waals terimi MM'den geliyor.
5. **Kaynak tahminleri:** Mertebe formülleri kullanılıyor. Makaledeki "double factorization" vb. ile gerçek sayılar daha düşük.
6. **İstatistik:** `test` ve `demo` ön ayarlarında örnekleme kısa, hata çubukları büyük. ΔG'nin deneye ne
   kadar yakın çıktığı bu ölçekte şans eseri olabilir.

## Gerçek kuantum donanımında küçük bir deneme (isteğe bağlı)

06. notebook, `results/<ön_ayar>/H_cas22_4qubit.json` dosyasına 4 qubitlik bir Hamiltoniyen yazıyor. Qiskit ile:

```python
import json
from qiskit.quantum_info import SparsePauliOp
pl = json.load(open("results/demo/H_cas22_4qubit.json"))
H = SparsePauliOp.from_list([(lab, complex(re, im)) for lab, (re, im) in pl])
# H'yi qiskit-ibm-runtime'daki Estimator ile bir VQE ansatz'ına uygulayabilirsiniz (IBM Quantum ücretsiz plan).
```

## Testler

`qfe` paketinin temel kontrolleri birkaç saniyede çalışır. MD ya da uzun QM hesabı yapmazlar:

```bash
python -m pytest tests
```

Kontrol edilenler:
- qubit Hamiltoniyeninin özdeğeri = PySCF CASCI
- VQE tam çözüme ulaşıyor; QPE tepesi doğru enerjide; Trotter hatası adım sayısıyla azalıyor
- Zwanzig/kümülant analitik sonuçla aynı; standart durum düzeltmesi doğru
- havuz ya da QM ayarı değişince eski QM etiketleri yeniden kullanılmıyor

## Sorun giderme

- **`CUDA` platformu görünmüyor:** WSL'de `nvidia-smi` çalışıyor mu? Ortamı `cuda-version` sürücünüzle uyumlu
  olacak şekilde yeniden kurun: `mamba install -c conda-forge openmm "cuda-version=12"`.
- **İndirme hatası / 429:** GitHub kısa süreli sınır koymuş olabilir. Kod 5 kez yeniden deniyor; biraz bekleyip
  tekrar çalıştırın.
- **DFT yavaş:** `export OMP_NUM_THREADS=<çekirdek sayısı>` ayarlayın. `config.py` içinde `qm_cutoff_A`
  değerini (nokta yük kesmesi) küçültebilirsiniz.
- **macOS'ta `No compatible OpenCL platform`:** Kod artık çalışmayan platformları atlayıp CPU'ya geçiyor.
- **Bellek:** `standart` ön ayarındaki CAS(14,14) örtüşme taraması birkaç GB RAM isteyebilir.

## Kaynaklar

- Günther, J. et al. *J. Chem. Theory Comput.* 2026, 22, 4329–4345. FreeQuantum kodu: github.com/Quantum-for-Life/pipeline
- Mobley, D. L.; Heinzelmann, G.; Henriksen, N. M.; Gilson, M. K. *Predicting binding free energies: Frontiers and benchmarks* — github.com/MobleyLab/benchmarksets
- Rekharsky, M. V. et al. *J. Phys. Chem. B* 1997, 101, 87–100 (deneysel veri)
- Henriksen, N. M.; Gilson, M. K. *J. Chem. Theory Comput.* 2017, 13, 4253 (CD kuvvet alanı karşılaştırması)
- Yazılım: OpenMM, openmmtools, pymbar, ParmEd, PySCF, OpenFermion, scikit-learn
