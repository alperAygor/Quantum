# Kuantum Triyaj

Fotodinamik terapi (PDT) fotosensitizörleri için kuantum-hibrit hesaplama platformu prototipi.

> **English summary.** *Quantum Triage* is a prototype workflow that decides **which molecules actually need a
> quantum computer**. A cheap classical triage step measures how multi-configurational a molecule's excited state is;
> only those that need it get a small quantum active space (CAS(4e,4o), solved with VQE), and dynamic correlation is
> added classically (NEVPT2). Applied to porphyrin-based photosensitizers for photodynamic therapy, it reaches a 0.07 eV
> error in the S1 excitation energy versus experiment (TD-DFT/B3LYP: 0.39 eV). The circuits were run on IBM's 156-qubit
> `ibm_kingston`; a compact 5-qubit encoding developed here halves the hardware error of the original 8-qubit circuit.
> No quantum advantage is claimed: these 8-qubit problems are exactly solvable on a laptop. Documentation is in Turkish.

**Fikir:** Kuantum hesaplama kıt, pahalı ve gürültülü bir kaynak. Her molekülü kuantumla hesaplamak yerine
bir *triyaj* adımı, molekülün elektron yapısının çok-konfigürasyonlu olup olmadığını ölçer. Kuantum kaynağı
yalnızca gereken moleküle ve yalnızca kritik bölgeye (aktif uzay) harcanır; gerisi klasik yöntemle hesaplanır.

Triyaj kuralı: S1 halinde ikinci elektron konfigürasyonunun ağırlığı **> 0,2 ise "kuantum gerekli"**.

## Aşamalar ve sonuçlar

| Aşama | Soru | Sonuç | Klasör |
|---|---|---|---|
| 1. β-siklodekstrin + siklopentanol | İş akışı doğru mu, bu sistemde kuantum gerekli mi? | ΔG° = −2,7 ± 0,6 kcal/mol (deney −3,05); kuantum düzeltmesi ≈ 0 → **kuantum gereksiz** | [`asama1_bcd/`](asama1_bcd/) |
| 2. Porfirin ailesi (PDT) | Hangi fotosensitizörde kuantum gerekli? | S1 hatası 0,07 eV (TD-DFT: 0,39 eV); 5 molekül "kuantum gerekli" | [`asama2_porfirin/`](asama2_porfirin/) |
| IBM Quantum | Devre gerçek donanımda ne kadar doğru? | Kompakt 5 kübitlik devre ile porfin S1 hatası +0,41 eV (8 kübitlik ilk devre: +0,93 eV) | [`ibm/`](ibm/) |

### IBM donanım sonuçları (`ibm_kingston`, S1 uyarılma enerjisi, NEVPT2 öncesi)

| Molekül (triyaj sınıfı) | Yöntem | Kübit | 2-kübit kapı (S0/S1) | Tam çözüm (CASCI) | Ham fark | Düzeltilmiş fark |
|---|---|---|---|---|---|---|
| Porfin (kuantum gerekli) | `ibm/ibm_vqe.py` (ADAPT, Jordan–Wigner) | 8 | 264 / 170 | 3,617 eV | −2,48 eV | +0,93 ± 0,19 eV |
| Porfin (kuantum gerekli) | `ibm/ibm_vqe_kompakt.py` | 5 | 60 / 60 | 3,617 eV | −0,90 eV | **+0,41 ± 0,14 eV** |
| Klorin (klasik yeterli) | `ibm/ibm_vqe_kompakt.py` (3072 atış) | 5 | 60 / 60 | 3,453 eV | −1,39 eV | **−0,34 ± 0,21 eV** |

**Kompakt yöntem:** Elektron sayısı ve singlet spin simetrisi sabit olduğu için 8 kübitin 256 durumundan yalnızca
21'i fiziksel. Bu 21 durum 5 kübite yerleştirilir. Hesaplanan şey aynıdır (parametreler klasikte, enerji donanımda
ölçülür), ama devre ~4 kat kısalır ve donanım hatası yaklaşık yarıya iner. Ham veri ve iş kimlikleri
[`ibm/sonuclar/`](ibm/sonuclar/) klasöründe.

Bugünkü donanım hatası (~0,3–0,4 eV) hâlâ kimyasal doğruluk hedefinin (~0,1 eV) üstünde.
Ayrıntılı teknik rapor: [`docs/00_Master_Teknik_Rapor.md`](docs/00_Master_Teknik_Rapor.md).

## Klasör yapısı

```
├── README.md
├── docs/                      ← teknik rapor, veri yorumu, kaynakça, literatür, IBM planı, sunum metni
├── asama1_bcd/
│   ├── gromacs_vqe/           ← β-CD çalışması: GROMACS PMF, xTB, MFCC, VQE (src, data, molecules, visuals)
│   └── qfe_demo/              ← aynı sistemde FreeQuantum pipeline'ının eğitim simülasyonu (notebook'lar, testler)
├── asama2_porfirin/
│   ├── src/                   ← 03 porfin kuantum · 04 tarama v3 · 05 görselleştirme · 06 TD-DFT referansı
│   ├── data/                  ← CSV sonuçlar, ham JSON ve loglar
│   ├── molecules/             ← SMILES, PDB, XYZ, SDF, molden
│   └── visuals/               ← grafikler, ChimeraX görselleri ve .cxc betikleri, görsel envanteri
├── ibm/                       ← IBM Quantum betikleri, donanım sonuçları, sandbox testi
└── araclar/                   ← arşiv oluşturma betiği (tarihsel)
```

## Kurulum ve çalıştırma

### Eğitim simülasyonu (β-CD, notebook'lar)
Ayrıntılar: [`asama1_bcd/qfe_demo/README.md`](asama1_bcd/qfe_demo/README.md).
```bash
cd asama1_bcd/qfe_demo && mamba env create -f environment.yml && conda activate qfe
```
Testler birkaç saniyede çalışır:
```bash
python -m pytest tests
```

### IBM Quantum
Sürümler [`ibm/requirements.txt`](ibm/requirements.txt) dosyasında sabitlendi. numpy 2.5 / scipy 1.18 macOS'ta PySCF'i
rastgele çökertiyor ve sessizce yanlış Hamiltonyen üretiyordu.
```bash
python3 -m venv .venv-ibm && .venv-ibm/bin/pip install -r ibm/requirements.txt
```
IBM Quantum hesabını bir kez kaydedin. Anahtar ekranda görünmez ve yalnızca Qiskit'in kendi ayar dosyasına yazılır:
```bash
.venv-ibm/bin/python ibm/hesap_kaydet.py
```
Çalışma klasöründe molekül geometrisiyle çalıştırın (örnek: porfin):
```bash
mkdir -p ibm/calisma && cp asama2_porfirin/molecules/xyz/porfin.xyz ibm/calisma/ && cd ibm/calisma
```
```bash
../../.venv-ibm/bin/python ../ibm_vqe_kompakt.py hazirla porfin
```
```bash
../../.venv-ibm/bin/python ../ibm_vqe_kompakt.py gonder porfin --backend=ibm_kingston
```
```bash
../../.venv-ibm/bin/python ../ibm_vqe_kompakt.py sonuc porfin
```
`hazirla` ve `gonder --sahte` yereldir ve kota harcamaz. `gonder` göndermeden önce onay ister. Ücretsiz Open Plan
ayda 10 dakika QPU süresi verir; bir porfin ölçümü ~4,5 dakika sürer.

### Porfirin taraması
[`asama2_porfirin/src/04_tarama_v3/BENIOKU.txt`](asama2_porfirin/src/04_tarama_v3/BENIOKU.txt) ve teknik raporun 6. bölümü.

## Güvenlik ve gizlilik

- Depoda API anahtarı, token veya IBM hesap kimliği **yoktur**. `.gitignore`, IBM iş dosyalarını ve yerel ortamları dışarıda tutar.
- Kişisel bilgisayar yolları ve e-posta adresleri `<KULLANICI>` / `[iletisim e-postasi]` yer tutucularıyla değiştirildi.
- Telif hakkı korumalı makale PDF'leri depoya konmaz; kaynakça [`docs/References.md`](docs/References.md) dosyasında.

## Bilinen eksikler

- Bazı dosyalar arşive hiç alınamadı: GROMACS `.mdp` dosyaları, yörüngeler, `misafir_qm.py` ve bazı JSON'lar.
  Listesi [`asama2_porfirin/src/KAYIP_DOSYALAR.txt`](asama2_porfirin/src/KAYIP_DOSYALAR.txt) dosyasında.
- TD-DFT referansının tetrametilporfin ve Pd-porfin hesapları tamamlanmadı.
- Bakteriyoklorin (NIR bölgesi) için 4 orbitalli model yetersiz (0,49 eV hata); CAS(8,8) = 16 kübit gerekli.
- Tarama betikleri orijinal ortamda (`~/tarama`, WSL2) düz bir klasörde çalıştırıldı ve girdi JSON'larını çalışılan
  klasörde arıyor. Çıktı klasörleri ortam değişkenleriyle seçilir: `CIKTI_DIR` (matplotlib figürleri),
  `GORSEL_DIR` / `GORSEL_DIR_CHIMERAX` (ChimeraX dosyaları), `KOPYA_DIR` (`calistir_v4.sh` kopyası).
- `asama2_porfirin/visuals/chimerax_betikleri/*.cxc` dosyaları tarihsel çıktılardır. Yeniden üretmek için
  `gorsel_hq.py` ve `gercek_ilac.py` betiklerini çalıştırın.

## Dürüst çerçeve

Kuantum avantajı iddia edilmiyor. 8 kübitlik hesap klasik bilgisayarda da tam çözülebiliyor. TD-DFT'ye karşı
doğruluk kazancı çok-konfigürasyonlu yaklaşımdan (CAS + NEVPT2) geliyor ve NEVPT2 klasik bir yöntem.
Ayrıntı için teknik raporun 1.4. bölümüne bakın.

## Emeği geçenler

- **Sinan Ahmet Yalçınkaya** (Yıldız Teknik Üniversitesi): proje fikri, hesaplar, porfirin taraması, ilk IBM ölçümü,
  YTÜ Startup House başvurusu.
- **[@alperAygor](https://github.com/alperAygor):** deponun düzenlenmesi, hata düzeltmeleri ve testler, kompakt 5 kübitlik
  IBM yöntemi ve 21 Ekim 2026 donanım ölçümleri.

## Lisans

Henüz bir lisans seçilmedi. Lisans eklenene kadar tüm hakları sahiplerine aittir; kodu kullanmak veya
yeniden dağıtmak için izin alın.
