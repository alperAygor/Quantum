# Kuantum Triyaj

Fotodinamik terapi (PDT) fotosensitizörleri için kuantum-hibrit hesaplama platformu prototipi.
Hazırlayan: Sinan Ahmet Yalçınkaya (Yıldız Teknik Üniversitesi). YTÜ Startup House başvurusu için hazırlandı.

**Fikir:** Kuantum hesaplama kıt, pahalı ve gürültülü bir kaynak. Her molekülü kuantumla hesaplamak yerine
bir *triyaj* adımı molekülün elektron yapısının çok-konfigürasyonlu olup olmadığını ölçer. Kuantum kaynağı
yalnızca gereken moleküle ve yalnızca kritik bölgeye (aktif uzay) harcanır, gerisi klasik yöntemle hesaplanır.

Triyaj kuralı: S1 halinde ikinci elektron konfigürasyonunun ağırlığı **> 0,2 ise "kuantum gerekli"**.

## Aşamalar

| Aşama | Soru | Sonuç | Klasör |
|---|---|---|---|
| 1. β-siklodekstrin + siklopentanol | İş akışı doğru mu, bu sistemde kuantum gerekli mi? | ΔG° = −2,7 ± 0,6 kcal/mol (deney −3,05); kuantum düzeltmesi ≈ 0 → **kuantum gereksiz** | [`asama1_bcd/`](asama1_bcd/) |
| 2. Porfirin ailesi (PDT) | Hangi fotosensitizörde kuantum gerekli? | S1 hatası 0,07 eV (TD-DFT: 0,39 eV); 5 molekül "kuantum gerekli" | [`asama2_porfirin/`](asama2_porfirin/) |
| IBM Quantum | Devre gerçek donanımda çalışıyor mu? | porfin S1, 8 kübit, ibm_kingston; gürültü düzeltmeli +0,93 ± 0,19 eV sapma | [`ibm/`](ibm/) |

Ayrıntılı teknik rapor: [`docs/00_Master_Teknik_Rapor.md`](docs/00_Master_Teknik_Rapor.md).

## Klasör yapısı

```
kuantum-triyaj/
├── README.md                  ← bu dosya
├── docs/                      ← teknik rapor, veri yorumu, kaynakça, literatür, IBM planı, sunum metni
├── asama1_bcd/
│   ├── gromacs_vqe/           ← asıl β-CD çalışması: GROMACS PMF, xTB, MFCC, VQE (src, data, molecules, visuals)
│   └── qfe_demo/              ← aynı sistemde FreeQuantum pipeline'ının eğitim simülasyonu (notebook'lar)
├── asama2_porfirin/
│   ├── src/                   ← 03 porfin kuantum · 04 tarama v3 · 05 görselleştirme · 06 TD-DFT referansı
│   ├── data/                  ← CSV sonuçlar, ham JSON ve loglar
│   ├── molecules/             ← SMILES, PDB, XYZ, SDF, molden
│   └── visuals/               ← grafikler, ChimeraX görselleri ve .cxc betikleri, görsel envanteri
├── ibm/                       ← IBM Quantum betikleri ve sandbox test çıktıları
└── araclar/                   ← arşiv oluşturma betiği (tarihsel)
```

## Hızlı başlangıç

- **Eğitim simülasyonu (β-CD, notebook'lar):** [`asama1_bcd/qfe_demo/README.md`](asama1_bcd/qfe_demo/README.md).
  Conda ortamı `environment.yml` dosyasında. CPU'da `QFE_PRESET=test` ile yaklaşık bir saat sürer.
- **Porfirin taraması:** [`asama2_porfirin/src/04_tarama_v3/BENIOKU.txt`](asama2_porfirin/src/04_tarama_v3/BENIOKU.txt)
  ve teknik raporun 6. bölümü.
- **IBM Quantum:** [`ibm/BENIOKU_v5.txt`](ibm/BENIOKU_v5.txt). API anahtarı `hesap_kaydet.py` ile gizli
  girilir ve Qiskit'in kendi ayar dosyasında saklanır. Depoya hiçbir anahtar konmaz.

## Bilinen eksikler

- Bazı dosyalar arşive hiç alınamadı: GROMACS `.mdp` dosyaları, yörüngeler, `misafir_qm.py` ve bazı JSON'lar.
  Listesi [`asama2_porfirin/src/KAYIP_DOSYALAR.txt`](asama2_porfirin/src/KAYIP_DOSYALAR.txt) dosyasında.
- TD-DFT referansının tetrametilporfin ve Pd-porfin hesapları tamamlanmadı.
- Betikler orijinal ortamda (`~/tarama`, WSL2) düz bir klasörde çalıştırıldı ve girdi JSON'larını çalışılan klasörde
  arıyor (bkz. teknik rapor, 6. bölüm). Çıktı klasörleri ortam değişkenleriyle seçilir:
  `CIKTI_DIR` (matplotlib figürleri), `GORSEL_DIR` / `GORSEL_DIR_CHIMERAX` (ChimeraX dosyaları; WSL'de Windows yolu),
  `KOPYA_DIR` (`calistir_v4.sh` sonuçlarının kopyası). Varsayılan: çalışılan klasör.
- `asama2_porfirin/visuals/chimerax_betikleri/*.cxc` dosyaları orijinal Windows yollarını içeren tarihsel çıktılardır.
  Yeniden üretmek için `gorsel_hq.py` ve `gercek_ilac.py` betiklerini çalıştırın.

## Dürüst çerçeve

Kuantum avantajı iddia edilmiyor. 8 kübitlik hesap klasik bilgisayarda da tam çözülebiliyor. TD-DFT'ye karşı
doğruluk kazancı çok-konfigürasyonlu yaklaşımdan (CAS + NEVPT2) geliyor. Ayrıntı için teknik raporun 1.4. bölümüne bakın.
