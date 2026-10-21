# Kuantum Triyaj — Proje Arşivi ve Master Teknik Rapor

**Proje:** Molecule / Kuantum Triyaj (PDT fotosensitizörleri için kuantum-hibrit hesaplama platformu) **Hazırlayan:** Sinan Ahmet Yalçınkaya (Yıldız Teknik Üniversitesi), Claude ile birlikte **Arşiv tarihi:** 20 Ekim 2026 **Durum:** YTÜ Startup House başvurusu için prototip tamamlandı. TD-DFT referans hesabının son iki molekülü (tetrametilporfin, Pd-porfin) kullanıcının bilgisayarında çalışmaya devam ediyor.

---

## 1\. Proje özeti ve amaç

### 1.1 Ne yaptık (kısa)

Proje iki aşamada ilerledi:

1. **β-siklodekstrin (β-CD) \+ siklopentanol konak–misafir kompleksi.** Klasik moleküler dinamik (GROMACS) ve şemsiye örnekleme (umbrella sampling) ile bağlanma serbest enerjisi hesaplandı ve deneyle doğrulandı. Ardından aynı sistemde kuantum iş akışı kuruldu: PySCF → aktif uzay → Jordan–Wigner → VQE. Bu iş akışı doğrulandı ve parçalama (fragmentasyon) yaklaşımı (MFCC) test edildi. Sonuç: bu sistemde bağlanma dispersiyon güdümlü; kuantum düzeltmesi pratikte sıfır, yani **"kuantum gereksiz"**.  
2. **Fotodinamik terapi (PDT) fotosensitizörleri: porfirin ailesi.** "Kuantum hangi molekülde gerçekten gerekli?" sorusunu yanıtlamak için, ışığı soğuran bölgeyi (Gouterman'ın 4 orbitali \= 8 kübit) VQE ile çözen ve kalan elektron korelasyonunu klasik NEVPT2 ile ekleyen hibrit bir iş akışı kuruldu. Bu iş akışı 8 porfirin türevine uygulandı ve deneyle karşılaştırıldı. Ardından üç klinik ilacın (Temoporfin, Talaporfin, Protoporfirin IX) gerçek yapıları bu çekirdeklerin üzerine oturtuldu. Sektör standardı TD-DFT ile karşılaştırma yapıldı ve VQE devresi IBM'in 156 kübitlik gerçek kuantum bilgisayarında (ibm\_kingston) çalıştırıldı.

### 1.2 Temel hipotez

> Kuantum hesaplama kaynağı kıt, pahalı ve gürültülüdür. Her molekülü kuantumla hesaplamak yerine, bir **triyaj** adımı molekülün elektron yapısının çok-konfigürasyonlu olup olmadığını ölçer. Kuantum kaynağını yalnızca gereken moleküle ve yalnızca kritik bölgeye (aktif uzay) harcar; geri kalanını klasik yöntemle hesaplar.

Triyaj kuralı: Kuantum bölgesinde S1 uyarılmış halinde **ikinci elektron konfigürasyonunun ağırlığı \> 0,2 ise "KUANTUM GEREKLİ"**, değilse "klasik yeterli".

### 1.3 Temel çıktılar

| Çıktı | Sonuç |
| :---- | :---- |
| β-CD PMF bağlanma enerjisi | ΔG° \= −2,7 ± 0,6 kcal/mol (deney −3,05) |
| MFCC parçalama hatası | −0,18 ± 0,13 kcal/mol (%2) |
| β-CD kuantum düzeltmesi | \+0,0009 kcal/mol → kuantum gereksiz (triyaj doğru) |
| VQE \= CASCI doğruluğu | 10⁻⁸–10⁻¹¹ mertebesinde fark (simülatörde) |
| Porfirin S1, deneye ortalama mutlak hata | **0,11 eV** (doğrudan deneyli 5 molekül), 0,18 eV (7 molekül) |
| Porfin \+ klorin: TD-DFT vs biz | TD-DFT B3LYP 0,39 eV, biz **0,07 eV** |
| Klinik ilaç çekirdekleri | Bizim hesap klinik ışığa 26–38 nm, TD-DFT 107–141 nm uzak |
| IBM gerçek donanım (porfin S1, 8 kübit) | Ham −2,48 eV, gürültü düzeltmeli \+0,93 ± 0,19 eV (tam çözüme göre) |
| Startup sunumu | 17 slaytlık "Kuantum Triyaj – YTÜ Startup House Başvuru Sunumu" (Claude Slides artifact) |

&nbsp;

### 1.4 Dürüst çerçeve (iddia etmediklerimiz)

- **Kuantum avantajı iddia edilmiyor.** 8 kübitlik hesap klasik bilgisayarda da tam çözülüyor; TD-DFT'ye karşı doğruluk üstünlüğü çok-konfigürasyonlu yaklaşımdan (CAS \+ NEVPT2) geliyor, NEVPT2 de klasik bir yöntem.  
- **İlacın tamamı hesaplanmadı.** Yalnızca fotoaktif çekirdek hesaplandı; yan gruplar ve çözücü yok.  
- **NIR bölgesi zayıf nokta.** Bakteriyoklorinde (padeliporfin/Tookad çekirdeği) 4 orbitalli model 0,49 eV hata veriyor; CAS(8,8) \= 16 kübit gerekli (betik hazır: `src/06_referans_v4/cas_buyuk.py`).  
- **T1 (triplet) güvenilir değil.** ΔSCF porfinde 2,00 eV (deney 1,58); TDA-B3LYP 1,83 eV ile daha iyi.  
- **Gerçek donanımda gürültü baskın.** Sinyalin yalnızca %41–53'ü korunuyor.

---

## 2\. Proje zaman çizelgesi (aşamalar)

| Aşama | İçerik | Ana dosyalar |
| :---- | :---- | :---- |
| A | β-CD \+ siklopentanol: topoloji (ATB), 10 ns MD, 37 pencereli umbrella sampling, WHAM | `src/01_bcd_md_analiz/`, `data/raw/` |
| B | Bağlanma modları (901 kare), temsilci kareler | `temsilci_kare.py`, `molecules/bcd/` |
| C | Kuantum doğrulama: misafir VQE (8 kübit) \= CASCI | `src/02_bcd_kuantum/vqe_misafir.py` |
| D | GFN2-xTB tam kompleks etkileşim, D4 dispersiyon | `xtb_etkilesim.py` |
| E | MFCC parçalama, bölgede HF/MP2/B3LYP-D3, bölge VQE \+ hibrit | `fragman_xtb.py`, `bolge_qm.py`, `bolge_vqe.py` |
| F | Porfin: Gouterman CAS(4,4) VQD \+ NEVPT2, literatür doğrulaması | `src/03_porfin_kuantum/` |
| G | 8 porfirin taraması (v1 → v2 → v3; v3 nihai) | `src/04_tarama_v3/` |
| H | Görselleştirme: cube/molden → ChimeraX; gerçek ilaçların hizalanması | `src/05_gorsellestirme/` |
| I | Literatür (Xanadu, Algorithmiq) ve klinik dalga boylarının doğrulanması | `docs/` |
| J | Klasik referans: TD-DFT B3LYP; büyük aktif uzay betiği | `src/06_referans_v4/` |
| K | IBM Quantum: kompakt ADAPT devresi, ibm\_kingston'da ölçüm | `src/07_ibm_v5/` |
| L | YTÜ Startup House sunumu (17 slayt) | `docs/Sunum_metni_ve_konusmaci_notlari.md` |

&nbsp;

---

## 3\. Yöntem ve iş akışı

```
                 SMILES / yapı
                      │
          RDKit (3B gömme, MMFF) ──► GFN2-xTB (tblite + ASE BFGS) ──► geometri (.xyz)
                      │
          PySCF B3LYP/def2-SVP (DF) ──► molekül orbitalleri, HOMO–LUMO, metal payı
                      │
          ┌───────────┴─────────────┐
   Kuantum bölgesi seçimi        Klasik referanslar
   (Gouterman 4 orbital          ├─ TD-DFT B3LYP (RPA singlet, TDA triplet)
    = CAS(4e,4o) = 8 kübit)      └─ ΔSCF T1 (UKS-B3LYP)
          │
   Jordan–Wigner → 361 Pauli terimli kübit Hamiltonyeni
          │
   VQE (S0) + VQD (S1, S2, T1)   [PennyLane, Givens uyarılma ansatzı, S² cezası, örtüşme cezası β=3]
          │
   Konfigürasyon ağırlıkları ──► TRİYAJ (2. ağırlık > 0,2 → kuantum gerekli)
          │
   Geçiş dipolleri (kübit dipol operatörü) ──► osilatör gücü f
          │
   SC-NEVPT2 (S0, S1) ──► dinamik korelasyon düzeltmesi ──► S1 (eV, nm)
          │
   Deney + literatür + klinik ışık ile karşılaştırma ──► figürler (matplotlib) + 3B görseller (ChimeraX)
          │
   (isteğe bağlı) IBM Quantum: kompakt ADAPT devresi → EstimatorV2 (TREX, twirling, DD) → gürültü düzeltmesi
```

&nbsp;

**Hibrit enerji formülü (β-CD bölge çalışması):** E \= HF \+ MP2(tüm MO) \+ \[VQE − MP2\]\_aktif uzay **Hibrit uyarılma (porfirin):** S1 \= \[E(S1) − E(S0)\]\_CASCI/VQD \+ \[ΔE\_NEVPT2(S1) − ΔE\_NEVPT2(S0)\] **Donanım sonucunun NEVPT2'ye eklenmesi:** S1\_donanım \= S1\_NEVPT2 \+ (ΔE\_donanım − ΔE\_CASCI)

---

## 4\. Moleküler görselleştirme ve yazılım altyapısı

### 4.1 Ortam

- **İşletim sistemi:** Windows 11 üzerinde WSL2 Ubuntu (16 thread, WSL belleği 7 GB, Windows RAM 15,6 GB). `/tmp` 3,8 GB'lık tmpfs olduğu için `PYSCF_TMPDIR=~/pyscf_tmp` zorunlu.  
- **Python ortamı:** `~/qenv` (Python 3.14 sanal ortamı).  
- **Bulut denemeleri:** Google Cloud ve Azure for Students kullanılmadı. Sonuç olarak tüm hesaplar kullanıcının bilgisayarında yapıldı.

### 4.2 Yazılımlar, kütüphaneler ve kullanım amaçları

| Araç | Sürüm | Kullanım amacı | Nerede |
| :---- | :---- | :---- | :---- |
| GROMACS | 2026.1 | MD, umbrella sampling, `gmx wham` | Aşama A |
| GROMOS 54A7 / ATB | molid 345755 (β-CD), 859 (siklopentanol) | Kuvvet alanı ve topoloji | Aşama A |
| PySCF | 2.14 | HF, DFT (B3LYP, CAM-B3LYP), MP2, CASCI, SC-NEVPT2, TD-DFT/TDA, cubegen, molden | Aşama C–K |
| PennyLane | 0.45.1 | Kübit Hamiltonyeni (`qml.qchem`), VQE/VQD, Givens kapıları, OpenQASM dışa aktarımı | C, E, F, G, K |
| tblite (GFN2-xTB) | 0.7 | Geometri optimizasyonu, etkileşim enerjileri (ALPB çözücü) | D, E, G |
| dftd3 / dftd4 | 1.6 / 4.3 | Dispersiyon düzeltmesi, D4 ayrıştırması | D, E |
| RDKit | 2026.9.1 | SMILES → 3B, MMFF, alt yapı eşleştirme (ilaç ↔ çekirdek) | G, H |
| ASE | 3.29 | BFGS optimizasyonu (xTB hesaplayıcısıyla) | G |
| Qiskit | 2.6.0 | Devre transpile, SparsePauliOp, QWC gruplama | K |
| qiskit-ibm-runtime | 0.50.0 | IBM Quantum hesabı, EstimatorV2, FakeKingston | K |
| qiskit-aer | 0.17.2 | Gürültülü simülasyon (sandbox'ta) | K |
| NumPy / SciPy | – | Doğrusal cebir, BFGS optimizasyonu | Tümü |
| matplotlib | – | Tüm grafikler (PMF, tarama, F1–F3, paneller) | Tümü |
| UCSF ChimeraX | masaüstü | 3B moleküler görseller (yüzey, top-çubuk, orbital izoyüzeyleri) | H |
| Avogadro 2 | masaüstü | İlk orbital görselleştirmeleri (molden → MO yüzeyleri) | F |
| IBM Quantum Platform | Open Plan | ibm\_kingston (156 kübit Heron r2) gerçek donanım | K |
| Claude Slides (artifact) | – | Başvuru sunumu | L |

&nbsp;

### 4.3 Kurulum

```sh
# WSL Ubuntu
sudo apt-get update -y && sudo apt-get install -y python3-venv python3-pip
python3 -m venv ~/qenv && source ~/qenv/bin/activate
pip install --upgrade pip
pip install numpy scipy pyscf pennylane tblite rdkit ase h5py matplotlib dftd3 dftd4
pip install "qiskit==2.6.0" "qiskit-ibm-runtime==0.50.0" "qiskit-aer==0.17.2"
mkdir -p ~/pyscf_tmp && echo 'export PYSCF_TMPDIR=~/pyscf_tmp' >> ~/.bashrc
```

&nbsp;

**Bilinen sorun:** pyscf, qiskit/pennylane'den sonra import edilirse `libgomp ... cannot allocate memory in static TLS block` hatası veriyor. Çözüm: pyscf'i en başta import etmek (`ibm_vqe.py` bu şekilde düzeltildi). qiskit\_aer kullanıcı ortamında aynı nedenle yüklenemiyor (gürültü simülasyonu ve `--sahte` modu yalnızca sandbox'ta çalıştı). Alternatif: `LD_PRELOAD=$(ls ~/qenv/lib/python3.14/site-packages/pyscf/lib/libgomp-*.so*)`.

### 4.4 Görselleştirme iş akışı

1. `gorsel_hq.py`: PySCF ile molden ve küp dosyaları üretir.  
   - Küpler: HOMO-1, HOMO, LUMO, LUMO+1 ve S1−S0 fark yoğunluğu (`_fark.cube`).  
   - Izgara aralığı 0,1 Å, kenar payı 7,0 bohr.  
   - Her görünüm için bir ChimeraX `.cxc` betiği yazar: HOMO, LUMO, fark, video.  
2. ChimeraX'te `open .../hepsi.cxc` veya `molekul_hepsi.cxc` çalıştırılır.  
   - Çıktı: 3000×2000 piksel, supersample 4\.  
3. `gercek_ilac.py`: Klinik ilacın SMILES'ı RDKit ile 3B'ye gömülür. Ardından:  
   - Jenerik çekirdek SMARTS'ı ile alt yapı eşleşmesi yapılır.  
   - İlaç, Kabsch hizalamasıyla bizim xTB çekirdeğimizin üzerine oturtulur.  
   - Sonuç PDB \+ CONECT olarak yazılır, ChimeraX betikleri üretilir.  
4. Ayrıntılı parametreler ve tüm şekil açıklamaları: `visuals/Gorsel_Envanteri_ve_Rehberi.md`.

---

## 5\. Sonuç özetleri

### 5.1 β-CD \+ siklopentanol

| Büyüklük | Değer |
| :---- | :---- |
| ΔG° (PMF/WHAM) | −2,7 ± 0,6 kcal/mol |
| ΔG° deney (ITC, Rekharsky 1997\) | −3,05 kcal/mol |
| Mod A / Mod B / geçiş | %53,5 / %38,7 / %7,8 |
| E\_int gaz / su (GFN2-xTB, 31 kare) | −9,34 ± 1,63 / −6,33 ± 1,45 kcal/mol |
| D4 dispersiyon payı | \~%96 |
| MFCC hatası | −0,18 ± 0,13 kcal/mol |
| Bölge (39 atom, CP): HF / B3LYP / MP2 / B3LYP-D3 / xTB | \+1,39 / \+0,73 / −1,32 / −2,83 / −2,80 kcal/mol |
| Kuantum düzeltmesi (VQE − MP2, aktif uzay) | \+0,0009 kcal/mol |

&nbsp;

### 5.2 Porfirin taraması (v3)

Tam tablo için bkz. `data/v3_tarama_sonuclari.csv` ve `data/Veriler_ve_Kuresel_Muadil_Analizi.md`. Özet:

- 5 doğrudan deneyli molekülde S1 ortalama mutlak hatası (MAE) **0,111 eV**; 7 molekülde 0,184 eV.  
- Triyaj:  
  - **Kuantum gerekli:** porfin, Zn-, Mg-, Pd-porfin, tetrametilporfin (ikinci konfigürasyon ağırlığı 0,36–0,47).  
  - **Klasik yeterli:** klorin, porfirazin, bakteriyoklorin (0,10–0,17).  
- Osilatör gücü klinik tasarım tercihini yeniden üretiyor: klorin f \= 0,579, porfin f \= 0,033. Klinik ilaçların klorin tabanlı olmasının nedeni güçlü kırmızı soğurma.

### 5.3 Klasik referans ve gerçek donanım

- TD-DFT B3LYP:  
  - porfin 2,348 eV (deney 1,98, \+0,37)  
  - klorin 2,369 eV (deney 1,95, \+0,42)  
  - Zn-, Mg-porfin, porfirazin, bakteriyoklorin kullanıcının bilgisayarında hesaplandı (JSON dosyaları `~/tarama`'da); tetrametilporfin ve Pd-porfin çalışıyor.  
- IBM ibm\_kingston, porfin:  
  - CASCI 3,617; ideal VQE 3,616  
  - ham 1,141; düzeltilmiş 4,546 ± 0,19 eV  
  - gürültü faktörü f \= 0,41 / 0,53  
  - İş kimliği: dbbl9ns24mrc73chvt8g

---

## 6\. Nasıl yeniden üretilir

```sh
source ~/qenv/bin/activate && cd ~/tarama
# 1) Geometri + kuantum triyaj taraması (8 molekül, ~20 dk/molekül, 16 thread)
bash calistir.sh 1            # geo.py + tarama_qm.py + ozet.py → *_v3_sonuc.json, tarama_sonuclari.csv
# 2) Görseller
python3 gorsel_hq.py porfin klorin ...    # molden + cube + .cxc → ChimeraX
python3 gercek_ilac.py                    # Temoporfin, Talaporfin, PpIX hizalama
# 3) Klasik referans + figürler
FONK=b3lyp nohup bash calistir_v4.sh > v4.log 2>&1 &
python3 final_figurler.py                 # F1, F2, F3, final_tablo.csv
# 4) Büyük aktif uzay (bakteriyoklorin)
python3 cas_buyuk.py bakteriyoklorin 8 8
# 5) IBM Quantum
python3 hesap_kaydet.py                   # API anahtarı (gizli giriş)
QMEM=2000 OMP_NUM_THREADS=4 python3 ibm_vqe.py hazirla porfin
python3 ibm_vqe.py gonder porfin          # onay sorar
python3 ibm_vqe.py sonuc porfin
```

&nbsp;

---

## 7\. Dizin rehberi

| Klasör | İçerik |
| :---- | :---- |
| `docs/` | Kaynakça, aşama özetleri, literatür karşılaştırması, IBM planı, sunum metni |
| `data/` | Tüm sayısal sonuçlar (CSV), veri yorumu ve küresel muadil analizi, ham çıktılar |
| `molecules/` | SMILES listesi, PDB (xTB geometrileri \+ hizalanmış ilaçlar), XYZ, SDF, molden, β-CD topolojileri |
| `visuals/` | Grafikler (yeniden çizilmiş), ChimeraX görselleri (standart adlarla), .cxc betikleri, görsel envanteri |
| `src/` | Tüm betikler, aşamalara göre numaralı klasörlerde, kısaltılmadan |

&nbsp;

## 8\. Açık işler

1. TD-DFT: tetrametilporfin ve Pd-porfin'in bitmesi; tüm `*_tddft.json` dosyalarının arşive eklenmesi ve F1'in 6–8 molekülle yeniden çizilmesi.  
2. Bakteriyoklorin için CAS(8,8)+NEVPT2 (16 kübit).  
3. Tetrametilporfin geometrisinin düzlemsel başlangıçla yeniden optimize edilmesi.  
4. T1 için TDA veya daha büyük aktif uzay.  
5. IBM: ZNE (`--zne`), daha sığ devre, kübit seçimi; kalan Open Plan süresinin kontrolü.  
6. Yan grup ve çözücü etkisi (tam ilaç hesabı).  
7. Müşteri görüşmeleri (5 PDT laboratuvarı), web demosu.