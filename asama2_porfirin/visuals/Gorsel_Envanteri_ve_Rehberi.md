# Görsel Envanteri ve Rehberi

Bu dosya arşivdeki tüm görselleri standart adlarla listeler. Her görsel için şunlar verilmiştir: (a) neyi gösterdiği ve bilimsel açıklaması, (b) üretildiği araç ve görüntüleme parametreleri, (c) rapor veya sunumda kullanılacak hazır şekil açıklaması (caption).

**Klasörler:**
- `visuals/grafikler/`: Python/matplotlib grafikleri (fig01–fig20).
- `visuals/chimerax/`: UCSF ChimeraX ile üretilmiş 3B moleküler görseller (fig21 ve sonrası; 3000×2000 piksel, supersample 4).
- `visuals/chimerax_betikleri/`: Görselleri yeniden üreten `.cxc` betikleri. Dosya adlarındaki `__` işareti klasör ayırıcısıdır.
- `visuals/sunum_kirpilmis/`: Sunum için kırpılmış 1200×1200 kopyalar.

## Ortak görüntüleme parametreleri (ChimeraX)

| Görsel türü | Stil | Renk kodları | Işık ve kalite | Diğer |
|---|---|---|---|---|
| Uzay dolgu (molekül) | `style sphere` | `color byelement`; C #8a8a8a, N mavi, O kırmızı, H beyaz, metal element rengi | `lighting soft`, `shadows true intensity 0.6`, `multiShadow 1024`, `material specular 0.4 exponent 40`, `graphics quality 6` | `graphics silhouettes true width 2 color #2a2a2a`; `turn x -35` (ilaçlarda -40), `turn z 15`; beyaz arka plan |
| Top-çubuk | `style ball`, `size stickRadius 0.14`, `size ballScale 0.25` | aynı | aynı | aynı |
| Orbital (HOMO/LUMO) | `style ball`, `stickRadius 0.11`, `ballScale 0.22` (ilaçlarda 0.13 / 0.22) | C #909090, N #3050f8, H #f4f4f4; izoyüzey: negatif faz #2f6fd0 (mavi), pozitif faz #f0782f (turuncu) | `lighting full`, `shadows 0.5`, `multiShadow 512`, `graphics quality 5`, `camera ortho` | `volume ... style surface level ±0.02 transparency 0.1` (ilaçlarda 0.25), `surfaceSmoothing true smoothingIterations 4`; küp ızgarası 0,1 Å, kenar 7 bohr (PySCF cubegen) |
| Fark yoğunluğu (S1−S0) | Orbital ile aynı | turuncu (#f0782f) = elektronun geldiği yer (+), mavi (#2f6fd0) = ayrıldığı yer (−) | aynı | Seviye otomatik: ±0,12 × max\|Δρ\| (porfinde ±0,00058) |

Görsellerin üst ve alt kenarlarındaki etiketler `2dlabels` komutuyla eklendi. Etiketlerde molekül adı, S1 (eV ve nm), f ve triyaj kararı bulunur.

## A. Python/matplotlib grafikleri (`visuals/grafikler/`)

### `fig01_PMF_bCD_siklopentanol.png`
- **(a) Ne gösteriyor:** β-CD + siklopentanol için ortalama kuvvet potansiyeli (PMF) W(ξ). Alt panelde 37 şemsiye penceresinin histogramları (örtüşme kontrolü). Bağlı durum ξ ≈ −0,05 nm'de −17,7 kJ/mol.
- **(b) Araç/parametre:** GROMACS 2026.1 `gmx wham` (-b 200) → `plot_pmf.py` (matplotlib). Mavi #2a78d6 çizgi, kesikli çizgiler su platoları.
- **(c) Şekil açıklaması:** Şekil 1. β-siklodekstrin–siklopentanol kompleksinin CD ekseni boyunca ortalama kuvvet potansiyeli (37 pencere × 2 ns, 300 K, WHAM). Hesaplanan ΔG° = −2,7 ± 0,6 kcal/mol; deney −3,05 kcal/mol.

### `fig02_porfin_gouterman_orbitaller.png`
- **(a) Ne gösteriyor:** Porfinin Gouterman dört orbitali (HOMO-1, HOMO, LUMO, LUMO+1) ve enerjileri.
- **(b) Araç/parametre:** PySCF B3LYP/def2-SVP → molden → matplotlib/ilk görselleştirme (porfin_qm.py / ciz.py).
- **(c) Şekil açıklaması:** Şekil 2. Porfinin Gouterman orbitalleri (B3LYP/def2-SVP). HOMO-1/HOMO (−5,53/−5,40 eV) ve LUMO/LUMO+1 (−2,46/−2,40 eV) neredeyse dejenere; bu, güçlü konfigürasyon karışımının kaynağıdır.

### `fig03_porfin_orbital_secimi.png`
- **(a) Ne gösteriyor:** Porfinin 406 orbitallik enerji merdiveni (sol) ve seçilen 4 orbitallik kuantum bölgesi (sağ, yakınlaştırma).
- **(b) Araç/parametre:** matplotlib (`merdiven.py`). Turuncu = boş orbitaller, mavi = dolu orbitaller.
- **(c) Şekil açıklaması:** Şekil 3. Kuantum bölgesinin seçimi: 406 orbitalden ışığı belirleyen 4 Gouterman orbitali (8 kübit) kuantum bilgisayara, kalan 402 orbital klasik NEVPT2'ye verilir.

### `fig04_tarama_v3_sonuclari.png`
- **(a) Ne gösteriyor:** 8 porfirin türevi için hesap–deney karşılaştırması (dumbbell), osilatör gücü f ve triyaj kararı.
- **(b) Araç/parametre:** matplotlib (tarama sonuç betiği), dataviz paleti (mavi #2a78d6, turuncu #eb6834).
- **(c) Şekil açıklaması:** Şekil 4. Kuantum triyaj taraması (v3): CAS(4,4)-VQD+NEVPT2 ile S1 enerjileri, deney değerleri ve triyaj sınıflandırması.

### `fig05_porfin_panel_2_3.png`
- **(a) Ne gösteriyor:** Porfin için enerji/soğurma paneli ve triyaj yorumu (sunum paneli 2–3).
- **(b) Araç/parametre:** matplotlib (`panel.py`).
- **(c) Şekil açıklaması:** Şekil 5. Porfin: kuantum bölgesi sonuçları (S1 = 1,87 eV, f = 0,033) ve triyaj kararı (kuantum gerekli).

### `fig06_F1_yontem_karsilastirma.png`
- **(a) Ne gösteriyor:** Sol: 8 molekülde deney, bizim yöntem ve TD-DFT B3LYP S1 değerleri. Sağ: ortak moleküllerde (porfin, klorin) ortalama mutlak hata (0,07 vs 0,39 eV).
- **(b) Araç/parametre:** `final_figurler.py` (matplotlib, 220 dpi). Siyah elmas = deney, mavi daire = bizim yöntem, turuncu kare = TD-DFT.
- **(c) Şekil açıklaması:** Şekil 6. Kuantum-hibrit yöntem (VQE 8 kübit + NEVPT2) ile sektör standardı TD-DFT B3LYP'nin deneye göre karşılaştırması. Porfin ve klorinde ortalama mutlak hata: 0,07 eV ve 0,39 eV.

### `fig07_F2_triyaj_dogrulama.png`
- **(a) Ne gösteriyor:** x ekseni: S1'de ikinci konfigürasyonun ağırlığı (triyaj ölçütü); y ekseni: deneye göre mutlak hata. Eşik 0,2.
- **(b) Araç/parametre:** `final_figurler.py`.
- **(c) Şekil açıklaması:** Şekil 7. Triyaj ölçütü ile hesap hatası arasındaki ilişki. Eşiğin (0,2) sağı 'kuantum gerekli' bölgesidir.

### `fig08_F3_klinik_ilaclar.png`
- **(a) Ne gösteriyor:** Dört klinik PDT ilacı için klinikte kullanılan ışık (siyah), bizim çekirdek hesabı (mavi) ve TD-DFT çekirdek hesabı (turuncu); 700–900 nm doku geçirgenlik penceresi.
- **(b) Araç/parametre:** `final_figurler.py`; x ekseni 500–800 nm.
- **(c) Şekil açıklaması:** Şekil 8. Klinik ilaçların fotoaktif çekirdekleri: bizim hesabımız klinik ışığa 26–38 nm, TD-DFT 107–141 nm uzaklıkta. Yakın kızılötesi ilaç padeliporfinde (753 nm) 4 orbitalli model yetersiz.

### `fig09_ibm_kingston_donanim.png`
- **(a) Ne gösteriyor:** Porfin S1 için tam çözüm, ideal VQE, gürültü simülasyonu ve gerçek IBM Kingston donanım sonuçları (ham ve düzeltilmiş).
- **(b) Araç/parametre:** matplotlib (build_archive.py); veriler `data/ibm_donanim_porfin.csv`.
- **(c) Şekil açıklaması:** Şekil 9. VQE devresinin IBM Kingston (156 kübit) gerçek kuantum işlemcisindeki ilk ölçümü. Ham sonuç 1,14 eV, gürültü düzeltmeli 4,55 ± 0,19 eV; tam çözüm 3,62 eV. Bugünkü donanımda gürültü baskındır.

### `fig10_panel_Bakteriyoklorin_panel_2_3.png`
- **(a) Ne gösteriyor:** Bakteriyoklorin için enerji, soğurma şiddeti (f) ve triyaj paneli.
- **(b) Araç/parametre:** matplotlib (`paneller.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Bakteriyoklorin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig11_panel_ILAC_Protoporfirin_IX_panel.png`
- **(a) Ne gösteriyor:** Protoporfirin_IX için enerji, soğurma şiddeti (f) ve triyaj paneli (klinik ilaç paneli; eski sürüm — PpIX panelindeki ~1 eV klasik hata ifadesi sonradan geri çekildi, sunumda kullanmayın).
- **(b) Araç/parametre:** matplotlib (`ilac_panel.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Protoporfirin_IX: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig12_panel_ILAC_Talaporfin_panel.png`
- **(a) Ne gösteriyor:** Talaporfin için enerji, soğurma şiddeti (f) ve triyaj paneli (klinik ilaç paneli; eski sürüm — PpIX panelindeki ~1 eV klasik hata ifadesi sonradan geri çekildi, sunumda kullanmayın).
- **(b) Araç/parametre:** matplotlib (`ilac_panel.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Talaporfin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig13_panel_ILAC_Temoporfin_panel.png`
- **(a) Ne gösteriyor:** Temoporfin için enerji, soğurma şiddeti (f) ve triyaj paneli (klinik ilaç paneli; eski sürüm — PpIX panelindeki ~1 eV klasik hata ifadesi sonradan geri çekildi, sunumda kullanmayın).
- **(b) Araç/parametre:** matplotlib (`ilac_panel.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Temoporfin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig14_panel_Klorin_panel_2_3.png`
- **(a) Ne gösteriyor:** Klorin için enerji, soğurma şiddeti (f) ve triyaj paneli.
- **(b) Araç/parametre:** matplotlib (`paneller.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Klorin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig15_panel_Mg-porfin_panel_2_3.png`
- **(a) Ne gösteriyor:** Mg-porfin için enerji, soğurma şiddeti (f) ve triyaj paneli.
- **(b) Araç/parametre:** matplotlib (`paneller.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Mg-porfin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig16_panel_Pd-porfin_panel_2_3.png`
- **(a) Ne gösteriyor:** Pd-porfin için enerji, soğurma şiddeti (f) ve triyaj paneli.
- **(b) Araç/parametre:** matplotlib (`paneller.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Pd-porfin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig17_panel_Porfin_panel_2_3.png`
- **(a) Ne gösteriyor:** Porfin için enerji, soğurma şiddeti (f) ve triyaj paneli.
- **(b) Araç/parametre:** matplotlib (`paneller.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Porfin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig18_panel_Porfirazin_panel_2_3.png`
- **(a) Ne gösteriyor:** Porfirazin için enerji, soğurma şiddeti (f) ve triyaj paneli.
- **(b) Araç/parametre:** matplotlib (`paneller.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Porfirazin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig19_panel_Tetrametilporfin_panel_2_3.png`
- **(a) Ne gösteriyor:** Tetrametilporfin için enerji, soğurma şiddeti (f) ve triyaj paneli.
- **(b) Araç/parametre:** matplotlib (`paneller.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Tetrametilporfin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

### `fig20_panel_Zn-porfin_panel_2_3.png`
- **(a) Ne gösteriyor:** Zn-porfin için enerji, soğurma şiddeti (f) ve triyaj paneli.
- **(b) Araç/parametre:** matplotlib (`paneller.py`), dataviz paleti.
- **(c) Şekil açıklaması:** Zn-porfin: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).

## B. ChimeraX 3B moleküler görseller (`visuals/chimerax/`)

### `fig21_molekul_porfin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Porfin molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı (C20H14N4, 38 atom); atomlar van der Waals yarıçaplarıyla gösterilmiş.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Porfin (C20H14N4, 38 atom) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).

### `fig22_molekul_porfin_topcubuk.png`
- **(a) Ne gösteriyor:** Porfin molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.
- **(b) Araç/parametre:** UCSF ChimeraX; top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Porfin — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).

### `fig23_molekul_porfin_HOMO.png`
- **(a) Ne gösteriyor:** Porfin HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Porfin HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = 1.87 eV (663 nm; deney 1.98 eV), f = 0.033, triyaj: KUANTUM GEREKLİ.

### `fig24_molekul_porfin_LUMO.png`
- **(a) Ne gösteriyor:** Porfin LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Porfin LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları 0.544 / 0.435.

### `fig25_molekul_porfin_S1_S0_fark_yogunlugu.png`
- **(a) Ne gösteriyor:** Porfin için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler.
- **(b) Araç/parametre:** UCSF ChimeraX; S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Porfin S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = 1.87 eV, f = 0.033.

### `fig26_molekul_Zn-porfin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Zn-porfin molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı (C20H12N4Zn, 37 atom); atomlar van der Waals yarıçaplarıyla gösterilmiş.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Zn-porfin (C20H12N4Zn, 37 atom) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).

### `fig27_molekul_Zn-porfin_topcubuk.png`
- **(a) Ne gösteriyor:** Zn-porfin molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.
- **(b) Araç/parametre:** UCSF ChimeraX; top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Zn-porfin — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).

### `fig28_molekul_Zn-porfin_HOMO.png`
- **(a) Ne gösteriyor:** Zn-porfin HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Zn-porfin HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = 2.07 eV (599 nm; deney 2.18 eV), f = 0.079, triyaj: KUANTUM GEREKLİ.

### `fig29_molekul_Zn-porfin_LUMO.png`
- **(a) Ne gösteriyor:** Zn-porfin LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Zn-porfin LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları 0.616 / 0.367.

### `fig30_molekul_Zn-porfin_S1_S0_fark_yogunlugu.png`
- **(a) Ne gösteriyor:** Zn-porfin için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler.
- **(b) Araç/parametre:** UCSF ChimeraX; S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Zn-porfin S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = 2.07 eV, f = 0.079.

### `fig31_molekul_Mg-porfin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Mg-porfin molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı (C20H12N4Mg, 37 atom); atomlar van der Waals yarıçaplarıyla gösterilmiş.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Mg-porfin (C20H12N4Mg, 37 atom) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).

### `fig32_molekul_Mg-porfin_topcubuk.png`
- **(a) Ne gösteriyor:** Mg-porfin molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.
- **(b) Araç/parametre:** UCSF ChimeraX; top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Mg-porfin — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).

### `fig33_molekul_Mg-porfin_HOMO.png`
- **(a) Ne gösteriyor:** Mg-porfin HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Mg-porfin HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = 2.04 eV (608 nm; deney 2.17 eV), f = 0.073, triyaj: KUANTUM GEREKLİ.

### `fig34_molekul_Mg-porfin_LUMO.png`
- **(a) Ne gösteriyor:** Mg-porfin LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Mg-porfin LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları 0.611 / 0.373.

### `fig35_molekul_Mg-porfin_S1_S0_fark_yogunlugu.png`
- **(a) Ne gösteriyor:** Mg-porfin için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler.
- **(b) Araç/parametre:** UCSF ChimeraX; S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Mg-porfin S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = 2.04 eV, f = 0.073.

### `fig36_molekul_Pd-porfin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Pd-porfin molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı (C20H12N4Pd, 37 atom); atomlar van der Waals yarıçaplarıyla gösterilmiş.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Pd-porfin (C20H12N4Pd, 37 atom) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).

### `fig37_molekul_Pd-porfin_topcubuk.png`
- **(a) Ne gösteriyor:** Pd-porfin molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.
- **(b) Araç/parametre:** UCSF ChimeraX; top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Pd-porfin — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).

### `fig38_molekul_Pd-porfin_HOMO.png`
- **(a) Ne gösteriyor:** Pd-porfin HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Pd-porfin HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = 2.14 eV (579 nm; deney değeri yok), f = 0.085, triyaj: KUANTUM GEREKLİ.

### `fig39_molekul_Pd-porfin_LUMO.png`
- **(a) Ne gösteriyor:** Pd-porfin LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Pd-porfin LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları 0.629 / 0.355.

### `fig40_molekul_Pd-porfin_S1_S0_fark_yogunlugu.png`
- **(a) Ne gösteriyor:** Pd-porfin için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler.
- **(b) Araç/parametre:** UCSF ChimeraX; S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Pd-porfin S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = 2.14 eV, f = 0.085.

### `fig41_molekul_tetrametilporfin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Tetrametilporfin molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı (C24H22N4, 50 atom); atomlar van der Waals yarıçaplarıyla gösterilmiş.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Tetrametilporfin (C24H22N4, 50 atom) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).

### `fig42_molekul_tetrametilporfin_topcubuk.png`
- **(a) Ne gösteriyor:** Tetrametilporfin molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.
- **(b) Araç/parametre:** UCSF ChimeraX; top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Tetrametilporfin — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).

### `fig43_molekul_tetrametilporfin_HOMO.png`
- **(a) Ne gösteriyor:** Tetrametilporfin HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Tetrametilporfin HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = 1.64 eV (756 nm; deney 1.88 eV), f = 0.001, triyaj: KUANTUM GEREKLİ.

### `fig44_molekul_tetrametilporfin_LUMO.png`
- **(a) Ne gösteriyor:** Tetrametilporfin LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Tetrametilporfin LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları 0.507 / 0.47.

### `fig45_molekul_tetrametilporfin_S1_S0_fark_yogunlugu.png`
- **(a) Ne gösteriyor:** Tetrametilporfin için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler.
- **(b) Araç/parametre:** UCSF ChimeraX; S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Tetrametilporfin S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = 1.64 eV, f = 0.001.

### `fig46_molekul_klorin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Klorin molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı (C20H16N4, 40 atom); atomlar van der Waals yarıçaplarıyla gösterilmiş.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Klorin (C20H16N4, 40 atom) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).

### `fig47_molekul_klorin_topcubuk.png`
- **(a) Ne gösteriyor:** Klorin molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.
- **(b) Araç/parametre:** UCSF ChimeraX; top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Klorin — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).

### `fig48_molekul_klorin_HOMO.png`
- **(a) Ne gösteriyor:** Klorin HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Klorin HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = 1.98 eV (626 nm; deney 1.95 eV), f = 0.579, triyaj: klasik yeterli.

### `fig49_molekul_klorin_LUMO.png`
- **(a) Ne gösteriyor:** Klorin LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Klorin LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları 0.743 / 0.098.

### `fig50_molekul_klorin_S1_S0_fark_yogunlugu.png`
- **(a) Ne gösteriyor:** Klorin için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler.
- **(b) Araç/parametre:** UCSF ChimeraX; S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Klorin S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = 1.98 eV, f = 0.579.

### `fig51_molekul_porfirazin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Porfirazin molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı (C16H10N8, 34 atom); atomlar van der Waals yarıçaplarıyla gösterilmiş.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Porfirazin (C16H10N8, 34 atom) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).

### `fig52_molekul_porfirazin_topcubuk.png`
- **(a) Ne gösteriyor:** Porfirazin molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.
- **(b) Araç/parametre:** UCSF ChimeraX; top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Porfirazin — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).

### `fig53_molekul_porfirazin_HOMO.png`
- **(a) Ne gösteriyor:** Porfirazin HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Porfirazin HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = 2.24 eV (553 nm; deney 2.065 eV), f = 0.633, triyaj: klasik yeterli.

### `fig54_molekul_porfirazin_LUMO.png`
- **(a) Ne gösteriyor:** Porfirazin LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Porfirazin LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları 0.788 / 0.139.

### `fig55_molekul_porfirazin_S1_S0_fark_yogunlugu.png`
- **(a) Ne gösteriyor:** Porfirazin için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler.
- **(b) Araç/parametre:** UCSF ChimeraX; S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Porfirazin S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = 2.24 eV, f = 0.633.

### `fig56_molekul_bakteriyoklorin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Bakteriyoklorin molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı (C20H18N4, 42 atom); atomlar van der Waals yarıçaplarıyla gösterilmiş.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Bakteriyoklorin (C20H18N4, 42 atom) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).

### `fig57_molekul_bakteriyoklorin_topcubuk.png`
- **(a) Ne gösteriyor:** Bakteriyoklorin molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.
- **(b) Araç/parametre:** UCSF ChimeraX; top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (molekul_hepsi.cxc ile üretildi).
- **(c) Şekil açıklaması:** Bakteriyoklorin — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).

### `fig58_molekul_bakteriyoklorin_HOMO.png`
- **(a) Ne gösteriyor:** Bakteriyoklorin HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Bakteriyoklorin HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = 2.23 eV (556 nm; deney 1.74 eV), f = 0.422, triyaj: klasik yeterli (model yetersiz).

### `fig59_molekul_bakteriyoklorin_LUMO.png`
- **(a) Ne gösteriyor:** Bakteriyoklorin LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Bakteriyoklorin LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları 0.748 / 0.168.

### `fig60_molekul_bakteriyoklorin_S1_S0_fark_yogunlugu.png`
- **(a) Ne gösteriyor:** Bakteriyoklorin için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler.
- **(b) Araç/parametre:** UCSF ChimeraX; S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (gorsel_hq.py ile üretildi).
- **(c) Şekil açıklaması:** Bakteriyoklorin S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = 2.23 eV, f = 0.422.

### `fig61_ilac_Temoporfin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Klinik ilaç Temoporfin (mTHPC, Foscan®) (C44H32N4O4, 84 atom); RDKit ile 3B'ye gömülüp klorin çekirdeğimize Kabsch ile hizalandı (RMSD 0.4 Å). Uzay dolgu.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Temoporfin (mTHPC, Foscan®) — AB'de onaylı; baş-boyun kanseri PDT; 652 nm ışık. Yapı RDKit + MMFF, klorin çekirdeğine hizalanmış (RMSD 0.4 Å).

### `fig62_ilac_Temoporfin_HOMO_kuantum_bolgesi.png`
- **(a) Ne gösteriyor:** Temoporfin (mTHPC, Foscan®) üzerinde, klorin çekirdeğinden hesaplanan HOMO izoyüzeyi. Renkli bulutlar kuantum bölgesini (8 kübit) gerçek ilacın fotoaktif halkası üzerinde gösterir.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Temoporfin (mTHPC, Foscan®) ve hesapladığımız kuantum bölgesi (HOMO, klorin çekirdeği, B3LYP/def2-SVP, ±0,02 a.u.). Yan gruplar hesaba dahil değildir.

### `fig63_ilac_Temoporfin_LUMO_kuantum_bolgesi.png`
- **(a) Ne gösteriyor:** Temoporfin (mTHPC, Foscan®) üzerinde klorin çekirdeğinin LUMO'su: ışıkla elektronun gittiği orbital.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Temoporfin (mTHPC, Foscan®) — klorin çekirdeğinin LUMO'su (±0,02 a.u.).

### `fig64_ilac_Talaporfin_uzaydolgu.png`
- **(a) Ne gösteriyor:** Klinik ilaç Talaporfin sodyum (NPe6, Laserphyrin®) (93 atom (nötr asit formu)); RDKit ile 3B'ye gömülüp klorin çekirdeğimize Kabsch ile hizalandı (RMSD 0.24 Å). Uzay dolgu.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Talaporfin sodyum (NPe6, Laserphyrin®) — Japonya'da onaylı; akciğer kanseri PDT; 664 nm ışık. Yapı RDKit + MMFF, klorin çekirdeğine hizalanmış (RMSD 0.24 Å).

### `fig65_ilac_Talaporfin_HOMO_kuantum_bolgesi.png`
- **(a) Ne gösteriyor:** Talaporfin sodyum (NPe6, Laserphyrin®) üzerinde, klorin çekirdeğinden hesaplanan HOMO izoyüzeyi. Renkli bulutlar kuantum bölgesini (8 kübit) gerçek ilacın fotoaktif halkası üzerinde gösterir.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Talaporfin sodyum (NPe6, Laserphyrin®) ve hesapladığımız kuantum bölgesi (HOMO, klorin çekirdeği, B3LYP/def2-SVP, ±0,02 a.u.). Yan gruplar hesaba dahil değildir.

### `fig66_ilac_Talaporfin_LUMO_kuantum_bolgesi.png`
- **(a) Ne gösteriyor:** Talaporfin sodyum (NPe6, Laserphyrin®) üzerinde klorin çekirdeğinin LUMO'su: ışıkla elektronun gittiği orbital.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Talaporfin sodyum (NPe6, Laserphyrin®) — klorin çekirdeğinin LUMO'su (±0,02 a.u.).

### `fig67_ilac_Protoporfirin-IX_uzaydolgu.png`
- **(a) Ne gösteriyor:** Klinik ilaç Protoporfirin IX (PpIX) (C34H34N4O4, 76 atom); RDKit ile 3B'ye gömülüp porfin çekirdeğimize Kabsch ile hizalandı (RMSD 0.22 Å). Uzay dolgu.
- **(b) Araç/parametre:** UCSF ChimeraX; uzay dolgu (van der Waals küreleri); ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Protoporfirin IX (PpIX) — 5-ALA/MAL ön ilacından vücutta oluşur; Ameluz® + 635 nm ışık. Yapı RDKit + MMFF, porfin çekirdeğine hizalanmış (RMSD 0.22 Å).

### `fig68_ilac_Protoporfirin-IX_HOMO_kuantum_bolgesi.png`
- **(a) Ne gösteriyor:** Protoporfirin IX (PpIX) üzerinde, porfin çekirdeğinden hesaplanan HOMO izoyüzeyi. Renkli bulutlar kuantum bölgesini (8 kübit) gerçek ilacın fotoaktif halkası üzerinde gösterir.
- **(b) Araç/parametre:** UCSF ChimeraX; HOMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Protoporfirin IX (PpIX) ve hesapladığımız kuantum bölgesi (HOMO, porfin çekirdeği, B3LYP/def2-SVP, ±0,02 a.u.). Yan gruplar hesaba dahil değildir.

### `fig69_ilac_Protoporfirin-IX_LUMO_kuantum_bolgesi.png`
- **(a) Ne gösteriyor:** Protoporfirin IX (PpIX) üzerinde porfin çekirdeğinin LUMO'su: ışıkla elektronun gittiği orbital.
- **(b) Araç/parametre:** UCSF ChimeraX; LUMO izoyüzeyi (±0,02) + top-çubuk iskelet; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında (ilaclar_hepsi.cxc / gercek_ilac.py ile üretildi).
- **(c) Şekil açıklaması:** Protoporfirin IX (PpIX) — porfin çekirdeğinin LUMO'su (±0,02 a.u.).

## C. Sohbette kullanılan ama arşive dahil edilmeyen dış görseller

- Ekran görüntüleri: terminal çıktıları, IBM Quantum Platform ekranları, Avogadro ekranları. Bilgi amaçlıdır; sonuçları `data/` klasöründe sayısal olarak yer alır.
- RCSB 4V6X (ribozom) örnek görseli: yalnızca görsel kalite tartışması için kullanıldı.
- Elektrostatik yüzeyli (Gasteiger + coulombic) Temoporfin denemesi: soluk ve beyaz çıktığı için kullanılmadı. Yeniden üretmek için: `addcharge #1 method gasteiger ; surface #1 ; coulombic #1 range -2,2 ; transparency #1 0 target s`.
