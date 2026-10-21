# Elde Edilen Veriler ve Küresel Muadil Analizi

Tüm sayılar bu projede hesaplanmıştır; dış kaynaklardan gelen sayılar ayrıca belirtilmiştir. Enerjiler eV, bağlanma enerjileri kcal/mol cinsindendir. CSV karşılıkları `data/` klasöründedir.

---

## A. β-siklodekstrin + siklopentanol (konak–misafir)

### A1. Serbest enerji ve bağlanma modları (`data/bcd_sonuclar.csv`)
| Büyüklük | Değer | Yorum |
|---|---|---|
| ΔG°(PMF, WHAM) | −2,7 ± 0,6 kcal/mol | Deneyle (−3,05) 0,33 kcal/mol fark. Literatürdeki RMSE aralığının (0,9–1,8) içinde. |
| Deney (ITC, Rekharsky 1997) | −3,05 kcal/mol (ΔH ≈ −1,1; TΔS ≈ +2,0) | Bağlanma entropi güdümlü. |
| Mod A (OH dar/primer kenarda) | %53,5 | Halka her iki modda da boşluk merkezinde. |
| Mod B (OH geniş/sekonder kenarda) | %38,7 | ΔG(A→B) ≈ +0,2 kcal/mol, yani modlar neredeyse eş enerjili. |
| Geçiş | %7,8 | |

### A2. Etkileşim enerjileri (GFN2-xTB, 31 MD karesi)
| Küme | E_int gaz | E_int su (ALPB) | D4 dispersiyon |
|---|---|---|---|
| Mod A (n=16) | −9,70 ± 1,62 | −6,78 ± 1,42 | −8,94 ± 0,86 |
| Mod B (n=13) | −8,67 ± 1,46 | −5,69 ± 1,35 | −9,06 ± 0,80 |
| Tümü (n=31) | −9,34 ± 1,63 | −6,33 ± 1,45 | −9,00 ± 0,84 |

**Yorum:** Gaz fazı etkileşiminin yaklaşık %96'sı dispersiyon; bağlanma dispersiyon güdümlü.

### A3. Parçalama (fragmentasyon) doğrulaması
| Şema | Toplam | Hata (tam −9,20'ye göre) |
|---|---|---|
| 7 tekli glukoz toplamı | −12,39 | −3,19 ± 0,22 (H kapakları çift sayılıyor) |
| **MFCC (örtüşen ikililer − tekliler)** | −9,38 | **−0,18 ± 0,13 (%2)** |

### A4. Bölgede ab initio (misafir + glukoz 2, 39 atom, def2-SVP, counterpoise)
| Yöntem | E_int (kcal/mol) |
|---|---|
| HF | +1,39 |
| B3LYP | +0,73 |
| MP2 | −1,32 |
| B3LYP-D3(BJ) | −2,83 |
| GFN2-xTB | −2,80 |

**Yorum:** Korelasyon içermeyen yöntemler itme, korelasyonlu yöntemler çekim veriyor (MP2 − HF = −2,71). Bağlanma bir elektron korelasyonu etkisi. B3LYP-D3 ≈ xTB çıktığı için xTB sonuçları bağımsız olarak doğrulanmış oldu.

### A5. Bölge VQE ve hibrit sonuç
- VQE korelasyon katkısı +0,0032 kcal/mol; CASCI ile farkı 2·10⁻⁸.
- Aynı aktif uzayda MP2 +0,0023. Kuantum düzeltmesi +0,0009 kcal/mol.
- Hibrit −1,320 vs MP2 −1,321 kcal/mol.
- Tüm bölgeyi VQE ile çözmek 696 kübit, tam kompleks ~3100 kübit gerektirir.

**Triyaj kararı:** β-CD için "klasik yeterli". Kuantum düzeltmesi kimyasal doğruluk eşiğinin (1 kcal/mol) yaklaşık 1000 kat altında.

---

## B. Porfin: kuantum bölgesinin literatürle doğrulanması

| Hal | CAS(4,4) = VQD | +NEVPT2 | f | Baskın konfigürasyonlar |
|---|---|---|---|---|
| S0 | 0 | 0 | – | HF 0,93 |
| S1 (Qx) | 3,62 | **1,87** | 0,033 | H→L 0,54; H-1→L+1 0,44 |
| S2 (Qy) | 3,74 | – | 0,068 | H-1→L 0,59; H→L+1 0,39 |
| S3 (By) | 5,35 | – | 2,13 | H→L+1 0,47; H-1→L 0,31 |
| S4 (Bx) | 5,52 | – | 2,76 | H→L 0,45; H-1→L+1 0,37 |
| T1 | 2,71 | 1,95 (güvenilmez) | – | H→L 0,70 |

**Literatür karşılaştırması (Qx bandı; JCTC 2019, 15, 1679, Tablo 3):**

| CASSCF | CASSCF/NEVPT2 | CASPT2(4,4) | MRCI | Deney | Bu çalışma (VQD+NEVPT2) |
|---|---|---|---|---|---|
| 3,42 | 1,96 | 1,70 | 1,90 | 1,98–2,02 | **1,87** |

**Yorum:** 8 kübitlik VQD sonucu ile NEVPT2 birlikte deneyden 0,11–0,15 eV uzakta. Bu fark, literatürdeki CASPT2 sonucundan (0,28 eV fark) daha küçük ve literatür NEVPT2 sonucuna 0,09 eV yakın.

---

## C. 8 porfirin türevi: tarama v3 (`data/v3_tarama_sonuclari.csv`)

Yöntem: B3LYP/def2-SVP orbitalleri + CAS(4,4) = 8 kübit VQD + SC-NEVPT2; geometri GFN2-xTB.

| Molekül | S1 CAS | S1 +NEVPT2 | Deney S1 | Fark | f(S1) | 2. konf. ağırlığı | T1 ΔSCF | Triyaj |
|---|---|---|---|---|---|---|---|---|
| Porfin | 3,62 | 1,87 | 1,98 | −0,11 | 0,033 | 0,435 | 2,00 (deney 1,58) | Kuantum gerekli |
| Zn-porfin | 3,75 | 2,07 | 2,18 | −0,11 | 0,079 | 0,367 | 2,08 ⚠ | Kuantum gerekli |
| Mg-porfin | 3,72 | 2,04 | 2,17 | −0,13 | 0,073 | 0,373 | 2,05 ⚠ | Kuantum gerekli |
| Pd-porfin | 3,80 | 2,14 | – | – | 0,085 | 0,355 | 2,16 ⚠ | Kuantum gerekli |
| Tetrametilporfin* | 3,42 | 1,64 | ~1,88 | −0,24 | 0,001 | 0,470 | 1,58 | Kuantum gerekli |
| Klorin | 3,45 | 1,98 | ~1,95 | +0,03 | 0,579 | 0,098 | 1,53 | Klasik yeterli |
| Porfirazin | 3,74 | 2,24 | 2,065 | +0,18 | 0,633 | 0,139 | 1,65 | Klasik yeterli |
| Bakteriyoklorin* | 3,59 | 2,23 | 1,74 | +0,49 | 0,422 | 0,168 | 1,65 | Klasik yeterli (model yetersiz) |

\* Deney değeri analog veya vekil molekülden alındı.

| Metrik | Değer |
|---|---|
| MAE, doğrudan deneyli 5 molekül | **0,111 eV** |
| MAE, deneyli 7 molekül | 0,184 eV |
| VQD − CASCI farkı (en büyük) | ≤ 10⁻⁸ eV |
| Süre | ~20 dk/molekül (16 thread), ~18 dk porfin |

**Bilimsel yorum:**
1. **Metal etkisi:** Metal eklenince Q bandı maviye kayıyor (porfin 1,87 → Zn 2,07, Mg 2,04). Deneydeki yönle uyumlu.
2. **Simetri bozulması:** Klorinde ve porfirazinde simetri bozuluyor ve güçlü Q bandı ortaya çıkıyor (f = 0,58–0,63). Simetrik porfinlerde ise Q bandı zayıf (f < 0,1). Bu, Gouterman modeliyle ve klinik tercihle uyumlu.
3. **Triyaj ve Gouterman:** Simetrik ve metalli porfinlerde iki geçiş neredeyse eş enerjili ve güçlü karışıyor, bu yüzden "kuantum gerekli" çıkıyor. Simetrisi bozulmuş halkalarda tek konfigürasyon baskın.
4. **Sınırlar:**
   - Bakteriyoklorinde 4 orbital yetersiz.
   - Tetrametilporfin xTB geometrisi çok bükük (RMS z = 0,48 Å).
   - ΔSCF ile hesaplanan T1, simetrik moleküllerde S1'e eşit çıkıyor; güvenilmez.

---

## D. Klasik referans: TD-DFT B3LYP/def2-SVP (`data/tddft_b3lyp.csv`)

| Molekül | TD-DFT S1 | f | TDA T1 | Deney S1 | TD-DFT farkı | Bizim fark |
|---|---|---|---|---|---|---|
| Porfin | 2,348 | 0,001 | 1,832 (deney 1,58) | 1,98 | +0,37 | −0,11 |
| Klorin | 2,369 | 0,174 | 1,540 | 1,95 | +0,42 | +0,03 |
| Zn-porfin, Mg-porfin, porfirazin, bakteriyoklorin | Hesaplandı, değerler henüz arşive aktarılmadı (`~/tarama/*_tddft.json`) | | | | | |
| Tetrametilporfin, Pd-porfin | Hesap sürüyor | | | | | |

**Ortak sette (porfin + klorin) ortalama mutlak hata:** TD-DFT 0,39 eV, bizim yöntem 0,07 eV (yaklaşık 5,6 kat daha yakın).
**T1:** TDA-B3LYP porfinde 1,83 eV veriyor (deney 1,58). Bizim ΔSCF değerimiz (2,00) bundan daha kötü; T1 için TDA'ya geçilmeli.
**Maliyet:** TD-DFT kullanıcının WSL ortamında molekül başına yaklaşık 1,7 saat sürdü (RPA, 4 kök, grid 2, DF). Bu süre, aynı makinede v3 iş akışından (~20 dk) daha uzun. Yavaşlığın büyük olasılıkla nedeni, eşzamanlı çalışan iki süreç ve RPA Davidson yakınsaması. Ölçeklemede bu bir optimizasyon konusu.

---

## E. Klinik ilaçlar: çekirdek hesabı ve klinik ışık (`data/klinik_ilaclar.csv`)

| İlaç | Çekirdek | Klinik ışık | Bizim çekirdek S1 | Fark | TD-DFT çekirdek | Fark | Hizalama RMSD | Triyaj |
|---|---|---|---|---|---|---|---|---|
| Temoporfin (Foscan) | Klorin | 652 nm (EMA) | 1,98 eV = 626 nm | 26 nm | 2,369 eV = 523 nm | 129 nm | 0,40 Å | Klasik yeterli |
| Talaporfin (Laserphyrin) | Klorin | 664 nm (664–667) | 626 nm | 38 nm | 523 nm | 141 nm | 0,24 Å | Klasik yeterli |
| Protoporfirin IX (5-ALA, Ameluz) | Porfin | 635 nm (FDA) | 1,87 eV = 663 nm | 28 nm | 2,348 eV = 528 nm | 107 nm | 0,22 Å | Kuantum gerekli |
| Padeliporfin (Tookad) | Bakteriyoklorin | 753 nm (SmPC) | 2,23 eV = 556 nm | 197 nm | – | – | – | Model yetersiz |

**Yorum:**
- Yan grupları ve çözücüyü içermeyen çekirdek hesabı, klinik ışığa 26–38 nm yakın; TD-DFT 107–141 nm uzak.
- Yakın kızılötesi (NIR, 700–900 nm) yeni nesil ilaçlarda model yetersiz. Bu, platformun en değerli pazar bölgesi ve aynı zamanda bir sonraki teknik hedef (CAS(8,8) = 16 kübit).
- Temoporfin C44H32N4O4 (84 atom), Talaporfin 93 atom, PpIX 76 atom.

---

## F. Gerçek kuantum donanımı: IBM ibm_kingston (`data/ibm_donanim_porfin.csv`)

| Hesap | S1 (eV) | Tam çözümden fark |
|---|---|---|
| Tam çözüm (CASCI, klasik) | 3,617 | – |
| VQE, ideal simülatör | 3,616 | −0,001 |
| VQE, FakeKingston gürültü simülasyonu (ham) | 2,825 | −0,792 |
| VQE, FakeKingston (düzeltilmiş) | 4,095 | +0,478 |
| **VQE, ibm_kingston (ham)** | **1,141** | **−2,477** |
| **VQE, ibm_kingston (düzeltilmiş)** | **4,546** | **+0,928 ± 0,190** |

**Ayrıntılar:**
- İş kimliği: dbbl9ns24mrc73chvt8g. Donanım: ibm_kingston (156 kübit). Tarih: 20 Ekim 2026, 14:47. Plan: Open Plan, Job modu.
- Devre: kompakt ADAPT, S0 için 8 ve S1 için 5 çift uyarma kapısı.
  - Transpile sonrası (FakeKingston) 264 / 170 CZ kapısı; devre derinliği yaklaşık 470 / 305.
- Ölçüm: |c| ≥ 1e-3 olan 121 Pauli terimi, 52 QWC grubu, 4 devre × 4096 atış.
- Gürültü azaltma:
  - TREX (resilience_level 1), kapı twirling, dinamik dekuplaj.
  - Açıları sıfırlanmış gürültü tahmin devresiyle depolarizan ölçekleme.
- Gürültü faktörü f: S0 0,41, S1 0,53.

**Yorum:**
- Gerçek donanım, gürültü modelinin öngördüğünden daha gürültülü çıktı.
- Düzeltme hatayı 2,5 eV'den 0,93 eV'ye indiriyor. Ama bu, kimyasal doğruluk hedefinin (~0,1 eV) yaklaşık 10 katı.
- Bugünkü hata düzeltmesiz donanımda 8 kübitlik bölge bile güvenilir değil.
- Bu sonuç alanın genel durumuyla tutarlı: Q4Bio'nun 5 M$'lık "kanıtlanabilir kuantum avantajı" ödülü kimseye verilmedi.

---

## G. Küresel muadil karşılaştırması

### G1. Akademik ve endüstriyel öncü çalışmalar
| Kriter | **Kuantum Triyaj (bu proje)** | Algorithmiq ΔADAPT-VQE (arXiv:2404.16149, 2024) | Algorithmiq Q4Bio / TLD1433 (2026) | Xanadu (arXiv:2512.15889, 2025) | JCTC 2019 (klasik çok-referanslı) |
|---|---|---|---|---|---|
| Molekül ailesi | Porfirin, klorin, metal-porfin, porfirazin, bakteriyoklorin + 3 klinik ilaç çekirdeği | 6 BODIPY | Ru(II) kompleksi TLD1433 (Faz II, mesane kanseri) | 4 BODIPY sınıfı (orijinal, triazolil aza-, bromlu aza-, Pt-BODIPY) | Porfin ve konjuge sistemler |
| Hesaplanan özellik | S1, f, T1, triyaj | S0→S1 dikey uyarılma | S0 temel hal enerjisi | Toplam soğurma (700–850 nm), ISC hızı (kaynak tahmini) | Q ve B bantları |
| Kuantum kaynağı | 8 kübit (16 kübit betiği hazır) | 8–16 kübit (UCCSD(4,4)–(8,8)) | 52–100 kübit | 180–350 mantıksal kübit, 10⁷–10⁹ Toffoli | – |
| Donanım | Simülatör + ibm_kingston (ilk ölçüm) | Gürültüsüz simülatör | ibm_cleveland (156 kübit), gerçek | Yok (hataya dayanıklı gelecek donanım) | Klasik |
| Deneye göre doğruluk | MAE 0,11 eV (5 mol.), 0,18 eV (7 mol.) | ΔUCCSD(8,8) 0,131 eV; klasik ΔCCSD 0,054 eV | Deneyle karşılaştırma yok | Yok | NEVPT2 porfin Qx 1,96 (deney 1,98–2,02) |
| Klasik referansa göre | TD-DFT 0,39 → 0,07 eV (2 mol.) | TD-DFT 0,20–0,39; EOM-CCSD 0,50 | Sabit bağ boyutunda DMRG'den daha düşük enerji | TD-DFT hatası 0,3–0,6 eV (alıntı) | – |
| Gerçek donanım sonucu | Ham −2,48, düzeltilmiş +0,93 eV | Yok | DMRG'yi aynı bağ boyutunda geçti; sınırsız kaynaklı klasik DMRG hâlâ daha iyi | Yok | – |
| Triyaj / seçici kuantum | **Var (temel özellik)** | Yok | Kavramsal ("kuantum kıt") | Yok | – |
| Finansman / durum | Öğrenci prototipi | Wellcome Leap destekli | 2 M$ Q4Bio ödülü (Nisan 2026) | ~500 M$ birleşme sermayesi beklentisi, U. Alberta ortaklığı (Ağustos 2026) | Akademik |

### G2. Yazılım ve hizmet alternatifleri (yöntem ve maliyet açısından)
| Alternatif | Tür | Bizimle ilişkisi |
|---|---|---|
| Gaussian, Q-Chem, ORCA, Schrödinger (Jaguar) | Klasik kuantum kimyası paketleri | TD-DFT/CASPT2/NEVPT2 sunarlar. Biz açık kaynak PySCF ile aynı sınıf yöntemleri kullanıyoruz. Avantajımız triyaj ve porfirin odaklı otomasyon. Lisans ücretli olanlar var; ORCA akademik kullanımda ücretsiz. |
| Qiskit Nature, PennyLane qchem, OpenFermion | Açık kaynak kuantum kimyası kütüphaneleri | Altyapı olarak PennyLane ve Qiskit'i kullanıyoruz. Bunlar uçtan uca bir ilaç tarama ürünü değil. |
| Algorithmiq Aurora, Quantinuum InQuanto, Xanadu çerçevesi | Ticari ve araştırma kuantum kimyası platformları | Doğrudan rakip sınıfı. Odak molekül aileleri BODIPY ve Ru kompleksleri; porfirin ailesi ve klinik ilaç çekirdekleri bizim ayrışma alanımız. |
| IBM Quantum (donanım) | QPU erişimi | Open Plan ücretsiz (ayda 10 dk). Ücretli plan "Pay-As-You-Go" dakikası 96 $'dan başlıyor. Ölçümümüz tahminen 3–6 QPU dakikası sürdü. |

### G3. Neredeyiz? (değerlendirme)
| Eksen | Konum | Gerekçe |
|---|---|---|
| **Doğruluk** | Güçlü (porfirin ailesinde) | 0,11 eV MAE; Algorithmiq'in BODIPY'de simülatörde elde ettiği 0,13 eV ile aynı sınıf. Molekül setleri farklı, doğrudan kıyas değil. TD-DFT'yi iki molekülde ~5 kat geçiyor. |
| **Kuantum donanımı** | Başlangıç | Tek ölçüm, 8 kübit, gürültü baskın. Algorithmiq 100 kübite kadar donanım + DMRG hibriti ile çok ileride. |
| **Ölçek** | Zayıf | 8 molekül; NIR (bakteriyoklorin) ve tam ilaç (yan grup + çözücü) eksik. T1/ISC güvenilmez. |
| **Yöntem özgünlüğü** | Orta–güçlü | Triyaj + porfirin odağı literatürde boş bir niş. Kullanılan yöntemlerin kendisi (CAS, NEVPT2, VQE) standart. |
| **Maliyet** | Çok düşük | Açık kaynak yazılım, kişisel bilgisayar (16 thread), ücretsiz IBM Open Plan. |
| **Ticari olgunluk** | Erken | Müşteri görüşmesi yok, web arayüzü yok, tek kişilik ekip. |

**Sonuç:**
- Bilimsel doğruluk, porfirin ailesi için uluslararası öncü çalışmalarla aynı mertebede (~0,1 eV). Bu doğruluğun kaynağı klasik çok-referanslı düzeltme; bugün kuantum donanımından gelmiyor.
- Rakiplerden ayrıştığımız iki nokta var: (1) klinikteki PDT ilaçlarının çoğunu oluşturan porfirin ailesine odaklanmak, (2) kuantum kaynağını seçici kullanan triyaj.
- En büyük açıklar: NIR ilaçları, tam molekül hesabı, donanım gürültüsü ve ürünleşme.
