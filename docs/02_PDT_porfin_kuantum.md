# PDT yönü — porfin üzerinde kuantum triyaj testi

## Amaç
β-CD/siklopentanol sisteminde kuantum düzeltmesi ≈ 0 çıktı. Triyaj fikrini kanıtlamak için zıt örnek: kuantumun gerekli olduğu bir sistem. Seçilen sistem serbest baz porfin (C20H14N4): bütün porfirin tipi PDT fotosensitizörlerinin iskeleti.

## Yöntem (porfin_geo.py, porfin_qm.py — WSL ~/pdt, ~/qenv)
- Geometri: RDKit → GFN2-xTB optimizasyonu (ASE BFGS, fmax 0.02 eV/Å).
- Orbitaller: B3LYP/def2-SVP (DF, grid level 1), 406 baz fonksiyonu, 162 elektron.
- Kuantum bölgesi: Gouterman 4 orbitali (HOMO-1 = MO 80, HOMO = 81, LUMO = 82, LUMO+1 = 83) → CAS(4e,4o) = 8 kübit, 361 Pauli terimi.
- VQE (taban hal) + VQD (uyarılmış haller): HF + 2 katman tekli/çiftli Givens uyarılma ansatzı; ⟨(S²−s)²⟩ spin cezası, örtüşme cezası β=3.
- Geçiş dipolleri: dipol operatörü de kübitlere çevrildi, ⟨S0|μ|Si⟩ kuantum durumlarından hesaplandı.
- NEVPT2 (strongly contracted, S0, S1, T1; durum-özgü, B3LYP orbitalleri üzerinde): `export PYSCF_TMPDIR=~/pyscf_tmp; nohup python3 porfin_qm.py nevpt2 > nevpt2.log &`.
  - WSL'de /tmp RAM'de duran 3.8 GB'lık tmpfs → ilk deneme "No space left on device" ile durdu (cvcv ara dosyası ~5 GB). PYSCF_TMPDIR = ~/pyscf_tmp ile çözüldü.
  - Süre (kullanıcının makinesi, 16 thread): B3LYP 70 s, CASCI 108 s, VQD ~390 s, NEVPT2 (3 hal) ~600 s; toplam 1103 s ≈ 18 dk/molekül.
- Görsel: porfin_orbitaller.molden (Avogadro 2.0: Analyze → Create Surfaces → Molecular Orbital, MO 80–83, izodeğer 0.02–0.03), porfin_gouterman_orbitaller.png.

## Sonuçlar (eV; VQD = CASCI, fark ≤ 10⁻⁹ eV)
| Hal | CAS(4,4) = VQD | +NEVPT2 | f | Baskın konfigürasyonlar |
|---|---|---|---|---|
| S0 | 0 | 0 | – | HF 0.93 |
| S1 (Qx) | 3.62 | **1.87** | 0.033 (x) | H→L 0.54, H-1→L+1 0.44 |
| S2 (Qy) | 3.74 | – | 0.068 (y) | H-1→L 0.59, H→L+1 0.39 |
| S3 (By) | 5.35 | – | 2.13 (y) | H→L+1 0.47, H-1→L 0.31 |
| S4 (Bx) | 5.52 | – | 2.76 (x) | H→L 0.45, H-1→L+1 0.37 |
| T1 | 2.71 | **1.95 (güvenilmez, aşağıya bak)** | – | H→L 0.70 |
| T2 | 2.91 | – | – | H-1→L+1 0.63 |

## Literatür karşılaştırması
JCTC 2019, 15(3), 1679–1689, doi:10.1021/acs.jctc.8b01196 ("Constructing Molecular π-Orbital Active Spaces for Multireference Calculations of Conjugated Systems"; ücretsiz tam metin PMC6526033), Tablo 3:
| Bant | CASSCF | CASSCF/NEVPT2 | CASPT2 (4e,4o) | MRCI | Deney | **Bu çalışma CAS / +NEVPT2** |
|---|---|---|---|---|---|---|
| Qx | 3.42 | 1.96 | 1.70 | 1.90 | 1.98–2.02 | **3.62 / 1.87** |
| Qy | 3.64 | 2.37 | 2.26 | 2.36 | 2.33–2.42 | 3.74 / – |
| B | 5.42 | 3.06 | 2.91 | 3.04 | 3.13–3.33 | 5.35 / – |
| B | 5.43 | 3.09 | 3.04 | 3.06 | 3.13–3.33 | 5.52 / – |
- **S1 (Qx) doğrulandı:** 8 kübit (VQD) + NEVPT2 = 1.87 eV; deneyden 0.11–0.15 eV, literatür NEVPT2'den 0.09 eV, CASPT2'den (1.70) daha iyi.
- Dört orbital modeli: M. Gouterman, J. Mol. Spectrosc. 6, 138–163 (1961). Q zayıf / B güçlü örüntüsü kuantum devresinden yeniden üretildi.
- HOMO-1/HOMO (−5.53/−5.40 eV) ve LUMO/LUMO+1 (−2.46/−2.40 eV) neredeyse dejenere → güçlü konfigürasyon karışımı → triyaj "kuantum gerekli".
- HOMO'nun yoğunluğu meso karbonlarda → meso sübstitüsyonu (ör. temoporfin/mTHPC) Q bandını ayarlar.
- Deneysel T1 değeri: henüz doğrulanmış kaynak yok (aday: "Porphyrin free base phosphorescence", J. Mol. Spectrosc.). Singlet oksijen ¹Δg 0.98 eV.

## T1 sorunu (dürüst değerlendirme)
- NEVPT2 sonrası T1 (1.95) > S1 (1.87): porfirinlerde fiziksel olarak yanlış sıralama (T1 her zaman S1'in altında olmalı).
- NEVPT2 düzeltmesi S1 için −1.75 eV, T1 için yalnız −0.75 eV → singlet ve triplet dengesiz işleniyor.
- Olası nedenler: orbitaller tek bir singlet B3LYP hesabından (durum-ortalamalı CASSCF değil); CAS(4,4) çok küçük, triplet için dinamik korelasyon dengesi bozuk.
- Çözüm adayları: (a) S0/S1/T1 üzerinden durum-ortalamalı (SA-)CASSCF orbitalleri + NEVPT2; (b) daha büyük π aktif uzayı (ör. JCTC 2019 PiOS yöntemi) → daha fazla kübit.
- Sonuç: "T1 > 0.98 eV → singlet oksijen üretebilir" kararı büyük olasılıkla doğru ama bu hesapla henüz kanıtlanmış sayılmaz.

## Sonraki adımlar (kullanıcının onayladığı plan)
1. T1 düzeltmesi: SA-CASSCF(4,4) orbitalleri ile NEVPT2'yi tekrar et (S1–T1 sıralamasını kontrol et).
2. **Toplu triyaj script'i:** molekül listesi (SMILES) → her molekül için triyaj kararı, S1 (+NEVPT2), T1, f → CSV.
   - Hedef: 10–20 fotosensitizör (porfin, klorin, bakteriyoklorin, Zn-porfirin, temoporfin iskeleti vb.).
3. **Ölçekleme:** önce kullanıcının bilgisayarında 2–3 molekül (≈18 dk/molekül → 20 molekül ≈ 6 saat sıralı); sonra TRUBA veya kiralık makinede paralel.
4. WSL belleğini artır (.wslconfig memory=…, `wsl --shutdown`).
5. Tüm proje için akış şeması + kaynaklar + kavram sözlüğü içeren rapor dokümanı.
6. β-CD + fotosensitizör: MD/PMF (taşıyıcı) + temsilci karelerde uyarılmış hal.
