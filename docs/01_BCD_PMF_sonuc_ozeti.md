# β-CD + siklopentanol — sonuç özeti

## 1. Kurulum
- GROMACS 2026.1, GROMOS 54A7 (ATB): β-CD molid 345755 (AOYB, AM1 yükleri, doğal stereokimya), siklopentanol molid 859 (G108). SPC su, PME, 1.4 nm kesme, 300 K.
- 10 ns serbest MD (~/bcd_cpo): misafir ayrışmadı.
- Umbrella sampling (~/bcd_pmf): 37 pencere, ξ = CD ekseni boyunca z, k=1000 kJ/mol/nm², silindirik düz-dip kısıt R=0.4 nm, 7 glikozidik O'ya posre, WHAM (-b 200). ΔG° = −RT ln K, K = πR²∫e^(−W/kT)dξ · 0.6022.
- Kuantum ortamı: WSL ~/qenv (Python 3.14): pyscf 2.14, pennylane 0.45, tblite 0.7, dftd3, dftd4.

## 2. PMF vs deney (kapandı)
- **ΔG°(PMF) = −2.7 ± 0.6 kcal/mol**, deney (ITC, Rekharsky 1997) **−3.05 kcal/mol** (ΔH ≈ −1.1, TΔS ≈ +2.0 → entropi güdümlü). Fark 0.33 kcal/mol; literatür RMSE (0.9–1.8) içinde.
- Sınırlama: merkez pencerelerde yönelim kilitlenmesi (uç platolar 4–4.6 kJ/mol farklı). Pencere uzatma = isteğe bağlı gelecek çalışma.

## 3. Bağlanma modları (temsilci_kare.py, 901 kare)
- Mod A (OH dar/primer kenarda): %53.5; Mod B (OH geniş/sekonder kenarda): %38.7; geçiş %7.8.
- Halka her iki modda da boşluk merkezinde (h ≈ ±0.05 nm, r ≈ 0.04 nm); yalnız OH ucu taraf değiştiriyor.
- Popülasyondan ΔG(A→B) ≈ +0.2 kcal/mol → modlar neredeyse eş enerjili.
- Temsilci kareler: A t=9990 ps, B t=3530 ps (mod_A/B .pdb/.xyz).

## 4. Kuantum iş akışı doğrulaması — misafir tek başına (misafir_qm.py, vqe_misafir.py)
- RHF/def2-SVP (134 baz fonk.), HOMO–LUMO ≈ 16 eV; MP2 korelasyonu −571 kcal/mol; CASCI(4e,4o) yalnız −0.05.
- VQE (8 kübit, Jordan–Wigner, 26 parametreli tekli/çiftli uyarılma ansatzı, BFGS): kübit Hamiltonyeni HF'yi birebir verdi; VQE − CASCI = 0.00000 kcal/mol. 12 kübit default.qubit+backprop ile >10 dk.

## 5. Tam kompleks etkileşim enerjisi (xtb_etkilesim.py, GFN2-xTB, 31 MD karesi)
| | E_int gaz | E_int su (ALPB) | D4 dispersiyon |
|---|---|---|---|
| Mod A (n=16) | −9.70 ± 1.62 | −6.78 ± 1.42 | −8.94 ± 0.86 |
| Mod B (n=13) | −8.67 ± 1.46 | −5.69 ± 1.35 | −9.06 ± 0.80 |
| Tümü (n=31) | −9.34 ± 1.63 | −6.33 ± 1.45 | −9.00 ± 0.84 |
- Gaz fazı etkileşimin ~%96'sı dispersiyon → **bağlanma dispersiyon güdümlü**. Çözücü ~3 kcal/mol zayıflatıyor.
- A, B'den ~1 kcal/mol güçlü (≈1.8σ); fark elektrostatik kaynaklı.

## 6. Bölge (parçalama) modeli (fragman_xtb.py)
- CD, itp bağ grafiğinden 7 glukoz birimine ayrıldı (glikozidik bağlar kesildi, H kapakları).
- Şema 1 (7 tekli glukoz toplamı): −12.39 vs tam −9.20 → hata −3.19 ± 0.22 (H kapaklarından çift sayım).
- **Şema 2 (MFCC: örtüşen ikililer − tekliler): −9.38, hata −0.18 ± 0.13 kcal/mol (%2) → bölge modeli doğrulandı.**
- 7 glukozun katkısı istatistiksel olarak eşit (−1.5…−1.9 ± 0.5–1.0) → sıcak nokta yok, simetrik dispersiyon sarması; tek bölge temsilci.

## 7. Bölgede ab initio (bolge_qm.py; mod A, misafir + glukoz 2, 39 atom, def2-SVP, counterpoise)
| Yöntem | E_int (kcal/mol) |
|---|---|
| HF | +1.39 (itme) |
| B3LYP | +0.73 (itme) |
| MP2 | −1.32 |
| B3LYP-D3(BJ) | −2.83 |
| GFN2-xTB | −2.80 |
- Korelasyonsuz yöntemler itme, korelasyonlu olanlar çekim veriyor: MP2 − HF = −2.71 → **bağlanma bir elektron korelasyonu (dispersiyon) etkisi**.
- B3LYP-D3 ≈ xTB → xTB tabanlı tam kompleks ve parçalama sonuçları bağımsız olarak doğrulandı.
- MP2/def2-SVP küçük baz nedeniyle dispersiyonu az tahmin ediyor.

## 8. Bölgede VQE + hibrit birleştirme (bolge_vqe.py)
- Hibrit formül: E_int = HF + MP2(tüm 348 MO) + [VQE − MP2]_aktif uzay.
- Aktif uzay: her parçadan 1 dolu + 1 boş orbital (dimer 4e,4o = 8 kübit; monomerler 4 kübit).
- VQE korelasyon katkısı +0.0032 kcal/mol = CASCI (fark 2·10⁻⁸); aynı uzayda MP2 +0.0023; **kuantum düzeltmesi +0.0009 kcal/mol**; tam korelasyon etkileşiminin ~%0.1'i (gürültü düzeyi).
- Hibrit −1.320 vs MP2 −1.321 kcal/mol.
- Tüm bölgeyi VQE ile çözmek: 696 kübit (def2-SVP); tam kompleks ~3100 kübit.

## 9. Genel sonuç
- MD/PMF deneyle uyumlu (−2.7 vs −3.05 kcal/mol).
- Bağlanma, halka boyunca simetrik dağılmış dispersiyonla sağlanıyor; ab initio olarak kanıtlandı.
- Kuantum iş akışı (PySCF → aktif uzay → Jordan–Wigner → VQE) ve bölge+hibrit birleştirme çerçevesi kuruldu ve doğrulandı.
- Bu sistemde kuantum düzeltmesi kimyasal doğruluk eşiğinin (1 kcal/mol) ~1000 kat altında → klasik MP2/DFT-D yeterli; VQE'nin katkı yapabilmesi için yüzlerce kübit veya güçlü korelasyonlu bir alt sistem gerekir.
