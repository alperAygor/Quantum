# Fotosensitizör kuantum triyaj taraması — v3 sonuçları ve literatür karşılaştırması

## Yöntem (v3, kuantum_triyaj_paketi_v3)
- Geometri: RDKit → GFN2-xTB (ASE BFGS). Metal porfinler: porfin iskeletine metal yerleştirilip xTB ile optimize.
- Orbitaller: B3LYP/def2-SVP (DF, grid 1); Pd için def2-SVP ECP.
- Kuantum bölgesi: Gouterman 4 orbitali → CAS(4e,4o) = 8 kübit; VQE/VQD (S0, S1, S2 singlet; T1 triplet), geçiş dipolleri kübit operatöründen.
- Dinamik korelasyon: SC-NEVPT2 (S0, S1).
- T1: dikey ΔSCF (UKS-B3LYP triplet − RKS singlet).
- Triyaj kuralı: S1'de ikinci konfigürasyon ağırlığı > 0.2 → "KUANTUM GEREKLİ".
- Süre: ~20 dk/molekül (16 thread, tek tek; iki paralel iş WSL 7 GB belleği aşıp "Killed" oldu). 8 molekül tamamlandı (2026-10-18).
- v2 (SA-CASSCF S0+S1+T1 ortak orbitaller) terk edildi: porfin S1 2.79 eV, klorin S1 3.33 eV ve f yanlış yönde → küçük aktif uzayda hal karışması.

## Sonuçlar ve deney (eV)
| Molekül | S1 CAS | S1 +NEVPT2 | Deney (Q, S1) | Fark | f(S1) | T1 ΔSCF | Triyaj |
|---|---|---|---|---|---|---|---|
| Porfin | 3.62 | 1.87 | 1.98 (Qx, buhar; Edwards 1971) | −0.11 | 0.033 | 2.00 (deney 1.58, Gouterman & Khalil 1974) ❌ | Kuantum gerekli |
| Klorin | 3.45 | 1.98 | ~1.95 (n-alkan kristal, 4.2 K; Lenth/Baumert) | +0.03 | 0.579 | 1.53 | Klasik yeterli |
| Zn-porfin | 3.75 | 2.07 | 2.18–2.23 (çözelti; Baerends 2002), ZnOEP gaz 2.16 | −0.1…−0.16 | 0.079 | 2.08 ⚠ | Kuantum gerekli |
| Mg-porfin | 3.72 | 2.04 | 2.17 (DMF; Devillers 2013), ~2.2 (matris) | −0.13…−0.16 | 0.073 | 2.05 ⚠ | Kuantum gerekli |
| Pd-porfin | 3.80 | 2.14 | deney değeri bulunamadı (sübstitüe edilmemiş PdP) | – | 0.085 | 2.16 ⚠ | Kuantum gerekli |
| Porfirazin | 3.74 | 2.24 | 2.065 (Qx, Ne matris 8 K; Henchy 2014) | +0.18 | 0.633 | 1.65 | Klasik yeterli |
| Bakteriyoklorin | 3.59 | 2.23 | 1.74 (sübstitüe analog B6, toluen; Reddy/Lindsey 2017) | +0.49 ❌ | 0.422 | 1.65 | Klasik yeterli |
| Tetrametilporfin | 3.42 | 1.64 | ~1.88 (meso-tetraalkil vekil, CH2Cl2; Paliteiro 2005) | −0.24 ⚠ | 0.001 | 1.58 | Kuantum gerekli |

- Ortalama mutlak hata (deneyi olan 7 molekül): ~0.19 eV; bakteriyoklorin hariç ~0.14 eV.
- Trendler: metal → mavi kayma (doğru yön); klorin ve porfirazinde güçlü Q bandı (doğru); triyaj simetrik/metalli porfirinleri "kuantum gerekli", simetrisi bozulmuş klorin/porfirazin/bakteriyoklorini "klasik yeterli" sınıfladı (Gouterman modeliyle uyumlu).
- Sınırlar: bakteriyoklorin (4 orbital yetersiz → daha büyük aktif uzay); tetrametilporfin (xTB geometrisi çok bükük, RMS z 0.48 Å, f≈0 şüpheli); ΔSCF T1 simetrik moleküllerde güvenilmez (porfin 2.00 vs deney 1.58; Zn/Mg/Pd'de T1≈S1).
- Grafik: tarama_sonuclari.png (hesap–deney dumbbell + f(S1) ve triyaj).

## Kaynaklar
- Edwards, Dolphin, Gouterman, Adler, J. Mol. Spectrosc. 38, 16 (1971) — porfin buhar spektrumu (Šeda, Burda, Leszczynski, J. Comput. Chem. 26, 294 (2005), Tablo 5–6 üzerinden).
- JCTC 15, 1679 (2019), doi:10.1021/acs.jctc.8b01196, Tablo 3.
- Gouterman & Khalil, J. Mol. Spectrosc. 53, 88 (1974) — porfin fosforesansı, T1 ≈ 785 nm (1.58 eV).
- Lenth, Baumert et al. — klorin, n-alkan, 4.2 K, ~635 nm (künye doğrulanmalı).
- Reddy, …, Lindsey, Molecules 22, 634 (2017), Tablo 2 (bakteriyoklorin B6).
- Devillers et al., C. R. Chimie 16, 540 (2013) — ZnP, MgP DMF.
- Baerends et al., Coord. Chem. Rev. 230, 5 (2002), Tablo 3–4.
- Henchy, McCaffrey et al., J. Chem. Phys. 141, 124303 (2014), Tablo I — porfirazin.
- Paliteiro & Sobral, Electrochim. Acta 50, 2445 (2005) — meso-tetraalkilporfirinler.

## Sonraki adımlar
1. Bakteriyoklorin ve porfin için büyük aktif uzay (8e,8o) klasik kontrolü.
2. Tetrametilporfin geometrisini düzelt (düzlemsel başlangıç).
3. T1 için daha iyi yöntem (TDA-B3LYP triplet veya CAS(8,8)+NEVPT2).
4. Rapor; görseller (molden/cube → Avogadro/ChimeraX); IBM Quantum donanım ölçümü; startup sunumu.
