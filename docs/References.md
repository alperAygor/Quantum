# Kaynakça ve Literatür

Her kaynağın yanında projede ne için kullanıldığı yazılıdır.

## 1. Bilimsel makaleler

### 1.1 PDT + kuantum hesaplama (rakip ve öncü çalışmalar)
1. **Zhou, Y.; Casares, P. A. M.; Dhawan, D.; Loaiza, I.; Jahangiri, S.; Lang, R. A.; Arrazola, J. M.; Fomichev, S. (Xanadu).** *Quantum Algorithms for Photoreactivity in Cancer-Targeted Photosensitizers.* arXiv:2512.15889 [quant-ph], 17 Aralık 2025.
   - **Kullanım:** Rakip analizi.
     - BODIPY türevleri; terapötik pencere 700–850 nm (s.3).
     - Kaynak tahminleri: 180–350 mantıksal kübit ve 10⁷–10⁹ Toffoli (Tablo I).
     - TD-DFT'nin 17 BODIPY'deki 0,3–0,6 eV hatası (s.7, literatürden alıntı).
     - Hesaplanmış enerji ve deney karşılaştırması içermiyor.
   - Kullanıcının arkadaşı tarafından incelendi (dosya: "Enes Makale.pdf").
2. **Nykänen, A.; Thiessen, L.; Borrelli, E.-M.; Krishna, V.; Knecht, S.; Pavošević, F. (Algorithmiq, Cleveland Clinic).** *ΔADAPT-VQE: Toward Accurate Calculation of Excitation Energies on Quantum Computers for BODIPY Molecules With Application in Photodynamic Therapy.* arXiv:2404.16149 [physics.chem-ph], 24 Nisan 2024.
   - **Kullanım:** Doğruluk kıyası. 6 BODIPY için deneye göre ortalama mutlak hata (Şekil 2, s.11–12):
     - TD-DFT: BLYP 0,199, B3LYP 0,368, CAM-B3LYP 0,392 eV
     - EOM-CCSD 0,498 eV, ΔHF 0,496 eV
     - Klasik ΔCCSD 0,054 eV (en iyi), ΔUCCSD(8,8) 0,131 eV
   - Gürültüsüz simülatörde yapıldı (s.10). Aktif uzayların klasik olarak çözülebildiği kabul ediliyor (s.13). Hedef < 0,1 eV (s.14). Terapötik pencere 750–900 nm alınmış (s.4).
   - Dosya: "3.Makale.pdf".
3. **Q4Bio Team Algorithmiq.** *Photodynamic Therapy as a Driver for Useful Quantum Advantage in Chemistry: Perspective from the Q4Bio program.* 2026.
   - **Kullanım:** Gerçek donanım kıyası.
     - TLD1433 (Ruvidar, Ru(II), Faz II mesane kanseri, Şekil 3).
     - ibm_cleveland (156 kübit), 52 kübitlik aktif uzayda QB-DMRG; DMRG'den sabit bağ boyutunda daha düşük enerji (Şekil 8). 100 kübite kadar ölçekleme (Şekil 9, s.16).
     - Sınırsız kaynaklı klasik DMRG hâlâ daha güçlü (s.16). Soğurma ve ISC gelecek çalışma (s.17).
     - "Yararlı kuantum avantajı" için 3 şart: bugünkü donanımda çalışabilirlik, biyolojik olarak anlamlı sistem, en iyi klasik yönteme karşı doğrulama.
     - "20 yılı aşkın sürede çok az temelde yeni fotosensitizer" ifadesi (s.3).
   - Dosya: "Q4Bio_Perspective_Paper.pdf".

### 1.2 Porfirin fotofiziği ve deney değerleri
4. **Gouterman, M.** *Spectra of porphyrins.* J. Mol. Spectrosc. **6**, 138–163 (1961).
   - **Kullanım:** Dört orbital modeli. Kuantum bölgesinin (HOMO-1, HOMO, LUMO, LUMO+1) seçim gerekçesi.
5. **JCTC 2019, 15(3), 1679–1689**, doi:10.1021/acs.jctc.8b01196. *Constructing Molecular π-Orbital Active Spaces for Multireference Calculations of Conjugated Systems* (ücretsiz tam metin: PMC6526033).
   - **Kullanım:** Porfin Q/B bantları, Tablo 3. CASSCF 3,42; NEVPT2 1,96; CASPT2 1,70; MRCI 1,90; deney 1,98–2,02 eV. Bizim 1,87 eV ile karşılaştırma.
6. **Edwards, L.; Dolphin, D. H.; Gouterman, M.; Adler, A. D.** J. Mol. Spectrosc. **38**, 16 (1971).
   - **Kullanım:** Porfin buhar spektrumu, Qx 1,98 eV. Šeda, Burda, Leszczynski, J. Comput. Chem. **26**, 294 (2005), Tablo 5–6 üzerinden.
7. **Gouterman, M.; Khalil, G.-E.** J. Mol. Spectrosc. **53**, 88 (1974).
   - **Kullanım:** Porfin fosforesansı, T1 ≈ 785 nm (1,58 eV).
8. **Lenth, Baumert vd.**
   - **Kullanım:** Klorin, n-alkan, 4,2 K, ~635 nm (1,95 eV).
   - Not: Künye arşivde tam doğrulanmadı.
9. **Reddy, K. R.; …; Lindsey, J. S.** Molecules **22**, 634 (2017), Tablo 2.
   - **Kullanım:** Bakteriyoklorin analoğu B6, toluende 1,74 eV (vekil deney değeri).
10. **Devillers, C. H. vd.** C. R. Chimie **16**, 540 (2013).
    - **Kullanım:** ZnP ve MgP Q bandı, DMF'de 2,17 eV (MgP).
11. **Baerends, E. J. vd.** Coord. Chem. Rev. **230**, 5 (2002), Tablo 3–4.
    - **Kullanım:** Zn-porfin Q bandı 2,18–2,23 eV, ZnOEP gaz 2,16 eV.
12. **Henchy, C.; McCaffrey, J. G. vd.** J. Chem. Phys. **141**, 124303 (2014), Tablo I.
    - **Kullanım:** Porfirazin Qx 2,065 eV (Ne matris, 8 K).
13. **Paliteiro, C.; Sobral, A.** Electrochim. Acta **50**, 2445 (2005).
    - **Kullanım:** Meso-tetraalkilporfirinler, tetrametilporfin için vekil deney değeri ~1,88 eV.

### 1.3 Siklodekstrin
14. **Rekharsky, M. V.; Inoue, Y. vd. (1997).** Siklodekstrin–alkol kompleksleşmesi, ITC.
    - **Kullanım:** β-CD + siklopentanol için ΔG° = −3,05 kcal/mol, ΔH ≈ −1,1, TΔS ≈ +2,0. PMF doğrulaması.
    - Not: Tam künye arşivde ayrıca doğrulanmalı.

## 2. Klinik ilaç bilgileri (dalga boyları)
15. **EMA — Foscan (temoporfin) ürün bilgisi.**
    - **Kullanım:** Aktivasyon ışığı 652 nm.
16. **Wikipedia — Talaporfin.** https://en.wikipedia.org/wiki/Talaporfin
    - **Kullanım:** "664–667 nm" kırmızı ışık (Laserphyrin, Japonya).
17. **FDA — Ameluz (aminolevulinik asit) etiketi, 2017.** https://www.accessdata.fda.gov/drugsatfda_docs/label/2017/208081s002lbl.pdf
    - **Kullanım:** BF-RhodoLED "yaklaşık 635 nm" kırmızı ışık; PpIX aktivasyonu.
18. **emc — Tookad (padeliporfin) SmPC.** https://medicines.org.uk/emc/product/9366/smpc
    - **Kullanım:** "753 nm" lazer; yakın kızılötesi, bakteriyoklorin çekirdeği.
19. **PubChem ve Wikipedia ilaç kayıtları:** Temoporfin CID 60751, Talaporfin CID 5486799, Protoporfirin IX CID 4971.
    - **Kullanım:** SMILES ve yapılar. PpIX SMILES'ı elle yazıldı.

## 3. Haber ve sektör kaynakları
20. **Quantum Computing Report — Xanadu outlines quantum framework for photodynamic cancer therapy research** (19 Aralık 2025). https://quantumcomputingreport.com/xanadu-outlines-quantum-framework-for-photodynamic-cancer-therapy-research/
    - **Kullanım:** Xanadu'nun PDT çerçevesi ve Crane Harbor birleşmesi (~500 M$).
21. **Indian Pharma Post — Xanadu and University of Alberta team up** (18 Ağustos 2026). https://www.indianpharmapost.com/digitisation/xanadu-and-university-of-alberta-team-up-to-use-quantum-computing-in-cancer-drug-discovery-21433
    - **Kullanım:** Xanadu–U. Alberta (Alex Brown) ortaklığı.
22. **Quantum Computing Report — Algorithmiq wins $2 million Q4Bio prize** (16 Nisan 2026). https://quantumcomputingreport.com/algorithmiq-wins-2-million-q4bio-prize-for-quantum-driven-cancer-therapy-simulation/
    - **Kullanım:** 2 M$ ödül; 5 M$ "kanıtlanabilir kuantum avantajı" ödülü verilmedi.
23. **StartupResearcher — Algorithmiq Wins $2M Wellcome Leap Prize.** https://www.startupresearcher.com/news/algorithmiq-wins-usd2-million-wellcome-leap-prize-for-quantum-cancer-research
    - **Kullanım:** Aynı haberin ikinci kaynağı.
24. **DHA — YTÜ Startup House dünyanın en iyi kuluçka merkezleri arasında ilk 10'da** (25 Ocak 2023). https://www.dha.com.tr/teknoloji/ytu-startup-house-dunyanin-en-iyi-kulucka-merkezleri-arasinda-ilk-10da-2195645
    - **Kullanım:** Başvuru stratejisi.
25. **DHA — YTÜ Yıldız Teknopark ve Microsoft'tan girişimcilere destek** (21 Aralık 2022). https://www.dha.com.tr/ekonomi/ytu-yildiz-teknopark-ve-microsofttan-girisimcilere-destek-2179676
    - **Kullanım:** Startup House destekleri (ofis, mentorluk, Azure kredileri).

## 4. IBM Quantum
26. **IBM Quantum — Plans overview.** https://quantum.cloud.ibm.com/docs/guides/plans-overview
    - **Kullanım:** Open Plan: ücretsiz, 28 günlük pencerede 10 dk; 16 Mart 2026'dan itibaren +180 dk promosyon.
27. **IBM Quantum — Pricing.** https://www.ibm.com/quantum/pricing
    - **Kullanım:** Pay-As-You-Go dakikası 96 $'dan başlıyor.
28. **Quantum Computing Report — IBM offers special promotion to Open Plan users** (Mart 2026). https://quantumcomputingreport.com/ibm-offers-special-promotion-to-open-plan-users/
    - **Kullanım:** ibm_kingston (Heron r2, 156 kübit) Open Plan'de erişilebilir.
29. **IBM — Cloud simulators and Lab are now retired** (15 Mayıs 2024). https://quantum.cloud.ibm.com/announcements/en/product-updates/2024-05-15-2024-sunset-final-lab-simulators
    - **Kullanım:** IBM sitesinde artık kod çalıştırılamıyor; yerel ortam veya Colab gerekiyor.
30. **IBM — Qiskit notebook environments** (blog). https://www.ibm.com/quantum/blog/qiskit-notebook-environments
    - **Kullanım:** Alternatif not defteri ortamları.

## 5. Yazılımlar ve veri tabanları
| Kaynak | Kullanım |
|---|---|
| PySCF — https://pyscf.org | DFT, CASCI, NEVPT2, TD-DFT, cubegen |
| PennyLane — https://pennylane.ai | VQE/VQD, qchem |
| Qiskit / qiskit-ibm-runtime — https://www.ibm.com/quantum/qiskit | Donanım erişimi |
| tblite (GFN2-xTB) — https://github.com/tblite/tblite | Geometri, etkileşim |
| RDKit — https://www.rdkit.org | SMILES, 3B, alt yapı eşleşmesi |
| ASE — https://wiki.fysik.dtu.dk/ase | Optimizasyon |
| GROMACS — https://www.gromacs.org | MD, umbrella sampling, WHAM |
| ATB (Automated Topology Builder) — https://atb.uq.edu.au | β-CD (molid 345755) ve siklopentanol (molid 859) topolojileri |
| UCSF ChimeraX — https://www.cgl.ucsf.edu/chimerax | 3B görseller |
| Avogadro 2 — https://www.openchemistry.org | İlk orbital görselleri |
| PubChem — https://pubchem.ncbi.nlm.nih.gov | İlaç yapıları |
| RCSB PDB — 4V6X (ribozom, görsel örnek), 1N5U (insan serum albümini + hem; önerildi, kullanılmadı) | Görsel kalite tartışması |
