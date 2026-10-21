# PDT + kuantum: literatür karşılaştırması (20 Ekim 2026)

Üç makale okundu ve doğrulandı (sayfa numaraları PDF'e göre).

## 1. Xanadu: "Quantum Algorithms for Photoreactivity in Cancer-Targeted Photosensitizers" (arXiv:2512.15889, 17 Ara 2025; dosya: Enes_Makale.pdf)
- Molekül: BODIPY türevleri (orijinal, triazolil aza-, bromlu aza-, Pt-BODIPY).
- Hesaplanan değer YOK; sadece hataya dayanıklı (fault-tolerant) kuantum bilgisayarlar için kaynak tahmini var: 180–350 mantıksal kübit, 10^7–10^9 Toffoli (Tablo I, s.3).
- Hedefledikleri gözlemlenebilirler: terapötik penceredeki (700–850 nm, s.3) toplam soğurma ve ISC hızı.
- 0.3–0.6 eV TD-DFT hatası literatürden alıntı, kendi sonuçları değil (s.7).

## 2. Algorithmiq: "ΔADAPT-VQE … BODIPY … PDT" (arXiv:2404.16149, Nis 2024; dosya: 3.Makale.pdf)
- 6 BODIPY, S0→S1 dikey uyarılma enerjisi, deneyle karşılaştırma (Şekil 2, s.11).
- Ortalama mutlak hatalar:
  - TD-DFT: BLYP 0.199, B3LYP 0.368, CAM-B3LYP 0.392 eV (s.11).
  - EOM-CCSD: 0.498 eV.
  - ΔHF: 0.496 eV.
  - **ΔCCSD (klasik): 0.054 eV, en iyi sonuç.**
  - ΔUCCSD(8,8) (kuantum): 0.131 eV (s.12).
- Hesaplar gürültüsüz simülatörde yapıldı (s.10). Kullanılan aktif uzayların klasik yöntemlerle de çözülebildiği makalede kabul ediliyor (s.13).
- Hedef < 0.1 eV (s.14). Terapötik pencere 750–900 nm alınmış (s.4).

## 3. Algorithmiq Q4Bio perspektif makalesi (2026; dosya: Q4Bio_Perspective_Paper.pdf)
- Molekül: TLD1433 (Ruvidar, Ru(II) kompleksi), mesane kanseri için Faz II klinik denemede (Şekil 3).
- Gerçek donanımda çalıştırılmış: ibm_cleveland (156 kübit), 52 kübitlik aktif uzay. QB-DMRG, aynı bağ boyutunda klasik DMRG'den daha düşük S0 enerjisi buldu (Şekil 8). 100 kübite kadar karşılaştırma Şekil 9'da (s.16).
- Kendi kabulleri: kaynak sınırı olmayan klasik DMRG hâlâ daha güçlü (s.16). Soğurma ve ISC gelecek çalışma (s.17).
- Makalenin "yararlı kuantum avantajı" için koyduğu 3 şart: bugünkü donanımda çalışabilirlik, biyolojik olarak anlamlı sistem, en iyi klasik yöntemlere karşı doğrulama.
- Haberler (QCR, Nis 2026): Algorithmiq 2 milyon $ Q4Bio ödülünü aldı. Kanıtlanabilir kuantum avantajı için konan 5 milyon $'lık büyük ödül kimseye verilmedi.

## Bizim için sonuçlar
- "Klasik yöntem PpIX'te ~1 eV hata yapıyor" iddiası adil değil. Bu karşılaştırma tek konfigürasyonlu tahminle yapıldı; adil karşılaştırma TD-DFT ile olmalı.
- Bizim doğruluğumuz (ortalama hata ~0.14–0.19 eV) NEVPT2'den geliyor, yani klasik bir yöntemden. VQE, CASCI(4,4) sonucunu birebir tekrarlıyor ama tek başına doğruluk avantajı sağlamıyor.
- Üç makalenin hiçbiri porfirin, klorin ya da bakteriyoklorin çalışmamış. Klinikte onaylı PDT ilaçlarının çoğu ise bu aileden. Konumlanma fırsatı burada.
- Yeni nesil hedef NIR (yakın kızılötesi, 700–900 nm) bölgesi, yani bakteriyoklorinler (ör. padeliporfin ~753 nm). Bizim yöntemin başarısız olduğu molekül de tam olarak bakteriyoklorin; (8e,8o) aktif uzay gerekli.
- Yapılacaklar:
  1. 8 molekül için TD-DFT referansı (B3LYP ve CAM-B3LYP).
  2. Bakteriyoklorin CAS(8,8).
  3. IBM donanımında çalıştırma.
  4. Sunumda kuantum avantajı iddia etmemek; Q4Bio'nun 3 şartını kontrol listesi olarak kullanmak.
