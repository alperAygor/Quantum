# YTÜ Startup House başvurusu: v4 durumu (20 Ekim 2026)

## Sunum
- "Kuantum Triyaj – YTÜ Startup House Başvuru Sunumu" (Slides artifact), 16 slayt.
- Akış:
  1. Problem
  2. Neden kuantum
  3. Dört adımlı çözüm
  4. Bölge seçimi
  5. Üç klinik ilaç
  6. 8 molekülün sonuçları
  7. TD-DFT karşılaştırması (F1, yeni)
  8. Triyaj
  9. Klinik ışıkla karşılaştırma (F3, TD-DFT eklendi)
  10. Rekabet
  11. İddialar
  12. İş modeli
  13. Yol haritası
  14. Beklenti ve ekip
  15. Kaynaklar
- Doldurulacak yer tutucular: [bölüm], [Ekip arkadaşı · rolü].

## TD-DFT referansı
- Yöntem: B3LYP/def2-SVP; singletler RPA ile, tripletler TDA ile. Kullanıcının bilgisayarında hesaplandı.

| Molekül | TD-DFT S1 (eV) | f | T1 (eV) | Deney S1 (eV) | TD-DFT farkı | Bizim fark |
|---|---|---|---|---|---|---|
| Porfin | 2,348 | 0,001 | 1,832 (deney 1,58) | 1,98 | +0,37 | −0,11 |
| Klorin | 2,369 | 0,174 | 1,540 | 1,95 | +0,42 | +0,03 |

- Aynı iki molekülde ortalama mutlak hata: TD-DFT 0,39 eV, bizim yöntem 0,07 eV.
- Çekirdekler dalga boyuna çevrildiğinde TD-DFT klorini 523 nm, porfini 528 nm buluyor. Klinik ışık 635–664 nm; fark 107–141 nm. Bizim farkımız 26–38 nm.
- T1 için TDA-B3LYP, bizim ΔSCF yaklaşımımızdan daha iyi: porfinde 1,83 eV'ye karşı 2,00 eV (deney 1,58). Yol haritasına alındı.
- Süre: her molekül yaklaşık 1,7 saat sürüyor. Bu yüzden yalnızca B3LYP çalıştırılıyor (FONK=b3lyp). Kalan moleküller sırayla devam ediyor: Zn, Mg, porfirazin, bakteriyoklorin, tetrametilporfin, Pd.
- Dürüst çerçeve: TD-DFT'ye karşı üstünlük çok-konfigürasyonlu yaklaşımdan (CAS + NEVPT2) geliyor. Bu ölçekte o hesap klasik bilgisayarda da yapılabiliyor.

## Kod
- v4:
  - tddft_referans.py
  - cas_buyuk.py
  - final_figurler.py: MAE artık yalnızca tüm yöntemlerin hesaplandığı ortak molekül setinde hesaplanıyor; F3 grafiğinin x ekseni 500–800 nm.
  - calistir_v4.sh
- v5: ibm_vqe.py ve hesap_kaydet.py. Ayrıntılar IBM_donanim_plani.md'de.
- IBM hesabı açıldı: Open Plan, open-instance, 10 dk; görülen donanımlar ibm_kingston, ibm_fez, ibm_marrakesh. Kayıt yapıldı.

## Doğrulanan sayılar
- Bizim yöntem: 0,111 eV (5 doğrudan deney), 0,184 eV (7 molekül).
- Klinik ışık:
  - Temoporfin: 652 nm (EMA)
  - Talaporfin: 664–667 nm
  - PpIX/Ameluz: 635 nm (FDA)
  - Tookad: 753 nm (SmPC)
- Temoporfin C44H32N4O4, 84 atom.
