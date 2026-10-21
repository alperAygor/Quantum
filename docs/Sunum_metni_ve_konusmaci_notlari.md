# Sunum Metni ve Konuşmacı Notları

Kaynak: Claude Slides artifact'i "Kuantum Triyaj – YTÜ Startup House Başvuru Sunumu" (17 slayt, 20 Ekim 2026).
Not: Ekip slaytındaki yer tutucular ([bölüm], [Ekip arkadaşı · rolü]) kullanıcı tarafından editörde doldurulmuş olabilir; buradaki metin Claude'un son yazdığı sürümdür.

## Slayt 1 — `cover`

- YTÜ Startup House başvurusu
- Kuantum Triyaj
- Fotodinamik terapi ilaçları için kuantum-hibrit hesaplama platformu
- Klasik hesap her molekülde · kuantum yalnızca gereken bölgede
- Sinan Ahmet Yalçınkaya · Ekim 2026
- Temoporfin (Foscan) ve hesapladığımız kuantum bölgesi

**Konuşmacı notu:** Açılış: PDT ilaçlarının rengini ve ışığa tepkisini hesaplayan, kuantum bilgisayarı sadece gereken molekülde ve sadece kritik bölgede kullanan bir platform geliştiriyoruz. Görseldeki ilaç Avrupa'da onaylı temoporfin; renkli bulutlar bizim kuantumla hesapladığımız bölge.

## Slayt 2 — `problem`

- Problem
- PDT ilaçları hâlâ büyük ölçüde deneme-yanılmayla tasarlanıyor
- 20+ yıl
- Bu sürede rutin kliniğe giren, temelde yeni fotosensitizer sayısı çok az.
- 630–670 nm
- Klinikteki porfirin ilaçlarının ışığı. Derin tümör için 700–900 nm gerekiyor.
- 0,2–0,4 eV
- Sektör standardı TD-DFT'nin hatası. 700 nm'de bu, yaklaşık 100 nm kayma demek.
- Kaynaklar: Q4Bio perspektif makalesi (Algorithmiq, 2026) s.3 · ilaç ürün bilgileri (EMA, FDA) · ΔADAPT-VQE, arXiv:2404.16149, s.11 (6 BODIPY)

**Konuşmacı notu:** Üç sorun: (1) Yeni ilaç çıkmıyor, Algorithmiq'in Q4Bio makalesi 20 yılı aşkın sürede çok az temelde yeni fotosensitizerin rutin kliniğe girdiğini söylüyor. (2) Mevcut porfirin ilaçları 630–670 nm kırmızı ışıkla çalışıyor; ışık derine gitmiyor, doku penceresi 700–900 nm. (3) Yeni molekülü bilgisayarda tasarlamak için kullanılan TD-DFT 0,2–0,4 eV hata yapıyor; bu, 700 nm civarında yaklaşık 100 nm kayma: molekülün pencereye girip girmediğini söyleyemiyor.

## Slayt 3 — `neden`

- Neden kuantum
- Porfirinlerde ışık, iki elektron yolunu aynı anda kullanıyor
- 406 orbitalden 4'ü ışığı belirliyor (Gouterman modeli).
- Bu 4 orbital = 8 kübitlik kuantum bölgesi .
- Porfinde uyarılmış hal iki konfigürasyonun karışımı: %54 + %44 . Tek resimli yöntemler burada zorlanıyor.
- Kalan 402 orbital klasik olarak düzeltiliyor (NEVPT2).

**Konuşmacı notu:** Porfinin hesapta 406 orbitali var. Işığın rengini belirleyen sadece ortadaki 4 orbital. Porfinde ışık soğurulduğunda elektron iki farklı geçişi neredeyse eşit oranda (yüzde 54 ve 44) kullanıyor. Klasik tek-referanslı yöntemler bu karışımı iyi tarif edemiyor; kuantum bilgisayar ise bu çok-konfigürasyonlu durumu doğal olarak temsil ediyor. Biz bu 4 orbitali 8 kübitlik kuantum bölgesi olarak alıyoruz, gerisini klasik yöntemle düzeltiyoruz.

## Slayt 4 — `cozum`

- Çözüm
- Dört adımlı iş akışı
- 1
- Klasik ön hesap
- Geometri (GFN2-xTB) ve elektron yapısı (DFT). Her molekül için.
- 2
- Triyaj
- Kuantum bölgesindeki konfigürasyon karışımı ölçülür. İkinci ağırlık 0,2'yi geçerse: kuantum gerekli.
- 3
- Kuantum bölgesi
- VQE/VQD, 8 kübit. Bugün simülatörde; devre parametreleri IBM Quantum donanımına hazır.
- 4
- Düzeltme ve rapor
- NEVPT2 düzeltmesi. Dalga boyu, soğurma şiddeti (f), literatür karşılaştırması.
- Kuantum kaynağı kıt ve pahalı. Platform onu yalnızca gereken moleküle ve yalnızca kritik bölgeye harcıyor.

**Konuşmacı notu:** Akış: önce her molekül ucuz klasik yöntemlerle hesaplanıyor. Sonra triyaj adımı, kuantum bölgesindeki elektron konfigürasyonlarının ne kadar karıştığına bakıyor. Karışım yüksekse molekül kuantum hesabına gidiyor; düşükse klasik sonuç yeterli. Kuantum hesabı sadece 8 kübitlik kritik bölgede yapılıyor, geri kalan etkiler klasik NEVPT2 ile ekleniyor. Çıktı: ilacın hangi dalga boyunda ışık soğurduğu ve ne kadar güçlü soğurduğu, literatürle karşılaştırmalı.

## Slayt 5 — `parca`

- Bölge seçimi
- Bütün ilaçtan kuantum bölgesine
- 1 · Bütün ilaç
- Temoporfin, 84 atom
- 2 · Fotoaktif çekirdek
- Klorin halkası: ışığı soğuran kısım
- 3 · Kuantum bölgesi
- 8 kübit
- 4 orbital, VQE ile çözülür. Sonuç: 1,98 eV = 626 nm.
- Klinikte kullanılan ışık: 652 nm.

**Konuşmacı notu:** Fragmentasyon mantığı: ilacın tamamı 84 atom. Işığı soğuran kısım ortadaki klorin halkası; yan gruplar (hidroksifenil) halkaya neredeyse dik durduğu için rengi az etkiliyor. Ortadaki görselde renkli bulutlar bizim kuantumla hesapladığımız orbital. Bu bölge 4 orbital, yani 8 kübit. Sonuç 626 nm; klinikte 652 nm ışık kullanılıyor, fark 26 nm. Önemli not: hesap çekirdek üzerinde, yan gruplar ve çözücü henüz yok.

## Slayt 6 — `ilaclar`

- Gerçek ilaçlar
- Klinikte kullanılan üç ilaç, iki çekirdek
- Temoporfin
- Foscan · AB onaylı / Çekirdek: klorin · Işık: 652 nm
- Talaporfin
- Laserphyrin · Japonya onaylı / Çekirdek: klorin · Işık: 664 nm
- Protoporfirin IX
- 5-ALA'dan vücutta oluşur / Çekirdek: porfin · Işık: 635 nm
- İlaç çekirdekleri hesapladığımız modellere 0,22–0,40 Å içinde oturuyor · Işık: EMA (Foscan), Talaporfin 664–667 nm, FDA Ameluz etiketi (635 nm)

**Konuşmacı notu:** Üç gerçek ilaç: Temoporfin (Avrupa'da onaylı, baş-boyun kanseri), Talaporfin (Japonya'da onaylı) ve Protoporfirin IX (5-ALA kremi verildiğinde vücutta oluşan fotosensitizer; ABD'de Ameluz ile 635 nm ışıkla kullanılıyor). İlk ikisinin çekirdeği klorin, üçüncüsünün porfin. Gerçek ilaç yapılarını bizim hesapladığımız çekirdeğin üstüne oturttuğumuzda atomlar ortalama 0,22–0,40 Å kayıyor: ilacın içindeki halka, bizim hesapladığımız halkayla neredeyse aynı.

## Slayt 7 — `sonuclar`

- Doğrulama
- 8 porfirin molekülü, deneyle karşılaştırma
- Molekül
- Hesap (eV)
- Deney (eV)
- Fark
- Triyaj
- Porfin
- 1,87
- 1,98
- −0,11
- Kuantum gerekli
- Zn-porfin
- 2,07
- 2,18
- −0,11
- Kuantum gerekli
- Mg-porfin
- 2,04
- 2,17
- −0,13
- Kuantum gerekli
- Klorin
- 1,98
- 1,95
- +0,03
- Klasik yeterli
- Porfirazin
- 2,24
- 2,07
- +0,18
- Klasik yeterli
- Bakteriyoklorin *
- 2,23
- 1,74
- +0,49
- Model yetersiz
- Tetrametilporfin *
- 1,64
- 1,88
- −0,24
- Kuantum gerekli
- Pd-porfin
- 2,14
- –
- –
- Kuantum gerekli
- 0,11 eV
- ortalama hata, deneyi doğrudan ölçülmüş 5 molekül
- 0,18 eV
- deney değeri olan 7 molekülün tamamında
- S1 (Q bandı). Hesap: B3LYP orbitalleri + VQE (8 kübit) + NEVPT2, def2-SVP. * deney değeri analog/vekil molekülden. Deney kaynakları: Edwards 1971; Baerends 2002; Devillers 2013; Henchy 2014; Reddy 2017

**Konuşmacı notu:** Sekiz porfirin türevinde ilk uyarılmış halin enerjisini deneyle karşılaştırdık. Deneyi doğrudan ölçülmüş beş molekülde ortalama hata 0,11 eV; analog/vekil deney değerleri dahil yedi molekülde 0,18 eV. En büyük hata bakteriyoklorinde: 4 orbitalli kuantum bölgesi bu molekül için yetersiz, 16 kübitlik hesap çalışıyor. Tetrametilporfinin geometrisi şüpheli. Bunları saklamıyoruz; bir sonraki slayttaki triyaj da bu sınırları gösteriyor.

## Slayt 8 — `yontem`

- Klasik yöntemle karşılaştırma
- Porfin ve klorinde TD-DFT'den 5 kat yakın

**Konuşmacı notu:** Sektör standardı TD-DFT B3LYP'yi aynı geometri ve aynı baz setiyle kendi bilgisayarımızda hesapladık. Porfinde TD-DFT 2,35 eV veriyor, deney 1,98: 0,37 eV fazla. Klorinde 2,37'ye karşı 1,95: 0,42 eV fazla. Bizim yöntem aynı iki molekülde ortalama 0,07 eV hata yapıyor. Diğer moleküllerin TD-DFT hesabı sürüyor; şimdilik karşılaştırma iki molekülde. Dürüst not: farkın kaynağı çok-konfigürasyonlu yaklaşım (kuantum bölgesi + NEVPT2); bu ölçekte o hesap klasik bilgisayarda da yapılabiliyor. Ayrıca triplet hal için TD-DFT bizim eski yöntemimizden daha iyi çıktı (porfin T1: TD-DFT 1,83, bizim ΔSCF 2,00, deney 1,58); bunu yol haritasına aldık.

## Slayt 9 — `donanim`

- Gerçek kuantum donanımı
- IBM Kingston'da ilk ölçüm: gürültü baskın
- Porfin, S1 (8 kübitlik kuantum bölgesi)
- Hesap
- S1 (eV)
- Tam çözümden fark
- Tam çözüm (CASCI, klasik)
- 3,62
- –
- VQE, ideal simülatör
- 3,62
- 0,00
- VQE, ibm_kingston, ham
- 1,14
- −2,48
- VQE, ibm_kingston, gürültü düzeltmeli
- 4,55
- +0,93 ± 0,19
- %41–53
- Sinyalin donanım gürültüsünden sağ çıkan kısmı
- Bugün 8 kübitlik bölge bile hata düzeltmesiz donanımda kimyasal doğruluğa ulaşmıyor.
- Triyajın değeri burada: kuantum kaynağını yalnızca gereken moleküle harcamak ve donanım olgunlaştıkça aynı iş akışını büyütmek.
- İş dbbl9ns24mrc73chvt8g · ibm_kingston (156 kübit) · 20 Ekim 2026 · 4 devre × 52 ölçüm grubu × 4096 atış · TREX, twirling, dinamik dekuplaj, gürültü tahmin devresi

**Konuşmacı notu:** Aynı VQE devresini IBM'in 156 kübitlik Kingston işlemcisinde çalıştırdık. İdeal simülatörde devre tam çözümü (3,62 eV) birebir veriyor. Gerçek donanımda ham sonuç 1,14 eV, yani gürültü sonucu tamamen bozuyor. Gürültü tahmin devresiyle düzeltince 4,55 eV, tam çözümden 0,93 eV uzak; istatistik belirsizlik yaklaşık 0,19 eV. Gürültü faktörü 0,41 ve 0,53: sinyalin yarısından fazlası kayboluyor. Bu sonuç, alandaki büyük oyuncuların da söylediğini doğruluyor: bugün kuantum avantajı yok. Sonraki adımlar: daha sığ devre, sıfır-gürültü ekstrapolasyonu (ZNE) ve en iyi kübitlerin seçimi.

## Slayt 10 — `triyaj`

- Triyaj
- Kuantum hangi molekülde gerekli?
- Kuantum gerekli
- İkinci konfigürasyon ağırlığı 0,36–0,47
- Porfin · Zn-porfin · Mg-porfin · Pd-porfin · Tetrametilporfin
- Simetrik halka: iki geçiş eş enerjili, güçlü karışım.
- Klasik yeterli
- İkinci konfigürasyon ağırlığı 0,10–0,17
- Klorin · Porfirazin · (Bakteriyoklorin)
- Simetri bozuk: tek geçiş baskın. Bakteriyoklorin büyük aktif uzayla yeniden test ediliyor.
- Farklı alan, aynı kural
- β-siklodekstrin + siklopentanol (ilaç taşıyıcı)
- Klasik MD: ΔG −2,7 ± 0,6 / Deney: −3,05 kcal/mol
- Etkileşimin %96'sı dispersiyon. Triyaj: kuantum gereksiz.

**Konuşmacı notu:** Triyaj kuralı: kuantum bölgesindeki uyarılmış halde ikinci elektron konfigürasyonunun ağırlığı 0,2'yi geçerse molekül kuantum hesabına gidiyor. Simetrik porfirinlerde ağırlık 0,36–0,47; klorin ve porfirazinde 0,10–0,17. Bu, 1960'lardan beri bilinen Gouterman modeliyle uyumlu. Aynı yaklaşımı ilaç taşıyıcı siklodekstrin sisteminde de denedik: klasik moleküler dinamik deneyi hata payı içinde yakaladı, etkileşim neredeyse tamamen dispersiyon; triyaj doğru şekilde kuantum gereksiz dedi. Not: bakteriyoklorinde kural yanlış sınıflandırdı, çünkü 4 orbital bu molekül için yetersiz; 16 kübitlik hesap sürüyor.

## Slayt 11 — `klinik`

- Klinik ışıkla karşılaştırma
- Biz 26–38 nm, TD-DFT 107–141 nm uzakta

**Konuşmacı notu:** Siyah elmaslar klinikte kullanılan ışık, mavi noktalar bizim hesabımız, turuncu kareler TD-DFT. TD-DFT klorin çekirdeğini 523 nm, porfini 528 nm buluyor: klinik ışıktan 107–141 nm uzak. Temoporfin 626 nm vs 652 nm, Talaporfin 626 vs 664, PpIX 663 vs 635: üçü de 26–38 nm içinde, ve bunlar yan grupları ve çözücüyü içermeyen çekirdek hesapları. Dördüncü ilaç padeliporfin (Tookad, 753 nm) yeni nesil yakın kızılötesi ilaçların örneği; çekirdeği bakteriyoklorin ve 4 orbitalli modelimiz burada yetersiz kalıyor. En değerli pazar bu bölge; büyük aktif uzay hesabı tam bu boşluğu kapatmak için.

## Slayt 12 — `rekabet`

- Konumlanma
- Liderler BODIPY'de, porfirinler açık alan
- Xanadu (Kanada)
- Algorithmiq (Finlandiya)
- Kuantum Triyaj
- Molekül ailesi
- BODIPY
- BODIPY, Ru kompleksi (TLD1433)
- Porfirin, klorin: klinik ilaçların çoğu
- Yaklaşım
- Gelecekteki donanım: 180–350 mantıksal kübit
- 8–16 kübit VQE; IBM'de 52–100 kübit + DMRG
- 8 kübitlik bölge + triyaj + NEVPT2
- Deneyle doğrulama
- Yok (kaynak tahmini)
- BODIPY: 0,13 eV (simülatör)
- Porfirinler: 0,11 eV (simülatör)
- Durum
- U. Alberta ortaklığı, 2026
- Q4Bio 2 M$ ödülü, 2026
- Prototip; IBM Kingston'da ilk çalıştırma
- Kanıtlanabilir kuantum avantajı için konan 5 M$'lık büyük ödül kimseye verilmedi : avantaj henüz kimsede yok.
- Kaynaklar: arXiv:2512.15889 (Xanadu) · arXiv:2404.16149 ve Q4Bio perspektifi (Algorithmiq) · Quantum Computing Report, Nisan 2026

**Konuşmacı notu:** Bu alanda iki güçlü oyuncu var. Xanadu, gelecekteki hataya dayanıklı kuantum bilgisayarlar için BODIPY ilaçlarında kaynak tahmini yaptı; hesaplanmış enerji yok. Algorithmiq, BODIPY'de simülatörde deneye 0,13 eV yaklaştı ve bir rutenyum ilacında IBM donanımında 100 kübite kadar çalıştı; 2 milyon dolarlık Q4Bio ödülünü aldı. Ama kuantum avantajı ödülü olan 5 milyon dolar kimseye verilmedi. İkisi de porfirin ailesine bakmıyor; klinikteki PDT ilaçlarının çoğu ise bu aileden. Biz oradayız. Dürüst not: Algorithmiq'in makalesinde klasik ΔCCSD 0,05 eV ile kuantum yöntemini de geçiyor.

## Slayt 13 — `iddia`

- Dürüst çerçeve
- Ne iddia ediyoruz, ne etmiyoruz
- İddia ediyoruz
- Porfirin ailesinde çalışan iş akışı: 5 molekülde deneye 0,11 eV; porfin ve klorinde TD-DFT 0,39 eV, biz 0,07 eV
- 8 kübitlik VQE simülatörde tam çözümle aynı; IBM Kingston'da gerçek donanımda çalıştı
- Triyaj, kuantumun gerektiği molekülleri tutarlı ayırıyor
- Üç klinik ilacın çekirdeğinde klinik ışığa 26–38 nm
- İddia etmiyoruz
- Kuantum avantajı: 8 kübit klasikte de çözülüyor; gerçek donanımda gürültü hatası ~0,9 eV
- İlacın tamamı: yan gruplar ve çözücü henüz yok
- Yakın kızılötesi ilaçlar: bakteriyoklorinde 0,49 eV hata
- Triplet hal ve ISC: henüz güvenilir değil

**Konuşmacı notu:** Bu slayt jürinin soracağı soruları baştan cevaplıyor. Bugün kuantum bilgisayar bize doğruluk avantajı vermiyor; 8 kübitlik hesap dizüstü bilgisayarda da yapılabiliyor ve iyi sonucun önemli kısmı klasik NEVPT2 düzeltmesinden geliyor. Bizim değerimiz, kuantum donanımı büyüdükçe aynı iş akışının daha büyük kuantum bölgesiyle çalışacak olması ve bugünden hangi molekülün kuantuma ihtiyaç duyduğunu söyleyebilmesi. Alanın liderleri de avantajın henüz kanıtlanmadığını söylüyor.

## Slayt 14 — `model`

- İş modeli
- Sentezden önce: hangi ışık, ne kadar güçlü, kuantuma değer mi?
- PDT ve fotokimya laboratuvarları
- Sentezden önce aday eleme, yayın için hesap desteği.
- İlaç ve biyoteknoloji Ar-Ge
- Fotosensitizer optimizasyonu, yakın kızılötesi yeni adaylar.
- Komşu pazarlar
- Boya duyarlı güneş pilleri, fotokataliz, OLED: aynı molekül ailesi, aynı soru.
- Analiz raporu
- Molekül başına ücret
- Web platformu
- Abonelik: SMILES gir, rapor al
- Ortak Ar-Ge
- Proje ortaklıkları ve hibeler

**Konuşmacı notu:** Müşteri soruyu sentezden önce soruyor: bu molekül hangi dalga boyunda ışık soğurur, ne kadar güçlü soğurur, ve pahalı kuantum hesabına değer mi? Üç müşteri grubu: akademik PDT ve fotokimya laboratuvarları, ilaç ve biyoteknoloji Ar-Ge ekipleri, ve aynı molekül ailesini kullanan güneş pili, fotokataliz ve OLED alanları. Gelir: molekül başına analiz raporu, abonelikli web platformu ve ortak Ar-Ge projeleri. Fiyatlandırma ilk müşteri görüşmelerinden sonra belirlenecek.

## Slayt 15 — `yol`

- Yol haritası
- 12 ayda prototipten ilk müşteriye
- Bu ay
- TD-DFT referansını tamamla (ilk 2 molekül bitti)
- 16 kübitlik bakteriyoklorin hesabı
- IBM'de gürültü azaltma: ZNE, daha sığ devre
- 1–3 ay
- 5 PDT laboratuvarıyla görüşme
- Yan grup ve çözücü etkisi
- Triplet hal için daha iyi yöntem
- 3–6 ay
- Web demo: SMILES gir, rapor al
- 20+ molekül tarama
- İlk pilot kullanıcı
- 6–12 ay
- TÜBİTAK BiGG başvurusu
- İlk ücretli analizler
- Bilimsel yayın

**Konuşmacı notu:** Bu ay: klasik TD-DFT referansı (şu an çalışıyor), bakteriyoklorin için 16 kübitlik büyük aktif uzay hesabı ve IBM Quantum donanımında porfin VQE denemesi; parametreler hazır. 1–3 ay: beş PDT laboratuvarıyla müşteri görüşmesi, yan grup ve çözücü etkilerinin eklenmesi, triplet hal için daha iyi yöntem. 3–6 ay: web demo ve 20'den fazla molekül tarama, ilk pilot kullanıcı. 6–12 ay: TÜBİTAK BiGG başvurusu, ilk ücretli analizler ve bilimsel yayın.

## Slayt 16 — `talep`

- Ekip ve beklenti
- Startup House'tan beklentimiz
- Mentorluk
- İlaç Ar-Ge, iş geliştirme, fiyatlandırma
- Hesaplama altyapısı
- Bulut kredisi, kuantum donanım erişimi
- Sektör bağlantıları
- Pilot kullanıcı olacak laboratuvar ve firmalar
- Ofis ve yatırımcı ağı
- Demo günleri, ilk yatırım görüşmeleri
- Ekip
- Sinan Ahmet Yalçınkaya / Kurucu · [bölüm], YTÜ
- Kuantum programlama (PennyLane), moleküler simülasyon (GROMACS, PySCF), biyomedikal optik
- [Ekip arkadaşı · rolü]
- Aranan: kimyager / ilaç Ar-Ge danışmanı

**Konuşmacı notu:** Startup House'tan dört şey istiyoruz: ilaç Ar-Ge ve iş geliştirme mentorluğu, hesaplama altyapısı (bulut kredisi ve kuantum donanım erişimi), pilot kullanıcı olacak laboratuvar ve firmalarla bağlantı, ve ofis ile yatırımcı ağı. Ekipte en büyük eksiğimiz bir kimyager veya ilaç Ar-Ge danışmanı; bunu açıkça arıyoruz.

## Slayt 17 — `kaynak`

- Kuantum Triyaj
- Doğru molekül, sentezden önce.
- [iletisim e-postasi]
- Kaynaklar
- Zhou vd. (Xanadu), arXiv:2512.15889, 2025
- Nykänen vd. (Algorithmiq), arXiv:2404.16149, 2024
- Algorithmiq Q4Bio perspektifi, 2026; QCR, Nisan 2026
- Edwards vd., J. Mol. Spectrosc. 38, 16 (1971)
- Baerends vd., Coord. Chem. Rev. 230, 5 (2002)
- Devillers vd., C. R. Chimie 16, 540 (2013)
- Henchy vd., J. Chem. Phys. 141, 124303 (2014)
- Reddy vd., Molecules 22, 634 (2017)
- Ürün bilgileri: Foscan, Tookad (EMA); Ameluz (FDA)

**Konuşmacı notu:** Kapanış: Kuantum Triyaj, ışıkla çalışan ilaçlarda doğru molekülü sentezden önce bulmak için. Teşekkürler, sorularınızı bekliyorum.
