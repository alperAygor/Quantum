# IBM kuantum donanımında VQE: plan, ön test ve İLK SONUÇ (20 Ekim 2026)

## İlk gerçek donanım sonucu (porfin, S1, 8 kübit)
- İş kimliği: dbbl9ns24mrc73chvt8g. Donanım: ibm_kingston (156 kübit). Gönderim: 20 Ekim 2026, 14:47. Instance: open-instance, Job modu.
- Ayarlar: 4 devre × 52 ölçüm grubu × 4096 atış. Gürültü azaltma: TREX, kapı twirling, dinamik dekuplaj, gürültü tahmin devresi.

| Hesap | S1 (eV) | Tam çözümden fark |
|---|---|---|
| Tam çözüm (CASCI, klasik) | 3,617 | – |
| VQE, ideal simülatör | 3,616 | −0,001 |
| VQE, ibm_kingston, ham | 1,141 | −2,477 |
| VQE, ibm_kingston, düzeltilmiş | 4,546 | +0,928 |

- İstatistik belirsizlik yaklaşık 0,19 eV.
- Gürültü faktörü f: S0 için 0,41, S1 için 0,53. Sinyalin yarısından fazlası gürültüde kayboluyor.
- Gerçek donanım, FakeKingston simülasyonundan daha gürültülü çıktı:

| | Simülasyon tahmini | Gerçek donanım |
|---|---|---|
| Ham fark | −0,79 | −2,48 |
| Düzeltilmiş fark | +0,48 | +0,93 |

- Yorum: Hata düzeltmesi olmayan bugünkü donanımda 8 kübitlik bölge bile kimyasal doğruluğa ulaşmıyor. Sunuma "IBM Kingston'da ilk ölçüm: gürültü baskın" slaytı eklendi; kuantum avantajı iddia edilmiyor.
- Sonraki adımlar (kalan Open Plan süresine göre):
  - --zne seçeneği
  - Daha sığ devre: daha az kapı, qubit tapering
  - En iyi kübitlerin seçilmesi
  - Ham sonuçta parçacık sayısı ile post-seleksiyon

## Kurulum notları (kullanıcının WSL ortamı)
- Hata: libgomp "cannot allocate memory in static TLS block".
  - Sebep: pyscf, qiskit/pennylane'den sonra yüklenince kendi libgomp kütüphanesini yükleyemiyor.
  - Çözüm: ibm_vqe.py'de pyscf, numpy'den hemen sonra import ediliyor (sed ile düzeltildi).
- qiskit_aer da aynı nedenle yüklenemiyor. Bu yüzden kullanıcının bilgisayarında gürültü simülasyonu ve --sahte modu çalışmıyor; gerçek donanım modu çalışıyor.
- Web: quantum.cloud.ibm.com/workloads (cloud.ibm.com değil). Girişte MFA kaydı gerekiyor.

## Paket ve yöntem
- v5 paketi: ibm_vqe.py (komutları: hazirla, gonder, sonuc; seçenekleri: --sahte, --zne, --backend=, --shots=), hesap_kaydet.py.
- Devre: kompakt ADAPT, S0 için 8 ve S1 için 5 çift uyarma kapısı. Parametreler simülatörde optimize ediliyor; donanımda yalnızca enerji ölçülüyor.
- Gözlenebilir: |c| ≥ 1e-3 olan 121 Pauli terimi, 52 QWC grubu. Kesme hatası yaklaşık 0,014 eV.
- Gürültü tahmin devresi: tüm açılar 0 alınıyor, depolarizan ölçekleme yapılıyor:
  E_duz = c0 + (E − c0)·(E_ref_ideal − c0)/(E_ref_ölçülen − c0)
