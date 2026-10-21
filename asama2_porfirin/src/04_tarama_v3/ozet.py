# Tum *_sonuc.json dosyalarini tek tabloda topla
import json, glob, csv
rows = [json.load(open(f)) for f in sorted(glob.glob("*_v3_sonuc.json"))]
keys = ["molekul", "atom", "triyaj", "S1_konfig_agirlik", "S1_CAS_eV", "S1_NEVPT2_eV", "T1_CAS_eV", "T1_dSCF_eV",
        "f_S1", "f_S2", "singlet_oksijen", "T1_gt_S1_hata", "VQD_CASCI_maxfark_eV", "metal_payi", "sure_s"]
with open("tarama_sonuclari.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(keys); [wr.writerow([r.get(k, "") for k in keys]) for r in rows]
print(f"{'molekul':17s} {'triyaj':16s} {'S1 CAS':>7s} {'S1 +NEV':>8s} {'T1 dSCF':>8s} {'f(S1)':>6s} {'1O2':>4s}")
for r in rows:
    g = lambda k, f="{:8.2f}": f.format(r[k]) if isinstance(r.get(k), (int, float)) else f"{'-':>8s}"
    print(f"{r['molekul']:17s} {r['triyaj']:16s} {r['S1_CAS_eV']:7.2f} {g('S1_NEVPT2_eV')} {g('T1_dSCF_eV')} {r['f_S1']:6.3f} {str(r.get('singlet_oksijen','-')):>4s}")
print("Kaydedildi: tarama_sonuclari.csv")
