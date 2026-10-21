# Final figurler + ozet tablo. ~/tarama klasorunde calistir:  python3 final_figurler.py
# Okur: *_v3_sonuc.json (bizim: CAS(4,4)/VQE + NEVPT2), *_tddft.json (klasik referans), *_cas88.json (opsiyonel)
# Yazar: F1_yontem_karsilastirma.png, F2_triyaj_dogrulama.png, F3_klinik_ilaclar.png, final_tablo.csv, final_ozet.txt
import json, glob, os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3de"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})

AD = {"porfin": "Porfin", "klorin": "Klorin", "bakteriyoklorin": "Bakteriyoklorin", "porfirazin": "Porfirazin",
      "Zn-porfin": "Zn-porfin", "Mg-porfin": "Mg-porfin", "Pd-porfin": "Pd-porfin", "tetrametilporfin": "Tetrametilporfin"}
# Deneysel S1 (Q bandi, eV). dogrudan=True: ayni molekulun olcumu; False: analog/vekil molekul
DENEY = {"porfin": (1.98, True, "Edwards 1971, buhar"), "klorin": (1.95, True, "n-alkan, 4.2 K"),
         "Zn-porfin": (2.18, True, "Baerends 2002"), "Mg-porfin": (2.17, True, "Devillers 2013"),
         "porfirazin": (2.065, True, "Henchy 2014, Ne matris"),
         "bakteriyoklorin": (1.74, False, "analog B6, Reddy 2017"),
         "tetrametilporfin": (1.88, False, "meso-tetraalkil vekil")}
# Klinik ilaclar: (isik nm, cekirdek modeli, kaynak)
KLINIK = [("Temoporfin\n(Foscan)", 652, "klorin", "EMA SmPC"),
          ("Talaporfin\n(Laserphyrin)", 664, "klorin", "664–667 nm"),
          ("PpIX\n(5-ALA, Ameluz)", 635, "porfin", "FDA etiketi"),
          ("Padeliporfin\n(Tookad)", 753, "bakteriyoklorin", "EMA SmPC")]
nm = lambda e: 1239.84 / e

V3 = {json.load(open(f))["molekul"]: json.load(open(f)) for f in glob.glob("*_v3_sonuc.json")}
TD = {json.load(open(f))["molekul"]: json.load(open(f)) for f in glob.glob("*_tddft.json")}
C8 = {json.load(open(f))["molekul"]: json.load(open(f)) for f in glob.glob("*_cas88.json")}
if not V3: raise SystemExit("*_v3_sonuc.json bulunamadi - ~/tarama klasorunde calistir")
ORDER = [m for m in ["porfin", "Zn-porfin", "Mg-porfin", "Pd-porfin", "tetrametilporfin", "klorin", "porfirazin", "bakteriyoklorin"] if m in V3]

def td(m, xc, k="S1_eV"):
    return TD.get(m, {}).get(xc, {}).get(k)
YONTEM = [("Bizim: VQE(8 kübit)+NEVPT2", BLUE, "o", lambda m: V3[m].get("S1_NEVPT2_eV")),
          ("TD-DFT B3LYP", ORANGE, "s", lambda m: td(m, "b3lyp")),
          ("TD-DFT CAM-B3LYP", AQUA, "^", lambda m: td(m, "camb3lyp"))]
YONTEM = [y for y in YONTEM if any(y[3](m) is not None for m in ORDER)]

# ---------- MAE ----------
ORTAK = [m for m in ORDER if m in DENEY and all(y[3](m) is not None for y in YONTEM)]   # adil kiyas: tum yontemlerin oldugu molekuller
def mae(fn, sadece_dogrudan):
    err = [abs(fn(m) - DENEY[m][0]) for m in ORTAK if (DENEY[m][1] or not sadece_dogrudan)]
    return (float(np.mean(err)), len(err)) if err else (None, 0)
ozet = []
for ad, _, _, fn in YONTEM:
    a, na = mae(fn, False); b, nb = mae(fn, True)
    ozet.append((ad, a, na, b, nb))

# ---------- F1: yontem karsilastirmasi ----------
fig, (ax, bx) = plt.subplots(1, 2, figsize=(14, 6.2), gridspec_kw={"width_ratios": [2.2, 1]})
ys = np.arange(len(ORDER))[::-1]
for y, m in zip(ys, ORDER):
    ax.axhline(y, color=GRID, lw=1, zorder=0)
    if m in DENEY:
        ax.scatter(DENEY[m][0], y, marker="D", s=90, color=INK, zorder=4, edgecolor="white", linewidth=2)
    for k, (ad, c, mk, fn) in enumerate(YONTEM):
        v = fn(m)
        if v is not None:
            ax.scatter(v, y + (k - 1) * 0.18, marker=mk, s=80, color=c, zorder=3, edgecolor="white", linewidth=2)
lab = [AD[m] + ("" if m not in DENEY else ("" if DENEY[m][1] else " *")) + ("" if m in DENEY else " (deney yok)") for m in ORDER]
ax.set_yticks(ys); ax.set_yticklabels(lab); ax.set_xlabel("S$_1$ uyarılma enerjisi (eV)")
ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
h = [plt.Line2D([], [], marker="D", ls="", color=INK, ms=8, label="Deney")] + \
    [plt.Line2D([], [], marker=mk, ls="", color=c, ms=8, label=ad) for ad, c, mk, _ in YONTEM]
ax.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False, fontsize=10)
ax.text(0, -0.13, "* deney değeri analog/vekil molekülden", transform=ax.transAxes, color=MUTED, fontsize=9)
# MAE paneli (adil: tum yontemlerin hesaplandigi ortak molekul seti)
iki = any(o[2] != o[4] for o in ozet)
w = 0.38 if iki else 0.6; x = np.arange(len(ozet))
for i, (ad, a, na, b, nb) in enumerate(ozet):
    c = YONTEM[i][1]
    if iki and a is not None:
        bx.bar(i - w / 2 - 0.01, a, w, color=c, alpha=0.45); bx.text(i - w / 2, a + 0.008, f"{a:.2f}", ha="center", fontsize=10, color=INK)
    if b is not None:
        xb = i + w / 2 + 0.01 if iki else i
        bx.bar(xb, b, w, color=c); bx.text(xb, b + 0.008, f"{b:.2f} eV", ha="center", fontsize=10, color=INK)
bx.set_xticks(x); bx.set_xticklabels([o[0].replace(": ", ":\n").replace("TD-DFT ", "TD-DFT\n") for o in ozet], fontsize=9)
bx.set_ylabel("Deneye göre ortalama mutlak hata (eV)"); bx.grid(axis="y", color=GRID); bx.set_axisbelow(True)
bx.set_title("Aynı moleküllerde (n=%d): %s" % (len(ORTAK), ", ".join(AD[m] for m in ORTAK)), fontsize=10, color=INK, pad=12)
if iki:
    n_all = max(o[2] for o in ozet); n_dir = max(o[4] for o in ozet)
    bx.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=MUTED, alpha=0.45, label=f"tüm deneyli (n={n_all})"),
                       plt.Rectangle((0, 0), 1, 1, color=MUTED, label=f"doğrudan deney (n={n_dir})")], loc="upper left", frameon=False, fontsize=9)
fig.suptitle("Porfirin ailesi: S$_1$ (Q bandı) — kuantum-hibrit yöntem, klasik TD-DFT ve deney", fontsize=14, color=INK, y=1.02)
fig.tight_layout(); fig.savefig("F1_yontem_karsilastirma.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ---------- F2: triyaj dogrulamasi ----------
fig, ax = plt.subplots(figsize=(10, 6))
ax.axvspan(0.2, 0.5, color=BLUE, alpha=0.06); ax.axvline(0.2, color=MUTED, ls="--", lw=1.2)
ax.text(0.205, 0.97, "KUANTUM GEREKLİ bölgesi →", transform=ax.get_xaxis_transform(), color=INK, fontsize=10, va="top")
ax.text(0.195, 0.97, "← klasik yeterli", transform=ax.get_xaxis_transform(), color=INK, fontsize=10, va="top", ha="right")
for m in ORDER:
    if m not in DENEY: continue
    xw = V3[m]["S1_konfig_agirlik"][1]
    eb = abs(V3[m]["S1_NEVPT2_eV"] - DENEY[m][0]); et = td(m, "b3lyp")
    if et is not None:
        et = abs(et - DENEY[m][0]); ax.plot([xw, xw], [eb, et], color=GRID, lw=2, zorder=1)
        ax.scatter(xw, et, marker="s", s=80, color=ORANGE, edgecolor="white", linewidth=2, zorder=3)
    ax.scatter(xw, eb, marker="o", s=80, color=BLUE, edgecolor="white", linewidth=2, zorder=3)
    ax.annotate(AD[m] + ("" if DENEY[m][1] else " *"), (xw, max(eb, et or 0)), xytext=(6, 6), textcoords="offset points", fontsize=9, color=INK)
ax.set_xlim(0.05, 0.5); ax.set_ylim(0, None)
ax.set_xlabel("S$_1$'de ikinci elektron konfigürasyonunun ağırlığı (8 kübitlik VQE/CASCI'den)")
ax.set_ylabel("Deneye göre mutlak hata (eV)"); ax.grid(color=GRID); ax.set_axisbelow(True)
h = [plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=8, label="Bizim (VQE + NEVPT2)")]
if TD: h.append(plt.Line2D([], [], marker="s", ls="", color=ORANGE, ms=8, label="TD-DFT B3LYP"))
ax.legend(handles=h, loc="upper left", bbox_to_anchor=(0, 0.92), frameon=False)
ax.set_title("Triyaj: konfigürasyon karışımı arttıkça tek-referanslı klasik yöntem zorlanıyor mu?", fontsize=13, color=INK)
fig.text(0.01, -0.02, "* deney değeri analog/vekil molekülden", color=MUTED, fontsize=9)
fig.tight_layout(); fig.savefig("F2_triyaj_dogrulama.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ---------- F3: klinik ilaclar (cekirdek modeli) ----------
fig, ax = plt.subplots(figsize=(12, 5.6))
ax.axvspan(700, 900, color=AQUA, alpha=0.10)
ax.text(705, -0.75, "doku geçirgenlik penceresi (700–900 nm)", color=INK, fontsize=9, va="bottom")
ys = np.arange(len(KLINIK))[::-1]
for y, (ad, lam, core, kay) in zip(ys, KLINIK):
    ax.axhline(y, color=GRID, lw=1, zorder=0)
    ax.scatter(lam, y, marker="D", s=100, color=INK, zorder=4, edgecolor="white", linewidth=2)
    ax.annotate(f"{lam} nm", (lam, y), xytext=(0, 10), textcoords="offset points", ha="center", fontsize=9, color=INK)
    if core in V3:
        v = nm(V3[core]["S1_NEVPT2_eV"]); ax.scatter(v, y - 0.12, marker="o", s=90, color=BLUE, zorder=3, edgecolor="white", linewidth=2)
        tri = "kuantum gerekli" if V3[core]["triyaj"].startswith("KUANTUM") else "klasik yeterli"
        ek = " (4 orbitalli model yetersiz: 16 kübitlik hesap gerekli)" if core == "bakteriyoklorin" and core not in C8 else ""
        ax.annotate(f"çekirdek: {AD[core].lower()} · triyaj: {tri}{ek}", (505, y - 0.36), fontsize=8.5, color=MUTED)
    if td(core, "b3lyp") is not None:
        ax.scatter(nm(td(core, "b3lyp")), y + 0.12, marker="s", s=80, color=ORANGE, zorder=3, edgecolor="white", linewidth=2)
    if core in C8 and C8[core].get("S1_NEVPT2_eV"):
        ax.scatter(nm(C8[core]["S1_NEVPT2_eV"]), y - 0.12, marker="o", s=90, facecolor="white", edgecolor=BLUE, linewidth=2, zorder=3)
ax.set_yticks(ys); ax.set_yticklabels([k[0] for k in KLINIK]); ax.set_xlim(500, 800); ax.set_ylim(-0.8, len(KLINIK) - 0.5)
ax.set_xlabel("Dalga boyu (nm)"); ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
h = [plt.Line2D([], [], marker="D", ls="", color=INK, ms=8, label="Klinikte kullanılan ışık"),
     plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=8, label="Bizim (çekirdek, VQE+NEVPT2)")]
if TD: h.append(plt.Line2D([], [], marker="s", ls="", color=ORANGE, ms=8, label="TD-DFT B3LYP (çekirdek)"))
if C8: h.append(plt.Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=BLUE, ms=8, label="Bizim, büyük aktif uzay (16 kübit)"))
ax.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=10)
fig.text(0.01, -0.03, "Hesaplar ilacın fotoaktif çekirdeği üzerinde (yan gruplar ve çözücü yok). Işık kaynakları: EMA/FDA ürün bilgileri.", color=MUTED, fontsize=9)
fig.tight_layout(); fig.savefig("F3_klinik_ilaclar.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ---------- tablo + ozet ----------
with open("final_tablo.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["molekul", "deney_S1", "deney_tip", "bizim_S1", "TD_B3LYP_S1", "TD_CAM_S1", "f_S1_bizim", "f_S1_TDB3LYP",
                 "ikinci_konfig_agirlik", "triyaj", "T1_dSCF", "T1_TDA_B3LYP", "T1_TDA_CAM", "CAS88_S1_NEVPT2", "CAS88_agirlik"])
    for m in ORDER:
        d = DENEY.get(m, (None, None, ""))
        wr.writerow([m, d[0], "dogrudan" if d[1] else ("analog/vekil" if d[1] is False else "yok"), V3[m].get("S1_NEVPT2_eV"),
                     td(m, "b3lyp"), td(m, "camb3lyp"), V3[m].get("f_S1"), td(m, "b3lyp", "f_S1"), V3[m]["S1_konfig_agirlik"][1],
                     V3[m]["triyaj"], V3[m].get("T1_dSCF_eV"), td(m, "b3lyp", "T1_eV"), td(m, "camb3lyp", "T1_eV"),
                     C8.get(m, {}).get("S1_NEVPT2_eV"), C8.get(m, {}).get("S1_konfig_agirlik")])
L = ["ORTALAMA MUTLAK HATA (eV), deneye gore"]
for ad, a, na, b, nb in ozet:
    L.append(f"  {ad:32s} tum: {a if a is None else round(a,3)} (n={na})   dogrudan deney: {b if b is None else round(b,3)} (n={nb})")
L.append("\nT1 (eV) - porfin deney 1.58 (fosforesans)")
for m in ORDER:
    L.append(f"  {AD[m]:17s} dSCF {V3[m].get('T1_dSCF_eV', 0):.2f}   TDA-B3LYP {td(m,'b3lyp','T1_eV') or float('nan'):.2f}   TDA-CAM {td(m,'camb3lyp','T1_eV') or float('nan'):.2f}")
if C8:
    L.append("\nBUYUK AKTIF UZAY (8e,8o) = 16 kubit")
    for m, r in C8.items():
        L.append(f"  {AD.get(m,m):17s} S1 CAS {r.get('S1_CAS_eV')}  +NEVPT2 {r.get('S1_NEVPT2_eV')}  agirliklar {r.get('S1_konfig_agirlik')}  f {r.get('f_S1')}  {r.get('triyaj','')}  deney {DENEY.get(m,('-',))[0]}")
open("final_ozet.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
print("\nYazildi: F1_yontem_karsilastirma.png F2_triyaj_dogrulama.png F3_klinik_ilaclar.png final_tablo.csv final_ozet.txt")
