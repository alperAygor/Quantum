import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
BLUE, ORANGE, INK, MUTED, GRID, BG, WARN = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a64", "#e4e3df", "#fcfcfb", "#9a5b00"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
# ad, S1 CAS, S1 NEV, S1 deney, deney notu, f, agirliklar, T1 (deger, etiket) veya None, uyari
M = [
 ("Porfin", 3.62, 1.87, 1.98, "gaz fazı (Edwards 1971)", 0.033, (0.544, 0.435), (1.58, "T1  deney (fosforesans)"), ""),
 ("Klorin", 3.45, 1.98, 1.95, "n-alkan, 4.2 K", 0.579, (0.743, 0.098), (1.53, "T1  hesap (ΔSCF, yaklaşık)"), ""),
 ("Bakteriyoklorin", 3.59, 2.23, 1.74, "sübstitüe analog, toluen", 0.422, (0.748, 0.168), (1.65, "T1  hesap (ΔSCF, yaklaşık)"),
  "Sınır: 4 orbital bu molekül için yetersiz; hesap deneyden 0.49 eV uzak. Daha büyük aktif uzay (daha fazla orbital) gerekiyor."),
 ("Zn-porfin", 3.75, 2.07, 2.18, "çözelti (Baerends 2002)", 0.079, (0.616, 0.367), None, ""),
 ("Mg-porfin", 3.72, 2.04, 2.17, "DMF (Devillers 2013)", 0.073, (0.611, 0.373), None, ""),
 ("Pd-porfin", 3.80, 2.14, None, "deney bulunamadı", 0.085, (0.629, 0.355), None, ""),
 ("Porfirazin", 3.74, 2.24, 2.07, "Ne matrisi, 8 K", 0.633, (0.788, 0.139), (1.65, "T1  hesap (ΔSCF, yaklaşık)"), ""),
 ("Tetrametilporfin", 3.42, 1.64, 1.88, "tetraalkil vekil, CH2Cl2", 0.001, (0.507, 0.470), (1.58, "T1  hesap (ΔSCF, yaklaşık)"),
  "Sınır: xTB geometrisi çok bükük çıktı; sonuç şüpheli, geometri düzeltilecek."),
]
for (n, cas, nev, exp, expn, f, (w1, w2), t1, warn) in M:
    q = w2 > 0.2
    fig = plt.figure(figsize=(13, 6.3), facecolor=BG)
    a = fig.add_axes([0.05, 0.1, 0.42, 0.76]); a.set_facecolor(BG)
    lv = [("S0  taban hal", 0.0, INK, 0.1), ("S1  hesap (VQE+NEVPT2)", nev, BLUE, 0.1)]
    if exp: lv.append(("S1  deney", exp, ORANGE, 0.78))
    if t1: lv.append((t1[1], t1[0], MUTED, 0.78))
    close = exp and t1 and abs(exp - t1[0]) < 0.3
    for lab, e, c, x in lv:
        a.plot([x, x + 0.35], [e, e], color=c, lw=3, solid_capstyle="round")
        dy = (0.12 if lab.startswith("S1  deney") else -0.12 if lab.startswith("T1") else 0) if close else 0
        a.text(x + 0.37, e + dy, f"{lab}\n{e:.2f} eV" + (f"  ({1239.84/e:.0f} nm)" if e else ""), va="center", fontsize=9.5, color=INK)
    a.annotate("", xy=(0.2, nev - 0.03), xytext=(0.2, 0.03), arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2))
    a.text(0.12, nev / 2, "ışık\nsoğurma", ha="right", va="center", fontsize=9.5, color=BLUE)
    if t1 and exp and t1[0] < exp:
        a.annotate("", xy=(1.05, t1[0] + 0.03), xytext=(1.05, exp - 0.03), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5, ls="--"))
    a.axhspan(0, 0.98, xmin=0.80, xmax=0.92, color=GRID)
    a.text(1.46, 0.49, "¹O₂ için\ngereken\n0.98 eV", ha="center", va="center", fontsize=8.5, color=MUTED)
    a.set_xlim(0, 1.7); a.set_ylim(-0.15, 2.5); a.set_xticks([]); a.set_ylabel("Enerji (eV)", color=MUTED)
    for s in ("top", "right", "bottom"): a.spines[s].set_visible(False)
    a.spines["left"].set_color(GRID); a.tick_params(colors=MUTED)
    a.set_title("2 · Kuantum hesabın söylediği: enerji seviyeleri", loc="left", color=INK, fontsize=12.5)
    b = fig.add_axes([0.53, 0.1, 0.44, 0.76]); b.axis("off"); b.set_xlim(0, 1); b.set_ylim(0, 1)
    b.set_title("3 · Yorum: kuantum gerekli mi?", loc="left", color=INK, fontsize=12.5)
    b.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0,rounding_size=0.03", fc="white", ec=GRID, lw=1.2))
    b.text(0.05, 0.9, "Işık soğurunca elektron hangi yoldan gidiyor?", fontsize=11, color=INK)
    b.barh(0.79, w1 * 0.9, left=0.05, height=0.07, color=BLUE)
    b.barh(0.79, w2 * 0.9, left=0.05 + w1 * 0.9 + 0.006, height=0.07, color=ORANGE)
    b.text(0.06, 0.71, f"1. yol  %{w1*100:.0f}", fontsize=9.5, color=INK)
    b.text(0.05 + w1 * 0.9 + 0.01, 0.71, f"2. yol  %{w2*100:.0f}", fontsize=9.5, color=INK)
    b.text(0.05, 0.60, "İki yol neredeyse eşit → tek bir resimle anlatılamaz" if q else "Tek yol baskın → klasik yöntemler yeterli", fontsize=10, color=MUTED)
    b.text(0.05, 0.47, "KUANTUM GEREKLİ" if q else "KLASİK YETERLİ", fontsize=20, color=INK if q else MUTED, fontweight="bold")
    rows = [("Kuantum bölgesi", "4 orbital = 8 kübit (VQE, kesin çözüm)"),
            ("Sadece kuantum bölgesi", f"{cas:.2f} eV"),
            ("+ klasik düzeltme", f"{nev:.2f} eV (kalan ~400 orbital, NEVPT2)"),
            ("Deney", f"{exp:.2f} eV  ({expn}),  fark {nev-exp:+.2f} eV" if exp else expn),
            ("Soğurma şiddeti f", f"{f:.3f}" + ("  (güçlü)" if f > 0.3 else "  (zayıf)"))]
    for i, (k, v) in enumerate(rows):
        b.text(0.05, 0.36 - i * 0.065, k, fontsize=9.3, color=MUTED); b.text(0.38, 0.36 - i * 0.065, v, fontsize=9.3, color=INK)
    if warn: fig.text(0.05, 0.025, "⚠ " + warn, fontsize=9.5, color=WARN)
    fig.suptitle(n, x=0.05, ha="left", fontsize=16, color=INK, fontweight="bold")
    plt.savefig(f"/mnt/user-data/outputs/paneller/{n}_panel_2_3.png", dpi=150, facecolor=BG); plt.close()
print("ok")
