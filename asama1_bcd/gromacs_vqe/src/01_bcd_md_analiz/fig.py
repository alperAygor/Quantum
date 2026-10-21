import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = [  # ad, hesap S1+NEVPT2, deney (None), f, triyaj, not
 ("Porfin",            1.87, 1.98, 0.033, "Kuantum gerekli", ""),
 ("Klorin",            1.98, 1.95, 0.579, "Klasik yeterli", ""),
 ("Zn-porfin",         2.07, 2.18, 0.079, "Kuantum gerekli", ""),
 ("Mg-porfin",         2.04, 2.17, 0.073, "Kuantum gerekli", ""),
 ("Pd-porfin",         2.14, None, 0.085, "Kuantum gerekli", "deney yok"),
 ("Porfirazin",        2.24, 2.07, 0.633, "Klasik yeterli", ""),
 ("Tetrametilporfin*", 1.64, 1.88, 0.001, "Kuantum gerekli", "geometri şüpheli"),
 ("Bakteriyoklorin*",  2.23, 1.74, 0.422, "Klasik yeterli", "4 orbital yetersiz"),
]
D = D[::-1]
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a64", "#e4e3df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5})
fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 5.6), gridspec_kw={"width_ratios": [1.35, 1]}, facecolor="#fcfcfb")
y = range(len(D))
for i, (n, c, e, f, t, nt) in enumerate(D):
    if e is not None: a.plot([c, e], [i, i], color=GRID, lw=2.5, zorder=1, solid_capstyle="round")
    a.scatter(c, i, s=70, color=BLUE, zorder=3, edgecolor="#fcfcfb", linewidth=2)
    if e is not None: a.scatter(e, i, s=70, color=ORANGE, zorder=3, edgecolor="#fcfcfb", linewidth=2)
    if nt: a.text(2.42, i, nt, va="center", fontsize=8.5, color=MUTED)
a.set_yticks(list(y)); a.set_yticklabels([d[0] for d in D], color=INK)
a.set_xlim(1.55, 2.8); a.set_xlabel("S1 enerjisi (eV)", color=MUTED)
a.set_title("Q bandı: hesap (8 kübit VQE + NEVPT2) ve deney", loc="left", color=INK, fontsize=12)
a.scatter([], [], color=BLUE, s=50, label="Hesap"); a.scatter([], [], color=ORANGE, s=50, label="Deney")
a.legend(loc="upper right", frameon=False, bbox_to_anchor=(1.0, 1.0), fontsize=9.5)
a.text(1.55, -1.25, "Ortalama mutlak hata: 0.19 eV (7 molekül), bakteriyoklorin hariç 0.14 eV.   * yöntemin sınırı",
       fontsize=8.5, color=MUTED)
for i, (n, c, e, f, t, nt) in enumerate(D):
    b.barh(i, f, height=0.55, color=BLUE)
    b.text(f + 0.015, i, f"{f:.3f}", va="center", fontsize=9, color=INK)
    b.text(0.98, i, t, va="center", ha="right", fontsize=9, color=INK if t.startswith("Kuantum") else MUTED,
           fontweight="bold" if t.startswith("Kuantum") else "normal", transform=b.get_yaxis_transform())
b.set_yticks(list(y)); b.set_yticklabels([]); b.set_xlim(0, 1.0)
b.set_xlabel("Soğurma şiddeti f(S1)", color=MUTED)
b.set_title("Soğurma şiddeti ve triyaj kararı", loc="left", color=INK, fontsize=12)
for ax in (a, b):
    ax.set_facecolor("#fcfcfb"); ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID); ax.tick_params(colors=MUTED, length=0)
plt.tight_layout(); plt.savefig("/mnt/user-data/outputs/tarama_sonuclari.png", dpi=150, facecolor="#fcfcfb")
