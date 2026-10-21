import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from pyscf.tools import molden
mol, e, C, occ = molden.load("pdt/porfin_orbitaller.molden")[:4]
e = np.array(e) * 27.2114; nocc = int((np.array(occ) > 0).sum())
BLUE, ORANGE, INK, MUTED, GRID, BG = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a64", "#c9c8c2", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 6.5), gridspec_kw={"width_ratios": [1, 1.25]}, facecolor=BG)
for ax in (a, b): ax.set_facecolor(BG)
# sol: tum orbitaller (sadece -20..+10 eV araligi gorunur)
act = list(range(nocc - 2, nocc + 2))
for i, ei in enumerate(e):
    if i in act: continue
    a.plot([0.15, 0.85], [ei, ei], color=GRID, lw=0.8)
a.axhspan(e[nocc - 1] + 0.01, e[nocc] - 0.01, color="#f3f2ee")
a.text(0.5, (e[nocc - 1] + e[nocc]) / 2, "dolu ↓   boş ↑", ha="center", va="center", fontsize=9.5, color=MUTED)
for i in act: a.plot([0.05, 0.95], [e[i], e[i]], color=BLUE if i < nocc else ORANGE, lw=3)
a.set_ylim(-20, 8); a.set_xlim(0, 1); a.set_xticks([]); a.set_ylabel("Orbital enerjisi (eV)", color=MUTED)
for s in ("top", "right", "bottom"): a.spines[s].set_visible(False)
a.spines["left"].set_color(GRID); a.tick_params(colors=MUTED)
a.set_title(f"Porfin: {len(e)} orbitalin tamamı", loc="left", color=INK, fontsize=12.5)
a.text(0.0, -21.8, f"{nocc} dolu + {len(e) - nocc} boş orbital  (görünen aralık: −20…+8 eV)", fontsize=9, color=MUTED)
# sag: yakinlastirma - secilen 4 orbital
names = ["HOMO-1", "HOMO", "LUMO", "LUMO+1"]
for k, i in enumerate(act):
    c = BLUE if i < nocc else ORANGE
    b.plot([0.1, 0.45], [e[i], e[i]], color=c, lw=4, solid_capstyle="round")
    b.text(0.48, e[i] + (0.09 if k in (1, 3) else -0.09), f"{names[k]}   {e[i]:.2f} eV", va="center", fontsize=10.5, color=INK)
for i in range(nocc - 6, nocc + 6):
    if i not in act and -8 < e[i] < 0: b.plot([0.1, 0.45], [e[i], e[i]], color=GRID, lw=1)
b.annotate("", xy=(0.28, e[nocc] - 0.05), xytext=(0.28, e[nocc - 1] + 0.05), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.5))
b.text(0.30, (e[nocc] + e[nocc - 1]) / 2, "ışık", fontsize=10, color=INK, va="center")
b.text(0.48, -3.95, "İki dolu ve iki boş orbital\nneredeyse eş enerjili\n→ elektron iki yolu birden kullanıyor\n→ KUANTUM BÖLGESİ: 4 orbital = 8 kübit",
       fontsize=10.5, color=INK, va="center", linespacing=1.5)
b.text(0.48, -7.0, f"Geri kalan {len(e) - 4} orbital → klasik hesap (NEVPT2)", fontsize=10.5, color=MUTED, va="center")
b.set_ylim(-7.8, -1.4); b.set_xlim(0, 1.25); b.set_xticks([]); b.set_yticks([])
for s in ("top", "right", "bottom", "left"): b.spines[s].set_visible(False)
b.set_title("Yakınlaştırma: seçilen kuantum bölgesi", loc="left", color=INK, fontsize=12.5)
plt.tight_layout(); plt.savefig("/mnt/user-data/outputs/porfin_orbital_secimi.png", dpi=150, facecolor=BG)
print(nocc, len(e), e[act])
