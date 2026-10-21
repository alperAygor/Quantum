import os  # CIKTI_DIR: figurlerin yazilacagi klasor (varsayilan: calisilan klasor)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
BLUE, ORANGE, INK, MUTED, GRID, BG = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a64", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
fig = plt.figure(figsize=(13, 6), facecolor=BG)
# --- Panel 2: enerji diyagramı
a = fig.add_axes([0.05, 0.1, 0.42, 0.78]); a.set_facecolor(BG)
lv = [("S0  taban hal", 0.0, INK, 0.1), ("S1  hesap (VQE+NEVPT2)", 1.87, BLUE, 0.1),
      ("S1  deney", 1.98, ORANGE, 0.78), ("T1  deney (fosforesans)", 1.58, MUTED, 0.78)]
for lab, e, c, x in lv:
    a.plot([x, x + 0.35], [e, e], color=c, lw=3, solid_capstyle="round")
    a.text(x + 0.37, e, f"{lab}\n{e:.2f} eV" + (f"  ({1239.84/e:.0f} nm)" if e else ""), va="center", fontsize=9.5, color=INK)
a.annotate("", xy=(0.2, 1.85), xytext=(0.2, 0.03), arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2))
a.text(0.12, 0.95, "ışık\nsoğurma", ha="right", va="center", fontsize=9.5, color=BLUE)
a.annotate("", xy=(1.05, 1.6), xytext=(1.05, 1.96), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5, ls="--"))
a.text(0.80, 1.78, "", fontsize=8)
a.axhspan(0, 0.98, xmin=0.80, xmax=0.92, color=GRID)
a.text(1.46, 0.49, "¹O₂ için\ngereken\n0.98 eV", ha="center", va="center", fontsize=8.5, color=MUTED)
a.set_xlim(0, 1.7); a.set_ylim(-0.15, 2.35); a.set_xticks([])
a.set_ylabel("Enerji (eV)", color=MUTED)
for s in ("top", "right", "bottom"): a.spines[s].set_visible(False)
a.spines["left"].set_color(GRID); a.tick_params(colors=MUTED)
a.set_title("2 · Kuantum hesabın söylediği: enerji seviyeleri", loc="left", color=INK, fontsize=12.5)
# --- Panel 3: triyaj kartı
b = fig.add_axes([0.53, 0.1, 0.44, 0.78]); b.axis("off")
b.set_title("3 · Yorum: kuantum gerekli mi?", loc="left", color=INK, fontsize=12.5)
b.add_patch(FancyBboxPatch((0.0, 0.0), 1.0, 1.0, boxstyle="round,pad=0.0,rounding_size=0.03", fc="white", ec=GRID, lw=1.2))
b.text(0.05, 0.9, "Işık soğurunca elektron hangi yoldan gidiyor?", fontsize=11, color=INK)
b.barh(0.78, 0.54*0.9, left=0.05, height=0.07, color=BLUE)
b.barh(0.78, 0.44*0.9, left=0.05 + 0.54*0.9 + 0.005, height=0.07, color=ORANGE)
b.text(0.06, 0.70, "HOMO-1 → LUMO+1   %54", fontsize=9.5, color=INK)
b.text(0.55, 0.70, "HOMO → LUMO   %44", fontsize=9.5, color=INK)
b.text(0.05, 0.58, "İki yol neredeyse eşit → tek bir resimle anlatılamaz", fontsize=10, color=MUTED)
b.text(0.05, 0.44, "KUANTUM GEREKLİ", fontsize=20, color=INK, fontweight="bold")
rows = [("Kuantum bölgesi", "4 orbital = 8 kübit (VQE, kesin çözüm)"),
        ("Klasik kısım", "kalan ~400 orbital (NEVPT2)"),
        ("Deneyle fark", "−0.11 eV (1.87 / 1.98 eV)"),
        ("Klasik tek-resim tahmini", "2.94 eV (0.96 eV hatalı)")]
for i, (k, v) in enumerate(rows):
    b.text(0.05, 0.32 - i*0.075, k, fontsize=9.5, color=MUTED); b.text(0.40, 0.32 - i*0.075, v, fontsize=9.5, color=INK)
b.set_xlim(0, 1); b.set_ylim(0, 1)
fig.suptitle("Porfin", x=0.05, ha="left", fontsize=16, color=INK, fontweight="bold")
plt.savefig(os.environ.get("CIKTI_DIR", ".") + "/porfin_panel_2_3.png", dpi=150, facecolor=BG)
