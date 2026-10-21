import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
BLUE, ORANGE, INK, MUTED, GRID, BG, WARN = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a64", "#e4e3df", "#fcfcfb", "#9a5b00"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
D = [  # ad, klinik, cekirdek, CAS, NEV, f, w1, w2, klinik nm, klinik not, arti metni, uyari
 ("Temoporfin (Foscan)", "baş-boyun kanseri, AB", "klorin", 3.45, 1.98, 0.579, 0.743, 0.098, 652, "klinik ışık 652 nm (EMA)",
  ["Tek yol baskın: klasik yöntem (TDDFT vb.) yeterli",
   "Kazanç: kuantum kaynağı harcamadan doğru tahmin",
   "Güçlü kırmızı soğurma (f yüksek) → klinikte tercih sebebi"], ""),
 ("Talaporfin (Laserphyrin)", "akciğer kanseri, Japonya", "klorin", 3.45, 1.98, 0.579, 0.743, 0.098, 664, "klinik lazer 664 nm (doğrulanacak)",
  ["Tek yol baskın: klasik yöntem yeterli",
   "Kazanç: hızlı ve ucuz tarama",
   "Klorin ailesi: güçlü kırmızı soğurma"], ""),
 ("Protoporfirin IX (5-ALA)", "floresan cerrahi ve PDT", "porfin", 3.62, 1.87, 0.033, 0.544, 0.435, 635, "klinik ışık ~635 nm (doğrulanacak)",
  ["İki yol neredeyse eşit: tek resimli klasik tahmin 2.94 eV (≈1 eV hatalı)",
   "Hibrit (8 kübit VQE + NEVPT2): 1.87 eV → klinik değere ~0.08 eV",
   "Zayıf Q bandı (f düşük) → gerçekte de ince ayar gereken ajan"],
  "Not: kenar grupları simetriyi bozar; gerçek PpIX'te karışım çekirdek modelden daha zayıf olabilir."),
]
for (ad, klin, core, cas, nev, f, w1, w2, nm, knot, arti, warn) in D:
    ek = 1239.84 / nm; q = w2 > 0.2
    fig = plt.figure(figsize=(13.5, 6.6), facecolor=BG)
    # Panel 2: soğurma
    a = fig.add_axes([0.20, 0.12, 0.29, 0.72]); a.set_facecolor(BG)
    a.barh([1], [nev], color=BLUE, height=0.45); a.barh([0], [ek], color=ORANGE, height=0.45)
    a.text(nev + 0.03, 1, f"{nev:.2f} eV  ({1239.84/nev:.0f} nm)", va="center", fontsize=10.5, color=INK)
    a.text(ek + 0.03, 0, f"{ek:.2f} eV  ({nm} nm)", va="center", fontsize=10.5, color=INK)
    a.set_yticks([0, 1]); a.set_yticklabels([f"Klinik / literatür\n({knot})", f"Bizim hesap\n({core} çekirdeği, VQE+NEVPT2)"], fontsize=9.5, color=INK)
    a.set_xlim(0, 3.2); a.set_xlabel("Soğurulan ışığın enerjisi (eV)", color=MUTED)
    a.text(-1.55, -0.85, f"Fark: {nev-ek:+.2f} eV      Soğurma şiddeti f = {f:.3f}  " + ("(güçlü)" if f > 0.3 else "(zayıf)"), fontsize=10.5, color=INK)
    a.text(-1.55, -1.15, "Doku penceresi: kırmızı/yakın kızılötesi (~600–850 nm) ışık dokuya daha derin işler", fontsize=9, color=MUTED)
    a.set_ylim(-1.35, 1.6)
    for s in ("top", "right", "left"): a.spines[s].set_visible(False)
    a.spines["bottom"].set_color(GRID); a.tick_params(colors=MUTED, length=0); a.grid(axis="x", color=GRID, lw=0.8); a.set_axisbelow(True)
    a.set_title("2 · Işık soğurma: hesap ve klinik", loc="left", color=INK, fontsize=12.5)
    # Panel 3: yorum
    b = fig.add_axes([0.55, 0.12, 0.42, 0.72]); b.axis("off"); b.set_xlim(0, 1); b.set_ylim(0, 1)
    b.set_title("3 · Kuantum gerekli mi? Ne kazandırır?", loc="left", color=INK, fontsize=12.5)
    b.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0,rounding_size=0.03", fc="white", ec=GRID, lw=1.2))
    b.text(0.05, 0.9, "Işık soğurunca elektron hangi yoldan gidiyor?", fontsize=11, color=INK)
    b.barh(0.79, w1 * 0.9, left=0.05, height=0.07, color=BLUE)
    b.barh(0.79, w2 * 0.9, left=0.05 + w1 * 0.9 + 0.006, height=0.07, color=ORANGE)
    b.text(0.06, 0.71, f"1. yol  %{w1*100:.0f}", fontsize=9.5, color=INK)
    b.text(0.05 + w1 * 0.9 + 0.01, 0.71, f"2. yol  %{w2*100:.0f}", fontsize=9.5, color=INK)
    b.text(0.05, 0.56, "KUANTUM GEREKLİ" if q else "KLASİK YETERLİ", fontsize=20, color=INK if q else MUTED, fontweight="bold")
    for i, t in enumerate(arti): b.text(0.05, 0.43 - i * 0.085, "• " + t, fontsize=9.8, color=INK)
    b.text(0.05, 0.12, "Kuantum bölgesi: 4 orbital = 8 kübit  |  kalan ~400 orbital: klasik (NEVPT2)", fontsize=9, color=MUTED)
    if warn: fig.text(0.05, 0.025, warn, fontsize=9.3, color=WARN)
    fig.suptitle(f"{ad}  —  {klin}", x=0.05, ha="left", fontsize=15.5, color=INK, fontweight="bold")
    fn = ad.split(" (")[0].replace(" ", "_")
    plt.savefig(f"/mnt/user-data/outputs/paneller/ILAC_{fn}_panel.png", dpi=150, facecolor=BG); plt.close()
print("ok")
