"""Notebook'ları üretir:  python tools/make_notebooks.py"""
from pathlib import Path

import nbformat as nbf

import os
OUT = Path(os.environ.get("NB_OUT", Path(__file__).resolve().parent.parent / "notebooks"))
OUT.mkdir(exist_ok=True)


def md(s):
    return nbf.v4.new_markdown_cell(s.strip("\n"))


def code(s):
    return nbf.v4.new_code_cell(s.strip("\n"))


SETUP = '''
import sys, os
sys.path.insert(0, os.path.abspath(".."))
from qfe import config, store
from qfe.plotting import C, LEG, LEG_TR
import numpy as np, matplotlib.pyplot as plt
# Ön ayar: "test" (hızlı kontrol), "demo" (ekibe gösterim), "standart" (daha güvenilir)
config.set_preset(os.environ.get("QFE_PRESET", "demo"))
print("Aktif ön ayar:", config.preset_name())
'''


def save(name, cells):
    nb = nbf.v4.new_notebook()
    nb["cells"] = cells
    nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nbf.write(nb, OUT / name)
    print("yazıldı:", OUT / name)


# =============================================================================
# 00
# =============================================================================
save("00_genel_bakis_ve_kurulum.ipynb", [
    md(r"""
# 00 · Genel bakış: FreeQuantum pipeline'ını laptopta simüle etmek

**Referans makale:** Günther et al., *How to Use Quantum Computers for Biomolecular Free Energies*,
J. Chem. Theory Comput. 2026, 22, 4329–4345 (FreeQuantum pipeline'ı).

## Makale ne yapıyor?
Bir ilacın proteine **bağlanma serbest enerjisini** ($\Delta G_\text{bağlanma}$) hesaplamak istiyoruz.
Tam kuantum mekaniksel hesap binlerce atom için imkânsız. Bu yüzden katmanlı bir model kullanılıyor:

1. **MM:** Klasik kuvvet alanı ile MD ve alkimik serbest enerji (FEP).
2. **QM/MM:** Ligand (kuantum bölgesi) DFT ile, gerisi MM ile hesaplanıyor. Bu veriyle bir ML potansiyeli (ML1) eğitiliyor.
3. **QM/QM/MM:** Ligandın içindeki küçük bir *kuantum çekirdeği* çok yüksek doğrulukla (CASCI/NEVPT2, CCSD(T)) hesaplanıyor.
   ML1, transfer öğrenme ile ML2'ye dönüştürülüyor.
4. **Kuantum bilgisayar:** Gelecekte 3. adımdaki kuantum çekirdeği enerjileri bir kuantum bilgisayarda
   (kuantum faz tahmini, QPE) hesaplanacak. Makale bunun için kaynak tahmini veriyor.

## Bu demoda ne yapıyoruz?
Aynı mantığı **küçültülmüş, laptopta çalışan** bir sistemde baştan sona simüle ediyoruz. Kullandığımız sistem
**β-siklodekstrin (konakçı) + siklopentanol (misafir)**. Bu ilaç benzeri bir "bağlanma" problemi ve
gerçek deneysel bağlanma serbest enerjisi biliniyor:

> ΔG_deney = −3.050 ± 0.010 kcal/mol = **−12.76 kJ/mol** (Rekharsky et al., *J. Phys. Chem. B* 1997, 101, 87).

Girdi dosyaları da gerçek, yayımlanmış bir benchmark setinden geliyor: MobleyLab *benchmarksets*, CD Set 2.

| Makale (FreeQuantum) | Bu demo | Bu demoda nerede? | Gerçek projede nerede? |
|---|---|---|---|
| GRP78 proteini + NKP1339 (Ru ilacı) | β-siklodekstrin + siklopentanol | — | — |
| MM-FEP (Amber ff99SB-ILDN), MBAR | MM-FEP (q4md-CD + GAFF), MBAR | laptop GPU | GPU kümesi |
| QM/MM, PBE-D3/def2-SVP, ~2000 nokta | QM/MM, PBE/def2-SVP, onlarca nokta | laptop CPU | **HPC (CPU)** |
| ML1 + aktif öğrenme + NEQ anahtarlama | Gauss süreci + aktif öğrenme + Zwanzig | laptop | GPU (ML eğitimi) |
| Kuantum çekirdek: DMRG-CASCI + NEVPT2 | CASCI(4,4) + NEVPT2 | laptop | **HPC** → gelecekte **kuantum bilgisayar** |
| ML2 (transfer öğrenme) | GP düzeltme modeli | laptop | GPU |
| QPE kaynak tahmini (Şekil 4) | QPE/VQE/Trotter simülasyonu + ölçekleme | laptop (≤ ~24 qubit) | **Hata düzeltmeli kuantum bilgisayar** |

## Notebook sırası
| # | Notebook | Makaledeki karşılığı |
|---|---|---|
| 01 | Veri ve sistem hazırlama | Şekil 2, gri kutu |
| 02 | MM serbest enerjisi (FEP + MBAR) | Bölüm 2.1, Denklem 3–5 |
| 03 | QM bölgesi ve QM/MM enerjileri | Bölüm 2.3, Denklem 6 |
| 04 | ML1 + aktif öğrenme → ΔG(MM+ML1) | Bölüm 2.2, Şekil 6–7 |
| 05 | Kuantum çekirdeği + ML2 → ΔG(MM+ML2) | Bölüm 2.3–2.4, Denklem 7 |
| 06 | Kuantum algoritmaları (VQE, QPE, Trotter) | Bölüm 3, Şekil 4d |
| 07 | Kaynak tahmini, hata yayılımı, yol haritası | Bölüm 3, Şekil 4a–c |

> **Dürüst not:** Bu bir *eğitim simülasyonu*. Amaç, her adımın ne yaptığını ve maliyetini görmek.
> Kısa örnekleme nedeniyle sayılar istatistiksel olarak gürültülü olacak. Hangi adımın büyütülünce
> HPC veya kuantum bilgisayar gerektirdiğini 07. notebook'ta ölçülen sürelerden çıkaracağız.
"""),
    md("## Ortam kontrolü"),
    code(SETUP),
    code('''
import importlib, platform
print("Python", platform.python_version(), "|", platform.platform())
for m in ["openmm", "openmmtools", "pymbar", "parmed", "pyscf", "openfermion", "sklearn", "scipy", "numpy"]:
    try:
        mod = importlib.import_module(m)
        print(f"  ✓ {m:12s} {getattr(mod, '__version__', '')}")
    except Exception as e:
        print(f"  ✗ {m:12s} EKSİK -> {e}")
'''),
    code('''
import openmm as mm
plats = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
print("OpenMM platformları:", plats)
if "CUDA" in plats:
    print("NVIDIA GPU (CUDA) kullanılacak ✓")
elif "OpenCL" in plats:
    print("OpenCL kullanılacak (GPU olabilir)")
else:
    print("Sadece CPU var. 'test' ön ayarıyla çalışın; 'demo' çok uzun sürebilir.")
'''),
    md("""
### Termodinamik döngü (çift ayrıştırma / double decoupling)

```
  misafir(su) + konakçı(su)  ──ΔG_bağlanma──►  kompleks(su)
        │                                          │
  ΔG_s  │ misafiri "kapat"                  ΔG_c   │ misafiri "kapat" (kısıtlı)
        ▼                                          ▼
  misafir(gaz benzeri) + konakçı(su) ◄── kT·ln(V0/V_site) ── kapalı misafir konakçıda
```
ΔG_bağlanma = ΔG_s − ΔG_c + kT ln(V0/V_site). Bu, makaledeki Denklem 5'in aynısı.
"""),
])

# =============================================================================
# 01
# =============================================================================
save("01_veri_ve_sistem_hazirlama.ipynb", [
    md("""
# 01 · Gerçek veri ve MM sistemlerinin hazırlanması
*Makale: Şekil 2 "Structure Preparation and Equilibration (MM)".*

1. Yayımlanmış benchmark girdilerini indiriyoruz: solvatize edilmiş bCD–siklopentanol kompleksi
   (q4md-CD + GAFF kuvvet alanı, TIP3P su, 1500 su molekülü).
2. **Kompleks bacağını** dosyadan olduğu gibi alıyoruz.
3. **Çözücü bacağı** için aynı kutudan konakçıyı silip misafiri saf suda bırakıyoruz.
4. Her iki sistemi NPT (298 K, 1 bar) koşullarında dengeliyoruz.
"""),
    code(SETUP),
    code('''
from qfe import mmsys
with store.timed("01_indirme"):
    files = mmsys.download_inputs()
'''),
    code('''
# Deneysel veri (benchmark README'sinden; Rekharsky et al. 1997)
dG, sem = config.exp_dG_kj()
print(f"Deneysel ΔG_bağlanma (bCD + siklopentanol, 298 K): {dG:.2f} ± {sem:.2f} kJ/mol "
      f"({config.SYSTEM['exp_dG_kcal']} kcal/mol)")
store.save_json("results.json", {"exp_dG_kJ": dG, "exp_dG_err_kJ": sem})
'''),
    code('''
legs = mmsys.build_legs(files)
for leg, (_, _, summ) in legs.items():
    print(f"{leg:8s}: {summ['atoms']} atom, rezidüler={summ['residues']}, net yük={summ['net_charge']}")
'''),
    md("""
### Moleküle bakış (isteğe bağlı)
`py3Dmol` kuruluysa kompleksin 3B görüntüsü aşağıda çıkar (su molekülleri gizli).
"""),
    code('''
try:
    import py3Dmol, parmed as pmd
    s = pmd.load_file(str(config.DATA / "complex.prmtop"), xyz=str(config.DATA / "complex.rst7"))
    sub = s[":MGO,MOL"]
    pdb = config.DATA / "complex_dry.pdb"; sub.save(str(pdb), overwrite=True)
    v = py3Dmol.view(width=520, height=380)
    v.addModel(open(pdb).read(), "pdb")
    v.setStyle({"resn": "MGO"}, {"stick": {"radius": 0.12, "color": "lightgray"}})
    v.setStyle({"resn": "MOL"}, {"stick": {"radius": 0.22, "colorscheme": "orangeCarbon"}, "sphere": {"scale": 0.25}})
    v.zoomTo(); v.show()
except ImportError:
    print("py3Dmol yok (pip install py3Dmol) — görselleştirme atlandı.")
'''),
    md("### NPT dengeleme"),
    code('''
eq = {}
for leg in ["complex", "solvent"]:
    with store.timed(f"01_dengeleme_{leg}") as rec:
        eq[leg] = mmsys.equilibrate(leg)
        rec["ns_per_day"] = eq[leg]["ns_per_day"]
store.save_json("results.json", {"md_ns_per_day": eq["complex"]["ns_per_day"]})

fig, ax = plt.subplots(figsize=(6.5, 3.4))
for leg in eq:
    v = np.array(eq[leg]["volumes_nm3"])
    ax.plot(np.linspace(0, config.P()["equil_npt_ps"], len(v)), v, color=LEG[leg], label=LEG_TR[leg])
ax.set_xlabel("süre (ps)"); ax.set_ylabel("kutu hacmi (nm³)")
ax.set_title("NPT dengeleme: kutu hacmi")
ax.legend(); plt.tight_layout(); plt.show()
'''),
    md("""
**Ne görmeliyiz?** Çözücü bacağında konakçı silindiği için boş bir oyuk kalıyor. Barostat kutuyu
küçülterek bu oyuğu suyla dolduruyor, bu yüzden hacim düşüp sabitlenmeli. Hacim hâlâ düşüyorsa
`equil_npt_ps` değerini artırın.

**Maliyet notu:** Hücrenin çıktısındaki *ns/gün* değeri MD hızınız. 07. notebook bunu kullanarak
gerçek bir protein sistemi için gereken GPU süresini tahmin edecek.
"""),
])

# =============================================================================
# 02
# =============================================================================
save("02_MM_serbest_enerji_FEP.ipynb", [
    md(r"""
# 02 · MM serbest enerjisi: alkimik FEP + MBAR
*Makale: Bölüm 2.1, Denklem 3–5; Şekil 2 kırmızı üst kutu.*

Misafirin çevresiyle etkileşimini bir $\lambda$ parametresiyle adım adım "kapatıyoruz":

$$H_\lambda = H_\text{çevre} + H_\text{misafir} + \lambda\,H_\text{etkileşim}\qquad(\text{makale Denklem 4})$$

Önce elektrostatik etkileşimi ($\lambda_e$: 1→0), sonra van der Waals etkileşimini ($\lambda_s$: 1→0)
kapatıyoruz. vdW kapatılırken, atomlar üst üste binince enerjinin sonsuza gitmemesi için "soft-core"
potansiyel kullanılıyor. Her $\lambda$ durumunda MD yapıyoruz ve **her örnekte tüm durumların enerjisini**
kaydediyoruz. Sonra makaledeki gibi **MBAR** ile durumlar arası serbest enerji farklarını hesaplıyoruz.
"""),
    code(SETUP),
    code('''
from qfe import fep, mmsys
states = fep.lambda_schedule()
print(f"{len(states)} λ durumu (λ_elektrostatik, λ_sterik):")
print(states)
p = config.P()
print(f"Her durumda {p['fep_equil_ps']} ps dengeleme + {p['fep_prod_ps']} ps üretim, "
      f"{p['fep_sample_ps']} ps'de bir örnek")
'''),
    code('''
raw = {}
for leg in ["complex", "solvent"]:
    with store.timed(f"02_fep_{leg}") as rec:
        raw[leg] = fep.run_leg(leg)
        rec["simulated_ps"] = len(states) * (p["fep_equil_ps"] + p["fep_prod_ps"])
    store.save_npz(f"fep_{leg}.npz", **raw[leg])
'''),
    code('''
an = {leg: fep.mbar_analysis(raw[leg]["u_kln"]) for leg in raw}
res = fep.combine(an["complex"], an["solvent"])
expG, expE = config.exp_dG_kj()
print(f"ΔG_c (kompleksde kapatma) = {res['dG_complex_decouple']:7.2f} kJ/mol")
print(f"ΔG_s (suda kapatma)       = {res['dG_solvent_decouple']:7.2f} kJ/mol")
print(f"Standart durum düzeltmesi = {res['standard_state_corr']:7.2f} kJ/mol  (V_site={res['v_site_nm3']:.3f} nm³)")
print(f"------------------------------------------")
print(f"ΔG_bağlanma (MM)          = {res['dG_bind']:7.2f} ± {res['dG_bind_err']:.2f} kJ/mol")
print(f"Deney                     = {expG:7.2f} ± {expE:.2f} kJ/mol")
store.save_json("results.json", {"dG_MM": res})
'''),
    code('''
fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
for leg in ["complex", "solvent"]:
    f = an[leg]["f_kJ"]; df = an[leg]["df_kJ"]
    axs[0].errorbar(range(len(f)), f, yerr=df, color=LEG[leg], marker="o", ms=4, label=LEG_TR[leg])
axs[0].set_xlabel("λ durumu (0 = tam etkileşim)"); axs[0].set_ylabel("G(λ) − G(0)  [kJ/mol]")
axs[0].set_title("Serbest enerji profili"); axs[0].legend(fontsize=8)
for ax, leg in zip(axs[1:], ["complex", "solvent"]):
    O = an[leg]["overlap"]
    if O is not None:
        im = ax.imshow(O, cmap="Blues", vmin=0, vmax=max(0.3, O.max()))
        ax.set_title(f"MBAR örtüşme matrisi — {leg}"); ax.set_xlabel("durum"); ax.set_ylabel("durum")
        ax.grid(False); plt.colorbar(im, ax=ax, fraction=0.046)
plt.tight_layout(); plt.show()
'''),
    md("""
**Nasıl okunur?**
- **Örtüşme matrisi:** Komşu λ durumları arasındaki örtüşme (köşegenin hemen yanındaki değerler)
  ~0.03'ten büyük olmalı. Değilse o aralığa yeni bir λ durumu eklenmeli. Makale de adımları bu ölçüte göre seçiyor.
- **MM sonucu ile deney arasındaki fark:** Bu fark kuvvet alanının hatası + örnekleme hatasından oluşuyor.
  Makaledeki fikir, kuvvet alanı hatasını kuantum mekaniğiyle (QM) düzeltmek.

### Kuantum hesapları için anlık görüntü havuzu
Tam etkileşimli (λ=1) durumdan yapılar topluyoruz. 03–05. notebook'larda QM/MM ve ML adımları
bu yapılar üzerinde çalışacak (makalede ~2000 anlık görüntü).
"""),
    code('''
for leg in ["complex", "solvent"]:
    with store.timed(f"02_havuz_{leg}") as rec:
        X, B = fep.sample_endstate_pool(leg)
        rec["frames"] = len(X)
    store.save_npz(f"pool_{leg}.npz", X=X, B=B)
    print(leg, X.shape)
'''),
])

# =============================================================================
# 03
# =============================================================================
save("03_QM_bolgesi_ve_QMMM.ipynb", [
    md(r"""
# 03 · Kuantum bölgesi ve QM/MM enerjileri
*Makale: Bölüm 2.3, Denklem 6; Şekil 2 turuncu/turkuaz kutular.*

**Kuantum bölgesi = misafir molekül** (makalede de QM bölgesi ligand). Çevredeki her atom (konakçı, su,
iyonlar) QM hesabına bir **nokta yük** olarak giriyor. Buna elektrostatik gömme deniyor:

$$H_\text{QM bölgesi} = H_\text{elektronik, QM} + V_\text{etkileşim, MM}\qquad(\text{makale Denklem 6})$$

Her anlık görüntü $x$ için QM ile MM arasındaki enerji farkını iki parçaya ayırıyoruz:

$$\Delta U_\text{etk}(x) = \big[E_\text{QM}(\text{misafir}\,|\,\text{yükler}) - E_\text{QM}(\text{misafir, vakum})\big] - E^\text{MM}_\text{Coulomb}(\text{misafir–çevre})$$
$$\Delta U_\text{iç}(x) = E_\text{QM}(\text{misafir, vakum}) - E^\text{MM}_\text{iç}(\text{misafir})$$

- $\Delta U_\text{etk}$: çevreyle **etkileşimin** kuantum düzeltmesi (polarizasyon, yük dağılımının doğru tarifi).
  Bağlanmayı doğrudan etkileyen kısım bu.
- $\Delta U_\text{iç}$: misafirin **kendi iç enerjisinin** düzeltmesi (DFT ile GAFF'ın bağ ve açı terimleri farklı).

Misafir–çevre van der Waals terimi her iki seviyede de MM'den geliyor ve sadeleşiyor.
QM seviyesi olarak **PBE/def2-SVP** kullanıyoruz (makale: PBE-D3(BJ)/def2-SVP).
"""),
    code(SETUP),
    code('''
from qfe import qmmm
from pyscf import lib
lib.num_threads(os.cpu_count())
pools = {leg: (lambda d: (d["X"], d["B"]))(store.load_npz(f"pool_{leg}.npz")) for leg in ["complex", "solvent"]}
regions = {leg: qmmm.QMRegion(leg) for leg in pools}
for leg, r in regions.items():
    fr = r.frame(pools[leg][0][0], pools[leg][1][0])
    print(f"{leg:8s}: QM atomu={len(fr['qm_xyz'])}, nokta yük={len(fr['pc_q'])}, "
          f"nokta yüklerin toplamı={fr['pc_q'].sum():+.3f} e  (kesme {config.P()['qm_cutoff_A']} Å)")
print("Misafir:", "".join(regions["complex"].guest_symbols()))
'''),
    md("""
### Tek bir QM/MM hesabı: neler hesaplanıyor ve maliyeti nedir?
Aşağıdaki hücre, bir anlık görüntü için gömülü ve vakum DFT hesaplarının süresini ölçüyor. Bu süre,
07. notebook'taki HPC tahmininin temel girdisi.
"""),
    code('''
leg, i = "complex", 0
fr = regions[leg].frame(pools[leg][0][i], pools[leg][1][i])
sym = regions[leg].guest_symbols()
with store.timed("03_tek_DFT"):
    e_emb = qmmm.energy_dft(fr, sym)
e_vac = qmmm.energy_dft(fr, sym, embed=False)
mol = qmmm._mol(fr, sym, config.P()["qm_basis"])
e_mm, e_int, e_coul = regions[leg].mm_qm_part(fr)
H = config.HARTREE_TO_KJ
print(f"E_DFT(gömülü) = {e_emb:.6f} Ha, E_DFT(vakum) = {e_vac:.6f} Ha  ({mol.nao} baz fonksiyonu)")
print(f"QM etkileşim enerjisi (gömülü − vakum) : {(e_emb - e_vac)*H:8.2f} kJ/mol")
print(f"MM Coulomb etkileşimi (sabit yükler)   : {e_coul:8.2f} kJ/mol")
print(f"→ ΔU_etk = {(e_emb - e_vac)*H - e_coul:6.2f} kJ/mol  (çoğunlukla polarizasyon: sabit MM yükleri polarize olamaz)")
store.save_json("results.json", {"qm_nbasis": int(mol.nao), "qm_natoms": len(fr["qm_xyz"]),
                                 "n_threads": os.cpu_count()})
'''),
    md("""
### Birkaç anlık görüntüde iki parçanın dağılımı
Her bacaktan rastgele birkaç yapı hesaplıyoruz. Bu etiketler önbelleğe alınıyor ve 04. notebook'ta
tekrar kullanılıyor, yani boşa gitmiyorlar.
"""),
    code('''
lab = qmmm.Labeler(regions, pools)
lab.load_dict(store.load_json("labels.json", default={}))
n_probe = max(4, config.P()["al_validation"] // 2)
rng = np.random.default_rng(1)
probe = {leg: rng.choice(len(pools[leg][0]), n_probe, replace=False) for leg in pools}
with store.timed("03_DFT_ornekleri") as rec:
    etk = {leg: np.array([lab.dU_dft(leg, i) for i in probe[leg]]) for leg in pools}
    ic = {leg: np.array([lab.dU_intra(leg, i) for i in probe[leg]]) for leg in pools}
    rec["n"] = 2 * n_probe
store.save_json("labels.json", lab.to_dict())
store.save_json("probe.json", {leg: [int(i) for i in probe[leg]] for leg in probe})

fig, axs = plt.subplots(1, 2, figsize=(11, 3.4))
off_ic = np.mean(np.concatenate(list(ic.values())))
for ax, data, off, title in [(axs[0], ic, off_ic, "ΔU_iç (iç enerji düzeltmesi)"),
                             (axs[1], etk, 0.0, "ΔU_etk (etkileşim düzeltmesi)")]:
    for k, leg in enumerate(data):
        ax.scatter(data[leg] - off, np.full(len(data[leg]), k), s=40, color=LEG[leg],
                   label=f"{LEG_TR[leg]}: ort {data[leg].mean()-off:+.1f}, σ={data[leg].std():.1f}")
        ax.errorbar(data[leg].mean() - off, k + 0.25, xerr=data[leg].std() / np.sqrt(len(data[leg])),
                    fmt="|", color=C["ink"], capsize=4, ms=12)
    ax.set_yticks([0, 1], ["kompleks", "çözücü"]); ax.set_ylim(-0.7, 1.9)
    ax.set_xlabel("kJ/mol"); ax.set_title(title); ax.legend(fontsize=8, loc="upper left")
plt.tight_layout(); plt.show()
se = lambda a: a.std(ddof=1) / np.sqrt(len(a))
print(f"ΔU_iç  ortalama farkı (kompleks − çözücü): {ic['complex'].mean()-ic['solvent'].mean():+.2f} ± "
      f"{np.hypot(se(ic['complex']), se(ic['solvent'])):.2f} kJ/mol")
print(f"ΔU_etk ortalama farkı (kompleks − çözücü): {etk['complex'].mean()-etk['solvent'].mean():+.2f} ± "
      f"{np.hypot(se(etk['complex']), se(etk['solvent'])):.2f} kJ/mol")
print(f"Ortalama DFT süresi (gömülü + vakum): {np.mean(lab.timings['dft']):.1f} s / yapı")
'''),
    md(r"""
**Bu grafik neden önemli?**

- **$\Delta U_\text{iç}$ çok geniş dağılıyor** (σ ≈ birkaç $kT$), ama **ortalaması iki bacakta aynı**
  (fark hata çubuğu içinde). Molekül aynı ve konformasyonları benzer olduğu için bu kısım bağlanma
  farkında sadeleşiyor. Biz de bu sadeleşmeyi varsayıp ML ve serbest enerji düzeltmesini
  $\Delta U_\text{etk}$ üzerinden yapacağız. (Bu varsayım esnek ya da konformasyonu bağlanınca değişen
  ligandlarda geçerli olmaz. Makale bu yüzden tam enerjiyi kullanıyor ve binlerce QM noktasıyla çalışıyor.)
- **$\Delta U_\text{etk}$ dar dağılıyor ve iki bacakta farklı.** Bağlanmaya giren kuantum düzeltmesi bu.
- Serbest enerji düzeltmesi $\Delta A = -kT \ln\langle e^{-\Delta U/kT}\rangle_\text{MM}$ ile hesaplanıyor.
  Dağılım $kT \approx 2.5$ kJ/mol'den çok genişse ortalama zor yakınsıyor ve çok sayıda QM hesabı gerekiyor.
  Makaledeki ML ve NEQ katmanlarının nedeni bu.
"""),
])

# =============================================================================
# 04
# =============================================================================
save("04_ML1_aktif_ogrenme.ipynb", [
    md(r"""
# 04 · ML1: aktif öğrenme ile QM/MM düzeltmesi
*Makale: Bölüm 2.2; Şekil 2 "Train ML1", Şekil 6–7.*

Pahalı DFT hesabını her yapı için yapmak yerine, $\Delta U_\text{etk}(x)$'i tahmin eden bir ML modeli eğitiyoruz:
- **Tanımlayıcılar:** her misafir atomunda çevrenin oluşturduğu elektrostatik potansiyel ve elektrik
  alanın karesi. Polarizasyon enerjisi yaklaşık $-\tfrac12\alpha|E|^2$ olduğu için bunlar fiziksel olarak anlamlı.
- **Model:** Gauss süreci (GP). Tahminle birlikte bir **belirsizlik** de veriyor.
- **Aktif öğrenme:** Her turda modelin *en emin olmadığı* yapıları DFT ile hesaplayıp eğitim kümesine
  ekliyoruz. Makale de aynı stratejiyi kullanıyor: "snapshots with high energy uncertainty".

Sonra ML1'i havuzdaki **tüm** yapılara uygulayıp serbest enerji düzeltmesini hesaplıyoruz:

$$\Delta G_\text{bağ}^\text{MM+ML1} = \Delta G_\text{bağ}^\text{MM} + \Delta A_\text{kompleks}^{\text{MM}\to\text{QM}} - \Delta A_\text{çözücü}^{\text{MM}\to\text{QM}}$$

(Burada $\Delta U = \Delta U_\text{etk}$; iç kısmın sadeleştiğini 03'te kontrol etmiştik.)

> **Basitleştirme:** Makale, kuvvet de veren bir ML potansiyeliyle MD yapıyor ve 10 ps'lik NEQ anahtarlama
> kullanıyor. Biz tek adımlı pertürbasyon (Zwanzig) kullanıyoruz. Bu, NEQ'nun "anlık anahtarlama" limiti.
"""),
    code(SETUP),
    code('''
from qfe import qmmm, ml
from pyscf import lib
lib.num_threads(os.cpu_count())
pools = {leg: (lambda d: (d["X"], d["B"]))(store.load_npz(f"pool_{leg}.npz")) for leg in ["complex", "solvent"]}
regions = {leg: qmmm.QMRegion(leg) for leg in pools}
lab = qmmm.Labeler(regions, pools)
lab.load_dict(store.load_json("labels.json", default={}))
print("Önbellekteki DFT etiketi:", len(lab.cache_dft))
with store.timed("04_tanimlayicilar"):
    feats = {leg: ml.descriptors(regions[leg], *pools[leg]) for leg in pools}
print({leg: f.shape for leg, f in feats.items()})
'''),
    code('''
p = config.P()
with store.timed("04_aktif_ogrenme") as rec:
    model1, train, val, hist = ml.active_learning(
        lab, feats, n_init=p["al_initial"], batch=p["al_batch"], rounds=p["al_rounds"], n_val=p["al_validation"])
    rec["n_dft"] = len(lab.cache_dft)
store.save_json("labels.json", lab.to_dict())
store.save_json("ml1.json", {"train": train, "val": val, "history": hist})
'''),
    code('''
fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.6))
axs[0].plot([h["n_train"] for h in hist], [h["rmse_val"] for h in hist], marker="o", color=C["blue"], label="doğrulama RMSE")
axs[0].axhline(hist[-1]["std_target"], color=C["muted"], ls="--", lw=1.2, label="hedefin std'si (model yokmuş gibi)")
axs[0].axhline(config.KT, color=C["aqua"], ls=":", lw=1.5, label="kT (2.5 kJ/mol)")
axs[0].set_xlabel("DFT eğitim noktası sayısı"); axs[0].set_ylabel("kJ/mol"); axs[0].set_title("Aktif öğrenme eğrisi")
axs[0].legend(fontsize=8)
Fv = np.array([feats[l][i] for l, i in val]); yv = np.array([lab.dU_dft(l, i) for l, i in val])
pv, sv = model1.predict(Fv, return_std=True)
for leg in ["complex", "solvent"]:
    m = np.array([l == leg for l, _ in val])
    axs[1].errorbar(yv[m] - yv.mean(), pv[m] - yv.mean(), yerr=sv[m], fmt="o", color=LEG[leg], ms=5, label=LEG_TR[leg])
lim = np.abs(np.concatenate([yv, pv]) - yv.mean()).max() * 1.15
axs[1].plot([-lim, lim], [-lim, lim], color=C["muted"], lw=1)
axs[1].set_xlabel("DFT ΔU (kJ/mol)"); axs[1].set_ylabel("ML1 tahmini (kJ/mol)"); axs[1].set_title("ML1: doğrulama kümesi")
axs[1].legend(fontsize=8); plt.tight_layout(); plt.show()
'''),
    md("### ML1 ile tüm havuzda serbest enerji düzeltmesi"),
    code('''
dU_ml1 = {leg: model1.predict(feats[leg]) for leg in pools}
corr1 = ml.binding_correction(dU_ml1["complex"], dU_ml1["solvent"])
r = store.load_json("results.json")
dG_mm, err_mm = r["dG_MM"]["dG_bind"], r["dG_MM"]["dG_bind_err"]
print(f"σ(ΔU) kompleks = {corr1['complex']['sigma']:.2f} kJ/mol ({corr1['complex']['sigma_over_kT']:.1f} kT), "
      f"çözücü = {corr1['solvent']['sigma']:.2f} kJ/mol")
for key, name in [("exp", "üstel ortalama (Zwanzig)"), ("cumulant", "2. mertebe kümülant")]:
    print(f"ΔΔG düzeltme [{name:24s}] = {corr1[key]:+6.2f} ± {corr1[key+'_err']:.2f} kJ/mol"
          f"  →  ΔG(MM+ML1) = {dG_mm + corr1[key]:6.2f} kJ/mol")
store.save_json("results.json", {"corr_ML1": corr1, "dG_MM_ML1": dG_mm + corr1["cumulant"],
                                 "dG_MM_ML1_err": float(np.hypot(err_mm, corr1["cumulant_err"]))})
store.save_npz("ml1_predictions.npz", complex=dU_ml1["complex"], solvent=dU_ml1["solvent"])

fig, ax = plt.subplots(figsize=(6.5, 3.4))
off = np.mean(np.concatenate(list(dU_ml1.values())))
bins = np.linspace(-4 * corr1["complex"]["sigma"], 4 * corr1["complex"]["sigma"], 35)
for leg in dU_ml1:
    ax.hist(dU_ml1[leg] - off, bins=bins, color=LEG[leg], alpha=0.55, label=LEG_TR[leg], edgecolor=C["surface"], linewidth=1.5)
ax.set_xlabel("ML1 ΔU − ortalama (kJ/mol)"); ax.set_ylabel("yapı sayısı")
ax.set_title("ML1 ile tüm havuz (makale Şekil 7a/b'nin karşılığı)"); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
'''),
    md("""
### Çapraz kontrol: ML olmadan, doğrudan DFT ortalaması
Doğrulama kümesi ve 03'teki yapılar havuzdan **rastgele** seçildi. Bu yüzden yalnızca onların DFT
değerleriyle de (ML'siz) bir tahmin yapılabilir. Daha az yapıya dayandığı için daha gürültülü olacak.
ML1, bu tahmini çok daha fazla yapıya genelleştiriyor.
"""),
    code('''
probe = store.load_json("probe.json", default={})
rand = {leg: sorted(set([i for l, i in val if l == leg] + probe.get(leg, []))) for leg in pools}
direct = {leg: np.array([lab.dU_dft(leg, i) for i in rand[leg]]) for leg in pools}
cd = ml.binding_correction(direct["complex"], direct["solvent"])
print(f"Doğrudan DFT ({len(rand['complex'])}+{len(rand['solvent'])} rastgele yapı): ΔΔG = {cd['cumulant']:+.2f} ± {cd['cumulant_err']:.2f} kJ/mol")
print(f"ML1 ({len(pools['complex'][0])}+{len(pools['solvent'][0])} yapı)               : ΔΔG = {corr1['cumulant']:+.2f} ± {corr1['cumulant_err']:.2f} kJ/mol")
store.save_json("results.json", {"corr_direct_DFT": cd})
'''),
    md(r"""
**Yorum:**
- **Üstel ortalama ile kümülant tahmini birbirine yakınsa** sonuç güvenilir. Uzaksa ($\sigma \gg kT$)
  örnekleme yetersiz demektir. Gerçek projede bu durumda daha fazla yapı ve NEQ anahtarlama gerekir.
- **ML1 neden mükemmel değil?** Öğrenme eğrisinde doğrulama hatası hedefin std'sinin altına iniyor ama
  sıfıra yaklaşmıyor. ΔU_etk'nin bir kısmı, kısa menzilli yük nüfuzu ve uzun menzilli elektrostatik gibi
  basit tanımlayıcılarımızın tam yakalayamadığı etkilerden geliyor. Makale bu yüzden binlerce noktayla
  eğitilmiş sinir ağı potansiyelleri kullanıyor. Bizim ölçeğimizde ML1 en azından ortalama kaymayı ve
  doğrusal eğilimi öğreniyor; çapraz kontrol de bunu doğruluyor.
- **Kaç DFT hesabı yaptık?** Yalnızca onlarca. Makale ise ~2000 hesap yapıyor ve QM bölgesi bir ağır metal
  içerdiği için her biri çok daha pahalı. **İlk HPC ihtiyacı burada ortaya çıkıyor** (bkz. 07).
"""),
])

# =============================================================================
# 05
# =============================================================================
save("05_kuantum_cekirdegi_ML2.ipynb", [
    md(r"""
# 05 · Kuantum çekirdeği (QM/QM/MM) ve ML2
*Makale: Bölüm 2.3–2.4, Denklem 7; Şekil 2 mor kutu ("with Traditional or Quantum Computer").*

DFT'nin hatası sistematik olarak düzeltilemiyor. Makale bu yüzden QM bölgesinin içinde daha küçük bir
**kuantum çekirdeği** tanımlıyor ve onu dalga fonksiyonu yöntemleriyle çözüyor:

- **CASCI(n, m):** Seçilen $m$ orbitalde $n$ elektron için Schrödinger denkleminin **tam** çözümü.
  Hilbert uzayı üstel büyüdüğü için klasik bilgisayarı zorlayan kısım bu. **Gelecekte bu kısmı kuantum
  bilgisayar (QPE) yapacak.**
- **NEVPT2:** Aktif uzay dışındaki "dinamik" korelasyonu ekleyen klasik pertürbasyon düzeltmesi.

Bizim misafirimiz küçük olduğu için kuantum çekirdeği bütün misafirden oluşuyor ve aktif uzay CAS(4,4).
Makalede ise çekirdek Ru d-orbitalleri; orada DFT içine Huzinaga gömmesi kullanılıyor.

**ML2 (transfer öğrenme):** Yüksek seviye hesap çok pahalı olduğu için **az** sayıda yapıda yapılıyor.
ML2, ML1'in üzerine yalnızca *farkı* öğreniyor:
$$\Delta U^\text{yüksek}_\text{etk}(x) \approx \text{ML1}(x) + \text{ML2}(x),\qquad \text{ML2} \leftarrow \Delta U_\text{etk}^\text{CASCI+NEVPT2} - \Delta U_\text{etk}^\text{DFT}$$

Yani ML2, kuantum çekirdeği seviyesinin misafir–çevre etkileşimini DFT'den ne kadar farklı gördüğünü öğreniyor.
Her etiket için hem gömülü hem vakum yüksek seviye hesabı yapılıyor.
"""),
    code(SETUP),
    code('''
from qfe import qmmm, ml
from pyscf import lib
lib.num_threads(os.cpu_count())
pools = {leg: (lambda d: (d["X"], d["B"]))(store.load_npz(f"pool_{leg}.npz")) for leg in ["complex", "solvent"]}
regions = {leg: qmmm.QMRegion(leg) for leg in pools}
lab = qmmm.Labeler(regions, pools)
lab.load_dict(store.load_json("labels.json"))
feats = {leg: ml.descriptors(regions[leg], *pools[leg]) for leg in pools}
ml1 = store.load_json("ml1.json")
train1 = [tuple(t) for t in ml1["train"]]
model1 = ml.GPModel().fit(np.array([feats[l][i] for l, i in train1]), np.array([lab.dU_dft(l, i) for l, i in train1]))
print("ML1 yeniden kuruldu,", len(train1), "eğitim noktası")
'''),
    md("### Tek bir yapıda seviyeleri karşılaştıralım"),
    code('''
leg, i = train1[0]
fr = regions[leg].frame(pools[leg][0][i], pools[leg][1][i])
with store.timed("05_tek_yuksek_seviye"):
    h, mf, mc = qmmm.energy_high(fr, regions[leg].guest_symbols(), return_objects=True)
print(f"HF           : {h['e_hf']:.6f} Ha")
print(f"CASCI{tuple(config.P()['cas'])} : {h['e_casci']:.6f} Ha   (statik korelasyon {1000*(h['e_casci']-h['e_hf']):.3f} mHa)")
print(f"+ NEVPT2     : {h['e_total']:.6f} Ha   (dinamik korelasyon {1000*h['e_nevpt2_corr']:.1f} mHa)")
occ = mf.mo_occ; nocc = int((occ > 0).sum()); ncas = config.P()["cas"][0]
print(f"Aktif orbitaller: HOMO-{ncas//2-1} … LUMO+{ncas//2-1}, enerjiler (Ha):",
      np.round(mf.mo_energy[nocc-ncas//2:nocc+ncas//2], 3))
'''),
    md("""
**Gözlem:** Siklopentanol kapalı kabuklu, "kolay" bir molekül. Statik korelasyon (CASCI − HF)
çok küçük, enerjinin büyük kısmı dinamik korelasyondan geliyor. Bu, makaledeki **Şekil 1'in sol alt
köşesi**. Bu tür sistemlerde kuantum avantajı beklenmiyor. Makalenin Ru kompleksi ise açık kabuklu
(spin dubleti) ve çok referanslı, yani kuantum bilgisayarın asıl fark yaratacağı türden bir sistem.
Burada amacımız avantaj göstermek değil, **iş akışının nasıl kurulduğunu** göstermek.

### ML2 eğitim verisi
"""),
    code('''
p = config.P()
rng = np.random.default_rng(2)
ml2_set = []
for leg in pools:
    cand = [t for t in train1 if t[0] == leg]
    sel = rng.choice(len(cand), size=min(len(cand), p["ml2_points"] // 2), replace=False)
    ml2_set += [cand[j] for j in sel]
with store.timed("05_yuksek_seviye_etiketler") as rec:
    delta = np.array([lab.high_minus_dft(l, i) for l, i in ml2_set])
    rec["n"] = len(ml2_set)
store.save_json("labels.json", lab.to_dict())
store.save_json("ml2.json", {"train": ml2_set})
print(f"{len(ml2_set)} yapı, ortalama yüksek seviye süresi (gömülü + vakum) {np.mean(lab.timings['high']):.1f} s")
print(f"δ = ΔU_etk(yüksek) − ΔU_etk(DFT): ortalama {delta.mean():+.2f}, std {delta.std():.2f} kJ/mol")
'''),
    code('''
model2 = ml.GPModel(seed=1).fit(np.array([feats[l][i] for l, i in ml2_set]), delta)
dU_ml2 = {leg: model1.predict(feats[leg]) + model2.predict(feats[leg]) for leg in pools}
corr2 = ml.binding_correction(dU_ml2["complex"], dU_ml2["solvent"])
r = store.load_json("results.json")
dG_mm, err_mm = r["dG_MM"]["dG_bind"], r["dG_MM"]["dG_bind_err"]
for key, name in [("exp", "Zwanzig"), ("cumulant", "kümülant")]:
    print(f"ΔΔG(MM→QM/QM/MM) [{name:9s}] = {corr2[key]:+6.2f} ± {corr2[key+'_err']:.2f}  →  "
          f"ΔG(MM+ML2) = {dG_mm + corr2[key]:6.2f} kJ/mol")
store.save_json("results.json", {"corr_ML2": corr2, "dG_MM_ML2": dG_mm + corr2["cumulant"],
                                 "dG_MM_ML2_err": float(np.hypot(err_mm, corr2["cumulant_err"]))})
fig, ax = plt.subplots(figsize=(6.2, 3.3))
for leg in pools:
    m = np.array([l == leg for l, _ in ml2_set])
    ax.scatter(delta[m] - delta.mean(), np.full(m.sum(), 0 if leg == "complex" else 1), color=LEG[leg], s=40, label=LEG_TR[leg])
ax.set_yticks([0, 1], ["kompleks", "çözücü"]); ax.set_ylim(-0.7, 1.7)
ax.set_xlabel("δ − ortalama (kJ/mol)")
ax.set_title("Kuantum çekirdeği düzeltmesinin dağılımı (makale Şekil 7, alt sıra)"); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
'''),
    md("""
**Dikkat:** Makalede en kritik nokta, δ'nın dağılımının **dar** olması (Şekil 7, alt sıra). Dar olduğunda
ML2 onu az veriyle öğrenebiliyor. Makale bunun için ~10.000 QM/QM/MM enerjisi kullanıyor; bunlar
gelecekte kuantum bilgisayardan gelecek. Bir sonraki notebook'ta bu enerjilerin kuantum algoritmasıyla
nasıl hesaplanacağını simüle ediyoruz.
"""),
])

# =============================================================================
# 06
# =============================================================================
save("06_kuantum_algoritmalari.ipynb", [
    md(r"""
# 06 · Kuantum çekirdeğini kuantum algoritmalarıyla çözmek (simülatörde)
*Makale: Bölüm 3, Şekil 4d.*

05'teki CASCI enerjisini bu sefer bir **kuantum bilgisayarın yapacağı gibi** hesaplıyoruz:

1. Aktif uzay Hamiltoniyenini fermiyon operatörlerine, sonra **Jordan–Wigner** dönüşümüyle **qubit**
   operatörlerine çeviriyoruz: $H = \sum_i h_i P_i$ ($P_i$: Pauli dizileri). $N$ orbital için $2N$ qubit gerekiyor.
2. **VQE:** Yakın vadeli, değişkenli yöntem. Makale büyük sistemlerde sorunlu olduğunu belirtiyor.
3. **QPE (kuantum faz tahmini):** Makalenin hedeflediği, hata düzeltmeli bilgisayar algoritması.
   Bir başlangıç durumu $|\psi\rangle$ hazırlanıyor, $e^{-iHt}$ uygulanıyor ve faz ölçülüyor.
   Doğru enerjiyi bulma olasılığı $|\langle\psi|\Psi_0\rangle|^2$ örtüşmesine bağlı.
4. **Trotter:** Gerçek donanımda $e^{-iHt}$ tam uygulanamıyor; Pauli terimlerinin üstellerinin çarpımına bölünüyor.

Hepsi klasik **durum vektörü simülasyonu**. Bu yaklaşım ~25–30 qubit'e kadar işliyor.
"""),
    code(SETUP),
    code('''
from qfe import qmmm, quantum
from qfe.plotting import SERIES
pools = {leg: (lambda d: (d["X"], d["B"]))(store.load_npz(f"pool_{leg}.npz")) for leg in ["complex"]}
reg = qmmm.QMRegion("complex")
fr = reg.frame(pools["complex"][0][0], pools["complex"][1][0])
mf = qmmm.hf_in_field(fr, reg.guest_symbols())
ncas, nel = config.P()["cas"]
nq = 2 * ncas
h1, h2, ecore = quantum.active_space_integrals(mf, ncas, nel)
qop = quantum.qubit_hamiltonian(h1, h2, ecore)
print(f"CAS({nel},{ncas}) → {nq} qubit, {len(qop.terms)} Pauli terimi")
for term, c in list(qop.terms.items())[:6]:
    print("  ", "".join(f"{p}{q}" for q, p in term) or "I", f"{c.real:+.5f}")
'''),
    md("### 1) Klasik referans: tam köşegenleştirme = CASCI olmalı"),
    code('''
from pyscf import mcscf
HN, idx, hf_pos, Hs = quantum.number_sector(qop, nq, nel)
w, v = quantum.exact_ground(HN)
mc = mcscf.CASCI(mf, ncas, nel); mc.verbose = 0; e_casci = mc.kernel()[0]
print(f"{nel}-elektron alt uzayının boyutu: {len(HN)} (tam Hilbert uzayı 2^{nq} = {2**nq})")
print(f"Qubit Hamiltoniyeni özdeğeri : {w[0]:.8f} Ha")
print(f"PySCF CASCI                  : {e_casci:.8f} Ha   fark = {abs(w[0]-e_casci):.1e}")
print(f"HF başlangıç durumu örtüşmesi |<HF|Ψ0>|² = {abs(v[hf_pos,0])**2:.5f}")
'''),
    md("### 2) VQE (UCCSD ansatz) — gürültüsüz ve atış gürültülü"),
    code('''
with store.timed("06_vqe"):
    vq = quantum.vqe_uccsd(qop, nq, nel)
    vq_noisy = quantum.vqe_uccsd(qop, nq, nel, shots=2000, maxiter=300)
print(f"VQE (gürültüsüz) : {vq['energy']:.8f} Ha, hata {1000*(vq['energy']-w[0]):.4f} mHa, {vq['n_params']} parametre, {vq['n_evals']} enerji ölçümü")
print(f"VQE (2000 atış)  : {vq_noisy['energy']:.8f} Ha, hata {1000*(vq_noisy['energy']-w[0]):.4f} mHa, {vq_noisy['n_evals']} enerji ölçümü")
fig, ax = plt.subplots(figsize=(6.5, 3.4))
ax.plot(1000 * (vq["trace"] - w[0]), color=C["blue"], label="gürültüsüz")
ax.plot(1000 * (vq_noisy["trace"] - w[0]), color=C["orange"], lw=1, alpha=0.8, label="2000 atış/ölçüm")
ax.axhline(1.6, color=C["muted"], ls="--", lw=1.2, label="kimyasal doğruluk (1.6 mHa)")
ax.set_yscale("symlog", linthresh=1e-3); ax.set_xlabel("enerji değerlendirmesi"); ax.set_ylabel("E − E_tam (mHa)")
ax.set_title("VQE yakınsaması"); ax.legend(fontsize=8); plt.tight_layout(); plt.show()
'''),
    md("""
**Not:** Gerçek donanımda her "enerji değerlendirmesi" binlerce devre çalıştırması (atış) demek.
Atış gürültüsü, optimizasyonu ve son doğruluğu sınırlıyor. Makalenin VQE yerine QPE'yi tercih
etmesinin bir nedeni de bu.

### 3) Kuantum faz tahmini (QPE) — devre simülasyonu
"""),
    code('''
psi_hf = np.zeros(len(HN)); psi_hf[hf_pos] = 1.0
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
errs = []
for k, bits in enumerate([6, 8, 10, 12]):
    q = quantum.qpe_simulate(HN, psi_hf, bits)
    e_est, _ = quantum.qpe_shots(q, 200, seed=k)
    errs.append((bits, abs(e_est - w[0]), q["resolution"]))
    if bits in (8, 12):
        m = (q["energies"] > w[0] - 0.03) & (q["energies"] < w[0] + 0.03)
        axs[0].plot(1000 * (q["energies"][m] - w[0]), q["prob"][m], color=SERIES[k % 4], lw=1.4,
                    marker="o", ms=3, label=f"{bits} yardımcı qubit")
axs[0].axvline(0, color=C["muted"], ls="--", lw=1)
axs[0].set_xlabel("ölçülen E − E₀ (mHa)"); axs[0].set_ylabel("olasılık"); axs[0].set_title("QPE ölçüm dağılımı (HF başlangıç, yakın plan)")
axs[0].legend(fontsize=8)
b, e, r_ = zip(*errs)
axs[1].semilogy(b, np.array(e) * 1000, marker="o", color=C["blue"], label="200 atışla tahmin hatası")
axs[1].semilogy(b, np.array(r_) * 1000, marker="s", ls="--", color=C["muted"], label="çözünürlük (1 bit)")
axs[1].axhline(1.6, color=C["aqua"], ls=":", label="kimyasal doğruluk")
axs[1].set_xlabel("yardımcı qubit (bit) sayısı"); axs[1].set_ylabel("mHa"); axs[1].set_title("QPE hassasiyeti")
axs[1].legend(fontsize=8); plt.tight_layout(); plt.show()
'''),
    md("""
Her ek bit hassasiyeti ~2 kat artırıyor ama gereken evrim süresini ($U^{2^k}$) de 2 kat uzatıyor.
Makaledeki kapı sayılarının ($10^{9}$–$10^{13}$ Toffoli) bu kadar büyük olmasının nedeni bu.

### 4) Trotter hatası: $e^{-iHt}$'yi kapılarla yaklaşık uygulamak
"""),
    code('''
psi0 = np.zeros(2 ** nq, complex); psi0[idx] = v[:, 0]
t_total = 10.0
steps = [2, 4, 8, 16, 32, 64, 128]
with store.timed("06_trotter"):
    tr1 = quantum.trotter_phase_error(qop, nq, psi0, w[0], t_total, steps, order=1)
fig, ax = plt.subplots(figsize=(6.3, 3.5))
ax.loglog([t_total / t["steps"] for t in tr1], [1000 * t["error"] + 1e-12 for t in tr1], marker="o", color=C["blue"], label="1. mertebe Trotter")
for t in tr1[::2]:
    ax.annotate(f"{t['n_exponentials']:,} kapı", (t_total / t["steps"], 1000 * t["error"]), textcoords="offset points",
                xytext=(6, -10), fontsize=8, color=C["ink2"])
ax.axhline(1.6, color=C["muted"], ls="--", lw=1.2, label="kimyasal doğruluk")
ax.set_xlabel("Trotter adım boyu δt (atomik birim)"); ax.set_ylabel("enerji hatası (mHa)")
ax.set_title(f"Trotter hatası (t = {t_total:.0f} a.u. evrim)"); ax.legend(fontsize=8); plt.tight_layout(); plt.show()
'''),
    md("""
Adım küçüldükçe hata hızla düşüyor ama **kapı sayısı** ters orantılı olarak artıyor. Bu küçük ve zayıf
korelasyonlu sistemde hata zaten küçük. Aktif uzay büyüdükçe Hamiltoniyenin normu λ büyüyor ve aynı
doğruluk için çok daha küçük adım gerekiyor (bkz. 07, kaynak tahmini). Makale bu Trotter hatasını
36 qubit'e kadar GPU statevector simülasyonuyla ölçüp daha büyük sistemlere ekstrapole ediyor (ref. 30).
"""),
    md("### 5) Başlangıç durumu örtüşmesi — aktif uzay büyüdükçe (makale Şekil 4d)"),
    code('''
with store.timed("06_ortusme_taramasi"):
    ov = quantum.hf_overlap_scan(mf, config.P()["overlap_scan"])
store.save_json("quantum.json", {"overlap_scan": ov})
fig, ax = plt.subplots(figsize=(6.3, 3.4))
ax.plot([o["n_orb"] for o in ov], [o["overlap_hf"] for o in ov], marker="o", color=C["green"], label="HF durumu")
ax.set_xlabel("aktif orbital sayısı"); ax.set_ylabel("|⟨HF|Ψ₀⟩|²"); ax.set_title("QPE başlangıç durumu kalitesi")
ax.legend(fontsize=8); plt.tight_layout(); plt.show()
for o in ov:
    print(f"CAS({o['n_orb']},{o['n_orb']}): {o['n_qubits']:2d} qubit, {o['n_determinants']:>9,d} determinant, "
          f"örtüşme {o['overlap_hf']:.4f}, klasik süre {o['seconds']:.2f} s")
'''),
    md("""
**Yorum:**
- Örtüşme, aktif uzay büyüdükçe düşüyor. Makalede Ru sistemi için 30 orbitalde HF örtüşmesi ~0.73'e iniyor;
  daha iyi başlangıç durumlarıyla (MPS, Slater determinant toplamı) ~0.9–0.97'ye çıkıyor.
- **Determinant sayısı üstel büyüyor.** ~18 orbitalden sonra klasik tam çözüm (FCI) laptopta imkânsız hale
  geliyor. Bu noktadan sonra **DMRG (HPC)** ya da **kuantum bilgisayar** gerekiyor.

### 6) Gerçek kuantum donanımına köprü
Aşağıdaki hücre Hamiltoniyeni Qiskit'in okuyabileceği bir JSON dosyasına yazıyor. IBM Quantum'un ücretsiz
planında küçük bir örnek (ör. CAS(2,2) → 4 qubit) için tahmin yapılabilir (bkz. README).
"""),
    code('''
import json
small = quantum.qubit_hamiltonian(*quantum.active_space_integrals(mf, 2, 2))
for name, op, n in [("H_cas22_4qubit.json", small, 4), (f"H_cas{nel}{ncas}_{nq}qubit.json", qop, nq)]:
    pl = [(lab, [c.real, c.imag]) for lab, c in quantum.pauli_list(op, n)]
    (store.path(name)).write_text(json.dumps(pl))
    print("yazıldı:", store.path(name), f"({len(pl)} terim)")
'''),
])

# =============================================================================
# 07
# =============================================================================
save("07_kaynak_tahmini_ve_yol_haritasi.ipynb", [
    md("""
# 07 · Kaynak tahmini, hata yayılımı ve yol haritası
*Makale: Bölüm 3, Şekil 4a–c ve Sonuçlar.*

Bu notebook üç soruyu yanıtlıyor:
1. Aktif uzay büyüdükçe kuantum bilgisayar için kaç qubit ve kaç kapı gerekiyor?
2. Kuantum enerjilerinin hatası ΔG'yi ne kadar etkiliyor? (Kuantum bilgisayarın ne kadar *doğru* olması gerekiyor?)
3. Ölçtüğümüz sürelere göre gerçek projede **nerede HPC, nerede kuantum bilgisayar** gerekecek?
"""),
    code(SETUP),
    code('''
from qfe import qmmm, quantum, ml
pools = {leg: (lambda d: (d["X"], d["B"]))(store.load_npz(f"pool_{leg}.npz")) for leg in ["complex", "solvent"]}
regions = {leg: qmmm.QMRegion(leg) for leg in pools}
fr = regions["complex"].frame(pools["complex"][0][0], pools["complex"][1][0])
mf = qmmm.hf_in_field(fr, regions["complex"].guest_symbols())
with store.timed("07_kaynak_tahmini"):
    rs = quantum.resource_estimates(mf, config.P()["resource_scan"])
store.save_json("quantum.json", {"resources": rs})
for r in rs:
    print(f"{r['n_orb']:2d} orbital | {r['n_qubits']:2d} qubit | {r['n_pauli_terms']:6d} Pauli terimi | λ={r['lambda_ha']:7.1f} Ha | "
          f"qubitization ~{r['qubitization_toffoli_naive']:.1e} Toffoli | qDRIFT ~{r['qdrift_exponentials']:.1e} üstel")
'''),
    code('''
n = np.array([r["n_orb"] for r in rs])
fig, axs = plt.subplots(1, 2, figsize=(11, 3.7))
axs[0].semilogy(n, [r["n_pauli_terms"] for r in rs], marker="o", color=C["blue"], label="Pauli terimi (L)")
axs[0].semilogy(n, [r["lambda_ha"] for r in rs], marker="s", color=C["orange"], label="1-norm λ (Ha)")
axs[0].set_xlabel("aktif orbital"); axs[0].set_title("Hamiltoniyen büyüklüğü"); axs[0].legend(fontsize=8)
axs[1].semilogy(n, [r["qubitization_toffoli_naive"] for r in rs], marker="o", color=C["green"], label="qubitization (naif SELECT)")
axs[1].semilogy(n, [r["qdrift_exponentials"] for r in rs], marker="^", color=C["blue"], label="qDRIFT")
# Makaledeki (Şekil 4b) Ru sistemi için okunan kaba değerler, kıyas için:
axs[1].scatter([30, 60], [1e10, 1e11], marker="*", s=140, color=C["red"], zorder=5, label="makale Ru (qubitization/DF, Şekil 4b)")
axs[1].set_xlabel("aktif orbital"); axs[1].set_ylabel("kapı sayısı / devre"); axs[1].set_title("QPE maliyeti (ε = 1.6 mHa)")
axs[1].legend(fontsize=8); plt.tight_layout(); plt.show()
'''),
    md("""
**Okuma:** Bizim formüllerimiz yalnızca **mertebe** veriyor. Makaledeki "double factorization" gibi teknikler
sayıları birkaç mertebe azaltıyor. Yine de mesaj aynı: ~30–60 orbitallik ilginç bir çekirdek için
**~100–1000 mantıksal qubit ve 10¹⁰ civarı kapı** gerekiyor. Bu ancak hata düzeltmeli, hızlı
(kapı başına ≲10⁻⁷ s) bir makineyle mümkün (makale Bölüm 3).

### Kuantum enerjilerinin hatası ΔG'ye nasıl yansıyor?
Yüksek seviye enerjilere (gelecekte kuantum bilgisayardan gelecek) σ büyüklüğünde rastgele hata
ekliyoruz, ML2'yi yeniden eğitiyoruz ve ΔG'nin ne kadar saptığına bakıyoruz.
"""),
    code('''
lab = qmmm.Labeler(regions, pools); lab.load_dict(store.load_json("labels.json"))
feats = {leg: ml.descriptors(regions[leg], *pools[leg]) for leg in pools}
train1 = [tuple(t) for t in store.load_json("ml1.json")["train"]]
ml2_set = [tuple(t) for t in store.load_json("ml2.json")["train"]]
model1 = ml.GPModel().fit(np.array([feats[l][i] for l, i in train1]), np.array([lab.dU_dft(l, i) for l, i in train1]))
base1 = {leg: model1.predict(feats[leg]) for leg in pools}
delta = np.array([lab.high_minus_dft(l, i) for l, i in ml2_set])
F2 = np.array([feats[l][i] for l, i in ml2_set])

def dG_corr(d):
    m2 = ml.GPModel(seed=1).fit(F2, d)
    return ml.binding_correction(base1["complex"] + m2.predict(feats["complex"]),
                                 base1["solvent"] + m2.predict(feats["solvent"]))["cumulant"]

ref = dG_corr(delta)
sigmas = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]
rng = np.random.default_rng(0)
spread = []
with store.timed("07_hata_yayilimi"):
    for s in sigmas:
        vals = [dG_corr(delta + rng.normal(0, s, len(delta))) for _ in range(8 if s > 0 else 1)]
        spread.append(np.sqrt(np.mean((np.array(vals) - ref) ** 2)))
fig, ax = plt.subplots(figsize=(6.3, 3.5))
ax.plot(sigmas, spread, marker="o", color=C["violet"], label="ΔG sapması (RMS)")
ax.axvline(1.6e-3 * config.HARTREE_TO_KJ, color=C["muted"], ls="--", lw=1.2, label="kimyasal doğruluk ≈ 4.2 kJ/mol")
ax.axhline(config.KT, color=C["aqua"], ls=":", label="kT")
ax.set_xlabel("kuantum enerjisi başına rastgele hata σ (kJ/mol)"); ax.set_ylabel("ΔG hatası (kJ/mol)")
ax.set_title("Kuantum hesap hatası → serbest enerji hatası"); ax.legend(fontsize=8); plt.tight_layout(); plt.show()
store.save_json("results.json", {"noise_propagation": {"sigma": sigmas, "dG_rms": spread}})
'''),
    md("""
**Yorum:** Rastgele hatalar ML ve ortalama alma sayesinde kısmen sönümleniyor. **Sistematik** hatalar
(örneğin yetersiz aktif uzay ya da kötü başlangıç durumu) ise sönümlenmiyor. Makalenin "kimyasal doğruluk
(~1 kJ/mol bölgesi)" hedefi bu yüzden konuluyor. Eğitim noktası azsa ML2 hatayı büyütebiliyor da.

### Sonuçların özeti
"""),
    code('''
r = store.load_json("results.json")
rows = [("MM (FEP)", r["dG_MM"]["dG_bind"], r["dG_MM"]["dG_bind_err"]),
        ("MM + ML1 (QM/MM, DFT)", r.get("dG_MM_ML1"), r.get("dG_MM_ML1_err")),
        ("MM + ML2 (QM/QM/MM)", r.get("dG_MM_ML2"), r.get("dG_MM_ML2_err")),
        ("Deney", r["exp_dG_kJ"], r["exp_dG_err_kJ"])]
for name, v, e in rows:
    print(f"{name:24s} {v:7.2f} ± {e:.2f} kJ/mol")
fig, ax = plt.subplots(figsize=(7.5, 3.2))
cols = [C["blue"], C["blue"], C["blue"], C["ink"]]
for k, (name, v, e) in enumerate(rows):
    ax.errorbar(v, k, xerr=e, fmt="o", color=cols[k], ms=7, capsize=4)
ax.axvline(r["exp_dG_kJ"], color=C["muted"], ls="--", lw=1)
ax.set_yticks(range(len(rows)), [x[0] for x in rows]); ax.invert_yaxis()
ax.set_xlabel("ΔG_bağlanma (kJ/mol)"); ax.set_title("Katman katman ΔG (makale Şekil 3b karşılığı)")
plt.tight_layout(); plt.show()
'''),
    md("""
### Nerede HPC, nerede kuantum bilgisayar? (ölçülen sürelerden ölçekleme)
Bu demoda ölçtüğümüz süreleri, makaledeki gerçek sistemin boyutuna ölçekliyoruz. DFT maliyeti
yaklaşık olarak (baz fonksiyonu)³ ile büyüyor; Ru ilacı NKP1339'un baz fonksiyonu sayısını PySCF'le sayıyoruz.
"""),
    code('''
from pyscf import gto
t = store.load_json("timings.json")
lab_t = store.load_json("labels.json")["timings"]
t_dft, t_high = float(np.mean(lab_t["dft"])), float(np.mean(lab_t["high"]))
nb_ours = r["qm_nbasis"]
# NKP1339 anyonu [RuCl4(indazol)2]^-  : C14 H12 N4 Cl4 Ru  (baz sayısı geometriden bağımsız)
atoms = ["Ru"] + ["Cl"] * 4 + ["N"] * 4 + ["C"] * 14 + ["H"] * 12
mol = gto.M(atom=[(a, (0, 0, 5.0 * k)) for k, a in enumerate(atoms)], basis="def2-svp",
            ecp={"Ru": "def2-svp"}, charge=-1, spin=1, verbose=0)
nb_ru = mol.nao
scale = (nb_ru / nb_ours) ** 3
ns_day = r.get("md_ns_per_day", float("nan"))
fep_ns_paper = 2 * 20 * 5      # 2 bacak × ~20 λ × ~5 ns (tipik protein-ligand FEP)
tab = [
    ("MM-FEP (protein, ~50k atom)", f"{fep_ns_paper} ns MD", f"{fep_ns_paper/ max(ns_day,1e-9) * 10:.1f} gün (bu makinede, 10× büyük sistem)", "GPU kümesi"),
    ("QM/MM DFT, ~2000 yapı", f"{t_dft:.0f} s/yapı ({r.get('n_threads', 1)} çekirdek)", f"~{2000 * t_dft * scale * r.get('n_threads', 1) / 3600:,.0f} çekirdek-saat (Ru ölçekli)", "HPC (CPU)"),
    ("QM/QM/MM, ~10000 yapı", f"{t_high:.0f} s/yapı ({r.get('n_threads', 1)} çekirdek)", f"≥{10000 * t_high * scale * r.get('n_threads', 1) / 3600:,.0f} çekirdek-saat (DMRG çok daha pahalı)", "HPC → kuantum bilgisayar"),
    ("Kuantum çekirdek ≥ 30 orbital", "klasik FCI imkânsız", "—", "Hata düzeltmeli kuantum bilgisayar"),
]
print(f"Baz fonksiyonu: bizim misafir {nb_ours}, NKP1339 {nb_ru}  → DFT maliyet çarpanı ~{scale:.0f}×\\n")
print(f"{'Adım':34s} | {'Ölçülen':22s} | {'Gerçek projeye ölçekli':44s} | Nerede")
print("-" * 130)
for row in tab:
    print(f"{row[0]:34s} | {row[1]:22s} | {row[2]:44s} | {row[3]}")
'''),
    md("""
## Ekip için yol haritası (özet)

| Aşama | Laptop (bu demo) | HPC gerektiğinde | Kuantum bilgisayar gerektiğinde |
|---|---|---|---|
| Sistem hazırlama, MM-FEP | ✓ küçük konakçı–misafir | Protein–ligand, çok sayıda ligand | — |
| QM/MM etiketleri (DFT) | Onlarca, ≤ ~30 atomluk QM bölgesi | Binlerce, metal içeren bölge | — |
| ML1/ML2 eğitimi | GP, yüzlerce nokta | Sinir ağı potansiyeli, GPU | — |
| Kuantum çekirdek (CASCI) | ≤ ~16 orbital (tam FCI) | DMRG ile ~30–60 orbital | ≥ 30–60 orbital, binlerce yapı, her biri ≤ 20 dk |
| QPE/VQE denemeleri | Simülatör ≤ ~24 qubit | GPU statevector ≤ ~36 qubit (makale ref. 30) | Gerçek donanım: küçük demo (IBM, 4–8 qubit); asıl hedef hata düzeltmeli makine |

**Sonraki adımlar için öneriler:**
1. `demo` ön ayarını GPU'da çalıştırıp sonuçları birlikte yorumlayın.
2. Misafiri değiştirin (ör. guest 5 sikloheptanol, guest 4 trans-4-metilsiklohekzanol) ve deneyle korelasyona bakın.
3. Açık kabuklu bir metal kompleksine geçmek: asıl kuantum avantajı sorusu orada başlıyor (makalenin Ru örneği).
4. HPC başvurusunda (UHeM öğrenci projesi) 03–05 adımlarının ölçeklenmiş versiyonunu proje olarak tanımlayın.
"""),
])
