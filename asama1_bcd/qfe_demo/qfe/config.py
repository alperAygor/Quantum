"""
Ayar dosyası — tüm notebook'lar buradan okur.

Ön ayarı (preset) değiştirmek için iki yol var:
  1) Ortam değişkeni:  export QFE_PRESET=demo
  2) Notebook'ta:      from qfe import config; config.set_preset("demo")

Ön ayarlar:
  "test"     : Kodun çalıştığını görmek için (CPU'da ~30-40 dk). Sonuçlar istatistiksel olarak anlamsız.
  "demo"     : Laptop + NVIDIA GPU'da ekibe gösterim için (~1-2 saat). Kabaca anlamlı sonuçlar.
  "standart" : Daha güvenilir sonuçlar (GPU'da birkaç saat + QM için gece boyu).
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
DATA.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

# ----------------------------------------------------------------------------
# Fiziksel sabitler
# ----------------------------------------------------------------------------
T_KELVIN = 298.15
KB_KJ = 0.0083144626          # kJ/(mol K)
KT = KB_KJ * T_KELVIN         # kJ/mol  (~2.479)
HARTREE_TO_KJ = 2625.499639
KCAL_TO_KJ = 4.184
V0_NM3 = 1.66054              # standart durum hacmi (1 M) nm^3

# ----------------------------------------------------------------------------
# Sistem: beta-siklodekstrin (bCD) + siklopentanol (MobleyLab CD Set 2, guest 3)
# Deneysel veri: Rekharsky et al., J. Phys. Chem. B 1997, 101, 87-100
# ----------------------------------------------------------------------------
SYSTEM = dict(
    name="bCD-siklopentanol",
    host_resname="MGO",
    guest_resname="MOL",
    guest_smiles="OC1CCCC1",
    orientation="s",                       # 's': OH ikincil (geniş) yüze bakıyor
    exp_dG_kcal=-3.050, exp_dG_sem_kcal=0.010,
    source_base="https://raw.githubusercontent.com/MobleyLab/benchmarksets/master/input_files/cd-set2",
)

PRESETS = {
    "test": dict(
        platform_preference=["CUDA", "OpenCL", "CPU"],
        equil_npt_ps=20,
        lambda_elec=[1.0, 0.75, 0.5, 0.25, 0.0],
        lambda_ster=[1.0, 0.85, 0.7, 0.55, 0.42, 0.3, 0.2, 0.1, 0.0],
        fep_equil_ps=4, fep_prod_ps=12, fep_sample_ps=0.4,
        pool_frames=120, pool_stride_ps=0.5,
        qm_cutoff_A=14.0, qm_basis="def2-svp", qm_xc="pbe",
        al_initial=8, al_batch=4, al_rounds=3, al_validation=6,
        ml2_points=6,
        cas=(4, 4),
        overlap_scan=[2, 4, 6, 8, 10],
        resource_scan=[2, 4, 6, 8],
    ),
    "demo": dict(
        platform_preference=["CUDA", "OpenCL", "CPU"],
        equil_npt_ps=100,
        lambda_elec=[1.0, 0.75, 0.5, 0.25, 0.0],
        lambda_ster=[1.0, 0.85, 0.7, 0.55, 0.42, 0.3, 0.2, 0.1, 0.0],
        fep_equil_ps=20, fep_prod_ps=200, fep_sample_ps=1.0,
        pool_frames=400, pool_stride_ps=2.0,
        qm_cutoff_A=16.0, qm_basis="def2-svp", qm_xc="pbe",
        al_initial=16, al_batch=8, al_rounds=5, al_validation=16,
        ml2_points=16,
        cas=(4, 4),
        overlap_scan=[2, 4, 6, 8, 10, 12],
        resource_scan=[2, 4, 6, 8, 10],
    ),
    "standart": dict(
        platform_preference=["CUDA", "OpenCL", "CPU"],
        equil_npt_ps=500,
        lambda_elec=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0],
        lambda_ster=[1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.42, 0.35, 0.28, 0.21, 0.14, 0.07, 0.0],
        fep_equil_ps=100, fep_prod_ps=1000, fep_sample_ps=1.0,
        pool_frames=1000, pool_stride_ps=2.0,
        qm_cutoff_A=16.0, qm_basis="def2-svp", qm_xc="pbe",
        al_initial=32, al_batch=16, al_rounds=8, al_validation=32,
        ml2_points=40,
        cas=(6, 6),
        overlap_scan=[2, 4, 6, 8, 10, 12, 14],
        resource_scan=[2, 4, 6, 8, 10],
    ),
}

_preset_name = os.environ.get("QFE_PRESET", "demo")


def set_preset(name: str):
    global _preset_name
    if name not in PRESETS:
        raise ValueError(f"Bilinmeyen ön ayar: {name}. Seçenekler: {list(PRESETS)}")
    _preset_name = name
    os.environ["QFE_PRESET"] = name


def preset_name() -> str:
    return _preset_name


def P() -> dict:
    """Aktif ön ayarın parametreleri."""
    return PRESETS[_preset_name]


def exp_dG_kj():
    return SYSTEM["exp_dG_kcal"] * KCAL_TO_KJ, SYSTEM["exp_dG_sem_kcal"] * KCAL_TO_KJ
