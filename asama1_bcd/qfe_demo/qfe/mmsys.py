"""Aşama 1-2: gerçek veriyi indirme, MM sistemlerini (kompleks / çözücü bacağı) kurma.

Makaledeki karşılığı: 'Structure Preparation and Equilibration (MM)' kutusu (Şekil 2, gri).
"""
import time
import urllib.request
from pathlib import Path

import numpy as np
_trapz = getattr(np, "trapezoid", None) or np.trapz
import openmm as mm
import openmm.app as app
import openmm.unit as u

from .config import DATA, SYSTEM, T_KELVIN, P
from . import store

DT_PS = 0.004  # integratör adımı (ps); HMR ile 4 fs


def n_steps(ps):
    """ps cinsinden süreyi adım sayısına çevirir. int(0.7/0.004) = 174 gibi kayan nokta
    kesmelerine karşı yuvarlanır."""
    return max(1, int(round(ps / DT_PS)))


# ----------------------------------------------------------------------------
# Veri indirme
# ----------------------------------------------------------------------------
def _download(url, dest, retries=5):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    for k in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
            if len(data) < 1000:
                raise IOError(f"Beklenenden küçük dosya ({len(data)} bayt): {data[:80]!r}")
            dest.write_bytes(data)
            return dest
        except Exception as e:  # GitHub bazen 429 (çok fazla istek) döndürür
            wait = 5 * (k + 1)
            print(f"  indirme denemesi {k+1} başarısız ({e}); {wait}s sonra tekrar...")
            time.sleep(wait)
    raise RuntimeError(f"İndirilemedi: {url}")


def download_inputs():
    """MobleyLab benchmarksets (CD Set 2) içinden bCD + siklopentanol (guest-3) dosyaları."""
    base = SYSTEM["source_base"]
    o = SYSTEM["orientation"]
    files = {
        "prmtop": (f"{base}/prmtop-rst7/bcd-3-{o}.prmtop", DATA / f"bcd-3-{o}.prmtop"),
        "rst7": (f"{base}/prmtop-rst7/bcd-3-{o}.rst7", DATA / f"bcd-3-{o}.rst7"),
        "guest_sdf": (f"{base}/sdf/guest-3.sdf", DATA / "guest-3.sdf"),
        "host_pdb": (f"{base}/pdb/host-bcd.pdb", DATA / "host-bcd.pdb"),
        "deney_tablosu": (f"{base}/README.md", DATA / "README_cdset2.md"),
    }
    out = {}
    for k, (url, dest) in files.items():
        out[k] = _download(url, dest)
        print(f"  {k:10s} -> {out[k].name} ({out[k].stat().st_size/1024:.0f} kB)")
    return out


# ----------------------------------------------------------------------------
# Platform seçimi
# ----------------------------------------------------------------------------
def _platform_works(plat, props):
    """Platform listede olsa da cihaz olmayabilir (ör. macOS'ta OpenCL). Küçük bir Context ile dener."""
    s = mm.System()
    s.addParticle(1.0)
    try:
        mm.Context(s, mm.VerletIntegrator(0.001), plat, props)
        return True
    except Exception:
        return False


def get_platform(verbose=True):
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    for pref in P()["platform_preference"]:
        if pref in names:
            plat = mm.Platform.getPlatformByName(pref)
            props = {}
            if pref in ("CUDA", "OpenCL"):
                props = {"Precision": "mixed"}
            if not _platform_works(plat, props):
                continue
            if verbose:
                print(f"OpenMM platformu: {pref}  (mevcut: {names})")
            return plat, props
    raise RuntimeError("Uygun OpenMM platformu yok")


# ----------------------------------------------------------------------------
# Sistem kurulumu
# ----------------------------------------------------------------------------
def build_legs(files):
    """Kompleks bacağı = dosyadaki çözülmüş bCD+misafir.
    Çözücü bacağı = aynı kutudan konakçı (bCD) silinerek elde edilir; boşluk NPT
    dengelemesinde suyla dolar. (Ek bir parametrizasyon aracı gerektirmeyen pratik bir hile.)"""
    import parmed as pmd

    s = pmd.load_file(str(files["prmtop"]), xyz=str(files["rst7"]))
    complex_prm = DATA / "complex.prmtop"
    complex_crd = DATA / "complex.rst7"
    s.save(str(complex_prm), overwrite=True)
    s.save(str(complex_crd), overwrite=True)

    solv = pmd.load_file(str(files["prmtop"]), xyz=str(files["rst7"]))
    solv.strip(f":{SYSTEM['host_resname']}")
    solv_prm = DATA / "solvent.prmtop"
    solv_crd = DATA / "solvent.rst7"
    solv.save(str(solv_prm), overwrite=True)
    solv.save(str(solv_crd), overwrite=True)

    def summary(st):
        from collections import Counter
        c = Counter(r.name for r in st.residues)
        return dict(atoms=len(st.atoms), residues=dict(c),
                    net_charge=round(sum(a.charge for a in st.atoms), 4))

    return {"complex": (complex_prm, complex_crd, summary(s)),
            "solvent": (solv_prm, solv_crd, summary(solv))}


def load_leg(leg):
    """ParmEd yapısı döner (.topology, .positions, .box_vectors, .createSystem).
    Not: Sistem ParmEd ile kuruluyor; konakçı silinmiş prmtop'u OpenMM'in kendi
    okuyucusu gereksiz yere yavaş CustomNonbondedForce ile kuruyor."""
    import parmed as pmd
    return pmd.load_file(str(DATA / f"{leg}.prmtop"), xyz=str(DATA / f"{leg}.rst7"))


def guest_atoms(topology):
    return [a.index for a in topology.atoms() if a.residue.name == SYSTEM["guest_resname"]]


def host_atoms(topology):
    return [a.index for a in topology.atoms() if a.residue.name == SYSTEM["host_resname"]]


def build_system(leg, struct):
    """Bacağın MM sistemi; kompleks bacağında konakçı-misafir kısıtı da eklenir."""
    system = make_system(struct)
    if leg == "complex":
        add_flat_bottom_restraint(system, struct.topology)
    return system


def make_context(system, positions, box, verbose=False):
    """Langevin integratörü + en uygun platformla Context kurar, kutu ve konumları yükler.
    Dengeleme, FEP ve anlık görüntü havuzu aynı integratör ayarlarını buradan alır."""
    integ = mm.LangevinMiddleIntegrator(T_KELVIN * u.kelvin, 1.0 / u.picosecond, DT_PS * u.picosecond)
    plat, props = get_platform(verbose=verbose)
    ctx = mm.Context(system, integ, plat, props)
    ctx.setPeriodicBoxVectors(*box)
    ctx.setPositions(positions)
    return ctx, integ


def make_system(struct, hmr=True):
    """PME, 9 Å kesme, H-bağ kısıtları, hidrojen kütle yeniden dağıtımı (4 fs adım)."""
    return struct.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=0.9 * u.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        hydrogenMass=(1.5 * u.amu if hmr else None),
    )


# Düz tabanlı (flat-bottom) konakçı-misafir kütle merkezi kısıtı
RESTRAINT_R0_NM = 0.35
RESTRAINT_K = 4184.0  # kJ/mol/nm^2  (= 10 kcal/mol/Å^2)


def add_flat_bottom_restraint(system, topology):
    """Misafir, konakçının boşluğundan uzaklaşırsa geri çeken kuvvet.
    Misafir tamamen 'kapalı' (decoupled) olduğunda kaybolup gitmesin diye gerekir;
    standart durum düzeltmesiyle analitik olarak geri alınır."""
    host = [i for i in host_atoms(topology)
            if system.getParticleMass(i).value_in_unit(u.amu) > 2.0]
    guest = [i for i in guest_atoms(topology)
             if system.getParticleMass(i).value_in_unit(u.amu) > 2.0]
    f = mm.CustomCentroidBondForce(2, "step(d-r0)*0.5*k*(d-r0)^2; d=distance(g1,g2)")
    f.addGlobalParameter("r0", RESTRAINT_R0_NM)
    f.addGlobalParameter("k", RESTRAINT_K)
    f.addGroup(host)
    f.addGroup(guest)
    f.addBond([0, 1], [])
    f.setUsesPeriodicBoundaryConditions(True)
    f.setForceGroup(5)
    system.addForce(f)
    return f


def restraint_standard_state_correction():
    """kT ln(V0 / V_site), V_site = ∫ 4πr² exp(-βU(r)) dr  (kJ/mol).
    ΔG_bağlanma'ya EKLENİR (bkz. fep.combine)."""
    from .config import KT, V0_NM3
    r = np.linspace(0, RESTRAINT_R0_NM + 2.0, 200001)
    U = np.where(r > RESTRAINT_R0_NM, 0.5 * RESTRAINT_K * (r - RESTRAINT_R0_NM) ** 2, 0.0)
    v_site = _trapz(4 * np.pi * r ** 2 * np.exp(-U / KT), r)
    return KT * np.log(V0_NM3 / v_site), v_site


def equilibrate(leg, ps=None, report=True):
    """Kısa enerji minimizasyonu + NPT dengeleme. Dengelenmiş koordinatları
    results/<ön_ayar>/ altına kaydeder (her ön ayarın kendi dengelemesi olur)."""
    ps = ps if ps is not None else P()["equil_npt_ps"]
    st0 = load_leg(leg)
    system = build_system(leg, st0)
    system.addForce(mm.MonteCarloBarostat(1.0 * u.bar, T_KELVIN * u.kelvin, 25))
    ctx, integ = make_context(system, st0.positions, st0.box_vectors, verbose=report)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 500)
    ctx.setVelocitiesToTemperature(T_KELVIN * u.kelvin)
    nsteps = n_steps(ps)
    vols, t0 = [], time.time()
    chunk = max(1, nsteps // 20)
    done = 0
    while done < nsteps:
        k = min(chunk, nsteps - done)
        integ.step(k)
        done += k
        st = ctx.getState()
        vols.append(st.getPeriodicBoxVolume().value_in_unit(u.nanometer ** 3))
    wall = time.time() - t0
    st = ctx.getState(getPositions=True)
    np.savez(_equil_path(leg),
             positions=st.getPositions(asNumpy=True).value_in_unit(u.nanometer),
             box=st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer))
    ns_day = (ps / 1000) / (wall / 86400) if wall > 0 else float("nan")
    if report:
        print(f"{leg}: {ps} ps NPT, {wall:.1f} s  (~{ns_day:.0f} ns/gün)")
    return dict(volumes_nm3=vols, wall_s=wall, ns_per_day=ns_day)


def _equil_path(leg):
    return store.path(f"{leg}_equil.npz")


def load_equilibrated(leg):
    p = _equil_path(leg)
    if not p.exists():
        raise FileNotFoundError(f"{p} yok — bu ön ayar için 01. notebook'u (dengeleme) çalıştırdınız mı?")
    d = np.load(p)
    return d["positions"], d["box"]
