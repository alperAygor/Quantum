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
def get_platform(verbose=True):
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    for pref in P()["platform_preference"]:
        if pref in names:
            plat = mm.Platform.getPlatformByName(pref)
            props = {}
            if pref in ("CUDA", "OpenCL"):
                props = {"Precision": "mixed"}
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
    """Kısa enerji minimizasyonu + NPT dengeleme. Dengelenmiş koordinatları kaydeder."""
    ps = ps if ps is not None else P()["equil_npt_ps"]
    st0 = load_leg(leg)
    system = make_system(st0)
    if leg == "complex":
        add_flat_bottom_restraint(system, st0.topology)
    system.addForce(mm.MonteCarloBarostat(1.0 * u.bar, T_KELVIN * u.kelvin, 25))
    integ = mm.LangevinMiddleIntegrator(T_KELVIN * u.kelvin, 1.0 / u.picosecond, 0.004 * u.picosecond)
    plat, props = get_platform(verbose=report)
    sim = app.Simulation(st0.topology, system, integ, plat, props)
    sim.context.setPeriodicBoxVectors(*st0.box_vectors)
    sim.context.setPositions(st0.positions)
    sim.minimizeEnergy(maxIterations=500)
    sim.context.setVelocitiesToTemperature(T_KELVIN * u.kelvin)
    nsteps = int(ps / 0.004)
    vols, t0 = [], time.time()
    chunk = max(1, nsteps // 20)
    for _ in range(0, nsteps, chunk):
        sim.step(chunk)
        st = sim.context.getState()
        vols.append(st.getPeriodicBoxVolume().value_in_unit(u.nanometer ** 3))
    wall = time.time() - t0
    st = sim.context.getState(getPositions=True)
    np.savez(DATA / f"{leg}_equil.npz",
             positions=st.getPositions(asNumpy=True).value_in_unit(u.nanometer),
             box=st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer))
    ns_day = (ps / 1000) / (wall / 86400) if wall > 0 else float("nan")
    if report:
        print(f"{leg}: {ps} ps NPT, {wall:.1f} s  (~{ns_day:.0f} ns/gün)")
    return dict(volumes_nm3=vols, wall_s=wall, ns_per_day=ns_day)


def load_equilibrated(leg):
    d = np.load(DATA / f"{leg}_equil.npz")
    return d["positions"], d["box"]
