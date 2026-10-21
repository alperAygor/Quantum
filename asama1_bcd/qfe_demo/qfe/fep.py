"""Aşama 2: Alkimik serbest enerji pertürbasyonu (FEP) + MBAR.

Makaledeki karşılığı: 'Free Energy Perturbation (MM)' (Şekil 2, kırmızı üst kutu), Denklem 3-5.

Çift ayrıştırma (double decoupling):
  kompleks bacağı : misafirin çevreyle etkileşimini kapat (konakçı içinde, kısıtlı)  -> ΔG_c
  çözücü bacağı   : misafirin çevreyle etkileşimini kapat (saf suda)                -> ΔG_s
  ΔG_bağlanma = ΔG_s - ΔG_c + kT ln(V0/V_site)
"""
import time

import numpy as np
import openmm as mm
import openmm.unit as u

from .config import KT, T_KELVIN, P
from . import mmsys


def lambda_schedule():
    """Önce elektrostatik (λ_e: 1→0), sonra van der Waals (λ_s: 1→0) kapatılır."""
    p = P()
    states = [(le, 1.0) for le in p["lambda_elec"]]
    states += [(0.0, ls) for ls in p["lambda_ster"][1:]]
    return states


def make_alchemical(leg):
    from openmmtools import alchemy

    st = mmsys.load_leg(leg)
    system = mmsys.make_system(st)
    if leg == "complex":
        mmsys.add_flat_bottom_restraint(system, st.topology)
    guest = mmsys.guest_atoms(st.topology)
    factory = alchemy.AbsoluteAlchemicalFactory(disable_alchemical_dispersion_correction=True)
    region = alchemy.AlchemicalRegion(alchemical_atoms=guest,
                                      annihilate_electrostatics=True,
                                      annihilate_sterics=False)
    asys = factory.create_alchemical_system(system, region)
    return st, asys


def _set_lambda(ctx, le, ls):
    ctx.setParameter("lambda_electrostatics", le)
    ctx.setParameter("lambda_sterics", ls)


def run_leg(leg, progress=True):
    """Her λ durumunda MD yapar; her örnekte TÜM durumlardaki indirgenmiş enerjiyi
    (u = βU) hesaplar -> MBAR için u_kln matrisi."""
    p = P()
    states = lambda_schedule()
    K = len(states)
    _, asys = make_alchemical(leg)
    pos, box = mmsys.load_equilibrated(leg)
    beta = 1.0 / KT
    integ = mm.LangevinMiddleIntegrator(T_KELVIN * u.kelvin, 1.0 / u.picosecond, 0.004 * u.picosecond)
    plat, props = mmsys.get_platform(verbose=progress)
    ctx = mm.Context(asys, integ, plat, props)
    ctx.setPeriodicBoxVectors(*box)
    ctx.setPositions(pos)

    steps_eq = int(p["fep_equil_ps"] / 0.004)
    steps_sample = max(1, int(p["fep_sample_ps"] / 0.004))
    n_samples = int(p["fep_prod_ps"] / p["fep_sample_ps"])
    u_kln = np.zeros((K, K, n_samples))
    t0 = time.time()
    for k, (le, ls) in enumerate(states):
        _set_lambda(ctx, le, ls)
        if k == 0:
            mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 200)
            ctx.setVelocitiesToTemperature(T_KELVIN * u.kelvin)
        integ.step(steps_eq)
        for n in range(n_samples):
            integ.step(steps_sample)
            for l, (le2, ls2) in enumerate(states):
                _set_lambda(ctx, le2, ls2)
                e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                u_kln[k, l, n] = beta * e
            _set_lambda(ctx, le, ls)
        if progress:
            print(f"  [{leg}] durum {k+1:2d}/{K}  λ_e={le:.2f} λ_s={ls:.2f}   "
                  f"geçen süre {time.time()-t0:6.0f} s")
    return dict(states=np.array(states), u_kln=u_kln, wall_s=time.time() - t0)


def mbar_analysis(u_kln):
    """MBAR (makalede kullanılan çok durumlu Bennett kabul oranı) ile serbest enerjiler."""
    from pymbar import MBAR, timeseries

    K, _, N = u_kln.shape
    # Her durum için korelasyonsuz örnekleri seç
    idx_list = []
    for k in range(K):
        try:
            g = timeseries.statistical_inefficiency(u_kln[k, k, :])
            idx = timeseries.subsample_correlated_data(u_kln[k, k, :], g=g)
        except Exception:
            idx = np.arange(N)
        if len(idx) < 3:
            idx = np.arange(N)
        idx_list.append(np.array(idx))
    N_k = np.array([len(i) for i in idx_list])
    u_kn = np.concatenate([u_kln[k][:, idx_list[k]] for k in range(K)], axis=1)
    mbar = MBAR(u_kn, N_k, solver_protocol="robust")
    res = mbar.compute_free_energy_differences()
    f = res["Delta_f"][0] * KT
    df = res["dDelta_f"][0] * KT
    try:
        O = mbar.compute_overlap()["matrix"]
    except Exception:
        O = None
    return dict(f_kJ=f, df_kJ=df, N_k=N_k, overlap=O)


def combine(res_complex, res_solvent):
    """ΔG_bağlanma = ΔG_s - ΔG_c + kT ln(V0/V_site).
    ΔG_x = f(kapalı) - f(açık) (MBAR dizisinin son elemanı)."""
    dGc, ddGc = res_complex["f_kJ"][-1], res_complex["df_kJ"][-1]
    dGs, ddGs = res_solvent["f_kJ"][-1], res_solvent["df_kJ"][-1]
    corr, v_site = mmsys.restraint_standard_state_correction()
    dG = dGs - dGc + corr
    err = float(np.sqrt(ddGc ** 2 + ddGs ** 2))
    return dict(dG_complex_decouple=float(dGc), dG_solvent_decouple=float(dGs),
                standard_state_corr=float(corr), v_site_nm3=float(v_site),
                dG_bind=float(dG), dG_bind_err=err)


def sample_endstate_pool(leg, n_frames=None, stride_ps=None, progress=True):
    """Tam etkileşimli (λ=1) durumdan konfigürasyon havuzu.
    Bunlar QM/MM hesapları ve ML için 'anlık görüntüler' (snapshots) olur.
    Makalede ~2000 anlık görüntü (FEP yörüngesinden) alınır."""
    p = P()
    n_frames = n_frames or p["pool_frames"]
    stride_ps = stride_ps or p["pool_stride_ps"]
    st = mmsys.load_leg(leg)
    system = mmsys.make_system(st)
    if leg == "complex":
        mmsys.add_flat_bottom_restraint(system, st.topology)
    pos, box = mmsys.load_equilibrated(leg)
    integ = mm.LangevinMiddleIntegrator(T_KELVIN * u.kelvin, 1.0 / u.picosecond, 0.004 * u.picosecond)
    plat, props = mmsys.get_platform(verbose=False)
    ctx = mm.Context(system, integ, plat, props)
    ctx.setPeriodicBoxVectors(*box)
    ctx.setPositions(pos)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 200)
    ctx.setVelocitiesToTemperature(T_KELVIN * u.kelvin)
    integ.step(int(10 / 0.004))
    steps = max(1, int(stride_ps / 0.004))
    X = np.zeros((n_frames, system.getNumParticles(), 3), dtype=np.float32)
    B = np.zeros((n_frames, 3), dtype=np.float32)
    t0 = time.time()
    for n in range(n_frames):
        integ.step(steps)
        st = ctx.getState(getPositions=True)   # periyodik sarma YOK: moleküller bütün kalır
        X[n] = st.getPositions(asNumpy=True).value_in_unit(u.nanometer)
        bv = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer)
        B[n] = np.diag(bv)
    if progress:
        print(f"  [{leg}] {n_frames} anlık görüntü, {time.time()-t0:.0f} s")
    return X, B
