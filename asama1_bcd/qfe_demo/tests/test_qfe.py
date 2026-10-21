"""qfe paketinin hızlı birim testleri (birkaç saniye; MD veya uzun QM hesabı yok).

Çalıştırma (qfe_demo klasöründe):  python -m pytest tests
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qfe import config  # noqa: E402

KT = config.KT


# ----------------------------------------------------------------------
# Serbest enerji yardımcıları
# ----------------------------------------------------------------------
def test_zwanzig_gaussian_matches_analytic():
    """Gauss dağılımlı ΔU için ΔA = μ − σ²/(2kT) analitik olarak doğrudur."""
    from qfe import ml

    rng = np.random.default_rng(0)
    mu, sigma = 3.0, 1.5
    z = ml.zwanzig(rng.normal(mu, sigma, 20000), n_boot=20)
    exact = mu - sigma ** 2 / (2 * KT)
    assert z["cumulant"] == pytest.approx(exact, abs=0.05)
    assert z["exp"] == pytest.approx(exact, abs=0.1)


def test_binding_correction_ignores_common_offset():
    """QM ve MM'nin farklı sıfır noktası (iki bacakta aynı sabit) sonucu değiştirmemeli."""
    from qfe import ml

    rng = np.random.default_rng(1)
    c, s = rng.normal(-2.0, 1.0, 500), rng.normal(0.5, 1.0, 500)
    a = ml.binding_correction(c, s)
    b = ml.binding_correction(c + 1234.5, s + 1234.5)
    assert a["exp"] == pytest.approx(b["exp"], abs=1e-9)
    assert a["cumulant"] == pytest.approx(b["cumulant"], abs=1e-9)


def test_restraint_standard_state_correction():
    """Düz tabanlı kısıtın hacmi yaklaşık kesme küresinin hacmidir; düzeltme pozitif olmalı
    (misafiri boşlukta tutmak serbest enerjiye ceza ekler)."""
    pytest.importorskip("openmm")
    from qfe import mmsys

    corr, v_site = mmsys.restraint_standard_state_correction()
    r0 = mmsys.RESTRAINT_R0_NM
    width = np.sqrt(KT / mmsys.RESTRAINT_K)
    assert 4 / 3 * np.pi * r0 ** 3 < v_site < 4 / 3 * np.pi * (r0 + 3 * width) ** 3
    assert corr == pytest.approx(KT * np.log(config.V0_NM3 / v_site))
    assert corr > 0


def test_n_steps_rounds_instead_of_truncating():
    pytest.importorskip("openmm")
    from qfe import mmsys

    assert int(0.7 / mmsys.DT_PS) == 174  # eski davranış: bir adım eksik
    assert mmsys.n_steps(0.7) == 175
    assert mmsys.n_steps(0.0) == 1


# ----------------------------------------------------------------------
# Kuantum algoritmaları (LiH, CAS(2e, 3o) → 6 qubit)
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def lih():
    pytest.importorskip("pyscf")
    pytest.importorskip("openfermion")
    from pyscf import gto, mcscf, scf
    from qfe import quantum

    mol = gto.M(atom="Li 0 0 0; H 0 0 1.6", basis="sto-3g", verbose=0)
    mf = scf.RHF(mol).run()
    ncas, nel = 3, 2
    h1, h2, ecore = quantum.active_space_integrals(mf, ncas, nel)
    qop = quantum.qubit_hamiltonian(h1, h2, ecore)
    mc = mcscf.CASCI(mf, ncas, nel)
    mc.verbose = 0
    return dict(qop=qop, nq=2 * ncas, nel=nel, e_casci=mc.kernel()[0])


def test_qubit_hamiltonian_matches_casci(lih):
    from qfe import quantum

    HN, idx, hf_pos, _ = quantum.number_sector(lih["qop"], lih["nq"], lih["nel"])
    w, v = quantum.exact_ground(HN)
    assert w[0] == pytest.approx(lih["e_casci"], abs=1e-8)
    assert abs(v[hf_pos, 0]) ** 2 > 0.9  # LiH'de HF iyi bir başlangıç durumu


def test_vqe_reaches_exact_energy(lih):
    from qfe import quantum

    HN, *_ = quantum.number_sector(lih["qop"], lih["nq"], lih["nel"])
    e0 = np.linalg.eigvalsh(HN)[0]
    res = quantum.vqe_uccsd(lih["qop"], lih["nq"], lih["nel"])
    assert res["energy"] == pytest.approx(e0, abs=1e-6)


def test_qpe_peak_at_ground_energy(lih):
    from qfe import quantum

    HN, *_ = quantum.number_sector(lih["qop"], lih["nq"], lih["nel"])
    w, v = np.linalg.eigh(HN)
    q = quantum.qpe_simulate(HN, v[:, 0], n_bits=8)
    e_peak = q["energies"][np.argmax(q["prob"])]
    assert abs(e_peak - w[0]) <= q["resolution"]


def test_trotter_error_decreases_with_steps(lih):
    from qfe import quantum

    HN, idx, hf_pos, Hs = quantum.number_sector(lih["qop"], lih["nq"], lih["nel"])
    w, v = np.linalg.eigh(HN)
    psi0 = np.zeros(Hs.shape[0], complex)
    psi0[idx] = v[:, 0]
    out = quantum.trotter_phase_error(lih["qop"], lih["nq"], psi0, w[0], t=1.0, steps_list=[1, 4, 16])
    errs = [o["error"] for o in out]
    assert errs[0] > errs[1] > errs[2]


def test_pauli_list_is_little_endian():
    pytest.importorskip("openfermion")
    from openfermion import QubitOperator
    from qfe import quantum

    pl = quantum.pauli_list(QubitOperator("Z0 X2", 0.5), 3)
    assert pl == [("XIZ", 0.5 + 0j)]  # Qiskit: qubit 0 en sağdaki karakter


# ----------------------------------------------------------------------
# QM etiket önbelleği
# ----------------------------------------------------------------------
def test_label_cache_rejected_for_other_pool_or_settings():
    """Havuz ya da QM ayarı değişince eski etiketler yeniden kullanılmamalı."""
    from qfe import qmmm

    config.set_preset("test")
    rng = np.random.default_rng(0)
    pools = {"complex": (rng.random((3, 5, 3), dtype=np.float32), np.ones((3, 3), np.float32))}
    lab = qmmm.Labeler(regions={}, pools=pools)
    lab.cache_dft[("complex", 0)] = dict(emb=1.0, vac=0.0, mm_int=0.0, mm_coul=0.0)
    saved = lab.to_dict()

    same = qmmm.Labeler(regions={}, pools=pools)
    same.load_dict(saved)
    assert ("complex", 0) in same.cache_dft

    moved = {"complex": (pools["complex"][0] + 0.01, pools["complex"][1])}
    other = qmmm.Labeler(regions={}, pools=moved)
    other.load_dict(saved)
    assert not other.cache_dft

    old = dict(config.PRESETS["test"])
    try:
        config.PRESETS["test"]["qm_cutoff_A"] = 10.0
        cut = qmmm.Labeler(regions={}, pools=pools)
        cut.load_dict(saved)
        assert not cut.cache_dft
    finally:
        config.PRESETS["test"].clear()
        config.PRESETS["test"].update(old)
