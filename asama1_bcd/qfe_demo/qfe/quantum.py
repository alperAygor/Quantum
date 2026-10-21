"""Aşama 6-7: Kuantum çekirdeğini kuantum algoritmalarıyla (simülatörde) çözmek ve
gerçek bir kuantum bilgisayar için kaynak tahmini.

Makaledeki karşılığı: Bölüm 3 ve Şekil 4
  - aktif uzay Hamiltoniyeni -> qubit Hamiltoniyeni (Jordan-Wigner)
  - kuantum faz tahmini (QPE), Trotter hatası, başlangıç durumu örtüşmesi (Şekil 4d)
  - qubit / kapı sayısı tahminleri (Şekil 4a-c)

Burada HER ŞEY klasik bilgisayarda durum vektörü (statevector) simülasyonudur.
Bu yöntem ~25-30 qubit'e kadar mümkündür; makaledeki ilgi çekici aktif uzaylar
(30-60 orbital = 60-120 qubit) için gerçek, hata düzeltmeli kuantum bilgisayar gerekir.
"""
import numpy as np
from scipy.sparse.linalg import expm_multiply, eigsh


# ----------------------------------------------------------------------
# Aktif uzay Hamiltoniyeni
# ----------------------------------------------------------------------
def active_space_integrals(mf, ncas, nelecas):
    """HF orbitallerinden aktif uzay integralleri: h1 (ncas×ncas), h2 (kimyacı (pq|rs)), E_çekirdek.
    E_çekirdek = nükleer itme + donmuş (aktif olmayan) elektronların enerjisi + nokta yük etkisi."""
    from pyscf import mcscf, ao2mo

    mc = mcscf.CASCI(mf, ncas, nelecas)
    h1, ecore = mc.get_h1eff()
    h2 = ao2mo.restore(1, mc.get_h2eff(), ncas)
    return h1, h2, ecore


def fermion_hamiltonian(h1, h2, ecore):
    from openfermion import InteractionOperator
    from openfermion.chem.molecular_data import spinorb_from_spatial

    one_so, two_so = spinorb_from_spatial(h1, np.asarray(h2.transpose(0, 2, 3, 1), order="C"))
    return InteractionOperator(ecore, one_so, 0.5 * two_so)


def qubit_hamiltonian(h1, h2, ecore):
    from openfermion import jordan_wigner

    return jordan_wigner(fermion_hamiltonian(h1, h2, ecore))


def pauli_list(qop, n_qubits):
    """Qiskit uyumlu (etiket, katsayı) listesi. Qiskit'te kullanım:
        from qiskit.quantum_info import SparsePauliOp
        H = SparsePauliOp.from_list(pauli_list(qop, n))
    Qiskit küçük-uçlu (little-endian) sıralama kullanır: qubit 0 en sağdaki karakterdir."""
    out = []
    for term, c in qop.terms.items():
        lab = ["I"] * n_qubits
        for q, p in term:
            lab[n_qubits - 1 - q] = p
        out.append(("".join(lab), complex(c)))
    return out


# ----------------------------------------------------------------------
# Tam çözüm (klasik referans)
# ----------------------------------------------------------------------
def number_sector(qop, n_qubits, n_elec):
    """Hamiltoniyeni N-elektron alt uzayına indirger (yoğun matris) + HF durumu indeksi."""
    from openfermion import get_sparse_operator, jw_number_indices

    Hs = get_sparse_operator(qop, n_qubits=n_qubits).tocsr()
    idx = np.array(jw_number_indices(n_elec, n_qubits))
    HN = Hs[idx][:, idx].toarray()
    hf_full = sum(1 << (n_qubits - 1 - q) for q in range(n_elec))  # OpenFermion: qubit 0 = en anlamlı bit
    hf_pos = int(np.where(idx == hf_full)[0][0])
    return HN, idx, hf_pos, Hs


def exact_ground(HN):
    w, v = np.linalg.eigh(HN)
    return w, v


# ----------------------------------------------------------------------
# VQE (UCCSD), durum vektörü simülasyonu
# ----------------------------------------------------------------------
def vqe_uccsd(qop, n_qubits, n_elec, maxiter=200, shots=None, seed=0):
    """Değişkenli kuantum özçözücü. shots verilirse her enerji ölçümüne
    istatistiksel (atış) gürültüsü eklenir (gerçek donanıma yakın)."""
    from openfermion import (get_sparse_operator, jordan_wigner, uccsd_singlet_generator,
                             uccsd_singlet_paramsize)
    from scipy.optimize import minimize

    rng = np.random.default_rng(seed)
    H = get_sparse_operator(qop, n_qubits=n_qubits).tocsr()
    hf = np.zeros(2 ** n_qubits, complex)
    hf[sum(1 << (n_qubits - 1 - q) for q in range(n_elec))] = 1.0
    npar = uccsd_singlet_paramsize(n_qubits, n_elec)
    gen_cache = {}

    def state(theta):
        G = get_sparse_operator(jordan_wigner(uccsd_singlet_generator(theta, n_qubits, n_elec)),
                                n_qubits=n_qubits)
        return expm_multiply(G, hf)

    # Var(H) ~ atış gürültüsü ölçeği için
    trace = []

    def energy(theta):
        psi = state(theta)
        e = float(np.real(np.vdot(psi, H @ psi)))
        if shots:
            var = float(np.real(np.vdot(psi, H @ (H @ psi)))) - e ** 2
            e = e + rng.normal(0, np.sqrt(max(var, 1e-12) / shots))
        trace.append(e)
        return e

    x0 = np.full(npar, 1e-3)   # tam sıfırda sayısal gradyan sorunlu; küçük başlangıç
    e_ref = float(np.real(np.vdot(hf, H @ hf)))
    # Korelasyon enerjisi küçük (mHa mertebesi) olduğu için optimizasyonu mHa biriminde yapıyoruz
    obj = lambda th: (energy(th) - e_ref) * 1e3
    if shots:
        res = minimize(obj, x0, method="COBYLA", options={"maxiter": maxiter, "rhobeg": 0.05})
    else:
        res = minimize(obj, x0, method="BFGS", options={"maxiter": maxiter, "gtol": 1e-7})
    psi = state(res.x)
    e_final = float(np.real(np.vdot(psi, H @ psi)))
    return dict(energy=e_final, n_params=npar, trace=np.array(trace), theta=res.x, n_evals=len(trace))


# ----------------------------------------------------------------------
# Kuantum faz tahmini (QPE) — devre seviyesinde durum vektörü simülasyonu
# ----------------------------------------------------------------------
def qpe_simulate(HN, psi_in, n_bits, e_window=None):
    """Standart QPE devresi: m yardımcı (ancilla) qubit üzerinde Hadamard,
    kontrollü U^(2^k) kapıları, ters kuantum Fourier dönüşümü, ölçüm.
    U = exp(-i (H - E_kaydırma) t).  Dönen: her ölçüm sonucunun olasılığı ve karşılık gelen enerji.
    """
    from scipy.linalg import expm

    w = np.linalg.eigvalsh(HN)
    emin, emax = (e_window if e_window else (w.min(), w.max()))
    width = (emax - emin) * 1.05 + 1e-9
    e_shift = emin - 0.025 * (emax - emin)
    t = 2 * np.pi / width                     # tüm spektrum [0, 2π) faz aralığına sığar (örtüşme yok)
    U = expm(-1j * (HN - e_shift * np.eye(len(HN))) * t)
    M = 2 ** n_bits
    A = np.empty((M, len(psi_in)), complex)
    v = psi_in.astype(complex)
    for j in range(M):                        # kontrollü-U^j uygulanmış dal
        A[j] = v
        v = U @ v
    A /= np.sqrt(M)
    # ters QFT (yardımcı yazmaç üzerinde): e^{+iφ2πj} faz -> k = φ M
    Aout = np.fft.ifft(A, axis=0) * np.sqrt(M)
    # U|E> = e^{-i(E-Es)t}|E>  =>  tepe k = (E-Es) t M / 2π  =>  E = Es + (k/M)·2π/t
    prob = np.sum(np.abs(Aout) ** 2, axis=1)
    phi = np.arange(M) / M
    energies = e_shift + phi * 2 * np.pi / t
    return dict(prob=prob / prob.sum(), energies=energies, t=t, e_shift=e_shift,
                resolution=2 * np.pi / t / M)


def qpe_shots(qpe, n_shots, seed=0):
    """Sonlu sayıda ölçüm. Tahmin: en sık gözlenen sonuç ve iki komşusunun
    sayımlarla ağırlıklı ortalaması (Fejér çekirdeğinin yan loblarını dengeler)."""
    rng = np.random.default_rng(seed)
    M = len(qpe["prob"])
    k = rng.choice(M, size=n_shots, p=qpe["prob"])
    counts = np.bincount(k, minlength=M)
    m = int(np.argmax(counts))
    nb = [j for j in (m - 1, m, m + 1) if 0 <= j < M]
    w = counts[nb].astype(float)
    e = float(np.sum(w * qpe["energies"][nb]) / w.sum())
    return e, counts


# ----------------------------------------------------------------------
# Trotter hatası
# ----------------------------------------------------------------------
def trotter_phase_error(qop, n_qubits, psi0, e0, t, steps_list, order=1):
    """Gerçek kuantum donanımında exp(-iHt) tam uygulanamaz; Pauli terimlerinin
    üstellerinin çarpımına bölünür (Trotter). Hata adım sayısıyla azalır."""
    from openfermion import get_sparse_operator, QubitOperator

    terms = [(QubitOperator(term, 1.0), float(np.real(c))) for term, c in qop.terms.items() if term]
    const = float(np.real(qop.terms.get((), 0.0)))
    mats = [(get_sparse_operator(op, n_qubits=n_qubits).tocsr(), c) for op, c in terms]
    out = []
    for r in steps_list:
        dt = t / r
        psi = psi0.astype(complex).copy()
        seq = mats if order == 1 else mats + mats[::-1]
        fac = 1.0 if order == 1 else 0.5
        for _ in range(r):
            for P, c in seq:
                a = c * dt * fac
                psi = np.cos(a) * psi - 1j * np.sin(a) * (P @ psi)
        amp = np.vdot(psi0, psi) * np.exp(-1j * const * t)
        e_eff = e0 - np.angle(amp * np.exp(1j * e0 * t)) / t
        out.append(dict(steps=r, energy=float(e_eff), error=float(abs(e_eff - e0)),
                        n_exponentials=r * len(seq)))
    return out


# ----------------------------------------------------------------------
# Başlangıç durumu örtüşmesi ve kaynak tahmini
# ----------------------------------------------------------------------
def hf_overlap_scan(mf, sizes):
    """CAS(n,n) için |<HF|Ψ0>|² — makale Şekil 4d'nin küçük ölçekli karşılığı.
    QPE'nin başarı olasılığı bu örtüşmeyle orantılıdır."""
    import time
    from pyscf import mcscf

    out = []
    for n in sizes:
        t0 = time.time()
        mc = mcscf.CASCI(mf, n, n)
        mc.verbose = 0
        mc.fcisolver.max_cycle = 200
        mc.kernel()
        ci = np.asarray(mc.ci)
        out.append(dict(n_orb=n, n_qubits=2 * n, overlap_hf=float(ci[0, 0] ** 2),
                        n_determinants=int(ci.size), seconds=time.time() - t0,
                        e_casci=float(mc.e_tot)))
    return out


def resource_estimates(mf, sizes, eps_ha=1.6e-3):
    """Aktif uzay büyüdükçe qubit Hamiltoniyeninin terim sayısı L, 1-normu λ ve
    kaba kaynak tahminleri. (Makale Şekil 4: qDRIFT, kısmen rastgele Trotter, qubitization.)

    KABA (mertebe) formüller — sadece ölçeklenmeyi göstermek için:
      qubitization: kuantum yürüyüş adımı ≈ π λ / (2 ε);  Toffoli ≈ adım × L (naif SELECT)
      qDRIFT      : üstel sayısı ≈ (π λ / ε)²  ; 1 yardımcı qubit
    Makaledeki 'double factorization' gibi gelişmiş yöntemler bu sayıları ciddi azaltır.
    """
    import time

    out = []
    for n in sizes:
        t0 = time.time()
        h1, h2, ecore = active_space_integrals(mf, n, n)
        qop = qubit_hamiltonian(h1, h2, ecore)
        coeffs = np.array([abs(c) for term, c in qop.terms.items() if term])
        L, lam = len(coeffs), float(coeffs.sum())
        walk = np.pi * lam / (2 * eps_ha)
        out.append(dict(n_orb=n, n_qubits=2 * n, n_pauli_terms=L, lambda_ha=lam,
                        qubitization_walk_steps=walk,
                        qubitization_toffoli_naive=walk * L,
                        qubitization_logical_qubits=2 * n + int(np.ceil(np.log2(L))) + 10,
                        qdrift_exponentials=(np.pi * lam / eps_ha) ** 2,
                        qdrift_logical_qubits=2 * n + 1,
                        seconds=time.time() - t0))
    return out
