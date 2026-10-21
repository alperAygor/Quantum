# Kompakt VQE degerlendirmesi: ayni CAS(4e,4o) problemi, 8 yerine 5 kubit, ~4 kat daha sig devre.
#
#   python3 ibm_vqe_kompakt.py hazirla porfin            # ucretsiz, yerel: devreler + gurultulu simulasyon
#   python3 ibm_vqe_kompakt.py gonder  porfin --sahte    # ayni akis FakeKingston uzerinde (kredi harcamaz)
#   python3 ibm_vqe_kompakt.py gonder  porfin            # IBM donanimina gonderir (onay sorar)
#   python3 ibm_vqe_kompakt.py sonuc   porfin [--sahte]
#   secenekler: --backend=ibm_kingston  --shots=4096
#
# ibm_vqe.py'den farki (ayni "parametreler klasikte, enerji donanimda olculur" yaklasimi):
#  * Elektron sayisi (2 alfa + 2 beta) ve singlet spin simetrisi sabit oldugu icin 8 kubitin 256 durumundan
#    yalnizca 21'i fiziksel: |a,a> ve (|a,b>+|b,a>)/sqrt2. Bu 21 durum 5 kubite (32 durum) yerlestirilir.
#    CASCI S0 ve S1 bu alt uzayda tam olarak temsil edilir (fark < 1e-10 eV).
#  * Durum hazirlama: gercek genlikler icin Mottonen agaci = 30 CNOT (ibm_vqe: S0 264, S1 170 iki-kubit kapi
#    donanim uzerinde; burada donanimda ~60).
#  * 21 durumun 32'lik tabana yerlesimi, |c| >= 1e-3 kesmesinde enerji kaymasi < 0.005 eV olacak ve
#    qubit-wise olcum grubu sayisi en az olacak sekilde aranir.
#  * S0 ve S1 devreleri ayni CNOT iskeletine sahip; acilar 0'a baglaninca ikisi de |00000> = HF hazirlar.
#    Bu yuzden tek bir ortak gurultu-referans devresi yeter: 4 yerine 3 devre.
#  * Gurultu duzeltmesi ibm_vqe.py ile ayni: E_duz = c0 + (E_olc - c0) * (E_ref_ideal - c0) / (E_ref_olc - c0)
import os, sys, json, time, pickle
os.environ.setdefault("PYSCF_TMPDIR", os.path.abspath("pyscf_tmp")); os.makedirs(os.environ["PYSCF_TMPDIR"], exist_ok=True)
import numpy as np
from pyscf import gto, scf, dft, mcscf, fci, ao2mo  # qiskit oncesi yuklenmeli (libgomp TLS hatasi)
from pyscf.fci import cistring

EV = 27.211386
KMD, NAME = sys.argv[1], sys.argv[2]
SAHTE = "--sahte" in sys.argv
ARG = {a.split("=")[0]: a.split("=")[1] for a in sys.argv[3:] if "=" in a}
SHOTS = int(ARG.get("--shots", 4096)); BACKEND = ARG.get("--backend")
ESIK, NQ, D = 1e-3, 5, 32
PKL = f"{NAME}_kompakt_hazir.pkl"; JOBF = f"{NAME}_kompakt_is{'_sahte' if SAHTE else ''}.json"
# Karsilastirma icin onceki donanim sonucu (ibm_vqe.py, docs/00_Master_Teknik_Rapor.md 5.3)
ONCEKI = {"porfin": dict(donanim="ibm_kingston", duzeltilmis_fark_eV=0.93, hata_eV=0.19)}

def log(s): print(s, flush=True)

# ---------------------------------------------------------------- Hamiltonyen ve 21 durumluk alt uzay
def hamiltonyen():
    lines = [l for l in open(f"{NAME}.xyz").read().split("\n")[2:] if l.strip()]
    ecp = "def2-svp" if any(l.split()[0] == "Pd" for l in lines) else None
    mol = gto.M(atom=f"{NAME}.xyz", basis="def2-svp", ecp=ecp, verbose=0, max_memory=int(os.environ.get("QMEM", "8000")))
    ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
    if not ks.converged: sys.exit("B3LYP yakinsamadi; Hamiltonyen guvenilir degil")
    mf = scf.RHF(mol).density_fit(); mo = ks.mo_coeff
    cS = mcscf.CASCI(mf, 4, 4); cS.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0)
    cS.fcisolver.nroots = 3; cS.kernel(mo)
    h1, ec = cS.get_h1eff(mo); eri = ao2mo.restore(1, cS.get_h2eff(mo), 4)
    na = cistring.num_strings(4, 2)
    h2e = fci.direct_spin1.absorb_h1e(h1, eri, 4, (2, 2), 0.5)
    Hdet = np.zeros((na * na, na * na))
    for k in range(na * na):
        e = np.zeros(na * na); e[k] = 1
        Hdet[:, k] = fci.direct_spin1.contract_2e(h2e, e.reshape(na, na), 4, (2, 2)).ravel()
    Hdet += ec * np.eye(na * na)
    basis = []
    for a in range(na):
        for b in range(a, na):
            v = np.zeros((na, na))
            if a == b: v[a, a] = 1
            else: v[a, b] = v[b, a] = 1 / np.sqrt(2)
            basis.append(v.ravel())
    B = np.array(basis).T                       # 36 x 21; sutun 0 = HF determinanti
    Hs = B.T @ Hdet @ B
    psi = [B.T @ np.asarray(c).ravel() for c in cS.ci[:2]]
    return Hs, psi, np.array(cS.e_tot)

def gozlenebilir(Hs, perm):
    from qiskit.quantum_info import SparsePauliOp
    H5 = np.diag(np.full(D, np.trace(Hs) / len(Hs)))   # kullanilmayan durumlar: ortalama enerji
    H5[np.ix_(perm, perm)] = Hs
    op = SparsePauliOp.from_operator(H5); op = SparsePauliOp(op.paulis, op.coeffs.real).simplify(atol=1e-12)
    keep = [i for i, (p, c) in enumerate(zip(op.paulis, op.coeffs)) if abs(c) >= ESIK or not (p.x.any() or p.z.any())]
    return SparsePauliOp(op.paulis[keep], op.coeffs[keep])

def yerlesim_ara(Hs, psi, ref_eV, n=300, seed=0):
    """21 durumun 32'lik tabandaki yerleri (HF = 0 sabit): kayma < 0.005 eV, en az olcum grubu."""
    rng = np.random.default_rng(seed); best = None
    for k in range(n):
        perm = np.arange(len(Hs)) if k == 0 else np.r_[0, 1 + rng.permutation(D - 1)[:len(Hs) - 1]]
        obs = gozlenebilir(Hs, perm); M = obs.to_matrix().real
        v = [np.zeros(D) for _ in psi]
        for s in range(2): v[s][perm] = psi[s]
        kayma = (v[1] @ M @ v[1] - v[0] @ M @ v[0]) * EV - ref_eV
        ng = len(obs.group_commuting(qubit_wise=True))
        if abs(kayma) < 0.005 and (best is None or ng < best[1]): best = (perm, ng, kayma)
    if best is None: sys.exit("Uygun yerlesim bulunamadi (kayma < 0.005 eV)")
    return best

# ---------------------------------------------------------------- devre
def devre(vec):
    """Gercek genlikli durum: Mottonen agaci (duzey k: 2^k RY + 2^k CX). Tum acilar Parameter."""
    from qiskit import QuantumCircuit
    from qiskit.circuit import Parameter
    qc = QuantumCircuit(NQ); vals = {}
    def ry(theta, q):
        P = Parameter(f"a{len(vals):03d}"); vals[P.name] = float(theta); qc.ry(P, q)
    for k in range(NQ):
        t, N = NQ - 1 - k, 2 ** k
        blk = vec.reshape(2 ** k, 2, 2 ** (NQ - 1 - k))
        if k == NQ - 1: a0, a1 = blk[:, 0, 0], blk[:, 1, 0]
        else: a0, a1 = np.linalg.norm(blk[:, 0, :], axis=1), np.linalg.norm(blk[:, 1, :], axis=1)
        alpha = 2 * np.arctan2(a1, a0)
        if k == 0: ry(alpha[0], t); continue
        ctrl = [t + 1 + l for l in range(k)]
        gray = [i ^ (i >> 1) for i in range(N)]
        M = np.array([[(-1) ** bin(j & gray[i]).count("1") for j in range(N)] for i in range(N)])
        theta = M @ alpha / N
        for i in range(N):
            ry(theta[i], t)
            qc.cx(ctrl[(gray[i] ^ gray[(i + 1) % N]).bit_length() - 1], t)
    return qc, vals

def bagla(qc, vals, sifir=False):
    return qc.assign_parameters({p: (0.0 if sifir else vals[p.name]) for p in qc.parameters})

def duzelt(E, Eref, Eref_id, c0):
    f = (Eref - c0) / (Eref_id - c0)
    return c0 + (E - c0) / f, f

def pubs_kur(D_, be):
    from qiskit import transpile
    tq = {nm: transpile(D_["q"][nm], be, optimization_level=3, seed_transpiler=1) for nm in ("S0", "S1")}
    pubs = [(bagla(tq[nm], D_["v"][nm]), D_["obs"].apply_layout(tq[nm].layout)) for nm in ("S0", "S1")]
    pubs.append((bagla(tq["S0"], D_["v"]["S0"], True), D_["obs"].apply_layout(tq["S0"].layout)))
    return pubs, tq

# ---------------------------------------------------------------- komutlar
if KMD == "hazirla":
    from qiskit.quantum_info import Statevector
    t0 = time.time()
    log(f"[{NAME}] 1) Hamiltonyen (B3LYP orbitalleri, CAS(4,4)) ve 21 durumluk singlet alt uzayi ...")
    Hs, psi, eS = hamiltonyen()
    ref = (eS[1] - eS[0]) * EV
    alt = (psi[1] @ Hs @ psi[1] - psi[0] @ Hs @ psi[0]) * EV
    log(f"   CASCI S1 = {ref:.4f} eV ; alt uzayda {alt:.4f} eV ; |<HF|S0>|^2 = {psi[0][0]**2:.3f}")
    log(f"[{NAME}] 2) 5 kubite yerlesim aramasi (esik {ESIK:g}) ...")
    perm, ng, kayma = yerlesim_ara(Hs, psi, ref)
    obs = gozlenebilir(Hs, perm); c0 = float(np.trace(obs.to_matrix()).real / D)
    vec = {}
    for s, nm in enumerate(("S0", "S1")):
        v = np.zeros(D); v[perm] = psi[s]; vec[nm] = v
    q, v = {}, {}
    for nm in ("S0", "S1"): q[nm], v[nm] = devre(vec[nm])
    ideal = {nm: Statevector(bagla(q[nm], v[nm])).expectation_value(obs).real for nm in ("S0", "S1")}
    ideal["ref"] = Statevector(bagla(q["S0"], v["S0"], True)).expectation_value(obs).real
    R = {"molekul": NAME, "yontem": "kompakt (5 kubit)", "CASCI_S1_eV": ref, "VQE_ideal_kesik_H_S1_eV": (ideal["S1"] - ideal["S0"]) * EV,
         "kesme_kaymasi_eV": kayma, "pauli_terim": len(obs), "olcum_grubu": ng, "cx": q["S0"].count_ops().get("cx", 0), "esik": ESIK}
    log(f"   {len(obs)} Pauli terimi, {ng} olcum grubu ; kesme kaymasi {kayma:+.4f} eV ; devre basina {R['cx']} CNOT")
    log(f"   ideal S1 = {R['VQE_ideal_kesik_H_S1_eV']:.4f} eV (CASCI {ref:.4f})")
    try:
        from qiskit_ibm_runtime.fake_provider import FakeKingston
        from qiskit_aer import AerSimulator
        from qiskit_aer.primitives import EstimatorV2 as AerEstimator
        be = FakeKingston()
        est = AerEstimator(options=dict(backend_options=dict(noise_model=AerSimulator.from_backend(be).options.noise_model), default_precision=0.0))
        D_ = dict(q=q, v=v, obs=obs)
        pubs, tq = pubs_kur(D_, be)
        r = est.run(pubs).result(); E0, E1, Er = (float(x.data.evs) for x in r)
        d0, _ = duzelt(E0, Er, ideal["ref"], c0); d1, _ = duzelt(E1, Er, ideal["ref"], c0)
        R.update(cz_S0=int(tq["S0"].count_ops().get("cz", 0)), cz_S1=int(tq["S1"].count_ops().get("cz", 0)),
                 sim_ham_S1_eV=(E1 - E0) * EV, sim_duzeltilmis_S1_eV=(d1 - d0) * EV)
        log(f"[{NAME}] 3) FakeKingston kapi gurultusu (atis gurultusu yok): 2-kubit kapi S0 {R['cz_S0']}, S1 {R['cz_S1']}")
        log(f"   S1 ham {R['sim_ham_S1_eV']:.3f} eV, duzeltilmis {R['sim_duzeltilmis_S1_eV']:.3f} eV  (CASCI {ref:.3f})")
    except ImportError as e:
        log(f"   (gurultulu simulasyon atlandi: {e})")
    pickle.dump(dict(q=q, v=v, obs=obs, c0=c0, ideal=ideal, R=R), open(PKL, "wb"))
    json.dump(R, open(f"{NAME}_kompakt_hazir.json", "w"), indent=1)
    log(f"[{NAME}] Hazir -> {PKL}  ({time.time()-t0:.0f} s)")
    log(f"   Donanim: 3 devre x {ng} grup x {SHOTS} atis ~ {3*ng*SHOTS/1e6:.2f} M atis")

elif KMD == "gonder":
    from qiskit_ibm_runtime import EstimatorV2 as Estimator
    D_ = pickle.load(open(PKL, "rb"))
    if SAHTE:
        from qiskit_ibm_runtime.fake_provider import FakeKingston
        be = FakeKingston(); log("SAHTE mod: FakeKingston (yerel, kredi harcamaz)")
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        svc = QiskitRuntimeService()
        be = svc.backend(BACKEND) if BACKEND else svc.least_busy(operational=True, simulator=False, min_num_qubits=NQ)
        log(f"Donanim: {be.name}  (bekleyen is: {be.status().pending_jobs})")
        try: log(f"Kalan Open Plan suresi: {svc.usage().get('usage_remaining_seconds')} s")
        except Exception: pass
    pubs, tq = pubs_kur(D_, be)
    for nm in ("S0", "S1"):
        log(f"  {nm}: 2-kubit kapi {tq[nm].count_ops().get('cz', 0) + tq[nm].count_ops().get('ecr', 0)}, derinlik {tq[nm].depth()}")
    ng = len(D_["obs"].group_commuting(qubit_wise=True))
    log(f"  Gonderilecek: 3 devre x {ng} olcum grubu x {SHOTS} atis (TREX + kapi twirling + dinamik dekuplaj).")
    if not SAHTE and input("  Gondereyim mi? (e/h): ").strip().lower() != "e":
        sys.exit("Iptal edildi.")
    est = Estimator(mode=be)
    est.options.default_shots = SHOTS
    if not SAHTE:
        est.options.resilience_level = 1
        est.options.twirling.enable_gates = True
        est.options.dynamical_decoupling.enable = True
    job = est.run(pubs)
    jid = job.job_id()
    json.dump({"job_id": jid, "backend": be.name, "meta": ["S0", "S1", "ref"], "shots": SHOTS, "zaman": time.ctime()}, open(JOBF, "w"))
    log(f"Is gonderildi: {jid}  -> {JOBF}")
    if SAHTE:
        r = job.result(); json.dump({"evs": [float(x.data.evs) for x in r], "stds": [float(x.data.stds) for x in r]}, open(f"{NAME}_kompakt_ham_sahte.json", "w"))
        log(f"SAHTE sonuc hazir: python3 ibm_vqe_kompakt.py sonuc {NAME} --sahte")
    else:
        log(f"Sirada bekleyebilir. Durum: https://quantum.cloud.ibm.com/workloads ; bitince: python3 ibm_vqe_kompakt.py sonuc {NAME}")

elif KMD == "sonuc":
    D_ = pickle.load(open(PKL, "rb")); J = json.load(open(JOBF)); R = dict(D_["R"])
    if SAHTE:
        H = json.load(open(f"{NAME}_kompakt_ham_sahte.json")); evs, stds = H["evs"], H["stds"]
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        job = QiskitRuntimeService().job(J["job_id"]); log(f"Is durumu: {job.status()}")
        r = job.result(); evs = [float(x.data.evs) for x in r]; stds = [float(x.data.stds) for x in r]
        try: R["QPU_saniye"] = job.usage()
        except Exception: pass
    E = dict(zip(J["meta"], evs)); S = dict(zip(J["meta"], stds)); c0, ideal = D_["c0"], D_["ideal"]
    d0, f = duzelt(E["S0"], E["ref"], ideal["ref"], c0); d1, _ = duzelt(E["S1"], E["ref"], ideal["ref"], c0)
    # hata yayilimi: S0, S1 ve ortak referansin istatistik hatalari (dogrusal yaklasim)
    g = (ideal["ref"] - c0) / (E["ref"] - c0)
    dS1 = (E["S1"] - E["S0"]) * g
    hata = np.sqrt((g * S["S0"]) ** 2 + (g * S["S1"]) ** 2 + (dS1 / (E["ref"] - c0) * S["ref"]) ** 2) * EV
    R.update(donanim=J["backend"], job_id=J["job_id"], atis=J["shots"],
             donanim_ham_S1_eV=(E["S1"] - E["S0"]) * EV, donanim_duzeltilmis_S1_eV=(d1 - d0) * EV,
             istatistik_hata_eV=float(hata), gurultu_f=float(f))
    json.dump(R, open(f"{NAME}_kompakt_sonuc{'_sahte' if SAHTE else ''}.json", "w"), indent=1)
    log(f"\n=== {NAME}: S1 uyarilma enerjisi, kompakt 5 kubit (NEVPT2 oncesi) ===")
    for k, ad in [("CASCI_S1_eV", "Tam cozum (CASCI, klasik)"), ("VQE_ideal_kesik_H_S1_eV", "Ideal simulator"),
                  ("sim_ham_S1_eV", "FakeKingston kapi gurultusu (ham)"), ("sim_duzeltilmis_S1_eV", "FakeKingston kapi gurultusu (duz.)"),
                  ("donanim_ham_S1_eV", f"{J['backend']} (ham)"), ("donanim_duzeltilmis_S1_eV", f"{J['backend']} (duzeltilmis)")]:
        if k in R: log(f"  {ad:38s} {R[k]:6.3f} eV   (tam cozumden fark {R[k]-R['CASCI_S1_eV']:+.3f})")
    log(f"  istatistik belirsizlik ~ {R['istatistik_hata_eV']:.3f} eV ; gurultu faktoru f = {f:.2f}")
    if NAME in ONCEKI and not SAHTE:
        o = ONCEKI[NAME]
        log(f"  Onceki yontem (ibm_vqe.py, {o['donanim']}): tam cozumden fark +{o['duzeltilmis_fark_eV']:.2f} +- {o['hata_eV']:.2f} eV")
else:
    sys.exit("komut: hazirla | gonder | sonuc")
