# Gercek kuantum bilgisayarda (IBM Quantum) VQE: S0 ve S1 enerjileri -> uyarilma enerjisi
#
#   python3 ibm_vqe.py hazirla porfin     # ucretsiz, yerel: kompakt devre + gurultulu simulasyon tahmini
#   python3 ibm_vqe.py gonder  porfin     # IBM donanimina is gonderir (Open Plan suresinden yer)
#   python3 ibm_vqe.py sonuc   porfin     # sonucu alir, gurultu duzeltmesi yapar, karsilastirir
#   python3 ibm_vqe.py gonder  porfin --sahte   # ayni akisi FakeKingston uzerinde yerel dener (kredi harcamaz)
#   secenekler: --backend=ibm_kingston  --shots=4096  --zne (gurultu ekstrapolasyonu, ~3 kat sure)
#
# Yontem:
#  * Kuantum bolgesi: CAS(4e,4o) = 8 kubit, B3LYP orbitalleri (v3 taramasiyla ayni).
#  * Devre: ADAPT ile secilen az sayida Givens kapisi (S0: ~8, S1: ~5 cift uyarma) -> tum-kapi VQE'den cok daha sig.
#    Parametreler klasik simulatorde optimize edilir; donanimda enerji olculur (standart "donanimda VQE degerlendirmesi").
#  * Gozlenebilir: |c| >= 1e-3 Pauli terimleri (52 olcum grubu; kesme hatasi ~0.01 eV, ideal referansta da ayni H kullanilir).
#  * Gurultu: TREX (olcum), kapi twirling, dinamik dekuplaj + gurultu tahmin devresi (ayni devre, acilar 0):
#    E_duz = c0 + (E_olc - c0) * (E_ref_ideal - c0) / (E_ref_olc - c0)   [depolarizan gurultu varsayimi]
import os, sys, json, time, pickle
os.makedirs(os.path.expanduser("~/pyscf_tmp"), exist_ok=True)
os.environ.setdefault("PYSCF_TMPDIR", os.path.expanduser("~/pyscf_tmp"))
import numpy as np
from pyscf import gto, scf, dft, mcscf, fci, ao2mo  # qiskit/pennylane oncesi yuklenmeli (libgomp TLS hatasi)

EV = 27.211386
KMD, NAME = sys.argv[1], sys.argv[2]
SAHTE = "--sahte" in sys.argv
ARG = {a.split("=")[0]: a.split("=")[1] for a in sys.argv[3:] if "=" in a}
ZNE = "--zne" in sys.argv
SHOTS = int(ARG.get("--shots", 2048 if ZNE else 4096)); BACKEND = ARG.get("--backend")
ESIK = float(ARG.get("--esik", 1e-3))
PKL = f"{NAME}_ibm_hazir.pkl"; JOBF = f"{NAME}_ibm_is{'_sahte' if SAHTE else ''}.json"

def log(s): print(s, flush=True)

# ---------------------------------------------------------------- Hamiltonyen
def hamiltonyen():
    from pyscf import gto, scf, dft, mcscf, fci, ao2mo
    lines = [l for l in open(f"{NAME}.xyz").read().split("\n")[2:] if l.strip()]
    ecp = "def2-svp" if any(l.split()[0] == "Pd" for l in lines) else None
    mol = gto.M(atom=f"{NAME}.xyz", basis="def2-svp", ecp=ecp, verbose=0, max_memory=int(os.environ.get("QMEM", "6000")))
    ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
    mf = scf.RHF(mol).density_fit(); mo = ks.mo_coeff
    cS = mcscf.CASCI(mf, 4, 4); cS.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0)
    cS.fcisolver.nroots = 3; cS.kernel(mo)
    h1, ec = cS.get_h1eff(mo); eri = ao2mo.restore(1, cS.get_h2eff(mo), 4)
    return h1, float(ec), eri, np.array(cS.e_tot)

def gozlenebilir(h1, ec, eri, esik):
    import pennylane as qml
    from qiskit.quantum_info import SparsePauliOp
    H = qml.qchem.qubit_observable(qml.qchem.fermionic_observable(np.array([ec]), h1, np.swapaxes(eri, 1, 3)))
    ident, L = 0.0, []
    for pw, c in H.pauli_rep.items():
        if len(pw) == 0: ident += float(np.real(c)); continue
        if abs(c) >= esik: L.append(("".join(pw[w] for w in pw), [int(w) for w in pw], float(np.real(c))))
    return SparsePauliOp.from_sparse_list(L + [("", [], ident)], num_qubits=8).simplify(), ident, H

# ---------------------------------------------------------------- kompakt ADAPT
def adapt(Hm, M, ref, hedef, maxop=14, tol=1.6e-3):
    import pennylane as qml
    from pennylane import numpy as pnp
    from scipy.optimize import minimize
    sing, dbl = qml.qchem.excitations(4, 8)
    pool = [("D", x) for x in dbl] + [("S", x) for x in sing]
    dev = qml.device("default.qubit", wires=8)
    def make(ops):
        @qml.qnode(dev, diff_method="backprop")
        def c(t):
            qml.BasisState(np.array(ref), wires=range(8))
            for k, (kind, wr) in enumerate(ops):
                (qml.DoubleExcitation if kind == "D" else qml.SingleExcitation)(t[k], wires=wr)
            return qml.state()
        return c
    ops, th = [], np.zeros(0)
    for it in range(maxop):
        gr = []
        for op in pool:
            c2 = make(ops + [op])
            e = lambda x: float(np.real(np.conj(c2(pnp.array(np.r_[th, x]))) @ M @ c2(pnp.array(np.r_[th, x]))))
            gr.append(abs(e(1e-3) - e(-1e-3)) / 2e-3)
        ops.append(pool[int(np.argmax(gr))]); th = np.r_[th, 0.0]; c = make(ops)
        f = lambda x: pnp.real(pnp.conj(c(x)) @ M @ c(x))
        g = qml.grad(f)
        r = minimize(lambda x: float(f(pnp.array(x, requires_grad=True))), th,
                     jac=lambda x: np.array(g(pnp.array(x, requires_grad=True))), method="BFGS")
        th = r.x; psi = np.array(c(pnp.array(th))); E = float(np.real(np.conj(psi) @ Hm @ psi))
        log(f"   kapi {it+1:2d}: {ops[-1]}  hedeften fark {(E-hedef)*1000:7.2f} mHa")
        if E - hedef < tol: break
    return ops, th, psi, E

def devre(ops, th, ref):
    """PennyLane -> Qiskit; tum sayisal acilar Parameter yapilir (gurultu tahmin devresi icin acilar 0'a baglanir)."""
    import pennylane as qml
    from qiskit import QuantumCircuit
    from qiskit.circuit import Parameter
    dev = qml.device("default.qubit", wires=8)
    @qml.qnode(dev)
    def c():
        qml.BasisState(np.array(ref), wires=range(8))
        for k, (kind, wr) in enumerate(ops):
            (qml.DoubleExcitation if kind == "D" else qml.SingleExcitation)(th[k], wires=wr)
        return qml.state()
    qc = QuantumCircuit.from_qasm_str(qml.to_openqasm(c, measure_all=False)())
    new = QuantumCircuit(8); vals = {}
    for inst in qc.data:
        op = inst.operation
        if op.params and all(isinstance(p, (float, int)) for p in op.params):
            op = op.copy(); ps = []
            for p in op.params:
                P = Parameter(f"a{len(vals):03d}"); vals[P.name] = float(p); ps.append(P)
            op.params = ps
        new.append(op, inst.qubits)
    return new, vals

def bagla(qc, vals, sifir=False):
    return qc.assign_parameters({p: (0.0 if sifir else vals[p.name]) for p in qc.parameters})

def duzelt(E, Eref, Eref_id, c0):
    f = (Eref - c0) / (Eref_id - c0)
    return c0 + (E - c0) / f, f

# ---------------------------------------------------------------- komutlar
if KMD == "hazirla":
    import pennylane as qml
    from qiskit.quantum_info import Statevector
    t0 = time.time()
    log(f"[{NAME}] 1) Hamiltonyen (B3LYP orbitalleri, CAS(4,4)) ...")
    h1, ec, eri, eS = hamiltonyen()
    obs_tam, c0, Hpl = gozlenebilir(h1, ec, eri, 0.0); obs, _, _ = gozlenebilir(h1, ec, eri, ESIK)
    Hm = qml.matrix(Hpl, wire_order=range(8)).real; S2m = qml.matrix(qml.qchem.spin2(4, 8), wire_order=range(8)).real
    log(f"   CASCI S1 = {(eS[1]-eS[0])*EV:.3f} eV ; olculecek Pauli terimi {len(obs)}, olcum grubu {len(obs.group_commuting(qubit_wise=True))}")
    log(f"[{NAME}] 2) Kompakt ADAPT-VQE: S0 (HF referansi)")
    hf = [1, 1, 1, 1, 0, 0, 0, 0]
    ops0, th0, psi0, E0 = adapt(Hm, Hm, hf, eS[0])
    log(f"[{NAME}]    S1 (VQD: HOMO->LUMO referansi, S^2 ve S0'a dikeylik cezasi)")
    ref1 = [1, 1, 1, 0, 0, 1, 0, 0]
    ops1, th1, psi1, E1 = adapt(Hm, Hm + S2m @ S2m + 3.0 * np.outer(psi0, psi0), ref1, eS[1], maxop=16)
    q0, v0 = devre(ops0, th0, hf); q1, v1 = devre(ops1, th1, ref1)
    ideal = {}
    for nm, q, v in (("S0", q0, v0), ("S1", q1, v1)):
        ideal[nm] = Statevector(bagla(q, v)).expectation_value(obs).real
        ideal[nm + "_ref"] = Statevector(bagla(q, v, True)).expectation_value(obs).real
    R = {"molekul": NAME, "CASCI_S1_eV": (eS[1]-eS[0])*EV, "VQE_ideal_S1_eV": (E1-E0)*EV,
         "VQE_ideal_kesik_H_S1_eV": (ideal["S1"]-ideal["S0"])*EV, "kapi_S0": len(ops0), "kapi_S1": len(ops1),
         "pauli_terim": len(obs), "esik": ESIK}
    log(f"   ideal VQE S1 = {R['VQE_ideal_S1_eV']:.3f} eV (CASCI {R['CASCI_S1_eV']:.3f}); kesik H ile {R['VQE_ideal_kesik_H_S1_eV']:.3f} eV")
    # 3) Gurultulu simulasyon tahmini (FakeKingston = ibm_kingston gurultu modeli)
    try:
        from qiskit import transpile
        from qiskit_ibm_runtime.fake_provider import FakeKingston
        from qiskit_aer import AerSimulator
        from qiskit_aer.primitives import EstimatorV2 as AerEstimator
        be = FakeKingston(); est = AerEstimator(options=dict(backend_options=dict(noise_model=AerSimulator.from_backend(be).options.noise_model), default_precision=0.004))
        sim = {}
        for nm, q, v in (("S0", q0, v0), ("S1", q1, v1)):
            t = transpile(q, be, optimization_level=3, seed_transpiler=1); ob = obs.apply_layout(t.layout)
            r = est.run([(bagla(t, v), ob), (bagla(t, v, True), ob)]).result()
            E, Er = float(r[0].data.evs), float(r[1].data.evs)
            sim[nm] = (E, duzelt(E, Er, ideal[nm + "_ref"], c0)[0]); R[f"cz_{nm}"] = int(t.count_ops().get("cz", 0))
        R["sim_ham_S1_eV"] = (sim["S1"][0]-sim["S0"][0])*EV; R["sim_duzeltilmis_S1_eV"] = (sim["S1"][1]-sim["S0"][1])*EV
        log(f"[{NAME}] 3) FakeKingston gurultulu simulasyon: 2-kubit kapi S0 {R['cz_S0']}, S1 {R['cz_S1']}")
        log(f"   S1 ham {R['sim_ham_S1_eV']:.2f} eV, duzeltilmis {R['sim_duzeltilmis_S1_eV']:.2f} eV  (ideal {R['VQE_ideal_kesik_H_S1_eV']:.2f})")
    except ImportError as e:
        log(f"   (gurultulu simulasyon atlandi: {e})")
    pickle.dump(dict(q0=q0, v0=v0, q1=q1, v1=v1, obs=obs, c0=c0, ideal=ideal, R=R), open(PKL, "wb"))
    json.dump(R, open(f"{NAME}_ibm_hazir.json", "w"), indent=1)
    ng = len(obs.group_commuting(qubit_wise=True))
    log(f"[{NAME}] Hazir -> {PKL}  ({time.time()-t0:.0f} s)")
    log(f"   Donanim maliyeti tahmini: 4 devre x {ng} grup x {SHOTS} atis ~ {4*ng*SHOTS/1e6:.2f} M atis ~ 3-6 dk QPU")

elif KMD == "gonder":
    from qiskit import transpile
    from qiskit_ibm_runtime import EstimatorV2 as Estimator
    D = pickle.load(open(PKL, "rb"))
    if SAHTE:
        from qiskit_ibm_runtime.fake_provider import FakeKingston
        be = FakeKingston(); log("SAHTE mod: FakeKingston (yerel, kredi harcamaz)")
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        svc = QiskitRuntimeService()
        be = svc.backend(BACKEND) if BACKEND else svc.least_busy(operational=True, simulator=False, min_num_qubits=8)
        log(f"Donanim: {be.name}  (bekleyen is: {be.status().pending_jobs})")
    pubs, meta = [], []
    for nm in ("S0", "S1"):
        q, v = D["q" + nm[-1]], D["v" + nm[-1]]
        t = transpile(q, be, optimization_level=3, seed_transpiler=1); ob = D["obs"].apply_layout(t.layout)
        pubs += [(bagla(t, v), ob), (bagla(t, v, True), ob)]; meta += [nm, nm + "_ref"]
        log(f"  {nm}: 2-kubit kapi {t.count_ops().get('cz', 0) + t.count_ops().get('ecr', 0)}, derinlik {t.depth()}")
    ng = len(D["obs"].group_commuting(qubit_wise=True))
    log(f"  Gonderilecek: 4 devre x {ng} olcum grubu x {SHOTS} atis{' + ZNE' if ZNE else ''}. Tahmini QPU suresi {'6-9' if ZNE else '3-6'} dk (Open Plan: 10 dk/ay).")
    if not SAHTE and input("  Gondereyim mi? (e/h): ").strip().lower() != "e":
        sys.exit("Iptal edildi.")
    est = Estimator(mode=be)
    est.options.default_shots = SHOTS
    if not SAHTE:
        est.options.resilience_level = 2 if ZNE else 1   # 1: TREX olcum duzeltmesi ; 2: + ZNE (sure ~3 kat)
        est.options.twirling.enable_gates = True
        est.options.dynamical_decoupling.enable = True
    job = est.run(pubs)
    jid = job.job_id(); json.dump({"job_id": jid, "backend": be.name, "meta": meta, "shots": SHOTS, "zaman": time.ctime()}, open(JOBF, "w"))
    log(f"Is gonderildi: {jid}  -> {JOBF}")
    if SAHTE:
        r = job.result(); json.dump({"evs": [float(x.data.evs) for x in r], "stds": [float(x.data.stds) for x in r]}, open(f"{NAME}_ibm_ham_sahte.json", "w"))
        log("SAHTE sonuc hazir: python3 ibm_vqe.py sonuc " + NAME + " --sahte")
    else:
        log("Sirada bekleyebilir. Durum: https://quantum.cloud.ibm.com/workloads  ; bitince: python3 ibm_vqe.py sonuc " + NAME)

elif KMD == "sonuc":
    D = pickle.load(open(PKL, "rb")); J = json.load(open(JOBF)); R = dict(D["R"])
    if SAHTE:
        H = json.load(open(f"{NAME}_ibm_ham_sahte.json")); evs, stds = H["evs"], H["stds"]
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        job = QiskitRuntimeService().job(J["job_id"]); log(f"Is durumu: {job.status()}")
        r = job.result(); evs = [float(x.data.evs) for x in r]; stds = [float(x.data.stds) for x in r]
        try: R["QPU_saniye"] = job.usage()
        except Exception: pass
    E = dict(zip(J["meta"], evs)); S = dict(zip(J["meta"], stds)); c0, ideal = D["c0"], D["ideal"]
    d0, f0 = duzelt(E["S0"], E["S0_ref"], ideal["S0_ref"], c0); d1, f1 = duzelt(E["S1"], E["S1_ref"], ideal["S1_ref"], c0)
    R.update(donanim=J["backend"], job_id=J["job_id"], atis=J["shots"],
             donanim_ham_S1_eV=(E["S1"]-E["S0"])*EV, donanim_duzeltilmis_S1_eV=(d1-d0)*EV,
             istatistik_hata_eV=float(np.hypot(S["S0"], S["S1"]) * EV), gurultu_f=[f0, f1])
    json.dump(R, open(f"{NAME}_ibm_sonuc{'_sahte' if SAHTE else ''}.json", "w"), indent=1)
    log(f"\n=== {NAME}: S1 uyarilma enerjisi (8 kubitlik kuantum bolgesi, NEVPT2 oncesi) ===")
    for k, ad in [("CASCI_S1_eV", "Tam cozum (CASCI, klasik)"), ("VQE_ideal_kesik_H_S1_eV", "VQE, ideal simulator"),
                  ("sim_ham_S1_eV", "VQE, gurultulu sim (ham)"), ("sim_duzeltilmis_S1_eV", "VQE, gurultulu sim (duzeltilmis)"),
                  ("donanim_ham_S1_eV", f"VQE, {J['backend']} (ham)"), ("donanim_duzeltilmis_S1_eV", f"VQE, {J['backend']} (duzeltilmis)")]:
        if k in R: log(f"  {ad:38s} {R[k]:6.3f} eV   (tam cozumden fark {R[k]-R['CASCI_S1_eV']:+.3f})")
    log(f"  istatistik belirsizlik ~ {R['istatistik_hata_eV']:.3f} eV ; gurultu faktoru f(S0,S1) = {f0:.2f}, {f1:.2f}")
    log("  Not: donanim farki (dE_donanim - dE_CASCI) dogrudan NEVPT2'li S1'e eklenebilir: S1_donanim = S1_NEVPT2 + (dE_donanim - dE_CASCI)")
else:
    sys.exit("komut: hazirla | gonder | sonuc")
