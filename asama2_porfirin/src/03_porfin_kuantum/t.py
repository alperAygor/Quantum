# Porfin: kuantum bolgesi (Gouterman 4 orbitali = 8 kubit) VQE/VQD + klasik NEVPT2 duzeltmesi
import sys, time
NEV = len(sys.argv) > 1 and sys.argv[1] == 'nevpt2'   # istege bagli, uzun surer
import numpy as np
from pyscf import gto, scf, dft, mcscf, ao2mo, mrpt, fci
from pyscf.fci import cistring
import pennylane as qml
from pennylane import numpy as pnp
from scipy.optimize import minimize
EV = 27.211386; t0 = time.time()
mol = gto.M(atom="porfin.xyz", basis="def2-svp", verbose=0, max_memory=4000)
print(f"Porfin: {mol.natm} atom, {mol.nelectron} elektron, {mol.nao} baz fonksiyonu", flush=True)

# --- 1) Orbitaller: B3LYP
ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
nocc = mol.nelectron // 2; e = ks.mo_energy * EV
print(f"B3LYP: HOMO-1 {e[nocc-2]:.2f}  HOMO {e[nocc-1]:.2f}  LUMO {e[nocc]:.2f}  LUMO+1 {e[nocc+1]:.2f} eV  [{time.time()-t0:.0f} s]", flush=True)
mf = scf.RHF(mol).density_fit()

# --- 2) CASCI(4e,4o): singletler ve tripletler (klasik kesin cozum = VQE hedefi)
def casci(ss, nroots):
    c = mcscf.CASCI(mf, 4, 4); c.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=ss)
    c.fcisolver.nroots = nroots; c.kernel(ks.mo_coeff); return c
cS, cT = casci(0, 5), casci(2, 2)
eS, eT = np.array(cS.e_tot), np.array(cT.e_tot)
def configs(ci):
    """Uzaysal konfigurasyon agirliklari (alfa+beta doluluklari, determinant ciftleri birlestirilir)."""
    sa = cistring.make_strings(range(4), 2); occ = {}
    for i, a in enumerate(sa):
        for j, b in enumerate(sa):
            key = tuple(((a >> k) & 1) + ((b >> k) & 1) for k in range(4))
            occ[key] = occ.get(key, 0) + ci[i, j] ** 2
    return sorted(occ.items(), key=lambda x: -x[1])
print(f"CASCI bitti  [{time.time()-t0:.0f} s]", flush=True)

# --- 3) Kubit Hamiltonyeni + VQE/VQD (spin cezali)
h1, ecore = cS.get_h1eff(ks.mo_coeff); eri = ao2mo.restore(1, cS.get_h2eff(ks.mo_coeff), 4)
H = qml.qchem.qubit_observable(qml.qchem.fermionic_observable(np.array([ecore]), h1, np.swapaxes(eri, 1, 3)))
Hm = qml.matrix(H, wire_order=range(8)).real
S2m = qml.matrix(qml.qchem.spin2(4, 8), wire_order=range(8)).real
print(f"Kubit Hamiltonyeni: 8 kubit, {len(H.operands)} Pauli terimi", flush=True)
hf_state = qml.qchem.hf_state(4, 8); singles, doubles = qml.qchem.excitations(4, 8)
NL = 2; npar = NL * (len(singles) + len(doubles))
dev = qml.device("default.qubit", wires=8)
@qml.qnode(dev, diff_method="backprop")
def circuit(t):
    qml.BasisState(hf_state, wires=range(8)); k = 0
    for _ in range(NL):
        for d in doubles: qml.DoubleExcitation(t[k], wires=d); k += 1
        for s in singles: qml.SingleExcitation(t[k], wires=s); k += 1
    return qml.state()
rng = np.random.default_rng(7)
def vqd(nstate, s2target):
    P = (S2m - s2target * np.eye(256)) @ (S2m - s2target * np.eye(256))
    found = []
    for k in range(nstate):
        def cost(t):
            psi = circuit(t); c = pnp.real(pnp.conj(psi) @ (Hm + 1.0 * P) @ psi)
            for p in found: c = c + 3.0 * pnp.abs(pnp.conj(p) @ psi) ** 2
            return c
        g = qml.grad(cost); best = None
        for trial in range(2):
            x0 = np.zeros(npar) if (k == 0 and trial == 0 and s2target == 0) else rng.normal(0, 0.5, npar)
            r = minimize(lambda x: float(cost(pnp.array(x, requires_grad=True))), x0,
                         jac=lambda x: np.array(g(pnp.array(x, requires_grad=True))), method="BFGS", options={"gtol": 1e-7})
            if best is None or r.fun < best.fun: best = r
        found.append(np.array(circuit(best.x)))
    return (np.array([np.real(np.conj(p) @ Hm @ p) for p in found]),
            np.array([np.real(np.conj(p) @ S2m @ p) for p in found]), found)
vS, sS, psiS = vqd(5, 0); vT, sT, psiT = vqd(2, 2)
print(f"VQE/VQD bitti  [{time.time()-t0:.0f} s]", flush=True)
# gecis dipolleri: dipol operatoru de kubitlere cevrilir, <S0|mu|Si> kuantum durumlarindan hesaplanir
mo_act = ks.mo_coeff[:, cS.ncore:cS.ncore + 4]
dact = np.einsum("xij,ip,jq->xpq", mol.intor("int1e_r"), mo_act, mo_act)
Dm = [qml.matrix(qml.qchem.qubit_observable(qml.qchem.fermionic_observable(np.array([0.0]), dact[x])), wire_order=range(8)).real for x in range(3)]
fosc, pol = [0.0], [""]
for i in range(1, 5):
    mu = -np.array([np.real(np.conj(psiS[0]) @ Dm[x] @ psiS[i]) for x in range(3)])
    fosc.append(2 / 3 * (vS[i] - vS[0]) * mu @ mu); pol.append("xyz"[int(np.argmax(np.abs(mu)))])

# --- 4) Klasik dinamik korelasyon duzeltmesi: NEVPT2 (her hal icin)
def nev(c, n): return np.array([mrpt.NEVPT(c, root=i).kernel() for i in range(n)])
if NEV:
    nS = np.concatenate([nev(cS, 2), np.full(3, np.nan)]); nT = np.append(nev(cT, 1), np.nan); print(f"NEVPT2 bitti  [{time.time()-t0:.0f} s]", flush=True)
else:
    nS, nT = np.full(5, np.nan), np.full(2, np.nan)

# --- 5) Rapor
lbl = ["H-1", "H", "L", "L+1"]
def desc(cfg):
    out = []
    for occ, w in cfg[:2]:
        frm = [lbl[i] for i in range(2) for _ in range(2 - occ[i])]
        to = [lbl[i] for i in range(2, 4) for _ in range(occ[i])]
        out.append(("HF" if not frm else ",".join(frm) + "->" + ",".join(to)) + f" {w:.2f}")
    return " | ".join(out)
E0q, E0n = vS[0], eS[0] + nS[0]
print("\nPORFIN UYARILMIS HALLERI (eV, taban hale gore)")
print(f"{'hal':>4} {'VQE/VQD':>8} {'CASCI':>8} {'fark':>8} {'+NEVPT2':>8} {'f (siddet)':>11}   baskin konfigurasyonlar (agirlik)")
names = ["S0", "S1", "S2", "S3", "S4"]
for i in range(5):
    print(f"{names[i]:>4} {(vS[i]-vS[0])*EV:8.3f} {(eS[i]-eS[0])*EV:8.3f} {(vS[i]-eS[i])*EV:8.1e} {(eS[i]+nS[i]-E0n)*EV:8.3f} {fosc[i]:8.3f} {pol[i]:1s}   {desc(configs(cS.ci[i]))}")
for i in range(2):
    print(f"{'T'+str(i+1):>4} {(vT[i]-vS[0])*EV:8.3f} {(eT[i]-eS[0])*EV:8.3f} {(vT[i]-eT[i])*EV:8.1e} {(eT[i]+nT[i]-E0n)*EV:8.3f} {'-':>8}     {desc(configs(cT.ci[i]))}")
print(f"\nVQD spin kontrolu <S^2>: singletler {np.round(sS,3)}, tripletler {np.round(sT,3)}")
print(f"Klasik tek-determinant tahmini (B3LYP HOMO-LUMO): {e[nocc]-e[nocc-1]:.2f} eV")
print("Deney (gaz fazi): Qx ~1.98, Qy ~2.42 (zayif), Soret/B ~3.33 eV (cok siddetli)")
T1 = (eT[0] + nT[0] - E0n) * EV if NEV else (vT[0] - vS[0]) * EV
print(f"PDT: T1 ({'CAS+NEVPT2' if NEV else 'VQD, sadece aktif uzay'}) = {T1:.2f} eV; singlet oksijen (1O2) icin T1 > 0.98 eV -> {'EVET' if T1 > 0.98 else 'HAYIR'}")
print(f"Toplam sure: {time.time()-t0:.0f} s")
