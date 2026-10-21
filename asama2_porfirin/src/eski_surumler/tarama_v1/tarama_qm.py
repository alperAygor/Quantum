# Kuantum triyaj: tek molekul. Kullanim: python3 tarama_qm.py porfin   (porfin.xyz okunur)
import os, sys, time, json
os.makedirs(os.path.expanduser("~/pyscf_tmp"), exist_ok=True)
os.environ.setdefault("PYSCF_TMPDIR", os.path.expanduser("~/pyscf_tmp"))
import numpy as np
from pyscf import gto, scf, dft, mcscf, fci, mrpt, ao2mo
from pyscf.fci import cistring
import pennylane as qml
from pennylane import numpy as pnp
from scipy.optimize import minimize
NAME = sys.argv[1]; NEV = "nonev" not in sys.argv
EV = 27.211386; t0 = time.time(); R = {"molekul": NAME}
def log(s): print(f"[{NAME}] {s}  [{time.time()-t0:.0f} s]", flush=True)
ecp = "def2-svp" if any(l.split()[0] == "Pd" for l in open(f"{NAME}.xyz").read().split("\n")[2:] if l.strip()) else None
mol = gto.M(atom=f"{NAME}.xyz", basis="def2-svp", ecp=ecp, verbose=0, max_memory=int(os.environ.get("QMEM", "6000")))
R.update(atom=mol.natm, elektron=mol.nelectron, baz=mol.nao)

# 1) Baslangic orbitalleri (B3LYP) ve aktif orbitallerdeki metal payi
ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
nocc = mol.nelectron // 2; act = list(range(nocc - 2, nocc + 2))
R["HOMO_LUMO_eV"] = float((ks.mo_energy[nocc] - ks.mo_energy[nocc-1]) * EV)
metal = [i for i in range(mol.natm) if mol.atom_symbol(i) in ("Zn", "Mg", "Pd")]
if metal:
    S = mol.intor("int1e_ovlp"); C = ks.mo_coeff
    ao = [k for k, lab in enumerate(mol.ao_labels(fmt=False)) if lab[0] in metal]
    R["metal_payi"] = [round(float((C[ao, i] * (S @ C[:, i])[ao]).sum()), 3) for i in act]
log("B3LYP bitti")

# 2) SA-CASSCF(4,4): S0, S1 (singlet) + T1 (triplet) ortak orbitaller
mf = scf.RHF(mol).density_fit()
sS = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0); sS.nroots = 2
sT = fci.direct_spin1.FCI(mol); sT.spin = 2; sT.nroots = 1
mc = mcscf.CASSCF(mf, 4, 4); mcscf.addons.state_average_mix_(mc, [sS, sT], [1/3, 1/3, 1/3])
mc.max_cycle_macro = 50; mc.kernel(ks.mo_coeff); mo = mc.mo_coeff
R["SA_CASSCF_yakinsadi"] = bool(mc.converged); log(f"SA-CASSCF bitti (yakinsadi: {mc.converged})")

# 3) CASCI referanslari
cS = mcscf.CASCI(mf, 4, 4); cS.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0)
cS.fcisolver.nroots = 3; cS.kernel(mo)
cT = mcscf.CASCI(mf, 4, (3, 1)); cT.kernel(mo)
eS, eT = np.array(cS.e_tot), float(cT.e_tot)
def cfg(ci):
    sa = cistring.make_strings(range(4), 2); occ = {}
    for i, a in enumerate(sa):
        for j, b in enumerate(sa):
            k = tuple(((a >> q) & 1) + ((b >> q) & 1) for q in range(4)); occ[k] = occ.get(k, 0) + ci[i, j] ** 2
    return sorted(occ.values(), reverse=True)
w = cfg(cS.ci[1]); R["S1_konfig_agirlik"] = [round(float(x), 3) for x in w[:2]]
R["S0_HF_agirlik"] = round(float(cfg(cS.ci[0])[0]), 3)
log("CASCI bitti")

# 4) Kuantum: 8 kubit VQE/VQD (S0, S1, S2 singlet; T1 triplet) + gecis dipolleri
h1, ec = cS.get_h1eff(mo); eri = ao2mo.restore(1, cS.get_h2eff(mo), 4)
H = qml.qchem.qubit_observable(qml.qchem.fermionic_observable(np.array([ec]), h1, np.swapaxes(eri, 1, 3)))
Hm = qml.matrix(H, wire_order=range(8)).real; S2m = qml.matrix(qml.qchem.spin2(4, 8), wire_order=range(8)).real
hf = qml.qchem.hf_state(4, 8); sing, dbl = qml.qchem.excitations(4, 8); NL = 2; npar = NL * (len(sing) + len(dbl))
dev = qml.device("default.qubit", wires=8)
@qml.qnode(dev, diff_method="backprop")
def circ(t):
    qml.BasisState(hf, wires=range(8)); k = 0
    for _ in range(NL):
        for d in dbl: qml.DoubleExcitation(t[k], wires=d); k += 1
        for s in sing: qml.SingleExcitation(t[k], wires=s); k += 1
    return qml.state()
rng = np.random.default_rng(7)
def vqd(n, s2):
    P = (S2m - s2 * np.eye(256)) @ (S2m - s2 * np.eye(256)); found = []; pars = []
    for k in range(n):
        def cost(t):
            psi = circ(t); c = pnp.real(pnp.conj(psi) @ (Hm + P) @ psi)
            for p in found: c = c + 3.0 * pnp.abs(pnp.conj(p) @ psi) ** 2
            return c
        g = qml.grad(cost); best = None
        for tr in range(2):
            x0 = np.zeros(npar) if (k == 0 and tr == 0 and s2 == 0) else rng.normal(0, 0.5, npar)
            r = minimize(lambda x: float(cost(pnp.array(x, requires_grad=True))), x0,
                         jac=lambda x: np.array(g(pnp.array(x, requires_grad=True))), method="BFGS", options={"gtol": 1e-7})
            if best is None or r.fun < best.fun: best = r
        found.append(np.array(circ(best.x))); pars.append(best.x.tolist())
    return np.array([np.real(np.conj(p) @ Hm @ p) for p in found]), found, pars
vS, psiS, parS = vqd(3, 0); vT, _, parT = vqd(1, 2)
R["VQD_CASCI_maxfark_eV"] = float(max(np.abs(vS - eS).max(), abs(vT[0] - eT)) * EV)
mo_act = mo[:, cS.ncore:cS.ncore + 4]; dact = np.einsum("xij,ip,jq->xpq", mol.intor("int1e_r"), mo_act, mo_act)
Dm = [qml.matrix(qml.qchem.qubit_observable(qml.qchem.fermionic_observable(np.array([0.0]), dact[x])), wire_order=range(8)).real for x in range(3)]
for i in (1, 2):
    mu = -np.array([np.real(np.conj(psiS[0]) @ Dm[x] @ psiS[i]) for x in range(3)])
    R[f"f_S{i}"] = float(2 / 3 * (vS[i] - vS[0]) * mu @ mu)
R["S1_CAS_eV"] = float((vS[1] - vS[0]) * EV); R["S2_CAS_eV"] = float((vS[2] - vS[0]) * EV); R["T1_CAS_eV"] = float((vT[0] - vS[0]) * EV)
json.dump({"S": parS, "T": parT}, open(f"{NAME}_vqe_parametreleri.json", "w"))   # gercek kuantum donanimi icin
log("VQE/VQD bitti")

# 5) NEVPT2 (klasik dinamik korelasyon duzeltmesi)
if NEV:
    n0 = mrpt.NEVPT(cS, root=0).kernel(); n1 = mrpt.NEVPT(cS, root=1).kernel(); nT = mrpt.NEVPT(cT).kernel()
    R["S1_NEVPT2_eV"] = float((eS[1] + n1 - eS[0] - n0) * EV); R["T1_NEVPT2_eV"] = float((eT + nT - eS[0] - n0) * EV)
    log("NEVPT2 bitti")

# 6) Triyaj karari
mix = R["S1_konfig_agirlik"][1] if len(R["S1_konfig_agirlik"]) > 1 else 0
R["triyaj"] = "KUANTUM GEREKLI" if mix > 0.2 else "klasik yeterli"
if NEV: R["T1_gt_S1_hata"] = R["T1_NEVPT2_eV"] > R["S1_NEVPT2_eV"]; R["singlet_oksijen"] = R["T1_NEVPT2_eV"] > 0.98
R["sure_s"] = round(time.time() - t0)
json.dump(R, open(f"{NAME}_sonuc.json", "w"), indent=1); log("SONUC: " + json.dumps(R))
