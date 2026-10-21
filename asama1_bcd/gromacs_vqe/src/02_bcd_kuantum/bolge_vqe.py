import sys, time
import numpy as np
from scipy.optimize import minimize
from pyscf import gto, scf, mcscf, ao2mo
from pyscf.mp.dfmp2 import DFMP2
import pennylane as qml
from pennylane import numpy as pnp
FN = sys.argv[1]                       # ornek: mod_A_bolge_G2.xyz
NG = 16                                # son 16 atom misafir
BASIS = "def2-svp"; HA = 627.5095
L = [l.split() for l in open(FN).read().split("\n")[2:] if l.strip()]
nA = len(L) - NG                       # A = glukoz, B = misafir
def mol(gA, gB):
    at = [("ghost-" if (gA and i < nA) or (gB and i >= nA) else "") + l[0] + " " + " ".join(l[1:4])
          for i, l in enumerate(L)]
    return gto.M(atom="; ".join(at), basis=BASIS, verbose=0, max_memory=4000)

def frag_pop(m, C):
    """Her MO'nun A (glukoz) atomlarinin baz fonksiyonlari uzerindeki Mulliken payi."""
    S = m.intor("int1e_ovlp"); SC = S @ C
    aoA = [i for i, lab in enumerate(m.ao_labels(fmt=False)) if lab[0] < nA]
    return (C[aoA] * SC[aoA]).sum(0)

def vqe(H, ne, nq):
    hf_state = qml.qchem.hf_state(ne, nq)
    singles, doubles = qml.qchem.excitations(ne, nq)
    dev = qml.device("default.qubit", wires=nq)
    @qml.qnode(dev, diff_method="backprop")
    def energy(t):
        qml.BasisState(hf_state, wires=range(nq))
        for i, d in enumerate(doubles): qml.DoubleExcitation(t[i], wires=d)
        for j, s in enumerate(singles): qml.SingleExcitation(t[len(doubles) + j], wires=s)
        return qml.expval(H)
    g = qml.grad(energy)
    f = lambda x: float(energy(pnp.array(x, requires_grad=True)))
    jac = lambda x: np.array(g(pnp.array(x, requires_grad=True)))
    r = minimize(f, np.zeros(len(singles) + len(doubles)), jac=jac, method="BFGS", options={"gtol": 1e-8})
    return r.fun, r.nit

out = {}; t0 = time.time()
for tag, gA, gB in [("AB", 0, 0), ("A", 0, 1), ("B", 1, 0)]:
    m = mol(gA, gB); mf = scf.RHF(m).density_fit().run()
    nocc = m.nelectron // 2; nmo = len(mf.mo_energy)
    ptf = DFMP2(mf); ptf.kernel(with_t2=False)          # tam MP2 (klasik referans)
    # VQE aktif uzayi: her parcadan 1 dolu + 1 bos orbital (dimerde 4e,4o / monomerde 2e,2o)
    pA = frag_pop(m, mf.mo_coeff)
    occ, vir = list(range(nocc)), list(range(nocc, nmo))
    sel = []
    for side in ([1] if tag == "A" else [0] if tag == "B" else [1, 0]):   # 1: A tarafi, 0: B tarafi
        isA = lambda i: pA[i] > 0.5 if side == 1 else pA[i] < 0.5
        sel.append(max(i for i in occ if isA(i))); sel.append(min(i for i in vir if isA(i)))
    sel = sorted(sel); ncas = len(sel); ne = ncas
    cas = mcscf.CASCI(mf, ncas, ne)
    mo = cas.sort_mo([i + 1 for i in sel]); cas.kernel(mo)
    h1, ecore = cas.get_h1eff(mo); eri = ao2mo.restore(1, cas.get_h2eff(mo), ncas)
    H = qml.qchem.qubit_observable(qml.qchem.fermionic_observable(np.array([ecore]), h1, np.swapaxes(eri, 1, 3)))
    ev, nit = vqe(H, ne, 2 * ncas)
    frz = [i for i in range(nmo) if i not in sel]
    pt = DFMP2(mf, frozen=frz); pt.kernel(with_t2=False)
    out[tag] = dict(hf=mf.e_tot, mp2=ptf.e_corr, cas=cas.e_tot - mf.e_tot, vqe=ev - mf.e_tot,
                    mp2act=pt.e_corr, sel=sel, pop=[round(float(pA[i]), 2) for i in sel], nq=2 * ncas, nit=nit)
    print(f"  {tag:2s}: {nmo} MO, aktif MO {[s+1 for s in sel]} (A payi {out[tag]['pop']}), "
          f"{2*ncas} kubit VQE {nit} iter  [{time.time()-t0:.0f} s]", flush=True)

I = lambda k: (out["AB"][k] - out["A"][k] - out["B"][k]) * HA
eHF, full = I("hf"), I("mp2")
print(f"\nmisafir + glukoz, {BASIS}, counterpoise  (kcal/mol)")
print(f"  HF etkilesimi                    : {eHF:8.3f}")
print(f"  MP2 korelasyon katkisi (tum {len(out['AB']['pop']) and nmo} MO) : {full:8.3f}")
print(f"  MP2 toplam (klasik)              : {eHF+full:8.3f}")
print("\n  Kuantum bolgesi (her parcadan 1 dolu + 1 bos orbital; dimer 8 kubit):")
print(f"    VQE korelasyon katkisi         : {I('vqe'):8.4f}")
print(f"    CASCI (kesin)                  : {I('cas'):8.4f}   (VQE - CASCI = {I('vqe')-I('cas'):.1e})")
print(f"    MP2, ayni 4 orbitalde          : {I('mp2act'):8.4f}")
dq = I("vqe") - I("mp2act")
print(f"    Kuantum duzeltmesi (VQE - MP2) : {dq:+8.4f}")
print(f"    Tam MP2 korelasyonunun yuzdesi : %{100*I('vqe')/full:.2f}")
print(f"\n  HIBRIT  E_int = HF + MP2(tam) + kuantum duzeltmesi = {eHF+full+dq:8.3f}")
print(f"  Tum bolge VQE ile yapilsaydi: {2*(nmo - 0)} kubit (def2-SVP, 348 MO)")
