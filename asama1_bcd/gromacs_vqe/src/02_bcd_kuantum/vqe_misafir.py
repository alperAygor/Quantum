import sys, time
import numpy as np
from scipy.optimize import minimize
from pyscf import gto, scf, mcscf, ao2mo
import pennylane as qml
from pennylane import numpy as pnp

fn = sys.argv[1] if len(sys.argv) > 1 else "mod_A.xyz"
NE = int(sys.argv[2]) if len(sys.argv) > 2 else 4   # aktif elektron
NO = int(sys.argv[3]) if len(sys.argv) > 3 else 4   # aktif orbital
HA = 627.5095                                         # Ha -> kcal/mol

# 1) Klasik kisim: HF + CASCI (referans)
lines = [l for l in open(fn).read().split("\n")[2:] if l.strip()]
mol = gto.M(atom="\n".join(lines[-16:]), basis="def2-svp", verbose=0)
mf = scf.RHF(mol).run()
cas = mcscf.CASCI(mf, NO, NE).run()
h1, ecore = cas.get_h1eff()
eri = ao2mo.restore(1, cas.get_h2eff(), NO)

# 2) Aktif uzay Hamiltonyeni -> kubit Hamiltonyeni (Jordan-Wigner)
two = np.swapaxes(eri, 1, 3)
Hf = qml.qchem.fermionic_observable(np.array([ecore]), h1, two)
H = qml.qchem.qubit_observable(Hf)
nq = 2 * NO
print(f"Aktif uzay: {NE} elektron, {NO} orbital -> {nq} kubit, {len(H.operands)} Pauli terimi")

# 3) Ansatz: HF durumu + tum tekli/ciftli uyarilmalar (UCCSD benzeri)
hf_state = qml.qchem.hf_state(NE, nq)
singles, doubles = qml.qchem.excitations(NE, nq)
print(f"Parametre sayisi: {len(singles)} tekli + {len(doubles)} ciftli = {len(singles)+len(doubles)}")
dev = qml.device("default.qubit", wires=nq)

@qml.qnode(dev, diff_method="backprop")
def energy(theta):
    qml.BasisState(hf_state, wires=range(nq))
    for i, d in enumerate(doubles):
        qml.DoubleExcitation(theta[i], wires=d)
    for j, s in enumerate(singles):
        qml.SingleExcitation(theta[len(doubles) + j], wires=s)
    return qml.expval(H)

grad = qml.grad(energy)
hist = []
def f(x):
    e = float(energy(pnp.array(x, requires_grad=True))); hist.append(e); return e
def g(x):
    return np.array(grad(pnp.array(x, requires_grad=True)))

t0 = time.time()
x0 = np.zeros(len(singles) + len(doubles))
e0 = f(x0)
res = minimize(f, x0, jac=g, method="BFGS", options={"gtol": 1e-8})
e_vqe = res.fun

print(f"\nHF (PySCF)            : {mf.e_tot:.8f} Ha")
print(f"HF (kubit, theta=0)   : {e0:.8f} Ha   <- ayni olmali (Hamiltonyen kontrolu)")
print(f"VQE                   : {e_vqe:.8f} Ha   ({res.nit} iterasyon, {time.time()-t0:.1f} s)")
print(f"CASCI (klasik kesin)  : {cas.e_tot:.8f} Ha")
print(f"VQE - CASCI           : {(e_vqe-cas.e_tot)*HA:.5f} kcal/mol  (kimyasal dogruluk < 1)")
print(f"Yakalanan korelasyon  : VQE {(e_vqe-mf.e_tot)*HA:.3f}  |  CASCI {(cas.e_tot-mf.e_tot)*HA:.3f} kcal/mol")
np.savetxt(fn.replace(".xyz", f"_vqe_{NE}e{NO}o.txt"), hist)
