import os  # CIKTI_DIR: figurlerin yazilacagi klasor (varsayilan: calisilan klasor)
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pyscf import gto, dft
from pyscf.tools import molden
d = molden.load("porfin_orbitaller.molden"); mol, en, C = d[0], d[1], d[2]
n = mol.nelectron // 2
g = np.linspace(-6.5, 6.5, 160); X, Y = np.meshgrid(g, g)
pts = np.c_[X.ravel(), Y.ravel(), np.full(X.size, 0.8)] / 0.52917721
ao = dft.numint.eval_ao(mol, pts)
xyz = mol.atom_coords() * 0.52917721; sym = [mol.atom_symbol(i) for i in range(mol.natm)]
fig, axs = plt.subplots(2, 2, figsize=(10, 10))
for ax, (k, nm) in zip(axs.ravel(), [(n-2, "HOMO-1"), (n-1, "HOMO"), (n, "LUMO"), (n+1, "LUMO+1")]):
    v = (ao @ C[:, k]).reshape(X.shape); m = np.abs(v).max()
    ax.contourf(X, Y, v, levels=np.linspace(-m, m, 21), cmap="RdBu_r")
    for i in range(mol.natm):
        for j in range(i):
            if np.linalg.norm(xyz[i]-xyz[j]) < 1.6 and "H" not in (sym[i], sym[j]):
                ax.plot(*zip(xyz[i][:2], xyz[j][:2]), "k-", lw=1)
    for s, c in zip(sym, xyz):
        if s == "N": ax.plot(c[0], c[1], "o", color="tab:blue", ms=7)
        if s == "H" and np.linalg.norm(c[:2]) < 1.5: ax.plot(c[0], c[1], "o", color="0.3", ms=4)
    ax.set_title(f"{nm}  (MO {k+1}, {en[k]*27.2114:.2f} eV)"); ax.set_aspect("equal"); ax.axis("off")
fig.suptitle("Porfin: Gouterman'in 4 orbitali (molekul duzleminin 0.8 A ustu)\nkirmizi/mavi = dalga fonksiyonunun isareti", fontsize=12)
plt.tight_layout(); plt.savefig(os.environ.get("CIKTI_DIR", ".") + "/porfin_gouterman_orbitaller.png", dpi=110)
print("ok")
