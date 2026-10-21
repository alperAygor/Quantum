# Molekul seti -> GFN2-xTB geometrileri (xyz)
import sys, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
P = "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2"
SET = {
 "porfin":          ("smi", P),
 "klorin":          ("smi", "C1CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2"),
 "bakteriyoklorin": ("smi", "C1CC2=CC3=NC(=CC4=CC=C(N4)C=C5CCC(=N5)C=C1N2)C=C3"),
 "tetrametilporfin":("smi", "CC1=C2C=CC(=C(C)C3=CC=C(N3)C(C)=C4C=CC(=N4)C(C)=C5C=CC1=N5)N2"),
 "porfirazin":      ("smi", "C1=CC2=NC3=CC=C(N3)N=C4C=CC(=N4)N=C5C=CC(=N5)N=C1N2"),
 "Zn-porfin":       ("metal", "Zn"),
 "Mg-porfin":       ("metal", "Mg"),
 "Pd-porfin":       ("metal", "Pd"),
}
def opt(sym, X, charge=0):
    at = Atoms(sym, positions=X); at.calc = TBLite(method="GFN2-xTB", charge=charge, verbosity=0)
    BFGS(at, logfile=None).run(fmax=0.02, steps=800)
    X = at.get_positions(); X -= X.mean(0); _, _, vt = np.linalg.svd(X); return at.get_chemical_symbols(), X @ vt.T
def from_smi(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=3); AllChem.MMFFOptimizeMolecule(m)
    return [a.GetSymbol() for a in m.GetAtoms()], m.GetConformer().GetPositions()
def write(name, sym, X):
    with open(f"{name}.xyz", "w") as f:
        f.write(f"{len(sym)}\n{name}\n")
        for s, x in zip(sym, X): f.write(f"{s:2s} {x[0]:12.6f} {x[1]:12.6f} {x[2]:12.6f}\n")
names = sys.argv[1:] or list(SET)
base = None
for n in names:
    kind, val = SET[n]
    if kind == "smi":
        sym, X = opt(*from_smi(val))
    else:
        if base is None: base = opt(*from_smi(P))
        sym0, X0 = base
        # ic NH hidrojenlerini cikar, metali merkeze koy
        Nidx = [i for i, s in enumerate(sym0) if s == "N"]
        drop = [i for i, s in enumerate(sym0) if s == "H" and min(np.linalg.norm(X0[i]-X0[j]) for j in Nidx) < 1.15]
        keep = [i for i in range(len(sym0)) if i not in drop]
        sym, X = opt([sym0[i] for i in keep] + [val], np.vstack([X0[keep], X0[Nidx].mean(0)]))
    write(n, sym, X)
    print(f"{n:17s} {len(sym):3d} atom  duzlemsellik {np.sqrt((X[:,2]**2).mean()):.3f} A", flush=True)
