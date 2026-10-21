# Porfin (serbest baz, C20H14N4): RDKit ile kur, GFN2-xTB ile optimize et
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
smi = "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2"
m = Chem.AddHs(Chem.MolFromSmiles(smi))
print("Formul:", Chem.rdMolDescriptors.CalcMolFormula(m))
AllChem.EmbedMolecule(m, randomSeed=3); AllChem.MMFFOptimizeMolecule(m)
at = Atoms([a.GetSymbol() for a in m.GetAtoms()], positions=m.GetConformer().GetPositions())
at.calc = TBLite(method="GFN2-xTB", verbosity=0)
BFGS(at, logfile=None).run(fmax=0.02, steps=500)
X = at.get_positions(); X -= X.mean(0)
_, _, vt = np.linalg.svd(X); X = X @ vt.T           # molekul duzlemi = xy, z = normal
print(f"Duzlemsellik (z RMS): {np.sqrt((X[:,2]**2).mean()):.3f} A")
with open("porfin.xyz", "w") as f:
    f.write(f"{len(at)}\nporfin GFN2-xTB\n")
    for s, x in zip(at.get_chemical_symbols(), X): f.write(f"{s:2s} {x[0]:12.6f} {x[1]:12.6f} {x[2]:12.6f}\n")
print("Kaydedildi: porfin.xyz")
