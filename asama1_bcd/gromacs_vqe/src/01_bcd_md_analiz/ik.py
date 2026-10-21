import sys
from rdkit import Chem
for f in sys.argv[1:]:
    m=Chem.MolFromPDBFile(f, removeHs=False, proximityBonding=False)
    Chem.AssignStereochemistryFrom3D(m)
    print(f.split('/')[-1][9:], Chem.MolToInchiKey(m))
