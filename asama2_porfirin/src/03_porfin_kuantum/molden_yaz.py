from pyscf import gto, dft
from pyscf.tools import molden
mol = gto.M(atom="porfin.xyz", basis="def2-svp", verbose=0, max_memory=3000)
ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
n = mol.nelectron // 2
molden.from_scf(ks, "porfin_orbitaller.molden")
print("HOMO no:", n, " LUMO no:", n + 1, " (Avogadro'da 1'den sayilir)")
for k, nm in [(n-2, "HOMO-1"), (n-1, "HOMO"), (n, "LUMO"), (n+1, "LUMO+1")]:
    print(f"  MO {k+1}: {nm:7s} {ks.mo_energy[k]*27.2114:7.2f} eV")
