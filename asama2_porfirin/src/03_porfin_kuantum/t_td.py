import time
from pyscf import gto, dft, tddft
t0=time.time()
mol = gto.M(atom="porfin.xyz", basis="def2-svp", verbose=0, max_memory=4000)
ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
print("ks", time.time()-t0, flush=True)
td = tddft.TDA(ks); td.nstates = 4; td.verbose=4; td.kernel(); print("tda", time.time()-t0, td.e*27.2114, flush=True)
