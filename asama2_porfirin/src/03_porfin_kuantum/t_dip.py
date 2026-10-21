import numpy as np
from pyscf import gto, scf, dft, mcscf, fci
EV=27.211386
mol = gto.M(atom="porfin.xyz", basis="def2-svp", verbose=0, max_memory=4000)
ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
mf = scf.RHF(mol).density_fit()
c = mcscf.CASCI(mf, 4, 4); c.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0)
c.fcisolver.nroots = 5; c.kernel(ks.mo_coeff)
mo = ks.mo_coeff[:, c.ncore:c.ncore+4]
d = np.einsum("xij,ip,jq->xpq", mol.intor("int1e_r"), mo, mo)
for i in range(1,5):
    t = c.fcisolver.trans_rdm1(c.ci[0], c.ci[i], 4, 4)
    mu = -np.einsum("xpq,pq->x", d, t); dE = c.e_tot[i]-c.e_tot[0]
    print(i, round(dE*EV,3), np.round(mu,3), "f=", round(2/3*dE*(mu@mu),4))
