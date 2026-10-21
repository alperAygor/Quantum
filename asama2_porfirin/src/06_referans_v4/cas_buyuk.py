# Buyuk aktif uzay kontrolu: CAS(8e,8o) = 16 kubit esdegeri, B3LYP orbitalleri + SC-NEVPT2
# Kullanim: python3 cas_buyuk.py bakteriyoklorin        (varsayilan 8 elektron, 8 orbital)
#           python3 cas_buyuk.py porfin 8 8
# CASCI(8,8) sonucu, ideal (gurultusuz) 16 kubitlik VQE'nin ulasacagi limittir.
import os, sys, time, json
os.makedirs(os.path.expanduser("~/pyscf_tmp"), exist_ok=True)
os.environ.setdefault("PYSCF_TMPDIR", os.path.expanduser("~/pyscf_tmp"))
import numpy as np
from pyscf import gto, scf, dft, mcscf, fci, mrpt
from pyscf.fci import cistring

NAME = sys.argv[1]; NE = int(sys.argv[2]) if len(sys.argv) > 2 else 8; NO = int(sys.argv[3]) if len(sys.argv) > 3 else 8
BASIS = os.environ.get("BAZ", "def2-svp")
EV = 27.211386; t0 = time.time(); R = {"molekul": NAME, "aktif_uzay": f"({NE}e,{NO}o)", "kubit": 2 * NO}
def log(s): print(f"[{NAME} CAS({NE},{NO})] {s}  [{time.time()-t0:.0f} s]", flush=True)

lines = [l for l in open(f"{NAME}.xyz").read().split("\n")[2:] if l.strip()]
ecp = BASIS if any(l.split()[0] == "Pd" for l in lines) else None
mol = gto.M(atom=f"{NAME}.xyz", basis=BASIS, ecp=ecp, verbose=0, max_memory=int(os.environ.get("QMEM", "6000")))

ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run(); log("B3LYP bitti")
mf = scf.RHF(mol).density_fit(); mo = ks.mo_coeff

cS = mcscf.CASCI(mf, NO, NE); cS.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0)
cS.fcisolver.nroots = 3; cS.kernel(mo); eS = np.array(cS.e_tot)
log("CASCI bitti")

def cfg(ci):
    sa = cistring.make_strings(range(NO), NE // 2); occ = {}
    for i, a in enumerate(sa):
        for j, b in enumerate(sa):
            k = tuple(((a >> q) & 1) + ((b >> q) & 1) for q in range(NO)); occ[k] = occ.get(k, 0) + ci[i, j] ** 2
    return sorted(occ.values(), reverse=True)
R["S1_konfig_agirlik"] = [round(float(x), 3) for x in cfg(cS.ci[1])[:3]]
R["S0_HF_agirlik"] = round(float(cfg(cS.ci[0])[0]), 3)

# Gecis dipolu ve osilator gucu (aktif uzay icinden)
mo_act = mo[:, cS.ncore:cS.ncore + NO]
dip = np.einsum("xij,ip,jq->xpq", mol.intor("int1e_r"), mo_act, mo_act)
for i in (1, 2):
    t = cS.fcisolver.trans_rdm1(cS.ci[0], cS.ci[i], NO, NE)
    mu = -np.einsum("xpq,pq->x", dip, t)
    R[f"f_S{i}"] = round(float(2 / 3 * (eS[i] - eS[0]) * mu @ mu), 4)
R["S1_CAS_eV"] = round(float((eS[1] - eS[0]) * EV), 4); R["S2_CAS_eV"] = round(float((eS[2] - eS[0]) * EV), 4)
json.dump(R, open(f"{NAME}_cas{NE}{NO}.json", "w"), indent=1)

n0 = mrpt.NEVPT(cS, root=0).kernel(); log("NEVPT2 S0 bitti")
n1 = mrpt.NEVPT(cS, root=1).kernel(); log("NEVPT2 S1 bitti")
R["S1_NEVPT2_eV"] = round(float((eS[1] + n1 - eS[0] - n0) * EV), 4)
mix = R["S1_konfig_agirlik"][1]
R["triyaj"] = "KUANTUM GEREKLI" if mix > 0.2 else "klasik yeterli"
R["sure_s"] = round(time.time() - t0)
json.dump(R, open(f"{NAME}_cas{NE}{NO}.json", "w"), indent=1)
log("SONUC: " + json.dumps(R))
