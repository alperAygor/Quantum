# Klasik referans: TD-DFT (B3LYP ve CAM-B3LYP) / def2-SVP
# Kullanim: python3 tddft_referans.py porfin      (porfin.xyz okunur, porfin_tddft.json yazilir)
# Singletler: TDDFT (RPA), 4 hal.  Tripletler: TDA, 2 hal (T1 icin Delta-SCF'den daha guvenilir).
import os, sys, time, json
os.makedirs(os.path.expanduser("~/pyscf_tmp"), exist_ok=True)
os.environ.setdefault("PYSCF_TMPDIR", os.path.expanduser("~/pyscf_tmp"))
import numpy as np
from pyscf import gto, dft

NAME = sys.argv[1]
BASIS = os.environ.get("BAZ", "def2-svp")
FONK = os.environ.get("FONK", "b3lyp,camb3lyp").split(",")
EV = 27.211386; t0 = time.time()
OUT = f"{NAME}_tddft.json"
R = json.load(open(OUT)) if os.path.exists(OUT) else {"molekul": NAME, "baz": BASIS}
def log(s): print(f"[{NAME}] {s}  [{time.time()-t0:.0f} s]", flush=True)

lines = [l for l in open(f"{NAME}.xyz").read().split("\n")[2:] if l.strip()]
ecp = BASIS if any(l.split()[0] == "Pd" for l in lines) else None
mol = gto.M(atom=f"{NAME}.xyz", basis=BASIS, ecp=ecp, verbose=0,
            max_memory=int(os.environ.get("QMEM", "6000")))

for xc in FONK:
    if xc in R:
        log(f"{xc} zaten var, atlaniyor"); continue
    ks = dft.RKS(mol, xc=xc).density_fit(); ks.grids.level = 2; ks.conv_tol = 1e-9; ks.run()
    log(f"{xc} SCF bitti (yakinsadi={ks.converged})")
    # Singletler (TDDFT/RPA)
    td = ks.TDDFT(); td.nstates = 4; td.conv_tol = 1e-6; td.kernel()
    eS = (np.array(td.e) * EV).tolist(); fS = np.array(td.oscillator_strength()).tolist()
    # Tripletler (TDA; RPA tripletlerinde kararsizlik sorunu var)
    tt = ks.TDA(); tt.singlet = False; tt.nstates = 2; tt.conv_tol = 1e-6; tt.kernel()
    eT = (np.array(tt.e) * EV).tolist()
    parlak = next((e for e, f in zip(eS, fS) if f > 0.1), None)
    R[xc] = {"S_eV": [round(x, 4) for x in eS], "f": [round(x, 4) for x in fS],
             "S1_eV": round(eS[0], 4), "f_S1": round(fS[0], 4), "ilk_parlak_eV": parlak,
             "T1_eV": round(eT[0], 4), "yakinsadi": bool(np.all(td.converged)) and bool(np.all(tt.converged))}
    json.dump(R, open(OUT, "w"), indent=1)
    log(f"{xc}: S1 {eS[0]:.3f} eV (f={fS[0]:.3f}), S2 {eS[1]:.3f} eV (f={fS[1]:.3f}), T1 {eT[0]:.3f} eV")

R["sure_s"] = round(time.time() - t0); json.dump(R, open(OUT, "w"), indent=1)
log("TD-DFT bitti -> " + OUT)
