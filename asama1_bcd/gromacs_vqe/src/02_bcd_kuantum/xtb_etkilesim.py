import sys, time
import numpy as np
from tblite.interface import Calculator
from dftd4.interface import DampingParam, DispersionModel
GRO = sys.argv[1] if len(sys.argv) > 1 else "view.gro"
STRIDE = int(sys.argv[2]) if len(sys.argv) > 2 else 30
PRIM = {"O35","O30","O25","O20","O15","O10","O5"}
GLYC = [21,42,63,84,105,126,138]
ZMAP = {"H": 1, "C": 6, "O": 8}
B, HA = 1/0.52917721, 627.5095
D4P = DampingParam(s6=1.0, s8=2.7, a1=0.52, a2=5.0, s9=5.0)

def frames(fn):
    with open(fn) as f:
        while True:
            title = f.readline()
            if not title.strip(): break
            n = int(f.readline())
            atoms = [f.readline() for _ in range(n)]
            f.readline()
            yield title.strip(), atoms

def xtb(Z, X, solv):
    c = Calculator("GFN2-xTB", Z, X * B); c.set("verbosity", 0)
    if solv: c.add("alpb-solvation", "water")
    return c.singlepoint().get("energy")

def disp(Z, X):
    return DispersionModel(Z, X * B).get_dispersion(D4P, grad=False)["energy"]

def eint(fun, Zc, Xc, Zh, Xh, Zg, Xg, *a):
    return (fun(Zc, Xc, *a) - fun(Zh, Xh, *a) - fun(Zg, Xg, *a)) * HA

rows = []; setup = None; t0 = time.time()
print(f"{'kare':>5} {'t(ps)':>8} {'mod':>4} {'cos':>6} {'E_gaz':>8} {'E_su':>8} {'E_disp':>8}  (kcal/mol)")
for k, (title, atoms) in enumerate(frames(GRO)):
    if k % STRIDE: continue
    xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in atoms]) * 10.0
    if setup is None:
        resn = [a[5:10].strip() for a in atoms]; name = [a[10:15].strip() for a in atoms]
        cd = [i for i, r in enumerate(resn) if r == "AOYB"]; gu = [i for i, r in enumerate(resn) if r == "G108"]
        prim = [i for i in cd if name[i] in PRIM]; gly = [cd[j-1] for j in GLYC]
        iO = [i for i in gu if name[i] == "O1"][0]; iC = [i for i in gu if name[i] == "C4"][0]
        Z = np.array([ZMAP[n[0]] for n in name]); keep = cd + gu; setup = True
    c = xyz[gly].mean(0); _, _, vt = np.linalg.svd(xyz[gly] - c); n = vt[2]
    if np.dot(xyz[prim].mean(0) - c, n) < 0: n = -n
    v = xyz[iO] - xyz[iC]; cos = float(np.dot(v, n) / np.linalg.norm(v))
    mod = "A" if cos > 0.5 else ("B" if cos < -0.5 else "ara")
    args = (Z[keep], xyz[keep], Z[cd], xyz[cd], Z[gu], xyz[gu])
    try:
        eg = eint(xtb, *args, False); ew = eint(xtb, *args, True); ed = eint(disp, *args)
    except Exception as ex:
        print(f"{k:5d}  hata: {ex}"); continue
    t = title.split("t=")[-1].split()[0] if "t=" in title else str(k)
    rows.append((k, float(t), mod, cos, eg, ew, ed))
    print(f"{k:5d} {float(t):8.0f} {mod:>4} {cos:+6.2f} {eg:8.2f} {ew:8.2f} {ed:8.2f}   [{time.time()-t0:.0f} s]", flush=True)

print("\nOZET (ortalama +- std, kcal/mol)")
for m in ["A", "B", "ara", "tum"]:
    s = [r for r in rows if r[2] == m or m == "tum"]
    if not s: continue
    a = np.array([r[4:] for r in s])
    print(f"  {m:>3} (n={len(s):2d}): gaz {a[:,0].mean():7.2f} +- {a[:,0].std():4.2f} | su(ALPB) {a[:,1].mean():7.2f} +- {a[:,1].std():4.2f}"
          f" | dispersiyon {a[:,2].mean():7.2f} +- {a[:,2].std():4.2f}")
with open("xtb_etkilesim.csv", "w") as f:
    f.write("kare,t_ps,mod,cos,E_gaz,E_su,E_disp\n")
    for r in rows: f.write(",".join(str(x) for x in r) + "\n")
print("Kaydedildi: xtb_etkilesim.csv")
