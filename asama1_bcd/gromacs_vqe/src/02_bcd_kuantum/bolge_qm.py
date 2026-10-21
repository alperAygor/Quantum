import sys, time
import numpy as np
from tblite.interface import Calculator
from pyscf import gto, scf, dft
from pyscf.mp.dfmp2 import DFMP2
from dftd3.interface import RationalDampingParam, DispersionModel
XYZ = sys.argv[1] if len(sys.argv) > 1 else "mod_A.xyz"
ITP = sys.argv[2] if len(sys.argv) > 2 else "bcd_345755.itp"
BASIS = sys.argv[3] if len(sys.argv) > 3 else "def2-svp"
ZMAP = {"H": 1, "C": 6, "O": 8}
B, HA = 1/0.52917721, 627.5095

# --- 1) itp'den bag grafigi -> 7 glukoz birimi + kesilen glikozidik baglar
sec = None; names = []; bonds = []
for line in open(ITP):
    s = line.split(";")[0].strip()
    if not s or s.startswith("#"): continue
    if s.startswith("["): sec = s.strip("[] ").lower(); continue
    p = s.split()
    if sec == "atoms": names.append(p[4])
    elif sec == "bonds": bonds.append((int(p[0]) - 1, int(p[1]) - 1))
el = [n[0] for n in names]; N = len(names)
nb = [set() for _ in range(N)]
for a, b in bonds: nb[a].add(b); nb[b].add(a)
nO = lambda i: sum(el[j] == "O" for j in nb[i])
cand = set()
for o in range(N):
    if el[o] == "O" and sum(el[j] == "C" for j in nb[o]) == 2:
        for c in nb[o]:
            if el[c] == "C" and nO(c) == 2: cand.add((min(o, c), max(o, c)))
kept = [bd for bd in bonds if (min(bd), max(bd)) not in cand]
par = list(range(N))
def root(i):
    while par[i] != i: par[i] = par[par[i]]; i = par[i]
    return i
for a, b in kept: par[root(a)] = root(b)
roots = sorted({root(i) for i in range(N)}, key=lambda r: min(i for i in range(N) if root(i) == r))
frag = [[i for i in range(N) if root(i) == r] for r in roots]
fid = {i: k for k, f in enumerate(frag) for i in f}
cuts = [(a, b) for a, b in cand if fid[a] != fid[b]]
print(f"CD atom: {N}, birim sayisi: {len(frag)} (atom: {[len(f) for f in frag]}), kesilen glikozidik bag: {len(cuts)}")
assert len(frag) == 7 and len(cuts) == 7, "Parcalama beklenen gibi degil!"

adj = {k: set() for k in range(7)}
for a, b in cuts: adj[fid[a]].add(fid[b]); adj[fid[b]].add(fid[a])
order = [0]
while len(order) < 7: order.append([j for j in adj[order[-1]] if j not in order][0])
pairs = [(order[i], order[(i + 1) % 7]) for i in range(7)]
def capped(ks, X):
    ks = set(ks) if isinstance(ks, (list, tuple, set)) else {ks}
    idx = [i for k in ks for i in frag[k]]; Z = [ZMAP[el[i]] for i in idx]; P = [X[i] for i in idx]
    for a, b in cuts:
        for u, v in ((a, b), (b, a)):
            if fid[u] in ks and fid[v] not in ks:
                d = 1.09 if el[u] == "C" else 0.96
                w = X[v] - X[u]; P.append(X[u] + d * w / np.linalg.norm(w)); Z.append(1)
    return np.array(Z), np.array(P)


SYM = {1: "H", 6: "C", 8: "O"}
def xtb(Z, X):
    c = Calculator("GFN2-xTB", Z, X * B); c.set("verbosity", 0)
    return c.singlepoint().get("energy")

# --- 1) temsilci kareyi oku: ilk N atom CD (itp sirasi), son 16 atom misafir
L = [l.split() for l in open(XYZ).read().split("\n")[2:] if l.strip()]
X = np.array([[float(v) for v in l[1:4]] for l in L])
Xc, Xg = X[:N], X[N:]
Zg = np.array([ZMAP[l[0]] for l in L[N:]])
print(f"Kare: {XYZ}  (CD {len(Xc)} atom, misafir {len(Xg)} atom)")

# --- 2) xTB ile 7 bolgeyi tara, en guclu etkilesen glukozu sec
eg = xtb(Zg, Xg); scan = []
for k in range(7):
    Zf, Xf = capped(k, Xc)
    scan.append((xtb(np.concatenate([Zf, Zg]), np.vstack([Xf, Xg])) - xtb(Zf, Xf) - eg) * HA)
print("xTB bolge taramasi (kcal/mol):", " ".join(f"G{k+1}:{e:5.2f}" for k, e in enumerate(scan)))
K = int(np.argmin(scan)); Zf, Xf = capped(K, Xc)
print(f"Secilen bolge: misafir + glukoz {K+1}  ({len(Zf)} + {len(Zg)} = {len(Zf)+len(Zg)} atom)")
out = XYZ.replace(".xyz", f"_bolge_G{K+1}.xyz")
with open(out, "w") as f:
    f.write(f"{len(Zf)+len(Zg)}\nglukoz {K+1}: ilk {len(Zf)} atom | misafir: son {len(Zg)} atom\n")
    for z, x in zip(np.concatenate([Zf, Zg]), np.vstack([Xf, Xg])):
        f.write(f"{SYM[z]:2s} {x[0]:12.5f} {x[1]:12.5f} {x[2]:12.5f}\n")
print("Bolge geometrisi:", out)

# --- 3) counterpoise duzeltmeli etkilesim enerjileri (PySCF)
ZA, XA, ZB, XB = Zf, Xf, Zg, Xg
def mol(ghostA, ghostB):
    at = [("ghost-" if ghostA else "") + SYM[z] + " %f %f %f" % tuple(x) for z, x in zip(ZA, XA)]
    at += [("ghost-" if ghostB else "") + SYM[z] + " %f %f %f" % tuple(x) for z, x in zip(ZB, XB)]
    return gto.M(atom="; ".join(at), basis=BASIS, verbose=0, max_memory=4000)
res = {}; t0 = time.time()
for tag, gA, gB in [("AB", 0, 0), ("A", 0, 1), ("B", 1, 0)]:
    m = mol(gA, gB)
    hf = scf.RHF(m).density_fit().run()
    pt = DFMP2(hf); pt.kernel(with_t2=False)
    ks = dft.RKS(m, xc="b3lyp").density_fit().run()
    res[tag] = (hf.e_tot, pt.e_tot, ks.e_tot)
    print(f"  {tag:2s}: {m.nao} baz fonk., {m.nelectron} elektron  [{time.time()-t0:.0f} s]", flush=True)
def d3(Z, Xa):
    return DispersionModel(Z, Xa * B).get_dispersion(RationalDampingParam(method="b3lyp"), grad=False)["energy"]
ed3 = (d3(np.concatenate([ZA, ZB]), np.vstack([XA, XB])) - d3(ZA, XA) - d3(ZB, XB)) * HA
e = {i: (res["AB"][i] - res["A"][i] - res["B"][i]) * HA for i in range(3)}
print(f"\nETKILESIM ENERJISI, misafir + glukoz {K+1}  ({BASIS}, counterpoise; kcal/mol)")
print(f"  HF            : {e[0]:7.2f}   (korelasyon/dispersiyon YOK)")
print(f"  MP2           : {e[1]:7.2f}")
print(f"  B3LYP         : {e[2]:7.2f}")
print(f"  B3LYP-D3(BJ)  : {e[2]+ed3:7.2f}   (D3 katkisi {ed3:.2f})")
print(f"  GFN2-xTB      : {scan[K]:7.2f}")
print(f"  MP2 - HF      : {e[1]-e[0]:7.2f}   <- elektron korelasyonunun bagladigi kisim (VQE hedefi)")
