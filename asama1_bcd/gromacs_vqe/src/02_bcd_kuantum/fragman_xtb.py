import sys, time
import numpy as np
from tblite.interface import Calculator
GRO = sys.argv[1] if len(sys.argv) > 1 else "view.gro"
ITP = sys.argv[2] if len(sys.argv) > 2 else "bcd_345755.itp"
STRIDE = int(sys.argv[3]) if len(sys.argv) > 3 else 90
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

def xtb(Z, X):
    c = Calculator("GFN2-xTB", Z, X * B); c.set("verbosity", 0)
    return c.singlepoint().get("energy")

def frames(fn):
    with open(fn) as f:
        while True:
            title = f.readline()
            if not title.strip(): break
            n = int(f.readline()); atoms = [f.readline() for _ in range(n)]; f.readline()
            yield title.strip(), atoms

# --- 2) her karede: tam kompleks E_int vs 7 parcanin toplami
res = []; t0 = time.time(); setup = None
for k, (title, atoms) in enumerate(frames(GRO)):
    if k % STRIDE: continue
    xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in atoms]) * 10.0
    if setup is None:
        resn = [a[5:10].strip() for a in atoms]; gname = [a[10:15].strip() for a in atoms]
        cd = [i for i, r in enumerate(resn) if r == "AOYB"]; gu = [i for i, r in enumerate(resn) if r == "G108"]
        assert len(cd) == N, "gro ve itp atom sayisi uyusmuyor"
        Zg = np.array([ZMAP[gname[i][0]] for i in gu]); Zc = np.array([ZMAP[e] for e in el]); setup = True
    Xc, Xg = xyz[cd], xyz[gu]
    eg = xtb(Zg, Xg)
    full = (xtb(np.concatenate([Zc, Zg]), np.vstack([Xc, Xg])) - xtb(Zc, Xc) - eg) * HA
    def ei(ks):
        Zf, Xf = capped(ks, Xc)
        return (xtb(np.concatenate([Zf, Zg]), np.vstack([Xf, Xg])) - xtb(Zf, Xf) - eg) * HA
    parts = [ei(f) for f in range(7)]
    dimer = [ei(p) for p in pairs]
    ovl = sum(dimer) - sum(parts)
    res.append([full] + parts + [ovl])
    print(f"kare {k:4d}: tam {full:7.2f} | tekli toplam {sum(parts):7.2f} | ikili-ortusmeli {ovl:7.2f}  [{time.time()-t0:.0f} s]", flush=True)

R = np.array(res)
print("\nOZET (kcal/mol, ortalama +- std)")
print(f"  Tam kompleks E_int     : {R[:,0].mean():7.2f} +- {R[:,0].std():.2f}")
S1 = R[:,1:8].sum(1); S2 = R[:,8]
print(f"  Sema 1 (7 tekli glukoz)       : {S1.mean():7.2f} +- {S1.std():.2f}   hata {(S1-R[:,0]).mean():+6.2f} +- {(S1-R[:,0]).std():.2f}")
print(f"  Sema 2 (ikili - tekli, MFCC)  : {S2.mean():7.2f} +- {S2.std():.2f}   hata {(S2-R[:,0]).mean():+6.2f} +- {(S2-R[:,0]).std():.2f}")
print("  Halka sirasi (glukoz no):", [o+1 for o in order])
for f in range(7):
    print(f"  Glukoz {f+1}: {R[:,f+1].mean():6.2f} +- {R[:,f+1].std():.2f}")
np.savetxt("fragman_xtb.txt", R, header="tam g1 g2 g3 g4 g5 g6 g7 sema2 (kcal/mol)")
