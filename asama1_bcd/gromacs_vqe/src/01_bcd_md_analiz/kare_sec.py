import math
def read_itp(f):
    atoms, bonds, sec = [], [], None
    for line in open(f):
        s = line.split(';')[0].strip()
        if not s or s.startswith('#'): continue
        if s.startswith('['): sec = s.strip('[] ').lower(); continue
        p = s.split()
        if sec == 'atoms': atoms.append((p[4], p[1]))
        elif sec == 'bonds': bonds.append((int(p[0]) - 1, int(p[1]) - 1))
    return atoms, bonds
def frames(f):
    L = open(f).read().splitlines(); i = 0
    while i < len(L):
        t = float(L[i].split('t=')[1].split()[0]) if 't=' in L[i] else 0.0
        n = int(L[i + 1])
        yield t, [[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in L[i + 2:i + 2 + n]]
        i += n + 3
at, bd = read_itp('bcd_345755.itp'); n = len(at)
gat, _ = read_itp('cyclopentanol_aa.itp')
nb = {i: set() for i in range(n)}
for i, j in bd: nb[i].add(j); nb[j].add(i)
nH = lambda i: sum(at[j][1] == 'HOpt' for j in nb[i])
C6 = [i for i in range(n) if at[i][1] == 'COpt' and nH(i) == 2]
C5 = [next(j for j in nb[c] if at[j][1] == 'COpt') for c in C6]
gly = [i for i in range(n) if at[i][1] == 'OE' and not any(j in C5 for j in nb[i])]
heavy_cd = [i for i in range(n) if not at[i][1].startswith('H')]
heavy_g = [n + i for i, (nm, tp) in enumerate(gat) if not tp.startswith('H')]
assert len(C6) == 7 and len(gly) == 7
sub = lambda a, b: [a[k] - b[k] for k in range(3)]
dot = lambda a, b: sum(a[k] * b[k] for k in range(3))
crs = lambda a, b: [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
mean = lambda ps: [sum(p[k] for p in ps) / len(ps) for k in range(3)]
rows = []
for t, x in frames('view.gro'):
    G = [x[i] for i in gly]; c = mean(G)
    nrm = [0, 0, 0]
    for k in range(7):
        v = crs(sub(G[k], c), sub(G[(k + 1) % 7], c)); nrm = [nrm[m] + v[m] for m in range(3)]
    side = [dot(sub(x[c6], c), nrm) for c6 in C6]
    ters = min(sum(s > 0 for s in side), sum(s < 0 for s in side))
    r = [math.dist(g, c) for g in G]
    d = math.dist(mean([x[i] for i in heavy_g]), mean([x[i] for i in heavy_cd]))
    rows.append((t, max(r) / min(r), d, ters))
med = sorted(e for _, e, _, _ in rows)[len(rows) // 2]
aday = [r for r in rows if r[3] == 0]
en = min(aday, key=lambda r: abs(r[1] - med) / 0.05 + r[2] / 0.05)
print(f"Medyan eliptiklik: {med:.2f}")
print(f"Secilen kare: t = {en[0]:.0f} ps | eliptiklik {en[1]:.2f} | misafir-CD mesafesi {en[2]:.3f} nm | ters glukoz {en[3]}")
