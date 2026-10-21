"""beta-CD glukoz 'flip' kontrolu (Henriksen & Gilson 2017, JCTC 13, 4253'teki uyari).
Her karede: glikozidik O duzlemine gore her glukozun C6'si hangi tarafta?
Normal beta-CD: 7/7 ayni tarafta. Ters donmus glukoz -> kavite bozulur."""
import math, sys

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
nb = {i: set() for i in range(n)}
for i, j in bd: nb[i].add(j); nb[j].add(i)
nH = lambda i: sum(at[j][1] == 'HOpt' for j in nb[i])
C6 = [i for i in range(n) if at[i][1] == 'COpt' and nH(i) == 2]
C5 = [next(j for j in nb[c] if at[j][1] == 'COpt') for c in C6]
gly = [i for i in range(n) if at[i][1] == 'OE' and not any(j in C5 for j in nb[i])]
assert len(C6) == 7 and len(gly) == 7

sub = lambda a, b: [a[k] - b[k] for k in range(3)]
dot = lambda a, b: sum(a[k] * b[k] for k in range(3))
crs = lambda a, b: [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
mean = lambda ps: [sum(p[k] for p in ps) / len(ps) for k in range(3)]

kotu, eliptik, toplam, ornek = 0, [], 0, []
for t, x in frames(sys.argv[1] if len(sys.argv) > 1 else 'view.gro'):
    G = [x[i] for i in gly]; c = mean(G)
    nrm = [0, 0, 0]                                   # duzlem normali: ardisik O'larin capraz carpimlari
    for k in range(7):
        v = crs(sub(G[k], c), sub(G[(k + 1) % 7], c)); nrm = [nrm[m] + v[m] for m in range(3)]
    side = [dot(sub(x[c6], c), nrm) for c6 in C6]
    cogunluk = 1 if sum(s > 0 for s in side) >= 4 else -1
    ters = sum(1 for s in side if s * cogunluk < 0)
    r = [math.sqrt(dot(sub(g, c), sub(g, c))) for g in G]
    eliptik.append(max(r) / min(r)); toplam += 1
    if ters: kotu += 1; ornek.append((t, ters))
print(f"Kare: {toplam}")
print(f"En az bir glukozu ters donmus kare: {kotu}  (%{100*kotu/toplam:.1f})")
if ornek: print("  ornekler (t ps, ters glukoz sayisi):", ornek[:10])
eliptik.sort()
print(f"Halka eliptikligi (max/min glikozidik O yaricapi): medyan {eliptik[len(eliptik)//2]:.2f}, en kotu {eliptik[-1]:.2f}  (ideal 1.00)")
