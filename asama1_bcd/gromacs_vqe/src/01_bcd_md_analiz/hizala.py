"""Kompleksi CD ekseni (dar->genis kenar) +z olacak sekilde dondurur,
CD agir atom merkezini kutunun (Lx/2, Ly/2, Lz/2) noktasina koyar."""
import math, sys
fin, fout = sys.argv[1], sys.argv[2]
Lx, Ly, Lz = (float(v) for v in sys.argv[3:6])

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

at, bd = read_itp('bcd_345755.itp'); n = len(at)
nb = {i: set() for i in range(n)}
for i, j in bd: nb[i].add(j); nb[j].add(i)
nH = lambda i: sum(at[j][1] == 'HOpt' for j in nb[i])
OH = [i for i in range(n) if at[i][1] == 'OAlc']
prim = [o for o in OH if any(at[c][1] == 'COpt' and nH(c) == 2 for c in nb[o])]
sec = [o for o in OH if o not in prim]
assert len(prim) == 7 and len(sec) == 14

L = open(fin).read().splitlines(); N = int(L[1]); A = L[2:2 + N]
X = [[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in A]
mean = lambda ps: [sum(p[k] for p in ps) / len(ps) for k in range(3)]
mp, ms = mean([X[i] for i in prim]), mean([X[i] for i in sec])
u = [ms[k] - mp[k] for k in range(3)]; l = math.sqrt(sum(v * v for v in u)); u = [v / l for v in u]
c = mean([X[i] for i in range(n) if not at[i][1].startswith('H')])
# u -> z donusu (Rodrigues)
z = [0, 0, 1]; k = [u[1], -u[0], 0.0]; s = math.sqrt(k[0]**2 + k[1]**2); cth = u[2]
def rot(v):
    if s < 1e-9: return v if cth > 0 else [v[0], -v[1], -v[2]]
    kk = [k[0] / s, k[1] / s, 0.0]; th = math.atan2(s, cth)
    cr = [kk[1]*v[2]-kk[2]*v[1], kk[2]*v[0]-kk[0]*v[2], kk[0]*v[1]-kk[1]*v[0]]
    d = sum(kk[i]*v[i] for i in range(3))
    return [v[i]*math.cos(th) + cr[i]*math.sin(th) + kk[i]*d*(1-math.cos(th)) for i in range(3)]
Y = [rot([x[i] - c[i] for i in range(3)]) for x in X]
Y = [[y[0] + Lx/2, y[1] + Ly/2, y[2] + Lz/2] for y in Y]
gz = mean(Y[n:]); 
with open(fout, 'w') as f:
    f.write('CD ekseni +z, PMF icin\n%5d\n' % N)
    for a, y in zip(A, Y): f.write(a[:20] + '%8.3f%8.3f%8.3f\n' % tuple(y))
    f.write('%10.5f%10.5f%10.5f\n' % (Lx, Ly, Lz))
zs = [y[2] for y in Y[:n]]
print(f"CD z araligi: {min(zs):.2f} - {max(zs):.2f} nm (merkez {Lz/2:.2f})")
print(f"Misafir merkezi: z = {gz[2]:.3f} nm  (CD merkezine gore {gz[2]-Lz/2:+.3f} nm)")
print("Yazildi:", fout)
