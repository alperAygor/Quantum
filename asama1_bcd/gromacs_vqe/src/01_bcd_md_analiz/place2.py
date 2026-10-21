import math
def read(f):
    L = open(f).read().splitlines(); n = int(L[1]); return L[2:2+n], L[2+n]
P = lambda l: [float(l[20:28]), float(l[28:36]), float(l[36:44])]
N = lambda l: l[10:15].strip()
sub  = lambda a,b: [a[i]-b[i] for i in range(3)]
add  = lambda a,b: [a[i]+b[i] for i in range(3)]
mul  = lambda a,s: [x*s for x in a]
dot  = lambda a,b: sum(a[i]*b[i] for i in range(3))
cross= lambda a,b: [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
unit = lambda a: mul(a, 1/math.sqrt(dot(a,a)))
mean = lambda ps: [sum(p[i] for p in ps)/len(ps) for i in range(3)]
def rot(v,k,t):
    c,s = math.cos(t), math.sin(t)
    return add(add(mul(v,c), mul(cross(k,v),s)), mul(k, dot(k,v)*(1-c)))

YON = +1   # +1: OH genis (ikincil) kenara, -1: OH dar (birincil) kenara

cd, box = read('cd_c.gro'); g, _ = read('guest_c.gro')
cdx = [P(l) for l in cd]; cdn = [N(l) for l in cd]
# molid 345755 (AOYB) atom adlari
prim = {'O35','O30','O25','O20','O15','O10','O5'}
sec  = {'O32','O31','O27','O26','O22','O21','O17','O16','O12','O11','O7','O6','O2','O1'}
assert sum(n in prim for n in cdn)==7 and sum(n in sec for n in cdn)==14, "CD atom adlari eslesmedi!"
mp = mean([x for x,n in zip(cdx,cdn) if n in prim])
ms = mean([x for x,n in zip(cdx,cdn) if n in sec])
u  = unit(sub(ms, mp))
c  = mean([x for x,n in zip(cdx,cdn) if not n.startswith('H')])

gx = [P(l) for l in g]; gn = [N(l) for l in g]
iO = [i for i,n in enumerate(gn) if n.startswith('O')][0]
iC = min((i for i,n in enumerate(gn) if n.startswith('C')), key=lambda i: math.dist(gx[i], gx[iO]))
gc = mean([x for x,n in zip(gx,gn) if not n.startswith('H')])
rel = [sub(x, gc) for x in gx]
a = unit(sub(gx[iO], gx[iC])); b = mul(u, YON)
k = cross(a, b); s = math.sqrt(dot(k,k)); ang = math.atan2(s, dot(a,b))
if s > 1e-6: k = unit(k); rel = [rot(r, k, ang) for r in rel]

best = None
for deg in range(0, 360, 10):
    rr = [rot(r, u, math.radians(deg)) for r in rel]
    for sh in [x*0.05 for x in range(-3, 4)]:
        pos = [add(add(c, mul(u, sh)), r) for r in rr]
        d = min(math.dist(p, q) for p in pos for q in cdx)
        if best is None or d > best[0]: best = (d, deg, sh, pos)
d, deg, sh, pos = best
print(f"Kenarlar arasi mesafe (dar-genis O ort.): {math.dist(mp, ms):.3f} nm")
print(f"Secilen konum: donme {deg} derece, eksen boyunca kayma {sh:+.2f} nm")
print(f"Misafir-CD en kisa atom mesafesi        : {d:.3f} nm")
print(f"Misafir O'nun eksen uzerindeki konumu   : {dot(sub(pos[iO], c), u):+.3f} nm  (+ = genis kenar)")
out = [l[:20] + '%8.3f%8.3f%8.3f' % tuple(p) + l[44:] for l, p in zip(g, pos)]
with open('complex.gro', 'w') as o:
    o.write('bCD-cyclopentanol complex\n%5d\n' % (len(cd) + len(g)))
    o.write('\n'.join(cd + out) + '\n' + box + '\n')
print("complex.gro yazildi:", len(cd) + len(g), "atom")
