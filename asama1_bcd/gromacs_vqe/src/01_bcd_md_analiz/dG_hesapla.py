"""WHAM profilinden standart baglanma serbest enerjisi.
K = pi R^2 * integral_site exp(-W/kT) dxi * 0.6022 [1/M] ; dG = -RT ln K"""
import math, sys
f = sys.argv[1] if len(sys.argv) > 1 else 'profil.xvg'
T, R = 300.0, 0.4                       # K, silindir yaricapi (nm)
kT = 0.0083145 * T                      # kJ/mol
xs, ws = [], []
for l in open(f):
    if l[0] in '#@' or not l.strip(): continue
    a, b = l.split()[:2]; xs.append(float(a)); ws.append(float(b))
# su (bulk) referansi: elde olan uc(lar)daki duz bolge
lo = [w for x, w in zip(xs, ws) if x <= -1.5]
hi = [w for x, w in zip(xs, ws) if x >= 1.5]
ref_parts = [p for p in (lo, hi) if p]
bulk = sum(sum(p) / len(p) for p in ref_parts) / len(ref_parts)
W = [w - bulk for w in ws]
i0 = min(range(len(W)), key=lambda i: W[i])
print(f"Bulk referansi: {'dar kenar' if lo else ''}{' + ' if lo and hi else ''}{'genis kenar' if hi else ''} tarafi")
if lo and hi: print(f"  iki ucun farki: {sum(hi)/len(hi) - sum(lo)/len(lo):+.2f} kJ/mol (0'a yakin olmali)")
print(f"Kuyu dibi: xi = {xs[i0]:+.2f} nm, W = {W[i0]:.2f} kJ/mol = {W[i0]/4.184:.2f} kcal/mol")
for xc in (0.8, 1.0, 1.2):
    s = 0.0
    for k in range(len(xs) - 1):
        xm = 0.5 * (xs[k] + xs[k + 1])
        if abs(xm) <= xc:
            s += 0.5 * (math.exp(-W[k] / kT) + math.exp(-W[k + 1] / kT)) * (xs[k + 1] - xs[k])
    K = math.pi * R * R * s * 0.6022
    dG = -kT * math.log(K)
    print(f"  kuyu sinirı |xi|<={xc:.1f} nm:  K = {K:7.1f} 1/M   dG = {dG:6.2f} kJ/mol = {dG/4.184:5.2f} kcal/mol")
print("Deney: dG = -12.76 kJ/mol = -3.05 kcal/mol (K ~ 172 1/M)")
