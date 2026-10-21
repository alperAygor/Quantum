"""beta-CD (AOYB) + siklopentanol (G108) yorunge analizi.
Girdi : view.gro (trjconv ile suyu cikarilmis, CD merkezlenmis cok kareli gro)
        bcd_345755.itp, cyclopentanol_aa.itp (atom tipleri ve baglar buradan okunur)
Cikti : ekrana ozet + analiz.csv (kare kare degerler)"""
import math, sys

def read_itp(f):
    atoms, bonds, sec = [], [], None
    for line in open(f):
        s = line.split(';')[0].strip()
        if not s: continue
        if s.startswith('['): sec = s.strip('[] ').lower(); continue
        p = s.split()
        if sec == 'atoms': atoms.append((p[4], p[1]))          # (ad, tip)
        elif sec == 'bonds': bonds.append((int(p[0]) - 1, int(p[1]) - 1))
    return atoms, bonds

def frames(f):
    L = open(f).read().splitlines(); i = 0
    while i < len(L):
        t = float(L[i].split('t=')[1].split()[0]) if 't=' in L[i] else None
        n = int(L[i + 1]); at = L[i + 2:i + 2 + n]
        yield t, [(a[5:10].strip(), a[10:15].strip(),
                   [float(a[20:28]), float(a[28:36]), float(a[36:44])]) for a in at]
        i += n + 3

sub = lambda a, b: [a[k] - b[k] for k in range(3)]
dot = lambda a, b: sum(a[k] * b[k] for k in range(3))
nrm = lambda a: math.sqrt(dot(a, a))
mean = lambda ps: [sum(p[k] for p in ps) / len(ps) for k in range(3)]
def angle(a, b, c):  # a-b-c acisi (derece)
    u, v = sub(a, b), sub(c, b)
    return math.degrees(math.acos(max(-1, min(1, dot(u, v) / (nrm(u) * nrm(v))))))

cd_at, cd_b = read_itp('bcd_345755.itp')
g_at, g_b = read_itp('cyclopentanol_aa.itp')
ncd = len(cd_at)

def hydroxyls(atoms, bonds, off=0):   # (O indeksi, H indeksi) ciftleri
    pairs = []
    for i, j in bonds:
        for o, h in ((i, j), (j, i)):
            if atoms[o][1].startswith('O') and atoms[h][1].startswith('H'):
                pairs.append((o + off, h + off))
    return pairs

cd_OH = hydroxyls(cd_at, cd_b)                     # CD hidroksilleri (verici)
cd_acc = [i for i, (n, t) in enumerate(cd_at) if t.startswith('O')]   # tum CD oksijenleri (alici)
g_OH = hydroxyls(g_at, g_b, ncd)                   # misafir OH
gO = g_OH[0][0]
# birincil O: CH2 karbonuna (2 HOpt) bagli OAlc ; ikincil: diger OAlc
nH = {i: 0 for i in range(ncd)}
for i, j in cd_b:
    if cd_at[j][1] == 'HOpt': nH[i] += 1
    if cd_at[i][1] == 'HOpt': nH[j] += 1
prim, sec = [], []
for o, h in cd_OH:
    c = [j for a, b in cd_b for j in ((b,) if a == o else (a,) if b == o else ()) if cd_at[j][1] == 'COpt'][0]
    (prim if nH[c] == 2 else sec).append(o)
assert len(prim) == 7 and len(sec) == 14, (len(prim), len(sec))

DHA_R, DHA_ANG = 0.35, 30.0     # H-bagi olcutu: O...O <= 0.35 nm, H-O...O acisi <= 30 derece
rows = []
for t, at in frames(sys.argv[1] if len(sys.argv) > 1 else 'view.gro'):
    x = [a[2] for a in at]
    mp, ms = mean([x[i] for i in prim]), mean([x[i] for i in sec])
    u = sub(ms, mp); u = [v / nrm(u) for v in u]                  # dar -> genis kenar ekseni
    c = mean([x[i] for i in range(ncd) if not cd_at[i][1].startswith('H')])
    gc = mean([x[i] for i in range(ncd, len(x)) if not g_at[i - ncd][1].startswith('H')])
    d = sub(gc, c); h = dot(d, u); r = math.sqrt(max(0, dot(d, d) - h * h))
    zO = dot(sub(x[gO], c), u)
    hb_gd = sum(1 for a in cd_acc for o, hh in g_OH          # misafir OH -> CD O
                if nrm(sub(x[a], x[o])) <= DHA_R and angle(x[hh], x[o], x[a]) <= DHA_ANG)
    hb_cd = sum(1 for o, hh in cd_OH                         # CD OH -> misafir O
                if nrm(sub(x[gO], x[o])) <= DHA_R and angle(x[hh], x[o], x[gO]) <= DHA_ANG)
    rows.append((t, nrm(d), h, r, zO, hb_gd, hb_cd))

with open('analiz.csv', 'w') as f:
    f.write('t_ps,mesafe_nm,eksen_h_nm,radyal_r_nm,O_eksen_nm,HB_misafir_verici,HB_CD_verici\n')
    for rr in rows: f.write(','.join(f'{v:.4f}' if isinstance(v, float) else str(v) for v in rr) + '\n')

N = len(rows); col = lambda k: [rr[k] for rr in rows]
avg = lambda v: sum(v) / len(v)
sd = lambda v: math.sqrt(sum((a - avg(v)) ** 2 for a in v) / len(v))
print(f"Kare sayisi: {N}  ({rows[0][0]:.0f} - {rows[-1][0]:.0f} ps)")
print(f"CD-misafir mesafesi      : {avg(col(1)):.3f} +/- {sd(col(1)):.3f} nm  (max {max(col(1)):.3f})")
print(f"  eksen boyunca (h)      : {avg(col(2)):+.3f} +/- {sd(col(2)):.3f} nm  (+ = genis kenar)")
print(f"  eksenden uzaklik (r)   : {avg(col(3)):.3f} +/- {sd(col(3)):.3f} nm")
zo = col(4)
print(f"Misafir O eksen konumu   : {avg(zo):+.3f} +/- {sd(zo):.3f} nm")
print(f"  O genis tarafta (>0)   : %{100 * sum(z > 0 for z in zo) / N:.0f} kare")
print(f"H-bagi misafir OH -> CD  : ort {avg(col(5)):.2f}, karelerin %{100 * sum(v > 0 for v in col(5)) / N:.0f}'inde")
print(f"H-bagi CD OH -> misafir O: ort {avg(col(6)):.2f}, karelerin %{100 * sum(v > 0 for v in col(6)) / N:.0f}'inde")
print("Kare kare degerler: analiz.csv")
