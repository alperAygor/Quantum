"""Cekme yorungelerinden umbrella pencerelerinin baslangic karelerini secer.
xi < 0 -> pull_minus, xi >= 0 -> pull_plus. Kareler 1 ps aralikla kayitli."""
def oku(f):
    d = []
    for l in open(f):
        if l[0] in '#@': continue
        p = l.split(); t = float(p[0])
        if abs(t - round(t)) < 1e-6: d.append((round(t), float(p[1])))
    return d
veri = {'plus': oku('pull_plus_pullx.xvg'), 'minus': oku('pull_minus_pullx.xvg')}
hedefler = [round(-1.8 + 0.1 * i, 1) for i in range(37)]
with open('pencereler.txt', 'w') as f:
    for i, h in enumerate(hedefler):
        yon = 'minus' if h < 0 else 'plus'
        t, xi = min(veri[yon], key=lambda r: abs(r[1] - h))
        f.write(f"{i:02d} {h:+.1f} {yon} {t} {xi:+.3f}\n")
        flag = '' if abs(xi - h) < 0.05 else '   <-- hedefe uzak!'
        print(f"pencere {i:02d}: hedef {h:+.1f} nm | {yon:5s} t={t:4d} ps | xi={xi:+.3f}{flag}")
