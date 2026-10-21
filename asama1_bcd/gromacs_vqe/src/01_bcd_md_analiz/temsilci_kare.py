#!/usr/bin/env python3
"""
temsilci_kare.py  -  view.gro'dan iki baglanma modu icin temsilci kare secer.

Mod A: siklopentanol OH'i DAR kenara (primer OH tarafi) bakiyor   (cos > +0.5)
Mod B: siklopentanol OH'i GENIS kenara (sekonder OH tarafi) bakiyor (cos < -0.5)

Her mod icin, o modun medyan geometrisine (h, r, cos) en yakin kare secilir
ve PDB + XYZ (Angstrom) olarak yazilir. XYZ dosyalari dogrudan PySCF'e girer.

Kullanim:  python3 temsilci_kare.py view.gro
"""
import sys
import numpy as np

GRO = sys.argv[1] if len(sys.argv) > 1 else "view.gro"
PRIM = {"O35", "O30", "O25", "O20", "O15", "O10", "O5"}   # primer O (345755)
GLYC = [21, 42, 63, 84, 105, 126, 138]                     # glikozidik O, CD icinde 1-tabanli


def frames(fn):
    with open(fn) as f:
        while True:
            title = f.readline()
            if not title.strip():
                break
            n = int(f.readline())
            atoms = [f.readline() for _ in range(n)]
            f.readline()  # kutu satiri
            yield title.strip(), atoms


recs = []
setup = None
for k, (title, atoms) in enumerate(frames(GRO)):
    xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in atoms])
    if setup is None:
        resn = [a[5:10].strip() for a in atoms]
        name = [a[10:15].strip() for a in atoms]
        cd = [i for i, r in enumerate(resn) if r == "AOYB"]
        gu = [i for i, r in enumerate(resn) if r == "G108"]
        prim = [i for i in cd if name[i] in PRIM]
        gly = [cd[j - 1] for j in GLYC]
        iO = [i for i in gu if name[i] == "O1"][0]
        iC = [i for i in gu if name[i] == "C4"][0]
        keep = cd + gu
        setup = True
        print(f"CD atom: {len(cd)}  misafir atom: {len(gu)}  primer O: {len(prim)}")

    c = xyz[gly].mean(0)
    _, _, vt = np.linalg.svd(xyz[gly] - c)
    n = vt[2]
    if np.dot(xyz[prim].mean(0) - c, n) < 0:   # ekseni dar (primer) kenara yonlendir
        n = -n
    g = xyz[gu].mean(0) - c
    h = float(np.dot(g, n))
    r = float(np.linalg.norm(g - h * n))
    v = xyz[iO] - xyz[iC]
    cos = float(np.dot(v, n) / np.linalg.norm(v))
    t = title.split("t=")[-1].split()[0] if "t=" in title else str(k)
    recs.append((k, t, h, r, cos, xyz - c))

H = np.array([x[2] for x in recs])
R = np.array([x[3] for x in recs])
C = np.array([x[4] for x in recs])
N = len(recs)
print(f"Toplam kare: {N}")


def write_files(tag, rec):
    k, t, h, r, cos, xyz = rec
    with open(f"mod_{tag}.pdb", "w") as f:
        for j, i in enumerate(keep):
            rn, resid = ("BCD", 1) if resn[i] == "AOYB" else ("CPO", 2)
            x, y, z = xyz[i] * 10.0
            el = name[i][0]
            f.write("HETATM%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
                    % (j + 1, name[i], rn, resid, x, y, z, el))
        f.write("END\n")
    with open(f"mod_{tag}.xyz", "w") as f:
        f.write(f"{len(keep)}\nmod {tag}  t={t} ps  h={h:.3f} nm  r={r:.3f} nm  cos={cos:+.2f}\n")
        for i in keep:
            x, y, z = xyz[i] * 10.0
            f.write(f"{name[i][0]:2s} {x:12.5f} {y:12.5f} {z:12.5f}\n")


for tag, mask, lab in [("A", C > 0.5, "OH dar kenara (primer)"),
                       ("B", C < -0.5, "OH genis kenara (sekonder)")]:
    idx = np.where(mask)[0]
    if len(idx) == 0:
        print(f"Mod {tag}: kare yok!")
        continue
    med = np.array([np.median(H[idx]), np.median(R[idx]), np.median(C[idx])])
    sd = np.array([H[idx].std(), R[idx].std(), C[idx].std()]) + 1e-6
    X = np.stack([H[idx], R[idx], C[idx]], 1)
    best = idx[np.argmin((((X - med) / sd) ** 2).sum(1))]
    rec = recs[best]
    write_files(tag, rec)
    print(f"\nMod {tag} ({lab}): {len(idx)} kare = %{100*len(idx)/N:.1f}")
    print(f"  medyan  h={med[0]:+.3f} nm  r={med[1]:.3f} nm  cos={med[2]:+.2f}")
    print(f"  secilen t={rec[1]} ps  h={rec[2]:+.3f}  r={rec[3]:.3f}  cos={rec[4]:+.2f}"
          f"  -> mod_{tag}.pdb / mod_{tag}.xyz")

print(f"\nAra (|cos|<=0.5): {np.sum(np.abs(C) <= 0.5)} kare")
