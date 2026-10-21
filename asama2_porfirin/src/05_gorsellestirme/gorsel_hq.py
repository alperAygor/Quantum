# Yuksek kaliteli ChimeraX gorselleri: ince izgara (0.1 A), genis kutu, film kalitesinde isik.
# Kullanim: python3 gorsel_hq.py [molekul ...]
import os, sys, glob, json
import numpy as np
from pyscf import gto, scf, dft, mcscf, fci
from pyscf.tools import cubegen, molden
OUT_WIN = "C:/Users/<KULLANICI>/Desktop/PDT_gorsel_HQ"
OUT = "/mnt/c/Users/<KULLANICI>/Desktop/PDT_gorsel_HQ"
RES, MARGIN = 0.1 / 0.529177, 7.0          # izgara araligi 0.1 A (bohr cinsinden), kenar payi 7 bohr (~3.7 A)
os.makedirs(OUT, exist_ok=True)
names = sys.argv[1:] or [f.replace("_v3_sonuc.json", "") for f in sorted(glob.glob("*_v3_sonuc.json"))]
ELEM = {"C": "#909090", "N": "#3050f8", "H": "#f4f4f4", "O": "#ff0d0d", "Zn": "#7d80b0", "Mg": "#8aff00", "Pd": "#006985"}
STYLE = """set bgColor white
style ball
size stickRadius 0.11
size ballScale 0.22
{elem}
lighting full
lighting shadows true intensity 0.5
lighting multiShadow 512
lighting depthCue false
material specular 0.3 exponent 30
graphics silhouettes true width 2 color #303030
graphics quality 5
camera ortho
"""
master = []
for n in names:
    R = json.load(open(f"{n}_v3_sonuc.json")) if os.path.exists(f"{n}_v3_sonuc.json") else {}
    d, dw = f"{OUT}/{n}", f"{OUT_WIN}/{n}"; os.makedirs(d, exist_ok=True)
    L = [l.split() for l in open(f"{n}.xyz").read().split("\n")[2:] if l.strip()]
    ecp = "def2-svp" if any(l[0] == "Pd" for l in L) else None
    mol = gto.M(atom=f"{n}.xyz", basis="def2-svp", ecp=ecp, verbose=0, max_memory=5000)
    with open(f"{d}/{n}.pdb", "w") as f:
        for i, l in enumerate(L):
            f.write("HETATM%5d %-4s MOL A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (i + 1, (l[0] + str(i + 1))[:4], *map(float, l[1:4]), l[0]))
        f.write("END\n")
    ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
    molden.from_scf(ks, f"{d}/{n}.molden")
    nocc = mol.nelectron // 2
    orbs = [(nocc - 2, "HOMO-1"), (nocc - 1, "HOMO"), (nocc, "LUMO"), (nocc + 1, "LUMO+1")]
    for k, lab in orbs:
        cubegen.orbital(mol, f"{d}/{n}_{lab}.cube", ks.mo_coeff[:, k], resolution=RES, margin=MARGIN)
    cas = mcscf.CASCI(scf.RHF(mol).density_fit(), 4, 4)
    cas.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0); cas.fcisolver.nroots = 2
    cas.kernel(ks.mo_coeff)
    g0 = cas.fcisolver.make_rdm1(cas.ci[0], 4, 4); g1 = cas.fcisolver.make_rdm1(cas.ci[1], 4, 4)
    Ca = ks.mo_coeff[:, cas.ncore:cas.ncore + 4]
    dd = cubegen.density(mol, f"{d}/{n}_fark.cube", Ca @ (g1 - g0) @ Ca.T, resolution=RES, margin=MARGIN)
    flev = float(np.round(0.12 * np.abs(dd).max(), 5))           # fark haritasi icin otomatik esik
    s1, f1, tr = R.get("S1_NEVPT2_eV"), R.get("f_S1"), R.get("triyaj", "")
    title = f"{n}   |   S1 = {s1:.2f} eV ({1239.84 / s1:.0f} nm)   |   f = {f1:.3f}   |   {tr}" if s1 else n
    elem = "\n".join(f"color {e} {c} target a" for e, c in ELEM.items() if any(l[0] == e for l in L))
    views = [("HOMO", 0.02, "HOMO: isik gelmeden once elektronun bulundugu orbital"),
             ("LUMO", 0.02, "LUMO: isikla elektronun siçradigi orbital"),
             ("fark", flev, "S0 -> S1: turuncu = elektronun geldigi yer, mavi = elektronun ayrildigi yer")]
    for lab, lev, cap in views:
        cxc = f"""close session
open {dw}/{n}.pdb
{STYLE.format(elem=elem)}
open {dw}/{n}_{lab}.cube
volume #2 style surface level {-lev} color #2f6fd0 level {lev} color #f0782f transparency 0.1 surfaceSmoothing true smoothingIterations 4 squareMesh false
view
turn x -35
turn z 15
view pad 0.15
2dlabels text "{title}" xpos 0.02 ypos 0.95 size 20 color black bold true
2dlabels text "{cap}" xpos 0.02 ypos 0.03 size 15 color dimgray
save {dw}/{n}_{lab}.png width 3000 height 2000 supersample 4 transparentBackground false
"""
        open(f"{d}/{n}_{lab}.cxc", "w").write(cxc)
        master.append(f"open {dw}/{n}/{n}_{lab}.cxc")
    vid = open(f"{d}/{n}_fark.cxc").read().split("save ")[0] + f"movie record supersample 3 size 1920,1080\nturn y 1 360\nwait 360\nmovie encode {dw}/{n}_video.mp4 quality highest framerate 30\n"
    open(f"{d}/{n}_video.cxc", "w").write(vid)
    print(f"{n}: tamam (fark esigi {flev})", flush=True)
open(f"{OUT}/hepsi.cxc", "w").write("\n".join(master) + "\n")
print(f"\nBITTI. ChimeraX'te:  open {OUT_WIN}/hepsi.cxc")
