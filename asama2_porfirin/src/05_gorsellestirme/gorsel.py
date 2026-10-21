# Her molekul icin ChimeraX/Avogadro gorselleri: orbitaller + "elektron nereden nereye" (S1 - S0 yogunluk farki)
# Kullanim: python3 gorsel.py [molekul ...]     (varsayilan: *_v3_sonuc.json olan tum molekuller)
import os, sys, glob, json
import numpy as np
from pyscf import gto, scf, dft, mcscf, fci
from pyscf.tools import cubegen, molden
OUT_WIN = "C:/Users/<KULLANICI>/Desktop/PDT_gorsel"            # ChimeraX'in gorecegi Windows yolu
OUT = "/mnt/c/Users/<KULLANICI>/Desktop/PDT_gorsel"            # ayni klasorun WSL yolu
os.makedirs(OUT, exist_ok=True)
names = sys.argv[1:] or [f.replace("_v3_sonuc.json", "") for f in sorted(glob.glob("*_v3_sonuc.json"))]
master = []
for n in names:
    R = json.load(open(f"{n}_v3_sonuc.json")) if os.path.exists(f"{n}_v3_sonuc.json") else {}
    d = f"{OUT}/{n}"; os.makedirs(d, exist_ok=True); dw = f"{OUT_WIN}/{n}"
    L = [l.split() for l in open(f"{n}.xyz").read().split("\n")[2:] if l.strip()]
    ecp = "def2-svp" if any(l[0] == "Pd" for l in L) else None
    mol = gto.M(atom=f"{n}.xyz", basis="def2-svp", ecp=ecp, verbose=0, max_memory=5000)
    # PDB (ChimeraX icin)
    with open(f"{d}/{n}.pdb", "w") as f:
        for i, l in enumerate(L):
            f.write("HETATM%5d %-4s MOL A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (i + 1, l[0] + str(i + 1)[-2:], *map(float, l[1:4]), l[0]))
        f.write("END\n")
    ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
    molden.from_scf(ks, f"{d}/{n}.molden")                  # Avogadro icin tum orbitaller
    nocc = mol.nelectron // 2
    for k, lab in [(nocc - 2, "HOMO-1"), (nocc - 1, "HOMO"), (nocc, "LUMO"), (nocc + 1, "LUMO+1")]:
        cubegen.orbital(mol, f"{d}/{n}_{lab}.cube", ks.mo_coeff[:, k], nx=60, ny=60, nz=40)
    # S0 ve S1 (CAS 4,4, singlet) -> yogunluk farki
    cas = mcscf.CASCI(scf.RHF(mol).density_fit(), 4, 4)
    cas.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0); cas.fcisolver.nroots = 2
    cas.kernel(ks.mo_coeff)
    g0 = cas.fcisolver.make_rdm1(cas.ci[0], 4, 4); g1 = cas.fcisolver.make_rdm1(cas.ci[1], 4, 4)
    Ca = ks.mo_coeff[:, cas.ncore:cas.ncore + 4]
    cubegen.density(mol, f"{d}/{n}_fark.cube", Ca @ (g1 - g0) @ Ca.T, nx=60, ny=60, nz=40)
    s1 = R.get("S1_NEVPT2_eV"); f1 = R.get("f_S1"); tr = R.get("triyaj", "")
    nm = f"{1239.84 / s1:.0f} nm" if s1 else ""
    etiket = f"{n}:  S1 = {s1:.2f} eV ({nm}),  f = {f1:.3f},  {tr}" if s1 else n
    # ChimeraX komut dosyasi
    cxc = f"""close session
open {dw}/{n}.pdb
style stick
color byelement
color C gray
open {dw}/{n}_fark.cube
volume #2 style surface level -0.0015 color #2a78d6 level 0.0015 color #eb6834 transparency 0.15
lighting soft
graphics silhouettes true
set bgColor white
view
2dlabels text "{etiket}" xpos 0.02 ypos 0.94 size 18 color black
2dlabels text "turuncu: elektronun geldigi yer   |   mavi: elektronun ayrildigi yer (S0 -> S1)" xpos 0.02 ypos 0.03 size 14 color dimgray
save {dw}/{n}_fark.png width 1600 height 1000 supersample 3
"""
    open(f"{d}/{n}_goster.cxc", "w").write(cxc)
    master.append(f"open {dw}/{n}/{n}_goster.cxc")
    # donen video
    open(f"{d}/{n}_video.cxc", "w").write(cxc.replace(f"save {dw}/{n}_fark.png width 1600 height 1000 supersample 3",
        f"movie record supersample 3\nturn y 2 180\nwait 180\nmovie encode {dw}/{n}_video.mp4 quality high"))
    print(f"{n}: tamam -> {d}", flush=True)
open(f"{OUT}/hepsi_png.cxc", "w").write("\n".join(master) + "\n")
print(f"\nChimeraX'te: open {OUT_WIN}/hepsi_png.cxc   (tum molekullerin PNG'lerini kaydeder)")
