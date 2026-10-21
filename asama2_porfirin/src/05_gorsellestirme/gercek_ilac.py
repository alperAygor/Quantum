# Gercek klinik PDT ilaclarini 3B olarak kur, bizim hesapladigimiz cekirdege (kuantum bolgesi) hizala,
# ve ChimeraX'te ilacin uzerinde bizim orbital/elektron farki sonuclarimizi gosteren dosyalar uret.
import os, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
OUT = "/mnt/c/Users/<KULLANICI>/Desktop/PDT_gorsel_HQ"; W = "C:/Users/<KULLANICI>/Desktop/PDT_gorsel_HQ"
CORE = {"porfin": "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2",
        "klorin": "C1CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2"}
ILAC = [  # ad, klinik bilgi, cekirdek, SMILES (Wikipedia/PubChem)
 ("Temoporfin", "Foscan (AB) - bas-boyun kanseri PDT", "klorin",
  "C1CC2=NC1=C(C3=CC=C(N3)C(=C4C=CC(=N4)C(=C5C=CC(=C2C6=CC(=CC=C6)O)N5)C7=CC(=CC=C7)O)C8=CC(=CC=C8)O)C9=CC(=CC=C9)O"),
 ("Talaporfin", "Laserphyrin (Japonya) - akciger kanseri PDT", "klorin",
  "CCC1=C2C=C3C(=C(C(=N3)C(=C4[C@H]([C@@H](C(=N4)C=C5C(=C(C(=CC(=C1C)N2)N5)C=C)C)C)CCC(=O)O)CC(=O)N[C@@H](CC(=O)O)C(=O)O)C(=O)O)C"),
 ("Protoporfirin-IX", "5-ALA ile vucutta olusan fotosensitizor (floresan cerrahi, PDT)", "porfin",
  "C=CC1=C(C)C2=CC3=NC(=CC4=C(C)C(CCC(=O)O)=C(N4)C=C5N=C(C=C1N2)C(C)=C5CCC(=O)O)C(C)=C3C=C"),
]
def generic(m):
    q = Chem.RWMol()
    for a in m.GetAtoms():
        q.AddAtom(Chem.AtomFromSmarts(f"[#{a.GetAtomicNum()}]"))
    for b in m.GetBonds():
        q.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
        q.ReplaceBond(q.GetNumBonds() - 1, Chem.BondFromSmarts("~"))
    return q.GetMol()
def kabsch(P, Q):   # P'yi Q'ya oturtan R,t
    pc, qc = P.mean(0), Q.mean(0); H = (P - pc).T @ (Q - qc); U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); D = np.diag([1, 1, d]); R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T
for ad, bilgi, core, smi in ILAC:
    cm = Chem.MolFromSmiles(CORE[core]); nheavy = cm.GetNumAtoms()
    X = np.array([[float(v) for v in l.split()[1:4]] for l in open(f"{core}.xyz").read().split("\n")[2:] if l.strip()])[:nheavy]
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=11); AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    Y = m.GetConformer().GetPositions()
    best = None
    for mt in m.GetSubstructMatches(generic(cm), uniquify=False, useChirality=False, maxMatches=2000):
        P = Y[list(mt)]; R, t = kabsch(P, X); rmsd = np.sqrt((((P @ R.T + t) - X) ** 2).sum(1).mean())
        if best is None or rmsd < best[0]: best = (rmsd, R, t)
    rmsd, R, t = best; Yt = Y @ R.T + t
    d, dw = f"{OUT}/{core}", f"{W}/{core}"; os.makedirs(d, exist_ok=True)
    pdb = f"{d}/{ad}.pdb"
    with open(pdb, "w") as f:
        for i, (a, p) in enumerate(zip(m.GetAtoms(), Yt)):
            f.write("HETATM%5d %-4s DRG A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (i + 1, (a.GetSymbol() + str(i + 1))[:4], *p, a.GetSymbol()))
        for b in m.GetBonds():
            f.write("CONECT%5d%5d\n" % (b.GetBeginAtomIdx() + 1, b.GetEndAtomIdx() + 1))
        f.write("END\n")
    for lab, lev, cap in [("HOMO", 0.02, "renkli bulutlar: bizim kuantum bolgemiz (4 orbital = 8 kubit) - HOMO"),
                          ("LUMO", 0.02, "renkli bulutlar: bizim kuantum bolgemiz - LUMO (isikla elektronun gittigi orbital)")]:
        cxc = f"""close session
open {dw}/{ad}.pdb
set bgColor white
style ball
size stickRadius 0.13
size ballScale 0.22
color byelement
color C #9a9a9a
lighting soft
lighting shadows true intensity 0.5
lighting multiShadow 1024
material specular 0.4 exponent 40
graphics silhouettes true width 2 color #2a2a2a
graphics quality 6
open {dw}/{core}_{lab}.cube
volume #2 style surface level -{lev} color #2f6fd0 level {lev} color #f0782f transparency 0.25 surfaceSmoothing true smoothingIterations 4
view
turn x -40
turn z 20
view pad 0.1
2dlabels text "{ad}  -  {bilgi}" xpos 0.02 ypos 0.94 size 22 color black bold true
2dlabels text "{cap}" xpos 0.02 ypos 0.03 size 15 color dimgray
save {dw}/{ad}_{lab}.png width 3000 height 2000 supersample 4
"""
        open(f"{d}/{ad}_{lab}.cxc", "w").write(cxc)
    sade = f"""close session
open {dw}/{ad}.pdb
set bgColor white
style sphere
color byelement
color C #8a8a8a
lighting soft
lighting shadows true intensity 0.6
lighting multiShadow 1024
material specular 0.4 exponent 40
graphics silhouettes true width 2 color #2a2a2a
turn x -40
turn z 20
view
view pad 0.1
2dlabels text "{ad}  -  {bilgi}" xpos 0.02 ypos 0.94 size 22 color black bold true
save {dw}/{ad}_molekul.png width 3000 height 2000 supersample 4
"""
    open(f"{d}/{ad}_molekul.cxc", "w").write(sade)
    print(f"{ad}: {m.GetNumAtoms()} atom, cekirdek ({core}) hizalama RMSD {rmsd:.2f} A -> {d}")
with open(f"{OUT}/ilaclar_hepsi.cxc", "w") as f:
    for ad, _, core, _ in ILAC:
        for v in ("molekul", "HOMO", "LUMO"): f.write(f"open {W}/{core}/{ad}_{v}.cxc\n")
print(f"\nChimeraX'te:  open {W}/ilaclar_hepsi.cxc")
