import sys, numpy as np
from collections import defaultdict
def load(f):
    xyz={};name={};nb=defaultdict(set)
    for l in open(f):
        if l.startswith(('ATOM','HETATM')):
            i=int(l[6:11]); name[i]=l[12:16].strip(); xyz[i]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])/10
        if l.startswith('CONECT'):
            fl=[int(l[k:k+5]) for k in range(6,len(l.rstrip()),5) if l[k:k+5].strip()]
            for j in fl[1:]: nb[fl[0]].add(j); nb[j].add(fl[0])
    return xyz,name,nb
def analyze(f):
    xyz,name,nb=load(f); el=lambda i:name[i][0]
    nH=lambda i:sum(el(j)=='H' for j in nb[i])
    C6=[i for i in xyz if el(i)=='C' and nH(i)==2]
    C5=[next(j for j in nb[c] if el(j)=='C') for c in C6]
    Os=[i for i in xyz if el(i)=='O' and sum(el(j)=='C' for j in nb[i])==2]
    ringO=[o for o in Os if any(j in C5 for j in nb[o])]
    gly=[o for o in Os if o not in ringO]
    O6=[next(j for j in nb[c] if el(j)=='O') for c in C6]
    O23=[i for i in xyz if el(i)=='O' and nH(i)==1 and i not in O6]
    G=np.array([xyz[o] for o in gly]); cen=G.mean(0)
    _,_,vt=np.linalg.svd(G-cen); n=vt[2]
    side=np.array([(xyz[c]-cen)@n for c in C6])
    if (side>0).sum()<(side<0).sum(): n=-n; side=-side   # n: C6 cogunluk tarafi
    rim=np.mean([xyz[o] for o in O23],0)-np.mean([xyz[o] for o in O6],0)
    # kavite: eksene en yakin atom (|h|<0.35 nm)
    rad=[]
    for i,p in xyz.items():
        d=p-cen; h=d@n
        if abs(h)<0.35: rad.append((np.linalg.norm(d-h*n),name[i]))
    rmin=min(rad)
    r=np.linalg.norm(G-cen,axis=1)
    print(f"{f.split('/')[-1][9:]:16s} gly_O={len(gly)} C6_ayni_taraf={int((side>0).sum())}/7  "
          f"rim_mes={np.linalg.norm(rim):.3f}nm  glyO_r={r.min():.2f}-{r.max():.2f}  "
          f"glyO_duzlem_disi_max={np.abs((G-cen)@n).max():.3f}  kavite_min_r={rmin[0]:.3f}({rmin[1]})")
for f in sys.argv[1:]: analyze(f)
