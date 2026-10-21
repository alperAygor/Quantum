import sys, numpy as np
from collections import defaultdict
xyz={};name={};nb=defaultdict(set)
for l in open(sys.argv[1]):
    if l.startswith(('ATOM','HETATM')):
        i=int(l[6:11]); name[i]=l[12:16].strip(); xyz[i]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    if l.startswith('CONECT'):
        f=[int(l[k:k+5]) for k in range(6,len(l.rstrip()),5) if l[k:k+5].strip()]
        for j in f[1:]: nb[f[0]].add(j); nb[j].add(f[0])
el=lambda i: name[i][0]
heavy=lambda i: [j for j in nb[i] if el(j)!='H']
def key(j, frm):
    return (el(j), len(heavy(j)), tuple(sorted(el(k)+str(len(heavy(k))) for k in heavy(j) if k!=frm)), sum(1 for k in nb[j] if el(k)=='H'))
res=defaultdict(list)
for i in xyz:
    if el(i)=='C' and len(nb[i])==4 and sum(el(j)=='H' for j in nb[i])==1:
        hv=sorted(heavy(i), key=lambda j:key(j,i))
        keys=[key(j,i) for j in hv]
        a,b,c=[xyz[j]-xyz[i] for j in hv]
        res[tuple(keys)].append((name[i], int(np.sign(np.dot(a,np.cross(b,c))))))
for k,v in res.items():
    print(len(v), [s for _,s in v], [n for n,_ in v])
