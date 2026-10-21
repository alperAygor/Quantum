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
def key(j, frm, d):
    nh=sum(el(k)=='H' for k in nb[j])
    if d==0: return (el(j), nh)
    return (el(j), nh, tuple(sorted(key(k,j,d-1) for k in nb[j] if k!=frm and el(k)!='H')))
res=defaultdict(list)
for i in xyz:
    if el(i)=='C' and len(nb[i])==4 and sum(el(j)=='H' for j in nb[i])==1:
        hv=[j for j in nb[i] if el(j)!='H']
        ks=[key(j,i,4) for j in hv]
        assert len(set(ks))==3, (name[i])
        hv=[j for _,j in sorted(zip(ks,hv))]
        a,b,c=[xyz[j]-xyz[i] for j in hv]
        res[tuple(sorted(ks))].append((name[i], int(np.sign(np.dot(a,np.cross(b,c))))))
for k in sorted(res): v=res[k]; print(sorted(set(s for _,s in v)))
