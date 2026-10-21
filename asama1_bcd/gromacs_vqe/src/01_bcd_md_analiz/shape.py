import math,sys
import numpy as np
at={}
for l in open(sys.argv[1]):
    if l.startswith(('ATOM','HETATM')):
        at[l[12:16].strip()]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])/10
gly=['O5','O10','O15','O20','O25','O30','O35']
c6=['C42','C6','C12','C18','C24','C30','C36']
G=np.array([at[n] for n in gly]); cen=G.mean(0)
_,s,vt=np.linalg.svd(G-cen); n=vt[2]
print("glycosidic O radii (nm):", np.round(np.linalg.norm(G-cen,axis=1),3))
print("glycosidic O out-of-plane (nm):", np.round((G-cen)@n,3))
print("C6 side wrt glyco plane (nm):", np.round((np.array([at[x] for x in c6])-cen)@n,3))
