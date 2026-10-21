import math,sys
f=sys.argv[1]
at=[]
for l in open(f):
    if l.startswith(('ATOM','HETATM')):
        at.append((l[12:16].strip(),[float(l[30:38])/10,float(l[38:46])/10,float(l[46:54])/10]))
prim={'O34','O4','O9','O14','O19','O24','O29'}
sec={'O33','O32','O3','O2','O8','O7','O13','O12','O18','O17','O23','O22','O28','O27'}
mean=lambda ps:[sum(p[i] for p in ps)/len(ps) for i in range(3)]
mp=mean([x for n,x in at if n in prim]); ms=mean([x for n,x in at if n in sec])
print("rim distance nm", round(math.dist(mp,ms),3))
# per-glucose: projection of O6 and O2/O3 along axis
u=[ms[i]-mp[i] for i in range(3)]; L=math.sqrt(sum(v*v for v in u)); u=[v/L for v in u]
c=mean([x for n,x in at if not n.startswith('H')])
proj=lambda x: sum((x[i]-c[i])*u[i] for i in range(3))
print("primary O proj:", [round(proj(x),2) for n,x in at if n in prim])
print("secondary O proj:", [round(proj(x),2) for n,x in at if n in sec])
