# Porfin T1 duzeltmesi: S0, S1, T1 icin ortak (durum-ortalamali) CASSCF(4,4) orbitalleri + NEVPT2
import os, sys, time
os.makedirs(os.path.expanduser("~/pyscf_tmp"), exist_ok=True)
os.environ.setdefault("PYSCF_TMPDIR", os.path.expanduser("~/pyscf_tmp"))
import numpy as np
from pyscf import gto, scf, dft, mcscf, fci, mrpt
NEV = not (len(sys.argv) > 1 and sys.argv[1] == "nonev")
EV = 27.211386; t0 = time.time()
mol = gto.M(atom="porfin.xyz", basis="def2-svp", verbose=0, max_memory=5000)
ks = dft.RKS(mol, xc="b3lyp").density_fit(); ks.grids.level = 1; ks.run()
print(f"B3LYP (baslangic orbitalleri)  [{time.time()-t0:.0f} s]", flush=True)
mf = scf.RHF(mol).density_fit()

# 1) SA-CASSCF(4,4): S0 + S1 (singlet) ve T1 (triplet) esit agirlikli ortak orbitaller
sS = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0); sS.nroots = 2
sT = fci.direct_spin1.FCI(mol); sT.spin = 2; sT.nroots = 1
mc = mcscf.CASSCF(mf, 4, 4)
mcscf.addons.state_average_mix_(mc, [sS, sT], [1/3, 1/3, 1/3])
mc.max_cycle_macro = 50; mc.verbose = 4
mc.kernel(ks.mo_coeff)
print(f"SA-CASSCF yakinsadi: {mc.converged}  [{time.time()-t0:.0f} s]", flush=True)
mo = mc.mo_coeff

# 2) Ortak orbitallerle durum-ozgu CASCI: S0,S1 (singlet), T1 (yuksek spin triplet)
cS = mcscf.CASCI(mf, 4, 4); cS.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0)
cS.fcisolver.nroots = 2; cS.kernel(mo)
cT = mcscf.CASCI(mf, 4, (3, 1)); cT.kernel(mo)
eS0, eS1, eT1 = cS.e_tot[0], cS.e_tot[1], cT.e_tot
print(f"CASCI (SA orbitalleri): S1 {(eS1-eS0)*EV:.3f} eV, T1 {(eT1-eS0)*EV:.3f} eV  [{time.time()-t0:.0f} s]", flush=True)

# 3) NEVPT2
if NEV:
    n0 = mrpt.NEVPT(cS, root=0).kernel(); print(f"  NEVPT2 S0  [{time.time()-t0:.0f} s]", flush=True)
    n1 = mrpt.NEVPT(cS, root=1).kernel(); print(f"  NEVPT2 S1  [{time.time()-t0:.0f} s]", flush=True)
    nT = mrpt.NEVPT(cT).kernel();         print(f"  NEVPT2 T1  [{time.time()-t0:.0f} s]", flush=True)
else:
    n0 = n1 = nT = np.nan
S1n, T1n = (eS1 + n1 - eS0 - n0) * EV, (eT1 + nT - eS0 - n0) * EV
print("\nPORFIN, SA-CASSCF(4,4) orbitalleri (eV)")
print(f"{'':6s} {'CAS':>7s} {'+NEVPT2':>8s}   Deney / Literatur")
print(f"{'S1':6s} {(eS1-eS0)*EV:7.3f} {S1n:8.3f}   deney 1.98-2.02, lit. NEVPT2 1.96")
print(f"{'T1':6s} {(eT1-eS0)*EV:7.3f} {T1n:8.3f}")
if NEV:
    print(f"S1 - T1 = {S1n-T1n:.2f} eV  ->  sira {'DOGRU (T1 < S1)' if T1n < S1n else 'HALA YANLIS (T1 > S1)'}")
    print(f"PDT: T1 = {T1n:.2f} eV; singlet oksijen icin T1 > 0.98 eV -> {'EVET' if T1n > 0.98 else 'HAYIR'}")
print(f"Toplam sure: {time.time()-t0:.0f} s")
