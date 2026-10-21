"""Aşama 3 ve 5: QM/MM ve QM/QM/MM enerji motorları.

Makaledeki karşılığı:
  - 'Determine Quantum Region' + 'Create ML1 Training Data (QM/MM)'  -> energy_dft  (PBE/MM)
  - 'Determine Quantum Cores' + 'Create ML2 Training Data (QM/QM/MM)' -> energy_high (CASCI+NEVPT2/MM)

QM bölgesi = misafir molekül (siklopentanol). Çevre (bCD + su + iyonlar) nokta yükler olarak
QM Hamiltoniyenine girer (elektrostatik gömme, Denklem 6).

Yeniden ağırlıklandırma için QM - MM enerji farkını ikiye ayırıyoruz:
  ΔU(x)        = ΔU_etk(x) + ΔU_iç(x)
  ΔU_etk(x)    = [E_QM(misafir | nokta yükler) - E_QM(misafir, vakum)] - E_MM,Coulomb(misafir-çevre)
                 -> çevreyle ETKİLEŞİMİN kuantum düzeltmesi (polarizasyon vb.)
  ΔU_iç(x)     = E_QM(misafir, vakum) - E_MM,iç(misafir)
                 -> misafirin İÇ enerjisinin düzeltmesi
Misafir-çevre van der Waals terimi her iki seviyede de MM'den geldiği için sadeleşir.
03. notebook'ta gösterildiği gibi ΔU_iç büyük dalgalanır ama ortalaması iki bacakta aynıdır
(aynı molekül, benzer konformasyonlar) ve bağlanma farkında sadeleşir. Bu yüzden ML ve serbest
enerji düzeltmesi ΔU_etk üzerinden yapılır; ΔU_iç'in sadeleştiği 03'te ayrıca kontrol edilir.
"""
import time

import numpy as np

from .config import DATA, SYSTEM, HARTREE_TO_KJ, P

ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O", 11: "Na", 17: "Cl"}


class QMRegion:
    """Bir bacağın (kompleks/çözücü) topolojisinden QM bölgesi ve nokta yük bilgisi."""

    def __init__(self, leg):
        import parmed as pmd

        self.leg = leg
        st = pmd.load_file(str(DATA / f"{leg}.prmtop"))
        top = st.topology
        self.charges = np.array([a.charge for a in st.atoms])
        self.elements = [a.element.atomic_number for a in top.atoms()]
        self.guest = np.array([a.index for a in top.atoms()
                               if a.residue.name == SYSTEM["guest_resname"]])
        self.residues = [np.array([a.index for a in r.atoms()]) for r in top.residues()
                         if r.name != SYSTEM["guest_resname"]]
        self.res_first = np.array([r[0] for r in self.residues])
        self._guest_mm = None

    # ------------------------------------------------------------------
    def guest_symbols(self):
        return [ELEMENTS[self.elements[i]] for i in self.guest]

    def frame(self, X, box, cutoff_A=None):
        """Bir anlık görüntüden QM koordinatları (Å) ve nokta yükleri (Å, e) çıkarır.
        Periyodik kutuda en yakın görüntü (minimum image) kullanılır; moleküller bölünmez."""
        cutoff = (cutoff_A or P()["qm_cutoff_A"]) / 10.0  # nm
        X = np.asarray(X, dtype=float)
        box = np.asarray(box, dtype=float)
        g = X[self.guest]
        c = g.mean(0)
        pc_xyz, pc_q = [], []
        for res in self.residues:
            r = X[res]
            shift = -box * np.round((r[0] - c) / box)
            r = r + shift
            d = np.min(np.linalg.norm(r - c, axis=1))
            # İyonlar (net yüklü) her zaman dahil: kesme sınırından girip çıkan bir iyon
            # enerjide birkaç kJ/mol'lük yapay sıçramalar yaratır.
            if d < cutoff or len(res) == 1:
                pc_xyz.append(r)
                pc_q.append(self.charges[res])
        pc_xyz = np.concatenate(pc_xyz) * 10.0
        pc_q = np.concatenate(pc_q)
        return dict(qm_xyz=g * 10.0, pc_xyz=pc_xyz, pc_q=pc_q)

    # ------------------------------------------------------------------
    def _guest_mm_context(self):
        if self._guest_mm is None:
            import parmed as pmd
            import openmm as mm
            import openmm.app as app

            s = pmd.load_file(str(DATA / f"{self.leg}.prmtop"))
            keep = f":{SYSTEM['guest_resname']}"
            s.strip(f"!{keep}")
            sysg = s.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
            integ = mm.VerletIntegrator(0.001)
            ctx = mm.Context(sysg, integ, mm.Platform.getPlatformByName("Reference"))
            self._guest_mm = (ctx, integ)
        return self._guest_mm[0]

    def mm_qm_part(self, fr):
        """E_MM,iç(misafir) + Coulomb(misafir-nokta yükler)  [kJ/mol]"""
        import openmm.unit as u

        ctx = self._guest_mm_context()
        ctx.setPositions(fr["qm_xyz"] / 10.0)
        e_int = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        qg = self.charges[self.guest]
        d = np.linalg.norm(fr["qm_xyz"][:, None, :] - fr["pc_xyz"][None, :, :], axis=2) / 10.0  # nm
        e_coul = 138.935458 * np.sum(qg[:, None] * fr["pc_q"][None, :] / d)
        return e_int + e_coul, e_int, e_coul


# ----------------------------------------------------------------------
# PySCF enerji motorları
# ----------------------------------------------------------------------
def _mol(fr, symbols, basis):
    from pyscf import gto

    return gto.M(atom=[(s, x) for s, x in zip(symbols, fr["qm_xyz"])],
                 basis=basis, charge=0, spin=0, verbose=0)


def _embed(mf, fr, embed):
    from pyscf import qmmm

    return qmmm.mm_charge(mf, fr["pc_xyz"], fr["pc_q"]) if embed else mf


def energy_dft(fr, symbols, xc=None, basis=None, embed=True):
    """QM/MM seviyesi (makalede PBE-D3/def2-SVP; burada PBE/def2-SVP). Hartree döner.
    embed=False: aynı geometride vakum hesabı."""
    from pyscf import dft

    p = P()
    mol = _mol(fr, symbols, basis or p["qm_basis"])
    mf = _embed(dft.RKS(mol, xc=xc or p["qm_xc"]).density_fit(), fr, embed)
    mf.conv_tol = 1e-8
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError("DFT yakınsamadı")
    return e


def hf_in_field(fr, symbols, basis=None, embed=True):
    from pyscf import scf

    mol = _mol(fr, symbols, basis or P()["qm_basis"])
    mf = _embed(scf.RHF(mol).density_fit(), fr, embed)
    mf.conv_tol = 1e-9
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("HF yakınsamadı")
    return mf


def energy_high(fr, symbols, cas=None, basis=None, embed=True, return_objects=False):
    """'Kuantum çekirdeği' seviyesi: CASCI(n,m) + NEVPT2, nokta yük alanında.

    - CASCI kısmı = aktif uzayda tam köşegenleştirme (FCI). Gelecekte bu adımı
      bir kuantum bilgisayar (kuantum faz tahmini, QPE) yapacak.
    - NEVPT2 = aktif uzay dışındaki dinamik korelasyon (klasik pertürbasyon kuramı).
    Makalede: DMRG-CASCI + NEVPT2 (Ru d-orbitalleri). Burada küçük bir organik molekül
    olduğu için CASCI'yi PySCF tam olarak çözebiliyor.
    """
    from pyscf import mcscf
    from pyscf.mrpt import nevpt2

    ncas, nelecas = cas or P()["cas"]
    mf = hf_in_field(fr, symbols, basis, embed=embed)
    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.verbose = 0
    e_cas = mc.kernel()[0]
    e_pt2 = nevpt2.NEVPT(mc).set(verbose=0).kernel()
    out = dict(e_hf=mf.e_tot, e_casci=e_cas, e_nevpt2_corr=e_pt2, e_total=e_cas + e_pt2)
    if return_objects:
        return out, mf, mc
    return out


# ----------------------------------------------------------------------
class Labeler:
    """Bir anlık görüntü için tüm etiketleri (DFT, yüksek seviye, MM) hesaplar ve önbelleğe alır.
    (Makaledeki 'QM/MM Calculations' modülünün basit karşılığı.)
    Her etiket = gömülü + vakum olmak üzere iki QM hesabı."""

    def __init__(self, regions, pools):
        self.regions = regions        # {"complex": QMRegion, "solvent": QMRegion}
        self.pools = pools            # {"complex": (X, B), "solvent": (X, B)}
        self.cache_dft = {}           # key -> dict(emb, vac, mm_int, mm_coul)   [kJ/mol]
        self.cache_high = {}          # key -> dict(emb_*, vac_*)                [kJ/mol]
        self.timings = {"dft": [], "high": []}
        self._sig = None

    def signature(self):
        """Etiketlerin hangi yapılara ve QM ayarlarına ait olduğunu gösteren parmak izi.
        Havuz yeniden üretilirse (02 tekrar çalışırsa) ya da qm_cutoff_A / baz / fonksiyonel /
        aktif uzay değişirse değişir; böylece eski etiketler yeni yapılarla karışmaz."""
        if self._sig is None:
            import hashlib

            h = hashlib.sha1()
            for leg in sorted(self.pools):
                X, B = self.pools[leg]
                h.update(leg.encode())
                h.update(np.ascontiguousarray(X).tobytes())
                h.update(np.ascontiguousarray(B).tobytes())
            p = P()
            h.update(repr([p["qm_cutoff_A"], p["qm_basis"], p["qm_xc"], tuple(p["cas"])]).encode())
            self._sig = h.hexdigest()[:16]
        return self._sig

    def _fr(self, leg, i):
        X, B = self.pools[leg]
        return self.regions[leg].frame(X[i], B[i])

    def _dft(self, leg, i):
        key = (leg, int(i))
        if key not in self.cache_dft:
            fr = self._fr(leg, i)
            sym = self.regions[leg].guest_symbols()
            t = time.time()
            emb = energy_dft(fr, sym) * HARTREE_TO_KJ
            vac = energy_dft(fr, sym, embed=False) * HARTREE_TO_KJ
            self.timings["dft"].append(time.time() - t)
            _, e_int, e_coul = self.regions[leg].mm_qm_part(fr)
            self.cache_dft[key] = dict(emb=emb, vac=vac, mm_int=e_int, mm_coul=e_coul)
        return self.cache_dft[key]

    def dU_dft(self, leg, i):
        """ΔU_etk (QM/MM - MM), etkileşim kısmı [kJ/mol] — ML1'in hedefi."""
        d = self._dft(leg, i)
        return d["emb"] - d["vac"] - d["mm_coul"]

    def dU_intra(self, leg, i):
        """ΔU_iç = E_QM(vakum) - E_MM,iç [kJ/mol] — tanı amaçlı."""
        d = self._dft(leg, i)
        return d["vac"] - d["mm_int"]

    def _high(self, leg, i):
        key = (leg, int(i))
        if key not in self.cache_high:
            fr = self._fr(leg, i)
            sym = self.regions[leg].guest_symbols()
            t = time.time()
            he = energy_high(fr, sym)
            hv = energy_high(fr, sym, embed=False)
            self.timings["high"].append(time.time() - t)
            rec = {f"emb_{k}": v * HARTREE_TO_KJ for k, v in he.items()}
            rec.update({f"vac_{k}": v * HARTREE_TO_KJ for k, v in hv.items()})
            self.cache_high[key] = rec
        return self.cache_high[key]

    def high_minus_dft(self, leg, i):
        """δ = ΔU_etk(yüksek seviye) - ΔU_etk(DFT) [kJ/mol] — ML2'nin hedefi.
        Yani kuantum çekirdeği seviyesinin, misafir-çevre etkileşimini DFT'den ne kadar farklı gördüğü."""
        d = self._dft(leg, i)
        h = self._high(leg, i)
        return (h["emb_e_total"] - h["vac_e_total"]) - (d["emb"] - d["vac"])

    def to_dict(self):
        return dict(
            signature=self.signature(),
            dft=[[k[0], k[1], v] for k, v in self.cache_dft.items()],
            high=[[k[0], k[1], v] for k, v in self.cache_high.items()],
            timings=self.timings,
        )

    def load_dict(self, d):
        if d and d.get("signature") != self.signature():
            print("Uyarı: kayıtlı QM etiketleri başka bir havuza veya QM ayarına ait; "
                  "kullanılmıyor, yeniden hesaplanacak.")
            return
        for leg, i, v in d.get("dft", []):
            self.cache_dft[(leg, int(i))] = v
        for leg, i, v in d.get("high", []):
            self.cache_high[(leg, int(i))] = v
        for k in ("dft", "high"):
            self.timings[k] += d.get("timings", {}).get(k, [])
