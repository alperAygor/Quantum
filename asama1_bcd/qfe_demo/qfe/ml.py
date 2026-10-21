"""Aşama 4-5: ML potansiyelleri (ML1, ML2), aktif öğrenme ve serbest enerji yeniden ağırlıklandırma.

Makaledeki karşılığı:
  - ML1: QM/MM verisiyle eğitilen ML kuvvet alanı, aktif öğrenme döngüsü (Bölüm 2.2)
  - ML2: QM/QM/MM verisiyle transfer öğrenme (Bölüm 2.2)
  - NEQ anahtarlama ile MM -> ML/MM serbest enerji düzeltmesi (Şekil 6)

BU DEMODAKİ BASİTLEŞTİRME (ekibe açıkça söyleyin):
  Makale, kuvvetleri de veren tam bir ML potansiyeliyle MD yapar ve 10 ps'lik
  dengesiz (NEQ) anahtarlama kullanır. Biz burada ML modelini yalnızca ENERJİ FARKI
  (ΔU = E_QM - E_MM) tahmincisi olarak kullanıyoruz ve tek adımlı serbest enerji
  pertürbasyonu (Zwanzig; NEQ'nun 'anlık anahtarlama' limiti) uyguluyoruz.
  Model: Gauss süreci (GP) — tahminle birlikte belirsizlik de verdiği için aktif
  öğrenmeyi göstermek için ideal.
"""
import numpy as np

from .config import KT


# ----------------------------------------------------------------------
# Tanımlayıcılar (descriptors)
# ----------------------------------------------------------------------
def descriptors(region, X, B):
    """Her anlık görüntü için çevrenin misafir üzerindeki elektrostatik etkisini özetleyen
    fiziksel tanımlayıcılar (ML1'in hedefi olan ΔU_etk büyük ölçüde polarizasyondan gelir):
      - her misafir atomundaki elektrostatik potansiyel  φ_i = Σ q_j / r_ij
      - her misafir atomundaki elektrik alanın karesi   |E_i|²  (polarizasyon enerjisi ~ -½ α |E|²)
    """
    feats = []
    for x, b in zip(X, B):
        fr = region.frame(x, b)
        g = fr["qm_xyz"] / 0.529177    # bohr
        pc = fr["pc_xyz"] / 0.529177
        d = g[:, None, :] - pc[None]
        r = np.linalg.norm(d, axis=2)
        esp = (fr["pc_q"][None, :] / r).sum(1)
        E = (fr["pc_q"][None, :, None] * d / r[..., None] ** 3).sum(1)
        feats.append(np.concatenate([esp, (E ** 2).sum(1)]))
    return np.array(feats)


class GPModel:
    """Standartlaştırma + Gauss süreci regresyonu.
    Çekirdek = doğrusal (DotProduct) + RBF + beyaz gürültü: az veride doğrusal eğilimi,
    çok veride doğrusal olmayan ayrıntıyı yakalar; WhiteKernel etiketlerdeki gürültüyü emer."""

    def __init__(self, seed=0):
        self.seed = seed

    def fit(self, F, y):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel, DotProduct, WhiteKernel
        from sklearn.preprocessing import StandardScaler

        self.sc = StandardScaler().fit(F)
        k = (ConstantKernel(1.0) * RBF(10.0, (1e-1, 1e3))
             + ConstantKernel(0.1, (1e-4, 1e2)) * DotProduct(1.0, (1e-3, 1e3))
             + WhiteKernel(0.5))
        self.gp = GaussianProcessRegressor(kernel=k, normalize_y=True,
                                           n_restarts_optimizer=2, random_state=self.seed)
        self.gp.fit(self.sc.transform(F), np.asarray(y))
        return self

    def predict(self, F, return_std=False):
        m, s = self.gp.predict(self.sc.transform(F), return_std=True)
        return (m, s) if return_std else m


# ----------------------------------------------------------------------
# Serbest enerji yeniden ağırlıklandırma (MM -> QM)
# ----------------------------------------------------------------------
def zwanzig(dU, n_boot=500, seed=0):
    """ΔA = -kT ln < exp(-ΔU/kT) >_MM   ve ikinci mertebe kümülant yaklaşımı
       ΔA ≈ <ΔU> - σ²/(2kT). Bootstrap ile hata çubuğu."""
    rng = np.random.default_rng(seed)
    dU = np.asarray(dU, float)

    def exp_avg(x):
        m = x.min()
        return m - KT * np.log(np.mean(np.exp(-(x - m) / KT)))

    def cumulant(x):
        return x.mean() - x.var() / (2 * KT)

    boots_e, boots_c = [], []
    for _ in range(n_boot):
        s = rng.choice(dU, size=len(dU), replace=True)
        boots_e.append(exp_avg(s))
        boots_c.append(cumulant(s))
    return dict(exp=exp_avg(dU), exp_err=float(np.std(boots_e)),
                cumulant=cumulant(dU), cumulant_err=float(np.std(boots_c)),
                sigma=float(dU.std()), sigma_over_kT=float(dU.std() / KT))


def binding_correction(dU_complex, dU_solvent):
    """ΔΔG = ΔA_kompleks(MM→QM) - ΔA_çözücü(MM→QM).
    ΔU'daki sabit kayma (QM ve MM enerjilerinin farklı sıfır noktası) iki bacakta aynı
    olduğu için farkta sadeleşir."""
    off = np.mean(np.concatenate([dU_complex, dU_solvent]))
    zc = zwanzig(np.asarray(dU_complex) - off)
    zs = zwanzig(np.asarray(dU_solvent) - off)
    out = {}
    for key in ("exp", "cumulant"):
        out[key] = zc[key] - zs[key]
        out[key + "_err"] = float(np.hypot(zc[key + "_err"], zs[key + "_err"]))
    out["complex"] = zc
    out["solvent"] = zs
    return out


# ----------------------------------------------------------------------
# Aktif öğrenme döngüsü (ML1)
# ----------------------------------------------------------------------
def active_learning(labeler, feats, legs=("complex", "solvent"), n_init=8, batch=4,
                    rounds=3, n_val=6, seed=0, log=print):
    """Havuzdaki anlık görüntülerden, modelin EN EMİN OLMADIĞI yapıları seçip DFT ile
    etiketler. Her turda doğrulama kümesindeki hata izlenir.
    feats: {leg: (N_leg, d)} ; dönen model her iki bacak için ortak (makaledeki gibi)."""
    rng = np.random.default_rng(seed)
    pool = [(leg, i) for leg in legs for i in range(len(feats[leg]))]
    # doğrulama kümesi (rastgele, eğitimde kullanılmaz)
    val = []
    for leg in legs:
        idx = rng.choice(len(feats[leg]), size=n_val // len(legs), replace=False)
        val += [(leg, int(i)) for i in idx]
    remaining = [p for p in pool if p not in val]
    # başlangıç eğitim kümesi (rastgele)
    train = []
    for leg in legs:
        cand = [p for p in remaining if p[0] == leg]
        sel = rng.choice(len(cand), size=n_init // len(legs), replace=False)
        train += [cand[j] for j in sel]

    def F(lst):
        return np.array([feats[l][i] for l, i in lst])

    def Y(lst):
        return np.array([labeler.dU_dft(l, i) for l, i in lst])

    y_val = Y(val)
    history = []
    model = None
    for r in range(rounds + 1):
        y_tr = Y(train)
        model = GPModel(seed=seed).fit(F(train), y_tr)
        pv = model.predict(F(val))
        rmse = float(np.sqrt(np.mean((pv - y_val) ** 2)))
        history.append(dict(round=r, n_train=len(train), rmse_val=rmse,
                            std_target=float(np.std(np.concatenate([y_tr, y_val])))))
        log(f"  tur {r}: eğitim={len(train):3d}  doğrulama RMSE={rmse:6.2f} kJ/mol "
            f"(hedefin std'si {history[-1]['std_target']:.2f})")
        if r == rounds:
            break
        cand = [p for p in remaining if p not in train]
        _, s = model.predict(F(cand), return_std=True)
        pick = np.argsort(-s)[:batch]
        train += [cand[j] for j in pick]
    return model, train, val, history
