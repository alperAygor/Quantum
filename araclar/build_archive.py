#!/usr/bin/env python3
# Proje_Arsivi olusturucu: tum dosyalari klasorlere ayirir, veri CSV'lerini, molekul dosyalarini,
# grafikleri ve gorsel envanterini uretir, sonra zip'ler.
import os, re, json, csv, shutil, glob, tarfile, html
import numpy as np

ROOT = "/home/claude/Proje_Arsivi"
SP = "/tmp/claude-0/-home-claude/7fd9d8ae-bc1f-52fd-9764-fda5ed1b2771/scratchpad"
OUT = "/mnt/user-data/outputs"
UP = "/mnt/user-data/uploads/PDT_gorsel_HQ"
SRC = "/home/claude/arsiv_src"
if os.path.exists(ROOT): shutil.rmtree(ROOT)
D = lambda *p: os.makedirs(os.path.join(ROOT, *p), exist_ok=True) or os.path.join(ROOT, *p)
def cp(src, *dst):
    d = os.path.join(ROOT, *dst); os.makedirs(os.path.dirname(d), exist_ok=True); shutil.copy2(src, d)

# ------------------------------------------------------------------ docs
D("docs")
cp(f"{SRC}/README.md", "README.md")
cp(f"{SRC}/References.md", "docs", "References.md")
DOCS = {"01_BCD_PMF_sonuc_ozeti.md": "PMF_sonuc_ozeti.md", "02_PDT_porfin_kuantum.md": "PDT_porfin_kuantum.md",
        "03_Tarama_v3_sonuclar.md": "Tarama_v3_sonuclar.md", "04_Literatur_karsilastirma_PDT_kuantum.md": "Literatur_karsilastirma_PDT_kuantum.md",
        "05_Basvuru_v4_durum.md": "Basvuru_v4_durum.md", "06_IBM_donanim_plani.md": "IBM_donanim_plani.md"}
for dst, src in DOCS.items():
    cp(f"{SRC}/project_docs/{src}", "docs", dst)

# Sunum metni: deck slaytlarindan metin + konusmaci notlari
deck = f"{SP}/deck/project"
order = json.load(open(f"{deck}/deck.json"))["order"]
L = ["# Sunum Metni ve Konuşmacı Notları", "",
     "Kaynak: Claude Slides artifact'i \"Kuantum Triyaj – YTÜ Startup House Başvuru Sunumu\" (17 slayt, 20 Ekim 2026).",
     "Not: Ekip slaytındaki yer tutucular ([bölüm], [Ekip arkadaşı · rolü]) kullanıcı tarafından editörde doldurulmuş olabilir; buradaki metin Claude'un son yazdığı sürümdür.", ""]
for i, sid in enumerate(order, 1):
    h = open(f"{deck}/slides/{sid}.html", encoding="utf-8").read()
    note = re.search(r"<aside>(.*?)</aside>", h, re.S); body = re.sub(r"<aside>.*?</aside>", "", h, flags=re.S)
    body = re.sub(r"<br\s*/?>", " / ", body)
    body = re.sub(r"</(p|h1|h2|h3|li|tr|th|td)>", r"</\1>\n", body)
    txt = [html.unescape(re.sub(r"<[^>]+>", " ", l)).strip() for l in body.split("\n")]
    txt = [re.sub(r"\s+", " ", t) for t in txt if t.strip()]
    L += [f"## Slayt {i} — `{sid}`", ""] + [f"- {t}" for t in txt] + [""]
    if note: L += ["**Konuşmacı notu:** " + html.unescape(note.group(1).strip()), ""]
open(os.path.join(ROOT, "docs", "Sunum_metni_ve_konusmaci_notlari.md"), "w", encoding="utf-8").write("\n".join(L))

# ------------------------------------------------------------------ data
D("data", "raw")
cp(f"{SRC}/Veriler_ve_Kuresel_Muadil_Analizi.md", "data", "Veriler_ve_Kuresel_Muadil_Analizi.md")
V3 = [  # molekul, S1_CAS, S1_NEVPT2, deney, deney_tipi, deney_kaynak, f_S1, w1, w2, triyaj, T1_dSCF
 ("porfin", 3.62, 1.87, 1.98, "dogrudan", "Edwards 1971 (buhar, Qx)", 0.033, 0.544, 0.435, "KUANTUM GEREKLI", 2.00),
 ("klorin", 3.45, 1.98, 1.95, "dogrudan (~)", "Lenth/Baumert (n-alkan, 4.2 K)", 0.579, 0.743, 0.098, "klasik yeterli", 1.53),
 ("bakteriyoklorin", 3.59, 2.23, 1.74, "analog", "Reddy/Lindsey 2017, B6 (toluen)", 0.422, 0.748, 0.168, "klasik yeterli (model yetersiz)", 1.65),
 ("Zn-porfin", 3.75, 2.07, 2.18, "dogrudan", "Baerends 2002 (cozelti)", 0.079, 0.616, 0.367, "KUANTUM GEREKLI", 2.08),
 ("Mg-porfin", 3.72, 2.04, 2.17, "dogrudan", "Devillers 2013 (DMF)", 0.073, 0.611, 0.373, "KUANTUM GEREKLI", 2.05),
 ("Pd-porfin", 3.80, 2.14, None, "yok", "-", 0.085, 0.629, 0.355, "KUANTUM GEREKLI", 2.16),
 ("porfirazin", 3.74, 2.24, 2.065, "dogrudan", "Henchy 2014 (Ne matris, 8 K)", 0.633, 0.788, 0.139, "klasik yeterli", 1.65),
 ("tetrametilporfin", 3.42, 1.64, 1.88, "vekil", "Paliteiro 2005 (meso-tetraalkil, CH2Cl2)", 0.001, 0.507, 0.470, "KUANTUM GEREKLI (geometri supheli)", 1.58)]
with open(os.path.join(ROOT, "data", "v3_tarama_sonuclari.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["molekul", "S1_CAS44_eV", "S1_CAS44_NEVPT2_eV", "deney_S1_eV", "deney_tipi", "deney_kaynak", "fark_eV",
                                   "f_S1", "konfig_agirlik_1", "konfig_agirlik_2", "triyaj", "T1_dSCF_eV", "S1_nm"])
    for r in V3:
        w.writerow([r[0], r[1], r[2], r[3] if r[3] else "", r[4], r[5], "" if r[3] is None else round(r[2] - r[3], 3), r[6], r[7], r[8], r[9], r[10], round(1239.84 / r[2], 1)])
TD = [("porfin", 2.348, 0.001, 1.832, 1.98, "tamamlandi"), ("klorin", 2.369, 0.174, 1.540, 1.95, "tamamlandi"),
      ("Zn-porfin", None, None, None, 2.18, "tamamlandi; deger ~/tarama/Zn-porfin_tddft.json icinde, arsive aktarilmadi"),
      ("Mg-porfin", None, None, None, 2.17, "tamamlandi; deger ~/tarama/Mg-porfin_tddft.json icinde, arsive aktarilmadi"),
      ("porfirazin", None, None, None, 2.065, "tamamlandi; deger ~/tarama/porfirazin_tddft.json icinde, arsive aktarilmadi"),
      ("bakteriyoklorin", None, None, None, 1.74, "tamamlandi; deger ~/tarama/bakteriyoklorin_tddft.json icinde, arsive aktarilmadi"),
      ("tetrametilporfin", None, None, None, 1.88, "hesap suruyor"), ("Pd-porfin", None, None, None, None, "hesap suruyor")]
with open(os.path.join(ROOT, "data", "tddft_b3lyp.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["molekul", "TDDFT_B3LYP_S1_eV", "f_S1", "TDA_B3LYP_T1_eV", "deney_S1_eV", "fark_eV", "S1_nm", "durum"])
    for m, s, fs, t, e, du in TD:
        w.writerow([m, s or "", fs if fs is not None else "", t or "", e or "", "" if (s is None or e is None) else round(s - e, 3), "" if s is None else round(1239.84 / s, 1), du])
with open(os.path.join(ROOT, "data", "ibm_donanim_porfin.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["hesap", "S1_eV", "tam_cozumden_fark_eV", "not"])
    for r in [("Tam cozum (CASCI, klasik)", 3.617, 0.0, "referans"), ("VQE ideal simulator", 3.616, -0.001, "statevector"),
              ("VQE FakeKingston gurultu sim. ham", 2.825, -0.792, "sandbox, qiskit-aer"), ("VQE FakeKingston duzeltilmis", 4.095, 0.478, "sandbox"),
              ("VQE ibm_kingston ham", 1.141, -2.477, "is dbbl9ns24mrc73chvt8g, 20.10.2026 14:47"),
              ("VQE ibm_kingston duzeltilmis", 4.546, 0.928, "istatistik belirsizlik 0.190 eV; gurultu faktoru f(S0,S1)=0.41,0.53")]:
        w.writerow(r)
with open(os.path.join(ROOT, "data", "klinik_ilaclar.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["ilac", "ticari_ad", "cekirdek", "atom_sayisi", "klinik_isik_nm", "isik_kaynagi", "hizalama_RMSD_A",
                                   "bizim_cekirdek_S1_eV", "bizim_nm", "fark_nm", "TDDFT_cekirdek_S1_eV", "TDDFT_nm", "TDDFT_fark_nm", "triyaj"])
    for r in [("Temoporfin", "Foscan", "klorin", 84, 652, "EMA urun bilgisi", 0.40, 1.98, 2.369, "klasik yeterli"),
              ("Talaporfin", "Laserphyrin", "klorin", 93, 664, "Wikipedia (664-667 nm)", 0.24, 1.98, 2.369, "klasik yeterli"),
              ("Protoporfirin IX", "5-ALA / Ameluz", "porfin", 76, 635, "FDA Ameluz etiketi 2017", 0.22, 1.87, 2.348, "KUANTUM GEREKLI"),
              ("Padeliporfin", "Tookad", "bakteriyoklorin", "", 753, "emc SmPC", "", 2.23, None, "model yetersiz")]:
        nm = 1239.84 / r[7]; tn = None if r[8] is None else 1239.84 / r[8]
        w.writerow([r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], round(nm, 1), round(abs(r[4] - nm), 1),
                    r[8] or "", "" if tn is None else round(tn, 1), "" if tn is None else round(abs(r[4] - tn), 1), r[9]])
with open(os.path.join(ROOT, "data", "bcd_sonuclar.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["buyukluk", "deger", "birim", "not"])
    for r in [("dG_PMF_WHAM", "-2.7 +/- 0.6", "kcal/mol", "37 pencere x 2 ns, 300 K"), ("dG_deney_ITC", "-3.05", "kcal/mol", "Rekharsky 1997; dH -1.1, TdS +2.0"),
              ("mod_A_orani", "53.5", "%", "OH dar/primer kenar"), ("mod_B_orani", "38.7", "%", "OH genis/sekonder kenar"), ("gecis_orani", "7.8", "%", ""),
              ("Eint_gaz_tumu", "-9.34 +/- 1.63", "kcal/mol", "GFN2-xTB, 31 kare"), ("Eint_su_tumu", "-6.33 +/- 1.45", "kcal/mol", "ALPB"),
              ("D4_dispersiyon_tumu", "-9.00 +/- 0.84", "kcal/mol", "~%96 dispersiyon"), ("MFCC_hata", "-0.18 +/- 0.13", "kcal/mol", "tam -9.20'ye gore"),
              ("bolge_HF", "+1.39", "kcal/mol", "39 atom, def2-SVP, CP"), ("bolge_B3LYP", "+0.73", "kcal/mol", ""), ("bolge_MP2", "-1.32", "kcal/mol", ""),
              ("bolge_B3LYP_D3BJ", "-2.83", "kcal/mol", ""), ("bolge_GFN2xTB", "-2.80", "kcal/mol", ""),
              ("kuantum_duzeltmesi_VQE_eksi_MP2", "+0.0009", "kcal/mol", "aktif uzay 4e,4o"), ("hibrit_Eint", "-1.320", "kcal/mol", "MP2 -1.321")]:
        w.writerow(r)
with open(os.path.join(ROOT, "data", "porfin_uyarilmis_haller.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["hal", "CAS44_VQD_eV", "NEVPT2_eV", "f", "baskin_konfigurasyonlar"])
    for r in [("S0", 0, 0, "", "HF 0.93"), ("S1 (Qx)", 3.62, 1.87, 0.033, "H->L 0.54; H-1->L+1 0.44"), ("S2 (Qy)", 3.74, "", 0.068, "H-1->L 0.59; H->L+1 0.39"),
              ("S3 (By)", 5.35, "", 2.13, "H->L+1 0.47; H-1->L 0.31"), ("S4 (Bx)", 5.52, "", 2.76, "H->L 0.45; H-1->L+1 0.37"),
              ("T1", 2.71, "1.95 (guvenilmez)", "", "H->L 0.70"), ("T2", 2.91, "", "", "H-1->L+1 0.63")]:
        w.writerow(r)
with open(os.path.join(ROOT, "data", "literatur_kiyas.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["calisma", "molekul_ailesi", "yontem", "kubit", "donanim", "deneye_MAE_eV", "not"])
    for r in [("Bu proje", "porfirin (8)", "CAS(4,4) VQD + NEVPT2", 8, "simulator + ibm_kingston", "0.111 (5 mol) / 0.184 (7 mol)", "TD-DFT B3LYP 0.39 vs 0.07 (porfin+klorin)"),
              ("Algorithmiq dADAPT 2024", "BODIPY (6)", "dUCCSD(8,8)", 16, "gurultusuz simulator", "0.131", "klasik dCCSD 0.054; TD-DFT 0.199/0.368/0.392; EOM-CCSD 0.498"),
              ("Algorithmiq Q4Bio 2026", "TLD1433 Ru(II)", "ADAPT-VMPE + QB-DMRG", "52-100", "ibm_cleveland", "-", "S0 enerjisi; deney kiyasi yok"),
              ("Xanadu 2025", "BODIPY (4 sinif)", "QSP esik projeksiyonu, ISC proxy", "180-350 mantiksal", "yok (FTQC)", "-", "kaynak tahmini"),
              ("JCTC 2019", "porfin", "CASSCF/NEVPT2, CASPT2, MRCI", "-", "klasik", "Qx: NEVPT2 1.96, CASPT2 1.70, MRCI 1.90 (deney 1.98-2.02)", "")]:
        w.writerow(r)
for f_ in glob.glob("/home/claude/figtest/*_v3_sonuc.json"): cp(f_, "data", "raw", "v3_json", os.path.basename(f_))
for f_ in glob.glob("/home/claude/figtest/*_tddft.json"): cp(f_, "data", "raw", "tddft_json", os.path.basename(f_))
cp("/home/claude/figtest/final_tablo.csv", "data", "final_tablo_F1F3.csv")
cp("/home/claude/figtest/final_ozet.txt", "data", "final_ozet.txt")
for f_ in ["porfin_ibm_hazir.json", "porfin_ibm_sonuc_sahte.json", "hazirla.log", "gonder.log", "sonuc.log"]:
    if os.path.exists(f"/home/claude/ibm/{f_}"): cp(f"/home/claude/ibm/{f_}", "data", "raw", "ibm_sandbox_testi", f_)
for f_, n in [(f"{SP}/pmftest/profil.xvg", "pmf_profil.xvg"), (f"{SP}/histo.xvg", "umbrella_histogramlari.xvg"), (f"{SP}/profil.xvg", "pmf_profil_tam.xvg"),
              (f"{SP}/analiz.csv", "md_konak_misafir_analiz.csv"), (f"{SP}/pmftest/pencereler.txt", "umbrella_pencereler.txt"),
              (f"{SP}/cdtest/fragman_xtb.txt", "fragman_xtb_cikti.txt"), (f"{SP}/pdt/log3.txt", "porfin_qm_log.txt"),
              (f"{SP}/pdt/t1log.txt", "porfin_T1_log.txt"), (f"{SP}/test_vqe_4e4o.txt", "vqe_4e4o_test.txt")]:
    if os.path.exists(f_): cp(f_, "data", "raw", n)

# ------------------------------------------------------------------ molecules
from rdkit import Chem
from rdkit.Chem import AllChem
D("molecules", "pdb"); D("molecules", "xyz"); D("molecules", "sdf"); D("molecules", "orbitals"); D("molecules", "bcd_siklodekstrin")
SMI = [("porfin", "cekirdek", "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2", "geo.py"),
       ("klorin", "cekirdek", "C1CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2", "geo.py"),
       ("bakteriyoklorin", "cekirdek", "C1CC2=CC3=NC(=CC4=CC=C(N4)C=C5CCC(=N5)C=C1N2)C=C3", "geo.py"),
       ("tetrametilporfin", "cekirdek", "CC1=C2C=CC(=C(C)C3=CC=C(N3)C(C)=C4C=CC(=N4)C(C)=C5C=CC1=N5)N2", "geo.py"),
       ("porfirazin", "cekirdek", "C1=CC2=NC3=CC=C(N3)N=C4C=CC(=N4)N=C5C=CC(=N5)N=C1N2", "geo.py"),
       ("Zn-porfin", "metal cekirdek", "porfin iskeleti, ic NH hidrojenleri cikarilip merkeze Zn (geo.py)", "geo.py"),
       ("Mg-porfin", "metal cekirdek", "porfin iskeleti, ic NH hidrojenleri cikarilip merkeze Mg (geo.py)", "geo.py"),
       ("Pd-porfin", "metal cekirdek", "porfin iskeleti, ic NH hidrojenleri cikarilip merkeze Pd (geo.py)", "geo.py"),
       ("Temoporfin", "klinik ilac (Foscan), PubChem CID 60751", "C1CC2=NC1=C(C3=CC=C(N3)C(=C4C=CC(=N4)C(=C5C=CC(=C2C6=CC(=CC=C6)O)N5)C7=CC(=CC=C7)O)C8=CC(=CC=C8)O)C9=CC(=CC=C9)O", "gercek_ilac.py"),
       ("Talaporfin", "klinik ilac (Laserphyrin), PubChem CID 5486799", "CCC1=C2C=C3C(=C(C(=N3)C(=C4[C@H]([C@@H](C(=N4)C=C5C(=C(C(=CC(=C1C)N2)N5)C=C)C)C)CCC(=O)O)CC(=O)N[C@@H](CC(=O)O)C(=O)O)C(=O)O)C", "gercek_ilac.py"),
       ("Protoporfirin-IX", "klinik (5-ALA ile olusur), PubChem CID 4971; SMILES elle yazildi", "C=CC1=C(C)C2=CC3=NC(=CC4=C(C)C(CCC(=O)O)=C(N4)C=C5N=C(C=C1N2)C(C)=C5CCC(=O)O)C(C)=C3C=C", "gercek_ilac.py"),
       ("siklopentanol", "misafir (ATB molid 859)", "OC1CCCC1", "GROMACS/ATB")]
with open(os.path.join(ROOT, "molecules", "SMILES_listesi.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["molekul", "tur", "SMILES_veya_tanim", "uretildigi_betik"]); [w.writerow(r) for r in SMI]
EL = {"C": 6, "N": 7, "O": 8, "H": 1, "Zn": 30, "Mg": 12, "Pd": 46}
def pdb_atoms(fn):
    at = []
    for l in open(fn):
        if l.startswith(("ATOM", "HETATM")):
            el = l[76:78].strip() or re.sub(r"[^A-Za-z]", "", l[12:16])[:2].strip()
            el = el.capitalize(); at.append((el, float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return at
pdbs = sorted(glob.glob(f"{UP}/*/*.pdb"))
for p in pdbs:
    n = os.path.basename(p)[:-4]; cp(p, "molecules", "pdb", os.path.basename(p))
    at = pdb_atoms(p)
    with open(os.path.join(ROOT, "molecules", "xyz", n + ".xyz"), "w") as f:
        f.write(f"{len(at)}\n{n} (GFN2-xTB geometrisi{'; ilac cekirdege Kabsch ile hizalandi' if n in ('Temoporfin','Talaporfin','Protoporfirin-IX') else ''})\n")
        for el, x, y, z in at: f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    m = Chem.MolFromPDBFile(p, removeHs=False, sanitize=False, proximityBonding=True)
    if m is not None:
        try:
            Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass
        m.SetProp("_Name", n)
        w_ = Chem.SDWriter(os.path.join(ROOT, "molecules", "sdf", n + ".sdf")); w_.SetKekulize(False); w_.write(m); w_.close()
for n, t, s, _ in SMI:
    if "iskeleti" in s or n in [os.path.basename(p)[:-4] for p in pdbs]: continue
    m = Chem.AddHs(Chem.MolFromSmiles(s)); AllChem.EmbedMolecule(m, randomSeed=7); AllChem.MMFFOptimizeMolecule(m); m.SetProp("_Name", n + " (RDKit+MMFF)")
    Chem.MolToMolFile(m, os.path.join(ROOT, "molecules", "sdf", n + "_rdkit.sdf"))
cp(f"{OUT}/porfin.xyz", "molecules", "xyz", "porfin_sandbox_xtb.xyz")
cp(f"{OUT}/porfin_orbitaller.molden", "molecules", "orbitals", "porfin_orbitaller.molden")
for f_, n in [(f"{SP}/bcd_345755.itp", "bcd_345755.itp"), (f"{SP}/cyclopentanol_aa.itp", "siklopentanol_aa.itp"), (f"{SP}/complex.gro", "kompleks_baslangic.gro"),
              (f"{SP}/view.gro", "md_yorunge_view.gro"), (f"{SP}/cdtest/mod_A.xyz", "mod_A_temsilci_kare.xyz"), (f"{SP}/cdtest/mod_A_bolge_G1.xyz", "mod_A_bolge_glukoz.xyz"),
              (f"{SP}/guest_c.gro", "misafir.gro"), (f"{SP}/cd_c.gro", "beta_CD.gro")]:
    if os.path.exists(f_): cp(f_, "molecules", "bcd_siklodekstrin", n)

# ------------------------------------------------------------------ visuals
D("visuals", "grafikler"); D("visuals", "chimerax"); D("visuals", "chimerax_betikleri"); D("visuals", "sunum_kirpilmis")
GRAF = [("fig01_PMF_bCD_siklopentanol.png", f"{OUT}/PMF_bCD_siklopentanol.png"),
        ("fig02_porfin_gouterman_orbitaller.png", f"{OUT}/porfin_gouterman_orbitaller.png"),
        ("fig03_porfin_orbital_secimi.png", f"{OUT}/porfin_orbital_secimi.png"),
        ("fig04_tarama_v3_sonuclari.png", f"{OUT}/tarama_sonuclari.png"),
        ("fig05_porfin_panel_2_3.png", f"{OUT}/porfin_panel_2_3.png"),
        ("fig06_F1_yontem_karsilastirma.png", "/home/claude/figtest/F1_yontem_karsilastirma.png"),
        ("fig07_F2_triyaj_dogrulama.png", "/home/claude/figtest/F2_triyaj_dogrulama.png"),
        ("fig08_F3_klinik_ilaclar.png", "/home/claude/figtest/F3_klinik_ilaclar.png")]
for n, s in GRAF: cp(s, "visuals", "grafikler", n)
PAN = sorted(glob.glob(f"{OUT}/paneller/*.png"))
for i, p in enumerate(PAN, 10): cp(p, "visuals", "grafikler", f"fig{i:02d}_panel_{os.path.basename(p)}")
for p in glob.glob("/home/claude/deckimg/*.jpg"): cp(p, "visuals", "sunum_kirpilmis", os.path.basename(p))
for p in glob.glob(f"{UP}/**/*.cxc", recursive=True):
    rel = os.path.relpath(p, UP).replace(os.sep, "__"); cp(p, "visuals", "chimerax_betikleri", rel)

# IBM donanim grafigi (fig09)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6a64", "#e4e3de"
rows = [("Tam çözüm (CASCI)", 3.617, INK), ("VQE, ideal simülatör", 3.616, BLUE), ("Gürültü sim., ham", 2.825, "#9db8e0"),
        ("Gürültü sim., düzeltilmiş", 4.095, "#9db8e0"), ("ibm_kingston, ham", 1.141, ORANGE), ("ibm_kingston, düzeltilmiş", 4.546, ORANGE)]
fig, ax = plt.subplots(figsize=(11, 5.2))
for i, (n, v, c) in enumerate(rows[::-1]):
    ax.barh(i, v, color=c, height=0.6); ax.text(v + (0.27 if "düzeltilmiş" in n and "kingston" in n else 0.05), i, f"{v:.2f} eV", va="center", fontsize=11, color=INK)
ax.errorbar(4.546, 0, xerr=0.19, color=INK, capsize=5)
ax.axvline(3.617, color=INK, ls="--", lw=1); ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=11)
ax.set_xlabel("Porfin S1 uyarılma enerjisi, 8 kübitlik kuantum bölgesi (eV)"); ax.set_xlim(0, 5.6)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
ax.set_title("IBM Kingston (156 kübit) gerçek donanım ölçümü — iş dbbl9ns24mrc73chvt8g, 20.10.2026", fontsize=12, color=INK)
fig.text(0.01, -0.02, "Kesikli çizgi: tam çözüm. Hata çubuğu: istatistik belirsizlik ±0,19 eV. Gürültü faktörü f(S0,S1) = 0,41 / 0,53.", fontsize=9, color=MUTED)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "visuals", "grafikler", "fig09_ibm_kingston_donanim.png"), dpi=250, bbox_inches="tight"); plt.close(fig)

# ChimeraX goruntuleri: standart adlar
MOL = {"porfin": ("Porfin", 1.87, 0.033, "KUANTUM GEREKLİ", 1.98, (0.544, 0.435)),
       "Zn-porfin": ("Zn-porfin", 2.07, 0.079, "KUANTUM GEREKLİ", 2.18, (0.616, 0.367)),
       "Mg-porfin": ("Mg-porfin", 2.04, 0.073, "KUANTUM GEREKLİ", 2.17, (0.611, 0.373)),
       "Pd-porfin": ("Pd-porfin", 2.14, 0.085, "KUANTUM GEREKLİ", None, (0.629, 0.355)),
       "tetrametilporfin": ("Tetrametilporfin", 1.64, 0.001, "KUANTUM GEREKLİ", 1.88, (0.507, 0.47)),
       "klorin": ("Klorin", 1.98, 0.579, "klasik yeterli", 1.95, (0.743, 0.098)),
       "porfirazin": ("Porfirazin", 2.24, 0.633, "klasik yeterli", 2.065, (0.788, 0.139)),
       "bakteriyoklorin": ("Bakteriyoklorin", 2.23, 0.422, "klasik yeterli (model yetersiz)", 1.74, (0.748, 0.168))}
ATOMS = {"porfin": "C20H14N4, 38 atom", "Zn-porfin": "C20H12N4Zn, 37 atom", "Mg-porfin": "C20H12N4Mg, 37 atom", "Pd-porfin": "C20H12N4Pd, 37 atom",
         "tetrametilporfin": "C24H22N4, 50 atom", "klorin": "C20H16N4, 40 atom", "porfirazin": "C16H10N8, 34 atom", "bakteriyoklorin": "C20H18N4, 42 atom"}
TUR = [("molekul", "uzaydolgu"), ("molekul_cubuk", "topcubuk"), ("HOMO", "HOMO"), ("LUMO", "LUMO"), ("fark", "S1_S0_fark_yogunlugu")]
inv = []; k = 21
for key in ["porfin", "Zn-porfin", "Mg-porfin", "Pd-porfin", "tetrametilporfin", "klorin", "porfirazin", "bakteriyoklorin"]:
    for t, ad in TUR:
        src = f"{UP}/Gorseller/{key}_{t}.png"
        if not os.path.exists(src): continue
        n = f"fig{k:02d}_molekul_{key}_{ad}.png"; cp(src, "visuals", "chimerax", n); inv.append((n, key, t, None)); k += 1
for drug, core, files in [("Temoporfin", "klorin", [("klorin", "Temoporfin_molekul"), ("klorin", "Temoporfin_HOMO"), ("klorin", "Temoporfin_LUMO")]),
                          ("Talaporfin", "klorin", [("klorin", "Talaporfin_molekul"), ("klorin", "Talaporfin_HOMO"), ("klorin", "Talaporfin_LUMO")]),
                          ("Protoporfirin-IX", "porfin", [("porfin", "Protoporfirin-IX_molekul"), ("porfin", "Protoporfirin-IX_HOMO"), ("porfin", "Protoporfirin-IX_LUMO")])]:
    for d, f_ in files:
        src = f"{UP}/{d}/{f_}.png"
        if not os.path.exists(src): continue
        t = f_.split("_")[-1]; n = f"fig{k:02d}_ilac_{drug}_{ {'molekul': 'uzaydolgu', 'HOMO': 'HOMO_kuantum_bolgesi', 'LUMO': 'LUMO_kuantum_bolgesi'}[t] }.png".replace(" ", "")
        cp(src, "visuals", "chimerax", n); inv.append((n, drug, t, core)); k += 1

# ------------------------------------------------------------------ Gorsel envanteri
DRUGINFO = {"Temoporfin": ("Temoporfin (mTHPC, Foscan®)", "C44H32N4O4, 84 atom", "AB'de onaylı; baş-boyun kanseri PDT; 652 nm ışık", 0.40),
            "Talaporfin": ("Talaporfin sodyum (NPe6, Laserphyrin®)", "93 atom (nötr asit formu)", "Japonya'da onaylı; akciğer kanseri PDT; 664 nm ışık", 0.24),
            "Protoporfirin-IX": ("Protoporfirin IX (PpIX)", "C34H34N4O4, 76 atom", "5-ALA/MAL ön ilacından vücutta oluşur; Ameluz® + 635 nm ışık", 0.22)}
G = ["# Görsel Envanteri ve Rehberi", "",
     "Bu dosya arşivdeki tüm görselleri standart adlarla listeler. Her görsel için şunlar verilmiştir: (a) neyi gösterdiği ve bilimsel açıklaması, (b) üretildiği araç ve görüntüleme parametreleri, (c) rapor veya sunumda kullanılacak hazır şekil açıklaması (caption).", "",
     "**Klasörler:**",
     "- `visuals/grafikler/`: Python/matplotlib grafikleri (fig01–fig20).",
     "- `visuals/chimerax/`: UCSF ChimeraX ile üretilmiş 3B moleküler görseller (fig21 ve sonrası; 3000×2000 piksel, supersample 4).",
     "- `visuals/chimerax_betikleri/`: Görselleri yeniden üreten `.cxc` betikleri. Dosya adlarındaki `__` işareti klasör ayırıcısıdır.",
     "- `visuals/sunum_kirpilmis/`: Sunum için kırpılmış 1200×1200 kopyalar.", "",
     "## Ortak görüntüleme parametreleri (ChimeraX)", "",
     "| Görsel türü | Stil | Renk kodları | Işık ve kalite | Diğer |",
     "|---|---|---|---|---|",
     "| Uzay dolgu (molekül) | `style sphere` | `color byelement`; C #8a8a8a, N mavi, O kırmızı, H beyaz, metal element rengi | `lighting soft`, `shadows true intensity 0.6`, `multiShadow 1024`, `material specular 0.4 exponent 40`, `graphics quality 6` | `graphics silhouettes true width 2 color #2a2a2a`; `turn x -35` (ilaçlarda -40), `turn z 15`; beyaz arka plan |",
     "| Top-çubuk | `style ball`, `size stickRadius 0.14`, `size ballScale 0.25` | aynı | aynı | aynı |",
     "| Orbital (HOMO/LUMO) | `style ball`, `stickRadius 0.11`, `ballScale 0.22` (ilaçlarda 0.13 / 0.22) | C #909090, N #3050f8, H #f4f4f4; izoyüzey: negatif faz #2f6fd0 (mavi), pozitif faz #f0782f (turuncu) | `lighting full`, `shadows 0.5`, `multiShadow 512`, `graphics quality 5`, `camera ortho` | `volume ... style surface level ±0.02 transparency 0.1` (ilaçlarda 0.25), `surfaceSmoothing true smoothingIterations 4`; küp ızgarası 0,1 Å, kenar 7 bohr (PySCF cubegen) |",
     "| Fark yoğunluğu (S1−S0) | Orbital ile aynı | turuncu (#f0782f) = elektronun geldiği yer (+), mavi (#2f6fd0) = ayrıldığı yer (−) | aynı | Seviye otomatik: ±0,12 × max\\|Δρ\\| (porfinde ±0,00058) |",
     "", "Görsellerin üst ve alt kenarlarındaki etiketler `2dlabels` komutuyla eklendi. Etiketlerde molekül adı, S1 (eV ve nm), f ve triyaj kararı bulunur.", "",
     "## A. Python/matplotlib grafikleri (`visuals/grafikler/`)", ""]
GDESC = {
 "fig01_PMF_bCD_siklopentanol.png": ("β-CD + siklopentanol için ortalama kuvvet potansiyeli (PMF) W(ξ). Alt panelde 37 şemsiye penceresinin histogramları (örtüşme kontrolü). Bağlı durum ξ ≈ −0,05 nm'de −17,7 kJ/mol.", "GROMACS 2026.1 `gmx wham` (-b 200) → `plot_pmf.py` (matplotlib). Mavi #2a78d6 çizgi, kesikli çizgiler su platoları.", "Şekil 1. β-siklodekstrin–siklopentanol kompleksinin CD ekseni boyunca ortalama kuvvet potansiyeli (37 pencere × 2 ns, 300 K, WHAM). Hesaplanan ΔG° = −2,7 ± 0,6 kcal/mol; deney −3,05 kcal/mol."),
 "fig02_porfin_gouterman_orbitaller.png": ("Porfinin Gouterman dört orbitali (HOMO-1, HOMO, LUMO, LUMO+1) ve enerjileri.", "PySCF B3LYP/def2-SVP → molden → matplotlib/ilk görselleştirme (porfin_qm.py / ciz.py).", "Şekil 2. Porfinin Gouterman orbitalleri (B3LYP/def2-SVP). HOMO-1/HOMO (−5,53/−5,40 eV) ve LUMO/LUMO+1 (−2,46/−2,40 eV) neredeyse dejenere; bu, güçlü konfigürasyon karışımının kaynağıdır."),
 "fig03_porfin_orbital_secimi.png": ("Porfinin 406 orbitallik enerji merdiveni (sol) ve seçilen 4 orbitallik kuantum bölgesi (sağ, yakınlaştırma).", "matplotlib (`merdiven.py`). Turuncu = boş orbitaller, mavi = dolu orbitaller.", "Şekil 3. Kuantum bölgesinin seçimi: 406 orbitalden ışığı belirleyen 4 Gouterman orbitali (8 kübit) kuantum bilgisayara, kalan 402 orbital klasik NEVPT2'ye verilir."),
 "fig04_tarama_v3_sonuclari.png": ("8 porfirin türevi için hesap–deney karşılaştırması (dumbbell), osilatör gücü f ve triyaj kararı.", "matplotlib (tarama sonuç betiği), dataviz paleti (mavi #2a78d6, turuncu #eb6834).", "Şekil 4. Kuantum triyaj taraması (v3): CAS(4,4)-VQD+NEVPT2 ile S1 enerjileri, deney değerleri ve triyaj sınıflandırması."),
 "fig05_porfin_panel_2_3.png": ("Porfin için enerji/soğurma paneli ve triyaj yorumu (sunum paneli 2–3).", "matplotlib (`panel.py`).", "Şekil 5. Porfin: kuantum bölgesi sonuçları (S1 = 1,87 eV, f = 0,033) ve triyaj kararı (kuantum gerekli)."),
 "fig06_F1_yontem_karsilastirma.png": ("Sol: 8 molekülde deney, bizim yöntem ve TD-DFT B3LYP S1 değerleri. Sağ: ortak moleküllerde (porfin, klorin) ortalama mutlak hata (0,07 vs 0,39 eV).", "`final_figurler.py` (matplotlib, 220 dpi). Siyah elmas = deney, mavi daire = bizim yöntem, turuncu kare = TD-DFT.", "Şekil 6. Kuantum-hibrit yöntem (VQE 8 kübit + NEVPT2) ile sektör standardı TD-DFT B3LYP'nin deneye göre karşılaştırması. Porfin ve klorinde ortalama mutlak hata: 0,07 eV ve 0,39 eV."),
 "fig07_F2_triyaj_dogrulama.png": ("x ekseni: S1'de ikinci konfigürasyonun ağırlığı (triyaj ölçütü); y ekseni: deneye göre mutlak hata. Eşik 0,2.", "`final_figurler.py`.", "Şekil 7. Triyaj ölçütü ile hesap hatası arasındaki ilişki. Eşiğin (0,2) sağı 'kuantum gerekli' bölgesidir."),
 "fig08_F3_klinik_ilaclar.png": ("Dört klinik PDT ilacı için klinikte kullanılan ışık (siyah), bizim çekirdek hesabı (mavi) ve TD-DFT çekirdek hesabı (turuncu); 700–900 nm doku geçirgenlik penceresi.", "`final_figurler.py`; x ekseni 500–800 nm.", "Şekil 8. Klinik ilaçların fotoaktif çekirdekleri: bizim hesabımız klinik ışığa 26–38 nm, TD-DFT 107–141 nm uzaklıkta. Yakın kızılötesi ilaç padeliporfinde (753 nm) 4 orbitalli model yetersiz."),
 "fig09_ibm_kingston_donanim.png": ("Porfin S1 için tam çözüm, ideal VQE, gürültü simülasyonu ve gerçek IBM Kingston donanım sonuçları (ham ve düzeltilmiş).", "matplotlib (build_archive.py); veriler `data/ibm_donanim_porfin.csv`.", "Şekil 9. VQE devresinin IBM Kingston (156 kübit) gerçek kuantum işlemcisindeki ilk ölçümü. Ham sonuç 1,14 eV, gürültü düzeltmeli 4,55 ± 0,19 eV; tam çözüm 3,62 eV. Bugünkü donanımda gürültü baskındır.")}
for n, _ in GRAF + [("fig09_ibm_kingston_donanim.png", None)]:
    a, b, c = GDESC[n]; G += [f"### `{n}`", f"- **(a) Ne gösteriyor:** {a}", f"- **(b) Araç/parametre:** {b}", f"- **(c) Şekil açıklaması:** {c}", ""]
for i, p in enumerate(PAN, 10):
    nm = os.path.basename(p); n = f"fig{i:02d}_panel_{nm}"; mol = nm.replace("_panel_2_3.png", "").replace("ILAC_", "").replace("_panel.png", "")
    isilac = nm.startswith("ILAC_")
    G += [f"### `{n}`",
          f"- **(a) Ne gösteriyor:** {mol} için enerji, soğurma şiddeti (f) ve triyaj paneli{' (klinik ilaç paneli; eski sürüm — PpIX panelindeki ~1 eV klasik hata ifadesi sonradan geri çekildi, sunumda kullanmayın)' if isilac else ''}.",
          f"- **(b) Araç/parametre:** matplotlib (`{'ilac_panel.py' if isilac else 'paneller.py'}`), dataviz paleti.",
          f"- **(c) Şekil açıklaması:** {mol}: CAS(4,4)-VQD+NEVPT2 sonuçlarının görsel özeti (S1, f, triyaj).", ""]
G += ["## B. ChimeraX 3B moleküler görseller (`visuals/chimerax/`)", ""]
TA = {"molekul": "uzay dolgu (van der Waals küreleri)", "molekul_cubuk": "top-çubuk", "HOMO": "HOMO izoyüzeyi (±0,02) + top-çubuk iskelet",
      "LUMO": "LUMO izoyüzeyi (±0,02) + top-çubuk iskelet", "fark": "S1−S0 elektron yoğunluk farkı izoyüzeyi + top-çubuk"}
for n, key, t, core in inv:
    if core is None:
        nm, s1, f_, tri, ex, w_ = MOL[key]; lam = 1239.84 / s1
        exs = f"deney {ex} eV" if ex else "deney değeri yok"
        a = {"molekul": f"{nm} molekülünün GFN2-xTB ile optimize edilmiş 3B yapısı ({ATOMS[key]}); atomlar van der Waals yarıçaplarıyla gösterilmiş.",
             "molekul_cubuk": f"{nm} molekülünün bağ yapısı (top-çubuk); düzlemsel π-halkası ve iç N/metal merkezi görünür.",
             "HOMO": f"{nm} HOMO'su: ışık soğurulmadan önce elektronun bulunduğu en yüksek dolu orbital. Kuantum bölgesinin (4 orbital = 8 kübit) parçasıdır. Mavi ve turuncu dalga fonksiyonunun zıt fazlarıdır.",
             "LUMO": f"{nm} LUMO'su: ışıkla elektronun sıçradığı en düşük boş orbital; kuantum bölgesinin parçası.",
             "fark": f"{nm} için S0→S1 geçişinde elektron yoğunluğunun nereden nereye kaydığı (kuantum bölgesi CASCI yoğunluk matrislerinden). Turuncu = yoğunluğun arttığı, mavi = azaldığı bölgeler."}[t]
        c = {"molekul": f"{nm} ({ATOMS[key]}) — GFN2-xTB geometrisi, uzay dolgu modeli (UCSF ChimeraX).",
             "molekul_cubuk": f"{nm} — top-çubuk modeli (GFN2-xTB geometrisi, UCSF ChimeraX).",
             "HOMO": f"{nm} HOMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). S1 = {s1:.2f} eV ({lam:.0f} nm; {exs}), f = {f_}, triyaj: {tri}.",
             "LUMO": f"{nm} LUMO'su (B3LYP/def2-SVP, izoyüzey ±0,02 a.u.). Kuantum bölgesi: CAS(4e,4o) = 8 kübit; konfigürasyon ağırlıkları {w_[0]} / {w_[1]}.",
             "fark": f"{nm} S0→S1 elektron yoğunluk farkı (turuncu: artış, mavi: azalma). S1 = {s1:.2f} eV, f = {f_}."}[t]
    else:
        nm, at, info, rmsd = DRUGINFO[key]
        a = {"molekul": f"Klinik ilaç {nm} ({at}); RDKit ile 3B'ye gömülüp {core} çekirdeğimize Kabsch ile hizalandı (RMSD {rmsd} Å). Uzay dolgu.",
             "HOMO": f"{nm} üzerinde, {core} çekirdeğinden hesaplanan HOMO izoyüzeyi. Renkli bulutlar kuantum bölgesini (8 kübit) gerçek ilacın fotoaktif halkası üzerinde gösterir.",
             "LUMO": f"{nm} üzerinde {core} çekirdeğinin LUMO'su: ışıkla elektronun gittiği orbital."}[t]
        c = {"molekul": f"{nm} — {info}. Yapı RDKit + MMFF, {core} çekirdeğine hizalanmış (RMSD {rmsd} Å).",
             "HOMO": f"{nm} ve hesapladığımız kuantum bölgesi (HOMO, {core} çekirdeği, B3LYP/def2-SVP, ±0,02 a.u.). Yan gruplar hesaba dahil değildir.",
             "LUMO": f"{nm} — {core} çekirdeğinin LUMO'su (±0,02 a.u.)."}[t]
    ttype = TA.get(t, "")
    b = f"UCSF ChimeraX; {ttype}; ortak parametreler yukarıdaki tabloda; betik `visuals/chimerax_betikleri/` altında ({'gorsel_hq.py' if core is None and t in ('HOMO', 'LUMO', 'fark') else ('molekul_hepsi.cxc' if core is None else 'ilaclar_hepsi.cxc / gercek_ilac.py')} ile üretildi)."
    G += [f"### `{n}`", f"- **(a) Ne gösteriyor:** {a}", f"- **(b) Araç/parametre:** {b}", f"- **(c) Şekil açıklaması:** {c}", ""]
G += ["## C. Sohbette kullanılan ama arşive dahil edilmeyen dış görseller", "",
      "- Ekran görüntüleri: terminal çıktıları, IBM Quantum Platform ekranları, Avogadro ekranları. Bilgi amaçlıdır; sonuçları `data/` klasöründe sayısal olarak yer alır.",
      "- RCSB 4V6X (ribozom) örnek görseli: yalnızca görsel kalite tartışması için kullanıldı.",
      "- Elektrostatik yüzeyli (Gasteiger + coulombic) Temoporfin denemesi: soluk ve beyaz çıktığı için kullanılmadı. Yeniden üretmek için: `addcharge #1 method gasteiger ; surface #1 ; coulombic #1 range -2,2 ; transparency #1 0 target s`.", ""]
open(os.path.join(ROOT, "visuals", "Gorsel_Envanteri_ve_Rehberi.md"), "w", encoding="utf-8").write("\n".join(G))

# ------------------------------------------------------------------ src
S = {"01_bcd_md_analiz": [f"{OUT}/konak_misafir.py", f"{OUT}/temsilci_kare.py", f"{SP}/plot_pmf.py", f"{SP}/analyze.py", f"{SP}/kare_sec.py",
                          f"{SP}/flip_kontrol.py", f"{SP}/hizala.py", f"{SP}/chir.py", f"{SP}/chir2.py", f"{SP}/place2.py", f"{SP}/rim.py", f"{SP}/shape.py",
                          f"{SP}/ik.py", f"{SP}/fig.py", f"{SP}/pmftest/dG_hesapla.py", f"{SP}/pmftest/pencere_sec.py"],
     "02_bcd_kuantum": [f"{SP}/vqe_misafir.py", f"{SP}/xtb_etkilesim.py", f"{SP}/cdtest/fragman_xtb.py", f"{SP}/cdtest/bolge_qm.py", f"{SP}/cdtest/bolge_vqe.py"],
     "03_porfin_kuantum": [f"{SP}/pdt/porfin_geo.py", f"{SP}/pdt/porfin_qm.py", f"{SP}/pdt/porfin_t1.py", f"{SP}/pdt/molden_yaz.py", f"{SP}/pdt/ciz.py",
                           f"{SP}/pdt/t.py", f"{SP}/pdt/t_dip.py", f"{SP}/pdt/t_td.py"],
     "04_tarama_v3": ["/home/claude/v3x/geo.py", "/home/claude/v3x/tarama_qm.py", "/home/claude/v3x/ozet.py", "/home/claude/v3x/calistir.sh",
                      "/home/claude/v3x/kurulum.sh", "/home/claude/v3x/BENIOKU.txt"],
     "05_gorsellestirme": [f"{OUT}/gorsel.py", f"{OUT}/gorsel_hq.py", f"{OUT}/gercek_ilac.py", f"{SP}/ilac_panel.py", f"{SP}/panel.py", f"{SP}/paneller.py",
                           f"{SP}/merdiven.py", f"{SP}/tarama/gi_test.py"],
     "06_referans_v4": ["/home/claude/v4/tddft_referans.py", "/home/claude/v4/cas_buyuk.py", "/home/claude/v4/final_figurler.py",
                        "/home/claude/v4/calistir_v4.sh", "/home/claude/v4/BENIOKU_v4.txt"],
     "07_ibm_v5": ["/home/claude/v5/ibm_vqe.py", "/home/claude/v5/hesap_kaydet.py", "/home/claude/v5/BENIOKU_v5.txt"],
     "08_arsiv_araclari": [f"{SRC}/build_archive.py"]}
missing = []
for d, files in S.items():
    for f_ in files:
        if os.path.exists(f_): cp(f_, "src", d, os.path.basename(f_))
        else: missing.append(f_)
for v in ["kuantum_triyaj_paketi.tar.gz", "kuantum_triyaj_paketi_v2.tar.gz"]:
    dd = D("src", "eski_surumler", v.replace(".tar.gz", "").replace("kuantum_triyaj_paketi", "tarama") or "tarama_v1")
    with tarfile.open(f"{OUT}/{v}") as t: t.extractall(dd)
os.rename(os.path.join(ROOT, "src", "eski_surumler", "tarama"), os.path.join(ROOT, "src", "eski_surumler", "tarama_v1")) if os.path.exists(os.path.join(ROOT, "src", "eski_surumler", "tarama")) else None
open(os.path.join(ROOT, "src", "KAYIP_DOSYALAR.txt"), "w").write(
    "Arsiv olusturulurken bulunamayan betikler (sohbetin ilk bolumunde kullanicinin bilgisayarinda yazilip sandbox'a hic gelmemis olabilir):\n" + "\n".join(missing) +
    "\n\nAyrica kullanicinin bilgisayarinda olup arsive alinamayanlar: misafir_qm.py (beta-CD misafir RHF/MP2), GROMACS .mdp dosyalari ve tum yorunge (.xtc/.trr) dosyalari, ~/tarama/*_tddft.json (porfin ve klorin disindakiler), *_vqe_parametreleri.json, ibm_vqe.py'nin kullanici bilgisayarinda olusturdugu porfin_ibm_hazir.pkl / porfin_ibm_is.json.\n")

# ------------------------------------------------------------------ dizin agaci
def tree(p, pre=""):
    out = []; items = sorted(os.listdir(p), key=lambda s: (not os.path.isdir(os.path.join(p, s)), s.lower()))
    for i, it in enumerate(items):
        last = i == len(items) - 1; fp = os.path.join(p, it)
        out.append(pre + ("└── " if last else "├── ") + it + ("/" if os.path.isdir(fp) else ""))
        if os.path.isdir(fp): out += tree(fp, pre + ("    " if last else "│   "))
    return out
T = ["Proje_Arsivi/"] + tree(ROOT)
open(os.path.join(ROOT, "docs", "Dizin_Agaci.txt"), "w", encoding="utf-8").write("\n".join(T) + "\n")
print("\n".join(T)); print("EKSIK:", missing)
