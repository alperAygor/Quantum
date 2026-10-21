"""Ortak grafik stili (tutarlı renkler, ince işaretler, sade eksenler)."""
import matplotlib as mpl
import matplotlib.pyplot as plt

# Kategorik renkler sabit sırada (kompleks hep mavi, çözücü hep turuncu ...)
C = dict(
    blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", yellow="#eda100",
    magenta="#e87ba4", green="#008300", violet="#6250d6", red="#e34948",
    ink="#0b0b0b", ink2="#52514e", muted="#8a8984", grid="#e6e5e0", surface="#fcfcfb",
)
LEG = {"complex": C["blue"], "solvent": C["orange"]}
LEG_TR = {"complex": "kompleks (bCD + misafir)", "solvent": "çözücü (sadece misafir)"}
SERIES = [C["blue"], C["orange"], C["aqua"], C["yellow"], C["magenta"], C["green"], C["violet"], C["red"]]


def setup():
    mpl.rcParams.update({
        "figure.facecolor": C["surface"], "axes.facecolor": C["surface"],
        "axes.edgecolor": C["muted"], "axes.labelcolor": C["ink"], "text.color": C["ink"],
        "xtick.color": C["ink2"], "ytick.color": C["ink2"],
        "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 2.0, "lines.markersize": 6,
        "axes.prop_cycle": mpl.cycler(color=SERIES),
        "figure.dpi": 110, "font.size": 10.5, "legend.frameon": False,
        "axes.titleweight": "bold", "axes.titlesize": 11.5,
    })


setup()
