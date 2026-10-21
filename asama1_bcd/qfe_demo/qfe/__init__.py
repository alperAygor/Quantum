"""qfe — FreeQuantum pipeline'ının (Günther et al., JCTC 2026) laptop ölçekli eğitim amaçlı simülasyonu."""
import warnings, logging
warnings.filterwarnings("ignore", message=".*JAX.*")
for _n in ("pymbar", "openmmtools", "openmmtools.alchemy"):
    logging.getLogger(_n).setLevel(logging.ERROR)
warnings.filterwarnings("ignore", module="sklearn")
from . import config  # noqa
