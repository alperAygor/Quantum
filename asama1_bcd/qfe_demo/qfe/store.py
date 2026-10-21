"""Aşamalar arası veri aktarımı: her aşama sonuçlarını results/ altına yazar,
sonraki aşama oradan okur. (Makaledeki MongoDB veritabanının basit karşılığı.)"""
import json
import time
from contextlib import contextmanager

import numpy as np

from .config import RESULTS, preset_name


def _path(name):
    d = RESULTS / preset_name()
    d.mkdir(parents=True, exist_ok=True)
    return d / name


def save_json(name, obj):
    p = _path(name)
    old = {}
    if p.exists():
        old = json.loads(p.read_text())
    old.update(obj)
    p.write_text(json.dumps(old, indent=2, ensure_ascii=False, default=float))
    return p


def load_json(name, default=None):
    p = _path(name)
    if not p.exists():
        if default is not None:
            return default
        raise FileNotFoundError(f"{p} yok — önceki notebook'u çalıştırdınız mı?")
    return json.loads(p.read_text())


def save_npz(name, **arrays):
    p = _path(name)
    np.savez_compressed(p, **arrays)
    return p


def load_npz(name):
    p = _path(name)
    if not p.exists():
        raise FileNotFoundError(f"{p} yok — önceki notebook'u çalıştırdınız mı?")
    return dict(np.load(p, allow_pickle=True))


def path(name):
    return _path(name)


@contextmanager
def timed(label):
    """Süre ölçer ve results/<preset>/timings.json içine yazar.
    07. notebook bu süreleri 'nerede HPC gerekir' tablosunda kullanır."""
    t0 = time.time()
    rec = {}
    yield rec
    dt = time.time() - t0
    rec["seconds"] = dt
    timings = load_json("timings.json", default={})
    timings[label] = rec
    save_json("timings.json", timings)
    print(f"[süre] {label}: {dt:.1f} s")
