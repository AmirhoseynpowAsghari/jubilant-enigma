# -*- coding: utf-8 -*-
"""eqrender.py -- render LaTeX-style math to transparent PNGs (matplotlib mathtext)."""
import hashlib
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pymupdf

plt.rcParams["mathtext.fontset"] = "dejavuserif"
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eqcache")
os.makedirs(CACHE, exist_ok=True)


def eq_img(latex, pt):
    key = hashlib.md5((latex + f"{pt:.1f}").encode()).hexdigest()[:12]
    path = f"{CACHE}/{key}.png"
    if not os.path.exists(path):
        fig = plt.figure(figsize=(0.5, 0.5))
        fig.text(0.5, 0.5, "$" + latex + "$", fontsize=pt, ha="center", va="center")
        fig.savefig(path, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.015)
        plt.close(fig)
    px = pymupdf.Pixmap(path)
    return path, px.width / 300.0 * 72.0, px.height / 300.0 * 72.0
