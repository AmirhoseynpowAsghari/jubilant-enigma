"""
build_pdf.py -- assembles the corrected version of arXiv:2607.02812.

* original 18 pages verbatim (text untouched);
* the seven wrong figures (pp. 11-15) replaced by the recomputed ones;
* errata notes inserted in the free band at the top of pp. 11-15;
* Errata & Addendum appended as pages 19-24.
"""
import shutil
import os
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
VEND = os.path.join(HERE, "vendor")
SRC = os.path.join(VEND, "arXiv_2607.02812_original.pdf")  # original paper, 18 pp (vendored; was misnamed src.tar.gz, F-16)
OUT = os.path.join(HERE, "results", "arXiv_2607.02812_corrected.pdf")
FIG = os.path.join(VEND, "figs") + os.sep          # recomputed replacement figures (vendored)
ETA1 = os.path.join(VEND, "eta1") + os.sep         # master-set addendum figures (vendored)
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
for _need in (SRC, os.path.join(VEND, "figs"), os.path.join(VEND, "eta1")):
    if not os.path.exists(_need):
        raise SystemExit(f"missing vendored input: {_need}")
from matplotlib import font_manager as _fm
# Prefer an installed system font: PyMuPDF rejects some otherwise-valid copies bundled
# inside matplotlib wheels. Resolve dynamically; never hard-code a host-specific path.
_sys_ttf = sorted(_fm.findSystemFonts(fontext="ttf"))
TTF = next((p for p in _sys_ttf if os.path.basename(p).lower() == "dejavusans.ttf"),
           _fm.findfont("DejaVu Sans", fallback_to_default=False))

doc = pymupdf.open(SRC)

# ---------------------------------------------------------------- replacements
# (page index, old xref, new file, rect)
REPL = [
    (10, 241, "fig_group_A.png",  (26.2, 228.8, 582.7, 748.5)),
    (11, 244, "fig_group_B.png",  (15.8, 169.5, 570.0, 756.0)),
    (12, 247, "fig_group_C.png",  (27.8, 169.5, 584.2, 615.0)),
    (13, 251, "fig_group_D.png",  (24.0, 169.5, 598.5, 564.0)),
    (13, 250, "fig_scan_U95.png", (84.8, 639.3, 523.3, 788.0)),
    (14, 256, "fig_scan_U125.png", (27.8, 211.8, 557.2, 474.3)),
    (14, 257, "fig_scan_Um75.png", (72.0, 548.9, 540.0, 788.0)),
]
for pno, xref, fn, rect in REPL:
    page = doc[pno]
    page.delete_image(xref)
    page.insert_image(pymupdf.Rect(*rect), filename=FIG + fn, keep_proportion=True)

# ---------------------------------------------------------------- errata notes
NOTE_STD = ("ERRATA NOTE (E1-E5; Addendum pp. 19-24).  Figure replaced by a recomputation:  "
            "W/W_N = S_sc/S_N,  S_sc = \u03a3 over helicity branches of "
            "\u222b dE \u03bd\u03c3(E) \u03bd\u03c3\u2032(E+\u0127\u03c9) \u00b7 \u00bd[1 + (\u03be\u03be\u2032 \u2212 \u0394\u0304\u0394\u0304\u2032)/(EE\u2032)] \u00b7 [f(E) \u2212 f(E+\u0127\u03c9)],  "
            "\u0127\u03c9 = 0.7t, \u03b7 = +1, scattering channel only; stage-1 inputs (Tc, \u0394_eff, \u03bc) unchanged.  "
            "The corrected ratio is < 1 for all T < Tc and saturates below 1 at Tc\u207b.")
NOTE_SHORT = ("ERRATA NOTE (E1-E5; Addendum pp. 19-24).  Figure replaced by a recomputation: "
              "W/W_N = S_sc/S_N with S_sc = \u03a3_{\u03c3,\u03c3\u2032} \u222b dE \u03bd\u03c3(E)\u03bd\u03c3\u2032(E+\u0127\u03c9) "
              "\u00b7 \u00bd[1+(\u03be\u03be\u2032\u2212\u0394\u0304\u0394\u0304\u2032)/(EE\u2032)] \u00b7 [f(E)\u2212f(E+\u0127\u03c9)], \u0127\u03c9=0.7t, \u03b7=+1, scattering channel; "
              "stage-1 inputs unchanged.  Corrected ratio < 1 for all T < Tc, saturating below 1 at Tc\u207b.")
NOTES = {
    10: (NOTE_SHORT, (72.0, 183.0, 540.0, 226.0)),
    11: (NOTE_STD + "   (E6: the Group B title should read \u0394_t = 4.5.)", (72.0, 84.0, 540.0, 165.0)),
    12: (NOTE_STD, (72.0, 108.0, 540.0, 165.0)),
    13: (NOTE_STD + "   Upper figure: Group D regenerated consistently from (U=95t, \u0394_t=4.5t, n=1.885); "
                    "the original auxiliary panels quoted Tc values belonging to Group A, U=\u22124t (E6).  "
                    "Lower figure: recomputed V_SO scan at fixed T/Tc.", (72.0, 66.0, 540.0, 166.0)),
    14: (NOTE_STD + "   Relaxed policy: residual-verified U=\u221275t roots may lie far outside the band; flagged mu_window.inside=false and not interpreted as Fermi-surface states (E7/E12).",
         (72.0, 84.0, 540.0, 168.0)),
}
for pno, (txt, rect) in NOTES.items():
    page = doc[pno]
    rc = page.insert_textbox(pymupdf.Rect(*rect), txt, fontsize=7.6, fontname="djv",
                             fontfile=TTF, color=(0.55, 0.0, 0.0), lineheight=1.25)
    assert rc >= 0, (pno, rc)

# ---------------------------------------------------------------- addendum (regenerated text)
import addendum_text_en as AT
T_ = AT.build()
W, H = 612.0, 792.0
M = (62.0, 62.0, W - 62.0, H - 56.0)


def new_page():
    return doc.new_page(width=W, height=H)


def put(page, txt, rect=M, size=9.3, color=(0, 0, 0)):
    rc = page.insert_textbox(pymupdf.Rect(*rect), txt, fontsize=size, fontname="djv",
                             fontfile=TTF, color=color, lineheight=1.32)
    assert rc >= 0, rc
    return rc


def flow(blocks, size):
    """Place text blocks on as many new pages as needed (no block is split across pages)."""
    page = new_page()
    cur = page
    for blk in blocks:
        rc = cur.insert_textbox(pymupdf.Rect(*M), blk, fontsize=size, fontname="djv",
                                fontfile=TTF, lineheight=1.32)
        if rc < 0:                       # does not fit: continue on a fresh page
            cur = new_page()
            rc = cur.insert_textbox(pymupdf.Rect(*M), blk, fontsize=size, fontname="djv",
                                    fontfile=TTF, lineheight=1.32)
            assert rc >= 0, "a single block exceeds one page"
    return cur


flow([T_["title"], T_["summary"]], 8.2)
flow([T_["table1"] + "\n" + T_["table2"] + "\n" + T_["table3"], T_["statements"],
      T_["validation"], T_["reproduction"], T_["reading"]], 8.4)

for fn, cap in (("fig_WWN_vs_T.png",
                 "Figure A1 (master set U=95t, Dt=4.5t, n=1.885; regenerated).  W/W_N vs T/Tc for V_SO = 0, 0.4082, 0.8, "
                 "scattering channel, eta=+1: variant (a) (signed xi xi', primary) against (b) |xi||xi'|, (c) the paper's "
                 "combination and (d) the form without xi xi' (not a Bogoliubov matrix element; E13)."),
                ("fig_branch_decomposition.png",
                 "Figure A2 (regenerated).  Branch decomposition (sigma,sigma') of the corrected scattering sum for the "
                 "master set at V_SO = 0.4082, variant (a)."),
                ("fig_limits.png",
                 "Figure A3 (regenerated).  Limits lim_{T->Tc-} W/W_N of the four variants versus V_SO (master set).  "
                 "The omega->0 plateau 2 f(Delta) is verified separately (T2)."),
):
    p = new_page()
    pix = pymupdf.Pixmap(ETA1 + fn)
    ar = pix.width / pix.height
    w = 512.0
    h = w / ar
    p.insert_image(pymupdf.Rect(50, 80, 50 + w, 80 + h), filename=ETA1 + fn)
    put(p, cap, rect=(62, 90 + h, 550, 740), size=9.0)

doc.save(OUT, garbage=4, deflate=True)
import shutil
shutil.copyfile(OUT, os.path.join(HERE, "arXiv_2607.02812_corrected.pdf"))
print("saved", OUT, "pages:", doc.page_count)
