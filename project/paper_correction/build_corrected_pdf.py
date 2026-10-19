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
SRC = os.path.join(VEND, "src.tar.gz")            # original LaTeX-built source (vendored input)
OUT = os.path.join(HERE, "results", "arXiv_2607.02812_corrected.pdf")
FIG = os.path.join(VEND, "figs") + os.sep          # recomputed replacement figures (vendored)
ETA1 = os.path.join(VEND, "eta1") + os.sep         # master-set addendum figures (vendored)
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
for _need in (SRC, os.path.join(VEND, "figs"), os.path.join(VEND, "eta1")):
    if not os.path.exists(_need):
        raise SystemExit(f"missing vendored input: {_need}")
TTF = "/usr/local/lib/python3.13/site-packages/matplotlib/mpl-data/fonts/ttf/DejaVuSans.ttf"

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
    14: (NOTE_STD + "   For U=\u221275t the stage-1 maximum-gap root has \u03bc \u2248 35t, outside the band "
                    "(Tc \u2248 12.7t): the set lies outside the model's physical window (E7).",
         (72.0, 84.0, 540.0, 168.0)),
}
for pno, (txt, rect) in NOTES.items():
    page = doc[pno]
    rc = page.insert_textbox(pymupdf.Rect(*rect), txt, fontsize=7.6, fontname="djv",
                             fontfile=TTF, color=(0.55, 0.0, 0.0), lineheight=1.25)
    assert rc >= 0, (pno, rc)

# ---------------------------------------------------------------- addendum
W, H = 612.0, 792.0
M = (62.0, 62.0, W - 62.0, H - 56.0)


def new_page():
    return doc.new_page(width=W, height=H)


def put(page, txt, rect=M, size=9.3, color=(0, 0, 0)):
    rc = page.insert_textbox(pymupdf.Rect(*rect), txt, fontsize=size, fontname="djv",
                             fontfile=TTF, color=color, lineheight=1.32)
    assert rc >= 0, rc
    return rc


TITLE = ("ERRATA AND ADDENDUM  (pages 19-24)\n"
         "to arXiv:2607.02812 - \"Theoretical investigation of power absorbed by superconductors "
         "with Rashba spin-orbit coupling in the presence external time-dependent perturbation\"  "
         "(R. Afzali, B. Takrim, A. H. Asghari)\n"
         "Prepared 8 October 2026.  The original text of the paper is reproduced verbatim; every "
         "correction in this file is additive: (i) the seven figures on pages 11-15 are replaced by "
         "recomputations (each tagged [corrected]); (ii) errata notes are inserted on pages 11-15; "
         "(iii) this addendum.  High-resolution corrected figures and all reproduction scripts "
         "accompany this file (paper_fix/figs/, paper_fix/*.py, wwn_eta1/).\n")

p = new_page()
put(p, TITLE + "\nSummary of corrections\n"
    "E1  Coherence factor and energy conservation.  The working ratio (equation on p. 11) combines "
    "the coherence factor with the roots \u221a(E\u00b2\u2212\u0394\u00b2) (i.e. |\u03be|) and with the term "
    "+\u03b7 \u0394\u0304\u0394\u0304\u2032/EE\u2032.  The correct golden-rule combination for \u03b7 = +1 (Case I, ultrasound), keeping the "
    "sign of \u03be and imposing energy conservation E\u2032 = E + \u0127\u03c9 on the quasiparticle energies, is\n"
    "        C_I(E,\u03be; E\u2032,\u03be\u2032) = \u00bd [ 1 + (\u03be\u03be\u2032 \u2212 \u0394\u0304\u0394\u0304\u2032) / (E E\u2032) ] ,\n"
    "with the sum over the two Rashba helicity branches made explicit.  Replacing \u03be\u03be\u2032 by |\u03be||\u03be\u2032| "
    "(variants (b), (c) of the Addendum) drives the ratio above 1 in the \u0394\u21920 limit (up to 1.97), "
    "which is unphysical for a scattering-only channel.\n"
    "E2  All tabulated and plotted W/W_N values (flat, 0.90-1.00) are wrong.  Corrected: W/W_N < 1 for "
    "every 0 < T < Tc and every parameter set, strongly suppressed at low T and rising monotonically "
    "towards Tc; e.g. master set U=95t, \u0394_t=4.5t, n=1.885: lim W/W_N = 0.4024 / 0.5086 / 0.7616 at "
    "V_SO = 0 / 0.4082 / 0.8 (Table 2).\n"
    "E3  Limit at Tc.  W/W_N does not tend to 1 at Tc\u207b in this model: it saturates at R(\u0127\u03c9,Tc), the "
    "band-structure value of the \u0394=0 ratio, because the helicity coherence factor is a normal-state "
    "(band) property, not a pairing property.  Adding the pair-breaking channel (E+E\u2032=\u0127\u03c9, occupation "
    "1\u2212f(E)\u2212f(E\u2032), coherence factor 1\u2212C_I) is examined as an OUTLOOK only: the scattering channel alone returns R(hbar*omega,T_c) at T_c-; restoration of W/W_N -> 1 would involve the pair-breaking channel, outside the quantitative result (REPORT App. A2).\n"
    "E4  Pair-breaking threshold.  The validity cut used in the paper, 2\u0394_eff(T) > \u0127\u03c9, is not the "
    "pair-breaking threshold.  The true kinematic threshold is 2 E_min(T) = \u0127\u03c9 with "
    "E_min = |\u0394\u00b0 \u2212 b\u03bc|/\u221a(1+b\u00b2), b = \u0394\u02e2/4t, the minimum of the quasiparticle energy.  Inside the domain "
    "2\u0394_eff > \u0127\u03c9 the pair-breaking channel is identically zero, so the scattering-only result is complete "
    "there and the cut is stricter than needed.\n"
    "E5  Normal-state reference and DOS.  W and W_N must use the same DOS N(\u03b5), normalised to "
    "\u222bN d\u03b5 = 2 and with exact band support |\u03b5| \u2264 4t\u221a(1+V_SO\u00b2/2t\u00b2).  The original stage-2 DOS leaked "
    "Gaussian weight above the band edge and was mis-normalised (\u222bN d\u03b5 = 0.643); negative W and W_N "
    "were masked by abs().\n"
    "E6  Parameter misstatements.  Group B's title reads \u0394_t = 4; the printed Tc and \u0394_eff reproduce only "
    "for \u0394_t = 4.5.  Group D's auxiliary panels (Tc = 0.718/0.688/0.590, \u0394\u2080 = 1.327/1.262/1.072) do not "
    "belong to (U=95t, n=1.885): they coincide with Group A's U=\u22124t row; Group D is regenerated "
    "consistently from (U=95t, \u0394_t=4.5t, n=1.885).\n"
    "E7  U=\u221275t scan.  For U=\u221275t the stage-1 maximum-gap root (the paper solver's own selection rule) "
    "is a strong-coupling solution with \u03bc \u2248 35t, far outside the band (Tc \u2248 12.7t, \u0394_eff(0) \u2248 17.8t).  The "
    "corrected attenuation of that root is plotted as stated, but the set lies outside the physical "
    "window of the tight-binding model and should not carry physical claims.\n"
    "E8  Gap shape.  Stage 2 uses the rigid master gap shape phi(T) = (1-(T/Tc)^3.434)^0.553 for all "
    "V_SO; a per-V_SO refit (3.674/0.556 at V_SO = 0) changes the reported W/W_N by 1.3-3.4 % "
    "(measured, results/gap_refit_check.json).\n"
    "E9  Citability.  U=-75t outputs are excluded from the citable presentation figures "
    "(SOCfigures/non_citable/) and retained only as a historical record.",
    size=8.6)

p = new_page()
put(p, "1  Corrected formalism (derivation before numbers)\n"
    "Model.  Bands \u03b5_\u00b1(k) = \u22122t(cos kx + cos ky) \u2213 2 V_SO |sin k|;  HHM gap \u0394\u0304(\u03b5) = \u0394\u00b0 \u2212 (\u0394\u02e2/4t)\u03b5;  "
    "quasiparticle energy E_\u03c3(\u03b5) = \u221a((\u03b5\u2212\u03bc)\u00b2 + \u0394\u0304(\u03b5)\u00b2).  Branch inversion (exact): with b = \u0394\u02e2/4t, "
    "\u03ba = \u221a(1+b\u00b2), \u03b5* = (\u03bc + b\u0394\u00b0)/(1+b\u00b2), \u2113 = \u221a(E\u00b2\u2212E_min\u00b2), E_min = |\u0394\u00b0\u2212b\u03bc|/\u03ba, the band energies at "
    "quasiparticle energy E on branch \u03c3 are \u03b5_\u03c3 = \u03b5* + \u03c3\u2113/\u03ba, and the integration measure becomes "
    "\u03bd_\u03c3(E) dE = N(\u03b5_\u03c3(E)) d\u2113/\u03ba, which removes the BCS 1/\u221a(E\u00b2\u2212E_min\u00b2) singularity exactly.\n"
    "Ultrasound (\u03b7 = +1, Case I).  The T-matrix element carries u u\u2032 \u2212 \u03b7 v v\u2032; squaring and averaging "
    "gives the coherence factor C_I = \u00bd[1 + (\u03be\u03be\u2032 \u2212 \u0394\u0304\u0394\u0304\u2032)/(EE\u2032)] with \u03be = \u03b5\u2212\u03bc signed.\n"
    "Scattering channel (the in-window answer):\n"
    "   W = \u03a3_{\u03c3,\u03c3\u2032} \u222b dE \u03bd_\u03c3(E) \u03bd_\u03c3\u2032(E+\u0127\u03c9) C_I [f(E) \u2212 f(E+\u0127\u03c9)],\n"
    "with per-pair upper limits from the band support: E \u2264 min(E_max(\u03c3), E_max(\u03c3\u2032)\u2212\u0127\u03c9), "
    "E_max(s) = \u221a(\u03ba\u00b2(W_band \u2212 s\u03b5*)\u00b2 + E_min\u00b2).  At \u0127\u03c9 = 0.7t the (+1,+1) and (\u22121,+1) pairs vanish "
    "identically (E_max(+1) < \u0127\u03c9).\n"
    "Normal reference:  W_N = \u222b dE N(E) N(E+\u0127\u03c9) [f_N(E) \u2212 f_N(E+\u0127\u03c9)]  with the same N(\u03b5) and \u03bc_N from "
    "the normal-state number equation.\n"
    "Pair-breaking channel (reported separately; zero inside the window 2\u0394_eff > \u0127\u03c9):\n"
    "   W_p = \u03a3_{\u03c3,\u03c3\u2032} \u222b dE \u03bd_\u03c3(E) \u03bd_\u03c3\u2032(\u0127\u03c9\u2212E) (1\u2212C_I) [1 \u2212 f(E) \u2212 f(\u0127\u03c9\u2212E)],  nonzero only for \u0127\u03c9 > 2E_min.\n"
    "Elementary checks: 0 \u2264 C_I \u2264 1 pointwise; at \u0394=0, C_I \u2192 branch projector (1 if \u03c3=\u03c3\u2032 else 0); "
    "M_scat\u00b2 + M_pair\u00b2 = 1 to 1e-16; at \u0394=0 the two-channel sum equals S_N to 1e-6 (a \u0394\u21920 limit where the sheets coincide); NO finite-\u0394 identity is claimed (REPORT App. A2).\n\n"
    "2  Corrected results\n"
    "Table 1 - groups A-D (scattering channel, \u03b7=+1, signed \u03be\u03be\u2032, \u0127\u03c9=0.7t):\n"
    "   set                 Tc/t     \u0394_eff(0)/t   min W/W_N   lim T\u2192Tc\u207b\n"
    "   A  U=\u2212t             0.0123    0.0214      0.000       0.0198\n"
    "   A  U=\u22122t            0.1818    0.3206      0.000       0.2949\n"
    "   A  U=\u22123t            0.4388    0.7969      0.000       0.5567\n"
    "   A  U=\u22124t            0.7184    1.3273      0.000       0.6731\n"
    "   B  U=75t            0.5009    0.6272      0.000       0.6717\n"
    "   B  U=115t           0.1399    0.2840      0.000       0.1837\n"
    "   B  U=125t           0.0861    0.2070      0.000       0.1063\n"
    "   C  \u0394_t=3.5          0.4094    0.6057      0.000       0.4652\n"
    "   C  \u0394_t=4.0          0.6720    0.8414      0.000       0.7930\n"
    "   C  \u0394_t=4.5          0.9536    1.0934      0.000       0.9890\n"
    "   D  V_SO=0           0.2861    0.4444      0.000       0.3999\n"
    "   D  V_SO=0.4         0.3770    0.5447      0.000       0.5019\n"
    "   D  V_SO=0.8         0.6494    0.8542      0.000       0.7602\n"
    "Table 2 - master set (U=95t, \u0394_t=4.5t, n=1.885), V_SO scan of the limit:\n"
    "   V_SO    Tc/t     \u0394_eff(0)/t   lim W/W_N   R(0.7,Tc) flat-DOS\n"
    "   0.000   0.2861   0.4444       0.40241     0.49867\n"
    "   0.4082  0.3808   0.5489       0.50856     0.59355\n"
    "   0.800   0.6494   0.8542       0.76163     0.74264\n"
    "Table 3 - V_SO scans at fixed T/Tc (range of W/W_N over V_SO=0\u21920.8 at T/Tc=0.9):\n"
    "   U=95t:   0.253 \u2192 0.435      U=125t:  0.069 \u2192 0.253      U=\u221275t:  0.489 \u2192 0.496  (E7)\n"
    "Statements.  (i) W/W_N < 1 strictly, for every set and every 0 < T < Tc.  (ii) The Tc\u207b limit "
    "equals R(\u0127\u03c9,Tc) < 1 (E3), not 1; the original flat 0.90-1.00 curves and their jump to exactly 1 at "
    "Tc are artefacts of E1/E5.  (iii) Inside the paper's window 2\u0394_eff > \u0127\u03c9 the pair-breaking channel is "
    "identically zero (E4); for the master set the true threshold 2E_min=\u0127\u03c9 sits at T/Tc = 0.820 / 0.912 / "
    "0.975 (V_SO = 0 / 0.4082 / 0.8) whereas the paper's cut 2\u0394_eff=\u0127\u03c9 sits at 0.740 / 0.840 / 0.936.  "
    "(iv) Including the pair-breaking channel above threshold, W/W_N rises to exactly 1 at Tc\u207b "
    "(overshoot 1.09 / 1.34 / 1.99 near T/Tc \u2248 0.96-0.98 for the master set).",
    size=8.6)

p = new_page()
put(p, "3  Validation of the recomputation\n"
    "T0  Quadrature: the \u2113-substitution equals the E-integral with BOTH Jacobians (E/\u03ba\u2113)(E\u2032/\u03ba\u2113\u2032); relative difference \u22643.9e-4 (3.9/1.8/0.8e-4 at \u0127\u03c9=0.05/0.01/0.002). The committed test originally carried a single Jacobian (0.30-0.35 residual); fixed.\n"
    "T1  \u0394\u21920 closed form: C_I collapses to the branch projector and W/W_N \u2192 R(\u0127\u03c9,T) = "
    "(2/\u0127\u03c9)\u222b\u2080^\u221e [f(E)\u2212f(E+\u0127\u03c9)] dE; verified to 3.1e-6.  Sum rule at \u0394=0: (S_sc+S_pair)/S_N = 100.0000%; above threshold the band-reference total exceeds 1 (App. A2).\n"
    "T2  \u03c9\u21920 limit: the exact golden-rule ratio plateaus at 2f(\u0394) = Levy 1963 Eq. (8); verified to 2e-6 "
    "(T5 reproduces Levy's s-wave ratio to 3e-6).\n"
    "T3  Independent brute-force lattice sum over (k,k\u2032) pairs with \u03c3 = sign(\u03b5\u2212\u03b5*): agrees with the "
    "energy-space pipeline per variant: (a) 1.0-1.3%, (b),(c) \u22640.7%, (d) 1.5-45% (variant (d) is not used "
    "for results); the (d) discrepancy is the gap-edge regularisation of the \u03be\u03be\u2032-free kernel.\n"
    "T4  0 \u2264 C_I \u2264 1 pointwise on a structured grid (gap edge, mixed signs); ratio-level bound at \u0394=0: (a),(d)=R\u22641 while (b),(c)=2R>1 (1.9002 at \u0127\u03c9=0.2, T=1).\n"
    "Defects of the original stage-2 code, for the record: \u222bN d\u03b5 = 0.6430 (mis-normalisation); Gaussian "
    "leak above the band edge (g(4.2) = 3.6e-2); negative W and W_N hidden by abs(); wrong pair-breaking "
    "threshold 2\u0394_eff; T < 1e-4 shortcut returning 0.\n\n"
    "4  Reproduction\n"
    "Stage 1 (Tc, \u0394\u00b0, \u0394\u02e2, \u03bc for every parameter set): paper_fix/stage1.py - verbatim replica of the "
    "repository solver (project/sc_params, branch master); reproduces every Tc and \u0394_eff printed in the "
    "paper to 4 decimals (Groups A-C) and the repository's own scan (Group D, master set).\n"
    "Stage 2 (W/W_N): wwn_eta1/model.py + wwn_eta1/energy_space.py (validated pipeline, tests T0-T4 in "
    "wwn_eta1/validate.py); production runs paper_fix/run_groups.py (Groups A-D) and "
    "paper_fix/run_scans.py (V_SO scans); master-set tables and figures in wwn_eta1/REPORT.md, "
    "wwn_eta1/results.{json,csv}.\n"
    "Figures: paper_fix/make_paper_figs.py (this file's replaced figures); this PDF: paper_fix/build_pdf.py.\n\n"
    "5  Reading guide to the replaced figures\n"
    "Pages 11-13 (Groups A-C) and page 14 top (Group D): six-panel composites in the original layout.  "
    "The two W/W_N panels are the corrected ones; the four auxiliary panels (\u03bc(T), \u0394(T), universal gap "
    "curve) are recomputed from the same stage-1 solution for internal consistency.  Note the qualitative "
    "change: exponential suppression at low T, monotonic rise, saturation below 1 at Tc\u207b.\n"
    "Page 14 bottom and page 15: W/W_N vs V_SO at fixed T/Tc = 0.3, 0.5, 0.7, 0.9.  The corrected curves "
    "increase with both V_SO and T/Tc and stay below 1; the non-monotonic bumps of the originals (and "
    "their 0.93-0.97 range) are gone.\n"
    "Pages 22-24: master-set detail figures - W/W_N vs T/Tc for the three V_SO values with all coherence-"
    "factor variants; branch decomposition of the scattering sum; the \u0394\u21920 and T\u2192Tc limits.",
    size=8.8)

for fn, cap in (("fig_WWN_vs_T.png",
                 "Figure A1 (master set U=95t, \u0394_t=4.5t, n=1.885).  Corrected W/W_N vs T/Tc for "
                 "V_SO = 0, 0.4082, 0.8: scattering channel only, \u03b7=+1, signed \u03be\u03be\u2032 (variant (a)); "
                 "compared with | \u03be || \u03be\u2032 | (b), the paper's combination (c) and the textbook Case-I form "
                 "without \u03be\u03be\u2032 (d).  Only (a) stays below 1 and has a consistent \u0394\u21920 limit."),
                ("fig_branch_decomposition.png",
                 "Figure A2.  Branch decomposition (\u03c3,\u03c3\u2032) of the corrected scattering sum for the master "
                 "set: the (\u22121,\u22121) intra-band channel dominates; (+1,+1) and (\u22121,+1) vanish at "
                 "\u0127\u03c9 = 0.7t because E_max(+1) < \u0127\u03c9."),
                ("fig_limits.png",
                 "Figure A3.  Limits of the corrected ratio: \u0394\u21920 (T\u2192Tc\u207b) saturation at R(\u0127\u03c9,Tc) and the "
                 "\u03c9\u21920 plateau 0.647\u00b72f(\u0394); adding the pair-breaking channel above the true threshold "
                 "2E_min = \u0127\u03c9 the pair channel is an outlook only; the quantitative result is scattering-only (App. A2).")):
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
