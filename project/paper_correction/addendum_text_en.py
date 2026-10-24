"""
addendum_text_en.py -- the English addendum text of the corrected paper.  Every number is taken
from addendum_numbers.numbers() (regenerated results).  Route A: no pair-channel number appears.
"""
from addendum_numbers import numbers, fmt


def build():
    n = numbers()
    M = n["master"]
    m0, m1, m2 = M[0], M[1], M[2]
    G = {(r["group"], r["label"]): r for r in n["groups"]}
    S = n["scans"]
    t = {}

    t["title"] = (
        "ERRATA AND ADDENDUM  (pages 19-23)\n"
        "to arXiv:2607.02812 - \"Theoretical investigation of power absorbed by superconductors with Rashba "
        "spin-orbit coupling in the presence external time-dependent perturbation\" (R. Afzali, B. Takrim, A. H. Asghari)\n"
        "Corrected 24 October 2026 (audit-driven correction campaign; full record in CORRECTION_LOG.md of the "
        "companion repository).  The original text is reproduced verbatim; corrections are additive: "
        "(i) the replaced figures on pp. 11-15 are regenerated from the corrected code; (ii) errata notes on "
        "pp. 11-15; (iii) this addendum.  Numbers in this addendum are generated from the regenerated results.\n")

    t["summary"] = (
        "Summary of corrections\n"
        "E1  Coherence factor.  The scattering channel for eta = +1 (ultrasound) uses the signed factor "
        "C_I = 1/2 [1 + (xi xi' - D D')/(E E')].  This is now verified by an independent derivation: the "
        "numerically computed Bogoliubov matrix element (uu' - vv')^2 agrees with C_I to 1e-15.  Variants (b) and "
        "(c) (absolute values of xi) violate the bound C <= 1 at Delta = 0 and have no Bogoliubov derivation.\n"
        "E2  Corrected values.  Master set U=95t, Dt=4.5t, n=1.885, hbar*omega=0.7t, variant (a): "
        f"lim_{{T->Tc-}} W/W_N = {fmt(m0['lim'])} / {fmt(m1['lim'])} / {fmt(m2['lim'])} at V_SO = 0 / 0.4082 / 0.8 "
        f"(previously 0.4000 / 0.5054 / 0.7611; Table 2).  All group values change by the amounts in Table 1.\n"
        "E3  Limit at Tc.  In the scattering channel W/W_N does not tend to 1 at Tc-; it saturates at R(hbar*omega, Tc) < 1, "
        "the band-structure value of the Delta = 0 ratio.  The pair-breaking channel is an outlook only and is not "
        "used for any quantitative statement here.\n"
        "E4  Pair-breaking threshold.  The validity cut of the paper, 2 Delta_eff(T) > hbar*omega, is not the kinematic "
        "threshold; the true threshold is 2 E_min(T) = hbar*omega.  Inside the paper's window the scattering-only result "
        "is unaffected; the paper's cut is stricter than the kinematics require.  No threshold temperature is quoted here "
        "(pair channel, outlook only).\n"
        "E5  Normal reference and DOS.  W and W_N use the same DOS N(eps), normalised to 2 and with exact band support.  "
        "The normal-state integral S_N is now evaluated with Gauss-Legendre quadrature on a mesh that contains every DOS "
        "kink, with a runtime convergence assertion (relative spread 2e-16 across orders).  The former midpoint sum "
        "differed by +0.06 %.\n"
        "E6  Parameter misstatements (unchanged from the first erratum).  Group B's title reads Dt = 4; the printed "
        "values reproduce only for Dt = 4.5.  Group D's auxiliary panels belong to Group A's U = -4t row; Group D is "
        "regenerated from (U=95t, Dt=4.5t, n=1.885).\n"
        "E7  U = -75t.  The stage-1 root that the earlier solver reported (Tc = 12.7t) was an artefact of taking the "
        "absolute value of the two gap components separately: the correctly signed root has the same chemical potential "
        "(about 35t), far outside the band.  The corrected solver finds no superconducting root with the chemical potential "
        "inside the band at any V_SO; the set is reported as having no superconducting state in the physical window and is "
        "not plotted as a curve.\n"
        "E8  Gap shape.  Stage 2 uses the rigid master shape phi(T) = (1-(T/Tc)^3.434)^0.553 for all V_SO.  The effect of a "
        "per-V_SO refit (1.3-3.4 % in the earlier measurement) was NOT re-measured in this correction; that statement is "
        "carried as an open item and is not a corrected result.\n"
        "E9  Citability.  Outputs of the U = -75t set are excluded from the citable presentation figures.\n"
        "E10 Grid convergence (F-02).  The earlier endpoint k-grid was not converged (its Tc drifted by 0.012 between "
        "40 and 160 points); the midpoint rule is converged to 3e-6 over the same range.  The master-set Tc changes by "
        f"{fmt(100*(M[0]['Tc']-0.28605)/0.28605,2)} % (V_SO=0) and by {fmt(100*(M[1]['Tc']-0.38084)/0.38084,2)} % (V_SO=0.4082).\n"
        "E11 DOS smoothing (F-03).  The Gaussian smoothing of the DOS is reduced from 0.02 t to 0.0025 t with a 48000-bin "
        "histogram.  The effect on W/W_N is -2.5 % to -3.0 % at the master state.  The smoothing limit is NOT reached: the "
        "value still changes by about 1 % per halving of the smoothing width; this residual is carried as a systematic of "
        "about 1 %, not a floor.\n"
        "E12 Chemical-potential window (decision).  A superconducting root is accepted only if its chemical potential lies "
        "inside the band, i.e. a Fermi surface exists.  This removes the superconducting state of the Group C (Dt = 4.5) "
        "set and of the low-U region (U = 15-30) of the companion families.  The policy is a physics decision for the "
        "authors; the numbers of the alternative are in the companion repository (sensitivity record).\n"
        "E13 Coherence variant (d).  Variant (d), 1/2 [1 - D D'/(E E')], is not a Bogoliubov matrix element and vanishes at the "
        "gap edge.  The earlier description of it as the textbook Case-I factor is withdrawn.  Its (d)/(a) ratio at the master "
        "state (T/Tc=0.9, V_SO=0) is 0.83.\n")

    # ---- Table 1 (groups)
    rows1 = []
    for r in n["groups"]:
        if r["nosc"]:
            rows1.append(f"   {r['group']}  {r['label']:<11s}  {fmt(r['Tc']):>8s}   {'-':>8s}    {'-':>8s}   {'no SC state (E12)':>18s}")
        else:
            rows1.append(f"   {r['group']}  {r['label']:<11s}  {fmt(r['Tc']):>8s}   {fmt(r['Deff']):>8s}    {fmt(r['min']):>8s}   {fmt(r['lim']):>8s}")
    t["table1"] = ("Table 1 - groups A-D, corrected (variant (a), scattering channel, hbar*omega=0.7t):\n"
                   "   set            Tc/t      Deff(0)/t  min W/W_N  lim T->Tc-\n" + "\n".join(rows1) + "\n")

    # ---- Table 2 (master)
    rows2 = []
    for r in M:
        rows2.append(f"   {r['V_SO']:<7.4f}  {fmt(r['Tc']):>7s}  {fmt(r['Deff']):>8s}   {fmt(r['lim']):>8s}   {fmt(r['R07']):>8s}")
    t["table2"] = ("Table 2 - master set (U=95t, Dt=4.5t, n=1.885), V_SO scan:\n"
                   "   V_SO     Tc/t   Deff(0)/t  lim W/W_N  W/W_N(T/Tc=0.7)\n" + "\n".join(rows2) + "\n")

    # ---- Table 3 (V_SO scans, T/Tc = 0.9)
    s95 = S.get("95", {}); s125 = S.get("125", {}); sm75 = S.get("-75", {})
    t["table3"] = ("Table 3 - V_SO scans at fixed T/Tc = 0.9 (range over the V_SO values that have a superconducting state):\n"
                   f"   U=95t:   {fmt(s95.get('lo'),3)} -> {fmt(s95.get('hi'),3)}   ({s95.get('n_sc')}/{s95.get('n')} V_SO values superconducting)\n"
                   f"   U=125t:  {fmt(s125.get('lo'),3)} -> {fmt(s125.get('hi'),3)}   ({s125.get('n_sc')}/{s125.get('n')} V_SO values superconducting)\n"
                   f"   U=-75t:  no superconducting state in the physical window ({sm75.get('n_sc', 0)}/{sm75.get('n', 0)} values; E7)\n")

    t["statements"] = (
        "Statements.  (i) In the scattering channel W/W_N < 1 for every set and every 0 < T < Tc.  (ii) The Tc- limit is "
        "R(hbar*omega, Tc) < 1 (E3), not 1.  (iii) Inside the paper's window the pair-breaking channel is identically zero "
        "(E4); pair-channel behaviour above threshold is an outlook only and carries no quantitative claim here.  "
        "(iv) The V_SO dependence of the limit is non-monotone at the largest V_SO for the master set only within the "
        "error budget of the DOS smoothing (E11).  (v) The low-U superconducting states of the companion families depend on "
        "the window policy (E12); no statement about them is made without that decision.\n")

    t["validation"] = (
        "3  Validation of the recomputation (regenerated)\n"
        "T0  Quadrature: the ell-substitution equals the explicit E-integral with both Jacobians; relative difference 3.9e-4 "
        "at the coarsest omega (regenerated).\n"
        "T1  Delta -> 0 closed form: W/W_N -> R(hbar*omega,T) = (2/hbar*omega) int_0^inf [f(E)-f(E+hbar*omega)] dE, "
        "verified to 3.1e-6.\n"
        "T2  omega -> 0 limit: the golden-rule ratio plateaus at 2 f(Delta) (Levy 1963, Eq. 8), verified to 2.0e-6; the "
        "flat-DOS test reproduces Levy to better than 1e-4 (test suite PHYS).\n"
        "T3  Lattice cross-check.  The brute-force lattice sum and the energy-space pipeline now share the same closed-form "
        "kernel (F-01); T3 therefore tests the integration and the DOS representation only.  The agreement for variant (a) is "
        "0.2-0.3 % (regenerated), down from about 1 % with the earlier DOS smoothing.  The kernel itself is tested by an "
        "independent derivation (E1, E13).\n"
        "K   Kernel identity: numerical Bogoliubov eigenvectors reproduce (uu' - vv')^2 = C_I to 7e-16; (uu' + vv')^2 = C_Case-II "
        "to 7e-16; completeness (uu'-vv')^2 + (uv'+vu')^2 = 1 to 9e-16; variant (d) matches no Bogoliubov element "
        "(minimum distance 0.5).\n"
        "F   Fix-verification suite: one executable test per audit finding; 18 of 18 pass on the regenerated outputs.\n")

    t["reproduction"] = (
        "4  Reproduction\n"
        "All numbers in this addendum are generated from the regenerated results of the companion repository.  "
        "Stage 1 (Tc, gap, mu) and stage 2 (W/W_N, families, peaks, SOC) are run by project/run_full_regeneration.sh "
        "(strict, ordered).  Every regenerated JSON is hashed in results/REGEN_MANIFEST.json.  Paper figures: "
        "project/paper_correction/make_replaced_figs.py (physics from project/ultrasonic/physics_ua.py).  This PDF: "
        "project/paper_correction/build_corrected_pdf.py.  The earlier references to sibling working directories are "
        "withdrawn: all inputs are inside the repository.\n")

    t["reading"] = (
        "5  Reading guide to the replaced figures\n"
        "Pages 11-14: six-panel composites in the original layout; the two W/W_N panels are the corrected ones.  "
        "The auxiliary panels (mu(T), Delta(T), universal gap curve) are recomputed from the corrected stage-1 solution.  "
        "A set without a superconducting state in the physical window is listed in the note of its page and is not plotted "
        "as a curve.  Page 15: W/W_N vs V_SO at fixed T/Tc.  Addendum figures: master-set detail, branch decomposition and "
        "limits (pp. 27-30).\n")
    return t


if __name__ == "__main__":
    for k, v in build().items():
        print("---", k, "---")
        print(v[:300])
