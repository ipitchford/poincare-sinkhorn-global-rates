#!/usr/bin/env python3
"""Floating-point checks for the numerical statements added in the revision.

Checked here (all are sanity checks of displayed formulas, not proofs):
  A. the refined expansion of the fixed-cusp deficit, eq. (44), and its
     constants a, b (closed forms against quadrature, 30 digits);
  B. the power-family constant c_a = 2^(1-a) Gamma(a+1) eta(a) (Remark 7.1);
  C. numerical estimates of the true Poincare constants of the uniform
     source, the fixed cusp, the Proposition 7 family and two power
     sources of Remark 7.1 (Sturm-Liouville shooting; Remark 7.2), compared with certified bounds;
  D. the factor between the Theorem 2 guarantee and the exact local deficit
     of the Proposition 7 family;
  E. an illustration of Theorem 3 / Corollary 4: Laplace source on R,
     quadratic cost, three target atoms (one-step ratios and the stopping
     certificate (33) of Remark 5.1);
  F. the local ratio (E(S psi)-E(psi)) / (lambda KL) -> 2 - D_lambda for the
     uniform two-atom example (Remark 5.2).
Every check appends to a failure list; the programme exits with status 1 if
any check fails (no assert statements, so it also fails under python -O).
Requires NumPy, SciPy and mpmath.
"""
import argparse
import json
import math
import platform
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

FAIL = []
LOG = []


def check(cond, label, value=None):
    if not cond:
        FAIL.append([label, value])


def say(s):
    LOG.append(s)
    print(s)


# ---------------------------------------------------------------- A
def part_a(report):
    mp.mp.dps = 30
    l2 = mp.log(2)
    sech2 = lambda u: mp.sech(u) ** 2
    brk = [0, 0.5, 1, 2, 4, 8, 16, 32, 64, mp.inf]
    a_quad = mp.quad(lambda u: u * mp.log(u) * sech2(u), brk)
    b_quad = mp.quad(lambda u: u * mp.log(u) ** 2 * sech2(u), brk)
    l_quad = mp.quad(lambda u: u * sech2(u), brk)
    g1 = mp.stieltjes(1)
    a_cf = l2 * (1 - mp.mpf(3) / 2 * l2)
    b_cf = l2 * (mp.pi ** 2 / 6 - mp.euler ** 2 - 2 * g1 - 3 * l2 + mp.mpf(7) / 3 * l2 ** 2)
    say(f"A. a = {mp.nstr(a_quad, 20)} (closed form {mp.nstr(a_cf, 20)})")
    say(f"   b = {mp.nstr(b_quad, 20)} (closed form {mp.nstr(b_cf, 20)})")
    check(abs(a_quad - a_cf) < mp.mpf("1e-25"), "a closed form", float(a_quad - a_cf))
    check(abs(b_quad - b_cf) < mp.mpf("1e-25"), "b closed form", float(b_quad - b_cf))
    check(abs(l_quad - l2) < mp.mpf("1e-25"), "int u sech^2 = log 2", float(l_quad - l2))
    rows = []
    for lam in ["0.1", "0.05", "0.01", "0.001", "1e-4", "1e-6"]:
        lam = mp.mpf(lam)
        L = mp.log(mp.e / lam)
        up = 1 / lam
        pts = [p for p in [0, 0.5, 1, 2, 4, 8, 16, 32, 64] if p < up] + [up]
        # D_lambda = int_{-1}^{1} sech^2(x/lambda) w(x) dx with x = lambda u (exact substitution)
        d_quad = mp.mpf(4) / 5 * lam ** 2 * mp.quad(
            lambda u: u * mp.log(mp.e / (lam * u)) ** 2 * sech2(u), pts)
        d_ref = mp.mpf(4) / 5 * lam ** 2 * (l2 * L ** 2 - 2 * a_cf * L + b_cf)
        tail = 69 / lam * mp.e ** (-2 / lam)
        lead = mp.mpf(4) / 5 * l2 * lam ** 2 * mp.log(1 / lam) ** 2
        rows.append({"lambda": float(lam), "D_quadrature": mp.nstr(d_quad, 15),
                     "D_refined": mp.nstr(d_ref, 15), "abs_diff": float(abs(d_quad - d_ref)),
                     "tail_bound": float(tail), "ratio_to_leading_term": float(d_quad / lead)})
        check(abs(d_quad - d_ref) <= tail + mp.mpf("1e-25") * d_quad, "refined expansion", float(lam))
        say(f"   lambda={mp.nstr(lam, 3):>7}: D={mp.nstr(d_quad, 13)}, refined={mp.nstr(d_ref, 13)}, "
            f"|diff|={float(abs(d_quad - d_ref)):.2e} <= {float(tail):.2e}, D/leading={float(d_quad / lead):.4f}")
    for lam, val in [("0.1", "0.0647035371183"), ("0.01", "0.00179454249623"), ("0.001", "0.0000353004893444")]:
        r = [x for x in rows if abs(x["lambda"] - float(lam)) < 1e-15][0]
        check(abs(float(r["D_quadrature"]) / float(val) - 1) < 1e-11, "paper D values", lam)
    report["A_refined_expansion"] = {"a": mp.nstr(a_quad, 25), "b": mp.nstr(b_quad, 25), "rows": rows}


# ---------------------------------------------------------------- B
def part_b(report):
    mp.mp.dps = 30
    out = []
    for a in ["0", "0.25", "0.5", "0.9", "1"]:
        a = mp.mpf(a)
        quad = mp.quad(lambda u: u ** a * mp.sech(u) ** 2, [0, 0.5, 1, 2, 4, 8, 16, 32, 64, mp.inf])
        eta = mp.altzeta(a)
        cf = 2 ** (1 - a) * mp.gamma(a + 1) * eta
        out.append({"a": float(a), "c_a": mp.nstr(quad, 20), "closed_form": mp.nstr(cf, 20),
                    "deficit_coefficient_2(a+1)c_a": mp.nstr(2 * (a + 1) * quad, 12)})
        check(abs(quad - cf) < mp.mpf("1e-24"), "power-family constant", float(a))
        say(f"B. a={float(a):.2f}: c_a={mp.nstr(quad, 16)} (closed form {mp.nstr(cf, 16)}), "
            f"local deficit ~ {mp.nstr(2 * (a + 1) * quad, 8)} lambda^(1+a)")
    report["B_power_constants"] = out


# ---------------------------------------------------------------- C
def first_root(mismatch, k_lo, k_hi):
    ks = np.geomspace(k_lo, k_hi, 120)
    vals = [mismatch(k) for k in ks]
    for i in range(len(ks) - 1):
        if vals[i] == 0:
            return ks[i]
        if vals[i] * vals[i + 1] < 0:
            return brentq(mismatch, ks[i], ks[i + 1], xtol=1e-14, rtol=1e-13)
    raise RuntimeError("no eigenvalue bracket")


def sl_problem(weight_scale, S, parity, rtol, breaks=()):
    """-(f')' = kappa W(s) f on (0,S) in natural-scale form y1'=y2, y2'=-kappa W y1.

    weight_scale(s) returns (p(s), W(s)) for the system y1' = y2/p, y2' = -kappa W y1.
    parity 'odd': y1(0)=0, y2(0)=1; 'even': y1(0)=1, y2(0)=0. Mismatch y2(S)."""
    knots = [0.0] + [b for b in breaks if 0 < b < S] + [S]

    def mismatch(k):
        y = np.array([0.0, 1.0]) if parity == "odd" else np.array([1.0, 0.0])
        for s0, s1 in zip(knots[:-1], knots[1:]):
            def rhs(s, v):
                p, W = weight_scale(s)
                return [v[1] / p, -k * W * v[0]]
            sol = solve_ivp(rhs, (s0, s1), y, method="DOP853", rtol=rtol, atol=1e-14 * rtol / 1e-12)
            y = sol.y[:, -1]
        return y[1]
    return mismatch


def poincare_estimate(weight_scale, S, breaks=(), k_lo=1e-3, k_hi=1e4):
    res = {}
    for rtol in (1e-10, 1e-12):
        ko = first_root(sl_problem(weight_scale, S, "odd", rtol, breaks), k_lo, k_hi)
        ke = first_root(sl_problem(weight_scale, S, "even", rtol, breaks), k_lo, k_hi)
        res[rtol] = (1 / ko, 1 / ke)
    co, ce = res[1e-12]
    drift = max(abs(res[1e-10][0] - co), abs(res[1e-10][1] - ce))
    return max(co, ce), co, ce, drift


def part_c(report):
    out = {}
    # uniform on (-1,1): exact C_P = 4/pi^2
    cp, co, ce, dr = poincare_estimate(lambda s: (1.0, 1.0), 1.0)
    out["uniform"] = {"C_num": cp, "odd": co, "even": ce, "exact": 4 / math.pi ** 2}
    check(abs(cp - 4 / math.pi ** 2) < 1e-9, "uniform sanity", cp)
    say(f"C. uniform: C_P(num)={cp:.10f} (exact 4/pi^2={4 / math.pi ** 2:.10f})")
    # fixed cusp in natural scale s = R(x) = 5/(2 log(e/x)), W(s) = (25/4) e^(2-5/s)/s^4
    def cusp(s):
        if s <= 0.0:
            return 1.0, 0.0
        return 1.0, 6.25 * math.exp(2.0 - 5.0 / s) / s ** 4
    cp, co, ce, dr = poincare_estimate(cusp, 2.5)
    out["cusp"] = {"C_num": cp, "odd": co, "even": ce, "certified": 0.75, "solver_drift": dr}
    check(cp <= 0.75 and dr < 1e-7, "cusp below certificate", cp)
    say(f"   fixed cusp (31): C_P(num)={cp:.6f} (odd {co:.6f}, even {ce:.6f}); certified 3/4")
    # power family w = ((a+1)/2)|x|^a, natural scale
    for a in (0.5, 0.9):
        S = 2.0 / ((a + 1) * (1 - a))
        def pw(s, a=a):
            x = ((a + 1) * (1 - a) * max(s, 0.0) / 2.0) ** (1.0 / (1 - a))
            return 1.0, ((a + 1) / 2.0) ** 2 * x ** (2 * a)
        cp, co, ce, dr = poincare_estimate(pw, S)
        bound = 1 / (2 * (1 - a))
        out[f"power_a={a}"] = {"C_num": cp, "odd": co, "even": ce, "certified": bound, "solver_drift": dr}
        check(cp <= bound and dr < 1e-7, "power family below certificate", [a, cp])
        say(f"   power a={a}: C_P(num)={cp:.6f}; certified 1/(2(1-a))={bound:.4f}")
    # Proposition 7 family, x-coordinates, p = W = rho_lambda (normalisation irrelevant)
    fam = []
    for lam in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005):
        xs = lam * math.acosh(1 / lam)
        rho = lambda x, lam=lam: min(lam * math.cosh(x / lam), 1.0)
        cp, co, ce, dr = poincare_estimate(lambda s, rho=rho: (rho(s), rho(s)), 1.0, breaks=(xs,))
        I = 1 - xs + lam * math.sqrt(1 - lam ** 2)
        J = math.acos(lam) + 1 - xs
        D = (lam ** 2 * math.acos(lam) + lam * (math.tanh(1 / lam) - math.sqrt(1 - lam ** 2))) / I
        delta = 1 - (1 - D) ** 2
        g_cert = lam ** 2 / (1 + math.pi / 2 + lam ** 2)
        g_num = lam ** 2 / (cp + lam ** 2)
        row = {"lambda": lam, "C_num": cp, "odd": co, "even": ce, "certificate_I_times_J": I * J,
               "local_deficit": delta, "deficit_over_pi_lambda2": delta / (math.pi * lam ** 2),
               "factor_vs_certified_guarantee": delta / g_cert, "factor_vs_numerical_guarantee": delta / g_num,
               "solver_drift": dr}
        fam.append(row)
        check(cp <= I * J + 1e-9 and I * J <= 1 + math.pi / 2 and dr < 1e-7, "Prop 7 family certificate", lam)
        check(delta <= g_num * 20 and delta >= g_num, "guarantee consistent with exact deficit", lam)
        say(f"   Prop 7, lambda={lam}: C_P(num)={cp:.5f}, I*J={I * J:.5f}, deficit/(pi lam^2)={delta / (math.pi * lam ** 2):.4f}, "
            f"exact/guaranteed deficit = {delta / g_cert:.3f} (certified C_P), {delta / g_num:.3f} (numerical C_P)")
    out["prop7_family"] = fam
    report["C_poincare_numerical"] = out


# ---------------------------------------------------------------- E
def part_e(report):
    # Laplace source (1/2)e^{-|x|} on R, quadratic cost; reduced cost -xy + y^2/2
    gx, gw = np.polynomial.legendre.leggauss(20)
    T, h = 60.0, 0.25
    edges = np.arange(0.0, T + h / 2, h)
    xs, ws = [], []
    for a0, a1 in zip(edges[:-1], edges[1:]):
        xs.append(0.5 * (a1 - a0) * gx + 0.5 * (a1 + a0))
        ws.append(0.5 * (a1 - a0) * gw)
    xp, wp = np.concatenate(xs), np.concatenate(ws)
    x = np.concatenate([-xp[::-1], xp])
    w = np.concatenate([wp[::-1], wp]) * 0.5 * np.exp(-np.abs(np.concatenate([-xp[::-1], xp])))
    mass = float(w.sum())
    check(abs(mass - 1) < 1e-13, "Laplace quadrature mass", mass)
    y = np.array([-1.0, 0.3, 1.0])
    nu = np.array([0.2, 0.5, 0.3])
    R, CP = 1.0, 4.0
    rows = []

    def quantities(psi, lam):
        logk = np.log(nu)[None, :] + (psi[None, :] + x[:, None] * y[None, :] - 0.5 * y[None, :] ** 2) / lam
        m = logk.max(axis=1, keepdims=True)
        logZ = (m + np.log(np.exp(logk - m).sum(axis=1, keepdims=True)))[:, 0]
        E = float(nu @ psi - lam * (w @ logZ))
        r = (w[:, None] * np.exp(logk - logZ[:, None])).sum(axis=0)
        return E, r

    for lam in (1.0, 0.5, 0.3):
        psi = np.zeros(3)
        Es = []
        for t in range(60000):
            E, r = quantities(psi, lam)
            Es.append(E)
            new = psi - lam * np.log(r / nu)
            new -= new.mean()
            if np.max(np.abs(new - psi)) < 1e-15:
                break
            psi = new
        Estar = max(Es)
        q = CP * R ** 2 / (CP * R ** 2 + lam ** 2)
        A_over_lam = 1 + CP * R ** 2 / lam ** 2
        for start in ("zero", "skewed"):
            psi = np.zeros(3) if start == "zero" else np.array([3.0, -2.0, 0.5])
            worst_ratio, worst_cert = 0.0, 0.0
            E, r = quantities(psi, lam)
            for t in range(400):
                gap = Estar - E
                if gap < 1e-9:
                    break
                psi2 = psi - lam * np.log(r / nu)
                E2, r2 = quantities(psi2, lam)
                gap2 = Estar - E2
                worst_ratio = max(worst_ratio, gap2 / gap)
                worst_cert = max(worst_cert, gap / (A_over_lam * (E2 - E)))
                psi, E, r = psi2, E2, r2
            rows.append({"lambda": lam, "start": start, "bound_q": q, "worst_one_step_ratio": worst_ratio,
                         "worst_certificate_ratio": worst_cert, "iterations_to_converge": len(Es)})
            check(worst_ratio <= q + 1e-9, "Theorem 3 one-step bound (Laplace)", [lam, start, worst_ratio])
            check(worst_cert <= 1 + 1e-9, "stopping certificate (Laplace)", [lam, start, worst_cert])
            say(f"E. Laplace/quadratic, lambda={lam}, start={start}: worst gap ratio {worst_ratio:.5f} <= bound {q:.5f}; "
                f"worst gap/(A/lam * increase) = {worst_cert:.4f} <= 1")
    report["E_laplace_quadratic"] = {"quadrature_mass": mass, "rows": rows}


# ---------------------------------------------------------------- F
def part_f(report):
    mp.mp.dps = 40
    out = []
    for lam in (mp.mpf("0.5"), mp.mpf("0.1")):
        D = lam * mp.tanh(1 / lam)
        th = mp.mpf("1e-6")
        Elc = lambda t: lam / 2 * mp.quad(lambda x: mp.log(mp.cosh((x + t) / lam)), [-1, -t, 0, 1])
        m = lam / 2 * (mp.log(mp.cosh((1 + th) / lam)) - mp.log(mp.cosh((th - 1) / lam)))
        th2 = th - lam * mp.atanh(m)
        increase = Elc(th) - Elc(th2)          # E(theta) = -lambda E log cosh((X+theta)/lambda)
        kl = -mp.log(1 - m ** 2) / 2
        ratio = increase / (lam * kl)
        out.append({"lambda": float(lam), "ratio": float(ratio), "two_minus_D": float(2 - D)})
        check(abs(ratio - (2 - D)) < 1e-4, "local ascent ratio 2-D", float(lam))
        say(f"F. uniform two-atom, lambda={float(lam)}: (E(S psi)-E(psi))/(lambda KL) = {float(ratio):.6f}, 2-D = {float(2 - D):.6f}")
    report["F_local_ascent_ratio"] = out


def run(outdir):
    t0 = time.time()
    outdir.mkdir(parents=True, exist_ok=True)
    report = {"python": sys.version.split()[0], "platform": platform.platform(), "numpy": np.__version__,
              "scipy": scipy.__version__, "mpmath": mp.__version__}
    for part in (part_a, part_b, part_c, part_e, part_f):
        try:
            part(report)
        except Exception as exc:  # any numerical breakdown is a failure, not a pass
            FAIL.append([part.__name__, repr(exc)])
    report["failures"] = FAIL
    report["passed"] = not FAIL
    report["elapsed_seconds"] = time.time() - t0
    say(f"Failures: {len(FAIL)}. Elapsed: {report['elapsed_seconds']:.1f}s.")
    for f in FAIL:
        say(f"FAILED: {f}")
    (outdir / "revision_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    (outdir / "revision_checks.log").write_text("\n".join(LOG) + "\n")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "checks")
    raise SystemExit(run(ap.parse_args().output_dir))
